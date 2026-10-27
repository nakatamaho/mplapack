/*
 * Copyright (c) 2026
 *      Nakata, Maho
 *      All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

/*
 * Host side of the MPFR complex Level 1/2/3 BLAS (and Rgemmtr) on CUDA:
 * the routines of ../../common/blas_cuda/Cblas_bridge_body.h with the device
 * types of Cblas_mpfr_dev_types.h, where every GPU thread runs the text of
 * the reference routine on its part of the output.  The entry points
 * <R>_mpfr_cuda (include/mpblas_mpfr_cuda.h) return true when they have
 * computed the result on the GPU, and false, leaving every
 * operand untouched, when the call is left to the CPU:
 *
 *   - invalid arguments and quick returns;
 *   - MPLAPACK_MPFR_CUDA=0 or MPLAPACK_MPFR_CUDA_L3=0; Level 1 and 2 unless
 *     MPLAPACK_MPFR_CUDA_L12=1; rounding mode other than MPFR_RNDN; fewer
 *     multiply-adds than MPLAPACK_MPFR_CUDA_MIN_MNK;
 *   - a default precision p (mpfr_class, both components of mpc_class) above
 *     4096 bits (16384 with MPLAPACK_MPFR_CUDA_RT_LARGE).  p = 256/512/768/1024/2048 uses cu_freal (zeros are +0),
 *     other precisions (or MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1) the MPFR/MPC
 *     functions themselves (libmpc_cuda; MPLAPACK_MPFR_CUDA_RUNTIME=0 turns
 *     them off), whose results are the reference ones signed zeros included;
 *   - an element of an operand (the whole rectangle the routine may read)
 *     or a scalar of another precision, not finite, or with an exponent
 *     outside the bound of mpfr_pack_cuda.h, in which no intermediate
 *     overflows or underflows; a result outside the exponent range (Ctrsm:
 *     outside the bound, as its results are factors of later operations);
 *   - Ctrsm with diag = "N" and a zero diagonal element of A.
 *
 * The results equal the reference ones bit for bit (cu_freal: up to the
 * sign of zeros).  Elements the routine does not change keep their value and
 * sign.
 */

#include <mpblas.h>
#include <cstdio>
#include <deque>
#include <type_traits>
#include <vector>
#include "mpfr_pack_cuda.h"
#include "Cblas_mpfr_dev_types.h"
#include "../../common/blas_cuda/Cblas_bridge_body.h"

namespace {

using namespace mplapack_mpfr_cuda_host;
using mplapack_blas::gpu_operand;
using mplapack_blas_dev::mcomplex;
using mplapack_blas_dev::mpfr_fixed;
using mplapack_blas_dev::mreal;

bool level_enabled(bool level3)
{
    static const bool l3 = env_long("MPLAPACK_MPFR_CUDA_L3", 1) != 0;
    static const bool l12 = env_long("MPLAPACK_MPFR_CUDA_L12", 0) != 0;
    return cuda_enabled() && (level3 ? l3 : l12);
}

bool verbose() { return env_long("MPLAPACK_MPFR_CUDA_VERBOSE", 0) > 0; }

template <int PB> bool same(const real_t<PB> &a, const real_t<PB> &b)
{
    if (a.is_zero() || b.is_zero())
        return a.is_zero() && b.is_zero();
    if (a.sign != b.sign || a.exp != b.exp)
        return false;
    for (int i = 0; i < real_t<PB>::N; i++)
        if (a.m[i] != b.m[i])
            return false;
    return true;
}

// The two engines: cu_freal<PB> (fixed precision PB; zeros are +0) and the
// MPFR/MPC functions at the runtime precision p (mpfr_rt<NL>, p <= 64 NL).
template <int PB> struct fixed_engine {
    typedef mpfr_fixed<PB> Tag;
    typedef mreal<PB> DR;
    typedef mcomplex<PB> DC;
    long p;
    explicit fixed_engine(long) : p(PB) {}
    // device stack: 16 KB per 512 bits (twice what 512 bits need on an A100)
    static size_t stack_bytes() { return (size_t)16384 * ((PB + 511) / 512); }
    bool pack1(mpfr_srcptr x, DR &d, mpfr_exp_t bound) const { return pack<PB>(x, d.v, bound); }
    void unpack1(const DR &d, mpfr_ptr x) const { unpack<PB>(d.v, x); }
    bool same1(const DR &a, const DR &b) const { return same<PB>(a.v, b.v); }
    bool zero1(const DR &d) const { return d.v.is_zero(); }
    long exp1(const DR &d) const { return d.v.exp; }
    template <class J> int run(const J &j, const gpu_operand<typename J::E> *ops, int nops) const { return mplapack_blas_dev::cgpu_run(j, ops, nops); }
};

template <int NL> struct rt_engine {
    typedef mplapack_blas_dev::mpfr_rt<NL> Tag;
    typedef mplapack_blas_dev::rreal<NL> DR;
    typedef mplapack_blas_dev::rcomplex<NL> DC;
    long p;
    int nl;
    explicit rt_engine(long prec) : p(prec), nl((int)((prec + 63) / 64)) {}
    bool pack1(mpfr_srcptr x, DR &d, mpfr_exp_t bound) const
    {
        if (mpfr_get_prec(x) != p)
            return false;
        const int kind = mpfr_custom_get_kind(x);
        d.kind = kind;
        d.exp = 0;
        if (kind == MPFR_ZERO_KIND || kind == -MPFR_ZERO_KIND)
            return true;
        if (kind != MPFR_REGULAR_KIND && kind != -MPFR_REGULAR_KIND)
            return false;
        const mpfr_exp_t e = mpfr_get_exp(x);
        if (e > bound || e < -bound)
            return false;
        d.exp = (long)e;
        const mp_limb_t *s = static_cast<const mp_limb_t *>(mpfr_custom_get_significand(x));
        for (int t = 0; t < nl; t++)
            d.m[t] = s[t];
        return true;
    }
    void unpack1(const DR &d, mpfr_ptr x) const
    {
        mpfr_t v;
        mpfr_custom_init_set(v, d.kind, (mpfr_exp_t)d.exp, (mpfr_prec_t)p, const_cast<mp_limb_t *>(reinterpret_cast<const mp_limb_t *>(d.m)));
        mpfr_set(x, v, MPFR_RNDN); // exact: x has precision p
    }
    bool zero1(const DR &d) const { return d.kind == MPFR_ZERO_KIND || d.kind == -MPFR_ZERO_KIND; }
    bool same1(const DR &a, const DR &b) const
    {
        if (a.kind != b.kind)
            return false;
        if (zero1(a))
            return true;
        if (a.exp != b.exp)
            return false;
        for (int t = 0; t < nl; t++)
            if (a.m[t] != b.m[t])
                return false;
        return true;
    }
    long exp1(const DR &d) const { return d.exp; }
    template <class J> int run(const J &j, const gpu_operand<typename J::E> *ops, int nops) const { return mplapack_blas_dev::cgpu_run_rt(j, ops, nops, p); }
    // device stack of the cu_mpfr / cu_mpc functions, as the real runtime kernels
    static size_t stack_bytes() { return (size_t)65536; }
};

// The element-type side of Cblas_bridge_body.h: operands are packed into
// device elements with the layout (leading dimension) of the host ones; the
// written ones are unpacked after the run.
template <class Eng> struct pack_ctx {
    typedef typename Eng::DR DR;
    typedef typename Eng::DC DC;
    template <class E> struct held {
        std::vector<E> buf, orig;
        void *host;
        long ld, rows, cols;
        bool download;
    };
    const Eng eng;
    std::deque<held<DC> > hc;
    std::deque<held<DR> > hr;
    bool ok;
    const mpfr_exp_t bound;

    explicit pack_ctx(long prec) : eng(prec), ok(true), bound(input_exponent_bound()) {}

    bool eligible(bool level3, double work) const
    {
        return level_enabled(level3) && work >= (double)min_mnk() && mpfr_get_default_rounding_mode() == MPFR_RNDN && mpfrxx::default_mpc_real_precision_bits() == eng.p && mpfrxx::default_mpc_imag_precision_bits() == eng.p;
    }

    bool pk(const mpfr_class &s, DR &d) const { return eng.pack1(s.get_mpfr_t(), d, bound); }
    bool pk(const mpc_class &s, DC &d) const { return eng.pack1(mpc_realref(s.mpc_data()), d.re_, bound) && eng.pack1(mpc_imagref(s.mpc_data()), d.im_, bound); }
    void sc(DR &d, const mpfr_class &s) { ok = ok && pk(s, d); }
    void sc(DC &d, const mpc_class &s) { ok = ok && pk(s, d); }

    template <class E> std::deque<held<E> > &store()
    {
        if constexpr (std::is_same<E, DC>::value)
            return hc;
        else
            return hr;
    }

    // rows x cols elements of x with leading dimension (vectors: increment) ld > 0
    template <class E, class T> gpu_operand<E> op(T *x, long ld, long rows, long cols, bool download)
    {
        std::deque<held<E> > &s = store<E>();
        s.push_back(held<E>());
        held<E> &h = s.back();
        h.host = x, h.ld = ld, h.rows = rows, h.cols = cols, h.download = download;
        const long count = (cols - 1) * ld + rows;
        h.buf.resize((size_t)count);
        for (long j = 0; j < cols && ok; j++)
            for (long i = 0; i < rows && ok; i++)
                ok = pk(x[i + j * ld], h.buf[(size_t)(i + j * ld)]);
        if (download)
            h.orig = h.buf;
        gpu_operand<E> o = {h.buf.data(), ld, 0, count, download, rows, cols, mplapack_blas::REGION_ALL};
        return o;
    }

    // The result r of an element may be written back.
    bool fits1(const DR &r, bool factor) const
    {
        if (eng.zero1(r))
            return true;
        const long e = eng.exp1(r);
        if (factor)
            return e >= -bound && e <= bound;
        return e >= mpfr_get_emin() && e <= mpfr_get_emax();
    }
    bool fits(const DR &r, bool factor) const { return fits1(r, factor); }
    bool fits(const DC &r, bool factor) const { return fits1(r.re_, factor) && fits1(r.im_, factor); }
    bool same_e(const DR &a, const DR &b) const { return eng.same1(a, b); }
    bool same_e(const DC &a, const DC &b) const { return eng.same1(a.re_, b.re_) && eng.same1(a.im_, b.im_); }
    void put(const DR &r, mpfr_class &x) const { eng.unpack1(r, x.get_mpfr_t()); }
    void put(const DC &r, mpc_class &x) const
    {
        eng.unpack1(r.re_, mpc_realref(x.mpc_data()));
        eng.unpack1(r.im_, mpc_imagref(x.mpc_data()));
    }

    template <class E, class T> bool check(const std::deque<held<E> > &s, bool factor) const
    {
        for (size_t q = 0; q < s.size(); q++) {
            const held<E> &h = s[q];
            if (!h.download)
                continue;
            for (long j = 0; j < h.cols; j++)
                for (long i = 0; i < h.rows; i++) {
                    const size_t t = (size_t)(i + j * h.ld);
                    if (!same_e(h.buf[t], h.orig[t]) && !fits(h.buf[t], factor))
                        return false;
                }
        }
        return true;
    }

    template <class E, class T> void write(const std::deque<held<E> > &s) const
    {
        for (size_t q = 0; q < s.size(); q++) {
            const held<E> &h = s[q];
            if (!h.download)
                continue;
            T *x = static_cast<T *>(h.host);
            for (long j = 0; j < h.cols; j++)
                for (long i = 0; i < h.rows; i++) {
                    const size_t t = (size_t)(i + j * h.ld);
                    if (!same_e(h.buf[t], h.orig[t]))
                        put(h.buf[t], x[i + j * h.ld]);
                }
        }
    }

    template <class J> bool run_checked(const J &j, const gpu_operand<typename J::E> *ops, int nops, bool factor)
    {
        if (!ok)
            return false;
        int rc = mplapack_blas_dev::mpfr_cblas_stack(Eng::stack_bytes());
        if (rc == 0)
            rc = eng.run(j, ops, nops);
        if (rc != 0) {
            if (verbose())
                std::fprintf(stderr, "mplapack mpfr cuda (complex BLAS): GPU error %d, computing on the CPU\n", rc);
            return false;
        }
        if (!check<DC, mpc_class>(hc, factor) || !check<DR, mpfr_class>(hr, factor))
            return false;
        write<DC, mpc_class>(hc);
        write<DR, mpfr_class>(hr);
        return true;
    }
    template <class J> bool run(const J &j, const gpu_operand<typename J::E> *ops, int nops) { return run_checked(j, ops, nops, false); }
    bool run(const mplapack_blas_dev::job_Ctrsm<typename Eng::Tag> &j, const gpu_operand<DC> *ops, int nops)
    {
        if (mplapack_blas_dev::lsame(j.diag, "N") && !hc.empty()) {
            // the divisors a(k, k): a zero one gives Inf/NaN
            const held<DC> &A = hc[0];
            for (long k = 0; k < A.cols && k < A.rows; k++)
                if (eng.zero1(A.buf[(size_t)(k + k * A.ld)].re_) && eng.zero1(A.buf[(size_t)(k + k * A.ld)].im_))
                    return false;
        }
        return run_checked(j, ops, nops, true);
    }
};

} // namespace

// The engine for the default precision p: cu_freal for the precisions of
// MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS (unless MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1),
// otherwise the MPFR/MPC functions with the smallest limb count >= p/64
// (p <= 4096, with MPLAPACK_MPFR_CUDA_RT_LARGE p <= 16384; larger precisions
// stay on the CPU).
#define MPLAPACK_MPFR_CBLAS_FIXED(PB, F, ...)                                     \
    case PB: {                                                                    \
        pack_ctx<fixed_engine<PB> > ctx(PB);                                      \
        return mplapack_cblas_gpu::F<mpfr_fixed<PB> >(ctx, __VA_ARGS__);          \
    }
#define MPLAPACK_MPFR_CBLAS_RT(NL, F, ...)                                        \
    if (p <= 64 * NL) {                                                           \
        pack_ctx<rt_engine<NL> > ctx(p);                                          \
        return mplapack_cblas_gpu::F<mplapack_blas_dev::mpfr_rt<NL> >(ctx, __VA_ARGS__); \
    }
#ifdef MPLAPACK_MPFR_CUDA_RT_LARGE
#define MPLAPACK_MPFR_CBLAS_RT_LARGE(F, ...) MPLAPACK_MPFR_CBLAS_RT(128, F, __VA_ARGS__) MPLAPACK_MPFR_CBLAS_RT(256, F, __VA_ARGS__)
#else
#define MPLAPACK_MPFR_CBLAS_RT_LARGE(F, ...)
#endif
#define MPLAPACK_MPFR_CBLAS_RUN(F, ...)                                           \
    ([&]() -> bool {                                                              \
        std::lock_guard<std::mutex> lock(gpu_mutex());                            \
        const long p = (long)mpfrxx::default_precision_bits();                    \
        if (!force_runtime())                                                     \
            switch (p) {                                                          \
                MPLAPACK_MPFR_CBLAS_FIXED(256, F, __VA_ARGS__)                    \
                MPLAPACK_MPFR_CBLAS_FIXED(512, F, __VA_ARGS__)                    \
                MPLAPACK_MPFR_CBLAS_FIXED(768, F, __VA_ARGS__)                    \
                MPLAPACK_MPFR_CBLAS_FIXED(1024, F, __VA_ARGS__)                   \
                MPLAPACK_MPFR_CBLAS_FIXED(2048, F, __VA_ARGS__)                   \
            default:                                                              \
                break;                                                            \
            }                                                                     \
        if (!runtime_enabled())                                                   \
            return false;                                                         \
        MPLAPACK_MPFR_CBLAS_RT(2, F, __VA_ARGS__)                                 \
        MPLAPACK_MPFR_CBLAS_RT(4, F, __VA_ARGS__)                                 \
        MPLAPACK_MPFR_CBLAS_RT(8, F, __VA_ARGS__)                                 \
        MPLAPACK_MPFR_CBLAS_RT(16, F, __VA_ARGS__)                                \
        MPLAPACK_MPFR_CBLAS_RT(32, F, __VA_ARGS__)                                \
        MPLAPACK_MPFR_CBLAS_RT(64, F, __VA_ARGS__)                                \
        MPLAPACK_MPFR_CBLAS_RT_LARGE(F, __VA_ARGS__)                              \
        return false;                                                             \
    }())
MPLAPACK_CBLAS_GPU_FUNCTIONS(MPLAPACK_MPFR_CBLAS_RUN)
