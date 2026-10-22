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
 * the reference routine on its part of the output.  R*_gpu return true when
 * they have computed the result on the GPU, and false, leaving every
 * operand untouched, when the call is left to the CPU:
 *
 *   - invalid arguments and quick returns;
 *   - MPLAPACK_MPFR_CUDA=0 or MPLAPACK_MPFR_CUDA_L3=0; Level 1 and 2 unless
 *     MPLAPACK_MPFR_CUDA_L12=1; rounding mode other than MPFR_RNDN; fewer
 *     multiply-adds than MPLAPACK_MPFR_CUDA_MIN_MNK;
 *   - a default precision (mpfr_class, both components of mpc_class) other
 *     than one of 256/512/768/1024/2048 bits (cu_freal; there is no runtime
 *     precision engine for the complex routines);
 *   - an element of an operand (the whole rectangle the routine may read)
 *     or a scalar of another precision, not finite, or with an exponent
 *     outside the bound of mpfr_pack_cuda.h, in which no intermediate
 *     overflows or underflows; a result outside the exponent range (Ctrsm:
 *     outside the bound, as its results are factors of later operations);
 *   - Ctrsm with diag = "N" (complex division is not computed on the GPU).
 *
 * The results equal the reference ones bit for bit, except that zeros the
 * GPU computes are +0.  Elements the routine does not change keep their
 * value and sign.
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

// The element-type side of Cblas_bridge_body.h for precision PB: operands
// are packed into device elements with the layout (leading dimension) of
// the host ones; the written ones are unpacked after the run.
template <int PB> struct pack_ctx {
    typedef mreal<PB> DR;
    typedef mcomplex<PB> DC;
    template <class E> struct held {
        std::vector<E> buf, orig;
        void *host;
        long ld, rows, cols;
        bool download;
    };
    std::deque<held<DC> > hc;
    std::deque<held<DR> > hr;
    bool ok;
    const mpfr_exp_t bound;

    pack_ctx() : ok(true), bound(input_exponent_bound()) {}

    bool eligible(bool level3, double work) const
    {
        return level_enabled(level3) && work >= (double)min_mnk() && mpfr_get_default_rounding_mode() == MPFR_RNDN && mpfrxx::default_mpc_real_precision_bits() == PB && mpfrxx::default_mpc_imag_precision_bits() == PB;
    }

    bool pk(const mpfr_class &s, DR &d) const { return pack<PB>(s.get_mpfr_t(), d.v, bound); }
    bool pk(const mpc_class &s, DC &d) const { return pack<PB>(mpc_realref(s.mpc_data()), d.re_.v, bound) && pack<PB>(mpc_imagref(s.mpc_data()), d.im_.v, bound); }
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
        gpu_operand<E> o = {h.buf.data(), ld, 0, count, download};
        return o;
    }

    // The result r of an element may be written back.
    bool fits(const real_t<PB> &r, bool factor) const
    {
        if (r.is_zero())
            return true;
        if (factor)
            return r.exp >= -bound && r.exp <= bound;
        return r.exp >= mpfr_get_emin() && r.exp <= mpfr_get_emax();
    }
    bool fits(const DR &r, bool factor) const { return fits(r.v, factor); }
    bool fits(const DC &r, bool factor) const { return fits(r.re_.v, factor) && fits(r.im_.v, factor); }
    static bool same_e(const DR &a, const DR &b) { return same<PB>(a.v, b.v); }
    static bool same_e(const DC &a, const DC &b) { return same<PB>(a.re_.v, b.re_.v) && same<PB>(a.im_.v, b.im_.v); }
    static void put(const DR &r, mpfr_class &x) { unpack<PB>(r.v, x.get_mpfr_t()); }
    static void put(const DC &r, mpc_class &x)
    {
        unpack<PB>(r.re_.v, mpc_realref(x.mpc_data()));
        unpack<PB>(r.im_.v, mpc_imagref(x.mpc_data()));
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
        const int rc = mplapack_blas_dev::cgpu_run(j, ops, nops);
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
    bool run(const mplapack_blas_dev::job_Ctrsm<mpfr_fixed<PB> > &j, const gpu_operand<DC> *ops, int nops)
    {
        if (mplapack_blas_dev::lsame(j.diag, "N"))
            return false;
        return run_checked(j, ops, nops, true);
    }
};

} // namespace

#define MPLAPACK_MPFR_CBLAS_CASE(PB, F, ...) \
    case PB: {                               \
        pack_ctx<PB> ctx;                    \
        return mplapack_cblas_gpu::F<mpfr_fixed<PB> >(ctx, __VA_ARGS__); \
    }
#define MPLAPACK_MPFR_CBLAS_RUN(F, ...)                         \
    ([&]() -> bool {                                            \
        switch (mpfrxx::default_precision_bits()) {             \
            MPLAPACK_MPFR_CBLAS_CASE(256, F, __VA_ARGS__)       \
            MPLAPACK_MPFR_CBLAS_CASE(512, F, __VA_ARGS__)       \
            MPLAPACK_MPFR_CBLAS_CASE(768, F, __VA_ARGS__)       \
            MPLAPACK_MPFR_CBLAS_CASE(1024, F, __VA_ARGS__)      \
            MPLAPACK_MPFR_CBLAS_CASE(2048, F, __VA_ARGS__)      \
        default:                                                \
            return false;                                       \
        }                                                       \
    }())
MPLAPACK_CBLAS_GPU_FUNCTIONS(MPLAPACK_MPFR_CBLAS_RUN)
