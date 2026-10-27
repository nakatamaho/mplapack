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
 * Host side of the MPFR Level 1/2/3 BLAS on CUDA (kernels in
 * blas3_kernel_cuda.h, lanes of ../../common/Rblas_lanes_common.h).
 * The entry points R*_mpfr_cuda() (include/mpblas_mpfr_cuda.h) return true
 * when they have computed the result on the GPU, and false, leaving every
 * operand untouched, when the call is left to the CPU:
 *
 *   - invalid arguments, quick returns, alpha == 0 (Level 2/3) and zero
 *     increments (the reference reports or handles them);
 *   - MPLAPACK_MPFR_CUDA=0 or MPLAPACK_MPFR_CUDA_L3=0; Level 1 and 2 unless
 *     MPLAPACK_MPFR_CUDA_L12=1 (they are bound by the transfers; meant for
 *     debugging); rounding mode other than MPFR_RNDN; fewer multiply-adds
 *     than MPLAPACK_MPFR_CUDA_MIN_MNK;
 *   - an operand or scalar of another precision than the MPFR default
 *     precision p, or not finite, among the referenced elements;
 *   - all routines but Rtrsm: an input exponent outside the bound that keeps
 *     every intermediate inside the exponent range, or a result outside it;
 *   - Rtrsm (always the cu_mpfr engine: cu_freal has no correctly rounded
 *     division): an MPFR exponent range other than the device default
 *     [1-2^30, 2^30-1], so that overflow and underflow happen exactly as on
 *     the host.
 *
 * p = 256/512/768/1024/2048 uses the cu_freal engine (zeros are +0), other
 * precisions (or MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1) the cu_mpfr engine,
 * whose results equal the reference bit for bit, signed zeros included.
 */

#include <mpblas_mpfr_cuda.h>
#include "mpfr_pack_cuda.h"
#include "blas3_kernel_cuda.h"
#include "../../common/Rblas_calls_common.h"

namespace {

using namespace mplapack_mpfr_cuda_host;
using mplapack_blas::blas_call;
using mplapack_blas::blas_operand;

typedef blas_operand<mpfr_class> operand;

bool l3_enabled()
{
    static const bool enabled = env_long("MPLAPACK_MPFR_CUDA_L3", 1) != 0;
    return cuda_enabled() && enabled;
}

bool l12_enabled()
{
    static const bool enabled = env_long("MPLAPACK_MPFR_CUDA_L12", 0) != 0;
    return cuda_enabled() && enabled;
}

// The operands as the lanes see them: A, B (Rtrmv: a copy of the original
// x), and the result O; A is written too by Rswap and Rrot.
struct gpu_call {
    const blas_call<mpfr_class> *c;
    operand A, B, O;
    bool read_o, write_a;
};

gpu_call gpu_operands(const blas_call<mpfr_class> &c)
{
    gpu_call g;
    g.c = &c;
    g.A = c.a;
    g.B = c.b_is_copy_of_c ? c.c : c.b;
    g.B.region = mplapack_blas::REGION_ALL;
    g.O = c.c;
    g.read_o = c.c_read;
    g.write_a = c.a_written;
    return g;
}

bool present(const operand &o) { return o.X != NULL && o.rows > 0 && o.cols > 0; }

template <int PB> bool fixed_pack(const operand &o, bool read, mpfr_prec_t prec, std::vector<real_t<PB>> &out, mpfr_exp_t bound)
{
    out.resize((size_t)o.rows * o.cols + 1);
    for (long j = 0; j < o.cols; j++)
        for (long i = 0; i < o.rows; i++) {
            real_t<PB> &r = out[i + j * o.rows];
            r.set_zero();
            if (!o.ref(i, j))
                continue;
            mpfr_srcptr x = o.X[i + j * o.ld].get_mpfr_t();
            if (!read) {
                if (mpfr_get_prec(x) != prec)
                    return false;
            } else if (!pack<PB>(x, r, bound)) {
                return false;
            }
        }
    return true;
}

template <int PB> bool fixed_fits(const operand &o, const std::vector<real_t<PB>> &p)
{
    for (long j = 0; j < o.cols; j++)
        for (long i = 0; i < o.rows; i++)
            if (o.ref(i, j) && !fits<PB>(p[i + j * o.rows]))
                return false;
    return true;
}

template <int PB> void fixed_unpack(const operand &o, const std::vector<real_t<PB>> &p)
{
    for (long j = 0; j < o.cols; j++)
        for (long i = 0; i < o.rows; i++)
            if (o.ref(i, j))
                unpack<PB>(p[i + j * o.rows], const_cast<mpfr_class &>(o.X[i + j * o.ld]).get_mpfr_t());
}

bool rt_pack_operand(const operand &o, bool read, mpfr_prec_t prec, int nl, rt_storage &st, mpfr_exp_t bound)
{
    st.resize((size_t)o.rows * o.cols, nl);
    for (long j = 0; j < o.cols; j++)
        for (long i = 0; i < o.rows; i++) {
            const size_t idx = (size_t)(i + j * o.rows);
            st.kind[idx] = (signed char)MPFR_ZERO_KIND;
            if (!o.ref(i, j))
                continue;
            mpfr_srcptr x = o.X[i + j * o.ld].get_mpfr_t();
            if (!read) {
                if (mpfr_get_prec(x) != prec)
                    return false;
            } else if (!rt_pack_value(x, prec, nl, st, idx, bound)) {
                return false;
            }
        }
    return true;
}

bool rt_operand_fits(const operand &o, const rt_storage &st)
{
    for (long j = 0; j < o.cols; j++)
        for (long i = 0; i < o.rows; i++)
            if (o.ref(i, j) && !rt_fits(st, (size_t)(i + j * o.rows)))
                return false;
    return true;
}

// Any MPFR value (also Inf and NaN, possible in Rtrsm).
void rt_unpack_any(const rt_storage &st, size_t idx, mpfr_prec_t prec, int nl, mpfr_ptr x)
{
    const int kind = st.kind[idx];
    if (kind == MPFR_NAN_KIND || kind == -MPFR_NAN_KIND)
        mpfr_set_nan(x);
    else if (kind == MPFR_INF_KIND || kind == -MPFR_INF_KIND)
        mpfr_set_inf(x, kind > 0 ? 1 : -1);
    else
        rt_unpack(st, idx, prec, nl, x);
}

void rt_unpack_operand(const operand &o, const rt_storage &st, mpfr_prec_t prec, int nl)
{
    for (long j = 0; j < o.cols; j++)
        for (long i = 0; i < o.rows; i++)
            if (o.ref(i, j))
                rt_unpack_any(st, (size_t)(i + j * o.rows), prec, nl, const_cast<mpfr_class &>(o.X[i + j * o.ld]).get_mpfr_t());
}

// The shape with the leading dimensions of the packed operands.
mplapack_blas::blas_shape packed_shape(const gpu_call &g)
{
    mplapack_blas::blas_shape s = g.c->s;
    s.lda = present(g.A) ? g.A.rows : 1;
    s.ldb = present(g.B) ? g.B.rows : 1;
    s.ldc = g.O.rows;
    return s;
}

size_t count(const operand &o) { return present(o) ? (size_t)o.rows * o.cols : 0; }

template <int PB> bool run_fixed(const gpu_call &g, mpfr_prec_t prec)
{
    typedef real_t<PB> F;
    const mpfr_exp_t bound = input_exponent_bound();
    F a, b;
    if (!pack<PB>(g.c->alpha.get_mpfr_t(), a, bound) || !pack<PB>(g.c->beta.get_mpfr_t(), b, bound))
        return false;
    std::vector<F> pA(1), pB(1), pO;
    if ((present(g.A) && !fixed_pack<PB>(g.A, true, prec, pA, bound)) || (present(g.B) && !fixed_pack<PB>(g.B, true, prec, pB, bound)) || !fixed_pack<PB>(g.O, g.read_o, prec, pO, bound))
        return false;
    if (mplapack_mpfr_cuda::blas3_fixed<PB>(packed_shape(g), a, b, pA.data(), count(g.A), pB.data(), count(g.B), pO.data(), count(g.O)) != 0)
        return false;
    if (!fixed_fits<PB>(g.O, pO) || (g.write_a && !fixed_fits<PB>(g.A, pA)))
        return false;
    fixed_unpack<PB>(g.O, pO);
    if (g.write_a)
        fixed_unpack<PB>(g.A, pA);
    return true;
}

bool run_runtime(const gpu_call &g, mpfr_prec_t prec)
{
    const int nl = (int)((prec + 63) / 64);
    const bool trsm = g.c->s.op == mplapack_mpfr_cuda::L3_TRSM;
    if (trsm && (mpfr_get_emax() != mplapack_mpfr_cuda::RT_DEVICE_EMAX || mpfr_get_emin() != -mplapack_mpfr_cuda::RT_DEVICE_EMAX))
        return false;
    const mpfr_exp_t bound = trsm ? (mpfr_exp_t)mplapack_mpfr_cuda::RT_DEVICE_EMAX : rt_input_exponent_bound();
    rt_storage sa, sb, pA, pB, pO;
    sa.resize(1, nl);
    sb.resize(1, nl);
    pA.resize(1, nl);
    pB.resize(1, nl);
    if (!rt_pack_value(g.c->alpha.get_mpfr_t(), prec, nl, sa, 0, bound) || !rt_pack_value(g.c->beta.get_mpfr_t(), prec, nl, sb, 0, bound))
        return false;
    if ((present(g.A) && !rt_pack_operand(g.A, true, prec, nl, pA, bound)) || (present(g.B) && !rt_pack_operand(g.B, true, prec, nl, pB, bound)) || !rt_pack_operand(g.O, g.read_o, prec, nl, pO, bound))
        return false;
    if (mplapack_mpfr_cuda::blas3_rt(packed_shape(g), (long)prec, sa.pack, sb.pack, pA.pack, count(g.A), pB.pack, count(g.B), pO.pack, count(g.O)) != 0)
        return false;
    if (!trsm && (!rt_operand_fits(g.O, pO) || (g.write_a && !rt_operand_fits(g.A, pA))))
        return false;
    rt_unpack_operand(g.O, pO, prec, nl);
    if (g.write_a)
        rt_unpack_operand(g.A, pA, prec, nl);
    return true;
}

bool dispatch(const blas_call<mpfr_class> &c)
{
    const bool level3 = c.s.op <= mplapack_mpfr_cuda::L3_TRSM;
    if (!(level3 ? l3_enabled() : l12_enabled()) || c.work < (double)min_mnk())
        return false;
    if (mpfr_get_default_rounding_mode() != MPFR_RNDN)
        return false;
    std::lock_guard<std::mutex> lock(gpu_mutex());
    const gpu_call g = gpu_operands(c);
    const mpfr_prec_t prec = mpfrxx::default_precision_bits();
    const bool fixed_ok = c.s.op != mplapack_mpfr_cuda::L3_TRSM;
    if (fixed_ok && !force_runtime()) {
        switch (prec) {
#define MPLAPACK_MPFR_CUDA_L3_DISPATCH(PB) \
    case PB:                               \
        return run_fixed<PB>(g, prec);
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(MPLAPACK_MPFR_CUDA_L3_DISPATCH)
#undef MPLAPACK_MPFR_CUDA_L3_DISPATCH
        default:
            break;
        }
    }
    if (!runtime_enabled())
        return false;
    return run_runtime(g, prec);
}

bool ls(const char *a, const char *b) { return Mlsame_mpfr(a, b) != 0; }

mpfr_class *w(const mpfr_class *p) { return const_cast<mpfr_class *>(p); }

} // namespace

bool Rsyrk_mpfr_cuda(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class const beta, mpfr_class *c, mplapackint const ldc)
{
    blas_call<mpfr_class> x;
    return mplapack_blas::make_syrk<mpfr_class>(x, ls, false, uplo, trans, n, k, alpha, w(a), lda, (mpfr_class *)NULL, 1, beta, c, ldc) && dispatch(x);
}

bool Rsyr2k_mpfr_cuda(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *c, mplapackint const ldc)
{
    blas_call<mpfr_class> x;
    return mplapack_blas::make_syrk<mpfr_class>(x, ls, true, uplo, trans, n, k, alpha, w(a), lda, w(b), ldb, beta, c, ldc) && dispatch(x);
}

bool Rsymm_mpfr_cuda(const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *c, mplapackint const ldc)
{
    blas_call<mpfr_class> x;
    return mplapack_blas::make_symm<mpfr_class>(x, ls, side, uplo, m, n, alpha, w(a), lda, w(b), ldb, beta, c, ldc) && dispatch(x);
}

bool Rtrmm_mpfr_cuda(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb)
{
    blas_call<mpfr_class> x;
    return mplapack_blas::make_trxm<mpfr_class>(x, ls, false, side, uplo, transa, diag, m, n, alpha, w(a), lda, b, ldb) && dispatch(x);
}

bool Rtrsm_mpfr_cuda(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb)
{
    blas_call<mpfr_class> x;
    return mplapack_blas::make_trxm<mpfr_class>(x, ls, true, side, uplo, transa, diag, m, n, alpha, w(a), lda, b, ldb) && dispatch(x);
}

// ---- Level 1 and 2 (MPLAPACK_MPFR_CUDA_L12=1) ----

bool Rgemv_mpfr_cuda(const char *trans, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy)
{
    blas_call<mpfr_class> c;
    if (!mplapack_blas::make_gemv<mpfr_class>(c, ls, trans, m, n, alpha, a, lda, x, incx, beta, y, incy))
        return false;
    return dispatch(c);
}

bool Rsymv_mpfr_cuda(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy)
{
    blas_call<mpfr_class> c;
    if (!mplapack_blas::make_symv<mpfr_class>(c, ls, uplo, n, alpha, a, lda, x, incx, beta, y, incy))
        return false;
    return dispatch(c);
}

bool Rtrmv_mpfr_cuda(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx)
{
    blas_call<mpfr_class> c;
    if (!mplapack_blas::make_trmv<mpfr_class>(c, ls, uplo, trans, diag, n, a, lda, x, incx))
        return false;
    return dispatch(c);
}

bool Rger_mpfr_cuda(mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *a, mplapackint const lda)
{
    blas_call<mpfr_class> c;
    if (!mplapack_blas::make_rank<mpfr_class>(c, ls, 0, "U", m, n, alpha, x, incx, y, incy, a, lda))
        return false;
    return dispatch(c);
}

bool Rsyr_mpfr_cuda(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *a, mplapackint const lda)
{
    blas_call<mpfr_class> c;
    if (!mplapack_blas::make_rank<mpfr_class>(c, ls, 1, uplo, n, n, alpha, x, incx, (mpfr_class *)0, 1, a, lda))
        return false;
    return dispatch(c);
}

bool Rsyr2_mpfr_cuda(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *a, mplapackint const lda)
{
    blas_call<mpfr_class> c;
    if (!mplapack_blas::make_rank<mpfr_class>(c, ls, 2, uplo, n, n, alpha, x, incx, y, incy, a, lda))
        return false;
    return dispatch(c);
}

bool Rscal_mpfr_cuda(mplapackint const n, mpfr_class const da, mpfr_class *dx, mplapackint const incx)
{
    blas_call<mpfr_class> c;
    if (!mplapack_blas::make_level1<mpfr_class>(c, mplapack_blas::L1_SCAL, n, da, mpfr_class(0.0), dx, incx, (mpfr_class *)0, 1))
        return false;
    return dispatch(c);
}

bool Rswap_mpfr_cuda(mplapackint const n, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy)
{
    blas_call<mpfr_class> c;
    if (!mplapack_blas::make_level1<mpfr_class>(c, mplapack_blas::L1_SWAP, n, mpfr_class(0.0), mpfr_class(0.0), dx, incx, dy, incy))
        return false;
    return dispatch(c);
}

bool Rrot_mpfr_cuda(mplapackint const n, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy, mpfr_class const c_, mpfr_class const s)
{
    blas_call<mpfr_class> call;
    if (!mplapack_blas::make_level1<mpfr_class>(call, mplapack_blas::L1_ROT, n, c_, s, dx, incx, dy, incy))
        return false;
    return dispatch(call);
}
