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
 * Host side of the MPFR Level 3 BLAS on CUDA (Rsyrk, Rsyr2k, Rsymm, Rtrmm,
 * Rtrsm; kernels in blas3_kernel_cuda.h).  R*_mpfr_cuda() returns true when
 * it has computed the result on the GPU, and false, leaving every operand
 * untouched, when the call is left to the reference implementation:
 *
 *   - invalid arguments, quick returns and alpha == 0 (the reference reports
 *     or handles them);
 *   - MPLAPACK_MPFR_CUDA=0 or MPLAPACK_MPFR_CUDA_L3=0, rounding mode other
 *     than MPFR_RNDN, size below MPLAPACK_MPFR_CUDA_MIN_MNK (n*n*k for
 *     Rsyrk/Rsyr2k, m*m*n or m*n*n for the others);
 *   - an operand of another precision than the MPFR default precision p, or
 *     not finite, among the referenced elements;
 *   - Rsyrk, Rsyr2k, Rsymm, Rtrmm: an input exponent outside the bound that
 *     keeps every intermediate inside the exponent range, or a result outside
 *     it;
 *   - Rtrsm (always the cu_mpfr engine: cu_freal has no correctly rounded
 *     division): an MPFR exponent range other than the device default
 *     [1-2^30, 2^30-1], so that overflow and underflow happen exactly as on
 *     the host.
 *
 * p = 256/512/768/1024/2048 uses the cu_freal engine (zeros are +0), other
 * precisions (or MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1) the cu_mpfr engine,
 * whose results equal the reference bit for bit, signed zeros included.
 */

#include "mpfr_pack_cuda.h"
#include "blas3_kernel_cuda.h"

namespace {

using namespace mplapack_mpfr_cuda_host;
using mplapack_mpfr_cuda::l3_shape;

bool l3_enabled()
{
    static const bool enabled = env_long("MPLAPACK_MPFR_CUDA_L3", 1) != 0;
    return cuda_enabled() && enabled;
}

// An operand matrix: rows x cols starting at X with leading dimension ld;
// ref(i, j) tells whether element (i, j) is referenced.
struct operand {
    const mpfr_class *X;
    mplapackint ld;
    long rows, cols;
    int region; // 0: all, 1: upper (i <= j), 2: lower (i >= j); +4: diagonal excluded
    bool ref(long i, long j) const
    {
        if ((region & 4) && i == j)
            return false;
        switch (region & 3) {
        case 1:
            return i <= j;
        case 2:
            return i >= j;
        default:
            return true;
        }
    }
};

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

// The call: shape, scalars, operands A, B (may be empty) and the output
// operand O (C, or B for Rtrmm/Rtrsm), written in its referenced region.
struct l3_call {
    l3_shape s;
    const mpfr_class *alpha, *beta;
    operand A, B, O;
    bool read_o;            // O is read
    mpfr_class *out;        // == O.X
};

template <int PB> bool run_fixed(l3_call &c, mpfr_prec_t prec)
{
    typedef real_t<PB> F;
    const mpfr_exp_t bound = input_exponent_bound();
    F a, b;
    b.set_zero();
    if (!pack<PB>(c.alpha->get_mpfr_t(), a, bound))
        return false;
    if (c.beta && !pack<PB>(c.beta->get_mpfr_t(), b, bound))
        return false;
    std::vector<F> pA, pB, pO;
    if (!fixed_pack<PB>(c.A, true, prec, pA, bound) || (c.B.X && !fixed_pack<PB>(c.B, true, prec, pB, bound)) || !fixed_pack<PB>(c.O, c.read_o, prec, pO, bound))
        return false;
    if (!c.B.X)
        pB.resize(1);
    c.s.lda = c.A.rows, c.s.ldb = c.B.X ? c.B.rows : 1, c.s.ldc = c.O.rows;
    if (mplapack_mpfr_cuda::blas3_fixed<PB>(c.s, a, b, pA.data(), (size_t)c.A.rows * c.A.cols, pB.data(), c.B.X ? (size_t)c.B.rows * c.B.cols : 0, pO.data(), (size_t)c.O.rows * c.O.cols) != 0)
        return false;
    for (long j = 0; j < c.O.cols; j++)
        for (long i = 0; i < c.O.rows; i++)
            if (c.O.ref(i, j) && !fits<PB>(pO[i + j * c.O.rows]))
                return false;
    for (long j = 0; j < c.O.cols; j++)
        for (long i = 0; i < c.O.rows; i++)
            if (c.O.ref(i, j))
                unpack<PB>(pO[i + j * c.O.rows], c.out[i + j * c.O.ld].get_mpfr_t());
    return true;
}

bool run_runtime(l3_call &c, mpfr_prec_t prec)
{
    const int nl = (int)((prec + 63) / 64);
    const bool trsm = c.s.op == mplapack_mpfr_cuda::L3_TRSM;
    if (trsm && (mpfr_get_emax() != mplapack_mpfr_cuda::RT_DEVICE_EMAX || mpfr_get_emin() != -mplapack_mpfr_cuda::RT_DEVICE_EMAX))
        return false;
    const mpfr_exp_t bound = trsm ? (mpfr_exp_t)mplapack_mpfr_cuda::RT_DEVICE_EMAX : rt_input_exponent_bound();
    rt_storage sa, sb, pA, pB, pO;
    sa.resize(1, nl);
    sb.resize(1, nl);
    sb.kind[0] = (signed char)MPFR_ZERO_KIND;
    if (!rt_pack_value(c.alpha->get_mpfr_t(), prec, nl, sa, 0, bound))
        return false;
    if (c.beta && !rt_pack_value(c.beta->get_mpfr_t(), prec, nl, sb, 0, bound))
        return false;
    if (!rt_pack_operand(c.A, true, prec, nl, pA, bound) || (c.B.X && !rt_pack_operand(c.B, true, prec, nl, pB, bound)) || !rt_pack_operand(c.O, c.read_o, prec, nl, pO, bound))
        return false;
    if (!c.B.X)
        pB.resize(1, nl);
    c.s.lda = c.A.rows, c.s.ldb = c.B.X ? c.B.rows : 1, c.s.ldc = c.O.rows;
    if (mplapack_mpfr_cuda::blas3_rt(c.s, (long)prec, sa.pack, sb.pack, pA.pack, (size_t)c.A.rows * c.A.cols, pB.pack, c.B.X ? (size_t)c.B.rows * c.B.cols : 0, pO.pack, (size_t)c.O.rows * c.O.cols) != 0)
        return false;
    if (!trsm)
        for (long j = 0; j < c.O.cols; j++)
            for (long i = 0; i < c.O.rows; i++)
                if (c.O.ref(i, j) && !rt_fits(pO, (size_t)(i + j * c.O.rows)))
                    return false;
    for (long j = 0; j < c.O.cols; j++)
        for (long i = 0; i < c.O.rows; i++)
            if (c.O.ref(i, j))
                rt_unpack_any(pO, (size_t)(i + j * c.O.rows), prec, nl, c.out[i + j * c.O.ld].get_mpfr_t());
    return true;
}

bool dispatch(l3_call &c, double work)
{
    if (!l3_enabled() || work < (double)min_mnk())
        return false;
    if (mpfr_get_default_rounding_mode() != MPFR_RNDN)
        return false;
    const mpfr_prec_t prec = mpfrxx::default_precision_bits();
    const bool fixed_ok = c.s.op != mplapack_mpfr_cuda::L3_TRSM;
    if (fixed_ok && !force_runtime()) {
        switch (prec) {
#define MPLAPACK_MPFR_CUDA_L3_DISPATCH(PB) \
    case PB:                               \
        return run_fixed<PB>(c, prec);
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(MPLAPACK_MPFR_CUDA_L3_DISPATCH)
#undef MPLAPACK_MPFR_CUDA_L3_DISPATCH
        default:
            break;
        }
    }
    if (!runtime_enabled())
        return false;
    return run_runtime(c, prec);
}

bool is(const char *a, const char *b) { return Mlsame_mpfr(a, b) != 0; }

l3_call new_call(int op)
{
    l3_call c;
    c.s.op = op;
    c.s.left = c.s.upper = c.s.trans = c.s.nounit = 0;
    c.s.m = c.s.n = c.s.k = 0;
    c.s.lda = c.s.ldb = c.s.ldc = 1;
    c.s.alpha_one = c.s.beta_zero = c.s.beta_one = 0;
    c.alpha = c.beta = NULL;
    c.A.X = c.B.X = c.O.X = NULL;
    c.A.region = c.B.region = c.O.region = 0;
    c.read_o = true;
    c.out = NULL;
    return c;
}

} // namespace

bool Rsyrk_mpfr_cuda(const char *uplo, const char *trans, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, const mpfr_class &beta, mpfr_class *c, mplapackint ldc)
{
    const bool upper = is(uplo, "U"), notr = is(trans, "N");
    if ((!upper && !is(uplo, "L")) || (!notr && !is(trans, "T") && !is(trans, "C")) || n <= 0 || k <= 0)
        return false;
    const mplapackint nrowa = notr ? n : k;
    if (lda < std::max((mplapackint)1, nrowa) || ldc < std::max((mplapackint)1, n) || alpha == 0.0)
        return false;
    l3_call x = new_call(mplapack_mpfr_cuda::L3_SYRK);
    x.s.upper = upper, x.s.trans = !notr, x.s.n = n, x.s.k = k;
    x.s.beta_zero = (beta == 0.0), x.s.beta_one = (beta == 1.0);
    x.alpha = &alpha, x.beta = &beta;
    x.A = {a, lda, (long)nrowa, notr ? (long)k : (long)n, 0};
    x.O = {c, ldc, (long)n, (long)n, upper ? 1 : 2};
    x.read_o = !x.s.beta_zero;
    x.out = c;
    return dispatch(x, (double)n * n * k);
}

bool Rsyr2k_mpfr_cuda(const char *uplo, const char *trans, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, const mpfr_class *b, mplapackint ldb, const mpfr_class &beta, mpfr_class *c, mplapackint ldc)
{
    const bool upper = is(uplo, "U"), notr = is(trans, "N");
    if ((!upper && !is(uplo, "L")) || (!notr && !is(trans, "T") && !is(trans, "C")) || n <= 0 || k <= 0)
        return false;
    const mplapackint nrowa = notr ? n : k;
    if (lda < std::max((mplapackint)1, nrowa) || ldb < std::max((mplapackint)1, nrowa) || ldc < std::max((mplapackint)1, n) || alpha == 0.0)
        return false;
    l3_call x = new_call(mplapack_mpfr_cuda::L3_SYR2K);
    x.s.upper = upper, x.s.trans = !notr, x.s.n = n, x.s.k = k;
    x.s.beta_zero = (beta == 0.0), x.s.beta_one = (beta == 1.0);
    x.alpha = &alpha, x.beta = &beta;
    x.A = {a, lda, (long)nrowa, notr ? (long)k : (long)n, 0};
    x.B = {b, ldb, (long)nrowa, notr ? (long)k : (long)n, 0};
    x.O = {c, ldc, (long)n, (long)n, upper ? 1 : 2};
    x.read_o = !x.s.beta_zero;
    x.out = c;
    return dispatch(x, 2.0 * n * n * k);
}

bool Rsymm_mpfr_cuda(const char *side, const char *uplo, mplapackint m, mplapackint n, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, const mpfr_class *b, mplapackint ldb, const mpfr_class &beta, mpfr_class *c, mplapackint ldc)
{
    const bool left = is(side, "L"), upper = is(uplo, "U");
    if ((!left && !is(side, "R")) || (!upper && !is(uplo, "L")) || m <= 0 || n <= 0)
        return false;
    const mplapackint nrowa = left ? m : n;
    if (lda < std::max((mplapackint)1, nrowa) || ldb < std::max((mplapackint)1, m) || ldc < std::max((mplapackint)1, m) || alpha == 0.0)
        return false;
    l3_call x = new_call(mplapack_mpfr_cuda::L3_SYMM);
    x.s.left = left, x.s.upper = upper, x.s.m = m, x.s.n = n;
    x.s.beta_zero = (beta == 0.0), x.s.beta_one = (beta == 1.0);
    x.alpha = &alpha, x.beta = &beta;
    x.A = {a, lda, (long)nrowa, (long)nrowa, upper ? 1 : 2};
    x.B = {b, ldb, (long)m, (long)n, 0};
    x.O = {c, ldc, (long)m, (long)n, 0};
    x.read_o = !x.s.beta_zero;
    x.out = c;
    return dispatch(x, (double)m * n * nrowa);
}

static bool trxm(int op, const char *side, const char *uplo, const char *transa, const char *diag, mplapackint m, mplapackint n, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, mpfr_class *b, mplapackint ldb)
{
    const bool left = is(side, "L"), upper = is(uplo, "U"), nota = is(transa, "N"), nounit = is(diag, "N");
    if ((!left && !is(side, "R")) || (!upper && !is(uplo, "L")) || (!nota && !is(transa, "T") && !is(transa, "C")) || (!nounit && !is(diag, "U")) || m <= 0 || n <= 0)
        return false;
    const mplapackint nrowa = left ? m : n;
    if (lda < std::max((mplapackint)1, nrowa) || ldb < std::max((mplapackint)1, m) || alpha == 0.0)
        return false;
    l3_call x = new_call(op);
    x.s.left = left, x.s.upper = upper, x.s.trans = !nota, x.s.nounit = nounit, x.s.m = m, x.s.n = n;
    x.s.alpha_one = (alpha == 1.0);
    x.alpha = &alpha;
    x.A = {a, lda, (long)nrowa, (long)nrowa, (upper ? 1 : 2) | (nounit ? 0 : 4)};
    x.O = {b, ldb, (long)m, (long)n, 0};
    x.out = b;
    return dispatch(x, (double)m * n * nrowa);
}

bool Rtrmm_mpfr_cuda(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint m, mplapackint n, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, mpfr_class *b, mplapackint ldb)
{
    return trxm(mplapack_mpfr_cuda::L3_TRMM, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

bool Rtrsm_mpfr_cuda(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint m, mplapackint n, const mpfr_class &alpha, const mpfr_class *a, mplapackint lda, mpfr_class *b, mplapackint ldb)
{
    return trxm(mplapack_mpfr_cuda::L3_TRSM, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}
