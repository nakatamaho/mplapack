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
 * Bridge between mpfr_class matrices and the CUDA Rgemm.
 *
 * Rgemm_mpfr_cuda() returns true when it has computed C on the GPU, and false
 * (leaving C untouched) when the call is not eligible, so that Rgemm falls back
 * to the CPU code.  A call is eligible when
 *   - the MPFR default precision, alpha, beta and every referenced element of
 *     A, B and C have the same precision p,
 *   - the MPFR default rounding mode is MPFR_RNDN,
 *   - every referenced value is finite (C is not read when beta == 0),
 *   - exponents are small enough that no intermediate can leave the MPFR
 *     exponent range, and every result fits in it,
 *   - m*n*k >= MPLAPACK_MPFR_CUDA_MIN_MNK (default 32768), and
 *   - MPLAPACK_MPFR_CUDA is not set to 0.
 *
 * p = 256, 512, 768, 1024 or 2048 uses the fixed-precision kernels (cu_fp::cu_freal<p>,
 * Rgemm_device_cuda.cu); the result equals the CPU Rgemm except that a zero
 * is always +0.  Any other p uses the runtime-precision kernels (cu_mpfr,
 * Rgemm_rt_cuda.cu), whose result equals the CPU Rgemm exactly.
 * MPLAPACK_MPFR_CUDA_FORCE_RUNTIME=1 sends those calls to the
 * runtime-precision kernels too; MPLAPACK_MPFR_CUDA_RUNTIME=0 disables them.
 *
 * MPLAPACK_MPFR_CUDA_WINOGRAD_CUTOFF=N (N > 0) computes fixed-precision calls
 * with min(m, n, k) > N by the Winograd variant of Strassen's algorithm
 * (Rgemm_winograd_cuda.h), recursing down to blocks of size N.  The result
 * then differs from the CPU Rgemm within a normwise error bound.  It is off
 * by default.
 */

#include <mpblas_mpfr.h>
#include <cstdlib>
#include <vector>
#include "Rgemm_kernel_cuda.h"
#include "Rgemm_rt_cuda.h"
#include "Rgemm_winograd_cuda.h"

static_assert(GMP_NUMB_BITS == 64 && sizeof(mp_limb_t) == sizeof(cu_fp::cu_limb), "mpfr cuda Rgemm requires 64-bit GMP limbs");
static_assert(MPFR_ZERO_KIND == mplapack_mpfr_cuda::RT_ZERO_KIND && MPFR_REGULAR_KIND == mplapack_mpfr_cuda::RT_REGULAR_KIND, "MPFR custom kinds differ");

namespace {

using mplapack_mpfr_cuda::gemm_shape;
using mplapack_mpfr_cuda::real_t;

long env_long(const char *name, long fallback)
{
    const char *v = std::getenv(name);
    if (v == NULL || v[0] == '\0')
        return fallback;
    char *end = NULL;
    long r = std::strtol(v, &end, 10);
    return (end != v && *end == '\0') ? r : fallback;
}

bool cuda_enabled()
{
    static const bool enabled = env_long("MPLAPACK_MPFR_CUDA", 1) != 0;
    return enabled;
}

long min_mnk()
{
    static const long v = env_long("MPLAPACK_MPFR_CUDA_MIN_MNK", 32768);
    return v;
}

bool runtime_enabled()
{
    static const bool enabled = env_long("MPLAPACK_MPFR_CUDA_RUNTIME", 1) != 0;
    return enabled;
}

long winograd_cutoff()
{
    static const long v = env_long("MPLAPACK_MPFR_CUDA_WINOGRAD_CUTOFF", 0);
    return v;
}

bool force_runtime()
{
    static const bool forced = env_long("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME", 0) != 0;
    return forced;
}

// Bound on |exponent| of the inputs: alpha*op(B)*op(A) multiplies three
// values and sums at most k terms, so 4*bound + 64 must stay inside the range.
mpfr_exp_t input_exponent_bound()
{
    mpfr_exp_t emax = mpfr_get_emax();
    mpfr_exp_t emin = -mpfr_get_emin();
    mpfr_exp_t e = emax < emin ? emax : emin;
    return (e - 128) / 4;
}

template <int PB> bool pack(mpfr_srcptr x, real_t<PB> &r, mpfr_exp_t bound)
{
    if (mpfr_get_prec(x) != PB)
        return false;
    if (mpfr_zero_p(x)) {
        r.set_zero();
        return true;
    }
    if (!mpfr_regular_p(x))
        return false;
    mpfr_exp_t e = mpfr_get_exp(x);
    if (e > bound || e < -bound)
        return false;
    r.sign = mpfr_signbit(x) ? -1 : 1;
    r.exp = (long)e;
    const mp_limb_t *d = static_cast<const mp_limb_t *>(mpfr_custom_get_significand(x));
    for (int i = 0; i < real_t<PB>::N; i++)
        r.m[i] = d[i];
    return true;
}

template <int PB> bool fits(const real_t<PB> &r) { return r.is_zero() || (r.exp >= mpfr_get_emin() && r.exp <= mpfr_get_emax()); }

template <int PB> void unpack(const real_t<PB> &r, mpfr_ptr x)
{
    if (r.is_zero()) {
        mpfr_set_zero(x, 1);
        return;
    }
    mp_limb_t limbs[real_t<PB>::N];
    for (int i = 0; i < real_t<PB>::N; i++)
        limbs[i] = r.m[i];
    mpfr_t v;
    mpfr_custom_init_set(v, r.sign > 0 ? MPFR_REGULAR_KIND : -MPFR_REGULAR_KIND, (mpfr_exp_t)r.exp, PB, limbs);
    mpfr_set(x, v, MPFR_RNDN); // exact: x has precision PB
}

// Packs rows x cols elements of X (leading dimension ld) into out (leading dimension rows).
template <int PB> bool pack_matrix(const mpfr_class *X, mplapackint ld, long rows, long cols, std::vector<real_t<PB>> &out, mpfr_exp_t bound)
{
    out.resize((size_t)rows * cols);
    for (long j = 0; j < cols; j++)
        for (long i = 0; i < rows; i++)
            if (!pack<PB>(X[i + j * ld].get_mpfr_t(), out[i + j * rows], bound))
                return false;
    return true;
}

// Packs op(X) (rows x cols) with leading dimension rows; trans: X is stored transposed.
template <int PB> bool pack_op_matrix(const mpfr_class *X, mplapackint ld, bool trans, long rows, long cols, std::vector<real_t<PB>> &out, mpfr_exp_t bound)
{
    out.resize((size_t)rows * cols);
    for (long j = 0; j < cols; j++)
        for (long i = 0; i < rows; i++)
            if (!pack<PB>((trans ? X[j + i * ld] : X[i + j * ld]).get_mpfr_t(), out[i + j * rows], bound))
                return false;
    return true;
}

template <int PB> bool run(bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *A, mplapackint lda, const mpfr_class *B, mplapackint ldb, const mpfr_class &beta, mpfr_class *C, mplapackint ldc)
{
    typedef real_t<PB> F;
    const mpfr_exp_t bound = input_exponent_bound();
    gemm_shape s;
    s.nota = nota;
    s.notb = notb;
    s.m = m;
    s.n = n;
    s.k = k;
    s.lda = nota ? m : k;
    s.ldb = notb ? k : n;
    s.ldc = m;
    s.beta_zero = (beta == 0.0);
    s.beta_one = (beta == 1.0);

    F a, b;
    if (!pack<PB>(alpha.get_mpfr_t(), a, bound) || !pack<PB>(beta.get_mpfr_t(), b, bound))
        return false;
    const long cutoff = winograd_cutoff();
    const bool use_winograd = cutoff > 0 && m > cutoff && n > cutoff && k > cutoff;
    std::vector<F> pA, pB, pC;
    if (use_winograd) {
        // op(A) and op(B) explicitly: m x k and k x n
        if (!pack_op_matrix<PB>(A, lda, !nota, m, k, pA, bound) || !pack_op_matrix<PB>(B, ldb, !notb, k, n, pB, bound))
            return false;
    } else if (!pack_matrix<PB>(A, lda, s.lda, nota ? k : m, pA, bound) || !pack_matrix<PB>(B, ldb, s.ldb, notb ? n : k, pB, bound)) {
        return false;
    }
    if (s.beta_zero) {
        // C is overwritten without being read; only its precision matters.
        for (long j = 0; j < n; j++)
            for (long i = 0; i < m; i++)
                if (mpfr_get_prec(C[i + j * ldc].get_mpfr_t()) != PB)
                    return false;
        pC.resize((size_t)m * n);
    } else if (!pack_matrix<PB>(C, ldc, m, n, pC, bound)) {
        return false;
    }

    const int rc = use_winograd ? mplapack_mpfr_cuda::gemm_winograd<PB>(s, a, b, pA.data(), pB.data(), pC.data(), cutoff) : mplapack_mpfr_cuda::gemm<PB>(s, a, b, pA.data(), pB.data(), pC.data());
    if (rc != 0)
        return false;
    for (size_t t = 0; t < pC.size(); t++)
        if (!fits<PB>(pC[t]))
            return false;
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++)
            unpack<PB>(pC[i + j * m], C[i + j * ldc].get_mpfr_t());
    return true;
}

// ---- runtime precision (cu_mpfr) ----

// The device MPFR uses its default exponent range, which may be narrower
// than the host one.
mpfr_exp_t rt_input_exponent_bound()
{
    mpfr_exp_t e = input_exponent_bound();
    mpfr_exp_t d = (mplapack_mpfr_cuda::RT_DEVICE_EMAX - 128) / 4;
    return e < d ? e : d;
}

struct rt_storage {
    std::vector<unsigned long long> limbs;
    std::vector<long> exp;
    std::vector<signed char> kind;
    mplapack_mpfr_cuda::rt_pack pack;
    void resize(size_t count, int nl)
    {
        limbs.assign(count * nl + 1, 0);
        exp.assign(count + 1, 0);
        kind.assign(count + 1, 0);
        pack.limbs = limbs.data();
        pack.exp = exp.data();
        pack.kind = kind.data();
    }
};

bool rt_pack_value(mpfr_srcptr x, mpfr_prec_t prec, int nl, rt_storage &st, size_t idx, mpfr_exp_t bound)
{
    if (mpfr_get_prec(x) != prec)
        return false;
    int kind = mpfr_custom_get_kind(x);
    if (kind == MPFR_ZERO_KIND || kind == -MPFR_ZERO_KIND) {
        st.kind[idx] = (signed char)kind;
        return true;
    }
    if (kind != MPFR_REGULAR_KIND && kind != -MPFR_REGULAR_KIND)
        return false;
    mpfr_exp_t e = mpfr_get_exp(x);
    if (e > bound || e < -bound)
        return false;
    st.kind[idx] = (signed char)kind;
    st.exp[idx] = (long)e;
    const mp_limb_t *d = static_cast<const mp_limb_t *>(mpfr_custom_get_significand(x));
    for (int t = 0; t < nl; t++)
        st.limbs[idx * nl + t] = d[t];
    return true;
}

bool rt_pack_matrix(const mpfr_class *X, mplapackint ld, long rows, long cols, mpfr_prec_t prec, int nl, rt_storage &st, mpfr_exp_t bound)
{
    st.resize((size_t)rows * cols, nl);
    for (long j = 0; j < cols; j++)
        for (long i = 0; i < rows; i++)
            if (!rt_pack_value(X[i + j * ld].get_mpfr_t(), prec, nl, st, (size_t)(i + j * rows), bound))
                return false;
    return true;
}

bool rt_fits(const rt_storage &st, size_t idx)
{
    int kind = st.kind[idx];
    if (kind == MPFR_ZERO_KIND || kind == -MPFR_ZERO_KIND)
        return true;
    if (kind != MPFR_REGULAR_KIND && kind != -MPFR_REGULAR_KIND)
        return false;
    return st.exp[idx] >= mpfr_get_emin() && st.exp[idx] <= mpfr_get_emax();
}

void rt_unpack(const rt_storage &st, size_t idx, mpfr_prec_t prec, int nl, mpfr_ptr x)
{
    int kind = st.kind[idx];
    if (kind == MPFR_ZERO_KIND || kind == -MPFR_ZERO_KIND) {
        mpfr_set_zero(x, kind > 0 ? 1 : -1);
        return;
    }
    mpfr_t v;
    mpfr_custom_init_set(v, kind, (mpfr_exp_t)st.exp[idx], prec, const_cast<unsigned long long *>(&st.limbs[idx * nl]));
    mpfr_set(x, v, MPFR_RNDN); // exact: x has precision prec
}

bool run_rt(mpfr_prec_t prec, bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *A, mplapackint lda, const mpfr_class *B, mplapackint ldb, const mpfr_class &beta, mpfr_class *C, mplapackint ldc)
{
    const int nl = (int)((prec + 63) / 64);
    const mpfr_exp_t bound = rt_input_exponent_bound();
    gemm_shape s;
    s.nota = nota;
    s.notb = notb;
    s.m = m;
    s.n = n;
    s.k = k;
    s.lda = nota ? m : k;
    s.ldb = notb ? k : n;
    s.ldc = m;
    s.beta_zero = (beta == 0.0);
    s.beta_one = (beta == 1.0);

    rt_storage sa, sb, pA, pB, pC;
    sa.resize(1, nl);
    sb.resize(1, nl);
    if (!rt_pack_value(alpha.get_mpfr_t(), prec, nl, sa, 0, bound) || !rt_pack_value(beta.get_mpfr_t(), prec, nl, sb, 0, bound))
        return false;
    const long ncola = nota ? k : m, ncolb = notb ? n : k;
    if (!rt_pack_matrix(A, lda, s.lda, ncola, prec, nl, pA, bound) || !rt_pack_matrix(B, ldb, s.ldb, ncolb, prec, nl, pB, bound))
        return false;
    if (s.beta_zero) {
        for (long j = 0; j < n; j++)
            for (long i = 0; i < m; i++)
                if (mpfr_get_prec(C[i + j * ldc].get_mpfr_t()) != prec)
                    return false;
        pC.resize((size_t)m * n, nl);
    } else if (!rt_pack_matrix(C, ldc, m, n, prec, nl, pC, bound)) {
        return false;
    }

    if (mplapack_mpfr_cuda::gemm_rt(s, (long)prec, sa.pack, sb.pack, pA.pack, (size_t)s.lda * ncola, pB.pack, (size_t)s.ldb * ncolb, pC.pack) != 0)
        return false;
    for (size_t t = 0; t < (size_t)m * n; t++)
        if (!rt_fits(pC, t))
            return false;
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++)
            rt_unpack(pC, (size_t)(i + j * m), prec, nl, C[i + j * ldc].get_mpfr_t());
    return true;
}

} // namespace

bool Rgemm_mpfr_cuda(bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *A, mplapackint lda, const mpfr_class *B, mplapackint ldb, const mpfr_class &beta, mpfr_class *C, mplapackint ldc)
{
    if (!cuda_enabled() || (double)m * n * k < (double)min_mnk())
        return false;
    if (mpfr_get_default_rounding_mode() != MPFR_RNDN)
        return false;
    // The CPU code rounds temporaries to the default precision.
    const mpfr_prec_t prec = mpfrxx::default_precision_bits();
    if (!force_runtime()) {
        switch (prec) {
#define MPLAPACK_MPFR_CUDA_DISPATCH(PB) \
    case PB:                            \
        return run<PB>(nota, notb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
            MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(MPLAPACK_MPFR_CUDA_DISPATCH)
#undef MPLAPACK_MPFR_CUDA_DISPATCH
        default:
            break;
        }
    }
    if (!runtime_enabled())
        return false;
    return run_rt(prec, nota, notb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
