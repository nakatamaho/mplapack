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
 * Bridge between mpfr_class matrices and the fixed-precision CUDA Rgemm.
 *
 * Rgemm_mpfr_cuda() returns true when it has computed C on the GPU, and false
 * (leaving C untouched) when the call is not eligible, so that Rgemm falls back
 * to the CPU code.  A call is eligible when
 *   - the MPFR default precision, alpha, beta and every referenced element of
 *     A, B and C have the same precision, 512 or 1024 bits,
 *   - every referenced value is finite (C is not read when beta == 0),
 *   - exponents are small enough that no intermediate can leave the MPFR
 *     exponent range, and every result fits in it,
 *   - m*n*k >= MPLAPACK_MPFR_CUDA_MIN_MNK (default 32768), and
 *   - MPLAPACK_MPFR_CUDA is not set to 0.
 * Under these conditions the result is identical to the CPU Rgemm, except
 * that a zero result is always +0 (cu_freal has no signed zero).
 */

#include <mpblas_mpfr.h>
#include <cstdlib>
#include <vector>
#include "Rgemm_kernel_cuda.h"

static_assert(GMP_NUMB_BITS == 64 && sizeof(mp_limb_t) == sizeof(cu_fp::cu_limb), "mpfr cuda Rgemm requires 64-bit GMP limbs");

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
    std::vector<F> pA, pB, pC;
    if (!pack_matrix<PB>(A, lda, s.lda, nota ? k : m, pA, bound) || !pack_matrix<PB>(B, ldb, s.ldb, notb ? n : k, pB, bound))
        return false;
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

    if (mplapack_mpfr_cuda::gemm<PB>(s, a, b, pA.data(), pB.data(), pC.data()) != 0)
        return false;
    for (size_t t = 0; t < pC.size(); t++)
        if (!fits<PB>(pC[t]))
            return false;
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++)
            unpack<PB>(pC[i + j * m], C[i + j * ldc].get_mpfr_t());
    return true;
}

} // namespace

bool Rgemm_mpfr_cuda(bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, const mpfr_class &alpha, const mpfr_class *A, mplapackint lda, const mpfr_class *B, mplapackint ldb, const mpfr_class &beta, mpfr_class *C, mplapackint ldc)
{
    if (!cuda_enabled() || (double)m * n * k < (double)min_mnk())
        return false;
    // The CPU code rounds temporaries to the default precision.
    switch (mpfrxx::default_precision_bits()) {
    case 512:
        return run<512>(nota, notb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    case 1024:
        return run<1024>(nota, notb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    default:
        return false;
    }
}
