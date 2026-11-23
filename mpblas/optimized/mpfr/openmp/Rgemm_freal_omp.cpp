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

// Rgemm_freal_cpu (Rgemm_freal.h): the checks that make the fixed-precision
// kernels give the result of the MPFR operations bit for bit, and the
// choice of kernel.

#include <cstdlib>
#include <mpblas_mpfr.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Rgemm_freal.h"

namespace {

long freal_env(const char *name, long def)
{
    const char *s = std::getenv(name);
    return (s && *s) ? std::atol(s) : def;
}

bool freal_precision(long p)
{
    switch (p) {
#define MPFREAL_CASE(P) case P:
        MPFREAL_PRECISIONS(MPFREAL_CASE)
#undef MPFREAL_CASE
        return true;
    default:
        return false;
    }
}

// Every operand: precision p, a regular number or zero; [lo, hi] the range
// of the exponents of the nonzero ones.
struct operand_scan {
    mpfr_prec_t p;
    mpfr_exp_t lo, hi;
    bool ok;
    explicit operand_scan(mpfr_prec_t p_) : p(p_), lo(0), hi(0), ok(true) {}
    void add(mpfr_srcptr x)
    {
        if (mpfr_get_prec(x) != p) {
            ok = false;
            return;
        }
        if (mpfr_zero_p(x))
            return;
        if (!mpfr_regular_p(x)) {
            ok = false;
            return;
        }
        const mpfr_exp_t e = mpfr_get_exp(x);
        lo = e < lo ? e : lo;
        hi = e > hi ? e : hi;
    }
    void matrix(const mpfr_class *X, long rows, long cols, long ld)
    {
        for (long j = 0; j < cols && ok; j++)
            for (long i = 0; i < rows && ok; i++)
                add(X[i + j * ld].get_mpfr_t());
    }
};

} // namespace

bool Rgemm_freal_cpu(bool nota, bool notb, long m, long n, long k, const mpfr_class &alpha, const mpfr_class *A, long lda, const mpfr_class *B, long ldb,
                     const mpfr_class &beta, mpfr_class *C, long ldc)
{
    static const long mode = freal_env("MPLAPACK_MPFR_CPU_FREAL", 1);
    static const long min_work = freal_env("MPLAPACK_MPFR_CPU_FREAL_MIN", 512);
    if (mode == 0 || m < 1 || n < 1 || k < 1 || (double)m * n * k < (double)min_work)
        return false;
    const mpfr_prec_t p = mpfrxx::default_precision_bits();
    if (!freal_precision(p) || mpfr_get_default_rounding_mode() != MPFR_RNDN)
        return false;
    const int beta_kind = mpfr_zero_p(beta.get_mpfr_t()) ? 0 : mpfr_cmp_ui(beta.get_mpfr_t(), 1) == 0 ? 1 : 2;
    operand_scan s(p);
    s.add(alpha.get_mpfr_t());
    s.add(beta.get_mpfr_t());
    s.matrix(A, nota ? m : k, nota ? k : m, lda);
    s.matrix(B, notb ? k : n, notb ? n : k, ldb);
    if (beta_kind == 0) { // C is only written: the precision is what matters
        for (long j = 0; j < n && s.ok; j++)
            for (long i = 0; i < m && s.ok; i++)
                s.ok = mpfr_get_prec(C[i + j * ldc].get_mpfr_t()) == p;
    } else
        s.matrix(C, m, n, ldc);
    if (!s.ok)
        return false;
    // Every intermediate stays in [emin, emax]: a product of three factors
    // (alpha op(B) A) is within [3 lo - 3, 3 hi], a sum of k + 1 of them at
    // most ~log2(k) + 1 higher, and a nonzero sum or difference of such
    // numbers at most 2p below its smallest operand.
    const double top = 3.0 * (double)s.hi + 128.0, bottom = 3.0 * (double)s.lo - 2.0 * (double)p - 128.0;
    if (top > (double)mpfr_get_emax() || bottom < (double)mpfr_get_emin())
        return false;
#ifdef _OPENMP
    const bool parallel = omp_get_max_threads() >= 2 && (double)m * n * k >= (double)freal_env("MPLAPACK_BLAS_OMP_MIN", 32768) && !omp_in_parallel();
#else
    const bool parallel = false;
#endif
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    static const bool x86 = mode != 2 && Rgemm_freal_x86_available() && __builtin_cpu_supports("bmi2") && __builtin_cpu_supports("adx");
#else
    static const bool x86 = false;
#endif
    if (x86)
        return Rgemm_freal_x86((int)p, nota, notb, m, n, k, alpha, A, lda, B, ldb, beta, beta_kind, C, ldc, parallel);
    return Rgemm_freal_generic((int)p, nota, notb, m, n, k, alpha, A, lda, B, ldb, beta, beta_kind, C, ldc, parallel);
}
