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

// The fixed-precision Rgemm kernels (openmp/Rgemm_freal_*.cpp) against the
// OpenMP loops openmp/Rgemm_{NN,NT,TN,TT}_omp.cpp, bit for bit (signed
// zeros and precisions included): every precision with a kernel, shapes that
// cross the tiles and the k panels, random full-precision operands mixed with
// zeros of both signs, short significands (exact ties) and spread exponents;
// both kernels (x86-64 BMI2/ADX where the CPU has it, portable), with and
// without OpenMP; and the cases Rgemm_freal_cpu must leave to the loops.

#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>
#include <mpblas_mpfr.h>
#include "../openmp/Rgemm_freal.h"

void Rgemm_NN_omp(mplapackint m, mplapackint n, mplapackint k, mpfr_class alpha, mpfr_class *A, mplapackint lda, mpfr_class *B, mplapackint ldb, mpfr_class beta, mpfr_class *C, mplapackint ldc);
void Rgemm_NT_omp(mplapackint m, mplapackint n, mplapackint k, mpfr_class alpha, mpfr_class *A, mplapackint lda, mpfr_class *B, mplapackint ldb, mpfr_class beta, mpfr_class *C, mplapackint ldc);
void Rgemm_TN_omp(mplapackint m, mplapackint n, mplapackint k, mpfr_class alpha, mpfr_class *A, mplapackint lda, mpfr_class *B, mplapackint ldb, mpfr_class beta, mpfr_class *C, mplapackint ldc);
void Rgemm_TT_omp(mplapackint m, mplapackint n, mplapackint k, mpfr_class alpha, mpfr_class *A, mplapackint lda, mpfr_class *B, mplapackint ldb, mpfr_class beta, mpfr_class *C, mplapackint ldc);

static std::mt19937_64 gen(2026);
static long checks = 0, failures = 0;

// kind 0: random full precision; 1: also zeros (+0, -0), short significands
// (exact halfway sums) and exponents in [-40, 40]
static mpfr_class value(int kind)
{
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    mpfr_class x = u(gen);
    const int p = (int)mpfr_get_prec(x.get_mpfr_t());
    for (int t = 1; 53 * t < p + 53; t++)
        x += mpfr_class(u(gen)) * mpfr_class(std::ldexp(1.0, -53 * t));
    if (kind == 1) {
        const unsigned r = gen() % 8;
        if (r == 0)
            x = 0.0;
        else if (r == 1)
            x = -mpfr_class(0.0);
        else if (r == 2)
            x = (double)((long)(gen() % 7) - 3); // small integers: exact ties
        if (r >= 3)
            mpfr_mul_2si(x.get_mpfr_t(), x.get_mpfr_t(), (long)(gen() % 81) - 40, MPFR_RNDN);
    }
    return x;
}

static bool same(const mpfr_class &a, const mpfr_class &b)
{
    return mpfr_equal_p(a.get_mpfr_t(), b.get_mpfr_t()) && mpfr_signbit(a.get_mpfr_t()) == mpfr_signbit(b.get_mpfr_t()) &&
           mpfr_get_prec(a.get_mpfr_t()) == mpfr_get_prec(b.get_mpfr_t());
}

static void compare(const std::string &what, const std::vector<mpfr_class> &R, const std::vector<mpfr_class> &G)
{
    checks++;
    for (size_t t = 0; t < R.size(); t++)
        if (!same(R[t], G[t])) {
            failures++;
            std::printf("FAIL %s: element %ld differs\n", what.c_str(), (long)t);
            return;
        }
}

static void reference(bool na, bool nb, long m, long n, long k, const mpfr_class &alpha, std::vector<mpfr_class> &A, long lda, std::vector<mpfr_class> &B, long ldb,
                      const mpfr_class &beta, std::vector<mpfr_class> &C, long ldc)
{
    if (na && nb)
        Rgemm_NN_omp(m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);
    else if (na)
        Rgemm_NT_omp(m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);
    else if (nb)
        Rgemm_TN_omp(m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);
    else
        Rgemm_TT_omp(m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, C.data(), ldc);
}

static void run_precision(int prec)
{
    mpfrxx::set_default_precision_bits(prec);
    const bool x86 = Rgemm_freal_x86_available()
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
                     && __builtin_cpu_supports("bmi2") && __builtin_cpu_supports("adx")
#endif
        ;
    // tiles of 32 x 32; k panels of 102 steps at 512 bits, 30 at 2048 bits
    const long shapes[][3] = {{1, 1, 9}, {7, 5, 3}, {33, 31, 17}, {33, 3, 105}, {3, 33, 105}, {65, 2, 31}, {5, 66, 33}};
    for (auto &z : shapes)
        for (int form = 0; form < 4; form++)
            for (int kind = 0; kind < 2; kind++)
                for (int ak = 0; ak < 2; ak++)
                    for (int bk = 0; bk < 3; bk++) {
                        const long m = z[0], n = z[1], k = z[2];
                        const bool na = form < 2, nb = form % 2 == 0;
                        const long lda = (na ? m : k) + 2, ldb = (nb ? k : n) + 1, ldc = m + 3;
                        std::vector<mpfr_class> A(lda * (na ? k : m)), B(ldb * (nb ? n : k)), C(ldc * n);
                        for (auto &x : A) x = value(kind);
                        for (auto &x : B) x = value(kind);
                        for (auto &x : C) x = value(kind);
                        mpfr_class alpha = ak == 0 ? mpfr_class(1.0) : value(0), beta = bk == 0 ? mpfr_class(0.0) : bk == 1 ? mpfr_class(1.0) : value(kind);
                        if (alpha == 0.0)
                            alpha = 3.0;
                        const int beta_kind = beta == 0.0 ? 0 : beta == 1.0 ? 1 : 2;
                        std::vector<mpfr_class> R(C);
                        reference(na, nb, m, n, k, alpha, A, lda, B, ldb, beta, R, ldc);
                        const std::string tag = "prec " + std::to_string(prec) + " " + (na ? "N" : "T") + (nb ? "N" : "T") + " m=" + std::to_string(m) +
                                                " n=" + std::to_string(n) + " k=" + std::to_string(k) + " kind=" + std::to_string(kind) + " alpha=" +
                                                std::to_string(ak) + " beta=" + std::to_string(bk);
                        for (int par = 0; par < 2; par++) {
                            std::vector<mpfr_class> G(C);
                            if (!Rgemm_freal_generic(prec, na, nb, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, beta_kind, G.data(), ldc, par != 0)) {
                                failures++;
                                std::printf("FAIL %s: no portable kernel\n", tag.c_str());
                            }
                            compare(tag + " portable" + (par ? " omp" : ""), R, G);
                            if (x86) {
                                G = C;
                                Rgemm_freal_x86(prec, na, nb, m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, beta_kind, G.data(), ldc, par != 0);
                                compare(tag + " x86" + (par ? " omp" : ""), R, G);
                            }
                        }
                        // and through Rgemm
                        std::vector<mpfr_class> G(C);
                        Rgemm(na ? "N" : "T", nb ? "N" : "T", m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, G.data(), ldc);
                        compare(tag + " Rgemm", R, G);
                    }
}

// Rgemm_freal_cpu must return false (C untouched) for these
static void declined()
{
    mpfrxx::set_default_precision_bits(512);
    const long m = 20, n = 20, k = 20;
    std::vector<mpfr_class> A(m * k), B(k * n), C(m * n);
    for (auto &x : A) x = value(0);
    for (auto &x : B) x = value(0);
    for (auto &x : C) x = value(0);
    const mpfr_class alpha = value(0), beta = value(0);
    auto expect_false = [&](const char *what) {
        std::vector<mpfr_class> G(C);
        checks++;
        if (Rgemm_freal_cpu(true, true, m, n, k, alpha, A.data(), m, B.data(), k, beta, G.data(), m)) {
            failures++;
            std::printf("FAIL declined: %s was computed\n", what);
        }
    };
    {
        checks++;
        std::vector<mpfr_class> G(C);
        if (!Rgemm_freal_cpu(true, true, m, n, k, alpha, A.data(), m, B.data(), k, beta, G.data(), m)) {
            failures++;
            std::printf("FAIL declined: the plain call was not computed\n");
        }
    }
    const mpfr_class a0 = A[7];
    mpfr_set_nan(A[7].get_mpfr_t());
    expect_false("NaN in A");
    mpfr_set_inf(A[7].get_mpfr_t(), -1);
    expect_false("Inf in A");
    mpfr_set_prec(A[7].get_mpfr_t(), 256);
    mpfr_set_ui(A[7].get_mpfr_t(), 1, MPFR_RNDN);
    expect_false("a 256-bit element");
    A[7] = a0;
    mpfr_mul_2si(B[3].get_mpfr_t(), B[3].get_mpfr_t(), mpfr_get_emax() / 2, MPFR_RNDN);
    expect_false("a huge exponent");
    mpfr_div_2si(B[3].get_mpfr_t(), B[3].get_mpfr_t(), mpfr_get_emax() / 2, MPFR_RNDN);
    mpfr_set_default_rounding_mode(MPFR_RNDZ);
    expect_false("MPFR_RNDZ");
    mpfr_set_default_rounding_mode(MPFR_RNDN);
    mpfrxx::set_default_precision_bits(300);
    {
        std::vector<mpfr_class> A3(m * k), B3(k * n), C3(m * n);
        for (auto &x : A3) x = value(0);
        for (auto &x : B3) x = value(0);
        for (auto &x : C3) x = value(0);
        const mpfr_class al = value(0), be = value(0);
        checks++;
        if (Rgemm_freal_cpu(true, true, m, n, k, al, A3.data(), m, B3.data(), k, be, C3.data(), m)) {
            failures++;
            std::printf("FAIL declined: 300 bits was computed\n");
        }
    }
}

int main()
{
    const int precs[] = {128, 192, 256, 448, 512, 768, 1024, 2048};
    for (int p : precs)
        run_precision(p);
    declined();
    std::printf("Rgemm_freal: %ld checks, %ld failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
