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


// Smoke test of the opt-in Winograd variant of Rgemm and Cgemm of every
// optimized library (<P>_GEMM_WINOGRAD_CUTOFF, <P>_CGEMM_WINOGRAD_CUTOFF;
// ../Rgemm_cpu_common.h): all transposes (and conjugates), odd sizes,
// leading dimensions larger than needed, beta in {0, 1, general}, cutoffs
// 1, 2, 3, 8.  Small integer matrices give the conventional result exactly;
// random matrices agree with it within 2^(20-p) (k |alpha| max|A| max|B| +
// |beta| max|C|), p the precision in bits (64 for MPFR and GMP).

#include <mpblas.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <string>
#include <vector>

namespace {

#if defined(MPLAPACK_BUILD_WITH_DOUBLE)
const char *PREFIX = "MPLAPACK_DOUBLE";
const int PREC = 53;
#elif defined(MPLAPACK_BUILD_WITH_BINARY80)
const char *PREFIX = "MPLAPACK_BINARY80";
const int PREC = 64;
#elif defined(MPLAPACK_BUILD_WITH_BINARY128)
const char *PREFIX = "MPLAPACK_BINARY128";
const int PREC = 113;
#elif defined(MPLAPACK_BUILD_WITH_DD)
const char *PREFIX = "MPLAPACK_DD";
const int PREC = 106;
#elif defined(MPLAPACK_BUILD_WITH_TD)
const char *PREFIX = "MPLAPACK_TD";
const int PREC = 159;
#elif defined(MPLAPACK_BUILD_WITH_QD)
const char *PREFIX = "MPLAPACK_QD";
const int PREC = 212;
#elif defined(MPLAPACK_BUILD_WITH_MPFR)
const char *PREFIX = "MPLAPACK_MPFR";
const int PREC = 64;
#elif defined(MPLAPACK_BUILD_WITH_GMP)
const char *PREFIX = "MPLAPACK_GMP";
const int PREC = 64;
#endif

long failures = 0, checks = 0, differ = 0; // differ: random cases off the conventional result (the variant ran)
std::mt19937_64 rng(20261023);
double u11() { return std::uniform_real_distribution<double>(-1.0, 1.0)(rng); }
REAL rnd(bool integers)
{
    if (integers)
        return REAL((double)((long)(rng() % 129) - 64));
    REAL x = u11();
    for (int t = 1; t < 4; t++)
        x += REAL(u11() * std::ldexp(1.0, -53 * t));
    return x;
}
template <class T> T make(bool integers);
template <> REAL make<REAL>(bool integers) { return rnd(integers); }
template <> COMPLEX make<COMPLEX>(bool integers)
{
    const REAL re = rnd(integers), im = rnd(integers);
    return COMPLEX(re, im);
}

REAL absr(const REAL &x)
{
    REAL a = x;
    if (a < REAL(0.0))
        a = -a;
    return a;
}
// the larger of the absolute parts
REAL mag(const REAL &x) { return absr(x); }
REAL mag(const COMPLEX &x)
{
    const REAL a = absr(REAL(x.real())), b = absr(REAL(x.imag()));
    return a > b ? a : b;
}

void setenv_s(const char *name, const char *value)
{
    const std::string v = std::string(PREFIX) + name;
    if (value)
        setenv(v.c_str(), value, 1);
    else
        unsetenv(v.c_str());
}

void gemm(const char *ta, const char *tb, INTEGER m, INTEGER n, INTEGER k, REAL alpha, REAL *A, INTEGER lda, REAL *B, INTEGER ldb, REAL beta, REAL *C, INTEGER ldc) { Rgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc); }
void gemm(const char *ta, const char *tb, INTEGER m, INTEGER n, INTEGER k, COMPLEX alpha, COMPLEX *A, INTEGER lda, COMPLEX *B, INTEGER ldb, COMPLEX beta, COMPLEX *C, INTEGER ldc) { Cgemm(ta, tb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc); }

template <class T> void run(const char *var, const std::vector<std::string> &trans)
{
    const long sizes[][3] = {{9, 7, 11}, {16, 16, 16}, {37, 29, 33}};
    const char *cutoffs[] = {"1", "2", "3", "8"};
    for (const std::string &ta : trans)
        for (const std::string &tb : trans)
            for (const auto &sz : sizes)
                for (int bm = 0; bm < 3; bm++)
                    for (int integers = 0; integers < 2; integers++) {
                        const long m = sz[0], n = sz[1], k = sz[2];
                        const long nra = ta == "N" ? m : k, nca = ta == "N" ? k : m, nrb = tb == "N" ? k : n, ncb = tb == "N" ? n : k;
                        const long lda = nra + 2, ldb = nrb + 1, ldc = m + 3;
                        std::vector<T> A(lda * nca), B(ldb * ncb), C(ldc * n);
                        for (auto &v : A)
                            v = make<T>(integers);
                        for (auto &v : B)
                            v = make<T>(integers);
                        for (auto &v : C)
                            v = make<T>(integers);
                        const T alpha = integers ? make<T>(true) : make<T>(false);
                        const T beta = bm == 0 ? T(REAL(0.0)) : (bm == 1 ? T(REAL(1.0)) : make<T>(integers));
                        setenv_s(var, NULL);
                        std::vector<T> ref(C);
                        gemm(ta.c_str(), tb.c_str(), m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, ref.data(), ldc);
                        REAL ma(0.0), mb(0.0), mc(0.0);
                        for (auto &v : A)
                            if (mag(v) > ma)
                                ma = mag(v);
                        for (auto &v : B)
                            if (mag(v) > mb)
                                mb = mag(v);
                        for (auto &v : C)
                            if (mag(v) > mc)
                                mc = mag(v);
                        const REAL bound = REAL(std::ldexp(1.0, 20 - PREC)) * (REAL((double)(2 * k)) * mag(alpha) * ma * mb + mag(beta) * mc);
                        for (const char *cut : cutoffs) {
                            setenv_s(var, cut);
                            std::vector<T> got(C);
                            gemm(ta.c_str(), tb.c_str(), m, n, k, alpha, A.data(), lda, B.data(), ldb, beta, got.data(), ldc);
                            bool ok = true, dif = false;
                            for (long j = 0; j < n && ok; j++)
                                for (long i = 0; i < ldc && ok; i++) {
                                    const T d = got[i + j * ldc] - ref[i + j * ldc];
                                    if (i >= m || integers)
                                        ok = mag(d) == REAL(0.0); // rows past m untouched, integers exact
                                    else {
                                        ok = mag(d) <= bound;
                                        dif = dif || mag(d) != REAL(0.0);
                                    }
                                }
                            checks++;
                            differ += dif;
                            if (!ok) {
                                failures++;
                                std::printf("FAIL %s %s%s %ldx%ldx%ld beta%d %s cutoff %s\n", var, ta.c_str(), tb.c_str(), m, n, k, bm, integers ? "integers" : "random", cut);
                            }
                        }
                        setenv_s(var, NULL);
                    }
}

} // namespace

int main()
{
    run<REAL>("_GEMM_WINOGRAD_CUTOFF", {"N", "T"});
    run<COMPLEX>("_CGEMM_WINOGRAD_CUTOFF", {"N", "T", "C"});
    if (differ == 0) {
        failures++;
        std::printf("FAIL no random case differs from the conventional result: the variant did not run\n");
    }
    std::printf("Winograd Rgemm/Cgemm: %ld checks, %ld failures (%ld random cases off the conventional result)\n", checks, failures, differ);
    return failures ? 1 : 0;
}
