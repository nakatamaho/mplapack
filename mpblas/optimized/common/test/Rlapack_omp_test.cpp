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

// Smoke test of the dedicated factorizations of the optimized libraries
// (mplapack/optimized/common/{R,C}{getrf,potrf}.cpp): the look-ahead loop
// against the reference routines (*_reference_cpu), bit for bit, info and
// pivots included; MPLAPACK_LAPACK_FAST=1 runs and its info agrees.

#include <mpblas.h>
#include <mplapack.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <string>
#include <vector>

void Rgetrf_reference_cpu(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info);
void Cgetrf_reference_cpu(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info);
void Rpotrf_reference_cpu(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, INTEGER &info);
void Cpotrf_reference_cpu(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER &info);

namespace {

long failures = 0, checks = 0;
std::mt19937_64 rng(20261023);
double u11() { return std::uniform_real_distribution<double>(-1.0, 1.0)(rng); }
REAL rnd()
{
    REAL x = u11();
    for (int t = 1; t < 6; t++)
        x += REAL(u11() * std::ldexp(1.0, -53 * t));
    return x;
}

bool same(const REAL &a, const REAL &b)
{
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    if (mpfr_nan_p(a.get_mpfr_t()) || mpfr_nan_p(b.get_mpfr_t()))
        return mpfr_nan_p(a.get_mpfr_t()) && mpfr_nan_p(b.get_mpfr_t());
    return mpfr_equal_p(a.get_mpfr_t(), b.get_mpfr_t()) && mpfr_signbit(a.get_mpfr_t()) == mpfr_signbit(b.get_mpfr_t());
#elif defined(MPLAPACK_BUILD_WITH_GMP)
    return mpf_cmp(a.get_mpf_t(), b.get_mpf_t()) == 0;
#elif defined(MPLAPACK_BUILD_WITH_BINARY80)
    return (a == b && std::signbit((long double)a) == std::signbit((long double)b)) || (a != a && b != b);
#else
    return std::memcmp(&a, &b, sizeof(REAL)) == 0;
#endif
}
bool same(const COMPLEX &a, const COMPLEX &b) { return same(REAL(a.real()), REAL(b.real())) && same(REAL(a.imag()), REAL(b.imag())); }

template <class T> void compare(const std::string &what, const std::vector<T> &ref, const std::vector<T> &got, INTEGER iref, INTEGER igot, const std::vector<INTEGER> &pref, const std::vector<INTEGER> &pgot)
{
    checks++;
    if (iref != igot || pref != pgot) {
        failures++;
        std::printf("FAIL %s: info %d/%d or pivots differ\n", what.c_str(), (int)iref, (int)igot);
        return;
    }
    for (size_t t = 0; t < ref.size(); t++)
        if (!same(ref[t], got[t])) {
            failures++;
            std::printf("FAIL %s: element %ld differs\n", what.c_str(), (long)t);
            return;
        }
}

} // namespace

int main()
{
    // the block size of iMlaenv is 64: several blocks, odd sizes
    const long sizes[][2] = {{70, 70}, {200, 130}, {130, 200}, {257, 257}};
    for (int fast = 0; fast < 2; fast++) {
        setenv("MPLAPACK_LAPACK_FAST", fast ? "1" : "0", 1);
        setenv("MPLAPACK_LAPACK_LOOKAHEAD", "2", 1); // the task form at every size
        for (auto &z : sizes) {
            const long m = z[0], n = z[1], lda = m + 3;
            {
                std::vector<REAL> A(lda * n);
                for (auto &x : A)
                    x = rnd();
                if (m > 100)
                    A[5 + 3 * lda] = 0.0; // a zero on the way
                std::vector<REAL> R(A), G(A);
                std::vector<INTEGER> pr(std::min(m, n)), pg(std::min(m, n));
                INTEGER ir = 0, ig = 0;
                Rgetrf_reference_cpu(m, n, R.data(), lda, pr.data(), ir);
                Rgetrf(m, n, G.data(), lda, pg.data(), ig);
                if (!fast)
                    compare("Rgetrf " + std::to_string(m) + "x" + std::to_string(n), R, G, ir, ig, pr, pg);
                else if (ir != ig) {
                    checks++, failures++;
                    std::printf("FAIL Rgetrf fast info\n");
                }
            }
            {
                std::vector<COMPLEX> A(lda * n);
                for (auto &x : A)
                    x = COMPLEX(rnd(), rnd());
                std::vector<COMPLEX> R(A), G(A);
                std::vector<INTEGER> pr(std::min(m, n)), pg(std::min(m, n));
                INTEGER ir = 0, ig = 0;
                Cgetrf_reference_cpu(m, n, R.data(), lda, pr.data(), ir);
                Cgetrf(m, n, G.data(), lda, pg.data(), ig);
                if (!fast)
                    compare("Cgetrf " + std::to_string(m) + "x" + std::to_string(n), R, G, ir, ig, pr, pg);
            }
            if (m != n)
                continue;
            for (const char *uplo : {"U", "L"})
                for (int bad = 0; bad < 2; bad++) {
                    // SPD (bad: not positive definite from row 100 on)
                    std::vector<REAL> A(lda * n);
                    std::vector<COMPLEX> Z(lda * n);
                    for (auto &x : A)
                        x = rnd();
                    for (auto &x : Z)
                        x = COMPLEX(rnd(), rnd());
                    std::vector<REAL> S(lda * n);
                    std::vector<COMPLEX> H(lda * n);
                    for (long j = 0; j < n; j++)
                        for (long i = 0; i < n; i++) {
                            S[i + j * lda] = A[i + j * lda] + A[j + i * lda];
                            H[i + j * lda] = Z[i + j * lda] + conj(Z[j + i * lda]);
                        }
                    for (long i = 0; i < n; i++) {
                        const REAL d = (bad && i >= 100) ? REAL(-1.0) : REAL(2.0 * n);
                        S[i + i * lda] = S[i + i * lda] + d;
                        H[i + i * lda] = COMPLEX(H[i + i * lda].real() + d, REAL(0.0));
                    }
                    std::vector<INTEGER> none;
                    std::vector<REAL> R(S), G(S);
                    INTEGER ir = 0, ig = 0;
                    Rpotrf_reference_cpu(uplo, n, R.data(), lda, ir);
                    Rpotrf(uplo, n, G.data(), lda, ig);
                    if (!fast)
                        compare(std::string("Rpotrf ") + uplo + " " + std::to_string(n) + (bad ? " bad" : ""), R, G, ir, ig, none, none);
                    else if ((ir == 0) != (ig == 0)) {
                        checks++, failures++;
                        std::printf("FAIL Rpotrf fast info %d %d\n", (int)ir, (int)ig);
                    }
                    std::vector<COMPLEX> RC(H), GC(H);
                    ir = ig = 0;
                    Cpotrf_reference_cpu(uplo, n, RC.data(), lda, ir);
                    Cpotrf(uplo, n, GC.data(), lda, ig);
                    if (!fast)
                        compare(std::string("Cpotrf ") + uplo + " " + std::to_string(n) + (bad ? " bad" : ""), RC, GC, ir, ig, none, none);
                }
        }
    }
    std::printf("dedicated LAPACK: %ld checks, %ld failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
