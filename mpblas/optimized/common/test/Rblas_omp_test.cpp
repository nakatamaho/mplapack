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
 * Test of the OpenMP lane BLAS of the optimized libraries (common/blas/*):
 * every routine must give the result of the reference implementation
 * (R*_reference_cpu, the unchanged reference code) bit for bit, for every
 * option, alpha/beta in {1, 0, general}, odd sizes, leading dimensions
 * larger than needed and increments 1, 2 and -2.  Compiled once per type
 * with -DMPLAPACK_BUILD_WITH_<TYPE> and linked with libmplapack_<type>_opt.
 */

#include <mpblas.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <string>
#include <vector>

void Rsyrk_reference_cpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL const beta, REAL *c, INTEGER const ldc);
void Rsyr2k_reference_cpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const beta, REAL *c, INTEGER const ldc);
void Rsymm_reference_cpu(const char *side, const char *uplo, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const beta, REAL *c, INTEGER const ldc);
void Rtrmm_reference_cpu(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb);
void Rtrsm_reference_cpu(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb);
void Rgemv_reference_cpu(const char *trans, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy);
void Rsymv_reference_cpu(const char *uplo, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy);
void Rtrmv_reference_cpu(const char *uplo, const char *trans, const char *diag, INTEGER const n, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx);
void Rger_reference_cpu(INTEGER const m, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *y, INTEGER const incy, REAL *a, INTEGER const lda);
void Rsyr_reference_cpu(const char *uplo, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *a, INTEGER const lda);
void Rsyr2_reference_cpu(const char *uplo, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *y, INTEGER const incy, REAL *a, INTEGER const lda);
void Rscal_reference_cpu(INTEGER const n, REAL const da, REAL *dx, INTEGER const incx);
void Rswap_reference_cpu(INTEGER const n, REAL *dx, INTEGER const incx, REAL *dy, INTEGER const incy);
void Rgemmtr_reference_cpu(const char *uplo, const char *transa, const char *transb, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const beta, REAL *c, INTEGER const ldc);
void Rrot_reference_cpu(INTEGER const n, REAL *dx, INTEGER const incx, REAL *dy, INTEGER const incy, REAL const c, REAL const s);

void Rlaswp(INTEGER const n, REAL *a, INTEGER const lda, INTEGER const k1, INTEGER const k2, INTEGER *ipiv, INTEGER const incx);
void Rlaswp_reference_cpu(INTEGER const n, REAL *a, INTEGER const lda, INTEGER const k1, INTEGER const k2, INTEGER *ipiv, INTEGER const incx);

namespace {

std::mt19937_64 rng(20261022);
double u11() { return std::uniform_real_distribution<double>(-1.0, 1.0)(rng); }
long failures = 0, checks = 0;

REAL rnd()
{
    if (rng() % 9 == 0)
        return REAL(0.0);
    // random bits down to 2^-636, so that MPFR/GMP products and sums round
    REAL x = REAL(u11());
    for (int t = 1; t < 12; t++)
        x += REAL(u11() * std::ldexp(1.0, -53 * t));
    return x;
}

void fill(std::vector<REAL> &v)
{
    for (auto &x : v)
        x = rnd();
}

// well-conditioned triangle for Rtrsm
void fill_tri(std::vector<REAL> &A, long n, long lda)
{
    fill(A);
    for (long j = 0; j < n; j++)
        for (long i = 0; i < n; i++)
            if (i == j)
                A[i + j * lda] = REAL(1.0) + REAL(std::fabs(u11()));
            else
                A[i + j * lda] = A[i + j * lda] / REAL(2.0 * n);
}

bool same(const REAL &a, const REAL &b)
{
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    if (mpfr_nan_p(a.get_mpfr_t()) || mpfr_nan_p(b.get_mpfr_t()))
        return mpfr_nan_p(a.get_mpfr_t()) && mpfr_nan_p(b.get_mpfr_t());
    return mpfr_equal_p(a.get_mpfr_t(), b.get_mpfr_t()) && mpfr_signbit(a.get_mpfr_t()) == mpfr_signbit(b.get_mpfr_t()) && mpfr_get_prec(a.get_mpfr_t()) == mpfr_get_prec(b.get_mpfr_t());
#elif defined(MPLAPACK_BUILD_WITH_GMP)
    return mpf_cmp(a.get_mpf_t(), b.get_mpf_t()) == 0;
#elif defined(MPLAPACK_BUILD_WITH_BINARY80)
    return (a == b && std::signbit((long double)a) == std::signbit((long double)b)) || (a != a && b != b);
#else
    return std::memcmp(&a, &b, sizeof(REAL)) == 0;
#endif
}

void compare(const std::string &what, const std::vector<REAL> &ref, const std::vector<REAL> &got)
{
    checks++;
    for (size_t t = 0; t < ref.size(); t++)
        if (!same(ref[t], got[t])) {
            failures++;
            std::printf("FAIL %s: element %ld differs\n", what.c_str(), (long)t);
            return;
        }
}

REAL scalar(int kind)
{
    if (kind == 0)
        return REAL(1.0);
    if (kind == 1)
        return REAL(0.0);
    REAL x;
    do
        x = rnd();
    while (x == 0.0);
    return x;
}

std::string tag(const char *r, std::initializer_list<const char *> opts, std::initializer_list<long> nums)
{
    std::string s = r;
    for (const char *o : opts)
        s += std::string(" ") + o;
    for (long v : nums)
        s += " " + std::to_string(v);
    return s;
}

} // namespace

int main()
{
    setenv("MPLAPACK_BLAS_OMP_MIN", "1", 1);
#ifdef RBLAS_TEST_GPU
    // Rblas_gpu_test.cpp: every call, Level 1 and 2 included, on the GPU
    rblas_test_gpu_setup();
    if (!rblas_test_have_gpu()) {
        std::printf("no CUDA device: skipped\n");
        return 77;
    }
#endif
    const char *ul[] = {"U", "L"}, *tr[] = {"N", "T"}, *sd[] = {"L", "R"}, *dg[] = {"N", "U"};
    const long sizes[][2] = {{1, 1}, {7, 5}, {19, 23}};
    const long incs[] = {1, 2, -2};
    for (auto &z : sizes) {
        const long m = z[0], n = z[1];
        for (int ak = 0; ak < 3; ak++)
            for (int bk = 0; bk < 3; bk++) {
                const REAL alpha = scalar(ak), beta = scalar(bk);
                // Level 3
                for (const char *u : ul)
                    for (const char *t : tr) {
                        const bool nt = t[0] == 'N';
                        const long lda = (nt ? n : m) + 2, ca = nt ? m : n;
                        std::vector<REAL> A(lda * ca), B(lda * ca), C((n + 1) * n);
                        fill(A), fill(B), fill(C);
                        std::vector<REAL> R(C), G(C);
                        Rsyrk_reference_cpu(u, t, n, m, alpha, A.data(), lda, beta, R.data(), n + 1);
                        Rsyrk(u, t, n, m, alpha, A.data(), lda, beta, G.data(), n + 1);
                        compare(tag("Rsyrk", {u, t}, {n, m, ak, bk}), R, G);
                        R = C, G = C;
                        Rsyr2k_reference_cpu(u, t, n, m, alpha, A.data(), lda, B.data(), lda, beta, R.data(), n + 1);
                        Rsyr2k(u, t, n, m, alpha, A.data(), lda, B.data(), lda, beta, G.data(), n + 1);
                        compare(tag("Rsyr2k", {u, t}, {n, m, ak, bk}), R, G);
                        for (const char *t2 : tr) {
                            const long rb = t2[0] == 'N' ? m : n, cb = t2[0] == 'N' ? n : m;
                            std::vector<REAL> B2((rb + 1) * cb);
                            fill(B2);
                            R = C, G = C;
                            Rgemmtr_reference_cpu(u, t, t2, n, m, alpha, A.data(), lda, B2.data(), rb + 1, beta, R.data(), n + 1);
                            Rgemmtr(u, t, t2, n, m, alpha, A.data(), lda, B2.data(), rb + 1, beta, G.data(), n + 1);
                            compare(tag("Rgemmtr", {u, t, t2}, {n, m, ak, bk}), R, G);
                        }
                    }
                for (const char *s : sd)
                    for (const char *u : ul) {
                        const long na = s[0] == 'L' ? m : n, lda = na + 1;
                        std::vector<REAL> A(lda * na), B((m + 2) * n), C((m + 3) * n);
                        fill(A), fill(B), fill(C);
                        std::vector<REAL> R(C), G(C);
                        Rsymm_reference_cpu(s, u, m, n, alpha, A.data(), lda, B.data(), m + 2, beta, R.data(), m + 3);
                        Rsymm(s, u, m, n, alpha, A.data(), lda, B.data(), m + 2, beta, G.data(), m + 3);
                        compare(tag("Rsymm", {s, u}, {m, n, ak, bk}), R, G);
                        if (bk != 0)
                            continue;
                        for (const char *t : tr)
                            for (const char *d : dg) {
                                std::vector<REAL> T(lda * na);
                                fill_tri(T, na, lda);
                                R = B, G = B;
                                Rtrmm_reference_cpu(s, u, t, d, m, n, alpha, T.data(), lda, R.data(), m + 2);
                                Rtrmm(s, u, t, d, m, n, alpha, T.data(), lda, G.data(), m + 2);
                                compare(tag("Rtrmm", {s, u, t, d}, {m, n, ak}), R, G);
                                R = B, G = B;
                                Rtrsm_reference_cpu(s, u, t, d, m, n, alpha, T.data(), lda, R.data(), m + 2);
                                Rtrsm(s, u, t, d, m, n, alpha, T.data(), lda, G.data(), m + 2);
                                compare(tag("Rtrsm", {s, u, t, d}, {m, n, ak}), R, G);
                            }
                    }
                // Level 2 and 1
                for (long ix : incs)
                    for (long iy : incs) {
                        const long lx = std::labs(ix), ly = std::labs(iy);
                        std::vector<REAL> A((m + 2) * std::max(m, n)), X(lx * std::max(m, n) + 1), Y(ly * std::max(m, n) + 1);
                        fill(A), fill(X), fill(Y);
                        for (const char *t : tr) {
                            std::vector<REAL> R(Y), G(Y);
                            Rgemv_reference_cpu(t, m, n, alpha, A.data(), m + 2, X.data(), ix, beta, R.data(), iy);
                            Rgemv(t, m, n, alpha, A.data(), m + 2, X.data(), ix, beta, G.data(), iy);
                            compare(tag("Rgemv", {t}, {m, n, ix, iy, ak, bk}), R, G);
                        }
                        const long nn = std::min(m, n);
                        for (const char *u : ul) {
                            std::vector<REAL> R(Y), G(Y);
                            Rsymv_reference_cpu(u, nn, alpha, A.data(), m + 2, X.data(), ix, beta, R.data(), iy);
                            Rsymv(u, nn, alpha, A.data(), m + 2, X.data(), ix, beta, G.data(), iy);
                            compare(tag("Rsymv", {u}, {nn, ix, iy, ak, bk}), R, G);
                            for (const char *t : tr)
                                for (const char *d : dg) {
                                    R = X, G = X;
                                    Rtrmv_reference_cpu(u, t, d, nn, A.data(), m + 2, R.data(), ix);
                                    Rtrmv(u, t, d, nn, A.data(), m + 2, G.data(), ix);
                                    compare(tag("Rtrmv", {u, t, d}, {nn, ix}), R, G);
                                }
                            R = A, G = A;
                            Rsyr_reference_cpu(u, nn, alpha, X.data(), ix, R.data(), m + 2);
                            Rsyr(u, nn, alpha, X.data(), ix, G.data(), m + 2);
                            compare(tag("Rsyr", {u}, {nn, ix, ak}), R, G);
                            R = A, G = A;
                            Rsyr2_reference_cpu(u, nn, alpha, X.data(), ix, Y.data(), iy, R.data(), m + 2);
                            Rsyr2(u, nn, alpha, X.data(), ix, Y.data(), iy, G.data(), m + 2);
                            compare(tag("Rsyr2", {u}, {nn, ix, iy, ak}), R, G);
                        }
                        {
                            std::vector<REAL> R(A), G(A);
                            Rger_reference_cpu(m, n, alpha, X.data(), ix, Y.data(), iy, R.data(), m + 2);
                            Rger(m, n, alpha, X.data(), ix, Y.data(), iy, G.data(), m + 2);
                            compare(tag("Rger", {}, {m, n, ix, iy, ak}), R, G);
                        }
                        {
                            std::vector<REAL> R(X), G(X);
                            if (ix > 0) {
                                Rscal_reference_cpu(n, alpha, R.data(), ix);
                                Rscal(n, alpha, G.data(), ix);
                                compare(tag("Rscal", {}, {n, ix, ak}), R, G);
                            }
                            std::vector<REAL> RX(X), RY(Y), GX(X), GY(Y);
                            Rswap_reference_cpu(n, RX.data(), ix, RY.data(), iy);
                            Rswap(n, GX.data(), ix, GY.data(), iy);
                            compare(tag("Rswap x", {}, {n, ix, iy}), RX, GX);
                            compare(tag("Rswap y", {}, {n, ix, iy}), RY, GY);
                            RX = X, RY = Y, GX = X, GY = Y;
                            Rrot_reference_cpu(n, RX.data(), ix, RY.data(), iy, alpha, beta);
                            Rrot(n, GX.data(), ix, GY.data(), iy, alpha, beta);
                            compare(tag("Rrot x", {}, {n, ix, iy, ak, bk}), RX, GX);
                            compare(tag("Rrot y", {}, {n, ix, iy, ak, bk}), RY, GY);
                        }
                    }
            }
        // Rlaswp (LAPACK): row interchanges of every column
        for (long incx : {1L, -1L, 2L}) {
            const long rows = m + 5, k1 = 2, k2 = std::min(rows, k1 + m);
            std::vector<REAL> A((rows + 1) * (n + 33));
            fill(A);
            std::vector<INTEGER> ipiv((size_t)(std::labs(incx) * (k2 + 1)));
            for (auto &v : ipiv)
                v = (INTEGER)(k1 + rng() % (rows - k1 + 1));
            std::vector<REAL> R(A), G(A);
            Rlaswp_reference_cpu(n + 33, R.data(), rows + 1, k1, k2, ipiv.data(), incx);
            Rlaswp(n + 33, G.data(), rows + 1, k1, k2, ipiv.data(), incx);
            compare(tag("Rlaswp", {}, {m, n, incx}), R, G);
        }
    }
#ifdef RBLAS_TEST_GPU
    std::printf("lane BLAS (GPU): %ld checks, %ld failures\n", checks, failures);
#else
    std::printf("lane BLAS (OpenMP): %ld checks, %ld failures\n", checks, failures);
#endif
    return failures == 0 ? 0 : 1;
}
