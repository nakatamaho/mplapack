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
#include <mplapack_arithmetic_params.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include <string>
#include <vector>

// The complex BLAS of the optimized libraries (../blas/C*.cpp) against the
// reference routines, bit for bit, with every call run in parallel.

void Cgemm_reference_cpu(const char *transa, const char *transb, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc);
void Csymm_reference_cpu(const char *side, const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc);
void Chemm_reference_cpu(const char *side, const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc);
void Ctrmm_reference_cpu(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb);
void Ctrsm_reference_cpu(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb);
void Cgemv_reference_cpu(const char *trans, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy);
void Cgeru_reference_cpu(INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda);
void Cgerc_reference_cpu(INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda);
void Caxpy_reference_cpu(INTEGER const n, COMPLEX const za, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy);
void Ccopy_reference_cpu(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy);
void Cscal_reference_cpu(INTEGER const n, COMPLEX const za, COMPLEX *zx, INTEGER const incx);
void CRscal_reference_cpu(INTEGER const n, REAL const da, COMPLEX *zx, INTEGER const incx);
void Cswap_reference_cpu(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy);
void CRrot_reference_cpu(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy, REAL const c, REAL const s);

void Csyrk_reference_cpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX const beta, COMPLEX *c, INTEGER const ldc);
void Cherk_reference_cpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, COMPLEX *a, INTEGER const lda, REAL const beta, COMPLEX *c, INTEGER const ldc);
void Csyr2k_reference_cpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc);
void Cher2k_reference_cpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL const beta, COMPLEX *c, INTEGER const ldc);
void Cgemmtr_reference_cpu(const char *uplo, const char *transa, const char *transb, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc);

void Claswp(INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER const k1, INTEGER const k2, INTEGER *ipiv, INTEGER const incx);
void Claswp_reference_cpu(INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER const k1, INTEGER const k2, INTEGER *ipiv, INTEGER const incx);

void Ctrsv_reference_cpu(const char *uplo, const char *trans, const char *diag, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx);
void Cher_reference_cpu(const char *uplo, INTEGER const n, REAL const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *a, INTEGER const lda);
void Cher2_reference_cpu(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda);
void Chemv_reference_cpu(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy);
void Chbmv_reference_cpu(const char *uplo, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy);
void Chpmv_reference_cpu(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *ap, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy);
void Cgbmv_reference_cpu(const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy);
void Ctrmv_reference_cpu(const char *uplo, const char *trans, const char *diag, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx);
void Ctbmv_reference_cpu(const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx);
void Ctpmv_reference_cpu(const char *uplo, const char *trans, const char *diag, INTEGER const n, COMPLEX *ap, COMPLEX *x, INTEGER const incx);
void Chpr_reference_cpu(const char *uplo, INTEGER const n, REAL const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *ap);
void Chpr2_reference_cpu(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *ap);
void Ctbsv_reference_cpu(const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx);
void Ctpsv_reference_cpu(const char *uplo, const char *trans, const char *diag, INTEGER const n, COMPLEX *ap, COMPLEX *x, INTEGER const incx);

COMPLEX Cdotu_reference_cpu(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy);
COMPLEX Cdotc_reference_cpu(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy);
REAL RCasum_reference_cpu(INTEGER const n, COMPLEX *zx, INTEGER const incx);
REAL RCnrm2_reference_cpu(INTEGER const n, COMPLEX *x, INTEGER const incx);
INTEGER iCamax_reference_cpu(INTEGER const n, COMPLEX *zx, INTEGER const incx);

namespace {

std::mt19937_64 rng(20261023);
double u11() { return std::uniform_real_distribution<double>(-1.0, 1.0)(rng); }
long failures = 0, checks = 0, zero_signs = 0;

REAL rnd()
{
    if (rng() % 9 == 0)
        return REAL(0.0);
    REAL x = REAL(u11());
    for (int t = 1; t < 12; t++)
        x += REAL(u11() * std::ldexp(1.0, -53 * t));
    return x;
}

COMPLEX crnd() { return COMPLEX(rnd(), rnd()); }

void fill(std::vector<COMPLEX> &v)
{
    for (auto &x : v)
        x = crnd();
}

bool same(const REAL &a, const REAL &b)
{
#if defined(MPLAPACK_BUILD_WITH_MPFR)
    if (mpfr_nan_p(a.get_mpfr_t()) || mpfr_nan_p(b.get_mpfr_t()))
        return mpfr_nan_p(a.get_mpfr_t()) && mpfr_nan_p(b.get_mpfr_t());
#ifdef RBLAS_TEST_GPU
    // The MPFR complex GPU routines compute +0 for every zero
    // (mpfr/cuda/Cblas_mpfr_bridge_cuda.cpp): counted, not failures.
    if (mpfr_zero_p(a.get_mpfr_t()) && mpfr_zero_p(b.get_mpfr_t()) && mpfr_signbit(a.get_mpfr_t()) != mpfr_signbit(b.get_mpfr_t())) {
        zero_signs++;
        return mpfr_get_prec(a.get_mpfr_t()) == mpfr_get_prec(b.get_mpfr_t());
    }
#endif
    return mpfr_equal_p(a.get_mpfr_t(), b.get_mpfr_t()) && mpfr_signbit(a.get_mpfr_t()) == mpfr_signbit(b.get_mpfr_t()) && mpfr_get_prec(a.get_mpfr_t()) == mpfr_get_prec(b.get_mpfr_t());
#elif defined(MPLAPACK_BUILD_WITH_GMP)
    return mpf_cmp(a.get_mpf_t(), b.get_mpf_t()) == 0;
#elif defined(MPLAPACK_BUILD_WITH_BINARY80)
    return (a == b && std::signbit((long double)a) == std::signbit((long double)b)) || (a != a && b != b);
#else
    return std::memcmp(&a, &b, sizeof(REAL)) == 0;
#endif
}

bool csame(const COMPLEX &a, const COMPLEX &b) { return same(a.real(), b.real()) && same(a.imag(), b.imag()); }

void compare(const std::string &what, const std::vector<COMPLEX> &ref, const std::vector<COMPLEX> &got)
{
    checks++;
    for (size_t t = 0; t < ref.size(); t++)
        if (!csame(ref[t], got[t])) {
            failures++;
            std::printf("FAIL %s: element %ld differs\n", what.c_str(), (long)t);
            return;
        }
}

COMPLEX scalar(int kind)
{
    if (kind == 0)
        return COMPLEX(REAL(1.0), REAL(0.0));
    if (kind == 1)
        return COMPLEX(REAL(0.0), REAL(0.0));
    return COMPLEX(rnd() + REAL(0.5), rnd());
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
    // the small-block engine (Rblas_block_common.h) for every type and size
    setenv("MPLAPACK_BLAS_BLOCK_MIN", "1", 1);
    setenv("MPLAPACK_BLAS_BLOCK", "2", 1);
#ifdef RBLAS_TEST_GPU
    // Cblas_gpu_test.cpp: every call, Level 1 and 2 included, on the GPU
    rblas_test_gpu_setup();
    if (!rblas_test_have_gpu()) {
        std::printf("no CUDA device: skipped\n");
        return 77;
    }
#endif
    const char *ul[] = {"U", "L"}, *tr[] = {"N", "T", "C"}, *sd[] = {"L", "R"}, *dg[] = {"N", "U"};
    // {150, 70, 40}: several blocks of the engine (alpha, beta general only)
    const long sizes[][3] = {{1, 1, 1}, {7, 5, 3}, {19, 23, 6}, {150, 70, 40}};
    const long incs[] = {1, 2, -2};
    for (auto &z : sizes) {
        const long m = z[0], n = z[1], k = z[2];
        for (int ak = 0; ak < 3; ak++)
            for (int bk = 0; bk < 3; bk++) {
                if (m > 100 && (ak != 2 || bk == 1))
                    continue;
                const COMPLEX alpha = scalar(ak), beta = scalar(bk);
                for (const char *ta : tr)
                    for (const char *tb : tr) {
                        const long ra = ta[0] == 'N' ? m : k, ca = ta[0] == 'N' ? k : m;
                        const long rb = tb[0] == 'N' ? k : n, cb = tb[0] == 'N' ? n : k;
                        std::vector<COMPLEX> A((ra + 1) * ca), B((rb + 2) * cb), C((m + 3) * n);
                        fill(A), fill(B), fill(C);
                        std::vector<COMPLEX> R(C), G(C);
                        Cgemm_reference_cpu(ta, tb, m, n, k, alpha, A.data(), ra + 1, B.data(), rb + 2, beta, R.data(), m + 3);
                        Cgemm(ta, tb, m, n, k, alpha, A.data(), ra + 1, B.data(), rb + 2, beta, G.data(), m + 3);
                        compare(tag("Cgemm", {ta, tb}, {m, n, k, ak, bk}), R, G);
                    }
                for (const char *u : ul)
                    for (const char *t : tr) {
                        const bool nt = t[0] == 'N';
                        const long lda = (nt ? n : k) + 2, ca = nt ? k : n;
                        std::vector<COMPLEX> A(lda * ca), B(lda * ca), C((n + 1) * n);
                        fill(A), fill(B), fill(C);
                        std::vector<COMPLEX> R(C), G(C);
                        if (t[0] != 'C') {
                            Csyrk_reference_cpu(u, t, n, k, alpha, A.data(), lda, beta, R.data(), n + 1);
                            Csyrk(u, t, n, k, alpha, A.data(), lda, beta, G.data(), n + 1);
                            compare(tag("Csyrk", {u, t}, {n, k, ak, bk}), R, G);
                            R = C, G = C;
                            Csyr2k_reference_cpu(u, t, n, k, alpha, A.data(), lda, B.data(), lda, beta, R.data(), n + 1);
                            Csyr2k(u, t, n, k, alpha, A.data(), lda, B.data(), lda, beta, G.data(), n + 1);
                            compare(tag("Csyr2k", {u, t}, {n, k, ak, bk}), R, G);
                        }
                        if (t[0] != 'T') {
                            R = C, G = C;
                            Cherk_reference_cpu(u, t, n, k, alpha.real(), A.data(), lda, beta.real(), R.data(), n + 1);
                            Cherk(u, t, n, k, alpha.real(), A.data(), lda, beta.real(), G.data(), n + 1);
                            compare(tag("Cherk", {u, t}, {n, k, ak, bk}), R, G);
                            R = C, G = C;
                            Cher2k_reference_cpu(u, t, n, k, alpha, A.data(), lda, B.data(), lda, beta.real(), R.data(), n + 1);
                            Cher2k(u, t, n, k, alpha, A.data(), lda, B.data(), lda, beta.real(), G.data(), n + 1);
                            compare(tag("Cher2k", {u, t}, {n, k, ak, bk}), R, G);
                        }
                        for (const char *t2 : tr) {
                            const long rb = t2[0] == 'N' ? k : n, cb = t2[0] == 'N' ? n : k;
                            std::vector<COMPLEX> B2((rb + 1) * cb);
                            fill(B2);
                            R = C, G = C;
                            Cgemmtr_reference_cpu(u, t, t2, n, k, alpha, A.data(), lda, B2.data(), rb + 1, beta, R.data(), n + 1);
                            Cgemmtr(u, t, t2, n, k, alpha, A.data(), lda, B2.data(), rb + 1, beta, G.data(), n + 1);
                            compare(tag("Cgemmtr", {u, t, t2}, {n, k, ak, bk}), R, G);
                        }
                    }
                for (const char *s : sd)
                    for (const char *u : ul) {
                        const long na = s[0] == 'L' ? m : n, lda = na + 1;
                        std::vector<COMPLEX> A(lda * na), B((m + 2) * n), C((m + 3) * n);
                        fill(A), fill(B), fill(C);
                        std::vector<COMPLEX> R(C), G(C);
                        Csymm_reference_cpu(s, u, m, n, alpha, A.data(), lda, B.data(), m + 2, beta, R.data(), m + 3);
                        Csymm(s, u, m, n, alpha, A.data(), lda, B.data(), m + 2, beta, G.data(), m + 3);
                        compare(tag("Csymm", {s, u}, {m, n, ak, bk}), R, G);
                        R = C, G = C;
                        Chemm_reference_cpu(s, u, m, n, alpha, A.data(), lda, B.data(), m + 2, beta, R.data(), m + 3);
                        Chemm(s, u, m, n, alpha, A.data(), lda, B.data(), m + 2, beta, G.data(), m + 3);
                        compare(tag("Chemm", {s, u}, {m, n, ak, bk}), R, G);
                        if (bk != 0)
                            continue;
                        for (const char *t : tr)
                            for (const char *d : dg) {
                                std::vector<COMPLEX> T(A);
                                for (long i = 0; i < na; i++)
                                    T[i + i * lda] = T[i + i * lda] + COMPLEX(REAL(2.0 + na), REAL(0.0));
                                R = B, G = B;
                                Ctrmm_reference_cpu(s, u, t, d, m, n, alpha, T.data(), lda, R.data(), m + 2);
                                Ctrmm(s, u, t, d, m, n, alpha, T.data(), lda, G.data(), m + 2);
                                compare(tag("Ctrmm", {s, u, t, d}, {m, n, ak}), R, G);
                                R = B, G = B;
                                Ctrsm_reference_cpu(s, u, t, d, m, n, alpha, T.data(), lda, R.data(), m + 2);
                                Ctrsm(s, u, t, d, m, n, alpha, T.data(), lda, G.data(), m + 2);
                                compare(tag("Ctrsm", {s, u, t, d}, {m, n, ak}), R, G);
                            }
                    }
                for (long ix : incs)
                    for (long iy : incs) {
                        const long lx = std::labs(ix), ly = std::labs(iy), mx = std::max(m, n);
                        std::vector<COMPLEX> A((m + 2) * mx), X(lx * mx + 1), Y(ly * mx + 1);
                        fill(A), fill(X), fill(Y);
                        for (const char *t : tr) {
                            std::vector<COMPLEX> R(Y), G(Y);
                            Cgemv_reference_cpu(t, m, n, alpha, A.data(), m + 2, X.data(), ix, beta, R.data(), iy);
                            Cgemv(t, m, n, alpha, A.data(), m + 2, X.data(), ix, beta, G.data(), iy);
                            compare(tag("Cgemv", {t}, {m, n, ix, iy, ak, bk}), R, G);
                        }
                        const long nn = std::min(m, n);
                        for (const char *u : ul) {
                            std::vector<COMPLEX> R(A), G(A);
                            Cher_reference_cpu(u, nn, alpha.real(), X.data(), ix, R.data(), m + 2);
                            Cher(u, nn, alpha.real(), X.data(), ix, G.data(), m + 2);
                            compare(tag("Cher", {u}, {nn, ix, ak}), R, G);
                            R = A, G = A;
                            Cher2_reference_cpu(u, nn, alpha, X.data(), ix, Y.data(), iy, R.data(), m + 2);
                            Cher2(u, nn, alpha, X.data(), ix, Y.data(), iy, G.data(), m + 2);
                            compare(tag("Cher2", {u}, {nn, ix, iy, ak}), R, G);
                            for (const char *t : tr)
                                for (const char *d : dg) {
                                    std::vector<COMPLEX> T(A);
                                    for (long q = 0; q < nn; q++)
                                        T[q + q * (m + 2)] = T[q + q * (m + 2)] + COMPLEX(REAL(2.0 + nn), REAL(0.0));
                                    std::vector<COMPLEX> RX(X), GX(X);
                                    Ctrsv_reference_cpu(u, t, d, nn, T.data(), m + 2, RX.data(), ix);
                                    Ctrsv(u, t, d, nn, T.data(), m + 2, GX.data(), ix);
                                    compare(tag("Ctrsv", {u, t, d}, {nn, ix}), RX, GX);
                                }
                        }
                        // Hermitian, band and packed Level 2 (Chemv ... Chpr2)
                        for (const char *u : ul) {
                            std::vector<COMPLEX> RY(Y), GY(Y);
                            Chemv_reference_cpu(u, nn, alpha, A.data(), m + 2, X.data(), ix, beta, RY.data(), iy);
                            Chemv(u, nn, alpha, A.data(), m + 2, X.data(), ix, beta, GY.data(), iy);
                            compare(tag("Chemv", {u}, {nn, ix, iy, ak, bk}), RY, GY);
                            std::vector<COMPLEX> AP((size_t)(nn * (nn + 1) / 2 + 1));
                            fill(AP);
                            RY = Y, GY = Y;
                            Chpmv_reference_cpu(u, nn, alpha, AP.data(), X.data(), ix, beta, RY.data(), iy);
                            Chpmv(u, nn, alpha, AP.data(), X.data(), ix, beta, GY.data(), iy);
                            compare(tag("Chpmv", {u}, {nn, ix, iy, ak, bk}), RY, GY);
                            std::vector<COMPLEX> RP(AP), GP(AP);
                            Chpr_reference_cpu(u, nn, alpha.real(), X.data(), ix, RP.data());
                            Chpr(u, nn, alpha.real(), X.data(), ix, GP.data());
                            compare(tag("Chpr", {u}, {nn, ix, ak}), RP, GP);
                            RP = AP, GP = AP;
                            Chpr2_reference_cpu(u, nn, alpha, X.data(), ix, Y.data(), iy, RP.data());
                            Chpr2(u, nn, alpha, X.data(), ix, Y.data(), iy, GP.data());
                            compare(tag("Chpr2", {u}, {nn, ix, iy, ak}), RP, GP);
                            for (long kb : {0L, 1L, 3L}) {
                                std::vector<COMPLEX> AB((size_t)((kb + 2) * mx));
                                fill(AB);
                                RY = Y, GY = Y;
                                Chbmv_reference_cpu(u, nn, kb, alpha, AB.data(), kb + 2, X.data(), ix, beta, RY.data(), iy);
                                Chbmv(u, nn, kb, alpha, AB.data(), kb + 2, X.data(), ix, beta, GY.data(), iy);
                                compare(tag("Chbmv", {u}, {nn, kb, ix, iy, ak, bk}), RY, GY);
                                for (const char *t : tr)
                                    for (const char *d : dg) {
                                        std::vector<COMPLEX> RX(X), GX(X);
                                        Ctbmv_reference_cpu(u, t, d, nn, kb, AB.data(), kb + 2, RX.data(), ix);
                                        Ctbmv(u, t, d, nn, kb, AB.data(), kb + 2, GX.data(), ix);
                                        compare(tag("Ctbmv", {u, t, d}, {nn, kb, ix}), RX, GX);
                                        std::vector<COMPLEX> TB(AB);
                                        for (long q = 0; q < nn; q++)
                                            TB[(u[0] == 'U' ? kb : 0) + q * (kb + 2)] = TB[(u[0] == 'U' ? kb : 0) + q * (kb + 2)] + COMPLEX(REAL(2.0 + kb), REAL(0.0));
                                        RX = X, GX = X;
                                        Ctbsv_reference_cpu(u, t, d, nn, kb, TB.data(), kb + 2, RX.data(), ix);
                                        Ctbsv(u, t, d, nn, kb, TB.data(), kb + 2, GX.data(), ix);
                                        compare(tag("Ctbsv", {u, t, d}, {nn, kb, ix}), RX, GX);
                                    }
                            }
                            for (const char *t : tr)
                                for (const char *d : dg) {
                                    std::vector<COMPLEX> RX(X), GX(X);
                                    Ctrmv_reference_cpu(u, t, d, nn, A.data(), m + 2, RX.data(), ix);
                                    Ctrmv(u, t, d, nn, A.data(), m + 2, GX.data(), ix);
                                    compare(tag("Ctrmv", {u, t, d}, {nn, ix}), RX, GX);
                                    RX = X, GX = X;
                                    Ctpmv_reference_cpu(u, t, d, nn, AP.data(), RX.data(), ix);
                                    Ctpmv(u, t, d, nn, AP.data(), GX.data(), ix);
                                    compare(tag("Ctpmv", {u, t, d}, {nn, ix}), RX, GX);
                                    std::vector<COMPLEX> TP(AP);
                                    for (long q = 0; q < nn; q++) {
                                        const long s = u[0] == 'U' ? q * (q + 1) / 2 + q : q * nn - q * (q - 1) / 2;
                                        TP[s] = TP[s] + COMPLEX(REAL(2.0 + nn), REAL(0.0));
                                    }
                                    RX = X, GX = X;
                                    Ctpsv_reference_cpu(u, t, d, nn, TP.data(), RX.data(), ix);
                                    Ctpsv(u, t, d, nn, TP.data(), GX.data(), ix);
                                    compare(tag("Ctpsv", {u, t, d}, {nn, ix}), RX, GX);
                                }
                        }
                        for (long kl : {0L, 2L})
                            for (long ku : {0L, 1L, 4L}) {
                                std::vector<COMPLEX> AB((size_t)((kl + ku + 2) * mx));
                                fill(AB);
                                for (const char *t : tr) {
                                    std::vector<COMPLEX> RY(Y), GY(Y);
                                    Cgbmv_reference_cpu(t, m, n, kl, ku, alpha, AB.data(), kl + ku + 2, X.data(), ix, beta, RY.data(), iy);
                                    Cgbmv(t, m, n, kl, ku, alpha, AB.data(), kl + ku + 2, X.data(), ix, beta, GY.data(), iy);
                                    compare(tag("Cgbmv", {t}, {m, n, kl, ku, ix, iy, ak, bk}), RY, GY);
                                }
                            }
                        {
                            std::vector<COMPLEX> R(A), G(A);
                            Cgeru_reference_cpu(m, n, alpha, X.data(), ix, Y.data(), iy, R.data(), m + 2);
                            Cgeru(m, n, alpha, X.data(), ix, Y.data(), iy, G.data(), m + 2);
                            compare(tag("Cgeru", {}, {m, n, ix, iy, ak}), R, G);
                            R = A, G = A;
                            Cgerc_reference_cpu(m, n, alpha, X.data(), ix, Y.data(), iy, R.data(), m + 2);
                            Cgerc(m, n, alpha, X.data(), ix, Y.data(), iy, G.data(), m + 2);
                            compare(tag("Cgerc", {}, {m, n, ix, iy, ak}), R, G);
                        }
                        {
                            std::vector<COMPLEX> RX(X), RY(Y), GX(X), GY(Y);
                            Caxpy_reference_cpu(n, alpha, RX.data(), ix, RY.data(), iy);
                            Caxpy(n, alpha, GX.data(), ix, GY.data(), iy);
                            compare(tag("Caxpy", {}, {n, ix, iy, ak}), RY, GY);
                            RX = X, RY = Y, GX = X, GY = Y;
                            Ccopy_reference_cpu(n, RX.data(), ix, RY.data(), iy);
                            Ccopy(n, GX.data(), ix, GY.data(), iy);
                            compare(tag("Ccopy", {}, {n, ix, iy}), RY, GY);
                            RX = X, RY = Y, GX = X, GY = Y;
                            Cswap_reference_cpu(n, RX.data(), ix, RY.data(), iy);
                            Cswap(n, GX.data(), ix, GY.data(), iy);
                            compare(tag("Cswap x", {}, {n, ix, iy}), RX, GX);
                            compare(tag("Cswap y", {}, {n, ix, iy}), RY, GY);
                            RX = X, RY = Y, GX = X, GY = Y;
                            CRrot_reference_cpu(n, RX.data(), ix, RY.data(), iy, alpha.real(), beta.real());
                            CRrot(n, GX.data(), ix, GY.data(), iy, alpha.real(), beta.real());
                            compare(tag("CRrot x", {}, {n, ix, iy, ak, bk}), RX, GX);
                            compare(tag("CRrot y", {}, {n, ix, iy, ak, bk}), RY, GY);
                            if (ix > 0) {
                                RX = X, GX = X;
                                Cscal_reference_cpu(n, alpha, RX.data(), ix);
                                Cscal(n, alpha, GX.data(), ix);
                                compare(tag("Cscal", {}, {n, ix, ak}), RX, GX);
                                RX = X, GX = X;
                                CRscal_reference_cpu(n, alpha.real(), RX.data(), ix);
                                CRscal(n, alpha.real(), GX.data(), ix);
                                compare(tag("CRscal", {}, {n, ix, ak}), RX, GX);
                            }
                        }
                    }
            }
        // Claswp (LAPACK): row interchanges of every column
        for (long incx : {1L, -1L, 2L}) {
            const long rows = m + 5, k1 = 2, k2 = std::min(rows, k1 + m);
            std::vector<COMPLEX> A((rows + 1) * (n + 33));
            fill(A);
            std::vector<INTEGER> ipiv((size_t)(std::labs(incx) * (k2 + 1)));
            for (auto &v : ipiv)
                v = (INTEGER)(k1 + rng() % (rows - k1 + 1));
            std::vector<COMPLEX> R(A), G(A);
            Claswp_reference_cpu(n + 33, R.data(), rows + 1, k1, k2, ipiv.data(), incx);
            Claswp(n + 33, G.data(), rows + 1, k1, k2, ipiv.data(), incx);
            compare(tag("Claswp", {}, {m, n, incx}), R, G);
        }
    }
    // Level 1 reductions (Cdotu, Cdotc, RCasum, RCnrm2, iCamax)
    {
        const auto &ap = mplapack::get_arithmetic_params<REAL>();
        const auto bp = mplapack::make_blue_scaling_params(ap);
        for (long len : {1L, 5L, 12L, 23L, 97L, 1000L})
            for (int kind = 0; kind < 3; kind++) {
                std::vector<COMPLEX> X((size_t)len), Y((size_t)len);
                fill(X), fill(Y);
                if (kind == 1)
                    for (long i = 0; i < len; i += 7)
                        X[(size_t)i] = COMPLEX((i % 2 ? bp.tbig * REAL(4.0) : bp.tsml * REAL(0.25)) * X[(size_t)i].real(), X[(size_t)i].imag());
                if (kind == 2)
                    for (long i = 1; i < len; i += 3)
                        X[(size_t)i] = COMPLEX(REAL(i % 2 ? 2.0 : -2.0), REAL(1.0));
                std::vector<COMPLEX> R(1), G(1);
                R[0] = Cdotu_reference_cpu(len, X.data(), 1, Y.data(), 1), G[0] = Cdotu(len, X.data(), 1, Y.data(), 1);
                compare(tag("Cdotu", {}, {len, kind}), R, G);
                R[0] = Cdotc_reference_cpu(len, X.data(), 1, Y.data(), 1), G[0] = Cdotc(len, X.data(), 1, Y.data(), 1);
                compare(tag("Cdotc", {}, {len, kind}), R, G);
                R[0] = COMPLEX(RCasum_reference_cpu(len, X.data(), 1), REAL(0.0)), G[0] = COMPLEX(RCasum(len, X.data(), 1), REAL(0.0));
                compare(tag("RCasum", {}, {len, kind}), R, G);
                R[0] = COMPLEX(RCnrm2_reference_cpu(len, X.data(), 1), REAL(0.0)), G[0] = COMPLEX(RCnrm2(len, X.data(), 1), REAL(0.0));
                compare(tag("RCnrm2", {}, {len, kind}), R, G);
                checks++;
                if (iCamax_reference_cpu(len, X.data(), 1) != iCamax(len, X.data(), 1)) {
                    failures++;
                    std::printf("FAIL iCamax %ld %d\n", len, kind);
                }
            }
    }

#if defined(MPLAPACK_BUILD_WITH_DOUBLE) && !defined(RBLAS_TEST_GPU)
    // complex products whose parts are both NaN: the C99 recovery
    // of the operator, which the kernel of std::complex<double> reproduces
    for (const char *ta : tr) {
        const long m = 37, n = 21, k = 19;
        std::vector<COMPLEX> A(m * k + k * m), B(k * n + n * k), C(m * n);
        fill(A), fill(B), fill(C);
        // (inf, nan)*(x, y): both parts NaN, the recovery gives infinities
        A[3] = COMPLEX(INFINITY, NAN), A[40] = COMPLEX(NAN, -INFINITY);
        for (size_t t = 0; t < B.size(); t += 5)
            B[t] = COMPLEX(0.0, 0.0);
        std::vector<COMPLEX> R(C), G(C);
        const long lda = ta[0] == 'N' ? m : k;
        const COMPLEX alpha(1.5, -0.5), beta(0.5, 0.25);
        Cgemm_reference_cpu(ta, "N", m, n, k, alpha, A.data(), lda, B.data(), k, beta, R.data(), m);
        Cgemm(ta, "N", m, n, k, alpha, A.data(), lda, B.data(), k, beta, G.data(), m);
        compare(tag("Cgemm nan", {ta}, {m, n, k}), R, G);
    }
#endif
#ifdef RBLAS_TEST_GPU
    std::printf("complex BLAS (GPU): %ld checks, %ld failures, %ld zeros differing only in sign\n", checks, failures, zero_signs);
#else
    std::printf("complex BLAS (OpenMP): %ld checks, %ld failures\n", checks, failures);
#endif
    return failures == 0 ? 0 : 1;
}
