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
 * dd/td/qd complex Level 2/3 BLAS on the CPU with AVX2: the column and
 * row-block forms of Rblas_simd_common.h for COMPLEX = qd3_complex<REAL>
 * (real part, then imaginary part; libQD3's complex multiplication
 * (ar*br - ai*bi, ar*bi + ai*br) and componentwise addition, on the words
 * of the branch-free arithmetic).  Covered: Cgemm with transa = N, Csyrk and
 * Cherk with trans = N, Ctrmm and Ctrsm with side L and transa = N or with
 * side R, Cgemv with trans = N, Cgeru, Cgerc.  The loops are those of the
 * reference routines, so the results are the reference results bit for bit;
 * the other cases run the reference code in parallel (Rblas_split_common.h).
 */

#ifndef _MPLAPACK_CBLAS_SIMD_COMMON_H_
#define _MPLAPACK_CBLAS_SIMD_COMMON_H_

#include "Rblas_simd_common.h"

#ifdef MPLAPACK_BK_TYPE

namespace mplapack_blas {

#define MPLAPACK_CBK_KERNELS(T, W)                                                                                        \
    /* y = y + t*x, or y = y - t*x (sub); complex t, x, y */                                                               \
    MPLAPACK_GEMM_CLONES static void cbk_axpy_##T(long n, double *__restrict y, const double *t, const double *__restrict x, int sub) \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(tr);                                                                                        \
        MPLAPACK_BK_##T##_DECL(ti);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(tr, t);                                                                                     \
        MPLAPACK_BK_##T##_LOAD(ti, t + W);                                                                                 \
        MPLAPACK_GEMM_SIMD                                                                                                 \
        for (long k = 0; k < n; k++) {                                                                                     \
            MPLAPACK_BK_##T##_DECL(xr);                                                                                    \
            MPLAPACK_BK_##T##_DECL(xi);                                                                                    \
            MPLAPACK_BK_##T##_DECL(yr);                                                                                    \
            MPLAPACK_BK_##T##_DECL(yi);                                                                                    \
            MPLAPACK_BK_##T##_DECL(a);                                                                                     \
            MPLAPACK_BK_##T##_DECL(b);                                                                                     \
            MPLAPACK_BK_##T##_LOAD(xr, x + 2 * W * k);                                                                     \
            MPLAPACK_BK_##T##_LOAD(xi, x + 2 * W * k + W);                                                                 \
            MPLAPACK_BK_##T##_LOAD(yr, y + 2 * W * k);                                                                     \
            MPLAPACK_BK_##T##_LOAD(yi, y + 2 * W * k + W);                                                                 \
            MPLAPACK_BK_##T##_MUL(a, tr, xr); /* re = tr*xr - ti*xi */                                                     \
            MPLAPACK_BK_##T##_MUL(b, ti, xi);                                                                              \
            MPLAPACK_BK_##T##_NEG(b);                                                                                      \
            MPLAPACK_BK_##T##_ADD(a, a, b);                                                                                \
            if (sub)                                                                                                       \
                MPLAPACK_BK_##T##_NEG(a);                                                                                  \
            MPLAPACK_BK_##T##_ADD(yr, yr, a);                                                                              \
            MPLAPACK_BK_##T##_MUL(a, tr, xi); /* im = tr*xi + ti*xr */                                                     \
            MPLAPACK_BK_##T##_MUL(b, ti, xr);                                                                              \
            MPLAPACK_BK_##T##_ADD(a, a, b);                                                                                \
            if (sub)                                                                                                       \
                MPLAPACK_BK_##T##_NEG(a);                                                                                  \
            MPLAPACK_BK_##T##_ADD(yi, yi, a);                                                                              \
            MPLAPACK_BK_##T##_STORE(y + 2 * W * k, yr);                                                                    \
            MPLAPACK_BK_##T##_STORE(y + 2 * W * k + W, yi);                                                                \
        }                                                                                                                  \
    }                                                                                                                      \
    /* p(k) = a(k)*x(k) (cj: conj(a(k))*x(k)), the product of qd3_complex */                                                                    \
    MPLAPACK_GEMM_CLONES static void cbk_vmul_##T(long n, double *__restrict p, const double *__restrict a, const double *__restrict x, int cj) \
    {                                                                                                                                           \
        MPLAPACK_GEMM_SIMD                                                                                                                      \
        for (long k = 0; k < n; k++) {                                                                                                          \
            MPLAPACK_BK_##T##_DECL(ar);                                                                                                         \
            MPLAPACK_BK_##T##_DECL(ai);                                                                                                         \
            MPLAPACK_BK_##T##_DECL(xr);                                                                                                         \
            MPLAPACK_BK_##T##_DECL(xi);                                                                                                         \
            MPLAPACK_BK_##T##_DECL(u);                                                                                                          \
            MPLAPACK_BK_##T##_DECL(v);                                                                                                          \
            MPLAPACK_BK_##T##_LOAD(ar, a + 2 * W * k);                                                                                          \
            MPLAPACK_BK_##T##_LOAD(ai, a + 2 * W * k + W);                                                                                      \
            if (cj)                                                                                                                             \
                MPLAPACK_BK_##T##_NEG(ai);                                                                                                      \
            MPLAPACK_BK_##T##_LOAD(xr, x + 2 * W * k);                                                                                          \
            MPLAPACK_BK_##T##_LOAD(xi, x + 2 * W * k + W);                                                                                      \
            MPLAPACK_BK_##T##_MUL(u, ar, xr); /* re = ar*xr - ai*xi */                                                                          \
            MPLAPACK_BK_##T##_MUL(v, ai, xi);                                                                                                   \
            MPLAPACK_BK_##T##_NEG(v);                                                                                                           \
            MPLAPACK_BK_##T##_ADD(u, u, v);                                                                                                     \
            MPLAPACK_BK_##T##_STORE(p + 2 * W * k, u);                                                                                          \
            MPLAPACK_BK_##T##_MUL(u, ar, xi); /* im = ar*xi + ai*xr */                                                                          \
            MPLAPACK_BK_##T##_MUL(v, ai, xr);                                                                                                   \
            MPLAPACK_BK_##T##_ADD(u, u, v);                                                                                                     \
            MPLAPACK_BK_##T##_STORE(p + 2 * W * k + W, u);                                                                                      \
        }                                                                                                                                       \
    }                                                                                                                                           \
    /* y = t*y; complex t */                                                                                               \
    MPLAPACK_GEMM_CLONES static void cbk_scale_##T(long n, double *__restrict y, const double *t)                          \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(tr);                                                                                        \
        MPLAPACK_BK_##T##_DECL(ti);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(tr, t);                                                                                     \
        MPLAPACK_BK_##T##_LOAD(ti, t + W);                                                                                 \
        MPLAPACK_GEMM_SIMD                                                                                                 \
        for (long k = 0; k < n; k++) {                                                                                     \
            MPLAPACK_BK_##T##_DECL(yr);                                                                                    \
            MPLAPACK_BK_##T##_DECL(yi);                                                                                    \
            MPLAPACK_BK_##T##_DECL(a);                                                                                     \
            MPLAPACK_BK_##T##_DECL(b);                                                                                     \
            MPLAPACK_BK_##T##_DECL(c);                                                                                     \
            MPLAPACK_BK_##T##_LOAD(yr, y + 2 * W * k);                                                                     \
            MPLAPACK_BK_##T##_LOAD(yi, y + 2 * W * k + W);                                                                 \
            MPLAPACK_BK_##T##_MUL(a, tr, yr);                                                                              \
            MPLAPACK_BK_##T##_MUL(b, ti, yi);                                                                              \
            MPLAPACK_BK_##T##_NEG(b);                                                                                      \
            MPLAPACK_BK_##T##_ADD(c, a, b);                                                                                \
            MPLAPACK_BK_##T##_MUL(a, tr, yi);                                                                              \
            MPLAPACK_BK_##T##_MUL(b, ti, yr);                                                                              \
            MPLAPACK_BK_##T##_ADD(a, a, b);                                                                                \
            MPLAPACK_BK_##T##_STORE(y + 2 * W * k, c);                                                                     \
            MPLAPACK_BK_##T##_STORE(y + 2 * W * k + W, a);                                                                 \
        }                                                                                                                  \
    }                                                                                                                      \
    /* y = r*y; real r */                                                                                                  \
    MPLAPACK_GEMM_CLONES static void cbk_rscale_##T(long n, double *__restrict y, const double *r)                         \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(rr);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(rr, r);                                                                                     \
        MPLAPACK_GEMM_SIMD                                                                                                 \
        for (long k = 0; k < 2 * n; k++) {                                                                                 \
            MPLAPACK_BK_##T##_DECL(v);                                                                                     \
            MPLAPACK_BK_##T##_LOAD(v, y + W * k);                                                                          \
            MPLAPACK_BK_##T##_MUL(v, v, rr);                                                                               \
            MPLAPACK_BK_##T##_STORE(y + W * k, v);                                                                         \
        }                                                                                                                  \
    }

#if defined(MPLAPACK_BUILD_WITH_DD)
MPLAPACK_CBK_KERNELS(DD, 2)
#elif defined(MPLAPACK_BUILD_WITH_TD)
MPLAPACK_CBK_KERNELS(TD, 3)
#elif defined(MPLAPACK_BUILD_WITH_QD)
MPLAPACK_CBK_KERNELS(QD, 4)
#endif

#define MPLAPACK_CBLAS_SIMD_AVAILABLE 1

// y(k) = y(k) + t*x(k) (sub: y(k) - t*x(k)), element strides
template <class C> void cbk_axpy(long n, C *y, long sy, const C &t, const C *x, long sx, bool sub)
{
    if (n <= 0)
        return;
    if (sy == 1 && sx == 1) {
        MPLAPACK_BK_FN(cbk_axpy_)(n, bk_w(y), bk_w(&t), bk_w(x), sub ? 1 : 0);
        return;
    }
    for (long k = 0; k < n; k++) {
        if (sub)
            y[k * sy] = y[k * sy] - t * x[k * sx];
        else
            y[k * sy] += t * x[k * sx];
    }
}

// p(k) = a(k)*x(k) (cj: conj(a(k))*x(k)), unit strides
template <class C> void cbk_vmul(long n, C *p, const C *a, const C *x, bool cj)
{
    if (n > 0)
        MPLAPACK_BK_FN(cbk_vmul_)(n, bk_w(p), bk_w(a), bk_w(x), cj ? 1 : 0);
}

template <class C> void cbk_scale(long n, C *y, long sy, const C &t)
{
    if (n <= 0)
        return;
    if (sy == 1) {
        MPLAPACK_BK_FN(cbk_scale_)(n, bk_w(y), bk_w(&t));
        return;
    }
    for (long k = 0; k < n; k++)
        y[k * sy] = t * y[k * sy];
}

template <class C, class R> void cbk_rscale(long n, C *y, const R &r)
{
    if (n > 0)
        MPLAPACK_BK_FN(cbk_rscale_)(n, bk_w(y), bk_w(&r));
}

template <class C> void cbk_zero(long n, C *y, long sy)
{
    for (long k = 0; k < n; k++)
        y[k * sy] = C(0.0, 0.0);
}

// Runs f(u) for the units u < units: in parallel, or on one thread.
template <class F> bool cbk_run(double work, long units, F f)
{
    if (units <= 0 || blas_env_long("MPLAPACK_BLAS_SIMD", 1) == 0 || work < (double)blas_env_long("MPLAPACK_BLAS_OMP_MIN", 32768))
        return false;
#ifdef _OPENMP
    if (blas_env_long("MPLAPACK_BLAS_OMP", 1) != 0 && omp_get_max_threads() >= 2) {
#pragma omp parallel for schedule(dynamic, 1)
        for (long u = 0; u < units; u++)
            f(u);
        return true;
    }
#endif
    for (long u = 0; u < units; u++)
        f(u);
    return true;
}

static const long cbk_rows = 64; // rows per unit of the row-block forms

// ---- Cgemm, transa = N: column j of C ----
template <class C> bool cgemm_simd(bool notb, bool conjb, long m, long n, long k, const C &alpha, const C *a, long lda, const C *b, long ldb, const C &beta, C *c, long ldc)
{
    const C zero(0.0, 0.0), one(1.0, 0.0);
    if (alpha == zero) // the reference only scales
        return false;
    return cbk_run(4.0 * m * n * k, n, [&](long u) {
        const long j = u;
        C *cj = c + j * ldc;
        if (beta == zero)
            cbk_zero(m, cj, 1);
        else if (beta != one)
            cbk_scale(m, cj, 1, beta);
        for (long l = 0; l < k; l++) {
            const C temp = notb ? alpha * b[l + j * ldb] : (conjb ? alpha * conj(b[j + l * ldb]) : alpha * b[j + l * ldb]);
            cbk_axpy(m, cj, 1, temp, a + l * lda, 1, false);
        }
    });
}

// ---- Csyrk / Cherk, trans = N: column j of the triangle of C ----
template <class C> bool csyrk_simd(bool upper, long n, long k, const C &alpha, const C *a, long lda, const C &beta, C *c, long ldc)
{
    const C zero(0.0, 0.0), one(1.0, 0.0);
    if (alpha == zero)
        return false;
    return cbk_run(4.0 * n * n * k, n, [&](long u) {
        const long j = u, r0 = upper ? 0 : j, cnt = upper ? j + 1 : n - j;
        C *cj = c + r0 + j * ldc;
        if (beta == zero)
            cbk_zero(cnt, cj, 1);
        else if (beta != one)
            cbk_scale(cnt, cj, 1, beta);
        for (long l = 0; l < k; l++)
            if (a[j + l * lda] != zero) {
                const C temp = alpha * a[j + l * lda];
                cbk_axpy(cnt, cj, 1, temp, a + r0 + l * lda, 1, false);
            }
    });
}

template <class C, class R> bool cherk_simd(bool upper, long n, long k, const R &alpha, const C *a, long lda, const R &beta, C *c, long ldc)
{
    const R zero = 0.0, one = 1.0;
    if (alpha == zero || (k == 0 && beta == one)) // the reference only scales, or returns
        return false;
    return cbk_run(4.0 * n * n * k, n, [&](long u) {
        const long j = u;
        const long r0 = upper ? 0 : j + 1, cnt = upper ? j : n - j - 1; // off-diagonal rows
        C *col = c + r0 + j * ldc;
        C &d = c[j + j * ldc];
        if (beta == zero) {
            cbk_zero(cnt, col, 1);
            d = zero;
        } else if (beta != one) {
            cbk_rscale(cnt, col, beta);
            d = beta * d.real();
        } else {
            d = d.real();
        }
        for (long l = 0; l < k; l++)
            if (a[j + l * lda] != C(zero)) {
                const C temp = alpha * conj(a[j + l * lda]);
                cbk_axpy(cnt, col, 1, temp, a + r0 + l * lda, 1, false);
                d = d.real() + (temp * a[j + l * lda]).real();
            }
    });
}

// ---- Ctrmm / Ctrsm (B is b) ----
template <class C> bool ctrxm_simd(bool solve, bool left, bool upper, bool notr, bool noconj, bool nounit, long m, long n, const C &alpha, const C *a, long lda, C *b, long ldb)
{
    const C zero(0.0, 0.0), one(1.0, 0.0);
    if (alpha == zero) // the reference sets B to zero
        return false;
    if (left) {
        if (!notr)
            return false;
        return cbk_run(4.0 * m * m * n, n, [&](long u) { // column j
            const long j = u;
            C *bj = b + j * ldb;
            if (solve) {
                if (alpha != one)
                    cbk_scale(m, bj, 1, alpha);
                for (long kk = 0; kk < m; kk++) {
                    const long k = upper ? m - 1 - kk : kk;
                    if (bj[k] != zero) {
                        if (nounit)
                            bj[k] = bj[k] / a[k + k * lda];
                        const C bk = bj[k];
                        if (upper)
                            cbk_axpy(k, bj, 1, bk, a + k * lda, 1, true);
                        else
                            cbk_axpy(m - k - 1, bj + k + 1, 1, bk, a + k + 1 + k * lda, 1, true);
                    }
                }
            } else {
                for (long kk = 0; kk < m; kk++) {
                    const long k = upper ? kk : m - 1 - kk;
                    if (bj[k] != zero) {
                        C temp = alpha * bj[k];
                        if (upper) {
                            cbk_axpy(k, bj, 1, temp, a + k * lda, 1, false);
                            if (nounit)
                                temp = temp * a[k + k * lda];
                            bj[k] = temp;
                        } else {
                            bj[k] = temp;
                            if (nounit)
                                bj[k] = bj[k] * a[k + k * lda];
                            cbk_axpy(m - k - 1, bj + k + 1, 1, temp, a + k + 1 + k * lda, 1, false);
                        }
                    }
                }
            }
        });
    }
    const long units = (m + cbk_rows - 1) / cbk_rows;
    return cbk_run(4.0 * m * n * n, units, [&](long u) { // rows [i0, i0 + cnt)
        const long i0 = u * cbk_rows, cnt = std::min(cbk_rows, m - i0);
        auto B = [&](long j) { return b + i0 + j * ldb; };
        auto op = [&](const C &x) { return noconj ? x : conj(x); };
        if (solve) {
            if (notr) {
                for (long jj = 0; jj < n; jj++) {
                    const long j = upper ? jj : n - 1 - jj;
                    if (alpha != one)
                        cbk_scale(cnt, B(j), 1, alpha);
                    const long k0 = upper ? 0 : j + 1, k1 = upper ? j : n;
                    for (long k = k0; k < k1; k++)
                        if (a[k + j * lda] != zero)
                            cbk_axpy(cnt, B(j), 1, a[k + j * lda], B(k), 1, true);
                    if (nounit) {
                        const C temp = one / a[j + j * lda];
                        cbk_scale(cnt, B(j), 1, temp);
                    }
                }
            } else {
                for (long kk = 0; kk < n; kk++) {
                    const long k = upper ? n - 1 - kk : kk;
                    if (nounit) {
                        const C temp = one / op(a[k + k * lda]);
                        cbk_scale(cnt, B(k), 1, temp);
                    }
                    const long j0 = upper ? 0 : k + 1, j1 = upper ? k : n;
                    for (long j = j0; j < j1; j++)
                        if (a[j + k * lda] != zero) {
                            const C temp = op(a[j + k * lda]);
                            cbk_axpy(cnt, B(j), 1, temp, B(k), 1, true);
                        }
                    if (alpha != one)
                        cbk_scale(cnt, B(k), 1, alpha);
                }
            }
        } else {
            if (notr) {
                for (long jj = 0; jj < n; jj++) {
                    const long j = upper ? n - 1 - jj : jj;
                    C temp = alpha;
                    if (nounit)
                        temp = temp * a[j + j * lda];
                    cbk_scale(cnt, B(j), 1, temp);
                    const long k0 = upper ? 0 : j + 1, k1 = upper ? j : n;
                    for (long k = k0; k < k1; k++)
                        if (a[k + j * lda] != zero) {
                            const C t = alpha * a[k + j * lda];
                            cbk_axpy(cnt, B(j), 1, t, B(k), 1, false);
                        }
                }
            } else {
                for (long kk = 0; kk < n; kk++) {
                    const long k = upper ? kk : n - 1 - kk;
                    const long j0 = upper ? 0 : k + 1, j1 = upper ? k : n;
                    for (long j = j0; j < j1; j++)
                        if (a[j + k * lda] != zero) {
                            const C t = noconj ? alpha * a[j + k * lda] : alpha * conj(a[j + k * lda]);
                            cbk_axpy(cnt, B(j), 1, t, B(k), 1, false);
                        }
                    C temp = alpha;
                    if (nounit)
                        temp = noconj ? temp * a[k + k * lda] : temp * conj(a[k + k * lda]);
                    if (temp != one)
                        cbk_scale(cnt, B(k), 1, temp);
                }
            }
        }
    });
}

// ---- Cgemv, trans = N: rows [i0, i0 + cnt) of y ----
template <class C> bool cgemv_simd(long m, long n, const C &alpha, const C *a, long lda, const C *x, long incx, const C &beta, C *y, long incy)
{
    const C zero(0.0, 0.0), one(1.0, 0.0);
    const long units = (m + cbk_rows - 1) / cbk_rows;
    const long kx = incx > 0 ? 0 : -(n - 1) * incx; // x(1) is x[kx]
    return cbk_run(4.0 * m * n, units, [&](long u) {
        const long i0 = u * cbk_rows, cnt = std::min(cbk_rows, m - i0);
        C *yb = y + i0 * incy;
        if (beta != one) {
            if (beta == zero)
                cbk_zero(cnt, yb, incy);
            else
                cbk_scale(cnt, yb, incy, beta);
        }
        if (alpha == zero)
            return;
        for (long j = 0; j < n; j++) {
            const C temp = alpha * x[kx + j * incx];
            cbk_axpy(cnt, yb, incy, temp, a + i0 + j * lda, 1, false);
        }
    });
}

// ---- Cgeru / Cgerc: column j of A ----
template <class C> bool cger_simd(bool conjy, long m, long n, const C &alpha, const C *x, long incx, const C *y, long incy, C *a, long lda)
{
    const C zero(0.0, 0.0);
    if (alpha == zero)
        return false;
    const long kx = incx > 0 ? 0 : -(m - 1) * incx, ky = incy > 0 ? 0 : -(n - 1) * incy;
    return cbk_run(4.0 * m * n, n, [&](long u) {
        const long j = u;
        const C &yj = y[ky + j * incy];
        if (yj != zero) {
            const C temp = conjy ? alpha * conj(yj) : alpha * yj;
            // a(i,j) += x(i)*temp: the complex product is commutative bit for bit
            cbk_axpy(m, a + j * lda, 1, temp, x + kx, incx, false);
        }
    });
}

} // namespace mplapack_blas

#endif // MPLAPACK_BK_TYPE

#endif
