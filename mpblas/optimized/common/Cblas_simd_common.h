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
 * of the branch-free arithmetic).  Covered: Cgemm (transa = T, C: over a
 * transposed copy), Csyrk and
 * Cherk (trans = N: columns; trans = T, C: over a transposed copy), Ctrmm and
 * Ctrsm with side L (transa = N: columns; T, C: blocks of columns of a
 * transposed copy) or with
 * side R, Cgemv with trans = N, Cgeru, Cgerc.  The column forms take groups
 * of columns on packed planes (bk_planes).  The loops are those of the
 * reference routines, so the results are the reference results bit for bit;
 * the other cases run the reference code in parallel (Rblas_split_common.h).
 */

#ifndef _MPLAPACK_CBLAS_SIMD_COMMON_H_
#define _MPLAPACK_CBLAS_SIMD_COMMON_H_

#include <vector>
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

// cbk_axpy on planes (bk_planes): word q of the real part of element k at
// y[q * ys + k], of the imaginary part at y[(W + q) * ys + k].
#define MPLAPACK_CBK_PLANE_KERNELS(T, W)                                                                                  \
    MPLAPACK_GEMM_CLONES static void cbk_axpyp_##T(long n, double *y, long ys, const double *t, const double *x, long xs, int sub) \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(tr);                                                                                        \
        MPLAPACK_BK_##T##_DECL(ti);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(tr, t);                                                                                     \
        MPLAPACK_BK_##T##_LOAD(ti, t + W);                                                                                 \
        MPLAPACK_BK_##T##_PTRS(YR, y, ys);                                                                                 \
        MPLAPACK_BK_##T##_PTRS(YI, y + W * ys, ys);                                                                        \
        MPLAPACK_BK_##T##_CPTRS(XR, x, xs);                                                                                \
        MPLAPACK_BK_##T##_CPTRS(XI, x + W * xs, xs);                                                                       \
        MPLAPACK_GEMM_SIMD                                                                                                 \
        for (long k = 0; k < n; k++) {                                                                                     \
            MPLAPACK_BK_##T##_DECL(xr);                                                                                    \
            MPLAPACK_BK_##T##_DECL(xi);                                                                                    \
            MPLAPACK_BK_##T##_DECL(yr);                                                                                    \
            MPLAPACK_BK_##T##_DECL(yi);                                                                                    \
            MPLAPACK_BK_##T##_DECL(a);                                                                                     \
            MPLAPACK_BK_##T##_DECL(b);                                                                                     \
            MPLAPACK_BK_##T##_LOADI(xr, XR, k);                                                                            \
            MPLAPACK_BK_##T##_LOADI(xi, XI, k);                                                                            \
            MPLAPACK_BK_##T##_LOADI(yr, YR, k);                                                                            \
            MPLAPACK_BK_##T##_LOADI(yi, YI, k);                                                                            \
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
            MPLAPACK_BK_##T##_STOREI(YR, k, yr);                                                                           \
            MPLAPACK_BK_##T##_STOREI(YI, k, yi);                                                                           \
        }                                                                                                                  \
    }

#if defined(MPLAPACK_BUILD_WITH_DD)
MPLAPACK_CBK_KERNELS(DD, 2)
MPLAPACK_CBK_PLANE_KERNELS(DD, 2)
#elif defined(MPLAPACK_BUILD_WITH_TD)
MPLAPACK_CBK_KERNELS(TD, 3)
MPLAPACK_CBK_PLANE_KERNELS(TD, 3)
#elif defined(MPLAPACK_BUILD_WITH_QD)
MPLAPACK_CBK_KERNELS(QD, 4)
MPLAPACK_CBK_PLANE_KERNELS(QD, 4)
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

// Whether cbk_run runs the units.
inline bool cbk_go(double work, long units) { return units > 0 && blas_env_long("MPLAPACK_BLAS_SIMD", 1) != 0 && work >= (double)blas_env_long("MPLAPACK_BLAS_OMP_MIN", 32768); }

// Runs f(u) for the units u < units: in parallel, or on one thread.
template <class F> bool cbk_run(double work, long units, F f)
{
    if (!cbk_go(work, units))
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

// Columns per unit of the column forms: the columns of a group are
// independent and each goes through its own operations in the reference
// order, the steps of the group interleaved so that a column of A is used
// by all of them while it is in the cache.
inline long cbk_cols() { return std::max(1L, blas_env_long("MPLAPACK_BLAS_SIMD_COLS", 8)); }

// y(k) = y(k) + t*x(k) (sub: y(k) - t*x(k)) on planes (bk_planes)
template <class C> void cbk_axpyp(long n, double *y, long ys, const C &t, const double *x, long xs, bool sub)
{
    if (n > 0)
        MPLAPACK_BK_FN(cbk_axpyp_)(n, y, ys, bk_w(&t), x, xs, sub ? 1 : 0);
}

// P (rows x cols) := f(i, l)
template <class C, class F> void cbk_pack(bk_planes<C> &P, long rows, long cols, F f)
{
    P.resize(rows, cols);
    for (long l = 0; l < cols; l++)
        for (long i = 0; i < rows; i++)
            P.set(i, l, f(i, l));
}

// The column forms below take groups of columns (cbk_cols); the matrices
// read by the loops over the rows are packed into planes once (bk_planes),
// and a unit packs its columns of the result.

// ---- Cgemm, transa = N: columns [g0, g1) of C ----
template <class C> bool cgemm_simd(bool notb, bool conjb, long m, long n, long k, const C &alpha, const C *a, long lda, const C *b, long ldb, const C &beta, C *c, long ldc)
{
    const C zero(0.0, 0.0), one(1.0, 0.0);
    const long cb = cbk_cols(), units = (n + cb - 1) / cb;
    if (alpha == zero || !cbk_go(4.0 * m * n * k, units)) // alpha == 0: the reference only scales
        return false;
    bk_planes<C> Ap;
    cbk_pack(Ap, m, k, [&](long i, long l) { return a[i + l * lda]; });
    return cbk_run(4.0 * m * n * k, units, [&](long u) {
        const long g0 = u * cb, g1 = std::min(g0 + cb, n);
        bk_planes<C> Cp;
        Cp.resize(m, g1 - g0);
        for (long j = g0; j < g1; j++) {
            C *cj = c + j * ldc;
            if (beta == zero)
                cbk_zero(m, cj, 1);
            else if (beta != one)
                cbk_scale(m, cj, 1, beta);
            for (long i = 0; i < m; i++)
                Cp.set(i, j - g0, cj[i]);
        }
        for (long l = 0; l < k; l++)
            for (long j = g0; j < g1; j++) {
                const C temp = notb ? alpha * b[l + j * ldb] : (conjb ? alpha * conj(b[j + l * ldb]) : alpha * b[j + l * ldb]);
                cbk_axpyp(m, Cp.at(0, j - g0), Cp.ps, temp, Ap.at(0, l), Ap.ps, false);
            }
        for (long j = g0; j < g1; j++)
            for (long i = 0; i < m; i++)
                c[i + j * ldc] = Cp.get(i, j - g0);
    });
}

// ---- Cgemm, transa = T or C: columns [g0, g1) of C, the dot products of
// the rows accumulated at once over the transposed op(A) (m x k) ----
template <class C> bool cgemm_t_simd(bool conja, bool notb, bool conjb, long m, long n, long k, const C &alpha, const C *a, long lda, const C *b, long ldb, const C &beta, C *c, long ldc)
{
    const C zero(0.0, 0.0);
    const long cb = cbk_cols(), units = (n + cb - 1) / cb;
    if (alpha == zero || !cbk_go(4.0 * m * n * k, units)) // alpha == 0: the reference only scales
        return false;
    bk_planes<C> At;
    cbk_pack(At, m, k, [&](long i, long l) { return conja ? conj(a[l + i * lda]) : a[l + i * lda]; });
    return cbk_run(4.0 * m * n * k, units, [&](long u) {
        const long g0 = u * cb, g1 = std::min(g0 + cb, n);
        bk_planes<C> tv;
        cbk_pack(tv, m, g1 - g0, [&](long, long) { return zero; }); // temp = zero
        for (long l = 0; l < k; l++)
            for (long j = g0; j < g1; j++) {
                // temp += op(a)(l,i)*op(b)(l,j)
                const C t = notb ? b[l + j * ldb] : (conjb ? conj(b[j + l * ldb]) : b[j + l * ldb]);
                cbk_axpyp(m, tv.at(0, j - g0), tv.ps, t, At.at(0, l), At.ps, false);
            }
        for (long j = g0; j < g1; j++)
            for (long i = 0; i < m; i++) {
                C &cij = c[i + j * ldc];
                const C t = tv.get(i, j - g0);
                if (beta == zero)
                    cij = alpha * t;
                else
                    cij = alpha * t + beta * cij;
            }
    });
}

// ---- Csyrk / Cherk, trans = N: columns [g0, g1) of the triangle of C ----
template <class C> bool csyrk_simd(bool upper, long n, long k, const C &alpha, const C *a, long lda, const C &beta, C *c, long ldc)
{
    const C zero(0.0, 0.0), one(1.0, 0.0);
    const long cb = cbk_cols(), units = (n + cb - 1) / cb;
    if (alpha == zero || !cbk_go(4.0 * n * n * k, units))
        return false;
    bk_planes<C> Ap;
    cbk_pack(Ap, n, k, [&](long i, long l) { return a[i + l * lda]; });
    return cbk_run(4.0 * n * n * k, units, [&](long u) {
        const long g0 = u * cb, g1 = std::min(g0 + cb, n);
        bk_planes<C> Cp;
        Cp.resize(n, g1 - g0);
        for (long j = g0; j < g1; j++) {
            const long r0 = upper ? 0 : j, cnt = upper ? j + 1 : n - j;
            C *cj = c + r0 + j * ldc;
            if (beta == zero)
                cbk_zero(cnt, cj, 1);
            else if (beta != one)
                cbk_scale(cnt, cj, 1, beta);
            for (long t = 0; t < cnt; t++)
                Cp.set(r0 + t, j - g0, cj[t]);
        }
        for (long l = 0; l < k; l++)
            for (long j = g0; j < g1; j++)
                if (a[j + l * lda] != zero) {
                    const long r0 = upper ? 0 : j, cnt = upper ? j + 1 : n - j;
                    const C temp = alpha * a[j + l * lda];
                    cbk_axpyp(cnt, Cp.at(r0, j - g0), Cp.ps, temp, Ap.at(r0, l), Ap.ps, false);
                }
        for (long j = g0; j < g1; j++) {
            const long r0 = upper ? 0 : j, cnt = upper ? j + 1 : n - j;
            for (long t = 0; t < cnt; t++)
                c[r0 + t + j * ldc] = Cp.get(r0 + t, j - g0);
        }
    });
}

template <class C, class R> bool cherk_simd(bool upper, long n, long k, const R &alpha, const C *a, long lda, const R &beta, C *c, long ldc)
{
    const R zero = 0.0, one = 1.0;
    const long cb = cbk_cols(), units = (n + cb - 1) / cb;
    if (alpha == zero || (k == 0 && beta == one) || !cbk_go(4.0 * n * n * k, units)) // the reference only scales, or returns
        return false;
    bk_planes<C> Ap;
    cbk_pack(Ap, n, k, [&](long i, long l) { return a[i + l * lda]; });
    return cbk_run(4.0 * n * n * k, units, [&](long u) {
        const long g0 = u * cb, g1 = std::min(g0 + cb, n);
        bk_planes<C> Cp; // the off-diagonal rows; the diagonal stays in c
        Cp.resize(n, g1 - g0);
        for (long j = g0; j < g1; j++) {
            const long r0 = upper ? 0 : j + 1, cnt = upper ? j : n - j - 1;
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
            for (long t = 0; t < cnt; t++)
                Cp.set(r0 + t, j - g0, col[t]);
        }
        for (long l = 0; l < k; l++)
            for (long j = g0; j < g1; j++)
                if (a[j + l * lda] != C(zero)) {
                    const long r0 = upper ? 0 : j + 1, cnt = upper ? j : n - j - 1;
                    C &d = c[j + j * ldc];
                    const C temp = alpha * conj(a[j + l * lda]);
                    cbk_axpyp(cnt, Cp.at(r0, j - g0), Cp.ps, temp, Ap.at(r0, l), Ap.ps, false);
                    d = d.real() + (temp * a[j + l * lda]).real();
                }
        for (long j = g0; j < g1; j++) {
            const long r0 = upper ? 0 : j + 1, cnt = upper ? j : n - j - 1;
            for (long t = 0; t < cnt; t++)
                c[r0 + t + j * ldc] = Cp.get(r0 + t, j - g0);
        }
    });
}

// ---- Csyrk trans = T / Cherk trans = C: columns [g0, g1), the dot
// products of the rows accumulated at once over the transposed A (Cherk:
// conj(A)), n x k ----
template <class C> bool csyrk_t_simd(bool upper, long n, long k, const C &alpha, const C *a, long lda, const C &beta, C *c, long ldc)
{
    const C zero(0.0, 0.0);
    const long cb = cbk_cols(), units = (n + cb - 1) / cb;
    if (alpha == zero || !cbk_go(4.0 * n * n * k, units)) // alpha == 0: the reference only scales
        return false;
    bk_planes<C> At;
    cbk_pack(At, n, k, [&](long i, long l) { return a[l + i * lda]; });
    return cbk_run(4.0 * n * n * k, units, [&](long u) {
        const long g0 = u * cb, g1 = std::min(g0 + cb, n);
        bk_planes<C> tv;
        cbk_pack(tv, n, g1 - g0, [&](long, long) { return zero; }); // temp = zero
        for (long l = 0; l < k; l++)
            for (long j = g0; j < g1; j++) {
                const long r0 = upper ? 0 : j, cnt = upper ? j + 1 : n - j;
                cbk_axpyp(cnt, tv.at(r0, j - g0), tv.ps, a[l + j * lda], At.at(r0, l), At.ps, false); // temp += a(l,i)*a(l,j)
            }
        for (long j = g0; j < g1; j++) {
            const long r0 = upper ? 0 : j, cnt = upper ? j + 1 : n - j;
            for (long t = 0; t < cnt; t++) {
                C &cij = c[(r0 + t) + j * ldc];
                const C tt = tv.get(r0 + t, j - g0);
                if (beta == zero)
                    cij = alpha * tt;
                else
                    cij = alpha * tt + beta * cij;
            }
        }
    });
}

template <class C, class R> bool cherk_c_simd(bool upper, long n, long k, const R &alpha, const C *a, long lda, const R &beta, C *c, long ldc)
{
    const R zero = 0.0, one = 1.0;
    const long cb = cbk_cols(), units = (n + cb - 1) / cb;
    if (alpha == zero || (k == 0 && beta == one) || !cbk_go(4.0 * n * n * k, units)) // the reference only scales, or returns
        return false;
    bk_planes<C> At;
    cbk_pack(At, n, k, [&](long i, long l) { return conj(a[l + i * lda]); });
    return cbk_run(4.0 * n * n * k, units, [&](long u) {
        const long g0 = u * cb, g1 = std::min(g0 + cb, n);
        bk_planes<C> tv; // temp = zero, off-diagonal rows
        cbk_pack(tv, n, g1 - g0, [&](long, long) { return C(zero); });
        for (long l = 0; l < k; l++)
            for (long j = g0; j < g1; j++) {
                const long r0 = upper ? 0 : j + 1, cnt = upper ? j : n - j - 1;
                cbk_axpyp(cnt, tv.at(r0, j - g0), tv.ps, a[l + j * lda], At.at(r0, l), At.ps, false); // temp += conj(a(l,i))*a(l,j)
            }
        for (long j = g0; j < g1; j++) {
            const long r0 = upper ? 0 : j + 1, cnt = upper ? j : n - j - 1;
            for (long t = 0; t < cnt; t++) {
                C &cij = c[(r0 + t) + j * ldc];
                const C tt = tv.get(r0 + t, j - g0);
                if (beta == zero)
                    cij = alpha * tt;
                else
                    cij = alpha * tt + beta * cij;
            }
            R rtemp = zero;
            for (long l = 0; l < k; l++)
                rtemp += (conj(a[l + j * lda]) * a[l + j * lda]).real();
            C &d = c[j + j * ldc];
            if (beta == zero)
                d = alpha * rtemp;
            else
                d = alpha * rtemp + beta * d.real();
        }
    });
}

// ---- Ctrmm / Ctrsm, side = L, transa = T or C: blocks of columns of B, the
// rows processed one after the other as in the reference, each over the
// columns of the block at once (Bt: the block transposed) ----
template <class C> void ctrxm_lt_unit(bool solve, bool upper, bool noconj, bool nounit, long m, long n, const C &alpha, const C *a, long lda, C *b, long ldb, long u)
{
    const long j0 = u * cbk_rows, w = std::min(cbk_rows, n - j0);
    std::vector<C> Bt((size_t)(w * m)), tv((size_t)w);
    for (long i = 0; i < m; i++)
        for (long t = 0; t < w; t++)
            Bt[(size_t)(t + i * w)] = b[i + (j0 + t) * ldb];
    for (long ii = 0; ii < m; ii++) {
        const long i = (solve == upper) ? ii : m - 1 - ii;
        C *row = Bt.data() + i * w;
        for (long t = 0; t < w; t++)
            tv[(size_t)t] = row[t];
        const long k0 = upper ? 0 : i + 1, k1 = upper ? i : m;
        if (solve) {
            // temp = alpha*b(i,j); temp = temp - a(k,i)*b(k,j); temp = temp/a(i,i)
            cbk_scale(w, tv.data(), 1, alpha);
            for (long k = k0; k < k1; k++)
                cbk_axpy(w, tv.data(), 1, noconj ? a[k + i * lda] : conj(a[k + i * lda]), Bt.data() + k * w, 1, true);
            if (nounit)
                for (long t = 0; t < w; t++)
                    tv[(size_t)t] = noconj ? tv[(size_t)t] / a[i + i * lda] : tv[(size_t)t] / conj(a[i + i * lda]);
        } else {
            // temp = b(i,j); temp = temp*a(i,i); temp += a(k,i)*b(k,j); b(i,j) = alpha*temp
            if (nounit)
                cbk_scale(w, tv.data(), 1, noconj ? a[i + i * lda] : conj(a[i + i * lda]));
            for (long k = k0; k < k1; k++)
                cbk_axpy(w, tv.data(), 1, noconj ? a[k + i * lda] : conj(a[k + i * lda]), Bt.data() + k * w, 1, false);
            cbk_scale(w, tv.data(), 1, alpha);
        }
        for (long t = 0; t < w; t++)
            row[t] = tv[(size_t)t];
    }
    for (long i = 0; i < m; i++)
        for (long t = 0; t < w; t++)
            b[i + (j0 + t) * ldb] = Bt[(size_t)(t + i * w)];
}

// ---- Ctrmm / Ctrsm (B is b) ----
template <class C> bool ctrxm_simd(bool solve, bool left, bool upper, bool notr, bool noconj, bool nounit, long m, long n, const C &alpha, const C *a, long lda, C *b, long ldb)
{
    const C zero(0.0, 0.0), one(1.0, 0.0);
    if (alpha == zero) // the reference sets B to zero
        return false;
    if (left) {
        if (!notr)
            return cbk_run(4.0 * m * m * n, (n + cbk_rows - 1) / cbk_rows, [&](long u) { ctrxm_lt_unit(solve, upper, noconj, nounit, m, n, alpha, a, lda, b, ldb, u); });
        const long cb = cbk_cols(), units = (n + cb - 1) / cb;
        if (!cbk_go(4.0 * m * m * n, units))
            return false;
        bk_planes<C> Ap;
        cbk_pack(Ap, m, m, [&](long i, long l) { return a[i + l * lda]; });
        return cbk_run(4.0 * m * m * n, units, [&](long u) { // columns [g0, g1) of B, on planes
            const long g0 = u * cb, g1 = std::min(g0 + cb, n);
            bk_planes<C> Bp;
            Bp.resize(m, g1 - g0);
            for (long j = g0; j < g1; j++) {
                if (solve && alpha != one)
                    cbk_scale(m, b + j * ldb, 1, alpha);
                for (long i = 0; i < m; i++)
                    Bp.set(i, j - g0, b[i + j * ldb]);
            }
            for (long kk = 0; kk < m; kk++) {
                for (long j = g0; j < g1; j++) {
                    const long t = j - g0;
                    const long k = upper == solve ? m - 1 - kk : kk; // solve: upper from the bottom
                    C bk = Bp.get(k, t);
                    if (bk == zero)
                        continue;
                    if (solve) {
                        if (nounit) {
                            bk = bk / a[k + k * lda];
                            Bp.set(k, t, bk);
                        }
                        if (upper)
                            cbk_axpyp(k, Bp.at(0, t), Bp.ps, bk, Ap.at(0, k), Ap.ps, true);
                        else
                            cbk_axpyp(m - k - 1, Bp.at(k + 1, t), Bp.ps, bk, Ap.at(k + 1, k), Ap.ps, true);
                    } else {
                        C temp = alpha * bk;
                        if (upper) {
                            cbk_axpyp(k, Bp.at(0, t), Bp.ps, temp, Ap.at(0, k), Ap.ps, false);
                            if (nounit)
                                temp = temp * a[k + k * lda];
                            Bp.set(k, t, temp);
                        } else {
                            Bp.set(k, t, nounit ? temp * a[k + k * lda] : temp);
                            cbk_axpyp(m - k - 1, Bp.at(k + 1, t), Bp.ps, temp, Ap.at(k + 1, k), Ap.ps, false);
                        }
                    }
                }
            }
            for (long j = g0; j < g1; j++)
                for (long i = 0; i < m; i++)
                    b[i + j * ldb] = Bp.get(i, j - g0);
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
