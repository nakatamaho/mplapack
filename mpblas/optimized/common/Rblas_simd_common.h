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
 * dd/td/qd Level 1/2/3 BLAS on the CPU with AVX2 (Rblas_cpu_common.h):
 * instead of one lane per element, a unit of work is a column of the result
 * (Rsyrk/Rsyr2k with trans = N, Rsymm/Rtrmm/Rtrsm with side L, Rger, Rsyr,
 * Rsyr2) or a block of rows (Rsymm/Rtrmm/Rtrsm with side R, Rgemv with
 * trans = N, Level 1), and the loops over the rows of a column are
 * vectorized: x86-64-v3 (AVX2 + FMA) and baseline variants
 * (Rgemm_simd_common.h), with libQD3's branch-free arithmetic on the words
 * (Rgemm_qd_ops_common.h).  The loops are those of the reference routines;
 * every element goes through the operations of the reference, in the same
 * order, so the result is the reference result bit for bit.  The other
 * routines and cases keep the lanes of Rblas_lanes_common.h.
 *
 *   MPLAPACK_BLAS_SIMD=0   use the lanes instead
 */

#ifndef _MPLAPACK_RBLAS_SIMD_COMMON_H_
#define _MPLAPACK_RBLAS_SIMD_COMMON_H_

#include "Rgemm_qd_ops_common.h"
#include "Rgemm_simd_common.h"
#include "Rblas_calls_common.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace mplapack_blas {

// ---- vectorized loops on unit-stride columns of words ----
// The words of an element are scalar variables, so that the compiler
// vectorizes the loops over the elements.

#define MPLAPACK_BK_DD_DECL(v) double v##0, v##1
#define MPLAPACK_BK_DD_LOAD(v, p) (v##0 = (p)[0], v##1 = (p)[1])
#define MPLAPACK_BK_DD_STORE(p, v) ((p)[0] = v##0, (p)[1] = v##1)
#define MPLAPACK_BK_DD_NEG(v) (v##0 = -v##0, v##1 = -v##1)
#define MPLAPACK_BK_DD_ZERO(p) ((p)[0] = 0.0, (p)[1] = 0.0)
#define MPLAPACK_BK_DD_MUL(r, a, b) mplapack_gemm::dd_bf_mul_w(a##0, a##1, b##0, b##1, r##0, r##1)
#define MPLAPACK_BK_DD_ADD(r, a, b) mplapack_gemm::dd_bf_add_w(a##0, a##1, b##0, b##1, r##0, r##1)

#define MPLAPACK_BK_TD_DECL(v) double v##0, v##1, v##2
#define MPLAPACK_BK_TD_LOAD(v, p) (v##0 = (p)[0], v##1 = (p)[1], v##2 = (p)[2])
#define MPLAPACK_BK_TD_STORE(p, v) ((p)[0] = v##0, (p)[1] = v##1, (p)[2] = v##2)
#define MPLAPACK_BK_TD_NEG(v) (v##0 = -v##0, v##1 = -v##1, v##2 = -v##2)
#define MPLAPACK_BK_TD_ZERO(p) ((p)[0] = 0.0, (p)[1] = 0.0, (p)[2] = 0.0)
#define MPLAPACK_BK_TD_MUL(r, a, b) mplapack_gemm::td_bf_mul_w(a##0, a##1, a##2, b##0, b##1, b##2, r##0, r##1, r##2)
#define MPLAPACK_BK_TD_ADD(r, a, b) mplapack_gemm::td_bf_add_w(a##0, a##1, a##2, b##0, b##1, b##2, r##0, r##1, r##2)

#define MPLAPACK_BK_QD_DECL(v) double v##0, v##1, v##2, v##3
#define MPLAPACK_BK_QD_LOAD(v, p) (v##0 = (p)[0], v##1 = (p)[1], v##2 = (p)[2], v##3 = (p)[3])
#define MPLAPACK_BK_QD_STORE(p, v) ((p)[0] = v##0, (p)[1] = v##1, (p)[2] = v##2, (p)[3] = v##3)
#define MPLAPACK_BK_QD_NEG(v) (v##0 = -v##0, v##1 = -v##1, v##2 = -v##2, v##3 = -v##3)
#define MPLAPACK_BK_QD_ZERO(p) ((p)[0] = 0.0, (p)[1] = 0.0, (p)[2] = 0.0, (p)[3] = 0.0)
#define MPLAPACK_BK_QD_MUL(r, a, b) mplapack_gemm::qd_bf_mul_w(a##0, a##1, a##2, a##3, b##0, b##1, b##2, b##3, r##0, r##1, r##2, r##3)
#define MPLAPACK_BK_QD_ADD(r, a, b) mplapack_gemm::qd_bf_add_w(a##0, a##1, a##2, a##3, b##0, b##1, b##2, b##3, r##0, r##1, r##2, r##3)

// T: the word macros above (DD, TD, QD); W: words per element.
#define MPLAPACK_BK_KERNELS(T, W)                                                                                         \
    /* y = y + t*x, or y = y - t*x (sub) */                                                                                \
    MPLAPACK_GEMM_CLONES static void bk_axpy_##T(long n, double *__restrict y, const double *t, const double *__restrict x, int sub) \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(tt);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(tt, t);                                                                                     \
        if (sub) {                                                                                                         \
            MPLAPACK_GEMM_SIMD                                                                                             \
            for (long k = 0; k < n; k++) {                                                                                 \
                MPLAPACK_BK_##T##_DECL(xv);                                                                                \
                MPLAPACK_BK_##T##_DECL(yv);                                                                                \
                MPLAPACK_BK_##T##_DECL(p);                                                                                 \
                MPLAPACK_BK_##T##_LOAD(xv, x + W * k);                                                                     \
                MPLAPACK_BK_##T##_LOAD(yv, y + W * k);                                                                     \
                MPLAPACK_BK_##T##_MUL(p, tt, xv);                                                                          \
                MPLAPACK_BK_##T##_NEG(p);                                                                                  \
                MPLAPACK_BK_##T##_ADD(yv, yv, p);                                                                          \
                MPLAPACK_BK_##T##_STORE(y + W * k, yv);                                                                    \
            }                                                                                                              \
        } else {                                                                                                           \
            MPLAPACK_GEMM_SIMD                                                                                             \
            for (long k = 0; k < n; k++) {                                                                                 \
                MPLAPACK_BK_##T##_DECL(xv);                                                                                \
                MPLAPACK_BK_##T##_DECL(yv);                                                                                \
                MPLAPACK_BK_##T##_DECL(p);                                                                                 \
                MPLAPACK_BK_##T##_LOAD(xv, x + W * k);                                                                     \
                MPLAPACK_BK_##T##_LOAD(yv, y + W * k);                                                                     \
                MPLAPACK_BK_##T##_MUL(p, tt, xv);                                                                          \
                MPLAPACK_BK_##T##_ADD(yv, yv, p);                                                                          \
                MPLAPACK_BK_##T##_STORE(y + W * k, yv);                                                                    \
            }                                                                                                              \
        }                                                                                                                  \
    }                                                                                                                      \
    /* y = y + ((x1*t1) + (x2*t2)) */                                                                                      \
    MPLAPACK_GEMM_CLONES static void bk_axpy2_##T(long n, double *__restrict y, const double *__restrict x1, const double *t1, const double *__restrict x2, const double *t2) \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(s1);                                                                                        \
        MPLAPACK_BK_##T##_DECL(s2);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(s1, t1);                                                                                    \
        MPLAPACK_BK_##T##_LOAD(s2, t2);                                                                                    \
        MPLAPACK_GEMM_SIMD                                                                                                 \
        for (long k = 0; k < n; k++) {                                                                                     \
            MPLAPACK_BK_##T##_DECL(a);                                                                                     \
            MPLAPACK_BK_##T##_DECL(b);                                                                                     \
            MPLAPACK_BK_##T##_DECL(yv);                                                                                    \
            MPLAPACK_BK_##T##_DECL(p);                                                                                     \
            MPLAPACK_BK_##T##_DECL(q);                                                                                     \
            MPLAPACK_BK_##T##_LOAD(a, x1 + W * k);                                                                         \
            MPLAPACK_BK_##T##_LOAD(b, x2 + W * k);                                                                         \
            MPLAPACK_BK_##T##_LOAD(yv, y + W * k);                                                                         \
            MPLAPACK_BK_##T##_MUL(p, a, s1);                                                                               \
            MPLAPACK_BK_##T##_MUL(q, b, s2);                                                                               \
            MPLAPACK_BK_##T##_ADD(p, p, q);                                                                                \
            MPLAPACK_BK_##T##_ADD(yv, yv, p);                                                                              \
            MPLAPACK_BK_##T##_STORE(y + W * k, yv);                                                                        \
        }                                                                                                                  \
    }                                                                                                                      \
    /* y = t*y */                                                                                                          \
    MPLAPACK_GEMM_CLONES static void bk_scale_##T(long n, double *__restrict y, const double *t)                           \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(tt);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(tt, t);                                                                                     \
        MPLAPACK_GEMM_SIMD                                                                                                 \
        for (long k = 0; k < n; k++) {                                                                                     \
            MPLAPACK_BK_##T##_DECL(yv);                                                                                    \
            MPLAPACK_BK_##T##_LOAD(yv, y + W * k);                                                                         \
            MPLAPACK_BK_##T##_MUL(yv, tt, yv);                                                                             \
            MPLAPACK_BK_##T##_STORE(y + W * k, yv);                                                                        \
        }                                                                                                                  \
    }                                                                                                                      \
    /* y = t*x, or y = (be*y) + (t*x) (beta) */                                                                            \
    MPLAPACK_GEMM_CLONES static void bk_symmr_##T(long n, double *__restrict y, const double *be, const double *t, const double *__restrict x, int beta) \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(tt);                                                                                        \
        MPLAPACK_BK_##T##_DECL(bb);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(tt, t);                                                                                     \
        MPLAPACK_BK_##T##_LOAD(bb, be);                                                                                    \
        if (beta) {                                                                                                        \
            MPLAPACK_GEMM_SIMD                                                                                             \
            for (long k = 0; k < n; k++) {                                                                                 \
                MPLAPACK_BK_##T##_DECL(xv);                                                                                \
                MPLAPACK_BK_##T##_DECL(yv);                                                                                \
                MPLAPACK_BK_##T##_DECL(p);                                                                                 \
                MPLAPACK_BK_##T##_LOAD(xv, x + W * k);                                                                     \
                MPLAPACK_BK_##T##_LOAD(yv, y + W * k);                                                                     \
                MPLAPACK_BK_##T##_MUL(yv, bb, yv);                                                                         \
                MPLAPACK_BK_##T##_MUL(p, tt, xv);                                                                          \
                MPLAPACK_BK_##T##_ADD(yv, yv, p);                                                                          \
                MPLAPACK_BK_##T##_STORE(y + W * k, yv);                                                                    \
            }                                                                                                              \
        } else {                                                                                                           \
            MPLAPACK_GEMM_SIMD                                                                                             \
            for (long k = 0; k < n; k++) {                                                                                 \
                MPLAPACK_BK_##T##_DECL(xv);                                                                                \
                MPLAPACK_BK_##T##_DECL(p);                                                                                 \
                MPLAPACK_BK_##T##_LOAD(xv, x + W * k);                                                                     \
                MPLAPACK_BK_##T##_MUL(p, tt, xv);                                                                          \
                MPLAPACK_BK_##T##_STORE(y + W * k, p);                                                                     \
            }                                                                                                              \
        }                                                                                                                  \
    }                                                                                                                      \
    /* x' = (c*x) + (s*y), y' = (c*y) - (s*x) */                                                                           \
    MPLAPACK_GEMM_CLONES static void bk_rot_##T(long n, double *__restrict x, double *__restrict y, const double *c, const double *s) \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(cc);                                                                                        \
        MPLAPACK_BK_##T##_DECL(ss);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(cc, c);                                                                                     \
        MPLAPACK_BK_##T##_LOAD(ss, s);                                                                                     \
        MPLAPACK_GEMM_SIMD                                                                                                 \
        for (long k = 0; k < n; k++) {                                                                                     \
            MPLAPACK_BK_##T##_DECL(xv);                                                                                    \
            MPLAPACK_BK_##T##_DECL(yv);                                                                                    \
            MPLAPACK_BK_##T##_DECL(u);                                                                                     \
            MPLAPACK_BK_##T##_DECL(p);                                                                                     \
            MPLAPACK_BK_##T##_DECL(q);                                                                                     \
            MPLAPACK_BK_##T##_LOAD(xv, x + W * k);                                                                         \
            MPLAPACK_BK_##T##_LOAD(yv, y + W * k);                                                                         \
            MPLAPACK_BK_##T##_MUL(u, cc, xv);                                                                              \
            MPLAPACK_BK_##T##_MUL(p, ss, yv);                                                                              \
            MPLAPACK_BK_##T##_ADD(u, u, p);                                                                                \
            MPLAPACK_BK_##T##_MUL(p, cc, yv);                                                                              \
            MPLAPACK_BK_##T##_MUL(q, ss, xv);                                                                              \
            MPLAPACK_BK_##T##_NEG(q);                                                                                      \
            MPLAPACK_BK_##T##_ADD(p, p, q);                                                                                \
            MPLAPACK_BK_##T##_STORE(y + W * k, p);                                                                         \
            MPLAPACK_BK_##T##_STORE(x + W * k, u);                                                                         \
        }                                                                                                                  \
    }

#if defined(MPLAPACK_BUILD_WITH_DD)
MPLAPACK_BK_KERNELS(DD, 2)
#define MPLAPACK_BK_TYPE DD
#elif defined(MPLAPACK_BUILD_WITH_TD)
MPLAPACK_BK_KERNELS(TD, 3)
#define MPLAPACK_BK_TYPE TD
#elif defined(MPLAPACK_BUILD_WITH_QD)
MPLAPACK_BK_KERNELS(QD, 4)
#define MPLAPACK_BK_TYPE QD
#endif

#ifdef MPLAPACK_BK_TYPE
#define MPLAPACK_BK_CAT2(a, b) a##b
#define MPLAPACK_BK_CAT(a, b) MPLAPACK_BK_CAT2(a, b)
#define MPLAPACK_BK_FN(name) MPLAPACK_BK_CAT(name, MPLAPACK_BK_TYPE)

#define MPLAPACK_BLAS_SIMD_AVAILABLE 1

// ---- strided forms: the same operations with the operators of T ----

template <class T> inline const double *bk_w(const T *p) { return reinterpret_cast<const double *>(p); }
template <class T> inline double *bk_w(T *p) { return reinterpret_cast<double *>(p); }

// y(k) = y(k) + t*x(k) (sub: y(k) - t*x(k)), k < n, with element strides
template <class T> void bk_axpy(long n, T *y, long sy, const T &t, const T *x, long sx, bool sub)
{
    if (n <= 0)
        return;
    if (sy == 1 && sx == 1) {
        MPLAPACK_BK_FN(bk_axpy_)(n, bk_w(y), bk_w(&t), bk_w(x), sub ? 1 : 0);
        return;
    }
    for (long k = 0; k < n; k++) {
        const T p = t * x[k * sx];
        y[k * sy] = sub ? y[k * sy] - p : y[k * sy] + p;
    }
}

// y(k) = y(k) + ((x1(k)*t1) + (x2(k)*t2))
template <class T> void bk_axpy2(long n, T *y, long sy, const T *x1, long s1, const T &t1, const T *x2, long s2, const T &t2)
{
    if (n <= 0)
        return;
    if (sy == 1 && s1 == 1 && s2 == 1) {
        MPLAPACK_BK_FN(bk_axpy2_)(n, bk_w(y), bk_w(x1), bk_w(&t1), bk_w(x2), bk_w(&t2));
        return;
    }
    for (long k = 0; k < n; k++) {
        const T p = x1[k * s1] * t1;
        const T q = x2[k * s2] * t2;
        const T r = p + q;
        y[k * sy] = y[k * sy] + r;
    }
}

// y(k) = t*y(k)
template <class T> void bk_scale(long n, T *y, long sy, const T &t)
{
    if (n <= 0)
        return;
    if (sy == 1) {
        MPLAPACK_BK_FN(bk_scale_)(n, bk_w(y), bk_w(&t));
        return;
    }
    for (long k = 0; k < n; k++)
        y[k * sy] = t * y[k * sy];
}

template <class T> void bk_zero(long n, T *y, long sy)
{
    for (long k = 0; k < n; k++)
        y[k * sy] = 0.0;
}

// y = beta*y as the first step of the routines (unchanged when beta == 1)
template <class T> void bk_beta(long n, T *y, long sy, const blas_call<T> &c)
{
    if (c.s.beta_zero)
        bk_zero(n, y, sy);
    else if (!c.s.beta_one)
        bk_scale(n, y, sy, c.beta);
}

// y = t*x (beta == 0) or y = (be*y) + (t*x); unit strides
template <class T> void bk_symmr(long n, T *y, const T &be, const T &t, const T *x, bool beta)
{
    if (n > 0)
        MPLAPACK_BK_FN(bk_symmr_)(n, bk_w(y), bk_w(&be), bk_w(&t), bk_w(x), beta ? 1 : 0);
}

template <class T> void bk_rot(long n, T *x, long sx, T *y, long sy, const T &c, const T &s)
{
    if (n <= 0)
        return;
    if (sx == 1 && sy == 1) {
        MPLAPACK_BK_FN(bk_rot_)(n, bk_w(x), bk_w(y), bk_w(&c), bk_w(&s));
        return;
    }
    for (long k = 0; k < n; k++) {
        T u = c * x[k * sx];
        T p = s * y[k * sy];
        u = u + p;
        p = c * y[k * sy];
        const T q = s * x[k * sx];
        y[k * sy] = p - q;
        x[k * sx] = u;
    }
}

// ---- units of work ----

template <class T> struct bk_mat {
    T *p;
    long ld;
    T &operator()(long i, long j) const { return p[i + j * ld]; }
    T *col(long i, long j) const { return p + i + j * ld; }
};

template <class T> inline bk_mat<T> bk_of(const blas_operand<T> &o)
{
    bk_mat<T> m = {o.X, o.ld};
    return m;
}

// Number of units of a call in this mode, 0 if the call keeps the lanes.
template <class T> long bk_units(const blas_call<T> &c, long rb)
{
    const blas_shape &s = c.s;
    switch (s.op) {
    case L3_SYRK:
    case L3_SYR2K:
        return s.trans ? 0 : s.n;
    case L3_SYMM:
        return s.left ? s.n : (s.m + rb - 1) / rb;
    case L3_TRMM:
    case L3_TRSM:
        if (s.left)
            return s.trans ? 0 : s.n;
        return (s.m + rb - 1) / rb;
    case L2_GEMV:
        return s.trans ? 0 : (s.m + rb - 1) / rb;
    case L2_GER:
    case L2_SYR:
    case L2_SYR2:
        return s.n;
    case L1_SCAL:
    case L1_SWAP:
    case L1_ROT:
        return (s.n + rb - 1) / rb;
    default:
        return 0;
    }
}

template <class T> void bk_unit(const blas_call<T> &c, long u, long rb)
{
    const blas_shape &s = c.s;
    const bk_mat<T> A = bk_of(c.a), B = bk_of(c.b), C = bk_of(c.c);
    const T &al = c.alpha;
    const long m = s.m, n = s.n;
    switch (s.op) {
    case L3_SYRK:
    case L3_SYR2K: { // column j, trans = N
        const long j = u, r0 = s.upper ? 0 : j, r1 = s.upper ? j + 1 : n;
        T *cj = C.col(r0, j);
        bk_beta(r1 - r0, cj, 1, c);
        for (long l = 0; l < s.k; l++) {
            if (s.op == L3_SYRK) {
                if (A(j, l) == 0.0)
                    continue;
                const T t = al * A(j, l); // temp = alpha * a(j,l)
                bk_axpy(r1 - r0, cj, 1, t, A.col(r0, l), 1, false);
            } else {
                if (A(j, l) == 0.0 && B(j, l) == 0.0)
                    continue;
                const T t1 = al * B(j, l), t2 = al * A(j, l);
                bk_axpy2(r1 - r0, cj, 1, A.col(r0, l), 1, t1, B.col(r0, l), 1, t2);
            }
        }
        break;
    }
    case L3_SYMM:
        if (s.left) { // column j
            const long j = u;
            for (long ii = 0; ii < m; ii++) {
                const long i = s.upper ? ii : m - 1 - ii;
                const T t1 = al * B(i, j); // temp1 = alpha * b(i,j)
                T t2 = 0.0;
                const long k0 = s.upper ? 0 : i + 1, k1 = s.upper ? i : m;
                bk_axpy(k1 - k0, C.col(k0, j), 1, t1, A.col(k0, i), 1, false); // c(k,j) += temp1*a(k,i)
                for (long k = k0; k < k1; k++) {
                    const T p = B(k, j) * A(k, i); // temp2 += b(k,j)*a(k,i)
                    t2 = t2 + p;
                }
                T v;
                if (s.beta_zero) {
                    v = t1 * A(i, i);
                } else {
                    v = c.beta * C(i, j);
                    const T p = t1 * A(i, i);
                    v = v + p;
                }
                const T p = al * t2;
                C(i, j) = v + p;
            }
        } else { // rows [i0, i1)
            const long i0 = u * rb, cnt = std::min(rb, m - i0);
            for (long j = 0; j < n; j++) {
                const T t1 = al * A(j, j); // temp1 = alpha * a(j,j)
                bk_symmr(cnt, C.col(i0, j), c.beta, t1, B.col(i0, j), !s.beta_zero);
                for (long k = 0; k < j; k++) {
                    const T t = al * (s.upper ? A(k, j) : A(j, k));
                    bk_axpy(cnt, C.col(i0, j), 1, t, B.col(i0, k), 1, false);
                }
                for (long k = j + 1; k < n; k++) {
                    const T t = al * (s.upper ? A(j, k) : A(k, j));
                    bk_axpy(cnt, C.col(i0, j), 1, t, B.col(i0, k), 1, false);
                }
            }
        }
        break;
    case L3_TRMM: // B is C
        if (s.left) { // column j, transa = N
            const long j = u;
            for (long kk = 0; kk < m; kk++) {
                const long k = s.upper ? kk : m - 1 - kk;
                if (C(k, j) == 0.0)
                    continue;
                T t = al * C(k, j); // temp = alpha * b(k,j)
                if (s.upper) {
                    bk_axpy(k, C.col(0, j), 1, t, A.col(0, k), 1, false);
                    if (s.nounit)
                        t = t * A(k, k);
                    C(k, j) = t;
                } else {
                    if (s.nounit)
                        C(k, j) = t * A(k, k);
                    else
                        C(k, j) = t;
                    bk_axpy(m - k - 1, C.col(k + 1, j), 1, t, A.col(k + 1, k), 1, false);
                }
            }
        } else { // rows [i0, i1)
            const long i0 = u * rb, cnt = std::min(rb, m - i0);
            if (!s.trans) {
                for (long jj = 0; jj < n; jj++) {
                    const long j = s.upper ? n - 1 - jj : jj;
                    T t = al; // temp = alpha [* a(j,j)]
                    if (s.nounit)
                        t = t * A(j, j);
                    bk_scale(cnt, C.col(i0, j), 1, t);
                    const long k0 = s.upper ? 0 : j + 1, k1 = s.upper ? j : n;
                    for (long k = k0; k < k1; k++)
                        if (!(A(k, j) == 0.0)) {
                            const T tk = al * A(k, j);
                            bk_axpy(cnt, C.col(i0, j), 1, tk, C.col(i0, k), 1, false);
                        }
                }
            } else {
                for (long kk = 0; kk < n; kk++) {
                    const long k = s.upper ? kk : n - 1 - kk;
                    const long j0 = s.upper ? 0 : k + 1, j1 = s.upper ? k : n;
                    for (long j = j0; j < j1; j++)
                        if (!(A(j, k) == 0.0)) {
                            const T tj = al * A(j, k);
                            bk_axpy(cnt, C.col(i0, j), 1, tj, C.col(i0, k), 1, false);
                        }
                    T t = al;
                    if (s.nounit)
                        t = t * A(k, k);
                    if (!(t == 1.0))
                        bk_scale(cnt, C.col(i0, k), 1, t);
                }
            }
        }
        break;
    case L3_TRSM: // B is C
        if (s.left) { // column j, transa = N
            const long j = u;
            if (!s.alpha_one)
                bk_scale(m, C.col(0, j), 1, al);
            for (long kk = 0; kk < m; kk++) {
                const long k = s.upper ? m - 1 - kk : kk;
                if (C(k, j) == 0.0)
                    continue;
                if (s.nounit)
                    C(k, j) = C(k, j) / A(k, k);
                const T bk = C(k, j);
                if (s.upper)
                    bk_axpy(k, C.col(0, j), 1, bk, A.col(0, k), 1, true);
                else
                    bk_axpy(m - k - 1, C.col(k + 1, j), 1, bk, A.col(k + 1, k), 1, true);
            }
        } else { // rows [i0, i1)
            const long i0 = u * rb, cnt = std::min(rb, m - i0);
            const T one = 1.0;
            if (!s.trans) {
                for (long jj = 0; jj < n; jj++) {
                    const long j = s.upper ? jj : n - 1 - jj;
                    if (!s.alpha_one)
                        bk_scale(cnt, C.col(i0, j), 1, al);
                    const long k0 = s.upper ? 0 : j + 1, k1 = s.upper ? j : n;
                    for (long k = k0; k < k1; k++)
                        if (!(A(k, j) == 0.0))
                            bk_axpy(cnt, C.col(i0, j), 1, A(k, j), C.col(i0, k), 1, true);
                    if (s.nounit) {
                        const T t = one / A(j, j);
                        bk_scale(cnt, C.col(i0, j), 1, t);
                    }
                }
            } else {
                for (long kk = 0; kk < n; kk++) {
                    const long k = s.upper ? n - 1 - kk : kk;
                    if (s.nounit) {
                        const T t = one / A(k, k);
                        bk_scale(cnt, C.col(i0, k), 1, t);
                    }
                    const long j0 = s.upper ? 0 : k + 1, j1 = s.upper ? k : n;
                    for (long j = j0; j < j1; j++)
                        if (!(A(j, k) == 0.0))
                            bk_axpy(cnt, C.col(i0, j), 1, A(j, k), C.col(i0, k), 1, true);
                    if (!s.alpha_one)
                        bk_scale(cnt, C.col(i0, k), 1, al);
                }
            }
        }
        break;
    case L2_GEMV: { // rows [i0, i1) of y, trans = N; x is b, y is c
        const long i0 = u * rb, cnt = std::min(rb, m - i0);
        T *y = c.c.X + i0 * c.c.ld;
        bk_beta(cnt, y, c.c.ld, c);
        for (long j = 0; j < n; j++) {
            const T t = al * c.b.X[j * c.b.ld]; // temp = alpha * x(j)
            bk_axpy(cnt, y, c.c.ld, t, A.col(i0, j), 1, false);
        }
        break;
    }
    case L2_GER:
    case L2_SYR:
    case L2_SYR2: { // column j of A (c); x is a, y is b
        const long j = u;
        const long r0 = s.op == L2_GER ? 0 : (s.upper ? 0 : j);
        const long r1 = s.op == L2_GER ? m : (s.upper ? j + 1 : n);
        const T *x = c.a.X, *yv = c.b.X;
        const long ix = c.a.ld, iy = c.b.ld;
        if (s.op == L2_SYR2) {
            if (x[j * ix] == 0.0 && yv[j * iy] == 0.0)
                break;
            const T t1 = al * yv[j * iy], t2 = al * x[j * ix];
            bk_axpy2(r1 - r0, C.col(r0, j), 1, x + r0 * ix, ix, t1, yv + r0 * iy, iy, t2);
        } else {
            const T &w = s.op == L2_GER ? yv[j * iy] : x[j * ix];
            if (w == 0.0)
                break;
            const T t = al * w;
            // a(i,j) += x(i)*temp: the branch-free multiplication is commutative
            bk_axpy(r1 - r0, C.col(r0, j), 1, t, x + r0 * ix, ix, false);
        }
        break;
    }
    case L1_SCAL:
    case L1_SWAP:
    case L1_ROT: { // elements [i0, i1); x is a, y (Rscal: x) is c
        const long i0 = u * rb, cnt = std::min(rb, n - i0);
        T *y = c.c.X + i0 * c.c.ld;
        if (s.op == L1_SCAL) {
            bk_scale(cnt, y, c.c.ld, al);
        } else {
            T *x = c.a.X + i0 * c.a.ld;
            if (s.op == L1_ROT) {
                bk_rot(cnt, x, c.a.ld, y, c.c.ld, al, c.beta);
            } else {
                for (long k = 0; k < cnt; k++) {
                    const T t = x[k * c.a.ld];
                    x[k * c.a.ld] = y[k * c.c.ld];
                    y[k * c.c.ld] = t;
                }
            }
        }
        break;
    }
    default:
        break;
    }
}

// Runs the call in this mode; false when it keeps the lanes.
template <class T> bool blas_simd_cpu(blas_call<T> &c, bool parallel)
{
    static const long rb = 64;
    if (blas_env_long("MPLAPACK_BLAS_SIMD", 1) == 0)
        return false;
    const long units = bk_units(c, rb);
    if (units <= 0)
        return false;
#ifdef _OPENMP
    if (parallel) {
#pragma omp parallel for schedule(dynamic, 1)
        for (long u = 0; u < units; u++)
            bk_unit(c, u, rb);
        return true;
    }
#endif
    (void)parallel;
    for (long u = 0; u < units; u++)
        bk_unit(c, u, rb);
    return true;
}

#endif // MPLAPACK_BK_TYPE

} // namespace mplapack_blas

#endif
