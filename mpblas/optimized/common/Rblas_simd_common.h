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
 * instead of one lane per element, a unit of work is a group of columns of
 * the result (Rsyrk/Rsyr2k, Rsymm, Rtrmm/Rtrsm with side L and transa = N;
 * bk_units), a column (Rger, Rsyr, Rsyr2), a block of rows (Rsymm/Rtrmm/Rtrsm with side R, Rgemv with
 * trans = N, Level 1) or a block of columns (Rtrmm/Rtrsm with side L and
 * transa = T, on a transposed copy), and the loops over the rows of a column
 * (the columns of a row) are vectorized; Rsyrk/Rsyr2k with trans = T
 * accumulate the dot products of a column over a transposed copy of A, B: x86-64-v3 (AVX2 + FMA) and baseline variants
 * (Rgemm_simd_common.h), with libQD3's branch-free arithmetic on the words
 * (Rgemm_qd_ops_common.h).  The loops are those of the reference routines;
 * every element goes through the operations of the reference, in the same
 * order, so the result is the reference result bit for bit.  The other
 * routines and cases keep the lanes of Rblas_lanes_common.h.
 *
 *   MPLAPACK_BLAS_SIMD=0        use the lanes instead
 *   MPLAPACK_BLAS_SIMD_COLS=N   columns per group (8)
 *   MPLAPACK_BLAS_SIMD_PLANES=0 column groups on the interleaved words
 */

#ifndef _MPLAPACK_RBLAS_SIMD_COMMON_H_
#define _MPLAPACK_RBLAS_SIMD_COMMON_H_

#include <algorithm>
#include <vector>
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

// The same on planes (structure of arrays): word q of element k at P##q[k].
#define MPLAPACK_BK_DD_PTRS(P, p, s) double *__restrict P##0 = (p), *__restrict P##1 = (p) + (s)
#define MPLAPACK_BK_DD_CPTRS(P, p, s) const double *__restrict P##0 = (p), *__restrict P##1 = (p) + (s)
#define MPLAPACK_BK_DD_LOADI(v, P, k) (v##0 = P##0[k], v##1 = P##1[k])
#define MPLAPACK_BK_DD_STOREI(P, k, v) (P##0[k] = v##0, P##1[k] = v##1)
#define MPLAPACK_BK_TD_PTRS(P, p, s) double *__restrict P##0 = (p), *__restrict P##1 = (p) + (s), *__restrict P##2 = (p) + 2 * (s)
#define MPLAPACK_BK_TD_CPTRS(P, p, s) const double *__restrict P##0 = (p), *__restrict P##1 = (p) + (s), *__restrict P##2 = (p) + 2 * (s)
#define MPLAPACK_BK_TD_LOADI(v, P, k) (v##0 = P##0[k], v##1 = P##1[k], v##2 = P##2[k])
#define MPLAPACK_BK_TD_STOREI(P, k, v) (P##0[k] = v##0, P##1[k] = v##1, P##2[k] = v##2)
#define MPLAPACK_BK_QD_PTRS(P, p, s) double *__restrict P##0 = (p), *__restrict P##1 = (p) + (s), *__restrict P##2 = (p) + 2 * (s), *__restrict P##3 = (p) + 3 * (s)
#define MPLAPACK_BK_QD_CPTRS(P, p, s) const double *__restrict P##0 = (p), *__restrict P##1 = (p) + (s), *__restrict P##2 = (p) + 2 * (s), *__restrict P##3 = (p) + 3 * (s)
#define MPLAPACK_BK_QD_LOADI(v, P, k) (v##0 = P##0[k], v##1 = P##1[k], v##2 = P##2[k], v##3 = P##3[k])
#define MPLAPACK_BK_QD_STOREI(P, k, v) (P##0[k] = v##0, P##1[k] = v##1, P##2[k] = v##2, P##3[k] = v##3)

// Kernels on planes: y with plane stride ys, x with plane stride xs.
#define MPLAPACK_BK_PLANE_KERNELS(T)                                                                                       \
    /* y = y + t*x, or y = y - t*x (sub) */                                                                                \
    MPLAPACK_GEMM_CLONES static void bk_axpyp_##T(long n, double *y, long ys, const double *t, const double *x, long xs, int sub) \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(tt);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(tt, t);                                                                                     \
        MPLAPACK_BK_##T##_PTRS(Y, y, ys);                                                                                  \
        MPLAPACK_BK_##T##_CPTRS(X, x, xs);                                                                                 \
        if (sub) {                                                                                                         \
            MPLAPACK_GEMM_SIMD                                                                                             \
            for (long k = 0; k < n; k++) {                                                                                 \
                MPLAPACK_BK_##T##_DECL(xv);                                                                                \
                MPLAPACK_BK_##T##_DECL(yv);                                                                                \
                MPLAPACK_BK_##T##_DECL(p);                                                                                 \
                MPLAPACK_BK_##T##_LOADI(xv, X, k);                                                                         \
                MPLAPACK_BK_##T##_LOADI(yv, Y, k);                                                                         \
                MPLAPACK_BK_##T##_MUL(p, tt, xv);                                                                          \
                MPLAPACK_BK_##T##_NEG(p);                                                                                  \
                MPLAPACK_BK_##T##_ADD(yv, yv, p);                                                                          \
                MPLAPACK_BK_##T##_STOREI(Y, k, yv);                                                                        \
            }                                                                                                              \
        } else {                                                                                                           \
            MPLAPACK_GEMM_SIMD                                                                                             \
            for (long k = 0; k < n; k++) {                                                                                 \
                MPLAPACK_BK_##T##_DECL(xv);                                                                                \
                MPLAPACK_BK_##T##_DECL(yv);                                                                                \
                MPLAPACK_BK_##T##_DECL(p);                                                                                 \
                MPLAPACK_BK_##T##_LOADI(xv, X, k);                                                                         \
                MPLAPACK_BK_##T##_LOADI(yv, Y, k);                                                                         \
                MPLAPACK_BK_##T##_MUL(p, tt, xv);                                                                          \
                MPLAPACK_BK_##T##_ADD(yv, yv, p);                                                                          \
                MPLAPACK_BK_##T##_STOREI(Y, k, yv);                                                                        \
            }                                                                                                              \
        }                                                                                                                  \
    }                                                                                                                      \
    /* y = y + ((x1*t1) + (x2*t2)) */                                                                                      \
    MPLAPACK_GEMM_CLONES static void bk_axpy2p_##T(long n, double *y, long ys, const double *x1, const double *t1, const double *x2, const double *t2, long xs) \
    {                                                                                                                      \
        MPLAPACK_BK_##T##_DECL(s1);                                                                                        \
        MPLAPACK_BK_##T##_DECL(s2);                                                                                        \
        MPLAPACK_BK_##T##_LOAD(s1, t1);                                                                                    \
        MPLAPACK_BK_##T##_LOAD(s2, t2);                                                                                    \
        MPLAPACK_BK_##T##_PTRS(Y, y, ys);                                                                                  \
        MPLAPACK_BK_##T##_CPTRS(X, x1, xs);                                                                                \
        MPLAPACK_BK_##T##_CPTRS(Z, x2, xs);                                                                                \
        MPLAPACK_GEMM_SIMD                                                                                                 \
        for (long k = 0; k < n; k++) {                                                                                     \
            MPLAPACK_BK_##T##_DECL(a);                                                                                     \
            MPLAPACK_BK_##T##_DECL(b);                                                                                     \
            MPLAPACK_BK_##T##_DECL(yv);                                                                                    \
            MPLAPACK_BK_##T##_DECL(p);                                                                                     \
            MPLAPACK_BK_##T##_DECL(q);                                                                                     \
            MPLAPACK_BK_##T##_LOADI(a, X, k);                                                                              \
            MPLAPACK_BK_##T##_LOADI(b, Z, k);                                                                              \
            MPLAPACK_BK_##T##_LOADI(yv, Y, k);                                                                             \
            MPLAPACK_BK_##T##_MUL(p, a, s1);                                                                               \
            MPLAPACK_BK_##T##_MUL(q, b, s2);                                                                               \
            MPLAPACK_BK_##T##_ADD(p, p, q);                                                                                \
            MPLAPACK_BK_##T##_ADD(yv, yv, p);                                                                              \
            MPLAPACK_BK_##T##_STOREI(Y, k, yv);                                                                            \
        }                                                                                                                  \
    }

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
    /* p(k) = a(k)*x(k) */                                                                                                             \
    MPLAPACK_GEMM_CLONES static void bk_vmul_##T(long n, double *__restrict p, const double *__restrict a, const double *__restrict x) \
    {                                                                                                                                  \
        MPLAPACK_GEMM_SIMD                                                                                                             \
        for (long k = 0; k < n; k++) {                                                                                                 \
            MPLAPACK_BK_##T##_DECL(av);                                                                                                \
            MPLAPACK_BK_##T##_DECL(xv);                                                                                                \
            MPLAPACK_BK_##T##_DECL(r);                                                                                                 \
            MPLAPACK_BK_##T##_LOAD(av, a + W * k);                                                                                     \
            MPLAPACK_BK_##T##_LOAD(xv, x + W * k);                                                                                     \
            MPLAPACK_BK_##T##_MUL(r, av, xv);                                                                                          \
            MPLAPACK_BK_##T##_STORE(p + W * k, r);                                                                                     \
        }                                                                                                                              \
    }                                                                                                                                  \
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
MPLAPACK_BK_PLANE_KERNELS(DD)
#define MPLAPACK_BK_TYPE DD
#elif defined(MPLAPACK_BUILD_WITH_TD)
MPLAPACK_BK_KERNELS(TD, 3)
MPLAPACK_BK_PLANE_KERNELS(TD)
#define MPLAPACK_BK_TYPE TD
#elif defined(MPLAPACK_BUILD_WITH_QD)
MPLAPACK_BK_KERNELS(QD, 4)
MPLAPACK_BK_PLANE_KERNELS(QD)
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

// p(k) = a(k)*x(k), unit strides
template <class T> void bk_vmul(long n, T *p, const T *a, const T *x)
{
    if (n > 0)
        MPLAPACK_BK_FN(bk_vmul_)(n, bk_w(p), bk_w(a), bk_w(x));
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

// ---- packed planes ----
// A matrix of T stored word by word: word q of element (i, j) is
// v[q * ps + i + j * ld].  The packed copies of A (and B) are shared by the
// units; a unit packs its group of columns of the result, so that the loops
// over the rows run on unit-stride planes.

template <class T> struct bk_planes {
    static const int W = (int)(sizeof(T) / sizeof(double));
    std::vector<double> v;
    long ld, ps;
    bk_planes() : ld(0), ps(0) {}
    double *at(long i, long j) { return v.data() + i + j * ld; }
    const double *at(long i, long j) const { return v.data() + i + j * ld; }
    T get(long i, long j) const
    {
        T x;
        double *w = bk_w(&x);
        for (int q = 0; q < W; q++)
            w[q] = v[(size_t)(q * ps + i + j * ld)];
        return x;
    }
    void set(long i, long j, const T &x)
    {
        const double *w = bk_w(&x);
        for (int q = 0; q < W; q++)
            v[(size_t)(q * ps + i + j * ld)] = w[q];
    }
    void resize(long rows, long cols)
    {
        ld = rows, ps = rows * cols;
        v.resize((size_t)(W * ps));
    }
};

// y(k) = y(k) + t*x(k) (sub: y(k) - t*x(k)), k < n, on planes
template <class T> void bk_axpyp(long n, double *y, long ys, const T &t, const double *x, long xs, bool sub)
{
    if (n > 0)
        MPLAPACK_BK_FN(bk_axpyp_)(n, y, ys, bk_w(&t), x, xs, sub ? 1 : 0);
}

// y(k) = y(k) + ((x1(k)*t1) + (x2(k)*t2)), on planes
template <class T> void bk_axpy2p(long n, double *y, long ys, const double *x1, const T &t1, const double *x2, const T &t2, long xs)
{
    if (n > 0)
        MPLAPACK_BK_FN(bk_axpy2p_)(n, y, ys, x1, bk_w(&t1), x2, bk_w(&t2), xs);
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
// Units over columns take groups of cb columns: the columns of a group are
// independent and each goes through its own operations in the reference
// order, the steps of the group interleaved so that a column of A is used
// by the cb columns while it is in the cache.
template <class T> long bk_units(const blas_call<T> &c, long rb, long cb)
{
    const blas_shape &s = c.s;
    switch (s.op) {
    case L3_SYRK:
    case L3_SYR2K:
        return (s.n + cb - 1) / cb; // groups of columns of C
    case L3_SYMM:
        return s.left ? (s.n + cb - 1) / cb : (s.m + rb - 1) / rb;
    case L3_TRMM:
    case L3_TRSM:
        if (s.left) // groups of columns of B; transa = T: blocks of rb columns
            return s.trans ? (s.n + rb - 1) / rb : (s.n + cb - 1) / cb;
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

// syrk/syr2k with trans = T: c(i,j) from the sums temp(1,2) of column j
template <class T> void bk_syrk_t_finish(const blas_call<T> &c, long j, const std::vector<T> &t1, const std::vector<T> &t2)
{
    const blas_shape &s = c.s;
    const bk_mat<T> C = bk_of(c.c);
    const T &al = c.alpha;
    const long r0 = s.upper ? 0 : j, cnt = (long)t1.size();
    for (long t = 0; t < cnt; t++) {
        T &cij = C(r0 + t, j);
        if (s.op == L3_SYRK) {
            if (s.beta_zero)
                cij = al * t1[(size_t)t];
            else
                cij = al * t1[(size_t)t] + c.beta * cij;
        } else {
            if (s.beta_zero)
                cij = al * t1[(size_t)t] + al * t2[(size_t)t];
            else
                cij = c.beta * cij + al * t1[(size_t)t] + al * t2[(size_t)t];
        }
    }
}

// syrk/syr2k with trans = T: c(i,j) is the sum over l of a(l,i)*a(l,j)
// (syr2k: a(l,i)*b(l,j) and b(l,i)*a(l,j)), accumulated for the rows of
// column j at once over the transposed copies At, Bt (n x k, leading
// dimension n) of A and B, for the columns [j0, j1); the final combination
// is the reference one.
template <class T> void bk_syrk_t(const blas_call<T> &c, long j0, long j1, const T *At, const T *Bt)
{
    const blas_shape &s = c.s;
    const bk_mat<T> A = bk_of(c.a), B = bk_of(c.b);
    const long n = s.n, g = j1 - j0;
    // temp(1,2) = zero for the rows [r0, r1) of every column of the group
    std::vector<std::vector<T> > t1v((size_t)g), t2v((size_t)g);
    for (long j = j0; j < j1; j++) {
        const long cnt = s.upper ? j + 1 : n - j;
        t1v[(size_t)(j - j0)].resize((size_t)cnt);
        if (s.op == L3_SYR2K)
            t2v[(size_t)(j - j0)].resize((size_t)cnt);
    }
    for (long l = 0; l < s.k; l++)
        for (long j = j0; j < j1; j++) {
            const long r0 = s.upper ? 0 : j, cnt = s.upper ? j + 1 : n - j;
            T *t1 = t1v[(size_t)(j - j0)].data();
            if (s.op == L3_SYRK) {
                bk_axpy(cnt, t1, 1, A(l, j), At + r0 + l * n, 1, false); // temp += a(l,i)*a(l,j)
            } else {
                bk_axpy(cnt, t1, 1, B(l, j), At + r0 + l * n, 1, false);                          // temp1 += a(l,i)*b(l,j)
                bk_axpy(cnt, t2v[(size_t)(j - j0)].data(), 1, A(l, j), Bt + r0 + l * n, 1, false); // temp2 += b(l,i)*a(l,j)
            }
        }
    for (long j = j0; j < j1; j++)
        bk_syrk_t_finish(c, j, t1v[(size_t)(j - j0)], t2v[(size_t)(j - j0)]);
}

// trmm/trsm with side = L and transa = T on the columns [j0, j0 + w) of B:
// the rows of the block are processed one after the other as in the
// reference, each over the w columns at once (Bt: the block transposed,
// w x m, so that a row is contiguous).
template <class T> void bk_trxm_lt(const blas_call<T> &c, long u, long rb)
{
    const blas_shape &s = c.s;
    const bk_mat<T> A = bk_of(c.a), C = bk_of(c.c);
    const T &al = c.alpha;
    const long m = s.m, j0 = u * rb, w = std::min(rb, s.n - j0);
    std::vector<T> Bt((size_t)(w * m)), tv((size_t)w);
    for (long i = 0; i < m; i++)
        for (long t = 0; t < w; t++)
            Bt[(size_t)(t + i * w)] = C(i, j0 + t);
    for (long ii = 0; ii < m; ii++) {
        T *row;
        if (s.op == L3_TRSM) {
            // temp = alpha*b(i,j); temp = temp - a(k,i)*b(k,j) over the solved
            // rows k; temp = temp/a(i,i)
            const long i = s.upper ? ii : m - 1 - ii;
            row = Bt.data() + i * w;
            for (long t = 0; t < w; t++)
                tv[(size_t)t] = row[t];
            bk_scale(w, tv.data(), 1, al);
            const long k0 = s.upper ? 0 : i + 1, k1 = s.upper ? i : m;
            for (long k = k0; k < k1; k++)
                bk_axpy(w, tv.data(), 1, A(k, i), Bt.data() + k * w, 1, true);
            if (s.nounit)
                for (long t = 0; t < w; t++)
                    tv[(size_t)t] = tv[(size_t)t] / A(i, i);
        } else {
            // temp = b(i,j); temp = temp*a(i,i); temp += a(k,i)*b(k,j) over the
            // rows k not yet processed; b(i,j) = alpha*temp
            const long i = s.upper ? m - 1 - ii : ii;
            row = Bt.data() + i * w;
            for (long t = 0; t < w; t++)
                tv[(size_t)t] = row[t];
            if (s.nounit)
                bk_scale(w, tv.data(), 1, A(i, i));
            const long k0 = s.upper ? 0 : i + 1, k1 = s.upper ? i : m;
            for (long k = k0; k < k1; k++)
                bk_axpy(w, tv.data(), 1, A(k, i), Bt.data() + k * w, 1, false);
            bk_scale(w, tv.data(), 1, al);
        }
        for (long t = 0; t < w; t++)
            row[t] = tv[(size_t)t];
    }
    for (long i = 0; i < m; i++)
        for (long t = 0; t < w; t++)
            C(i, j0 + t) = Bt[(size_t)(t + i * w)];
}

// The column groups of Rsyrk/Rsyr2k (trans = N) and Rtrmm/Rtrsm (side L,
// transa = N) on planes: Ap, Bp are A and B packed, the columns [g0, g1) of
// C (B for trmm/trsm) are packed here and unpacked at the end.  The
// operations are those of bk_unit.
template <class T> void bk_unit_planes(const blas_call<T> &c, long g0, long g1, const bk_planes<T> &Ap, const bk_planes<T> &Bp)
{
    const blas_shape &s = c.s;
    const bk_mat<T> A = bk_of(c.a), B = bk_of(c.b), C = bk_of(c.c);
    const T &al = c.alpha;
    const long g = g1 - g0, rows = s.op == L3_SYRK || s.op == L3_SYR2K ? s.n : s.m;
    bk_planes<T> Cp;
    Cp.resize(rows, g);
    for (long j = g0; j < g1; j++) {
        const long r0 = s.op == L3_TRMM || s.op == L3_TRSM || s.upper ? 0 : j;
        const long r1 = s.op == L3_TRMM || s.op == L3_TRSM || !s.upper ? rows : j + 1;
        if (s.op == L3_SYRK || s.op == L3_SYR2K)
            bk_beta(r1 - r0, C.col(r0, j), 1, c);
        else if (s.op == L3_TRSM && !s.alpha_one)
            bk_scale(rows, C.col(0, j), 1, al);
        for (long i = r0; i < r1; i++)
            Cp.set(i, j - g0, C(i, j));
    }
    switch (s.op) {
    case L3_SYRK:
    case L3_SYR2K:
        for (long l = 0; l < s.k; l++)
            for (long j = g0; j < g1; j++) {
                const long r0 = s.upper ? 0 : j, r1 = s.upper ? j + 1 : rows;
                double *cj = Cp.at(r0, j - g0);
                if (s.op == L3_SYRK) {
                    if (A(j, l) == 0.0)
                        continue;
                    const T t = al * A(j, l); // temp = alpha * a(j,l)
                    bk_axpyp(r1 - r0, cj, Cp.ps, t, Ap.at(r0, l), Ap.ps, false);
                } else {
                    if (A(j, l) == 0.0 && B(j, l) == 0.0)
                        continue;
                    const T t1 = al * B(j, l), t2 = al * A(j, l);
                    bk_axpy2p(r1 - r0, cj, Cp.ps, Ap.at(r0, l), t1, Bp.at(r0, l), t2, Ap.ps);
                }
            }
        break;
    case L3_TRMM: {
        const long m = rows;
        for (long kk = 0; kk < m; kk++) {
            for (long j = g0; j < g1; j++) {
                const long k = s.upper ? kk : m - 1 - kk, t = j - g0;
                const T bkj = Cp.get(k, t);
                if (bkj == 0.0)
                    continue;
                T tk = al * bkj; // temp = alpha * b(k,j)
                if (s.upper) {
                    bk_axpyp(k, Cp.at(0, t), Cp.ps, tk, Ap.at(0, k), Ap.ps, false);
                    if (s.nounit)
                        tk = tk * A(k, k);
                    Cp.set(k, t, tk);
                } else {
                    if (s.nounit)
                        Cp.set(k, t, tk * A(k, k));
                    else
                        Cp.set(k, t, tk);
                    bk_axpyp(m - k - 1, Cp.at(k + 1, t), Cp.ps, tk, Ap.at(k + 1, k), Ap.ps, false);
                }
            }
        }
        break;
    }
    case L3_TRSM: {
        const long m = rows;
        for (long kk = 0; kk < m; kk++) {
            for (long j = g0; j < g1; j++) {
                const long k = s.upper ? m - 1 - kk : kk, t = j - g0;
                T bk = Cp.get(k, t);
                if (bk == 0.0)
                    continue;
                if (s.nounit) {
                    bk = bk / A(k, k);
                    Cp.set(k, t, bk);
                }
                if (s.upper)
                    bk_axpyp(k, Cp.at(0, t), Cp.ps, bk, Ap.at(0, k), Ap.ps, true);
                else
                    bk_axpyp(m - k - 1, Cp.at(k + 1, t), Cp.ps, bk, Ap.at(k + 1, k), Ap.ps, true);
            }
        }
        break;
    }
    default:
        break;
    }
    for (long j = g0; j < g1; j++) {
        const long r0 = s.op == L3_TRMM || s.op == L3_TRSM || s.upper ? 0 : j;
        const long r1 = s.op == L3_TRMM || s.op == L3_TRSM || !s.upper ? rows : j + 1;
        for (long i = r0; i < r1; i++)
            C(i, j) = Cp.get(i, j - g0);
    }
}

// The calls that run on planes (bk_unit_planes).
inline bool bk_on_planes(const blas_shape &s)
{
    return ((s.op == L3_SYRK || s.op == L3_SYR2K) && !s.trans) || ((s.op == L3_TRMM || s.op == L3_TRSM) && s.left && !s.trans);
}

template <class T> void bk_unit(const blas_call<T> &c, long u, long rb, long cb, const T *At = 0, const T *Bt = 0)
{
    const blas_shape &s = c.s;
    const bk_mat<T> A = bk_of(c.a), B = bk_of(c.b), C = bk_of(c.c);
    const T &al = c.alpha;
    const long m = s.m, n = s.n;
    // the group of columns [g0, g1) of a unit over columns
    const long g0 = u * cb, g1 = std::min(g0 + cb, n);
    if ((s.op == L3_SYRK || s.op == L3_SYR2K) && s.trans) {
        bk_syrk_t(c, g0, g1, At, Bt);
        return;
    }
    if ((s.op == L3_TRMM || s.op == L3_TRSM) && s.left && s.trans) {
        bk_trxm_lt(c, u, rb);
        return;
    }
    switch (s.op) {
    case L3_SYRK:
    case L3_SYR2K: { // columns [g0, g1), trans = N
        for (long j = g0; j < g1; j++) {
            const long r0 = s.upper ? 0 : j, r1 = s.upper ? j + 1 : n;
            bk_beta(r1 - r0, C.col(r0, j), 1, c);
        }
        for (long l = 0; l < s.k; l++)
            for (long j = g0; j < g1; j++) {
                const long r0 = s.upper ? 0 : j, r1 = s.upper ? j + 1 : n;
                T *cj = C.col(r0, j);
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
        if (s.left) { // columns [g0, g1)
            for (long ii = 0; ii < m; ii++) {
                for (long j = g0; j < g1; j++) {
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
        if (s.left) { // columns [g0, g1), transa = N
            for (long kk = 0; kk < m; kk++) {
                for (long j = g0; j < g1; j++) {
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
        if (s.left) { // columns [g0, g1), transa = N
            if (!s.alpha_one)
                for (long j = g0; j < g1; j++)
                    bk_scale(m, C.col(0, j), 1, al);
            for (long kk = 0; kk < m; kk++) {
                for (long j = g0; j < g1; j++) {
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
    const long cb = std::max(1L, blas_env_long("MPLAPACK_BLAS_SIMD_COLS", 8));
    const long units = bk_units(c, rb, cb);
    if (units <= 0)
        return false;
    // syrk/syr2k with trans = T: A and B transposed (n x k), shared by the units
    std::vector<T> At, Bt;
    if ((c.s.op == L3_SYRK || c.s.op == L3_SYR2K) && c.s.trans) {
        const long n = c.s.n, k = c.s.k;
        const bk_mat<T> A = bk_of(c.a), B = bk_of(c.b);
        At.resize((size_t)(n * k));
        if (c.s.op == L3_SYR2K)
            Bt.resize((size_t)(n * k));
        for (long l = 0; l < k; l++)
            for (long i = 0; i < n; i++) {
                At[(size_t)(i + l * n)] = A(l, i);
                if (c.s.op == L3_SYR2K)
                    Bt[(size_t)(i + l * n)] = B(l, i);
            }
    }
    const T *at = At.empty() ? 0 : At.data(), *bt = Bt.empty() ? 0 : Bt.data();
    // the column forms on planes: A (and B) packed, shared by the units
    bk_planes<T> Ap, Bp;
    const bool planes = bk_on_planes(c.s) && blas_env_long("MPLAPACK_BLAS_SIMD_PLANES", 1) != 0;
    if (planes) {
        const bool tr = c.s.op == L3_TRMM || c.s.op == L3_TRSM;
        const long rows = tr ? c.s.m : c.s.n, cols = tr ? c.s.m : c.s.k;
        const bk_mat<T> A = bk_of(c.a), B = bk_of(c.b);
        Ap.resize(rows, cols);
        if (c.s.op == L3_SYR2K)
            Bp.resize(rows, cols);
        for (long l = 0; l < cols; l++)
            for (long i = 0; i < rows; i++) {
                Ap.set(i, l, A(i, l));
                if (c.s.op == L3_SYR2K)
                    Bp.set(i, l, B(i, l));
            }
    }
#ifdef _OPENMP
    if (parallel) {
#pragma omp parallel for schedule(dynamic, 1)
        for (long u = 0; u < units; u++) {
            if (planes)
                bk_unit_planes(c, u * cb, std::min(u * cb + cb, c.s.n), Ap, Bp);
            else
                bk_unit(c, u, rb, cb, at, bt);
        }
        return true;
    }
#endif
    (void)parallel;
    for (long u = 0; u < units; u++) {
        if (planes)
            bk_unit_planes(c, u * cb, std::min(u * cb + cb, c.s.n), Ap, Bp);
        else
            bk_unit(c, u, rb, cb, at, bt);
    }
    return true;
}

#endif // MPLAPACK_BK_TYPE

} // namespace mplapack_blas

#endif
