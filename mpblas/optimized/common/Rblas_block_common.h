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
 * Small-block engine for the Level 3 updates of every CPU type.
 *
 * An update is, for every element (i, j) of an m x n block of C and for
 * l = 0, 1, ..., k - 1 in this order,
 *
 *   BLK_ADD:  c(i,j) = c(i,j) + coef(l,j) * x(i,l)
 *   BLK_SUB:  c(i,j) = c(i,j) - coef(l,j) * x(i,l)
 *   BLK_DOT:  c(i,j) = c(i,j) + x(i,l) * coef(l,j)
 *   BLK_SUBDOT:  c(i,j) = c(i,j) - x(i,l) * coef(l,j)
 *   BLK_ADD2: c(i,j) = c(i,j) + (x(i,l) * coef(l,j) + x(i,l+1) * coef(l+1,j))
 *             for the pairs l = 0, 2, 4, ... (Rsyr2k; k even, a pair skipped
 *             as a whole by skip(l, j))
 *
 * with the steps (l, j) for which skip(l, j) holds left out, as the
 * reference routines leave out the zero coefficients.  Only the order in
 * which the elements are visited changes; the chain of operations on one
 * element is that of the reference, so the result is the same bit for bit.
 *
 * The blocking is that of Goto and van de Geijn (see also
 * github.com/nakatamaho/dgemm_tutorial): NC columns of C, KC steps of l and
 * MC rows; for every (NC, MC) block and every KC chunk, x is packed into
 * micro-panels of MR rows and coef into micro-panels of NR columns, in the
 * order in which the micro-kernel reads them, in page-aligned per-thread
 * buffers (the panels are contiguous and touch as few pages, hence TLB
 * entries, as possible); the micro-kernel keeps an MR x NR tile of C in
 * registers over the KC steps.  The (NC, MC) blocks are distributed over
 * the threads.
 *
 * Kernels (K):
 *   value kernel (T trivially copyable: the REAL of the double, binary80 and
 *     binary128 libraries, mplapack_binary80_t and mplapack_binary128_t):
 *     the panels hold the values, the tile is an array of T;
 *   pointer kernel (MPFR, GMP): the panels hold pointers to the values, the
 *     kernel works on C in place (a copy would change the precision);
 *   plane kernels (dd, td, qd): the panels hold the words plane by plane and
 *     the tile is held word by word, with libQD3's branch-free operations.
 *
 *   MPLAPACK_BLAS_BLOCK=0       do not use the engine
 *   MPLAPACK_BLAS_BLOCK=2       use it for every type (default: blk_profitable)
 *   MPLAPACK_BLAS_BLOCK_MIN=N   smallest m*n*k (32768)
 */

#ifndef _MPLAPACK_RBLAS_BLOCK_COMMON_H_
#define _MPLAPACK_RBLAS_BLOCK_COMMON_H_

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <type_traits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Rgemm_simd_common.h"
#include <mpblas.h> // REAL, COMPLEX of the library
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "Rgemm_qd_ops_common.h"
#endif

namespace mplapack_blas {

enum { BLK_ADD = 0, BLK_SUB = 1, BLK_DOT = 2, BLK_SUBDOT = 3, BLK_ADD2 = 4 };

// Page-aligned scratch memory of one thread.
struct blk_buffer {
    void *p;
    size_t n;
    blk_buffer() : p(0), n(0) {}
    ~blk_buffer() { std::free(p); }
    void *get(size_t bytes)
    {
        if (bytes > n) {
            std::free(p);
            n = (bytes + 4095) / 4096 * 4096;
            p = std::aligned_alloc(4096, n);
        }
        return p;
    }
};

// ---- value kernel ----

// The AVX2 variant must not contract the products and sums into FMAs (the
// reference loops are compiled for the baseline, without FMA).
// The multiversioned kernels are static: their clones and resolvers are
// global symbols otherwise, defined in every object that includes this file.
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#define MPLAPACK_BLK_NO_CONTRACT __attribute__((optimize("fp-contract=off")))
#else
#define MPLAPACK_BLK_NO_CONTRACT
#endif

template <class T, int MR, int NR, int MODE, bool HM> MPLAPACK_GEMM_CLONES MPLAPACK_BLK_NO_CONTRACT static void blk_value_micro(long kc, const T *A, const T *B, const unsigned char *mask, T *C, long ldc, int mr, int nr)
{
    T acc[NR][MR];
    for (int jj = 0; jj < NR; jj++)
        for (int ii = 0; ii < MR; ii++)
            acc[jj][ii] = (ii < mr && jj < nr) ? C[ii + jj * ldc] : T(0.0);
    if (MODE == BLK_ADD2) {
        for (long l = 0; l + 1 < kc; l += 2) {
            const T *a1 = A + l * MR, *a2 = a1 + MR, *b1 = B + l * NR, *b2 = b1 + NR;
            for (int jj = 0; jj < NR; jj++) {
                if (HM && mask[l * NR + jj])
                    continue;
                const T c1 = b1[jj], c2 = b2[jj];
                MPLAPACK_GEMM_SIMD
                for (int ii = 0; ii < MR; ii++)
                    acc[jj][ii] += a1[ii] * c1 + a2[ii] * c2;
            }
        }
    } else
    for (long l = 0; l < kc; l++) {
        const T *a = A + l * MR, *b = B + l * NR;
        for (int jj = 0; jj < NR; jj++) {
            if (HM && mask[l * NR + jj])
                continue;
            const T bv = b[jj];
            if (MODE == BLK_ADD) {
                MPLAPACK_GEMM_SIMD
                for (int ii = 0; ii < MR; ii++)
                    acc[jj][ii] += bv * a[ii];
            } else if (MODE == BLK_SUB) {
                MPLAPACK_GEMM_SIMD
                for (int ii = 0; ii < MR; ii++)
                    acc[jj][ii] = acc[jj][ii] - bv * a[ii];
            } else if (MODE == BLK_DOT) {
                MPLAPACK_GEMM_SIMD
                for (int ii = 0; ii < MR; ii++)
                    acc[jj][ii] += a[ii] * bv;
            } else {
                MPLAPACK_GEMM_SIMD
                for (int ii = 0; ii < MR; ii++)
                    acc[jj][ii] = acc[jj][ii] - a[ii] * bv;
            }
        }
    }
    for (int jj = 0; jj < nr; jj++)
        for (int ii = 0; ii < mr; ii++)
            C[ii + jj * ldc] = acc[jj][ii];
}


template <class T, int MR_, int NR_, long MC_, long KC_, long NC_> struct blk_value_kernel {
    typedef T slot;
    static const int W = 1, MR = MR_, NR = NR_;
    static const long MC = MC_, KC = KC_, NC = NC_;
    static void put(slot *p, long, const T &v) { *p = v; }
    static void pad(slot *p, long) { *p = T(0.0); }
    // the tile of C at (C, ldc), mr x nr valid; A: kc x MR, B: kc x NR
    template <int MODE> static void micro(long kc, const slot *A, const slot *B, const unsigned char *mask, T *C, long ldc, int mr, int nr)
    {
        if (mask)
            blk_value_micro<T, MR, NR, MODE, true>(kc, A, B, mask, C, ldc, mr, nr);
        else
            blk_value_micro<T, MR, NR, MODE, false>(kc, A, B, mask, C, ldc, mr, nr);
    }
};

// ---- pointer kernel ----

template <class T, int MR_, int NR_, long MC_, long KC_, long NC_> struct blk_pointer_kernel {
    typedef const T *slot;
    static const int W = 1, MR = MR_, NR = NR_;
    static const long MC = MC_, KC = KC_, NC = NC_;
    static void put(slot *p, long, const T &v) { *p = &v; }
    static void pad(slot *p, long) { *p = 0; }
    template <int MODE> static void micro(long kc, const slot *A, const slot *B, const unsigned char *mask, T *C, long ldc, int mr, int nr)
    {
        if (MODE == BLK_ADD2) {
            for (long l = 0; l + 1 < kc; l += 2) {
                const slot *a1 = A + l * MR, *a2 = a1 + MR, *b1 = B + l * NR, *b2 = b1 + NR;
                for (int jj = 0; jj < nr; jj++) {
                    if (mask && mask[l * NR + jj])
                        continue;
                    T *c = C + jj * ldc;
                    for (int ii = 0; ii < mr; ii++)
                        c[ii] += *a1[ii] * *b1[jj] + *a2[ii] * *b2[jj];
                }
            }
            return;
        }
        for (long l = 0; l < kc; l++) {
            const slot *a = A + l * MR, *b = B + l * NR;
            for (int jj = 0; jj < nr; jj++) {
                if (mask && mask[l * NR + jj])
                    continue;
                const T &bv = *b[jj];
                T *c = C + jj * ldc;
                for (int ii = 0; ii < mr; ii++) {
                    if (MODE == BLK_ADD)
                        c[ii] += bv * *a[ii];
                    else if (MODE == BLK_SUB)
                        c[ii] = c[ii] - bv * *a[ii];
                    else if (MODE == BLK_DOT)
                        c[ii] += *a[ii] * bv;
                    else
                        c[ii] = c[ii] - *a[ii] * bv;
                }
            }
        }
    }
};

// ---- COMPLEX of the REAL library on planes ----

#if defined(MPLAPACK_BUILD_WITH_DOUBLE)

// The COMPLEX multiplication (std::complex of REAL) as GCC computes it: (a,b)*(c,d) =
// (a*c - b*d, a*d + b*c), and, when both parts are NaN, the C99 Annex G
// recovery (__muldc3); the kernel computes the first form vectorized and
// recomputes the lanes whose parts are both NaN with the operator.  Planes:
// real part of row ii of step l at A[l * 2 * MR + ii], imaginary part at
// A[l * 2 * MR + MR + ii]; B likewise with NR.

template <int MR, int NR, int MODE, bool HM> MPLAPACK_GEMM_CLONES MPLAPACK_BLK_NO_CONTRACT static void blk_zplane_micro(long kc, const REAL *A, const REAL *B, const unsigned char *mask, REAL *C, long ldc, int mr, int nr)
{
    REAL cr[NR][MR], ci[NR][MR];
    for (int jj = 0; jj < NR; jj++)
        for (int ii = 0; ii < MR; ii++) {
            const bool in = ii < mr && jj < nr;
            cr[jj][ii] = in ? C[(ii + jj * ldc) * 2] : 0.0;
            ci[jj][ii] = in ? C[(ii + jj * ldc) * 2 + 1] : 0.0;
        }
    if (MODE == BLK_ADD2) {
        // acc += (x1*c1 + x2*c2), each product as std::complex computes it
        for (long l = 0; l + 1 < kc; l += 2) {
            const REAL *a1 = A + l * 2 * MR, *a2 = a1 + 2 * MR, *b1 = B + l * 2 * NR, *b2 = b1 + 2 * NR;
            for (int jj = 0; jj < NR; jj++) {
                if (HM && mask[l * NR + jj])
                    continue;
                const REAL c1r = b1[jj], c1i = b1[NR + jj], c2r = b2[jj], c2i = b2[NR + jj];
                REAL pr[MR], pi[MR], qr[MR], qi[MR];
                int nan = 0;
                MPLAPACK_GEMM_SIMD
                for (int ii = 0; ii < MR; ii++) {
                    pr[ii] = a1[ii] * c1r - a1[MR + ii] * c1i;
                    pi[ii] = a1[ii] * c1i + a1[MR + ii] * c1r;
                    qr[ii] = a2[ii] * c2r - a2[MR + ii] * c2i;
                    qi[ii] = a2[ii] * c2i + a2[MR + ii] * c2r;
                    nan |= ((pr[ii] != pr[ii]) & (pi[ii] != pi[ii])) | ((qr[ii] != qr[ii]) & (qi[ii] != qi[ii]));
                }
                if (nan)
                    for (int ii = 0; ii < MR; ii++) {
                        if (pr[ii] != pr[ii] && pi[ii] != pi[ii]) {
                            const COMPLEX p = COMPLEX(a1[ii], a1[MR + ii]) * COMPLEX(c1r, c1i);
                            pr[ii] = p.real(), pi[ii] = p.imag();
                        }
                        if (qr[ii] != qr[ii] && qi[ii] != qi[ii]) {
                            const COMPLEX q = COMPLEX(a2[ii], a2[MR + ii]) * COMPLEX(c2r, c2i);
                            qr[ii] = q.real(), qi[ii] = q.imag();
                        }
                    }
                MPLAPACK_GEMM_SIMD
                for (int ii = 0; ii < MR; ii++)
                    cr[jj][ii] += pr[ii] + qr[ii], ci[jj][ii] += pi[ii] + qi[ii];
            }
        }
    } else
    for (long l = 0; l < kc; l++) {
        const REAL *a = A + l * 2 * MR, *b = B + l * 2 * NR;
        for (int jj = 0; jj < NR; jj++) {
            if (HM && mask[l * NR + jj])
                continue;
            const REAL br = b[jj], bi = b[NR + jj];
            REAL pr[MR], pi[MR];
            int nan = 0;
            MPLAPACK_GEMM_SIMD
            for (int ii = 0; ii < MR; ii++) {
                const REAL ar = a[ii], ai = a[MR + ii];
                if (MODE == BLK_DOT || MODE == BLK_SUBDOT) { // x*coef
                    pr[ii] = ar * br - ai * bi;
                    pi[ii] = ar * bi + ai * br;
                } else { // coef*x
                    pr[ii] = br * ar - bi * ai;
                    pi[ii] = br * ai + bi * ar;
                }
                nan |= (pr[ii] != pr[ii]) & (pi[ii] != pi[ii]);
            }
            if (nan)
                for (int ii = 0; ii < MR; ii++)
                    if (pr[ii] != pr[ii] && pi[ii] != pi[ii]) {
                        const COMPLEX x(a[ii], a[MR + ii]), t(br, bi);
                        const COMPLEX p = (MODE == BLK_DOT || MODE == BLK_SUBDOT) ? x * t : t * x;
                        pr[ii] = p.real(), pi[ii] = p.imag();
                    }
            if (MODE == BLK_SUB || MODE == BLK_SUBDOT) {
                MPLAPACK_GEMM_SIMD
                for (int ii = 0; ii < MR; ii++)
                    cr[jj][ii] -= pr[ii], ci[jj][ii] -= pi[ii];
            } else {
                MPLAPACK_GEMM_SIMD
                for (int ii = 0; ii < MR; ii++)
                    cr[jj][ii] += pr[ii], ci[jj][ii] += pi[ii];
            }
        }
    }
    for (int jj = 0; jj < nr; jj++)
        for (int ii = 0; ii < mr; ii++) {
            C[(ii + jj * ldc) * 2] = cr[jj][ii];
            C[(ii + jj * ldc) * 2 + 1] = ci[jj][ii];
        }
}

template <int MR_, int NR_, long MC_, long KC_, long NC_> struct blk_zplane_kernel {
    typedef REAL slot;
    typedef COMPLEX T;
    static_assert(sizeof(COMPLEX) == 2 * sizeof(REAL), "COMPLEX is two REAL");
    static const int W = 2, MR = MR_, NR = NR_;
    static const long MC = MC_, KC = KC_, NC = NC_;
    static void put(slot *p, long stride, const T &v) { p[0] = v.real(), p[stride] = v.imag(); }
    static void pad(slot *p, long stride) { p[0] = 0.0, p[stride] = 0.0; }
    template <int MODE> static void micro(long kc, const slot *A, const slot *B, const unsigned char *mask, T *C, long ldc, int mr, int nr)
    {
        if (mask)
            blk_zplane_micro<MR, NR, MODE, true>(kc, A, B, mask, reinterpret_cast<REAL *>(C), ldc, mr, nr);
        else
            blk_zplane_micro<MR, NR, MODE, false>(kc, A, B, mask, reinterpret_cast<REAL *>(C), ldc, mr, nr);
    }
};

#endif

// ---- plane kernels (dd, td, qd) ----

#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)

// Micro-kernel on planes: word q of row ii of step l of A at A[(l * W + q) * MR
// + ii], of column jj of B at B[(l * W + q) * NR + jj]; C holds elements of W
// words.  The tile is held word by word (acc[q][jj][ii]); the loop over the
// rows of the tile is vectorized.
#define MPLAPACK_BLK_PLANE_STEP(W, DECL, LOAD, STORE, MUL, ADD, NEG)                                                     \
    for (long l = 0; l < kc; l++) {                                                                                      \
        const double *a = A + l * W * MR, *b = B + l * W * NR;                                                           \
        for (int jj = 0; jj < NR; jj++) {                                                                                \
            if (HM && mask[l * NR + jj])                                                                                 \
                continue;                                                                                                \
            DECL(bv);                                                                                                    \
            LOAD(bv, b + jj, NR);                                                                                        \
            MPLAPACK_GEMM_SIMD                                                                                           \
            for (int ii = 0; ii < MR; ii++) {                                                                            \
                DECL(av);                                                                                                \
                DECL(cv);                                                                                                \
                DECL(p);                                                                                                 \
                LOAD(av, a + ii, MR);                                                                                    \
                LOAD(cv, &acc[0][jj][ii], NR * MR);                                                                      \
                if (MODE == BLK_DOT || MODE == BLK_SUBDOT)                                                              \
                    MUL(p, av, bv);                                                                                      \
                else                                                                                                     \
                    MUL(p, bv, av);                                                                                      \
                if (MODE == BLK_SUB || MODE == BLK_SUBDOT)                                                              \
                    NEG(p);                                                                                              \
                ADD(cv, cv, p);                                                                                          \
                STORE(&acc[0][jj][ii], NR * MR, cv);                                                                     \
            }                                                                                                            \
        }                                                                                                                \
    }

// The pairs of BLK_ADD2: acc += (x1*c1 + x2*c2).
#define MPLAPACK_BLK_PLANE_PAIR(W, DECL, LOAD, STORE, MUL, ADD)                                                          \
    for (long l = 0; l + 1 < kc; l += 2) {                                                                               \
        const double *a1 = A + l * W * MR, *a2 = a1 + W * MR, *b1 = B + l * W * NR, *b2 = b1 + W * NR;                   \
        for (int jj = 0; jj < NR; jj++) {                                                                                \
            if (HM && mask[l * NR + jj])                                                                                 \
                continue;                                                                                                \
            DECL(c1);                                                                                                    \
            DECL(c2);                                                                                                    \
            LOAD(c1, b1 + jj, NR);                                                                                       \
            LOAD(c2, b2 + jj, NR);                                                                                       \
            MPLAPACK_GEMM_SIMD                                                                                           \
            for (int ii = 0; ii < MR; ii++) {                                                                            \
                DECL(x1);                                                                                                \
                DECL(x2);                                                                                                \
                DECL(cv);                                                                                                \
                DECL(p);                                                                                                 \
                DECL(q);                                                                                                 \
                LOAD(x1, a1 + ii, MR);                                                                                   \
                LOAD(x2, a2 + ii, MR);                                                                                   \
                LOAD(cv, &acc[0][jj][ii], NR * MR);                                                                      \
                MUL(p, x1, c1);                                                                                          \
                MUL(q, x2, c2);                                                                                          \
                ADD(p, p, q);                                                                                            \
                ADD(cv, cv, p);                                                                                          \
                STORE(&acc[0][jj][ii], NR * MR, cv);                                                                     \
            }                                                                                                            \
        }                                                                                                                \
    }

#define MPLAPACK_BLK_DD_DECL(v) double v##0, v##1
#define MPLAPACK_BLK_DD_LOAD(v, p, s) (v##0 = (p)[0], v##1 = (p)[s])
#define MPLAPACK_BLK_DD_STORE(p, s, v) ((p)[0] = v##0, (p)[s] = v##1)
#define MPLAPACK_BLK_DD_MUL(r, a, b) mplapack_gemm::dd_bf_mul_w(a##0, a##1, b##0, b##1, r##0, r##1)
#define MPLAPACK_BLK_DD_ADD(r, a, b) mplapack_gemm::dd_bf_add_w(a##0, a##1, b##0, b##1, r##0, r##1)
#define MPLAPACK_BLK_DD_NEG(v) (v##0 = -v##0, v##1 = -v##1)
#define MPLAPACK_BLK_TD_DECL(v) double v##0, v##1, v##2
#define MPLAPACK_BLK_TD_LOAD(v, p, s) (v##0 = (p)[0], v##1 = (p)[s], v##2 = (p)[2 * (s)])
#define MPLAPACK_BLK_TD_STORE(p, s, v) ((p)[0] = v##0, (p)[s] = v##1, (p)[2 * (s)] = v##2)
#define MPLAPACK_BLK_TD_MUL(r, a, b) mplapack_gemm::td_bf_mul_w(a##0, a##1, a##2, b##0, b##1, b##2, r##0, r##1, r##2)
#define MPLAPACK_BLK_TD_ADD(r, a, b) mplapack_gemm::td_bf_add_w(a##0, a##1, a##2, b##0, b##1, b##2, r##0, r##1, r##2)
#define MPLAPACK_BLK_TD_NEG(v) (v##0 = -v##0, v##1 = -v##1, v##2 = -v##2)
#define MPLAPACK_BLK_QD_DECL(v) double v##0, v##1, v##2, v##3
#define MPLAPACK_BLK_QD_LOAD(v, p, s) (v##0 = (p)[0], v##1 = (p)[s], v##2 = (p)[2 * (s)], v##3 = (p)[3 * (s)])
#define MPLAPACK_BLK_QD_STORE(p, s, v) ((p)[0] = v##0, (p)[s] = v##1, (p)[2 * (s)] = v##2, (p)[3 * (s)] = v##3)
#define MPLAPACK_BLK_QD_MUL(r, a, b) mplapack_gemm::qd_bf_mul_w(a##0, a##1, a##2, a##3, b##0, b##1, b##2, b##3, r##0, r##1, r##2, r##3)
#define MPLAPACK_BLK_QD_ADD(r, a, b) mplapack_gemm::qd_bf_add_w(a##0, a##1, a##2, a##3, b##0, b##1, b##2, b##3, r##0, r##1, r##2, r##3)
#define MPLAPACK_BLK_QD_NEG(v) (v##0 = -v##0, v##1 = -v##1, v##2 = -v##2, v##3 = -v##3)

template <int W, int MR, int NR, int MODE, bool HM> struct blk_plane_steps;
#define MPLAPACK_BLK_PLANE_STEPS(W, T)                                                                                   \
    template <int MR, int NR, int MODE, bool HM> struct blk_plane_steps<W, MR, NR, MODE, HM> {                            \
        static MPLAPACK_GEMM_HD void run(long kc, const double *A, const double *B, const unsigned char *mask, double (*acc)[NR][MR]) \
        {                                                                                                                \
            if (MODE == BLK_ADD2) {                                                                                      \
                MPLAPACK_BLK_PLANE_PAIR(W, MPLAPACK_BLK_##T##_DECL, MPLAPACK_BLK_##T##_LOAD, MPLAPACK_BLK_##T##_STORE,   \
                                        MPLAPACK_BLK_##T##_MUL, MPLAPACK_BLK_##T##_ADD)                                  \
                return;                                                                                                  \
            }                                                                                                            \
            MPLAPACK_BLK_PLANE_STEP(W, MPLAPACK_BLK_##T##_DECL, MPLAPACK_BLK_##T##_LOAD, MPLAPACK_BLK_##T##_STORE,       \
                                    MPLAPACK_BLK_##T##_MUL, MPLAPACK_BLK_##T##_ADD, MPLAPACK_BLK_##T##_NEG)              \
        }                                                                                                                \
    };
MPLAPACK_BLK_PLANE_STEPS(2, DD)
MPLAPACK_BLK_PLANE_STEPS(3, TD)
MPLAPACK_BLK_PLANE_STEPS(4, QD)

template <int W, int MR, int NR, int MODE, bool HM> MPLAPACK_GEMM_CLONES static void blk_plane_micro(long kc, const double *A, const double *B, const unsigned char *mask, double *C, long ldc, int mr, int nr)
{
    double acc[W][NR][MR];
    for (int jj = 0; jj < NR; jj++)
        for (int ii = 0; ii < MR; ii++)
            for (int q = 0; q < W; q++)
                acc[q][jj][ii] = (ii < mr && jj < nr) ? C[(ii + jj * ldc) * W + q] : 0.0;
    blk_plane_steps<W, MR, NR, MODE, HM>::run(kc, A, B, mask, acc);
    for (int jj = 0; jj < nr; jj++)
        for (int ii = 0; ii < mr; ii++)
            for (int q = 0; q < W; q++)
                C[(ii + jj * ldc) * W + q] = acc[q][jj][ii];
}

template <class T, int W_, int MR_, int NR_, long MC_, long KC_, long NC_> struct blk_plane_kernel {
    typedef double slot;
    static const int W = W_, MR = MR_, NR = NR_;
    static const long MC = MC_, KC = KC_, NC = NC_;
    static void put(slot *p, long stride, const T &v)
    {
        const double *w = reinterpret_cast<const double *>(&v);
        for (int q = 0; q < W; q++)
            p[q * stride] = w[q];
    }
    static void pad(slot *p, long stride)
    {
        for (int q = 0; q < W; q++)
            p[q * stride] = 0.0;
    }
    template <int MODE> static void micro(long kc, const slot *A, const slot *B, const unsigned char *mask, T *C, long ldc, int mr, int nr)
    {
        if (mask)
            blk_plane_micro<W, MR, NR, MODE, true>(kc, A, B, mask, reinterpret_cast<double *>(C), ldc, mr, nr);
        else
            blk_plane_micro<W, MR, NR, MODE, false>(kc, A, B, mask, reinterpret_cast<double *>(C), ldc, mr, nr);
    }
};

#endif

// ---- the engine ----

struct blk_no_skip {
    bool operator()(long, long) const { return false; }
};

// All the (MC, NC) blocks of the update.
struct blk_all_tiles {
    bool operator()(long, long) const { return true; }
};

// Context: the arithmetic context of the caller, installed in every worker
// (blas_thread_context of Rblas_cpu_common.h).  x(i, l) and coef(l, j) return
// references that stay valid during the call.  With tile > 0, the blocks
// are tile x tile and only those (bi, bj) for which tiles(bi, bj) holds are
// updated (Rsyrk: the blocks off the diagonal of a triangle).
template <class K, class Context, class T, class FX, class FB, class FS, class FT = blk_all_tiles> void blk_update(int mode, long m, long n, long k, T *C, long ldc, FX x, FB coef, FS skip, bool use_skip, bool parallel, long tile = 0, FT tiles = FT())
{
    typedef typename K::slot S;
    const int MR = K::MR, NR = K::NR, W = K::W;
    if (m <= 0 || n <= 0 || k <= 0)
        return;
    int nthreads = 1;
#ifdef _OPENMP
    if (parallel)
        nthreads = omp_get_max_threads();
#endif
    long MC = tile > 0 ? tile : K::MC, NC = tile > 0 ? tile : K::NC;
    const long KC = K::KC;
    const long nbi = (m + MC - 1) / MC;
    // narrower column blocks when there are too few blocks for the threads
    while (tile <= 0 && NC > NR && nbi * ((n + NC - 1) / NC) < 2L * nthreads)
        NC = std::max((long)NR, (NC / 2 + NR - 1) / NR * NR);
    const long nbj = (n + NC - 1) / NC;
    const Context ctx;
#ifdef _OPENMP
#pragma omp parallel if (parallel && nthreads > 1)
#endif
    {
        const typename Context::scope sc(ctx);
        blk_buffer abuf, bbuf, mbuf;
        S *Ap = static_cast<S *>(abuf.get(sizeof(S) * (size_t)(W * ((MC + MR - 1) / MR * MR) * KC)));
        S *Bp = static_cast<S *>(bbuf.get(sizeof(S) * (size_t)(W * ((NC + NR - 1) / NR * NR) * KC)));
        unsigned char *Mp = use_skip ? static_cast<unsigned char *>(mbuf.get((size_t)(((NC + NR - 1) / NR * NR) * KC))) : 0;
        std::vector<unsigned char> anyskip((size_t)((NC + NR - 1) / NR));
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (long blk = 0; blk < nbi * nbj; blk++) {
            if (!tiles(blk % nbi, blk / nbi))
                continue;
            const long i0 = (blk % nbi) * MC, j0 = (blk / nbi) * NC;
            const long mc = std::min(MC, m - i0), nc = std::min(NC, n - j0);
            for (long p0 = 0; p0 < k; p0 += KC) {
                const long kc = std::min(KC, k - p0);
                // B: micro-panels of NR columns, kc x NR each; anyskip: whether
                // a micro-panel has steps to leave out (otherwise its kernel
                // runs without the mask)
                for (long jr = 0; jr < nc; jr += NR) {
                    S *bp = Bp + jr * W * kc;
                    unsigned char *mp = Mp ? Mp + jr * kc : 0;
                    unsigned char any = 0;
                    for (long l = 0; l < kc; l++)
                        for (int jj = 0; jj < NR; jj++) {
                            S *dst = bp + l * W * NR + jj;
                            if (jr + jj < nc) {
                                K::put(dst, NR, coef(p0 + l, j0 + jr + jj));
                                if (mp)
                                    any |= (mp[l * NR + jj] = skip(p0 + l, j0 + jr + jj) ? 1 : 0);
                            } else {
                                K::pad(dst, NR);
                                if (mp)
                                    mp[l * NR + jj] = 1;
                            }
                        }
                    if (Mp)
                        anyskip[jr / NR] = any;
                }
                // A: micro-panels of MR rows, kc x MR each
                for (long ir = 0; ir < mc; ir += MR) {
                    S *ap = Ap + ir * W * kc;
                    for (long l = 0; l < kc; l++)
                        for (int ii = 0; ii < MR; ii++) {
                            S *dst = ap + l * W * MR + ii;
                            if (ir + ii < mc)
                                K::put(dst, MR, x(i0 + ir + ii, p0 + l));
                            else
                                K::pad(dst, MR);
                        }
                }
                for (long jr = 0; jr < nc; jr += NR) {
                    const int nr = (int)std::min((long)NR, nc - jr);
                    const S *bp = Bp + jr * W * kc;
                    const unsigned char *mp = Mp && anyskip[jr / NR] ? Mp + jr * kc : 0;
                    for (long ir = 0; ir < mc; ir += MR) {
                        const int mr = (int)std::min((long)MR, mc - ir);
                        const S *ap = Ap + ir * W * kc;
                        T *ct = C + (i0 + ir) + (j0 + jr) * ldc;
                        switch (mode) {
                        case BLK_ADD:
                            K::template micro<BLK_ADD>(kc, ap, bp, mp, ct, ldc, mr, nr);
                            break;
                        case BLK_SUB:
                            K::template micro<BLK_SUB>(kc, ap, bp, mp, ct, ldc, mr, nr);
                            break;
                        case BLK_DOT:
                            K::template micro<BLK_DOT>(kc, ap, bp, mp, ct, ldc, mr, nr);
                            break;
                        case BLK_ADD2:
                            K::template micro<BLK_ADD2>(kc, ap, bp, mp, ct, ldc, mr, nr);
                            break;
                        default:
                            K::template micro<BLK_SUBDOT>(kc, ap, bp, mp, ct, ldc, mr, nr);
                            break;
                        }
                    }
                }
            }
        }
    }
}

// ---- the kernel of each type ----

template <class T, class Enable = void> struct blk_kernel_of {
    // MPFR, GMP and other types with non-trivial copies
    typedef blk_pointer_kernel<T, 4, 4, 64, 128, 256> type;
};
template <class T> struct blk_kernel_of<T, typename std::enable_if<std::is_trivially_copyable<T>::value>::type> {
    // mplapack_binary80_t, mplapack_binary128_t and their complex types
    typedef blk_value_kernel<T, 4, 4, 64, 128, 256> type;
};
#if defined(MPLAPACK_BUILD_WITH_DOUBLE)
template <> struct blk_kernel_of<REAL, void> {
    typedef blk_value_kernel<REAL, 16, 4, 128, 256, 1024> type;
};
template <> struct blk_kernel_of<COMPLEX, void> {
    typedef blk_zplane_kernel<16, 4, 64, 256, 512> type;
};
#endif
#if defined(MPLAPACK_BUILD_WITH_DD)
template <> struct blk_kernel_of<dd_real, void> {
    typedef blk_plane_kernel<dd_real, 2, 4, 4, 64, 256, 512> type;
};
#elif defined(MPLAPACK_BUILD_WITH_TD)
template <> struct blk_kernel_of<td_real, void> {
    typedef blk_plane_kernel<td_real, 3, 4, 2, 64, 256, 512> type;
};
#elif defined(MPLAPACK_BUILD_WITH_QD)
template <> struct blk_kernel_of<qd_real, void> {
    typedef blk_plane_kernel<qd_real, 4, 4, 2, 64, 128, 512> type;
};
#endif

// Whether the engine is faster than the lanes for T and so is used by
// default: double (vectorized products and sums) and dd, td, qd (branch-free
// words).  For mplapack_binary80_t, mplapack_binary128_t, MPFR and GMP the
// arithmetic dominates
// and the engine measured slower (0.7 to 1.0 times); MPLAPACK_BLAS_BLOCK=2
// uses it for every type.
template <class T> struct blk_profitable {
    static const bool value = false;
};
#if defined(MPLAPACK_BUILD_WITH_DOUBLE)
template <> struct blk_profitable<REAL> {
    static const bool value = true;
};
template <> struct blk_profitable<COMPLEX> {
    static const bool value = true;
};
#endif
#if defined(MPLAPACK_BUILD_WITH_DD)
template <> struct blk_profitable<dd_real> {
    static const bool value = true;
};
#elif defined(MPLAPACK_BUILD_WITH_TD)
template <> struct blk_profitable<td_real> {
    static const bool value = true;
};
#elif defined(MPLAPACK_BUILD_WITH_QD)
template <> struct blk_profitable<qd_real> {
    static const bool value = true;
};
#endif

// The same for Rtrsm and Rsyrk: double only.  For dd, td and qd the
// column forms on planes of Rblas_simd_common.h measured as fast or faster
// (0.76 to 1.05 times on the engine, slower on the small calls of Rpotrf).
template <class T> struct blk_profitable_l3 {
    static const bool value = false;
};
#if defined(MPLAPACK_BUILD_WITH_DOUBLE)
template <> struct blk_profitable_l3<REAL> {
    static const bool value = true;
};
template <> struct blk_profitable_l3<COMPLEX> {
    static const bool value = true;
};
#endif

inline long blk_env(const char *name, long dflt)
{
    const char *v = std::getenv(name);
    return (v == 0 || *v == '\0') ? dflt : std::strtol(v, 0, 10);
}

// Whether the calls run in parallel.
inline bool blk_parallel()
{
#ifdef _OPENMP
    return blk_env("MPLAPACK_BLAS_OMP", 1) != 0 && omp_get_max_threads() >= 2;
#else
    return false;
#endif
}

// f(j) for j < n, in parallel when par, with the caller's context.
template <class Context, class F> void blk_columns(long n, bool par, F f)
{
    const Context ctx;
#ifdef _OPENMP
#pragma omp parallel if (par)
#endif
    {
        const typename Context::scope sc(ctx);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (long j = 0; j < n; j++)
            f(j);
    }
    (void)par;
}

// ---- Rgemm ----
// C := alpha*op(A)*op(B) + beta*C with the operations of Rgemm_{NN,NT,TN,TT}_omp
// (openmp/ of every type): C := beta*C (C := 0 when beta == 0), then
//   op(A) = A:   c(i,j) += temp*a(i,l), temp = alpha*op(b)(l,j)
//   op(A) = A':  temp = 0; temp += a(l,i)*op(b)(l,j) over l; c(i,j) += alpha*temp
// Returns false (the caller runs its own code) for m*n*k < min_mnk (default:
// MPLAPACK_BLAS_BLOCK_MIN, 32768) or with MPLAPACK_BLAS_BLOCK=0.  alpha != 0
// and m, n > 0 are assumed.
template <class Context, class T> bool blk_gemm(bool nota, bool notb, long m, long n, long k, const T &alpha, const T *A, long lda, const T *B, long ldb, const T &beta, T *C, long ldc, double min_mnk = -1.0)
{
    typedef typename blk_kernel_of<T>::type K;
    if (min_mnk < 0.0)
        min_mnk = (double)blk_env("MPLAPACK_BLAS_BLOCK_MIN", 32768);
    const long use = blk_env("MPLAPACK_BLAS_BLOCK", 1);
    if (use == 0 || (use == 1 && !blk_profitable<T>::value) || (double)m * n * k < min_mnk)
        return false;
    const bool par = blk_parallel();
    blk_columns<Context>(n, par, [&](long j) {
        if (beta == 0.0) {
            for (long i = 0; i < m; i++)
                C[i + j * ldc] = 0.0;
        } else if (beta != 1.0) {
            for (long i = 0; i < m; i++)
                C[i + j * ldc] = beta * C[i + j * ldc];
        }
    });
    if (k <= 0)
        return true;
    if (nota) {
        std::vector<T> temp((size_t)(k * n)); // temp = alpha*op(b)(l,j)
        blk_columns<Context>(n, par, [&](long j) {
            for (long l = 0; l < k; l++)
                temp[(size_t)(l + j * k)] = alpha * (notb ? B[l + j * ldb] : B[j + l * ldb]);
        });
        blk_update<K, Context>(BLK_ADD, m, n, k, C, ldc, [&](long i, long l) -> const T & { return A[i + l * lda]; }, [&](long l, long j) -> const T & { return temp[(size_t)(l + j * k)]; }, blk_no_skip(), false, par);
    } else {
        std::vector<T> temp((size_t)(m * n), T(0.0)); // temp = 0
        blk_update<K, Context>(BLK_DOT, m, n, k, temp.data(), m, [&](long i, long l) -> const T & { return A[l + i * lda]; }, [&](long l, long j) -> const T & { return notb ? B[l + j * ldb] : B[j + l * ldb]; }, blk_no_skip(), false, par);
        blk_columns<Context>(n, par, [&](long j) {
            for (long i = 0; i < m; i++)
                C[i + j * ldc] += alpha * temp[(size_t)(i + j * m)];
        });
    }
    return true;
}

// ---- Cgemm ----
// C := alpha*op(A)*op(B) + beta*C with the operations of the reference Cgemm
// (mpblas/reference): op(A) = A: c := beta*c (zero when beta == 0), then
// c(i,j) += temp*a(i,l), temp = alpha*op(b)(l,j); op(A) = A' or A^H:
// temp = zero; temp += op(a)(l,i)*op(b)(l,j); c = alpha*temp [+ beta*c].
// op(b)(l,j) is b(l,j), b(j,l) or conj(b(j,l)).  Same conditions as blk_gemm.
template <class Context, class T> bool blk_cgemm(bool nota, bool conja, bool notb, bool conjb, long m, long n, long k, const T &alpha, const T *A, long lda, const T *B, long ldb, const T &beta, T *C, long ldc)
{
    typedef typename blk_kernel_of<T>::type K;
    const long use = blk_env("MPLAPACK_BLAS_BLOCK", 1);
    if (use == 0 || (use == 1 && !blk_profitable<T>::value) || (double)m * n * k < (double)blk_env("MPLAPACK_BLAS_BLOCK_MIN", 32768))
        return false;
    const bool par = blk_parallel();
    const T zero(0.0, 0.0), one(1.0, 0.0);
    // op(b)(l,j)
    const auto opb = [&](long l, long j) -> T { return notb ? B[l + j * ldb] : (conjb ? T(conj(B[j + l * ldb])) : B[j + l * ldb]); };
    if (nota) {
        blk_columns<Context>(n, par, [&](long j) {
            if (beta == zero) {
                for (long i = 0; i < m; i++)
                    C[i + j * ldc] = zero;
            } else if (beta != one) {
                for (long i = 0; i < m; i++)
                    C[i + j * ldc] = beta * C[i + j * ldc];
            }
        });
        if (k <= 0)
            return true;
        std::vector<T> temp((size_t)(k * n)); // temp = alpha*op(b)(l,j)
        blk_columns<Context>(n, par, [&](long j) {
            for (long l = 0; l < k; l++)
                temp[(size_t)(l + j * k)] = alpha * opb(l, j);
        });
        blk_update<K, Context>(BLK_ADD, m, n, k, C, ldc, [&](long i, long l) -> const T & { return A[i + l * lda]; }, [&](long l, long j) -> const T & { return temp[(size_t)(l + j * k)]; }, blk_no_skip(), false, par);
        return true;
    }
    // op(a)(l,i) and op(b)(l,j) as matrices (conjugates copied; copies keep
    // the precision)
    std::vector<T> ac, bc;
    if (conja) {
        ac.reserve((size_t)(k * m));
        for (long i = 0; i < m; i++)
            for (long l = 0; l < k; l++)
                ac.push_back(T(conj(A[l + i * lda])));
    }
    if (!notb) {
        bc.reserve((size_t)(k * n));
        for (long j = 0; j < n; j++)
            for (long l = 0; l < k; l++)
                bc.push_back(opb(l, j));
    }
    std::vector<T> temp((size_t)(m * n), zero); // temp = zero
    blk_update<K, Context>(
        BLK_DOT, m, n, k, temp.data(), m, [&](long i, long l) -> const T & { return conja ? ac[(size_t)(l + i * k)] : A[l + i * lda]; }, [&](long l, long j) -> const T & { return notb ? B[l + j * ldb] : bc[(size_t)(l + j * k)]; }, blk_no_skip(), false, par);
    blk_columns<Context>(n, par, [&](long j) {
        for (long i = 0; i < m; i++) {
            if (beta == zero)
                C[i + j * ldc] = alpha * temp[(size_t)(i + j * m)];
            else
                C[i + j * ldc] = alpha * temp[(size_t)(i + j * m)] + beta * C[i + j * ldc];
        }
    });
    return true;
}

// f(i0, i1) for blocks of 64 rows of [0, m), in parallel when par.
template <class Context, class F> void blk_rows(long m, bool par, F f)
{
    const long nb = (m + 63) / 64;
    blk_columns<Context>(nb, par, [&](long u) { f(u * 64, std::min(m, u * 64 + 64)); });
}

// Whether blk_trsm takes the case: the solves whose updates by the solved
// rows (columns) of a block can be applied at once in the reference order.
// Left, lower, transa = T and right, lower, transa = N add the terms of
// the block before those of the solved rows (columns) and keep the lanes.
inline bool blk_trsm_case(bool left, bool upper, bool trans) { return left ? (!trans || upper) : (trans || upper); }

// ---- Rtrsm ----
// B := alpha*inv(op(A))*B or alpha*B*inv(op(A)) with the operations of the
// reference Rtrsm (mpblas/reference) in their order, for the cases of
// blk_trsm_case.  The triangular index range is split in two recursively:
// the half the reference solves first, then the terms of its rows (columns)
// in the other half by the engine, in the reference order, then the other
// half; ranges of at most NB are solved by the reference loops.  Every
// element gets its terms in the order of the reference.  alpha != 0, m, n > 0.
template <class Context, class T> struct blk_trsm_run {
    typedef typename blk_kernel_of<T>::type K;
    bool left, upper, trans, nounit, par;
    long m, n, lda, ldb, NB;
    const T *A;
    T *B;
    std::vector<unsigned char> nz; // left, transa = N: b(k,j) != 0 before the division
    std::vector<T> temp;           // left, transa = T: the temp of the reference
    std::vector<unsigned char> az; // right: a(i,j) == 0 (the comparison is not cheap for every type)
    const T &a(long i, long j) const { return A[i + j * lda]; }
    T &b(long i, long j) const { return B[i + j * ldb]; }

    // the reference loops on the range [p0, p1), terms of the range only
    void base(long p0, long p1)
    {
        const T one(1.0);
        if (left && !trans) {
            blk_columns<Context>(n, par, [&](long j) {
                for (long t = 0; t < p1 - p0; t++) {
                    const long k = upper ? p1 - 1 - t : p0 + t;
                    unsigned char &f = nz[(size_t)(k + j * m)];
                    f = b(k, j) != T(0.0);
                    if (!f)
                        continue;
                    if (nounit)
                        b(k, j) = b(k, j) / a(k, k);
                    const long i0 = upper ? p0 : k + 1, i1 = upper ? k : p1;
                    for (long i = i0; i < i1; i++)
                        b(i, j) = b(i, j) - b(k, j) * a(i, k);
                }
            });
        } else if (left) { // upper
            blk_columns<Context>(n, par, [&](long j) {
                for (long i = p0; i < p1; i++) {
                    T &tp = temp[(size_t)(i + j * m)];
                    for (long k = p0; k < i; k++)
                        tp = tp - a(k, i) * b(k, j);
                    if (nounit)
                        tp = tp / a(i, i);
                    b(i, j) = tp;
                }
            });
        } else if (!trans) { // upper
            blk_rows<Context>(m, par, [&](long r0, long r1) {
                for (long j = p0; j < p1; j++) {
                    for (long k = p0; k < j; k++)
                        if (a(k, j) != T(0.0))
                            for (long i = r0; i < r1; i++)
                                b(i, j) = b(i, j) - a(k, j) * b(i, k);
                    if (nounit) {
                        const T tmp = one / a(j, j);
                        for (long i = r0; i < r1; i++)
                            b(i, j) = tmp * b(i, j);
                    }
                }
            });
        } else {
            blk_rows<Context>(m, par, [&](long r0, long r1) {
                for (long t = 0; t < p1 - p0; t++) {
                    const long k = upper ? p1 - 1 - t : p0 + t;
                    if (nounit) {
                        const T tmp = one / a(k, k);
                        for (long i = r0; i < r1; i++)
                            b(i, k) = tmp * b(i, k);
                    }
                    const long c0 = upper ? p0 : k + 1, c1 = upper ? k : p1;
                    for (long j = c0; j < c1; j++)
                        if (a(j, k) != T(0.0))
                            for (long i = r0; i < r1; i++)
                                b(i, j) = b(i, j) - a(j, k) * b(i, k);
                }
            });
        }
    }

    // the terms of the solved range [s0, s1) in the range [d0, d1)
    void update(long s0, long s1, long d0, long d1)
    {
        const long ns = s1 - s0, nd = d1 - d0;
        // the reference order of the solved indices: descending for the
        // right-looking upper forms, ascending otherwise
        const bool desc = upper && (left ? !trans : trans);
        if (left && !trans)
            blk_update<K, Context>(
                BLK_SUB, nd, n, ns, B + d0, ldb, [&](long i, long t) -> const T & { return a(d0 + i, desc ? s1 - 1 - t : s0 + t); }, [&](long t, long j) -> const T & { return b(desc ? s1 - 1 - t : s0 + t, j); }, [&](long t, long j) { return !nz[(size_t)((desc ? s1 - 1 - t : s0 + t) + j * m)]; }, true, par);
        else if (left)
            blk_update<K, Context>(
                BLK_SUBDOT, nd, n, ns, temp.data() + d0, m, [&](long i, long t) -> const T & { return a(s0 + t, d0 + i); }, [&](long t, long j) -> const T & { return b(s0 + t, j); }, blk_no_skip(), false, par);
        else if (!trans)
            blk_update<K, Context>(
                BLK_SUB, m, nd, ns, B + d0 * ldb, ldb, [&](long i, long t) -> const T & { return b(i, s0 + t); }, [&](long t, long j) -> const T & { return a(s0 + t, d0 + j); }, [&](long t, long j) { return az[(size_t)((s0 + t) + (d0 + j) * n)] != 0; }, true, par);
        else
            blk_update<K, Context>(
                BLK_SUB, m, nd, ns, B + d0 * ldb, ldb, [&](long i, long t) -> const T & { return b(i, desc ? s1 - 1 - t : s0 + t); }, [&](long t, long j) -> const T & { return a(d0 + j, desc ? s1 - 1 - t : s0 + t); },
                [&](long t, long j) { return az[(size_t)((d0 + j) + (desc ? s1 - 1 - t : s0 + t) * n)] != 0; }, true, par);
    }

    void solve(long p0, long p1)
    {
        if (p1 - p0 <= NB) {
            base(p0, p1);
            return;
        }
        const long mid = p0 + (p1 - p0) / 2;
        // the reference solves the high indices first in the right-looking
        // upper forms
        if (upper && (left ? !trans : trans)) {
            solve(mid, p1);
            update(mid, p1, p0, mid);
            solve(p0, mid);
        } else {
            solve(p0, mid);
            update(p0, mid, mid, p1);
            solve(mid, p1);
        }
    }
};

// Whether blk_trsm takes the call.
template <class T> bool blk_trsm_takes(bool left, bool upper, bool trans, long m, long n)
{
    const long use = blk_env("MPLAPACK_BLAS_BLOCK", 1);
    const long na = left ? m : n;
    return blk_trsm_case(left, upper, trans) && use != 0 && (use != 1 || blk_profitable_l3<T>::value) && (double)m * n * na / 2 >= (double)blk_env("MPLAPACK_BLAS_BLOCK_MIN", 32768);
}

// Complex Ctrsm with transa = C: the reference uses conj(a) where Rtrsm
// uses a (conj(a) != 0 exactly when a != 0), so it is blk_trsm with
// transa = T on a conjugated copy of A; Ctrsm with N or T is blk_trsm.
template <class T> struct blk_conj_copy {
    std::vector<T> v;
    blk_conj_copy(const T *A, long lda, long na)
    {
        v.reserve((size_t)(na * na));
        for (long j = 0; j < na; j++)
            for (long i = 0; i < na; i++)
                v.push_back(T(conj(A[i + j * lda]))); // copies keep the precision
    }
};

template <class Context, class T> bool blk_trsm(bool left, bool upper, bool trans, bool nounit, long m, long n, const T &alpha, const T *A, long lda, T *B, long ldb)
{
    const long na = left ? m : n;
    if (!blk_trsm_takes<T>(left, upper, trans, m, n))
        return false;
    blk_trsm_run<Context, T> r;
    r.left = left, r.upper = upper, r.trans = trans, r.nounit = nounit, r.par = blk_parallel();
    r.m = m, r.n = n, r.lda = lda, r.ldb = ldb, r.A = A, r.B = B;
    r.NB = blk_env("MPLAPACK_BLAS_BLOCK_NB", blk_kernel_of<T>::type::W > 1 ? 16 : 64);
    const bool alpha_one = alpha == T(1.0);
    if (left && trans) {
        // temp = alpha*b(i,j)
        r.temp.resize((size_t)(m * n));
        blk_columns<Context>(n, r.par, [&](long j) {
            for (long i = 0; i < m; i++)
                r.temp[(size_t)(i + j * m)] = alpha * B[i + j * ldb];
        });
    } else if (!(left == false && trans) && !alpha_one) {
        // left or right with transa = N: b := alpha*b first
        blk_columns<Context>(n, r.par, [&](long j) {
            for (long i = 0; i < m; i++)
                B[i + j * ldb] = alpha * B[i + j * ldb];
        });
    }
    if (left && !trans)
        r.nz.resize((size_t)(m * n));
    if (!left) {
        r.az.resize((size_t)(n * n));
        blk_columns<Context>(n, r.par, [&](long j) {
            for (long i = 0; i < n; i++)
                r.az[(size_t)(i + j * n)] = A[i + j * lda] == T(0.0);
        });
    }
    r.solve(0, na);
    if (!left && trans && !alpha_one) {
        // the reference scales column k by alpha after its last use; it is
        // not changed afterwards, so the scaling can wait until the end
        blk_columns<Context>(n, r.par, [&](long j) {
            for (long i = 0; i < m; i++)
                B[i + j * ldb] = alpha * B[i + j * ldb];
        });
    }
    return true;
}

// ---- Rtrmm ----
// B := alpha*op(A)*B or alpha*B*op(A) with the operations of the reference
// Rtrmm in their order.  In every form the terms of an element come from the
// original B (B0) in a fixed order of k:
//   left, N:  b(i) = [alpha*b0(i)][*a(i,i)], then += (alpha*b0(k))*a(i,k) for
//             the k after i (upper: ascending k > i; lower: descending k < i),
//             the zero b0(k) left out;
//   left, T:  temp = b0(i)[*a(i,i)]; temp += a(k,i)*b0(k) (upper: k < i,
//             lower: k > i, ascending); b(i) = alpha*temp;
//   right:    column j scaled first, then += (alpha*a)*b0(:,k) over its k
//             (N upper: k < j asc., N lower: k > j asc., T upper: k > j asc.,
//             T lower: k < j desc.), the zero a left out.
// The left forms split the rows recursively (the half whose rows are done
// first, then the terms of the other half by the engine); the right forms
// take tiles of columns.  alpha != 0, m, n > 0.
template <class Context, class T> struct blk_trmm_run {
    typedef typename blk_kernel_of<T>::type K;
    bool upper, trans, nounit, par;
    long m, n, lda, ldb, NB;
    T alpha;
    const T *A;
    T *B;
    std::vector<T> tc;             // left N: alpha*b0(k,j); left T: temp
    std::vector<unsigned char> nz; // left N: b0(k,j) != 0
    const T &a(long i, long j) const { return A[i + j * lda]; }
    T &b(long i, long j) const { return B[i + j * ldb]; }

    // left: the terms of the rows [p0, p1) from the rows [p0, p1)
    void base(long p0, long p1)
    {
        blk_columns<Context>(n, par, [&](long j) {
            if (!trans) {
                for (long t = 0; t < p1 - p0; t++) {
                    const long k = upper ? p0 + t : p1 - 1 - t;
                    if (!nz[(size_t)(k + j * m)])
                        continue;
                    const T &tp = tc[(size_t)(k + j * m)];
                    if (upper) {
                        for (long i = p0; i < k; i++)
                            b(i, j) += tp * a(i, k);
                        b(k, j) = nounit ? T(tp * a(k, k)) : tp;
                    } else {
                        b(k, j) = tp;
                        if (nounit)
                            b(k, j) = b(k, j) * a(k, k);
                        for (long i = k + 1; i < p1; i++)
                            b(i, j) += tp * a(i, k);
                    }
                }
            } else {
                for (long i = p0; i < p1; i++) {
                    T &tp = tc[(size_t)(i + j * m)];
                    const long k0 = upper ? p0 : i + 1, k1 = upper ? i : p1;
                    for (long k = k0; k < k1; k++)
                        tp += a(k, i) * b(k, j);
                }
            }
        });
    }

    // left: the terms of the rows [d0, d1) from the rows [s0, s1)
    void update(long s0, long s1, long d0, long d1)
    {
        const long ns = s1 - s0, nd = d1 - d0;
        if (!trans) {
            // upper: ascending k; lower: descending k
            blk_update<K, Context>(
                BLK_ADD, nd, n, ns, B + d0, ldb, [&](long i, long t) -> const T & { return a(d0 + i, upper ? s0 + t : s1 - 1 - t); }, [&](long t, long j) -> const T & { return tc[(size_t)((upper ? s0 + t : s1 - 1 - t) + j * m)]; },
                [&](long t, long j) { return !nz[(size_t)((upper ? s0 + t : s1 - 1 - t) + j * m)]; }, true, par);
        } else {
            blk_update<K, Context>(
                BLK_DOT, nd, n, ns, tc.data() + d0, m, [&](long i, long t) -> const T & { return a(s0 + t, d0 + i); }, [&](long t, long j) -> const T & { return b(s0 + t, j); }, blk_no_skip(), false, par);
        }
    }

    void solve(long p0, long p1)
    {
        if (p1 - p0 <= NB) {
            base(p0, p1);
            return;
        }
        const long mid = p0 + (p1 - p0) / 2;
        // N: rows done first are the low (upper) or high (lower) ones and
        // receive the terms of the other half afterwards; T: the rows whose
        // first terms are the other half's (upper: high rows from low ones)
        if (!trans) {
            if (upper) {
                solve(p0, mid);
                update(mid, p1, p0, mid);
                solve(mid, p1);
            } else {
                solve(mid, p1);
                update(p0, mid, mid, p1);
                solve(p0, mid);
            }
        } else {
            if (upper) {
                solve(p0, mid);
                update(p0, mid, mid, p1);
                solve(mid, p1);
            } else {
                solve(p0, mid);
                update(mid, p1, p0, mid);
                solve(mid, p1);
            }
        }
    }
};

template <class T> bool blk_trmm_takes(long m, long n, bool left)
{
    const long use = blk_env("MPLAPACK_BLAS_BLOCK", 1);
    const long na = left ? m : n;
    return use != 0 && (use != 1 || blk_profitable_l3<T>::value) && (double)m * n * na / 2 >= (double)blk_env("MPLAPACK_BLAS_BLOCK_MIN", 32768);
}

template <class Context, class T> bool blk_trmm(bool left, bool upper, bool trans, bool nounit, long m, long n, const T &alpha, const T *A, long lda, T *B, long ldb)
{
    typedef typename blk_kernel_of<T>::type K;
    if (!blk_trmm_takes<T>(m, n, left))
        return false;
    const bool par = blk_parallel();
    const long NB = blk_env("MPLAPACK_BLAS_BLOCK_NB", blk_kernel_of<T>::type::W > 1 ? 16 : 64);
    const T one(1.0), zero(0.0);
    if (left) {
        blk_trmm_run<Context, T> r;
        r.upper = upper, r.trans = trans, r.nounit = nounit, r.par = par;
        r.m = m, r.n = n, r.lda = lda, r.ldb = ldb, r.NB = NB, r.A = A, r.B = B;
        r.tc.resize((size_t)(m * n));
        if (!trans) {
            r.nz.resize((size_t)(m * n));
            blk_columns<Context>(n, par, [&](long j) {
                for (long k = 0; k < m; k++) {
                    unsigned char &f = r.nz[(size_t)(k + j * m)];
                    f = B[k + j * ldb] != zero;
                    if (f)
                        r.tc[(size_t)(k + j * m)] = alpha * B[k + j * ldb]; // temp = alpha*b(k,j)
                }
            });
        } else {
            blk_columns<Context>(n, par, [&](long j) {
                for (long i = 0; i < m; i++) {
                    T &tp = r.tc[(size_t)(i + j * m)];
                    tp = B[i + j * ldb]; // temp = b(i,j); temp = temp*a(i,i)
                    if (nounit)
                        tp = tp * A[i + i * lda];
                }
            });
        }
        r.solve(0, m);
        if (trans)
            blk_columns<Context>(n, par, [&](long j) {
                for (long i = 0; i < m; i++)
                    B[i + j * ldb] = alpha * r.tc[(size_t)(i + j * m)];
            });
        return true;
    }
    // right: the original columns as sources; the zero tests of a once
    std::vector<T> B0;
    B0.reserve((size_t)(m * n));
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++)
            B0.push_back(B[i + j * ldb]); // copies keep the precision
    std::vector<unsigned char> az((size_t)(n * n));
    std::vector<T> co((size_t)(n * n)); // alpha*a(k,j) (N) or alpha*a(j,k) (T), as coef(k, j)
    blk_columns<Context>(n, par, [&](long j) {
        for (long k = 0; k < n; k++) {
            const T &akj = trans ? A[j + k * lda] : A[k + j * lda];
            az[(size_t)(k + j * n)] = akj == zero;
            if (!az[(size_t)(k + j * n)])
                co[(size_t)(k + j * n)] = alpha * akj;
        }
    });
    // the k of column j: ascending (or descending: T lower) over [k0, k1)
    const bool desc = trans && !upper;
    const auto krange = [&](long j, long &k0, long &k1) {
        const bool below = trans ? upper : !upper; // k > j
        k0 = below ? j + 1 : 0, k1 = below ? n : j;
    };
    for (long j0 = 0; j0 < n; j0 += NB) {
        const long j1 = std::min(n, j0 + NB), jb = j1 - j0;
        // scaling: N: b(:,j) = (alpha[*a(j,j)])*b(:,j); T: the same if != one
        blk_columns<Context>(jb, par, [&](long t) {
            const long j = j0 + t;
            T tp = alpha;
            if (nounit)
                tp = tp * A[j + j * lda];
            if (!trans || tp != one)
                for (long i = 0; i < m; i++)
                    B[i + j * ldb] = tp * B[i + j * ldb];
        });
        // the terms of the tile in their order: the k outside the tile that
        // come before the tile's own k, the tile's k, the k after
        const auto outside = [&](long c0, long c1) { // k in [c0, c1), outside the tile
            if (c1 <= c0)
                return;
            blk_update<K, Context>(
                BLK_ADD, m, jb, c1 - c0, B + j0 * ldb, ldb, [&](long i, long t) -> const T & { return B0[(size_t)(i + (desc ? c1 - 1 - t : c0 + t) * m)]; },
                [&](long t, long j) -> const T & { return co[(size_t)((desc ? c1 - 1 - t : c0 + t) + (j0 + j) * n)]; }, [&](long t, long j) { return az[(size_t)((desc ? c1 - 1 - t : c0 + t) + (j0 + j) * n)] != 0; }, true, par);
        };
        const auto inside = [&]() {
            blk_rows<Context>(m, par, [&](long r0, long r1) {
                for (long j = j0; j < j1; j++) {
                    long k0, k1;
                    krange(j, k0, k1);
                    k0 = std::max(k0, j0), k1 = std::min(k1, j1);
                    for (long t = 0; t < k1 - k0; t++) {
                        const long k = desc ? k1 - 1 - t : k0 + t;
                        if (az[(size_t)(k + j * n)])
                            continue;
                        const T &tp = co[(size_t)(k + j * n)];
                        for (long i = r0; i < r1; i++)
                            B[i + j * ldb] += tp * B0[(size_t)(i + k * m)];
                    }
                }
            });
        };
        const bool below = trans ? upper : !upper; // the k of column j are > j
        if (below) {
            inside();          // k in (j, j1)
            outside(j1, n);    // k >= j1, ascending
        } else if (!desc) {
            outside(0, j0);    // k < j0, ascending
            inside();          // k in [j0, j)
        } else {
            inside();          // k in [j0, j), descending
            outside(0, j0);    // k < j0, descending
        }
    }
    return true;
}

// ---- Rsymm ----
// C := alpha*A*B + beta*C or alpha*B*A + beta*C, A symmetric, with the
// operations of the reference Rsymm in their order:
//   left, upper (lower): temp2(i) = sum over k < i (k > i), ascending, of
//     b(k)*a(k,i); c(i) = [beta*c(i) +] temp1*a(i,i) + alpha*temp2(i) with
//     temp1 = alpha*b(i); then c(i) += (alpha*b(i'))*a(i,i') for the later
//     rows i' > i ascending (i' < i descending);
//   right: c(:,j) = [beta*c(:,j) +] (alpha*a(j,j))*b(:,j), then
//     c(:,j) += (alpha*a(k,j))*b(:,k) for k != j ascending.
// Rows (columns) in tiles of NB: the terms from outside the tile on the
// engine, those from inside by the reference loops.  m, n > 0, alpha != 0.
// HERM: Chemm, which uses conj(a) where the symmetric routine uses the
// reflected element in temp2 (left) and temp1 (right), and real(a(i,i)).
template <bool HERM, class T> inline T blk_conj_if(const T &x)
{
    if constexpr (HERM)
        return T(conj(x));
    else
        return x;
}
template <class Context, class T, bool HERM = false> bool blk_symm(bool left, bool upper, long m, long n, const T &alpha, const T *A, long lda, const T *B, long ldb, const T &beta, T *C, long ldc)
{
    typedef typename blk_kernel_of<T>::type K;
    const long use = blk_env("MPLAPACK_BLAS_BLOCK", 1), na = left ? m : n;
    if (use == 0 || (use == 1 && !blk_profitable_l3<T>::value) || (double)m * n * na < (double)blk_env("MPLAPACK_BLAS_BLOCK_MIN", 32768))
        return false;
    const bool par = blk_parallel();
    const long NB = K::MC;
    const T zero(0.0);
    const bool beta_zero = beta == zero;
#define MPLAPACK_BLK_A(i, j) A[(i) + (j) * lda]
#define MPLAPACK_BLK_B(i, j) B[(i) + (j) * ldb]
#define MPLAPACK_BLK_C(i, j) C[(i) + (j) * ldc]
    // the symmetric a(k,i) as stored
    const auto as = [&](long k, long i) -> const T & { return (upper ? k <= i : k >= i) ? MPLAPACK_BLK_A(k, i) : MPLAPACK_BLK_A(i, k); };
    // HERM: conj(a) of the stored triangle, for temp2 (left) or temp1 (right)
    std::vector<T> ah;
    if (HERM) {
        ah.reserve((size_t)(na * na));
        for (long j = 0; j < na; j++)
            for (long i = 0; i < na; i++)
                ah.push_back(blk_conj_if<HERM>(MPLAPACK_BLK_A(i, j))); // copies keep the precision
    }
    // the element of temp2 (left): a(k,i) as stored, conjugated for HERM
    const auto a2 = [&](long k, long i) -> const T & {
        const long r = (upper ? k <= i : k >= i) ? k : i, c = (upper ? k <= i : k >= i) ? i : k;
        return HERM ? ah[(size_t)(r + c * na)] : MPLAPACK_BLK_A(r, c);
    };
    if (left) {
        std::vector<T> t1((size_t)(m * n)), t2((size_t)(m * n), zero); // temp1 = alpha*b, temp2 = zero
        blk_columns<Context>(n, par, [&](long j) {
            for (long i = 0; i < m; i++)
                t1[(size_t)(i + j * m)] = alpha * MPLAPACK_BLK_B(i, j);
        });
        // temp2(i) += b(k)*a(k,i) over the k before (upper: k < i) or after
        // (lower: k > i) the row, ascending: coef b(k,j), x a(k,i)
        for (long i0 = 0; i0 < m; i0 += NB) {
            const long i1 = std::min(m, i0 + NB), ib = i1 - i0;
            const auto outside = [&](long k0, long k1) {
                if (k1 > k0)
                    blk_update<K, Context>(
                        BLK_ADD, ib, n, k1 - k0, t2.data() + i0, m, [&](long i, long t) -> const T & { return a2(k0 + t, i0 + i); }, [&](long t, long j) -> const T & { return MPLAPACK_BLK_B(k0 + t, j); }, blk_no_skip(), false, par);
            };
            const auto inside = [&]() {
                blk_columns<Context>(n, par, [&](long j) {
                    for (long i = i0; i < i1; i++) {
                        T &tp = t2[(size_t)(i + j * m)];
                        const long k0 = upper ? i0 : i + 1, k1 = upper ? i : i1;
                        for (long k = k0; k < k1; k++)
                            tp += MPLAPACK_BLK_B(k, j) * a2(k, i);
                    }
                });
            };
            if (upper) {
                outside(0, i0);
                inside();
            } else {
                inside();
                outside(i1, m);
            }
        }
        // c(i) = [beta*c(i) +] temp1*a(i,i) + alpha*temp2
        blk_columns<Context>(n, par, [&](long j) {
            for (long i = 0; i < m; i++) {
                T &cij = MPLAPACK_BLK_C(i, j);
                if constexpr (HERM) {
                    if (beta_zero)
                        cij = t1[(size_t)(i + j * m)] * MPLAPACK_BLK_A(i, i).real() + alpha * t2[(size_t)(i + j * m)];
                    else
                        cij = beta * cij + t1[(size_t)(i + j * m)] * MPLAPACK_BLK_A(i, i).real() + alpha * t2[(size_t)(i + j * m)];
                } else {
                    if (beta_zero)
                        cij = t1[(size_t)(i + j * m)] * MPLAPACK_BLK_A(i, i) + alpha * t2[(size_t)(i + j * m)];
                    else
                        cij = beta * cij + t1[(size_t)(i + j * m)] * MPLAPACK_BLK_A(i, i) + alpha * t2[(size_t)(i + j * m)];
                }
            }
        });
        // c(i) += temp1(i')*a(i,i') over the later rows: upper i' > i
        // ascending, lower i' < i descending
        for (long i0 = 0; i0 < m; i0 += NB) {
            const long i1 = std::min(m, i0 + NB), ib = i1 - i0;
            const auto outside = [&](long k0, long k1) {
                if (k1 > k0)
                    blk_update<K, Context>(
                        BLK_ADD, ib, n, k1 - k0, C + i0, ldc, [&](long i, long t) -> const T & { return as(i0 + i, upper ? k0 + t : k1 - 1 - t); }, [&](long t, long j) -> const T & { return t1[(size_t)((upper ? k0 + t : k1 - 1 - t) + j * m)]; }, blk_no_skip(), false, par);
            };
            const auto inside = [&]() {
                blk_columns<Context>(n, par, [&](long j) {
                    for (long i = i0; i < i1; i++) {
                        T &cij = MPLAPACK_BLK_C(i, j);
                        if (upper) {
                            for (long k = i + 1; k < i1; k++)
                                cij += t1[(size_t)(k + j * m)] * as(i, k);
                        } else {
                            for (long k = i - 1; k >= i0; k--)
                                cij += t1[(size_t)(k + j * m)] * as(i, k);
                        }
                    }
                });
            };
            if (upper) {
                inside();
                outside(i1, m);
            } else {
                inside();
                outside(0, i0);
            }
        }
    } else {
        // the coefficients alpha*a(k,j) (k != j)
        std::vector<T> co((size_t)(n * n));
        blk_columns<Context>(n, par, [&](long j) {
            // temp1 = alpha*a(k,j) as stored, alpha*conj(a(j,k)) for the
            // reflected element of HERM
            for (long k = 0; k < n; k++)
                if (k != j)
                    co[(size_t)(k + j * n)] = alpha * ((upper ? k <= j : k >= j) ? MPLAPACK_BLK_A(k, j) : a2(k, j));
            T t1;
            if constexpr (HERM)
                t1 = alpha * MPLAPACK_BLK_A(j, j).real();
            else
                t1 = alpha * MPLAPACK_BLK_A(j, j);
            for (long i = 0; i < m; i++) {
                T &cij = MPLAPACK_BLK_C(i, j);
                if (beta_zero)
                    cij = t1 * MPLAPACK_BLK_B(i, j);
                else
                    cij = beta * cij + t1 * MPLAPACK_BLK_B(i, j);
            }
        });
        for (long j0 = 0; j0 < n; j0 += NB) {
            const long j1 = std::min(n, j0 + NB), jb = j1 - j0;
            const auto outside = [&](long k0, long k1) {
                if (k1 > k0)
                    blk_update<K, Context>(
                        BLK_ADD, m, jb, k1 - k0, C + j0 * ldc, ldc, [&](long i, long t) -> const T & { return MPLAPACK_BLK_B(i, k0 + t); }, [&](long t, long j) -> const T & { return co[(size_t)((k0 + t) + (j0 + j) * n)]; }, blk_no_skip(), false, par);
            };
            outside(0, j0);
            blk_rows<Context>(m, par, [&](long r0, long r1) {
                for (long j = j0; j < j1; j++)
                    for (long k = j0; k < j1; k++) {
                        if (k == j)
                            continue;
                        const T &tp = co[(size_t)(k + j * n)];
                        for (long i = r0; i < r1; i++)
                            MPLAPACK_BLK_C(i, j) += tp * MPLAPACK_BLK_B(i, k);
                    }
            });
            outside(j1, n);
        }
    }
#undef MPLAPACK_BLK_A
#undef MPLAPACK_BLK_B
#undef MPLAPACK_BLK_C
    return true;
}

// ---- Rsyrk ----
// C := alpha*A*A' + beta*C or alpha*A'*A + beta*C on the triangle, with the
// operations of the reference Rsyrk in their order: the triangle in tiles of
// NB x NB, the tiles off the diagonal on the engine (one parallel loop), the
// diagonal tiles by the reference loops.  The terms of an element are all
// added in one place, in the reference order.  alpha != 0, n > 0.
template <class Context, class T> bool blk_syrk(bool upper, bool trans, long n, long k, const T &alpha, const T *A, long lda, const T &beta, T *C, long ldc)
{
    typedef typename blk_kernel_of<T>::type K;
    const long use = blk_env("MPLAPACK_BLAS_BLOCK", 1);
    if (use == 0 || (use == 1 && !blk_profitable_l3<T>::value) || (double)n * n * k / 2 < (double)blk_env("MPLAPACK_BLAS_BLOCK_MIN", 32768))
        return false;
    const bool par = blk_parallel();
    const long NB = K::MC, nb = (n + NB - 1) / NB;
    if (use == 1 && par && nb < 4) // too few tiles for the threads: the column forms
        return false;
    std::vector<T> temp;           // trans = N: alpha*a(j,l); T: the dot products
    std::vector<unsigned char> az; // trans = N: a(j,l) == 0
#define MPLAPACK_BLK_A(i, j) A[(i) + (j) * lda]
#define MPLAPACK_BLK_C(i, j) C[(i) + (j) * ldc]
    const auto off = [&](long bi, long bj) { return upper ? bi < bj : bi > bj; };
    // the diagonal tile d, columns [p0, p1) with p0 = d*NB, on the engine as
    // well: trans = N on a copy of the tile, of which the triangle is written
    // back; trans = T into temp, whose other elements are not used
    const auto diag = [&](long d) {
        const long p0 = d * NB, p1 = std::min(n, p0 + NB), w = p1 - p0;
        if (!trans) {
            std::vector<T> D;
            D.reserve((size_t)(w * w));
            for (long j = p0; j < p1; j++)
                for (long i = p0; i < p1; i++)
                    D.push_back(MPLAPACK_BLK_C(i, j)); // copies keep the precision
            blk_update<K, Context>(
                BLK_ADD, w, w, k, D.data(), w, [&](long i, long l) -> const T & { return MPLAPACK_BLK_A(p0 + i, l); }, [&](long l, long j) -> const T & { return temp[(size_t)(l + (p0 + j) * k)]; }, [&](long l, long j) { return az[(size_t)(l + (p0 + j) * k)] != 0; }, true, false);
            for (long j = p0; j < p1; j++) {
                const long i0 = upper ? p0 : j, i1 = upper ? j + 1 : p1;
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = D[(size_t)((i - p0) + (j - p0) * w)];
            }
        } else {
            blk_update<K, Context>(
                BLK_DOT, w, w, k, temp.data() + p0 + p0 * n, n, [&](long i, long l) -> const T & { return MPLAPACK_BLK_A(l, p0 + i); }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_A(l, p0 + j); }, blk_no_skip(), false, false);
        }
    };
    if (!trans) {
        // column j: c := beta*c (0 when beta == 0); for l: if a(j,l) != 0:
        // temp = alpha*a(j,l); c(i,j) += temp*a(i,l)
        temp.resize((size_t)(k * n));
        az.resize((size_t)(k * n));
        blk_columns<Context>(n, par, [&](long j) {
            const long i0 = upper ? 0 : j, i1 = upper ? j + 1 : n;
            if (beta == T(0.0)) {
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = T(0.0);
            } else if (beta != T(1.0)) {
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = beta * MPLAPACK_BLK_C(i, j);
            }
            for (long l = 0; l < k; l++) {
                az[(size_t)(l + j * k)] = MPLAPACK_BLK_A(j, l) == T(0.0);
                if (!az[(size_t)(l + j * k)])
                    temp[(size_t)(l + j * k)] = alpha * MPLAPACK_BLK_A(j, l);
            }
        });
        blk_update<K, Context>(
            BLK_ADD, n, n, k, C, ldc, [&](long i, long l) -> const T & { return MPLAPACK_BLK_A(i, l); }, [&](long l, long j) -> const T & { return temp[(size_t)(l + j * k)]; }, [&](long l, long j) { return az[(size_t)(l + j * k)] != 0; }, true, par, NB, off);
        blk_columns<Context>(nb, par, diag);
    } else {
        // c(i,j): temp = 0; temp += a(l,i)*a(l,j); c = alpha*temp [+ beta*c]
        temp.assign((size_t)(n * n), T(0.0));
        blk_update<K, Context>(
            BLK_DOT, n, n, k, temp.data(), n, [&](long i, long l) -> const T & { return MPLAPACK_BLK_A(l, i); }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_A(l, j); }, blk_no_skip(), false, par, NB, off);
        blk_columns<Context>(nb, par, diag);
        blk_columns<Context>(n, par, [&](long j) {
            const long i0 = upper ? 0 : j, i1 = upper ? j + 1 : n;
            for (long i = i0; i < i1; i++) {
                if (beta == T(0.0))
                    MPLAPACK_BLK_C(i, j) = alpha * temp[(size_t)(i + j * n)];
                else
                    MPLAPACK_BLK_C(i, j) = alpha * temp[(size_t)(i + j * n)] + beta * MPLAPACK_BLK_C(i, j);
            }
        });
    }
#undef MPLAPACK_BLK_A
#undef MPLAPACK_BLK_C
    return true;
}

// ---- Rsyr2k ----
// C := alpha*A*B' + alpha*B*A' + beta*C (N) or alpha*A'*B + alpha*B'*A +
// beta*C (T) on the triangle, with the operations of the reference Rsyr2k in
// their order: N: c := beta*c; for l with a(j,l) or b(j,l) nonzero:
// c(i,j) += a(i,l)*temp1 + b(i,l)*temp2 (temp1 = alpha*b(j,l), temp2 =
// alpha*a(j,l)), the pairs of BLK_ADD2; T: temp1 += a(l,i)*b(l,j), temp2 +=
// b(l,i)*a(l,j), c = [beta*c +] alpha*temp1 + alpha*temp2.  The tiles as in
// blk_syrk.  alpha != 0, n > 0, k > 0.
template <class Context, class T> bool blk_syr2k(bool upper, bool trans, long n, long k, const T &alpha, const T *A, long lda, const T *B, long ldb, const T &beta, T *C, long ldc)
{
    typedef typename blk_kernel_of<T>::type K;
    const long use = blk_env("MPLAPACK_BLAS_BLOCK", 1);
    if (k <= 0 || use == 0 || (use == 1 && !blk_profitable_l3<T>::value) || (double)n * n * k < (double)blk_env("MPLAPACK_BLAS_BLOCK_MIN", 32768))
        return false;
    const bool par = blk_parallel();
    const long NB = K::MC, nb = (n + NB - 1) / NB;
    if (use == 1 && par && nb < 4) // too few tiles for the threads: the column forms
        return false;
    const T zero(0.0), one(1.0);
#define MPLAPACK_BLK_A(i, j) A[(i) + (j) * lda]
#define MPLAPACK_BLK_B(i, j) B[(i) + (j) * ldb]
#define MPLAPACK_BLK_C(i, j) C[(i) + (j) * ldc]
    const auto off = [&](long bi, long bj) { return upper ? bi < bj : bi > bj; };
    std::vector<T> co, t1, t2;     // N: temp1, temp2 as pairs coef(2l, j), coef(2l+1, j); T: the sums
    std::vector<unsigned char> az; // N: a(j,l) == 0 and b(j,l) == 0
    // N: step l2 of the pairs: x(i, 2l) = a(i,l), x(i, 2l+1) = b(i,l)
    const auto xn = [&](long r0) { return [=](long i, long l2) -> const T & { return (l2 & 1) ? MPLAPACK_BLK_B(r0 + i, l2 >> 1) : MPLAPACK_BLK_A(r0 + i, l2 >> 1); }; };
    if (!trans) {
        co.resize((size_t)(2 * k * n));
        az.resize((size_t)(k * n));
        blk_columns<Context>(n, par, [&](long j) {
            const long i0 = upper ? 0 : j, i1 = upper ? j + 1 : n;
            if (beta == zero) {
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = zero;
            } else if (beta != one) {
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = beta * MPLAPACK_BLK_C(i, j);
            }
            for (long l = 0; l < k; l++) {
                az[(size_t)(l + j * k)] = MPLAPACK_BLK_A(j, l) == zero && MPLAPACK_BLK_B(j, l) == zero;
                if (!az[(size_t)(l + j * k)]) {
                    co[(size_t)(2 * l + j * 2 * k)] = alpha * MPLAPACK_BLK_B(j, l);
                    co[(size_t)(2 * l + 1 + j * 2 * k)] = alpha * MPLAPACK_BLK_A(j, l);
                }
            }
        });
        const auto coef = [&](long c0) { return [=, &co](long l2, long j) -> const T & { return co[(size_t)(l2 + (c0 + j) * 2 * k)]; }; };
        const auto skip = [&](long c0) { return [=, &az](long l2, long j) { return az[(size_t)((l2 >> 1) + (c0 + j) * k)] != 0; }; };
        blk_update<K, Context>(BLK_ADD2, n, n, 2 * k, C, ldc, xn(0), coef(0), skip(0), true, par, NB, off);
        blk_columns<Context>(nb, par, [&](long d) {
            const long p0 = d * NB, p1 = std::min(n, p0 + NB), w = p1 - p0;
            std::vector<T> D;
            D.reserve((size_t)(w * w));
            for (long j = p0; j < p1; j++)
                for (long i = p0; i < p1; i++)
                    D.push_back(MPLAPACK_BLK_C(i, j)); // copies keep the precision
            blk_update<K, Context>(BLK_ADD2, w, w, 2 * k, D.data(), w, xn(p0), coef(p0), skip(p0), true, false);
            for (long j = p0; j < p1; j++) {
                const long i0 = upper ? p0 : j, i1 = upper ? j + 1 : p1;
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = D[(size_t)((i - p0) + (j - p0) * w)];
            }
        });
    } else {
        t1.assign((size_t)(n * n), zero);
        t2.assign((size_t)(n * n), zero);
        const auto dots = [&](long r0, long c0, long rows, long cols, long tile, bool p) {
            blk_update<K, Context>(
                BLK_DOT, rows, cols, k, t1.data() + r0 + c0 * n, n, [&](long i, long l) -> const T & { return MPLAPACK_BLK_A(l, r0 + i); }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_B(l, c0 + j); }, blk_no_skip(), false, p, tile, off);
            blk_update<K, Context>(
                BLK_DOT, rows, cols, k, t2.data() + r0 + c0 * n, n, [&](long i, long l) -> const T & { return MPLAPACK_BLK_B(l, r0 + i); }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_A(l, c0 + j); }, blk_no_skip(), false, p, tile, off);
        };
        dots(0, 0, n, n, NB, par);
        blk_columns<Context>(nb, par, [&](long d) {
            const long p0 = d * NB, w = std::min(n, p0 + NB) - p0;
            blk_update<K, Context>(
                BLK_DOT, w, w, k, t1.data() + p0 + p0 * n, n, [&](long i, long l) -> const T & { return MPLAPACK_BLK_A(l, p0 + i); }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_B(l, p0 + j); }, blk_no_skip(), false, false);
            blk_update<K, Context>(
                BLK_DOT, w, w, k, t2.data() + p0 + p0 * n, n, [&](long i, long l) -> const T & { return MPLAPACK_BLK_B(l, p0 + i); }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_A(l, p0 + j); }, blk_no_skip(), false, false);
        });
        blk_columns<Context>(n, par, [&](long j) {
            const long i0 = upper ? 0 : j, i1 = upper ? j + 1 : n;
            for (long i = i0; i < i1; i++) {
                const T &u = t1[(size_t)(i + j * n)], &v = t2[(size_t)(i + j * n)];
                if (beta == zero)
                    MPLAPACK_BLK_C(i, j) = alpha * u + alpha * v;
                else
                    MPLAPACK_BLK_C(i, j) = beta * MPLAPACK_BLK_C(i, j) + alpha * u + alpha * v;
            }
        });
    }
#undef MPLAPACK_BLK_A
#undef MPLAPACK_BLK_B
#undef MPLAPACK_BLK_C
    return true;
}

// ---- Cher2k ----
// The reference Cher2k (beta real): as blk_syr2k off the diagonal with
// N: temp1 = alpha*conj(b(j,l)), temp2 = conj(alpha*a(j,l)); C: temp1 +=
// conj(a(l,i))*b(l,j), temp2 += conj(b(l,i))*a(l,j), c = [beta*c +]
// alpha*temp1 + conj(alpha)*temp2; the diagonal with the reference
// expressions on the real parts.  alpha != 0, n > 0, k > 0.
template <class Context, class T, class R> bool blk_her2k(bool upper, bool trans, long n, long k, const T &alpha, const T *A, long lda, const T *B, long ldb, const R &beta, T *C, long ldc)
{
    typedef typename blk_kernel_of<T>::type K;
    const long use = blk_env("MPLAPACK_BLAS_BLOCK", 1);
    if (k <= 0 || use == 0 || (use == 1 && !blk_profitable_l3<T>::value) || (double)n * n * k < (double)blk_env("MPLAPACK_BLAS_BLOCK_MIN", 32768))
        return false;
    const bool par = blk_parallel();
    const long NB = K::MC, nb = (n + NB - 1) / NB;
    if (use == 1 && par && nb < 4) // too few tiles for the threads: the column forms
        return false;
    const R rzero(0.0), rone(1.0);
    const T zero(0.0, 0.0);
#define MPLAPACK_BLK_A(i, j) A[(i) + (j) * lda]
#define MPLAPACK_BLK_B(i, j) B[(i) + (j) * ldb]
#define MPLAPACK_BLK_C(i, j) C[(i) + (j) * ldc]
    const auto off = [&](long bi, long bj) { return upper ? bi < bj : bi > bj; };
    if (!trans) {
        std::vector<T> co((size_t)(2 * k * n)); // temp1, temp2 as the pairs
        std::vector<unsigned char> az((size_t)(k * n));
        blk_columns<Context>(n, par, [&](long j) {
            const long i0 = upper ? 0 : j + 1, i1 = upper ? j : n; // off the diagonal
            T &cjj = MPLAPACK_BLK_C(j, j);
            if (beta == rzero) {
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = zero;
                cjj = zero;
            } else if (beta != rone) {
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = beta * MPLAPACK_BLK_C(i, j);
                cjj = beta * cjj.real();
            } else {
                cjj = cjj.real();
            }
            for (long l = 0; l < k; l++) {
                az[(size_t)(l + j * k)] = MPLAPACK_BLK_A(j, l) == zero && MPLAPACK_BLK_B(j, l) == zero;
                if (az[(size_t)(l + j * k)])
                    continue;
                T &t1 = co[(size_t)(2 * l + j * 2 * k)], &t2 = co[(size_t)(2 * l + 1 + j * 2 * k)];
                t1 = alpha * conj(MPLAPACK_BLK_B(j, l));
                t2 = conj(alpha * MPLAPACK_BLK_A(j, l));
                cjj = cjj.real() + (MPLAPACK_BLK_A(j, l) * t1 + MPLAPACK_BLK_B(j, l) * t2).real();
            }
        });
        const auto xn = [&](long r0) { return [=](long i, long l2) -> const T & { return (l2 & 1) ? MPLAPACK_BLK_B(r0 + i, l2 >> 1) : MPLAPACK_BLK_A(r0 + i, l2 >> 1); }; };
        const auto coef = [&](long c0) { return [=, &co](long l2, long j) -> const T & { return co[(size_t)(l2 + (c0 + j) * 2 * k)]; }; };
        const auto skip = [&](long c0) { return [=, &az](long l2, long j) { return az[(size_t)((l2 >> 1) + (c0 + j) * k)] != 0; }; };
        blk_update<K, Context>(BLK_ADD2, n, n, 2 * k, C, ldc, xn(0), coef(0), skip(0), true, par, NB, off);
        blk_columns<Context>(nb, par, [&](long d) {
            const long p0 = d * NB, p1 = std::min(n, p0 + NB), w = p1 - p0;
            std::vector<T> D;
            D.reserve((size_t)(w * w));
            for (long j = p0; j < p1; j++)
                for (long i = p0; i < p1; i++)
                    D.push_back(MPLAPACK_BLK_C(i, j)); // copies keep the precision
            blk_update<K, Context>(BLK_ADD2, w, w, 2 * k, D.data(), w, xn(p0), coef(p0), skip(p0), true, false);
            for (long j = p0; j < p1; j++) {
                const long i0 = upper ? p0 : j + 1, i1 = upper ? j : p1; // off the diagonal
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = D[(size_t)((i - p0) + (j - p0) * w)];
            }
        });
    } else {
        // conj(a), conj(b) as copies (they keep the precision)
        std::vector<T> ac, bc, t1((size_t)(n * n), zero), t2((size_t)(n * n), zero);
        ac.reserve((size_t)(k * n)), bc.reserve((size_t)(k * n));
        for (long i = 0; i < n; i++)
            for (long l = 0; l < k; l++) {
                ac.push_back(T(conj(MPLAPACK_BLK_A(l, i))));
                bc.push_back(T(conj(MPLAPACK_BLK_B(l, i))));
            }
        const auto dots = [&](long p0, long w, bool tiles, bool p) {
            const long rows = tiles ? n : w;
            if (tiles) {
                blk_update<K, Context>(
                    BLK_DOT, rows, rows, k, t1.data(), n, [&](long i, long l) -> const T & { return ac[(size_t)(l + i * k)]; }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_B(l, j); }, blk_no_skip(), false, p, NB, off);
                blk_update<K, Context>(
                    BLK_DOT, rows, rows, k, t2.data(), n, [&](long i, long l) -> const T & { return bc[(size_t)(l + i * k)]; }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_A(l, j); }, blk_no_skip(), false, p, NB, off);
            } else {
                blk_update<K, Context>(
                    BLK_DOT, w, w, k, t1.data() + p0 + p0 * n, n, [&](long i, long l) -> const T & { return ac[(size_t)(l + (p0 + i) * k)]; }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_B(l, p0 + j); }, blk_no_skip(), false, p);
                blk_update<K, Context>(
                    BLK_DOT, w, w, k, t2.data() + p0 + p0 * n, n, [&](long i, long l) -> const T & { return bc[(size_t)(l + (p0 + i) * k)]; }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_A(l, p0 + j); }, blk_no_skip(), false, p);
            }
        };
        dots(0, n, true, par);
        blk_columns<Context>(nb, par, [&](long d) { dots(d * NB, std::min(n, d * NB + NB) - d * NB, false, false); });
        blk_columns<Context>(n, par, [&](long j) {
            const long i0 = upper ? 0 : j, i1 = upper ? j + 1 : n;
            for (long i = i0; i < i1; i++) {
                const T &u = t1[(size_t)(i + j * n)], &v = t2[(size_t)(i + j * n)];
                T &cij = MPLAPACK_BLK_C(i, j);
                if (i == j) {
                    if (beta == rzero)
                        cij = (alpha * u + conj(alpha) * v).real();
                    else
                        cij = beta * cij.real() + (alpha * u + conj(alpha) * v).real();
                } else {
                    if (beta == rzero)
                        cij = alpha * u + conj(alpha) * v;
                    else
                        cij = beta * cij + alpha * u + conj(alpha) * v;
                }
            }
        });
    }
#undef MPLAPACK_BLK_A
#undef MPLAPACK_BLK_B
#undef MPLAPACK_BLK_C
    return true;
}

// ---- Cherk ----
// C := alpha*A*A^H + beta*C or alpha*A^H*A + beta*C (alpha, beta real) on
// the triangle, with the operations of the reference Cherk in their order:
// as blk_syrk for the elements off the diagonal (trans = N: temp =
// alpha*conj(a(j,l)), c(i,j) += temp*a(i,l); C: temp += conj(a(l,i))*a(l,j)),
// and the reference expressions for the diagonal (N: c(j,j) =
// real(c(j,j)) + real(temp*a(j,l)); C: rtemp += real(conj(a(l,j))*a(l,j))).
// alpha != 0, n > 0, k > 0.
template <class Context, class T, class R> bool blk_herk(bool upper, bool trans, long n, long k, const R &alpha, const T *A, long lda, const R &beta, T *C, long ldc)
{
    typedef typename blk_kernel_of<T>::type K;
    const long use = blk_env("MPLAPACK_BLAS_BLOCK", 1);
    if (k <= 0 || use == 0 || (use == 1 && !blk_profitable_l3<T>::value) || (double)n * n * k / 2 < (double)blk_env("MPLAPACK_BLAS_BLOCK_MIN", 32768))
        return false;
    const bool par = blk_parallel();
    const long NB = K::MC, nb = (n + NB - 1) / NB;
    if (use == 1 && par && nb < 4) // too few tiles for the threads: the column forms
        return false;
    const R zero(0.0), one(1.0);
    const T czero(0.0, 0.0);
    std::vector<T> temp, ac;       // N: alpha*conj(a(j,l)); C: the dot products, conj(a(l,i))
    std::vector<unsigned char> az; // N: a(j,l) == 0
#define MPLAPACK_BLK_A(i, j) A[(i) + (j) * lda]
#define MPLAPACK_BLK_C(i, j) C[(i) + (j) * ldc]
    const auto off = [&](long bi, long bj) { return upper ? bi < bj : bi > bj; };
    // the elements off the diagonal of the diagonal tile d
    const auto diag = [&](long d) {
        const long p0 = d * NB, p1 = std::min(n, p0 + NB), w = p1 - p0;
        if (!trans) {
            std::vector<T> D;
            D.reserve((size_t)(w * w));
            for (long j = p0; j < p1; j++)
                for (long i = p0; i < p1; i++)
                    D.push_back(MPLAPACK_BLK_C(i, j)); // copies keep the precision
            blk_update<K, Context>(
                BLK_ADD, w, w, k, D.data(), w, [&](long i, long l) -> const T & { return MPLAPACK_BLK_A(p0 + i, l); }, [&](long l, long j) -> const T & { return temp[(size_t)(l + (p0 + j) * k)]; }, [&](long l, long j) { return az[(size_t)(l + (p0 + j) * k)] != 0; }, true, false);
            for (long j = p0; j < p1; j++) {
                const long i0 = upper ? p0 : j + 1, i1 = upper ? j : p1;
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = D[(size_t)((i - p0) + (j - p0) * w)];
            }
        } else {
            blk_update<K, Context>(
                BLK_DOT, w, w, k, temp.data() + p0 + p0 * n, n, [&](long i, long l) -> const T & { return ac[(size_t)(l + (p0 + i) * k)]; }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_A(l, p0 + j); }, blk_no_skip(), false, false);
        }
    };
    if (!trans) {
        temp.resize((size_t)(k * n));
        az.resize((size_t)(k * n));
        blk_columns<Context>(n, par, [&](long j) {
            const long i0 = upper ? 0 : j + 1, i1 = upper ? j : n; // off the diagonal
            T &cjj = MPLAPACK_BLK_C(j, j);
            if (beta == zero) {
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = czero;
                cjj = czero;
            } else if (beta != one) {
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = beta * MPLAPACK_BLK_C(i, j);
                cjj = beta * cjj.real();
            } else {
                cjj = cjj.real();
            }
            for (long l = 0; l < k; l++) {
                az[(size_t)(l + j * k)] = MPLAPACK_BLK_A(j, l) == czero;
                if (!az[(size_t)(l + j * k)])
                    temp[(size_t)(l + j * k)] = alpha * conj(MPLAPACK_BLK_A(j, l));
            }
            // the diagonal element
            for (long l = 0; l < k; l++)
                if (!az[(size_t)(l + j * k)])
                    cjj = cjj.real() + (temp[(size_t)(l + j * k)] * MPLAPACK_BLK_A(j, l)).real();
        });
        blk_update<K, Context>(
            BLK_ADD, n, n, k, C, ldc, [&](long i, long l) -> const T & { return MPLAPACK_BLK_A(i, l); }, [&](long l, long j) -> const T & { return temp[(size_t)(l + j * k)]; }, [&](long l, long j) { return az[(size_t)(l + j * k)] != 0; }, true, par, NB, off);
        blk_columns<Context>(nb, par, diag);
    } else {
        temp.assign((size_t)(n * n), czero); // temp = zero
        ac.reserve((size_t)(k * n));
        for (long i = 0; i < n; i++)
            for (long l = 0; l < k; l++)
                ac.push_back(T(conj(MPLAPACK_BLK_A(l, i))));
        blk_update<K, Context>(
            BLK_DOT, n, n, k, temp.data(), n, [&](long i, long l) -> const T & { return ac[(size_t)(l + i * k)]; }, [&](long l, long j) -> const T & { return MPLAPACK_BLK_A(l, j); }, blk_no_skip(), false, par, NB, off);
        blk_columns<Context>(nb, par, diag);
        blk_columns<Context>(n, par, [&](long j) {
            const long i0 = upper ? 0 : j + 1, i1 = upper ? j : n;
            for (long i = i0; i < i1; i++) {
                if (beta == zero)
                    MPLAPACK_BLK_C(i, j) = alpha * temp[(size_t)(i + j * n)];
                else
                    MPLAPACK_BLK_C(i, j) = alpha * temp[(size_t)(i + j * n)] + beta * MPLAPACK_BLK_C(i, j);
            }
            R rtemp = zero;
            for (long l = 0; l < k; l++)
                rtemp += (conj(MPLAPACK_BLK_A(l, j)) * MPLAPACK_BLK_A(l, j)).real();
            T &cjj = MPLAPACK_BLK_C(j, j);
            if (beta == zero)
                cjj = alpha * rtemp;
            else
                cjj = alpha * rtemp + beta * cjj.real();
        });
    }
#undef MPLAPACK_BLK_A
#undef MPLAPACK_BLK_C
    return true;
}

} // namespace mplapack_blas

#endif
