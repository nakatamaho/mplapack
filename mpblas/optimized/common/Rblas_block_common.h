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
 *   value kernel (T trivially copyable: double, long double, __float128):
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
#include <type_traits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Rgemm_simd_common.h"
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "Rgemm_qd_ops_common.h"
#endif

namespace mplapack_blas {

enum { BLK_ADD = 0, BLK_SUB = 1, BLK_DOT = 2, BLK_SUBDOT = 3 };

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
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#define MPLAPACK_BLK_NO_CONTRACT __attribute__((optimize("fp-contract=off")))
#else
#define MPLAPACK_BLK_NO_CONTRACT
#endif

template <class T, int MR, int NR, int MODE, bool HM> MPLAPACK_GEMM_CLONES MPLAPACK_BLK_NO_CONTRACT void blk_value_micro(long kc, const T *A, const T *B, const unsigned char *mask, T *C, long ldc, int mr, int nr)
{
    T acc[NR][MR];
    for (int jj = 0; jj < NR; jj++)
        for (int ii = 0; ii < MR; ii++)
            acc[jj][ii] = (ii < mr && jj < nr) ? C[ii + jj * ldc] : T(0.0);
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
            MPLAPACK_BLK_PLANE_STEP(W, MPLAPACK_BLK_##T##_DECL, MPLAPACK_BLK_##T##_LOAD, MPLAPACK_BLK_##T##_STORE,       \
                                    MPLAPACK_BLK_##T##_MUL, MPLAPACK_BLK_##T##_ADD, MPLAPACK_BLK_##T##_NEG)              \
        }                                                                                                                \
    };
MPLAPACK_BLK_PLANE_STEPS(2, DD)
MPLAPACK_BLK_PLANE_STEPS(3, TD)
MPLAPACK_BLK_PLANE_STEPS(4, QD)

template <int W, int MR, int NR, int MODE, bool HM> MPLAPACK_GEMM_CLONES void blk_plane_micro(long kc, const double *A, const double *B, const unsigned char *mask, double *C, long ldc, int mr, int nr)
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
    // long double, __float128, complex types
    typedef blk_value_kernel<T, 4, 4, 64, 128, 256> type;
};
template <> struct blk_kernel_of<double, void> {
    typedef blk_value_kernel<double, 16, 4, 128, 256, 1024> type;
};
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
// words).  For long double, __float128, MPFR and GMP the arithmetic dominates
// and the engine measured slower (0.7 to 1.0 times); MPLAPACK_BLAS_BLOCK=2
// uses it for every type.
template <class T> struct blk_profitable {
    static const bool value = false;
};
template <> struct blk_profitable<double> {
    static const bool value = true;
};
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
template <> struct blk_profitable_l3<double> {
    static const bool value = true;
};

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
                    f = b(k, j) != 0.0;
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
                        if (a(k, j) != 0.0)
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
                        if (a(j, k) != 0.0)
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

template <class Context, class T> bool blk_trsm(bool left, bool upper, bool trans, bool nounit, long m, long n, const T &alpha, const T *A, long lda, T *B, long ldb)
{
    const long use = blk_env("MPLAPACK_BLAS_BLOCK", 1);
    const long na = left ? m : n;
    if (!blk_trsm_case(left, upper, trans) || use == 0 || (use == 1 && !blk_profitable_l3<T>::value) || (double)m * n * na / 2 < (double)blk_env("MPLAPACK_BLAS_BLOCK_MIN", 32768))
        return false;
    blk_trsm_run<Context, T> r;
    r.left = left, r.upper = upper, r.trans = trans, r.nounit = nounit, r.par = blk_parallel();
    r.m = m, r.n = n, r.lda = lda, r.ldb = ldb, r.A = A, r.B = B;
    r.NB = blk_env("MPLAPACK_BLAS_BLOCK_NB", blk_kernel_of<T>::type::W > 1 ? 16 : 64);
    const bool alpha_one = alpha == 1.0;
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
                r.az[(size_t)(i + j * n)] = A[i + j * lda] == 0.0;
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
            if (beta == 0.0) {
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = 0.0;
            } else if (beta != 1.0) {
                for (long i = i0; i < i1; i++)
                    MPLAPACK_BLK_C(i, j) = beta * MPLAPACK_BLK_C(i, j);
            }
            for (long l = 0; l < k; l++) {
                az[(size_t)(l + j * k)] = MPLAPACK_BLK_A(j, l) == 0.0;
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
                if (beta == 0.0)
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

} // namespace mplapack_blas

#endif
