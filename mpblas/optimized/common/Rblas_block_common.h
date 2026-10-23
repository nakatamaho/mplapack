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

enum { BLK_ADD = 0, BLK_SUB = 1, BLK_DOT = 2 };

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

template <class T, int MR, int NR, int MODE> MPLAPACK_GEMM_CLONES MPLAPACK_BLK_NO_CONTRACT void blk_value_micro(long kc, const T *A, const T *B, const unsigned char *mask, T *C, long ldc, int mr, int nr)
{
    T acc[NR][MR];
    for (int jj = 0; jj < NR; jj++)
        for (int ii = 0; ii < MR; ii++)
            acc[jj][ii] = (ii < mr && jj < nr) ? C[ii + jj * ldc] : T(0.0);
    for (long l = 0; l < kc; l++) {
        const T *a = A + l * MR, *b = B + l * NR;
        for (int jj = 0; jj < NR; jj++) {
            if (mask && mask[l * NR + jj])
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
            } else {
                MPLAPACK_GEMM_SIMD
                for (int ii = 0; ii < MR; ii++)
                    acc[jj][ii] += a[ii] * bv;
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
    template <int MODE> static void micro(long kc, const slot *A, const slot *B, const unsigned char *mask, T *C, long ldc, int mr, int nr) { blk_value_micro<T, MR, NR, MODE>(kc, A, B, mask, C, ldc, mr, nr); }
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
                    else
                        c[ii] += *a[ii] * bv;
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
            if (mask && mask[l * NR + jj])                                                                               \
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
                if (MODE == BLK_DOT)                                                                                     \
                    MUL(p, av, bv);                                                                                      \
                else                                                                                                     \
                    MUL(p, bv, av);                                                                                      \
                if (MODE == BLK_SUB)                                                                                     \
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

template <int W, int MR, int NR, int MODE> struct blk_plane_steps;
#define MPLAPACK_BLK_PLANE_STEPS(W, T)                                                                                   \
    template <int MR, int NR, int MODE> struct blk_plane_steps<W, MR, NR, MODE> {                                        \
        static MPLAPACK_GEMM_HD void run(long kc, const double *A, const double *B, const unsigned char *mask, double (*acc)[NR][MR]) \
        {                                                                                                                \
            MPLAPACK_BLK_PLANE_STEP(W, MPLAPACK_BLK_##T##_DECL, MPLAPACK_BLK_##T##_LOAD, MPLAPACK_BLK_##T##_STORE,       \
                                    MPLAPACK_BLK_##T##_MUL, MPLAPACK_BLK_##T##_ADD, MPLAPACK_BLK_##T##_NEG)              \
        }                                                                                                                \
    };
MPLAPACK_BLK_PLANE_STEPS(2, DD)
MPLAPACK_BLK_PLANE_STEPS(3, TD)
MPLAPACK_BLK_PLANE_STEPS(4, QD)

template <int W, int MR, int NR, int MODE> MPLAPACK_GEMM_CLONES void blk_plane_micro(long kc, const double *A, const double *B, const unsigned char *mask, double *C, long ldc, int mr, int nr)
{
    double acc[W][NR][MR];
    for (int jj = 0; jj < NR; jj++)
        for (int ii = 0; ii < MR; ii++)
            for (int q = 0; q < W; q++)
                acc[q][jj][ii] = (ii < mr && jj < nr) ? C[(ii + jj * ldc) * W + q] : 0.0;
    blk_plane_steps<W, MR, NR, MODE>::run(kc, A, B, mask, acc);
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
    template <int MODE> static void micro(long kc, const slot *A, const slot *B, const unsigned char *mask, T *C, long ldc, int mr, int nr) { blk_plane_micro<W, MR, NR, MODE>(kc, A, B, mask, reinterpret_cast<double *>(C), ldc, mr, nr); }
};

#endif

// ---- the engine ----

struct blk_no_skip {
    bool operator()(long, long) const { return false; }
};

// Context: the arithmetic context of the caller, installed in every worker
// (blas_thread_context of Rblas_cpu_common.h).  x(i, l) and coef(l, j) return
// references that stay valid during the call.
template <class K, class Context, class T, class FX, class FB, class FS> void blk_update(int mode, long m, long n, long k, T *C, long ldc, FX x, FB coef, FS skip, bool use_skip, bool parallel)
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
    long MC = K::MC, NC = K::NC;
    const long KC = K::KC;
    const long nbi = (m + MC - 1) / MC;
    // narrower column blocks when there are too few blocks for the threads
    while (NC > NR && nbi * ((n + NC - 1) / NC) < 2L * nthreads)
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
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (long blk = 0; blk < nbi * nbj; blk++) {
            const long i0 = (blk % nbi) * MC, j0 = (blk / nbi) * NC;
            const long mc = std::min(MC, m - i0), nc = std::min(NC, n - j0);
            for (long p0 = 0; p0 < k; p0 += KC) {
                const long kc = std::min(KC, k - p0);
                // B: micro-panels of NR columns, kc x NR each
                for (long jr = 0; jr < nc; jr += NR) {
                    S *bp = Bp + jr * W * kc;
                    unsigned char *mp = Mp ? Mp + jr * kc : 0;
                    for (long l = 0; l < kc; l++)
                        for (int jj = 0; jj < NR; jj++) {
                            S *dst = bp + l * W * NR + jj;
                            if (jr + jj < nc) {
                                K::put(dst, NR, coef(p0 + l, j0 + jr + jj));
                                if (mp)
                                    mp[l * NR + jj] = skip(p0 + l, j0 + jr + jj) ? 1 : 0;
                            } else {
                                K::pad(dst, NR);
                                if (mp)
                                    mp[l * NR + jj] = 1;
                            }
                        }
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
                    const unsigned char *mp = Mp ? Mp + jr * kc : 0;
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
                        default:
                            K::template micro<BLK_DOT>(kc, ap, bp, mp, ct, ldc, mr, nr);
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

} // namespace mplapack_blas

#endif
