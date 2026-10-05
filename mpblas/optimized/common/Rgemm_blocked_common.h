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
 * Cache-blocked, OpenMP-parallel Rgemm engine shared by the dd, qd and
 * binary128 optimized libraries.
 *
 * The engine performs, for every element of C, exactly the operations of
 * openmp/Rgemm_{NN,NT,TN,TT}_omp.cpp in the same order, so its result is
 * bit-identical to theirs:
 *
 *   C := beta*C  (C := 0 when beta == 0, untouched when beta == 1), then
 *   op(A) == A  (NN, NT):  for l = 0..k-1:  c = c + (alpha*op(B)(l,j)) * A(i,l)
 *   op(A) == A' (TN, TT):  t = 0;  for l = 0..k-1:  t = t + A(l,i) * op(B)(l,j);
 *                          c = c + alpha*t
 *
 * Only the loops over i and j are blocked and reordered; the chain of
 * operations on one element of C is never split or reassociated.  A block of
 * C (MC x NC) is owned by one thread; for every block of k (KC) the panels of
 * op(A) and of op(B) (or alpha*op(B)) are packed, plane by plane (structure of
 * arrays), and a kernel updates the accumulators of the block.
 *
 * Kernel interface (K):
 *   typedef ... value_type;
 *   static const int planes;               // doubles per packed element
 *   static const long MC, NC, KC;          // block sizes
 *   static void pack(double *p, long stride, long e, const value_type &x);
 *   static value_type unpack(const double *p, long stride, long e);
 *   // acc(i, j) = acc(i, j) + b(l, j) * a(l, i)   (axpy == true)
 *   // acc(i, j) = acc(i, j) + a(l, i) * b(l, j)   (axpy == false)
 *   // for j < nc, l < kc, i < mc, in the order j, l, i; plane p of
 *   // acc(i, j) is acc[p * as + i + j * mc], of a(l, i) a[p * ps + i + l * mc],
 *   // of b(l, j) b[p * bs + l + j * kc].
 *   static void block(bool axpy, long mc, long nc, long kc, double *acc, long as,
 *                     const double *a, long ps, const double *b, long bs);
 */

#ifndef _MPLAPACK_RGEMM_BLOCKED_COMMON_H_
#define _MPLAPACK_RGEMM_BLOCKED_COMMON_H_

#include <cstdlib>
#include <cstring>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mplapack_gemm {

// C (m x n) := alpha*op(A)*op(B) + beta*C, column major.
struct shape {
    bool nota, notb;
    long m, n, k, lda, ldb, ldc;
};

// Integer value of an environment variable, or dflt when unset or empty.
inline long env_long(const char *name, long dflt) {
    const char *v = std::getenv(name);
    if (v == NULL || *v == '\0')
        return dflt;
    return std::strtol(v, NULL, 10);
}

template <class T> inline const T &opA(const shape &s, const T *A, long i, long l) { return s.nota ? A[i + l * s.lda] : A[l + i * s.lda]; }
template <class T> inline const T &opB(const shape &s, const T *B, long l, long j) { return s.notb ? B[l + j * s.ldb] : B[j + l * s.ldb]; }

// The value C(i,j) takes before the products are added (the first loop of
// Rgemm_*_omp).
template <class T> inline T scaled_c(const T &beta, const T &c) {
    if (beta == 0.0)
        return T(0.0);
    if (beta != 1.0)
        return beta * c;
    return c;
}

template <class K> void gemm_blocked(const shape &s, const typename K::value_type &alpha, const typename K::value_type &beta, const typename K::value_type *A, const typename K::value_type *B, typename K::value_type *C) {
    typedef typename K::value_type T;
    const int P = K::planes;
    const bool axpy = s.nota;
    const long m = s.m, n = s.n, k = s.k;
    if (m <= 0 || n <= 0)
        return;

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    // Smaller column blocks when there are too few blocks to keep every
    // thread busy.
    long MC = K::MC, NC = K::NC;
    const long nbi = (m + MC - 1) / MC;
    while (NC > 4 && nbi * ((n + NC - 1) / NC) < 2L * nthreads)
        NC /= 2;
    const long nbj = (n + NC - 1) / NC;
    const long KC = K::KC;
    const long kcmax = k < KC ? (k > 0 ? k : 1) : KC;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector<double> accbuf((size_t)P * MC * NC), abuf((size_t)P * MC * kcmax), bbuf((size_t)P * kcmax * NC);
        double *acc = accbuf.data(), *ap = abuf.data(), *bp = bbuf.data();
        const T zero(0.0);
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
        for (long blk = 0; blk < nbi * nbj; blk++) {
            const long bi = blk % nbi, bj = blk / nbi;
            const long i0 = bi * MC, j0 = bj * NC;
            const long mc = (m - i0 < MC) ? m - i0 : MC;
            const long nc = (n - j0 < NC) ? n - j0 : NC;
            const long as = mc * nc;

            for (long j = 0; j < nc; j++)
                for (long i = 0; i < mc; i++)
                    K::pack(acc, as, i + j * mc, axpy ? scaled_c(beta, C[(i0 + i) + (j0 + j) * s.ldc]) : zero);

            for (long p0 = 0; p0 < k; p0 += KC) {
                const long kc = (k - p0 < KC) ? k - p0 : KC;
                const long ps = mc * kc, bs = kc * nc;
                for (long l = 0; l < kc; l++)
                    for (long i = 0; i < mc; i++)
                        K::pack(ap, ps, i + l * mc, opA(s, A, i0 + i, p0 + l));
                for (long j = 0; j < nc; j++)
                    for (long l = 0; l < kc; l++) {
                        if (axpy)
                            K::pack(bp, bs, l + j * kc, alpha * opB(s, B, p0 + l, j0 + j));
                        else
                            K::pack(bp, bs, l + j * kc, opB(s, B, p0 + l, j0 + j));
                    }
                K::block(axpy, mc, nc, kc, acc, as, ap, ps, bp, bs);
            }

            for (long j = 0; j < nc; j++)
                for (long i = 0; i < mc; i++) {
                    T &c = C[(i0 + i) + (j0 + j) * s.ldc];
                    if (axpy) {
                        c = K::unpack(acc, as, i + j * mc);
                    } else {
                        T t = scaled_c(beta, c);
                        t += alpha * K::unpack(acc, as, i + j * mc);
                        c = t;
                    }
                }
        }
    }
}

// Kernel for any value type with the arithmetic operators, element by
// element (no SIMD): packed elements are stored as raw bytes, P doubles each.
template <class T, long MC_ = 32, long NC_ = 32, long KC_ = 128> struct scalar_kernel {
    typedef T value_type;
    static const int planes = (int)((sizeof(T) + sizeof(double) - 1) / sizeof(double));
    static const long MC = MC_, NC = NC_, KC = KC_;
    static void pack(double *p, long stride, long e, const T &x) {
        double w[planes];
        std::memset(w, 0, sizeof(w));
        std::memcpy(w, &x, sizeof(T));
        for (int q = 0; q < planes; q++)
            p[q * stride + e] = w[q];
    }
    static T unpack(const double *p, long stride, long e) {
        double w[planes];
        for (int q = 0; q < planes; q++)
            w[q] = p[q * stride + e];
        T x;
        std::memcpy(&x, w, sizeof(T));
        return x;
    }
    static void block(bool axpy, long mc, long nc, long kc, double *acc, long as, const double *a, long ps, const double *b, long bs) {
        std::vector<T> av((size_t)mc * kc), cv((size_t)mc);
        for (long e = 0; e < mc * kc; e++)
            av[e] = unpack(a, ps, e);
        for (long j = 0; j < nc; j++) {
            for (long i = 0; i < mc; i++)
                cv[i] = unpack(acc, as, i + j * mc);
            for (long l = 0; l < kc; l++) {
                const T bv = unpack(b, bs, l + j * kc);
                const T *al = &av[l * mc];
                if (axpy)
                    for (long i = 0; i < mc; i++)
                        cv[i] += bv * al[i];
                else
                    for (long i = 0; i < mc; i++)
                        cv[i] += al[i] * bv;
            }
            for (long i = 0; i < mc; i++)
                pack(acc, as, i + j * mc, cv[i]);
        }
    }
};

} // namespace mplapack_gemm

#endif
