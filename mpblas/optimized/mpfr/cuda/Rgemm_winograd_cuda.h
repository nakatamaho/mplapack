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
 * Winograd variant of Strassen's matrix multiplication for the
 * fixed-precision (cu_freal<PB>) MPFR Rgemm.
 *
 * The schedule is the one of mul_mpfmatrix_winograd_even() in BNCmatmul
 * (https://github.com/tkouya/bncmatmul, src/matmul_strassen_general_mpf.c):
 *
 *   S1 = A21 + A22   S2 = S1 - A11   S3 = A11 - A21   S4 = A12 - S2
 *   S5 = B12 - B11   S6 = B22 - S5   S7 = B22 - B12   S8 = S6 - B21
 *   M1 = S2 S6  M2 = A11 B11  M3 = A12 B21  M4 = S3 S7
 *   M5 = S1 S5  M6 = S4 B22   M7 = A22 S8
 *   T1 = M1 + M2   T2 = T1 + M4
 *   C11 = M2 + M3  C12 = (T1 + M5) + M6  C21 = T2 - M7  C22 = T2 + M5
 *
 * The recursion stops when m, n or k is at most `cutoff` and uses the
 * conventional product (temp = 0; temp += a*b in the order of k), as
 * mul_mpfmatrix_simple() does.  An odd dimension is padded with a zero row
 * or column at that level.
 *
 * The result is not the one of the conventional Rgemm (different operations
 * and rounding); its error is bounded normwise, not elementwise.
 *
 * The recursion is written once against a backend that owns the storage and
 * performs the elementwise operations and the leaf products, so the CUDA
 * implementation (Rgemm_device_cuda.cu) and the host one used by the tests
 * (winograd_host_backend below) run exactly the same operations.
 */

#ifndef _MPLAPACK_MPFR_RGEMM_WINOGRAD_CUDA_H_
#define _MPLAPACK_MPFR_RGEMM_WINOGRAD_CUDA_H_

#include <vector>
#include "Rgemm_kernel_cuda.h"

namespace mplapack_mpfr_cuda {

// ---- element operations shared by the backends ----

template <int PB> __host__ __device__ inline void wg_zero(real_t<PB> *X, long ldx, long i, long j) { X[i + j * ldx].set_zero(); }

template <int PB> __host__ __device__ inline void wg_copy(real_t<PB> *D, long ldd, const real_t<PB> *S, long lds, long i, long j) { D[i + j * ldd] = S[i + j * lds]; }

// Z = X + Y, or X - Y when sub != 0
template <int PB> __host__ __device__ inline void wg_addsub(real_t<PB> *Z, long ldz, const real_t<PB> *X, long ldx, const real_t<PB> *Y, long ldy, int sub, long i, long j) {
    Z[i + j * ldz] = sub ? X[i + j * ldx] - Y[i + j * ldy] : X[i + j * ldx] + Y[i + j * ldy];
}

// C(i,j) = sum_l A(i,l) * B(l,j), accumulated in the order of l
template <int PB> __host__ __device__ inline void wg_leaf(long k, real_t<PB> *C, long ldc, const real_t<PB> *A, long lda, const real_t<PB> *B, long ldb, long i, long j) {
    real_t<PB> temp;
    temp.set_zero();
    for (long l = 0; l < k; l++)
        temp = temp + A[i + l * lda] * B[l + j * ldb];
    C[i + j * ldc] = temp;
}

// C(i,j) = alpha*P(i,j) + beta*C(i,j)  (beta == 0: C is not read)
template <int PB> __host__ __device__ inline void wg_combine(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *P, real_t<PB> *C, long i, long j) {
    real_t<PB> c = alpha * P[i + j * s.m];
    if (!s.beta_zero)
        c = c + (s.beta_one ? C[i + j * s.ldc] : beta * C[i + j * s.ldc]);
    C[i + j * s.ldc] = c;
}

// ---- the recursion ----
//
// Backend interface:
//   real_t<PB> *alloc(size_t count);   // NULL on failure (sets failed())
//   void release(real_t<PB> *);
//   bool failed() const;
//   void zero(long r, long c, real_t<PB> *X, long ldx);
//   void copy(long r, long c, real_t<PB> *D, long ldd, const real_t<PB> *S, long lds);
//   void addsub(long r, long c, real_t<PB> *Z, long ldz, const real_t<PB> *X, long ldx, const real_t<PB> *Y, long ldy, int sub);
//   void leaf(long m, long n, long k, real_t<PB> *C, long ldc, const real_t<PB> *A, long lda, const real_t<PB> *B, long ldb);

// C (m x n) = A (m x k) * B (k x n); all column-major with leading dimensions.
template <int PB, class Backend> void winograd(Backend &be, long m, long n, long k, const real_t<PB> *A, long lda, const real_t<PB> *B, long ldb, real_t<PB> *C, long ldc, long cutoff) {
    typedef real_t<PB> F;
    if (be.failed())
        return;
    if (m <= cutoff || n <= cutoff || k <= cutoff) {
        be.leaf(m, n, k, C, ldc, A, lda, B, ldb);
        return;
    }
    if ((m | n | k) & 1) {
        // Pad the odd dimensions with zeros and continue on even sizes.
        const long m2 = m + (m & 1), n2 = n + (n & 1), k2 = k + (k & 1);
        F *Ap = be.alloc((size_t)m2 * k2), *Bp = be.alloc((size_t)k2 * n2), *Cp = be.alloc((size_t)m2 * n2);
        if (!be.failed()) {
            be.zero(m2, k2, Ap, m2);
            be.zero(k2, n2, Bp, k2);
            be.copy(m, k, Ap, m2, A, lda);
            be.copy(k, n, Bp, k2, B, ldb);
            winograd<PB>(be, m2, n2, k2, Ap, m2, Bp, k2, Cp, m2, cutoff);
            be.copy(m, n, C, ldc, Cp, m2);
        }
        be.release(Ap);
        be.release(Bp);
        be.release(Cp);
        return;
    }
    const long mh = m / 2, nh = n / 2, kh = k / 2;
    const F *A11 = A, *A21 = A + mh, *A12 = A + kh * lda, *A22 = A + mh + kh * lda;
    const F *B11 = B, *B21 = B + kh, *B12 = B + nh * ldb, *B22 = B + kh + nh * ldb;
    F *C11 = C, *C21 = C + mh, *C12 = C + nh * ldc, *C22 = C + mh + nh * ldc;

    F *S[8], *M[7], *T[2];
    for (int t = 0; t < 4; t++)
        S[t] = be.alloc((size_t)mh * kh);
    for (int t = 4; t < 8; t++)
        S[t] = be.alloc((size_t)kh * nh);
    for (int t = 0; t < 7; t++)
        M[t] = be.alloc((size_t)mh * nh);
    for (int t = 0; t < 2; t++)
        T[t] = be.alloc((size_t)mh * nh);

    if (!be.failed()) {
        be.addsub(mh, kh, S[0], mh, A21, lda, A22, lda, 0);  // S1 = A21 + A22
        be.addsub(mh, kh, S[1], mh, S[0], mh, A11, lda, 1);  // S2 = S1 - A11
        be.addsub(mh, kh, S[2], mh, A11, lda, A21, lda, 1);  // S3 = A11 - A21
        be.addsub(mh, kh, S[3], mh, A12, lda, S[1], mh, 1);  // S4 = A12 - S2
        be.addsub(kh, nh, S[4], kh, B12, ldb, B11, ldb, 1);  // S5 = B12 - B11
        be.addsub(kh, nh, S[5], kh, B22, ldb, S[4], kh, 1);  // S6 = B22 - S5
        be.addsub(kh, nh, S[6], kh, B22, ldb, B12, ldb, 1);  // S7 = B22 - B12
        be.addsub(kh, nh, S[7], kh, S[5], kh, B21, ldb, 1);  // S8 = S6 - B21

        winograd<PB>(be, mh, nh, kh, S[1], mh, S[5], kh, M[0], mh, cutoff); // M1 = S2 S6
        winograd<PB>(be, mh, nh, kh, A11, lda, B11, ldb, M[1], mh, cutoff); // M2 = A11 B11
        winograd<PB>(be, mh, nh, kh, A12, lda, B21, ldb, M[2], mh, cutoff); // M3 = A12 B21
        winograd<PB>(be, mh, nh, kh, S[2], mh, S[6], kh, M[3], mh, cutoff); // M4 = S3 S7
        winograd<PB>(be, mh, nh, kh, S[0], mh, S[4], kh, M[4], mh, cutoff); // M5 = S1 S5
        winograd<PB>(be, mh, nh, kh, S[3], mh, B22, ldb, M[5], mh, cutoff); // M6 = S4 B22
        winograd<PB>(be, mh, nh, kh, A22, lda, S[7], kh, M[6], mh, cutoff); // M7 = A22 S8

        be.addsub(mh, nh, T[0], mh, M[0], mh, M[1], mh, 0);  // T1 = M1 + M2
        be.addsub(mh, nh, T[1], mh, T[0], mh, M[3], mh, 0);  // T2 = T1 + M4
        be.addsub(mh, nh, C11, ldc, M[1], mh, M[2], mh, 0);  // C11 = M2 + M3
        be.addsub(mh, nh, C12, ldc, T[0], mh, M[4], mh, 0);  // C12 = T1 + M5
        be.addsub(mh, nh, C12, ldc, C12, ldc, M[5], mh, 0);  //     + M6
        be.addsub(mh, nh, C21, ldc, T[1], mh, M[6], mh, 1);  // C21 = T2 - M7
        be.addsub(mh, nh, C22, ldc, T[1], mh, M[4], mh, 0);  // C22 = T2 + M5
    }
    for (int t = 0; t < 8; t++)
        be.release(S[t]);
    for (int t = 0; t < 7; t++)
        be.release(M[t]);
    for (int t = 0; t < 2; t++)
        be.release(T[t]);
}

// ---- host backend (reference for the tests) ----

template <int PB> class winograd_host_backend {
  public:
    typedef real_t<PB> F;
    winograd_host_backend() {}
    ~winograd_host_backend() {}
    F *alloc(size_t count) { return new F[count ? count : 1]; }
    void release(F *p) { delete[] p; }
    bool failed() const { return false; }
    void zero(long r, long c, F *X, long ldx) {
        for (long j = 0; j < c; j++)
            for (long i = 0; i < r; i++)
                wg_zero<PB>(X, ldx, i, j);
    }
    void copy(long r, long c, F *D, long ldd, const F *S, long lds) {
        for (long j = 0; j < c; j++)
            for (long i = 0; i < r; i++)
                wg_copy<PB>(D, ldd, S, lds, i, j);
    }
    void addsub(long r, long c, F *Z, long ldz, const F *X, long ldx, const F *Y, long ldy, int sub) {
        for (long j = 0; j < c; j++)
            for (long i = 0; i < r; i++)
                wg_addsub<PB>(Z, ldz, X, ldx, Y, ldy, sub, i, j);
    }
    void leaf(long m, long n, long k, F *C, long ldc, const F *A, long lda, const F *B, long ldb) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (long j = 0; j < n; j++)
            for (long i = 0; i < m; i++)
                wg_leaf<PB>(k, C, ldc, A, lda, B, ldb, i, j);
    }
};

// C := alpha*A*B + beta*C on the host; A is m x k, B is k x n (s.m, s.n, s.k),
// C has leading dimension s.ldc.
template <int PB> void gemm_winograd_host(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, real_t<PB> *C, long cutoff) {
    winograd_host_backend<PB> be;
    std::vector<real_t<PB>> P((size_t)s.m * s.n + 1);
    winograd<PB>(be, s.m, s.n, s.k, A, s.m, B, s.k, P.data(), s.m, cutoff);
    for (long j = 0; j < s.n; j++)
        for (long i = 0; i < s.m; i++)
            wg_combine<PB>(s, alpha, beta, P.data(), C, i, j);
}

// Same on the GPU when available (Rgemm_device_cuda.cu); the host build
// (Rgemm_host_cuda.cpp) runs gemm_winograd_host.  A and B are op(A) and
// op(B) already (m x k and k x n).  Returns 0 on success.
template <int PB> int gemm_winograd(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, real_t<PB> *C, long cutoff);

} // namespace mplapack_mpfr_cuda

#endif
