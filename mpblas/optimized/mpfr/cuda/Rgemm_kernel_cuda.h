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
 * Fixed-precision (512/1024-bit) MPFR-compatible Rgemm kernels.
 *
 * Values are cu_fp::cu_freal<PB> from mpc_cuda (https://github.com/tkouya/mpc_cuda):
 * the significand layout is the MPFR one and + - * are bit-exact with MPFR
 * round-to-nearest-even.  The element routines below perform the same
 * operations, in the same order, as mpblas/optimized/mpfr/openmp/Rgemm_*_omp.cpp,
 * so the result equals the CPU result for operands of precision PB.
 *
 * This header must not include gmpfrxx_mkII or MPFR: it is compiled by nvcc.
 */

#ifndef _MPLAPACK_MPFR_RGEMM_KERNEL_CUDA_H_
#define _MPLAPACK_MPFR_RGEMM_KERNEL_CUDA_H_

#include <mpc_cuda/cu_freal.cuh>

namespace mplapack_mpfr_cuda {

template <int PB> using real_t = cu_fp::cu_freal<PB>;

// Operands are packed column-major: A is nrowa x ncola with lda = nrowa,
// B is nrowb x ncolb with ldb = nrowb, C is m x n with ldc = m.
struct gemm_shape {
    int nota;       // transa == "N"
    int notb;       // transb == "N"
    long m, n, k;
    long lda, ldb, ldc;
    int beta_zero;  // beta == 0
    int beta_one;   // beta == 1
};

// op(B)(l, j)
template <int PB> __host__ __device__ inline const real_t<PB> &gemm_opB(const gemm_shape &s, const real_t<PB> *B, long l, long j) {
    return s.notb ? B[l + j * s.ldb] : B[j + l * s.ldb];
}

// AB(l, j) = alpha * op(B)(l, j), used when op(A) = A (NN and NT).
// Mirrors "temp = alpha * B[...]" in Rgemm_NN_omp / Rgemm_NT_omp.
template <int PB> __host__ __device__ inline void gemm_alpha_opB(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> *B, real_t<PB> *AB, long l, long j) { AB[l + j * s.k] = alpha * gemm_opB<PB>(s, B, l, j); }

// C(i, j) for one output element.
//   op(A) = A  : C(i,j) = beta*C(i,j); for l: C(i,j) += AB(l,j) * A(i,l)
//   op(A) = A' : C(i,j) = beta*C(i,j); temp = 0; for l: temp += A(l,i) * op(B)(l,j);
//                C(i,j) += alpha * temp
template <int PB> __host__ __device__ inline void gemm_element(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, const real_t<PB> *AB, real_t<PB> *C, long i, long j) {
    real_t<PB> c;
    if (s.beta_zero) {
        c.set_zero();
    } else if (s.beta_one) {
        c = C[i + j * s.ldc];
    } else {
        c = beta * C[i + j * s.ldc];
    }
    if (s.nota) {
        for (long l = 0; l < s.k; l++) {
            c = c + AB[l + j * s.k] * A[i + l * s.lda];
        }
    } else {
        real_t<PB> temp;
        temp.set_zero();
        for (long l = 0; l < s.k; l++) {
            temp = temp + A[l + i * s.lda] * gemm_opB<PB>(s, B, l, j);
        }
        c = c + alpha * temp;
    }
    C[i + j * s.ldc] = c;
}

// Computes C := alpha*op(A)*op(B) + beta*C on packed operands.
// Returns 0 on success; any other value means nothing was computed and the
// caller must fall back to the CPU implementation.
template <int PB> int gemm(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, real_t<PB> *C);

// Number of usable devices (the host implementation reports 1).
int device_count();

} // namespace mplapack_mpfr_cuda

#endif
