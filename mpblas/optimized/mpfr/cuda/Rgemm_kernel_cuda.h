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
 * Fixed-precision (256/512/768/1024/2048-bit) MPFR-compatible Rgemm kernels.
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

// Precisions (bits) served by the fixed-precision kernels; X(PB) is expanded
// for each of them (explicit instantiations, dispatch).
#define MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(X) X(256) X(512) X(768) X(1024) X(2048)

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
// Engine form also serves device-resident LAPACK. packed_ab=false computes
// alpha*op(B) locally, retaining exactly the same rounding and operation order.
// Register slots 0..3 are private to this element.
template<class E> __host__ __device__ inline void gemm_engine_element(
    E &e, const gemm_shape &s, const typename E::M &al, const typename E::M &be,
    const typename E::M &A, const typename E::M &B, const typename E::M &AB,
    typename E::M &C, long i, long j, bool packed_ab) {
    auto c=e.reg(0), temp=e.reg(1), p=e.reg(2), ab=e.reg(3);
    if(s.beta_zero) e.zero(c);
    else if(s.beta_one) e.set(c,e.at(C,i,j));
    else e.mul(c,e.at(be,0,0),e.at(C,i,j));
    if(s.nota) {
        for(long l=0;l<s.k;++l) {
            if(packed_ab) e.set(ab,e.at(AB,l,j));
            else e.mul(ab,e.at(al,0,0),s.notb?e.at(B,l,j):e.at(B,j,l));
            e.mul(p,ab,e.at(A,i,l));
            e.add(c,c,p);
        }
    } else {
        e.zero(temp);
        for(long l=0;l<s.k;++l) {
            e.mul(p,e.at(A,l,i),s.notb?e.at(B,l,j):e.at(B,j,l));
            e.add(temp,temp,p);
        }
        e.mul(p,e.at(al,0,0),temp);
        e.add(c,c,p);
    }
    e.put(C,i,j,c);
}
template<int PB> struct gemm_fixed_engine {
    using V=real_t<PB>;
    struct M { const V *p; long ld; };
    __host__ __device__ V reg(int) const { V v; v.set_zero(); return v; }
    __host__ __device__ V at(const M &a,long i,long j) const { return a.p[i+j*a.ld]; }
    __host__ __device__ void put(M &a,long i,long j,const V &v) const { const_cast<V*>(a.p)[i+j*a.ld]=v; }
    __host__ __device__ void set(V &r,const V &a) const { r=a; }
    __host__ __device__ void zero(V &r) const { r.set_zero(); }
    __host__ __device__ void mul(V &r,const V &a,const V &b) const { r=a*b; }
    __host__ __device__ void add(V &r,const V &a,const V &b) const { r=a+b; }
};
template <int PB> __host__ __device__ inline void gemm_element(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, const real_t<PB> *AB, real_t<PB> *C, long i, long j) {
    gemm_fixed_engine<PB> e;
    using M=typename gemm_fixed_engine<PB>::M;
    M c{C,s.ldc};
    gemm_engine_element(e,s,M{&alpha,1},M{&beta,1},M{A,s.lda},M{B,s.ldb},M{AB,s.k},c,i,j,true);
}

// Computes C := alpha*op(A)*op(B) + beta*C on packed operands.
// Returns 0 on success; any other value means nothing was computed and the
// caller must fall back to the CPU implementation.
template <int PB> int gemm(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, real_t<PB> *C);

// Number of usable devices (the host implementation reports 1).
int device_count();

} // namespace mplapack_mpfr_cuda

#endif
