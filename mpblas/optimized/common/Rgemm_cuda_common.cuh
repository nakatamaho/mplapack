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
 * Rgemm on an NVIDIA GPU for double-double and quad-double, shared by
 * libmplapack_dd_opt_cuda and libmplapack_qd_opt_cuda.
 *
 * Algorithms (gemm_cuda's `algo`):
 *   GPU_NAIVE     one thread per element of C, operands read from global memory
 *   GPU_TILED     16 x 16 threads per block, 16 x 16 tiles of op(A) and of
 *                 alpha*op(B) / op(B) staged in shared memory
 *   GPU_WINOGRAD  Winograd/Strassen recursion (Rgemm_winograd_common.h), all
 *                 temporaries on the GPU, tiled leaves
 *   GPU_OZAKI     Ozaki scheme (Rgemm_ozaki_common.h): the slices are made on
 *                 the host, the exact slice products by cuBLAS DGEMM, the sums
 *                 in the target precision on the GPU
 *
 * NAIVE and TILED perform, for every element of C, the operations of
 * openmp/Rgemm_*_omp.cpp in the same order, with device functions that
 * reproduce libqd's (Rgemm_qd_ops_common.h): their result is bit-identical
 * to the CPU Rgemm.  OZAKI is bit-identical to the CPU Ozaki path
 * (gemm_ozaki_cpu) with the same slices.  WINOGRAD differs from the
 * conventional result in rounding (normwise error bound).
 *
 * When compiled without nvcc and with MPLAPACK_GEMM_CUDA_HOST_EMULATION
 * defined, the same element code runs in host loops and the DGEMM of the
 * Ozaki scheme is ozaki_dgemm; this is what the GPU-less tests run.
 *
 * Code including this file must be compiled without floating-point
 * contraction (nvcc -fmad=false, g++ -ffp-contract=off).
 */

#ifndef _MPLAPACK_RGEMM_CUDA_COMMON_CUH_
#define _MPLAPACK_RGEMM_CUDA_COMMON_CUH_

#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "Rgemm_cuda_ops_common.h"
#include "Rgemm_blocked_common.h"
#include "Rgemm_winograd_common.h"
#include "Rgemm_ozaki_common.h"

#if defined(__CUDACC__)
#include <cuda_runtime.h>
#include <cublas_v2.h>
#define MPLAPACK_GEMM_GLOBAL __global__
#elif defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
#define MPLAPACK_GEMM_GLOBAL
#else
#error "Rgemm_cuda_common.cuh: compile with nvcc, or define MPLAPACK_GEMM_CUDA_HOST_EMULATION"
#endif

namespace mplapack_gemm {

// ---- element code ----

template <class Ops> MPLAPACK_GEMM_HD typename Ops::F gpu_scaled_c(const gpu_scalars<typename Ops::F> &sc, const typename Ops::F &c) {
    if (sc.beta_zero)
        return Ops::from_d(0.0);
    if (!sc.beta_one)
        return Ops::mul(sc.beta, c);
    return c;
}

// Bs(l,j) = alpha*op(B)(l,j)  (k x n, leading dimension k), for op(A) == A
template <class Ops> MPLAPACK_GEMM_HD void gpu_alpha_opB(const shape &s, const typename Ops::F &alpha, const typename Ops::F *B, typename Ops::F *Bs, long l, long j) {
    Bs[l + j * s.k] = Ops::mul(alpha, s.notb ? B[l + j * s.ldb] : B[j + l * s.ldb]);
}

// The conventional element (i,j).  op(A) == A: Bs is alpha*op(B) (k x n).
// op(A) == A': Bs is B itself (leading dimension s.ldb, op given by s.notb).
template <class Ops> MPLAPACK_GEMM_HD void gpu_element(const shape &s, const gpu_scalars<typename Ops::F> &sc, const typename Ops::F *A, const typename Ops::F *Bs, typename Ops::F *C, long i, long j) {
    typedef typename Ops::F F;
    F &c = C[i + j * s.ldc];
    if (s.nota) {
        F acc = gpu_scaled_c<Ops>(sc, c);
        for (long l = 0; l < s.k; l++)
            acc = Ops::add(acc, Ops::mul(Bs[l + j * s.k], A[i + l * s.lda]));
        c = acc;
    } else {
        F t = Ops::from_d(0.0);
        for (long l = 0; l < s.k; l++)
            t = Ops::add(t, Ops::mul(A[l + i * s.lda], s.notb ? Bs[l + j * s.ldb] : Bs[j + l * s.ldb]));
        c = Ops::add(gpu_scaled_c<Ops>(sc, c), Ops::mul(sc.alpha, t));
    }
}

// Winograd leaf: C(i,j) = sum_l A(i,l) * B(l,j) (t = 0; t = t + a*b)
template <class Ops> MPLAPACK_GEMM_HD void gpu_leaf_element(long k, typename Ops::F *C, long ldc, const typename Ops::F *A, long lda, const typename Ops::F *B, long ldb, long i, long j) {
    typename Ops::F t = Ops::from_d(0.0);
    for (long l = 0; l < k; l++)
        t = Ops::add(t, Ops::mul(A[i + l * lda], B[l + j * ldb]));
    C[i + j * ldc] = t;
}

// C(i,j) := alpha*P(i,j) + beta*C(i,j) (as gemm_winograd_cpu)
template <class Ops> MPLAPACK_GEMM_HD void gpu_combine(const shape &s, const gpu_scalars<typename Ops::F> &sc, const typename Ops::F *P, typename Ops::F *C, long i, long j) {
    typedef typename Ops::F F;
    F &c = C[i + j * s.ldc];
    F r = Ops::mul(sc.alpha, P[i + j * s.m]);
    if (!sc.beta_zero)
        r = Ops::add(r, sc.beta_one ? c : Ops::mul(sc.beta, c));
    c = r;
}

// Ozaki: acc(i,j) += D(i,j) 2^-sh, then C(i,j) := alpha 2^(e_i+f_j) acc + beta*C (as ozaki_finish)
template <class Ops> MPLAPACK_GEMM_HD void gpu_ozaki_accumulate(long m, typename Ops::F *acc, const double *D, int sh, long i, long j) { acc[i + j * m] = Ops::add_d(acc[i + j * m], ldexp(D[i + j * m], -sh)); }

template <class Ops> MPLAPACK_GEMM_HD void gpu_ozaki_finish(const shape &s, const gpu_scalars<typename Ops::F> &sc, const typename Ops::F *acc, const int *ea, const int *fb, typename Ops::F *C, long i, long j) {
    typedef typename Ops::F F;
    F &c = C[i + j * s.ldc];
    F p = Ops::mul(Ops::mul(acc[i + j * s.m], Ops::from_d(ldexp(1.0, ea[i]))), Ops::from_d(ldexp(1.0, fb[j])));
    F r = Ops::mul(sc.alpha, p);
    if (!sc.beta_zero)
        r = Ops::add(r, sc.beta_one ? c : Ops::mul(sc.beta, c));
    c = r;
}

template <class Ops> MPLAPACK_GEMM_HD void gpu_addsub(typename Ops::F *Z, long ldz, const typename Ops::F *X, long ldx, const typename Ops::F *Y, long ldy, int sub, long i, long j) {
    Z[i + j * ldz] = sub ? Ops::sub(X[i + j * ldx], Y[i + j * ldy]) : Ops::add(X[i + j * ldx], Y[i + j * ldy]);
}

// The digit slices as doubles: dst[e] = src[e]
MPLAPACK_GEMM_HD void gpu_i32_to_double(const int32_t *src, double *dst, size_t e) { dst[e] = (double)src[e]; }

#if defined(__CUDACC__)
// ======================================================================
// Device implementation
// ======================================================================

static const int GPU_TILE = 16;

inline dim3 gpu_grid2(long r, long c) { return dim3((unsigned)((r + GPU_TILE - 1) / GPU_TILE), (unsigned)((c + GPU_TILE - 1) / GPU_TILE)); }
inline unsigned gpu_grid1(size_t count) { return (unsigned)((count + 255) / 256); }

template <class Ops> __global__ void k_alpha_opB(shape s, typename Ops::F alpha, const typename Ops::F *B, typename Ops::F *Bs) {
    const long l = blockIdx.x * (long)blockDim.x + threadIdx.x, j = blockIdx.y * (long)blockDim.y + threadIdx.y;
    if (l < s.k && j < s.n)
        gpu_alpha_opB<Ops>(s, alpha, B, Bs, l, j);
}

template <class Ops> __global__ void k_naive(shape s, gpu_scalars<typename Ops::F> sc, const typename Ops::F *A, const typename Ops::F *Bs, typename Ops::F *C) {
    const long i = blockIdx.x * (long)blockDim.x + threadIdx.x, j = blockIdx.y * (long)blockDim.y + threadIdx.y;
    if (i < s.m && j < s.n)
        gpu_element<Ops>(s, sc, A, Bs, C, i, j);
}

// Tiled conventional kernel.  mode 0: C = conventional Rgemm element
// (requires the shape's nota/notb semantics as gpu_element); mode 1: plain
// product C = A*B for the Winograd leaves (A m x k, B k x n, column major).
template <class Ops> __global__ void k_tiled(shape s, gpu_scalars<typename Ops::F> sc, const typename Ops::F *A, const typename Ops::F *Bs, typename Ops::F *C, int leaf) {
    typedef typename Ops::F F;
    __shared__ F As[GPU_TILE][GPU_TILE + 1]; // As[l][i]
    __shared__ F Bt[GPU_TILE][GPU_TILE + 1]; // Bt[j][l]
    const int ti = threadIdx.x, tj = threadIdx.y;
    const long i = blockIdx.x * (long)GPU_TILE + ti, j = blockIdx.y * (long)GPU_TILE + tj;
    const bool inside = i < s.m && j < s.n;
    const bool axpy = !leaf && s.nota;
    F acc;
    if (leaf)
        acc = Ops::from_d(0.0);
    else if (axpy)
        acc = inside ? gpu_scaled_c<Ops>(sc, C[i + j * s.ldc]) : Ops::from_d(0.0);
    else
        acc = Ops::from_d(0.0);
    for (long p0 = 0; p0 < s.k; p0 += GPU_TILE) {
        // load op(A)(i0 + ti, p0 + tj) and the B operand (p0 + ti, j0 + tj)
        {
            const long ai = blockIdx.x * (long)GPU_TILE + ti, al = p0 + tj;
            if (ai < s.m && al < s.k)
                As[tj][ti] = (leaf || s.nota) ? A[ai + al * s.lda] : A[al + ai * s.lda];
            const long bl = p0 + ti, bj = blockIdx.y * (long)GPU_TILE + tj;
            if (bl < s.k && bj < s.n) {
                if (leaf)
                    Bt[tj][ti] = Bs[bl + bj * s.ldb];
                else if (s.nota)
                    Bt[tj][ti] = Bs[bl + bj * s.k];
                else
                    Bt[tj][ti] = s.notb ? Bs[bl + bj * s.ldb] : Bs[bj + bl * s.ldb];
            }
        }
        __syncthreads();
        if (inside) {
            const long kc = (s.k - p0 < GPU_TILE) ? s.k - p0 : GPU_TILE;
            if (axpy)
                for (long l = 0; l < kc; l++)
                    acc = Ops::add(acc, Ops::mul(Bt[tj][l], As[l][ti]));
            else
                for (long l = 0; l < kc; l++)
                    acc = Ops::add(acc, Ops::mul(As[l][ti], Bt[tj][l]));
        }
        __syncthreads();
    }
    if (!inside)
        return;
    if (leaf || axpy)
        C[i + j * s.ldc] = acc;
    else
        C[i + j * s.ldc] = Ops::add(gpu_scaled_c<Ops>(sc, C[i + j * s.ldc]), Ops::mul(sc.alpha, acc));
}

template <class Ops> __global__ void k_zero(long r, long c, typename Ops::F *X, long ldx) {
    const long i = blockIdx.x * (long)blockDim.x + threadIdx.x, j = blockIdx.y * (long)blockDim.y + threadIdx.y;
    if (i < r && j < c)
        X[i + j * ldx] = Ops::from_d(0.0);
}
template <class Ops> __global__ void k_copy(long r, long c, typename Ops::F *D, long ldd, const typename Ops::F *S, long lds) {
    const long i = blockIdx.x * (long)blockDim.x + threadIdx.x, j = blockIdx.y * (long)blockDim.y + threadIdx.y;
    if (i < r && j < c)
        D[i + j * ldd] = S[i + j * lds];
}
template <class Ops> __global__ void k_addsub(long r, long c, typename Ops::F *Z, long ldz, const typename Ops::F *X, long ldx, const typename Ops::F *Y, long ldy, int sub) {
    const long i = blockIdx.x * (long)blockDim.x + threadIdx.x, j = blockIdx.y * (long)blockDim.y + threadIdx.y;
    if (i < r && j < c)
        gpu_addsub<Ops>(Z, ldz, X, ldx, Y, ldy, sub, i, j);
}
template <class Ops> __global__ void k_combine(shape s, gpu_scalars<typename Ops::F> sc, const typename Ops::F *P, typename Ops::F *C) {
    const long i = blockIdx.x * (long)blockDim.x + threadIdx.x, j = blockIdx.y * (long)blockDim.y + threadIdx.y;
    if (i < s.m && j < s.n)
        gpu_combine<Ops>(s, sc, P, C, i, j);
}
template <class Ops> __global__ void k_ozaki_accumulate(long m, long n, typename Ops::F *acc, const double *D, int sh) {
    const long i = blockIdx.x * (long)blockDim.x + threadIdx.x, j = blockIdx.y * (long)blockDim.y + threadIdx.y;
    if (i < m && j < n)
        gpu_ozaki_accumulate<Ops>(m, acc, D, sh, i, j);
}
template <class Ops> __global__ void k_ozaki_finish(shape s, gpu_scalars<typename Ops::F> sc, const typename Ops::F *acc, const int *ea, const int *fb, typename Ops::F *C) {
    const long i = blockIdx.x * (long)blockDim.x + threadIdx.x, j = blockIdx.y * (long)blockDim.y + threadIdx.y;
    if (i < s.m && j < s.n)
        gpu_ozaki_finish<Ops>(s, sc, acc, ea, fb, C, i, j);
}
__global__ inline void k_i32_to_double(const int32_t *src, double *dst, size_t count) {
    const size_t e = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
    if (e < count)
        gpu_i32_to_double(src, dst, e);
}

// Device memory owned by one call; every allocation and launch error is
// remembered and makes the call fail (the caller then runs on the CPU).
class gpu_session {
  public:
    gpu_session() : err_(cudaSuccess) {}
    ~gpu_session() {
        for (size_t t = 0; t < ptrs_.size(); t++)
            cudaFree(ptrs_[t]);
    }
    template <class X> X *alloc(size_t count) {
        void *p = NULL;
        if (err_ != cudaSuccess)
            return NULL;
        err_ = cudaMalloc(&p, (count ? count : 1) * sizeof(X));
        if (err_ != cudaSuccess)
            return NULL;
        ptrs_.push_back(p);
        return (X *)p;
    }
    void release(void *p) {
        for (size_t t = 0; t < ptrs_.size(); t++)
            if (ptrs_[t] == p) {
                cudaFree(p);
                ptrs_.erase(ptrs_.begin() + t);
                return;
            }
    }
    void check(cudaError_t e) {
        if (err_ == cudaSuccess && e != cudaSuccess)
            err_ = e;
    }
    void check_launch() { check(cudaGetLastError()); }
    bool ok() const { return err_ == cudaSuccess; }
    cudaError_t error() const { return err_; }

  private:
    cudaError_t err_;
    std::vector<void *> ptrs_;
};

template <class Ops> class winograd_gpu_backend {
  public:
    typedef typename Ops::F F;
    explicit winograd_gpu_backend(gpu_session &g) : g_(g) {}
    F *alloc(size_t count) { return g_.alloc<F>(count); }
    void release(F *p) { g_.release(p); }
    bool failed() const { return !g_.ok(); }
    void zero(long r, long c, F *X, long ldx) {
        if (r > 0 && c > 0)
            k_zero<Ops><<<gpu_grid2(r, c), dim3(GPU_TILE, GPU_TILE)>>>(r, c, X, ldx);
        g_.check_launch();
    }
    void copy(long r, long c, F *D, long ldd, const F *S, long lds) {
        if (r > 0 && c > 0)
            k_copy<Ops><<<gpu_grid2(r, c), dim3(GPU_TILE, GPU_TILE)>>>(r, c, D, ldd, S, lds);
        g_.check_launch();
    }
    void addsub(long r, long c, F *Z, long ldz, const F *X, long ldx, const F *Y, long ldy, int sub) {
        if (r > 0 && c > 0)
            k_addsub<Ops><<<gpu_grid2(r, c), dim3(GPU_TILE, GPU_TILE)>>>(r, c, Z, ldz, X, ldx, Y, ldy, sub);
        g_.check_launch();
    }
    void leaf(long m, long n, long k, F *C, long ldc, const F *A, long lda, const F *B, long ldb) {
        if (m <= 0 || n <= 0)
            return;
        shape s = {true, true, m, n, k, lda, ldb, ldc};
        gpu_scalars<F> sc;
        sc.beta_zero = true;
        sc.beta_one = false;
        k_tiled<Ops><<<gpu_grid2(m, n), dim3(GPU_TILE, GPU_TILE)>>>(s, sc, A, B, C, 1);
        g_.check_launch();
    }

  private:
    gpu_session &g_;
};

// Number of elements of a column-major r x c matrix with leading dimension ld
inline size_t gpu_span(long r, long c, long ld) { return (r <= 0 || c <= 0) ? 0 : (size_t)ld * (c - 1) + r; }

template <class Ops> int gemm_cuda(const shape &s, int algo, long cutoff, const ozaki_plan *plan, const gpu_scalars<typename Ops::F> &sc, const typename Ops::F *A, const typename Ops::F *B, typename Ops::F *C) {
    typedef typename Ops::F F;
    gpu_session g;
    const long m = s.m, n = s.n, k = s.k;
    const size_t nC = gpu_span(m, n, s.ldc);
    F *dC = g.alloc<F>(nC);
    if (g.ok())
        g.check(cudaMemcpy(dC, C, nC * sizeof(F), cudaMemcpyHostToDevice));

    if (algo == GPU_OZAKI) {
        const int S = plan->S;
        int32_t *di = g.alloc<int32_t>((size_t)S * (m * k > k * n ? m * k : k * n));
        double *dA = g.alloc<double>((size_t)S * m * k), *dB = g.alloc<double>((size_t)S * k * n), *dD = g.alloc<double>((size_t)m * n);
        F *dacc = g.alloc<F>((size_t)m * n);
        int *dea = g.alloc<int>(m), *dfb = g.alloc<int>(n);
        if (g.ok()) {
            g.check(cudaMemcpy(di, plan->ad.data(), (size_t)S * m * k * sizeof(int32_t), cudaMemcpyHostToDevice));
            if (S * m * k > 0)
                k_i32_to_double<<<gpu_grid1((size_t)S * m * k), 256>>>(di, dA, (size_t)S * m * k);
            g.check(cudaMemcpy(di, plan->bd.data(), (size_t)S * k * n * sizeof(int32_t), cudaMemcpyHostToDevice));
            if (S * k * n > 0)
                k_i32_to_double<<<gpu_grid1((size_t)S * k * n), 256>>>(di, dB, (size_t)S * k * n);
            g.check(cudaMemcpy(dea, plan->ea.data(), m * sizeof(int), cudaMemcpyHostToDevice));
            g.check(cudaMemcpy(dfb, plan->fb.data(), n * sizeof(int), cudaMemcpyHostToDevice));
            k_zero<Ops><<<gpu_grid2(m, n), dim3(GPU_TILE, GPU_TILE)>>>(m, n, dacc, m);
            g.check_launch();
        }
        cublasHandle_t h = NULL;
        if (g.ok()) {
            const cublasStatus_t st = cublasCreate(&h);
            if (st != CUBLAS_STATUS_SUCCESS)
                return (int)st;
            cublasSetMathMode(h, CUBLAS_DEFAULT_MATH);
        }
        int blas_err = 0;
        for (int gl = S - 1; gl >= 0 && g.ok() && !blas_err; gl--)
            for (int a = 0; a <= gl && g.ok() && !blas_err; a++) {
                const int b = gl - a;
                const double one = 1.0, zero = 0.0;
                const cublasStatus_t st = cublasDgemm(h, CUBLAS_OP_N, CUBLAS_OP_N, (int)m, (int)n, (int)k, &one, dA + (size_t)a * m * k, (int)m, dB + (size_t)b * k * n, (int)k, &zero, dD, (int)m);
                if (st != CUBLAS_STATUS_SUCCESS) {
                    blas_err = (int)st;
                    break;
                }
                k_ozaki_accumulate<Ops><<<gpu_grid2(m, n), dim3(GPU_TILE, GPU_TILE)>>>(m, n, dacc, dD, (gl + 2) * plan->beta);
                g.check_launch();
            }
        if (h)
            cublasDestroy(h);
        if (blas_err)
            return blas_err;
        if (g.ok()) {
            k_ozaki_finish<Ops><<<gpu_grid2(m, n), dim3(GPU_TILE, GPU_TILE)>>>(s, sc, dacc, dea, dfb, dC);
            g.check_launch();
        }
    } else if (algo == GPU_WINOGRAD) {
        // op(A) (m x k) and op(B) (k x n) as plain column-major matrices
        const size_t nA = gpu_span(s.nota ? m : k, s.nota ? k : m, s.lda), nB = gpu_span(s.notb ? k : n, s.notb ? n : k, s.ldb);
        F *dA = g.alloc<F>(nA), *dB = g.alloc<F>(nB), *oA = g.alloc<F>((size_t)m * k), *oB = g.alloc<F>((size_t)k * n), *dP = g.alloc<F>((size_t)m * n);
        if (g.ok()) {
            g.check(cudaMemcpy(dA, A, nA * sizeof(F), cudaMemcpyHostToDevice));
            g.check(cudaMemcpy(dB, B, nB * sizeof(F), cudaMemcpyHostToDevice));
        }
        if (g.ok()) {
            // transposition through the leaf kernel would change the result;
            // copy element by element instead
            std::vector<F> hA((size_t)m * k), hB((size_t)k * n);
            for (long l = 0; l < k; l++)
                for (long i = 0; i < m; i++)
                    hA[i + l * m] = s.nota ? A[i + l * s.lda] : A[l + i * s.lda];
            for (long j = 0; j < n; j++)
                for (long l = 0; l < k; l++)
                    hB[l + j * k] = s.notb ? B[l + j * s.ldb] : B[j + l * s.ldb];
            g.check(cudaMemcpy(oA, hA.data(), hA.size() * sizeof(F), cudaMemcpyHostToDevice));
            g.check(cudaMemcpy(oB, hB.data(), hB.size() * sizeof(F), cudaMemcpyHostToDevice));
            g.release(dA);
            g.release(dB);
        }
        if (g.ok()) {
            winograd_gpu_backend<Ops> be(g);
            winograd<F>(be, m, n, k, oA, m, oB, k, dP, m, cutoff);
        }
        if (g.ok()) {
            k_combine<Ops><<<gpu_grid2(m, n), dim3(GPU_TILE, GPU_TILE)>>>(s, sc, dP, dC);
            g.check_launch();
        }
    } else {
        const size_t nA = gpu_span(s.nota ? m : k, s.nota ? k : m, s.lda), nB = gpu_span(s.notb ? k : n, s.notb ? n : k, s.ldb);
        F *dA = g.alloc<F>(nA), *dB = g.alloc<F>(nB);
        F *dBs = s.nota ? g.alloc<F>((size_t)k * n) : dB;
        if (g.ok()) {
            g.check(cudaMemcpy(dA, A, nA * sizeof(F), cudaMemcpyHostToDevice));
            g.check(cudaMemcpy(dB, B, nB * sizeof(F), cudaMemcpyHostToDevice));
        }
        if (g.ok() && s.nota && k > 0) {
            k_alpha_opB<Ops><<<gpu_grid2(k, n), dim3(GPU_TILE, GPU_TILE)>>>(s, sc.alpha, dB, dBs);
            g.check_launch();
        }
        if (g.ok()) {
            if (algo == GPU_NAIVE)
                k_naive<Ops><<<gpu_grid2(m, n), dim3(GPU_TILE, GPU_TILE)>>>(s, sc, dA, dBs, dC);
            else
                k_tiled<Ops><<<gpu_grid2(m, n), dim3(GPU_TILE, GPU_TILE)>>>(s, sc, dA, dBs, dC, 0);
            g.check_launch();
        }
    }
    if (g.ok())
        g.check(cudaDeviceSynchronize());
    if (g.ok())
        g.check(cudaMemcpy(C, dC, nC * sizeof(F), cudaMemcpyDeviceToHost));
    return g.ok() ? 0 : (int)g.error();
}

template <class Ops> int gpu_device_count() {
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess)
        return 0;
    return count;
}

#else
// ======================================================================
// Host emulation (tests without a GPU): the same element code in loops
// ======================================================================

template <class Ops> class winograd_emul_backend {
  public:
    typedef typename Ops::F F;
    F *alloc(size_t count) { return new F[count ? count : 1]; }
    void release(F *p) { delete[] p; }
    bool failed() const { return false; }
    void zero(long r, long c, F *X, long ldx) {
        for (long j = 0; j < c; j++)
            for (long i = 0; i < r; i++)
                X[i + j * ldx] = Ops::from_d(0.0);
    }
    void copy(long r, long c, F *D, long ldd, const F *S, long lds) {
        for (long j = 0; j < c; j++)
            for (long i = 0; i < r; i++)
                D[i + j * ldd] = S[i + j * lds];
    }
    void addsub(long r, long c, F *Z, long ldz, const F *X, long ldx, const F *Y, long ldy, int sub) {
        for (long j = 0; j < c; j++)
            for (long i = 0; i < r; i++)
                gpu_addsub<Ops>(Z, ldz, X, ldx, Y, ldy, sub, i, j);
    }
    void leaf(long m, long n, long k, F *C, long ldc, const F *A, long lda, const F *B, long ldb) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (long j = 0; j < n; j++)
            for (long i = 0; i < m; i++)
                gpu_leaf_element<Ops>(k, C, ldc, A, lda, B, ldb, i, j);
    }
};

template <class Ops> int gemm_cuda(const shape &s, int algo, long cutoff, const ozaki_plan *plan, const gpu_scalars<typename Ops::F> &sc, const typename Ops::F *A, const typename Ops::F *B, typename Ops::F *C) {
    typedef typename Ops::F F;
    const long m = s.m, n = s.n, k = s.k;
    if (algo == GPU_OZAKI) {
        std::vector<F> acc((size_t)m * n, Ops::from_d(0.0));
        std::vector<double> dA((size_t)plan->S * m * k), dB((size_t)plan->S * k * n), D((size_t)m * n);
        for (size_t e = 0; e < dA.size(); e++)
            gpu_i32_to_double(plan->ad.data(), dA.data(), e);
        for (size_t e = 0; e < dB.size(); e++)
            gpu_i32_to_double(plan->bd.data(), dB.data(), e);
        for (int gl = plan->S - 1; gl >= 0; gl--)
            for (int a = 0; a <= gl; a++) {
                const int b = gl - a;
                // exact, so any order of operations gives this result
                for (long j = 0; j < n; j++)
                    for (long i = 0; i < m; i++) {
                        double t = 0.0;
                        for (long l = 0; l < k; l++)
                            t += dA[(size_t)a * m * k + i + l * m] * dB[(size_t)b * k * n + l + j * k];
                        D[i + j * m] = t;
                    }
                for (long j = 0; j < n; j++)
                    for (long i = 0; i < m; i++)
                        gpu_ozaki_accumulate<Ops>(m, acc.data(), D.data(), (gl + 2) * plan->beta, i, j);
            }
        for (long j = 0; j < n; j++)
            for (long i = 0; i < m; i++)
                gpu_ozaki_finish<Ops>(s, sc, acc.data(), plan->ea.data(), plan->fb.data(), C, i, j);
        return 0;
    }
    if (algo == GPU_WINOGRAD) {
        std::vector<F> oA((size_t)m * k), oB((size_t)k * n), P((size_t)m * n);
        for (long l = 0; l < k; l++)
            for (long i = 0; i < m; i++)
                oA[i + l * m] = s.nota ? A[i + l * s.lda] : A[l + i * s.lda];
        for (long j = 0; j < n; j++)
            for (long l = 0; l < k; l++)
                oB[l + j * k] = s.notb ? B[l + j * s.ldb] : B[j + l * s.ldb];
        winograd_emul_backend<Ops> be;
        winograd<F>(be, m, n, k, oA.data(), m, oB.data(), k, P.data(), m, cutoff);
        for (long j = 0; j < n; j++)
            for (long i = 0; i < m; i++)
                gpu_combine<Ops>(s, sc, P.data(), C, i, j);
        return 0;
    }
    std::vector<F> Bs;
    const F *bop = B;
    if (s.nota) {
        Bs.resize((size_t)k * n + 1);
        for (long j = 0; j < n; j++)
            for (long l = 0; l < k; l++)
                gpu_alpha_opB<Ops>(s, sc.alpha, B, Bs.data(), l, j);
        bop = Bs.data();
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (long j = 0; j < n; j++)
        for (long i = 0; i < m; i++)
            gpu_element<Ops>(s, sc, A, bop, C, i, j);
    return 0;
}

template <class Ops> int gpu_device_count() { return 1; }

#endif

} // namespace mplapack_gemm

#endif
