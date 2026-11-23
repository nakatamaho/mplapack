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
 * CUDA implementation of mplapack_mpfr_cuda::gemm<PB> (PB = 256, 512, 768, 1024, 2048).
 * Shared operand tiles with one ordered accumulator per output element;
 * host emulation executes the same tile phases.
 */

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include "Rgemm_kernel_cuda.h"
#include "Rgemm_tiled_cuda.h"
#include "Rgemm_winograd_cuda.h"
#include "../../common/cuda_multi_gpu.h"
#include "../../common/Rgemm_winograd_multi_gpu.cuh"
#include <memory>
#include <vector>

namespace mplapack_mpfr_cuda {

namespace {

const int BLOCK_X = 16;
const int BLOCK_Y = 8;
const int BLOCK_1D = 128;
const long MAX_GRID_Y = 65535;

template <int PB> __global__ void alpha_opB_kernel(gemm_shape s, real_t<PB> alpha, const real_t<PB> *B, real_t<PB> *AB) {
    long total = s.k * s.n;
    for (long t = (long)blockIdx.x * blockDim.x + threadIdx.x; t < total; t += (long)gridDim.x * blockDim.x) {
        gemm_alpha_opB<PB>(s, alpha, B, AB, t % s.k, t / s.k);
    }
}

template <int PB> __global__ void gemm_kernel(gemm_shape s, real_t<PB> alpha, real_t<PB> beta, const real_t<PB> *A, const real_t<PB> *B, const real_t<PB> *AB, real_t<PB> *C) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= s.m)
        return;
    for (long j = (long)blockIdx.y * blockDim.y + threadIdx.y; j < s.n; j += (long)gridDim.y * blockDim.y) {
        gemm_element<PB>(s, alpha, beta, A, B, AB, C, i, j);
    }
}

// Few output elements (m*n small): a block owns SMALL_E elements of C and,
// for every SMALL_KR steps of k, its SMALL_E*SMALL_KR threads form the
// products of those steps in parallel into shared memory; then the owner of
// each element adds them in order.  Every product and every sum is the one
// of gemm_engine_element, in the same order, so the result is the same; but
// the sequential chain of an element holds only the additions, from shared
// memory (with one thread per element the chain of products and sums, with
// its loads, leaves the device mostly idle).  Without operand reuse it loses
// to the tiled kernels once m*n fills the device (small_gemm_max_mn).
const int SMALL_E = 32;
const int SMALL_KR = 4;
template <int PB, bool Nota> __global__ void __launch_bounds__(SMALL_E *SMALL_KR) small_gemm_kernel(gemm_shape s, real_t<PB> alpha, real_t<PB> beta, const real_t<PB> *A, const real_t<PB> *B, const real_t<PB> *AB, real_t<PB> *C) {
    typedef real_t<PB> F;
    __shared__ __align__(16) unsigned char storage[sizeof(F) * SMALL_KR * SMALL_E];
    F *sh = reinterpret_cast<F *>(storage); // sh[l * SMALL_E + e]
    const int e = threadIdx.x % SMALL_E, l = threadIdx.x / SMALL_E;
    const long g = (long)blockIdx.x * SMALL_E + e;
    const bool in = g < s.m * s.n, owner = in && l == 0;
    const long i = in ? g % s.m : 0, j = in ? g / s.m : 0;
    F v;
    v.set_zero();
    // op(A) = A: C(i,j) = beta*C(i,j), then += AB(l,j) * A(i,l)
    if (Nota && owner) {
        if (s.beta_one)
            v = C[i + j * s.ldc];
        else if (!s.beta_zero)
            v = beta * C[i + j * s.ldc];
    }
    for (long q = 0; q < s.k; q += SMALL_KR) {
        const int count = s.k - q < SMALL_KR ? (int)(s.k - q) : SMALL_KR;
        if (in && l < count) {
            // no zero skip, as the reference
            if (Nota)
                sh[l * SMALL_E + e] = AB[(q + l) + j * s.k] * A[i + (q + l) * s.lda];
            else
                sh[l * SMALL_E + e] = A[(q + l) + i * s.lda] * gemm_opB<PB>(s, B, q + l, j);
        }
        __syncthreads();
        if (owner)
            for (int u = 0; u < count; ++u)
                v = v + sh[u * SMALL_E + e];
        __syncthreads();
    }
    if (!owner)
        return;
    if (Nota) {
        C[i + j * s.ldc] = v;
    } else {
        // op(A) = A': C(i,j) = beta*C(i,j) + alpha*temp
        const F z = alpha * v;
        F c0;
        if (s.beta_zero)
            c0.set_zero();
        else if (s.beta_one)
            c0 = C[i + j * s.ldc];
        else
            c0 = beta * C[i + j * s.ldc];
        C[i + j * s.ldc] = c0 + z;
    }
}
// Largest m*n for small_gemm_kernel: where it was faster than the tiled
// kernels on the V100 (m = n = k sweep: 1.3-2.1 times below these sizes);
// MPLAPACK_MPFR_CUDA_GEMM_SMALL_MN=N overrides it (0: never).
template <int PB> long small_gemm_max_mn() {
    if (const char *v = std::getenv("MPLAPACK_MPFR_CUDA_GEMM_SMALL_MN")) {
        char *end = nullptr;
        const long x = std::strtol(v, &end, 10);
        if (end != v && !*end && x >= 0)
            return x;
    }
    return PB == 256 ? 128L * 128 : PB == 512 ? 96L * 96 : PB == 1024 ? 112L * 112 : 0;
}

bool report_errors()
{
    const char *v = std::getenv("MPLAPACK_MPFR_CUDA_VERBOSE");
    return v != NULL && v[0] != '\0' && v[0] != '0';
}

// Logs the CUDA error (when MPLAPACK_MPFR_CUDA_VERBOSE is set) and returns true on failure.
bool failed(cudaError_t e, const char *what)
{
    if (e == cudaSuccess)
        return false;
    if (report_errors())
        std::fprintf(stderr, "mplapack mpfr cuda: %s: %s\n", what, cudaGetErrorString(e));
    return true;
}

class device_buffer {
  public:
    device_buffer() : p_(NULL) {}
    ~device_buffer()
    {
        if (p_)
            cudaFree(p_);
    }
    cudaError_t allocate(size_t bytes) { return cudaMalloc(&p_, bytes == 0 ? 1 : bytes); }
    template <typename T> T *get() const { return static_cast<T *>(p_); }

  private:
    device_buffer(const device_buffer &);
    device_buffer &operator=(const device_buffer &);
    void *p_;
};

// ---- Winograd (Rgemm_winograd_cuda.h) ----

enum wg_op { WG_ZERO, WG_COPY, WG_ADD, WG_SUB };

template <int PB> __global__ void wg_elementwise_kernel(int op, long r, long c, real_t<PB> *Z, long ldz, const real_t<PB> *X, long ldx, const real_t<PB> *Y, long ldy) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= r)
        return;
    for (long j = (long)blockIdx.y * blockDim.y + threadIdx.y; j < c; j += (long)gridDim.y * blockDim.y) {
        if (op == WG_ZERO)
            wg_zero<PB>(Z, ldz, i, j);
        else if (op == WG_COPY)
            wg_copy<PB>(Z, ldz, X, ldx, i, j);
        else
            wg_addsub<PB>(Z, ldz, X, ldx, Y, ldy, op == WG_SUB, i, j);
    }
}

template <int PB> __global__ void wg_leaf_kernel(long m, long n, long k, real_t<PB> *C, long ldc, const real_t<PB> *A, long lda, const real_t<PB> *B, long ldb) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= m)
        return;
    for (long j = (long)blockIdx.y * blockDim.y + threadIdx.y; j < n; j += (long)gridDim.y * blockDim.y)
        wg_leaf<PB>(k, C, ldc, A, lda, B, ldb, i, j);
}

template <int PB> __global__ void wg_combine_kernel(gemm_shape s, real_t<PB> alpha, real_t<PB> beta, const real_t<PB> *P, real_t<PB> *C) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= s.m)
        return;
    for (long j = (long)blockIdx.y * blockDim.y + threadIdx.y; j < s.n; j += (long)gridDim.y * blockDim.y)
        wg_combine<PB>(s, alpha, beta, P, C, i, j);
}

dim3 grid_for(long r, long c)
{
    long gx = (r + BLOCK_X - 1) / BLOCK_X, gy = (c + BLOCK_Y - 1) / BLOCK_Y;
    if (gy > MAX_GRID_Y)
        gy = MAX_GRID_Y;
    return dim3((unsigned)(gx ? gx : 1), (unsigned)(gy ? gy : 1));
}

// Device backend: all matrices live in device memory; operations are queued
// on the default stream and checked once at the end.
template <int PB> class winograd_device_backend {
  public:
    typedef real_t<PB> F;
    winograd_device_backend() : failed_(false) {}
    // Stream-ordered allocations from the device's pool: the recursion
    // allocates at every node, and cudaMalloc/cudaFree would serialize the
    // devices of a multi-device call.
    F *alloc(size_t count)
    {
        void *p = NULL;
        if (failed_ || mplapack_mpfr_cuda::failed(cudaMallocAsync(&p, sizeof(F) * (count ? count : 1), 0), "cudaMallocAsync (Winograd)")) {
            failed_ = true;
            return NULL;
        }
        return static_cast<F *>(p);
    }
    void release(F *p)
    {
        if (p)
            cudaFreeAsync(p, 0);
    }
    bool failed() const { return failed_; }
    void zero(long r, long c, F *X, long ldx) { elementwise(WG_ZERO, r, c, X, ldx, NULL, 0, NULL, 0); }
    void copy(long r, long c, F *D, long ldd, const F *S, long lds) { elementwise(WG_COPY, r, c, D, ldd, S, lds, NULL, 0); }
    void addsub(long r, long c, F *Z, long ldz, const F *X, long ldx, const F *Y, long ldy, int sub) { elementwise(sub ? WG_SUB : WG_ADD, r, c, Z, ldz, X, ldx, Y, ldy); }
    void leaf(long m, long n, long k, F *C, long ldc, const F *A, long lda, const F *B, long ldb)
    {
        if (failed_ || m == 0 || n == 0)
            return;
        wg_leaf_kernel<PB><<<grid_for(m, n), dim3(BLOCK_X, BLOCK_Y)>>>(m, n, k, C, ldc, A, lda, B, ldb);
        check_launch("wg_leaf_kernel");
    }
    void check_launch(const char *what)
    {
        if (mplapack_mpfr_cuda::failed(cudaGetLastError(), what))
            failed_ = true;
    }

  private:
    void elementwise(int op, long r, long c, F *Z, long ldz, const F *X, long ldx, const F *Y, long ldy)
    {
        if (failed_ || r == 0 || c == 0)
            return;
        wg_elementwise_kernel<PB><<<grid_for(r, c), dim3(BLOCK_X, BLOCK_Y)>>>(op, r, c, Z, ldz, X, ldx, Y, ldy);
        check_launch("wg_elementwise_kernel");
    }
    bool failed_;
};

} // namespace

// Winograd with the seven top-level products over G devices
// (../../common/Rgemm_winograd_multi_gpu.cuh); the same result as on one.
template <int PB> static int gemm_winograd_multi(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, real_t<PB> *C, long cutoff, int G)
{
    typedef real_t<PB> F;
    typedef winograd_device_backend<PB> Backend;
    using namespace mplapack_cuda_multi;
    device_restore restore;
    const size_t bytesA = sizeof(F) * (size_t)s.m * s.k, bytesB = sizeof(F) * (size_t)s.k * s.n, bytesC = sizeof(F) * (size_t)s.m * s.n;
    enable_peer_access(G);
    std::vector<std::unique_ptr<device_buffer> > dA((size_t)G), dB((size_t)G);
    std::vector<std::unique_ptr<Backend> > be((size_t)G);
    std::vector<Backend *> bp((size_t)G);
    std::vector<const F *> pA((size_t)G), pB((size_t)G);
    device_buffer dC, dP;
    int rc = 0;
    for (int d = 0; d < G && !rc; d++) {
        if (failed(cudaSetDevice(d), "cudaSetDevice")) {
            rc = 1;
            break;
        }
        dA[d].reset(new device_buffer);
        dB[d].reset(new device_buffer);
        be[d].reset(new Backend);
        bp[d] = be[d].get();
        if (failed(dA[d]->allocate(bytesA), "cudaMalloc A") || failed(dB[d]->allocate(bytesB), "cudaMalloc B"))
            rc = 2;
        else if (failed(cudaMemcpy(dA[d]->get<F>(), A, bytesA, cudaMemcpyHostToDevice), "copy A") || failed(cudaMemcpy(dB[d]->get<F>(), B, bytesB, cudaMemcpyHostToDevice), "copy B"))
            rc = 3;
        else
            pA[d] = dA[d]->get<F>(), pB[d] = dB[d]->get<F>();
    }
    if (!rc && failed(cudaSetDevice(0), "cudaSetDevice"))
        rc = 1;
    if (!rc && (failed(dC.allocate(bytesC), "cudaMalloc C") || failed(dP.allocate(bytesC), "cudaMalloc P")))
        rc = 2;
    if (!rc && !s.beta_zero && failed(cudaMemcpy(dC.get<F>(), C, bytesC, cudaMemcpyHostToDevice), "copy C"))
        rc = 3;
    if (!rc && mplapack_gemm::winograd_multi<F>(G, bp.data(), s.m, s.n, s.k, pA.data(), pB.data(), dP.get<F>(), [cutoff](Backend &b, long m, long n, long k, const F *a, long lda, const F *bb, long ldb, F *c, long ldc) { winograd<PB>(b, m, n, k, a, lda, bb, ldb, c, ldc, cutoff); }))
        rc = 4;
    if (!rc) {
        cudaSetDevice(0);
        wg_combine_kernel<PB><<<grid_for(s.m, s.n), dim3(BLOCK_X, BLOCK_Y)>>>(s, alpha, beta, dP.get<F>(), dC.get<F>());
        if (failed(cudaGetLastError(), "wg_combine_kernel launch") || failed(cudaDeviceSynchronize(), "Winograd kernels"))
            rc = 4;
        else if (failed(cudaMemcpy(C, dC.get<F>(), bytesC, cudaMemcpyDeviceToHost), "copy C back"))
            rc = 5;
    }
    for (int d = 0; d < G; d++) {
        cudaSetDevice(d);
        dA[d].reset();
        dB[d].reset();
    }
    cudaSetDevice(0);
    return rc;
}

template <int PB> int gemm_winograd(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, real_t<PB> *C, long cutoff)
{
    typedef real_t<PB> F;
    int ndev = 0;
    if (failed(cudaGetDeviceCount(&ndev), "cudaGetDeviceCount") || ndev == 0)
        return 1;
    {
        using namespace mplapack_cuda_multi;
        int G = devices(env_gpus("MPLAPACK_MPFR"));
        if (G > 7)
            G = 7;
        if (G > 1 && mplapack_gemm::winograd_multi_applies(s.m, s.n, s.k, cutoff, mplapack_gemm::winograd_multi_min("MPLAPACK_MPFR")))
            return gemm_winograd_multi<PB>(s, alpha, beta, A, B, C, cutoff, G);
    }
    const size_t bytesA = sizeof(F) * (size_t)s.m * s.k, bytesB = sizeof(F) * (size_t)s.k * s.n, bytesC = sizeof(F) * (size_t)s.m * s.n;
    device_buffer dA, dB, dC, dP;
    if (failed(dA.allocate(bytesA), "cudaMalloc A") || failed(dB.allocate(bytesB), "cudaMalloc B") || failed(dC.allocate(bytesC), "cudaMalloc C") || failed(dP.allocate(bytesC), "cudaMalloc P"))
        return 2;
    if (failed(cudaMemcpy(dA.get<F>(), A, bytesA, cudaMemcpyHostToDevice), "copy A") || failed(cudaMemcpy(dB.get<F>(), B, bytesB, cudaMemcpyHostToDevice), "copy B"))
        return 3;
    if (!s.beta_zero && failed(cudaMemcpy(dC.get<F>(), C, bytesC, cudaMemcpyHostToDevice), "copy C"))
        return 3;
    winograd_device_backend<PB> be;
    winograd<PB>(be, s.m, s.n, s.k, dA.get<F>(), s.m, dB.get<F>(), s.k, dP.get<F>(), s.m, cutoff);
    if (be.failed())
        return 4;
    wg_combine_kernel<PB><<<grid_for(s.m, s.n), dim3(BLOCK_X, BLOCK_Y)>>>(s, alpha, beta, dP.get<F>(), dC.get<F>());
    if (failed(cudaGetLastError(), "wg_combine_kernel launch") || failed(cudaDeviceSynchronize(), "Winograd kernels"))
        return 4;
    if (failed(cudaMemcpy(C, dC.get<F>(), bytesC, cudaMemcpyDeviceToHost), "copy C back"))
        return 5;
    return 0;
}

// One device (the current one): C = alpha op(A) op(B) + beta C for the
// shape s, whose leading dimensions are those of the packed device copies
// (rows of A, B and C).  The host operands have the leading dimensions
// hlda, hldb, hldc (a column block of a larger matrix).
template <int PB> static int gemm_one(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, long hlda, const real_t<PB> *B, long hldb, real_t<PB> *C, long hldc)
{
    typedef real_t<PB> F;
    const long nrowa = s.nota ? s.m : s.k, ncola = s.nota ? s.k : s.m;
    const long nrowb = s.notb ? s.k : s.n, ncolb = s.notb ? s.n : s.k;
    const size_t bytesA = sizeof(F) * (size_t)nrowa * ncola;
    const size_t bytesB = sizeof(F) * (size_t)nrowb * ncolb;
    const size_t bytesC = sizeof(F) * (size_t)s.m * s.n;
    // Volta (Rgemm_tiled_cuda.h): every form through the Lockstep layout,
    // which forms alpha*op(B) as it loads the tiles (one kernel, no AB).
    // Below its size threshold the default tiled layout, also one kernel for
    // NoTrans A: there the original two-kernel path runs at half the speed
    // of the tiled transposed form.
    bool lockstep = false, small = false;
    const bool few = s.k > 0 && s.m * s.n <= small_gemm_max_mn<PB>();
    if constexpr (PB == 256 || PB == 512 || PB == 768) {
        static_assert(value_gemm_lockstep_fits<F>(), "Lockstep tiles exceed static shared memory");
        if (s.k > 0 && value_gemm_lockstep()) {
            lockstep = value_gemm_lockstep_shape<F>(s);
            small = !lockstep;
        }
    }
    if (few)
        lockstep = small = false;
    const size_t bytesAB = s.nota && !lockstep && !small ? sizeof(F) * (size_t)s.k * s.n : 0;

    device_buffer dA, dB, dC, dAB;
    if (failed(dA.allocate(bytesA), "cudaMalloc A") || failed(dB.allocate(bytesB), "cudaMalloc B") || failed(dC.allocate(bytesC), "cudaMalloc C") || failed(dAB.allocate(bytesAB), "cudaMalloc alpha*op(B)"))
        return 2;
    if (failed(cudaMemcpy2D(dA.get<F>(), sizeof(F) * nrowa, A, sizeof(F) * hlda, sizeof(F) * nrowa, ncola, cudaMemcpyHostToDevice), "copy A") || failed(cudaMemcpy2D(dB.get<F>(), sizeof(F) * nrowb, B, sizeof(F) * hldb, sizeof(F) * nrowb, ncolb, cudaMemcpyHostToDevice), "copy B"))
        return 3;
    if (!s.beta_zero && failed(cudaMemcpy2D(dC.get<F>(), sizeof(F) * s.m, C, sizeof(F) * hldc, sizeof(F) * s.m, s.n, cudaMemcpyHostToDevice), "copy C"))
        return 3;

    if (s.nota && s.k > 0 && !lockstep && !small) {
        long total = s.k * s.n;
        long blocks = (total + BLOCK_1D - 1) / BLOCK_1D;
        if (blocks > 65535)
            blocks = 65535;
        alpha_opB_kernel<PB><<<(unsigned)blocks, BLOCK_1D>>>(s, alpha, dB.get<F>(), dAB.get<F>());
        if (failed(cudaGetLastError(), "alpha_opB_kernel launch"))
            return 4;
    }
    long gx = (s.m + BLOCK_X - 1) / BLOCK_X;
    long gy = (s.n + BLOCK_Y - 1) / BLOCK_Y;
    if (gy > MAX_GRID_Y)
        gy = MAX_GRID_Y;
    dim3 grid((unsigned)gx, (unsigned)gy), block(BLOCK_X, BLOCK_Y);
    cudaError_t launch_error=cudaSuccess;
    bool tiled=false;
    if(few) {
        const unsigned nblk=(unsigned)((s.m*s.n+SMALL_E-1)/SMALL_E);
        if(s.nota) small_gemm_kernel<PB,true><<<nblk,SMALL_E*SMALL_KR>>>(s,alpha,beta,dA.get<F>(),dB.get<F>(),dAB.get<F>(),dC.get<F>());
        else small_gemm_kernel<PB,false><<<nblk,SMALL_E*SMALL_KR>>>(s,alpha,beta,dA.get<F>(),dB.get<F>(),dAB.get<F>(),dC.get<F>());
        launch_error=cudaGetLastError();
        tiled=true;
    }
    if constexpr(PB==256 || PB==512 || PB==768) if(!tiled) {
        if(lockstep) {
            if(s.nota) launch_error=(cudaError_t)value_gemm_device(make_value_gemm_job<PB,true,true,false>(s,alpha,beta,dA.get<F>(),dB.get<F>(),dAB.get<F>(),dC.get<F>()),0);
            else launch_error=(cudaError_t)value_gemm_device(make_value_gemm_job<PB,false,true>(s,alpha,beta,dA.get<F>(),dB.get<F>(),dAB.get<F>(),dC.get<F>()),0);
            tiled=true;
        } else if(small) {
            if(s.nota) launch_error=(cudaError_t)value_gemm_device(make_value_gemm_job<PB,true,false,false>(s,alpha,beta,dA.get<F>(),dB.get<F>(),dAB.get<F>(),dC.get<F>()),0);
            else launch_error=(cudaError_t)value_gemm_device(make_value_gemm_job<PB,false>(s,alpha,beta,dA.get<F>(),dB.get<F>(),dAB.get<F>(),dC.get<F>()),0);
            tiled=true;
        }
    }
    if constexpr(PB==256 || PB==512 || PB==1024) if(!tiled) {
        // Shared loads alone do not improve 256/512-bit NoTrans; retain the
        // original kernel there and at the untuned precisions.
        if(s.k>0 && (PB==1024 || !s.nota)) {
            if constexpr(PB==1024) {
                if(s.nota) launch_error=(cudaError_t)value_gemm_device(make_value_gemm_job<PB,true>(s,alpha,beta,dA.get<F>(),dB.get<F>(),dAB.get<F>(),dC.get<F>()),0);
                else launch_error=(cudaError_t)value_gemm_device(make_value_gemm_job<PB,false>(s,alpha,beta,dA.get<F>(),dB.get<F>(),dAB.get<F>(),dC.get<F>()),0);
            } else launch_error=(cudaError_t)value_gemm_device(make_value_gemm_job<PB,false>(s,alpha,beta,dA.get<F>(),dB.get<F>(),dAB.get<F>(),dC.get<F>()),0);
            tiled=true;
        }
    }
    if(!tiled) {
        gemm_kernel<PB><<<grid, block>>>(s, alpha, beta, dA.get<F>(), dB.get<F>(), dAB.get<F>(), dC.get<F>());
        launch_error=cudaGetLastError();
    }
    if (failed(launch_error, "gemm_kernel launch") || failed(cudaDeviceSynchronize(), "gemm_kernel"))
        return 4;
    if (failed(cudaMemcpy2D(C, sizeof(F) * hldc, dC.get<F>(), sizeof(F) * s.m, sizeof(F) * s.m, s.n, cudaMemcpyDeviceToHost), "copy C back"))
        return 5;
    return 0;
}

// The columns of C (and op(B)) are cut over the devices
// (../../common/cuda_multi_gpu.h, MPLAPACK_MPFR_CUDA_NGPUS); every element is
// computed as on one device, so the result does not depend on their number.
// A block keeps at least 128 columns, so small calls stay on one device.
template <int PB> int gemm(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, real_t<PB> *C)
{
    using namespace mplapack_cuda_multi;
    const int ndev = devices(env_gpus("MPLAPACK_MPFR"));
    if (ndev == 0) {
        failed(cudaErrorNoDevice, "cudaGetDeviceCount");
        return 1;
    }
    const int nb = s.k > 0 ? blocks(s.n, ndev, 128) : 1;
    return for_column_blocks(s.n, nb, 16, [&](long j0, long j1) {
        gemm_shape t = s;
        t.n = j1 - j0;
        t.ldb = s.notb ? s.k : t.n;
        t.ldc = s.m;
        return gemm_one<PB>(t, alpha, beta, A, s.lda, s.notb ? B + j0 * s.ldb : B + j0, s.ldb, C + j0 * s.ldc, s.ldc);
    });
}

int device_count()
{
    int ndev = 0;
    if (cudaGetDeviceCount(&ndev) != cudaSuccess)
        return 0;
    return ndev;
}

#define MPLAPACK_MPFR_CUDA_INSTANTIATE_WINOGRAD(PB) template int gemm_winograd<PB>(const gemm_shape &, const real_t<PB> &, const real_t<PB> &, const real_t<PB> *, const real_t<PB> *, real_t<PB> *, long);
MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(MPLAPACK_MPFR_CUDA_INSTANTIATE_WINOGRAD)

#define MPLAPACK_MPFR_CUDA_INSTANTIATE_GEMM(PB) template int gemm<PB>(const gemm_shape &, const real_t<PB> &, const real_t<PB> &, const real_t<PB> *, const real_t<PB> *, real_t<PB> *);
MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(MPLAPACK_MPFR_CUDA_INSTANTIATE_GEMM)

} // namespace mplapack_mpfr_cuda
