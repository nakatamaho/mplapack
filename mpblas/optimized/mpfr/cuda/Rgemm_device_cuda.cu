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
 * CUDA implementation of mplapack_mpfr_cuda::gemm<PB> (PB = 512, 1024).
 * One thread per output element; the per-element arithmetic is in
 * Rgemm_kernel_cuda.h and is shared with the host implementation.
 */

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include "Rgemm_kernel_cuda.h"
#include "Rgemm_winograd_cuda.h"

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
    F *alloc(size_t count)
    {
        void *p = NULL;
        if (failed_ || mplapack_mpfr_cuda::failed(cudaMalloc(&p, sizeof(F) * (count ? count : 1)), "cudaMalloc (Winograd)")) {
            failed_ = true;
            return NULL;
        }
        return static_cast<F *>(p);
    }
    void release(F *p)
    {
        if (p)
            cudaFree(p);
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

template <int PB> int gemm_winograd(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, real_t<PB> *C, long cutoff)
{
    typedef real_t<PB> F;
    int ndev = 0;
    if (failed(cudaGetDeviceCount(&ndev), "cudaGetDeviceCount") || ndev == 0)
        return 1;
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

template <int PB> int gemm(const gemm_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, const real_t<PB> *B, real_t<PB> *C)
{
    typedef real_t<PB> F;
    int ndev = 0;
    if (failed(cudaGetDeviceCount(&ndev), "cudaGetDeviceCount") || ndev == 0)
        return 1;

    const long nrowa = s.nota ? s.m : s.k, ncola = s.nota ? s.k : s.m;
    const long nrowb = s.notb ? s.k : s.n, ncolb = s.notb ? s.n : s.k;
    const size_t bytesA = sizeof(F) * (size_t)nrowa * ncola;
    const size_t bytesB = sizeof(F) * (size_t)nrowb * ncolb;
    const size_t bytesC = sizeof(F) * (size_t)s.m * s.n;
    const size_t bytesAB = s.nota ? sizeof(F) * (size_t)s.k * s.n : 0;

    device_buffer dA, dB, dC, dAB;
    if (failed(dA.allocate(bytesA), "cudaMalloc A") || failed(dB.allocate(bytesB), "cudaMalloc B") || failed(dC.allocate(bytesC), "cudaMalloc C") || failed(dAB.allocate(bytesAB), "cudaMalloc alpha*op(B)"))
        return 2;
    if (failed(cudaMemcpy(dA.get<F>(), A, bytesA, cudaMemcpyHostToDevice), "copy A") || failed(cudaMemcpy(dB.get<F>(), B, bytesB, cudaMemcpyHostToDevice), "copy B"))
        return 3;
    if (!s.beta_zero && failed(cudaMemcpy(dC.get<F>(), C, bytesC, cudaMemcpyHostToDevice), "copy C"))
        return 3;

    if (s.nota && s.k > 0) {
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
    gemm_kernel<PB><<<grid, block>>>(s, alpha, beta, dA.get<F>(), dB.get<F>(), dAB.get<F>(), dC.get<F>());
    if (failed(cudaGetLastError(), "gemm_kernel launch") || failed(cudaDeviceSynchronize(), "gemm_kernel"))
        return 4;
    if (failed(cudaMemcpy(C, dC.get<F>(), bytesC, cudaMemcpyDeviceToHost), "copy C back"))
        return 5;
    return 0;
}

int device_count()
{
    int ndev = 0;
    if (cudaGetDeviceCount(&ndev) != cudaSuccess)
        return 0;
    return ndev;
}

template int gemm_winograd<512>(const gemm_shape &, const real_t<512> &, const real_t<512> &, const real_t<512> *, const real_t<512> *, real_t<512> *, long);
template int gemm_winograd<1024>(const gemm_shape &, const real_t<1024> &, const real_t<1024> &, const real_t<1024> *, const real_t<1024> *, real_t<1024> *, long);

template int gemm<512>(const gemm_shape &, const real_t<512> &, const real_t<512> &, const real_t<512> *, const real_t<512> *, real_t<512> *);
template int gemm<1024>(const gemm_shape &, const real_t<1024> &, const real_t<1024> &, const real_t<1024> *, const real_t<1024> *, real_t<1024> *);

} // namespace mplapack_mpfr_cuda
