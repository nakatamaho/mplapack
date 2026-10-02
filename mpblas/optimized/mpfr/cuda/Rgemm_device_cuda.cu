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

} // namespace

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

template int gemm<512>(const gemm_shape &, const real_t<512> &, const real_t<512> &, const real_t<512> *, const real_t<512> *, real_t<512> *);
template int gemm<1024>(const gemm_shape &, const real_t<1024> &, const real_t<1024> &, const real_t<1024> *, const real_t<1024> *, real_t<1024> *);

} // namespace mplapack_mpfr_cuda
