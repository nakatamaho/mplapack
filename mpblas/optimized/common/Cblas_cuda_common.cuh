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
 * dd/td/qd complex Level 1/2/3 BLAS (and Rgemmtr) on CUDA: every thread runs
 * the reference routine (blas_dev/<R>_dev.h, the reference text on the
 * device types of Rblas_dev_types.h) on its part of the output: an element
 * (Cgemm, Cgeru, Cgerc, Level 1), a column or row (Csymm, Chemm, Ctrmm,
 * Ctrsm, Cgemv), or a column of the triangle (Csyrk, Cherk, Csyr2k, Cher2k,
 * Cgemmtr, Rgemmtr, Cher, Cher2), so the results are the reference results
 * bit for bit.  The host side is ../common/blas_cuda/Cblas_bridge_cuda.cpp;
 * instantiated per type in <type>/cuda/Rgemm_gpu_cuda.cu.  With
 * MPLAPACK_GEMM_CUDA_HOST_EMULATION the threads run on the host.
 */

#ifndef _MPLAPACK_CBLAS_CUDA_COMMON_CUH_
#define _MPLAPACK_CBLAS_CUDA_COMMON_CUH_

#include <vector>
#include "blas_dev/blas_dev.h"
#include "Rblas_gpu_operand_common.h"

// MPLAPACK_GEMM_CUDA_HOST_EMULATION (dd/td/qd) or
// MPLAPACK_MPFR_CUDA_HOST_EMULATION (MPFR): the jobs run on the host.
#if defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION) || defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION)
#define MPLAPACK_CBLAS_CUDA_DEVICE 0
#elif defined(__CUDACC__)
#define MPLAPACK_CBLAS_CUDA_DEVICE 1
#include <cuda_runtime.h>
#else
#error "Cblas_cuda_common.cuh: compile with nvcc, or define MPLAPACK_GEMM_CUDA_HOST_EMULATION"
#endif

#include "Cblas_cuda_jobs_common.h"

namespace mplapack_blas_dev {

using mplapack_blas::gpu_operand;

#if MPLAPACK_CBLAS_CUDA_DEVICE
template <class Job> __global__ void k_job(Job job, long units)
{
    const long u = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (u < units)
        job(u);
}
#endif

// Uploads the operands (spans), runs job(u) for u < units, downloads the
// written ones.  ops[q].host is the logical element (0, 0); bind() gets the
// device addresses of these elements.
template <class E, class Job> int run_job(Job job, long units, const gpu_operand<E> *ops, int nops)
{
    E *dev[4] = {0, 0, 0, 0};
#if MPLAPACK_CBLAS_CUDA_DEVICE
    E *buf[4] = {0, 0, 0, 0};
    cudaError_t e = cudaSuccess;
    for (int q = 0; q < nops && e == cudaSuccess; q++) {
        const size_t count = (size_t)(ops[q].count ? ops[q].count : 1);
        e = cudaMalloc(&buf[q], count * sizeof(E));
        if (e == cudaSuccess && ops[q].count)
            e = cudaMemcpy(buf[q], ops[q].host + ops[q].lo, ops[q].count * sizeof(E), cudaMemcpyHostToDevice);
        dev[q] = buf[q] - ops[q].lo;
    }
    if (e == cudaSuccess && units > 0) {
        job.bind(dev);
        const int threads = 64;
        k_job<Job><<<(unsigned)((units + threads - 1) / threads), threads>>>(job, units);
        e = cudaGetLastError();
        if (e == cudaSuccess)
            e = cudaDeviceSynchronize();
    }
    for (int q = 0; q < nops && e == cudaSuccess; q++)
        if (ops[q].download && ops[q].count)
            e = cudaMemcpy(const_cast<E *>(ops[q].host) + ops[q].lo, buf[q], ops[q].count * sizeof(E), cudaMemcpyDeviceToHost);
    for (int q = 0; q < nops; q++)
        cudaFree(buf[q]);
    return (int)e;
#else
    std::vector<E> buf[4];
    for (int q = 0; q < nops; q++) {
        buf[q].assign(ops[q].host + ops[q].lo, ops[q].host + ops[q].lo + ops[q].count);
        if (buf[q].empty())
            buf[q].resize(1);
        dev[q] = buf[q].data() - ops[q].lo;
    }
    job.bind(dev);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (long u = 0; u < units; u++)
        job(u);
    for (int q = 0; q < nops; q++)
        if (ops[q].download)
            for (long t = 0; t < ops[q].count; t++)
                const_cast<E *>(ops[q].host)[ops[q].lo + t] = buf[q][(size_t)t];
    return 0;
#endif
}

template <class Job> int cgpu_run(Job job, const mplapack_blas::gpu_operand<typename Job::E> *ops, int nops) { return run_job<typename Job::E>(job, job.units(), ops, nops); }

} // namespace mplapack_blas_dev

#endif
