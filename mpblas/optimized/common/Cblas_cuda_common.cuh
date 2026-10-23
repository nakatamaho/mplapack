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
// The job is read from device memory: its scalars may exceed the size of the
// kernel parameters (MPFR runtime precision).
template <class Job> __global__ void k_job(const Job *job, long units)
{
    const long u = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (u < units)
        (*job)(u);
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
    Job *djob = 0;
    if (e == cudaSuccess && units > 0) {
        job.bind(dev);
        e = cudaMalloc(&djob, sizeof(Job));
        if (e == cudaSuccess)
            e = cudaMemcpy(djob, &job, sizeof(Job), cudaMemcpyHostToDevice);
        const int threads = 64;
        if (e == cudaSuccess)
            k_job<Job><<<(unsigned)((units + threads - 1) / threads), threads>>>(djob, units);
        if (e == cudaSuccess)
            e = cudaGetLastError();
        if (e == cudaSuccess)
            e = cudaDeviceSynchronize();
    }
    cudaFree(djob);
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

// ---- Cgemm tiled (dd/td/qd): a CG_TILE x CG_TILE block of threads computes a
// tile of C, the columns (rows) l of op(A) and op(B) in chunks of CG_TILE held
// in shared memory.  Every thread performs the operations of the reference
// Cgemm on its element in the reference order (transa = N: c = beta*c, then
// c += (alpha*op(b)(l,j))*a(i,l); otherwise temp += op(a)(l,i)*op(b)(l,j) and
// c = alpha*temp [+ beta*c]), so the result is the same bit for bit.  The
// phases are separate functions so that the host emulation runs them for all
// the threads of a block in turn.

static const int CG_TILE = 16;

template <class Ops> struct cgemm_tile {
    typedef typename cgpu_types<Ops>::DC DC;
    struct shared {
        DC a[CG_TILE][CG_TILE + 1]; // a[l][t]: op(A)(row i0 + t, l) as used
        DC b[CG_TILE][CG_TILE + 1]; // b[l][t]: op(B)(l, column j0 + t)
    };
    const job_Cgemm<Ops> &j;
    bool nota, conja, notb, conjb;
    MPLAPACK_DEV_HD explicit cgemm_tile(const job_Cgemm<Ops> &job) : j(job)
    {
        nota = lsame(j.ta, "N"), conja = lsame(j.ta, "C");
        notb = lsame(j.tb, "N"), conjb = lsame(j.tb, "C");
    }
    MPLAPACK_DEV_HD bool mine(long i, long jj) const { return i < j.m && jj < j.n; }
    MPLAPACK_DEV_HD void init(DC &acc, long i, long jj) const
    {
        const DC zero = DC(0.0, 0.0), one = DC(1.0, 0.0);
        if (!mine(i, jj))
            return;
        if (nota) {
            if (j.beta == zero)
                acc = zero;
            else if (j.beta != one)
                acc = j.beta * j.C[i + jj * j.ldc];
            else
                acc = j.C[i + jj * j.ldc];
        } else {
            acc = zero; // temp = zero
        }
    }
    MPLAPACK_DEV_HD void load(shared &sh, long i0, long j0, long p0, int ti, int tj) const
    {
        const long l = p0 + tj;
        if (l >= j.k)
            return;
        const long i = i0 + ti, jj = j0 + ti;
        if (i < j.m)
            sh.a[tj][ti] = nota ? j.A[i + l * j.lda] : (conja ? conj(j.A[l + i * j.lda]) : j.A[l + i * j.lda]);
        if (jj < j.n)
            sh.b[tj][ti] = notb ? j.B[l + jj * j.ldb] : (conjb ? conj(j.B[jj + l * j.ldb]) : j.B[jj + l * j.ldb]);
    }
    MPLAPACK_DEV_HD void step(const shared &sh, DC &acc, long i, long jj, long kc, int ti, int tj) const
    {
        if (!mine(i, jj))
            return;
        DC temp = 0.0;
        for (long l = 0; l < kc; l++) {
            if (nota) {
                temp = j.alpha * sh.b[l][tj];  // temp = alpha*op(b)(l,j)
                acc += temp * sh.a[l][ti];     // c(i,j) += temp*a(i,l)
            } else {
                acc += sh.a[l][ti] * sh.b[l][tj]; // temp += op(a)(l,i)*op(b)(l,j)
            }
        }
    }
    MPLAPACK_DEV_HD void finish(const DC &acc, long i, long jj) const
    {
        const DC zero = DC(0.0, 0.0);
        if (!mine(i, jj))
            return;
        DC &cij = j.C[i + jj * j.ldc];
        if (nota)
            cij = acc;
        else if (j.beta == zero)
            cij = j.alpha * acc;
        else
            cij = j.alpha * acc + j.beta * cij;
    }
};

// the tiled form: dd/td/qd elements (small enough for shared memory), alpha != 0
template <class Job> struct cgpu_is_cgemm {
    static const bool value = false;
};
template <class Ops> struct cgpu_is_cgemm<job_Cgemm<Ops> > {
    static const bool value = sizeof(typename cgpu_types<Ops>::DC) <= 64;
};

#if MPLAPACK_CBLAS_CUDA_DEVICE
template <class Ops> __global__ void k_cgemm_tiled(const job_Cgemm<Ops> *job)
{
    typedef cgemm_tile<Ops> T;
    __shared__ typename T::shared sh;
    const T t(*job);
    const int ti = threadIdx.x, tj = threadIdx.y;
    const long i0 = (long)blockIdx.x * CG_TILE, j0 = (long)blockIdx.y * CG_TILE;
    const long i = i0 + ti, jj = j0 + tj;
    typename T::DC acc;
    t.init(acc, i, jj);
    for (long p0 = 0; p0 < job->k; p0 += CG_TILE) {
        t.load(sh, i0, j0, p0, ti, tj);
        __syncthreads();
        const long kc = job->k - p0 < CG_TILE ? job->k - p0 : CG_TILE;
        t.step(sh, acc, i, jj, kc, ti, tj);
        __syncthreads();
    }
    t.finish(acc, i, jj);
}
#endif

template <class Ops> int cgemm_tiled_run(job_Cgemm<Ops> job, const gpu_operand<typename job_Cgemm<Ops>::E> *ops, int nops)
{
    typedef typename job_Cgemm<Ops>::E E;
    typedef cgemm_tile<Ops> T;
    E *dev[4] = {0, 0, 0, 0};
#if MPLAPACK_CBLAS_CUDA_DEVICE
    E *buf[4] = {0, 0, 0, 0};
    job_Cgemm<Ops> *djob = 0;
    cudaError_t e = cudaSuccess;
    for (int q = 0; q < nops && e == cudaSuccess; q++) {
        const size_t count = (size_t)(ops[q].count ? ops[q].count : 1);
        e = cudaMalloc(&buf[q], count * sizeof(E));
        if (e == cudaSuccess && ops[q].count)
            e = cudaMemcpy(buf[q], ops[q].host + ops[q].lo, ops[q].count * sizeof(E), cudaMemcpyHostToDevice);
        dev[q] = buf[q] - ops[q].lo;
    }
    if (e == cudaSuccess) {
        job.bind(dev);
        e = cudaMalloc(&djob, sizeof job);
        if (e == cudaSuccess)
            e = cudaMemcpy(djob, &job, sizeof job, cudaMemcpyHostToDevice);
        const dim3 grid((unsigned)((job.m + CG_TILE - 1) / CG_TILE), (unsigned)((job.n + CG_TILE - 1) / CG_TILE));
        if (e == cudaSuccess)
            k_cgemm_tiled<Ops><<<grid, dim3(CG_TILE, CG_TILE)>>>(djob);
        if (e == cudaSuccess)
            e = cudaGetLastError();
        if (e == cudaSuccess)
            e = cudaDeviceSynchronize();
    }
    for (int q = 0; q < nops && e == cudaSuccess; q++)
        if (ops[q].download && ops[q].count)
            e = cudaMemcpy(const_cast<E *>(ops[q].host) + ops[q].lo, buf[q], ops[q].count * sizeof(E), cudaMemcpyDeviceToHost);
    cudaFree(djob);
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
    const T t(job);
    const long nbi = (job.m + CG_TILE - 1) / CG_TILE, nbj = (job.n + CG_TILE - 1) / CG_TILE;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (long b = 0; b < nbi * nbj; b++) {
        const long i0 = (b % nbi) * CG_TILE, j0 = (b / nbi) * CG_TILE;
        std::vector<typename T::shared> shv(1);
        std::vector<typename T::DC> acc((size_t)(CG_TILE * CG_TILE));
        for (int tj = 0; tj < CG_TILE; tj++)
            for (int ti = 0; ti < CG_TILE; ti++)
                t.init(acc[(size_t)(ti + tj * CG_TILE)], i0 + ti, j0 + tj);
        for (long p0 = 0; p0 < job.k; p0 += CG_TILE) {
            for (int tj = 0; tj < CG_TILE; tj++)
                for (int ti = 0; ti < CG_TILE; ti++)
                    t.load(shv[0], i0, j0, p0, ti, tj);
            const long kc = job.k - p0 < CG_TILE ? job.k - p0 : CG_TILE;
            for (int tj = 0; tj < CG_TILE; tj++)
                for (int ti = 0; ti < CG_TILE; ti++)
                    t.step(shv[0], acc[(size_t)(ti + tj * CG_TILE)], i0 + ti, j0 + tj, kc, ti, tj);
        }
        for (int tj = 0; tj < CG_TILE; tj++)
            for (int ti = 0; ti < CG_TILE; ti++)
                t.finish(acc[(size_t)(ti + tj * CG_TILE)], i0 + ti, j0 + tj);
    }
    for (int q = 0; q < nops; q++)
        if (ops[q].download)
            for (long u = 0; u < ops[q].count; u++)
                const_cast<E *>(ops[q].host)[ops[q].lo + u] = buf[q][(size_t)u];
    return 0;
#endif
}

template <class Job> int cgpu_run(Job job, const mplapack_blas::gpu_operand<typename Job::E> *ops, int nops)
{
    if constexpr (cgpu_is_cgemm<Job>::value) {
        typedef typename Job::DC DC;
        if (job.m > 0 && job.n > 0 && job.k > 0 && !(job.alpha == DC(0.0, 0.0)))
            return cgemm_tiled_run(job, ops, nops);
    }
    return run_job<typename Job::E>(job, job.units(), ops, nops);
}

} // namespace mplapack_blas_dev

#endif
