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
 * Fixed-precision MPFR Level 3 BLAS on CUDA (see blas3_kernel_cuda.h): one
 * thread per lane.  With -DMPLAPACK_MPFR_CUDA_HOST_EMULATION the lanes run
 * on the host instead (tests without a GPU).
 */

#include <cstdio>
#include <cstdlib>
#include <vector>
#include "blas3_kernel_cuda.h"

#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
#include <cuda_runtime.h>
#endif

namespace mplapack_mpfr_cuda {

namespace {

template <int PB> struct fixed_args {
    l3_shape s;
    real_t<PB> *al, *be, *a, *b, *c;
};

template <int PB, bool TRSM> L3_HD bool fixed_lane(const fixed_args<PB> &f, long lane)
{
    typedef l3_fixed_engine<PB, TRSM> E;
    E e;
    blas_mats<E> x;
    x.al.p = f.al, x.al.ld = 1;
    x.be.p = f.be, x.be.ld = 1;
    x.a.p = f.a, x.a.ld = f.s.lda;
    x.b.p = f.b, x.b.ld = f.s.ldb;
    x.c.p = f.c, x.c.ld = f.s.ldc;
    if constexpr (TRSM)
        l3_trsm(e, f.s, x, lane);
    else
        blas_lane(e, f.s, x, lane);
    return e.ok;
}

#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
template <int PB, bool TRSM> __global__ void fixed_kernel(fixed_args<PB> f, long lanes, int *failed)
{
    const long lane = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (lane < lanes && !fixed_lane<PB, TRSM>(f, lane) && failed)
        atomicExch(failed, 1);
}

bool l3_failed(cudaError_t e, const char *what)
{
    if (e == cudaSuccess)
        return false;
    const char *v = std::getenv("MPLAPACK_MPFR_CUDA_VERBOSE");
    if (v != NULL && v[0] != '\0' && v[0] != '0')
        std::fprintf(stderr, "mplapack mpfr cuda (Level 3): %s: %s\n", what, cudaGetErrorString(e));
    return true;
}

template <class T> struct dbuf {
    T *p;
    dbuf() : p(NULL) {}
    ~dbuf() { cudaFree(p); }
    bool up(const T *h, size_t count)
    {
        if (l3_failed(cudaMalloc(&p, (count ? count : 1) * sizeof(T)), "cudaMalloc"))
            return false;
        return count == 0 || !l3_failed(cudaMemcpy(p, h, count * sizeof(T), cudaMemcpyHostToDevice), "copy to device");
    }

  private:
    dbuf(const dbuf &);
    dbuf &operator=(const dbuf &);
};
#endif

} // namespace

template <int PB> int blas3_fixed(const l3_shape &s, const real_t<PB> &alpha, const real_t<PB> &beta, const real_t<PB> *A, size_t countA, const real_t<PB> *B, size_t countB, real_t<PB> *C, size_t countC)
{
    const long lanes = blas_lanes(s);
#ifdef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    (void)countA, (void)countB, (void)countC;
    real_t<PB> al = alpha, be = beta;
    fixed_args<PB> f;
    f.s = s;
    f.al = &al, f.be = &be;
    f.a = const_cast<real_t<PB> *>(A), f.b = const_cast<real_t<PB> *>(B), f.c = C;
    int failed = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) reduction(|:failed)
#endif
    for (long lane = 0; lane < lanes; lane++)
        failed |= !(s.op == L3_TRSM ? fixed_lane<PB, true>(f, lane) : fixed_lane<PB, false>(f, lane));
    return failed;
#else
    int ndev = 0;
    if (l3_failed(cudaGetDeviceCount(&ndev), "cudaGetDeviceCount") || ndev == 0)
        return 1;
    // Only Rtrsm calls Newton. Leave other real BLAS launches' limits alone.
    if (s.op == L3_TRSM) {
        size_t cur = 0;
        if (l3_failed(cudaDeviceGetLimit(&cur, cudaLimitStackSize), "get stack") ||
            (cur < 65536 && l3_failed(cudaDeviceSetLimit(cudaLimitStackSize, 65536), "set stack")))
            return 2;
    }
    int failed = 0;
    dbuf<int> status;
    if (s.op == L3_TRSM && !status.up(&failed, 1))
        return 2;
    dbuf<real_t<PB>> dal, dbe, dA, dB, dC;
    if (!dal.up(&alpha, 1) || !dbe.up(&beta, 1) || !dA.up(A, countA) || !dB.up(B, countB) || !dC.up(C, countC))
        return 2;
    fixed_args<PB> f;
    f.s = s;
    f.al = dal.p, f.be = dbe.p, f.a = dA.p, f.b = dB.p, f.c = dC.p;
    const int threads = 64;
    if (lanes > 0) {
        if (s.op == L3_TRSM)
            fixed_kernel<PB, true><<<(unsigned)((lanes + threads - 1) / threads), threads>>>(f, lanes, status.p);
        else
            fixed_kernel<PB, false><<<(unsigned)((lanes + threads - 1) / threads), threads>>>(f, lanes, nullptr);
        if (l3_failed(cudaGetLastError(), "kernel launch") || l3_failed(cudaDeviceSynchronize(), "kernel"))
            return 3;
    }
    if (s.op == L3_TRSM &&
        (l3_failed(cudaMemcpy(&failed, status.p, sizeof(failed), cudaMemcpyDeviceToHost), "copy status") || failed))
        return 3;
    if (countC && l3_failed(cudaMemcpy(C, dC.p, countC * sizeof(real_t<PB>), cudaMemcpyDeviceToHost), "copy to host"))
        return 4;
    // Rswap and Rrot also write x (A), as the host emulation does in place
    if ((s.op == mplapack_blas::L1_SWAP || s.op == mplapack_blas::L1_ROT) && countA &&
        l3_failed(cudaMemcpy(const_cast<real_t<PB> *>(A), dA.p, countA * sizeof(real_t<PB>), cudaMemcpyDeviceToHost), "copy x to host"))
        return 4;
    return 0;
#endif
}

#define MPLAPACK_MPFR_CUDA_L3_INSTANTIATE(PB) template int blas3_fixed<PB>(const l3_shape &, const real_t<PB> &, const real_t<PB> &, const real_t<PB> *, size_t, const real_t<PB> *, size_t, real_t<PB> *, size_t);
MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(MPLAPACK_MPFR_CUDA_L3_INSTANTIATE)
#undef MPLAPACK_MPFR_CUDA_L3_INSTANTIATE

} // namespace mplapack_mpfr_cuda
