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
 * dd/td/qd Level 1/2/3 BLAS on CUDA: the lanes of Rblas_lanes_common.h, one
 * thread per lane, with the arithmetic of Rgemm_cuda_ops_common.h (libQD3's
 * branch-free operations), so that the results are bit-identical to the CPU
 * routines.  Instantiated per type in <type>/cuda/Rgemm_gpu_cuda.cu.  With
 * MPLAPACK_GEMM_CUDA_HOST_EMULATION and a C++ compiler the lanes run on the
 * host instead (tests without a GPU), on copies of the operands made as for
 * the GPU.
 */

#ifndef _MPLAPACK_RBLAS_CUDA_COMMON_CUH_
#define _MPLAPACK_RBLAS_CUDA_COMMON_CUH_

#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include "Rgemm_cuda_ops_common.h"
#include "Rblas_lanes_common.h"
#include "Rblas_cuda_bridge_common.h"

#if defined(__CUDACC__)
#include <cuda_runtime.h>
#elif !defined(MPLAPACK_GEMM_CUDA_HOST_EMULATION)
#error "Rblas_cuda_common.cuh: compile with nvcc, or define MPLAPACK_GEMM_CUDA_HOST_EMULATION"
#endif

namespace mplapack_blas {

// Engine of Rblas_lanes_common.h on the element type of Ops.
template <class Ops> struct ops_engine {
    typedef typename Ops::F V;
    struct M {
        V *p;
        long ld;
    };
    static const bool has_div = true;
    MPLAPACK_BLAS_HD V reg(int) const { return Ops::from_d(0.0); }
    MPLAPACK_BLAS_HD V at(const M &a, long i, long j) const { return a.p[i + j * a.ld]; }
    MPLAPACK_BLAS_HD void put(M &a, long i, long j, const V &v) const { a.p[i + j * a.ld] = v; }
    MPLAPACK_BLAS_HD void set(V &r, const V &a) const { r = a; }
    MPLAPACK_BLAS_HD void zero(V &r) const { r = Ops::from_d(0.0); }
    MPLAPACK_BLAS_HD void one(V &r) const { r = Ops::from_d(1.0); }
    MPLAPACK_BLAS_HD void add(V &r, const V &a, const V &b) const { r = Ops::add(a, b); }
    MPLAPACK_BLAS_HD void sub(V &r, const V &a, const V &b) const { r = Ops::sub(a, b); }
    MPLAPACK_BLAS_HD void mul(V &r, const V &a, const V &b) const { r = Ops::mul(a, b); }
    MPLAPACK_BLAS_HD void div(V &r, const V &a, const V &b) const { r = Ops::div(a, b); }
    // a == 0.0 and a == 1.0 of the C++ types
    MPLAPACK_BLAS_HD bool iszero(const V &a) const { return Ops::eq(a, Ops::from_d(0.0)); }
    MPLAPACK_BLAS_HD bool isone(const V &a) const { return Ops::eq(a, Ops::from_d(1.0)); }
};

template <class Ops> struct gpu_lane_args {
    blas_shape s;
    typename Ops::F *al, *be, *a, *b, *c;
    long lda, ldb, ldc; // of the operands (Rtrmv: b is a copy of c)
};

template <class Ops> MPLAPACK_BLAS_HD void gpu_lane(const gpu_lane_args<Ops> &g, long lane)
{
    typedef ops_engine<Ops> E;
    E e;
    blas_mats<E> x;
    x.al.p = g.al, x.al.ld = 1;
    x.be.p = g.be, x.be.ld = 1;
    x.a.p = g.a, x.a.ld = g.lda;
    x.b.p = g.b, x.b.ld = g.ldb;
    x.c.p = g.c, x.c.ld = g.ldc;
    blas_lane(e, g.s, x, lane);
}

#if defined(__CUDACC__)

template <class Ops> __global__ void k_blas_lanes(gpu_lane_args<Ops> g, long lanes)
{
    const long lane = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (lane < lanes)
        gpu_lane<Ops>(g, lane);
}

// Device copy of the span of an operand.
template <class F> class gpu_span {
  public:
    gpu_span() : d_(NULL) {}
    ~gpu_span() { cudaFree(d_); }
    cudaError_t up(const gpu_operand<F> &o)
    {
        o_ = o;
        cudaError_t e = cudaMalloc(&d_, (o.count ? o.count : 1) * sizeof(F));
        if (e == cudaSuccess && o.count)
            e = cudaMemcpy(d_, o.host + o.lo, o.count * sizeof(F), cudaMemcpyHostToDevice);
        return e;
    }
    cudaError_t down() { return o_.download && o_.count ? cudaMemcpy(const_cast<F *>(o_.host) + o_.lo, d_, o_.count * sizeof(F), cudaMemcpyDeviceToHost) : cudaSuccess; }
    F *at0() const { return d_ - o_.lo; } // logical element (0, 0)

  private:
    F *d_;
    gpu_operand<F> o_;
    gpu_span(const gpu_span &);
    gpu_span &operator=(const gpu_span &);
};

template <class Ops> int blas_gpu(const blas_shape &s, const typename Ops::F &alpha, const typename Ops::F &beta, const gpu_operand<typename Ops::F> &a, const gpu_operand<typename Ops::F> &b, const gpu_operand<typename Ops::F> &c)
{
    typedef typename Ops::F F;
    gpu_operand<F> sc[2];
    sc[0] = gpu_scalar_operand(&alpha);
    sc[1] = gpu_scalar_operand(&beta);
    gpu_span<F> dal, dbe, da, db, dc;
    cudaError_t e = dal.up(sc[0]);
    if (e == cudaSuccess)
        e = dbe.up(sc[1]);
    if (e == cudaSuccess)
        e = da.up(a);
    if (e == cudaSuccess)
        e = db.up(b);
    if (e == cudaSuccess)
        e = dc.up(c);
    if (e != cudaSuccess)
        return (int)e;
    gpu_lane_args<Ops> g;
    g.s = s;
    g.al = dal.at0(), g.be = dbe.at0(), g.a = da.at0(), g.b = db.at0(), g.c = dc.at0();
    g.lda = a.ld, g.ldb = b.ld, g.ldc = c.ld;
    const long lanes = blas_lanes(s);
    const int threads = 64;
    if (lanes > 0) {
        k_blas_lanes<Ops><<<(unsigned)((lanes + threads - 1) / threads), threads>>>(g, lanes);
        e = cudaGetLastError();
        if (e == cudaSuccess)
            e = cudaDeviceSynchronize();
    }
    if (e == cudaSuccess)
        e = dc.down();
    if (e == cudaSuccess)
        e = da.down();
    return (int)e;
}

#else

// Host emulation: the same lanes on copies of the operand spans.
template <class F> struct emul_span {
    std::vector<F> v;
    gpu_operand<F> o;
    void up(const gpu_operand<F> &op)
    {
        o = op;
        v.assign((size_t)(op.count ? op.count : 1), F());
        if (op.count)
            for (long t = 0; t < op.count; t++)
                v[(size_t)t] = op.host[op.lo + t];
    }
    void down()
    {
        if (o.download)
            for (long t = 0; t < o.count; t++)
                const_cast<F *>(o.host)[o.lo + t] = v[(size_t)t];
    }
    F *at0() { return v.data() - o.lo; }
};

template <class Ops> int blas_gpu(const blas_shape &s, const typename Ops::F &alpha, const typename Ops::F &beta, const gpu_operand<typename Ops::F> &a, const gpu_operand<typename Ops::F> &b, const gpu_operand<typename Ops::F> &c)
{
    typedef typename Ops::F F;
    emul_span<F> dal, dbe, da, db, dc;
    dal.up(gpu_scalar_operand(&alpha));
    dbe.up(gpu_scalar_operand(&beta));
    da.up(a), db.up(b), dc.up(c);
    gpu_lane_args<Ops> g;
    g.s = s;
    g.al = dal.at0(), g.be = dbe.at0(), g.a = da.at0(), g.b = db.at0(), g.c = dc.at0();
    g.lda = a.ld, g.ldb = b.ld, g.ldc = c.ld;
    const long lanes = blas_lanes(s);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (long lane = 0; lane < lanes; lane++)
        gpu_lane<Ops>(g, lane);
    dc.down();
    da.down();
    return 0;
}

#endif

} // namespace mplapack_blas

#endif
