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
 * Shared by the translation units of the MPFR complex BLAS device code
 * (Cblas_mpfr_device*_cuda.cu, split so that they compile in parallel): the
 * jobs of ../../common/Cblas_cuda_jobs_common.h and cgpu_run_rt.
 */

#ifndef _MPLAPACK_CBLAS_MPFR_DEVICE_CUDA_H_
#define _MPLAPACK_CBLAS_MPFR_DEVICE_CUDA_H_

#include "Rgemm_kernel_cuda.h"
#include "Cblas_mpfr_dev_types.h"
#include "../../common/Cblas_cuda_common.cuh"
#if MPLAPACK_CBLAS_CUDA_DEVICE
#include <cstdlib>
#include "mpc_arena_multi.cuh"
#endif

// X(R, TAG) for every job
#define MPLAPACK_MPFR_CBLAS_JOBS(X, TAG) \
    X(Cgemm, TAG) \
    X(Csymm, TAG) \
    X(Chemm, TAG) \
    X(Ctrmm, TAG) \
    X(Ctrsm, TAG) \
    X(Cgemv, TAG) \
    X(Cgeru, TAG) \
    X(Cgerc, TAG) \
    X(Caxpy, TAG) \
    X(Ccopy, TAG) \
    X(Cswap, TAG) \
    X(CRrot, TAG) \
    X(Cscal, TAG) \
    X(CRscal, TAG) \
    X(Csyrk, TAG) \
    X(Cherk, TAG) \
    X(Csyr2k, TAG) \
    X(Cher2k, TAG) \
    X(Cgemmtr, TAG) \
    X(Rgemmtr, TAG) \
    X(Cher, TAG) \
    X(Cher2, TAG) \
    X(Chemv, TAG) \
    X(Chbmv, TAG) \
    X(Rsbmv, TAG) \
    X(Chpmv, TAG) \
    X(Rspmv, TAG) \
    X(Rgbmv, TAG) \
    X(Cgbmv, TAG) \
    X(Ctrmv, TAG) \
    X(Rtbmv, TAG) \
    X(Ctbmv, TAG) \
    X(Rtpmv, TAG) \
    X(Ctpmv, TAG) \
    X(Rspr, TAG) \
    X(Rspr2, TAG) \
    X(Chpr, TAG) \
    X(Chpr2, TAG)

namespace mplapack_blas_dev {

// Fixed-precision complex BLAS uses value-semantic components; the same team
// schedule is valid with its MPC division and the existing scoped stack guard.
template<int PB> struct cgpu_is_ctrsm<job_Ctrsm<mpfr_fixed<PB>>> { static constexpr bool value=true; };
template<int PB> struct cgpu_is_herk<job_Cherk<mpfr_fixed<PB>>> { static constexpr bool value=true; };
template<int PB> struct herk_constants<job_Cherk<mpfr_fixed<PB>>> {
    static bool zero_alpha(const job_Cherk<mpfr_fixed<PB>> &j) { return j.alpha.v.is_zero(); }
};

// Host launch predicates inspect the value representation directly. Avoid
// instantiating noinline arithmetic helpers (and their weak exports) on host.
template<int PB> struct ctrsm_constants<job_Ctrsm<mpfr_fixed<PB>>> {
    using V=mcomplex<PB>;
    MPLAPACK_DEV_HD static bool zero(const V &v) { return v.re_.v.is_zero() && v.im_.v.is_zero(); }
    MPLAPACK_DEV_HD static bool one(const V &v) {
        const auto &r=v.re_.v;
        if(!v.im_.v.is_zero() || r.sign!=1 || r.exp!=1 || r.m[r.N-1]!=0x8000000000000000ULL) return false;
        for(int i=0;i<r.N-1;++i) if(r.m[i]) return false;
        return true;
    }
};

// Cgemm on several devices (cgpu_run: the columns of C cut over them).  The
// runtime precisions share one arena (mpfr_rt_arena_session) and their
// precision is a per-device variable, which cgpu_run_rt sets on every device
// the call may use.
} // namespace mplapack_blas_dev
template <int PB> struct mplapack_cuda_multi::prefix_of<mplapack_blas_dev::mpfr_fixed<PB> > {
    static constexpr const char *value = "MPLAPACK_MPFR";
};
template <int NL> struct mplapack_cuda_multi::prefix_of<mplapack_blas_dev::mpfr_rt<NL> > {
    static constexpr const char *value = "MPLAPACK_MPFR";
};
namespace mplapack_blas_dev {

#if MPLAPACK_CBLAS_CUDA_DEVICE
// The runtime precisions on a pool of threads with libmpc_cuda's arena
// (cgpu_arena), as Rgemm_rt_cuda.cu: MPLAPACK_MPFR_CUDA_RT_BLOCKS (256) x
// MPLAPACK_MPFR_CUDA_RT_THREADS (32) threads per device, 16 KB per 1024 bits
// each.  Every MPFR/MPC call returns its slab to where it was
// (Cblas_mpfr_dev_types.h), so a slab holds the scratch of one call; one
// that does not fit takes it from the device heap as without an arena.
struct mpfr_rt_arena_session {
    int error = 0;
    long pool = 0, threads = 32;
    char *saved_base = mpc_cuda_arena_base;
    size_t saved_slab = mpc_cuda_arena_slab;
    size_t *saved_top = mpc_cuda_arena_top;
    mplapack_mpfr_cuda::mpc_arena_multi multi;
    void *base1 = nullptr, *top1 = nullptr;
    bool installed = false;
    static long env(const char *name, long fallback) {
        const char *v = std::getenv(name);
        char *end = nullptr;
        const long x = v ? std::strtol(v, &end, 10) : 0;
        return v && end != v && !*end && x > 0 ? x : fallback;
    }
    explicit mpfr_rt_arena_session(int G) {
        if (G <= 0)
            return;
        threads = env("MPLAPACK_MPFR_CUDA_RT_THREADS", 32);
        pool = env("MPLAPACK_MPFR_CUDA_RT_BLOCKS", 256) * threads;
        const size_t slab = (size_t)16384 * (size_t)((mplapack_mpfr_rt_prec_h + 1023) / 1024);
        char *base;
        size_t *top;
        if (G == 1) {
            if ((error = (int)cudaMalloc(&base1, (size_t)pool * slab)) || (error = (int)cudaMalloc(&top1, (size_t)pool * sizeof(size_t))) || (error = (int)cudaMemset(top1, 0, (size_t)pool * sizeof(size_t))))
                return;
            base = (char *)base1;
            top = (size_t *)top1;
        } else {
            if ((error = multi.allocate(G, (size_t)pool, slab)))
                return;
            base = multi.base;
            top = multi.top;
        }
        mpc_cuda_arena_base = base;
        mpc_cuda_arena_slab = slab;
        mpc_cuda_arena_top = top;
        installed = true;
    }
    ~mpfr_rt_arena_session() {
        if (installed) {
            mpc_cuda_arena_base = saved_base;
            mpc_cuda_arena_slab = saved_slab;
            mpc_cuda_arena_top = saved_top;
        }
        cudaFree(base1);
        cudaFree(top1);
    }
};
template <int NL> struct cgpu_arena<rcomplex<NL> > {
    static constexpr bool pooled = true;
    typedef mpfr_rt_arena_session session;
};
template <int NL> struct cgpu_arena<rreal<NL> > {
    static constexpr bool pooled = true;
    typedef mpfr_rt_arena_session session;
};
#endif

// Runtime precision (mpfr_rt<NL>): the precision of the call is set before
// the launch (the variables are defined in Cblas_mpfr_device_cuda.cu).
template <class Job> int cgpu_run_rt(Job job, const mplapack_blas::gpu_operand<typename Job::E> *ops, int nops, long prec)
{
    mplapack_mpfr_rt_prec_h = prec;
#if MPLAPACK_CBLAS_CUDA_DEVICE
    cudaError_t e = cudaMemcpyToSymbol(mplapack_mpfr_rt_prec_d, &prec, sizeof prec);
    if (e != cudaSuccess)
        return (int)e;
    if constexpr (cgpu_cgemm_ops<Job>::value) {
        using namespace mplapack_cuda_multi;
        const int G = devices(env_gpus("MPLAPACK_MPFR"));
        if (G > 1) {
            device_restore restore;
            for (int d = 0; d < G; ++d)
                if ((e = cudaSetDevice(d)) != cudaSuccess || (e = cudaMemcpyToSymbol(mplapack_mpfr_rt_prec_d, &prec, sizeof prec)) != cudaSuccess)
                    return (int)e;
        }
    }
#endif
    return cgpu_run(job, ops, nops);
}

} // namespace mplapack_blas_dev

#define MPLAPACK_MPFR_CBLAS_FIXED_JOB(R, PB) template int cgpu_run<job_##R<mpfr_fixed<PB> > >(job_##R<mpfr_fixed<PB> >, const mplapack_blas::gpu_operand<job_##R<mpfr_fixed<PB> >::E> *, int);
#define MPLAPACK_MPFR_CBLAS_RT_JOB(R, NL) template int cgpu_run_rt<job_##R<mpfr_rt<NL> > >(job_##R<mpfr_rt<NL> >, const mplapack_blas::gpu_operand<job_##R<mpfr_rt<NL> >::E> *, int, long);

#endif
