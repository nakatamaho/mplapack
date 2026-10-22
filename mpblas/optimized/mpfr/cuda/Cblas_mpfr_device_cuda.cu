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
 * Device code of the MPFR complex Level 1/2/3 BLAS (and Rgemmtr): the jobs
 * of ../../common/Cblas_cuda_jobs_common.h with the cu_freal device types of
 * Cblas_mpfr_dev_types.h, for the precisions of
 * MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS, and with the runtime-precision types
 * (MPFR/MPC functions) for the limb counts of MPLAPACK_MPFR_CUDA_RT_TIERS.
 * Compile with -fmad=false.  With
 * MPLAPACK_MPFR_CUDA_HOST_EMULATION the jobs run on the host.
 */

#include "Rgemm_kernel_cuda.h"
#include "Cblas_mpfr_dev_types.h"
#include "../../common/Cblas_cuda_common.cuh"

namespace mplapack_blas_dev {

#define MPLAPACK_MPFR_CBLAS_JOB(R, PB) template int cgpu_run<job_##R<mpfr_fixed<PB> > >(job_##R<mpfr_fixed<PB> >, const mplapack_blas::gpu_operand<job_##R<mpfr_fixed<PB> >::E> *, int);
#define MPLAPACK_MPFR_CBLAS_INSTANTIATE(PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cgemm, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Csymm, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Chemm, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Ctrmm, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Ctrsm, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cgemv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cgeru, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cgerc, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Caxpy, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Ccopy, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cswap, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(CRrot, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cscal, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(CRscal, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Csyrk, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cherk, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Csyr2k, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cher2k, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cgemmtr, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Rgemmtr, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cher, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cher2, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Chemv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Chbmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Rsbmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Chpmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Rspmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Rgbmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Cgbmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Ctrmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Rtbmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Ctbmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Rtpmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Ctpmv, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Rspr, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Rspr2, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Chpr, PB) \
    MPLAPACK_MPFR_CBLAS_JOB(Chpr2, PB)
MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(MPLAPACK_MPFR_CBLAS_INSTANTIATE)
#undef MPLAPACK_MPFR_CBLAS_INSTANTIATE

// Runtime precision (mpfr_rt<NL>): the precision of the call is set before
// the launch.
long mplapack_mpfr_rt_prec_h = 0;
#if MPLAPACK_CBLAS_CUDA_DEVICE
__device__ long mplapack_mpfr_rt_prec_d = 0;
#endif

template <class Job> int cgpu_run_rt(Job job, const mplapack_blas::gpu_operand<typename Job::E> *ops, int nops, long prec)
{
    mplapack_mpfr_rt_prec_h = prec;
#if MPLAPACK_CBLAS_CUDA_DEVICE
    const cudaError_t e = cudaMemcpyToSymbol(mplapack_mpfr_rt_prec_d, &prec, sizeof prec);
    if (e != cudaSuccess)
        return (int)e;
#endif
    return cgpu_run(job, ops, nops);
}

#define MPLAPACK_MPFR_CBLAS_RT_JOB(R, NL) template int cgpu_run_rt<job_##R<mpfr_rt<NL> > >(job_##R<mpfr_rt<NL> >, const mplapack_blas::gpu_operand<job_##R<mpfr_rt<NL> >::E> *, int, long);
#define MPLAPACK_MPFR_CBLAS_RT_INSTANTIATE(NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cgemm, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Csymm, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Chemm, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Ctrmm, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Ctrsm, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cgemv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cgeru, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cgerc, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Caxpy, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Ccopy, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cswap, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(CRrot, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cscal, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(CRscal, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Csyrk, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cherk, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Csyr2k, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cher2k, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cgemmtr, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Rgemmtr, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cher, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cher2, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Chemv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Chbmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Rsbmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Chpmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Rspmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Rgbmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Cgbmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Ctrmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Rtbmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Ctbmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Rtpmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Ctpmv, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Rspr, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Rspr2, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Chpr, NL) \
    MPLAPACK_MPFR_CBLAS_RT_JOB(Chpr2, NL)
MPLAPACK_MPFR_CUDA_RT_TIERS(MPLAPACK_MPFR_CBLAS_RT_INSTANTIATE)
#undef MPLAPACK_MPFR_CBLAS_RT_INSTANTIATE
#undef MPLAPACK_MPFR_CBLAS_RT_JOB
#undef MPLAPACK_MPFR_CBLAS_JOB

} // namespace mplapack_blas_dev
