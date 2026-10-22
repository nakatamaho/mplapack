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

// Runtime precision (mpfr_rt<NL>): the precision of the call is set before
// the launch (the variables are defined in Cblas_mpfr_device_cuda.cu).
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

} // namespace mplapack_blas_dev

#define MPLAPACK_MPFR_CBLAS_FIXED_JOB(R, PB) template int cgpu_run<job_##R<mpfr_fixed<PB> > >(job_##R<mpfr_fixed<PB> >, const mplapack_blas::gpu_operand<job_##R<mpfr_fixed<PB> >::E> *, int);
#define MPLAPACK_MPFR_CBLAS_RT_JOB(R, NL) template int cgpu_run_rt<job_##R<mpfr_rt<NL> > >(job_##R<mpfr_rt<NL> >, const mplapack_blas::gpu_operand<job_##R<mpfr_rt<NL> >::E> *, int, long);

#endif
