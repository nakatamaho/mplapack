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
 * MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS.  Compile with -fmad=false.  With
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
#undef MPLAPACK_MPFR_CBLAS_JOB

} // namespace mplapack_blas_dev
