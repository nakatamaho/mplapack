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
 * Cblas_mpfr_dev_types.h for the precisions of
 * MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS; the runtime-precision types are in
 * Cblas_mpfr_device_rt*_cuda.cu.  Compile with -fmad=false.  With
 * MPLAPACK_MPFR_CUDA_HOST_EMULATION the jobs run on the host.
 */

#include "Cblas_mpfr_device_cuda.h"

namespace mplapack_blas_dev {

int mpfr_cblas_stack(size_t bytes)
{
#if MPLAPACK_CBLAS_CUDA_DEVICE
    size_t cur = 0;
    cudaError_t e = cudaDeviceGetLimit(&cur, cudaLimitStackSize);
    if (e == cudaSuccess && cur < bytes)
        e = cudaDeviceSetLimit(cudaLimitStackSize, bytes);
    return (int)e;
#else
    (void)bytes;
    return 0;
#endif
}

long mplapack_mpfr_rt_prec_h = 0;
#if MPLAPACK_CBLAS_CUDA_DEVICE
__device__ long mplapack_mpfr_rt_prec_d = 0;
#endif

#define MPLAPACK_MPFR_CBLAS_INSTANTIATE(PB) MPLAPACK_MPFR_CBLAS_JOBS(MPLAPACK_MPFR_CBLAS_FIXED_JOB, PB)
MPLAPACK_MPFR_CUDA_FIXED_PRECISIONS(MPLAPACK_MPFR_CBLAS_INSTANTIATE)
#undef MPLAPACK_MPFR_CBLAS_INSTANTIATE

} // namespace mplapack_blas_dev
