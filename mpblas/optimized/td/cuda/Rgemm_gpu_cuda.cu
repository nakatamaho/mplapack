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

// Device code of the td CUDA Rgemm (../../common/Rgemm_cuda_common.cuh) and
// Level 1/2/3 BLAS (../../common/Rblas_cuda_common.cuh, and for the complex
// routines ../../common/Cblas_cuda_common.cuh).
// Compile with -fmad=false.  With MPLAPACK_GEMM_CUDA_HOST_EMULATION and a
// C++ compiler (-x c++ -ffp-contract=off) the same code runs on the host.

#include "../../common/Rgemm_cuda_common.cuh"
#include "../../common/Rblas_cuda_common.cuh"
#include "../../common/Rblas_dev_types.h"
#include "../../common/Cblas_cuda_common.cuh"

namespace mplapack_gemm {
template int gemm_cuda<td_ops>(const shape &, int, long, const ozaki_plan *, const gpu_scalars<td_ops::F> &, const td_ops::F *, const td_ops::F *, td_ops::F *);
template int gpu_device_count<td_ops>();
} // namespace mplapack_gemm

namespace mplapack_blas {
template int blas_gpu<mplapack_gemm::td_ops>(const blas_shape &, const mplapack_gemm::td_ops::F &, const mplapack_gemm::td_ops::F &, const gpu_operand<mplapack_gemm::td_ops::F> &, const gpu_operand<mplapack_gemm::td_ops::F> &, const gpu_operand<mplapack_gemm::td_ops::F> &);
} // namespace mplapack_blas

namespace mplapack_blas_dev {
typedef mplapack_gemm::td_ops OPS;
template int cgpu_run<job_Cgemm<OPS> >(job_Cgemm<OPS>, const mplapack_blas::gpu_operand<job_Cgemm<OPS>::E> *, int);
template int cgpu_run<job_Csymm<OPS> >(job_Csymm<OPS>, const mplapack_blas::gpu_operand<job_Csymm<OPS>::E> *, int);
template int cgpu_run<job_Chemm<OPS> >(job_Chemm<OPS>, const mplapack_blas::gpu_operand<job_Chemm<OPS>::E> *, int);
template int cgpu_run<job_Ctrmm<OPS> >(job_Ctrmm<OPS>, const mplapack_blas::gpu_operand<job_Ctrmm<OPS>::E> *, int);
template int cgpu_run<job_Ctrsm<OPS> >(job_Ctrsm<OPS>, const mplapack_blas::gpu_operand<job_Ctrsm<OPS>::E> *, int);
template int cgpu_run<job_Cgemv<OPS> >(job_Cgemv<OPS>, const mplapack_blas::gpu_operand<job_Cgemv<OPS>::E> *, int);
template int cgpu_run<job_Cgeru<OPS> >(job_Cgeru<OPS>, const mplapack_blas::gpu_operand<job_Cgeru<OPS>::E> *, int);
template int cgpu_run<job_Cgerc<OPS> >(job_Cgerc<OPS>, const mplapack_blas::gpu_operand<job_Cgerc<OPS>::E> *, int);
template int cgpu_run<job_Caxpy<OPS> >(job_Caxpy<OPS>, const mplapack_blas::gpu_operand<job_Caxpy<OPS>::E> *, int);
template int cgpu_run<job_Ccopy<OPS> >(job_Ccopy<OPS>, const mplapack_blas::gpu_operand<job_Ccopy<OPS>::E> *, int);
template int cgpu_run<job_Cswap<OPS> >(job_Cswap<OPS>, const mplapack_blas::gpu_operand<job_Cswap<OPS>::E> *, int);
template int cgpu_run<job_CRrot<OPS> >(job_CRrot<OPS>, const mplapack_blas::gpu_operand<job_CRrot<OPS>::E> *, int);
template int cgpu_run<job_Cscal<OPS> >(job_Cscal<OPS>, const mplapack_blas::gpu_operand<job_Cscal<OPS>::E> *, int);
template int cgpu_run<job_CRscal<OPS> >(job_CRscal<OPS>, const mplapack_blas::gpu_operand<job_CRscal<OPS>::E> *, int);
template int cgpu_run<job_Csyrk<OPS> >(job_Csyrk<OPS>, const mplapack_blas::gpu_operand<job_Csyrk<OPS>::E> *, int);
template int cgpu_run<job_Cherk<OPS> >(job_Cherk<OPS>, const mplapack_blas::gpu_operand<job_Cherk<OPS>::E> *, int);
template int cgpu_run<job_Csyr2k<OPS> >(job_Csyr2k<OPS>, const mplapack_blas::gpu_operand<job_Csyr2k<OPS>::E> *, int);
template int cgpu_run<job_Cher2k<OPS> >(job_Cher2k<OPS>, const mplapack_blas::gpu_operand<job_Cher2k<OPS>::E> *, int);
template int cgpu_run<job_Cgemmtr<OPS> >(job_Cgemmtr<OPS>, const mplapack_blas::gpu_operand<job_Cgemmtr<OPS>::E> *, int);
template int cgpu_run<job_Rgemmtr<OPS> >(job_Rgemmtr<OPS>, const mplapack_blas::gpu_operand<job_Rgemmtr<OPS>::E> *, int);
template int cgpu_run<job_Cher<OPS> >(job_Cher<OPS>, const mplapack_blas::gpu_operand<job_Cher<OPS>::E> *, int);
template int cgpu_run<job_Cher2<OPS> >(job_Cher2<OPS>, const mplapack_blas::gpu_operand<job_Cher2<OPS>::E> *, int);
template int cgpu_run<job_Chemv<OPS> >(job_Chemv<OPS>, const mplapack_blas::gpu_operand<job_Chemv<OPS>::E> *, int);
template int cgpu_run<job_Chbmv<OPS> >(job_Chbmv<OPS>, const mplapack_blas::gpu_operand<job_Chbmv<OPS>::E> *, int);
template int cgpu_run<job_Rsbmv<OPS> >(job_Rsbmv<OPS>, const mplapack_blas::gpu_operand<job_Rsbmv<OPS>::E> *, int);
template int cgpu_run<job_Chpmv<OPS> >(job_Chpmv<OPS>, const mplapack_blas::gpu_operand<job_Chpmv<OPS>::E> *, int);
template int cgpu_run<job_Rspmv<OPS> >(job_Rspmv<OPS>, const mplapack_blas::gpu_operand<job_Rspmv<OPS>::E> *, int);
template int cgpu_run<job_Rgbmv<OPS> >(job_Rgbmv<OPS>, const mplapack_blas::gpu_operand<job_Rgbmv<OPS>::E> *, int);
template int cgpu_run<job_Cgbmv<OPS> >(job_Cgbmv<OPS>, const mplapack_blas::gpu_operand<job_Cgbmv<OPS>::E> *, int);
template int cgpu_run<job_Ctrmv<OPS> >(job_Ctrmv<OPS>, const mplapack_blas::gpu_operand<job_Ctrmv<OPS>::E> *, int);
template int cgpu_run<job_Rtbmv<OPS> >(job_Rtbmv<OPS>, const mplapack_blas::gpu_operand<job_Rtbmv<OPS>::E> *, int);
template int cgpu_run<job_Ctbmv<OPS> >(job_Ctbmv<OPS>, const mplapack_blas::gpu_operand<job_Ctbmv<OPS>::E> *, int);
template int cgpu_run<job_Rtpmv<OPS> >(job_Rtpmv<OPS>, const mplapack_blas::gpu_operand<job_Rtpmv<OPS>::E> *, int);
template int cgpu_run<job_Ctpmv<OPS> >(job_Ctpmv<OPS>, const mplapack_blas::gpu_operand<job_Ctpmv<OPS>::E> *, int);
template int cgpu_run<job_Rspr<OPS> >(job_Rspr<OPS>, const mplapack_blas::gpu_operand<job_Rspr<OPS>::E> *, int);
template int cgpu_run<job_Rspr2<OPS> >(job_Rspr2<OPS>, const mplapack_blas::gpu_operand<job_Rspr2<OPS>::E> *, int);
template int cgpu_run<job_Chpr<OPS> >(job_Chpr<OPS>, const mplapack_blas::gpu_operand<job_Chpr<OPS>::E> *, int);
template int cgpu_run<job_Chpr2<OPS> >(job_Chpr2<OPS>, const mplapack_blas::gpu_operand<job_Chpr2<OPS>::E> *, int);
} // namespace mplapack_blas_dev
