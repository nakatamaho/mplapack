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

// Rscal of libmplapack_mpfr_opt_cuda: with MPLAPACK_MPFR_CUDA_L12=1 eligible
// calls run on the GPU (blas3_bridge_cuda.cpp), the others on the CPU
// (../../common/blas/Rscal_omp.cpp: lane-parallel OpenMP, else the reference
// code).

#include <mpblas_mpfr_cuda.h>

bool Rscal_omp(mplapackint const n, mpfr_class const da, mpfr_class *dx, mplapackint const incx);
void Rscal_reference_cpu(mplapackint const n, mpfr_class const da, mpfr_class *dx, mplapackint const incx);

void Rscal(mplapackint const n, mpfr_class const da, mpfr_class *dx, mplapackint const incx)
{
    if (Rscal_mpfr_cuda(n, da, dx, incx))
        return;
    if (Rscal_omp(n, da, dx, incx))
        return;
    Rscal_reference_cpu(n, da, dx, incx);
}
