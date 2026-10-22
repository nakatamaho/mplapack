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

// Rsyr2 of libmplapack_mpfr_opt_cuda: with MPLAPACK_MPFR_CUDA_L12=1 eligible
// calls run on the GPU (blas3_bridge_cuda.cpp), the others on the CPU
// (../../common/blas/Rsyr2_omp.cpp: lane-parallel OpenMP, else the reference
// code).

#include <mpblas_mpfr.h>

bool Rsyr2_mpfr_cuda(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *a, mplapackint const lda);
bool Rsyr2_omp(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *a, mplapackint const lda);
void Rsyr2_reference_cpu(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *a, mplapackint const lda);

void Rsyr2(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *a, mplapackint const lda)
{
    if (Rsyr2_mpfr_cuda(uplo, n, alpha, x, incx, y, incy, a, lda))
        return;
    if (Rsyr2_omp(uplo, n, alpha, x, incx, y, incy, a, lda))
        return;
    Rsyr2_reference_cpu(uplo, n, alpha, x, incx, y, incy, a, lda);
}
