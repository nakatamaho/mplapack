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

// Rspmv of the CUDA libraries: eligible calls (Level 2: only with
// <P>_BLAS_CUDA_L12=1 / MPLAPACK_MPFR_CUDA_L12=1) run on the GPU
// (Cblas_bridge_body.h, ../blas_dev/l2_rows_dev.h), the others on the CPU
// (../blas/Rspmv_omp.cpp, then the reference code).  All of them give the
// reference result bit for bit.

#include <mpblas.h>

bool Rspmv_gpu(const char *uplo, INTEGER const n, REAL const alpha, REAL *ap, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy);
bool Rspmv_omp(const char *uplo, INTEGER const n, REAL const alpha, REAL *ap, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy);
void Rspmv_reference_cpu(const char *uplo, INTEGER const n, REAL const alpha, REAL *ap, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy);

void Rspmv(const char *uplo, INTEGER const n, REAL const alpha, REAL *ap, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy)
{
    if (Rspmv_gpu(uplo, n, alpha, ap, x, incx, beta, y, incy))
        return;
    if (Rspmv_omp(uplo, n, alpha, ap, x, incx, beta, y, incy))
        return;
    Rspmv_reference_cpu(uplo, n, alpha, ap, x, incx, beta, y, incy);
}
