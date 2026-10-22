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

// Rgbmv of the optimized libraries: large calls with unit increments run
// in parallel (Rgbmv_omp.cpp), with the reference result bit for bit; the
// others run the reference code.

#include <mpblas.h>

bool Rgbmv_omp(const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy);
void Rgbmv_reference_cpu(const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy);

void Rgbmv(const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy)
{
    if (Rgbmv_omp(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy))
        return;
    Rgbmv_reference_cpu(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy);
}
