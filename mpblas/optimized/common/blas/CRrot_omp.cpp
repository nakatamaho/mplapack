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

// The reference CRrot (mpblas/reference/CRrot.cpp, unchanged) under the name
// CRrot_reference_cpu, and CRrot_omp: the reference code on independent parts of
// the output in parallel.

#define CRrot CRrot_reference_cpu
#include "../../../reference/CRrot.cpp"
#undef CRrot

#include <algorithm>
#include "../Rblas_split_common.h"

bool CRrot_omp(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy, REAL const c, REAL const s)
{
    using namespace mplapack_blas;
    if (n <= 0 || !(incx > 0 && incy > 0))
        return false;
    return blas_split(16.0 * n, n, [&](long p0, long p1) { CRrot_reference_cpu(p1 - p0, zx + p0 * incx, incx, zy + p0 * incy, incy, c, s); });
}
