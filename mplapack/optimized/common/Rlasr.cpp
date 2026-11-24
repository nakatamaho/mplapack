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

// Rlasr of the optimized libraries: try freal before the reference routine.

#include <mpblas.h>
#include <mplapack.h>
#if defined(MPLAPACK_BUILD_WITH_MPFR) && !defined(MPLAPACK_FREAL_INSTANCE)
#include "../../../mpblas/optimized/mpfr/openmp/Rfreal_blas.h"
#endif

#define Rlasr Rlasr_reference_cpu
#include "../../reference/Rlasr.cpp"
#undef Rlasr

void Rlasr(const char *side, const char *pivot, const char *direct, INTEGER const m, INTEGER const n, REAL *c, REAL *s, REAL *a, INTEGER const lda) {
#if defined(MPLAPACK_BUILD_WITH_MPFR) && !defined(MPLAPACK_FREAL_INSTANCE)
    if (Rlasr_freal(side, pivot, direct, m, n, c, s, a, lda))
        return;
#endif
    Rlasr_reference_cpu(side, pivot, direct, m, n, c, s, a, lda);
}
