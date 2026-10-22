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
 * Parallel runs of a reference routine on independent parts of its output:
 * the columns of C in Cgemm, of B in Ctrmm/Ctrsm with side L, the rows with
 * side R, and so on.  Every element of the output is computed by the
 * reference code with the same operations in the same order as in one call,
 * so the result is the reference result bit for bit, for every type.
 *
 * blas_split(work, parts, run) calls run(p0, p1) on disjoint ranges covering
 * [0, parts) from OpenMP threads that carry the caller's MPFR/GMP precision
 * context (Rblas_cpu_common.h).  It returns false (the caller then makes one
 * reference call) when OpenMP is off (MPLAPACK_BLAS_OMP=0), there is one
 * thread, parts < 2 or work < MPLAPACK_BLAS_OMP_MIN.
 */

#ifndef _MPLAPACK_RBLAS_SPLIT_COMMON_H_
#define _MPLAPACK_RBLAS_SPLIT_COMMON_H_

#include <algorithm>
#include "Rblas_cpu_common.h"

namespace mplapack_blas {

template <class F> bool blas_split(double work, long parts, F run)
{
#ifndef _OPENMP
    (void)work, (void)parts, (void)run;
    return false;
#else
    if (blas_env_long("MPLAPACK_BLAS_OMP", 1) == 0 || omp_get_max_threads() < 2 || parts < 2)
        return false;
    if (work < (double)blas_env_long("MPLAPACK_BLAS_OMP_MIN", 32768))
        return false;
    const long chunks = std::min(parts, 4L * omp_get_max_threads());
    const blas_thread_context ctx;
#pragma omp parallel for schedule(dynamic, 1)
    for (long q = 0; q < chunks; q++) {
        const blas_thread_context::scope sc(ctx);
        const long p0 = parts * q / chunks, p1 = parts * (q + 1) / chunks;
        if (p1 > p0)
            run(p0, p1);
    }
    return true;
#endif
}

} // namespace mplapack_blas

#endif
