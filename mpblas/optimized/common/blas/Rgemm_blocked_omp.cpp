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

// Rgemm_blocked_omp of every optimized library (called by <type>/Rgemm.cpp):
// the small-block engine (../Rblas_block_common.h) with the operations of
// <type>/openmp/Rgemm_{NN,NT,TN,TT}_omp.cpp, so the result is theirs bit for
// bit; for dd, td, qd and binary128 also the opt-in Winograd and Ozaki
// variants (../Rgemm_cpu_common.h).  Returns false when the caller is to run
// its own code (small calls, MPLAPACK_BLAS_BLOCK=0).

#include <mpblas.h>
#include "../Rblas_cpu_common.h"

#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD) || defined(MPLAPACK_BUILD_WITH_BINARY128)
#include "../Rgemm_cpu_common.h"

namespace mplapack_gemm {
#if defined(MPLAPACK_BUILD_WITH_DD)
template <> struct ozaki_traits<dd_real> : qd_like_ozaki_traits<dd_real, 106, 2> {};
#define MPLAPACK_GEMM_PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
template <> struct ozaki_traits<td_real> : qd_like_ozaki_traits<td_real, 159, 3> {};
#define MPLAPACK_GEMM_PREFIX "MPLAPACK_TD"
#elif defined(MPLAPACK_BUILD_WITH_QD)
template <> struct ozaki_traits<qd_real> : qd_like_ozaki_traits<qd_real, 212, 4> {};
#define MPLAPACK_GEMM_PREFIX "MPLAPACK_QD"
#else
template <> struct ozaki_traits<REAL> {
    static const int precision = 113;
    static const int ncomp = 3;
    static double to_double(const REAL &x) { return (double)x; }
    static void add_double(REAL &acc, double d) { acc += (REAL)d; }
};
#define MPLAPACK_GEMM_PREFIX "MPLAPACK_BINARY128"
#endif
} // namespace mplapack_gemm
#endif

bool Rgemm_blocked_omp(bool nota, bool notb, INTEGER m, INTEGER n, INTEGER k, REAL alpha, REAL *A, INTEGER lda, REAL *B, INTEGER ldb, REAL beta, REAL *C, INTEGER ldc)
{
#ifdef MPLAPACK_GEMM_PREFIX
    const mplapack_gemm::shape s = {nota, notb, (long)m, (long)n, (long)k, (long)lda, (long)ldb, (long)ldc};
    return mplapack_gemm::gemm_cpu<REAL>(MPLAPACK_GEMM_PREFIX, s, alpha, beta, A, B, C);
#else
    return mplapack_blas::blk_gemm<mplapack_blas::blas_thread_context>(nota, notb, (long)m, (long)n, (long)k, alpha, A, (long)lda, B, (long)ldb, beta, C, (long)ldc);
#endif
}
