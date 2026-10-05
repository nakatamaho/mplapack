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

// Blocked SIMD/OpenMP Rgemm for dd_real, plus the opt-in Winograd and Ozaki
// variants (../../common/Rgemm_cpu_common.h).  Called by ../Rgemm.cpp.

#include <mpblas_dd.h>
#include "../../common/Rgemm_cpu_common.h"
#include "../../common/Rgemm_simd_kernels_common.h"

namespace mplapack_gemm {
template <> struct ozaki_traits<dd_real> : qd_like_ozaki_traits<dd_real, 106, 2> {};
} // namespace mplapack_gemm

bool Rgemm_blocked_omp(bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, dd_real alpha, dd_real *A, mplapackint lda, dd_real *B, mplapackint ldb, dd_real beta, dd_real *C, mplapackint ldc) {
    const mplapack_gemm::shape s = {nota, notb, (long)m, (long)n, (long)k, (long)lda, (long)ldb, (long)ldc};
#if MPLAPACK_GEMM_QD_OPS_MATCH_LIBQD
    return mplapack_gemm::gemm_cpu<mplapack_gemm::dd_simd_kernel<dd_real> >("MPLAPACK_DD", s, alpha, beta, A, B, C);
#else
    return mplapack_gemm::gemm_cpu<mplapack_gemm::scalar_kernel<dd_real, 64, 64, 256> >("MPLAPACK_DD", s, alpha, beta, A, B, C);
#endif
}
