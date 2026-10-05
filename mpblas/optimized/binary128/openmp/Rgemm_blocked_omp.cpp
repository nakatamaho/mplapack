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

// Blocked OpenMP Rgemm for binary128, plus the opt-in Winograd and Ozaki
// variants (../../common/Rgemm_cpu_common.h).  Called by ../Rgemm.cpp.

#include <mpblas_binary128.h>
#include "../../common/Rgemm_cpu_common.h"

namespace mplapack_gemm {
template <> struct ozaki_traits<mplapack_binary128_t> {
    static const int precision = 113;
    static const int ncomp = 3;
    static double to_double(const mplapack_binary128_t &x) { return (double)x; }
    static void add_double(mplapack_binary128_t &acc, double d) { acc += (mplapack_binary128_t)d; }
};
} // namespace mplapack_gemm

bool Rgemm_blocked_omp(bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, mplapack_binary128_t alpha, mplapack_binary128_t *A, mplapackint lda, mplapack_binary128_t *B, mplapackint ldb, mplapack_binary128_t beta, mplapack_binary128_t *C, mplapackint ldc) {
    const mplapack_gemm::shape s = {nota, notb, (long)m, (long)n, (long)k, (long)lda, (long)ldb, (long)ldc};
    return mplapack_gemm::gemm_cpu<mplapack_gemm::scalar_kernel<mplapack_binary128_t, 32, 32, 128> >("MPLAPACK_BINARY128", s, alpha, beta, A, B, C);
}
