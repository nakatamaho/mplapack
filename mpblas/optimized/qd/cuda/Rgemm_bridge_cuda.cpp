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

// Host side of the qd CUDA Rgemm (../../common/Rgemm_cuda_bridge_common.h).

#include <mpblas_qd.h>
#include "../../common/Rgemm_cuda_bridge_common.h"

namespace mplapack_gemm {
template <> struct ozaki_traits<qd_real> : qd_like_ozaki_traits<qd_real, 212, 4> {};
} // namespace mplapack_gemm

// Runs C := alpha*op(A)*op(B) + beta*C on the GPU; false when the call is
// left to the CPU.  alpha != 0 and m, n > 0 are assumed.
bool Rgemm_gpu(bool nota, bool notb, mplapackint m, mplapackint n, mplapackint k, qd_real alpha, qd_real *A, mplapackint lda, qd_real *B, mplapackint ldb, qd_real beta, qd_real *C, mplapackint ldc) {
    const mplapack_gemm::shape s = {nota, notb, (long)m, (long)n, (long)k, (long)lda, (long)ldb, (long)ldc};
    return mplapack_gemm::gemm_cuda_bridge<qd_real, mplapack_gemm::qd_ops>("MPLAPACK_QD", s, alpha, beta, A, B, C);
}
