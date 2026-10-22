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
 * Host side of the dd/td/qd Level 1/2/3 BLAS on CUDA (Rblas_cuda_common.cuh):
 * decides whether a call runs on the GPU, from environment variables read at
 * every call (<P> is MPLAPACK_DD, MPLAPACK_TD or MPLAPACK_QD):
 *
 *   <P>_BLAS_CUDA=0        never use the GPU for these routines
 *   <P>_BLAS_CUDA_MIN=N    smallest number of multiply-adds sent to the GPU
 *                          (32768)
 *   <P>_BLAS_CUDA_L12=1    also run Level 1 and 2 on the GPU (off: they are
 *                          bound by the transfers; meant for debugging)
 *   <P>_CUDA_VERBOSE=1     report CUDA errors that cause a CPU fallback
 *
 * Each operand is copied to the GPU as the span of memory it occupies; the
 * written operands are copied back.  Calls that do not run on the GPU, or
 * whose GPU run fails, return false and are computed on the CPU.
 */

#ifndef _MPLAPACK_RBLAS_CUDA_BRIDGE_COMMON_H_
#define _MPLAPACK_RBLAS_CUDA_BRIDGE_COMMON_H_

#include <cstdio>
#include <cstring>
#include <string>
#include "Rgemm_cuda_ops_common.h"
#include "Rblas_calls_common.h"
#include "Rblas_gpu_operand_common.h"

namespace mplapack_blas {

template <class F, class T> inline gpu_operand<F> gpu_operand_of(const blas_operand<T> &op, bool download)
{
    gpu_operand<F> o = {reinterpret_cast<const F *>(op.X), op.ld, 0, 0, download};
    if (op.X == 0 || op.rows <= 0 || op.cols <= 0)
        return o;
    const long last = (op.cols - 1) * op.ld; // offset of the last column
    o.lo = last < 0 ? last : 0;
    o.count = (last < 0 ? -last : last) + op.rows;
    return o;
}

// Defined in Rblas_cuda_common.cuh, instantiated in <type>/cuda/Rgemm_gpu_cuda.cu.
template <class Ops> int blas_gpu(const blas_shape &s, const typename Ops::F &alpha, const typename Ops::F &beta, const gpu_operand<typename Ops::F> &a, const gpu_operand<typename Ops::F> &b, const gpu_operand<typename Ops::F> &c);

template <class T, class Ops> bool blas_gpu_bridge(const char *prefix, const blas_call<T> &c)
{
    typedef typename Ops::F F;
    static_assert(sizeof(T) == sizeof(F), "element layout");
    const std::string P(prefix);
    if (mplapack_gemm::env_long((P + "_BLAS_CUDA").c_str(), 1) == 0)
        return false;
    if (c.s.op > L3_TRSM && mplapack_gemm::env_long((P + "_BLAS_CUDA_L12").c_str(), 0) == 0)
        return false;
    if (c.work < (double)mplapack_gemm::env_long((P + "_BLAS_CUDA_MIN").c_str(), 32768))
        return false;
    if (mplapack_gemm::gpu_device_count<Ops>() <= 0)
        return false;
    F alpha, beta;
    std::memcpy(&alpha, &c.alpha, sizeof(F));
    std::memcpy(&beta, &c.beta, sizeof(F));
    // Rtrmv: b is the original x, uploaded separately
    const gpu_operand<F> a = gpu_operand_of<F>(c.a, c.a_written);
    const gpu_operand<F> b = gpu_operand_of<F>(c.b_is_copy_of_c ? c.c : c.b, false);
    const gpu_operand<F> r = gpu_operand_of<F>(c.c, true);
    const int rc = blas_gpu<Ops>(c.s, alpha, beta, a, b, r);
    if (rc != 0 && mplapack_gemm::env_long((P + "_CUDA_VERBOSE").c_str(), 0) > 0)
        std::fprintf(stderr, "%s BLAS: GPU error %d, computing on the CPU\n", prefix, rc);
    return rc == 0;
}

} // namespace mplapack_blas

#endif
