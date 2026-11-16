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

// Host side of the dd/td/qd complex Level 1/2/3 BLAS (and Rgemmtr) on CUDA
// (../Cblas_cuda_common.cuh): the entry points <R>_<backend>_cuda return
// false when the call is left to the CPU.  Environment (<P> = MPLAPACK_DD, MPLAPACK_TD or MPLAPACK_QD), as
// for the real routines (../Rblas_cuda_bridge_common.h): <P>_BLAS_CUDA=0,
// <P>_BLAS_CUDA_MIN=N, <P>_BLAS_CUDA_L12=1 (Level 1 and 2 on the GPU),
// <P>_CUDA_VERBOSE=1.  Calls with alpha == 0 or quick returns run on the GPU
// too: every thread runs the reference code.

#include <mpblas.h>
#include <cstring>
#include <string>
#include "../Rblas_cuda_bridge_common.h"
#include "../Rblas_dev_types.h"
#include "Cblas_bridge_body.h"

#if defined(MPLAPACK_BUILD_WITH_DD)
typedef mplapack_gemm::dd_ops OPS;
#define PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
typedef mplapack_gemm::td_ops OPS;
#define PREFIX "MPLAPACK_TD"
#elif defined(MPLAPACK_BUILD_WITH_QD)
typedef mplapack_gemm::qd_ops OPS;
#define PREFIX "MPLAPACK_QD"
#else
#error "Cblas_bridge_cuda.cpp: dd, td or qd only"
#endif

using mplapack_blas::gpu_operand;

namespace {

typedef mplapack_blas_dev::cgpu_types<OPS>::DR DR;
typedef mplapack_blas_dev::cgpu_types<OPS>::DC DC;
static_assert(sizeof(DR) == sizeof(REAL) && sizeof(DC) == sizeof(COMPLEX), "element layout");

// The element-type side of Cblas_bridge_body.h: the device elements have
// the layout of the host ones, so operands are copied as they are.
struct copy_ctx {
    bool eligible(bool level3, double work) const
    {
        const std::string P(PREFIX);
        if (mplapack_gemm::env_long((P + "_BLAS_CUDA").c_str(), 1) == 0)
            return false;
        if (!level3 && mplapack_gemm::env_long((P + "_BLAS_CUDA_L12").c_str(), 0) == 0)
            return false;
        if (work < (double)mplapack_gemm::env_long((P + "_BLAS_CUDA_MIN").c_str(), 32768))
            return false;
        return mplapack_gemm::gpu_device_count<OPS>() > 0;
    }
    void sc(DC &d, const COMPLEX &s) const { std::memcpy(&d, &s, sizeof d); }
    void sc(DR &d, const REAL &s) const { std::memcpy(&d, &s, sizeof d); }
    // rows x cols elements of x with leading dimension (vectors: increment) ld > 0
    template <class E, class T> gpu_operand<E> op(T *x, long ld, long rows, long cols, bool download) const
    {
        gpu_operand<E> o = {reinterpret_cast<const E *>(x), ld, 0, (cols - 1) * ld + rows, download, rows, cols, mplapack_blas::REGION_ALL};
        return o;
    }
    template <class J> bool run(const J &j, const gpu_operand<typename J::E> *ops, int nops) const
    {
        const int rc = mplapack_blas_dev::cgpu_run(j, ops, nops);
        if (rc != 0 && mplapack_gemm::env_long(PREFIX "_CUDA_VERBOSE", 0) > 0)
            std::fprintf(stderr, "%s complex BLAS: GPU error %d, computing on the CPU\n", PREFIX, rc);
        return rc == 0;
    }
};

} // namespace

#define MPLAPACK_CBLAS_GPU_RUN(F, ...) ([&]() { copy_ctx ctx; return mplapack_cblas_gpu::F<OPS>(ctx, __VA_ARGS__); }())
MPLAPACK_CBLAS_GPU_FUNCTIONS(MPLAPACK_CBLAS_GPU_RUN)
