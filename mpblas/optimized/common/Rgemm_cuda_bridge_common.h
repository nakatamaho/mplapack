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
 * Host side of the dd/qd CUDA Rgemm: decides whether a call runs on the GPU
 * and with which algorithm, from environment variables read at every call
 * (<P> is MPLAPACK_DD or MPLAPACK_QD):
 *
 *   <P>_CUDA=0                    never use the GPU
 *   <P>_CUDA_MIN_MNK=N            smallest m*n*k sent to the GPU (32768)
 *   <P>_CUDA_KERNEL=naive|tiled   conventional kernel (tiled)
 *   <P>_CUDA_WINOGRAD_CUTOFF=N    Winograd when min(m,n,k) > N (0: off)
 *   <P>_CUDA_OZAKI=1              Ozaki scheme with cuBLAS DGEMM (off)
 *   <P>_CUDA_OZAKI_SLICES=S       number of Ozaki slices
 *   <P>_CUDA_VERBOSE=1            report CUDA errors that cause a CPU fallback
 *
 * Calls that do not run on the GPU, or whose GPU run fails, return false and
 * are computed by the CPU Rgemm.
 */

#ifndef _MPLAPACK_RGEMM_CUDA_BRIDGE_COMMON_H_
#define _MPLAPACK_RGEMM_CUDA_BRIDGE_COMMON_H_

#include <climits>
#include <cstdio>
#include <cstring>
#include <string>
#include "Rgemm_cuda_ops_common.h"

namespace mplapack_gemm {

template <class T, class Ops> bool gemm_cuda_bridge(const char *prefix, const shape &s, const T &alpha, const T &beta, const T *A, const T *B, T *C) {
    typedef typename Ops::F F;
    static_assert(sizeof(T) == sizeof(F), "element layout");
#if !MPLAPACK_GEMM_QD_OPS_MATCH_LIBQD
    // the device arithmetic reproduces libqd's IEEE add / accurate mul only
    (void)prefix, (void)s, (void)alpha, (void)beta, (void)A, (void)B, (void)C;
    return false;
#else
    const std::string P(prefix);
    if (env_long((P + "_CUDA").c_str(), 1) == 0)
        return false;
    if ((double)s.m * s.n * s.k < (double)env_long((P + "_CUDA_MIN_MNK").c_str(), 32768))
        return false;
    if (s.m > INT_MAX || s.n > INT_MAX || s.k > INT_MAX)
        return false;
    if (gpu_device_count<Ops>() <= 0)
        return false;

    gpu_scalars<F> sc;
    std::memcpy(&sc.alpha, &alpha, sizeof(F));
    std::memcpy(&sc.beta, &beta, sizeof(F));
    sc.beta_zero = (beta == 0.0);
    sc.beta_one = (beta == 1.0);
    const F *dA = reinterpret_cast<const F *>(A);
    const F *dB = reinterpret_cast<const F *>(B);
    F *dC = reinterpret_cast<F *>(C);
    const bool verbose = env_long((P + "_CUDA_VERBOSE").c_str(), 0) > 0;

    int algo = GPU_TILED;
    const char *kernel = std::getenv((P + "_CUDA_KERNEL").c_str());
    if (kernel && std::strcmp(kernel, "naive") == 0)
        algo = GPU_NAIVE;
    const long cutoff = env_long((P + "_CUDA_WINOGRAD_CUTOFF").c_str(), 0);
    if (cutoff > 0 && s.m > cutoff && s.n > cutoff && s.k > cutoff)
        algo = GPU_WINOGRAD;
    ozaki_plan plan;
    if (s.k > 0 && env_long((P + "_CUDA_OZAKI").c_str(), 0) > 0 && ozaki_split(s, A, B, (int)env_long((P + "_CUDA_OZAKI_SLICES").c_str(), 0), plan))
        algo = GPU_OZAKI;

    const int rc = gemm_cuda<Ops>(s, algo, cutoff, &plan, sc, dA, dB, dC);
    if (rc != 0 && verbose)
        std::fprintf(stderr, "%s Rgemm: GPU error %d, computing on the CPU\n", prefix, rc);
    return rc == 0;
#endif
}

} // namespace mplapack_gemm

#endif
