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
 * Element types, arithmetic and entry points of the dd/td/qd CUDA Rgemm
 * (Rgemm_cuda_common.cuh), usable from host C++ code that does not include
 * the CUDA headers.
 */

#ifndef _MPLAPACK_RGEMM_CUDA_OPS_COMMON_H_
#define _MPLAPACK_RGEMM_CUDA_OPS_COMMON_H_

#include "Rgemm_qd_ops_common.h"
#include "Rgemm_blocked_common.h"
#include "Rgemm_ozaki_common.h"

namespace mplapack_gemm {

enum { GPU_NAIVE = 0, GPU_TILED = 1, GPU_WINOGRAD = 2, GPU_OZAKI = 3 };

// Arithmetic of the element type F, as libQD3's operators with QD_BF:
// a + b = bf_add(a, b), a - b = bf_add(a, -b), a * b = bf_mul(a, b).
struct dd_ops {
    typedef dd2 F;
    static MPLAPACK_GEMM_HD F add(const F &a, const F &b) { return dd_bf_add(a, b); }
    static MPLAPACK_GEMM_HD F mul(const F &a, const F &b) { return dd_bf_mul(a, b); }
    static MPLAPACK_GEMM_HD F sub(const F &a, const F &b) {
        F nb;
        nb.x[0] = -b.x[0];
        nb.x[1] = -b.x[1];
        return dd_bf_add(a, nb);
    }
    static MPLAPACK_GEMM_HD F add_d(const F &a, double b) { return dd_add_d(a, b); }
    static MPLAPACK_GEMM_HD F div(const F &a, const F &b) { return dd_div(a, b); }
    static MPLAPACK_GEMM_HD F from_d(double d) {
        F r;
        r.x[0] = d;
        r.x[1] = 0.0;
        return r;
    }
    // a == b: qd::dd_comparison_eq
    static MPLAPACK_GEMM_HD bool eq(const F &a, const F &b) {
        if (qd_ops_isnan(a.x[0]) || qd_ops_isnan(a.x[1]) || qd_ops_isnan(b.x[0]) || qd_ops_isnan(b.x[1]))
            return false;
        if (MPLAPACK_GEMM_ISINF(a.x[0]) || MPLAPACK_GEMM_ISINF(b.x[0]))
            return a.x[0] == b.x[0];
        return sub(a, b).x[0] == 0.0;
    }
};

struct qd_ops {
    typedef qd4 F;
    static MPLAPACK_GEMM_HD F add(const F &a, const F &b) { return qd_bf_add(a, b); }
    static MPLAPACK_GEMM_HD F mul(const F &a, const F &b) { return qd_bf_mul(a, b); }
    static MPLAPACK_GEMM_HD F sub(const F &a, const F &b) {
        F nb;
        for (int q = 0; q < 4; q++)
            nb.x[q] = -b.x[q];
        return qd_bf_add(a, nb);
    }
    static MPLAPACK_GEMM_HD F add_d(const F &a, double b) { return qd_add_d(a, b); }
    static MPLAPACK_GEMM_HD F div(const F &a, const F &b) { return qd_div(a, b); }
    static MPLAPACK_GEMM_HD F from_d(double d) {
        F r;
        r.x[0] = d;
        r.x[1] = r.x[2] = r.x[3] = 0.0;
        return r;
    }
    // a == b: qd::comparison_eq
    static MPLAPACK_GEMM_HD bool eq(const F &a, const F &b) {
        for (int q = 0; q < 4; q++)
            if (qd_ops_isnan(a.x[q]) || qd_ops_isnan(b.x[q]))
                return false;
        if (MPLAPACK_GEMM_ISINF(a.x[0]) || MPLAPACK_GEMM_ISINF(b.x[0]))
            return a.x[0] == b.x[0];
        return sub(a, b).x[0] == 0.0;
    }
};

struct td_ops {
    typedef td3 F;
    static MPLAPACK_GEMM_HD F add(const F &a, const F &b) { return td_bf_add(a, b); }
    static MPLAPACK_GEMM_HD F mul(const F &a, const F &b) { return td_bf_mul(a, b); }
    static MPLAPACK_GEMM_HD F sub(const F &a, const F &b) { return td_bf_add(a, td_make(-b.x[0], -b.x[1], -b.x[2])); }
    static MPLAPACK_GEMM_HD F add_d(const F &a, double b) { return td_add_d(a, b); }
    static MPLAPACK_GEMM_HD F div(const F &a, const F &b) { return td_div(a, b); }
    static MPLAPACK_GEMM_HD F from_d(double d) { return td_make(d, 0.0, 0.0); }
    // a == b: td::comparison_eq
    static MPLAPACK_GEMM_HD bool eq(const F &a, const F &b) {
        for (int q = 0; q < 3; q++)
            if (qd_ops_isnan(a.x[q]) || qd_ops_isnan(b.x[q]))
                return false;
        if (MPLAPACK_GEMM_ISINF(a.x[0]) || MPLAPACK_GEMM_ISINF(b.x[0]))
            return a.x[0] == b.x[0];
        const F d = sub(a, b);
        return d.x[0] == 0.0 && d.x[1] == 0.0 && d.x[2] == 0.0;
    }
};

// Scalars of a call; beta_zero / beta_one are computed on the host with the
// comparison operators of the C++ type.
template <class F> struct gpu_scalars {
    F alpha, beta;
    bool beta_zero, beta_one;
};

// C := alpha*op(A)*op(B) + beta*C on the GPU.  A, B, C are host arrays laid
// out as in Rgemm; cutoff is the Winograd cutoff (GPU_WINOGRAD) and plan the
// Ozaki splitting (GPU_OZAKI).  Returns 0 on success, otherwise a
// cudaError_t / cublasStatus_t code, C unchanged.  Instantiated in the
// libraries' Rgemm_gpu_cuda.cu.
template <class Ops> int gemm_cuda(const shape &s, int algo, long cutoff, const ozaki_plan *plan, const gpu_scalars<typename Ops::F> &sc, const typename Ops::F *A, const typename Ops::F *B, typename Ops::F *C);

// Number of CUDA devices (0 when the runtime reports an error).
template <class Ops> int gpu_device_count();

} // namespace mplapack_gemm

#endif
