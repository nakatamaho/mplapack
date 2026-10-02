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
 * Runtime-precision MPFR Rgemm on CUDA, ported from mpc_cuda's
 * demos/matmul_mpfr.cu: the arithmetic is cu_mpfr (MPFR 4.2.2 compiled for
 * the device) and every operation is correctly rounded exactly like the CPU
 * MPFR, so the result is identical to the CPU Rgemm, signed zeros included.
 *
 * Operands are passed as packed arrays: for element e, limbs[e*nlimbs ..]
 * hold the MPFR significand, exp[e] the MPFR exponent and kind[e] the signed
 * MPFR custom kind (+-MPFR_ZERO_KIND or +-MPFR_REGULAR_KIND).
 *
 * This header must not include MPFR or gmpfrxx_mkII.
 */

#ifndef _MPLAPACK_MPFR_RGEMM_RT_CUDA_H_
#define _MPLAPACK_MPFR_RGEMM_RT_CUDA_H_

#include <cstddef>
#include "Rgemm_kernel_cuda.h"

namespace mplapack_mpfr_cuda {

struct rt_pack {
    unsigned long long *limbs;
    long *exp;
    signed char *kind;
};

// Values of the MPFR custom kinds, shared by MPFR and cu_mpfr.
const int RT_ZERO_KIND = 2;
const int RT_REGULAR_KIND = 3;

// Largest |exponent| the device MPFR accepts with its default range.
const long RT_DEVICE_EMAX = (1L << 30) - 1;

// C := alpha*op(A)*op(B) + beta*C with every value of precision prec.
// alpha and beta are packs of one element; A, B and C are packed column-major
// as described for gemm_shape.  Returns 0 on success; otherwise nothing was
// computed and the caller falls back to the CPU.
int gemm_rt(const gemm_shape &s, long prec, const rt_pack &alpha, const rt_pack &beta, const rt_pack &A, size_t countA, const rt_pack &B, size_t countB, rt_pack &C);

} // namespace mplapack_mpfr_cuda

#endif
