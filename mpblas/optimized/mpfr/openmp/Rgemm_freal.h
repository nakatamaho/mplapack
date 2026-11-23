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

// The fixed-precision MPFR Rgemm on the CPU (Rgemm_freal_*.cpp).

#ifndef MPLAPACK_RGEMM_FREAL_H
#define MPLAPACK_RGEMM_FREAL_H

#include <mpblas_mpfr.h>

// The precisions (bits) with a fixed-precision kernel.
#define MPFREAL_PRECISIONS(X) X(128) X(192) X(256) X(320) X(384) X(448) X(512) X(576) X(640) X(704) X(768) X(832) X(896) X(960) X(1024) X(1536) X(2048)

// C := alpha op(A) op(B) + beta C with cu_freal<prec> (Rgemm_freal_engine.h),
// bit for bit the result of openmp/Rgemm_{NN,NT,TN,TT}_omp.cpp.  Returns
// false (C untouched) when the call is left to them: MPLAPACK_MPFR_CPU_FREAL=0,
// a precision without a kernel, a rounding mode other than MPFR_RNDN, an
// operand of another precision, NaN or Inf, exponents that could leave
// [emin, emax], or a call below MPLAPACK_MPFR_CPU_FREAL_MIN (m n k, 512).
// The caller has done the quick returns (alpha == 0, k == 0).
// MPLAPACK_MPFR_CPU_FREAL=2 uses the portable kernel even where the x86-64
// BMI2/ADX one could run.
bool Rgemm_freal_cpu(bool nota, bool notb, long m, long n, long k, const mpfr_class &alpha, const mpfr_class *A, long lda, const mpfr_class *B, long ldb,
                     const mpfr_class &beta, mpfr_class *C, long ldc);

// The kernels (beta_kind: 0 beta == 0, 1 beta == 1, 2 otherwise); false
// for a precision without a kernel.
bool Rgemm_freal_x86(int prec, bool nota, bool notb, long m, long n, long k, const mpfr_class &alpha, const mpfr_class *A, long lda, const mpfr_class *B,
                     long ldb, const mpfr_class &beta, int beta_kind, mpfr_class *C, long ldc, bool parallel);
bool Rgemm_freal_generic(int prec, bool nota, bool notb, long m, long n, long k, const mpfr_class &alpha, const mpfr_class *A, long lda, const mpfr_class *B,
                         long ldb, const mpfr_class &beta, int beta_kind, mpfr_class *C, long ldc, bool parallel);
bool Rgemm_freal_x86_available();

#endif
