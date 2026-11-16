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
 * On x86-64 ELF targets, functions marked MPLAPACK_GEMM_CLONES are compiled
 * for AVX2+FMA (x86-64-v3) and for the baseline, and the variant is selected
 * at load time (GNU ifunc).  AVX-512 is deliberately not targeted.  Define
 * MPLAPACK_GEMM_NO_TARGET_CLONES to compile only for the target given on the
 * command line.  Loops marked MPLAPACK_GEMM_SIMD are vectorized.
 */

#ifndef _MPLAPACK_RGEMM_SIMD_COMMON_H_
#define _MPLAPACK_RGEMM_SIMD_COMMON_H_

#if defined(__x86_64__) && defined(__ELF__) && defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(__NVCC__) && !defined(MPLAPACK_GEMM_NO_TARGET_CLONES)
#define MPLAPACK_GEMM_CLONES __attribute__((target_clones("arch=x86-64-v3", "default")))
#else
#define MPLAPACK_GEMM_CLONES
#endif

#if defined(_OPENMP)
#define MPLAPACK_GEMM_SIMD _Pragma("omp simd")
#elif defined(__clang__)
#define MPLAPACK_GEMM_SIMD _Pragma("clang loop vectorize(enable)")
#elif defined(__GNUC__)
#define MPLAPACK_GEMM_SIMD _Pragma("GCC ivdep")
#else
#define MPLAPACK_GEMM_SIMD
#endif

#endif
