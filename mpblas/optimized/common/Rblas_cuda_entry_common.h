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
 * MPLAPACK_CUDA_ENTRY(R) names the CUDA entry point R_<backend>_cuda of the
 * library a source shared by several backends is compiled for (AGENTS.md,
 * "Accelerator entry points"), declared in include/mpblas_<backend>_cuda.h.
 */

#ifndef _MPLAPACK_RBLAS_CUDA_ENTRY_COMMON_H_
#define _MPLAPACK_RBLAS_CUDA_ENTRY_COMMON_H_

#if defined(MPLAPACK_BUILD_WITH_DD)
#include <mpblas_dd_cuda.h>
#define MPLAPACK_CUDA_ENTRY(r) r##_dd_cuda
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include <mpblas_td_cuda.h>
#define MPLAPACK_CUDA_ENTRY(r) r##_td_cuda
#elif defined(MPLAPACK_BUILD_WITH_QD)
#include <mpblas_qd_cuda.h>
#define MPLAPACK_CUDA_ENTRY(r) r##_qd_cuda
#elif defined(MPLAPACK_BUILD_WITH_MPFR)
#include <mpblas_mpfr_cuda.h>
#define MPLAPACK_CUDA_ENTRY(r) r##_mpfr_cuda
#else
#error "Rblas_cuda_entry_common.h: dd, td, qd or mpfr only"
#endif

#endif
