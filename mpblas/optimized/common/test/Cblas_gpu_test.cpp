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

// Cblas_omp_test.cpp for the CUDA complex Level 1/2/3 BLAS (dd/td/qd and
// MPFR: ../blas_cuda): every call is sent to the GPU and compared bit for bit with
// the reference routines.  With
// MPLAPACK_GEMM_CUDA_HOST_EMULATION (dd/td/qd) or
// MPLAPACK_MPFR_CUDA_HOST_EMULATION (MPFR) the GPU code runs on the host;
// otherwise link the CUDA library, and the test is skipped (77) without a
// device.

#include <cstdlib>

#define RBLAS_TEST_GPU 1

#if defined(MPLAPACK_BUILD_WITH_MPFR)
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
#include "../../mpfr/cuda/test/stack_policy_test.h"
#endif
namespace mplapack_mpfr_cuda {
int device_count();
}
static void rblas_test_gpu_setup()
{
    setenv("MPLAPACK_MPFR_CUDA_MIN_MNK", "0", 1);
    setenv("MPLAPACK_MPFR_CUDA_L12", "1", 1);
    // signed zeros as the reference (the cu_freal engine: =0 in the environment)
    setenv("MPLAPACK_MPFR_CUDA_FORCE_RUNTIME", "1", 0);
}
static bool rblas_test_have_gpu() {
    if(mplapack_mpfr_cuda::device_count()<=0) return false;
#ifndef MPLAPACK_MPFR_CUDA_HOST_EMULATION
    long p=mpfrxx::default_precision_bits();
    size_t required=mpfr_stack_test::fixed(p)?size_t(16384)*((p+511)/512):65536;
    mpfr_stack_test::policies("Cblas Ctrsm",required,[] {
        mpfr_stack_test::result<mpc_class> r;
        r.a.assign(15,mpc_class(19)); r.b.assign(5,mpc_class(23));
        for(long j=0;j<3;++j) for(long i=0;i<3;++i) r.a[i+j*5]=i==j?mpc_class(4):mpc_class(.125,.0625);
        r.computed=Ctrsm_mpfr_cuda("L","U","N","N",3,1,mpc_class(1),r.a.data(),5,r.b.data(),5);
        return r;
    });
#endif
    return true;
}
#else
#if defined(MPLAPACK_BUILD_WITH_DD)
#define RBLAS_TEST_PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#define RBLAS_TEST_PREFIX "MPLAPACK_TD"
#elif defined(MPLAPACK_BUILD_WITH_QD)
#define RBLAS_TEST_PREFIX "MPLAPACK_QD"
#else
#error "Rblas_gpu_test.cpp: dd, td, qd or mpfr only"
#endif
static void rblas_test_gpu_setup()
{
    setenv(RBLAS_TEST_PREFIX "_BLAS_CUDA_MIN", "0", 1);
    setenv(RBLAS_TEST_PREFIX "_BLAS_CUDA_L12", "1", 1);
}
#ifdef MPLAPACK_GEMM_CUDA_HOST_EMULATION
static bool rblas_test_have_gpu() { return true; }
#else
#include <cuda_runtime.h>
static bool rblas_test_have_gpu()
{
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}
#endif
#endif

#include "Cblas_omp_test.cpp"

#ifdef MPLAPACK_GEMM_CUDA_HOST_EMULATION
// the device code, compiled for the host
#if defined(MPLAPACK_BUILD_WITH_DD)
#include "../../dd/cuda/Rgemm_gpu_cuda.cu"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#include "../../td/cuda/Rgemm_gpu_cuda.cu"
#else
#include "../../qd/cuda/Rgemm_gpu_cuda.cu"
#endif
#endif
