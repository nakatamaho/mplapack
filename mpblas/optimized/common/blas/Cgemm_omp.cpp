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

// The reference Cgemm (mpblas/reference/Cgemm.cpp, unchanged) under the name
// Cgemm_reference_cpu, and Cgemm_omp: the reference code on independent parts of
// the output in parallel.

#define Cgemm Cgemm_reference_cpu
#include "../../../reference/Cgemm.cpp"
#undef Cgemm

#include <algorithm>
#include "../Rblas_split_common.h"
#include "../Rblas_cpu_common.h"
#if defined(MPLAPACK_BUILD_WITH_DD) || defined(MPLAPACK_BUILD_WITH_TD) || defined(MPLAPACK_BUILD_WITH_QD)
#include "../Cblas_simd_common.h"
#endif

#define MPLAPACK_GEMM_NO_OZAKI
#include "../Rgemm_cpu_common.h"
#if defined(MPLAPACK_BUILD_WITH_DOUBLE)
#define MPLAPACK_CGEMM_PREFIX "MPLAPACK_DOUBLE"
#elif defined(MPLAPACK_BUILD_WITH_BINARY80)
#define MPLAPACK_CGEMM_PREFIX "MPLAPACK_BINARY80"
#elif defined(MPLAPACK_BUILD_WITH_BINARY128)
#define MPLAPACK_CGEMM_PREFIX "MPLAPACK_BINARY128"
#elif defined(MPLAPACK_BUILD_WITH_DD)
#define MPLAPACK_CGEMM_PREFIX "MPLAPACK_DD"
#elif defined(MPLAPACK_BUILD_WITH_TD)
#define MPLAPACK_CGEMM_PREFIX "MPLAPACK_TD"
#elif defined(MPLAPACK_BUILD_WITH_QD)
#define MPLAPACK_CGEMM_PREFIX "MPLAPACK_QD"
#elif defined(MPLAPACK_BUILD_WITH_MPFR)
#define MPLAPACK_CGEMM_PREFIX "MPLAPACK_MPFR"
#elif defined(MPLAPACK_BUILD_WITH_GMP)
#define MPLAPACK_CGEMM_PREFIX "MPLAPACK_GMP"
#endif

bool Cgemm_omp(const char *transa, const char *transb, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc);

namespace {
bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }

// The leaves of the opt-in Winograd variant: C = A*B by Cgemm_omp (the
// leaves are at most the cutoff in one dimension, so it does not recurse).
struct cgemm_leaf {
    void product(long m, long n, long k, COMPLEX *C, long ldc, const COMPLEX *A, long lda, const COMPLEX *B, long ldb)
    {
        const COMPLEX one(1.0, 0.0), zero(0.0, 0.0);
        if ((double)m * n * k < 8192.0) { // small leaves serially (a fork costs more)
            Cgemm_reference_cpu("N", "N", (INTEGER)m, (INTEGER)n, (INTEGER)k, one, const_cast<COMPLEX *>(A), (INTEGER)lda, const_cast<COMPLEX *>(B), (INTEGER)ldb, zero, C, (INTEGER)ldc);
            return;
        }
        if (!Cgemm_omp("N", "N", (INTEGER)m, (INTEGER)n, (INTEGER)k, one, const_cast<COMPLEX *>(A), (INTEGER)lda, const_cast<COMPLEX *>(B), (INTEGER)ldb, zero, C, (INTEGER)ldc))
            Cgemm_reference_cpu("N", "N", (INTEGER)m, (INTEGER)n, (INTEGER)k, one, const_cast<COMPLEX *>(A), (INTEGER)lda, const_cast<COMPLEX *>(B), (INTEGER)ldb, zero, C, (INTEGER)ldc);
    }
};
} // namespace

bool Cgemm_omp(const char *transa, const char *transb, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    using namespace mplapack_blas;
    const bool nota = ls(transa, "N"), notb = ls(transb, "N");
    if ((!nota && !ls(transa, "T") && !ls(transa, "C")) || (!notb && !ls(transb, "T") && !ls(transb, "C")) || m <= 0 || n <= 0 || k < 0)
        return false;
    if (lda < std::max((INTEGER)1, nota ? m : k) || ldb < std::max((INTEGER)1, notb ? k : n) || ldc < std::max((INTEGER)1, m))
        return false;
    // the opt-in Winograd variant, <P>_CGEMM_WINOGRAD_CUTOFF=N (not the
    // reference result: its error is bounded normwise)
    if (!(alpha == COMPLEX(0.0, 0.0))) {
        const mplapack_gemm::shape s = {nota, notb, (long)m, (long)n, (long)k, (long)lda, (long)ldb, (long)ldc};
        cgemm_leaf leaf;
        if (mplapack_gemm::gemm_winograd_env<COMPLEX>(MPLAPACK_CGEMM_PREFIX, "_CGEMM_WINOGRAD_CUTOFF", s, alpha, beta, a, b, c, leaf, COMPLEX(0.0, 0.0), COMPLEX(1.0, 0.0), ls(transa, "C"), ls(transb, "C")))
            return true;
    }
    // the small-block engine where it is faster (../Rblas_block_common.h)
    if (!(alpha == COMPLEX(0.0, 0.0)) && blk_cgemm<blas_thread_context>(nota, ls(transa, "C"), notb, ls(transb, "C"), m, n, k, alpha, a, lda, b, ldb, beta, c, ldc))
        return true;
#ifdef MPLAPACK_CBLAS_SIMD_AVAILABLE
    if (nota ? cgemm_simd<COMPLEX>(notb, ls(transb, "C"), m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) : cgemm_t_simd<COMPLEX>(ls(transa, "C"), notb, ls(transb, "C"), m, n, k, alpha, a, lda, b, ldb, beta, c, ldc))
        return true;
#endif
    const double work = 4.0 * m * n * k;
    if (n >= m) // columns of C
        return blas_split(work, n, [&](long p0, long p1) { Cgemm_reference_cpu(transa, transb, m, p1 - p0, k, alpha, a, lda, notb ? b + p0 * ldb : b + p0, ldb, beta, c + p0 * ldc, ldc); });
    return blas_split(work, m, [&](long p0, long p1) { Cgemm_reference_cpu(transa, transb, p1 - p0, n, k, alpha, nota ? a + p0 : a + p0 * lda, lda, b, ldb, beta, c + p0, ldc); });
}
