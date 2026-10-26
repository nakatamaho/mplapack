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

// Host side of the Level 1/2/3 BLAS of the dd/td/qd CUDA libraries: the
// calls of ../blas/R*_omp.cpp, run on the GPU (../Rblas_cuda_bridge_common.h).
// The entry points R*_<backend>_cuda return false when the call is left to
// the CPU.

#include <mpblas.h>
#include "../Rblas_cuda_entry_common.h"
#include "../Rblas_cuda_bridge_common.h"

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
#error "Rblas_bridge_cuda.cpp: dd, td or qd only"
#endif

using namespace mplapack_blas;

namespace {
bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
} // namespace

bool MPLAPACK_CUDA_ENTRY(Rsyrk)(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL const beta, REAL *cc, INTEGER const ldc)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_syrk<REAL>(c, ls, false, uplo, trans, n, k, alpha, a, lda, (REAL *)0, 1, beta, cc, ldc))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rsyr2k)(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const beta, REAL *cc, INTEGER const ldc)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_syrk<REAL>(c, ls, true, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, cc, ldc))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rsymm)(const char *side, const char *uplo, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const beta, REAL *cc, INTEGER const ldc)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_symm<REAL>(c, ls, side, uplo, m, n, alpha, a, lda, b, ldb, beta, cc, ldc))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rtrmm)(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_trxm<REAL>(c, ls, false, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rtrsm)(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_trxm<REAL>(c, ls, true, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rgemv)(const char *trans, INTEGER const m, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_gemv<REAL>(c, ls, trans, m, n, alpha, a, lda, x, incx, beta, y, incy))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rsymv)(const char *uplo, INTEGER const n, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_symv<REAL>(c, ls, uplo, n, alpha, a, lda, x, incx, beta, y, incy))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rtrmv)(const char *uplo, const char *trans, const char *diag, INTEGER const n, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_trmv<REAL>(c, ls, uplo, trans, diag, n, a, lda, x, incx))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rger)(INTEGER const m, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *y, INTEGER const incy, REAL *a, INTEGER const lda)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_rank<REAL>(c, ls, 0, "U", m, n, alpha, x, incx, y, incy, a, lda))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rsyr)(const char *uplo, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *a, INTEGER const lda)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_rank<REAL>(c, ls, 1, uplo, n, n, alpha, x, incx, (REAL *)0, 1, a, lda))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rsyr2)(const char *uplo, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *y, INTEGER const incy, REAL *a, INTEGER const lda)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_rank<REAL>(c, ls, 2, uplo, n, n, alpha, x, incx, y, incy, a, lda))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rscal)(INTEGER const n, REAL const da, REAL *dx, INTEGER const incx)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_level1<REAL>(c, L1_SCAL, n, da, REAL(0.0), dx, incx, (REAL *)0, 1))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rswap)(INTEGER const n, REAL *dx, INTEGER const incx, REAL *dy, INTEGER const incy)
{
    using namespace mplapack_blas;
    blas_call<REAL> c;
    if (!make_level1<REAL>(c, L1_SWAP, n, REAL(0.0), REAL(0.0), dx, incx, dy, incy))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, c);
}

bool MPLAPACK_CUDA_ENTRY(Rrot)(INTEGER const n, REAL *dx, INTEGER const incx, REAL *dy, INTEGER const incy, REAL const c_, REAL const s)
{
    using namespace mplapack_blas;
    blas_call<REAL> call;
    if (!make_level1<REAL>(call, L1_ROT, n, c_, s, dx, incx, dy, incy))
        return false;
    return blas_gpu_bridge<REAL, OPS>(PREFIX, call);
}
