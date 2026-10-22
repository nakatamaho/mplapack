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
// (../Cblas_cuda_common.cuh): R*_gpu return false when the call is left to
// the CPU.  Environment (<P> = MPLAPACK_DD, MPLAPACK_TD or MPLAPACK_QD), as
// for the real routines (../Rblas_cuda_bridge_common.h): <P>_BLAS_CUDA=0,
// <P>_BLAS_CUDA_MIN=N, <P>_BLAS_CUDA_L12=1 (Level 1 and 2 on the GPU),
// <P>_CUDA_VERBOSE=1.  Calls with alpha == 0 or quick returns run on the GPU
// too: every thread runs the reference code.

#include <mpblas.h>
#include <cstring>
#include <string>
#include <algorithm>
#include "../Cblas_cuda_jobs_common.h"

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

using namespace mplapack_blas_dev;
using mplapack_blas::gpu_operand;

namespace {

typedef cgpu_types<OPS>::DR DR;
typedef cgpu_types<OPS>::DC DC;
static_assert(sizeof(DR) == sizeof(REAL) && sizeof(DC) == sizeof(COMPLEX), "element layout");

bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }

bool eligible(bool level3, double work)
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

bool finish(int rc)
{
    if (rc != 0 && mplapack_gemm::env_long(PREFIX "_CUDA_VERBOSE", 0) > 0)
        std::fprintf(stderr, "%s complex BLAS: GPU error %d, computing on the CPU\n", PREFIX, rc);
    return rc == 0;
}

void ch(char *d, const char *s) { d[0] = s[0], d[1] = '\0'; }
void sc(DC &d, const COMPLEX &s) { std::memcpy(&d, &s, sizeof d); }
void sc(DR &d, const REAL &s) { std::memcpy(&d, &s, sizeof d); }

// rows x cols elements of x with leading dimension (vectors: increment) ld > 0
template <class E, class T> gpu_operand<E> op(T *x, long ld, long rows, long cols, bool download)
{
    gpu_operand<E> o = {reinterpret_cast<const E *>(x), ld, 0, (cols - 1) * ld + rows, download};
    return o;
}

} // namespace

bool Cgemm_gpu(const char *transa, const char *transb, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool nota = ls(transa, "N"), notb = ls(transb, "N");
    if ((!nota && !ls(transa, "T") && !ls(transa, "C")) || (!notb && !ls(transb, "T") && !ls(transb, "C")) || m <= 0 || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, nota ? m : k) || ldb < std::max((INTEGER)1, notb ? k : n) || ldc < std::max((INTEGER)1, m))
        return false;
    if (!eligible(true, 4.0 * m * n * k))
        return false;
    typedef job_Cgemm<OPS> J;
    typedef J::E E;
    J j;
    ch(j.ta, transa); ch(j.tb, transb); j.m = m, j.n = n, j.k = k, j.lda = lda, j.ldb = ldb, j.ldc = ldc; sc(j.alpha, alpha); sc(j.beta, beta);
    const gpu_operand<E> ops[] = {op<E>(a, lda, nota ? m : k, nota ? k : m, false), op<E>(b, ldb, notb ? k : n, notb ? n : k, false), op<E>(c, ldc, m, n, true)};
    return finish(cgpu_run(j, ops, 3));
}

bool Csymm_gpu(const char *side, const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool left = ls(side, "L");
    if ((!left && !ls(side, "R")) || (!ls(uplo, "U") && !ls(uplo, "L")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, left ? m : n) || ldb < std::max((INTEGER)1, m) || ldc < std::max((INTEGER)1, m))
        return false;
    if (!eligible(true, 4.0 * m * n * (left ? m : n)))
        return false;
    typedef job_Csymm<OPS> J;
    typedef J::E E;
    J j;
    ch(j.side, side); ch(j.uplo, uplo); j.m = m, j.n = n, j.lda = lda, j.ldb = ldb, j.ldc = ldc; sc(j.alpha, alpha); sc(j.beta, beta);
    const gpu_operand<E> ops[] = {op<E>(a, lda, left ? m : n, left ? m : n, false), op<E>(b, ldb, m, n, false), op<E>(c, ldc, m, n, true)};
    return finish(cgpu_run(j, ops, 3));
}

bool Chemm_gpu(const char *side, const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool left = ls(side, "L");
    if ((!left && !ls(side, "R")) || (!ls(uplo, "U") && !ls(uplo, "L")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, left ? m : n) || ldb < std::max((INTEGER)1, m) || ldc < std::max((INTEGER)1, m))
        return false;
    if (!eligible(true, 4.0 * m * n * (left ? m : n)))
        return false;
    typedef job_Chemm<OPS> J;
    typedef J::E E;
    J j;
    ch(j.side, side); ch(j.uplo, uplo); j.m = m, j.n = n, j.lda = lda, j.ldb = ldb, j.ldc = ldc; sc(j.alpha, alpha); sc(j.beta, beta);
    const gpu_operand<E> ops[] = {op<E>(a, lda, left ? m : n, left ? m : n, false), op<E>(b, ldb, m, n, false), op<E>(c, ldc, m, n, true)};
    return finish(cgpu_run(j, ops, 3));
}

bool Ctrmm_gpu(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb)
{
    const bool left = ls(side, "L");
    if ((!left && !ls(side, "R")) || (!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(transa, "N") && !ls(transa, "T") && !ls(transa, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, left ? m : n) || ldb < std::max((INTEGER)1, m))
        return false;
    if (!eligible(true, 4.0 * m * n * (left ? m : n)))
        return false;
    typedef job_Ctrmm<OPS> J;
    typedef J::E E;
    J j;
    ch(j.side, side); ch(j.uplo, uplo); ch(j.ta, transa); ch(j.diag, diag); j.m = m, j.n = n, j.lda = lda, j.ldb = ldb; sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {op<E>(a, lda, left ? m : n, left ? m : n, false), op<E>(b, ldb, m, n, true)};
    return finish(cgpu_run(j, ops, 2));
}

bool Ctrsm_gpu(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb)
{
    const bool left = ls(side, "L");
    if ((!left && !ls(side, "R")) || (!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(transa, "N") && !ls(transa, "T") && !ls(transa, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, left ? m : n) || ldb < std::max((INTEGER)1, m))
        return false;
    if (!eligible(true, 4.0 * m * n * (left ? m : n)))
        return false;
    typedef job_Ctrsm<OPS> J;
    typedef J::E E;
    J j;
    ch(j.side, side); ch(j.uplo, uplo); ch(j.ta, transa); ch(j.diag, diag); j.m = m, j.n = n, j.lda = lda, j.ldb = ldb; sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {op<E>(a, lda, left ? m : n, left ? m : n, false), op<E>(b, ldb, m, n, true)};
    return finish(cgpu_run(j, ops, 2));
}

bool Cgemv_gpu(const char *trans, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy)
{
    const bool notr = ls(trans, "N");
    if ((!notr && !ls(trans, "T") && !ls(trans, "C")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, m) || incx <= 0 || incy <= 0)
        return false;
    if (!eligible(false, 4.0 * m * n))
        return false;
    typedef job_Cgemv<OPS> J;
    typedef J::E E;
    J j;
    ch(j.tr, trans); j.m = m, j.n = n, j.lda = lda, j.incx = incx, j.incy = incy; sc(j.alpha, alpha); sc(j.beta, beta);
    const gpu_operand<E> ops[] = {op<E>(a, lda, m, n, false), op<E>(x, incx, 1, notr ? n : m, false), op<E>(y, incy, 1, notr ? m : n, true)};
    return finish(cgpu_run(j, ops, 3));
}

bool Cgeru_gpu(INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda)
{
    
    if (m <= 0 || n <= 0 || incx <= 0 || incy <= 0 || lda < std::max((INTEGER)1, m))
        return false;
    if (!eligible(false, 4.0 * m * n))
        return false;
    typedef job_Cgeru<OPS> J;
    typedef J::E E;
    J j;
    j.m = m, j.n = n, j.incx = incx, j.incy = incy, j.lda = lda; sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {op<E>(x, incx, 1, m, false), op<E>(y, incy, 1, n, false), op<E>(a, lda, m, n, true)};
    return finish(cgpu_run(j, ops, 3));
}

bool Cgerc_gpu(INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda)
{
    
    if (m <= 0 || n <= 0 || incx <= 0 || incy <= 0 || lda < std::max((INTEGER)1, m))
        return false;
    if (!eligible(false, 4.0 * m * n))
        return false;
    typedef job_Cgerc<OPS> J;
    typedef J::E E;
    J j;
    j.m = m, j.n = n, j.incx = incx, j.incy = incy, j.lda = lda; sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {op<E>(x, incx, 1, m, false), op<E>(y, incy, 1, n, false), op<E>(a, lda, m, n, true)};
    return finish(cgpu_run(j, ops, 3));
}

bool Caxpy_gpu(INTEGER const n, COMPLEX const za, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy)
{
    
    if (n <= 0 || incx <= 0 || incy <= 0)
        return false;
    if (!eligible(false, 4.0 * n))
        return false;
    typedef job_Caxpy<OPS> J;
    typedef J::E E;
    J j;
    j.n = n, j.incx = incx, j.incy = incy; sc(j.alpha, za);
    const gpu_operand<E> ops[] = {op<E>(zx, incx, 1, n, false), op<E>(zy, incy, 1, n, true)};
    return finish(cgpu_run(j, ops, 2));
}

bool Ccopy_gpu(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy)
{
    
    if (n <= 0 || incx <= 0 || incy <= 0)
        return false;
    if (!eligible(false, 4.0 * n))
        return false;
    typedef job_Ccopy<OPS> J;
    typedef J::E E;
    J j;
    j.n = n, j.incx = incx, j.incy = incy;
    const gpu_operand<E> ops[] = {op<E>(zx, incx, 1, n, false), op<E>(zy, incy, 1, n, true)};
    return finish(cgpu_run(j, ops, 2));
}

bool Cswap_gpu(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy)
{
    
    if (n <= 0 || incx <= 0 || incy <= 0)
        return false;
    if (!eligible(false, 4.0 * n))
        return false;
    typedef job_Cswap<OPS> J;
    typedef J::E E;
    J j;
    j.n = n, j.incx = incx, j.incy = incy;
    const gpu_operand<E> ops[] = {op<E>(zx, incx, 1, n, true), op<E>(zy, incy, 1, n, true)};
    return finish(cgpu_run(j, ops, 2));
}

bool CRrot_gpu(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy, REAL const c, REAL const s)
{
    
    if (n <= 0 || incx <= 0 || incy <= 0)
        return false;
    if (!eligible(false, 4.0 * n))
        return false;
    typedef job_CRrot<OPS> J;
    typedef J::E E;
    J j;
    j.n = n, j.incx = incx, j.incy = incy; sc(j.c, c); sc(j.s, s);
    const gpu_operand<E> ops[] = {op<E>(zx, incx, 1, n, true), op<E>(zy, incy, 1, n, true)};
    return finish(cgpu_run(j, ops, 2));
}

bool Cscal_gpu(INTEGER const n, COMPLEX const za, COMPLEX *zx, INTEGER const incx)
{
    
    if (n <= 0 || incx <= 0)
        return false;
    if (!eligible(false, 4.0 * n))
        return false;
    typedef job_Cscal<OPS> J;
    typedef J::E E;
    J j;
    j.n = n, j.incx = incx; sc(j.alpha, za);
    const gpu_operand<E> ops[] = {op<E>(zx, incx, 1, n, true)};
    return finish(cgpu_run(j, ops, 1));
}

bool CRscal_gpu(INTEGER const n, REAL const da, COMPLEX *zx, INTEGER const incx)
{
    
    if (n <= 0 || incx <= 0)
        return false;
    if (!eligible(false, 4.0 * n))
        return false;
    typedef job_CRscal<OPS> J;
    typedef J::E E;
    J j;
    j.n = n, j.incx = incx; sc(j.da, da);
    const gpu_operand<E> ops[] = {op<E>(zx, incx, 1, n, true)};
    return finish(cgpu_run(j, ops, 1));
}

bool Csyrk_gpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool notr = ls(trans, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!notr && !ls(trans, "T")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, notr ? n : k) || ldc < std::max((INTEGER)1, n))
        return false;
    if (!eligible(true, 4.0 * n * n * k))
        return false;
    typedef job_Csyrk<OPS> J;
    typedef J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); j.n = n, j.k = k, j.lda = lda, j.ldc = ldc; sc(j.alpha, alpha); sc(j.beta, beta);
    const gpu_operand<E> ops[] = {op<E>(a, lda, notr ? n : k, notr ? k : n, false), op<E>(c, ldc, n, n, true)};
    return finish(cgpu_run(j, ops, 2));
}

bool Cherk_gpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, COMPLEX *a, INTEGER const lda, REAL const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool notr = ls(trans, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!notr && !ls(trans, "C")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, notr ? n : k) || ldc < std::max((INTEGER)1, n))
        return false;
    if (!eligible(true, 4.0 * n * n * k))
        return false;
    typedef job_Cherk<OPS> J;
    typedef J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); j.n = n, j.k = k, j.lda = lda, j.ldc = ldc; sc(j.alpha, alpha); sc(j.beta, beta);
    const gpu_operand<E> ops[] = {op<E>(a, lda, notr ? n : k, notr ? k : n, false), op<E>(c, ldc, n, n, true)};
    return finish(cgpu_run(j, ops, 2));
}

bool Csyr2k_gpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool notr = ls(trans, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!notr && !ls(trans, "T")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, notr ? n : k) || ldc < std::max((INTEGER)1, n) || ldb < std::max((INTEGER)1, notr ? n : k))
        return false;
    if (!eligible(true, 4.0 * n * n * k))
        return false;
    typedef job_Csyr2k<OPS> J;
    typedef J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); j.n = n, j.k = k, j.lda = lda, j.ldc = ldc; j.ldb = ldb; sc(j.alpha, alpha); sc(j.beta, beta);
    const gpu_operand<E> ops[] = {op<E>(a, lda, notr ? n : k, notr ? k : n, false), op<E>(b, ldb, notr ? n : k, notr ? k : n, false), op<E>(c, ldc, n, n, true)};
    return finish(cgpu_run(j, ops, 3));
}

bool Cher2k_gpu(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool notr = ls(trans, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!notr && !ls(trans, "C")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, notr ? n : k) || ldc < std::max((INTEGER)1, n) || ldb < std::max((INTEGER)1, notr ? n : k))
        return false;
    if (!eligible(true, 4.0 * n * n * k))
        return false;
    typedef job_Cher2k<OPS> J;
    typedef J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); j.n = n, j.k = k, j.lda = lda, j.ldc = ldc; j.ldb = ldb; sc(j.alpha, alpha); sc(j.beta, beta);
    const gpu_operand<E> ops[] = {op<E>(a, lda, notr ? n : k, notr ? k : n, false), op<E>(b, ldb, notr ? n : k, notr ? k : n, false), op<E>(c, ldc, n, n, true)};
    return finish(cgpu_run(j, ops, 3));
}

bool Cgemmtr_gpu(const char *uplo, const char *transa, const char *transb, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool nota = ls(transa, "N"), notb = ls(transb, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!nota && !ls(transa, "T") && !ls(transa, "C")) || (!notb && !ls(transb, "T") && !ls(transb, "C")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, nota ? n : k) || ldb < std::max((INTEGER)1, notb ? k : n) || ldc < std::max((INTEGER)1, n))
        return false;
    if (!eligible(true, 2.0 * n * n * k))
        return false;
    typedef job_Cgemmtr<OPS> J;
    typedef J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.ta, transa); ch(j.tb, transb); j.n = n, j.k = k, j.lda = lda, j.ldb = ldb, j.ldc = ldc; sc(j.alpha, alpha); sc(j.beta, beta);
    const gpu_operand<E> ops[] = {op<E>(a, lda, nota ? n : k, nota ? k : n, false), op<E>(b, ldb, notb ? k : n, notb ? n : k, false), op<E>(c, ldc, n, n, true)};
    return finish(cgpu_run(j, ops, 3));
}

bool Rgemmtr_gpu(const char *uplo, const char *transa, const char *transb, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const beta, REAL *c, INTEGER const ldc)
{
    const bool nota = ls(transa, "N"), notb = ls(transb, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!nota && !ls(transa, "T") && !ls(transa, "C")) || (!notb && !ls(transb, "T") && !ls(transb, "C")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, nota ? n : k) || ldb < std::max((INTEGER)1, notb ? k : n) || ldc < std::max((INTEGER)1, n))
        return false;
    if (!eligible(true, 2.0 * n * n * k))
        return false;
    typedef job_Rgemmtr<OPS> J;
    typedef J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.ta, transa); ch(j.tb, transb); j.n = n, j.k = k, j.lda = lda, j.ldb = ldb, j.ldc = ldc; sc(j.alpha, alpha); sc(j.beta, beta);
    const gpu_operand<E> ops[] = {op<E>(a, lda, nota ? n : k, nota ? k : n, false), op<E>(b, ldb, notb ? k : n, notb ? n : k, false), op<E>(c, ldc, n, n, true)};
    return finish(cgpu_run(j, ops, 3));
}

bool Cher_gpu(const char *uplo, INTEGER const n, REAL const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *a, INTEGER const lda)
{
    
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || lda < std::max((INTEGER)1, n))
        return false;
    if (!eligible(false, 4.0 * n * n))
        return false;
    typedef job_Cher<OPS> J;
    typedef J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n, j.incx = incx, j.lda = lda; sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {op<E>(x, incx, 1, n, false), op<E>(a, lda, n, n, true)};
    return finish(cgpu_run(j, ops, 2));
}

bool Cher2_gpu(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda)
{
    
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || lda < std::max((INTEGER)1, n))
        return false;
    if (!eligible(false, 8.0 * n * n))
        return false;
    typedef job_Cher2<OPS> J;
    typedef J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n, j.incx = incx, j.incy = incy, j.lda = lda; sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {op<E>(x, incx, 1, n, false), op<E>(y, incy, 1, n, false), op<E>(a, lda, n, n, true)};
    return finish(cgpu_run(j, ops, 3));
}
