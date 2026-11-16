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

// The complex Level 1/2/3 BLAS (and Rgemmtr) on CUDA, independent of the
// element type: <R>_gpu_t<Ops>(ctx, arguments) checks the arguments of the
// call (invalid ones and quick returns are left to the CPU), asks
// ctx.eligible(level3, work), fills the job of ../Cblas_cuda_jobs_common.h
// with ctx.sc() for the scalars and ctx.op<E>() for the operands, and returns
// ctx.run(job, operands, count): true when the GPU has computed the result.
// X is the element-type side: Cblas_bridge_cuda.cpp (dd/td/qd: memory
// copies) and mpfr/cuda/Cblas_mpfr_bridge_cuda.cpp (MPFR: packing).
// Include after <mpblas.h> and the device types.

#ifndef _MPLAPACK_CBLAS_BRIDGE_BODY_H_
#define _MPLAPACK_CBLAS_BRIDGE_BODY_H_

#include <algorithm>
#include "../Cblas_cuda_jobs_common.h"
#include "../Rblas_cuda_entry_common.h"

namespace mplapack_cblas_gpu {

using namespace mplapack_blas_dev;
using mplapack_blas::gpu_operand;

inline bool ls(const char *a, const char *b) { return Mlsame(a, b) != 0; }
inline void ch(char *d, const char *s) { d[0] = s[0], d[1] = '\0'; }

template <class Ops, class X> bool Cgemm_gpu_t(X &ctx, const char *transa, const char *transb, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool nota = ls(transa, "N"), notb = ls(transb, "N");
    if ((!nota && !ls(transa, "T") && !ls(transa, "C")) || (!notb && !ls(transb, "T") && !ls(transb, "C")) || m <= 0 || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, nota ? m : k) || ldb < std::max((INTEGER)1, notb ? k : n) || ldc < std::max((INTEGER)1, m))
        return false;
    if (!ctx.eligible(true, 4.0 * m * n * k))
        return false;
    typedef job_Cgemm<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.ta, transa); ch(j.tb, transb); j.m = m, j.n = n, j.k = k, j.lda = lda, j.ldb = ldb, j.ldc = ldc; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, nota ? m : k, nota ? k : m, false), ctx.template op<E>(b, ldb, notb ? k : n, notb ? n : k, false), ctx.template op<E>(c, ldc, m, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Csymm_gpu_t(X &ctx, const char *side, const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool left = ls(side, "L");
    if ((!left && !ls(side, "R")) || (!ls(uplo, "U") && !ls(uplo, "L")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, left ? m : n) || ldb < std::max((INTEGER)1, m) || ldc < std::max((INTEGER)1, m))
        return false;
    if (!ctx.eligible(true, 4.0 * m * n * (left ? m : n)))
        return false;
    typedef job_Csymm<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.side, side); ch(j.uplo, uplo); j.m = m, j.n = n, j.lda = lda, j.ldb = ldb, j.ldc = ldc; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, left ? m : n, left ? m : n, false), ctx.template op<E>(b, ldb, m, n, false), ctx.template op<E>(c, ldc, m, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Chemm_gpu_t(X &ctx, const char *side, const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool left = ls(side, "L");
    if ((!left && !ls(side, "R")) || (!ls(uplo, "U") && !ls(uplo, "L")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, left ? m : n) || ldb < std::max((INTEGER)1, m) || ldc < std::max((INTEGER)1, m))
        return false;
    if (!ctx.eligible(true, 4.0 * m * n * (left ? m : n)))
        return false;
    typedef job_Chemm<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.side, side); ch(j.uplo, uplo); j.m = m, j.n = n, j.lda = lda, j.ldb = ldb, j.ldc = ldc; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, left ? m : n, left ? m : n, false), ctx.template op<E>(b, ldb, m, n, false), ctx.template op<E>(c, ldc, m, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Ctrmm_gpu_t(X &ctx, const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb)
{
    const bool left = ls(side, "L");
    if ((!left && !ls(side, "R")) || (!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(transa, "N") && !ls(transa, "T") && !ls(transa, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, left ? m : n) || ldb < std::max((INTEGER)1, m))
        return false;
    if (!ctx.eligible(true, 4.0 * m * n * (left ? m : n)))
        return false;
    typedef job_Ctrmm<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.side, side); ch(j.uplo, uplo); ch(j.ta, transa); ch(j.diag, diag); j.m = m, j.n = n, j.lda = lda, j.ldb = ldb; ctx.sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, left ? m : n, left ? m : n, false), ctx.template op<E>(b, ldb, m, n, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool Ctrsm_gpu_t(X &ctx, const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb)
{
    const bool left = ls(side, "L");
    if ((!left && !ls(side, "R")) || (!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(transa, "N") && !ls(transa, "T") && !ls(transa, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, left ? m : n) || ldb < std::max((INTEGER)1, m))
        return false;
    if (!ctx.eligible(true, 4.0 * m * n * (left ? m : n)))
        return false;
    typedef job_Ctrsm<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.side, side); ch(j.uplo, uplo); ch(j.ta, transa); ch(j.diag, diag); j.m = m, j.n = n, j.lda = lda, j.ldb = ldb; ctx.sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, left ? m : n, left ? m : n, false), ctx.template op<E>(b, ldb, m, n, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool Cgemv_gpu_t(X &ctx, const char *trans, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy)
{
    const bool notr = ls(trans, "N");
    if ((!notr && !ls(trans, "T") && !ls(trans, "C")) || m <= 0 || n <= 0 || lda < std::max((INTEGER)1, m) || incx <= 0 || incy <= 0)
        return false;
    if (!ctx.eligible(false, 4.0 * m * n))
        return false;
    typedef job_Cgemv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.tr, trans); j.m = m, j.n = n, j.lda = lda, j.incx = incx, j.incy = incy; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, m, n, false), ctx.template op<E>(x, incx, 1, notr ? n : m, false), ctx.template op<E>(y, incy, 1, notr ? m : n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Cgeru_gpu_t(X &ctx, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda)
{
    if (m <= 0 || n <= 0 || incx <= 0 || incy <= 0 || lda < std::max((INTEGER)1, m))
        return false;
    if (!ctx.eligible(false, 4.0 * m * n))
        return false;
    typedef job_Cgeru<Ops> J;
    typedef typename J::E E;
    J j;
    j.m = m, j.n = n, j.incx = incx, j.incy = incy, j.lda = lda; ctx.sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {ctx.template op<E>(x, incx, 1, m, false), ctx.template op<E>(y, incy, 1, n, false), ctx.template op<E>(a, lda, m, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Cgerc_gpu_t(X &ctx, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda)
{
    if (m <= 0 || n <= 0 || incx <= 0 || incy <= 0 || lda < std::max((INTEGER)1, m))
        return false;
    if (!ctx.eligible(false, 4.0 * m * n))
        return false;
    typedef job_Cgerc<Ops> J;
    typedef typename J::E E;
    J j;
    j.m = m, j.n = n, j.incx = incx, j.incy = incy, j.lda = lda; ctx.sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {ctx.template op<E>(x, incx, 1, m, false), ctx.template op<E>(y, incy, 1, n, false), ctx.template op<E>(a, lda, m, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Caxpy_gpu_t(X &ctx, INTEGER const n, COMPLEX const za, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy)
{
    if (n <= 0 || incx <= 0 || incy <= 0)
        return false;
    if (!ctx.eligible(false, 4.0 * n))
        return false;
    typedef job_Caxpy<Ops> J;
    typedef typename J::E E;
    J j;
    j.n = n, j.incx = incx, j.incy = incy; ctx.sc(j.alpha, za);
    const gpu_operand<E> ops[] = {ctx.template op<E>(zx, incx, 1, n, false), ctx.template op<E>(zy, incy, 1, n, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool Ccopy_gpu_t(X &ctx, INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy)
{
    if (n <= 0 || incx <= 0 || incy <= 0)
        return false;
    if (!ctx.eligible(false, 4.0 * n))
        return false;
    typedef job_Ccopy<Ops> J;
    typedef typename J::E E;
    J j;
    j.n = n, j.incx = incx, j.incy = incy;
    const gpu_operand<E> ops[] = {ctx.template op<E>(zx, incx, 1, n, false), ctx.template op<E>(zy, incy, 1, n, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool Cswap_gpu_t(X &ctx, INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy)
{
    if (n <= 0 || incx <= 0 || incy <= 0)
        return false;
    if (!ctx.eligible(false, 4.0 * n))
        return false;
    typedef job_Cswap<Ops> J;
    typedef typename J::E E;
    J j;
    j.n = n, j.incx = incx, j.incy = incy;
    const gpu_operand<E> ops[] = {ctx.template op<E>(zx, incx, 1, n, true), ctx.template op<E>(zy, incy, 1, n, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool CRrot_gpu_t(X &ctx, INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy, REAL const c, REAL const s)
{
    if (n <= 0 || incx <= 0 || incy <= 0)
        return false;
    if (!ctx.eligible(false, 4.0 * n))
        return false;
    typedef job_CRrot<Ops> J;
    typedef typename J::E E;
    J j;
    j.n = n, j.incx = incx, j.incy = incy; ctx.sc(j.c, c); ctx.sc(j.s, s);
    const gpu_operand<E> ops[] = {ctx.template op<E>(zx, incx, 1, n, true), ctx.template op<E>(zy, incy, 1, n, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool Cscal_gpu_t(X &ctx, INTEGER const n, COMPLEX const za, COMPLEX *zx, INTEGER const incx)
{
    if (n <= 0 || incx <= 0)
        return false;
    if (!ctx.eligible(false, 4.0 * n))
        return false;
    typedef job_Cscal<Ops> J;
    typedef typename J::E E;
    J j;
    j.n = n, j.incx = incx; ctx.sc(j.alpha, za);
    const gpu_operand<E> ops[] = {ctx.template op<E>(zx, incx, 1, n, true)};
    return ctx.run(j, ops, 1);
}

template <class Ops, class X> bool CRscal_gpu_t(X &ctx, INTEGER const n, REAL const da, COMPLEX *zx, INTEGER const incx)
{
    if (n <= 0 || incx <= 0)
        return false;
    if (!ctx.eligible(false, 4.0 * n))
        return false;
    typedef job_CRscal<Ops> J;
    typedef typename J::E E;
    J j;
    j.n = n, j.incx = incx; ctx.sc(j.da, da);
    const gpu_operand<E> ops[] = {ctx.template op<E>(zx, incx, 1, n, true)};
    return ctx.run(j, ops, 1);
}

template <class Ops, class X> bool Csyrk_gpu_t(X &ctx, const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool notr = ls(trans, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!notr && !ls(trans, "T")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, notr ? n : k) || ldc < std::max((INTEGER)1, n))
        return false;
    if (!ctx.eligible(true, 4.0 * n * n * k))
        return false;
    typedef job_Csyrk<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); j.n = n, j.k = k, j.lda = lda, j.ldc = ldc; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, notr ? n : k, notr ? k : n, false), ctx.template op<E>(c, ldc, n, n, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool Cherk_gpu_t(X &ctx, const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, COMPLEX *a, INTEGER const lda, REAL const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool notr = ls(trans, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!notr && !ls(trans, "C")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, notr ? n : k) || ldc < std::max((INTEGER)1, n))
        return false;
    if (!ctx.eligible(true, 4.0 * n * n * k))
        return false;
    typedef job_Cherk<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); j.n = n, j.k = k, j.lda = lda, j.ldc = ldc; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, notr ? n : k, notr ? k : n, false), ctx.template op<E>(c, ldc, n, n, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool Csyr2k_gpu_t(X &ctx, const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool notr = ls(trans, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!notr && !ls(trans, "T")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, notr ? n : k) || ldc < std::max((INTEGER)1, n) || ldb < std::max((INTEGER)1, notr ? n : k))
        return false;
    if (!ctx.eligible(true, 4.0 * n * n * k))
        return false;
    typedef job_Csyr2k<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); j.n = n, j.k = k, j.lda = lda, j.ldc = ldc; j.ldb = ldb; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, notr ? n : k, notr ? k : n, false), ctx.template op<E>(b, ldb, notr ? n : k, notr ? k : n, false), ctx.template op<E>(c, ldc, n, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Cher2k_gpu_t(X &ctx, const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool notr = ls(trans, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!notr && !ls(trans, "C")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, notr ? n : k) || ldc < std::max((INTEGER)1, n) || ldb < std::max((INTEGER)1, notr ? n : k))
        return false;
    if (!ctx.eligible(true, 4.0 * n * n * k))
        return false;
    typedef job_Cher2k<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); j.n = n, j.k = k, j.lda = lda, j.ldc = ldc; j.ldb = ldb; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, notr ? n : k, notr ? k : n, false), ctx.template op<E>(b, ldb, notr ? n : k, notr ? k : n, false), ctx.template op<E>(c, ldc, n, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Cgemmtr_gpu_t(X &ctx, const char *uplo, const char *transa, const char *transb, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc)
{
    const bool nota = ls(transa, "N"), notb = ls(transb, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!nota && !ls(transa, "T") && !ls(transa, "C")) || (!notb && !ls(transb, "T") && !ls(transb, "C")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, nota ? n : k) || ldb < std::max((INTEGER)1, notb ? k : n) || ldc < std::max((INTEGER)1, n))
        return false;
    if (!ctx.eligible(true, 2.0 * n * n * k))
        return false;
    typedef job_Cgemmtr<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.ta, transa); ch(j.tb, transb); j.n = n, j.k = k, j.lda = lda, j.ldb = ldb, j.ldc = ldc; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, nota ? n : k, nota ? k : n, false), ctx.template op<E>(b, ldb, notb ? k : n, notb ? n : k, false), ctx.template op<E>(c, ldc, n, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Rgemmtr_gpu_t(X &ctx, const char *uplo, const char *transa, const char *transb, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const beta, REAL *c, INTEGER const ldc)
{
    const bool nota = ls(transa, "N"), notb = ls(transb, "N");
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!nota && !ls(transa, "T") && !ls(transa, "C")) || (!notb && !ls(transb, "T") && !ls(transb, "C")) || n <= 0 || k <= 0 || lda < std::max((INTEGER)1, nota ? n : k) || ldb < std::max((INTEGER)1, notb ? k : n) || ldc < std::max((INTEGER)1, n))
        return false;
    if (!ctx.eligible(true, 2.0 * n * n * k))
        return false;
    typedef job_Rgemmtr<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.ta, transa); ch(j.tb, transb); j.n = n, j.k = k, j.lda = lda, j.ldb = ldb, j.ldc = ldc; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, nota ? n : k, nota ? k : n, false), ctx.template op<E>(b, ldb, notb ? k : n, notb ? n : k, false), ctx.template op<E>(c, ldc, n, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Cher_gpu_t(X &ctx, const char *uplo, INTEGER const n, REAL const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *a, INTEGER const lda)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || lda < std::max((INTEGER)1, n))
        return false;
    if (!ctx.eligible(false, 4.0 * n * n))
        return false;
    typedef job_Cher<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n, j.incx = incx, j.lda = lda; ctx.sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {ctx.template op<E>(x, incx, 1, n, false), ctx.template op<E>(a, lda, n, n, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool Cher2_gpu_t(X &ctx, const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || lda < std::max((INTEGER)1, n))
        return false;
    if (!ctx.eligible(false, 8.0 * n * n))
        return false;
    typedef job_Cher2<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n, j.incx = incx, j.incy = incy, j.lda = lda; ctx.sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {ctx.template op<E>(x, incx, 1, n, false), ctx.template op<E>(y, incy, 1, n, false), ctx.template op<E>(a, lda, n, n, true)};
    return ctx.run(j, ops, 3);
}

// Level 2 range variants (Level 2: eligible(false, ...))

template <class Ops, class X> bool Chemv_gpu_t(X &ctx, const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || alpha == COMPLEX(0.0, 0.0) || lda < std::max((INTEGER)1, n))
        return false;
    if (!ctx.eligible(false, 4.0 * n * n))
        return false;
    typedef job_Chemv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n; j.lda = lda; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, n, n, false), ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(y, 1, 1, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Chbmv_gpu_t(X &ctx, const char *uplo, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || alpha == COMPLEX(0.0, 0.0) || k < 0 || lda < k + 1)
        return false;
    if (!ctx.eligible(false, 4.0 * n * (k + 1)))
        return false;
    typedef job_Chbmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n; j.k = k; j.lda = lda; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, k + 1, n, false), ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(y, 1, 1, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Rsbmv_gpu_t(X &ctx, const char *uplo, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || alpha == REAL(0.0) || k < 0 || lda < k + 1)
        return false;
    if (!ctx.eligible(false, 4.0 * n * (k + 1)))
        return false;
    typedef job_Rsbmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n; j.k = k; j.lda = lda; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, k + 1, n, false), ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(y, 1, 1, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Chpmv_gpu_t(X &ctx, const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *ap, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || alpha == COMPLEX(0.0, 0.0))
        return false;
    if (!ctx.eligible(false, 4.0 * n * n))
        return false;
    typedef job_Chpmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(ap, 1, 1, n * (n + 1) / 2, false), ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(y, 1, 1, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Rspmv_gpu_t(X &ctx, const char *uplo, INTEGER const n, REAL const alpha, REAL *ap, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || alpha == REAL(0.0))
        return false;
    if (!ctx.eligible(false, 4.0 * n * n))
        return false;
    typedef job_Rspmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(ap, 1, 1, n * (n + 1) / 2, false), ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(y, 1, 1, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Rgbmv_gpu_t(X &ctx, const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy)
{
    if ((!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || m <= 0 || n <= 0 || kl < 0 || ku < 0 || lda < kl + ku + 1 || incx != 1 || incy != 1 || alpha == REAL(0.0))
        return false;
    if (!ctx.eligible(false, 4.0 * (ls(trans, "N") ? m : n) * (kl + ku + 1)))
        return false;
    typedef job_Rgbmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.tr, trans); j.m = m; j.n = n; j.kl = kl; j.ku = ku; j.lda = lda; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, kl + ku + 1, n, false), ctx.template op<E>(x, 1, 1, ls(trans, "N") ? n : m, false), ctx.template op<E>(y, 1, 1, ls(trans, "N") ? m : n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Cgbmv_gpu_t(X &ctx, const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy)
{
    if ((!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || m <= 0 || n <= 0 || kl < 0 || ku < 0 || lda < kl + ku + 1 || incx != 1 || incy != 1 || alpha == COMPLEX(0.0, 0.0))
        return false;
    if (!ctx.eligible(false, 4.0 * (ls(trans, "N") ? m : n) * (kl + ku + 1)))
        return false;
    typedef job_Cgbmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.tr, trans); j.m = m; j.n = n; j.kl = kl; j.ku = ku; j.lda = lda; ctx.sc(j.alpha, alpha); ctx.sc(j.beta, beta);
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, kl + ku + 1, n, false), ctx.template op<E>(x, 1, 1, ls(trans, "N") ? n : m, false), ctx.template op<E>(y, 1, 1, ls(trans, "N") ? m : n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Ctrmv_gpu_t(X &ctx, const char *uplo, const char *trans, const char *diag, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || n <= 0 || incx != 1 || lda < std::max((INTEGER)1, n))
        return false;
    if (!ctx.eligible(false, 4.0 * n * n))
        return false;
    typedef job_Ctrmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); ch(j.diag, diag); j.n = n; j.lda = lda;
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, n, n, false), ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(x, 1, 1, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Rtbmv_gpu_t(X &ctx, const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const k, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || n <= 0 || incx != 1 || k < 0 || lda < k + 1)
        return false;
    if (!ctx.eligible(false, 4.0 * n * (k + 1)))
        return false;
    typedef job_Rtbmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); ch(j.diag, diag); j.n = n; j.k = k; j.lda = lda;
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, k + 1, n, false), ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(x, 1, 1, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Ctbmv_gpu_t(X &ctx, const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || n <= 0 || incx != 1 || k < 0 || lda < k + 1)
        return false;
    if (!ctx.eligible(false, 4.0 * n * (k + 1)))
        return false;
    typedef job_Ctbmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); ch(j.diag, diag); j.n = n; j.k = k; j.lda = lda;
    const gpu_operand<E> ops[] = {ctx.template op<E>(a, lda, k + 1, n, false), ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(x, 1, 1, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Rtpmv_gpu_t(X &ctx, const char *uplo, const char *trans, const char *diag, INTEGER const n, REAL *ap, REAL *x, INTEGER const incx)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || n <= 0 || incx != 1)
        return false;
    if (!ctx.eligible(false, 4.0 * n * n))
        return false;
    typedef job_Rtpmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); ch(j.diag, diag); j.n = n;
    const gpu_operand<E> ops[] = {ctx.template op<E>(ap, 1, 1, n * (n + 1) / 2, false), ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(x, 1, 1, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Ctpmv_gpu_t(X &ctx, const char *uplo, const char *trans, const char *diag, INTEGER const n, COMPLEX *ap, COMPLEX *x, INTEGER const incx)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || (!ls(trans, "N") && !ls(trans, "T") && !ls(trans, "C")) || (!ls(diag, "U") && !ls(diag, "N")) || n <= 0 || incx != 1)
        return false;
    if (!ctx.eligible(false, 4.0 * n * n))
        return false;
    typedef job_Ctpmv<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); ch(j.tr, trans); ch(j.diag, diag); j.n = n;
    const gpu_operand<E> ops[] = {ctx.template op<E>(ap, 1, 1, n * (n + 1) / 2, false), ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(x, 1, 1, n, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Rspr_gpu_t(X &ctx, const char *uplo, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *ap)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || alpha == REAL(0.0))
        return false;
    if (!ctx.eligible(false, 2.0 * n * n))
        return false;
    typedef job_Rspr<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n; ctx.sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(ap, 1, 1, n * (n + 1) / 2, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool Rspr2_gpu_t(X &ctx, const char *uplo, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *y, INTEGER const incy, REAL *ap)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || alpha == REAL(0.0))
        return false;
    if (!ctx.eligible(false, 4.0 * n * n))
        return false;
    typedef job_Rspr2<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n; ctx.sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(y, 1, 1, n, false), ctx.template op<E>(ap, 1, 1, n * (n + 1) / 2, true)};
    return ctx.run(j, ops, 3);
}

template <class Ops, class X> bool Chpr_gpu_t(X &ctx, const char *uplo, INTEGER const n, REAL const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *ap)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || alpha == REAL(0.0))
        return false;
    if (!ctx.eligible(false, 2.0 * n * n))
        return false;
    typedef job_Chpr<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n; ctx.sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(ap, 1, 1, n * (n + 1) / 2, true)};
    return ctx.run(j, ops, 2);
}

template <class Ops, class X> bool Chpr2_gpu_t(X &ctx, const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *ap)
{
    if ((!ls(uplo, "U") && !ls(uplo, "L")) || n <= 0 || incx != 1 || incy != 1 || alpha == COMPLEX(0.0, 0.0))
        return false;
    if (!ctx.eligible(false, 4.0 * n * n))
        return false;
    typedef job_Chpr2<Ops> J;
    typedef typename J::E E;
    J j;
    ch(j.uplo, uplo); j.n = n; ctx.sc(j.alpha, alpha);
    const gpu_operand<E> ops[] = {ctx.template op<E>(x, 1, 1, n, false), ctx.template op<E>(y, 1, 1, n, false), ctx.template op<E>(ap, 1, 1, n * (n + 1) / 2, true)};
    return ctx.run(j, ops, 3);
}

} // namespace mplapack_cblas_gpu

// The functions R_gpu(arguments) { return RUN(R_gpu_t, arguments); }, for a
// translation unit that defines the macro RUN(F, ...): it calls
// mplapack_cblas_gpu::F<Ops>(ctx, __VA_ARGS__) with its X ctx.
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cgemm(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cgemm)(const char *transa, const char *transb, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc) { return RUN(Cgemm_gpu_t, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Csymm(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Csymm)(const char *side, const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc) { return RUN(Csymm_gpu_t, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Chemm(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Chemm)(const char *side, const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc) { return RUN(Chemm_gpu_t, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Ctrmm(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Ctrmm)(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb) { return RUN(Ctrmm_gpu_t, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Ctrsm(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Ctrsm)(const char *side, const char *uplo, const char *transa, const char *diag, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb) { return RUN(Ctrsm_gpu_t, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cgemv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cgemv)(const char *trans, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy) { return RUN(Cgemv_gpu_t, trans, m, n, alpha, a, lda, x, incx, beta, y, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cgeru(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cgeru)(INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda) { return RUN(Cgeru_gpu_t, m, n, alpha, x, incx, y, incy, a, lda); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cgerc(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cgerc)(INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda) { return RUN(Cgerc_gpu_t, m, n, alpha, x, incx, y, incy, a, lda); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Caxpy(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Caxpy)(INTEGER const n, COMPLEX const za, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy) { return RUN(Caxpy_gpu_t, n, za, zx, incx, zy, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Ccopy(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Ccopy)(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy) { return RUN(Ccopy_gpu_t, n, zx, incx, zy, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cswap(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cswap)(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy) { return RUN(Cswap_gpu_t, n, zx, incx, zy, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_CRrot(RUN) \
    bool MPLAPACK_CUDA_ENTRY(CRrot)(INTEGER const n, COMPLEX *zx, INTEGER const incx, COMPLEX *zy, INTEGER const incy, REAL const c, REAL const s) { return RUN(CRrot_gpu_t, n, zx, incx, zy, incy, c, s); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cscal(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cscal)(INTEGER const n, COMPLEX const za, COMPLEX *zx, INTEGER const incx) { return RUN(Cscal_gpu_t, n, za, zx, incx); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_CRscal(RUN) \
    bool MPLAPACK_CUDA_ENTRY(CRscal)(INTEGER const n, REAL const da, COMPLEX *zx, INTEGER const incx) { return RUN(CRscal_gpu_t, n, da, zx, incx); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Csyrk(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Csyrk)(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX const beta, COMPLEX *c, INTEGER const ldc) { return RUN(Csyrk_gpu_t, uplo, trans, n, k, alpha, a, lda, beta, c, ldc); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cherk(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cherk)(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, COMPLEX *a, INTEGER const lda, REAL const beta, COMPLEX *c, INTEGER const ldc) { return RUN(Cherk_gpu_t, uplo, trans, n, k, alpha, a, lda, beta, c, ldc); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Csyr2k(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Csyr2k)(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc) { return RUN(Csyr2k_gpu_t, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cher2k(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cher2k)(const char *uplo, const char *trans, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL const beta, COMPLEX *c, INTEGER const ldc) { return RUN(Cher2k_gpu_t, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cgemmtr(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cgemmtr)(const char *uplo, const char *transa, const char *transb, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc) { return RUN(Cgemmtr_gpu_t, uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Rgemmtr(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Rgemmtr)(const char *uplo, const char *transa, const char *transb, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const beta, REAL *c, INTEGER const ldc) { return RUN(Rgemmtr_gpu_t, uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cher(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cher)(const char *uplo, INTEGER const n, REAL const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *a, INTEGER const lda) { return RUN(Cher_gpu_t, uplo, n, alpha, x, incx, a, lda); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cher2(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cher2)(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *a, INTEGER const lda) { return RUN(Cher2_gpu_t, uplo, n, alpha, x, incx, y, incy, a, lda); }

#define MPLAPACK_CBLAS_GPU_FUNCTION_Chemv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Chemv)(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy) { return RUN(Chemv_gpu_t, uplo, n, alpha, a, lda, x, incx, beta, y, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Chbmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Chbmv)(const char *uplo, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy) { return RUN(Chbmv_gpu_t, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Rsbmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Rsbmv)(const char *uplo, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy) { return RUN(Rsbmv_gpu_t, uplo, n, k, alpha, a, lda, x, incx, beta, y, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Chpmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Chpmv)(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *ap, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy) { return RUN(Chpmv_gpu_t, uplo, n, alpha, ap, x, incx, beta, y, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Rspmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Rspmv)(const char *uplo, INTEGER const n, REAL const alpha, REAL *ap, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy) { return RUN(Rspmv_gpu_t, uplo, n, alpha, ap, x, incx, beta, y, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Rgbmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Rgbmv)(const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, REAL const alpha, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx, REAL const beta, REAL *y, INTEGER const incy) { return RUN(Rgbmv_gpu_t, trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Cgbmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Cgbmv)(const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx, COMPLEX const beta, COMPLEX *y, INTEGER const incy) { return RUN(Cgbmv_gpu_t, trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Ctrmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Ctrmv)(const char *uplo, const char *trans, const char *diag, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx) { return RUN(Ctrmv_gpu_t, uplo, trans, diag, n, a, lda, x, incx); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Rtbmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Rtbmv)(const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const k, REAL *a, INTEGER const lda, REAL *x, INTEGER const incx) { return RUN(Rtbmv_gpu_t, uplo, trans, diag, n, k, a, lda, x, incx); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Ctbmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Ctbmv)(const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx) { return RUN(Ctbmv_gpu_t, uplo, trans, diag, n, k, a, lda, x, incx); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Rtpmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Rtpmv)(const char *uplo, const char *trans, const char *diag, INTEGER const n, REAL *ap, REAL *x, INTEGER const incx) { return RUN(Rtpmv_gpu_t, uplo, trans, diag, n, ap, x, incx); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Ctpmv(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Ctpmv)(const char *uplo, const char *trans, const char *diag, INTEGER const n, COMPLEX *ap, COMPLEX *x, INTEGER const incx) { return RUN(Ctpmv_gpu_t, uplo, trans, diag, n, ap, x, incx); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Rspr(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Rspr)(const char *uplo, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *ap) { return RUN(Rspr_gpu_t, uplo, n, alpha, x, incx, ap); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Rspr2(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Rspr2)(const char *uplo, INTEGER const n, REAL const alpha, REAL *x, INTEGER const incx, REAL *y, INTEGER const incy, REAL *ap) { return RUN(Rspr2_gpu_t, uplo, n, alpha, x, incx, y, incy, ap); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Chpr(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Chpr)(const char *uplo, INTEGER const n, REAL const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *ap) { return RUN(Chpr_gpu_t, uplo, n, alpha, x, incx, ap); }
#define MPLAPACK_CBLAS_GPU_FUNCTION_Chpr2(RUN) \
    bool MPLAPACK_CUDA_ENTRY(Chpr2)(const char *uplo, INTEGER const n, COMPLEX const alpha, COMPLEX *x, INTEGER const incx, COMPLEX *y, INTEGER const incy, COMPLEX *ap) { return RUN(Chpr2_gpu_t, uplo, n, alpha, x, incx, y, incy, ap); }

#define MPLAPACK_CBLAS_GPU_FUNCTIONS(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cgemm(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Csymm(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Chemm(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Ctrmm(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Ctrsm(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cgemv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cgeru(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cgerc(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Caxpy(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Ccopy(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cswap(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_CRrot(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cscal(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_CRscal(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Csyrk(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cherk(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Csyr2k(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cher2k(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cgemmtr(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Rgemmtr(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cher(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cher2(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Chemv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Chbmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Rsbmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Chpmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Rspmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Rgbmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Cgbmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Ctrmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Rtbmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Ctbmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Rtpmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Ctpmv(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Rspr(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Rspr2(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Chpr(RUN) \
    MPLAPACK_CBLAS_GPU_FUNCTION_Chpr2(RUN)

#endif
