/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpblas_dd.h,v 1.10 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPBLAS_DD_H_
#define _MPBLAS_DD_H_

#include <qd/dd_real.h>
#include <dd_complex.h>
#include <mplapack_config.h>
#include <mplapack_utils_dd.h>

bool Mlsame_dd(const char *a, const char *b);
dd_complex Cdotc(mplapackint const n, dd_complex *zx, mplapackint const incx, dd_complex *zy, mplapackint const incy);
dd_complex Cdotu(mplapackint const n, dd_complex *zx, mplapackint const incx, dd_complex *zy, mplapackint const incy);
dd_real RCabs1(dd_complex const z);
dd_real RCasum(mplapackint const n, dd_complex *zx, mplapackint const incx);
dd_real RCnrm2(mplapackint const n, dd_complex *x, mplapackint const incx);
dd_real Rasum(mplapackint const n, dd_real *dx, mplapackint const incx);
dd_real Rdot(mplapackint const n, dd_real *dx, mplapackint const incx, dd_real *dy, mplapackint const incy);
dd_real Rnrm2(mplapackint const n, dd_real *x, mplapackint const incx);
mplapackint iCamax(mplapackint const n, dd_complex *zx, mplapackint const incx);
mplapackint iRamax(mplapackint const n, dd_real *dx, mplapackint const incx);
void CRrot(mplapackint const n, dd_complex *zx, mplapackint const incx, dd_complex *zy, mplapackint const incy, dd_real const c, dd_real const s);
void CRscal(mplapackint const n, dd_real const da, dd_complex *zx, mplapackint const incx);
void Caxpy(mplapackint const n, dd_complex const za, dd_complex *zx, mplapackint const incx, dd_complex *zy, mplapackint const incy);
void Ccopy(mplapackint const n, dd_complex *zx, mplapackint const incx, dd_complex *zy, mplapackint const incy);
void Cgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *x, mplapackint const incx, dd_complex const beta, dd_complex *y, mplapackint const incy);
void Cgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *b, mplapackint const ldb, dd_complex const beta, dd_complex *c, mplapackint const ldc);
void Cgemv(const char *trans, mplapackint const m, mplapackint const n, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *x, mplapackint const incx, dd_complex const beta, dd_complex *y, mplapackint const incy);
void Cgerc(mplapackint const m, mplapackint const n, dd_complex const alpha, dd_complex *x, mplapackint const incx, dd_complex *y, mplapackint const incy, dd_complex *a, mplapackint const lda);
void Cgeru(mplapackint const m, mplapackint const n, dd_complex const alpha, dd_complex *x, mplapackint const incx, dd_complex *y, mplapackint const incy, dd_complex *a, mplapackint const lda);
void Chbmv(const char *uplo, mplapackint const n, mplapackint const k, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *x, mplapackint const incx, dd_complex const beta, dd_complex *y, mplapackint const incy);
void Chemm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *b, mplapackint const ldb, dd_complex const beta, dd_complex *c, mplapackint const ldc);
void Chemv(const char *uplo, mplapackint const n, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *x, mplapackint const incx, dd_complex const beta, dd_complex *y, mplapackint const incy);
void Cher(const char *uplo, mplapackint const n, dd_real const alpha, dd_complex *x, mplapackint const incx, dd_complex *a, mplapackint const lda);
void Cher2(const char *uplo, mplapackint const n, dd_complex const alpha, dd_complex *x, mplapackint const incx, dd_complex *y, mplapackint const incy, dd_complex *a, mplapackint const lda);
void Cher2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *b, mplapackint const ldb, dd_real const beta, dd_complex *c, mplapackint const ldc);
void Cherk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, dd_real const alpha, dd_complex *a, mplapackint const lda, dd_real const beta, dd_complex *c, mplapackint const ldc);
void Chpmv(const char *uplo, mplapackint const n, dd_complex const alpha, dd_complex *ap, dd_complex *x, mplapackint const incx, dd_complex const beta, dd_complex *y, mplapackint const incy);
void Chpr(const char *uplo, mplapackint const n, dd_real const alpha, dd_complex *x, mplapackint const incx, dd_complex *ap);
void Chpr2(const char *uplo, mplapackint const n, dd_complex const alpha, dd_complex *x, mplapackint const incx, dd_complex *y, mplapackint const incy, dd_complex *ap);
void Crotg(dd_complex &ca, dd_complex const cb, dd_real &c, dd_complex &s);
void Cscal(mplapackint const n, dd_complex const za, dd_complex *zx, mplapackint const incx);
void Cswap(mplapackint const n, dd_complex *zx, mplapackint const incx, dd_complex *zy, mplapackint const incy);
void Csymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *b, mplapackint const ldb, dd_complex const beta, dd_complex *c, mplapackint const ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *b, mplapackint const ldb, dd_complex const beta, dd_complex *c, mplapackint const ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex const beta, dd_complex *c, mplapackint const ldc);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, dd_complex *a, mplapackint const lda, dd_complex *x, mplapackint const incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, dd_complex *a, mplapackint const lda, dd_complex *x, mplapackint const incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, dd_complex *ap, dd_complex *x, mplapackint const incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, dd_complex *ap, dd_complex *x, mplapackint const incx);
void Ctrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *b, mplapackint const ldb);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, dd_complex *a, mplapackint const lda, dd_complex *x, mplapackint const incx);
void Ctrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, dd_complex const alpha, dd_complex *a, mplapackint const lda, dd_complex *b, mplapackint const ldb);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, dd_complex *a, mplapackint const lda, dd_complex *x, mplapackint const incx);
void Mxerbla_dd(const char *srname, int info);
void Raxpy(mplapackint const n, dd_real const da, dd_real *dx, mplapackint const incx, dd_real *dy, mplapackint const incy);
void Rcopy(mplapackint const n, dd_real *dx, mplapackint const incx, dd_real *dy, mplapackint const incy);
void Rgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx, dd_real const beta, dd_real *y, mplapackint const incy);
void Rgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb, dd_real const beta, dd_real *c, mplapackint const ldc);
void Rgemv(const char *trans, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx, dd_real const beta, dd_real *y, mplapackint const incy);
void Rger(mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *x, mplapackint const incx, dd_real *y, mplapackint const incy, dd_real *a, mplapackint const lda);
void Rrot(mplapackint const n, dd_real *dx, mplapackint const incx, dd_real *dy, mplapackint const incy, dd_real const c, dd_real const s);
void Rrotg(dd_real &da, dd_real &db, dd_real &c, dd_real &s);
void Rrotm(mplapackint const n, dd_real *dx, mplapackint const incx, dd_real *dy, mplapackint const incy, dd_real *dparam);
void Rrotmg(dd_real &dd1, dd_real &dd2, dd_real &dx1, dd_real const dy1, dd_real *dparam);
void Rsbmv(const char *uplo, mplapackint const n, mplapackint const k, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx, dd_real const beta, dd_real *y, mplapackint const incy);
void Rscal(mplapackint const n, dd_real const da, dd_real *dx, mplapackint const incx);
void Rspmv(const char *uplo, mplapackint const n, dd_real const alpha, dd_real *ap, dd_real *x, mplapackint const incx, dd_real const beta, dd_real *y, mplapackint const incy);
void Rspr(const char *uplo, mplapackint const n, dd_real const alpha, dd_real *x, mplapackint const incx, dd_real *ap);
void Rspr2(const char *uplo, mplapackint const n, dd_real const alpha, dd_real *x, mplapackint const incx, dd_real *y, mplapackint const incy, dd_real *ap);
void Rswap(mplapackint const n, dd_real *dx, mplapackint const incx, dd_real *dy, mplapackint const incy);
void Rsymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb, dd_real const beta, dd_real *c, mplapackint const ldc);
void Rsymv(const char *uplo, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx, dd_real const beta, dd_real *y, mplapackint const incy);
void Rsyr(const char *uplo, mplapackint const n, dd_real const alpha, dd_real *x, mplapackint const incx, dd_real *a, mplapackint const lda);
void Rsyr2(const char *uplo, mplapackint const n, dd_real const alpha, dd_real *x, mplapackint const incx, dd_real *y, mplapackint const incy, dd_real *a, mplapackint const lda);
void Rsyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb, dd_real const beta, dd_real *c, mplapackint const ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real const beta, dd_real *c, mplapackint const ldc);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, dd_real *ap, dd_real *x, mplapackint const incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, dd_real *ap, dd_real *x, mplapackint const incx);
void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, dd_real const alpha, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const incx);
#endif
