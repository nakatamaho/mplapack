/*
 * Copyright (c) 2008-2025
 *	Nakata, Maho
 * 	All rights reserved.
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

#ifndef _MPBLAS_TD_H_
#define _MPBLAS_TD_H_

#include <qd/td_real.h>
#include <qd/td_complex.h>
#include <mplapack_config.h>
#include <mplapack_utils_td.h>

bool Mlsame_td(const char *a, const char *b);
mplapackint iCamax(mplapackint const n, td_complex *zx, mplapackint const incx);
mplapackint iRamax(mplapackint const n, td_real *dx, mplapackint const incx);
td_complex Cdotc(mplapackint const n, td_complex *zx, mplapackint const incx, td_complex *zy, mplapackint const incy);
td_complex Cdotu(mplapackint const n, td_complex *zx, mplapackint const incx, td_complex *zy, mplapackint const incy);
td_real RCabs1(td_complex const z);
td_real RCasum(mplapackint const n, td_complex *zx, mplapackint const incx);
td_real RCnrm2(mplapackint const n, td_complex *x, mplapackint const incx);
td_real Rasum(mplapackint const n, td_real *dx, mplapackint const incx);
td_real Rdot(mplapackint const n, td_real *dx, mplapackint const incx, td_real *dy, mplapackint const incy);
td_real Rnrm2(mplapackint const n, td_real *x, mplapackint const incx);
void CRrot(mplapackint const n, td_complex *zx, mplapackint const incx, td_complex *zy, mplapackint const incy, td_real const c, td_real const s);
void CRscal(mplapackint const n, td_real const da, td_complex *zx, mplapackint const incx);
void Caxpy(mplapackint const n, td_complex const za, td_complex *zx, mplapackint const incx, td_complex *zy, mplapackint const incy);
void Ccopy(mplapackint const n, td_complex *zx, mplapackint const incx, td_complex *zy, mplapackint const incy);
void Cgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *x, mplapackint const incx, td_complex const beta, td_complex *y, mplapackint const incy);
void Cgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb, td_complex const beta, td_complex *c, mplapackint const ldc);
void Cgemmtr(const char *uplo, const char *transa, const char *transb, mplapackint const n, mplapackint const k, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb, td_complex const beta, td_complex *c, mplapackint const ldc);
void Cgemv(const char *trans, mplapackint const m, mplapackint const n, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *x, mplapackint const incx, td_complex const beta, td_complex *y, mplapackint const incy);
void Cgerc(mplapackint const m, mplapackint const n, td_complex const alpha, td_complex *x, mplapackint const incx, td_complex *y, mplapackint const incy, td_complex *a, mplapackint const lda);
void Cgeru(mplapackint const m, mplapackint const n, td_complex const alpha, td_complex *x, mplapackint const incx, td_complex *y, mplapackint const incy, td_complex *a, mplapackint const lda);
void Chbmv(const char *uplo, mplapackint const n, mplapackint const k, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *x, mplapackint const incx, td_complex const beta, td_complex *y, mplapackint const incy);
void Chemm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb, td_complex const beta, td_complex *c, mplapackint const ldc);
void Chemv(const char *uplo, mplapackint const n, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *x, mplapackint const incx, td_complex const beta, td_complex *y, mplapackint const incy);
void Cher(const char *uplo, mplapackint const n, td_real const alpha, td_complex *x, mplapackint const incx, td_complex *a, mplapackint const lda);
void Cher2(const char *uplo, mplapackint const n, td_complex const alpha, td_complex *x, mplapackint const incx, td_complex *y, mplapackint const incy, td_complex *a, mplapackint const lda);
void Cher2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb, td_real const beta, td_complex *c, mplapackint const ldc);
void Cherk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, td_real const alpha, td_complex *a, mplapackint const lda, td_real const beta, td_complex *c, mplapackint const ldc);
void Chpmv(const char *uplo, mplapackint const n, td_complex const alpha, td_complex *ap, td_complex *x, mplapackint const incx, td_complex const beta, td_complex *y, mplapackint const incy);
void Chpr(const char *uplo, mplapackint const n, td_real const alpha, td_complex *x, mplapackint const incx, td_complex *ap);
void Chpr2(const char *uplo, mplapackint const n, td_complex const alpha, td_complex *x, mplapackint const incx, td_complex *y, mplapackint const incy, td_complex *ap);
void Crotg(td_complex &a, td_complex const b, td_real &c, td_complex &s);
void Cscal(mplapackint const n, td_complex const za, td_complex *zx, mplapackint const incx);
void Cswap(mplapackint const n, td_complex *zx, mplapackint const incx, td_complex *zy, mplapackint const incy);
void Csymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb, td_complex const beta, td_complex *c, mplapackint const ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb, td_complex const beta, td_complex *c, mplapackint const ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex const beta, td_complex *c, mplapackint const ldc);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, td_complex *a, mplapackint const lda, td_complex *x, mplapackint const incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, td_complex *a, mplapackint const lda, td_complex *x, mplapackint const incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, td_complex *ap, td_complex *x, mplapackint const incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, td_complex *ap, td_complex *x, mplapackint const incx);
void Ctrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, td_complex *a, mplapackint const lda, td_complex *x, mplapackint const incx);
void Ctrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, td_complex const alpha, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, td_complex *a, mplapackint const lda, td_complex *x, mplapackint const incx);
void Mxerbla_td(const char *srname, int info);
void Raxpy(mplapackint const n, td_real const da, td_real *dx, mplapackint const incx, td_real *dy, mplapackint const incy);
void Rcopy(mplapackint const n, td_real *dx, mplapackint const incx, td_real *dy, mplapackint const incy);
void Rgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, td_real const alpha, td_real *a, mplapackint const lda, td_real *x, mplapackint const incx, td_real const beta, td_real *y, mplapackint const incy);
void Rgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, td_real const alpha, td_real *a, mplapackint const lda, td_real *b, mplapackint const ldb, td_real const beta, td_real *c, mplapackint const ldc);
void Rgemmtr(const char *uplo, const char *transa, const char *transb, mplapackint const n, mplapackint const k, td_real const alpha, td_real *a, mplapackint const lda, td_real *b, mplapackint const ldb, td_real const beta, td_real *c, mplapackint const ldc);
void Rgemv(const char *trans, mplapackint const m, mplapackint const n, td_real const alpha, td_real *a, mplapackint const lda, td_real *x, mplapackint const incx, td_real const beta, td_real *y, mplapackint const incy);
void Rger(mplapackint const m, mplapackint const n, td_real const alpha, td_real *x, mplapackint const incx, td_real *y, mplapackint const incy, td_real *a, mplapackint const lda);
void Rrot(mplapackint const n, td_real *dx, mplapackint const incx, td_real *dy, mplapackint const incy, td_real const c, td_real const s);
void Rrotg(td_real &a, td_real &b, td_real &c, td_real &s);
void Rrotm(mplapackint const n, td_real *dx, mplapackint const incx, td_real *dy, mplapackint const incy, td_real *dparam);
void Rrotmg(td_real &dd1, td_real &dd2, td_real &dx1, td_real const dy1, td_real *dparam);
void Rsbmv(const char *uplo, mplapackint const n, mplapackint const k, td_real const alpha, td_real *a, mplapackint const lda, td_real *x, mplapackint const incx, td_real const beta, td_real *y, mplapackint const incy);
void Rscal(mplapackint const n, td_real const da, td_real *dx, mplapackint const incx);
void Rspmv(const char *uplo, mplapackint const n, td_real const alpha, td_real *ap, td_real *x, mplapackint const incx, td_real const beta, td_real *y, mplapackint const incy);
void Rspr(const char *uplo, mplapackint const n, td_real const alpha, td_real *x, mplapackint const incx, td_real *ap);
void Rspr2(const char *uplo, mplapackint const n, td_real const alpha, td_real *x, mplapackint const incx, td_real *y, mplapackint const incy, td_real *ap);
void Rswap(mplapackint const n, td_real *dx, mplapackint const incx, td_real *dy, mplapackint const incy);
void Rsymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, td_real const alpha, td_real *a, mplapackint const lda, td_real *b, mplapackint const ldb, td_real const beta, td_real *c, mplapackint const ldc);
void Rsymv(const char *uplo, mplapackint const n, td_real const alpha, td_real *a, mplapackint const lda, td_real *x, mplapackint const incx, td_real const beta, td_real *y, mplapackint const incy);
void Rsyr(const char *uplo, mplapackint const n, td_real const alpha, td_real *x, mplapackint const incx, td_real *a, mplapackint const lda);
void Rsyr2(const char *uplo, mplapackint const n, td_real const alpha, td_real *x, mplapackint const incx, td_real *y, mplapackint const incy, td_real *a, mplapackint const lda);
void Rsyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, td_real const alpha, td_real *a, mplapackint const lda, td_real *b, mplapackint const ldb, td_real const beta, td_real *c, mplapackint const ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, td_real const alpha, td_real *a, mplapackint const lda, td_real const beta, td_real *c, mplapackint const ldc);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, td_real *a, mplapackint const lda, td_real *x, mplapackint const incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, td_real *a, mplapackint const lda, td_real *x, mplapackint const incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, td_real *ap, td_real *x, mplapackint const incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, td_real *ap, td_real *x, mplapackint const incx);
void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, td_real const alpha, td_real *a, mplapackint const lda, td_real *b, mplapackint const ldb);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, td_real *a, mplapackint const lda, td_real *x, mplapackint const incx);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, td_real const alpha, td_real *a, mplapackint const lda, td_real *b, mplapackint const ldb);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, td_real *a, mplapackint const lda, td_real *x, mplapackint const incx);
#endif
