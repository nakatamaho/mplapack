/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpblas_qd.h,v 1.12 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPBLAS_QD_H_
#define _MPBLAS_QD_H_

#include <qd/qd_real.h>
#include "qd_complex.h"
#include "mplapack_config.h"
#include "mplapack_utils_qd.h"

#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

bool Mlsame_qd(const char *a, const char *b);
mplapackint iCamax(mplapackint const n, qd_complex *zx, mplapackint const incx);
mplapackint iRamax(mplapackint const n, qd_real *dx, mplapackint const incx);
qd_complex Cdotc(mplapackint const n, qd_complex *zx, mplapackint const incx, qd_complex *zy, mplapackint const incy);
qd_complex Cdotu(mplapackint const n, qd_complex *zx, mplapackint const incx, qd_complex *zy, mplapackint const incy);
qd_real RCabs1(qd_complex const z);
qd_real RCasum(mplapackint const n, qd_complex *zx, mplapackint const incx);
qd_real RCnrm2(mplapackint const n, qd_complex *x, mplapackint const incx);
qd_real Rasum(mplapackint const n, qd_real *dx, mplapackint const incx);
qd_real Rdot(mplapackint const n, qd_real *dx, mplapackint const incx, qd_real *dy, mplapackint const incy);
qd_real Rnrm2(mplapackint const n, qd_real *x, mplapackint const incx);
void CRrot(mplapackint const n, qd_complex *zx, mplapackint const incx, qd_complex *zy, mplapackint const incy, qd_real const c, qd_real const s);
void CRscal(mplapackint const n, qd_real const da, qd_complex *zx, mplapackint const incx);
void Caxpy(mplapackint const n, qd_complex const za, qd_complex *zx, mplapackint const incx, qd_complex *zy, mplapackint const incy);
void Ccopy(mplapackint const n, qd_complex *zx, mplapackint const incx, qd_complex *zy, mplapackint const incy);
void Cgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *x, mplapackint const incx, qd_complex const beta, qd_complex *y, mplapackint const incy);
void Cgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *b, mplapackint const ldb, qd_complex const beta, qd_complex *c, mplapackint const ldc);
void Cgemv(const char *trans, mplapackint const m, mplapackint const n, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *x, mplapackint const incx, qd_complex const beta, qd_complex *y, mplapackint const incy);
void Cgerc(mplapackint const m, mplapackint const n, qd_complex const alpha, qd_complex *x, mplapackint const incx, qd_complex *y, mplapackint const incy, qd_complex *a, mplapackint const lda);
void Cgeru(mplapackint const m, mplapackint const n, qd_complex const alpha, qd_complex *x, mplapackint const incx, qd_complex *y, mplapackint const incy, qd_complex *a, mplapackint const lda);
void Chbmv(const char *uplo, mplapackint const n, mplapackint const k, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *x, mplapackint const incx, qd_complex const beta, qd_complex *y, mplapackint const incy);
void Chemm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *b, mplapackint const ldb, qd_complex const beta, qd_complex *c, mplapackint const ldc);
void Chemv(const char *uplo, mplapackint const n, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *x, mplapackint const incx, qd_complex const beta, qd_complex *y, mplapackint const incy);
void Cher(const char *uplo, mplapackint const n, qd_real const alpha, qd_complex *x, mplapackint const incx, qd_complex *a, mplapackint const lda);
void Cher2(const char *uplo, mplapackint const n, qd_complex const alpha, qd_complex *x, mplapackint const incx, qd_complex *y, mplapackint const incy, qd_complex *a, mplapackint const lda);
void Cher2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *b, mplapackint const ldb, qd_real const beta, qd_complex *c, mplapackint const ldc);
void Cherk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, qd_real const alpha, qd_complex *a, mplapackint const lda, qd_real const beta, qd_complex *c, mplapackint const ldc);
void Chpmv(const char *uplo, mplapackint const n, qd_complex const alpha, qd_complex *ap, qd_complex *x, mplapackint const incx, qd_complex const beta, qd_complex *y, mplapackint const incy);
void Chpr(const char *uplo, mplapackint const n, qd_real const alpha, qd_complex *x, mplapackint const incx, qd_complex *ap);
void Chpr2(const char *uplo, mplapackint const n, qd_complex const alpha, qd_complex *x, mplapackint const incx, qd_complex *y, mplapackint const incy, qd_complex *ap);
void Crotg(qd_complex &ca, qd_complex const cb, qd_real &c, qd_complex &s);
void Cscal(mplapackint const n, qd_complex const za, qd_complex *zx, mplapackint const incx);
void Cswap(mplapackint const n, qd_complex *zx, mplapackint const incx, qd_complex *zy, mplapackint const incy);
void Csymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *b, mplapackint const ldb, qd_complex const beta, qd_complex *c, mplapackint const ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *b, mplapackint const ldb, qd_complex const beta, qd_complex *c, mplapackint const ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex const beta, qd_complex *c, mplapackint const ldc);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, qd_complex *a, mplapackint const lda, qd_complex *x, mplapackint const incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, qd_complex *a, mplapackint const lda, qd_complex *x, mplapackint const incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, qd_complex *ap, qd_complex *x, mplapackint const incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, qd_complex *ap, qd_complex *x, mplapackint const incx);
void Ctrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *b, mplapackint const ldb);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, qd_complex *a, mplapackint const lda, qd_complex *x, mplapackint const incx);
void Ctrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, qd_complex const alpha, qd_complex *a, mplapackint const lda, qd_complex *b, mplapackint const ldb);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, qd_complex *a, mplapackint const lda, qd_complex *x, mplapackint const incx);
void Mxerbla_qd(const char *srname, int info);
void Raxpy(mplapackint const n, qd_real const da, qd_real *dx, mplapackint const incx, qd_real *dy, mplapackint const incy);
void Rcopy(mplapackint const n, qd_real *dx, mplapackint const incx, qd_real *dy, mplapackint const incy);
void Rgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, qd_real const alpha, qd_real *a, mplapackint const lda, qd_real *x, mplapackint const incx, qd_real const beta, qd_real *y, mplapackint const incy);
void Rgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, qd_real const alpha, qd_real *a, mplapackint const lda, qd_real *b, mplapackint const ldb, qd_real const beta, qd_real *c, mplapackint const ldc);
void Rgemv(const char *trans, mplapackint const m, mplapackint const n, qd_real const alpha, qd_real *a, mplapackint const lda, qd_real *x, mplapackint const incx, qd_real const beta, qd_real *y, mplapackint const incy);
void Rger(mplapackint const m, mplapackint const n, qd_real const alpha, qd_real *x, mplapackint const incx, qd_real *y, mplapackint const incy, qd_real *a, mplapackint const lda);
void Rrot(mplapackint const n, qd_real *dx, mplapackint const incx, qd_real *dy, mplapackint const incy, qd_real const c, qd_real const s);
void Rrotg(qd_real &da, qd_real &db, qd_real &c, qd_real &s);
void Rrotm(mplapackint const n, qd_real *dx, mplapackint const incx, qd_real *dy, mplapackint const incy, qd_real *dparam);
void Rrotmg(qd_real &dd1, qd_real &dd2, qd_real &dx1, qd_real const dy1, qd_real *dparam);
void Rsbmv(const char *uplo, mplapackint const n, mplapackint const k, qd_real const alpha, qd_real *a, mplapackint const lda, qd_real *x, mplapackint const incx, qd_real const beta, qd_real *y, mplapackint const incy);
void Rscal(mplapackint const n, qd_real const da, qd_real *dx, mplapackint const incx);
void Rspmv(const char *uplo, mplapackint const n, qd_real const alpha, qd_real *ap, qd_real *x, mplapackint const incx, qd_real const beta, qd_real *y, mplapackint const incy);
void Rspr(const char *uplo, mplapackint const n, qd_real const alpha, qd_real *x, mplapackint const incx, qd_real *ap);
void Rspr2(const char *uplo, mplapackint const n, qd_real const alpha, qd_real *x, mplapackint const incx, qd_real *y, mplapackint const incy, qd_real *ap);
void Rswap(mplapackint const n, qd_real *dx, mplapackint const incx, qd_real *dy, mplapackint const incy);
void Rsymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, qd_real const alpha, qd_real *a, mplapackint const lda, qd_real *b, mplapackint const ldb, qd_real const beta, qd_real *c, mplapackint const ldc);
void Rsymv(const char *uplo, mplapackint const n, qd_real const alpha, qd_real *a, mplapackint const lda, qd_real *x, mplapackint const incx, qd_real const beta, qd_real *y, mplapackint const incy);
void Rsyr(const char *uplo, mplapackint const n, qd_real const alpha, qd_real *x, mplapackint const incx, qd_real *a, mplapackint const lda);
void Rsyr2(const char *uplo, mplapackint const n, qd_real const alpha, qd_real *x, mplapackint const incx, qd_real *y, mplapackint const incy, qd_real *a, mplapackint const lda);
void Rsyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, qd_real const alpha, qd_real *a, mplapackint const lda, qd_real *b, mplapackint const ldb, qd_real const beta, qd_real *c, mplapackint const ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, qd_real const alpha, qd_real *a, mplapackint const lda, qd_real const beta, qd_real *c, mplapackint const ldc);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, qd_real *a, mplapackint const lda, qd_real *x, mplapackint const incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, qd_real *a, mplapackint const lda, qd_real *x, mplapackint const incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, qd_real *ap, qd_real *x, mplapackint const incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, qd_real *ap, qd_real *x, mplapackint const incx);
void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, qd_real const alpha, qd_real *a, mplapackint const lda, qd_real *b, mplapackint const ldb);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, qd_real *a, mplapackint const lda, qd_real *x, mplapackint const incx);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, qd_real const alpha, qd_real *a, mplapackint const lda, qd_real *b, mplapackint const ldb);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, qd_real *a, mplapackint const lda, qd_real *x, mplapackint const incx);
#endif
