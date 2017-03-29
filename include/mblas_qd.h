/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mblas_qd.h,v 1.12 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MBLAS_QD_H_
#define _MBLAS_QD_H_

#include <mpack_config.h>

#include <qd/qd_real.h>

#include <qd_complex.h>
#include <mutils_qd.h>

#if !defined __MPACK_ERRNO__
#define _MPACK_EXTERN_ extern
#else
#define _MPACK_EXTERN_
#endif
_MPACK_EXTERN_ int mpack_errno;

/* LEVEL 1 MBLAS */
void Rrot(mpackint n, qd_real * dx, mpackint incx, qd_real * dy, mpackint incy, qd_real c, qd_real s);
void Rrotm(mpackint n, qd_real * x, mpackint incx, qd_real * y, mpackint incy, qd_real * dparam);
void CRrot(mpackint n, qd_complex * cx, mpackint incx, qd_complex * cy, mpackint incy, qd_real c, qd_real s);
void Cswap(mpackint n, qd_complex * cx, mpackint incx, qd_complex * cy, mpackint incy);
void Rswap(mpackint n, qd_real * dx, mpackint incx, qd_real * dy, mpackint incy);
void Cscal(mpackint n, qd_complex ca, qd_complex * cx, mpackint incx);
void CRscal(mpackint n, qd_real sa, qd_complex * cx, mpackint incx);
void Rscal(mpackint n, qd_real ca, qd_real * cx, mpackint incx);
void Ccopy(mpackint n, qd_complex * cx, mpackint incx, qd_complex * cy, mpackint incy);
void Crotg(qd_complex * ca, qd_complex cb, qd_real * c, qd_complex * s);
void Rrotg(qd_real * da, qd_real * db, qd_real * c, qd_real * s);
void Rcopy(mpackint n, qd_real * dx, mpackint incx, qd_real * dy, mpackint incy);
void Raxpy(mpackint n, qd_real da, qd_real * dx, mpackint incx, qd_real * dy, mpackint incy);
void Caxpy(mpackint n, qd_complex ca, qd_complex * cx, mpackint incx, qd_complex * cy, mpackint incy);
qd_real Rdot(mpackint n, qd_real * dx, mpackint incx, qd_real * dy, mpackint incy);
qd_complex Cdotc(mpackint n, qd_complex * cx, mpackint incx, qd_complex * cy, mpackint incy);
qd_complex Cdotu(mpackint n, qd_complex * cx, mpackint incx, qd_complex * cy, mpackint incy);
qd_real RCnrm2(mpackint n, qd_complex * x, mpackint incx);
qd_real Rnrm2(mpackint n, qd_real * x, mpackint incx);
qd_real RCasum(mpackint n, qd_complex * zx, mpackint incx);
qd_real Rasum(mpackint n, qd_real * dx, mpackint incx);
mpackint iCamax(mpackint n, qd_complex * dx, mpackint incx);
mpackint iRamax(mpackint n, qd_real * dx, mpackint incx);
qd_real RCabs1(qd_complex z);
mpackint Mlsame_qd(const char *a, const char *b);
void Mxerbla_qd(const char *srname, int info);

/* LEVEL 2 MBLAS */
void Cgemv(const char *trans, mpackint m, mpackint n, qd_complex alpha,
	   qd_complex * A, mpackint lda, qd_complex * x, mpackint incx, qd_complex beta, qd_complex * y, mpackint incy);
void Rgemv(const char *trans, mpackint m, mpackint n, qd_real alpha,
	   qd_real * A, mpackint lda, qd_real * x, mpackint incx, qd_real beta, qd_real * y, mpackint incy);
void Cgbmv(const char *trans, mpackint m, mpackint n, mpackint kl, mpackint ku,
	   qd_complex alpha, qd_complex * A, mpackint lda, qd_complex * x, mpackint incx, qd_complex beta, qd_complex * y, mpackint incy);
void Rgbmv(const char *trans, mpackint m, mpackint n, mpackint kl, mpackint ku,
	   qd_real alpha, qd_real * A, mpackint lda, qd_real * x, mpackint incx, qd_real beta, qd_real * y, mpackint incy);
void Chemv(const char *uplo, mpackint n, qd_complex alpha, qd_complex * A,
	   mpackint lda, qd_complex * x, mpackint incx, qd_complex beta, qd_complex * y, mpackint incy);
void Chbmv(const char *uplo, mpackint n, mpackint k, qd_complex alpha,
	   qd_complex * A, mpackint lda, qd_complex * x, mpackint incx, qd_complex beta, qd_complex * y, mpackint incy);
void Chpmv(const char *uplo, mpackint n, qd_complex alpha, qd_complex * AP,
	   qd_complex * x, mpackint incx, qd_complex beta, qd_complex * y, mpackint incy);
void Rsymv(const char *uplo, mpackint n, qd_real alpha, qd_real * A,
	   mpackint lda, qd_real * x, mpackint incx, qd_real beta, qd_real * y, mpackint incy);
void Rsbmv(const char *uplo, mpackint n, mpackint k, qd_real alpha,
	   qd_real * A, mpackint lda, qd_real * x, mpackint incx, qd_real beta, qd_real * y, mpackint incy);
void Rspmv(const char *uplo, mpackint n, qd_real alpha, qd_real * AP, qd_real * x, mpackint incx, qd_real beta, qd_real * y, mpackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mpackint n, qd_complex * A, mpackint lda, qd_complex * x, mpackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mpackint n, qd_real * A, mpackint lda, qd_real * x, mpackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, qd_complex * A, mpackint lda, qd_complex * x, mpackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, qd_real * A, mpackint lda, qd_real * x, mpackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mpackint n, qd_complex * AP, qd_complex * x, mpackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mpackint n, qd_real * AP, qd_real * x, mpackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mpackint n, qd_complex * A, mpackint lda, qd_complex * x, mpackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mpackint n, qd_real * A, mpackint lda, qd_real * x, mpackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, qd_complex * A, mpackint lda, qd_complex * x, mpackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, qd_real * A, mpackint lda, qd_real * x, mpackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mpackint n, qd_complex * AP, qd_complex * x, mpackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mpackint n, qd_real * AP, qd_real * x, mpackint incx);
void Rger(mpackint m, mpackint n, qd_real alpha, qd_real * x, mpackint incx, qd_real * y, mpackint incy, qd_real * A, mpackint lda);
void Cgeru(mpackint m, mpackint n, qd_complex alpha, qd_complex * x, mpackint incx, qd_complex * y, mpackint incy, qd_complex * A, mpackint lda);
void Cgerc(mpackint m, mpackint n, qd_complex alpha, qd_complex * x, mpackint incx, qd_complex * y, mpackint incy, qd_complex * A, mpackint lda);
void Cher(const char *uplo, mpackint n, qd_real alpha, qd_complex * x, mpackint incx, qd_complex * A, mpackint lda);
void Chpr(const char *uplo, mpackint n, qd_real alpha, qd_complex * x, mpackint incx, qd_complex * AP);
void Cher2(const char *uplo, mpackint n, qd_complex alpha, qd_complex * x,
	   mpackint incx, qd_complex * y, mpackint incy, qd_complex * A, mpackint lda);
void Chpr2(const char *uplo, mpackint n, qd_complex alpha, qd_complex * x, mpackint incx, qd_complex * y, mpackint incy, qd_complex * AP);
void Rsyr(const char *uplo, mpackint n, qd_real alpha, qd_real * x, mpackint incx, qd_real * A, mpackint lda);
void Rspr(const char *uplo, mpackint n, qd_real alpha, qd_real * x, mpackint incx, qd_real * AP);
void Rsyr2(const char *uplo, mpackint n, qd_real alpha, qd_real * x, mpackint incx, qd_real * y, mpackint incy, qd_real * A, mpackint lda);
void Rspr2(const char *uplo, mpackint n, qd_real alpha, qd_real * x, mpackint incx, qd_real * y, mpackint incy, qd_real * AP);

/* LEVEL 3 MBLAS */
void Cgemm(const char *transa, const char *transb, mpackint m, mpackint n,
	   mpackint k, qd_complex alpha, qd_complex * A, mpackint lda, qd_complex * B, mpackint ldb, qd_complex beta, qd_complex * C, mpackint ldc);
void Rgemm(const char *transa, const char *transb, mpackint m, mpackint n,
	   mpackint k, qd_real alpha, qd_real * A, mpackint lda, qd_real * B, mpackint ldb, qd_real beta, qd_real * C, mpackint ldc);
void Csymm(const char *side, const char *uplo, mpackint m, mpackint n,
	   qd_complex alpha, qd_complex * A, mpackint lda, qd_complex * B, mpackint ldb, qd_complex beta, qd_complex * C, mpackint ldc);
void Rsymm(const char *side, const char *uplo, mpackint m, mpackint n,
	   qd_real alpha, qd_real * A, mpackint lda, qd_real * B, mpackint ldb, qd_real beta, qd_real * C, mpackint ldc);
void Chemm(const char *side, const char *uplo, mpackint m, mpackint n,
	   qd_complex alpha, qd_complex * A, mpackint lda, qd_complex * B, mpackint ldb, qd_complex beta, qd_complex * C, mpackint ldc);
void Csyrk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   qd_complex alpha, qd_complex * A, mpackint lda, qd_complex beta, qd_complex * C, mpackint ldc);
void Rsyrk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   qd_real alpha, qd_real * A, mpackint lda, qd_real beta, qd_real * C, mpackint ldc);
void Cherk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   qd_real alpha, qd_complex * A, mpackint lda, qd_real beta, qd_complex * C, mpackint ldc);
void Csyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    qd_complex alpha, qd_complex * A, mpackint lda, qd_complex * B, mpackint ldb, qd_complex beta, qd_complex * C, mpackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    qd_real alpha, qd_real * A, mpackint lda, qd_real * B, mpackint ldb, qd_real beta, qd_real * C, mpackint ldc);
void Cher2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    qd_complex alpha, qd_complex * A, mpackint lda, qd_complex * B, mpackint ldb, qd_real beta, qd_complex * C, mpackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, qd_complex alpha, qd_complex * A, mpackint lda, qd_complex * B, mpackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, qd_real alpha, qd_real * A, mpackint lda, qd_real * B, mpackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, qd_complex alpha, qd_complex * A, mpackint lda, qd_complex * B, mpackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, qd_real alpha, qd_real * A, mpackint lda, qd_real * B, mpackint ldb);

#endif
