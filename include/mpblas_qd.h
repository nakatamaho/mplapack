/*
 * Copyright (c) 2008-2010
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

#include "mplapack_config.h"

#include <qd/qd_real.h>

#include "qd_complex.h"
#include "mplapack_utils_qd.h"

#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

/* LEVEL 1 MPBLAS */
void Rrot(mplapackint n, qd_real * dx, mplapackint incx, qd_real * dy, mplapackint incy, qd_real c, qd_real s);
void Rrotm(mplapackint n, qd_real * x, mplapackint incx, qd_real * y, mplapackint incy, qd_real * dparam);
void CRrot(mplapackint n, qd_complex * cx, mplapackint incx, qd_complex * cy, mplapackint incy, qd_real c, qd_real s);
void Cswap(mplapackint n, qd_complex * cx, mplapackint incx, qd_complex * cy, mplapackint incy);
void Rswap(mplapackint n, qd_real * dx, mplapackint incx, qd_real * dy, mplapackint incy);
void Cscal(mplapackint n, qd_complex ca, qd_complex * cx, mplapackint incx);
void CRscal(mplapackint n, qd_real sa, qd_complex * cx, mplapackint incx);
void Rscal(mplapackint n, qd_real ca, qd_real * cx, mplapackint incx);
void Ccopy(mplapackint n, qd_complex * cx, mplapackint incx, qd_complex * cy, mplapackint incy);
void Crotg(qd_complex * ca, qd_complex cb, qd_real * c, qd_complex * s);
void Rrotg(qd_real * da, qd_real * db, qd_real * c, qd_real * s);
void Rcopy(mplapackint n, qd_real * dx, mplapackint incx, qd_real * dy, mplapackint incy);
void Raxpy(mplapackint n, qd_real da, qd_real * dx, mplapackint incx, qd_real * dy, mplapackint incy);
void Caxpy(mplapackint n, qd_complex ca, qd_complex * cx, mplapackint incx, qd_complex * cy, mplapackint incy);
qd_real Rdot(mplapackint n, qd_real * dx, mplapackint incx, qd_real * dy, mplapackint incy);
qd_complex Cdotc(mplapackint n, qd_complex * cx, mplapackint incx, qd_complex * cy, mplapackint incy);
qd_complex Cdotu(mplapackint n, qd_complex * cx, mplapackint incx, qd_complex * cy, mplapackint incy);
qd_real RCnrm2(mplapackint n, qd_complex * x, mplapackint incx);
qd_real Rnrm2(mplapackint n, qd_real * x, mplapackint incx);
qd_real RCasum(mplapackint n, qd_complex * zx, mplapackint incx);
qd_real Rasum(mplapackint n, qd_real * dx, mplapackint incx);
mplapackint iCamax(mplapackint n, qd_complex * dx, mplapackint incx);
mplapackint iRamax(mplapackint n, qd_real * dx, mplapackint incx);
qd_real RCabs1(qd_complex z);
mplapackint Mlsame_qd(const char *a, const char *b);
void Mxerbla_qd(const char *srname, int info);

/* LEVEL 2 MPBLAS */
void Cgemv(const char *trans, mplapackint m, mplapackint n, qd_complex alpha,
	   qd_complex * A, mplapackint lda, qd_complex * x, mplapackint incx, qd_complex beta, qd_complex * y, mplapackint incy);
void Rgemv(const char *trans, mplapackint m, mplapackint n, qd_real alpha,
	   qd_real * A, mplapackint lda, qd_real * x, mplapackint incx, qd_real beta, qd_real * y, mplapackint incy);
void Cgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   qd_complex alpha, qd_complex * A, mplapackint lda, qd_complex * x, mplapackint incx, qd_complex beta, qd_complex * y, mplapackint incy);
void Rgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   qd_real alpha, qd_real * A, mplapackint lda, qd_real * x, mplapackint incx, qd_real beta, qd_real * y, mplapackint incy);
void Chemv(const char *uplo, mplapackint n, qd_complex alpha, qd_complex * A,
	   mplapackint lda, qd_complex * x, mplapackint incx, qd_complex beta, qd_complex * y, mplapackint incy);
void Chbmv(const char *uplo, mplapackint n, mplapackint k, qd_complex alpha,
	   qd_complex * A, mplapackint lda, qd_complex * x, mplapackint incx, qd_complex beta, qd_complex * y, mplapackint incy);
void Chpmv(const char *uplo, mplapackint n, qd_complex alpha, qd_complex * AP,
	   qd_complex * x, mplapackint incx, qd_complex beta, qd_complex * y, mplapackint incy);
void Rsymv(const char *uplo, mplapackint n, qd_real alpha, qd_real * A,
	   mplapackint lda, qd_real * x, mplapackint incx, qd_real beta, qd_real * y, mplapackint incy);
void Rsbmv(const char *uplo, mplapackint n, mplapackint k, qd_real alpha,
	   qd_real * A, mplapackint lda, qd_real * x, mplapackint incx, qd_real beta, qd_real * y, mplapackint incy);
void Rspmv(const char *uplo, mplapackint n, qd_real alpha, qd_real * AP, qd_real * x, mplapackint incx, qd_real beta, qd_real * y, mplapackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, qd_complex * A, mplapackint lda, qd_complex * x, mplapackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, qd_real * A, mplapackint lda, qd_real * x, mplapackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, qd_complex * A, mplapackint lda, qd_complex * x, mplapackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, qd_real * A, mplapackint lda, qd_real * x, mplapackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, qd_complex * AP, qd_complex * x, mplapackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, qd_real * AP, qd_real * x, mplapackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, qd_complex * A, mplapackint lda, qd_complex * x, mplapackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, qd_real * A, mplapackint lda, qd_real * x, mplapackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, qd_complex * A, mplapackint lda, qd_complex * x, mplapackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, qd_real * A, mplapackint lda, qd_real * x, mplapackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, qd_complex * AP, qd_complex * x, mplapackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, qd_real * AP, qd_real * x, mplapackint incx);
void Rger(mplapackint m, mplapackint n, qd_real alpha, qd_real * x, mplapackint incx, qd_real * y, mplapackint incy, qd_real * A, mplapackint lda);
void Cgeru(mplapackint m, mplapackint n, qd_complex alpha, qd_complex * x, mplapackint incx, qd_complex * y, mplapackint incy, qd_complex * A, mplapackint lda);
void Cgerc(mplapackint m, mplapackint n, qd_complex alpha, qd_complex * x, mplapackint incx, qd_complex * y, mplapackint incy, qd_complex * A, mplapackint lda);
void Cher(const char *uplo, mplapackint n, qd_real alpha, qd_complex * x, mplapackint incx, qd_complex * A, mplapackint lda);
void Chpr(const char *uplo, mplapackint n, qd_real alpha, qd_complex * x, mplapackint incx, qd_complex * AP);
void Cher2(const char *uplo, mplapackint n, qd_complex alpha, qd_complex * x,
	   mplapackint incx, qd_complex * y, mplapackint incy, qd_complex * A, mplapackint lda);
void Chpr2(const char *uplo, mplapackint n, qd_complex alpha, qd_complex * x, mplapackint incx, qd_complex * y, mplapackint incy, qd_complex * AP);
void Rsyr(const char *uplo, mplapackint n, qd_real alpha, qd_real * x, mplapackint incx, qd_real * A, mplapackint lda);
void Rspr(const char *uplo, mplapackint n, qd_real alpha, qd_real * x, mplapackint incx, qd_real * AP);
void Rsyr2(const char *uplo, mplapackint n, qd_real alpha, qd_real * x, mplapackint incx, qd_real * y, mplapackint incy, qd_real * A, mplapackint lda);
void Rspr2(const char *uplo, mplapackint n, qd_real alpha, qd_real * x, mplapackint incx, qd_real * y, mplapackint incy, qd_real * AP);

/* LEVEL 3 MPBLAS */
void Cgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, qd_complex alpha, qd_complex * A, mplapackint lda, qd_complex * B, mplapackint ldb, qd_complex beta, qd_complex * C, mplapackint ldc);
void Rgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, qd_real alpha, qd_real * A, mplapackint lda, qd_real * B, mplapackint ldb, qd_real beta, qd_real * C, mplapackint ldc);
void Csymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   qd_complex alpha, qd_complex * A, mplapackint lda, qd_complex * B, mplapackint ldb, qd_complex beta, qd_complex * C, mplapackint ldc);
void Rsymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   qd_real alpha, qd_real * A, mplapackint lda, qd_real * B, mplapackint ldb, qd_real beta, qd_real * C, mplapackint ldc);
void Chemm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   qd_complex alpha, qd_complex * A, mplapackint lda, qd_complex * B, mplapackint ldb, qd_complex beta, qd_complex * C, mplapackint ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   qd_complex alpha, qd_complex * A, mplapackint lda, qd_complex beta, qd_complex * C, mplapackint ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   qd_real alpha, qd_real * A, mplapackint lda, qd_real beta, qd_real * C, mplapackint ldc);
void Cherk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   qd_real alpha, qd_complex * A, mplapackint lda, qd_real beta, qd_complex * C, mplapackint ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    qd_complex alpha, qd_complex * A, mplapackint lda, qd_complex * B, mplapackint ldb, qd_complex beta, qd_complex * C, mplapackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    qd_real alpha, qd_real * A, mplapackint lda, qd_real * B, mplapackint ldb, qd_real beta, qd_real * C, mplapackint ldc);
void Cher2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    qd_complex alpha, qd_complex * A, mplapackint lda, qd_complex * B, mplapackint ldb, qd_real beta, qd_complex * C, mplapackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, qd_complex alpha, qd_complex * A, mplapackint lda, qd_complex * B, mplapackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, qd_real alpha, qd_real * A, mplapackint lda, qd_real * B, mplapackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, qd_complex alpha, qd_complex * A, mplapackint lda, qd_complex * B, mplapackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, qd_real alpha, qd_real * A, mplapackint lda, qd_real * B, mplapackint ldb);

#endif
