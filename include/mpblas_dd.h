/*
 * Copyright (c) 2008-2010
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

#include "mplapack_config.h"
#include "qd/dd_real.h"

#include "dd_complex.h"
#include "mutils_dd.h"

#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

/* LEVEL 1 MPBLAS */
void Rrot(mplapackint n, dd_real * dx, mplapackint incx, dd_real * dy, mplapackint incy, dd_real c, dd_real s);
void Rrotm(mplapackint n, dd_real * x, mplapackint incx, dd_real * y, mplapackint incy, dd_real * dparam);
void CRrot(mplapackint n, dd_complex * cx, mplapackint incx, dd_complex * cy, mplapackint incy, dd_real c, dd_real s);
void Cswap(mplapackint n, dd_complex * cx, mplapackint incx, dd_complex * cy, mplapackint incy);
void Rswap(mplapackint n, dd_real * dx, mplapackint incx, dd_real * dy, mplapackint incy);
void Cscal(mplapackint n, dd_complex ca, dd_complex * cx, mplapackint incx);
void CRscal(mplapackint n, dd_real sa, dd_complex * cx, mplapackint incx);
void Rscal(mplapackint n, dd_real ca, dd_real * cx, mplapackint incx);
void Ccopy(mplapackint n, dd_complex * cx, mplapackint incx, dd_complex * cy, mplapackint incy);
void Crotg(dd_complex * ca, dd_complex cb, dd_real * c, dd_complex * s);
void Rrotg(dd_real * da, dd_real * db, dd_real * c, dd_real * s);
void Rcopy(mplapackint n, dd_real * dx, mplapackint incx, dd_real * dy, mplapackint incy);
void Raxpy(mplapackint n, dd_real da, dd_real * dx, mplapackint incx, dd_real * dy, mplapackint incy);
void Caxpy(mplapackint n, dd_complex ca, dd_complex * cx, mplapackint incx, dd_complex * cy, mplapackint incy);
dd_real Rdot(mplapackint n, dd_real * dx, mplapackint incx, dd_real * dy, mplapackint incy);
dd_complex Cdotc(mplapackint n, dd_complex * cx, mplapackint incx, dd_complex * cy, mplapackint incy);
dd_complex Cdotu(mplapackint n, dd_complex * cx, mplapackint incx, dd_complex * cy, mplapackint incy);
dd_real RCnrm2(mplapackint n, dd_complex * x, mplapackint incx);
dd_real Rnrm2(mplapackint n, dd_real * x, mplapackint incx);
dd_real RCasum(mplapackint n, dd_complex * zx, mplapackint incx);
dd_real Rasum(mplapackint n, dd_real * dx, mplapackint incx);
mplapackint iCamax(mplapackint n, dd_complex * dx, mplapackint incx);
mplapackint iRamax(mplapackint n, dd_real * dx, mplapackint incx);
dd_real RCabs1(dd_complex z);
mplapackint Mlsame_dd(const char *a, const char *b);
void Mxerbla_dd(const char *srname, int info);

/* LEVEL 2 MPBLAS */
void Cgemv(const char *trans, mplapackint m, mplapackint n, dd_complex alpha,
	   dd_complex * A, mplapackint lda, dd_complex * x, mplapackint incx, dd_complex beta, dd_complex * y, mplapackint incy);
void Rgemv(const char *trans, mplapackint m, mplapackint n, dd_real alpha,
	   dd_real * A, mplapackint lda, dd_real * x, mplapackint incx, dd_real beta, dd_real * y, mplapackint incy);
void Cgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   dd_complex alpha, dd_complex * A, mplapackint lda, dd_complex * x, mplapackint incx, dd_complex beta, dd_complex * y, mplapackint incy);
void Rgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   dd_real alpha, dd_real * A, mplapackint lda, dd_real * x, mplapackint incx, dd_real beta, dd_real * y, mplapackint incy);
void Chemv(const char *uplo, mplapackint n, dd_complex alpha, dd_complex * A,
	   mplapackint lda, dd_complex * x, mplapackint incx, dd_complex beta, dd_complex * y, mplapackint incy);
void Chbmv(const char *uplo, mplapackint n, mplapackint k, dd_complex alpha,
	   dd_complex * A, mplapackint lda, dd_complex * x, mplapackint incx, dd_complex beta, dd_complex * y, mplapackint incy);
void Chpmv(const char *uplo, mplapackint n, dd_complex alpha, dd_complex * AP,
	   dd_complex * x, mplapackint incx, dd_complex beta, dd_complex * y, mplapackint incy);
void Rsymv(const char *uplo, mplapackint n, dd_real alpha, dd_real * A,
	   mplapackint lda, dd_real * x, mplapackint incx, dd_real beta, dd_real * y, mplapackint incy);
void Rsbmv(const char *uplo, mplapackint n, mplapackint k, dd_real alpha,
	   dd_real * A, mplapackint lda, dd_real * x, mplapackint incx, dd_real beta, dd_real * y, mplapackint incy);
void Rspmv(const char *uplo, mplapackint n, dd_real alpha, dd_real * AP, dd_real * x, mplapackint incx, dd_real beta, dd_real * y, mplapackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, dd_complex * A, mplapackint lda, dd_complex * x, mplapackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, dd_real * A, mplapackint lda, dd_real * x, mplapackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, dd_complex * A, mplapackint lda, dd_complex * x, mplapackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, dd_real * A, mplapackint lda, dd_real * x, mplapackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, dd_complex * AP, dd_complex * x, mplapackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, dd_real * AP, dd_real * x, mplapackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, dd_complex * A, mplapackint lda, dd_complex * x, mplapackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, dd_real * A, mplapackint lda, dd_real * x, mplapackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, dd_complex * A, mplapackint lda, dd_complex * x, mplapackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, dd_real * A, mplapackint lda, dd_real * x, mplapackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, dd_complex * AP, dd_complex * x, mplapackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, dd_real * AP, dd_real * x, mplapackint incx);
void Rger(mplapackint m, mplapackint n, dd_real alpha, dd_real * x, mplapackint incx, dd_real * y, mplapackint incy, dd_real * A, mplapackint lda);
void Cgeru(mplapackint m, mplapackint n, dd_complex alpha, dd_complex * x, mplapackint incx, dd_complex * y, mplapackint incy, dd_complex * A, mplapackint lda);
void Cgerc(mplapackint m, mplapackint n, dd_complex alpha, dd_complex * x, mplapackint incx, dd_complex * y, mplapackint incy, dd_complex * A, mplapackint lda);
void Cher(const char *uplo, mplapackint n, dd_real alpha, dd_complex * x, mplapackint incx, dd_complex * A, mplapackint lda);
void Chpr(const char *uplo, mplapackint n, dd_real alpha, dd_complex * x, mplapackint incx, dd_complex * AP);
void Cher2(const char *uplo, mplapackint n, dd_complex alpha, dd_complex * x,
	   mplapackint incx, dd_complex * y, mplapackint incy, dd_complex * A, mplapackint lda);
void Chpr2(const char *uplo, mplapackint n, dd_complex alpha, dd_complex * x, mplapackint incx, dd_complex * y, mplapackint incy, dd_complex * AP);
void Rsyr(const char *uplo, mplapackint n, dd_real alpha, dd_real * x, mplapackint incx, dd_real * A, mplapackint lda);
void Rspr(const char *uplo, mplapackint n, dd_real alpha, dd_real * x, mplapackint incx, dd_real * AP);
void Rsyr2(const char *uplo, mplapackint n, dd_real alpha, dd_real * x, mplapackint incx, dd_real * y, mplapackint incy, dd_real * A, mplapackint lda);
void Rspr2(const char *uplo, mplapackint n, dd_real alpha, dd_real * x, mplapackint incx, dd_real * y, mplapackint incy, dd_real * AP);

/* LEVEL 3 MPBLAS */
void Cgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, dd_complex alpha, dd_complex * A, mplapackint lda, dd_complex * B, mplapackint ldb, dd_complex beta, dd_complex * C, mplapackint ldc);
void Rgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, dd_real alpha, dd_real * A, mplapackint lda, dd_real * B, mplapackint ldb, dd_real beta, dd_real * C, mplapackint ldc);
void Csymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   dd_complex alpha, dd_complex * A, mplapackint lda, dd_complex * B, mplapackint ldb, dd_complex beta, dd_complex * C, mplapackint ldc);
void Rsymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   dd_real alpha, dd_real * A, mplapackint lda, dd_real * B, mplapackint ldb, dd_real beta, dd_real * C, mplapackint ldc);
void Chemm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   dd_complex alpha, dd_complex * A, mplapackint lda, dd_complex * B, mplapackint ldb, dd_complex beta, dd_complex * C, mplapackint ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   dd_complex alpha, dd_complex * A, mplapackint lda, dd_complex beta, dd_complex * C, mplapackint ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   dd_real alpha, dd_real * A, mplapackint lda, dd_real beta, dd_real * C, mplapackint ldc);
void Cherk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   dd_real alpha, dd_complex * A, mplapackint lda, dd_real beta, dd_complex * C, mplapackint ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    dd_complex alpha, dd_complex * A, mplapackint lda, dd_complex * B, mplapackint ldb, dd_complex beta, dd_complex * C, mplapackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    dd_real alpha, dd_real * A, mplapackint lda, dd_real * B, mplapackint ldb, dd_real beta, dd_real * C, mplapackint ldc);
void Cher2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    dd_complex alpha, dd_complex * A, mplapackint lda, dd_complex * B, mplapackint ldb, dd_real beta, dd_complex * C, mplapackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, dd_complex alpha, dd_complex * A, mplapackint lda, dd_complex * B, mplapackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, dd_real alpha, dd_real * A, mplapackint lda, dd_real * B, mplapackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, dd_complex alpha, dd_complex * A, mplapackint lda, dd_complex * B, mplapackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, dd_real alpha, dd_real * A, mplapackint lda, dd_real * B, mplapackint ldb);

#endif
