/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpblas_mpfr.h,v 1.4 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPBLAS_MPFR_H_
#define _MPBLAS_MPFR_H_

#include "mplapack_config.h"
#include "mpreal.h"
#include "mpcomplex.h"
using namespace mpfr;
#include "mplapack_utils_mpfr.h"
#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

/* LEVEL 1 MPBLAS */
void Rrot(mplapackint n, mpreal * dx, mplapackint incx, mpreal * dy, mplapackint incy, mpreal c, mpreal s);
void Rrotm(mplapackint n, mpreal * x, mplapackint incx, mpreal * y, mplapackint incy, mpreal * dparam);
void CRrot(mplapackint n, mpcomplex * cx, mplapackint incx, mpcomplex * cy, mplapackint incy, mpreal c, mpreal s);
void Cswap(mplapackint n, mpcomplex * cx, mplapackint incx, mpcomplex * cy, mplapackint incy);
void Rswap(mplapackint n, mpreal * dx, mplapackint incx, mpreal * dy, mplapackint incy);
void Cscal(mplapackint n, mpcomplex ca, mpcomplex * cx, mplapackint incx);
void CRscal(mplapackint n, mpreal sa, mpcomplex * cx, mplapackint incx);
void Rscal(mplapackint n, mpreal ca, mpreal * cx, mplapackint incx);
void Ccopy(mplapackint n, mpcomplex * cx, mplapackint incx, mpcomplex * cy, mplapackint incy);
void Crotg(mpcomplex * ca, mpcomplex cb, mpreal * c, mpcomplex * s);
void Rrotg(mpreal * da, mpreal * db, mpreal * c, mpreal * s);
void Rcopy(mplapackint n, mpreal * dx, mplapackint incx, mpreal * dy, mplapackint incy);
void Raxpy(mplapackint n, mpreal da, mpreal * dx, mplapackint incx, mpreal * dy, mplapackint incy);
void Caxpy(mplapackint n, mpcomplex ca, mpcomplex * cx, mplapackint incx, mpcomplex * cy, mplapackint incy);
mpreal Rdot(mplapackint n, mpreal * dx, mplapackint incx, mpreal * dy, mplapackint incy);
mpcomplex Cdotc(mplapackint n, mpcomplex * cx, mplapackint incx, mpcomplex * cy, mplapackint incy);
mpcomplex Cdotu(mplapackint n, mpcomplex * cx, mplapackint incx, mpcomplex * cy, mplapackint incy);
mpreal RCnrm2(mplapackint n, mpcomplex * x, mplapackint incx);
mpreal Rnrm2(mplapackint n, mpreal * x, mplapackint incx);
mpreal RCasum(mplapackint n, mpcomplex * zx, mplapackint incx);
mpreal Rasum(mplapackint n, mpreal * dx, mplapackint incx);
mplapackint iCamax(mplapackint n, mpcomplex * dx, mplapackint incx);
mplapackint iRamax(mplapackint n, mpreal * dx, mplapackint incx);
mpreal RCabs1(mpcomplex z);
mplapackint Mlsame_mpfr(const char *a, const char *b);
void Mxerbla_mpfr(const char *srname, int info);

/* LEVEL 2 MPBLAS */
void Cgemv(const char *trans, mplapackint m, mplapackint n, mpcomplex alpha,
	   mpcomplex * A, mplapackint lda, mpcomplex * x, mplapackint incx, mpcomplex beta, mpcomplex * y, mplapackint incy);
void Rgemv(const char *trans, mplapackint m, mplapackint n, mpreal alpha,
	   mpreal * A, mplapackint lda, mpreal * x, mplapackint incx, mpreal beta, mpreal * y, mplapackint incy);
void Cgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   mpcomplex alpha, mpcomplex * A, mplapackint lda, mpcomplex * x, mplapackint incx, mpcomplex beta, mpcomplex * y, mplapackint incy);
void Rgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   mpreal alpha, mpreal * A, mplapackint lda, mpreal * x, mplapackint incx, mpreal beta, mpreal * y, mplapackint incy);
void Chemv(const char *uplo, mplapackint n, mpcomplex alpha, mpcomplex * A,
	   mplapackint lda, mpcomplex * x, mplapackint incx, mpcomplex beta, mpcomplex * y, mplapackint incy);
void Chbmv(const char *uplo, mplapackint n, mplapackint k, mpcomplex alpha,
	   mpcomplex * A, mplapackint lda, mpcomplex * x, mplapackint incx, mpcomplex beta, mpcomplex * y, mplapackint incy);
void Chpmv(const char *uplo, mplapackint n, mpcomplex alpha, mpcomplex * AP,
	   mpcomplex * x, mplapackint incx, mpcomplex beta, mpcomplex * y, mplapackint incy);
void Rsymv(const char *uplo, mplapackint n, mpreal alpha, mpreal * A,
	   mplapackint lda, mpreal * x, mplapackint incx, mpreal beta, mpreal * y, mplapackint incy);
void Rsbmv(const char *uplo, mplapackint n, mplapackint k, mpreal alpha,
	   mpreal * A, mplapackint lda, mpreal * x, mplapackint incx, mpreal beta, mpreal * y, mplapackint incy);
void Rspmv(const char *uplo, mplapackint n, mpreal alpha, mpreal * AP, mpreal * x, mplapackint incx, mpreal beta, mpreal * y, mplapackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpcomplex * A, mplapackint lda, mpcomplex * x, mplapackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpreal * A, mplapackint lda, mpreal * x, mplapackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, mpcomplex * A, mplapackint lda, mpcomplex * x, mplapackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, mpreal * A, mplapackint lda, mpreal * x, mplapackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpcomplex * AP, mpcomplex * x, mplapackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpreal * AP, mpreal * x, mplapackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpcomplex * A, mplapackint lda, mpcomplex * x, mplapackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpreal * A, mplapackint lda, mpreal * x, mplapackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, mpcomplex * A, mplapackint lda, mpcomplex * x, mplapackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, mpreal * A, mplapackint lda, mpreal * x, mplapackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpcomplex * AP, mpcomplex * x, mplapackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpreal * AP, mpreal * x, mplapackint incx);
void Rger(mplapackint m, mplapackint n, mpreal alpha, mpreal * x, mplapackint incx, mpreal * y, mplapackint incy, mpreal * A, mplapackint lda);
void Cgeru(mplapackint m, mplapackint n, mpcomplex alpha, mpcomplex * x, mplapackint incx, mpcomplex * y, mplapackint incy, mpcomplex * A, mplapackint lda);
void Cgerc(mplapackint m, mplapackint n, mpcomplex alpha, mpcomplex * x, mplapackint incx, mpcomplex * y, mplapackint incy, mpcomplex * A, mplapackint lda);
void Cher(const char *uplo, mplapackint n, mpreal alpha, mpcomplex * x, mplapackint incx, mpcomplex * A, mplapackint lda);
void Chpr(const char *uplo, mplapackint n, mpreal alpha, mpcomplex * x, mplapackint incx, mpcomplex * AP);
void Cher2(const char *uplo, mplapackint n, mpcomplex alpha, mpcomplex * x,
	   mplapackint incx, mpcomplex * y, mplapackint incy, mpcomplex * A, mplapackint lda);
void Chpr2(const char *uplo, mplapackint n, mpcomplex alpha, mpcomplex * x, mplapackint incx, mpcomplex * y, mplapackint incy, mpcomplex * AP);
void Rsyr(const char *uplo, mplapackint n, mpreal alpha, mpreal * x, mplapackint incx, mpreal * A, mplapackint lda);
void Rspr(const char *uplo, mplapackint n, mpreal alpha, mpreal * x, mplapackint incx, mpreal * AP);
void Rsyr2(const char *uplo, mplapackint n, mpreal alpha, mpreal * x, mplapackint incx, mpreal * y, mplapackint incy, mpreal * A, mplapackint lda);
void Rspr2(const char *uplo, mplapackint n, mpreal alpha, mpreal * x, mplapackint incx, mpreal * y, mplapackint incy, mpreal * AP);

/* LEVEL 3 MPBLAS */
void Cgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, mpcomplex alpha, mpcomplex * A, mplapackint lda, mpcomplex * B, mplapackint ldb, mpcomplex beta, mpcomplex * C, mplapackint ldc);
void Rgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, mpreal alpha, mpreal * A, mplapackint lda, mpreal * B, mplapackint ldb, mpreal beta, mpreal * C, mplapackint ldc);
void Csymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   mpcomplex alpha, mpcomplex * A, mplapackint lda, mpcomplex * B, mplapackint ldb, mpcomplex beta, mpcomplex * C, mplapackint ldc);
void Rsymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   mpreal alpha, mpreal * A, mplapackint lda, mpreal * B, mplapackint ldb, mpreal beta, mpreal * C, mplapackint ldc);
void Chemm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   mpcomplex alpha, mpcomplex * A, mplapackint lda, mpcomplex * B, mplapackint ldb, mpcomplex beta, mpcomplex * C, mplapackint ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   mpcomplex alpha, mpcomplex * A, mplapackint lda, mpcomplex beta, mpcomplex * C, mplapackint ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   mpreal alpha, mpreal * A, mplapackint lda, mpreal beta, mpreal * C, mplapackint ldc);
void Cherk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   mpreal alpha, mpcomplex * A, mplapackint lda, mpreal beta, mpcomplex * C, mplapackint ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    mpcomplex alpha, mpcomplex * A, mplapackint lda, mpcomplex * B, mplapackint ldb, mpcomplex beta, mpcomplex * C, mplapackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    mpreal alpha, mpreal * A, mplapackint lda, mpreal * B, mplapackint ldb, mpreal beta, mpreal * C, mplapackint ldc);
void Cher2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    mpcomplex alpha, mpcomplex * A, mplapackint lda, mpcomplex * B, mplapackint ldb, mpreal beta, mpcomplex * C, mplapackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, mpcomplex alpha, mpcomplex * A, mplapackint lda, mpcomplex * B, mplapackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, mpreal alpha, mpreal * A, mplapackint lda, mpreal * B, mplapackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, mpcomplex alpha, mpcomplex * A, mplapackint lda, mpcomplex * B, mplapackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, mpreal alpha, mpreal * A, mplapackint lda, mpreal * B, mplapackint ldb);

#endif
