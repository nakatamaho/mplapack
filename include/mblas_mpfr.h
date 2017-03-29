/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mblas_mpfr.h,v 1.4 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MBLAS_MPFR_H_
#define _MBLAS_MPFR_H_

#include <mpack_config.h>
#include <mpreal.h>
#include <mpcomplex.h>
using namespace mpfr;
#include <mutils_mpfr.h>
#if !defined __MPACK_ERRNO__
#define _MPACK_EXTERN_ extern
#else
#define _MPACK_EXTERN_
#endif
_MPACK_EXTERN_ int mpack_errno;

/* LEVEL 1 MBLAS */
void Rrot(mpackint n, mpreal * dx, mpackint incx, mpreal * dy, mpackint incy, mpreal c, mpreal s);
void Rrotm(mpackint n, mpreal * x, mpackint incx, mpreal * y, mpackint incy, mpreal * dparam);
void CRrot(mpackint n, mpcomplex * cx, mpackint incx, mpcomplex * cy, mpackint incy, mpreal c, mpreal s);
void Cswap(mpackint n, mpcomplex * cx, mpackint incx, mpcomplex * cy, mpackint incy);
void Rswap(mpackint n, mpreal * dx, mpackint incx, mpreal * dy, mpackint incy);
void Cscal(mpackint n, mpcomplex ca, mpcomplex * cx, mpackint incx);
void CRscal(mpackint n, mpreal sa, mpcomplex * cx, mpackint incx);
void Rscal(mpackint n, mpreal ca, mpreal * cx, mpackint incx);
void Ccopy(mpackint n, mpcomplex * cx, mpackint incx, mpcomplex * cy, mpackint incy);
void Crotg(mpcomplex * ca, mpcomplex cb, mpreal * c, mpcomplex * s);
void Rrotg(mpreal * da, mpreal * db, mpreal * c, mpreal * s);
void Rcopy(mpackint n, mpreal * dx, mpackint incx, mpreal * dy, mpackint incy);
void Raxpy(mpackint n, mpreal da, mpreal * dx, mpackint incx, mpreal * dy, mpackint incy);
void Caxpy(mpackint n, mpcomplex ca, mpcomplex * cx, mpackint incx, mpcomplex * cy, mpackint incy);
mpreal Rdot(mpackint n, mpreal * dx, mpackint incx, mpreal * dy, mpackint incy);
mpcomplex Cdotc(mpackint n, mpcomplex * cx, mpackint incx, mpcomplex * cy, mpackint incy);
mpcomplex Cdotu(mpackint n, mpcomplex * cx, mpackint incx, mpcomplex * cy, mpackint incy);
mpreal RCnrm2(mpackint n, mpcomplex * x, mpackint incx);
mpreal Rnrm2(mpackint n, mpreal * x, mpackint incx);
mpreal RCasum(mpackint n, mpcomplex * zx, mpackint incx);
mpreal Rasum(mpackint n, mpreal * dx, mpackint incx);
mpackint iCamax(mpackint n, mpcomplex * dx, mpackint incx);
mpackint iRamax(mpackint n, mpreal * dx, mpackint incx);
mpreal RCabs1(mpcomplex z);
mpackint Mlsame_mpfr(const char *a, const char *b);
void Mxerbla_mpfr(const char *srname, int info);

/* LEVEL 2 MBLAS */
void Cgemv(const char *trans, mpackint m, mpackint n, mpcomplex alpha,
	   mpcomplex * A, mpackint lda, mpcomplex * x, mpackint incx, mpcomplex beta, mpcomplex * y, mpackint incy);
void Rgemv(const char *trans, mpackint m, mpackint n, mpreal alpha,
	   mpreal * A, mpackint lda, mpreal * x, mpackint incx, mpreal beta, mpreal * y, mpackint incy);
void Cgbmv(const char *trans, mpackint m, mpackint n, mpackint kl, mpackint ku,
	   mpcomplex alpha, mpcomplex * A, mpackint lda, mpcomplex * x, mpackint incx, mpcomplex beta, mpcomplex * y, mpackint incy);
void Rgbmv(const char *trans, mpackint m, mpackint n, mpackint kl, mpackint ku,
	   mpreal alpha, mpreal * A, mpackint lda, mpreal * x, mpackint incx, mpreal beta, mpreal * y, mpackint incy);
void Chemv(const char *uplo, mpackint n, mpcomplex alpha, mpcomplex * A,
	   mpackint lda, mpcomplex * x, mpackint incx, mpcomplex beta, mpcomplex * y, mpackint incy);
void Chbmv(const char *uplo, mpackint n, mpackint k, mpcomplex alpha,
	   mpcomplex * A, mpackint lda, mpcomplex * x, mpackint incx, mpcomplex beta, mpcomplex * y, mpackint incy);
void Chpmv(const char *uplo, mpackint n, mpcomplex alpha, mpcomplex * AP,
	   mpcomplex * x, mpackint incx, mpcomplex beta, mpcomplex * y, mpackint incy);
void Rsymv(const char *uplo, mpackint n, mpreal alpha, mpreal * A,
	   mpackint lda, mpreal * x, mpackint incx, mpreal beta, mpreal * y, mpackint incy);
void Rsbmv(const char *uplo, mpackint n, mpackint k, mpreal alpha,
	   mpreal * A, mpackint lda, mpreal * x, mpackint incx, mpreal beta, mpreal * y, mpackint incy);
void Rspmv(const char *uplo, mpackint n, mpreal alpha, mpreal * AP, mpreal * x, mpackint incx, mpreal beta, mpreal * y, mpackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mpackint n, mpcomplex * A, mpackint lda, mpcomplex * x, mpackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mpackint n, mpreal * A, mpackint lda, mpreal * x, mpackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, mpcomplex * A, mpackint lda, mpcomplex * x, mpackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, mpreal * A, mpackint lda, mpreal * x, mpackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mpackint n, mpcomplex * AP, mpcomplex * x, mpackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mpackint n, mpreal * AP, mpreal * x, mpackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mpackint n, mpcomplex * A, mpackint lda, mpcomplex * x, mpackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mpackint n, mpreal * A, mpackint lda, mpreal * x, mpackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, mpcomplex * A, mpackint lda, mpcomplex * x, mpackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, mpreal * A, mpackint lda, mpreal * x, mpackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mpackint n, mpcomplex * AP, mpcomplex * x, mpackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mpackint n, mpreal * AP, mpreal * x, mpackint incx);
void Rger(mpackint m, mpackint n, mpreal alpha, mpreal * x, mpackint incx, mpreal * y, mpackint incy, mpreal * A, mpackint lda);
void Cgeru(mpackint m, mpackint n, mpcomplex alpha, mpcomplex * x, mpackint incx, mpcomplex * y, mpackint incy, mpcomplex * A, mpackint lda);
void Cgerc(mpackint m, mpackint n, mpcomplex alpha, mpcomplex * x, mpackint incx, mpcomplex * y, mpackint incy, mpcomplex * A, mpackint lda);
void Cher(const char *uplo, mpackint n, mpreal alpha, mpcomplex * x, mpackint incx, mpcomplex * A, mpackint lda);
void Chpr(const char *uplo, mpackint n, mpreal alpha, mpcomplex * x, mpackint incx, mpcomplex * AP);
void Cher2(const char *uplo, mpackint n, mpcomplex alpha, mpcomplex * x,
	   mpackint incx, mpcomplex * y, mpackint incy, mpcomplex * A, mpackint lda);
void Chpr2(const char *uplo, mpackint n, mpcomplex alpha, mpcomplex * x, mpackint incx, mpcomplex * y, mpackint incy, mpcomplex * AP);
void Rsyr(const char *uplo, mpackint n, mpreal alpha, mpreal * x, mpackint incx, mpreal * A, mpackint lda);
void Rspr(const char *uplo, mpackint n, mpreal alpha, mpreal * x, mpackint incx, mpreal * AP);
void Rsyr2(const char *uplo, mpackint n, mpreal alpha, mpreal * x, mpackint incx, mpreal * y, mpackint incy, mpreal * A, mpackint lda);
void Rspr2(const char *uplo, mpackint n, mpreal alpha, mpreal * x, mpackint incx, mpreal * y, mpackint incy, mpreal * AP);

/* LEVEL 3 MBLAS */
void Cgemm(const char *transa, const char *transb, mpackint m, mpackint n,
	   mpackint k, mpcomplex alpha, mpcomplex * A, mpackint lda, mpcomplex * B, mpackint ldb, mpcomplex beta, mpcomplex * C, mpackint ldc);
void Rgemm(const char *transa, const char *transb, mpackint m, mpackint n,
	   mpackint k, mpreal alpha, mpreal * A, mpackint lda, mpreal * B, mpackint ldb, mpreal beta, mpreal * C, mpackint ldc);
void Csymm(const char *side, const char *uplo, mpackint m, mpackint n,
	   mpcomplex alpha, mpcomplex * A, mpackint lda, mpcomplex * B, mpackint ldb, mpcomplex beta, mpcomplex * C, mpackint ldc);
void Rsymm(const char *side, const char *uplo, mpackint m, mpackint n,
	   mpreal alpha, mpreal * A, mpackint lda, mpreal * B, mpackint ldb, mpreal beta, mpreal * C, mpackint ldc);
void Chemm(const char *side, const char *uplo, mpackint m, mpackint n,
	   mpcomplex alpha, mpcomplex * A, mpackint lda, mpcomplex * B, mpackint ldb, mpcomplex beta, mpcomplex * C, mpackint ldc);
void Csyrk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   mpcomplex alpha, mpcomplex * A, mpackint lda, mpcomplex beta, mpcomplex * C, mpackint ldc);
void Rsyrk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   mpreal alpha, mpreal * A, mpackint lda, mpreal beta, mpreal * C, mpackint ldc);
void Cherk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   mpreal alpha, mpcomplex * A, mpackint lda, mpreal beta, mpcomplex * C, mpackint ldc);
void Csyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    mpcomplex alpha, mpcomplex * A, mpackint lda, mpcomplex * B, mpackint ldb, mpcomplex beta, mpcomplex * C, mpackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    mpreal alpha, mpreal * A, mpackint lda, mpreal * B, mpackint ldb, mpreal beta, mpreal * C, mpackint ldc);
void Cher2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    mpcomplex alpha, mpcomplex * A, mpackint lda, mpcomplex * B, mpackint ldb, mpreal beta, mpcomplex * C, mpackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, mpcomplex alpha, mpcomplex * A, mpackint lda, mpcomplex * B, mpackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, mpreal alpha, mpreal * A, mpackint lda, mpreal * B, mpackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, mpcomplex alpha, mpcomplex * A, mpackint lda, mpcomplex * B, mpackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, mpreal alpha, mpreal * A, mpackint lda, mpreal * B, mpackint ldb);

#endif
