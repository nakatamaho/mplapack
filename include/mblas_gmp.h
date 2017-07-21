/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mblas_gmp.h,v 1.12 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MBLAS_GMP_H_
#define _MBLAS_GMP_H_

#include "mpack_config.h"
#include "gmpxx.h"
#include "mpc_class.h"
#include "mutils_gmp.h"
#if !defined __MPACK_ERRNO__
#define _MPACK_EXTERN_ extern
#else
#define _MPACK_EXTERN_
#endif
_MPACK_EXTERN_ int mpack_errno;

/* LEVEL 1 MBLAS */
void Rrot(mpackint n, mpf_class * dx, mpackint incx, mpf_class * dy, mpackint incy, mpf_class c, mpf_class s);
void Rrotm(mpackint n, mpf_class * x, mpackint incx, mpf_class * y, mpackint incy, mpf_class * dparam);
void CRrot(mpackint n, mpc_class * cx, mpackint incx, mpc_class * cy, mpackint incy, mpf_class c, mpf_class s);
void Cswap(mpackint n, mpc_class * cx, mpackint incx, mpc_class * cy, mpackint incy);
void Rswap(mpackint n, mpf_class * dx, mpackint incx, mpf_class * dy, mpackint incy);
void Cscal(mpackint n, mpc_class ca, mpc_class * cx, mpackint incx);
void CRscal(mpackint n, mpf_class sa, mpc_class * cx, mpackint incx);
void Rscal(mpackint n, mpf_class ca, mpf_class * cx, mpackint incx);
void Ccopy(mpackint n, mpc_class * cx, mpackint incx, mpc_class * cy, mpackint incy);
void Crotg(mpc_class * ca, mpc_class cb, mpf_class * c, mpc_class * s);
void Rrotg(mpf_class * da, mpf_class * db, mpf_class * c, mpf_class * s);
void Rcopy(mpackint n, mpf_class * dx, mpackint incx, mpf_class * dy, mpackint incy);
void Raxpy(mpackint n, mpf_class da, mpf_class * dx, mpackint incx, mpf_class * dy, mpackint incy);
void Caxpy(mpackint n, mpc_class ca, mpc_class * cx, mpackint incx, mpc_class * cy, mpackint incy);
mpf_class Rdot(mpackint n, mpf_class * dx, mpackint incx, mpf_class * dy, mpackint incy);
mpc_class Cdotc(mpackint n, mpc_class * cx, mpackint incx, mpc_class * cy, mpackint incy);
mpc_class Cdotu(mpackint n, mpc_class * cx, mpackint incx, mpc_class * cy, mpackint incy);
mpf_class RCnrm2(mpackint n, mpc_class * x, mpackint incx);
mpf_class Rnrm2(mpackint n, mpf_class * x, mpackint incx);
mpf_class RCasum(mpackint n, mpc_class * zx, mpackint incx);
mpf_class Rasum(mpackint n, mpf_class * dx, mpackint incx);
mpackint iCamax(mpackint n, mpc_class * dx, mpackint incx);
mpackint iRamax(mpackint n, mpf_class * dx, mpackint incx);
mpf_class RCabs1(mpc_class z);
mpackint Mlsame_gmp(const char *a, const char *b);
void Mxerbla_gmp(const char *srname, int info);

/* LEVEL 2 MBLAS */
void Cgemv(const char *trans, mpackint m, mpackint n, mpc_class alpha,
	   mpc_class * A, mpackint lda, mpc_class * x, mpackint incx, mpc_class beta, mpc_class * y, mpackint incy);
void Rgemv(const char *trans, mpackint m, mpackint n, mpf_class alpha,
	   mpf_class * A, mpackint lda, mpf_class * x, mpackint incx, mpf_class beta, mpf_class * y, mpackint incy);
void Cgbmv(const char *trans, mpackint m, mpackint n, mpackint kl, mpackint ku,
	   mpc_class alpha, mpc_class * A, mpackint lda, mpc_class * x, mpackint incx, mpc_class beta, mpc_class * y, mpackint incy);
void Rgbmv(const char *trans, mpackint m, mpackint n, mpackint kl, mpackint ku,
	   mpf_class alpha, mpf_class * A, mpackint lda, mpf_class * x, mpackint incx, mpf_class beta, mpf_class * y, mpackint incy);
void Chemv(const char *uplo, mpackint n, mpc_class alpha, mpc_class * A,
	   mpackint lda, mpc_class * x, mpackint incx, mpc_class beta, mpc_class * y, mpackint incy);
void Chbmv(const char *uplo, mpackint n, mpackint k, mpc_class alpha,
	   mpc_class * A, mpackint lda, mpc_class * x, mpackint incx, mpc_class beta, mpc_class * y, mpackint incy);
void Chpmv(const char *uplo, mpackint n, mpc_class alpha, mpc_class * AP,
	   mpc_class * x, mpackint incx, mpc_class beta, mpc_class * y, mpackint incy);
void Rsymv(const char *uplo, mpackint n, mpf_class alpha, mpf_class * A,
	   mpackint lda, mpf_class * x, mpackint incx, mpf_class beta, mpf_class * y, mpackint incy);
void Rsbmv(const char *uplo, mpackint n, mpackint k, mpf_class alpha,
	   mpf_class * A, mpackint lda, mpf_class * x, mpackint incx, mpf_class beta, mpf_class * y, mpackint incy);
void Rspmv(const char *uplo, mpackint n, mpf_class alpha, mpf_class * AP, mpf_class * x, mpackint incx, mpf_class beta, mpf_class * y, mpackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mpackint n, mpc_class * A, mpackint lda, mpc_class * x, mpackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mpackint n, mpf_class * A, mpackint lda, mpf_class * x, mpackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, mpc_class * A, mpackint lda, mpc_class * x, mpackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, mpf_class * A, mpackint lda, mpf_class * x, mpackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mpackint n, mpc_class * AP, mpc_class * x, mpackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mpackint n, mpf_class * AP, mpf_class * x, mpackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mpackint n, mpc_class * A, mpackint lda, mpc_class * x, mpackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mpackint n, mpf_class * A, mpackint lda, mpf_class * x, mpackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, mpc_class * A, mpackint lda, mpc_class * x, mpackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, mpf_class * A, mpackint lda, mpf_class * x, mpackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mpackint n, mpc_class * AP, mpc_class * x, mpackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mpackint n, mpf_class * AP, mpf_class * x, mpackint incx);
void Rger(mpackint m, mpackint n, mpf_class alpha, mpf_class * x, mpackint incx, mpf_class * y, mpackint incy, mpf_class * A, mpackint lda);
void Cgeru(mpackint m, mpackint n, mpc_class alpha, mpc_class * x, mpackint incx, mpc_class * y, mpackint incy, mpc_class * A, mpackint lda);
void Cgerc(mpackint m, mpackint n, mpc_class alpha, mpc_class * x, mpackint incx, mpc_class * y, mpackint incy, mpc_class * A, mpackint lda);
void Cher(const char *uplo, mpackint n, mpf_class alpha, mpc_class * x, mpackint incx, mpc_class * A, mpackint lda);
void Chpr(const char *uplo, mpackint n, mpf_class alpha, mpc_class * x, mpackint incx, mpc_class * AP);
void Cher2(const char *uplo, mpackint n, mpc_class alpha, mpc_class * x,
	   mpackint incx, mpc_class * y, mpackint incy, mpc_class * A, mpackint lda);
void Chpr2(const char *uplo, mpackint n, mpc_class alpha, mpc_class * x, mpackint incx, mpc_class * y, mpackint incy, mpc_class * AP);
void Rsyr(const char *uplo, mpackint n, mpf_class alpha, mpf_class * x, mpackint incx, mpf_class * A, mpackint lda);
void Rspr(const char *uplo, mpackint n, mpf_class alpha, mpf_class * x, mpackint incx, mpf_class * AP);
void Rsyr2(const char *uplo, mpackint n, mpf_class alpha, mpf_class * x, mpackint incx, mpf_class * y, mpackint incy, mpf_class * A, mpackint lda);
void Rspr2(const char *uplo, mpackint n, mpf_class alpha, mpf_class * x, mpackint incx, mpf_class * y, mpackint incy, mpf_class * AP);

/* LEVEL 3 MBLAS */
void Cgemm(const char *transa, const char *transb, mpackint m, mpackint n,
	   mpackint k, mpc_class alpha, mpc_class * A, mpackint lda, mpc_class * B, mpackint ldb, mpc_class beta, mpc_class * C, mpackint ldc);
void Rgemm(const char *transa, const char *transb, mpackint m, mpackint n,
	   mpackint k, mpf_class alpha, mpf_class * A, mpackint lda, mpf_class * B, mpackint ldb, mpf_class beta, mpf_class * C, mpackint ldc);
void Csymm(const char *side, const char *uplo, mpackint m, mpackint n,
	   mpc_class alpha, mpc_class * A, mpackint lda, mpc_class * B, mpackint ldb, mpc_class beta, mpc_class * C, mpackint ldc);
void Rsymm(const char *side, const char *uplo, mpackint m, mpackint n,
	   mpf_class alpha, mpf_class * A, mpackint lda, mpf_class * B, mpackint ldb, mpf_class beta, mpf_class * C, mpackint ldc);
void Chemm(const char *side, const char *uplo, mpackint m, mpackint n,
	   mpc_class alpha, mpc_class * A, mpackint lda, mpc_class * B, mpackint ldb, mpc_class beta, mpc_class * C, mpackint ldc);
void Csyrk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   mpc_class alpha, mpc_class * A, mpackint lda, mpc_class beta, mpc_class * C, mpackint ldc);
void Rsyrk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   mpf_class alpha, mpf_class * A, mpackint lda, mpf_class beta, mpf_class * C, mpackint ldc);
void Cherk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   mpf_class alpha, mpc_class * A, mpackint lda, mpf_class beta, mpc_class * C, mpackint ldc);
void Csyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    mpc_class alpha, mpc_class * A, mpackint lda, mpc_class * B, mpackint ldb, mpc_class beta, mpc_class * C, mpackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    mpf_class alpha, mpf_class * A, mpackint lda, mpf_class * B, mpackint ldb, mpf_class beta, mpf_class * C, mpackint ldc);
void Cher2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    mpc_class alpha, mpc_class * A, mpackint lda, mpc_class * B, mpackint ldb, mpf_class beta, mpc_class * C, mpackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, mpc_class alpha, mpc_class * A, mpackint lda, mpc_class * B, mpackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, mpf_class alpha, mpf_class * A, mpackint lda, mpf_class * B, mpackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, mpc_class alpha, mpc_class * A, mpackint lda, mpc_class * B, mpackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, mpf_class alpha, mpf_class * A, mpackint lda, mpf_class * B, mpackint ldb);

#endif
