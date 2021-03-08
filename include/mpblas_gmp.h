/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpblas_gmp.h,v 1.12 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPBLAS_GMP_H_
#define _MPBLAS_GMP_H_

#include "mplapack_config.h"
#include "gmpxx.h"
#include "mpc_class.h"
#include "mutils_gmp.h"
#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

/* LEVEL 1 MPBLAS */
void Rrot(mplapackint n, mpf_class * dx, mplapackint incx, mpf_class * dy, mplapackint incy, mpf_class c, mpf_class s);
void Rrotm(mplapackint n, mpf_class * x, mplapackint incx, mpf_class * y, mplapackint incy, mpf_class * dparam);
void CRrot(mplapackint n, mpc_class * cx, mplapackint incx, mpc_class * cy, mplapackint incy, mpf_class c, mpf_class s);
void Cswap(mplapackint n, mpc_class * cx, mplapackint incx, mpc_class * cy, mplapackint incy);
void Rswap(mplapackint n, mpf_class * dx, mplapackint incx, mpf_class * dy, mplapackint incy);
void Cscal(mplapackint n, mpc_class ca, mpc_class * cx, mplapackint incx);
void CRscal(mplapackint n, mpf_class sa, mpc_class * cx, mplapackint incx);
void Rscal(mplapackint n, mpf_class ca, mpf_class * cx, mplapackint incx);
void Ccopy(mplapackint n, mpc_class * cx, mplapackint incx, mpc_class * cy, mplapackint incy);
void Crotg(mpc_class * ca, mpc_class cb, mpf_class * c, mpc_class * s);
void Rrotg(mpf_class * da, mpf_class * db, mpf_class * c, mpf_class * s);
void Rcopy(mplapackint n, mpf_class * dx, mplapackint incx, mpf_class * dy, mplapackint incy);
void Raxpy(mplapackint n, mpf_class da, mpf_class * dx, mplapackint incx, mpf_class * dy, mplapackint incy);
void Caxpy(mplapackint n, mpc_class ca, mpc_class * cx, mplapackint incx, mpc_class * cy, mplapackint incy);
mpf_class Rdot(mplapackint n, mpf_class * dx, mplapackint incx, mpf_class * dy, mplapackint incy);
mpc_class Cdotc(mplapackint n, mpc_class * cx, mplapackint incx, mpc_class * cy, mplapackint incy);
mpc_class Cdotu(mplapackint n, mpc_class * cx, mplapackint incx, mpc_class * cy, mplapackint incy);
mpf_class RCnrm2(mplapackint n, mpc_class * x, mplapackint incx);
mpf_class Rnrm2(mplapackint n, mpf_class * x, mplapackint incx);
mpf_class RCasum(mplapackint n, mpc_class * zx, mplapackint incx);
mpf_class Rasum(mplapackint n, mpf_class * dx, mplapackint incx);
mplapackint iCamax(mplapackint n, mpc_class * dx, mplapackint incx);
mplapackint iRamax(mplapackint n, mpf_class * dx, mplapackint incx);
mpf_class RCabs1(mpc_class z);
mplapackint Mlsame_gmp(const char *a, const char *b);
void Mxerbla_gmp(const char *srname, int info);

/* LEVEL 2 MPBLAS */
void Cgemv(const char *trans, mplapackint m, mplapackint n, mpc_class alpha,
	   mpc_class * A, mplapackint lda, mpc_class * x, mplapackint incx, mpc_class beta, mpc_class * y, mplapackint incy);
void Rgemv(const char *trans, mplapackint m, mplapackint n, mpf_class alpha,
	   mpf_class * A, mplapackint lda, mpf_class * x, mplapackint incx, mpf_class beta, mpf_class * y, mplapackint incy);
void Cgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   mpc_class alpha, mpc_class * A, mplapackint lda, mpc_class * x, mplapackint incx, mpc_class beta, mpc_class * y, mplapackint incy);
void Rgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   mpf_class alpha, mpf_class * A, mplapackint lda, mpf_class * x, mplapackint incx, mpf_class beta, mpf_class * y, mplapackint incy);
void Chemv(const char *uplo, mplapackint n, mpc_class alpha, mpc_class * A,
	   mplapackint lda, mpc_class * x, mplapackint incx, mpc_class beta, mpc_class * y, mplapackint incy);
void Chbmv(const char *uplo, mplapackint n, mplapackint k, mpc_class alpha,
	   mpc_class * A, mplapackint lda, mpc_class * x, mplapackint incx, mpc_class beta, mpc_class * y, mplapackint incy);
void Chpmv(const char *uplo, mplapackint n, mpc_class alpha, mpc_class * AP,
	   mpc_class * x, mplapackint incx, mpc_class beta, mpc_class * y, mplapackint incy);
void Rsymv(const char *uplo, mplapackint n, mpf_class alpha, mpf_class * A,
	   mplapackint lda, mpf_class * x, mplapackint incx, mpf_class beta, mpf_class * y, mplapackint incy);
void Rsbmv(const char *uplo, mplapackint n, mplapackint k, mpf_class alpha,
	   mpf_class * A, mplapackint lda, mpf_class * x, mplapackint incx, mpf_class beta, mpf_class * y, mplapackint incy);
void Rspmv(const char *uplo, mplapackint n, mpf_class alpha, mpf_class * AP, mpf_class * x, mplapackint incx, mpf_class beta, mpf_class * y, mplapackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpc_class * A, mplapackint lda, mpc_class * x, mplapackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpf_class * A, mplapackint lda, mpf_class * x, mplapackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, mpc_class * A, mplapackint lda, mpc_class * x, mplapackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, mpf_class * A, mplapackint lda, mpf_class * x, mplapackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpc_class * AP, mpc_class * x, mplapackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpf_class * AP, mpf_class * x, mplapackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpc_class * A, mplapackint lda, mpc_class * x, mplapackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpf_class * A, mplapackint lda, mpf_class * x, mplapackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, mpc_class * A, mplapackint lda, mpc_class * x, mplapackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, mpf_class * A, mplapackint lda, mpf_class * x, mplapackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpc_class * AP, mpc_class * x, mplapackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, mpf_class * AP, mpf_class * x, mplapackint incx);
void Rger(mplapackint m, mplapackint n, mpf_class alpha, mpf_class * x, mplapackint incx, mpf_class * y, mplapackint incy, mpf_class * A, mplapackint lda);
void Cgeru(mplapackint m, mplapackint n, mpc_class alpha, mpc_class * x, mplapackint incx, mpc_class * y, mplapackint incy, mpc_class * A, mplapackint lda);
void Cgerc(mplapackint m, mplapackint n, mpc_class alpha, mpc_class * x, mplapackint incx, mpc_class * y, mplapackint incy, mpc_class * A, mplapackint lda);
void Cher(const char *uplo, mplapackint n, mpf_class alpha, mpc_class * x, mplapackint incx, mpc_class * A, mplapackint lda);
void Chpr(const char *uplo, mplapackint n, mpf_class alpha, mpc_class * x, mplapackint incx, mpc_class * AP);
void Cher2(const char *uplo, mplapackint n, mpc_class alpha, mpc_class * x,
	   mplapackint incx, mpc_class * y, mplapackint incy, mpc_class * A, mplapackint lda);
void Chpr2(const char *uplo, mplapackint n, mpc_class alpha, mpc_class * x, mplapackint incx, mpc_class * y, mplapackint incy, mpc_class * AP);
void Rsyr(const char *uplo, mplapackint n, mpf_class alpha, mpf_class * x, mplapackint incx, mpf_class * A, mplapackint lda);
void Rspr(const char *uplo, mplapackint n, mpf_class alpha, mpf_class * x, mplapackint incx, mpf_class * AP);
void Rsyr2(const char *uplo, mplapackint n, mpf_class alpha, mpf_class * x, mplapackint incx, mpf_class * y, mplapackint incy, mpf_class * A, mplapackint lda);
void Rspr2(const char *uplo, mplapackint n, mpf_class alpha, mpf_class * x, mplapackint incx, mpf_class * y, mplapackint incy, mpf_class * AP);

/* LEVEL 3 MPBLAS */
void Cgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, mpc_class alpha, mpc_class * A, mplapackint lda, mpc_class * B, mplapackint ldb, mpc_class beta, mpc_class * C, mplapackint ldc);
void Rgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, mpf_class alpha, mpf_class * A, mplapackint lda, mpf_class * B, mplapackint ldb, mpf_class beta, mpf_class * C, mplapackint ldc);
void Csymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   mpc_class alpha, mpc_class * A, mplapackint lda, mpc_class * B, mplapackint ldb, mpc_class beta, mpc_class * C, mplapackint ldc);
void Rsymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   mpf_class alpha, mpf_class * A, mplapackint lda, mpf_class * B, mplapackint ldb, mpf_class beta, mpf_class * C, mplapackint ldc);
void Chemm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   mpc_class alpha, mpc_class * A, mplapackint lda, mpc_class * B, mplapackint ldb, mpc_class beta, mpc_class * C, mplapackint ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   mpc_class alpha, mpc_class * A, mplapackint lda, mpc_class beta, mpc_class * C, mplapackint ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   mpf_class alpha, mpf_class * A, mplapackint lda, mpf_class beta, mpf_class * C, mplapackint ldc);
void Cherk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   mpf_class alpha, mpc_class * A, mplapackint lda, mpf_class beta, mpc_class * C, mplapackint ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    mpc_class alpha, mpc_class * A, mplapackint lda, mpc_class * B, mplapackint ldb, mpc_class beta, mpc_class * C, mplapackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    mpf_class alpha, mpf_class * A, mplapackint lda, mpf_class * B, mplapackint ldb, mpf_class beta, mpf_class * C, mplapackint ldc);
void Cher2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    mpc_class alpha, mpc_class * A, mplapackint lda, mpc_class * B, mplapackint ldb, mpf_class beta, mpc_class * C, mplapackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, mpc_class alpha, mpc_class * A, mplapackint lda, mpc_class * B, mplapackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, mpf_class alpha, mpf_class * A, mplapackint lda, mpf_class * B, mplapackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, mpc_class alpha, mpc_class * A, mplapackint lda, mpc_class * B, mplapackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, mpf_class alpha, mpf_class * A, mplapackint lda, mpf_class * B, mplapackint ldb);

#endif
