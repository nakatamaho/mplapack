/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mblas_double.h,v 1.3 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MBLAS_DOUBLE_H_
#define _MBLAS_DOUBLE_H_

#include <mpack_config.h>
#include <mutils_double.h>
#if !defined __MPACK_ERRNO__
#define _MPACK_EXTERN_ extern
#else
#define _MPACK_EXTERN_
#endif
_MPACK_EXTERN_ int mpack_errno;

/* LEVEL 1 MBLAS */
void Rrot(mpackint n, double * dx, mpackint incx, double * dy, mpackint incy, double c, double s);
void Rrotm(mpackint n, double * x, mpackint incx, double * y, mpackint incy, double * dparam);
void CRrot(mpackint n, std::complex<double> * cx, mpackint incx, std::complex<double> * cy, mpackint incy, double c, double s);
void Cswap(mpackint n, std::complex<double> * cx, mpackint incx, std::complex<double> * cy, mpackint incy);
void Rswap(mpackint n, double * dx, mpackint incx, double * dy, mpackint incy);
void Cscal(mpackint n, std::complex<double> ca, std::complex<double> * cx, mpackint incx);
void CRscal(mpackint n, double sa, std::complex<double> * cx, mpackint incx);
void Rscal(mpackint n, double ca, double * cx, mpackint incx);
void Ccopy(mpackint n, std::complex<double> * cx, mpackint incx, std::complex<double> * cy, mpackint incy);
void Crotg(std::complex<double> * ca, std::complex<double> cb, double * c, std::complex<double> * s);
void Rrotg(double * da, double * db, double * c, double * s);
void Rcopy(mpackint n, double * dx, mpackint incx, double * dy, mpackint incy);
void Raxpy(mpackint n, double da, double * dx, mpackint incx, double * dy, mpackint incy);
void Caxpy(mpackint n, std::complex<double> ca, std::complex<double> * cx, mpackint incx, std::complex<double> * cy, mpackint incy);
double Rdot(mpackint n, double * dx, mpackint incx, double * dy, mpackint incy);
std::complex<double> Cdotc(mpackint n, std::complex<double> * cx, mpackint incx, std::complex<double> * cy, mpackint incy);
std::complex<double> Cdotu(mpackint n, std::complex<double> * cx, mpackint incx, std::complex<double> * cy, mpackint incy);
double RCnrm2(mpackint n, std::complex<double> * x, mpackint incx);
double Rnrm2(mpackint n, double * x, mpackint incx);
double RCasum(mpackint n, std::complex<double> * zx, mpackint incx);
double Rasum(mpackint n, double * dx, mpackint incx);
mpackint iCamax(mpackint n, std::complex<double> * dx, mpackint incx);
mpackint iRamax(mpackint n, double * dx, mpackint incx);
double RCabs1(std::complex<double> z);
mpackint Mlsame_double(const char *a, const char *b);
void Mxerbla_double(const char *srname, int info);

/* LEVEL 2 MBLAS */
void Cgemv(const char *trans, mpackint m, mpackint n, std::complex<double> alpha,
	   std::complex<double> * A, mpackint lda, std::complex<double> * x, mpackint incx, std::complex<double> beta, std::complex<double> * y, mpackint incy);
void Rgemv(const char *trans, mpackint m, mpackint n, double alpha,
	   double * A, mpackint lda, double * x, mpackint incx, double beta, double * y, mpackint incy);
void Cgbmv(const char *trans, mpackint m, mpackint n, mpackint kl, mpackint ku,
	   std::complex<double> alpha, std::complex<double> * A, mpackint lda, std::complex<double> * x, mpackint incx, std::complex<double> beta, std::complex<double> * y, mpackint incy);
void Rgbmv(const char *trans, mpackint m, mpackint n, mpackint kl, mpackint ku,
	   double alpha, double * A, mpackint lda, double * x, mpackint incx, double beta, double * y, mpackint incy);
void Chemv(const char *uplo, mpackint n, std::complex<double> alpha, std::complex<double> * A,
	   mpackint lda, std::complex<double> * x, mpackint incx, std::complex<double> beta, std::complex<double> * y, mpackint incy);
void Chbmv(const char *uplo, mpackint n, mpackint k, std::complex<double> alpha,
	   std::complex<double> * A, mpackint lda, std::complex<double> * x, mpackint incx, std::complex<double> beta, std::complex<double> * y, mpackint incy);
void Chpmv(const char *uplo, mpackint n, std::complex<double> alpha, std::complex<double> * AP,
	   std::complex<double> * x, mpackint incx, std::complex<double> beta, std::complex<double> * y, mpackint incy);
void Rsymv(const char *uplo, mpackint n, double alpha, double * A,
	   mpackint lda, double * x, mpackint incx, double beta, double * y, mpackint incy);
void Rsbmv(const char *uplo, mpackint n, mpackint k, double alpha,
	   double * A, mpackint lda, double * x, mpackint incx, double beta, double * y, mpackint incy);
void Rspmv(const char *uplo, mpackint n, double alpha, double * AP, double * x, mpackint incx, double beta, double * y, mpackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mpackint n, std::complex<double> * A, mpackint lda, std::complex<double> * x, mpackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mpackint n, double * A, mpackint lda, double * x, mpackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, std::complex<double> * A, mpackint lda, std::complex<double> * x, mpackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, double * A, mpackint lda, double * x, mpackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mpackint n, std::complex<double> * AP, std::complex<double> * x, mpackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mpackint n, double * AP, double * x, mpackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mpackint n, std::complex<double> * A, mpackint lda, std::complex<double> * x, mpackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mpackint n, double * A, mpackint lda, double * x, mpackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, std::complex<double> * A, mpackint lda, std::complex<double> * x, mpackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, double * A, mpackint lda, double * x, mpackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mpackint n, std::complex<double> * AP, std::complex<double> * x, mpackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mpackint n, double * AP, double * x, mpackint incx);
void Rger(mpackint m, mpackint n, double alpha, double * x, mpackint incx, double * y, mpackint incy, double * A, mpackint lda);
void Cgeru(mpackint m, mpackint n, std::complex<double> alpha, std::complex<double> * x, mpackint incx, std::complex<double> * y, mpackint incy, std::complex<double> * A, mpackint lda);
void Cgerc(mpackint m, mpackint n, std::complex<double> alpha, std::complex<double> * x, mpackint incx, std::complex<double> * y, mpackint incy, std::complex<double> * A, mpackint lda);
void Cher(const char *uplo, mpackint n, double alpha, std::complex<double> * x, mpackint incx, std::complex<double> * A, mpackint lda);
void Chpr(const char *uplo, mpackint n, double alpha, std::complex<double> * x, mpackint incx, std::complex<double> * AP);
void Cher2(const char *uplo, mpackint n, std::complex<double> alpha, std::complex<double> * x,
	   mpackint incx, std::complex<double> * y, mpackint incy, std::complex<double> * A, mpackint lda);
void Chpr2(const char *uplo, mpackint n, std::complex<double> alpha, std::complex<double> * x, mpackint incx, std::complex<double> * y, mpackint incy, std::complex<double> * AP);
void Rsyr(const char *uplo, mpackint n, double alpha, double * x, mpackint incx, double * A, mpackint lda);
void Rspr(const char *uplo, mpackint n, double alpha, double * x, mpackint incx, double * AP);
void Rsyr2(const char *uplo, mpackint n, double alpha, double * x, mpackint incx, double * y, mpackint incy, double * A, mpackint lda);
void Rspr2(const char *uplo, mpackint n, double alpha, double * x, mpackint incx, double * y, mpackint incy, double * AP);

/* LEVEL 3 MBLAS */
void Cgemm(const char *transa, const char *transb, mpackint m, mpackint n,
	   mpackint k, std::complex<double> alpha, std::complex<double> * A, mpackint lda, std::complex<double> * B, mpackint ldb, std::complex<double> beta, std::complex<double> * C, mpackint ldc);
void Rgemm(const char *transa, const char *transb, mpackint m, mpackint n,
	   mpackint k, double alpha, double * A, mpackint lda, double * B, mpackint ldb, double beta, double * C, mpackint ldc);
void Csymm(const char *side, const char *uplo, mpackint m, mpackint n,
	   std::complex<double> alpha, std::complex<double> * A, mpackint lda, std::complex<double> * B, mpackint ldb, std::complex<double> beta, std::complex<double> * C, mpackint ldc);
void Rsymm(const char *side, const char *uplo, mpackint m, mpackint n,
	   double alpha, double * A, mpackint lda, double * B, mpackint ldb, double beta, double * C, mpackint ldc);
void Chemm(const char *side, const char *uplo, mpackint m, mpackint n,
	   std::complex<double> alpha, std::complex<double> * A, mpackint lda, std::complex<double> * B, mpackint ldb, std::complex<double> beta, std::complex<double> * C, mpackint ldc);
void Csyrk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   std::complex<double> alpha, std::complex<double> * A, mpackint lda, std::complex<double> beta, std::complex<double> * C, mpackint ldc);
void Rsyrk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   double alpha, double * A, mpackint lda, double beta, double * C, mpackint ldc);
void Cherk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   double alpha, std::complex<double> * A, mpackint lda, double beta, std::complex<double> * C, mpackint ldc);
void Csyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    std::complex<double> alpha, std::complex<double> * A, mpackint lda, std::complex<double> * B, mpackint ldb, std::complex<double> beta, std::complex<double> * C, mpackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    double alpha, double * A, mpackint lda, double * B, mpackint ldb, double beta, double * C, mpackint ldc);
void Cher2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    std::complex<double> alpha, std::complex<double> * A, mpackint lda, std::complex<double> * B, mpackint ldb, double beta, std::complex<double> * C, mpackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, std::complex<double> alpha, std::complex<double> * A, mpackint lda, std::complex<double> * B, mpackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, double alpha, double * A, mpackint lda, double * B, mpackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, std::complex<double> alpha, std::complex<double> * A, mpackint lda, std::complex<double> * B, mpackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, double alpha, double * A, mpackint lda, double * B, mpackint ldb);

#endif
