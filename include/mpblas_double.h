/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpblas_double.h,v 1.3 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPBLAS_DOUBLE_H_
#define _MPBLAS_DOUBLE_H_

#include "mplapack_config.h"
#include "mutils_double.h"
#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

/* LEVEL 1 MPBLAS */
void Rrot(mplapackint n, double * dx, mplapackint incx, double * dy, mplapackint incy, double c, double s);
void Rrotm(mplapackint n, double * x, mplapackint incx, double * y, mplapackint incy, double * dparam);
void CRrot(mplapackint n, std::complex<double> * cx, mplapackint incx, std::complex<double> * cy, mplapackint incy, double c, double s);
void Cswap(mplapackint n, std::complex<double> * cx, mplapackint incx, std::complex<double> * cy, mplapackint incy);
void Rswap(mplapackint n, double * dx, mplapackint incx, double * dy, mplapackint incy);
void Cscal(mplapackint n, std::complex<double> ca, std::complex<double> * cx, mplapackint incx);
void CRscal(mplapackint n, double sa, std::complex<double> * cx, mplapackint incx);
void Rscal(mplapackint n, double ca, double * cx, mplapackint incx);
void Ccopy(mplapackint n, std::complex<double> * cx, mplapackint incx, std::complex<double> * cy, mplapackint incy);
void Crotg(std::complex<double> * ca, std::complex<double> cb, double * c, std::complex<double> * s);
void Rrotg(double * da, double * db, double * c, double * s);
void Rcopy(mplapackint n, double * dx, mplapackint incx, double * dy, mplapackint incy);
void Raxpy(mplapackint n, double da, double * dx, mplapackint incx, double * dy, mplapackint incy);
void Caxpy(mplapackint n, std::complex<double> ca, std::complex<double> * cx, mplapackint incx, std::complex<double> * cy, mplapackint incy);
double Rdot(mplapackint n, double * dx, mplapackint incx, double * dy, mplapackint incy);
std::complex<double> Cdotc(mplapackint n, std::complex<double> * cx, mplapackint incx, std::complex<double> * cy, mplapackint incy);
std::complex<double> Cdotu(mplapackint n, std::complex<double> * cx, mplapackint incx, std::complex<double> * cy, mplapackint incy);
double RCnrm2(mplapackint n, std::complex<double> * x, mplapackint incx);
double Rnrm2(mplapackint n, double * x, mplapackint incx);
double RCasum(mplapackint n, std::complex<double> * zx, mplapackint incx);
double Rasum(mplapackint n, double * dx, mplapackint incx);
mplapackint iCamax(mplapackint n, std::complex<double> * dx, mplapackint incx);
mplapackint iRamax(mplapackint n, double * dx, mplapackint incx);
double RCabs1(std::complex<double> z);
mplapackint Mlsame_double(const char *a, const char *b);
void Mxerbla_double(const char *srname, int info);

/* LEVEL 2 MPBLAS */
void Cgemv(const char *trans, mplapackint m, mplapackint n, std::complex<double> alpha,
	   std::complex<double> * A, mplapackint lda, std::complex<double> * x, mplapackint incx, std::complex<double> beta, std::complex<double> * y, mplapackint incy);
void Rgemv(const char *trans, mplapackint m, mplapackint n, double alpha,
	   double * A, mplapackint lda, double * x, mplapackint incx, double beta, double * y, mplapackint incy);
void Cgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   std::complex<double> alpha, std::complex<double> * A, mplapackint lda, std::complex<double> * x, mplapackint incx, std::complex<double> beta, std::complex<double> * y, mplapackint incy);
void Rgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   double alpha, double * A, mplapackint lda, double * x, mplapackint incx, double beta, double * y, mplapackint incy);
void Chemv(const char *uplo, mplapackint n, std::complex<double> alpha, std::complex<double> * A,
	   mplapackint lda, std::complex<double> * x, mplapackint incx, std::complex<double> beta, std::complex<double> * y, mplapackint incy);
void Chbmv(const char *uplo, mplapackint n, mplapackint k, std::complex<double> alpha,
	   std::complex<double> * A, mplapackint lda, std::complex<double> * x, mplapackint incx, std::complex<double> beta, std::complex<double> * y, mplapackint incy);
void Chpmv(const char *uplo, mplapackint n, std::complex<double> alpha, std::complex<double> * AP,
	   std::complex<double> * x, mplapackint incx, std::complex<double> beta, std::complex<double> * y, mplapackint incy);
void Rsymv(const char *uplo, mplapackint n, double alpha, double * A,
	   mplapackint lda, double * x, mplapackint incx, double beta, double * y, mplapackint incy);
void Rsbmv(const char *uplo, mplapackint n, mplapackint k, double alpha,
	   double * A, mplapackint lda, double * x, mplapackint incx, double beta, double * y, mplapackint incy);
void Rspmv(const char *uplo, mplapackint n, double alpha, double * AP, double * x, mplapackint incx, double beta, double * y, mplapackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<double> * A, mplapackint lda, std::complex<double> * x, mplapackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, double * A, mplapackint lda, double * x, mplapackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, std::complex<double> * A, mplapackint lda, std::complex<double> * x, mplapackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, double * A, mplapackint lda, double * x, mplapackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<double> * AP, std::complex<double> * x, mplapackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, double * AP, double * x, mplapackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<double> * A, mplapackint lda, std::complex<double> * x, mplapackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, double * A, mplapackint lda, double * x, mplapackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, std::complex<double> * A, mplapackint lda, std::complex<double> * x, mplapackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, double * A, mplapackint lda, double * x, mplapackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<double> * AP, std::complex<double> * x, mplapackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, double * AP, double * x, mplapackint incx);
void Rger(mplapackint m, mplapackint n, double alpha, double * x, mplapackint incx, double * y, mplapackint incy, double * A, mplapackint lda);
void Cgeru(mplapackint m, mplapackint n, std::complex<double> alpha, std::complex<double> * x, mplapackint incx, std::complex<double> * y, mplapackint incy, std::complex<double> * A, mplapackint lda);
void Cgerc(mplapackint m, mplapackint n, std::complex<double> alpha, std::complex<double> * x, mplapackint incx, std::complex<double> * y, mplapackint incy, std::complex<double> * A, mplapackint lda);
void Cher(const char *uplo, mplapackint n, double alpha, std::complex<double> * x, mplapackint incx, std::complex<double> * A, mplapackint lda);
void Chpr(const char *uplo, mplapackint n, double alpha, std::complex<double> * x, mplapackint incx, std::complex<double> * AP);
void Cher2(const char *uplo, mplapackint n, std::complex<double> alpha, std::complex<double> * x,
	   mplapackint incx, std::complex<double> * y, mplapackint incy, std::complex<double> * A, mplapackint lda);
void Chpr2(const char *uplo, mplapackint n, std::complex<double> alpha, std::complex<double> * x, mplapackint incx, std::complex<double> * y, mplapackint incy, std::complex<double> * AP);
void Rsyr(const char *uplo, mplapackint n, double alpha, double * x, mplapackint incx, double * A, mplapackint lda);
void Rspr(const char *uplo, mplapackint n, double alpha, double * x, mplapackint incx, double * AP);
void Rsyr2(const char *uplo, mplapackint n, double alpha, double * x, mplapackint incx, double * y, mplapackint incy, double * A, mplapackint lda);
void Rspr2(const char *uplo, mplapackint n, double alpha, double * x, mplapackint incx, double * y, mplapackint incy, double * AP);

/* LEVEL 3 MPBLAS */
void Cgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, std::complex<double> alpha, std::complex<double> * A, mplapackint lda, std::complex<double> * B, mplapackint ldb, std::complex<double> beta, std::complex<double> * C, mplapackint ldc);
void Rgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, double alpha, double * A, mplapackint lda, double * B, mplapackint ldb, double beta, double * C, mplapackint ldc);
void Csymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   std::complex<double> alpha, std::complex<double> * A, mplapackint lda, std::complex<double> * B, mplapackint ldb, std::complex<double> beta, std::complex<double> * C, mplapackint ldc);
void Rsymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   double alpha, double * A, mplapackint lda, double * B, mplapackint ldb, double beta, double * C, mplapackint ldc);
void Chemm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   std::complex<double> alpha, std::complex<double> * A, mplapackint lda, std::complex<double> * B, mplapackint ldb, std::complex<double> beta, std::complex<double> * C, mplapackint ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   std::complex<double> alpha, std::complex<double> * A, mplapackint lda, std::complex<double> beta, std::complex<double> * C, mplapackint ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   double alpha, double * A, mplapackint lda, double beta, double * C, mplapackint ldc);
void Cherk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   double alpha, std::complex<double> * A, mplapackint lda, double beta, std::complex<double> * C, mplapackint ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    std::complex<double> alpha, std::complex<double> * A, mplapackint lda, std::complex<double> * B, mplapackint ldb, std::complex<double> beta, std::complex<double> * C, mplapackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    double alpha, double * A, mplapackint lda, double * B, mplapackint ldb, double beta, double * C, mplapackint ldc);
void Cher2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    std::complex<double> alpha, std::complex<double> * A, mplapackint lda, std::complex<double> * B, mplapackint ldb, double beta, std::complex<double> * C, mplapackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, std::complex<double> alpha, std::complex<double> * A, mplapackint lda, std::complex<double> * B, mplapackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, double alpha, double * A, mplapackint lda, double * B, mplapackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, std::complex<double> alpha, std::complex<double> * A, mplapackint lda, std::complex<double> * B, mplapackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, double alpha, double * A, mplapackint lda, double * B, mplapackint ldb);

#endif
