/*
 * Copyright (c) 2012
 *	Nakata, Maho
 * 	All rights reserved.
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

#ifndef _MPBLAS_LONGDOUBLE_H_
#define _MPBLAS_LONGDOUBLE_H_

#include "mutils_longdouble.h"
#include "mplapack_config.h"

#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

/* LEVEL 1 MPBLAS */
void Rrot(mplapackint n, long double * dx, mplapackint incx, long double * dy, mplapackint incy, long double c, long double s);
void Rrotm(mplapackint n, long double * x, mplapackint incx, long double * y, mplapackint incy, long double * dparam);
void CRrot(mplapackint n, std::complex<long double> * cx, mplapackint incx, std::complex<long double> * cy, mplapackint incy, long double c, long double s);
void Cswap(mplapackint n, std::complex<long double> * cx, mplapackint incx, std::complex<long double> * cy, mplapackint incy);
void Rswap(mplapackint n, long double * dx, mplapackint incx, long double * dy, mplapackint incy);
void Cscal(mplapackint n, std::complex<long double> ca, std::complex<long double> * cx, mplapackint incx);
void CRscal(mplapackint n, long double sa, std::complex<long double> * cx, mplapackint incx);
void Rscal(mplapackint n, long double ca, long double * cx, mplapackint incx);
void Ccopy(mplapackint n, std::complex<long double> * cx, mplapackint incx, std::complex<long double> * cy, mplapackint incy);
void Crotg(std::complex<long double> * ca, std::complex<long double> cb, long double * c, std::complex<long double> * s);
void Rrotg(long double * da, long double * db, long double * c, long double * s);
void Rcopy(mplapackint n, long double * dx, mplapackint incx, long double * dy, mplapackint incy);
void Raxpy(mplapackint n, long double da, long double * dx, mplapackint incx, long double * dy, mplapackint incy);
void Caxpy(mplapackint n, std::complex<long double> ca, std::complex<long double> * cx, mplapackint incx, std::complex<long double> * cy, mplapackint incy);
long double Rdot(mplapackint n, long double * dx, mplapackint incx, long double * dy, mplapackint incy);
std::complex<long double> Cdotc(mplapackint n, std::complex<long double> * cx, mplapackint incx, std::complex<long double> * cy, mplapackint incy);
std::complex<long double> Cdotu(mplapackint n, std::complex<long double> * cx, mplapackint incx, std::complex<long double> * cy, mplapackint incy);
long double RCnrm2(mplapackint n, std::complex<long double> * x, mplapackint incx);
long double Rnrm2(mplapackint n, long double * x, mplapackint incx);
long double RCasum(mplapackint n, std::complex<long double> * zx, mplapackint incx);
long double Rasum(mplapackint n, long double * dx, mplapackint incx);
mplapackint iCamax(mplapackint n, std::complex<long double> * dx, mplapackint incx);
mplapackint iRamax(mplapackint n, long double * dx, mplapackint incx);
long double RCabs1(std::complex<long double> z);
mplapackint Mlsame_longdouble(const char *a, const char *b);
void Mxerbla_longdouble(const char *srname, int info);

/* LEVEL 2 MPBLAS */
void Cgemv(const char *trans, mplapackint m, mplapackint n, std::complex<long double> alpha,
	   std::complex<long double> * A, mplapackint lda, std::complex<long double> * x, mplapackint incx, std::complex<long double> beta, std::complex<long double> * y, mplapackint incy);
void Rgemv(const char *trans, mplapackint m, mplapackint n, long double alpha,
	   long double * A, mplapackint lda, long double * x, mplapackint incx, long double beta, long double * y, mplapackint incy);
void Cgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   std::complex<long double> alpha, std::complex<long double> * A, mplapackint lda, std::complex<long double> * x, mplapackint incx, std::complex<long double> beta, std::complex<long double> * y, mplapackint incy);
void Rgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   long double alpha, long double * A, mplapackint lda, long double * x, mplapackint incx, long double beta, long double * y, mplapackint incy);
void Chemv(const char *uplo, mplapackint n, std::complex<long double> alpha, std::complex<long double> * A,
	   mplapackint lda, std::complex<long double> * x, mplapackint incx, std::complex<long double> beta, std::complex<long double> * y, mplapackint incy);
void Chbmv(const char *uplo, mplapackint n, mplapackint k, std::complex<long double> alpha,
	   std::complex<long double> * A, mplapackint lda, std::complex<long double> * x, mplapackint incx, std::complex<long double> beta, std::complex<long double> * y, mplapackint incy);
void Chpmv(const char *uplo, mplapackint n, std::complex<long double> alpha, std::complex<long double> * AP,
	   std::complex<long double> * x, mplapackint incx, std::complex<long double> beta, std::complex<long double> * y, mplapackint incy);
void Rsymv(const char *uplo, mplapackint n, long double alpha, long double * A,
	   mplapackint lda, long double * x, mplapackint incx, long double beta, long double * y, mplapackint incy);
void Rsbmv(const char *uplo, mplapackint n, mplapackint k, long double alpha,
	   long double * A, mplapackint lda, long double * x, mplapackint incx, long double beta, long double * y, mplapackint incy);
void Rspmv(const char *uplo, mplapackint n, long double alpha, long double * AP, long double * x, mplapackint incx, long double beta, long double * y, mplapackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<long double> * A, mplapackint lda, std::complex<long double> * x, mplapackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, long double * A, mplapackint lda, long double * x, mplapackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, std::complex<long double> * A, mplapackint lda, std::complex<long double> * x, mplapackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, long double * A, mplapackint lda, long double * x, mplapackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<long double> * AP, std::complex<long double> * x, mplapackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, long double * AP, long double * x, mplapackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<long double> * A, mplapackint lda, std::complex<long double> * x, mplapackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, long double * A, mplapackint lda, long double * x, mplapackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, std::complex<long double> * A, mplapackint lda, std::complex<long double> * x, mplapackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, long double * A, mplapackint lda, long double * x, mplapackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<long double> * AP, std::complex<long double> * x, mplapackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, long double * AP, long double * x, mplapackint incx);
void Rger(mplapackint m, mplapackint n, long double alpha, long double * x, mplapackint incx, long double * y, mplapackint incy, long double * A, mplapackint lda);
void Cgeru(mplapackint m, mplapackint n, std::complex<long double> alpha, std::complex<long double> * x, mplapackint incx, std::complex<long double> * y, mplapackint incy, std::complex<long double> * A, mplapackint lda);
void Cgerc(mplapackint m, mplapackint n, std::complex<long double> alpha, std::complex<long double> * x, mplapackint incx, std::complex<long double> * y, mplapackint incy, std::complex<long double> * A, mplapackint lda);
void Cher(const char *uplo, mplapackint n, long double alpha, std::complex<long double> * x, mplapackint incx, std::complex<long double> * A, mplapackint lda);
void Chpr(const char *uplo, mplapackint n, long double alpha, std::complex<long double> * x, mplapackint incx, std::complex<long double> * AP);
void Cher2(const char *uplo, mplapackint n, std::complex<long double> alpha, std::complex<long double> * x,
	   mplapackint incx, std::complex<long double> * y, mplapackint incy, std::complex<long double> * A, mplapackint lda);
void Chpr2(const char *uplo, mplapackint n, std::complex<long double> alpha, std::complex<long double> * x, mplapackint incx, std::complex<long double> * y, mplapackint incy, std::complex<long double> * AP);
void Rsyr(const char *uplo, mplapackint n, long double alpha, long double * x, mplapackint incx, long double * A, mplapackint lda);
void Rspr(const char *uplo, mplapackint n, long double alpha, long double * x, mplapackint incx, long double * AP);
void Rsyr2(const char *uplo, mplapackint n, long double alpha, long double * x, mplapackint incx, long double * y, mplapackint incy, long double * A, mplapackint lda);
void Rspr2(const char *uplo, mplapackint n, long double alpha, long double * x, mplapackint incx, long double * y, mplapackint incy, long double * AP);

/* LEVEL 3 MPBLAS */
void Cgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, std::complex<long double> alpha, std::complex<long double> * A, mplapackint lda, std::complex<long double> * B, mplapackint ldb, std::complex<long double> beta, std::complex<long double> * C, mplapackint ldc);
void Rgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, long double alpha, long double * A, mplapackint lda, long double * B, mplapackint ldb, long double beta, long double * C, mplapackint ldc);
void Csymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   std::complex<long double> alpha, std::complex<long double> * A, mplapackint lda, std::complex<long double> * B, mplapackint ldb, std::complex<long double> beta, std::complex<long double> * C, mplapackint ldc);
void Rsymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   long double alpha, long double * A, mplapackint lda, long double * B, mplapackint ldb, long double beta, long double * C, mplapackint ldc);
void Chemm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   std::complex<long double> alpha, std::complex<long double> * A, mplapackint lda, std::complex<long double> * B, mplapackint ldb, std::complex<long double> beta, std::complex<long double> * C, mplapackint ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   std::complex<long double> alpha, std::complex<long double> * A, mplapackint lda, std::complex<long double> beta, std::complex<long double> * C, mplapackint ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   long double alpha, long double * A, mplapackint lda, long double beta, long double * C, mplapackint ldc);
void Cherk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   long double alpha, std::complex<long double> * A, mplapackint lda, long double beta, std::complex<long double> * C, mplapackint ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    std::complex<long double> alpha, std::complex<long double> * A, mplapackint lda, std::complex<long double> * B, mplapackint ldb, std::complex<long double> beta, std::complex<long double> * C, mplapackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    long double alpha, long double * A, mplapackint lda, long double * B, mplapackint ldb, long double beta, long double * C, mplapackint ldc);
void Cher2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    std::complex<long double> alpha, std::complex<long double> * A, mplapackint lda, std::complex<long double> * B, mplapackint ldb, long double beta, std::complex<long double> * C, mplapackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, std::complex<long double> alpha, std::complex<long double> * A, mplapackint lda, std::complex<long double> * B, mplapackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, long double alpha, long double * A, mplapackint lda, long double * B, mplapackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, std::complex<long double> alpha, std::complex<long double> * A, mplapackint lda, std::complex<long double> * B, mplapackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, long double alpha, long double * A, mplapackint lda, long double * B, mplapackint ldb);

#endif
