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

#ifndef _MPBLAS__FLOAT64X_H_
#define _MPBLAS__FLOAT64X_H_

#include "mplapack_utils_longdouble.h"
#include "mplapack_config.h"

#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

/* LEVEL 1 MPBLAS */
void Rrot(mplapackint n, _Float64x * dx, mplapackint incx, _Float64x * dy, mplapackint incy, _Float64x c, _Float64x s);
void Rrotm(mplapackint n, _Float64x * x, mplapackint incx, _Float64x * y, mplapackint incy, _Float64x * dparam);
void CRrot(mplapackint n, std::complex<_Float64x> * cx, mplapackint incx, std::complex<_Float64x> * cy, mplapackint incy, _Float64x c, _Float64x s);
void Cswap(mplapackint n, std::complex<_Float64x> * cx, mplapackint incx, std::complex<_Float64x> * cy, mplapackint incy);
void Rswap(mplapackint n, _Float64x * dx, mplapackint incx, _Float64x * dy, mplapackint incy);
void Cscal(mplapackint n, std::complex<_Float64x> ca, std::complex<_Float64x> * cx, mplapackint incx);
void CRscal(mplapackint n, _Float64x sa, std::complex<_Float64x> * cx, mplapackint incx);
void Rscal(mplapackint n, _Float64x ca, _Float64x * cx, mplapackint incx);
void Ccopy(mplapackint n, std::complex<_Float64x> * cx, mplapackint incx, std::complex<_Float64x> * cy, mplapackint incy);
void Crotg(std::complex<_Float64x> * ca, std::complex<_Float64x> cb, _Float64x * c, std::complex<_Float64x> * s);
void Rrotg(_Float64x * da, _Float64x * db, _Float64x * c, _Float64x * s);
void Rcopy(mplapackint n, _Float64x * dx, mplapackint incx, _Float64x * dy, mplapackint incy);
void Raxpy(mplapackint n, _Float64x da, _Float64x * dx, mplapackint incx, _Float64x * dy, mplapackint incy);
void Caxpy(mplapackint n, std::complex<_Float64x> ca, std::complex<_Float64x> * cx, mplapackint incx, std::complex<_Float64x> * cy, mplapackint incy);
_Float64x Rdot(mplapackint n, _Float64x * dx, mplapackint incx, _Float64x * dy, mplapackint incy);
std::complex<_Float64x> Cdotc(mplapackint n, std::complex<_Float64x> * cx, mplapackint incx, std::complex<_Float64x> * cy, mplapackint incy);
std::complex<_Float64x> Cdotu(mplapackint n, std::complex<_Float64x> * cx, mplapackint incx, std::complex<_Float64x> * cy, mplapackint incy);
_Float64x RCnrm2(mplapackint n, std::complex<_Float64x> * x, mplapackint incx);
_Float64x Rnrm2(mplapackint n, _Float64x * x, mplapackint incx);
_Float64x RCasum(mplapackint n, std::complex<_Float64x> * zx, mplapackint incx);
_Float64x Rasum(mplapackint n, _Float64x * dx, mplapackint incx);
mplapackint iCamax(mplapackint n, std::complex<_Float64x> * dx, mplapackint incx);
mplapackint iRamax(mplapackint n, _Float64x * dx, mplapackint incx);
_Float64x RCabs1(std::complex<_Float64x> z);
mplapackint Mlsame_longdouble(const char *a, const char *b);
void Mxerbla_longdouble(const char *srname, int info);

/* LEVEL 2 MPBLAS */
void Cgemv(const char *trans, mplapackint m, mplapackint n, std::complex<_Float64x> alpha,
	   std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * x, mplapackint incx, std::complex<_Float64x> beta, std::complex<_Float64x> * y, mplapackint incy);
void Rgemv(const char *trans, mplapackint m, mplapackint n, _Float64x alpha,
	   _Float64x * A, mplapackint lda, _Float64x * x, mplapackint incx, _Float64x beta, _Float64x * y, mplapackint incy);
void Cgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   std::complex<_Float64x> alpha, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * x, mplapackint incx, std::complex<_Float64x> beta, std::complex<_Float64x> * y, mplapackint incy);
void Rgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   _Float64x alpha, _Float64x * A, mplapackint lda, _Float64x * x, mplapackint incx, _Float64x beta, _Float64x * y, mplapackint incy);
void Chemv(const char *uplo, mplapackint n, std::complex<_Float64x> alpha, std::complex<_Float64x> * A,
	   mplapackint lda, std::complex<_Float64x> * x, mplapackint incx, std::complex<_Float64x> beta, std::complex<_Float64x> * y, mplapackint incy);
void Chbmv(const char *uplo, mplapackint n, mplapackint k, std::complex<_Float64x> alpha,
	   std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * x, mplapackint incx, std::complex<_Float64x> beta, std::complex<_Float64x> * y, mplapackint incy);
void Chpmv(const char *uplo, mplapackint n, std::complex<_Float64x> alpha, std::complex<_Float64x> * AP,
	   std::complex<_Float64x> * x, mplapackint incx, std::complex<_Float64x> beta, std::complex<_Float64x> * y, mplapackint incy);
void Rsymv(const char *uplo, mplapackint n, _Float64x alpha, _Float64x * A,
	   mplapackint lda, _Float64x * x, mplapackint incx, _Float64x beta, _Float64x * y, mplapackint incy);
void Rsbmv(const char *uplo, mplapackint n, mplapackint k, _Float64x alpha,
	   _Float64x * A, mplapackint lda, _Float64x * x, mplapackint incx, _Float64x beta, _Float64x * y, mplapackint incy);
void Rspmv(const char *uplo, mplapackint n, _Float64x alpha, _Float64x * AP, _Float64x * x, mplapackint incx, _Float64x beta, _Float64x * y, mplapackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * x, mplapackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, _Float64x * A, mplapackint lda, _Float64x * x, mplapackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * x, mplapackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, _Float64x * A, mplapackint lda, _Float64x * x, mplapackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<_Float64x> * AP, std::complex<_Float64x> * x, mplapackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, _Float64x * AP, _Float64x * x, mplapackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * x, mplapackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, _Float64x * A, mplapackint lda, _Float64x * x, mplapackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * x, mplapackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, _Float64x * A, mplapackint lda, _Float64x * x, mplapackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<_Float64x> * AP, std::complex<_Float64x> * x, mplapackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, _Float64x * AP, _Float64x * x, mplapackint incx);
void Rger(mplapackint m, mplapackint n, _Float64x alpha, _Float64x * x, mplapackint incx, _Float64x * y, mplapackint incy, _Float64x * A, mplapackint lda);
void Cgeru(mplapackint m, mplapackint n, std::complex<_Float64x> alpha, std::complex<_Float64x> * x, mplapackint incx, std::complex<_Float64x> * y, mplapackint incy, std::complex<_Float64x> * A, mplapackint lda);
void Cgerc(mplapackint m, mplapackint n, std::complex<_Float64x> alpha, std::complex<_Float64x> * x, mplapackint incx, std::complex<_Float64x> * y, mplapackint incy, std::complex<_Float64x> * A, mplapackint lda);
void Cher(const char *uplo, mplapackint n, _Float64x alpha, std::complex<_Float64x> * x, mplapackint incx, std::complex<_Float64x> * A, mplapackint lda);
void Chpr(const char *uplo, mplapackint n, _Float64x alpha, std::complex<_Float64x> * x, mplapackint incx, std::complex<_Float64x> * AP);
void Cher2(const char *uplo, mplapackint n, std::complex<_Float64x> alpha, std::complex<_Float64x> * x,
	   mplapackint incx, std::complex<_Float64x> * y, mplapackint incy, std::complex<_Float64x> * A, mplapackint lda);
void Chpr2(const char *uplo, mplapackint n, std::complex<_Float64x> alpha, std::complex<_Float64x> * x, mplapackint incx, std::complex<_Float64x> * y, mplapackint incy, std::complex<_Float64x> * AP);
void Rsyr(const char *uplo, mplapackint n, _Float64x alpha, _Float64x * x, mplapackint incx, _Float64x * A, mplapackint lda);
void Rspr(const char *uplo, mplapackint n, _Float64x alpha, _Float64x * x, mplapackint incx, _Float64x * AP);
void Rsyr2(const char *uplo, mplapackint n, _Float64x alpha, _Float64x * x, mplapackint incx, _Float64x * y, mplapackint incy, _Float64x * A, mplapackint lda);
void Rspr2(const char *uplo, mplapackint n, _Float64x alpha, _Float64x * x, mplapackint incx, _Float64x * y, mplapackint incy, _Float64x * AP);

/* LEVEL 3 MPBLAS */
void Cgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, std::complex<_Float64x> alpha, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * B, mplapackint ldb, std::complex<_Float64x> beta, std::complex<_Float64x> * C, mplapackint ldc);
void Rgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, _Float64x alpha, _Float64x * A, mplapackint lda, _Float64x * B, mplapackint ldb, _Float64x beta, _Float64x * C, mplapackint ldc);
void Csymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   std::complex<_Float64x> alpha, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * B, mplapackint ldb, std::complex<_Float64x> beta, std::complex<_Float64x> * C, mplapackint ldc);
void Rsymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   _Float64x alpha, _Float64x * A, mplapackint lda, _Float64x * B, mplapackint ldb, _Float64x beta, _Float64x * C, mplapackint ldc);
void Chemm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   std::complex<_Float64x> alpha, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * B, mplapackint ldb, std::complex<_Float64x> beta, std::complex<_Float64x> * C, mplapackint ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   std::complex<_Float64x> alpha, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> beta, std::complex<_Float64x> * C, mplapackint ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   _Float64x alpha, _Float64x * A, mplapackint lda, _Float64x beta, _Float64x * C, mplapackint ldc);
void Cherk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   _Float64x alpha, std::complex<_Float64x> * A, mplapackint lda, _Float64x beta, std::complex<_Float64x> * C, mplapackint ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    std::complex<_Float64x> alpha, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * B, mplapackint ldb, std::complex<_Float64x> beta, std::complex<_Float64x> * C, mplapackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    _Float64x alpha, _Float64x * A, mplapackint lda, _Float64x * B, mplapackint ldb, _Float64x beta, _Float64x * C, mplapackint ldc);
void Cher2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    std::complex<_Float64x> alpha, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * B, mplapackint ldb, _Float64x beta, std::complex<_Float64x> * C, mplapackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, std::complex<_Float64x> alpha, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * B, mplapackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, _Float64x alpha, _Float64x * A, mplapackint lda, _Float64x * B, mplapackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, std::complex<_Float64x> alpha, std::complex<_Float64x> * A, mplapackint lda, std::complex<_Float64x> * B, mplapackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, _Float64x alpha, _Float64x * A, mplapackint lda, _Float64x * B, mplapackint ldb);

#endif
