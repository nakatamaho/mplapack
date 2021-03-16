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

#ifndef _MPBLAS___FLOAT128_H_
#define _MPBLAS___FLOAT128_H_

#include "mutils___float128.h"
#include "mplapack_config.h"

#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

/* LEVEL 1 MPBLAS */
void Rrot(mplapackint n, __float128 * dx, mplapackint incx, __float128 * dy, mplapackint incy, __float128 c, __float128 s);
void Rrotm(mplapackint n, __float128 * x, mplapackint incx, __float128 * y, mplapackint incy, __float128 * dparam);
void CRrot(mplapackint n, std::complex<__float128> * cx, mplapackint incx, std::complex<__float128> * cy, mplapackint incy, __float128 c, __float128 s);
void Cswap(mplapackint n, std::complex<__float128> * cx, mplapackint incx, std::complex<__float128> * cy, mplapackint incy);
void Rswap(mplapackint n, __float128 * dx, mplapackint incx, __float128 * dy, mplapackint incy);
void Cscal(mplapackint n, std::complex<__float128> ca, std::complex<__float128> * cx, mplapackint incx);
void CRscal(mplapackint n, __float128 sa, std::complex<__float128> * cx, mplapackint incx);
void Rscal(mplapackint n, __float128 ca, __float128 * cx, mplapackint incx);
void Ccopy(mplapackint n, std::complex<__float128> * cx, mplapackint incx, std::complex<__float128> * cy, mplapackint incy);
void Crotg(std::complex<__float128> * ca, std::complex<__float128> cb, __float128 * c, std::complex<__float128> * s);
void Rrotg(__float128 * da, __float128 * db, __float128 * c, __float128 * s);
void Rcopy(mplapackint n, __float128 * dx, mplapackint incx, __float128 * dy, mplapackint incy);
void Raxpy(mplapackint n, __float128 da, __float128 * dx, mplapackint incx, __float128 * dy, mplapackint incy);
void Caxpy(mplapackint n, std::complex<__float128> ca, std::complex<__float128> * cx, mplapackint incx, std::complex<__float128> * cy, mplapackint incy);
__float128 Rdot(mplapackint n, __float128 * dx, mplapackint incx, __float128 * dy, mplapackint incy);
std::complex<__float128> Cdotc(mplapackint n, std::complex<__float128> * cx, mplapackint incx, std::complex<__float128> * cy, mplapackint incy);
std::complex<__float128> Cdotu(mplapackint n, std::complex<__float128> * cx, mplapackint incx, std::complex<__float128> * cy, mplapackint incy);
__float128 RCnrm2(mplapackint n, std::complex<__float128> * x, mplapackint incx);
__float128 Rnrm2(mplapackint n, __float128 * x, mplapackint incx);
__float128 RCasum(mplapackint n, std::complex<__float128> * zx, mplapackint incx);
__float128 Rasum(mplapackint n, __float128 * dx, mplapackint incx);
mplapackint iCamax(mplapackint n, std::complex<__float128> * dx, mplapackint incx);
mplapackint iRamax(mplapackint n, __float128 * dx, mplapackint incx);
__float128 RCabs1(std::complex<__float128> z);
mplapackint Mlsame___float128(const char *a, const char *b);
void Mxerbla___float128(const char *srname, int info);

/* LEVEL 2 MPBLAS */
void Cgemv(const char *trans, mplapackint m, mplapackint n, std::complex<__float128> alpha,
	   std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * x, mplapackint incx, std::complex<__float128> beta, std::complex<__float128> * y, mplapackint incy);
void Rgemv(const char *trans, mplapackint m, mplapackint n, __float128 alpha,
	   __float128 * A, mplapackint lda, __float128 * x, mplapackint incx, __float128 beta, __float128 * y, mplapackint incy);
void Cgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   std::complex<__float128> alpha, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * x, mplapackint incx, std::complex<__float128> beta, std::complex<__float128> * y, mplapackint incy);
void Rgbmv(const char *trans, mplapackint m, mplapackint n, mplapackint kl, mplapackint ku,
	   __float128 alpha, __float128 * A, mplapackint lda, __float128 * x, mplapackint incx, __float128 beta, __float128 * y, mplapackint incy);
void Chemv(const char *uplo, mplapackint n, std::complex<__float128> alpha, std::complex<__float128> * A,
	   mplapackint lda, std::complex<__float128> * x, mplapackint incx, std::complex<__float128> beta, std::complex<__float128> * y, mplapackint incy);
void Chbmv(const char *uplo, mplapackint n, mplapackint k, std::complex<__float128> alpha,
	   std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * x, mplapackint incx, std::complex<__float128> beta, std::complex<__float128> * y, mplapackint incy);
void Chpmv(const char *uplo, mplapackint n, std::complex<__float128> alpha, std::complex<__float128> * AP,
	   std::complex<__float128> * x, mplapackint incx, std::complex<__float128> beta, std::complex<__float128> * y, mplapackint incy);
void Rsymv(const char *uplo, mplapackint n, __float128 alpha, __float128 * A,
	   mplapackint lda, __float128 * x, mplapackint incx, __float128 beta, __float128 * y, mplapackint incy);
void Rsbmv(const char *uplo, mplapackint n, mplapackint k, __float128 alpha,
	   __float128 * A, mplapackint lda, __float128 * x, mplapackint incx, __float128 beta, __float128 * y, mplapackint incy);
void Rspmv(const char *uplo, mplapackint n, __float128 alpha, __float128 * AP, __float128 * x, mplapackint incx, __float128 beta, __float128 * y, mplapackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * x, mplapackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint n, __float128 * A, mplapackint lda, __float128 * x, mplapackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * x, mplapackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, __float128 * A, mplapackint lda, __float128 * x, mplapackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<__float128> * AP, std::complex<__float128> * x, mplapackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint n, __float128 * AP, __float128 * x, mplapackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * x, mplapackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint n, __float128 * A, mplapackint lda, __float128 * x, mplapackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * x, mplapackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint n,
	   mplapackint k, __float128 * A, mplapackint lda, __float128 * x, mplapackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, std::complex<__float128> * AP, std::complex<__float128> * x, mplapackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint n, __float128 * AP, __float128 * x, mplapackint incx);
void Rger(mplapackint m, mplapackint n, __float128 alpha, __float128 * x, mplapackint incx, __float128 * y, mplapackint incy, __float128 * A, mplapackint lda);
void Cgeru(mplapackint m, mplapackint n, std::complex<__float128> alpha, std::complex<__float128> * x, mplapackint incx, std::complex<__float128> * y, mplapackint incy, std::complex<__float128> * A, mplapackint lda);
void Cgerc(mplapackint m, mplapackint n, std::complex<__float128> alpha, std::complex<__float128> * x, mplapackint incx, std::complex<__float128> * y, mplapackint incy, std::complex<__float128> * A, mplapackint lda);
void Cher(const char *uplo, mplapackint n, __float128 alpha, std::complex<__float128> * x, mplapackint incx, std::complex<__float128> * A, mplapackint lda);
void Chpr(const char *uplo, mplapackint n, __float128 alpha, std::complex<__float128> * x, mplapackint incx, std::complex<__float128> * AP);
void Cher2(const char *uplo, mplapackint n, std::complex<__float128> alpha, std::complex<__float128> * x,
	   mplapackint incx, std::complex<__float128> * y, mplapackint incy, std::complex<__float128> * A, mplapackint lda);
void Chpr2(const char *uplo, mplapackint n, std::complex<__float128> alpha, std::complex<__float128> * x, mplapackint incx, std::complex<__float128> * y, mplapackint incy, std::complex<__float128> * AP);
void Rsyr(const char *uplo, mplapackint n, __float128 alpha, __float128 * x, mplapackint incx, __float128 * A, mplapackint lda);
void Rspr(const char *uplo, mplapackint n, __float128 alpha, __float128 * x, mplapackint incx, __float128 * AP);
void Rsyr2(const char *uplo, mplapackint n, __float128 alpha, __float128 * x, mplapackint incx, __float128 * y, mplapackint incy, __float128 * A, mplapackint lda);
void Rspr2(const char *uplo, mplapackint n, __float128 alpha, __float128 * x, mplapackint incx, __float128 * y, mplapackint incy, __float128 * AP);

/* LEVEL 3 MPBLAS */
void Cgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, std::complex<__float128> alpha, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * B, mplapackint ldb, std::complex<__float128> beta, std::complex<__float128> * C, mplapackint ldc);
void Rgemm(const char *transa, const char *transb, mplapackint m, mplapackint n,
	   mplapackint k, __float128 alpha, __float128 * A, mplapackint lda, __float128 * B, mplapackint ldb, __float128 beta, __float128 * C, mplapackint ldc);
void Csymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   std::complex<__float128> alpha, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * B, mplapackint ldb, std::complex<__float128> beta, std::complex<__float128> * C, mplapackint ldc);
void Rsymm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   __float128 alpha, __float128 * A, mplapackint lda, __float128 * B, mplapackint ldb, __float128 beta, __float128 * C, mplapackint ldc);
void Chemm(const char *side, const char *uplo, mplapackint m, mplapackint n,
	   std::complex<__float128> alpha, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * B, mplapackint ldb, std::complex<__float128> beta, std::complex<__float128> * C, mplapackint ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   std::complex<__float128> alpha, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> beta, std::complex<__float128> * C, mplapackint ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   __float128 alpha, __float128 * A, mplapackint lda, __float128 beta, __float128 * C, mplapackint ldc);
void Cherk(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	   __float128 alpha, std::complex<__float128> * A, mplapackint lda, __float128 beta, std::complex<__float128> * C, mplapackint ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    std::complex<__float128> alpha, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * B, mplapackint ldb, std::complex<__float128> beta, std::complex<__float128> * C, mplapackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    __float128 alpha, __float128 * A, mplapackint lda, __float128 * B, mplapackint ldb, __float128 beta, __float128 * C, mplapackint ldc);
void Cher2k(const char *uplo, const char *trans, mplapackint n, mplapackint k,
	    std::complex<__float128> alpha, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * B, mplapackint ldb, __float128 beta, std::complex<__float128> * C, mplapackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, std::complex<__float128> alpha, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * B, mplapackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, __float128 alpha, __float128 * A, mplapackint lda, __float128 * B, mplapackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, std::complex<__float128> alpha, std::complex<__float128> * A, mplapackint lda, std::complex<__float128> * B, mplapackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mplapackint m, mplapackint n, __float128 alpha, __float128 * A, mplapackint lda, __float128 * B, mplapackint ldb);

#endif
