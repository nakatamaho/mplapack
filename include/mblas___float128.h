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

#ifndef _MBLAS___FLOAT128_H_
#define _MBLAS___FLOAT128_H_

#include <mutils___float128.h>
#include <mpack_config.h>

#if !defined __MPACK_ERRNO__
#define _MPACK_EXTERN_ extern
#else
#define _MPACK_EXTERN_
#endif
_MPACK_EXTERN_ int mpack_errno;

/* LEVEL 1 MBLAS */
void Rrot(mpackint n, __float128 * dx, mpackint incx, __float128 * dy, mpackint incy, __float128 c, __float128 s);
void Rrotm(mpackint n, __float128 * x, mpackint incx, __float128 * y, mpackint incy, __float128 * dparam);
void CRrot(mpackint n, std::complex<__float128> * cx, mpackint incx, std::complex<__float128> * cy, mpackint incy, __float128 c, __float128 s);
void Cswap(mpackint n, std::complex<__float128> * cx, mpackint incx, std::complex<__float128> * cy, mpackint incy);
void Rswap(mpackint n, __float128 * dx, mpackint incx, __float128 * dy, mpackint incy);
void Cscal(mpackint n, std::complex<__float128> ca, std::complex<__float128> * cx, mpackint incx);
void CRscal(mpackint n, __float128 sa, std::complex<__float128> * cx, mpackint incx);
void Rscal(mpackint n, __float128 ca, __float128 * cx, mpackint incx);
void Ccopy(mpackint n, std::complex<__float128> * cx, mpackint incx, std::complex<__float128> * cy, mpackint incy);
void Crotg(std::complex<__float128> * ca, std::complex<__float128> cb, __float128 * c, std::complex<__float128> * s);
void Rrotg(__float128 * da, __float128 * db, __float128 * c, __float128 * s);
void Rcopy(mpackint n, __float128 * dx, mpackint incx, __float128 * dy, mpackint incy);
void Raxpy(mpackint n, __float128 da, __float128 * dx, mpackint incx, __float128 * dy, mpackint incy);
void Caxpy(mpackint n, std::complex<__float128> ca, std::complex<__float128> * cx, mpackint incx, std::complex<__float128> * cy, mpackint incy);
__float128 Rdot(mpackint n, __float128 * dx, mpackint incx, __float128 * dy, mpackint incy);
std::complex<__float128> Cdotc(mpackint n, std::complex<__float128> * cx, mpackint incx, std::complex<__float128> * cy, mpackint incy);
std::complex<__float128> Cdotu(mpackint n, std::complex<__float128> * cx, mpackint incx, std::complex<__float128> * cy, mpackint incy);
__float128 RCnrm2(mpackint n, std::complex<__float128> * x, mpackint incx);
__float128 Rnrm2(mpackint n, __float128 * x, mpackint incx);
__float128 RCasum(mpackint n, std::complex<__float128> * zx, mpackint incx);
__float128 Rasum(mpackint n, __float128 * dx, mpackint incx);
mpackint iCamax(mpackint n, std::complex<__float128> * dx, mpackint incx);
mpackint iRamax(mpackint n, __float128 * dx, mpackint incx);
__float128 RCabs1(std::complex<__float128> z);
mpackint Mlsame___float128(const char *a, const char *b);
void Mxerbla___float128(const char *srname, int info);

/* LEVEL 2 MBLAS */
void Cgemv(const char *trans, mpackint m, mpackint n, std::complex<__float128> alpha,
	   std::complex<__float128> * A, mpackint lda, std::complex<__float128> * x, mpackint incx, std::complex<__float128> beta, std::complex<__float128> * y, mpackint incy);
void Rgemv(const char *trans, mpackint m, mpackint n, __float128 alpha,
	   __float128 * A, mpackint lda, __float128 * x, mpackint incx, __float128 beta, __float128 * y, mpackint incy);
void Cgbmv(const char *trans, mpackint m, mpackint n, mpackint kl, mpackint ku,
	   std::complex<__float128> alpha, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * x, mpackint incx, std::complex<__float128> beta, std::complex<__float128> * y, mpackint incy);
void Rgbmv(const char *trans, mpackint m, mpackint n, mpackint kl, mpackint ku,
	   __float128 alpha, __float128 * A, mpackint lda, __float128 * x, mpackint incx, __float128 beta, __float128 * y, mpackint incy);
void Chemv(const char *uplo, mpackint n, std::complex<__float128> alpha, std::complex<__float128> * A,
	   mpackint lda, std::complex<__float128> * x, mpackint incx, std::complex<__float128> beta, std::complex<__float128> * y, mpackint incy);
void Chbmv(const char *uplo, mpackint n, mpackint k, std::complex<__float128> alpha,
	   std::complex<__float128> * A, mpackint lda, std::complex<__float128> * x, mpackint incx, std::complex<__float128> beta, std::complex<__float128> * y, mpackint incy);
void Chpmv(const char *uplo, mpackint n, std::complex<__float128> alpha, std::complex<__float128> * AP,
	   std::complex<__float128> * x, mpackint incx, std::complex<__float128> beta, std::complex<__float128> * y, mpackint incy);
void Rsymv(const char *uplo, mpackint n, __float128 alpha, __float128 * A,
	   mpackint lda, __float128 * x, mpackint incx, __float128 beta, __float128 * y, mpackint incy);
void Rsbmv(const char *uplo, mpackint n, mpackint k, __float128 alpha,
	   __float128 * A, mpackint lda, __float128 * x, mpackint incx, __float128 beta, __float128 * y, mpackint incy);
void Rspmv(const char *uplo, mpackint n, __float128 alpha, __float128 * AP, __float128 * x, mpackint incx, __float128 beta, __float128 * y, mpackint incy);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mpackint n, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * x, mpackint incx);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mpackint n, __float128 * A, mpackint lda, __float128 * x, mpackint incx);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * x, mpackint incx);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, __float128 * A, mpackint lda, __float128 * x, mpackint incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mpackint n, std::complex<__float128> * AP, std::complex<__float128> * x, mpackint incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mpackint n, __float128 * AP, __float128 * x, mpackint incx);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mpackint n, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * x, mpackint incx);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mpackint n, __float128 * A, mpackint lda, __float128 * x, mpackint incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * x, mpackint incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mpackint n,
	   mpackint k, __float128 * A, mpackint lda, __float128 * x, mpackint incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mpackint n, std::complex<__float128> * AP, std::complex<__float128> * x, mpackint incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mpackint n, __float128 * AP, __float128 * x, mpackint incx);
void Rger(mpackint m, mpackint n, __float128 alpha, __float128 * x, mpackint incx, __float128 * y, mpackint incy, __float128 * A, mpackint lda);
void Cgeru(mpackint m, mpackint n, std::complex<__float128> alpha, std::complex<__float128> * x, mpackint incx, std::complex<__float128> * y, mpackint incy, std::complex<__float128> * A, mpackint lda);
void Cgerc(mpackint m, mpackint n, std::complex<__float128> alpha, std::complex<__float128> * x, mpackint incx, std::complex<__float128> * y, mpackint incy, std::complex<__float128> * A, mpackint lda);
void Cher(const char *uplo, mpackint n, __float128 alpha, std::complex<__float128> * x, mpackint incx, std::complex<__float128> * A, mpackint lda);
void Chpr(const char *uplo, mpackint n, __float128 alpha, std::complex<__float128> * x, mpackint incx, std::complex<__float128> * AP);
void Cher2(const char *uplo, mpackint n, std::complex<__float128> alpha, std::complex<__float128> * x,
	   mpackint incx, std::complex<__float128> * y, mpackint incy, std::complex<__float128> * A, mpackint lda);
void Chpr2(const char *uplo, mpackint n, std::complex<__float128> alpha, std::complex<__float128> * x, mpackint incx, std::complex<__float128> * y, mpackint incy, std::complex<__float128> * AP);
void Rsyr(const char *uplo, mpackint n, __float128 alpha, __float128 * x, mpackint incx, __float128 * A, mpackint lda);
void Rspr(const char *uplo, mpackint n, __float128 alpha, __float128 * x, mpackint incx, __float128 * AP);
void Rsyr2(const char *uplo, mpackint n, __float128 alpha, __float128 * x, mpackint incx, __float128 * y, mpackint incy, __float128 * A, mpackint lda);
void Rspr2(const char *uplo, mpackint n, __float128 alpha, __float128 * x, mpackint incx, __float128 * y, mpackint incy, __float128 * AP);

/* LEVEL 3 MBLAS */
void Cgemm(const char *transa, const char *transb, mpackint m, mpackint n,
	   mpackint k, std::complex<__float128> alpha, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * B, mpackint ldb, std::complex<__float128> beta, std::complex<__float128> * C, mpackint ldc);
void Rgemm(const char *transa, const char *transb, mpackint m, mpackint n,
	   mpackint k, __float128 alpha, __float128 * A, mpackint lda, __float128 * B, mpackint ldb, __float128 beta, __float128 * C, mpackint ldc);
void Csymm(const char *side, const char *uplo, mpackint m, mpackint n,
	   std::complex<__float128> alpha, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * B, mpackint ldb, std::complex<__float128> beta, std::complex<__float128> * C, mpackint ldc);
void Rsymm(const char *side, const char *uplo, mpackint m, mpackint n,
	   __float128 alpha, __float128 * A, mpackint lda, __float128 * B, mpackint ldb, __float128 beta, __float128 * C, mpackint ldc);
void Chemm(const char *side, const char *uplo, mpackint m, mpackint n,
	   std::complex<__float128> alpha, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * B, mpackint ldb, std::complex<__float128> beta, std::complex<__float128> * C, mpackint ldc);
void Csyrk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   std::complex<__float128> alpha, std::complex<__float128> * A, mpackint lda, std::complex<__float128> beta, std::complex<__float128> * C, mpackint ldc);
void Rsyrk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   __float128 alpha, __float128 * A, mpackint lda, __float128 beta, __float128 * C, mpackint ldc);
void Cherk(const char *uplo, const char *trans, mpackint n, mpackint k,
	   __float128 alpha, std::complex<__float128> * A, mpackint lda, __float128 beta, std::complex<__float128> * C, mpackint ldc);
void Csyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    std::complex<__float128> alpha, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * B, mpackint ldb, std::complex<__float128> beta, std::complex<__float128> * C, mpackint ldc);
void Rsyr2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    __float128 alpha, __float128 * A, mpackint lda, __float128 * B, mpackint ldb, __float128 beta, __float128 * C, mpackint ldc);
void Cher2k(const char *uplo, const char *trans, mpackint n, mpackint k,
	    std::complex<__float128> alpha, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * B, mpackint ldb, __float128 beta, std::complex<__float128> * C, mpackint ldc);
void Ctrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, std::complex<__float128> alpha, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * B, mpackint ldb);
void Rtrmm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, __float128 alpha, __float128 * A, mpackint lda, __float128 * B, mpackint ldb);
void Ctrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, std::complex<__float128> alpha, std::complex<__float128> * A, mpackint lda, std::complex<__float128> * B, mpackint ldb);
void Rtrsm(const char *side, const char *uplo, const char *transa,
	   const char *diag, mpackint m, mpackint n, __float128 alpha, __float128 * A, mpackint lda, __float128 * B, mpackint ldb);

#endif
