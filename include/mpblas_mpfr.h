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

void CRrot(mplapackint const &n, mpcomplex *zx, mplapackint const &incx, mpcomplex *zy, mplapackint const &incy, mpreal const &c, mpreal const &s);
void CRscal(mplapackint const &n, mpreal const &da, mpcomplex *zx, mplapackint const &incx);
void Caxpy(mplapackint const &n, mpcomplex const &za, mpcomplex *zx, mplapackint const &incx, mpcomplex *zy, mplapackint const &incy);
void Ccopy(mplapackint const &n, mpcomplex *zx, mplapackint const &incx, mpcomplex *zy, mplapackint const &incy);
mpcomplex Cdotc(mplapackint const &n, mpcomplex *zx, mplapackint const &incx, mpcomplex *zy, mplapackint const &incy);
mpcomplex Cdotu(mplapackint const &n, mpcomplex *zx, mplapackint const &incx, mpcomplex *zy, mplapackint const &incy);
void Cgbmv(const char *trans, mplapackint const &m, mplapackint const &n, mplapackint const &kl, mplapackint const &ku, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *x, mplapackint const &incx, mpcomplex const &beta, mpcomplex *y, mplapackint const &incy);
void Cgemm(const char *transa, const char *transb, mplapackint const &m, mplapackint const &n, mplapackint const &k, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *b, mplapackint const &ldb, mpcomplex const &beta, mpcomplex *c, mplapackint const &ldc);
void Cgemv(const char *trans, mplapackint const &m, mplapackint const &n, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *x, mplapackint const &incx, mpcomplex const &beta, mpcomplex *y, mplapackint const &incy);
void Cgerc(mplapackint const &m, mplapackint const &n, mpcomplex const &alpha, mpcomplex *x, mplapackint const &incx, mpcomplex *y, mplapackint const &incy, mpcomplex *a, mplapackint const &lda);
void Cgeru(mplapackint const &m, mplapackint const &n, mpcomplex const &alpha, mpcomplex *x, mplapackint const &incx, mpcomplex *y, mplapackint const &incy, mpcomplex *a, mplapackint const &lda);
void Chbmv(const char *uplo, mplapackint const &n, mplapackint const &k, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *x, mplapackint const &incx, mpcomplex const &beta, mpcomplex *y, mplapackint const &incy);
void Chemm(const char *side, const char *uplo, mplapackint const &m, mplapackint const &n, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *b, mplapackint const &ldb, mpcomplex const &beta, mpcomplex *c, mplapackint const &ldc);
void Chemv(const char *uplo, mplapackint const &n, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *x, mplapackint const &incx, mpcomplex const &beta, mpcomplex *y, mplapackint const &incy);
void Cher(const char *uplo, mplapackint const &n, mpreal const &alpha, mpcomplex *x, mplapackint const &incx, mpcomplex *a, mplapackint const &lda);
void Cher2(const char *uplo, mplapackint const &n, mpcomplex const &alpha, mpcomplex *x, mplapackint const &incx, mpcomplex *y, mplapackint const &incy, mpcomplex *a, mplapackint const &lda);
void Cher2k(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *b, mplapackint const &ldb, mpreal const &beta, mpcomplex *c, mplapackint const &ldc);
void Cherk(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, mpreal const &alpha, mpcomplex *a, mplapackint const &lda, mpreal const &beta, mpcomplex *c, mplapackint const &ldc);
void Chpmv(const char *uplo, mplapackint const &n, mpcomplex const &alpha, mpcomplex *ap, mpcomplex *x, mplapackint const &incx, mpcomplex const &beta, mpcomplex *y, mplapackint const &incy);
void Chpr(const char *uplo, mplapackint const &n, mpreal const &alpha, mpcomplex *x, mplapackint const &incx, mpcomplex *ap);
void Chpr2(const char *uplo, mplapackint const &n, mpcomplex const &alpha, mpcomplex *x, mplapackint const &incx, mpcomplex *y, mplapackint const &incy, mpcomplex *ap);
void Crotg(mpcomplex &ca, mpcomplex const &cb, mpreal &c, mpcomplex &s);
void Cscal(mplapackint const &n, mpcomplex const &za, mpcomplex *zx, mplapackint const &incx);
void Cswap(mplapackint const &n, mpcomplex *zx, mplapackint const &incx, mpcomplex *zy, mplapackint const &incy);
void Csymm(const char *side, const char *uplo, mplapackint const &m, mplapackint const &n, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *b, mplapackint const &ldb, mpcomplex const &beta, mpcomplex *c, mplapackint const &ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *b, mplapackint const &ldb, mpcomplex const &beta, mpcomplex *c, mplapackint const &ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex const &beta, mpcomplex *c, mplapackint const &ldc);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mplapackint const &k, mpcomplex *a, mplapackint const &lda, mpcomplex *x, mplapackint const &incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mplapackint const &k, mpcomplex *a, mplapackint const &lda, mpcomplex *x, mplapackint const &incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mpcomplex *ap, mpcomplex *x, mplapackint const &incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mpcomplex *ap, mpcomplex *x, mplapackint const &incx);
void Ctrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const &m, mplapackint const &n, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *b, mplapackint const &ldb);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mpcomplex *a, mplapackint const &lda, mpcomplex *x, mplapackint const &incx);
void Ctrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const &m, mplapackint const &n, mpcomplex const &alpha, mpcomplex *a, mplapackint const &lda, mpcomplex *b, mplapackint const &ldb);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mpcomplex *a, mplapackint const &lda, mpcomplex *x, mplapackint const &incx);
mplapackint Mlsame_mpfr(const char *a, const char *b);
void Mxerbla_mpfr(const char *srname, int info);
mpreal RCabs1(mpcomplex const &z);
mpreal RCasum(mplapackint const &n, mpcomplex *zx, mplapackint const &incx);
mpreal RCnrm2(mplapackint const &n, mpcomplex *x, mplapackint const &incx);
mpreal Rasum(mplapackint const &n, mpreal *dx, mplapackint const &incx);
void Raxpy(mplapackint const &n, mpreal const &da, mpreal *dx, mplapackint const &incx, mpreal *dy, mplapackint const &incy);
void Rcopy(mplapackint const &n, mpreal *dx, mplapackint const &incx, mpreal *dy, mplapackint const &incy);
mpreal Rdot(mplapackint const &n, mpreal *dx, mplapackint const &incx, mpreal *dy, mplapackint const &incy);
void Rgbmv(const char *trans, mplapackint const &m, mplapackint const &n, mplapackint const &kl, mplapackint const &ku, mpreal const &alpha, mpreal *a, mplapackint const &lda, mpreal *x, mplapackint const &incx, mpreal const &beta, mpreal *y, mplapackint const &incy);
void Rgemm(const char *transa, const char *transb, mplapackint const &m, mplapackint const &n, mplapackint const &k, mpreal const &alpha, mpreal *a, mplapackint const &lda, mpreal *b, mplapackint const &ldb, mpreal const &beta, mpreal *c, mplapackint const &ldc);
void Rgemv(const char *trans, mplapackint const &m, mplapackint const &n, mpreal const &alpha, mpreal *a, mplapackint const &lda, mpreal *x, mplapackint const &incx, mpreal const &beta, mpreal *y, mplapackint const &incy);
void Rger(mplapackint const &m, mplapackint const &n, mpreal const &alpha, mpreal *x, mplapackint const &incx, mpreal *y, mplapackint const &incy, mpreal *a, mplapackint const &lda);
mpreal Rnrm2(mplapackint const &n, mpreal *x, mplapackint const &incx);
void Rrot(mplapackint const &n, mpreal *dx, mplapackint const &incx, mpreal *dy, mplapackint const &incy, mpreal const &c, mpreal const &s);
void Rrotg(mpreal &da, mpreal &db, mpreal &c, mpreal &s);
void Rrotm(mplapackint const &n, mpreal *dx, mplapackint const &incx, mpreal *dy, mplapackint const &incy, mpreal *dparam);
void Rrotmg(mpreal &dd1, mpreal &dd2, mpreal &dx1, mpreal const &dy1, mpreal *dparam);
void Rsbmv(const char *uplo, mplapackint const &n, mplapackint const &k, mpreal const &alpha, mpreal *a, mplapackint const &lda, mpreal *x, mplapackint const &incx, mpreal const &beta, mpreal *y, mplapackint const &incy);
void Rscal(mplapackint const &n, mpreal const &da, mpreal *dx, mplapackint const &incx);
void Rspmv(const char *uplo, mplapackint const &n, mpreal const &alpha, mpreal *ap, mpreal *x, mplapackint const &incx, mpreal const &beta, mpreal *y, mplapackint const &incy);
void Rspr(const char *uplo, mplapackint const &n, mpreal const &alpha, mpreal *x, mplapackint const &incx, mpreal *ap);
void Rspr2(const char *uplo, mplapackint const &n, mpreal const &alpha, mpreal *x, mplapackint const &incx, mpreal *y, mplapackint const &incy, mpreal *ap);
void Rswap(mplapackint const &n, mpreal *dx, mplapackint const &incx, mpreal *dy, mplapackint const &incy);
void Rsymm(const char *side, const char *uplo, mplapackint const &m, mplapackint const &n, mpreal const &alpha, mpreal *a, mplapackint const &lda, mpreal *b, mplapackint const &ldb, mpreal const &beta, mpreal *c, mplapackint const &ldc);
void Rsymv(const char *uplo, mplapackint const &n, mpreal const &alpha, mpreal *a, mplapackint const &lda, mpreal *x, mplapackint const &incx, mpreal const &beta, mpreal *y, mplapackint const &incy);
void Rsyr(const char *uplo, mplapackint const &n, mpreal const &alpha, mpreal *x, mplapackint const &incx, mpreal *a, mplapackint const &lda);
void Rsyr2(const char *uplo, mplapackint const &n, mpreal const &alpha, mpreal *x, mplapackint const &incx, mpreal *y, mplapackint const &incy, mpreal *a, mplapackint const &lda);
void Rsyr2k(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, mpreal const &alpha, mpreal *a, mplapackint const &lda, mpreal *b, mplapackint const &ldb, mpreal const &beta, mpreal *c, mplapackint const &ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, mpreal const &alpha, mpreal *a, mplapackint const &lda, mpreal const &beta, mpreal *c, mplapackint const &ldc);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mplapackint const &k, mpreal *a, mplapackint const &lda, mpreal *x, mplapackint const &incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mplapackint const &k, mpreal *a, mplapackint const &lda, mpreal *x, mplapackint const &incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mpreal *ap, mpreal *x, mplapackint const &incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mpreal *ap, mpreal *x, mplapackint const &incx);
void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const &m, mplapackint const &n, mpreal const &alpha, mpreal *a, mplapackint const &lda, mpreal *b, mplapackint const &ldb);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mpreal *a, mplapackint const &lda, mpreal *x, mplapackint const &incx);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const &m, mplapackint const &n, mpreal const &alpha, mpreal *a, mplapackint const &lda, mpreal *b, mplapackint const &ldb);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mpreal *a, mplapackint const &lda, mpreal *x, mplapackint const &incx);
mplapackint iCamax(mplapackint const &n, mpcomplex *zx, mplapackint const &incx);
mplapackint iRamax(mplapackint const &n, mpreal *dx, mplapackint const &incx);
#endif
