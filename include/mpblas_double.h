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
#include "mplapack_utils_double.h"
#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

void CRrot(mplapackint const &n, std::complex<double> *zx, mplapackint const &incx, std::complex<double> *zy, mplapackint const &incy, double const &c, double const &s);
void CRscal(mplapackint const &n, double const &da, std::complex<double> *zx, mplapackint const &incx);
void Caxpy(mplapackint const &n, std::complex<double> const &za, std::complex<double> *zx, mplapackint const &incx, std::complex<double> *zy, mplapackint const &incy);
void Ccopy(mplapackint const &n, std::complex<double> *zx, mplapackint const &incx, std::complex<double> *zy, mplapackint const &incy);
std::complex<double> Cdotc(mplapackint const &n, std::complex<double> *zx, mplapackint const &incx, std::complex<double> *zy, mplapackint const &incy);
std::complex<double> Cdotu(mplapackint const &n, std::complex<double> *zx, mplapackint const &incx, std::complex<double> *zy, mplapackint const &incy);
void Cgbmv(const char *trans, mplapackint const &m, mplapackint const &n, mplapackint const &kl, mplapackint const &ku, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *x, mplapackint const &incx, std::complex<double> const &beta, std::complex<double> *y, mplapackint const &incy);
void Cgemm(const char *transa, const char *transb, mplapackint const &m, mplapackint const &n, mplapackint const &k, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *b, mplapackint const &ldb, std::complex<double> const &beta, std::complex<double> *c, mplapackint const &ldc);
void Cgemv(const char *trans, mplapackint const &m, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *x, mplapackint const &incx, std::complex<double> const &beta, std::complex<double> *y, mplapackint const &incy);
void Cgerc(mplapackint const &m, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *x, mplapackint const &incx, std::complex<double> *y, mplapackint const &incy, std::complex<double> *a, mplapackint const &lda);
void Cgeru(mplapackint const &m, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *x, mplapackint const &incx, std::complex<double> *y, mplapackint const &incy, std::complex<double> *a, mplapackint const &lda);
void Chbmv(const char *uplo, mplapackint const &n, mplapackint const &k, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *x, mplapackint const &incx, std::complex<double> const &beta, std::complex<double> *y, mplapackint const &incy);
void Chemm(const char *side, const char *uplo, mplapackint const &m, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *b, mplapackint const &ldb, std::complex<double> const &beta, std::complex<double> *c, mplapackint const &ldc);
void Chemv(const char *uplo, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *x, mplapackint const &incx, std::complex<double> const &beta, std::complex<double> *y, mplapackint const &incy);
void Cher(const char *uplo, mplapackint const &n, double const &alpha, std::complex<double> *x, mplapackint const &incx, std::complex<double> *a, mplapackint const &lda);
void Cher2(const char *uplo, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *x, mplapackint const &incx, std::complex<double> *y, mplapackint const &incy, std::complex<double> *a, mplapackint const &lda);
void Cher2k(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *b, mplapackint const &ldb, double const &beta, std::complex<double> *c, mplapackint const &ldc);
void Cherk(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, double const &alpha, std::complex<double> *a, mplapackint const &lda, double const &beta, std::complex<double> *c, mplapackint const &ldc);
void Chpmv(const char *uplo, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *ap, std::complex<double> *x, mplapackint const &incx, std::complex<double> const &beta, std::complex<double> *y, mplapackint const &incy);
void Chpr(const char *uplo, mplapackint const &n, double const &alpha, std::complex<double> *x, mplapackint const &incx, std::complex<double> *ap);
void Chpr2(const char *uplo, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *x, mplapackint const &incx, std::complex<double> *y, mplapackint const &incy, std::complex<double> *ap);
void Crotg(std::complex<double> &ca, std::complex<double> const &cb, double &c, std::complex<double> &s);
void Cscal(mplapackint const &n, std::complex<double> const &za, std::complex<double> *zx, mplapackint const &incx);
void Cswap(mplapackint const &n, std::complex<double> *zx, mplapackint const &incx, std::complex<double> *zy, mplapackint const &incy);
void Csymm(const char *side, const char *uplo, mplapackint const &m, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *b, mplapackint const &ldb, std::complex<double> const &beta, std::complex<double> *c, mplapackint const &ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *b, mplapackint const &ldb, std::complex<double> const &beta, std::complex<double> *c, mplapackint const &ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> const &beta, std::complex<double> *c, mplapackint const &ldc);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mplapackint const &k, std::complex<double> *a, mplapackint const &lda, std::complex<double> *x, mplapackint const &incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mplapackint const &k, std::complex<double> *a, mplapackint const &lda, std::complex<double> *x, mplapackint const &incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, std::complex<double> *ap, std::complex<double> *x, mplapackint const &incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, std::complex<double> *ap, std::complex<double> *x, mplapackint const &incx);
void Ctrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const &m, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *b, mplapackint const &ldb);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, std::complex<double> *a, mplapackint const &lda, std::complex<double> *x, mplapackint const &incx);
void Ctrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const &m, mplapackint const &n, std::complex<double> const &alpha, std::complex<double> *a, mplapackint const &lda, std::complex<double> *b, mplapackint const &ldb);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, std::complex<double> *a, mplapackint const &lda, std::complex<double> *x, mplapackint const &incx);
mplapackint Mlsame_double(const char *a, const char *b);
void Mxerbla_double(const char *srname, int info);
double RCabs1(std::complex<double> const &z);
double RCasum(mplapackint const &n, std::complex<double> *zx, mplapackint const &incx);
double RCnrm2(mplapackint const &n, std::complex<double> *x, mplapackint const &incx);
double Rasum(mplapackint const &n, double *dx, mplapackint const &incx);
void Raxpy(mplapackint const &n, double const &da, double *dx, mplapackint const &incx, double *dy, mplapackint const &incy);
void Rcopy(mplapackint const &n, double *dx, mplapackint const &incx, double *dy, mplapackint const &incy);
double Rdot(mplapackint const &n, double *dx, mplapackint const &incx, double *dy, mplapackint const &incy);
void Rgbmv(const char *trans, mplapackint const &m, mplapackint const &n, mplapackint const &kl, mplapackint const &ku, double const &alpha, double *a, mplapackint const &lda, double *x, mplapackint const &incx, double const &beta, double *y, mplapackint const &incy);
void Rgemm(const char *transa, const char *transb, mplapackint const &m, mplapackint const &n, mplapackint const &k, double const &alpha, double *a, mplapackint const &lda, double *b, mplapackint const &ldb, double const &beta, double *c, mplapackint const &ldc);
void Rgemv(const char *trans, mplapackint const &m, mplapackint const &n, double const &alpha, double *a, mplapackint const &lda, double *x, mplapackint const &incx, double const &beta, double *y, mplapackint const &incy);
void Rger(mplapackint const &m, mplapackint const &n, double const &alpha, double *x, mplapackint const &incx, double *y, mplapackint const &incy, double *a, mplapackint const &lda);
double Rnrm2(mplapackint const &n, double *x, mplapackint const &incx);
void Rrot(mplapackint const &n, double *dx, mplapackint const &incx, double *dy, mplapackint const &incy, double const &c, double const &s);
void Rrotg(double &da, double &db, double &c, double &s);
void Rrotm(mplapackint const &n, double *dx, mplapackint const &incx, double *dy, mplapackint const &incy, double *dparam);
void Rrotmg(double &dd1, double &dd2, double &dx1, double const &dy1, double *dparam);
void Rsbmv(const char *uplo, mplapackint const &n, mplapackint const &k, double const &alpha, double *a, mplapackint const &lda, double *x, mplapackint const &incx, double const &beta, double *y, mplapackint const &incy);
void Rscal(mplapackint const &n, double const &da, double *dx, mplapackint const &incx);
double Rsdot(mplapackint const &n, arr_cref<float> sx, mplapackint const &incx, arr_cref<float> sy, mplapackint const &incy);
void Rspmv(const char *uplo, mplapackint const &n, double const &alpha, double *ap, double *x, mplapackint const &incx, double const &beta, double *y, mplapackint const &incy);
void Rspr(const char *uplo, mplapackint const &n, double const &alpha, double *x, mplapackint const &incx, double *ap);
void Rspr2(const char *uplo, mplapackint const &n, double const &alpha, double *x, mplapackint const &incx, double *y, mplapackint const &incy, double *ap);
void Rswap(mplapackint const &n, double *dx, mplapackint const &incx, double *dy, mplapackint const &incy);
void Rsymm(const char *side, const char *uplo, mplapackint const &m, mplapackint const &n, double const &alpha, double *a, mplapackint const &lda, double *b, mplapackint const &ldb, double const &beta, double *c, mplapackint const &ldc);
void Rsymv(const char *uplo, mplapackint const &n, double const &alpha, double *a, mplapackint const &lda, double *x, mplapackint const &incx, double const &beta, double *y, mplapackint const &incy);
void Rsyr(const char *uplo, mplapackint const &n, double const &alpha, double *x, mplapackint const &incx, double *a, mplapackint const &lda);
void Rsyr2(const char *uplo, mplapackint const &n, double const &alpha, double *x, mplapackint const &incx, double *y, mplapackint const &incy, double *a, mplapackint const &lda);
void Rsyr2k(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, double const &alpha, double *a, mplapackint const &lda, double *b, mplapackint const &ldb, double const &beta, double *c, mplapackint const &ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint const &n, mplapackint const &k, double const &alpha, double *a, mplapackint const &lda, double const &beta, double *c, mplapackint const &ldc);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mplapackint const &k, double *a, mplapackint const &lda, double *x, mplapackint const &incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, mplapackint const &k, double *a, mplapackint const &lda, double *x, mplapackint const &incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, double *ap, double *x, mplapackint const &incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, double *ap, double *x, mplapackint const &incx);
void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const &m, mplapackint const &n, double const &alpha, double *a, mplapackint const &lda, double *b, mplapackint const &ldb);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, double *a, mplapackint const &lda, double *x, mplapackint const &incx);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const &m, mplapackint const &n, double const &alpha, double *a, mplapackint const &lda, double *b, mplapackint const &ldb);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const &n, double *a, mplapackint const &lda, double *x, mplapackint const &incx);
mplapackint iCamax(mplapackint const &n, std::complex<double> *zx, mplapackint const &incx);
mplapackint iRamax(mplapackint const &n, double *dx, mplapackint const &incx);
#endif
