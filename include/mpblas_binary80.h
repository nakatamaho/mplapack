/*
 * Copyright (c) 2012-2025
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

#ifndef _MPBLAS_BINARY80_H_
#define _MPBLAS_BINARY80_H_

#include "mplapack_config.h"
#include <complex>
#include "mplapack_utils_binary80.h"

bool Mlsame_binary80(const char *a, const char *b);
mplapack_binary80_t RCabs1(std::complex<mplapack_binary80_t> const z);
mplapack_binary80_t RCasum(mplapackint const n, std::complex<mplapack_binary80_t> *zx, mplapackint const incx);
mplapack_binary80_t RCnrm2(mplapackint const n, std::complex<mplapack_binary80_t> *x, mplapackint const incx);
mplapack_binary80_t Rasum(mplapackint const n, mplapack_binary80_t *dx, mplapackint const incx);
mplapack_binary80_t Rdot(mplapackint const n, mplapack_binary80_t *dx, mplapackint const incx, mplapack_binary80_t *dy, mplapackint const incy);
mplapack_binary80_t Rnrm2(mplapackint const n, mplapack_binary80_t *x, mplapackint const incx);
mplapackint iCamax(mplapackint const n, std::complex<mplapack_binary80_t> *zx, mplapackint const incx);
mplapackint iRamax(mplapackint const n, mplapack_binary80_t *dx, mplapackint const incx);
std::complex<mplapack_binary80_t> Cdotc(mplapackint const n, std::complex<mplapack_binary80_t> *zx, mplapackint const incx, std::complex<mplapack_binary80_t> *zy, mplapackint const incy);
std::complex<mplapack_binary80_t> Cdotu(mplapackint const n, std::complex<mplapack_binary80_t> *zx, mplapackint const incx, std::complex<mplapack_binary80_t> *zy, mplapackint const incy);
void CRrot(mplapackint const n, std::complex<mplapack_binary80_t> *zx, mplapackint const incx, std::complex<mplapack_binary80_t> *zy, mplapackint const incy, mplapack_binary80_t const c, mplapack_binary80_t const s);
void CRscal(mplapackint const n, mplapack_binary80_t const da, std::complex<mplapack_binary80_t> *zx, mplapackint const incx);
void Caxpy(mplapackint const n, std::complex<mplapack_binary80_t> const za, std::complex<mplapack_binary80_t> *zx, mplapackint const incx, std::complex<mplapack_binary80_t> *zy, mplapackint const incy);
void Ccopy(mplapackint const n, std::complex<mplapack_binary80_t> *zx, mplapackint const incx, std::complex<mplapack_binary80_t> *zy, mplapackint const incy);
void Cgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *y, mplapackint const incy);
void Cgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *b, mplapackint const ldb, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *c, mplapackint const ldc);
void Cgemmtr(const char *uplo, const char *transa, const char *transb, mplapackint const n, mplapackint const k, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *b, mplapackint const ldb, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *c, mplapackint const ldc);
void Cgemv(const char *trans, mplapackint const m, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *y, mplapackint const incy);
void Cgerc(mplapackint const m, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> *y, mplapackint const incy, std::complex<mplapack_binary80_t> *a, mplapackint const lda);
void Cgeru(mplapackint const m, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> *y, mplapackint const incy, std::complex<mplapack_binary80_t> *a, mplapackint const lda);
void Chbmv(const char *uplo, mplapackint const n, mplapackint const k, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *y, mplapackint const incy);
void Chemm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *b, mplapackint const ldb, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *c, mplapackint const ldc);
void Chemv(const char *uplo, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *y, mplapackint const incy);
void Cher(const char *uplo, mplapackint const n, mplapack_binary80_t const alpha, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> *a, mplapackint const lda);
void Cher2(const char *uplo, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> *y, mplapackint const incy, std::complex<mplapack_binary80_t> *a, mplapackint const lda);
void Cher2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *b, mplapackint const ldb, mplapack_binary80_t const beta, std::complex<mplapack_binary80_t> *c, mplapackint const ldc);
void Cherk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mplapack_binary80_t const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, mplapack_binary80_t const beta, std::complex<mplapack_binary80_t> *c, mplapackint const ldc);
void Chpmv(const char *uplo, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *ap, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *y, mplapackint const incy);
void Chpr(const char *uplo, mplapackint const n, mplapack_binary80_t const alpha, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> *ap);
void Chpr2(const char *uplo, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *x, mplapackint const incx, std::complex<mplapack_binary80_t> *y, mplapackint const incy, std::complex<mplapack_binary80_t> *ap);
void Crotg(std::complex<mplapack_binary80_t> &a, std::complex<mplapack_binary80_t> const b, mplapack_binary80_t &c, std::complex<mplapack_binary80_t> &s);
void Cscal(mplapackint const n, std::complex<mplapack_binary80_t> const za, std::complex<mplapack_binary80_t> *zx, mplapackint const incx);
void Cswap(mplapackint const n, std::complex<mplapack_binary80_t> *zx, mplapackint const incx, std::complex<mplapack_binary80_t> *zy, mplapackint const incy);
void Csymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *b, mplapackint const ldb, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *c, mplapackint const ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *b, mplapackint const ldb, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *c, mplapackint const ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> const beta, std::complex<mplapack_binary80_t> *c, mplapackint const ldc);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *x, mplapackint const incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *x, mplapackint const incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, std::complex<mplapack_binary80_t> *ap, std::complex<mplapack_binary80_t> *x, mplapackint const incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, std::complex<mplapack_binary80_t> *ap, std::complex<mplapack_binary80_t> *x, mplapackint const incx);
void Ctrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *b, mplapackint const ldb);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *x, mplapackint const incx);
void Ctrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, std::complex<mplapack_binary80_t> const alpha, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *b, mplapackint const ldb);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, std::complex<mplapack_binary80_t> *a, mplapackint const lda, std::complex<mplapack_binary80_t> *x, mplapackint const incx);
void Mxerbla_binary80(const char *srname, int info);
void Raxpy(mplapackint const n, mplapack_binary80_t const da, mplapack_binary80_t *dx, mplapackint const incx, mplapack_binary80_t *dy, mplapackint const incy);
void Rcopy(mplapackint const n, mplapack_binary80_t *dx, mplapackint const incx, mplapack_binary80_t *dy, mplapackint const incy);
void Rgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *x, mplapackint const incx, mplapack_binary80_t const beta, mplapack_binary80_t *y, mplapackint const incy);
void Rgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *b, mplapackint const ldb, mplapack_binary80_t const beta, mplapack_binary80_t *c, mplapackint const ldc);
void Rgemmtr(const char *uplo, const char *transa, const char *transb, mplapackint const n, mplapackint const k, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *b, mplapackint const ldb, mplapack_binary80_t const beta, mplapack_binary80_t *c, mplapackint const ldc);
void Rgemv(const char *trans, mplapackint const m, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *x, mplapackint const incx, mplapack_binary80_t const beta, mplapack_binary80_t *y, mplapackint const incy);
void Rger(mplapackint const m, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *x, mplapackint const incx, mplapack_binary80_t *y, mplapackint const incy, mplapack_binary80_t *a, mplapackint const lda);
void Rrot(mplapackint const n, mplapack_binary80_t *dx, mplapackint const incx, mplapack_binary80_t *dy, mplapackint const incy, mplapack_binary80_t const c, mplapack_binary80_t const s);
void Rrotg(mplapack_binary80_t &a, mplapack_binary80_t &b, mplapack_binary80_t &c, mplapack_binary80_t &s);
void Rrotm(mplapackint const n, mplapack_binary80_t *dx, mplapackint const incx, mplapack_binary80_t *dy, mplapackint const incy, mplapack_binary80_t *dparam);
void Rrotmg(mplapack_binary80_t &dd1, mplapack_binary80_t &dd2, mplapack_binary80_t &dx1, mplapack_binary80_t const dy1, mplapack_binary80_t *dparam);
void Rsbmv(const char *uplo, mplapackint const n, mplapackint const k, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *x, mplapackint const incx, mplapack_binary80_t const beta, mplapack_binary80_t *y, mplapackint const incy);
void Rscal(mplapackint const n, mplapack_binary80_t const da, mplapack_binary80_t *dx, mplapackint const incx);
void Rspmv(const char *uplo, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *ap, mplapack_binary80_t *x, mplapackint const incx, mplapack_binary80_t const beta, mplapack_binary80_t *y, mplapackint const incy);
void Rspr(const char *uplo, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *x, mplapackint const incx, mplapack_binary80_t *ap);
void Rspr2(const char *uplo, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *x, mplapackint const incx, mplapack_binary80_t *y, mplapackint const incy, mplapack_binary80_t *ap);
void Rswap(mplapackint const n, mplapack_binary80_t *dx, mplapackint const incx, mplapack_binary80_t *dy, mplapackint const incy);
void Rsymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *b, mplapackint const ldb, mplapack_binary80_t const beta, mplapack_binary80_t *c, mplapackint const ldc);
void Rsymv(const char *uplo, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *x, mplapackint const incx, mplapack_binary80_t const beta, mplapack_binary80_t *y, mplapackint const incy);
void Rsyr(const char *uplo, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *x, mplapackint const incx, mplapack_binary80_t *a, mplapackint const lda);
void Rsyr2(const char *uplo, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *x, mplapackint const incx, mplapack_binary80_t *y, mplapackint const incy, mplapack_binary80_t *a, mplapackint const lda);
void Rsyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *b, mplapackint const ldb, mplapack_binary80_t const beta, mplapack_binary80_t *c, mplapackint const ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t const beta, mplapack_binary80_t *c, mplapackint const ldc);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *x, mplapackint const incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *x, mplapackint const incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapack_binary80_t *ap, mplapack_binary80_t *x, mplapackint const incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapack_binary80_t *ap, mplapack_binary80_t *x, mplapackint const incx);
void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *b, mplapackint const ldb);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *x, mplapackint const incx);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mplapack_binary80_t const alpha, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *b, mplapackint const ldb);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapack_binary80_t *a, mplapackint const lda, mplapack_binary80_t *x, mplapackint const incx);
#endif
