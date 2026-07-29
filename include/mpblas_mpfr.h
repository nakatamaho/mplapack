/*
 * Copyright (c) 2008-2025
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

#ifndef _MPBLAS_MPFR_H_
#define _MPBLAS_MPFR_H_

#define MPLAPACK_MPFR_DEFAULT_PRECISION 512

#include <mplapack_gmpfrxx_mkII_config.h>
#include <mpfrxx_mkII.h>
#include <mpcxx_mkII.h>
using namespace mpfrxx;
#include "mplapack_config.h"
#include "mplapack_utils_mpfr.h"

bool Mlsame_mpfr(const char *a, const char *b);
mpc_class Cdotc(mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy);
mpc_class Cdotu(mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy);
mplapackint iCamax(mplapackint const n, mpc_class *zx, mplapackint const incx);
mplapackint iRamax(mplapackint const n, mpfr_class *dx, mplapackint const incx);
mpfr_class RCabs1(mpc_class const z);
mpfr_class RCasum(mplapackint const n, mpc_class *zx, mplapackint const incx);
mpfr_class RCnrm2(mplapackint const n, mpc_class *x, mplapackint const incx);
mpfr_class Rasum(mplapackint const n, mpfr_class *dx, mplapackint const incx);
mpfr_class Rdot(mplapackint const n, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy);
mpfr_class Rnrm2(mplapackint const n, mpfr_class *x, mplapackint const incx);
void CRrot(mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy, mpfr_class const c, mpfr_class const s);
void CRscal(mplapackint const n, mpfr_class const da, mpc_class *zx, mplapackint const incx);
void Caxpy(mplapackint const n, mpc_class const za, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy);
void Ccopy(mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy);
void Cgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy);
void Cgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc);
void Cgemmtr(const char *uplo, const char *transa, const char *transb, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc);
void Cgemv(const char *trans, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy);
void Cgerc(mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *a, mplapackint const lda);
void Cgeru(mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *a, mplapackint const lda);
void Chbmv(const char *uplo, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy);
void Chemm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc);
void Chemv(const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy);
void Cher(const char *uplo, mplapackint const n, mpfr_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *a, mplapackint const lda);
void Cher2(const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *a, mplapackint const lda);
void Cher2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpfr_class const beta, mpc_class *c, mplapackint const ldc);
void Cherk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpc_class *a, mplapackint const lda, mpfr_class const beta, mpc_class *c, mplapackint const ldc);
void Chpmv(const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *ap, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy);
void Chpr(const char *uplo, mplapackint const n, mpfr_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *ap);
void Chpr2(const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *ap);
void Crotg(mpc_class &a, mpc_class const b, mpfr_class &c, mpc_class &s);
void Cscal(mplapackint const n, mpc_class const za, mpc_class *zx, mplapackint const incx);
void Cswap(mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy);
void Csymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc);
void Csyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc);
void Csyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class const beta, mpc_class *c, mplapackint const ldc);
void Ctbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx);
void Ctbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx);
void Ctpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpc_class *ap, mpc_class *x, mplapackint const incx);
void Ctpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpc_class *ap, mpc_class *x, mplapackint const incx);
void Ctrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb);
void Ctrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx);
void Ctrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb);
void Ctrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx);
void Mxerbla_mpfr(const char *srname, int info);
void Raxpy(mplapackint const n, mpfr_class const da, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy);
void Rcopy(mplapackint const n, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy);
void Rgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy);
void Rgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *c, mplapackint const ldc);
void Rgemmtr(const char *uplo, const char *transa, const char *transb, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *c, mplapackint const ldc);
void Rgemv(const char *trans, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy);
void Rger(mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *a, mplapackint const lda);
void Rrot(mplapackint const n, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy, mpfr_class const c, mpfr_class const s);
void Rrotg(mpfr_class &a, mpfr_class &b, mpfr_class &c, mpfr_class &s);
void Rrotm(mplapackint const n, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy, mpfr_class *dparam);
void Rrotmg(mpfr_class &dd1, mpfr_class &dd2, mpfr_class &dx1, mpfr_class const dy1, mpfr_class *dparam);
void Rsbmv(const char *uplo, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy);
void Rscal(mplapackint const n, mpfr_class const da, mpfr_class *dx, mplapackint const incx);
void Rspmv(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *ap, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy);
void Rspr(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *ap);
void Rspr2(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *ap);
void Rswap(mplapackint const n, mpfr_class *dx, mplapackint const incx, mpfr_class *dy, mplapackint const incy);
void Rsymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *c, mplapackint const ldc);
void Rsymv(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx, mpfr_class const beta, mpfr_class *y, mplapackint const incy);
void Rsyr(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *a, mplapackint const lda);
void Rsyr2(const char *uplo, mplapackint const n, mpfr_class const alpha, mpfr_class *x, mplapackint const incx, mpfr_class *y, mplapackint const incy, mpfr_class *a, mplapackint const lda);
void Rsyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class const beta, mpfr_class *c, mplapackint const ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class const beta, mpfr_class *c, mplapackint const ldc);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpfr_class *ap, mpfr_class *x, mplapackint const incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpfr_class *ap, mpfr_class *x, mplapackint const incx);
void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpfr_class const alpha, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const incx);
#endif
