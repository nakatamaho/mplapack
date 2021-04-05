/*
 * Copyright (c) 2008-2021
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

#ifndef _MPBLAS_GMP_H_
#define _MPBLAS_GMP_H_

#include "mplapack_config.h"
#include "gmpxx.h"
#include "mpc_class.h"
#include "mplapack_utils_gmp.h"
#if !defined __MPLAPACK_ERRNO__
#define _MPLAPACK_EXTERN_ extern
#else
#define _MPLAPACK_EXTERN_
#endif
_MPLAPACK_EXTERN_ int mplapack_errno;

mpc_class Cdotc(mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy);
mpc_class Cdotu(mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy);
mpf_class RCabs1(mpc_class const z);
mpf_class RCasum(mplapackint const n, mpc_class *zx, mplapackint const incx);
mpf_class RCnrm2(mplapackint const n, mpc_class *x, mplapackint const incx);
mpf_class Rasum(mplapackint const n, mpf_class *dx, mplapackint const incx);
mpf_class Rdot(mplapackint const n, mpf_class *dx, mplapackint const incx, mpf_class *dy, mplapackint const incy);
mpf_class Rnrm2(mplapackint const n, mpf_class *x, mplapackint const incx);
mplapackint Mlsame_gmp(const char *a, const char *b);
mplapackint iCamax(mplapackint const n, mpc_class *zx, mplapackint const incx);
mplapackint iRamax(mplapackint const n, mpf_class *dx, mplapackint const incx);
void CRrot(mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy, mpf_class const c, mpf_class const s);
void CRscal(mplapackint const n, mpf_class const da, mpc_class *zx, mplapackint const incx);
void Caxpy(mplapackint const n, mpc_class const za, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy);
void Ccopy(mplapackint const n, mpc_class *zx, mplapackint const incx, mpc_class *zy, mplapackint const incy);
void Cgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy);
void Cgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc);
void Cgemv(const char *trans, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy);
void Cgerc(mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *a, mplapackint const lda);
void Cgeru(mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *a, mplapackint const lda);
void Chbmv(const char *uplo, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy);
void Chemm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class const beta, mpc_class *c, mplapackint const ldc);
void Chemv(const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy);
void Cher(const char *uplo, mplapackint const n, mpf_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *a, mplapackint const lda);
void Cher2(const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *a, mplapackint const lda);
void Cher2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpc_class const alpha, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpf_class const beta, mpc_class *c, mplapackint const ldc);
void Cherk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpf_class const alpha, mpc_class *a, mplapackint const lda, mpf_class const beta, mpc_class *c, mplapackint const ldc);
void Chpmv(const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *ap, mpc_class *x, mplapackint const incx, mpc_class const beta, mpc_class *y, mplapackint const incy);
void Chpr(const char *uplo, mplapackint const n, mpf_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *ap);
void Chpr2(const char *uplo, mplapackint const n, mpc_class const alpha, mpc_class *x, mplapackint const incx, mpc_class *y, mplapackint const incy, mpc_class *ap);
void Crotg(mpc_class &ca, mpc_class const cb, mpf_class &c, mpc_class &s);
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
void Mxerbla_gmp(const char *srname, int info);
void Raxpy(mplapackint const n, mpf_class const da, mpf_class *dx, mplapackint const incx, mpf_class *dy, mplapackint const incy);
void Rcopy(mplapackint const n, mpf_class *dx, mplapackint const incx, mpf_class *dy, mplapackint const incy);
void Rgbmv(const char *trans, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx, mpf_class const beta, mpf_class *y, mplapackint const incy);
void Rgemm(const char *transa, const char *transb, mplapackint const m, mplapackint const n, mplapackint const k, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb, mpf_class const beta, mpf_class *c, mplapackint const ldc);
void Rgemv(const char *trans, mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx, mpf_class const beta, mpf_class *y, mplapackint const incy);
void Rger(mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class *x, mplapackint const incx, mpf_class *y, mplapackint const incy, mpf_class *a, mplapackint const lda);
void Rrot(mplapackint const n, mpf_class *dx, mplapackint const incx, mpf_class *dy, mplapackint const incy, mpf_class const c, mpf_class const s);
void Rrotg(mpf_class &da, mpf_class &db, mpf_class &c, mpf_class &s);
void Rrotm(mplapackint const n, mpf_class *dx, mplapackint const incx, mpf_class *dy, mplapackint const incy, mpf_class *dparam);
void Rrotmg(mpf_class &dd1, mpf_class &dd2, mpf_class &dx1, mpf_class const dy1, mpf_class *dparam);
void Rsbmv(const char *uplo, mplapackint const n, mplapackint const k, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx, mpf_class const beta, mpf_class *y, mplapackint const incy);
void Rscal(mplapackint const n, mpf_class const da, mpf_class *dx, mplapackint const incx);
void Rspmv(const char *uplo, mplapackint const n, mpf_class const alpha, mpf_class *ap, mpf_class *x, mplapackint const incx, mpf_class const beta, mpf_class *y, mplapackint const incy);
void Rspr(const char *uplo, mplapackint const n, mpf_class const alpha, mpf_class *x, mplapackint const incx, mpf_class *ap);
void Rspr2(const char *uplo, mplapackint const n, mpf_class const alpha, mpf_class *x, mplapackint const incx, mpf_class *y, mplapackint const incy, mpf_class *ap);
void Rswap(mplapackint const n, mpf_class *dx, mplapackint const incx, mpf_class *dy, mplapackint const incy);
void Rsymm(const char *side, const char *uplo, mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb, mpf_class const beta, mpf_class *c, mplapackint const ldc);
void Rsymv(const char *uplo, mplapackint const n, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx, mpf_class const beta, mpf_class *y, mplapackint const incy);
void Rsyr(const char *uplo, mplapackint const n, mpf_class const alpha, mpf_class *x, mplapackint const incx, mpf_class *a, mplapackint const lda);
void Rsyr2(const char *uplo, mplapackint const n, mpf_class const alpha, mpf_class *x, mplapackint const incx, mpf_class *y, mplapackint const incy, mpf_class *a, mplapackint const lda);
void Rsyr2k(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb, mpf_class const beta, mpf_class *c, mplapackint const ldc);
void Rsyrk(const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class const beta, mpf_class *c, mplapackint const ldc);
void Rtbmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx);
void Rtbsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const k, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx);
void Rtpmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpf_class *ap, mpf_class *x, mplapackint const incx);
void Rtpsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpf_class *ap, mpf_class *x, mplapackint const incx);
void Rtrmm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb);
void Rtrmv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx);
void Rtrsm(const char *side, const char *uplo, const char *transa, const char *diag, mplapackint const m, mplapackint const n, mpf_class const alpha, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb);
void Rtrsv(const char *uplo, const char *trans, const char *diag, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const incx);
#endif
