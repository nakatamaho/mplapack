/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_mpfr.h,v 1.5 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPLAPACK_MPFR_H_
#define _MPLAPACK_MPFR_H_

#include "mplapack_config.h"
#include "mpc_class.h"
#include "gmpxx.h"
#include "mpreal.h"
#include "mpcomplex.h"

using namespace mpfr;

bool Risnan(mpreal const din);
bool Rlaisnan(mpreal const din1, mpreal const din2);
mpcomplex Cladiv(mpcomplex const x, mpcomplex const y);
mplapackint Rlaneg(mplapackint const n, mpreal *d, mpreal *lld, mpreal const sigma, mpreal const, mplapackint const r);
mplapackint iCmax1(mplapackint const n, mpcomplex *zx, mplapackint const incx);
mplapackint iMladlc(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda);
mplapackint iMladlr(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda);
mplapackint iMlaenv_mpfr(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
mplapackint iMlazlc(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda);
mplapackint iMlazlr(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda);
mpreal Clangt(const char *norm, mplapackint const n, mpcomplex *dl, mpcomplex *d, mpcomplex *du);
mpreal Clanhe(const char *norm, const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpreal *work);
mpreal Clanht(const char *norm, mplapackint const n, mpreal *d, mpcomplex *e);
mpreal Clansy(const char *norm, const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpreal *work);
mpreal RCsum1(mplapackint const n, mpcomplex *cx, mplapackint const incx);
mpreal Rla_gerpvgrw(mplapackint const n, mplapackint const ncols, mpreal *a, mplapackint const lda, mpreal *af, mplapackint const ldaf);
mpreal Rladiv2(mpreal const &a, mpreal const &b, mpreal const &c, mpreal const &d, mpreal const &r, mpreal const &t);
mpreal Rladiv2(mpreal const &a, mpreal const &b, mpreal const &c, mpreal const &d, mpreal const &r, mpreal const &t);
mpreal Rlamch_mpfr(const char *cmach);
mpreal Rlangt(const char *norm, mplapackint const n, mpreal *dl, mpreal *d, mpreal *du);
mpreal Rlanst(const char *norm, mplapackint const n, mpreal *d, mpreal *e);
mpreal Rlansy(const char *norm, const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *work);
mpreal Rlapy2(mpreal const x, mpreal const y);
mpreal Rlapy3(mpreal const x, mpreal const y, mpreal const z);
mpreal abs1(mpcomplex ff);
mpreal abs1(mpcomplex ff);
mpreal abs1(mpcomplex ff);
mpreal abssq(mpcomplex ff);
mpreal abssq(mpcomplex ff);
void CRrscl(mplapackint const n, mpreal const sa, mpcomplex *sx, mplapackint const incx);
void Cgebak(const char *job, const char *side, mplapackint const n, mplapackint const ilo, mplapackint const ihi, mpreal *scale, mplapackint const m, mpcomplex *v, mplapackint const ldv, mplapackint &info);
void Cgelq(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *t, mplapackint const tsize, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Cgelqt(mplapackint const m, mplapackint const n, mplapackint const mb, mpcomplex *a, mplapackint const lda, mpcomplex *t, mplapackint const ldt, mpcomplex *work, mplapackint &info);
void Cgemlqt(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const mb, mpcomplex *v, mplapackint const ldv, mpcomplex *t, mplapackint const ldt, mpcomplex *c, mplapackint const ldc, mpcomplex *work, mplapackint &info);
void Cgemqrt(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const nb, mpcomplex *v, mplapackint const ldv, mpcomplex *t, mplapackint const ldt, mpcomplex *c, mplapackint const ldc, mpcomplex *work, mplapackint &info);
void Cgeqr(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *t, mplapackint const tsize, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Cgeqrt(mplapackint const m, mplapackint const n, mplapackint const nb, mpcomplex *a, mplapackint const lda, mpcomplex *t, mplapackint const ldt, mpcomplex *work, mplapackint &info);
void Cgesv(mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cgetc2(mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mplapackint *jpiv, mplapackint &info);
void Cgetf2(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info);
void Cgetrf(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info);
void Cgetrf2(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info);
void Cgetri(mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Cgetrs(const char *trans, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cggbak(const char *job, const char *side, mplapackint const n, mplapackint const ilo, mplapackint const ihi, mpreal *lscale, mpreal *rscale, mplapackint const m, mpcomplex *v, mplapackint const ldv, mplapackint &info);
void Cgtts2(mplapackint const itrans, mplapackint const n, mplapackint const nrhs, mpcomplex *dl, mpcomplex *d, mpcomplex *du, mpcomplex *du2, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb);
void Cheev(const char *jobz, const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpreal *w, mpcomplex *work, mplapackint const lwork, mpreal *rwork, mplapackint &info);
void Chegst(mplapackint const itype, const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Chesv(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Chesv_rook(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Chetd2(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpreal *d, mpreal *e, mpcomplex *tau, mplapackint &info);
void Chetrd(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpreal *d, mpreal *e, mpcomplex *tau, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Chetri(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *work, mplapackint &info);
void Chetri2(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Chetri_3(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *e, mplapackint *ipiv, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Chetri_rook(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *work, mplapackint &info);
void Chetrs(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Chetrs2(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mpcomplex *work, mplapackint &info);
void Chetrs_3(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mpcomplex *e, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Chetrs_aa(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Chetrs_rook(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Chfrk(const char *transr, const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpreal const alpha, mpcomplex *a, mplapackint const lda, mpreal const beta, mpcomplex *c);
void Chpsv(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *ap, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cla_wwaddw(mplapackint const n, mpcomplex *x, mpcomplex *y, mpcomplex *w);
void Clacgv(mplapackint const n, mpcomplex *x, mplapackint const incx);
void Clacp2(const char *uplo, mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb);
void Clacpy(const char *uplo, mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb);
void Clacrm(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mpcomplex *c, mplapackint const ldc, mpreal *rwork);
void Clacrt(mplapackint const n, mpcomplex *cx, mplapackint const incx, mpcomplex *cy, mplapackint const incy, mpcomplex const c, mpcomplex const s);
void Claesy(mpcomplex const a, mpcomplex const b, mpcomplex const c, mpcomplex &rt1, mpcomplex &rt2, mpcomplex &evscal, mpcomplex &cs1, mpcomplex &sn1);
void Claev2(mpcomplex const a, mpcomplex const b, mpcomplex const c, mpreal &rt1, mpreal &rt2, mpreal &cs1, mpcomplex &sn1);
void Clagtm(const char *trans, mplapackint const n, mplapackint const nrhs, mpreal const alpha, mpcomplex *dl, mpcomplex *d, mpcomplex *du, mpcomplex *x, mplapackint const ldx, mpreal const beta, mpcomplex *b, mplapackint const ldb);
void Clahef(const char *uplo, mplapackint const n, mplapackint const nb, mplapackint &kb, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *w, mplapackint const ldw, mplapackint &info);
void Clapmr(bool const forwrd, mplapackint const m, mplapackint const n, mpcomplex *x, mplapackint const ldx, mplapackint *k);
void Clapmt(bool const forwrd, mplapackint const m, mplapackint const n, mpcomplex *x, mplapackint const ldx, mplapackint *k);
void Claqgb(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpcomplex *ab, mplapackint const ldab, mpreal *r, mpreal *c, mpreal const rowcnd, mpreal const colcnd, mpreal const amax, char *equed);
void Claqge(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mpreal *r, mpreal *c, mpreal const rowcnd, mpreal const colcnd, mpreal const amax, char *equed);
void Claqhb(const char *uplo, mplapackint const n, mplapackint const kd, mpcomplex *ab, mplapackint const ldab, mpreal *s, mpreal const scond, mpreal const amax, char *equed);
void Claqhe(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpreal *s, mpreal const scond, mpreal const amax, char *equed);
void Claqhp(const char *uplo, mplapackint const n, mpcomplex *ap, mpreal *s, mpreal const scond, mpreal const amax, char *equed);
void Claqsb(const char *uplo, mplapackint const n, mplapackint const kd, mpcomplex *ab, mplapackint const ldab, mpreal *s, mpreal const scond, mpreal const amax, char *equed);
void Claqsp(const char *uplo, mplapackint const n, mpcomplex *ap, mpreal *s, mpreal const scond, mpreal const amax, char *equed);
void Claqsy(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpreal *s, mpreal const scond, mpreal const amax, char *equed);
void Clar2v(mplapackint const n, mpcomplex *x, mpcomplex *y, mpcomplex *z, mplapackint const incx, mpreal *c, mpcomplex *s, mplapackint const incc);
void Clarf(const char *side, mplapackint const m, mplapackint const n, mpcomplex *v, mplapackint const incv, mpcomplex const tau, mpcomplex *c, mplapackint const ldc, mpcomplex *work);
void Clarfb(const char *side, const char *trans, const char *direct, const char *storev, mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *v, mplapackint const ldv, mpcomplex *t, mplapackint const ldt, mpcomplex *c, mplapackint const ldc, mpcomplex *work, mplapackint const ldwork);
void Clarfb_gett(const char *ident, mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *t, mplapackint const ldt, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mpcomplex *work, mplapackint const ldwork);
void Clarfg(mplapackint const n, mpcomplex &alpha, mpcomplex *x, mplapackint const incx, mpcomplex &tau);
void Clarft(const char *direct, const char *storev, mplapackint const n, mplapackint const k, mpcomplex *v, mplapackint const ldv, mpcomplex *tau, mpcomplex *t, mplapackint const ldt);
void Clarfx(const char *side, mplapackint const m, mplapackint const n, mpcomplex *v, mpcomplex const tau, mpcomplex *c, mplapackint const ldc, mpcomplex *work);
void Clarfy(const char *uplo, mplapackint const n, mpcomplex *v, mplapackint const incv, mpcomplex const tau, mpcomplex *c, mplapackint const ldc, mpcomplex *work);
void Clarscl2(mplapackint const m, mplapackint const n, mpreal *d, mpcomplex *x, mplapackint const ldx);
void Clartg(mpcomplex const f, mpcomplex const g, mpreal &cs, mpcomplex &sn, mpcomplex &r);
void Clartv(mplapackint const n, mpcomplex *x, mplapackint const incx, mpcomplex *y, mplapackint const incy, mpreal *c, mpcomplex *s, mplapackint const incc);
void Clarz(const char *side, mplapackint const m, mplapackint const n, mplapackint const l, mpcomplex *v, mplapackint const incv, mpcomplex const tau, mpcomplex *c, mplapackint const ldc, mpcomplex *work);
void Clarzt(const char *direct, const char *storev, mplapackint const n, mplapackint const k, mpcomplex *v, mplapackint const ldv, mpcomplex *tau, mpcomplex *t, mplapackint const ldt);
void Clascl(const char *type, mplapackint const kl, mplapackint const ku, mpreal const cfrom, mpreal const cto, mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint &info);
void Clascl2(mplapackint const m, mplapackint const n, mpreal *d, mpcomplex *x, mplapackint const ldx);
void Claset(const char *uplo, mplapackint const m, mplapackint const n, mpcomplex const alpha, mpcomplex const beta, mpcomplex *a, mplapackint const lda);
void Clasr(const char *side, const char *pivot, const char *direct, mplapackint const m, mplapackint const n, mpreal *c, mpreal *s, mpcomplex *a, mplapackint const lda);
void Classq(mplapackint const n, mpcomplex *x, mplapackint const incx, mpreal &scale, mpreal &sumsq);
void Claswlq(mplapackint const m, mplapackint const n, mplapackint const mb, mplapackint const nb, mpcomplex *a, mplapackint const lda, mpcomplex *t, mplapackint const ldt, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Claswp(mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint const k1, mplapackint const k2, mplapackint *ipiv, mplapackint const incx);
void Clasyf(const char *uplo, mplapackint const n, mplapackint const nb, mplapackint &kb, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *w, mplapackint const ldw, mplapackint &info);
void Clatrd(const char *uplo, mplapackint const n, mplapackint const nb, mpcomplex *a, mplapackint const lda, mpreal *e, mpcomplex *tau, mpcomplex *w, mplapackint const ldw);
void Clatsqr(mplapackint const m, mplapackint const n, mplapackint const mb, mplapackint const nb, mpcomplex *a, mplapackint const lda, mpcomplex *t, mplapackint const ldt, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Claunhr_col_getrfnp(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *d, mplapackint &info);
void Clauum(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint &info);
void Cpbequ(const char *uplo, mplapackint const n, mplapackint const kd, mpcomplex *ab, mplapackint const ldab, mpreal *s, mpreal &scond, mpreal &amax, mplapackint &info);
void Cpbsv(const char *uplo, mplapackint const n, mplapackint const kd, mplapackint const nrhs, mpcomplex *ab, mplapackint const ldab, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cpbtrs(const char *uplo, mplapackint const n, mplapackint const kd, mplapackint const nrhs, mpcomplex *ab, mplapackint const ldab, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cpoequ(mplapackint const n, mpcomplex *a, mplapackint const lda, mpreal *s, mpreal &scond, mpreal &amax, mplapackint &info);
void Cposv(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cpotf2(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint &info);
void Cpotrf(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint &info);
void Cpotrf2(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint &info);
void Cpotri(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint &info);
void Cpotrs(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cppequ(const char *uplo, mplapackint const n, mpcomplex *ap, mpreal *s, mpreal &scond, mpreal &amax, mplapackint &info);
void Cppsv(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *ap, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cpptrs(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *ap, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cptcon(mplapackint const n, mpreal *d, mpcomplex *e, mpreal const anorm, mpreal &rcond, mpreal *rwork, mplapackint &info);
void Cptsv(mplapackint const n, mplapackint const nrhs, mpreal *d, mpcomplex *e, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cpttrf(mplapackint const n, mpreal *d, mpcomplex *e, mplapackint &info);
void Cpttrs(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *d, mpcomplex *e, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Cptts2(mplapackint const iuplo, mplapackint const n, mplapackint const nrhs, mpreal *d, mpcomplex *e, mpcomplex *b, mplapackint const ldb);
void Crot(mplapackint const n, mpcomplex *cx, mplapackint const incx, mpcomplex *cy, mplapackint const incy, mpreal const c, mpcomplex const s);
void Cspmv(const char *uplo, mplapackint const n, mpcomplex const alpha, mpcomplex *ap, mpcomplex *x, mplapackint const incx, mpcomplex const beta, mpcomplex *y, mplapackint const incy);
void Cspr(const char *uplo, mplapackint const n, mpcomplex const alpha, mpcomplex *x, mplapackint const incx, mpcomplex *ap);
void Cspsv(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *ap, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Csteqr(const char *compz, mplapackint const n, mpreal *d, mpreal *e, mpcomplex *z, mplapackint const ldz, mpreal *work, mplapackint &info);
void Csyconv(const char *uplo, const char *way, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *e, mplapackint &info);
void Csyconvf(const char *uplo, const char *way, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *e, mplapackint *ipiv, mplapackint &info);
void Csyconvf_rook(const char *uplo, const char *way, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *e, mplapackint *ipiv, mplapackint &info);
void Csymv(const char *uplo, mplapackint const n, mpcomplex const alpha, mpcomplex *a, mplapackint const lda, mpcomplex *x, mplapackint const incx, mpcomplex const beta, mpcomplex *y, mplapackint const incy);
void Csyr(const char *uplo, mplapackint const n, mpcomplex const alpha, mpcomplex *x, mplapackint const incx, mpcomplex *a, mplapackint const lda);
void Csytri2(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Csytri_3(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *e, mplapackint *ipiv, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Csytrs(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Csytrs2(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mpcomplex *work, mplapackint &info);
void Csytrs_3(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mpcomplex *e, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Csytrs_aa(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Csytrs_rook(const char *uplo, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mplapackint *ipiv, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Ctbtrs(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const kd, mplapackint const nrhs, mpcomplex *ab, mplapackint const ldab, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Ctgexc(bool const wantq, bool const wantz, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mpcomplex *q, mplapackint const ldq, mpcomplex *z, mplapackint const ldz, mplapackint const ifst, mplapackint &ilst, mplapackint &info);
void Ctplqt(mplapackint const m, mplapackint const n, mplapackint const l, mplapackint const mb, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mpcomplex *t, mplapackint const ldt, mpcomplex *work, mplapackint &info);
void Ctpmlqt(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const l, mplapackint const mb, mpcomplex *v, mplapackint const ldv, mpcomplex *t, mplapackint const ldt, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mpcomplex *work, mplapackint &info);
void Ctpmqrt(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const l, mplapackint const nb, mpcomplex *v, mplapackint const ldv, mpcomplex *t, mplapackint const ldt, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mpcomplex *work, mplapackint &info);
void Ctpqrt(mplapackint const m, mplapackint const n, mplapackint const l, mplapackint const nb, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mpcomplex *t, mplapackint const ldt, mpcomplex *work, mplapackint &info);
void Ctprfb(const char *side, const char *trans, const char *direct, const char *storev, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const l, mpcomplex *v, mplapackint const ldv, mpcomplex *t, mplapackint const ldt, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mpcomplex *work, mplapackint const ldwork);
void Ctptrs(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const nrhs, mpcomplex *ap, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Ctpttr(const char *uplo, mplapackint const n, mpcomplex *ap, mpcomplex *a, mplapackint const lda, mplapackint &info);
void Ctrti2(const char *uplo, const char *diag, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint &info);
void Ctrtri(const char *uplo, const char *diag, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint &info);
void Ctrtrs(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mplapackint &info);
void Ctrttp(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *ap, mplapackint &info);
void Cunbdb5(mplapackint const m1, mplapackint const m2, mplapackint const n, mpcomplex *x1, mplapackint const incx1, mpcomplex *x2, mplapackint const incx2, mpcomplex *q1, mplapackint const ldq1, mpcomplex *q2, mplapackint const ldq2, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Cunbdb6(mplapackint const m1, mplapackint const m2, mplapackint const n, mpcomplex *x1, mplapackint const incx1, mpcomplex *x2, mplapackint const incx2, mpcomplex *q1, mplapackint const ldq1, mpcomplex *q2, mplapackint const ldq2, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Cung2l(mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *work, mplapackint &info);
void Cung2r(mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *work, mplapackint &info);
void Cunghr(mplapackint const n, mplapackint const ilo, mplapackint const ihi, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Cungl2(mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *work, mplapackint &info);
void Cungql(mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Cungqr(mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Cungr2(mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *work, mplapackint &info);
void Cungtr(const char *uplo, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Cungtsqr(mplapackint const m, mplapackint const n, mplapackint const mb, mplapackint const nb, mpcomplex *a, mplapackint const lda, mpcomplex *t, mplapackint const ldt, mpcomplex *work, mplapackint const lwork, mplapackint &info);
void Cunhr_col(mplapackint const m, mplapackint const n, mplapackint const nb, mpcomplex *a, mplapackint const lda, mpcomplex *t, mplapackint const ldt, mpcomplex *d, mplapackint &info);
void Cunm2l(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *c, mplapackint const ldc, mpcomplex *work, mplapackint &info);
void Cunm2r(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *c, mplapackint const ldc, mpcomplex *work, mplapackint &info);
void Cunml2(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *c, mplapackint const ldc, mpcomplex *work, mplapackint &info);
void Cunmr2(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *c, mplapackint const ldc, mpcomplex *work, mplapackint &info);
void Cunmr3(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const l, mpcomplex *a, mplapackint const lda, mpcomplex *tau, mpcomplex *c, mplapackint const ldc, mpcomplex *work, mplapackint &info);
void Rbdsqr(const char *uplo, mplapackint const n, mplapackint const ncvt, mplapackint const nru, mplapackint const ncc, mpreal *d, mpreal *e, mpreal *vt, mplapackint const ldvt, mpreal *u, mplapackint const ldu, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint &info);
void Rcombssq(mpreal *v1, mpreal *v2);
void Rgebak(const char *job, const char *side, mplapackint const n, mplapackint const ilo, mplapackint const ihi, mpreal *scale, mplapackint const m, mpreal *v, mplapackint const ldv, mplapackint &info);
void Rgebd2(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *d, mpreal *e, mpreal *tauq, mpreal *taup, mpreal *work, mplapackint &info);
void Rgebrd(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *d, mpreal *e, mpreal *tauq, mpreal *taup, mpreal *work, mplapackint const lwork, int &info);
void Rgelq(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *t, mplapackint const tsize, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rgelq2(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint &info);
void Rgelqf(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rgelqt(mplapackint const m, mplapackint const n, mplapackint const mb, mpreal *a, mplapackint const lda, mpreal *t, mplapackint const ldt, mpreal *work, mplapackint &info);
void Rgemlqt(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const mb, mpreal *v, mplapackint const ldv, mpreal *t, mplapackint const ldt, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint &info);
void Rgemqrt(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const nb, mpreal *v, mplapackint const ldv, mpreal *t, mplapackint const ldt, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint &info);
void Rgeqr(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *t, mplapackint const tsize, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rgeqr2(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint &info);
void Rgeqrf(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rgeqrt(mplapackint const m, mplapackint const n, mplapackint const nb, mpreal *a, mplapackint const lda, mpreal *t, mplapackint const ldt, mpreal *work, mplapackint &info);
void Rgesv(mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rgesvd(const char *jobu, const char *jobvt, mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *s, mpreal *u, mplapackint const ldu, mpreal *vt, mplapackint const ldvt, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rgetc2(mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mplapackint *jpiv, mplapackint &info);
void Rgetf2(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info);
void Rgetrf(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info);
void Rgetrf2(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info);
void Rgetri(mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rgetrs(const char *trans, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rggbak(const char *job, const char *side, mplapackint const n, mplapackint const ilo, mplapackint const ihi, mpreal *lscale, mpreal *rscale, mplapackint const m, mpreal *v, mplapackint const ldv, mplapackint &info);
void Rgtsv(mplapackint const n, mplapackint const nrhs, mpreal *dl, mpreal *d, mpreal *du, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rgttrf(mplapackint const n, mpreal *dl, mpreal *d, mpreal *du, mpreal *du2, mplapackint *ipiv, mplapackint &info);
void Rgtts2(mplapackint const itrans, mplapackint const n, mplapackint const nrhs, mpreal *dl, mpreal *d, mpreal *du, mpreal *du2, mplapackint *ipiv, mpreal *b, mplapackint const ldb);
void Rla_wwaddw(mplapackint const n, mpreal *x, mpreal *y, mpreal *w);
void Rlabad(mpreal &small, mpreal &large);
void Rlabrd(mplapackint const m, mplapackint const n, mplapackint const nb, mpreal *a, mplapackint const lda, mpreal *d, mpreal *e, mpreal *tauq, mpreal *taup, mpreal *x, mplapackint const ldx, mpreal *y, mplapackint const ldy);
void Rlacpy(const char *uplo, mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb);
void Rladiv(mpreal const &a, mpreal const &b, mpreal const &c, mpreal const &d, mpreal &p, mpreal &q);
void Rladiv(mpreal const &a, mpreal const &b, mpreal const &c, mpreal const &d, mpreal &p, mpreal &q);
void Rladiv1(mpreal &a, mpreal const &b, mpreal const &c, mpreal const &d, mpreal &p, mpreal &q);
void Rladiv1(mpreal &a, mpreal const &b, mpreal const &c, mpreal const &d, mpreal &p, mpreal &q);
void Rlae2(mpreal const a, mpreal const b, mpreal const c, mpreal &rt1, mpreal &rt2);
void Rlaed5(mplapackint const i, mpreal *d, mpreal *z, mpreal *delta, mpreal const rho, mpreal &dlam);
void Rlaev2(mpreal const a, mpreal const b, mpreal const c, mpreal &rt1, mpreal &rt2, mpreal &cs1, mpreal &sn1);
void Rlagtf(mplapackint const n, mpreal *a, mpreal const lambda, mpreal *b, mpreal *c, mpreal const tol, mpreal *d, mplapackint *in, mplapackint &info);
void Rlagtm(const char *trans, mplapackint const n, mplapackint const nrhs, mpreal const alpha, mpreal *dl, mpreal *d, mpreal *du, mpreal *x, mplapackint const ldx, mpreal const beta, mpreal *b, mplapackint const ldb);
void Rlaic1(mplapackint const job, mplapackint const j, mpreal *x, mpreal const sest, mpreal *w, mpreal const gamma, mpreal &sestpr, mpreal &s, mpreal &c);
void Rlamrg(mplapackint const n1, mplapackint const n2, mpreal *a, mplapackint const dtrd1, mplapackint const dtrd2, mplapackint *index);
void Rlanv2(mpreal &a, mpreal &b, mpreal &c, mpreal &d, mpreal &rt1r, mpreal &rt1i, mpreal &rt2r, mpreal &rt2i, mpreal &cs, mpreal &sn);
void Rlaorhr_col_getrfnp(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *d, mplapackint &info);
void Rlaorhr_col_getrfnp2(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *d, mplapackint &info);
void Rlapmr(bool const forwrd, mplapackint const m, mplapackint const n, mpreal *x, mplapackint const ldx, mplapackint *k);
void Rlapmt(bool const forwrd, mplapackint const m, mplapackint const n, mpreal *x, mplapackint const ldx, mplapackint *k);
void Rlaqgb(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpreal *ab, mplapackint const ldab, mpreal *r, mpreal *c, mpreal const rowcnd, mpreal const colcnd, mpreal const amax, char *equed);
void Rlaqge(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *r, mpreal *c, mpreal const rowcnd, mpreal const colcnd, mpreal const amax, char *equed);
void Rlaqr1(mplapackint const n, mpreal *h, mplapackint const ldh, mpreal const sr1, mpreal const si1, mpreal const sr2, mpreal const si2, mpreal *v);
void Rlaqsb(const char *uplo, mplapackint const n, mplapackint const kd, mpreal *ab, mplapackint const ldab, mpreal *s, mpreal const scond, mpreal const amax, char *equed);
void Rlaqsp(const char *uplo, mplapackint const n, mpreal *ap, mpreal *s, mpreal const scond, mpreal const amax, char *equed);
void Rlaqsy(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *s, mpreal const scond, mpreal const amax, char *equed);
void Rlar1v(mplapackint const n, mplapackint const b1, mplapackint const bn, mpreal const lambda, mpreal *d, mpreal *l, mpreal *ld, mpreal *lld, mpreal const pivmin, mpreal const gaptol, mpreal *z, bool const wantnc, mplapackint &negcnt, mpreal &ztz, mpreal &mingma, mplapackint &r, mplapackint *isuppz, mpreal &nrminv, mpreal &resid, mpreal &rqcorr, mpreal *work);
void Rlar2v(mplapackint const n, mpreal *x, mpreal *y, mpreal *z, mplapackint const incx, mpreal *c, mpreal *s, mplapackint const incc);
void Rlarf(const char *side, mplapackint const m, mplapackint const n, mpreal *v, mplapackint const incv, mpreal const tau, mpreal *c, mplapackint const ldc, mpreal *work);
void Rlarfb(const char *side, const char *trans, const char *direct, const char *storev, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *v, mplapackint const ldv, mpreal *t, mplapackint const ldt, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint const ldwork);
void Rlarfb_gett(const char *ident, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *t, mplapackint const ldt, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint const ldwork);
void Rlarfg(mplapackint const n, mpreal &alpha, mpreal *x, mplapackint const incx, mpreal &tau);
void Rlarfgp(mplapackint const n, mpreal &alpha, mpreal *x, mplapackint const incx, mpreal &tau);
void Rlarft(const char *direct, const char *storev, mplapackint const n, mplapackint const k, mpreal *v, mplapackint const ldv, mpreal *tau, mpreal *t, mplapackint const ldt);
void Rlarfx(const char *side, mplapackint const m, mplapackint const n, mpreal *v, mpreal const tau, mpreal *c, mplapackint const ldc, mpreal *work);
void Rlarfy(const char *uplo, mplapackint const n, mpreal *v, mplapackint const incv, mpreal const tau, mpreal *c, mplapackint const ldc, mpreal *work);
void Rlargv(mplapackint const n, mpreal *x, mplapackint const incx, mpreal *y, mplapackint const incy, mpreal *c, mplapackint const incc);
void Rlarra(mplapackint const n, mpreal *d, mpreal *e, mpreal *e2, mpreal const spltol, mpreal const tnrm, mplapackint &nsplit, mplapackint *isplit, mplapackint &info);
void Rlarrc(const char *jobt, mplapackint const n, mpreal const vl, mpreal const vu, mpreal *d, mpreal *e, mpreal const, mplapackint &eigcnt, mplapackint &lcnt, mplapackint &rcnt, mplapackint &info);
void Rlarrj(mplapackint const n, mpreal *d, mpreal *e2, mplapackint const ifirst, mplapackint const ilast, mpreal const rtol, mplapackint const offset, mpreal *w, mpreal *werr, mpreal *work, mplapackint *iwork, mpreal const pivmin, mpreal const spdiam, mplapackint &info);
void Rlarrr(mplapackint const n, mpreal *d, mpreal *e, mplapackint &info);
void Rlarscl2(mplapackint const m, mplapackint const n, mpreal *d, mpreal *x, mplapackint const ldx);
void Rlartg(mpreal const f, mpreal const g, mpreal &cs, mpreal &sn, mpreal &r);
void Rlartgp(mpreal const f, mpreal const g, mpreal &cs, mpreal &sn, mpreal &r);
void Rlartgs(mpreal const x, mpreal const y, mpreal const sigma, mpreal cs, mpreal sn);
void Rlartv(mplapackint const n, mpreal *x, mplapackint const incx, mpreal *y, mplapackint const incy, mpreal *c, mpreal *s, mplapackint const incc);
void Rlarz(const char *side, mplapackint const m, mplapackint const n, mplapackint const l, mpreal *v, mplapackint const incv, mpreal const tau, mpreal *c, mplapackint const ldc, mpreal *work);
void Rlarzt(const char *direct, const char *storev, mplapackint const n, mplapackint const k, mpreal *v, mplapackint const ldv, mpreal *tau, mpreal *t, mplapackint const ldt);
void Rlas2(mpreal const f, mpreal const g, mpreal const h, mpreal &ssmin, mpreal &ssmax);
void Rlascl(const char *type, mplapackint const kl, mplapackint const ku, mpreal const cfrom, mpreal const cto, mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint &info);
void Rlascl2(mplapackint const m, mplapackint const n, mpreal *d, mpreal *x, mplapackint const ldx);
void Rlasd5(mplapackint const i, mpreal *d, mpreal *z, mpreal *delta, mpreal const rho, mpreal &dsigma, mpreal *work);
void Rlaset(const char *uplo, mplapackint const m, mplapackint const n, mpreal const alpha, mpreal const beta, mpreal *a, mplapackint const lda);
void Rlasq1(mplapackint const n, mpreal *d, mpreal *e, mpreal *work, mplapackint &info);
void Rlasq2(mplapackint const n, mpreal *z, mplapackint &info);
void Rlasq3(mplapackint const i0, mplapackint &n0, mpreal *z, mplapackint &pp, mpreal &dmin, mpreal &sigma, mpreal &desig, mpreal &qmax, mplapackint &nfail, mplapackint &iter, mplapackint &ndiv, bool const ieee, mplapackint &ttype, mpreal const dmin1, mpreal &dmin2, mpreal const dn, mpreal const dn1, mpreal const dn2, mpreal const g, mpreal &tau);
void Rlasq4(mplapackint const i0, mplapackint const n0, mpreal *z, mplapackint const pp, mplapackint const n0in, mpreal const dmin, mpreal const dmin1, mpreal const dmin2, mpreal const dn, mpreal const dn1, mpreal const dn2, mpreal &tau, mplapackint &ttype, mpreal &g);
void Rlasq6(mplapackint const i0, mplapackint const n0, mpreal *z, mplapackint const pp, mpreal &dmin, mpreal &dmin1, mpreal &dmin2, mpreal &dn, mpreal &dnm1, mpreal &dnm2);
void Rlasr(const char *side, const char *pivot, const char *direct, mplapackint const m, mplapackint const n, mpreal *c, mpreal *s, mpreal *a, mplapackint const lda);
void Rlasrt(const char *id, mplapackint const n, mpreal *d, mplapackint &info);
void Rlassq(mplapackint const n, mpreal *x, mplapackint const incx, mpreal &scale, mpreal &sumsq);
void Rlasv2(mpreal const f, mpreal const g, mpreal const h, mpreal &ssmin, mpreal &ssmax, mpreal &snr, mpreal &csr, mpreal &snl, mpreal &csl);
void Rlaswlq(mplapackint const m, mplapackint const n, mplapackint const mb, mplapackint const nb, mpreal *a, mplapackint const lda, mpreal *t, mplapackint const ldt, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rlaswp(mplapackint const n, mpreal *a, mplapackint const lda, mplapackint const k1, mplapackint const k2, mplapackint *ipiv, mplapackint const incx);
void Rlasyf(const char *uplo, mplapackint const n, mplapackint const nb, mplapackint &kb, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *w, mplapackint const ldw, mplapackint &info);
void Rlasyf_aa(const char *uplo, mplapackint const j1, mplapackint const m, mplapackint const nb, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *h, mplapackint const ldh, mpreal *work);
void Rlasyf_rk(const char *uplo, mplapackint const n, mplapackint const nb, mplapackint &kb, mpreal *a, mplapackint const lda, mpreal *e, mplapackint *ipiv, mpreal *w, mplapackint const ldw, mplapackint &info);
void Rlasyf_rook(const char *uplo, mplapackint const n, mplapackint const nb, mplapackint &kb, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *w, mplapackint const ldw, mplapackint &info);
void Rlatrd(const char *uplo, mplapackint const n, mplapackint const nb, mpreal *a, mplapackint const lda, mpreal *e, mpreal *tau, mpreal *w, mplapackint const ldw);
void Rlatsqr(mplapackint const m, mplapackint const n, mplapackint const mb, mplapackint const nb, mpreal *a, mplapackint const lda, mpreal *t, mplapackint const ldt, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rlauu2(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint &info);
void Rlauum(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint &info);
void Rorbdb5(mplapackint const m1, mplapackint const m2, mplapackint const n, mpreal *x1, mplapackint const incx1, mpreal *x2, mplapackint const incx2, mpreal *q1, mplapackint const ldq1, mpreal *q2, mplapackint const ldq2, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rorbdb6(mplapackint const m1, mplapackint const m2, mplapackint const n, mpreal *x1, mplapackint const incx1, mpreal *x2, mplapackint const incx2, mpreal *q1, mplapackint const ldq1, mpreal *q2, mplapackint const ldq2, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rorg2l(mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint &info);
void Rorg2r(mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint &info);
void Rorgbr(const char *vect, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rorghr(mplapackint const n, mplapackint const ilo, mplapackint const ihi, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rorgl2(mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint &info);
void Rorglq(mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rorgql(mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rorgqr(mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rorgtr(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rorhr_col(mplapackint const m, mplapackint const n, mplapackint const nb, mpreal *a, mplapackint const lda, mpreal *t, mplapackint const ldt, mpreal *d, mplapackint &info);
void Rorm2l(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint &info);
void Rorm2l(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint &info);
void Rorm2r(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint &info);
void Rorm2r(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint &info);
void Rormbr(const char *vect, const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rorml2(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint &info);
void Rormlq(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rormql(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rormqr(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mpreal *a, mplapackint const lda, mpreal *tau, mpreal *c, mplapackint const ldc, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rpbequ(const char *uplo, mplapackint const n, mplapackint const kd, mpreal *ab, mplapackint const ldab, mpreal *s, mpreal &scond, mpreal &amax, mplapackint &info);
void Rpbsv(const char *uplo, mplapackint const n, mplapackint const kd, mplapackint const nrhs, mpreal *ab, mplapackint const ldab, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rpbtrs(const char *uplo, mplapackint const n, mplapackint const kd, mplapackint const nrhs, mpreal *ab, mplapackint const ldab, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rpoequ(mplapackint const n, mpreal *a, mplapackint const lda, mpreal *s, mpreal &scond, mpreal &amax, mplapackint &info);
void Rpoequb(mplapackint const n, mpreal *a, mplapackint const lda, mpreal *s, mpreal &scond, mpreal &amax, mplapackint &info);
void Rposv(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rposvx(const char *fact, const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mpreal *af, mplapackint const ldaf, char *equed, mpreal *s, mpreal *b, mplapackint const ldb, mpreal *x, mplapackint const ldx, mpreal &rcond, mpreal *ferr, mpreal *berr, mpreal *work, mplapackint *iwork, mplapackint &info);
void Rpotf2(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint &info);
void Rpotrf(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint &info);
void Rpotrf2(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint &info);
void Rpotri(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint &info);
void Rpotrs(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rppequ(const char *uplo, mplapackint const n, mpreal *ap, mpreal *s, mpreal &scond, mpreal &amax, mplapackint &info);
void Rppsv(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *ap, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rpptrs(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *ap, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rptcon(mplapackint const n, mpreal *d, mpreal *e, mpreal const anorm, mpreal &rcond, mpreal *work, mplapackint &info);
void Rptsv(mplapackint const n, mplapackint const nrhs, mpreal *d, mpreal *e, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rptsvx(const char *fact, mplapackint const n, mplapackint const nrhs, mpreal *d, mpreal *e, mpreal *df, mpreal *ef, mpreal *b, mplapackint const ldb, mpreal *x, mplapackint const ldx, mpreal &rcond, mpreal *ferr, mpreal *berr, mpreal *work, mplapackint &info);
void Rpttrf(mplapackint const n, mpreal *d, mpreal *e, mplapackint &info);
void Rpttrs(mplapackint const n, mplapackint const nrhs, mpreal *d, mpreal *e, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rptts2(mplapackint const n, mplapackint const nrhs, mpreal *d, mpreal *e, mpreal *b, mplapackint const ldb);
void Rrscl(mplapackint const n, mpreal const sa, mpreal *sx, mplapackint const incx);
void Rsfrk(const char *transr, const char *uplo, const char *trans, mplapackint const n, mplapackint const k, mpreal const alpha, mpreal *a, mplapackint const lda, mpreal const beta, mpreal *c);
void Rspsv(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *ap, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rsteqr(const char *compz, mplapackint const n, mpreal *d, mpreal *e, mpreal *z, mplapackint const ldz, mpreal *work, mplapackint &info);
void Rsterf(mplapackint const n, mpreal *d, mpreal *e, mplapackint &info);
void Rstev(const char *jobz, mplapackint const n, mpreal *d, mpreal *e, mpreal *z, mplapackint const ldz, mpreal *work, mplapackint &info);
void Rstevd(const char *jobz, mplapackint const n, mpreal *d, mpreal *e, mpreal *z, mplapackint const ldz, mpreal *work, mplapackint const lwork, mplapackint *iwork, mplapackint const liwork, mplapackint &info);
void Rsyconv(const char *uplo, const char *way, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *e, mplapackint &info);
void Rsyconvf(const char *uplo, const char *way, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *e, mplapackint *ipiv, mplapackint &info);
void Rsyconvf_rook(const char *uplo, const char *way, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *e, mplapackint *ipiv, mplapackint &info);
void Rsyequb(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *s, mpreal &scond, mpreal &amax, mpreal *work, mplapackint &info);
void Rsyev(const char *jobz, const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *w, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rsyevd(const char *jobz, const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *w, mpreal *work, mplapackint const lwork, mplapackint *iwork, mplapackint const liwork, mplapackint &info);
void Rsygs2(mplapackint const itype, const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rsygst(mplapackint const itype, const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rsysv(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rsysv_aa(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rsysv_aa_2stage(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mpreal *tb, mplapackint const ltb, mplapackint *ipiv, mplapackint *ipiv2, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rsysv_rk(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mpreal *e, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rsysv_rook(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rsysvx(const char *fact, const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mpreal *af, mplapackint const ldaf, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mpreal *x, mplapackint const ldx, mpreal &rcond, mpreal *ferr, mpreal *berr, mpreal *work, mplapackint const lwork, mplapackint *iwork, mplapackint &info);
void Rsytd2(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *d, mpreal *e, mpreal *tau, mplapackint &info);
void Rsytf2(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info);
void Rsytf2_rk(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *e, mplapackint *ipiv, mplapackint &info);
void Rsytf2_rook(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mplapackint &info);
void Rsytrd(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *d, mpreal *e, mpreal *tau, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rsytri(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *work, mplapackint &info);
void Rsytri2(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rsytri_3(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *e, mplapackint *ipiv, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rsytri_rook(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *work, mplapackint &info);
void Rsytrs(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rsytrs2(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint &info);
void Rsytrs_3(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mpreal *e, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rsytrs_aa(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rsytrs_aa_2stage(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mpreal *tb, mplapackint const ltb, mplapackint *ipiv, mplapackint *ipiv2, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rsytrs_rook(const char *uplo, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mplapackint *ipiv, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rtbtrs(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const kd, mplapackint const nrhs, mpreal *ab, mplapackint const ldab, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rtgexc(bool const wantq, bool const wantz, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mpreal *q, mplapackint const ldq, mpreal *z, mplapackint const ldz, mplapackint &ifst, mplapackint &ilst, mpreal *work, mplapackint const lwork, mplapackint &info);
void Rtplqt(mplapackint const m, mplapackint const n, mplapackint const l, mplapackint const mb, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mpreal *t, mplapackint const ldt, mpreal *work, mplapackint &info);
void Rtpmlqt(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const l, mplapackint const mb, mpreal *v, mplapackint const ldv, mpreal *t, mplapackint const ldt, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint &info);
void Rtpmqrt(const char *side, const char *trans, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const l, mplapackint const nb, mpreal *v, mplapackint const ldv, mpreal *t, mplapackint const ldt, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint &info);
void Rtpqrt(mplapackint const m, mplapackint const n, mplapackint const l, mplapackint const nb, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mpreal *t, mplapackint const ldt, mpreal *work, mplapackint &info);
void Rtprfb(const char *side, const char *trans, const char *direct, const char *storev, mplapackint const m, mplapackint const n, mplapackint const k, mplapackint const l, mpreal *v, mplapackint const ldv, mpreal *t, mplapackint const ldt, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint const ldwork);
void Rtptrs(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const nrhs, mpreal *ap, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rtpttr(const char *uplo, mplapackint const n, mpreal *ap, mpreal *a, mplapackint const lda, mplapackint &info);
void Rtrexc(const char *compq, mplapackint const n, mpreal *t, mplapackint const ldt, mpreal *q, mplapackint const ldq, mplapackint &ifst, mplapackint &ilst, mpreal *work, mplapackint &info);
void Rtrti2(const char *uplo, const char *diag, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint &info);
void Rtrtri(const char *uplo, const char *diag, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint &info);
void Rtrtrs(const char *uplo, const char *trans, const char *diag, mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mplapackint &info);
void Rtrttp(const char *uplo, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *ap, mplapackint &info);
#endif
