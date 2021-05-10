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

#ifndef _MPLAPACK_MATGEN_MPFR_H_
#define _MPLAPACK_MATGEN_MPFR_H_

#include "mplapack_config.h"
#include "mpc_class.h"
#include "gmpxx.h"
#include "mpreal.h"
#include "mpcomplex.h"

using namespace mpfr;

mpcomplex Clarnd(mplapackint const idist, mplapackint *iseed);
mpcomplex Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, mpcomplex *d, mplapackint const igrade, mpcomplex *dl, mpcomplex *dr, mplapackint const ipvtng, mplapackint *iwork, mpreal const sparse);
mpcomplex Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, mpcomplex *d, mplapackint const igrade, mpcomplex *dl, mpcomplex *dr, mplapackint const ipvtng, mplapackint *iwork, mpreal const sparse);
mplapackint iMlaenv_mpfr(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
mplapackint iMlaenv_mpfr2stage(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
mpreal Rlamch_mpfr(const char *cmach);
mpreal Rlaran(mplapackint *iseed);
mpreal Rlarnd(mplapackint const idist, mplapackint *iseed);
mpreal Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, mpreal *d, mplapackint const igrade, mpreal *dl, mpreal *dr, mplapackint const ipvtng, mplapackint *iwork, mpreal const sparse);
mpreal Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, mpreal *d, mplapackint const igrade, mpreal *dl, mpreal *dr, mplapackint const ipvtng, mplapackint *iwork, mpreal const sparse);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpreal *d, mpcomplex *a, mplapackint const lda, mplapackint *iseed, mpcomplex *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, mpreal *d, mpcomplex *a, mplapackint const lda, mplapackint *iseed, mpcomplex *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, mpreal *d, mpcomplex *a, mplapackint const lda, mplapackint *iseed, mpcomplex *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, mpcomplex *a, mplapackint const lda, mpcomplex *x, mplapackint const ldx, mpcomplex *b, mplapackint const ldb, mpreal *work, mplapackint &info, const char *path);
void Clakf2(mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *b, mpcomplex *d, mpcomplex *e, mpcomplex *z, mplapackint const ldz);
void Clarge(mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *iseed, mpcomplex *work, mplapackint &info);
void Claror(const char *side, const char *init, mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mplapackint *iseed, mpcomplex *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, mpcomplex const c, mpcomplex const s, mpcomplex *a, mplapackint const lda, mpcomplex &xleft, mpcomplex &xright);
void Clatm1(mplapackint const mode, mpreal const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, mpcomplex *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *b, mplapackint const ldb, mpcomplex *c, mplapackint const ldc, mpcomplex *d, mplapackint const ldd, mpcomplex *e, mplapackint const lde, mpcomplex *f, mplapackint const ldf, mpcomplex *r, mplapackint const ldr, mpcomplex *l, mplapackint const ldl, mpreal const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, mpcomplex *a, mplapackint const lda, mpcomplex *b, mpcomplex *x, mplapackint const ldx, mpcomplex *y, mplapackint const ldy, mpcomplex const alpha, mpcomplex const beta, mpcomplex const wx, mpcomplex const wy, mpreal *s, mpreal *dif);
void Clatme(mplapackint const n, const char *dist, mplapackint *iseed, mpcomplex *d, mplapackint const mode, mpreal const cond, mpcomplex const dmax, const char *rsign, const char *upper, const char *sim, mpreal *ds, mplapackint const modes, mpreal const conds, mplapackint const kl, mplapackint const ku, mpreal const anorm, mpcomplex *a, mplapackint const lda, mpcomplex *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpcomplex *d, mplapackint const mode, mpreal const cond, mpcomplex const dmax, const char *rsign, const char *grade, mpcomplex *dl, mplapackint const model, mpreal const condl, mpcomplex *dr, mplapackint const moder, mpreal const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, mpreal const sparse, mpreal const anorm, const char *pack, mpcomplex *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpreal *d, mplapackint const mode, mpreal const cond, mpreal const dmax, mplapackint const kl, mplapackint const ku, const char *pack, mpcomplex *a, mplapackint const lda, mpcomplex *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpreal *d, mplapackint const mode, mpreal const cond, mpreal const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, mpcomplex *a, mplapackint const lda, mpcomplex *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpreal *d, mpreal *a, mplapackint const lda, mplapackint *iseed, mpreal *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, mpreal *d, mpreal *a, mplapackint const lda, mplapackint *iseed, mpreal *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, mpreal *a, mplapackint const lda, mpreal *x, mplapackint const ldx, mpreal *b, mplapackint const ldb, mpreal *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *b, mpreal *d, mpreal *e, mpreal *z, mplapackint const ldz);
void Rlarge(mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *iseed, mpreal *work, mplapackint &info);
void Rlaror(const char *side, const char *init, mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mplapackint *iseed, mpreal *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, mpreal const c, mpreal const s, mpreal *a, mplapackint const lda, mpreal &xleft, mpreal &xright);
void Rlatm1(mplapackint const mode, mpreal const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, mpreal *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *b, mplapackint const ldb, mpreal *c, mplapackint const ldc, mpreal *d, mplapackint const ldd, mpreal *e, mplapackint const lde, mpreal *f, mplapackint const ldf, mpreal *r, mplapackint const ldr, mpreal *l, mplapackint const ldl, mpreal const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *b, mpreal *x, mplapackint const ldx, mpreal *y, mplapackint const ldy, mpreal const alpha, mpreal const beta, mpreal const wx, mpreal const wy, mpreal *s, mpreal *dif);
void Rlatm7(mplapackint const mode, mpreal const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, mpreal *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatme(mplapackint const n, const char *dist, mplapackint *iseed, mpreal *d, mplapackint const mode, mpreal const cond, mpreal const dmax, const char *ei, const char *rsign, const char *upper, const char *sim, mpreal *ds, mplapackint const modes, mpreal const conds, mplapackint const kl, mplapackint const ku, mpreal const anorm, mpreal *a, mplapackint const lda, mpreal *work, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpreal *d, mplapackint const mode, mpreal const cond, mpreal const dmax, const char *rsign, const char *grade, mpreal *dl, mplapackint const model, mpreal const condl, mpreal *dr, mplapackint const moder, mpreal const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, mpreal const sparse, mpreal const anorm, const char *pack, mpreal *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpreal *d, mplapackint const mode, mpreal const cond, mpreal const dmax, mplapackint const kl, mplapackint const ku, const char *pack, mpreal *a, mplapackint const lda, mpreal *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpreal *d, mplapackint const mode, mpreal const cond, mpreal const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, mpreal *a, mplapackint const lda, mpreal *work, mplapackint &info);
#endif
