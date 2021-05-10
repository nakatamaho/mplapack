/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_dd.h,v 1.31 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPLAPACK_MATGEN_DD_H_
#define _MPLAPACK_MATGEN_DD_H_

#include "mplapack_config.h"
#include "qd/dd_real.h"
#include "dd_complex.h"

dd_complex Clarnd(mplapackint const idist, mplapackint *iseed);
dd_complex Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, dd_complex *d, mplapackint const igrade, dd_complex *dl, dd_complex *dr, mplapackint const ipvtng, mplapackint *iwork, dd_real const sparse);
dd_complex Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, dd_complex *d, mplapackint const igrade, dd_complex *dl, dd_complex *dr, mplapackint const ipvtng, mplapackint *iwork, dd_real const sparse);
dd_real Rlamch_dd(const char *cmach);
dd_real Rlaran(mplapackint *iseed);
dd_real Rlarnd(mplapackint const idist, mplapackint *iseed);
dd_real Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, dd_real *d, mplapackint const igrade, dd_real *dl, dd_real *dr, mplapackint const ipvtng, mplapackint *iwork, dd_real const sparse);
dd_real Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, dd_real *d, mplapackint const igrade, dd_real *dl, dd_real *dr, mplapackint const ipvtng, mplapackint *iwork, dd_real const sparse);
mplapackint iMlaenv_dd(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
mplapackint iMlaenv_dd2stage(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, dd_real *d, dd_complex *a, mplapackint const lda, mplapackint *iseed, dd_complex *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, dd_real *d, dd_complex *a, mplapackint const lda, mplapackint *iseed, dd_complex *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, dd_real *d, dd_complex *a, mplapackint const lda, mplapackint *iseed, dd_complex *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, dd_complex *a, mplapackint const lda, dd_complex *x, mplapackint const ldx, dd_complex *b, mplapackint const ldb, dd_real *work, mplapackint &info, const char *path);
void Clakf2(mplapackint const m, mplapackint const n, dd_complex *a, mplapackint const lda, dd_complex *b, dd_complex *d, dd_complex *e, dd_complex *z, mplapackint const ldz);
void Clarge(mplapackint const n, dd_complex *a, mplapackint const lda, mplapackint *iseed, dd_complex *work, mplapackint &info);
void Claror(const char *side, const char *init, mplapackint const m, mplapackint const n, dd_complex *a, mplapackint const lda, mplapackint *iseed, dd_complex *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, dd_complex const c, dd_complex const s, dd_complex *a, mplapackint const lda, dd_complex &xleft, dd_complex &xright);
void Clatm1(mplapackint const mode, dd_real const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, dd_complex *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, dd_complex *a, mplapackint const lda, dd_complex *b, mplapackint const ldb, dd_complex *c, mplapackint const ldc, dd_complex *d, mplapackint const ldd, dd_complex *e, mplapackint const lde, dd_complex *f, mplapackint const ldf, dd_complex *r, mplapackint const ldr, dd_complex *l, mplapackint const ldl, dd_real const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, dd_complex *a, mplapackint const lda, dd_complex *b, dd_complex *x, mplapackint const ldx, dd_complex *y, mplapackint const ldy, dd_complex const alpha, dd_complex const beta, dd_complex const wx, dd_complex const wy, dd_real *s, dd_real *dif);
void Clatme(mplapackint const n, const char *dist, mplapackint *iseed, dd_complex *d, mplapackint const mode, dd_real const cond, dd_complex const dmax, const char *rsign, const char *upper, const char *sim, dd_real *ds, mplapackint const modes, dd_real const conds, mplapackint const kl, mplapackint const ku, dd_real const anorm, dd_complex *a, mplapackint const lda, dd_complex *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, dd_complex *d, mplapackint const mode, dd_real const cond, dd_complex const dmax, const char *rsign, const char *grade, dd_complex *dl, mplapackint const model, dd_real const condl, dd_complex *dr, mplapackint const moder, dd_real const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, dd_real const sparse, dd_real const anorm, const char *pack, dd_complex *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, dd_real *d, mplapackint const mode, dd_real const cond, dd_real const dmax, mplapackint const kl, mplapackint const ku, const char *pack, dd_complex *a, mplapackint const lda, dd_complex *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, dd_real *d, mplapackint const mode, dd_real const cond, dd_real const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, dd_complex *a, mplapackint const lda, dd_complex *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, dd_real *d, dd_real *a, mplapackint const lda, mplapackint *iseed, dd_real *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, dd_real *d, dd_real *a, mplapackint const lda, mplapackint *iseed, dd_real *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, dd_real *a, mplapackint const lda, dd_real *x, mplapackint const ldx, dd_real *b, mplapackint const ldb, dd_real *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, dd_real *a, mplapackint const lda, dd_real *b, dd_real *d, dd_real *e, dd_real *z, mplapackint const ldz);
void Rlarge(mplapackint const n, dd_real *a, mplapackint const lda, mplapackint *iseed, dd_real *work, mplapackint &info);
void Rlaror(const char *side, const char *init, mplapackint const m, mplapackint const n, dd_real *a, mplapackint const lda, mplapackint *iseed, dd_real *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, dd_real const c, dd_real const s, dd_real *a, mplapackint const lda, dd_real &xleft, dd_real &xright);
void Rlatm1(mplapackint const mode, dd_real const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, dd_real *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, dd_real *a, mplapackint const lda, dd_real *b, mplapackint const ldb, dd_real *c, mplapackint const ldc, dd_real *d, mplapackint const ldd, dd_real *e, mplapackint const lde, dd_real *f, mplapackint const ldf, dd_real *r, mplapackint const ldr, dd_real *l, mplapackint const ldl, dd_real const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, dd_real *a, mplapackint const lda, dd_real *b, dd_real *x, mplapackint const ldx, dd_real *y, mplapackint const ldy, dd_real const alpha, dd_real const beta, dd_real const wx, dd_real const wy, dd_real *s, dd_real *dif);
void Rlatm7(mplapackint const mode, dd_real const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, dd_real *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatme(mplapackint const n, const char *dist, mplapackint *iseed, dd_real *d, mplapackint const mode, dd_real const cond, dd_real const dmax, const char *ei, const char *rsign, const char *upper, const char *sim, dd_real *ds, mplapackint const modes, dd_real const conds, mplapackint const kl, mplapackint const ku, dd_real const anorm, dd_real *a, mplapackint const lda, dd_real *work, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, dd_real *d, mplapackint const mode, dd_real const cond, dd_real const dmax, const char *rsign, const char *grade, dd_real *dl, mplapackint const model, dd_real const condl, dd_real *dr, mplapackint const moder, dd_real const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, dd_real const sparse, dd_real const anorm, const char *pack, dd_real *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, dd_real *d, mplapackint const mode, dd_real const cond, dd_real const dmax, mplapackint const kl, mplapackint const ku, const char *pack, dd_real *a, mplapackint const lda, dd_real *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, dd_real *d, mplapackint const mode, dd_real const cond, dd_real const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, dd_real *a, mplapackint const lda, dd_real *work, mplapackint &info);
#endif
