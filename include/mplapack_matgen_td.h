/*
 * Copyright (c) 2021-2025
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

#ifndef MPLAPACK_MATGEN_TD_H
#define MPLAPACK_MATGEN_TD_H

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include "mplapack_config.h"
#include "qd/td_real.h"
#include <qd/td_complex.h>

mplapackint iMlaenv_td(mplapackint const ispec, const char *name, const char *opts, mplapackint const n1, mplapackint const n2, mplapackint const n3, mplapackint const n4);
mplapackint iMlaenv_td2stage(mplapackint const ispec, const char *name, const char *opts, mplapackint const n1, mplapackint const n2, mplapackint const n3, mplapackint const n4);
td_complex Clarnd(mplapackint const idist, mplapackint (&iseed)[4]);
td_complex Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], td_complex *d, mplapackint const igrade, td_complex *dl, td_complex *dr, mplapackint const ipvtng, mplapackint *iwork, td_real const sparse);
td_complex Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], td_complex *d, mplapackint const igrade, td_complex *dl, td_complex *dr, mplapackint const ipvtng, mplapackint *iwork, td_real const sparse);
td_real Rlamch_td(const char *cmach);
td_real Rlaran(mplapackint (&iseed)[4]);
td_real Rlarnd(mplapackint const idist, mplapackint (&iseed)[4]);
td_real Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], td_real *d, mplapackint const igrade, td_real *dl, td_real *dr, mplapackint const ipvtng, mplapackint *iwork, td_real const sparse);
td_real Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], td_real *d, mplapackint const igrade, td_real *dl, td_real *dr, mplapackint const ipvtng, mplapackint *iwork, td_real const sparse);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, td_real *d, td_complex *a, mplapackint const lda, mplapackint (&iseed)[4], td_complex *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, td_real *d, td_complex *a, mplapackint const lda, mplapackint (&iseed)[4], td_complex *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, td_real *d, td_complex *a, mplapackint const lda, mplapackint (&iseed)[4], td_complex *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, td_complex *a, mplapackint const lda, td_complex *x, mplapackint const ldx, td_complex *b, mplapackint const ldb, td_real *work, mplapackint &info, fem::str_cref path);
void Clakf2(mplapackint const m, mplapackint const n, td_complex *a, mplapackint const lda, td_complex *b, td_complex *d, td_complex *e, td_complex *z, mplapackint const ldz);
void Clarge(mplapackint const n, td_complex *a, mplapackint const lda, mplapackint (&iseed)[4], td_complex *work, mplapackint &info);
void Claror(fem::str_cref side, fem::str_cref init, mplapackint const m, mplapackint const n, td_complex *a, mplapackint const lda, mplapackint (&iseed)[4], td_complex *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, td_complex const c, td_complex const s, td_complex *a, mplapackint const lda, td_complex &xleft, td_complex &xright);
void Clatm1(mplapackint const mode, td_real const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], td_complex *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, td_complex *a, mplapackint const lda, td_complex *b, mplapackint const ldb, td_complex *c, mplapackint const ldc, td_complex *d, mplapackint const ldd, td_complex *e, mplapackint const lde, td_complex *f, mplapackint const ldf, td_complex *r, mplapackint const ldr, td_complex *l, mplapackint const ldl, td_real const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, td_complex *a, mplapackint const lda, td_complex *b, td_complex *x, mplapackint const ldx, td_complex *y, mplapackint const ldy, td_complex const alpha, td_complex const beta, td_complex const wx, td_complex const wy, td_real *s, td_real *dif);
void Clatme(mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], td_complex *d, mplapackint const mode, td_real const cond, td_complex const dmax, fem::str_cref rsign, fem::str_cref upper, fem::str_cref sim, td_real *ds, mplapackint const modes, td_real const conds, mplapackint const kl, mplapackint const ku, td_real const anorm, td_complex *a, mplapackint const lda, td_complex *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, td_complex *d, mplapackint const mode, td_real const cond, td_complex const dmax, fem::str_cref rsign, fem::str_cref grade, td_complex *dl, mplapackint const model, td_real const condl, td_complex *dr, mplapackint const moder, td_real const condr, fem::str_cref pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, td_real const sparse, td_real const anorm, fem::str_cref pack, td_complex *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, td_real *d, mplapackint const mode, td_real const cond, td_real const dmax, mplapackint const kl, mplapackint const ku, fem::str_cref pack, td_complex *a, mplapackint const lda, td_complex *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, td_real *d, mplapackint const mode, td_real const cond, td_real const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, fem::str_cref pack, td_complex *a, mplapackint const lda, td_complex *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, td_real *d, td_real *a, mplapackint const lda, mplapackint (&iseed)[4], td_real *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, td_real *d, td_real *a, mplapackint const lda, mplapackint (&iseed)[4], td_real *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, td_real *a, mplapackint const lda, td_real *x, mplapackint const ldx, td_real *b, mplapackint const ldb, td_real *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, td_real *a, mplapackint const lda, td_real *b, td_real *d, td_real *e, td_real *z, mplapackint const ldz);
void Rlarge(mplapackint const n, td_real *a, mplapackint const lda, mplapackint (&iseed)[4], td_real *work, mplapackint &info);
void Rlaror(fem::str_cref side, fem::str_cref init, mplapackint const m, mplapackint const n, td_real *a, mplapackint const lda, mplapackint (&iseed)[4], td_real *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, td_real const c, td_real const s, td_real *a, mplapackint const lda, td_real &xleft, td_real &xright);
void Rlatm1(mplapackint const mode, td_real const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], td_real *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, td_real *a, mplapackint const lda, td_real *b, mplapackint const ldb, td_real *c, mplapackint const ldc, td_real *d, mplapackint const ldd, td_real *e, mplapackint const lde, td_real *f, mplapackint const ldf, td_real *r, mplapackint const ldr, td_real *l, mplapackint const ldl, td_real const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, td_real *a, mplapackint const lda, td_real *b, td_real *x, mplapackint const ldx, td_real *y, mplapackint const ldy, td_real const alpha, td_real const beta, td_real const wx, td_real const wy, td_real *s, td_real *dif);
void Rlatm7(mplapackint const mode, td_real const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], td_real *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatme(mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], td_real *d, mplapackint const mode, td_real const cond, td_real const dmax, const char *ei, fem::str_cref rsign, fem::str_cref upper, fem::str_cref sim, td_real *ds, mplapackint const modes, td_real const conds, mplapackint const kl, mplapackint const ku, td_real const anorm, td_real *a, mplapackint const lda, td_real *work, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, td_real *d, mplapackint const mode, td_real const cond, td_real const dmax, fem::str_cref rsign, fem::str_cref grade, td_real *dl, mplapackint const model, td_real const condl, td_real *dr, mplapackint const moder, td_real const condr, fem::str_cref pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, td_real const sparse, td_real const anorm, fem::str_cref pack, td_real *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, td_real *d, mplapackint const mode, td_real const cond, td_real const dmax, mplapackint const kl, mplapackint const ku, fem::str_cref pack, td_real *a, mplapackint const lda, td_real *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, td_real *d, mplapackint const mode, td_real const cond, td_real const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, fem::str_cref pack, td_real *a, mplapackint const lda, td_real *work, mplapackint &info);
#endif
