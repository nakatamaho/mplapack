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

#ifndef MPLAPACK_MATGEN_QD_H
#define MPLAPACK_MATGEN_QD_H

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include "mplapack_config.h"
#include <qd/qd_real.h>
#include <qd/qd_complex.h>

mplapackint iMlaenv_qd(mplapackint const ispec, const char *name, const char *opts, mplapackint const n1, mplapackint const n2, mplapackint const n3, mplapackint const n4);
mplapackint iMlaenv_qd2stage(mplapackint const ispec, const char *name, const char *opts, mplapackint const n1, mplapackint const n2, mplapackint const n3, mplapackint const n4);
qd_complex Clarnd(mplapackint const idist, mplapackint (&iseed)[4]);
qd_complex Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], qd_complex *d, mplapackint const igrade, qd_complex *dl, qd_complex *dr, mplapackint const ipvtng, mplapackint *iwork, qd_real const sparse);
qd_complex Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], qd_complex *d, mplapackint const igrade, qd_complex *dl, qd_complex *dr, mplapackint const ipvtng, mplapackint *iwork, qd_real const sparse);
qd_real Rlamch_qd(const char *cmach);
qd_real Rlaran(mplapackint (&iseed)[4]);
qd_real Rlarnd(mplapackint const idist, mplapackint (&iseed)[4]);
qd_real Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], qd_real *d, mplapackint const igrade, qd_real *dl, qd_real *dr, mplapackint const ipvtng, mplapackint *iwork, qd_real const sparse);
qd_real Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], qd_real *d, mplapackint const igrade, qd_real *dl, qd_real *dr, mplapackint const ipvtng, mplapackint *iwork, qd_real const sparse);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, qd_real *d, qd_complex *a, mplapackint const lda, mplapackint (&iseed)[4], qd_complex *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, qd_real *d, qd_complex *a, mplapackint const lda, mplapackint (&iseed)[4], qd_complex *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, qd_real *d, qd_complex *a, mplapackint const lda, mplapackint (&iseed)[4], qd_complex *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, qd_complex *a, mplapackint const lda, qd_complex *x, mplapackint const ldx, qd_complex *b, mplapackint const ldb, qd_real *work, mplapackint &info, fem::str_cref path);
void Clakf2(mplapackint const m, mplapackint const n, qd_complex *a, mplapackint const lda, qd_complex *b, qd_complex *d, qd_complex *e, qd_complex *z, mplapackint const ldz);
void Clarge(mplapackint const n, qd_complex *a, mplapackint const lda, mplapackint (&iseed)[4], qd_complex *work, mplapackint &info);
void Claror(fem::str_cref side, fem::str_cref init, mplapackint const m, mplapackint const n, qd_complex *a, mplapackint const lda, mplapackint (&iseed)[4], qd_complex *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, qd_complex const c, qd_complex const s, qd_complex *a, mplapackint const lda, qd_complex &xleft, qd_complex &xright);
void Clatm1(mplapackint const mode, qd_real const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], qd_complex *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, qd_complex *a, mplapackint const lda, qd_complex *b, mplapackint const ldb, qd_complex *c, mplapackint const ldc, qd_complex *d, mplapackint const ldd, qd_complex *e, mplapackint const lde, qd_complex *f, mplapackint const ldf, qd_complex *r, mplapackint const ldr, qd_complex *l, mplapackint const ldl, qd_real const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, qd_complex *a, mplapackint const lda, qd_complex *b, qd_complex *x, mplapackint const ldx, qd_complex *y, mplapackint const ldy, qd_complex const alpha, qd_complex const beta, qd_complex const wx, qd_complex const wy, qd_real *s, qd_real *dif);
void Clatme(mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], qd_complex *d, mplapackint const mode, qd_real const cond, qd_complex const dmax, fem::str_cref rsign, fem::str_cref upper, fem::str_cref sim, qd_real *ds, mplapackint const modes, qd_real const conds, mplapackint const kl, mplapackint const ku, qd_real const anorm, qd_complex *a, mplapackint const lda, qd_complex *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, qd_complex *d, mplapackint const mode, qd_real const cond, qd_complex const dmax, fem::str_cref rsign, fem::str_cref grade, qd_complex *dl, mplapackint const model, qd_real const condl, qd_complex *dr, mplapackint const moder, qd_real const condr, fem::str_cref pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, qd_real const sparse, qd_real const anorm, fem::str_cref pack, qd_complex *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, qd_real *d, mplapackint const mode, qd_real const cond, qd_real const dmax, mplapackint const kl, mplapackint const ku, fem::str_cref pack, qd_complex *a, mplapackint const lda, qd_complex *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, qd_real *d, mplapackint const mode, qd_real const cond, qd_real const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, fem::str_cref pack, qd_complex *a, mplapackint const lda, qd_complex *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, qd_real *d, qd_real *a, mplapackint const lda, mplapackint (&iseed)[4], qd_real *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, qd_real *d, qd_real *a, mplapackint const lda, mplapackint (&iseed)[4], qd_real *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, qd_real *a, mplapackint const lda, qd_real *x, mplapackint const ldx, qd_real *b, mplapackint const ldb, qd_real *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, qd_real *a, mplapackint const lda, qd_real *b, qd_real *d, qd_real *e, qd_real *z, mplapackint const ldz);
void Rlarge(mplapackint const n, qd_real *a, mplapackint const lda, mplapackint (&iseed)[4], qd_real *work, mplapackint &info);
void Rlaror(fem::str_cref side, fem::str_cref init, mplapackint const m, mplapackint const n, qd_real *a, mplapackint const lda, mplapackint (&iseed)[4], qd_real *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, qd_real const c, qd_real const s, qd_real *a, mplapackint const lda, qd_real &xleft, qd_real &xright);
void Rlatm1(mplapackint const mode, qd_real const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], qd_real *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, qd_real *a, mplapackint const lda, qd_real *b, mplapackint const ldb, qd_real *c, mplapackint const ldc, qd_real *d, mplapackint const ldd, qd_real *e, mplapackint const lde, qd_real *f, mplapackint const ldf, qd_real *r, mplapackint const ldr, qd_real *l, mplapackint const ldl, qd_real const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, qd_real *a, mplapackint const lda, qd_real *b, qd_real *x, mplapackint const ldx, qd_real *y, mplapackint const ldy, qd_real const alpha, qd_real const beta, qd_real const wx, qd_real const wy, qd_real *s, qd_real *dif);
void Rlatm7(mplapackint const mode, qd_real const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], qd_real *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatme(mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], qd_real *d, mplapackint const mode, qd_real const cond, qd_real const dmax, const char *ei, fem::str_cref rsign, fem::str_cref upper, fem::str_cref sim, qd_real *ds, mplapackint const modes, qd_real const conds, mplapackint const kl, mplapackint const ku, qd_real const anorm, qd_real *a, mplapackint const lda, qd_real *work, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, qd_real *d, mplapackint const mode, qd_real const cond, qd_real const dmax, fem::str_cref rsign, fem::str_cref grade, qd_real *dl, mplapackint const model, qd_real const condl, qd_real *dr, mplapackint const moder, qd_real const condr, fem::str_cref pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, qd_real const sparse, qd_real const anorm, fem::str_cref pack, qd_real *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, qd_real *d, mplapackint const mode, qd_real const cond, qd_real const dmax, mplapackint const kl, mplapackint const ku, fem::str_cref pack, qd_real *a, mplapackint const lda, qd_real *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, qd_real *d, mplapackint const mode, qd_real const cond, qd_real const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, fem::str_cref pack, qd_real *a, mplapackint const lda, qd_real *work, mplapackint &info);
#endif
