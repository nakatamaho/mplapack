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

#ifndef MPLAPACK_MATGEN_MPFR_H
#define MPLAPACK_MATGEN_MPFR_H

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include "mplapack_config.h"
#include <mplapack_gmpfrxx_mkII_config.h>
#include <mpfrxx_mkII.h>
#include <mpcxx_mkII.h>
using namespace mpfrxx;


mpc_class Clarnd(mplapackint const idist, mplapackint (&iseed)[4]);
mpc_class Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], mpc_class *d, mplapackint const igrade, mpc_class *dl, mpc_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpfr_class const sparse);
mpc_class Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], mpc_class *d, mplapackint const igrade, mpc_class *dl, mpc_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpfr_class const sparse);
mplapackint iMlaenv_mpfr(mplapackint const ispec, const char *name, const char *opts, mplapackint const n1, mplapackint const n2, mplapackint const n3, mplapackint const n4);
mplapackint iMlaenv_mpfr2stage(mplapackint const ispec, const char *name, const char *opts, mplapackint const n1, mplapackint const n2, mplapackint const n3, mplapackint const n4);
mpfr_class Rlamch_mpfr(const char *cmach);
mpfr_class Rlaran(mplapackint (&iseed)[4]);
mpfr_class Rlarnd(mplapackint const idist, mplapackint (&iseed)[4]);
mpfr_class Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], mpfr_class *d, mplapackint const igrade, mpfr_class *dl, mpfr_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpfr_class const sparse);
mpfr_class Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], mpfr_class *d, mplapackint const igrade, mpfr_class *dl, mpfr_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpfr_class const sparse);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpfr_class *d, mpc_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpc_class *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, mpfr_class *d, mpc_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpc_class *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, mpfr_class *d, mpc_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpc_class *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const ldx, mpc_class *b, mplapackint const ldb, mpfr_class *work, mplapackint &info, fem::str_cref path);
void Clakf2(mplapackint const m, mplapackint const n, mpc_class *a, mplapackint const lda, mpc_class *b, mpc_class *d, mpc_class *e, mpc_class *z, mplapackint const ldz);
void Clarge(mplapackint const n, mpc_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpc_class *work, mplapackint &info);
void Claror(fem::str_cref side, fem::str_cref init, mplapackint const m, mplapackint const n, mpc_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpc_class *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, mpc_class const c, mpc_class const s, mpc_class *a, mplapackint const lda, mpc_class &xleft, mpc_class &xright);
void Clatm1(mplapackint const mode, mpfr_class const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], mpc_class *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class *c, mplapackint const ldc, mpc_class *d, mplapackint const ldd, mpc_class *e, mplapackint const lde, mpc_class *f, mplapackint const ldf, mpc_class *r, mplapackint const ldr, mpc_class *l, mplapackint const ldl, mpfr_class const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, mpc_class *a, mplapackint const lda, mpc_class *b, mpc_class *x, mplapackint const ldx, mpc_class *y, mplapackint const ldy, mpc_class const alpha, mpc_class const beta, mpc_class const wx, mpc_class const wy, mpfr_class *s, mpfr_class *dif);
void Clatme(mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], mpc_class *d, mplapackint const mode, mpfr_class const cond, mpc_class const dmax, fem::str_cref rsign, fem::str_cref upper, fem::str_cref sim, mpfr_class *ds, mplapackint const modes, mpfr_class const conds, mplapackint const kl, mplapackint const ku, mpfr_class const anorm, mpc_class *a, mplapackint const lda, mpc_class *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpc_class *d, mplapackint const mode, mpfr_class const cond, mpc_class const dmax, fem::str_cref rsign, fem::str_cref grade, mpc_class *dl, mplapackint const model, mpfr_class const condl, mpc_class *dr, mplapackint const moder, mpfr_class const condr, fem::str_cref pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, mpfr_class const sparse, mpfr_class const anorm, fem::str_cref pack, mpc_class *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpfr_class *d, mplapackint const mode, mpfr_class const cond, mpfr_class const dmax, mplapackint const kl, mplapackint const ku, fem::str_cref pack, mpc_class *a, mplapackint const lda, mpc_class *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpfr_class *d, mplapackint const mode, mpfr_class const cond, mpfr_class const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, fem::str_cref pack, mpc_class *a, mplapackint const lda, mpc_class *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpfr_class *d, mpfr_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpfr_class *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, mpfr_class *d, mpfr_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpfr_class *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, mpfr_class *a, mplapackint const lda, mpfr_class *x, mplapackint const ldx, mpfr_class *b, mplapackint const ldb, mpfr_class *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, mpfr_class *a, mplapackint const lda, mpfr_class *b, mpfr_class *d, mpfr_class *e, mpfr_class *z, mplapackint const ldz);
void Rlarge(mplapackint const n, mpfr_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpfr_class *work, mplapackint &info);
void Rlaror(fem::str_cref side, fem::str_cref init, mplapackint const m, mplapackint const n, mpfr_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpfr_class *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, mpfr_class const c, mpfr_class const s, mpfr_class *a, mplapackint const lda, mpfr_class &xleft, mpfr_class &xright);
void Rlatm1(mplapackint const mode, mpfr_class const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], mpfr_class *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, mpfr_class *a, mplapackint const lda, mpfr_class *b, mplapackint const ldb, mpfr_class *c, mplapackint const ldc, mpfr_class *d, mplapackint const ldd, mpfr_class *e, mplapackint const lde, mpfr_class *f, mplapackint const ldf, mpfr_class *r, mplapackint const ldr, mpfr_class *l, mplapackint const ldl, mpfr_class const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, mpfr_class *a, mplapackint const lda, mpfr_class *b, mpfr_class *x, mplapackint const ldx, mpfr_class *y, mplapackint const ldy, mpfr_class const alpha, mpfr_class const beta, mpfr_class const wx, mpfr_class const wy, mpfr_class *s, mpfr_class *dif);
void Rlatm7(mplapackint const mode, mpfr_class const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], mpfr_class *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatme(mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], mpfr_class *d, mplapackint const mode, mpfr_class const cond, mpfr_class const dmax, const char *ei, fem::str_cref rsign, fem::str_cref upper, fem::str_cref sim, mpfr_class *ds, mplapackint const modes, mpfr_class const conds, mplapackint const kl, mplapackint const ku, mpfr_class const anorm, mpfr_class *a, mplapackint const lda, mpfr_class *work, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpfr_class *d, mplapackint const mode, mpfr_class const cond, mpfr_class const dmax, fem::str_cref rsign, fem::str_cref grade, mpfr_class *dl, mplapackint const model, mpfr_class const condl, mpfr_class *dr, mplapackint const moder, mpfr_class const condr, fem::str_cref pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, mpfr_class const sparse, mpfr_class const anorm, fem::str_cref pack, mpfr_class *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpfr_class *d, mplapackint const mode, mpfr_class const cond, mpfr_class const dmax, mplapackint const kl, mplapackint const ku, fem::str_cref pack, mpfr_class *a, mplapackint const lda, mpfr_class *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpfr_class *d, mplapackint const mode, mpfr_class const cond, mpfr_class const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, fem::str_cref pack, mpfr_class *a, mplapackint const lda, mpfr_class *work, mplapackint &info);
#endif
