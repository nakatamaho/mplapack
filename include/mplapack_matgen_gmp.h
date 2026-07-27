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

#ifndef _MPLAPACK_MATGEN_GMP_H_
#define _MPLAPACK_MATGEN_GMP_H_

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include "mplapack_config.h"
#include <mplapack_gmpfrxx_mkII_config.h>
#include <gmpxx_mkII.h>
using namespace gmpxx;

mpfc_class Clarnd(mplapackint const idist, mplapackint (&iseed)[4]);
mpfc_class Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], mpfc_class *d, mplapackint const igrade, mpfc_class *dl, mpfc_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpf_class const sparse);
mpfc_class Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], mpfc_class *d, mplapackint const igrade, mpfc_class *dl, mpfc_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpf_class const sparse);
mpf_class Rlamch_gmp(const char *cmach);
mpf_class Rlaran(mplapackint (&iseed)[4]);
mpf_class Rlarnd(mplapackint const idist, mplapackint (&iseed)[4]);
mpf_class Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], mpf_class *d, mplapackint const igrade, mpf_class *dl, mpf_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpf_class const sparse);
mpf_class Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], mpf_class *d, mplapackint const igrade, mpf_class *dl, mpf_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpf_class const sparse);
mplapackint iMlaenv_gmp(mplapackint const ispec, const char *name, const char *opts, mplapackint const n1, mplapackint const n2, mplapackint const n3, mplapackint const n4);
mplapackint iMlaenv_gmp2stage(mplapackint const ispec, const char *name, const char *opts, mplapackint const n1, mplapackint const n2, mplapackint const n3, mplapackint const n4);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpf_class *d, mpfc_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpfc_class *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, mpf_class *d, mpfc_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpfc_class *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, mpf_class *d, mpfc_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpfc_class *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, mpfc_class *a, mplapackint const lda, mpfc_class *x, mplapackint const ldx, mpfc_class *b, mplapackint const ldb, mpf_class *work, mplapackint &info, fem::str_cref path);
void Clakf2(mplapackint const m, mplapackint const n, mpfc_class *a, mplapackint const lda, mpfc_class *b, mpfc_class *d, mpfc_class *e, mpfc_class *z, mplapackint const ldz);
void Clarge(mplapackint const n, mpfc_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpfc_class *work, mplapackint &info);
void Claror(fem::str_cref side, fem::str_cref init, mplapackint const m, mplapackint const n, mpfc_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpfc_class *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, mpfc_class const c, mpfc_class const s, mpfc_class *a, mplapackint const lda, mpfc_class &xleft, mpfc_class &xright);
void Clatm1(mplapackint const mode, mpf_class const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], mpfc_class *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, mpfc_class *a, mplapackint const lda, mpfc_class *b, mplapackint const ldb, mpfc_class *c, mplapackint const ldc, mpfc_class *d, mplapackint const ldd, mpfc_class *e, mplapackint const lde, mpfc_class *f, mplapackint const ldf, mpfc_class *r, mplapackint const ldr, mpfc_class *l, mplapackint const ldl, mpf_class const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, mpfc_class *a, mplapackint const lda, mpfc_class *b, mpfc_class *x, mplapackint const ldx, mpfc_class *y, mplapackint const ldy, mpfc_class const alpha, mpfc_class const beta, mpfc_class const wx, mpfc_class const wy, mpf_class *s, mpf_class *dif);
void Clatme(mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], mpfc_class *d, mplapackint const mode, mpf_class const cond, mpfc_class const dmax, fem::str_cref rsign, fem::str_cref upper, fem::str_cref sim, mpf_class *ds, mplapackint const modes, mpf_class const conds, mplapackint const kl, mplapackint const ku, mpf_class const anorm, mpfc_class *a, mplapackint const lda, mpfc_class *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpfc_class *d, mplapackint const mode, mpf_class const cond, mpfc_class const dmax, fem::str_cref rsign, fem::str_cref grade, mpfc_class *dl, mplapackint const model, mpf_class const condl, mpfc_class *dr, mplapackint const moder, mpf_class const condr, fem::str_cref pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, mpf_class const sparse, mpf_class const anorm, fem::str_cref pack, mpfc_class *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, mplapackint const kl, mplapackint const ku, fem::str_cref pack, mpfc_class *a, mplapackint const lda, mpfc_class *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, fem::str_cref pack, mpfc_class *a, mplapackint const lda, mpfc_class *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpf_class *d, mpf_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpf_class *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, mpf_class *d, mpf_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpf_class *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const ldx, mpf_class *b, mplapackint const ldb, mpf_class *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *b, mpf_class *d, mpf_class *e, mpf_class *z, mplapackint const ldz);
void Rlarge(mplapackint const n, mpf_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpf_class *work, mplapackint &info);
void Rlaror(fem::str_cref side, fem::str_cref init, mplapackint const m, mplapackint const n, mpf_class *a, mplapackint const lda, mplapackint (&iseed)[4], mpf_class *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, mpf_class const c, mpf_class const s, mpf_class *a, mplapackint const lda, mpf_class &xleft, mpf_class &xright);
void Rlatm1(mplapackint const mode, mpf_class const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], mpf_class *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb, mpf_class *c, mplapackint const ldc, mpf_class *d, mplapackint const ldd, mpf_class *e, mplapackint const lde, mpf_class *f, mplapackint const ldf, mpf_class *r, mplapackint const ldr, mpf_class *l, mplapackint const ldl, mpf_class const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *b, mpf_class *x, mplapackint const ldx, mpf_class *y, mplapackint const ldy, mpf_class const alpha, mpf_class const beta, mpf_class const wx, mpf_class const wy, mpf_class *s, mpf_class *dif);
void Rlatm7(mplapackint const mode, mpf_class const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], mpf_class *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatme(mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, const char *ei, fem::str_cref rsign, fem::str_cref upper, fem::str_cref sim, mpf_class *ds, mplapackint const modes, mpf_class const conds, mplapackint const kl, mplapackint const ku, mpf_class const anorm, mpf_class *a, mplapackint const lda, mpf_class *work, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, fem::str_cref rsign, fem::str_cref grade, mpf_class *dl, mplapackint const model, mpf_class const condl, mpf_class *dr, mplapackint const moder, mpf_class const condr, fem::str_cref pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, mpf_class const sparse, mpf_class const anorm, fem::str_cref pack, mpf_class *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, mplapackint const kl, mplapackint const ku, fem::str_cref pack, mpf_class *a, mplapackint const lda, mpf_class *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, fem::str_cref pack, mpf_class *a, mplapackint const lda, mpf_class *work, mplapackint &info);
#endif
