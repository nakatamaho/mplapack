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

#ifndef MPLAPACK_MATGEN_BINARY128_H
#define MPLAPACK_MATGEN_BINARY128_H

#include "mplapack_config.h"

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

mplapack_binary128_t Rlamch_binary128(const char *cmach);
mplapack_binary128_t Rlaran(mplapackint (&iseed)[4]);
mplapack_binary128_t Rlarnd(mplapackint const idist, mplapackint (&iseed)[4]);
mplapack_binary128_t Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], mplapack_binary128_t *d, mplapackint const igrade, mplapack_binary128_t *dl, mplapack_binary128_t *dr, mplapackint const ipvtng, mplapackint *iwork, mplapack_binary128_t const sparse);
mplapack_binary128_t Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], mplapack_binary128_t *d, mplapackint const igrade, mplapack_binary128_t *dl, mplapack_binary128_t *dr, mplapackint const ipvtng, mplapackint *iwork, mplapack_binary128_t const sparse);
mplapackint iMlaenv_binary128(mplapackint const ispec, const char *name, const char *opts, mplapackint const n1, mplapackint const n2, mplapackint const n3, mplapackint const n4);
mplapackint iMlaenv_binary1282stage(mplapackint const ispec, const char *name, const char *opts, mplapackint const n1, mplapackint const n2, mplapackint const n3, mplapackint const n4);
std::complex<mplapack_binary128_t> Clarnd(mplapackint const idist, mplapackint (&iseed)[4]);
std::complex<mplapack_binary128_t> Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], std::complex<mplapack_binary128_t> *d, mplapackint const igrade, std::complex<mplapack_binary128_t> *dl, std::complex<mplapack_binary128_t> *dr, mplapackint const ipvtng, mplapackint *iwork, mplapack_binary128_t const sparse);
std::complex<mplapack_binary128_t> Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint (&iseed)[4], std::complex<mplapack_binary128_t> *d, mplapackint const igrade, std::complex<mplapack_binary128_t> *dl, std::complex<mplapack_binary128_t> *dr, mplapackint const ipvtng, mplapackint *iwork, mplapack_binary128_t const sparse);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mplapack_binary128_t *d, std::complex<mplapack_binary128_t> *a, mplapackint const lda, mplapackint (&iseed)[4], std::complex<mplapack_binary128_t> *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, mplapack_binary128_t *d, std::complex<mplapack_binary128_t> *a, mplapackint const lda, mplapackint (&iseed)[4], std::complex<mplapack_binary128_t> *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, mplapack_binary128_t *d, std::complex<mplapack_binary128_t> *a, mplapackint const lda, mplapackint (&iseed)[4], std::complex<mplapack_binary128_t> *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, std::complex<mplapack_binary128_t> *a, mplapackint const lda, std::complex<mplapack_binary128_t> *x, mplapackint const ldx, std::complex<mplapack_binary128_t> *b, mplapackint const ldb, mplapack_binary128_t *work, mplapackint &info, fem::str_cref path);
void Clakf2(mplapackint const m, mplapackint const n, std::complex<mplapack_binary128_t> *a, mplapackint const lda, std::complex<mplapack_binary128_t> *b, std::complex<mplapack_binary128_t> *d, std::complex<mplapack_binary128_t> *e, std::complex<mplapack_binary128_t> *z, mplapackint const ldz);
void Clarge(mplapackint const n, std::complex<mplapack_binary128_t> *a, mplapackint const lda, mplapackint (&iseed)[4], std::complex<mplapack_binary128_t> *work, mplapackint &info);
void Claror(fem::str_cref side, fem::str_cref init, mplapackint const m, mplapackint const n, std::complex<mplapack_binary128_t> *a, mplapackint const lda, mplapackint (&iseed)[4], std::complex<mplapack_binary128_t> *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, std::complex<mplapack_binary128_t> const c, std::complex<mplapack_binary128_t> const s, std::complex<mplapack_binary128_t> *a, mplapackint const lda, std::complex<mplapack_binary128_t> &xleft, std::complex<mplapack_binary128_t> &xright);
void Clatm1(mplapackint const mode, mplapack_binary128_t const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], std::complex<mplapack_binary128_t> *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, std::complex<mplapack_binary128_t> *a, mplapackint const lda, std::complex<mplapack_binary128_t> *b, mplapackint const ldb, std::complex<mplapack_binary128_t> *c, mplapackint const ldc, std::complex<mplapack_binary128_t> *d, mplapackint const ldd, std::complex<mplapack_binary128_t> *e, mplapackint const lde, std::complex<mplapack_binary128_t> *f, mplapackint const ldf, std::complex<mplapack_binary128_t> *r, mplapackint const ldr, std::complex<mplapack_binary128_t> *l, mplapackint const ldl, mplapack_binary128_t const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, std::complex<mplapack_binary128_t> *a, mplapackint const lda, std::complex<mplapack_binary128_t> *b, std::complex<mplapack_binary128_t> *x, mplapackint const ldx, std::complex<mplapack_binary128_t> *y, mplapackint const ldy, std::complex<mplapack_binary128_t> const alpha, std::complex<mplapack_binary128_t> const beta, std::complex<mplapack_binary128_t> const wx, std::complex<mplapack_binary128_t> const wy, mplapack_binary128_t *s, mplapack_binary128_t *dif);
void Clatme(mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], std::complex<mplapack_binary128_t> *d, mplapackint const mode, mplapack_binary128_t const cond, std::complex<mplapack_binary128_t> const dmax, fem::str_cref rsign, fem::str_cref upper, fem::str_cref sim, mplapack_binary128_t *ds, mplapackint const modes, mplapack_binary128_t const conds, mplapackint const kl, mplapackint const ku, mplapack_binary128_t const anorm, std::complex<mplapack_binary128_t> *a, mplapackint const lda, std::complex<mplapack_binary128_t> *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, std::complex<mplapack_binary128_t> *d, mplapackint const mode, mplapack_binary128_t const cond, std::complex<mplapack_binary128_t> const dmax, fem::str_cref rsign, fem::str_cref grade, std::complex<mplapack_binary128_t> *dl, mplapackint const model, mplapack_binary128_t const condl, std::complex<mplapack_binary128_t> *dr, mplapackint const moder, mplapack_binary128_t const condr, fem::str_cref pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, mplapack_binary128_t const sparse, mplapack_binary128_t const anorm, fem::str_cref pack, std::complex<mplapack_binary128_t> *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mplapack_binary128_t *d, mplapackint const mode, mplapack_binary128_t const cond, mplapack_binary128_t const dmax, mplapackint const kl, mplapackint const ku, fem::str_cref pack, std::complex<mplapack_binary128_t> *a, mplapackint const lda, std::complex<mplapack_binary128_t> *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mplapack_binary128_t *d, mplapackint const mode, mplapack_binary128_t const cond, mplapack_binary128_t const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, fem::str_cref pack, std::complex<mplapack_binary128_t> *a, mplapackint const lda, std::complex<mplapack_binary128_t> *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mplapack_binary128_t *d, mplapack_binary128_t *a, mplapackint const lda, mplapackint (&iseed)[4], mplapack_binary128_t *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, mplapack_binary128_t *d, mplapack_binary128_t *a, mplapackint const lda, mplapackint (&iseed)[4], mplapack_binary128_t *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, mplapack_binary128_t *a, mplapackint const lda, mplapack_binary128_t *x, mplapackint const ldx, mplapack_binary128_t *b, mplapackint const ldb, mplapack_binary128_t *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, mplapack_binary128_t *a, mplapackint const lda, mplapack_binary128_t *b, mplapack_binary128_t *d, mplapack_binary128_t *e, mplapack_binary128_t *z, mplapackint const ldz);
void Rlarge(mplapackint const n, mplapack_binary128_t *a, mplapackint const lda, mplapackint (&iseed)[4], mplapack_binary128_t *work, mplapackint &info);
void Rlaror(fem::str_cref side, fem::str_cref init, mplapackint const m, mplapackint const n, mplapack_binary128_t *a, mplapackint const lda, mplapackint (&iseed)[4], mplapack_binary128_t *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, mplapack_binary128_t const c, mplapack_binary128_t const s, mplapack_binary128_t *a, mplapackint const lda, mplapack_binary128_t &xleft, mplapack_binary128_t &xright);
void Rlatm1(mplapackint const mode, mplapack_binary128_t const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], mplapack_binary128_t *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, mplapack_binary128_t *a, mplapackint const lda, mplapack_binary128_t *b, mplapackint const ldb, mplapack_binary128_t *c, mplapackint const ldc, mplapack_binary128_t *d, mplapackint const ldd, mplapack_binary128_t *e, mplapackint const lde, mplapack_binary128_t *f, mplapackint const ldf, mplapack_binary128_t *r, mplapackint const ldr, mplapack_binary128_t *l, mplapackint const ldl, mplapack_binary128_t const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, mplapack_binary128_t *a, mplapackint const lda, mplapack_binary128_t *b, mplapack_binary128_t *x, mplapackint const ldx, mplapack_binary128_t *y, mplapackint const ldy, mplapack_binary128_t const alpha, mplapack_binary128_t const beta, mplapack_binary128_t const wx, mplapack_binary128_t const wy, mplapack_binary128_t *s, mplapack_binary128_t *dif);
void Rlatm7(mplapackint const mode, mplapack_binary128_t const cond, mplapackint const irsign, mplapackint const idist, mplapackint (&iseed)[4], mplapack_binary128_t *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatme(mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], mplapack_binary128_t *d, mplapackint const mode, mplapack_binary128_t const cond, mplapack_binary128_t const dmax, const char *ei, fem::str_cref rsign, fem::str_cref upper, fem::str_cref sim, mplapack_binary128_t *ds, mplapackint const modes, mplapack_binary128_t const conds, mplapackint const kl, mplapackint const ku, mplapack_binary128_t const anorm, mplapack_binary128_t *a, mplapackint const lda, mplapack_binary128_t *work, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mplapack_binary128_t *d, mplapackint const mode, mplapack_binary128_t const cond, mplapack_binary128_t const dmax, fem::str_cref rsign, fem::str_cref grade, mplapack_binary128_t *dl, mplapackint const model, mplapack_binary128_t const condl, mplapack_binary128_t *dr, mplapackint const moder, mplapack_binary128_t const condr, fem::str_cref pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, mplapack_binary128_t const sparse, mplapack_binary128_t const anorm, fem::str_cref pack, mplapack_binary128_t *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mplapack_binary128_t *d, mplapackint const mode, mplapack_binary128_t const cond, mplapack_binary128_t const dmax, mplapackint const kl, mplapackint const ku, fem::str_cref pack, mplapack_binary128_t *a, mplapackint const lda, mplapack_binary128_t *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, fem::str_cref dist, mplapackint (&iseed)[4], fem::str_cref sym, mplapack_binary128_t *d, mplapackint const mode, mplapack_binary128_t const cond, mplapack_binary128_t const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, fem::str_cref pack, mplapack_binary128_t *a, mplapackint const lda, mplapack_binary128_t *work, mplapackint &info);
#endif
