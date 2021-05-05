/*
 * Copyright (c) 2012
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_double.h,v 1.4 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPLAPACK_MATGEN__FLOAT64X_H_
#define _MPLAPACK_MATGEN__FLOAT64X_H_

#include "mplapack_config.h"

_Float64x Rlamch__Float64x(const char *cmach);
_Float64x Rlaran(mplapackint *iseed);
_Float64x Rlarnd(mplapackint const idist, mplapackint *iseed);
_Float64x Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, _Float64x *d, mplapackint const igrade, _Float64x *dl, _Float64x *dr, mplapackint const ipvtng, mplapackint *iwork, _Float64x const sparse);
_Float64x Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, _Float64x *d, mplapackint const igrade, _Float64x *dl, _Float64x *dr, mplapackint const ipvtng, mplapackint *iwork, _Float64x const sparse);
mplapackint iMlaenv__Float64x(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
mplapackint iMlaenv__Float64x2stage(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
std::complex<_Float64x> Clarnd(mplapackint const idist, mplapackint *iseed);
std::complex<_Float64x> Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, std::complex<_Float64x> *d, mplapackint const igrade, std::complex<_Float64x> *dl, std::complex<_Float64x> *dr, mplapackint const ipvtng, mplapackint *iwork, _Float64x const sparse);
std::complex<_Float64x> Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, std::complex<_Float64x> *d, mplapackint const igrade, std::complex<_Float64x> *dl, std::complex<_Float64x> *dr, mplapackint const ipvtng, mplapackint *iwork, _Float64x const sparse);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, _Float64x *d, std::complex<_Float64x> *a, mplapackint const lda, mplapackint *iseed, std::complex<_Float64x> *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, _Float64x *d, std::complex<_Float64x> *a, mplapackint const lda, mplapackint *iseed, std::complex<_Float64x> *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, _Float64x *d, std::complex<_Float64x> *a, mplapackint const lda, mplapackint *iseed, std::complex<_Float64x> *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, std::complex<_Float64x> *a, mplapackint const lda, std::complex<_Float64x> *x, mplapackint const ldx, std::complex<_Float64x> *b, mplapackint const ldb, _Float64x *work, mplapackint &info, const char *path);
void Clakf2(mplapackint const m, mplapackint const n, std::complex<_Float64x> *a, mplapackint const lda, std::complex<_Float64x> *b, std::complex<_Float64x> *d, std::complex<_Float64x> *e, std::complex<_Float64x> *z, mplapackint const ldz);
void Clarge(mplapackint const n, std::complex<_Float64x> *a, mplapackint const lda, mplapackint *iseed, std::complex<_Float64x> *work, mplapackint &info);
void Claror(const char *side, const char *init, mplapackint const m, mplapackint const n, std::complex<_Float64x> *a, mplapackint const lda, mplapackint *iseed, std::complex<_Float64x> *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, std::complex<_Float64x> const c, std::complex<_Float64x> const s, std::complex<_Float64x> *a, mplapackint const lda, std::complex<_Float64x> &xleft, std::complex<_Float64x> &xright);
void Clatm1(mplapackint const mode, _Float64x const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, std::complex<_Float64x> *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, std::complex<_Float64x> *a, mplapackint const lda, std::complex<_Float64x> *b, mplapackint const ldb, std::complex<_Float64x> *c, mplapackint const ldc, std::complex<_Float64x> *d, mplapackint const ldd, std::complex<_Float64x> *e, mplapackint const lde, std::complex<_Float64x> *f, mplapackint const ldf, std::complex<_Float64x> *r, mplapackint const ldr, std::complex<_Float64x> *l, mplapackint const ldl, _Float64x const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, std::complex<_Float64x> *a, mplapackint const lda, std::complex<_Float64x> *b, std::complex<_Float64x> *x, mplapackint const ldx, std::complex<_Float64x> *y, mplapackint const ldy, std::complex<_Float64x> const alpha, std::complex<_Float64x> const beta, std::complex<_Float64x> const wx, std::complex<_Float64x> const wy, _Float64x *s, _Float64x *dif);
void Clatme(mplapackint const n, const char *dist, mplapackint *iseed, std::complex<_Float64x> *d, mplapackint const mode, _Float64x const cond, std::complex<_Float64x> const dmax, const char *rsign, const char *upper, const char *sim, _Float64x *ds, mplapackint const modes, _Float64x const conds, mplapackint const kl, mplapackint const ku, _Float64x const anorm, std::complex<_Float64x> *a, mplapackint const lda, std::complex<_Float64x> *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, std::complex<_Float64x> *d, mplapackint const mode, _Float64x const cond, std::complex<_Float64x> const dmax, const char *rsign, const char *grade, std::complex<_Float64x> *dl, mplapackint const model, _Float64x const condl, std::complex<_Float64x> *dr, mplapackint const moder, _Float64x const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, _Float64x const sparse, _Float64x const anorm, const char *pack, std::complex<_Float64x> *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, _Float64x *d, mplapackint const mode, _Float64x const cond, _Float64x const dmax, mplapackint const kl, mplapackint const ku, const char *pack, std::complex<_Float64x> *a, mplapackint const lda, std::complex<_Float64x> *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, _Float64x *d, mplapackint const mode, _Float64x const cond, _Float64x const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, std::complex<_Float64x> *a, mplapackint const lda, std::complex<_Float64x> *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, _Float64x *d, _Float64x *a, mplapackint const lda, mplapackint *iseed, _Float64x *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, _Float64x *d, _Float64x *a, mplapackint const lda, mplapackint *iseed, _Float64x *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, _Float64x *a, mplapackint const lda, _Float64x *x, mplapackint const ldx, _Float64x *b, mplapackint const ldb, _Float64x *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, _Float64x *a, mplapackint const lda, _Float64x *b, _Float64x *d, _Float64x *e, _Float64x *z, mplapackint const ldz);
void Rlarge(mplapackint const n, _Float64x *a, mplapackint const lda, mplapackint *iseed, _Float64x *work, mplapackint &info);
void Rlaror(const char *side, const char *init, mplapackint const m, mplapackint const n, _Float64x *a, mplapackint const lda, mplapackint *iseed, _Float64x *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, _Float64x const c, _Float64x const s, _Float64x *a, mplapackint const lda, _Float64x &xleft, _Float64x &xright);
void Rlatm1(mplapackint const mode, _Float64x const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, _Float64x *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, _Float64x *a, mplapackint const lda, _Float64x *b, mplapackint const ldb, _Float64x *c, mplapackint const ldc, _Float64x *d, mplapackint const ldd, _Float64x *e, mplapackint const lde, _Float64x *f, mplapackint const ldf, _Float64x *r, mplapackint const ldr, _Float64x *l, mplapackint const ldl, _Float64x const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, _Float64x *a, mplapackint const lda, _Float64x *b, _Float64x *x, mplapackint const ldx, _Float64x *y, mplapackint const ldy, _Float64x const alpha, _Float64x const beta, _Float64x const wx, _Float64x const wy, _Float64x *s, _Float64x *dif);
void Rlatm7(mplapackint const mode, _Float64x const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, _Float64x *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, _Float64x *d, mplapackint const mode, _Float64x const cond, _Float64x const dmax, const char *rsign, const char *grade, _Float64x *dl, mplapackint const model, _Float64x const condl, _Float64x *dr, mplapackint const moder, _Float64x const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, _Float64x const sparse, _Float64x const anorm, const char *pack, _Float64x *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, _Float64x *d, mplapackint const mode, _Float64x const cond, _Float64x const dmax, mplapackint const kl, mplapackint const ku, const char *pack, _Float64x *a, mplapackint const lda, _Float64x *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, _Float64x *d, mplapackint const mode, _Float64x const cond, _Float64x const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, _Float64x *a, mplapackint const lda, _Float64x *work, mplapackint &info);
#endif
