/*
 * Copyright (c) 2012-2021
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

#ifndef _MPLAPACK_MATGEN__FLOAT128_H_
#define _MPLAPACK_MATGEN__FLOAT128_H_

#include "mplapack_config.h"

#if defined ___MPLAPACK_MATGEN_WANT_LIBQUADMATH___
#include "quadmath.h"
#endif

_Float128 Rlamch__Float128(const char *cmach);
_Float128 Rlaran(mplapackint *iseed);
_Float128 Rlarnd(mplapackint const idist, mplapackint *iseed);
_Float128 Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, _Float128 *d, mplapackint const igrade, _Float128 *dl, _Float128 *dr, mplapackint const ipvtng, mplapackint *iwork, _Float128 const sparse);
_Float128 Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, _Float128 *d, mplapackint const igrade, _Float128 *dl, _Float128 *dr, mplapackint const ipvtng, mplapackint *iwork, _Float128 const sparse);
mplapackint iMlaenv__Float128(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
mplapackint iMlaenv__Float1282stage(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
std::complex<_Float128> Clarnd(mplapackint const idist, mplapackint *iseed);
std::complex<_Float128> Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, std::complex<_Float128> *d, mplapackint const igrade, std::complex<_Float128> *dl, std::complex<_Float128> *dr, mplapackint const ipvtng, mplapackint *iwork, _Float128 const sparse);
std::complex<_Float128> Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, std::complex<_Float128> *d, mplapackint const igrade, std::complex<_Float128> *dl, std::complex<_Float128> *dr, mplapackint const ipvtng, mplapackint *iwork, _Float128 const sparse);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, _Float128 *d, std::complex<_Float128> *a, mplapackint const lda, mplapackint *iseed, std::complex<_Float128> *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, _Float128 *d, std::complex<_Float128> *a, mplapackint const lda, mplapackint *iseed, std::complex<_Float128> *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, _Float128 *d, std::complex<_Float128> *a, mplapackint const lda, mplapackint *iseed, std::complex<_Float128> *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, std::complex<_Float128> *a, mplapackint const lda, std::complex<_Float128> *x, mplapackint const ldx, std::complex<_Float128> *b, mplapackint const ldb, _Float128 *work, mplapackint &info, const char *path);
void Clakf2(mplapackint const m, mplapackint const n, std::complex<_Float128> *a, mplapackint const lda, std::complex<_Float128> *b, std::complex<_Float128> *d, std::complex<_Float128> *e, std::complex<_Float128> *z, mplapackint const ldz);
void Clarge(mplapackint const n, std::complex<_Float128> *a, mplapackint const lda, mplapackint *iseed, std::complex<_Float128> *work, mplapackint &info);
void Claror(const char *side, const char *init, mplapackint const m, mplapackint const n, std::complex<_Float128> *a, mplapackint const lda, mplapackint *iseed, std::complex<_Float128> *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, std::complex<_Float128> const c, std::complex<_Float128> const s, std::complex<_Float128> *a, mplapackint const lda, std::complex<_Float128> &xleft, std::complex<_Float128> &xright);
void Clatm1(mplapackint const mode, _Float128 const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, std::complex<_Float128> *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, std::complex<_Float128> *a, mplapackint const lda, std::complex<_Float128> *b, mplapackint const ldb, std::complex<_Float128> *c, mplapackint const ldc, std::complex<_Float128> *d, mplapackint const ldd, std::complex<_Float128> *e, mplapackint const lde, std::complex<_Float128> *f, mplapackint const ldf, std::complex<_Float128> *r, mplapackint const ldr, std::complex<_Float128> *l, mplapackint const ldl, _Float128 const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, std::complex<_Float128> *a, mplapackint const lda, std::complex<_Float128> *b, std::complex<_Float128> *x, mplapackint const ldx, std::complex<_Float128> *y, mplapackint const ldy, std::complex<_Float128> const alpha, std::complex<_Float128> const beta, std::complex<_Float128> const wx, std::complex<_Float128> const wy, _Float128 *s, _Float128 *dif);
void Clatme(mplapackint const n, const char *dist, mplapackint *iseed, std::complex<_Float128> *d, mplapackint const mode, _Float128 const cond, std::complex<_Float128> const dmax, const char *rsign, const char *upper, const char *sim, _Float128 *ds, mplapackint const modes, _Float128 const conds, mplapackint const kl, mplapackint const ku, _Float128 const anorm, std::complex<_Float128> *a, mplapackint const lda, std::complex<_Float128> *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, std::complex<_Float128> *d, mplapackint const mode, _Float128 const cond, std::complex<_Float128> const dmax, const char *rsign, const char *grade, std::complex<_Float128> *dl, mplapackint const model, _Float128 const condl, std::complex<_Float128> *dr, mplapackint const moder, _Float128 const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, _Float128 const sparse, _Float128 const anorm, const char *pack, std::complex<_Float128> *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, _Float128 *d, mplapackint const mode, _Float128 const cond, _Float128 const dmax, mplapackint const kl, mplapackint const ku, const char *pack, std::complex<_Float128> *a, mplapackint const lda, std::complex<_Float128> *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, _Float128 *d, mplapackint const mode, _Float128 const cond, _Float128 const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, std::complex<_Float128> *a, mplapackint const lda, std::complex<_Float128> *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, _Float128 *d, _Float128 *a, mplapackint const lda, mplapackint *iseed, _Float128 *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, _Float128 *d, _Float128 *a, mplapackint const lda, mplapackint *iseed, _Float128 *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, _Float128 *a, mplapackint const lda, _Float128 *x, mplapackint const ldx, _Float128 *b, mplapackint const ldb, _Float128 *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, _Float128 *a, mplapackint const lda, _Float128 *b, _Float128 *d, _Float128 *e, _Float128 *z, mplapackint const ldz);
void Rlarge(mplapackint const n, _Float128 *a, mplapackint const lda, mplapackint *iseed, _Float128 *work, mplapackint &info);
void Rlaror(const char *side, const char *init, mplapackint const m, mplapackint const n, _Float128 *a, mplapackint const lda, mplapackint *iseed, _Float128 *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, _Float128 const c, _Float128 const s, _Float128 *a, mplapackint const lda, _Float128 &xleft, _Float128 &xright);
void Rlatm1(mplapackint const mode, _Float128 const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, _Float128 *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, _Float128 *a, mplapackint const lda, _Float128 *b, mplapackint const ldb, _Float128 *c, mplapackint const ldc, _Float128 *d, mplapackint const ldd, _Float128 *e, mplapackint const lde, _Float128 *f, mplapackint const ldf, _Float128 *r, mplapackint const ldr, _Float128 *l, mplapackint const ldl, _Float128 const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, _Float128 *a, mplapackint const lda, _Float128 *b, _Float128 *x, mplapackint const ldx, _Float128 *y, mplapackint const ldy, _Float128 const alpha, _Float128 const beta, _Float128 const wx, _Float128 const wy, _Float128 *s, _Float128 *dif);
void Rlatm7(mplapackint const mode, _Float128 const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, _Float128 *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatme(mplapackint const n, const char *dist, mplapackint *iseed, _Float128 *d, mplapackint const mode, _Float128 const cond, _Float128 const dmax, const char *ei, const char *rsign, const char *upper, const char *sim, _Float128 *ds, mplapackint const modes, _Float128 const conds, mplapackint const kl, mplapackint const ku, _Float128 const anorm, _Float128 *a, mplapackint const lda, _Float128 *work, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, _Float128 *d, mplapackint const mode, _Float128 const cond, _Float128 const dmax, const char *rsign, const char *grade, _Float128 *dl, mplapackint const model, _Float128 const condl, _Float128 *dr, mplapackint const moder, _Float128 const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, _Float128 const sparse, _Float128 const anorm, const char *pack, _Float128 *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, _Float128 *d, mplapackint const mode, _Float128 const cond, _Float128 const dmax, mplapackint const kl, mplapackint const ku, const char *pack, _Float128 *a, mplapackint const lda, _Float128 *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, _Float128 *d, mplapackint const mode, _Float128 const cond, _Float128 const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, _Float128 *a, mplapackint const lda, _Float128 *work, mplapackint &info);
#endif
