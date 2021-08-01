/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_gmp.h,v 1.34 2010/08/07 03:15:46 nakatamaho Exp $
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

#include "mplapack_config.h"
#include "gmpxx.h"
#include "mpc_class.h"

mpc_class Clarnd(mplapackint const idist, mplapackint *iseed);
mpc_class Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, mpc_class *d, mplapackint const igrade, mpc_class *dl, mpc_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpf_class const sparse);
mpc_class Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, mpc_class *d, mplapackint const igrade, mpc_class *dl, mpc_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpf_class const sparse);
mpf_class Rlamch_gmp(const char *cmach);
mpf_class Rlaran(mplapackint *iseed);
mpf_class Rlarnd(mplapackint const idist, mplapackint *iseed);
mpf_class Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, mpf_class *d, mplapackint const igrade, mpf_class *dl, mpf_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpf_class const sparse);
mpf_class Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, mpf_class *d, mplapackint const igrade, mpf_class *dl, mpf_class *dr, mplapackint const ipvtng, mplapackint *iwork, mpf_class const sparse);
mplapackint iMlaenv_gmp(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
mplapackint iMlaenv_gmp2stage(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpf_class *d, mpc_class *a, mplapackint const lda, mplapackint *iseed, mpc_class *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, mpf_class *d, mpc_class *a, mplapackint const lda, mplapackint *iseed, mpc_class *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, mpf_class *d, mpc_class *a, mplapackint const lda, mplapackint *iseed, mpc_class *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, mpc_class *a, mplapackint const lda, mpc_class *x, mplapackint const ldx, mpc_class *b, mplapackint const ldb, mpf_class *work, mplapackint &info, const char *path);
void Clakf2(mplapackint const m, mplapackint const n, mpc_class *a, mplapackint const lda, mpc_class *b, mpc_class *d, mpc_class *e, mpc_class *z, mplapackint const ldz);
void Clarge(mplapackint const n, mpc_class *a, mplapackint const lda, mplapackint *iseed, mpc_class *work, mplapackint &info);
void Claror(const char *side, const char *init, mplapackint const m, mplapackint const n, mpc_class *a, mplapackint const lda, mplapackint *iseed, mpc_class *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, mpc_class const c, mpc_class const s, mpc_class *a, mplapackint const lda, mpc_class &xleft, mpc_class &xright);
void Clatm1(mplapackint const mode, mpf_class const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, mpc_class *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, mpc_class *a, mplapackint const lda, mpc_class *b, mplapackint const ldb, mpc_class *c, mplapackint const ldc, mpc_class *d, mplapackint const ldd, mpc_class *e, mplapackint const lde, mpc_class *f, mplapackint const ldf, mpc_class *r, mplapackint const ldr, mpc_class *l, mplapackint const ldl, mpf_class const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, mpc_class *a, mplapackint const lda, mpc_class *b, mpc_class *x, mplapackint const ldx, mpc_class *y, mplapackint const ldy, mpc_class const alpha, mpc_class const beta, mpc_class const wx, mpc_class const wy, mpf_class *s, mpf_class *dif);
void Clatme(mplapackint const n, const char *dist, mplapackint *iseed, mpc_class *d, mplapackint const mode, mpf_class const cond, mpc_class const dmax, const char *rsign, const char *upper, const char *sim, mpf_class *ds, mplapackint const modes, mpf_class const conds, mplapackint const kl, mplapackint const ku, mpf_class const anorm, mpc_class *a, mplapackint const lda, mpc_class *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpc_class *d, mplapackint const mode, mpf_class const cond, mpc_class const dmax, const char *rsign, const char *grade, mpc_class *dl, mplapackint const model, mpf_class const condl, mpc_class *dr, mplapackint const moder, mpf_class const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, mpf_class const sparse, mpf_class const anorm, const char *pack, mpc_class *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, mplapackint const kl, mplapackint const ku, const char *pack, mpc_class *a, mplapackint const lda, mpc_class *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, mpc_class *a, mplapackint const lda, mpc_class *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mpf_class *d, mpf_class *a, mplapackint const lda, mplapackint *iseed, mpf_class *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, mpf_class *d, mpf_class *a, mplapackint const lda, mplapackint *iseed, mpf_class *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, mpf_class *a, mplapackint const lda, mpf_class *x, mplapackint const ldx, mpf_class *b, mplapackint const ldb, mpf_class *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *b, mpf_class *d, mpf_class *e, mpf_class *z, mplapackint const ldz);
void Rlarge(mplapackint const n, mpf_class *a, mplapackint const lda, mplapackint *iseed, mpf_class *work, mplapackint &info);
void Rlaror(const char *side, const char *init, mplapackint const m, mplapackint const n, mpf_class *a, mplapackint const lda, mplapackint *iseed, mpf_class *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, mpf_class const c, mpf_class const s, mpf_class *a, mplapackint const lda, mpf_class &xleft, mpf_class &xright);
void Rlatm1(mplapackint const mode, mpf_class const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, mpf_class *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *b, mplapackint const ldb, mpf_class *c, mplapackint const ldc, mpf_class *d, mplapackint const ldd, mpf_class *e, mplapackint const lde, mpf_class *f, mplapackint const ldf, mpf_class *r, mplapackint const ldr, mpf_class *l, mplapackint const ldl, mpf_class const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, mpf_class *a, mplapackint const lda, mpf_class *b, mpf_class *x, mplapackint const ldx, mpf_class *y, mplapackint const ldy, mpf_class const alpha, mpf_class const beta, mpf_class const wx, mpf_class const wy, mpf_class *s, mpf_class *dif);
void Rlatm7(mplapackint const mode, mpf_class const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, mpf_class *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatme(mplapackint const n, const char *dist, mplapackint *iseed, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, const char *ei, const char *rsign, const char *upper, const char *sim, mpf_class *ds, mplapackint const modes, mpf_class const conds, mplapackint const kl, mplapackint const ku, mpf_class const anorm, mpf_class *a, mplapackint const lda, mpf_class *work, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, const char *rsign, const char *grade, mpf_class *dl, mplapackint const model, mpf_class const condl, mpf_class *dr, mplapackint const moder, mpf_class const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, mpf_class const sparse, mpf_class const anorm, const char *pack, mpf_class *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, mplapackint const kl, mplapackint const ku, const char *pack, mpf_class *a, mplapackint const lda, mpf_class *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, mpf_class *d, mplapackint const mode, mpf_class const cond, mpf_class const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, mpf_class *a, mplapackint const lda, mpf_class *work, mplapackint &info);
#endif
