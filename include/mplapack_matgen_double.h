/*
 * Copyright (c) 2008-2021
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

#ifndef _MPLAPACK_MATGEN_DOUBLE_H_
#define _MPLAPACK_MATGEN_DOUBLE_H_

#include "mplapack_config.h"

double Rlamch_double(const char *cmach);
double Rlaran(mplapackint *iseed);
double Rlarnd(mplapackint const idist, mplapackint *iseed);
double Rlatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, double *d, mplapackint const igrade, double *dl, double *dr, mplapackint const ipvtng, mplapackint *iwork, double const sparse);
double Rlatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, double *d, mplapackint const igrade, double *dl, double *dr, mplapackint const ipvtng, mplapackint *iwork, double const sparse);
mplapackint iMlaenv_double(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
mplapackint iMlaenv_double2stage(mplapackint ispec, const char *name, const char *opts, mplapackint n1, mplapackint n2, mplapackint n3, mplapackint n4);
std::complex<double> Clarnd(mplapackint const idist, mplapackint *iseed);
std::complex<double> Clatm2(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, std::complex<double> *d, mplapackint const igrade, std::complex<double> *dl, std::complex<double> *dr, mplapackint const ipvtng, mplapackint *iwork, double const sparse);
std::complex<double> Clatm3(mplapackint const m, mplapackint const n, mplapackint const i, mplapackint const j, mplapackint &isub, mplapackint &jsub, mplapackint const kl, mplapackint const ku, mplapackint const idist, mplapackint *iseed, std::complex<double> *d, mplapackint const igrade, std::complex<double> *dl, std::complex<double> *dr, mplapackint const ipvtng, mplapackint *iwork, double const sparse);
void Clagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, double *d, std::complex<double> *a, mplapackint const lda, mplapackint *iseed, std::complex<double> *work, mplapackint &info);
void Claghe(mplapackint const n, mplapackint const k, double *d, std::complex<double> *a, mplapackint const lda, mplapackint *iseed, std::complex<double> *work, mplapackint &info);
void Clagsy(mplapackint const n, mplapackint const k, double *d, std::complex<double> *a, mplapackint const lda, mplapackint *iseed, std::complex<double> *work, mplapackint &info);
void Clahilb(mplapackint const n, mplapackint const nrhs, std::complex<double> *a, mplapackint const lda, std::complex<double> *x, mplapackint const ldx, std::complex<double> *b, mplapackint const ldb, double *work, mplapackint &info, const char *path);
void Clakf2(mplapackint const m, mplapackint const n, std::complex<double> *a, mplapackint const lda, std::complex<double> *b, std::complex<double> *d, std::complex<double> *e, std::complex<double> *z, mplapackint const ldz);
void Clarge(mplapackint const n, std::complex<double> *a, mplapackint const lda, mplapackint *iseed, std::complex<double> *work, mplapackint &info);
void Claror(const char *side, const char *init, mplapackint const m, mplapackint const n, std::complex<double> *a, mplapackint const lda, mplapackint *iseed, std::complex<double> *x, mplapackint &info);
void Clarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, std::complex<double> const c, std::complex<double> const s, std::complex<double> *a, mplapackint const lda, std::complex<double> &xleft, std::complex<double> &xright);
void Clatm1(mplapackint const mode, double const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, std::complex<double> *d, mplapackint const n, mplapackint &info);
void Clatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, std::complex<double> *a, mplapackint const lda, std::complex<double> *b, mplapackint const ldb, std::complex<double> *c, mplapackint const ldc, std::complex<double> *d, mplapackint const ldd, std::complex<double> *e, mplapackint const lde, std::complex<double> *f, mplapackint const ldf, std::complex<double> *r, mplapackint const ldr, std::complex<double> *l, mplapackint const ldl, double const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Clatm6(mplapackint const type, mplapackint const n, std::complex<double> *a, mplapackint const lda, std::complex<double> *b, std::complex<double> *x, mplapackint const ldx, std::complex<double> *y, mplapackint const ldy, std::complex<double> const alpha, std::complex<double> const beta, std::complex<double> const wx, std::complex<double> const wy, double *s, double *dif);
void Clatme(mplapackint const n, const char *dist, mplapackint *iseed, std::complex<double> *d, mplapackint const mode, double const cond, std::complex<double> const dmax, const char *rsign, const char *upper, const char *sim, double *ds, mplapackint const modes, double const conds, mplapackint const kl, mplapackint const ku, double const anorm, std::complex<double> *a, mplapackint const lda, std::complex<double> *work, mplapackint &info);
void Clatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, std::complex<double> *d, mplapackint const mode, double const cond, std::complex<double> const dmax, const char *rsign, const char *grade, std::complex<double> *dl, mplapackint const model, double const condl, std::complex<double> *dr, mplapackint const moder, double const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, double const sparse, double const anorm, const char *pack, std::complex<double> *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Clatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, double *d, mplapackint const mode, double const cond, double const dmax, mplapackint const kl, mplapackint const ku, const char *pack, std::complex<double> *a, mplapackint const lda, std::complex<double> *work, mplapackint &info);
void Clatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, double *d, mplapackint const mode, double const cond, double const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, std::complex<double> *a, mplapackint const lda, std::complex<double> *work, mplapackint &info);
void Rlagge(mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, double *d, double *a, mplapackint const lda, mplapackint *iseed, double *work, mplapackint &info);
void Rlagsy(mplapackint const n, mplapackint const k, double *d, double *a, mplapackint const lda, mplapackint *iseed, double *work, mplapackint &info);
void Rlahilb(mplapackint const n, mplapackint const nrhs, double *a, mplapackint const lda, double *x, mplapackint const ldx, double *b, mplapackint const ldb, double *work, mplapackint &info);
void Rlakf2(mplapackint const m, mplapackint const n, double *a, mplapackint const lda, double *b, double *d, double *e, double *z, mplapackint const ldz);
void Rlarge(mplapackint const n, double *a, mplapackint const lda, mplapackint *iseed, double *work, mplapackint &info);
void Rlaror(const char *side, const char *init, mplapackint const m, mplapackint const n, double *a, mplapackint const lda, mplapackint *iseed, double *x, mplapackint &info);
void Rlarot(bool const lrows, bool const lleft, bool const lright, mplapackint const nl, double const c, double const s, double *a, mplapackint const lda, double &xleft, double &xright);
void Rlatm1(mplapackint const mode, double const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, double *d, mplapackint const n, mplapackint &info);
void Rlatm5(mplapackint const prtype, mplapackint const m, mplapackint const n, double *a, mplapackint const lda, double *b, mplapackint const ldb, double *c, mplapackint const ldc, double *d, mplapackint const ldd, double *e, mplapackint const lde, double *f, mplapackint const ldf, double *r, mplapackint const ldr, double *l, mplapackint const ldl, double const alpha, mplapackint &qblcka, mplapackint &qblckb);
void Rlatm6(mplapackint const type, mplapackint const n, double *a, mplapackint const lda, double *b, double *x, mplapackint const ldx, double *y, mplapackint const ldy, double const alpha, double const beta, double const wx, double const wy, double *s, double *dif);
void Rlatm7(mplapackint const mode, double const cond, mplapackint const irsign, mplapackint const idist, mplapackint *iseed, double *d, mplapackint const n, mplapackint const rank, mplapackint &info);
void Rlatme(mplapackint const n, const char *dist, mplapackint *iseed, double *d, mplapackint const mode, double const cond, double const dmax, const char *ei, const char *rsign, const char *upper, const char *sim, double *ds, mplapackint const modes, double const conds, mplapackint const kl, mplapackint const ku, double const anorm, double *a, mplapackint const lda, double *work, mplapackint &info);
void Rlatmr(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, double *d, mplapackint const mode, double const cond, double const dmax, const char *rsign, const char *grade, double *dl, mplapackint const model, double const condl, double *dr, mplapackint const moder, double const condr, const char *pivtng, mplapackint *ipivot, mplapackint const kl, mplapackint const ku, double const sparse, double const anorm, const char *pack, double *a, mplapackint const lda, mplapackint *iwork, mplapackint &info);
void Rlatms(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, double *d, mplapackint const mode, double const cond, double const dmax, mplapackint const kl, mplapackint const ku, const char *pack, double *a, mplapackint const lda, double *work, mplapackint &info);
void Rlatmt(mplapackint const m, mplapackint const n, const char *dist, mplapackint *iseed, const char *sym, double *d, mplapackint const mode, double const cond, double const dmax, mplapackint const rank, mplapackint const kl, mplapackint const ku, const char *pack, double *a, mplapackint const lda, double *work, mplapackint &info);
#endif
