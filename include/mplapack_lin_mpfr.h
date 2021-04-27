/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_mpfr.h,v 1.5 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPLAPACK_LIN_MPFR_H_
#define _MPLAPACK_LIN_MPFR_H_

#include "mplapack_config.h"
#include "mpc_class.h"
#include "gmpxx.h"
#include "mpreal.h"
#include "mpcomplex.h"

using namespace mpfr;

void Alaerh(const char *path, const char *subnam, mplapackint const info, mplapackint const infoe, const char *opts, mplapackint const m, mplapackint const n, mplapackint const kl, mplapackint const ku, mplapackint const n5, mplapackint const imat, mplapackint const nfail, mplapackint &nerrs, mplapackint const nout);
void Alareq(const char *path, mplapackint const nmats, bool *dotype, mplapackint const ntypes, mplapackint const nin, mplapackint const nout);
void Alasum(const char *type, mplapackint const nout, mplapackint const nfail, mplapackint const nrun, mplapackint const nerrs);
void Rchkaa(void);
void Rchkge(bool *dotype, mplapackint const nm, mplapackint *mval, mplapackint const nn, mplapackint *nval, mplapackint const nnb, mplapackint *nbval, mplapackint const nns, mplapackint *nsval, mpreal const thresh, bool const tsterr, mplapackint const nmax, mpreal *a, mpreal *afac, mpreal *ainv, mpreal *b, mpreal *x, mpreal *xact, mpreal *work, mpreal *rwork, mplapackint *iwork, mplapackint const nout);
void Rchklq(bool *dotype, mplapackint const nm, mplapackint *mval, mplapackint const nn, mplapackint *nval, mplapackint const nnb, mplapackint *nbval, mplapackint *nxval, mplapackint const nrhs, mpreal const thresh, bool const tsterr, mplapackint const nmax, mpreal *a, mpreal *af, mpreal *aq, mpreal *al, mpreal *ac, mpreal *b, mpreal *x, mpreal *xact, mpreal *tau, mpreal *work, mpreal *rwork, mplapackint const nout);
void Rerrge(const char *path, mplapackint const nunit);
void Rerrlq(const char *path, mplapackint const nunit);
void Rget01(mplapackint const m, mplapackint const n, mpreal *a, mplapackint const lda, mpreal *afac, mplapackint const ldafac, mplapackint *ipiv, mpreal *rwork, mpreal &resid);
void Rlatb4(const char *path, mplapackint const imat, mplapackint const m, mplapackint const n, char *type, mplapackint &kl, mplapackint &ku, mpreal &anorm, mplapackint &mode, mpreal &cndnum, char *dist);
void Rlqt01(mplapackint const m, mplapackint const n, mpreal *a, mpreal *af, mpreal *q, mpreal *l, mplapackint const lda, mpreal *tau, mpreal *work, mplapackint const lwork, mpreal *rwork, mpreal *result);
void chkxer(const char *srnamt, mplapackint const &infot, mplapackint const &nout, bool &lerr, bool &ok);
#endif
