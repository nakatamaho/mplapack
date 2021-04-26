/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_dd.h,v 1.31 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MPLAPACK_MATGEN_DD_H_
#define _MPLAPACK_MATGEN_DD_H_

#include "mplapack_config.h"
#include "qd/dd_real.h"
#include "dd_complex.h"

void Rchklq(bool *dotype, mplapackint const nm, mplapackint *mval, mplapackint const nn, mplapackint *nval, mplapackint const nnb, mplapackint *nbval, mplapackint *nxval, mplapackint const nrhs, dd_real const thresh, bool const tsterr, mplapackint const nmax, dd_real *a, dd_real *af, dd_real *aq, dd_real *al, dd_real *ac, dd_real *b, dd_real *x, dd_real *xact, dd_real *tau, dd_real *work, dd_real *rwork, mplapackint const nout);
void Rerrlq(const char *path, mplapackint const nunit);
void Rlatb4(const char *path, mplapackint const imat, mplapackint const m, mplapackint const n, char *type, mplapackint &kl, mplapackint &ku, dd_real &anorm, mplapackint &mode, dd_real &cndnum, char *dist);
void Rlqt01(mplapackint const m, mplapackint const n, dd_real *a, dd_real *af, dd_real *q, dd_real *l, mplapackint const lda, dd_real *tau, dd_real *work, mplapackint const lwork, dd_real *rwork, dd_real *result);
#endif
