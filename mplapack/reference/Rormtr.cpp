/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rormtr.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/

#include <mpblas.h>
#include <mplapack.h>
#include <string.h>

void
Rormtr(const char *side, const char *uplo, const char *trans, INTEGER m,
       INTEGER n, REAL * A, INTEGER lda, REAL * tau, REAL * c, INTEGER ldc, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i1, i2, nb, mi, ni, nq, nw;
    INTEGER left;
    INTEGER iinfo;
    INTEGER upper;
    INTEGER lwkopt;
    INTEGER lquery;
    char side_trans[3];
    REAL One = 1.0;

//Test the input arguments
    *info = 0;
    left = Mlsame(side, "L");
    upper = Mlsame(uplo, "U");
    lquery = lwork == -1;

//NQ is the order of Q and NW is the minimum dimension of WORK
    if (left) {
	nq = m;
	nw = n;
    } else {
	nq = n;
	nw = m;
    }
    if (!left && !Mlsame(side, "R")) {
	*info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
	*info = -2;
    } else if (!Mlsame(trans, "N") && !Mlsame(trans, "T")) {
	*info = -3;
    } else if (m < 0) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, nq)) {
	*info = -7;
    } else if (ldc < max((INTEGER) 1, m)) {
	*info = -10;
    } else if (lwork < max((INTEGER) 1, nw) && !lquery) {
	*info = -12;
    }

    if (*info == 0) {
	side_trans[0] = side[0];
	side_trans[1] = trans[0];
	side_trans[2] = '\0';
	if (upper) {
	    if (left) {
		nb = iMlaenv(1, "Rormql", side_trans, m - 1, n, m - 1, -1);
	    } else {
		nb = iMlaenv(1, "Rormql", side_trans, m, n - 1, n - 1, -1);
	    }
	} else {
	    if (left) {
		nb = iMlaenv(1, "Rormqr", side_trans, m - 1, n, m - 1, -1);
	    } else {
		nb = iMlaenv(1, "Rormqr", side_trans, m, n - 1, n - 1, -1);
	    }
	}
	lwkopt = max((INTEGER) 1, nw) * nb;
	work[1] = lwkopt;
    }

    if (*info != 0) {
	Mxerbla("Rormtr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0 || nq == 1) {
	work[1] = One;
	return;
    }

    if (left) {
	mi = m - 1;
	ni = n;
    } else {
	mi = m;
	ni = n - 1;
    }
    if (upper) {
//Q was determined by a call to DSYTRD with UPLO = 'U'
	Rormql(side, trans, mi, ni, nq - 1, &A[(lda << 1) + 1], lda, &tau[1], &c[0], ldc, &work[0], lwork, &iinfo);
    } else {

//Q was determined by a call to DSYTRD with UPLO = 'L'
	if (left) {
	    i1 = 2;
	    i2 = 1;
	} else {
	    i1 = 1;
	    i2 = 2;
	}
	Rormqr(side, trans, mi, ni, nq - 1, &A[lda + 2], lda, &tau[1], &c[i1 + i2 * ldc], ldc, &work[0], lwork, &iinfo);
    }
    work[1] = lwkopt;
    return;
}
