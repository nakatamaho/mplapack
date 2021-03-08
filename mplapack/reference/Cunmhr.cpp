/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cunmhr.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Cunmhr(const char *side, const char *trans, INTEGER m, INTEGER n,
       INTEGER ilo, INTEGER ihi, COMPLEX * A, INTEGER lda, COMPLEX * tau, COMPLEX * c, INTEGER ldc, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i1, i2, nb, mi, nh, ni, nq, nw;
    INTEGER left;
    INTEGER iinfo;
    INTEGER lwkopt;
    INTEGER lquery;
    char ch[3];

//Test the input arguments
    *info = 0;
    nh = ihi - ilo;
    left = Mlsame(side, "L");
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
    } else if (!Mlsame(trans, "N") && !Mlsame(trans, "C")) {
	*info = -2;
    } else if (m < 0) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (ilo < 1 || ilo > max((INTEGER) 1, nq)) {
	*info = -5;
    } else if (ihi < min(ilo, nq) || ihi > nq) {
	*info = -6;
    } else if (lda < max((INTEGER) 1, nq)) {
	*info = -8;
    } else if (ldc < max((INTEGER) 1, m)) {
	*info = -11;
    } else if (lwork < max((INTEGER) 1, nw) && !lquery) {
	*info = -13;
    }
    if (*info == 0) {
	if (left) {
/* Writing concatenation */
	    ch[0] = (*side);
	    ch[1] = (*trans);
	    ch[2] = '\0';
	    nb = iMlaenv(1, "Cunmqr", ch, nh, n, nh, -1);
	} else {
/* Writing concatenation */
	    ch[0] = (*side);
	    ch[1] = (*trans);
	    ch[2] = '\0';
	    nb = iMlaenv(1, "Cunmqr", ch, m, nh, nh, -1);
	}
	lwkopt = max((INTEGER) 1, nw) * nb;
	work[1] = lwkopt;
    }
    if (*info != 0) {
	Mxerbla("Cunmhr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0 || nh == 0) {
	work[1] = 1;
	return;
    }
    if (left) {
	mi = nh;
	ni = n;
	i1 = ilo + 1;
	i2 = 1;
    } else {
	mi = m;
	ni = nh;
	i1 = 1;
	i2 = ilo + 1;
    }
    Cunmqr(side, trans, mi, ni, nh, &A[ilo + 1 + ilo * lda], lda, &tau[ilo], &c[i1 + i2 * ldc], ldc, work, lwork, &iinfo);
    work[1] = lwkopt;
    return;
}
