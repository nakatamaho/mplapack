/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cunmql.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Cunmql(const char *side, const char *trans, INTEGER m, INTEGER n, INTEGER k, COMPLEX * A,
       INTEGER lda, COMPLEX * tau, COMPLEX * c, INTEGER ldc, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i;
    COMPLEX t[4160];
    INTEGER i1, i2, i3, ib, ic = 0, jc = 0, nb = 0, mi = 0, ni = 0, nq, nw, iws;
    INTEGER left;
    INTEGER nbmin, iinfo;
    INTEGER notran;
    INTEGER ldwork;
    INTEGER lwkopt;
    INTEGER lquery;
    INTEGER ldt = 0;
    char ch[3];

//Test the input arguments
    *info = 0;
    left = Mlsame(side, "L");
    notran = Mlsame(trans, "N");
    lquery = lwork == -1;
//NQ is the order of Q and NW is the minimum dimension of WORK
    if (left) {
	nq = m;
	nw = max((INTEGER) 1, n);
    } else {
	nq = n;
	nw = max((INTEGER) 1, m);
    }
    if (!left && !Mlsame(side, "R")) {
	*info = -1;
    } else if (!notran && !Mlsame(trans, "C")) {
	*info = -2;
    } else if (m < 0) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (k < 0 || k > nq) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, nq)) {
	*info = -7;
    } else if (ldc < max((INTEGER) 1, m)) {
	*info = -10;
    }
    if (*info == 0) {
//Determine the block size.  NB may be at most NBMAX, where
//NBMAX is used to define the local array T.
	ch[0] = (*side);
	ch[1] = (*trans);
	ch[2] = '\0';
	nb = max((INTEGER) 64, iMlaenv(1, "Cunmql", ch, m, n, k, -1));
	lwkopt = nw * nb;
	work[1] = lwkopt;
    }
    if (*info != 0) {
	Mxerbla("Cunmql", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	work[1] = 1;
	return;
    }
    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < k) {
	iws = nw * nb;
	if (lwork < iws) {
	    nb = lwork / ldwork;
	    ch[0] = (*side);
	    ch[1] = (*trans);
	    ch[2] = '\0';
	    nbmin = max((INTEGER) 2, iMlaenv(2, "Cunmql", ch, m, n, k, -1));
	}
    } else {
	iws = nw;
    }
    if (nb < nbmin || nb >= k) {
//Use unblocked code
	Cunm2l(side, trans, m, n, k, A, lda, tau, c, ldc, work, &iinfo);
    } else {
//Use blocked code
	if ((left && notran) || (!left && !notran)) {
	    i1 = 1;
	    i2 = k;
	    i3 = nb;
	} else {
	    i1 = (k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}
	if (left) {
	    ni = n;
	} else {
	    mi = m;
	}
	for (i = i1; i <= i2; i = i + i3) {
	    ib = min(nb, k - i + 1);
//Form the triangular factor of the block reflector
//H or H' is applied to C(1:m-k+i+ib-1,1:n)
	    Clarft("Backward", "Columnwise", nq - k + i + ib + 1, ib, &A[0 + i * lda], lda, &tau[i], t, ldt);
	    if (left) {
//H or H' is applied to C(i:m,1:n)
		mi = m - k + i + ib - 1;
	    } else {
//H or H' is applied to C(1:m,1:n-k+i+ib-1)
		ni = n - k + i + ib - 1;
		jc = i;
	    }
//Apply H or H'
	    Clarfb(side, trans, "Backward", "Columnwise", mi, ni, ib, &A[0 + i * lda], lda, t, ldt, &c[ic + jc * ldc], ldc, work, ldwork);
	}
    }
    work[1] = lwkopt;
    return;
}
