/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rormql.cpp,v 1.10 2010/08/07 04:48:33 nakatamaho Exp $ 
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
Rormql(const char *side, const char *trans, INTEGER m, INTEGER n, INTEGER k, REAL * A, INTEGER lda, REAL * tau, REAL * c, INTEGER ldc, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i;
    REAL t[4160] /* was [65][64] */ ;
    INTEGER i1, i2, i3, ib, nb = 0, mi = 0, ni = 0, nq, nw, iws;
    INTEGER left;
    INTEGER nbmin, iinfo;
    INTEGER notran;
    INTEGER ldwork, lwkopt;
    INTEGER lquery;
    char side_trans[3];

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
    } else if (!notran && !Mlsame(trans, "T")) {
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
	if (m == 0 || n == 0) {
	    lwkopt = 1;
	} else {
//Determine the block size.  NB may be at most NBMAX, where
//NBMAX is used to define the local array T.
	    side_trans[0] = side[0];
	    side_trans[1] = trans[0];
	    side_trans[2] = '\0';
	    nb = min((INTEGER) 64, iMlaenv(1, "Rormql", side_trans, m, n, k, -1));
	    lwkopt = nw * nb;
	}
	work[1] = lwkopt;
	if (lwork < nw && !lquery) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	Mxerbla("Rormql", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	return;
    }
    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < k) {
	iws = nw * nb;
	if (lwork < iws) {
	    nb = lwork / ldwork;
	    side_trans[0] = side[0];
	    side_trans[1] = trans[0];
	    side_trans[2] = '\0';
	    nbmin = max((INTEGER) 2, iMlaenv(2, "Rormql", side_trans, m, n, k, -1));
	}
    } else {
	iws = nw;
    }
    if (nb < nbmin || nb >= k) {
//Use unblocked code
	Rorm2l(side, trans, m, n, k, &A[0], lda, &tau[1], &c[0], ldc, &work[0], &iinfo);
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
	if (left)
	    ni = n;
	else
	    mi = m;

	for (i = i1; i <= i2; i += i3) {
	    ib = min(nb, k - i + 1);
//Form the triangular factor of the block reflector
//H = H(i+ib-1) . . . H(i+1) H(i)
	    Rlarft("Backward", "Columnwise", nq - k + i + ib - 1, ib, &A[i * lda], lda, &tau[i], t, 650);
	    if (left) {
//H or H' is applied to C(1:m-k+i+ib-1,1:n)
		mi = m - k + i + ib - 1;
	    } else {
//H or H' is applied to C(1:m,1:n-k+i+ib-1)
		ni = n - k + i + ib - 1;
	    }
//Apply H or H'
	    Rlarfb(side, trans, "Backward", "Columnwise", mi, ni, ib, &A[i * lda], lda, t, 65, &c[0], ldc, &work[0], ldwork);
	}
    }
    work[0] = lwkopt;
    return;
}
