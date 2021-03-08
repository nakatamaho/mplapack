/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cunmbr.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cunmbr(const char *vect, const char *side, const char *trans, INTEGER m,
	    INTEGER n, INTEGER k, COMPLEX * A, INTEGER lda, COMPLEX * tau, COMPLEX * c, INTEGER ldc, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i1, i2, nb, mi, ni, nq, nw;
    INTEGER left;
    INTEGER iinfo;
    INTEGER notran, applyq;
    char transt;
    INTEGER lwkopt;
    INTEGER lquery;
    char side_trans[3];

//Test the input arguments
    *info = 0;
    applyq = Mlsame(vect, "Q");
    left = Mlsame(side, "L");
    notran = Mlsame(trans, "N");
    lquery = lwork == -1;
//NQ is the order of Q or P and NW is the minimum dimension of WORK
    if (left) {
	nq = m;
	nw = n;
    } else {
	nq = n;
	nw = m;
    }
    if (m == 0 || n == 0) {
	nw = 0;
    }
    if (!applyq && !Mlsame(vect, "P")) {
	*info = -1;
    } else if (!left && !Mlsame(side, "R")) {
	*info = -2;
    } else if (!notran && !Mlsame(trans, "C")) {
	*info = -3;
    } else if (m < 0) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (k < 0) {
	*info = -6;
    } else {
	if ((applyq && lda < max((INTEGER) 1, nq)) || (!applyq && lda < max((INTEGER) 1, min(nq, k)))) {
	    *info = -8;
	} else if (ldc < max((INTEGER) 1, m)) {
	    *info = -11;
	} else if (lwork < max((INTEGER) 1, nw) && !lquery) {
	    *info = -13;
	}
    }
    if (*info == 0) {
	if (nw > 0) {
	    if (applyq) {
		if (left) {
		    side_trans[0] = side[0];
		    side_trans[1] = trans[0];
		    side_trans[2] = '\0';
		    nb = iMlaenv(1, "Cunmqr", side_trans, m - 1, n, m - 1, -1);
		} else {
		    side_trans[0] = side[0];
		    side_trans[1] = trans[0];
		    side_trans[2] = '\0';
		    nb = iMlaenv(1, "Cunmqr", side_trans, m, n - 1, n - 1, -1);
		}
	    } else {
		if (left) {
		    side_trans[0] = side[0];
		    side_trans[1] = trans[0];
		    side_trans[2] = '\0';
		    nb = iMlaenv(1, "Cunmlq", side_trans, m - 1, n, m - 1, -1);
		} else {
		    side_trans[0] = side[0];
		    side_trans[1] = trans[0];
		    side_trans[2] = '\0';
		    nb = iMlaenv(1, "Cunmlq", side_trans, m, n - 1, n - 2, -1);
		}
	    }
	    lwkopt = max((INTEGER) 1, nw * nb);
	} else {
	    lwkopt = 1;
	}
	work[1] = lwkopt;
    }
    if (*info != 0) {
	Mxerbla("Cunmber", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	return;
    }
    if (applyq) {
//Apply Q
	if (nq >= k) {
//Q was determined by a call to ZGEBRD with nq >= k
	    Cunmqr(side, trans, m, n, k, &A[0], lda, &tau[1], &c[0], ldc, &work[0], lwork, &iinfo);
	} else if (nq > 1) {
//Q was determined by a call to ZGEBRD with nq < k
	    if (left) {
		mi = m - 1;
		ni = n;
		i1 = 2;
		i2 = 1;
	    } else {
		mi = m;
		ni = n - 1;
		i1 = 1;
		i2 = 2;
	    }
	    Cunmqr(side, trans, mi, ni, nq - 1, &A[lda + 2], lda, &tau[1]
		   , &c[i1 + i2 * ldc], ldc, &work[0], lwork, &iinfo);
	}
    } else {
//Apply P
	if (notran) {
	    transt = 'C';
	} else {
	    transt = 'N';
	}
	if (nq > k) {
//P was determined by a call to ZGEBRD with nq > k
	    Cunmlq(side, &transt, m, n, k, &A[0], lda, &tau[1], &c[0], ldc, &work[0], lwork, &iinfo);
	} else if (nq > 1) {
//P was determined by a call to ZGEBRD with nq <= k
	    if (left) {
		mi = m - 1;
		ni = n;
		i1 = 2;
		i2 = 1;
	    } else {
		mi = m;
		ni = n - 1;
		i1 = 1;
		i2 = 2;
	    }
	    Cunmlq(side, &transt, mi, ni, nq - 1, &A[(lda * 2) + 1], lda, &tau[1], &c[nq - 1 + i2 * ldc], ldc, &work[0], lwork, &iinfo);
	}
    }
    work[1] = lwkopt;
    return;
}
