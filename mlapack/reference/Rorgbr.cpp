/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rorgbr.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#include <mblas.h>
#include <mlapack.h>

void Rorgbr(const char *vect, INTEGER m, INTEGER n, INTEGER k, REAL * A, INTEGER lda, REAL * tau, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, j, nb, mn;
    INTEGER iinfo;
    INTEGER wantq;
    INTEGER lwkopt;
    INTEGER lquery;
    REAL One = 1.0, Zero = 0.0;

//Test the input arguments
    *info = 0;
    wantq = Mlsame(vect, "Q");
    mn = min(m, n);
    lquery = lwork == -1;
    if (!wantq && !Mlsame(vect, "P")) {
	*info = -1;
    } else if (m < 0) {
	*info = -2;
    } else if (n < 0 || (wantq && (n > m || n < min(m, k))) || (!wantq && (m > n || m < min(n, k)))) {
	*info = -3;
    } else if (k < 0) {
	*info = -4;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -6;
    } else if (lwork < max((INTEGER) 1, mn) && !lquery) {
	*info = -9;
    }

    if (*info == 0) {
	if (wantq) {
	    nb = iMlaenv(1, "Rorgqr", " ", m, n, k, -1);
	} else {
	    nb = iMlaenv(1, "Rorgql", " ", m, n, k, -1);
	}
	lwkopt = max((INTEGER) 1, mn) * nb;
	work[1] = lwkopt;
    }
    if (*info != 0) {
	Mxerbla("Rorgbr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	work[1] = One;
	return;
    }

    if (wantq) {
//Form Q, determined by a call to DGEBRD to reduce an m-by-k
//matrix
	if (m >= k) {
//If m >= k, assume m >= n >= k
	    Rorgqr(m, n, k, &A[0], lda, &tau[1], &work[0], lwork, &iinfo);
	} else {
//If m < k, assume m = n
//Shift the vectors which define the elementary reflectors one
//column to the right, and set the first row and column of Q
//to those of the unit matrix
	    for (j = m; j >= 2; j--) {
		A[j * lda] = Zero;
		for (i = j + 1; i <= m; i++) {
		    A[i + j * lda] = A[i + (j - 1) * lda];
		}
	    }
	    A[lda + 1] = One;
	    for (i = 1; i < m; i++) {
		A[i + lda] = Zero;
	    }
	    if (m > 1) {
//Form Q(2:m,2:m)
		Rorgqr(m - 1, m - 1, m - 1, &A[(lda << 1) + 2], lda, &tau[1], &work[0], lwork, &iinfo);
	    }
	}
    } else {
//Form P', determined by a call to DGEBRD to reduce a k-by-n
//matrix
	if (k < n) {
//If k < n, assume k <= m <= n
	    Rorglq(m, n, k, &A[0], lda, &tau[1], &work[0], lwork, &iinfo);
	} else {
//If k >= n, assume m = n
//Shift the vectors which define the elementary reflectors one
//row downward, and set the first row and column of P' to
//those of the unit matrix
	    A[lda + 1] = One;
	    for (i = 1; i < n; i++) {
		A[i + lda] = Zero;
	    }
	    for (j = 2; j <= n; j++) {
		for (i = j - 1; i >= 2; i--) {
		    A[i + j * lda] = A[i - 1 + j * lda];
		}
		A[j * lda] = Zero;
	    }
	    if (n > 1) {
//Form P'(2:n,2:n)
		Rorglq(n - 1, n - 1, n - 1, &A[(lda << 1) + 2], lda, &tau[1], &work[0], lwork, &iinfo);
	    }
	}
    }
    work[1] = lwkopt;
    return;
}
