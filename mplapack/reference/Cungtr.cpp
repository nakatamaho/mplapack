/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cungtr.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cungtr(const char *uplo, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * tau, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER lquery, lwkopt, iinfo, upper, nb;
    INTEGER i, j;
    REAL Zero = 0.0, One = 1.0;

//Test the input arguments
    *info = 0;
    lquery = lwork == -1;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    } else {
	if (lwork < max((INTEGER) 1, n - 1) && !lquery) {
	    *info = -7;
	}
    }
    if (*info == 0) {
	if (upper) {
	    nb = iMlaenv(1, "Cungql", " ", n - 1, n - 1, n - 1, -1);
	} else {
	    nb = iMlaenv(1, "Cungqr", " ", n - 1, n - 1, n - 1, -1);
	}
	lwkopt = max((INTEGER) 1, n - 1) * nb;
	work[0] = lwkopt;
    }
    if (*info != 0) {
	Mxerbla("Cungtr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	work[0] = One;
	return;
    }
    if (upper) {
//Q was determined by a call to Chetrd with UPLO = 'U'
//Shift the vectors which define the elementary reflectors one
//column to the left, and set the last row and column of Q to
//those of the unit matrix
	for (j = 1; j <= n - 1; j++) {
	    for (i = 1; i <= j - 1; i++) {
		A[(i - 1) + (j - 1) * lda] = A[(i - 1) + j * lda];
	    }
	    A[(n - 1) + (j - 1) * lda] = Zero;
	}
	for (i = 1; i <= n - 1; i++) {
	    A[(i - 1) + (n - 1) * lda] = Zero;
	}
	A[(n - 1) + (n - 1) * lda] = One;
//Generate Q(1:n-1,1:n-1)
	Cungql(n - 1, n - 1, n - 1, A, lda, tau, work, lwork, &iinfo);
    } else {
//Q was determined by a call to ZHETRD with UPLO = 'L'.
//Shift the vectors which define the elementary reflectors one
//column to the right, and set the first row and column of Q to
//those of the unit matrix
	for (j = n; j >= 2; j--) {
	    A[0 + (j - 1) * lda] = Zero;
	    for (i = j + 1; i <= n; i++) {
		A[(i - 1) + (j - 1) * lda] = A[(i - 1) + (j - 2) * lda];
	    }
	}
	A[0 + 0 * lda] = One;
	for (i = 2; i <= n; i++) {
	    A[(i - 1) + 0 * lda] = Zero;
	}
	if (n > 1) {
//Generate Q(2:n,2:n)
	    Cungqr(n - 1, n - 1, n - 1, &A[1 + (1 * lda)], lda, tau, work, lwork, &iinfo);
	}
    }
    work[0] = lwkopt;
    return;
}
