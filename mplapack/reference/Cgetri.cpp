/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgetri.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgetri(INTEGER n, COMPLEX * A, INTEGER lda, INTEGER * ipiv, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, j, jb, nb, jj, jp, nn, iws, nbmin;
    INTEGER ldwork, lwkopt;
    INTEGER lquery;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters.
    *info = 0;
    nb = iMlaenv(1, "Cgetri", " ", n, -1, -1, -1);
    lwkopt = n * nb;
    work[0] = lwkopt;
    lquery = lwork == -1;
    if (n < 0) {
	*info = -1;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -3;
    } else if (lwork < max((INTEGER) 1, n) && !lquery) {
	*info = -6;
    }
    if (*info != 0) {
	Mxerbla("Cgetri", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Form inv(U). If INFO > 0 from Ctrtri, then U is singular,
//and the inverse is not computed.
    Ctrtri("Upper", "Non-unit", n, A, lda, info);
    if (*info > 0) {
	return;
    }
    nbmin = 2;
    ldwork = n;
    if (nb > 1 && nb < n) {
	iws = max(ldwork * nb, (INTEGER) 1);
	if (lwork < iws) {
	    nb = lwork / ldwork;
	    nbmin = max((INTEGER) 2, iMlaenv(2, "Cgetri", " ", n, -1, -1, -1));
	}
    } else {
	iws = n;
    }
//Solve the equation inv(A)*L = inv(U) for inv(A).
    if (nb < nbmin || nb >= n) {
//Use unblocked code.
	for (j = n; j >= 1; j--) {
//Copy current column of L to WORK and replace with zeros.
	    for (i = j + 1; i <= n; i++) {
		work[(i - 1)] = A[(i - 1) + (j - 1) * lda];
		A[(i - 1) + (j - 1) * lda] = Zero;
	    }
//Compute current column of inv(A).
	    if (j < n) {
		Cgemv("No transpose", n, n - j, (COMPLEX) - One, &A[0 + j * lda], lda, &work[j], 1, (COMPLEX) One, &A[0 + (j - 1) * lda], 1);
	    }
	}
    } else {
//Use blocked code.
	nn = ((n - 1) / nb) * nb + 1;
	for (j = nn; j >= 1; j = j - nb) {
	    jb = min(nb, n - j + 1);
//Copy current block column of L to WORK and replace with
//zeros.
	    for (jj = j; jj <= j + jb - 1; jj++) {
		for (i = jj + 1; i <= n; i++) {
		    work[(i + (jj - j) * ldwork) - 1] = A[(i - 1) + (jj - 1) * lda];
		    A[(i - 1) + (jj - 1) * lda] = Zero;
		}
	    }
//Compute current block column of inv(A).
	    if (j + jb <= n) {
		Cgemm("No transpose", "No transpose", n, jb, n - j - jb + 1,
		      (COMPLEX) - One, &A[0 + (j + jb - 1) * lda], lda, &work[j + jb - 1], ldwork, (COMPLEX) One, &A[0 + (j - 1) * lda], lda);
	    }
	    Ctrsm("Right", "Lower", "No transpose", "Unit", n, jb, (COMPLEX) One, &work[j - 1], ldwork, &A[0 + (j - 1) * lda], lda);
	}
    }
//Apply column interchanges.
    for (j = n - 1; j >= 1; j--) {
	jp = ipiv[j - 1];
	if (jp != j) {
	    Cswap(n, &A[0 + (j - 1) * lda], 1, &A[0 + (jp - 1) * lda], 1);
	}
    }
    work[0] = iws;
    return;
}
