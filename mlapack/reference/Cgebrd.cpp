/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgebrd.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgebrd(INTEGER m, INTEGER n, COMPLEX * A, INTEGER lda, REAL * d, REAL * e, COMPLEX * tauq, COMPLEX * taup, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, j, nb, nx;
    REAL ws;
    INTEGER nbmin, iinfo, minmn;
    INTEGER ldwrkx, ldwrky, lwkopt;
    INTEGER lquery;
    REAL One = 1.0;

//Test the input parameters
    *info = 0;
    nb = max((INTEGER) 1, iMlaenv((INTEGER) 1, "Cgebrd", " ", m, n, -1, -1));
    lwkopt = (m + n) * nb;
    work[1] = lwkopt;
    lquery = lwork == -1;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -4;
    } else {
	if (lwork < max(max((INTEGER) 1, m), n) && !lquery) {
	    *info = -10;
	}
    }
    if (*info < 0) {
	Mxerbla("Cgebrd", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    minmn = min(m, n);
    if (minmn == 0) {
	work[1] = 1;
	return;
    }
    ws = max(m, n);
    ldwrkx = m;
    ldwrky = n;
    if (nb > 1 && nb < minmn) {
//Set the crossover poINTEGER NX.
	nx = max(nb, iMlaenv(3, "Cgebrd", " ", m, n, -1, -1));
//Determine when to switch from blocked to unblocked code.
	if (nx < minmn) {
	    ws = ((m + n) * nb);
	    if (double (lwork) < ws) {
//Not enough work space for the optimal NB, consider using
//a smaller block size.
		nbmin = iMlaenv(2, "Cgebrd", " ", m, n, -1, -1);
		if (lwork >= (m + n) * nbmin) {
		    nb = lwork / (m + n);
		} else {
		    nb = 1;
		    nx = minmn;
		}
	    }
	}
    } else {
	nx = minmn;
    }
    for (i = 1; i <= minmn - nx; i += nb) {
//Reduce rows and columns i:i+ib-1 to bidiagonal form and return
//the matrices X and Y which are needed to update the unreduced
//part of the matrix
	Clabrd(m - i + 1, n - i + 1, nb, &A[i + i * lda], lda, &d[i], &e[i], &tauq[i], &taup[i], &work[0], ldwrkx, &work[ldwrkx * nb + 1], ldwrky);

//Update the trailing submatrix A(i+ib:m,i+ib:n), using
//an update of the form  A := A - V*Y' - X*U'
	Cgemm("No transpose", "Conjugate transpose", m - i - nb + 1, n - i - nb + 1, nb, (COMPLEX) - One, &A[i + nb + i * lda], lda,
	      &work[ldwrkx * nb + nb + 1], ldwrky, One, &A[i + nb + (i + nb) * lda], lda);
	Cgemm("No transpose", "No transpose", m - i - nb + 1, n - i - nb + 1, nb, (COMPLEX) - One, &work[nb + 1], ldwrkx,
	      &A[i + (i + nb) * lda], lda, One, &A[i + nb + (i + nb) * lda], lda);
//Copy diagonal and off-diagonal elements of B back into A
	if (m >= n) {
	    for (j = i; j <= i + nb - 1; j++) {
		A[j + j * lda] = d[j];
		A[j + (j + 1) * lda] = e[j];

	    }
	} else {
	    for (j = i; j <= i + nb - 1; j++) {
		A[j + j * lda] = d[j];
		A[j + 1 + j * lda] = e[j];
	    }
	}
    }
//Use unblocked code to reduce the remainder of the matrix
    Cgebd2(m - i + 1, n - i + 1, &A[i + i * lda], lda, &d[i], &e[i], &tauq[i], &taup[i], &work[0], &iinfo);
    work[1] = ws;
    return;
}
