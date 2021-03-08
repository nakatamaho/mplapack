/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rorgrq.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rorgrq(INTEGER m, INTEGER n, INTEGER k, REAL * A, INTEGER lda, REAL * tau, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, j, l, ib, nb = 0, ii, kk, nx, iws, nbmin, iinfo;
    INTEGER ldwork = 0, lwkopt;
    INTEGER lquery;
    REAL Zero = 0.0;

//Test the input arguments

    *info = 0;
    lquery = lwork == -1;
    if (m < 0) {
	*info = -1;
    } else if (n < m) {
	*info = -2;
    } else if (k < 0 || k > m) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -5;
    }

    if (*info == 0) {
	if (m <= 0) {
	    lwkopt = 1;
	} else {
	    nb = iMlaenv(1, "Rorgrq", " ", m, n, k, -1);
	    lwkopt = m * nb;
	}
	work[0] = lwkopt;
	if (lwork < max((INTEGER) 1, m) && !lquery) {
	    *info = -8;
	}
    }

    if (*info != 0) {
	Mxerbla("Rorgrq", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m <= 0) {
	return;
    }

    nbmin = 2;
    nx = 0;
    iws = m;
    if (nb > 1 && nb < k) {
//Determine when to cross over from blocked to unblocked code.
	nx = max((INTEGER) 0, iMlaenv(3, "Rorgrq", " ", m, n, k, -1));
	if (nx < k) {

//Determine if workspace is large enough for blocked code.
	    ldwork = m;
	    iws = ldwork * nb;
	    if (lwork < iws) {

//Not enough workspace to use optimal NB:  reduce NB and
//determine the minimum value of NB.
		nb = lwork / ldwork;
		nbmin = max((INTEGER) 2, iMlaenv(2, "Dorgrq", " ", m, n, k, -1));
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

//Use blocked code after the first block.
//The last kk rows are handled by the block method.
	kk = min(k, (k - nx + nb - 1) / nb * nb);

//Set A(1:m-kk,n-kk+1:n) to zero.
	for (j = n - kk + 1; j <= n; j++) {
	    for (i = 0; i < m - kk; i++) {
		A[i + j * lda] = Zero;

	    }

	}
    } else {
	kk = 0;
    }
//Use unblocked code for the first or only block.
    Rorgr2(m - kk, n - kk, k - kk, &A[0], lda, &tau[1], &work[0], &iinfo);

    if (kk > 0) {

//Use blocked code
	for (i = k - kk + 1; i < k; i += nb) {
	    ib = min(nb, k - i + 1);
	    ii = m - k + i;
	    if (ii > 1) {

//Form the triangular factor of the block reflector
//H = H(i+ib-1) . . . H(i+1) H(i)
		Rlarft("Backward", "Rowwise", n - k + i + ib - 1, ib, &A[ii + lda], lda, &tau[i], &work[0], ldwork);

//Apply H' to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
		Rlarfb("Right", "Transpose", "Backward", "Rowwise", ii - 1, n - k + i + ib - 1, ib, &A[ii + lda], lda, &work[0], ldwork, &A[0], lda, &work[ib + 1], ldwork);
	    }
//Apply H' to columns 1:n-k+i+ib-1 of current block
	    Rorgr2(ib, n - k + i + ib - 1, ib, &A[ii + lda], lda, &tau[i], &work[0], &iinfo);

//Set columns n-k+i+ib:n of current block to zero
	    for (l = n - k + i + ib; l <= n; l++) {
		for (j = ii; j <= ii + ib - 1; j++) {
		    A[j + l * lda] = Zero;
		}
	    }
	}
    }
    work[0] = iws;
    return;
}
