/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rorgqr.cpp,v 1.11 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rorgqr(INTEGER m, INTEGER n, INTEGER k, REAL * A, INTEGER lda, REAL * tau, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER nbmin, nx, iws, nb, lwkopt, lquery, ki = 0, kk;
    INTEGER i, j, l, iinfo, ldwork = 0, ib;
    REAL Zero = 0.0, One = 1.0;

//Test the input arguments
    *info = 0;
    nb = iMlaenv(1, "Rorgqr", " ", m, n, k, -1);
    lwkopt = max((INTEGER) 1, n) * nb;
    work[0] = lwkopt;
    lquery = lwork == -1;
    if (m < 0) {
	*info = -1;
    } else if (n < 0 || n > m) {
	*info = -2;
    } else if (k < 0 || k > n) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -5;
    } else if (lwork < max((INTEGER) 1, n) && !lquery) {
	*info = -8;
    }
    if (*info != 0) {
	Mxerbla("Rorgqr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n <= 0) {
	work[0] = One;
	return;
    }
    nbmin = 2;
    nx = 0;
    iws = n;
    if (nb > 1 && nb < k) {
//Determine when to cross over from blocked to unblocked code.
	nx = max((INTEGER) 0, iMlaenv(3, "Rorgqr", " ", m, n, k, -1));
	if (nx < k) {
//Determine if workspace is large enough for blocked code.
	    ldwork = n;
	    iws = ldwork * nb;
	    if (lwork < iws) {
//Not enough workspace to use optimal NB:  reduce NB and
//determine the minimum value of NB.
		nb = lwork / ldwork;
		nbmin = max((INTEGER) 2, iMlaenv(2, "Rorgqr", " ", m, n, k, -1));
	    }
	}
    }
    if (nb >= nbmin && nb < k && nx < k) {
//Use blocked code after the last block.
//The first kk columns are handled by the block method.
	ki = (k - nx - 1) / nb * nb;
	kk = min(k, ki + nb);
//Set A(1:kk,kk+1:n) to zero.
	for (j = kk + 1; j <= n; j++) {
	    for (i = 1; i <= kk; i++) {
		A[(i - 1) + (j - 1) * lda] = Zero;
	    }
	}
    } else {
	kk = 0;
    }
//Use unblocked code for the last or only block.
    if (kk < n) {
	Rorg2r(m - kk, n - kk, k - kk, &A[kk + kk * lda], lda, &tau[kk], &work[0], &iinfo);
    }
    if (kk > 0) {
//Use blocked code
	for (i = ki + 1; i >= 1; i = i - nb) {
	    ib = min(nb, k - i + 1);
	    if (i + ib <= n) {
//Form the triangular factor of the block reflector
//H = H(i) H(i+1) . . . H(i+ib-1)
		Rlarft("Forward", "Columnwise", m - i + 1, ib, &A[(i - 1) + (i - 1) * lda], lda, &tau[i - 1], work, ldwork);
//Apply H to A(i:m,i+ib:n) from the left
		Rlarfb("Left", "No transpose", "Forward", "Columnwise",
		       m - i + 1, n - i - ib + 1, ib, &A[(i - 1) + (i - 1) * lda], lda, work, ldwork, &A[(i - 1) + (i + ib - 1) * lda], lda, &work[ib], ldwork);
	    }
//Apply H to rows i:m of current block
	    Rorg2r(m - i + 1, ib, ib, &A[(i - 1) + (i - 1) * lda], lda, &tau[i - 1], work, &iinfo);
//Set rows 1:i-1 of current block to zero
	    for (j = i; j <= i + ib - 1; j++) {
		for (l = 1; l <= i - 1; l++) {
		    A[(l - 1) + (j - 1) * lda] = Zero;
		}
	    }
	}
    }
    work[0] = iws;
    return;
}
