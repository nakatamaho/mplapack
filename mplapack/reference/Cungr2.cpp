/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cungr2.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cungr2(INTEGER m, INTEGER n, INTEGER k, COMPLEX * A, INTEGER lda, COMPLEX * tau, COMPLEX * work, INTEGER * info)
{
    INTEGER i, j, l, ii;
    REAL One = 1.0, Zero = 0.0;

//Test the input arguments
    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if (n < m) {
	*info = -2;
    } else if (k < 0 || k > m) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -5;
    }
    if (*info != 0) {
	Mxerbla("Cungr2", -(*info));
	return;
    }
//Quick return if possible
    if (m <= 0)
	return;

    if (k < m) {
//Initialise rows 1:m-k to rows of the unit matrix
	for (j = 0; j < n; j++) {
	    for (l = 0; l < m - k; l++) {
		A[l + j * lda] = Zero;
	    }
	    if (j > n - m && j <= n - k) {
		A[m - n + j + j * lda] = One;
	    }

	}
    }
    for (i = 0; i < k; i++) {
	ii = m - k + i;
//Apply H(i)' to A(1:m-k+i,1:n-k+i) from the right
	Clacgv(n - m + ii - 1, &A[ii + lda], lda);
	A[ii + (n - m + ii) * lda] = One;
	Clarf("Right", ii - 1, n - m + ii, &A[ii + lda], lda, conj(tau[i]), &A[0], lda, &work[0]);

	Cscal(n - m + ii - 1, -tau[i], &A[ii + lda], lda);
	Clacgv(n - m + ii - 1, &A[ii + lda], lda);

	A[ii + (n - m + ii) * lda] = One - conj(tau[i]);
//Set A(m-k+i,n-k+i+1:n) to zero
	for (l = n - m + ii + 1; l <= n; l++) {
	    A[ii + l * lda] = Zero;
	}
    }
    return;
}
