/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgetrf.cpp,v 1.11 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgetrf(INTEGER m, INTEGER n, REAL * A, INTEGER lda, INTEGER * ipiv, INTEGER * info)
{
    INTEGER i, j, jb, nb, iinfo;
    REAL One = 1.0;

//Test the input parameters.
    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Rgetrf", -(*info));
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	return;
    }
//Determine the block size for this environment.
    nb = iMlaenv(1, "Rgetrf", " ", m, n, -1, -1);
    if (nb <= 1 || nb >= min(m, n)) {
//Use unblocked code.
	Rgetf2(m, n, A, lda, ipiv, info);
    } else {
//Use blocked code.
	for (j = 1; j <= min(m, n); j = j + nb) {
	    jb = min(min(m, n) - j + 1, nb);
//Factor diagonal and subdiagonal blocks and test for exact
//singularity.
	    Rgetf2(m - j + 1, jb, &A[(j - 1) + (j - 1) * lda], lda, &ipiv[j - 1], &iinfo);
//Adjust INFO and the pivot indices.
	    if (*info == 0 && iinfo > 0) {
		*info = iinfo + j - 1;
	    }
	    for (i = j; i <= min(m, j + jb - 1); i++) {
		ipiv[i - 1] = j - 1 + ipiv[i - 1];
	    }
//Apply interchanges to columns 1:J-one
	    Rlaswp(j - 1, A, lda, j, j + jb - 1, ipiv, 1);
	    if (j + jb <= n) {
//Apply interchanges to columns J+JB:N.
		Rlaswp(n - j - jb + 1, &A[0 + (j + jb - 1) * lda], lda, j, j + jb - 1, ipiv, 1);
//Compute block row of U.
		Rtrsm("Left", "Lower", "No transpose", "Unit", jb, n - j - jb + 1, One, &A[(j - 1) + (j - 1) * lda], lda, &A[(j - 1) + (j + jb - 1) * lda], lda);
		if (j + jb <= m) {
//Update trailing submatrix.
		    Rgemm("No transpose", "No transpose", m - j - jb + 1,
			  n - j - jb + 1, jb, -One, &A[(j + jb - 1) + (j - 1) * lda], lda, &A[(j - 1) + (j + jb - 1) * lda], lda, One, &A[(j + jb - 1) + (j + jb - 1) * lda], lda);
		}
	    }
	}
    }
    return;
}
