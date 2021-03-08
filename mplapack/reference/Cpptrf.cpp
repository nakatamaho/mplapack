/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cpptrf.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cpptrf(const char *uplo, INTEGER n, COMPLEX * ap, INTEGER * info)
{
    INTEGER j, jc, jj;
    REAL ajj;
    REAL Zero = 0.0, One = 1.0;
    INTEGER upper;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	Mxerbla("Cpptrf", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (upper) {
//Compute the Cholesky factorization A = U'*U.
	jj = 0;
	for (j = 0; j < n; j++) {
	    jc = jj + 1;
	    jj = jj + j;
//Compute elements 1:J-1 of column J.
	    if (j > 1) {
		Ctpsv("Upper", "Conjugate transpose", "Non-unit", j - 1, ap, &ap[jc], 1);
	    }
//Compute U(J,J) and test for non-positive-definiteness.
	    ajj = ap[jj].real() - Cdotc(j - 1, &ap[jc], 1, &ap[jc], 1).real();
	    if (ajj <= Zero) {
		ap[jj] = ajj;
		goto L30;
	    }
	    ap[jj] = sqrt(ajj);
	}
    } else {
//Compute the Cholesky factorization A = L*L'.
	jj = 0;
	for (j = 0; j < n; j++) {
//Compute L(J,J) and test for non-positive-definiteness.
	    ajj = ap[jj].real();
	    if (ajj <= Zero) {
		ap[jj] = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    ap[jj] = ajj;
//Compute elements J+1:N of column J and update the trailing
//submatrix.
	    if (j < n) {
		CRscal(n - j, One / ajj, &ap[jj + 1], 1);
		Chpr("Lower", n - j, -One, &ap[jj + 1], 1, &ap[jj + n - j + 1]);
		jj = jj + n - j + 1;
	    }
	}
    }
    goto L40;
  L30:
    *info = j;
  L40:
    return;
}
