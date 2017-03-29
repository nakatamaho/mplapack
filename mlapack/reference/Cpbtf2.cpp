/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cpbtf2.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cpbtf2(const char *uplo, INTEGER n, INTEGER kd, COMPLEX * AB, INTEGER ldab, INTEGER * info)
{
    INTEGER j, kn;
    REAL ajj;
    INTEGER kld;
    INTEGER upper;
    REAL Zero = 0.0, One = 1.0;

    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (kd < 0) {
	*info = -3;
    } else if (ldab < kd + 1) {
	*info = -5;
    }
    if (*info != 0) {
	Mxerbla("Cpbtf2", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    kld = max((INTEGER) 1, ldab - 1);
    if (upper) {
//Compute the Cholesky factorization A = U'*U.
	for (j = 0; j < n; j++) {
//Compute U(J,J) and test for non-positive-definiteness.
	    ajj = AB[kd + 1 + j * ldab].real();
	    if (ajj <= Zero) {
		AB[kd + 1 + j * ldab] = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    AB[kd + 1 + j * ldab] = ajj;
//Compute elements J+1:J+KN of row J and update the
//trailing submatrix within the band.
	    kn = min(kd, n - j);
	    if (kn > 0) {
		CRscal(kn, One / ajj, &AB[kd + (j + 1) * ldab], kld);
		Clacgv(kn, &AB[kd + (j + 1) * ldab], kld);
		Cher("Upper", kn, -One, &AB[kd + (j + 1) * ldab], kld, &AB[kd + 1 + (j + 1) * ldab], kld);
		Clacgv(kn, &AB[kd + (j + 1) * ldab], kld);
	    }
	}
    } else {
//Compute the Cholesky factorization A = L*L'.
	for (j = 0; j < n; j++) {
//Compute L(J,J) and test for non-positive-definiteness.
	    ajj = AB[j * ldab + 1].real();
	    if (ajj <= Zero) {
		AB[j * ldab + 1] = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    AB[j * ldab + 1] = ajj;
//Compute elements J+1:J+KN of column J and update the
//trailing submatrix within the band.
	    kn = min(kd, n - j);
	    if (kn > 0) {
		CRscal(kn, One / ajj, &AB[j * ldab + 2], 1);
		Cher("Lower", kn, -One, &AB[j * ldab + 2], 1, &AB[(j + 1) * ldab + 1], kld);
	    }
	}
    }
    return;
  L30:
    *info = j;
    return;
}
