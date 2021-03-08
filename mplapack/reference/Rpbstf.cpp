/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rpbstf.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rpbstf(const char *uplo, INTEGER n, INTEGER kd, REAL * AB, INTEGER ldab, INTEGER * info)
{
    INTEGER j, m, km;
    REAL ajj;
    INTEGER kld;
    INTEGER upper;
    REAL One = 1.0, Zero = 0.0;

//Test the input parameters.
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
	Mxerbla("Rpbstf", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;
    kld = max((INTEGER) 1, ldab - 1);
//Set the splitting point m.
    m = (n + kd) / 2;
    if (upper) {
//Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).
	for (j = n; j >= m + 1; j--) {
//Compute s(j,j) and test for non-positive-definiteness.
	    ajj = AB[kd + 1 + j * ldab];
	    if (ajj <= Zero) {
		goto L50;
	    }
	    ajj = sqrt(ajj);
	    AB[kd + 1 + j * ldab] = ajj;
	    km = min(j - 1, kd);
//Compute elements j-km:j-1 of the j-th column and update the
//the leading submatrix within the band.
	    Rscal(km, One / ajj, &AB[kd + 1 - km + j * ldab], 1);
	    Rsyr("Upper", km, -One, &AB[kd + 1 - km + j * ldab], 1, &AB[kd + 1 + (j - km) * ldab], kld);
	}
//Factorize the updated submatrix A(1:m,1:m) as U**T*U.
	for (j = 0; j < m; j++) {
//Compute s(j,j) and test for non-positive-definiteness.
	    ajj = AB[kd + 1 + j * ldab];
	    if (ajj <= Zero) {
		goto L50;
	    }
	    ajj = sqrt(ajj);
	    AB[kd + 1 + j * ldab] = ajj;
	    km = min(kd, m - j);
//Compute elements j+1:j+km of the j-th row and update the
//trailing submatrix within the band.
	    if (km > 0) {
		Rscal(km, One / ajj, &AB[kd + (j + 1) * ldab], kld);
		Rsyr("Upper", km, -One, &AB[kd + (j + 1) * ldab], kld, &AB[kd + 1 + (j + 1) * ldab], kld);
	    }
	}
    } else {
//Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).
	for (j = n; j >= m + 1; j--) {
//Compute s(j,j) and test for non-positive-definiteness.
	    ajj = AB[j * ldab + 1];
	    if (ajj <= Zero) {
		goto L50;
	    }
	    ajj = sqrt(ajj);
	    AB[j * ldab + 1] = ajj;
	    km = min(j - 1, kd);

//Compute elements j-km:j-1 of the j-th row and update the
//trailing submatrix within the band.
	    Rscal(km, One / ajj, &AB[km + 1 + (j - km) * ldab], kld);
	    Rsyr("Lower", km, -One, &AB[km + 1 + (j - km) * ldab], kld, &AB[(j - km) * ldab + 1], kld);
	}

//Factorize the updated submatrix A(1:m,1:m) as U**T*U.
	for (j = 0; j < m; j++) {
//Compute s(j,j) and test for non-positive-definiteness.
	    ajj = AB[j * ldab + 1];
	    if (ajj <= Zero) {
		goto L50;
	    }
	    ajj = sqrt(ajj);
	    AB[j * ldab + 1] = ajj;
	    km = min(kd, m - j);
//Compute elements j+1:j+km of the j-th column and update the
//trailing submatrix within the band.
	    if (km > 0) {
		Rscal(km, One / ajj, &AB[j * ldab + 2], 1);
		Rsyr("Lower", km, -One, &AB[j * ldab + 2], 1, &AB[(j + 1) * ldab + 1], kld);
	    }

	}
    }
    return;

  L50:
    *info = j;
    return;
}
