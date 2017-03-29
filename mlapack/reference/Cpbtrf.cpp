/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cpbtrf.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cpbtrf(const char *uplo, INTEGER n, INTEGER kd, COMPLEX * AB, INTEGER ldab, INTEGER * info)
{
    INTEGER i, j, i2, i3, ib, nb, ii, jj;
    COMPLEX work[1056];
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters.
    *info = 0;
    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (kd < 0) {
	*info = -3;
    } else if (ldab < kd + 1) {
	*info = -5;
    }
    if (*info != 0) {
	Mxerbla("Cpbtrf", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Determine the block size for this environment
    nb = iMlaenv(1, "Cpbtrf", uplo, n, kd, -1, -1);
//The block size must not exceed the semi-bandwidth KD, and must not
//exceed the limit set by the size of the local array WORK.
    nb = min(nb, (INTEGER) 32);
    if (nb <= 1 || nb > kd) {
//Use unblocked code
	Cpbtf2(uplo, n, kd, AB, ldab, info);
    } else {
//Use blocked code
	if (Mlsame(uplo, "U")) {
//Compute the Cholesky factorization of a Hermitian band
//matrix, given the upper triangle of the matrix in band
//storage.
//Zero the upper triangle of the work array.
	    for (j = 0; j < nb; j++) {
		for (i = 0; i < j - 1; i++) {
		    work[i + j * 33 - 34] = Zero;
		}
	    }
//Process the band matrix one diagonal block at a time.
	    for (i = 1; i <= n; i = i + nb) {
		ib = min(nb, n - i + 1);
//Factorize the diagonal block
		Cpotf2(uplo, ib, &AB[kd + 1 + i * ldab], ldab - 1, &ii);
		if (ii != 0) {
		    *info = i + ii - 1;
		    goto L150;
		}
		if (i + ib <= n) {
//Update the relevant part of the trailing submatrix.
//If A11 denotes the diagonal block which has just been
//factorized, then we need to update the remaining
//blocks in the diagram:
//   A11   A12   A13
//         A22   A23
//               A33
//The numbers of rows and columns in the partitioning
//are IB, I2, I3 respectively. The blocks A12, A22 and
//A23 are empty if IB = KD. The upper triangle of A13
//lies outside the band.
		    i2 = min(kd - ib, n - i - ib + 1);
		    i3 = min(ib, n - i - kd + 1);
		    if (i2 > 0) {
//Update A12
			Ctrsm("Left", "Upper", "Conjugate transpose", "Non-" "unit", ib, i2, One, &AB[kd + 1 + i * ldab], ldab - 1, &AB[kd + 1 - ib + (i + ib)
																	* ldab], ldab - 1);
//Update A22
			Cherk("Upper", "Conjugate transpose", i2, ib, -One, &AB[kd + 1 - ib + (i + ib) * ldab], ldab - 1, One, &AB[kd + 1 + (i + ib) * ldab], ldab - 1);
		    }
		    if (i3 > 0) {
//Copy the lower triangle of A13 into the work array.
			for (jj = 0; jj <= i3; jj++) {
			    for (ii = jj; ii <= ib; ii++) {
				work[ii + jj * 33 - 34] = AB[ii - jj + 1 + (jj + i + kd - 1) * ldab];
			    }
			}
//Update A13 (in the work array).
			Ctrsm("Left", "Upper", "Conjugate transpose", "Non-" "unit", ib, ldab - 1, One, &AB[kd + 1 + i * ldab], ldab - 1, work, 33);
//Update A23
			if (i2 > 0) {
			    Cgemm("Conjugate transpose", "No transpose", i2,
				  ldab - 1, ib, (COMPLEX) - One, &AB[kd + 1 - ib + (i + ib) * ldab], ldab - 1, work, 33, (COMPLEX) One, &AB[ib + 1 + (i + kd) * ldab], ldab - 1);
			}
//Update A33
			Cherk("Upper", "Conjugate transpose", ldab - 1, ib, -One, work, 33, One, &AB[kd + 1 + (i + kd) * ldab], ldab - 1);
//Copy the lower triangle of A13 back into place.
			for (jj = 0; jj <= i3; jj++) {
			    for (ii = jj; ii <= ib; ii++) {
				AB[ii - jj + 1 + (jj + i + kd - 1) * ldab] = work[ii + jj * 33 - 34];
			    }
			}
		    }
		}
	    }
	} else {
//Compute the Cholesky factorization of a Hermitian band
//matrix, given the lower triangle of the matrix in band
//storage.
//Zero the lower triangle of the work array.
	    for (j = 0; j < nb; j++) {
		for (i = j + 1; i <= nb; i++) {
		    work[i + j * 33 - 34] = Zero;
		}
	    }
//Process the band matrix one diagonal block at a time.
	    for (i = 1; i <= n; i = i + nb) {
		ib = min(nb, n - i + 1);
//Factorize the diagonal block
		Cpotf2(uplo, ib, &AB[i * ldab + 1], ldab - 1, &ii);
		if (ii != 0) {
		    *info = i + ii - 1;
		    goto L150;
		}
		if (i + ib <= n) {
//Update the relevant part of the trailing submatrix.
//If A11 denotes the diagonal block which has just been
//factorized, then we need to update the remaining
//blocks in the diagram:
//   A11
//   A21   A22
//   A31   A32   A33
//The numbers of rows and columns in the partitioning
//are IB, I2, I3 respectively. The blocks A21, A22 and
//A32 are empty if IB = KD. The lower triangle of A31
//lies outside the band.
		    i2 = min(kd - ib, n - i - ib + 1);
		    i3 = min(ib, n - i - kd + 1);
		    if (i2 > 0) {
//Update A21
			Ctrsm("Right", "Lower", "Conjugate transpose", "Non" "-unit", i2, ib, One, &AB[i * ldab + 1], ldab - 1, &AB[ib + 1 + i * ldab], ldab - 1);
//Update A22
			Cherk("Lower", "No transpose", i2, ib, -One, &AB[ib + 1 + i * ldab], ldab - 1, One, &AB[(i + ib) * ldab + 1], ldab - 1);
		    }
		    if (i3 > 0) {
//Copy the upper triangle of A31 into the work array.
			for (jj = 0; jj <= ib; jj++) {
			    for (ii = 0; ii <= min(jj, i3); ii++) {
				work[ii + jj * 33 - 34] = AB[kd + 1 - jj + ii + (jj + i - 1) * ldab];
			    }
			}
//Update A31 (in the work array).
			Ctrsm("Right", "Lower", "Conjugate transpose", "Non" "-unit", ldab - 1, ib, One, &AB[i * ldab + 1], ldab - 1, work, 33);
//Update A32
			if (i2 > 0) {
			    Cgemm("No transpose", "Conjugate transpose",
				  ldab - 1, i2, ib, (COMPLEX) - One, work, 33, &AB[ib + 1 + i * ldab], ldab - 1, (COMPLEX) One, &AB[kd + 1 - ib + (i + ib) * ldab], ldab - 1);
			}
//Update A33
			Cherk("Lower", "No transpose", ldab - 1, ib, -One, work, 33, One, &AB[(i + kd) * ldab + 1], ldab - 1);
//Copy the upper triangle of A31 back into place.
			for (jj = 0; jj <= ib; jj++) {
			    for (ii = 0; ii <= min(jj, i3); ii++) {
				AB[kd + 1 - jj + ii + (jj + i - 1) * ldab] = work[ii + jj * 33 - 34];
			    }
			}
		    }
		}
	    }
	}
    }
    return;
  L150:
    return;
}
