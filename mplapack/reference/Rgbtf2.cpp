/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgbtf2.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgbtf2(INTEGER m, INTEGER n, INTEGER kl, INTEGER ku, REAL * AB, INTEGER ldab, INTEGER * ipiv, INTEGER * info)
{
    INTEGER i, j, km, jp, ju, kv;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1;

    kv = ku + kl;
//Test the input parameters.
    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (kl < 0) {
	*info = -3;
    } else if (ku < 0) {
	*info = -4;
    } else if (ldab < kl + kv + 1) {
	*info = -6;
    }
    if (*info != 0) {
	Mxerbla("Rgbtf2", -(*info));
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0)
	return;

//Gaussian elimination with partial pivoting
//Set fill-in elements in columns KU+2 to KV to zero.
    for (j = ku + 2; j < min(kv, n); j++) {
	for (i = kv - j + 2; i < kl; i++) {
	    AB[i + j * ldab] = Zero;

	}
    }

//JU is the index of the last column affected by the current stage
//of the factorization.
    ju = 1;
    for (j = 0; j < min(m, n); j++) {

//Set fill-in elements in column J+KV to zero.
	if (j + kv <= n) {

	    for (i = 0; i < kl; i++) {
		AB[i + (j + kv) * ldab] = Zero;

	    }
	}
//Find pivot and test for singularity. KM is the number of
//subdiagonal elements in the current column.
	km = min(kl, m - j);
	jp = iRamax(km + 1, &AB[kv + 1 + j * ldab], 1);
	ipiv[j] = jp + j - 1;
	if (AB[kv + jp + j * ldab] != Zero) {
	    ju = max(ju, min(j + ku + jp - 1, n));
//Apply interchange to columns J to JU.
	    if (jp != 1) {
		Rswap(ju - j + 1, &AB[kv + jp + j * ldab], ldab - 1, &AB[kv + 1 + j * ldab], ldab - 1);
	    }
	    if (km > 0) {
//Compute multipliers.
		mtemp1 = One / AB[kv + 1 + j * ldab];
		Rscal(km, mtemp1, &AB[kv + 2 + j * ldab], 1);
//Update trailing submatrix within the band.

		if (ju > j) {
		    Rger(km, ju - j, -One, &AB[kv + 2 + j * ldab], 1, &AB[kv + (j + 1) * ldab], ldab - 1, &AB[kv + 1 + (j + 1) * ldab], ldab - 1);
		}
	    }
	} else {

//If pivot is zero, set INFO to the index of the pivot
//unless a zero pivot has already been found.
	    if (*info == 0) {
		*info = j;
	    }
	}

    }
    return;
}
