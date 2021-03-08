/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rspgst.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rspgst(INTEGER itype, const char *uplo, INTEGER n, REAL * ap, REAL * bp, INTEGER * info)
{
    INTEGER j, k, j1, k1, jj, kk;
    REAL ct, ajj;
    INTEGER j1j1;
    REAL akk;
    INTEGER k1k1;
    REAL bjj, bkk;
    REAL One = 1.0, Half = 0.5;
    REAL mtemp1;
    INTEGER upper;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (itype < 1 || itype > 3) {
	*info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    }
    if (*info != 0) {
	Mxerbla("Rspgst", -(*info));
	return;
    }

    if (itype == 1) {
	if (upper) {
//Compute inv(U')*A*inv(U)
//J1 and JJ are the indices of A(1,j) and A(j,j)

	    jj = 0;
	    for (j = 0; j < n; j++) {
		j1 = jj + 1;
		jj += j;
//Compute the j-th column of the upper triangle of A

		bjj = bp[jj];
		Rtpsv(uplo, "Transpose", "Nonunit", j, &bp[0], &ap[j1], 1);
		Rspmv(uplo, j - 1, -One, &ap[0], &bp[j1], 1, One, &ap[j1], 1);
		mtemp1 = One / bjj;
		Rscal(j - 1, mtemp1, &ap[j1], 1);
		ap[jj] = (ap[jj] - Rdot(j - 1, &ap[j1], 1, &bp[j1], 1)) / bjj;
	    }
	} else {
//Compute inv(L)*A*inv(L')
//KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)

	    kk = 0;
	    for (k = 0; k < n; k++) {
		k1k1 = kk + n - k + 1;
//Update the lower triangle of A(k:n,k:n)
		akk = ap[kk];
		bkk = bp[kk];
		akk /= bkk * bkk;
		ap[kk] = akk;
		if (k < n) {
		    mtemp1 = One / bkk;
		    Rscal(n - k, mtemp1, &ap[kk + 1], 1);
		    ct = akk * -Half;
		    Raxpy(n - k, ct, &bp[kk + 1], 1, &ap[kk + 1], 1);
		    Rspr2(uplo, n - k, -One, &ap[kk + 1], 1, &bp[kk + 1], 1, &ap[k1k1]);
		    Raxpy(n - k, ct, &bp[kk + 1], 1, &ap[kk + 1], 1);
		    Rtpsv(uplo, "No transpose", "Non-unit", n - k, &bp[k1k1], &ap[kk + 1], 1);
		}
		kk = k1k1;

	    }
	}
    } else {
	if (upper) {
//Compute U*A*U'
//K1 and KK are the indices of A(1,k) and A(k,k)
	    kk = 0;
	    for (k = 0; k < n; k++) {
		k1 = kk + 1;
		kk += k;
//Update the upper triangle of A(1:k,1:k)
		akk = ap[kk];
		bkk = bp[kk];

		Rtpmv(uplo, "No transpose", "Non-unit", k - 1, &bp[0], &ap[k1], 12);
		ct = akk * Half;

		Raxpy(k - 1, ct, &bp[k1], 1, &ap[k1], 1);
		Rspr2(uplo, k - 1, One, &ap[k1], 1, &bp[k1], 1, &ap[0]);
		Raxpy(k - 1, ct, &bp[k1], 1, &ap[k1], 1);
		Rscal(k - 1, bkk, &ap[k1], 1);
		ap[kk] = akk * (bkk * bkk);

	    }
	} else {
//Compute L'*A*L
//JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)

	    jj = 0;
	    for (j = 0; j < n; j++) {
		j1j1 = jj + n - j + 1;
//Compute the j-th column of the lower triangle of A
		ajj = ap[jj];
		bjj = bp[jj];

		ap[jj] = ajj * bjj + Rdot(n - j, &ap[jj + 1], 1, &bp[jj + 1], 1);
		Rscal(n - j, bjj, &ap[jj + 1], 1);
		Rspmv(uplo, n - j, One, &ap[j1j1], &bp[jj + 1], 1, One, &ap[jj + 1], 1);
		Rtpmv(uplo, "Transpose", "Non-unit", n - j + 1, &bp[jj], &ap[jj], 1);
		jj = j1j1;
	    }
	}
    }
    return;
}
