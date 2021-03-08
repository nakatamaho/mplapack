/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsytrf.cpp,v 1.11 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rsytrf(const char *uplo, INTEGER n, REAL * A, INTEGER lda, INTEGER * ipiv, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER j, k, kb = 0, nb, iws;
    INTEGER nbmin, iinfo;
    INTEGER upper;
    INTEGER ldwork, lwkopt;
    INTEGER lquery;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    lquery = lwork == -1;
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    } else if (lwork < 1 && !lquery) {
	*info = -7;
    }
    if (*info == 0) {
//Determine the block size
	nb = iMlaenv(1, "Rsytrf", uplo, n, -1, -1, -1);
	lwkopt = n * nb;
	work[0] = lwkopt;
    }
    if (*info != 0) {
	Mxerbla("Rsytrf", -(*info));
	return;
    } else if (lquery) {
	return;
    }
    nbmin = 2;
    ldwork = n;
    if (nb > 1 && nb < n) {
	iws = ldwork * nb;
	if (lwork < iws) {
	    nb = max(lwork / ldwork, (INTEGER) 1);
	    nbmin = max((INTEGER) 2, iMlaenv(2, "Rsytrf", uplo, n, -1, -1, -1));
	}
    } else {
	iws = 1;
    }
    if (nb < nbmin) {
	nb = n;
    }
    if (upper) {
//Factorize A as U*D*U' using the upper triangle of A
//K is the main loop index, decreasing from N to 1 in steps of
//KB, where KB is the number of columns factorized by DLASYF;
//KB is either NB or NB-1, or K for the last block
	k = n;
	while (1) {
//If K < 1, exit from loop
	    if (k < 1)
		break;
	    if (k > nb) {
//Factorize columns k-kb+1:k of A and use blocked code to
//update columns 1:k-kb
		Rlasyf(uplo, k, nb, kb, A, lda, ipiv, work, ldwork, &iinfo);
	    } else {
//Use unblocked code to factorize columns 1:k of A
		Rsytf2(uplo, k, A, lda, ipiv, &iinfo);
		kb = k;
	    }

//Set INFO on the first occurrence of a zero pivot
	    if (*info == 0 && iinfo > 0) {
		*info = iinfo;
	    }
//Decrease K and return to the start of the main loop
	    k = k - kb;
	}
    } else {
//Factorize A as L*D*L' using the lower triangle of A
//K is the main loop index, increasing from 1 to N in steps of
//KB, where KB is the number of columns factorized by DLASYF;
//KB is either NB or NB-1, or N-K+1 for the last block 
	k = 0;
	while (1) {
//If K > N, exit from loop
	    if (k > n)
		break;
	    if (k <= n - nb) {
//Factorize columns k:k+kb-1 of A and use blocked code to
// update columns k+kb:n
		Rlasyf(uplo, n - k + 1, nb, kb, &A[k + k * lda], lda, &ipiv[k], &work[0], ldwork, &iinfo);
	    } else {

//Use unblocked code to factorize columns k:n of A
		Rsytf2(uplo, n - k + 1, &A[k + k * lda], lda, &ipiv[k], &iinfo);
		kb = n - k + 1;
	    }
//Set INFO on the first occurrence of a zero pivot
	    if (*info == 0 && iinfo > 0) {
		*info = iinfo + k - 1;
	    }
//Adjust IPIV
	    for (j = k; j < k + kb - 1; j++) {
		if (ipiv[j] > 0) {
		    ipiv[j] = ipiv[j] + k - 1;
		} else {
		    ipiv[j] = ipiv[j] - k + 1;
		}
	    }
//Increase K and return to the start of the main loop
	    k += kb;
	}
    }
    work[0] = (double) lwkopt;	//cast from INTEGER to REAL
    return;
}
