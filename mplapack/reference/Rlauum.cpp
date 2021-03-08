/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlauum.cpp,v 1.9 2010/08/18 08:50:27 nakatamaho Exp $ 
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

void Rlauum(const char *uplo, INTEGER n, REAL * A, INTEGER lda, INTEGER * info)
{
    INTEGER i, ib, nb;
    INTEGER upper;
    REAL One = 1.0;

    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Rlauum", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;
//Determine the block size for this environment.
    nb = iMlaenv(1, "Rlauum", uplo, n, -1, -1, -1);
    if (nb <= 1 || nb >= n) {
//Use unblocked code
	Rlauu2(uplo, n, A, lda, info);
    } else {
//Use blocked code
	if (upper) {
//Compute the product U * U'.
	    for (i = 1; i <= n; i = i + nb) {
		ib = min(nb, n - i + 1);
		Rtrmm("Right", "Upper", "Transpose", "Non-unit", i - 1, ib, One, &A[(i - 1) + (i - 1) * lda], lda, &A[(i - 1) * lda], lda);
		Rlauu2("Upper", ib, &A[(i - 1) + (i - 1) * lda], lda, info);
		if (i + ib <= n) {
		    Rgemm("No transpose", "Transpose", i - 1, ib, n - i - ib + 1, One, &A[(i + ib - 1) * lda], lda, &A[(i - 1) + (i + ib - 1) * lda], lda, One, &A[(i - 1) * lda], lda);
		    Rsyrk("Upper", "No transpose", ib, n - i - ib + 1, One, &A[(i - 1) + (i + ib - 1) * lda], lda, One, &A[(i - 1) + (i - 1) * lda], lda);
		}
	    }
	} else {
//Compute the product L' * L.
	    for (i = 1; i <= n; i = i + nb) {
		ib = min(nb, n - i + 1);
		Rtrmm("Left", "Lower", "Transpose", "Non-unit", ib, i - 1, One, &A[(i - 1) + (i - 1) * lda], lda, &A[(i - 1)], lda); 
		Rlauu2("Lower", ib, &A[(i - 1) + (i - 1) * lda], lda, info);
		if (i + ib <= n) {
		    Rgemm("Transpose", "No transpose", ib, i - 1, n - i - ib + 1, One, &A[(i + ib - 1) + (i - 1) * lda], lda, &A[(i + ib - 1)], lda, One, &A[(i - 1)], lda);
		    Rsyrk("Lower", "Transpose", ib, n - i - ib + 1, One, &A[(i + ib - 1) + (i - 1 ) * lda], lda, One, &A[(i - 1 ) + (i - 1) * lda], lda);
		}
	    }
	}
    }
    return;
}
