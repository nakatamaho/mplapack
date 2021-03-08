/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clatrd.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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
#include <stdio.h>

void Clatrd(const char *uplo, INTEGER n, INTEGER nb, COMPLEX * A, INTEGER lda, REAL * e, COMPLEX * tau, COMPLEX * w, INTEGER ldw)
{
    INTEGER i, iw;
    COMPLEX alpha;
    REAL Zero = 0.0, Half = 0.5, One = 1.0;

//Quick return if possible
    if (n <= 0) {
	return;
    }
    if (Mlsame(uplo, "U")) {
//Reduce last NB columns of upper triangle
	for (i = n; i >= n - nb + 1; i--) {
	    iw = i - n + nb;
	    if (i < n) {
//Update A(1:i,i)
		A[(i - 1) + (i - 1) * lda] = A[(i - 1) + (i - 1) * lda].real();
		Clacgv(n - i, &w[(i - 1) + iw * ldw], ldw);
		Cgemv("No transpose", i, n - i, (COMPLEX) - One, &A[0 + i * lda], lda, &w[(i - 1) + iw * ldw], ldw, (COMPLEX) One, &A[0 + (i - 1) * lda], 1);
		Clacgv(n - i, &w[(i - 1) + iw * ldw], ldw);
		Clacgv(n - i, &A[(i - 1) + i * lda], lda);
		Cgemv("No transpose", i, n - i, (COMPLEX) - One, &w[0 + iw * ldw], ldw, &A[(i - 1) + i * lda], lda, (COMPLEX) One, &A[0 + (i - 1) * lda], 1);
		Clacgv(n - i, &A[(i - 1) + i * lda], lda);
		A[(i - 1) + (i - 1) * lda] = A[(i - 1) + (i - 1) * lda].real();
	    }
	    if (i > 1) {
//Generate elementary reflector H(i) to annihilate
//A(1:i-2,i)
		alpha = A[(i - 2) + (i - 1) * lda];
		Clarfg(i - 1, &alpha, &A[(i - 1) * lda], 1, &tau[i - 2]);
		e[i - 2] = alpha.real();
		A[(i - 2) + (i - 1) * lda] = One;
//Compute W(1:i-1,i)
		Chemv("Upper", i - 1, (COMPLEX) One, &A[0], lda, &A[0 + (i - 1) * lda], 1, (COMPLEX) Zero, &w[0 + (iw - 1) * ldw], 1);
		if (i < n) {
		    Cgemv("Conjugate transpose", i - 1, n - i, (COMPLEX) One, &w[0 + iw * ldw], ldw, &A[0 + (i - 1) * lda], 1, (COMPLEX) Zero, &w[i + (iw - 1) * ldw], 1);
		    Cgemv("No transpose", i - 1, n - i, (COMPLEX) - One, &A[0 + i * lda], lda, &w[i + (iw - 1) * ldw], 1, (COMPLEX) One, &w[0 + (iw - 1) * ldw], 1);
		    Cgemv("Conjugate transpose", i - 1, n - i, (COMPLEX) One, &A[0 + i * lda], lda, &A[0 + (i - 1) * lda], 1, (COMPLEX) Zero, &w[i + (iw - 1) * ldw], 1);
		    Cgemv("No transpose", i - 1, n - i, (COMPLEX) - One, &w[0 + iw * ldw], ldw, &w[i + (iw - 1) * ldw], 1, (COMPLEX) One, &w[0 + (iw - 1) * ldw], 1);
		}
		Cscal(i - 1, tau[i - 2], &w[0 + (iw - 1) * ldw], 1);
		alpha = -Half * tau[i - 2] * Cdotc(i - 1, &w[0 + (iw - 1) * ldw], 1, &A[0 + (i - 1) * lda], 1);
		Caxpy(i - 1, alpha, &A[0 + (i - 1) * lda], 1, &w[0 + (iw - 1) * ldw], 1);
	    }
	}
    } else {
//Reduce first NB columns of lower triangle
	for (i = 1; i <= nb; i++) {
//Update A(i:n,i)
	    A[(i - 1) + (i - 1) * lda] = A[(i - 1) + (i - 1) * lda].real();
	    Clacgv(i - 1, &w[i - 1], ldw);
	    Cgemv("No transpose", n - i + 1, i - 1, (COMPLEX) - One, &A[(i - 1) + 0 * lda], lda, &w[(i - 1) + 0 * ldw], ldw, (COMPLEX) One, &A[(i - 1) + (i - 1) * lda], 1);
	    Clacgv(i - 1, &w[i - 1], ldw);
	    Clacgv(i - 1, &A[i - 1], lda);
	    Cgemv("No transpose", n - i + 1, i - 1, (COMPLEX) - One, &w[(i - 1) + 0 * ldw], ldw, &A[(i - 1) + 0 * lda], lda, (COMPLEX) One, &A[(i - 1) + (i - 1) * lda], 1);
	    Clacgv(i - 1, &A[i - 1], lda);
	    A[(i - 1) + (i - 1) * lda] = A[(i - 1) + (i - 1) * lda].real();
	    if (i < n) {
//Generate elementary reflector H(i) to annihilate
//A(i+2:n,i)
		alpha = A[i + (i - 1) * lda];
		Clarfg(n - i, &alpha, &A[min(i + 2, n) - 1 + (i - 1) * lda], 1, &tau[i - 1]);
		e[i - 1] = alpha.real();
		A[i + (i - 1) * lda] = One;
//Compute W(i+1:n,i)
		Chemv("Lower", n - i, (COMPLEX) One, &A[i + i * lda], lda, &A[i + (i - 1) * lda], 1, (COMPLEX) Zero, &w[i + (i - 1) * ldw], 1);
		Cgemv("Conjugate transpose", n - i, i - 1, (COMPLEX) One, &w[i + 0 * ldw], ldw, &A[i + (i - 1) * lda], 1, (COMPLEX) Zero, &w[0 + (i - 1) * ldw], 1);
		Cgemv("No transpose", n - i, i - 1, (COMPLEX) - One, &A[i + 0 * lda], lda, &w[0 + (i - 1) * ldw], 1, (COMPLEX) One, &w[i + (i - 1) * ldw], 1);
		Cgemv("Conjugate transpose", n - i, i - 1, (COMPLEX) One, &A[i + 0 * lda], lda, &A[i + (i - 1) * lda], 1, (COMPLEX) Zero, &w[0 + (i - 1) * ldw], 1);
		Cgemv("No transpose", n - i, i - 1, -(COMPLEX) One, &w[i + 0 * ldw], ldw, &w[0 + (i - 1) * ldw], 1, (COMPLEX) One, &w[i + (i - 1) * ldw], 1);
		Cscal(n - i, tau[i - 1], &w[i + (i - 1) * ldw], 1);
		alpha = -Half * tau[i - 1] * Cdotc(n - i, &w[i + (i - 1) * ldw], 1, &A[i + (i - 1) * lda], 1);
		Caxpy(n - i, alpha, &A[i + (i - 1) * lda], 1, &w[i + (i - 1) * ldw], 1);
	    }
	}
    }
    return;
}
