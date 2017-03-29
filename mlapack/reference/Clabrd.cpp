/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clabrd.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Clabrd(INTEGER m, INTEGER n, INTEGER nb, COMPLEX * A, INTEGER lda, REAL * d, REAL * e, COMPLEX * tauq, COMPLEX * taup, COMPLEX * x, INTEGER ldx, COMPLEX * y, INTEGER ldy)
{
    INTEGER i;
    COMPLEX alpha;
    REAL One = 1.0, Zero = 0.0;

    if (m <= 0 || n <= 0) {
	return;
    }
    if (m >= n) {
//Reduce to upper bidiagonal form
	for (i = 0; i < nb; i++) {
//Update A(i:m,i)
	    Clacgv(i - 1, &y[i + ldy], ldy);
	    Cgemv("No transpose", m - i + 1, i - 1, (COMPLEX) - One, &A[i + lda], lda, &y[i + ldy], ldy, (COMPLEX) One, &A[i + i * lda], 1);
	    Clacgv(i - 1, &y[i + ldy], ldy);
	    Cgemv("No transpose", m - i + 1, i - 1, (COMPLEX) - One, &x[i + ldx], ldx, &A[i * lda], 1, (COMPLEX) One, &A[i + i * lda], 1);
//Generate reflection Q(i) to annihilate A(i+1:m,i)
	    alpha = A[i + i * lda];
	    Clarfg(m - i + 1, &alpha, &A[min(i + 1, m) + i * lda], 1, &tauq[i]);
	    d[i] = alpha.real();
	    if (i < n) {
		A[i + i * lda] = (COMPLEX) One;
//Compute Y(i+1:n,i)
		Cgemv("Conjugate transpose", m - i + 1, n - i, (COMPLEX) One, &A[i + (i + 1) * lda], lda, &A[i + i * lda], 1, (COMPLEX) Zero, &y[i + 1 + i * ldy], 1);
		Cgemv("Conjugate transpose", m - i + 1, i - 1, (COMPLEX) One, &A[i + lda], lda, &A[i + i * lda], 1, (COMPLEX) Zero, &y[i * ldy + 1], 1);
		Cgemv("No transpose", n - i, i - 1, (COMPLEX) - One, &y[i + 1 + ldy], ldy, &y[i * ldy + 1], 1, (COMPLEX) One, &y[i + 1 + i * ldy], 1);
		Cgemv("Conjugate transpose", m - i + 1, i - 1, (COMPLEX) One, &x[i + ldx], ldx, &A[i + i * lda], 1, (COMPLEX) Zero, &y[i * ldy + 1], 1);
		Cgemv("Conjugate transpose", i - 1, n - i, (COMPLEX) - One, &A[(i + 1) * lda], lda, &y[i * ldy + 1], 1, (COMPLEX) One, &y[i + 1 + i * ldy], 1);
		Cscal(n - i, tauq[i], &y[i + 1 + i * ldy], 1);
//Update A(i,i+1:n)
		Clacgv(n - i, &A[i + (i + 1) * lda], lda);
		Clacgv(i, &A[i + lda], lda);
		Cgemv("No transpose", n - i, i, (COMPLEX) - One, &y[i + 1 + ldy], ldy, &A[i + lda], lda, (COMPLEX) One, &A[i + (i + 1) * lda], lda);
		Clacgv(i, &A[i + lda], lda);
		Clacgv(i - 1, &x[i + ldx], ldx);
		Cgemv("Conjugate transpose", i - 1, n - i, (COMPLEX) - One, &A[(i + 1) * lda], lda, &x[i + ldx], ldx, (COMPLEX) One, &A[i + (i + 1) * lda], lda);
		Clacgv(i - 1, &x[i + ldx], ldx);
//Generate reflection P(i) to annihilate A(i,i+2:n)
		alpha = A[i + (i + 1) * lda];
		Clarfg(n - i, &alpha, &A[i + min(i + 2, n) * lda], lda, &taup[i]);
		e[i] = alpha.real();
		A[i + (i + 1) * lda] = (COMPLEX) One;
//Compute X(i+1:m,i)
		Cgemv("No transpose", m - i, n - i, (COMPLEX) One, &A[i + 1 + (i + 1) * lda], lda, &A[i + (i + 1) * lda], lda, (COMPLEX) Zero, &x[i + 1 + i * ldx], 1);
		Cgemv("Conjugate transpose", n - i, i, (COMPLEX) One, &y[i + 1 + ldy], ldy, &A[i + (i + 1) * lda], lda, (COMPLEX) Zero, &x[i * ldx + 1], 1);
		Cgemv("No transpose", m - i, i, (COMPLEX) - One, &A[i + 1 + lda], lda, &x[i * ldx + 1], 1, (COMPLEX) One, &x[i + 1 + i * ldx], 1);
		Cgemv("No transpose", i - 1, n - i, (COMPLEX) One, &A[(i + 1) * lda], lda, &A[i + (i + 1) * lda], lda, (COMPLEX) Zero, &x[i * ldx + 1], 1);
		Cgemv("No transpose", m - i, i - 1, (COMPLEX) - One, &x[i + 1 + ldx], ldx, &x[i * ldx + 1], 1, (COMPLEX) One, &x[i + 1 + i * ldx], 1);
		Cscal(m - i, taup[i], &x[i + 1 + i * ldx], 1);
		Clacgv(n - i, &A[i + (i + 1) * lda], lda);
	    }
	}
    } else {
//Reduce to lower bidiagonal form
	for (i = 0; i < nb; i++) {
//Update A(i,i:n)
	    Clacgv(n - i + 1, &A[i + i * lda], lda);
	    Clacgv(i - 1, &A[i + lda], lda);
	    Cgemv("No transpose", n - i + 1, i - 1, (COMPLEX) - One, &y[i + ldy], ldy, &A[i + lda], lda, (COMPLEX) One, &A[i + i * lda], lda);
	    Clacgv(i - 1, &A[i + lda], lda);
	    Clacgv(i - 1, &x[i + ldx], ldx);
	    Cgemv("Conjugate transpose", i - 1, n - i + 1, (COMPLEX) - One, &A[i * lda], lda, &x[i + ldx], ldx, (COMPLEX) One, &A[i + i * lda], lda);
	    Clacgv(i - 1, &x[i + ldx], ldx);
//Generate reflection P(i) to annihilate A(i,i+1:n)
	    alpha = A[i + i * lda];
	    Clarfg(n - i + 1, &alpha, &A[i + min(i + 1, n) * lda], lda, &taup[i]);
	    d[i] = alpha.real();
	    if (i < m) {
		A[i + i * lda] = (COMPLEX) One;
//Compute X(i+1:m,i)
		Cgemv("No transpose", m - i, n - i + 1, (COMPLEX) One, &A[i + 1 + i * lda], lda, &A[i + i * lda], lda, (COMPLEX) Zero, &x[i + 1 + i * ldx], 1);
		Cgemv("Conjugate transpose", n - i + 1, i - 1, (COMPLEX) One, &y[i + ldy], ldy, &A[i + i * lda], lda, (COMPLEX) Zero, &x[i * ldx + 1], 1);
		Cgemv("No transpose", m - i, i - 1, (COMPLEX) - One, &A[i + 1 + lda], lda, &x[i * ldx + 1], 1, (COMPLEX) One, &x[i + 1 + i * ldx], 1);
		Cgemv("No transpose", i - 1, n - i + 1, (COMPLEX) One, &A[i * lda], lda, &A[i + i * lda], lda, (COMPLEX) Zero, &x[i * ldx + 1], 1);
		Cgemv("No transpose", m - i, i - 1, (COMPLEX) - One, &x[i + 1 + ldx], ldx, &x[i * ldx + 1], 1, (COMPLEX) One, &x[i + 1 + i * ldx], 1);
		Cscal(m - i, taup[i], &x[i + 1 + i * ldx], 1);
		Clacgv(n - i + 1, &A[i + i * lda], lda);
//Update A(i+1:m,i)
		Clacgv(i - 1, &y[i + ldy], ldy);
		Cgemv("No transpose", m - i, i - 1, (COMPLEX) - One, &A[i + 1 + lda], lda, &y[i + ldy], ldy, (COMPLEX) One, &A[i + 1 + i * lda], 1);
		Clacgv(i - 1, &y[i + ldy], ldy);
		Cgemv("No transpose", m - i, i, (COMPLEX) - One, &x[i + 1 + ldx], ldx, &A[i * lda], 1, (COMPLEX) One, &A[i + 1 + i * lda], 1);
/*              Generate reflection Q(i) to annihilate A(i+2:m,i) */
		alpha = A[i + 1 + i * lda];
		Clarfg(m - i, &alpha, &A[min(i + 2, m) + i * lda], 1, &tauq[i]);
		e[i] = alpha.real();
		A[i + 1 + i * lda] = (COMPLEX) One;
//Compute Y(i+1:n,i)
		Cgemv("Conjugate transpose", m - i, n - i, (COMPLEX) One, &A[i + 1 + (i + 1) * lda], lda, &A[i + 1 + i * lda]
		      , 1, (COMPLEX) Zero, &y[i + 1 + i * ldy], 1);
		Cgemv("Conjugate transpose", m - i, i - 1, (COMPLEX) One, &A[i + 1 + lda], lda, &A[i + 1 + i * lda], 1, (COMPLEX) Zero, &y[i * ldy + 1], 1);
		Cgemv("No transpose", n - i, i - 1, (COMPLEX) - One, &y[i + 1 + ldy], ldy, &y[i * ldy + 1], 1, (COMPLEX) One, &y[i + 1 + i * ldy], 1);
		Cgemv("Conjugate transpose", m - i, i, (COMPLEX) One, &x[i + 1 + ldx], ldx, &A[i + 1 + i * lda], 1, (COMPLEX) Zero, &y[i * ldy + 1], 1);
		Cgemv("Conjugate transpose", i, n - i, (COMPLEX) - One, &A[(i + 1)
									   * lda], lda, &y[i * ldy + 1], 1, (COMPLEX) One, &y[i + 1 + i * ldy], 1);
		Cscal(n - i, tauq[i], &y[i + 1 + i * ldy], 1);
	    } else {
		Clacgv(n - i + 1, &A[i + i * lda], lda);
	    }
	}
    }
    return;
}
