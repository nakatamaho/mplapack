/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlabrd.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlabrd(INTEGER m, INTEGER n, INTEGER nb, REAL * A, INTEGER lda, REAL * d, REAL * e, REAL * tauq, REAL * taup, REAL * x, INTEGER ldx, REAL * y, INTEGER ldy)
{
    INTEGER i;
    REAL Zero = 0.0, One = 1.0;

//Quick return if possible
    if (m <= 0 || n <= 0) {
	return;
    }
    if (m >= n) {
//Reduce to upper bidiagonal form
	for (i = 0; i < nb; i++) {
//Update A(i:m,i)
	    Rgemv("No transpose", m - i + 1, i - 1, -One, &A[i + lda], lda, &y[i + ldy], ldy, One, &A[i + i * lda], 12);
	    Rgemv("No transpose", m - i + 1, i - 1, -One, &x[i + ldx], ldx, &A[i * lda], 1, One, &A[i + i * lda], 12);
//Generate reflection Q(i) to annihilate A(i+1:m,i)
	    Rlarfg(m - i + 1, &A[i + i * lda], &A[min(i + 1, m) + i * lda], 1, &tauq[i]);
	    d[i] = A[i + i * lda];
	    if (i < n) {
		A[i + i * lda] = One;
//Compute Y(i+1:n,i)
		Rgemv("Transpose", m - i + 1, n - i, One, &A[i + (i + 1) * lda], lda, &A[i + i * lda], 1, Zero, &y[i + 1 + i * ldy], 1);
		Rgemv("Transpose", m - i + 1, i - 1, One, &A[i + lda], lda, &A[i + i * lda], 1, Zero, &y[i * ldy + 1], 1);
		Rgemv("No transpose", n - i, i - 1, -One, &y[i + 1 + ldy], ldy, &y[i * ldy + 1], 1, One, &y[i + 1 + i * ldy], 12);
		Rgemv("Transpose", m - i + 1, i - 1, One, &x[i + ldx], ldx, &A[i + i * lda], 1, Zero, &y[i * ldy + 1], 1);
		Rgemv("Transpose", i - 1, n - i, -One, &A[(i + 1) * lda], lda, &y[i * ldy + 1], 1, One, &y[i + 1 + i * ldy], 1);
		Rscal(n - i, tauq[i], &y[i + 1 + i * ldy], 1);
//Update A(i,i+1:n)
		Rgemv("No transpose", n - i, i, -One, &y[i + 1 + ldy], ldy, &A[i + lda], lda, One, &A[i + (i + 1) * lda], lda);
		Rgemv("Transpose", i - 1, n - i, -One, &A[(i + 1) * lda], lda, &x[i + ldx], ldx, One, &A[i + (i + 1) * lda], lda);
//Generate reflection P(i) to annihilate A(i,i+2:n)
		Rlarfg(n - i, &A[i + (i + 1) * lda], &A[i + min(i + 2, n) * lda], lda, &taup[i]);
		e[i] = A[i + (i + 1) * lda];
		A[i + (i + 1) * lda] = One;
//Compute X(i+1:m,i)
		Rgemv("No transpose", m - i, n - i, One, &A[i + 1 + (i + 1) * lda], lda, &A[i + (i + 1) * lda], lda, Zero, &x[i + 1 + i * ldx], 12);
		Rgemv("Transpose", n - i, i, One, &y[i + 1 + ldy], ldy, &A[i + (i + 1) * lda], lda, Zero, &x[i * ldx + 1], 1);
		Rgemv("No transpose", m - i, i, -One, &A[i + 1 + lda], lda, &x[i * ldx + 1], 1, One, &x[i + 1 + i * ldx], 12);
		Rgemv("No transpose", i - 1, n - i, One, &A[(i + 1) * lda], lda, &A[i + (i + 1) * lda], lda, Zero, &x[i * ldx + 1], 12);
		Rgemv("No transpose", m - i, i - 1, -One, &x[i + 1 + ldx], ldx, &x[i * ldx + 1], 1, One, &x[i + 1 + i * ldx], 12);
		Rscal(m - i, taup[i], &x[i + 1 + i * ldx], 1);
	    }
	}
    } else {
//Reduce to lower bidiagonal form
	for (i = 0; i < nb; i++) {
//Update A(i,i:n)
	    Rgemv("No transpose", n - i + 1, i - 1, -One, &y[i + ldy], ldy, &A[i + lda], lda, One, &A[i + i * lda], lda);
	    Rgemv("Transpose", i - 1, n - i + 1, -One, &A[i * lda], lda, &x[i + ldx], ldx, One, &A[i + i * lda], lda);
//Generate reflection P(i) to annihilate A(i,i+1:n)
	    Rlarfg(n - i + 1, &A[i + i * lda], &A[i + min(i + 1, n) * lda], lda, &taup[i]);
	    d[i] = A[i + i * lda];
	    if (i < m) {
		A[i + i * lda] = One;
//Compute X(i+1:m,i)
		Rgemv("No transpose", m - i, n - i + 1, One, &A[i + 1 + i * lda], lda, &A[i + i * lda], lda, Zero, &x[i + 1 + i * ldx], 1);
		Rgemv("Transpose", n - i + 1, i - 1, One, &y[i + ldy], ldy, &A[i + i * lda], lda, Zero, &x[i * ldx + 1], 1);
		Rgemv("No transpose", m - i, i - 1, -One, &A[i + 1 + lda], lda, &x[i * ldx + 1], 1, One, &x[i + 1 + i * ldx], 1);
		Rgemv("No transpose", i - 1, n - i + 1, One, &A[i * lda], lda, &A[i + i * lda], lda, Zero, &x[i * ldx + 1], 1);
		Rgemv("No transpose", m - i, i - 1, -One, &x[i + 1 + ldx], ldx, &x[i * ldx + 1], 1, One, &x[i + 1 + i * ldx], 1);
		Rscal(m - i, taup[i], &x[i + 1 + i * ldx], 1);
//Update A(i+1:m,i)
		Rgemv("No transpose", m - i, i - 1, -One, &A[i + 1 + lda], lda, &y[i + ldy], ldy, One, &A[i + 1 + i * lda], 1);
		Rgemv("No transpose", m - i, i, -One, &x[i + 1 + ldx], ldx, &A[i * lda], 1, One, &A[i + 1 + i * lda], 1);
//Generate reflection Q(i) to annihilate A(i+2:m,i)
		Rlarfg(m - i, &A[i + 1 + i * lda], &A[min(i + 2, m) + i * lda], 1, &tauq[i]);
		e[i] = A[i + 1 + i * lda];
		A[i + 1 + i * lda] = One;
//Compute Y(i+1:n,i)
		Rgemv("Transpose", m - i, n - i, One, &A[i + 1 + (i + 1) * lda], lda, &A[i + 1 + i * lda], 1, Zero, &y[i + 1 + i * ldy], 1);
		Rgemv("Transpose", m - i, i - 1, One, &A[i + 1 + lda], lda, &A[i + 1 + i * lda], 1, Zero, &y[i * ldy + 1], 1);
		Rgemv("No transpose", n - i, i - 1, -One, &y[i + 1 + ldy], ldy, &y[i * ldy + 1], 1, One, &y[i + 1 + i * ldy], 1);
		Rgemv("Transpose", m - i, i, One, &x[i + 1 + ldx], ldx, &A[i + 1 + i * lda], 1, Zero, &y[i * ldy + 1], 1);
		Rgemv("Transpose", i, n - i, -One, &A[(i + 1) * lda], lda, &y[i * ldy + 1], 1, One, &y[i + 1 + i * ldy], 1);
		Rscal(n - i, tauq[i], &y[i + 1 + i * ldy], 1);
	    }
	}
    }
    return;
}
