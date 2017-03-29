/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlahrd.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlahrd(INTEGER n, INTEGER k, INTEGER nb, REAL * A, INTEGER lda, REAL * tau, REAL * t, INTEGER ldt, REAL * y, INTEGER ldy)
{
    INTEGER i;
    REAL ei = 0;
    REAL One = 1.0, Zero = 0.0;

//Quick return if possible
    if (n <= 1) {
	return;
    }
    for (i = 0; i < nb; i++) {
	if (i > 1) {
//Update A(1:n,i)
//Compute i-th column of A - Y * V'
	    Rgemv("No transpose", n, i - 1, -One, &y[0], ldy, &A[k + i - 1 + lda], lda, One, &A[i * lda], 1);
//Apply I - V * T' * V' to this column (call it b) from the
//left, using the last column of T as workspace
//Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
//         ( V2 )             ( b2 )
//where V1 is unit lower triangular
//w := V1' * b1
	    Rcopy(i - 1, &A[k + 1 + i * lda], 1, &t[nb * ldt + 1], 1);
	    Rtrmv("Lower", "Transpose", "Unit", i - 1, &A[k + 1 + lda], lda, &t[nb * ldt + 1], 1);
//w := w + V2'*b2
	    Rgemv("Transpose", n - k - i + 1, i - 1, One, &A[k + i + lda], lda, &A[k + i + i * lda], 1, One, &t[nb * ldt + 1], 1);
//w := T'*w
	    Rtrmv("Upper", "Transpose", "Non-unit", i - 1, &t[0], ldt, &t[nb * ldt + 1], 1);
//b2 := b2 - V2*w
	    Rgemv("No transpose", n - k - i + 1, i - 1, -One, &A[k + i + lda], lda, &t[nb * ldt + 1], 1, One, &A[k + i + i * lda], 1);
//b1 := b1 - V1*w
	    Rtrmv("Lower", "No transpose", "Unit", i - 1, &A[k + 1 + lda]
		  , lda, &t[nb * ldt + 1], 1);
	    Raxpy(i - 1, -One, &t[nb * ldt + 1], 1, &A[k + 1 + i * lda], 1);
	    A[k + i - 1 + (i - 1) * lda] = ei;
	}
//Generate the elementary reflector H(i) to annihilate
//A(k+i+1:n,i)
	Rlarfg(n - k - i + 1, &A[k + i + i * lda], &A[min(k + i + 1, n) + i * lda], 1, &tau[i]);
	ei = A[k + i + i * lda];
	A[k + i + i * lda] = One;
//Compute  Y(1:n,i)
	Rgemv("No transpose", n, n - k - i + 1, One, &A[(i + 1) * lda], lda, &A[k + i + i * lda], 1, Zero, &y[i * ldy + 1], 12);
	Rgemv("Transpose", n - k - i + 1, i - 1, One, &A[k + i + lda], lda, &A[k + i + i * lda], 1, Zero, &t[i * ldt + 1], 1);
	Rgemv("No transpose", n, i - 1, -One, &y[0], ldy, &t[i * ldt + 1], 1, One, &y[i * ldy + 1], 1);
	Rscal(n, tau[i], &y[i * ldy + 1], 1);
//Compute T(1:i,i)
	Rscal(i - 1, -tau[i], &t[i * ldt + 1], 1);
	Rtrmv("Upper", "No transpose", "Non-unit", i - 1, &t[0], ldt, &t[i * ldt + 1], 1);
	t[i + i * ldt] = tau[i];
    }
    A[k + nb + nb * lda] = ei;
    return;
}
