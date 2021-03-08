/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chetd2.cpp,v 1.11 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Chetd2(const char *uplo, INTEGER n, COMPLEX * A, INTEGER lda, REAL * d, REAL * e, COMPLEX * tau, INTEGER * info)
{
    INTEGER upper;
    INTEGER i;
    COMPLEX taui, alpha;
    REAL Zero = 0.0, Half = 0.5, One = 1.0;

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
	Mxerbla("Chetd2", -(*info));
	return;
    }
//Quick return if possible
    if (n <= 0) {
	return;
    }
    if (upper) {
//Reduce the upper triangle of A
	A[(n - 1) + (n - 1) * lda] = A[(n - 1) + (n - 1) * lda].real();
	for (i = n - 1; i >= 1; i--) {
//Generate elementary reflector H(i) = I - tau * v * v'
//to annihilate A(1:i-1,i+1)
	    alpha = A[(i - 1) + i * lda];
	    Clarfg(i, &alpha, &A[0 + i * lda], 1, &taui);
	    e[i - 1] = alpha.real();
	    if (taui != Zero) {
//Apply H(i) from both sides to A(1:i,1:i)
		A[(i - 1) + i * lda] = One;
//Compute x := tau * A * v  storing x in TAU(1:i)
		Chemv(uplo, i, taui, A, lda, &A[0 + i * lda], 1, Zero, tau, 1);
//Compute w := x - 1/2 * tau * (x'*v) * v
		alpha = -Half * taui * Cdotc(i, tau, 1, &A[0 + i * lda], 1);
		Caxpy(i, alpha, &A[0 + i * lda], 1, tau, 1);
//Apply the transformation as a rank-2 update
//A := A - v * w' - w * v'
		Cher2(uplo, i, -(COMPLEX) One, &A[0 + i * lda], 1, tau, 1, A, lda);
	    } else {
		A[(i - 1) + (i - 1) * lda] = A[(i - 1) + (i - 1) * lda].real();
	    }
	    A[(i - 1) + i * lda] = e[i - 1];
	    d[i] = A[i + i * lda].real();
	    tau[i - 1] = taui;
	}
	d[0] = A[0].real();
    } else {
//Reduce the lower triangle of A
	A[0] = A[0].real();
	for (i = 1; i <= n - 1; i++) {
//Generate elementary reflector H(i) = I - tau * v * v'
//to annihilate A(i+2:n,i)
	    alpha = A[i + (i - 1) * lda];
	    Clarfg(n - i, &alpha, &A[min(i + 2, n) - 1 + (i - 1) * lda], 1, &taui);
	    e[i - 1] = alpha.real();
	    if (taui != Zero) {
//Apply H(i) from both sides to A(i+1:n,i+1:n)
		A[i + (i - 1) * lda] = One;
//Compute x := tau * A * v  storing y in TAU(i:n-1)
		Chemv(uplo, n - i, taui, &A[i + i * lda], lda, &A[i + (i - 1) * lda], 1, Zero, &tau[i - 1], 1);
//Compute w := x - 1/2 * tau * (x'*v) * v
		alpha = -Half * taui * Cdotc(n - i, &tau[i - 1], 1, &A[i + (i - 1) * lda], 1);
		Caxpy(n - i, alpha, &A[i + (i - 1) * lda], 1, &tau[i - 1], 1);
//Apply the transformation as a rank-2 update:
//A := A - v * w' - w * v'
		Cher2(uplo, n - i, -(COMPLEX) One, &A[i + (i - 1) * lda], 1, &tau[i - 1], 1, &A[i + i * lda], lda);
	    } else {
		A[i + i * lda] = A[i + i * lda].real();
	    }
	    A[i + (i - 1) * lda] = e[i - 1];
	    d[i - 1] = A[(i - 1) + (i - 1) * lda].real();
	    tau[i - 1] = taui;
	}
	d[n - 1] = A[(n - 1) + (n - 1) * lda].real();
    }
    return;
}
