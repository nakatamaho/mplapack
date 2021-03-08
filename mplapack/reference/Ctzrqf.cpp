/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctzrqf.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctzrqf(INTEGER m, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * tau, INTEGER * info)
{
    INTEGER i, k, m1;
    COMPLEX alpha;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters.
    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if (n < m) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Ctzrqf", -(*info));
	return;
    }
//Perform the factorization.
    if (m == 0) {
	return;
    }
    if (m == n) {
	for (i = 0; i < n; i++) {
	    tau[i] = Zero;
	}
    } else {
	m1 = min(m + 1, n);
	for (k = m; k >= 1; k--) {
//Use a Householder reflection to zero the kth row of A.
//First set up the reflection.
	    A[k + k * lda] = conj(A[k + k * lda]);
	    Clacgv(n - m, &A[k + m1 * lda], lda);
	    alpha = A[k + k * lda];
	    Clarfg(n - m + 1, &alpha, &A[k + m1 * lda], lda, &tau[k]);
	    A[k + k * lda] = alpha;
	    tau[k] = conj(tau[k]);
	    if (tau[k] != Zero && k > 1) {
//We now perform the operation  A := A*P( k )'.
//Use the first ( k - 1 ) elements of TAU to store  a( k ),
//where  a( k ) consists of the first ( k - 1 ) elements of
//the  kth column  of  A.  Also  let  B  denote  the  first
//( k - 1 ) rows of the last ( n - m ) columns of A.
		Ccopy(k - 1, &A[k * lda], 1, &tau[1], 1);
//Form   w = a( k ) + B*z( k )  in TAU.
		Cgemv("No transpose", k - 1, n - m, One, &A[m1 * lda], lda, &A[k + m1 * lda], lda, One, tau, 1);
//Now form  a( k ) := a( k ) - conjg(tau)*w
//and       B      := B      - conjg(tau)*w*z( k )'.
		Caxpy(k - 1, -conj(tau[k]), &tau[1], 1, &A[k * lda], 1);
		Cgerc(k - 1, n - m, conj(tau[k]), &tau[1], 1, &A[k + m1 * lda], lda, &A[m1 * lda], lda);
	    }
	}
    }
    return;
}
