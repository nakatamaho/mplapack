/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgebd2.cpp,v 1.10 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgebd2(INTEGER m, INTEGER n, COMPLEX * A, INTEGER lda, REAL * d, REAL * e, COMPLEX * tauq, COMPLEX * taup, COMPLEX * work, INTEGER * info)
{
    INTEGER i;
    COMPLEX alpha;
    REAL Zero = 0.0, One = 1.0;

    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -4;
    }
    if (*info < 0) {
	Mxerbla("Cgebd", -(*info));
	return;
    }
    if (m >= n) {
//Reduce to upper bidiagonal form
	for (i = 0; i < n; i++) {
//Generate elementary reflector H(i) to annihilate A(i+1:m,i)
	    alpha = A[i + i * lda];
	    Clarfg(m - i + 1, &alpha, &A[min(i + 1, m) + i * lda], 1, &tauq[i]);
	    d[i] = alpha.real();
	    A[i + i * lda] = One;
//Apply H(i)' to A(i:m,i+1:n) from the left
	    if (i < n) {
		Clarf("Left", m - i + 1, n - i, &A[i + i * lda], 1, conj(tauq[i]), &A[i + (i + 1) * lda], lda, &work[0]);
	    }
	    A[i + i * lda] = d[i];
	    if (i < n) {
//Generate elementary reflector G(i) to annihilate
//A(i,i+2:n)
		Clacgv(n - i, &A[i + (i + 1) * lda], lda);
		alpha = A[i + (i + 1) * lda];
		Clarfg(n - i, &alpha, &A[i + min(i + 2, n) * lda], lda, &taup[i]);
		e[i] = alpha.real();
		A[i + (i + 1) * lda] = One;
//Apply G(i) to A(i+1:m,i+1:n) from the right
		Clarf("Right", m - i, n - i, &A[i + (i + 1) * lda], lda, taup[i], &A[i + 1 + (i + 1) * lda], lda, &work[0]);
		Clacgv(n - i, &A[i + (i + 1) * lda], lda);
		A[i + (i + 1) * lda] = e[i];
	    } else {
		taup[i] = Zero;
	    }
	}
    } else {
//Reduce to lower bidiagonal form
	for (i = 0; i < m; i++) {
//Generate elementary reflector G(i) to annihilate A(i,i+1:n)
	    Clacgv(n - i + 1, &A[i + i * lda], lda);
	    alpha = A[i + i * lda];
	    Clarfg(n - i + 1, &alpha, &A[i + min(i + 1, n) * lda], lda, &taup[i]);
	    d[i] = alpha.real();
	    A[i + i * lda] = One;
//Apply G(i) to A(i+1:m,i:n) from the right
	    if (i < m) {
		Clarf("Right", m - i, n - i + 1, &A[i + i * lda], lda, taup[i], &A[i + 1 + i * lda], lda, &work[0]);
	    }
	    Clacgv(n - i + 1, &A[i + i * lda], lda);
	    A[i + i * lda] = d[i];
	    if (i < m) {
//Generate elementary reflector H(i) to annihilate
//A(i+2:m,i)
		alpha = A[i + 1 + i * lda];
		Clarfg(m - i, &alpha, &A[min(i + 2, m) + i * lda], 1, &tauq[i]);
		e[i] = alpha.real();
		A[i + 1 + i * lda] = One;
//Apply H(i)' to A(i+1:m,i+1:n) from the left
		Clarf("Left", m - i, n - i, &A[i + 1 + i * lda], 1, conj(tauq[i]), &A[i + 1 + (i + 1) * lda], lda, &work[0]);
		A[i + 1 + i * lda] = e[i];
	    } else {
		tauq[i] = Zero;
	    }

	}
    }
    return;
}
