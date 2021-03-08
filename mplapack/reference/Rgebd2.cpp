/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgebd2.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgebd2(INTEGER m, INTEGER n, REAL * A, INTEGER lda, REAL * d, REAL * e, REAL * tauq, REAL * taup, REAL * work, INTEGER * info)
{
    INTEGER i;
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
	Mxerbla("Rgebd2", -(*info));
	return;
    }
    if (m >= n) {
//Reduce to upper bidiagonal form
	for (i = 0; i < n; i++) {
//Generate elementary reflector H(i) to annihilate A(i+1:m,i)
	    Rlarfg(m - i + 1, &A[i + i * lda], &A[min(i + 1, m) + i * lda], 1, &tauq[i]);
	    d[i] = A[i + i * lda];
	    A[i + i * lda] = One;
//Apply H(i) to A(i:m,i+1:n) from the left
	    if (i < n) {
		Rlarf("Left", m - i + 1, n - i, &A[i + i * lda], 1, tauq[0], &A[i + (i + 1) * lda], lda, &work[0]);
	    }
	    A[i + i * lda] = d[i];
	    if (i < n) {
//Generate elementary reflector G(i) to annihilate
//A(i,i+2:n)
		Rlarfg(n - i, &A[i + (i + 1) * lda], &A[i + min(i + 2, n) * lda], lda, &taup[i]);
		e[i] = A[i + (i + 1) * lda];
		A[i + (i + 1) * lda] = One;
//Apply G(i) to A(i+1:m,i+1:n) from the right
		Rlarf("Right", m - i, n - i, &A[i + (i + 1) * lda], lda, taup[i], &A[i + 1 + (i + 1) * lda], lda, &work[0]);
		A[i + (i + 1) * lda] = e[i];
	    } else {
		taup[i] = Zero;
	    }
	}
    } else {
//Reduce to lower bidiagonal form
	for (i = 0; i < m; i++) {
//Generate elementary reflector G(i) to annihilate A(i,i+1:n)
	    Rlarfg(n - i + 1, &A[i + i * lda], &A[i + min(i + 1, n) * lda], lda, &taup[i]);
	    d[i] = A[i + i * lda];
	    A[i + i * lda] = One;
//Apply G(i) to A(i+1:m,i:n) from the right
	    if (i < m) {
		Rlarf("Right", m - i, n - i + 1, &A[i + i * lda], lda, taup[i], &A[i + 1 + i * lda], lda, &work[0]);
	    }
	    A[i + i * lda] = d[i];
	    if (i < m) {
//Generate elementary reflector H(i) to annihilate
//A(i+2:m,i)
		Rlarfg(m - i, &A[i + 1 + i * lda], &A[min(i + 2, m) + i * lda], 1, &tauq[i]);
		e[i] = A[i + 1 + i * lda];
		A[i + 1 + i * lda] = One;
//Apply H(i) to A(i+1:m,i+1:n) from the left
		Rlarf("Left", m - i, n - i, &A[i + 1 + i * lda], 1, tauq[i], &A[i + 1 + (i + 1) * lda], lda, &work[0]);
		A[i + 1 + i * lda] = e[i];
	    } else {
		tauq[i] = Zero;
	    }
	}
    }
    return;
}
