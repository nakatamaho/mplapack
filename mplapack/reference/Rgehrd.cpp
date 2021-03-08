/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgehrd.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgehrd(INTEGER n, INTEGER ilo, INTEGER ihi, REAL * A, INTEGER lda, REAL * tau, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, j;
    REAL t[4160];
    INTEGER ib;
    REAL ei;
    INTEGER nb, nh, nx, iws;
    INTEGER nbmin, iinfo;
    INTEGER ldwork, lwkopt;
    INTEGER lquery;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters
    *info = 0;
    nb = max((INTEGER) 64, iMlaenv(1, "Rgehrd", " ", n, ilo, ihi, -1));
    lwkopt = n * nb;
    work[1] = (double) lwkopt;
    lquery = lwork == -1;
    if (n < 0) {
	*info = -1;
    } else if (ilo < 1 || ilo > max((INTEGER) 1, n)) {
	*info = -2;
    } else if (ihi < min(ilo, n) || ihi > n) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    } else if (lwork < max((INTEGER) 1, n) && !lquery) {
	*info = -8;
    }
    if (*info != 0) {
	Mxerbla("Rgehrd", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
    for (i = 0; i < ilo - 1; i++) {
	tau[i] = Zero;
    }
    for (i = max((INTEGER) 1, ihi); i <= n - 1; i++) {
	tau[i] = Zero;
    }
//Quick return if possible
    nh = ihi - ilo + 1;
    if (nh <= 1) {
	work[1] = One;
	return;
    }
//Determine the block size
    nb = max((INTEGER) 64, iMlaenv(1, "Rgehrd", " ", n, ilo, ihi, -1));
    nbmin = 2;
    iws = 1;
    if (nb > 1 && nb < nh) {
//Determine when to cross over from blocked to unblocked code
//(last block is always handled by unblocked code)
	nx = max(nb, iMlaenv(3, "Rgehrd", " ", n, ilo, ihi, -1));
	if (nx < nh) {
//Determine if workspace is large enough for blocked code
	    iws = n * nb;
	    if (lwork < iws) {
//Not enough workspace to use optimal NB:  determine the
//minimum value of NB, and reduce NB or force use of
//unblocked code
		nbmin = max((INTEGER) 2, iMlaenv(2, "Rgehrd", " ", n, ilo, ihi, -1));
		if (lwork >= n * nbmin) {
		    nb = lwork / n;
		} else {
		    nb = 1;
		}
	    }
	}
    }
    ldwork = n;
    if (nb < nbmin || nb >= nh) {
//Use unblocked code below
	i = ilo;
    } else {
//Use blocked code
	for (i = ilo; i <= ihi - 1 - nx; i = i + nb) {
	    ib = min(nb, ihi - i);
//Reduce columns i:i+ib-1 to Hessenberg form, returning the
//matrices V and T of the block reflector H = I - V*T*V'
//which performs the reduction, and also the matrix Y = A*V*T
	    Rlahr2(ihi, i, ib, &A[i * lda], lda, &tau[i], t, 65, work, ldwork);
//Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
//right, computing  A := A - Y * V'. V(i+ib,ib-1) must be set
//to 1
	    ei = A[i + ib + (i + ib - 1) * lda];
	    A[i + ib + (i + ib - 1) * lda] = One;
	    Rgemm("No transpose", "Transpose", ihi, ihi - i - ib + 1, ib, -One, work, ldwork, &A[i + ib + i * lda], lda, One, &A[(i + ib) * lda], lda);
	    A[i + ib + (i + ib - 1) * lda] = ei;
//Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
//right
	    Rtrmm("Right", "Lower", "Transpose", "Unit", i, ib - 1, One, &A[i + 1 + i * lda], lda, work, ldwork);
	    for (j = 0; j < ib - 2; j++) {
		Raxpy(i, -One, &work[ldwork * j + 1], 1, &A[(i + j + 1) * lda], 1);
	    }
//Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
//left
	    Rlarfb("Left", "Transpose", "Forward", "Columnwise", ihi - i, n - i - ib + 1, ib, &A[i + 1 + i * lda], lda, t, 65, &A[i + 1 + (i + ib) * lda], lda, work, ldwork);
	}
    }
//Use unblocked code to reduce the rest of the matrix
    Rgehd2(n, i, ihi, A, lda, tau, work, &iinfo);
    work[1] = (REAL) double (iws);
    return;
}
