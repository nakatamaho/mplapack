/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgglse.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgglse(INTEGER m, INTEGER n, INTEGER p, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, COMPLEX * c, COMPLEX * d,
	    COMPLEX * x, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER nb, mn, nr, nb1, nb2, nb3, nb4, lopt;
    INTEGER lwkmin, lwkopt;
    INTEGER lquery;
    REAL One = 1.0;
//Test the input parameters
    *info = 0;
    mn = min(m, n);
    lquery = lwork == -1;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (p < 0 || p > n || p < n - m) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, p)) {
	*info = -7;
    }
//Calculate workspace
    if (*info == 0) {
	if (n == 0) {
	    lwkmin = 1;
	    lwkopt = 1;
	} else {
	    nb1 = iMlaenv(1, "Cgeqrf", " ", m, n, -1, -1);
	    nb2 = iMlaenv(1, "Cgerqf", " ", m, n, -1, -1);
	    nb3 = iMlaenv(1, "Cunmqr", " ", m, n, p, -1);
	    nb4 = iMlaenv(1, "Cunmrq", " ", m, n, p, -1);
	    nb = max(max(max(nb1, nb2), nb3), nb4);
	    lwkmin = m + n + p;
	    lwkopt = p + mn + max(m, n) * nb;
	}
	work[1] = lwkopt;
	if (lwork < lwkmin && !lquery) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	Mxerbla("Cgglse", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Compute the GRQ factorization of matrices B and A:
//       B*Q' = (  0  T12 ) P   Z'*A*Q' = ( R11 R12 ) N-P
//                N-P  P                  (  0  R22 ) M+P-N
//                                          N-P  P

//where T12 and R11 are upper triangular, and Q and Z are
//unitary.
    Cggrqf(p, m, n, &B[0], ldb, &work[0], &A[0], lda, &work[p + 1], &work[p + mn + 1], lwork - p - mn, info);
    lopt = (INTEGER) cast2double(work[p + mn + 1].real());
//Update c = Z'*c = ( c1 ) N-P
//                  ( c2 ) M+P-N
    Cunmqr("Left", "Conjugate Transpose", m, 1, mn, &A[0], lda, &work[p + 1], &c[1], max((INTEGER) 1, m), &work[p + mn + 1], max((INTEGER) 1, m), info);
    lopt = max(lopt, (INTEGER) cast2double(work[p + mn + 1].real()));
//Solve T12*x2 = d for x2
    if (p > 0) {
	Ctrtrs("Upper", "No transpose", "Non-unit", p, 1, &B[(n - p + 1) * ldb + 1], ldb, &d[0], p, info);
	if (*info > 0) {
	    *info = 1;
	    return;
	}
//Put the solution in X
	Ccopy(p, &d[0], 1, &x[n - p + 1], 1);
//Update c1
	Cgemv("No transpose", n - p, p, (COMPLEX) - One, &A[(n - p + 1) * lda], lda, &d[0], 1, One, &c[1], 1);
    }
//Solve R11*x1 = c1 for x1
    if (n > p) {
	Ctrtrs("Upper", "No transpose", "Non-unit", n - p, 1, &A[0], lda, &c[1], n - p, info);
	if (*info > 0) {
	    *info = 2;
	    return;
	}
//Put the solutions in X
	Ccopy(n - p, &c[1], 1, &x[0], 1);
    }
//Compute the residual vector:
    if (m < n) {
	nr = m + p - n;
	if (nr > 0) {
	    Cgemv("No transpose", nr, n - m, (COMPLEX) - One, &A[n - p + 1 + (m + 1) * lda], lda, &d[nr + 1], 1, One, &c[n - p + 1], 1);
	}
    } else {
	nr = p;
    }
    if (nr > 0) {
	Ctrmv("Upper", "No transpose", "Non unit", nr, &A[n - p + 1 + (n - p + 1) * lda], lda, &d[0], 1);
	Caxpy(nr, (COMPLEX) - One, &d[0], 1, &c[n - p + 1], 1);
    }
//Backward transformation x = Q'*x
    Cunmrq("Left", "Conjugate Transpose", n, 1, p, &B[0], ldb, &work[0], &x[0], n, &work[p + mn + 1], lwork - p - mn, info);
    work[1] = p + mn + max(lopt, (INTEGER) cast2double(work[p + mn + 1].real()));
    return;
}
