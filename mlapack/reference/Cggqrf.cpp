/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cggqrf.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cggqrf(INTEGER n, INTEGER m, INTEGER p, COMPLEX * A, INTEGER lda, COMPLEX * taua, COMPLEX * B, INTEGER ldb, COMPLEX * taub, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER nb, nb1, nb2, nb3, lopt;
    INTEGER lwkopt;
    INTEGER lquery;

    *info = 0;
    nb1 = iMlaenv(1, "Cgeqrf", " ", n, m, -1, -1);
    nb2 = iMlaenv(1, "Cgerqf", " ", n, p, -1, -1);
    nb3 = iMlaenv(1, "Cunmqr", " ", n, m, p, -1);
    nb = max(max(nb1, nb2), nb3);
    lwkopt = max(max(n, m), p) * nb;
    work[1] = lwkopt;
    lquery = lwork == -1;
    if (n < 0) {
	*info = -1;
    } else if (m < 0) {
	*info = -2;
    } else if (p < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -8;
    } else {
	if (lwork < max(max(max((INTEGER) 1, n), m), p) && !lquery) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	Mxerbla("Cggrqf", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//QR factorization of N-by-M matrix A: A = Q*R
    Cgeqrf(n, m, &A[0], lda, taua, &work[0], lwork, info);
    lopt = (INTEGER) cast2double(work[1].real());
//Update B := Q'*B.
    Cunmqr("Left", "Conjugate Transpose", n, p, min(n, m), &A[0], lda, &taua[1], &B[0], ldb, &work[0], lwork, info);
    lopt = max(lopt, (INTEGER) cast2double(work[1].real()));
//RQ factorization of N-by-P matrix B: B = T*Z.
    Cgerqf(n, p, &B[0], ldb, &taub[1], &work[0], lwork, info);
    work[1] = max(lopt, (INTEGER) cast2double(work[1].real()));
    return;
}
