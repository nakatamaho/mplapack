/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rggrqf.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rggrqf(INTEGER m, INTEGER p, INTEGER n, REAL * A, INTEGER lda, REAL * taua, REAL * B, INTEGER ldb, REAL * taub, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER nb, nb1, nb2, nb3, lopt;
    INTEGER lwkopt;
    INTEGER lquery;

//Test the input parameters
    *info = 0;
    nb1 = iMlaenv(1, "Rgerqf", " ", m, n, -1, -1);
    nb2 = iMlaenv(1, "Rgeqrf", " ", p, n, -1, -1);
    nb3 = iMlaenv(1, "Rormrq", " ", m, n, p, -1);
    nb = max(max(nb1, nb2), nb3);
    lwkopt = max(max(n, m), p) * nb;
    work[0] = (REAL) double (lwkopt);
    lquery = lwork == -1;
    if (m < 0) {
	*info = -1;
    } else if (p < 0) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, p)) {
	*info = -8;
    } else {
	if (lwork < max(max(max((INTEGER) 1, m), p), n) && !lquery) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	Mxerbla("Rggrqf", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//RQ factorization of M-by-N matrix A: A = R*Q
    Rgerqf(m, n, &A[0], lda, &taua[1], &work[0], lwork, info);
    lopt = (INTEGER) cast2double(work[0]);
//Update B := B*Q'
    Rormrq("Right", "Transpose", p, n, min(m, n), &A[max((INTEGER) 1, m - n + 1) + lda], lda, &taua[1], &B[0], ldb, &work[0], lwork, info);
    lopt = max(lopt, (INTEGER) cast2double(work[0]));
//QR factorization of P-by-N matrix B: B = Z*T
    Rgeqrf(p, n, &B[0], ldb, &taub[0], &work[0], lwork, info);
    work[0] = (REAL) double (max(lopt, (INTEGER) cast2double(work[1])));
    return;
}
