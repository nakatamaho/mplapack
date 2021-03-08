/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsygv.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rsygv(INTEGER itype, const char *jobz, const char *uplo, INTEGER n, REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL * w, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER nb, neig;
    char trans;
    INTEGER upper;
    INTEGER wantz;
    INTEGER lwkmin;
    INTEGER lwkopt;
    INTEGER lquery;
    REAL One = 1.0;

//Test the input parameters.
    wantz = Mlsame(jobz, "V");
    upper = Mlsame(uplo, "U");
    lquery = lwork == -1;

    *info = 0;
    if (itype < 1 || itype > 3) {
	*info = -1;
    } else if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -2;
    } else if (!(upper || Mlsame(uplo, "L"))) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -8;
    }

    if (*info == 0) {
	lwkmin = max((INTEGER) 1, n * 3 - 1);
	nb = iMlaenv(1, "Rsytrd", uplo, n, -1, -1, -1);
	lwkopt = max(lwkmin, (nb + 2) * n);
	work[0] = lwkopt;
	if (lwork < lwkmin && !lquery) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	Mxerbla("Rsygv ", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0)
	return;

//Form a Cholesky factorization of B.
    Rpotrf(uplo, n, &B[0], ldb, info);
    if (*info != 0) {
	*info = n + *info;
	return;
    }
//Transform problem to standard eigenvalue problem and solve.
    Rsygst(itype, uplo, n, &A[0], lda, &B[0], ldb, info);
    Rsyev(jobz, uplo, n, &A[0], lda, &w[0], &work[0], lwork, info);
    if (wantz) {
//Backtransform eigenvectors to the original problem.
	neig = n;
	if (*info > 0) {
	    neig = *info - 1;
	}
	if (itype == 1 || itype == 2) {
//For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
//backtransform eigenvectors: x = inv(L)'*y or inv(U)*y
	    if (upper) {
		trans = 'N';
	    } else {
		trans = 'T';
	    }
	    Rtrsm("Left", uplo, &trans, "Non-unit", n, neig, One, &B[0], ldb, &A[0], lda);
	} else if (itype == 3) {
//For B*A*x=(lambda)*x;
//backtransform eigenvectors: x = L*y or U'*y
	    if (upper) {
		trans = 'T';
	    } else {
		trans = 'N';
	    }
	    Rtrmm("Left", uplo, &trans, "Non-unit", n, neig, One, &B[0], ldb, &A[0], lda);
	}
    }
    work[0] = lwkopt;
    return;
}
