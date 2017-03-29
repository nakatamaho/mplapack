/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chpgvd.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Chpgvd(INTEGER itype, const char *jobz, const char *uplo, INTEGER n, COMPLEX * ap, COMPLEX * bp, REAL * w, COMPLEX
       * z, INTEGER ldz, COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER lrwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER j, neig;
    INTEGER lwmin;
    char trans;
    INTEGER upper, wantz;
    INTEGER liwmin;
    INTEGER lrwmin;
    INTEGER lquery;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    wantz = Mlsame(jobz, "V");
    upper = Mlsame(uplo, "U");
    lquery = lwork == -1 || lrwork == -1 || liwork == -1;

    *info = 0;
    if (itype < 1 || itype > 3) {
	*info = -1;
    } else if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -2;
    } else if (!(upper || Mlsame(uplo, "L"))) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (ldz < 1 || (wantz && ldz < n)) {
	*info = -9;
    }
    if (*info == 0) {
	if (n <= 1) {
	    lwmin = 1;
	    liwmin = 1;
	    lrwmin = 1;
	} else {
	    if (wantz) {
		lwmin = n << 1;
		lrwmin = n * 5 + 1 + (n * n * 2);
		liwmin = n * 5 + 3;
	    } else {
		lwmin = n;
		lrwmin = n;
		liwmin = 1;
	    }
	}

	work[1] = lwmin;
	rwork[1] = (double) lrwmin;
	iwork[1] = liwmin;
	if (lwork < lwmin && !lquery) {
	    *info = -11;
	} else if (lrwork < lrwmin && !lquery) {
	    *info = -13;
	} else if (liwork < liwmin && !lquery) {
	    *info = -15;
	}
    }

    if (*info != 0) {

	Mxerbla("CHPGVD", -(*info));
	return;
    } else if (lquery) {
	return;
    }

/*     Quick return if possible */

    if (n == 0) {
	return;
    }

/*     Form a Cholesky factorization of B. */

    Cpptrf(uplo, n, &bp[1], info);
    if (*info != 0) {
	*info = n + *info;
	return;
    }
//Transform problem to standard eigenvalue problem and solve.
    Chpgst(&itype, uplo, n, &ap[1], &bp[1], info);
    Chpevd(jobz, uplo, n, &ap[1], &w[1], &z[0], ldz, &work[0], lwork, &rwork[1], lrwork, &iwork[1], liwork, info);
    mtemp1 = lwmin, mtemp2 = work[1].real();
    lwmin = (INTEGER) cast2double(max(mtemp1, mtemp2));
    mtemp1 = lrwmin;
    lrwmin = (INTEGER) cast2double(max(mtemp1, rwork[1]));
    mtemp1 = liwmin, mtemp2 = iwork[1];
    liwmin = (INTEGER) cast2double(max(mtemp1, mtemp2));
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
		trans = 'C';
	    }
	    for (j = 0; j < neig; j++) {
		Ctpsv(uplo, &trans, "Non-unit", n, &bp[1], &z[j * ldz + 1], 1);
	    }
	} else if (itype == 3) {
//For B*A*x=(lambda)*x;
//backtransform eigenvectors: x = L*y or U'*y
	    if (upper) {
		trans = 'C';
	    } else {
		trans = 'N';
	    }
	    for (j = 0; j < neig; j++) {
		Ctpmv(uplo, &trans, "Non-unit", n, &bp[1], &z[j * ldz + 1], 1);
	    }
	}
    }
    work[1] = lwmin;
    rwork[1] = (double) lrwmin;
    iwork[1] = liwmin;
    return;
}
