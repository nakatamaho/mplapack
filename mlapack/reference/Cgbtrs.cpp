/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgbtrs.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgbtrs(const char *trans, INTEGER n, INTEGER kl, INTEGER ku, INTEGER nrhs, COMPLEX * AB, INTEGER ldab, INTEGER * ipiv, COMPLEX * B, INTEGER ldb, INTEGER * info)
{
    INTEGER i, j, l, kd, lm;
    INTEGER lnoti;
    INTEGER notran;
    REAL One = 1.0;

//Test the input parameters.
    *info = 0;
    notran = Mlsame(trans, "N");
    if (!notran && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (kl < 0) {
	*info = -3;
    } else if (ku < 0) {
	*info = -4;
    } else if (nrhs < 0) {
	*info = -5;
    } else if (ldab < (kl * 2) + ku + 1) {
	*info = -7;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -10;
    }
    if (*info != 0) {
	Mxerbla("Cgbtrs", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0 || nrhs == 0) {
	return;
    }
    kd = ku + kl + 1;
    lnoti = kl > 0;
    if (notran) {
//Solve  A*X = B.
//Solve L*X = B, overwriting B with X.
//L is represented as a product of permutations and unit lower
//triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
//where each transformation L(i) is a rank-one modification of
//the identity matrix.
	if (lnoti) {
	    for (j = 0; j < n - 1; j++) {
		lm = min(kl, n - j);
		l = ipiv[j];
		if (l != j) {
		    Cswap(nrhs, &B[l + ldb], ldb, &B[j + ldb], ldb);
		}
		Cgeru(lm, nrhs, (COMPLEX) - One, &AB[kd + 1 + j * ldab], 1, &B[j + ldb], ldb, &B[j + 1 + ldb], ldb);
	    }
	}
	for (i = 0; i < nrhs; i++) {
//Solve U*X = B, overwriting B with X.
	    Ctbsv("Upper", "No transpose", "Non-unit", n, kl + ku, AB, ldab, &B[i * ldb + 1], 12);
	}
    } else if (Mlsame(trans, "T")) {
//Solve A**T * X = B.
	for (i = 0; i < nrhs; i++) {
//Solve U**T * X = B, overwriting B with X.
	    Ctbsv("Upper", "Transpose", "Non-unit", n, kl + ku, AB, ldab, &B[i * ldb + 1], 1);
	}
//Solve L**T * X = B, overwriting B with X.
	if (lnoti) {
	    for (j = n - 1; j >= 1; j--) {
		lm = min(kl, n - j);
		Cgemv("Transpose", lm, nrhs, (COMPLEX) - One, &B[j + 1 + ldb], ldb, &AB[kd + 1 + j * ldab], 1, (COMPLEX) One, &B[j + ldb], ldb);
		l = ipiv[j];
		if (l != j) {
		    Cswap(nrhs, &B[l + ldb], ldb, &B[j + ldb], ldb);
		}
	    }
	}
    } else {
//Solve A**H * X = B.
	for (i = 0; i < nrhs; i++) {
//Solve U**H * X = B, overwriting B with X.
	    Ctbsv("Upper", "Conjugate transpose", "Non-unit", n, kl + ku, AB, ldab, &B[i * ldb + 1], 1);
	}
//Solve L**H * X = B, overwriting B with X.
	if (lnoti) {
	    for (j = n - 1; j >= 1; j--) {
		lm = min(kl, n - j);
		Clacgv(nrhs, &B[j + ldb], ldb);
		Cgemv("Conjugate transpose", lm, nrhs, (COMPLEX) - One, &B[j + 1 + ldb], ldb, &AB[kd + 1 + j * ldab], 1, (COMPLEX) One, &B[j + ldb], ldb);
		Clacgv(nrhs, &B[j + ldb], ldb);
		l = ipiv[j];
		if (l != j) {
		    Cswap(nrhs, &B[l + ldb], ldb, &B[j + ldb], ldb);
		}
	    }
	}
    }
    return;
}
