/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cpbtrs.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cpbtrs(const char *uplo, INTEGER n, INTEGER kd, INTEGER nrhs, COMPLEX * AB, INTEGER ldab, COMPLEX * B, INTEGER ldb, INTEGER * info)
{
    INTEGER j;
    INTEGER upper;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (kd < 0) {
	*info = -3;
    } else if (nrhs < 0) {
	*info = -4;
    } else if (ldab < kd + 1) {
	*info = -6;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -8;
    }
    if (*info != 0) {
	Mxerbla("Cpbtrs", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0 || nrhs == 0) {
	return;
    }
    if (upper) {
//Solve A*X = B where A = U'*U.
	for (j = 0; j < nrhs; j++) {
//Solve U'*X = B, overwriting B with X.
	    Ctbsv("Upper", "Conjugate transpose", "Non-unit", n, kd, AB, ldab, &B[j * ldb + 1], 1);
//Solve U*X = B, overwriting B with X.
	    Ctbsv("Upper", "No transpose", "Non-unit", n, kd, AB, ldab, &B[j * ldb + 1], 1);
	}
    } else {
//Solve A*X = B where A = L*L'.
	for (j = 0; j < nrhs; j++) {
//Solve L*X = B, overwriting B with X.
	    Ctbsv("Lower", "No transpose", "Non-unit", n, kd, AB, ldab, &B[j * ldb + 1], 1);
//Solve L'*X = B, overwriting B with X.
	    Ctbsv("Lower", "Conjugate transpose", "Non-unit", n, kd, AB, ldab, &B[j * ldb + 1], 1);
	}
    }
    return;
}
