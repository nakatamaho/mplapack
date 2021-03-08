/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rptsvx.cpp,v 1.4 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rptsvx(const char *fact, INTEGER n, INTEGER nrhs,
	    REAL * d, REAL * e, REAL * df, REAL * ef, REAL * B, INTEGER ldb, REAL * x, INTEGER ldx, REAL * rcond, REAL * ferr, REAL * berr, REAL * work, INTEGER * info)
{
    REAL anorm;
    INTEGER nofact;
    REAL Zero = 0.0;

//Test the input parameters.
    *info = 0;
    nofact = Mlsame(fact, "N");
    if (!nofact && !Mlsame(fact, "F")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (nrhs < 0) {
	*info = -3;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -9;
    } else if (ldx < max((INTEGER) 1, n)) {
	*info = -11;
    }
    if (*info != 0) {
	Mxerbla("Rptsvx", -(*info));
	return;
    }
    if (nofact) {
//Compute the L*D*L' (or U'*D*U) factorization of A.
	Rcopy(n, &d[0], 1, &df[1], 1);
	if (n > 1) {
	    Rcopy(n - 1, &e[0], 1, &ef[1], 1);
	}
	Rpttrf(n, &df[1], &ef[1], info);
//Return if INFO is non-zero.
	if (*info > 0) {
	    *rcond = Zero;
	    return;
	}
    }
//Compute the norm of the matrix A.
    anorm = Rlanst("1", n, &d[0], &e[0]);
//Compute the reciprocal of the condition number of A.
    Rptcon(n, &df[1], &ef[1], anorm, rcond, &work[0], info);
//Compute the solution vectors X.
    Rlacpy("Full", n, nrhs, &B[0], ldb, &x[0], ldx);
    Rpttrs(n, nrhs, &df[1], &ef[1], &x[0], ldx, info);
//Use iterative refinement to improve the computed solutions and
//compute error bounds and backward error estimates for them.
    Rptrfs(n, nrhs, &d[0], &e[0], &df[1], &ef[1], &B[0], ldb, &x[0], ldx, &ferr[1], &berr[1], &work[0], info);
//Set INFO = N+1 if the matrix is singular to working precision.
    if (*rcond < Rlamch("Epsilon")) {
	*info = n + 1;
    }
    return;
}
