/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Claein.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Claein(INTEGER rightv, INTEGER noinit, INTEGER n, COMPLEX * h, INTEGER ldh, COMPLEX w, COMPLEX * v, COMPLEX * B, INTEGER ldb, REAL * rwork, REAL eps3, REAL smlnum, INTEGER * info)
{
    INTEGER i, j;
    COMPLEX x, ei, ej;
    INTEGER its, ierr;
    COMPLEX temp;
    REAL scale;
    char trans;
    REAL rtemp, rootn, vnorm;
    char normin;
    REAL nrmsml, growto;
    REAL mtemp1, mtemp2;
    REAL One = 1.0, Zero = 0.0, Tenth = 0.0000000001;

    *info = 0;
//GROWTO is the threshold used in the acceptance test for an
//eigenvector.

    rootn = (REAL) sqrt(n);
    growto = Tenth / rootn;

    mtemp1 = One, mtemp2 = eps3 * rootn;
    nrmsml = max(mtemp1, mtemp2);
    nrmsml = nrmsml * smlnum;

//Form B = H - W*I (except that the subdiagonal elements are not
//stored).
    for (j = 0; j < n; j++) {
	for (i = 0; i < j - 1; i++) {
	    B[i + j * ldb] = h[i + j * ldh];
	}
	B[j + j * ldb] = h[j + j * ldh] - w;
    }
    if (noinit) {
//Initialize V.
	for (i = 0; i < n; i++) {
	    v[i] = eps3;

	}
    } else {
//Scale supplied initial vector.
	vnorm = RCnrm2(n, &v[1], 1);
	mtemp1 = eps3 * rootn / max(vnorm, nrmsml);
	CRscal(n, mtemp1, &v[1], 1);
    }
    if (rightv) {
//LU decomposition with partial pivoting of B, replacing zero
//pivots by EPS3.
	for (i = 0; i < n - 1; i++) {
	    ei = h[i + 1 + i * ldh];
	    if (RCabs1(B[i + i * ldb]) < RCabs1(ei)) {
//Interchange rows and eliminate.
		x = Cladiv(B[i + i * ldb], ei);
		B[i + i * ldb] = ei;
		for (j = i + 1; j <= n; j++) {
		    temp = B[i + 1 + j * ldb];
		    B[i + 1 + j * ldb] = B[i + j * ldb] - x * temp;
		    B[i + j * ldb] = temp;
		}
	    } else {
//Eliminate without interchange.
		if (B[i + i * ldb] == Zero) {
		    B[i + i * ldb] = eps3;
		}
		x = Cladiv(ei, B[i + i * ldb]);
		if (x != Zero) {
		    for (j = i + 1; j <= n; j++) {
			B[i + 1 + j * ldb] = B[i + 1 + j * ldb] - x * B[i + j * ldb];
		    }
		}
	    }
	}
	if (B[n + n * ldb] == Zero) {
	    B[n + n * ldb] = eps3;
	}
	trans = 'N';
    } else {
//UL decomposition with partial pivoting of B, replacing zero
//pivots by EPS3.
	for (j = n; j >= 2; j--) {
	    ej = h[j + (j - 1) * ldh];
	    if (RCabs1(B[j + j * ldb]) < RCabs1(ej)) {
//Interchange columns and eliminate.
		x = Cladiv(B[j + j * ldb], ej);
		B[j + j * ldb] = ej;
		for (i = 0; i < j - 1; i++) {
		    temp = B[i + (j - 1) * ldb];
		    B[i + (j - 1) * ldb] = B[i + j * ldb] - x * temp;
		    B[i + j * ldb] = temp;

		}
	    } else {
//Eliminate without interchange.
		if (B[j + j * ldb] == Zero) {
		    B[j + j * ldb] = eps3;
		}
		x = Cladiv(ej, B[j + j * ldb]);
		if (x != Zero) {
		    for (i = 0; i < j - 1; i++) {
			B[i + (j - 1) * ldb] = B[i + (j - 1) * ldb] - x * B[i + j * ldb];
		    }
		}
	    }
	}
	if (B[ldb + 1] == Zero) {
	    B[ldb + 1] = eps3;
	}
	trans = 'C';
    }
    normin = 'N';
    for (its = 1; its <= n; its++) {
//Solve U*x = scale*v for a right eigenvector
//  or U'*x = scale*v for a left eigenvector,
//overwriting x on v.
	Clatrs("Upper", &trans, "Nonunit", &normin, n, &B[0], ldb, &v[1], &scale, &rwork[1], &ierr);
	normin = 'Y';
//Test for sufficient growth in the norm of v.
	vnorm = RCasum(n, &v[1], 1);
	if (vnorm >= growto * scale) {
	    goto L120;
	}
//Choose new orthogonal starting vector and try again.
	rtemp = eps3 / (rootn + One);
	v[1] = eps3;
	for (i = 1; i < n; i++) {
	    v[i] = rtemp;
	}
	v[n - its + 1] = v[n - its + 1] - eps3 * rootn;
    }
//Failure to find eigenvector in N iterations.
    *info = 1;
  L120:
//Normalize eigenvector.
    i = iCamax(n, &v[1], 1);
    CRscal(n, One / RCabs1(v[i]), &v[1], 1);
    return;
}
