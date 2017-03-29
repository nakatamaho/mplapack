/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasdq.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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
Rlasdq(const char *uplo, INTEGER sqre, INTEGER n, INTEGER ncvt, INTEGER nru, INTEGER ncc,
       REAL * d, REAL * e, REAL * vt, INTEGER ldvt, REAL * u, INTEGER ldu, REAL * c, INTEGER ldc, REAL * work, INTEGER * info)
{
    INTEGER i, j;
    REAL r, cs, sn;
    INTEGER np1, isub;
    REAL smin;
    INTEGER sqre1;
    INTEGER iuplo;
    INTEGER rotate;
    REAL Zero = 0.0;

//Test the input parameters.
    *info = 0;
    iuplo = 0;
    if (Mlsame(uplo, "U")) {
	iuplo = 1;
    }
    if (Mlsame(uplo, "L")) {
	iuplo = 2;
    }
    if (iuplo == 0) {
	*info = -1;
    } else if (sqre < 0 || sqre > 1) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (ncvt < 0) {
	*info = -4;
    } else if (nru < 0) {
	*info = -5;
    } else if (ncc < 0) {
	*info = -6;
    } else if ((ncvt == 0 && ldvt < 1) || (ncvt > 0 && ldvt < max((INTEGER) 1, n))) {
	*info = -10;
    } else if (ldu < max((INTEGER) 1, nru)) {
	*info = -12;
    } else if ((ncc == 0 && ldc < 1) || (ncc > 0 && ldc < max((INTEGER) 1, n))) {
	*info = -14;
    }
    if (*info != 0) {
	Mxerbla("Rlasdq", -(*info));
	return;
    }
    if (n == 0) {
	return;
    }
//ROTATE is true if any singular vectors desired, false otherwise
    rotate = ncvt > 0 || nru > 0 || ncc > 0;
    np1 = n + 1;
    sqre1 = sqre;

//If matrix non-square upper bidiagonal, rotate to be lower
//bidiagonal.  The rotations are on the right.
    if (iuplo == 1 && sqre1 == 1) {
	for (i = 0; i < n - 1; i++) {
	    Rlartg(d[i], e[i], &cs, &sn, &r);
	    d[i] = r;
	    e[i] = sn * d[i + 1];
	    d[i + 1] = cs * d[i + 1];
	    if (rotate) {
		work[i] = cs;
		work[n + i] = sn;
	    }
	}
	Rlartg(d[n], e[n], &cs, &sn, &r);
	d[n] = r;
	e[n] = Zero;
	if (rotate) {
	    work[n] = cs;
	    work[n + n] = sn;
	}
	iuplo = 2;
	sqre1 = 0;
//Update singular vectors if desired.
	if (ncvt > 0) {
	    Rlasr("L", "V", "F", np1, ncvt, &work[0], &work[np1], &vt[0], ldvt);
	}
    }
//If matrix lower bidiagonal, rotate to be upper bidiagonal
//by applying Givens rotations on the left.
    if (iuplo == 2) {
	for (i = 0; i < n - 1; i++) {
	    Rlartg(d[i], e[i], &cs, &sn, &r);
	    d[i] = r;
	    e[i] = sn * d[i + 1];
	    d[i + 1] = cs * d[i + 1];
	    if (rotate) {
		work[i] = cs;
		work[n + i] = sn;
	    }
	}
//If matrix (N+1)-by-N lower bidiagonal, one additional
//rotation is needed.
	if (sqre1 == 1) {
	    Rlartg(d[n], e[n], &cs, &sn, &r);
	    d[n] = r;
	    if (rotate) {
		work[n] = cs;
		work[n + n] = sn;
	    }
	}
//Update singular vectors if desired.
	if (nru > 0) {
	    if (sqre1 == 0) {
		Rlasr("R", "V", "F", nru, n, &work[0], &work[np1], &u[0], ldu);
	    } else {
		Rlasr("R", "V", "F", nru, np1, &work[0], &work[np1], &u[0], ldu);
	    }
	}
	if (ncc > 0) {
	    if (sqre1 == 0) {
		Rlasr("L", "V", "F", n, ncc, &work[0], &work[np1], &c[0], ldc);
	    } else {
		Rlasr("L", "V", "F", np1, ncc, &work[0], &work[np1], &c[0], ldc);
	    }
	}
    }
//Call DBDSQR to compute the SVD of the reduced real
//N-by-N upper bidiagonal matrix.
    Rbdsqr("U", n, ncvt, nru, ncc, &d[0], &e[0], &vt[0], ldvt, &u[0], ldu, &c[0], ldc, &work[0], info);

//Sort the singular values into ascending order (insertion sort on
//singular values, but only one transposition per singular vector)
    for (i = 0; i < n; i++) {
//Scan for smallest D(I).
	isub = i;
	smin = d[i];
	for (j = i + 1; j <= n; j++) {
	    if (d[j] < smin) {
		isub = j;
		smin = d[j];
	    }
	}
	if (isub != i) {
//Swap singular values and vectors.
	    d[isub] = d[i];
	    d[i] = smin;
	    if (ncvt > 0) {
		Rswap(ncvt, &vt[isub + ldvt], ldvt, &vt[i + ldvt], ldvt);
	    }
	    if (nru > 0) {
		Rswap(nru, &u[isub * ldu + 1], 1, &u[i * ldu + 1]
		      , 1);
	    }
	    if (ncc > 0) {
		Rswap(ncc, &c[isub + ldc], ldc, &c[i + ldc], ldc);
	    }
	}
    }
    return;
}
