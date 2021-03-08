/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rbdsdc.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rbdsdc(const char *uplo, const char *compq, INTEGER n, REAL * d, REAL * e, REAL * u, INTEGER ldu, REAL * vt, INTEGER ldvt,
	    REAL * q, INTEGER * iq, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j, k;
    REAL p, r;
    INTEGER z = 0, ic = 0, ii, kk;
    REAL cs;
    INTEGER is = 0, iu, iuplo;
    REAL sn;
    INTEGER nm1;
    REAL eps;
    INTEGER ivt = 0, difl = 0, difr = 0, ierr, perm, mlvl, sqre;
    INTEGER givcol = 0, poles = 0, start;
    INTEGER icompq, nsize;
    REAL orgnrm;
    INTEGER givnum = 0, givptr, qstart, smlsiz, wstart, smlszp;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;

//Test the input parameters.
    *info = 0;
    iuplo = 0;
    if (Mlsame(uplo, "U")) {
	iuplo = 1;
    }
    if (Mlsame(uplo, "L")) {
	iuplo = 2;
    }
    if (Mlsame(compq, "N")) {
	icompq = 0;
    } else if (Mlsame(compq, "P")) {
	icompq = 1;
    } else if (Mlsame(compq, "I")) {
	icompq = 2;
    } else {
	icompq = -1;
    }
    if (iuplo == 0) {
	*info = -1;
    } else if (icompq < 0) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (ldu < 1 || (icompq == 2 && ldu < n)) {
	*info = -7;
    } else if (ldvt < 1 || (icompq == 2 && ldvt < n)) {
	*info = -9;
    }
    if (*info != 0) {
	Mxerbla("Rbdsdc", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    smlsiz = iMlaenv(9, "Rbdsdc", " ", 0, 0, 0, 0);
    if (n == 1) {
	if (icompq == 1) {
	    q[1] = sign(One, d[0]);
	    q[smlsiz * n + 1] = One;
	} else if (icompq == 2) {
	    u[ldu + 1] = sign(One, d[0]);
	    vt[ldvt + 1] = One;
	}
	d[1] = abs(d[1]);
	return;
    }
    nm1 = n - 1;
/*     If matrix lower bidiagonal, rotate to be upper bidiagonal */
/*     by applying Givens rotations on the left */
    wstart = 1;
    qstart = 3;
    if (icompq == 1) {
	Rcopy(n, &d[0], 1, &q[1], 1);
	Rcopy(n - 1, &e[0], 1, &q[n + 1], 1);
    }
    if (iuplo == 2) {
	qstart = 5;
	wstart = (n << 1) - 1;
	for (i = 0; i < n - 1; i++) {
	    Rlartg(d[i], e[i], &cs, &sn, &r);
	    d[i] = r;
	    e[i] = sn * d[i + 1];
	    d[i + 1] = cs * d[i + 1];
	    if (icompq == 1) {
		q[i + (n << 1)] = cs;
		q[i + n * 3] = sn;
	    } else if (icompq == 2) {
		work[i] = cs;
		work[nm1 + i] = -sn;
	    }

	}
    }
//If ICOMPQ = 0, use DLASDQ to compute the singular values.
    if (icompq == 0) {
	Rlasdq("U", 0, n, 0, 0, 0, &d[0], &e[0], &vt[0], ldvt, &u[0], ldu, &u[0], ldu, &work[wstart], info);
	goto L40;
    }
//If N is smaller than the minimum divide size SMLSIZ, then solve
//the problem with another solver.
    if (n <= smlsiz) {
	if (icompq == 2) {
	    Rlaset("A", n, n, Zero, One, &u[0], ldu);
	    Rlaset("A", n, n, Zero, One, &vt[0], ldvt);
	    Rlasdq("U", 0, n, n, n, 0, &d[0], &e[0], &vt[0], ldvt, &u[0], ldu, &u[0], ldu, &work[wstart], info);
	} else if (icompq == 1) {
	    iu = 1;
	    ivt = iu + n;
	    Rlaset("A", n, n, Zero, One, &q[iu + (qstart - 1) * n], n);
	    Rlaset("A", n, n, Zero, One, &q[ivt + (qstart - 1) * n], n);
	    Rlasdq("U", 0, n, n, n, 0, &d[0], &e[0], &q[ivt + (qstart - 1) * n], n, &q[iu + (qstart - 1) * n], n, &q[iu + (qstart - 1) * n], n, &work[wstart], info);
	}
	goto L40;
    }
    if (icompq == 2) {
	Rlaset("A", n, n, Zero, One, &u[0], ldu);
	Rlaset("A", n, n, Zero, One, &vt[0], ldvt);
    }
//Scale.
    orgnrm = Rlanst("M", n, &d[0], &e[0]);
    if (orgnrm == Zero) {
	return;
    }
    Rlascl("G", 0, 0, orgnrm, One, n, 1, &d[0], n, &ierr);
    Rlascl("G", 0, 0, orgnrm, One, nm1, 1, &e[0], nm1, &ierr);
    eps = Rlamch("Epsilon");
    mlvl = (INTEGER) cast2double(log(double (n) / (double)(smlsiz + 1)) /log(Two)) + 1;
    smlszp = smlsiz + 1;
    if (icompq == 1) {
	iu = 1;
	ivt = smlsiz + 1;
	difl = ivt + smlszp;
	difr = difl + mlvl;
	z = difr + (mlvl << 1);
	ic = z + mlvl;
	is = ic + 1;
	poles = is + 1;
	givnum = poles + (mlvl << 1);
	k = 0;
	givptr = 2;
	perm = 3;
	givcol = perm + mlvl;
    }
    for (i = 0; i < n; i++) {
	if (abs(d[i]) < eps) {
	    d[i] = sign(eps, d[i]);
	}
    }
    start = 1;
    sqre = 0;
    for (i = 0; i < nm1; i++) {
	if (abs(e[i]) < eps || i == nm1) {
//Subproblem found. First determine its size and then
//apply divide and conquer on it.
	    if (i < nm1) {
//A subproblem with E(I) small for I < NM1
		nsize = i - start + 1;
	    } else if (abs(e[i]) >= eps) {
//A subproblem with E(NM1) not too small but I = NM1
		nsize = n - start + 1;
	    } else {
//A subproblem with E(NM1) small. This implies an
//1-by-1 subproblem at D(N). Solve this 1-by-1 problem
//first.
		nsize = i - start + 1;
		if (icompq == 2) {
		    u[n + n * ldu] = sign(One, d[n]);
		    vt[n + n * ldvt] = One;
		} else if (icompq == 1) {
		    q[n + (qstart - 1) * n] = sign(One, d[n]);
		    q[n + (smlsiz + qstart - 1) * n] = One;
		}
		d[n] = abs(d[n]);
	    }
	    if (icompq == 2) {
		Rlasd0(nsize, sqre, &d[start], &e[start], &u[start + start * ldu], ldu, &vt[start + start * ldvt], ldvt, smlsiz, &iwork[1], &work[wstart], info);
	    } else {
		Rlasda(icompq, smlsiz, nsize, sqre, &d[start], &e[start], &q[start + (iu + qstart - 2) * n], n,
		       &q[start + (ivt + qstart - 2) * n], &iq[start + k * n], &q[start + (difl + qstart - 2) * n],
		       &q[start + (difr + qstart - 2) * n], &q[start + (z + qstart - 2) * n],
		       &q[start + (poles + qstart - 2) * n], &iq[start + givptr * n], &iq[start + givcol * n], n,
		       &iq[start + perm * n], &q[start + (givnum + qstart - 2) * n], &q[start + (ic + qstart - 2) * n],
		       &q[start + (is + qstart - 2) * n], &work[wstart], &iwork[1], info);
		if (*info != 0) {
		    return;
		}
	    }
	    start = i + 1;
	}
    }
//Unscale
    Rlascl("G", 0, 0, One, orgnrm, n, 1, &d[0], n, &ierr);
  L40:
//Use Selection Sort to minimize swaps of singular vectors
    for (ii = 1; ii <= n; ii++) {
	i = ii - 1;
	kk = i;
	p = d[i];
	for (j = ii; j <= n; j++) {
	    if (d[j] > p) {
		kk = j;
		p = d[j];
	    }
	}
	if (kk != i) {
	    d[kk] = d[i];
	    d[i] = p;
	    if (icompq == 1) {
		iq[i] = kk;
	    } else if (icompq == 2) {
		Rswap(n, &u[i * ldu + 1], 1, &u[kk * ldu + 1], 1);
		Rswap(n, &vt[i + ldvt], ldvt, &vt[kk + ldvt], ldvt);
	    }
	} else if (icompq == 1) {
	    iq[i] = i;
	}
    }
//If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO
    if (icompq == 1) {
	if (iuplo == 1) {
	    iq[n] = 1;
	} else {
	    iq[n] = 0;
	}
    }
//If B is lower bidiagonal, update U by those Givens rotations
//which rotated B to be upper bidiagonal
    if (iuplo == 2 && icompq == 2) {
	Rlasr("L", "V", "B", n, n, &work[0], &work[n], &u[0], ldu);
    }
    return;
}
