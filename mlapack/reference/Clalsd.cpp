/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clalsd.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Clalsd(const char *uplo, INTEGER smlsiz, INTEGER n, INTEGER nrhs, REAL * d,
       REAL * e, COMPLEX * B, INTEGER ldb, REAL rcond, INTEGER * rank, COMPLEX * work, REAL * rwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER c, i, j, k;
    REAL r;
    INTEGER s, u, z;
    REAL cs;
    INTEGER bx;
    REAL sn;
    INTEGER st, vt, nm1, st1;
    REAL eps;
    INTEGER iwk;
    REAL tol;
    INTEGER difl, difr;
    REAL rcnd;
    INTEGER jcol, irwb, perm, nsub, nlvl, sqre, bxst, jrow, irwu, jimag;
    INTEGER jreal, irwib, poles, sizei, irwrb, nsize;
    INTEGER irwvt, icmpq1, icmpq2;
    INTEGER givcol;
    REAL orgnrm;
    INTEGER givnum, givptr, nrwork, irwwrk, smlszp;
    REAL Zero = 0.0, One = 1.0;

    *info = 0;

    if (n < 0) {
	*info = -3;
    } else if (nrhs < 1) {
	*info = -4;
    } else if (ldb < 1 || ldb < n) {
	*info = -8;
    }
    if (*info != 0) {
	Mxerbla("Clalsd", -(*info));
	return;
    }

    eps = Rlamch("Epsilon");
//Set up the tolerance.
    if (rcond <= Zero || rcond >= One) {
	rcnd = eps;
    } else {
	rcnd = rcond;
    }

    *rank = 0;
//Quick return if possible.
    if (n == 0) {
	return;
    } else if (n == 1) {
	if (d[1] == Zero) {
	    Claset("A", 1, nrhs, Zero, Zero, &B[0], ldb);
	} else {
	    *rank = 0;
	    Clascl("G", 0, 0, d[0], One, 1, nrhs, &B[0], ldb, info);
	    d[1] = abs(d[1]);
	}
	return;
    }
//Rotate the matrix if it is lower bidiagonal.
    if (Mlsame(uplo, "L")) {
	for (i = 0; i < n - 1; i++) {
	    Rlartg(d[i], e[i], &cs, &sn, &r);
	    d[i] = r;
	    e[i] = sn * d[i + 1];
	    d[i + 1] = cs * d[i + 1];
	    if (nrhs == 1) {
		CRrot(1, &B[i + ldb], 1, &B[i + 1 + ldb], 1, cs, sn);
	    } else {
		rwork[(i * 2) - 1] = cs;
		rwork[i * 2] = sn;
	    }

	}
	if (nrhs > 1) {
	    for (i = 0; i < nrhs; i++) {
		for (j = 0; j < n - 1; j++) {
		    cs = rwork[(j * 2) - 1];
		    sn = rwork[j * 2];
		    CRrot(1, &B[j + i * ldb], 1, &B[j + 1 + i * ldb], 1, cs, sn);
		}
	    }
	}
    }
//Scale.
    nm1 = n - 1;
    orgnrm = Rlanst("M", n, &d[0], &e[0]);
    if (orgnrm == Zero) {
	Claset("A", n, nrhs, Zero, Zero, &B[0], ldb);
	return;
    }

    Rlascl("G", 0, 0, orgnrm, One, n, 1, &d[0], n, info);
    Rlascl("G", 0, 0, orgnrm, One, nm1, 1, &e[0], nm1, info);
//If N is smaller than the minimum divide size SMLSIZ, then solve
//the problem with another solver.
    if (n <= smlsiz) {
	irwu = 1;
	irwvt = irwu + n * n;
	irwwrk = irwvt + n * n;
	irwrb = irwwrk;
	irwib = irwrb + n * nrhs;
	irwb = irwib + n * nrhs;
	Rlaset("A", n, n, Zero, One, &rwork[irwu], n);
	Rlaset("A", n, n, Zero, One, &rwork[irwvt], n);
	Rlasdq("U", 0, n, n, n, 0, &d[0], &e[0], &rwork[irwvt], n, &rwork[irwu], n, &rwork[irwwrk], 1, &rwork[irwwrk], info);
	if (*info != 0) {
	    return;
	}
//In the real version, B is passed to DLASDQ and multiplied
//internally by Q'. Here B is complex and that product is
//computed below in two steps (real and imaginary parts).
	j = irwb - 1;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = 1; jrow <= n; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].real();
	    }
	}
	Rgemm("T", "N", n, nrhs, n, One, &rwork[irwu], n, &rwork[irwb], n, Zero, &rwork[irwrb], n);
	j = irwb - 1;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = 1; jrow <= n; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].imag();
	    }
	}
	Rgemm("T", "N", n, nrhs, n, One, &rwork[irwu], n, &rwork[irwb], n, Zero, &rwork[irwib], n);
	jreal = irwrb - 1;
	jimag = irwib - 1;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = 1; jrow <= n; jrow++) {
		jreal++;
		jimag++;
		B[jrow + jcol * ldb] = rwork[jreal];
	    }
	}

	tol = rcnd * abs(d[iRamax(n, &d[0], 1)]);
	for (i = 0; i < n; i++) {
	    if (d[i] <= tol) {
		Claset("A", 1, nrhs, Zero, Zero, &B[i + ldb], ldb);
	    } else {
		Clascl("G", 0, 0, d[i], One, 1, nrhs, &B[i + ldb], ldb, info);
		++(*rank);
	    }
	}
//Since B is complex, the following call to DGEMM is performed
//in two steps (real and imaginary parts). That is for V * B
//(in the real version of the code V' is stored in WORK).
//CALL DGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO,
//            WORK( NWORK ), N )
	j = irwb - 1;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = 1; jrow <= n; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].real();
	    }
	}
	Rgemm("T", "N", n, nrhs, n, One, &rwork[irwvt], n, &rwork[irwb], n, Zero, &rwork[irwrb], n);
	j = irwb - 1;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = 1; jrow <= n; jrow++) {
		j++;
		rwork[j] = B[jrow + jcol * ldb].imag();
	    }
	}
	Rgemm("T", "N", n, nrhs, n, One, &rwork[irwvt], n, &rwork[irwb], n, Zero, &rwork[irwib], n);
	jreal = irwrb - 1;
	jimag = irwib - 1;
	for (jcol = 0; jcol <= nrhs; jcol++) {
	    for (jrow = 1; jrow <= n; jrow++) {
		jreal++;
		jimag++;
		B[jrow + jcol * ldb] = rwork[jreal];
	    }
	}
//Unscale.
	Rlascl("G", 0, 0, One, orgnrm, n, 1, &d[0], n, info);
	Rlasrt("D", n, &d[0], info);
	Clascl("G", 0, 0, orgnrm, One, n, nrhs, &B[0], ldb, info);
	return;
    }
//Book-keeping and setting up some constants.
    nlvl = (INTEGER) (log((double) n / (double) (smlsiz + 1)) / log((double) 2)) + 1;

    smlszp = smlsiz + 1;

    u = 1;
    vt = smlsiz * n + 1;
    difl = vt + smlszp * n;
    difr = difl + nlvl * n;
    z = difr + (nlvl * n * 2);
    c = z + nlvl * n;
    s = c + n;
    poles = s + n;
    givnum = poles + (nlvl * 2) * n;
    nrwork = givnum + (nlvl * 2) * n;
    bx = 1;

    irwrb = nrwork;
    irwib = irwrb + smlsiz * nrhs;
    irwb = irwib + smlsiz * nrhs;

    sizei = n + 1;
    k = sizei + n;
    givptr = k + n;
    perm = givptr + n;
    givcol = perm + nlvl * n;
    iwk = givcol + (nlvl * n * 2);

    st = 1;
    sqre = 0;
    icmpq1 = 1;
    icmpq2 = 0;
    nsub = 0;

    for (i = 0; i < n; i++) {
	if (abs(d[i]) < eps) {
	    d[i] = sign(eps, d[i]);
	}
    }

    for (i = 0; i < nm1; i++) {
	if (abs(e[i]) < eps || i == nm1) {
	    ++nsub;
	    iwork[nsub] = st;
//Subproblem found. First determine its size and then
//apply divide and conquer on it.
	    if (i < nm1) {
//A subproblem with E(I) small for I < NM1
		nsize = i - st + 1;
		iwork[sizei + nsub - 1] = nsize;
	    } else if (abs(e[i]) >= eps) {
//A subproblem with E(NM1) not too small but I = NM1
		nsize = n - st + 1;
		iwork[sizei + nsub - 1] = nsize;
	    } else {
//A subproblem with E(NM1) small. This implies an
//1-by-1 subproblem at D(N), which is not solved
//explicitly.
		nsize = i - st + 1;
		iwork[sizei + nsub - 1] = nsize;
		++nsub;
		iwork[nsub] = n;
		iwork[sizei + nsub - 1] = 1;
		Ccopy(nrhs, &B[n + ldb], ldb, &work[bx + nm1], n);
	    }
	    st1 = st - 1;
	    if (nsize == 1) {
//This is a 1-by-1 subproblem and is not solved
//explicitly.
		Ccopy(nrhs, &B[st + ldb], ldb, &work[bx + st1], n);
	    } else if (nsize <= smlsiz) {
//This is a small subproblem and is solved by DLASDQ.
		Rlaset("A", nsize, nsize, Zero, One, &rwork[vt + st1], n);
		Rlaset("A", nsize, nsize, Zero, One, &rwork[u + st1], n);
		Rlasdq("U", 0, nsize, nsize, nsize, 0, &d[st], &e[st], &rwork[vt + st1], n, &rwork[u + st1], n, &rwork[nrwork], 1, &rwork[nrwork], info);
		if (*info != 0) {
		    return;
		}
//In the real version, B is passed to DLASDQ and multiplied
//internally by Q'. Here B is complex and that product is
//computed below in two steps (real and imaginary parts).
		j = irwb - 1;
		for (jcol = 0; jcol <= nrhs; jcol++) {
		    for (jrow = st; jrow <= st + nsize - 1; jrow++) {
			j++;
			rwork[j] = B[jrow + jcol * ldb].real();
		    }
		}
		Rgemm("T", "N", nsize, nrhs, nsize, One, &rwork[u + st1]
		      , n, &rwork[irwb], nsize, Zero, &rwork[irwrb], nsize);
		j = irwb - 1;
		for (jcol = 0; jcol <= nrhs; jcol++) {
		    for (jrow = st; jrow <= st + nsize - 1; jrow++) {
			j++;
			rwork[j] = B[jrow + jcol * ldb].imag();
		    }
		}
		Rgemm("T", "N", nsize, nrhs, nsize, One, &rwork[u + st1]
		      , n, &rwork[irwb], nsize, Zero, &rwork[irwib], nsize);
		jreal = irwrb - 1;
		jimag = irwib - 1;
		for (jcol = 0; jcol <= nrhs; jcol++) {
		    for (jrow = st; jrow <= st + nsize - 1; jrow++) {
			jreal++;
			jimag++;
			B[jrow + jcol * ldb] = rwork[jreal];
		    }
		}
		Clacpy("A", nsize, nrhs, &B[st + ldb], ldb, &work[bx + st1], n);
	    } else {
//A large problem. Solve it using divide and conquer.
		Rlasda(icmpq1, smlsiz, nsize, sqre, &d[st], &e[st],
		       &rwork[u + st1], n, &rwork[vt + st1], &iwork[k + st1],
		       &rwork[difl + st1], &rwork[difr + st1], &rwork[z + st1],
		       &rwork[poles + st1], &iwork[givptr + st1],
		       &iwork[givcol + st1], n, &iwork[perm + st1], &rwork[givnum + st1], &rwork[c + st1], &rwork[s + st1], &rwork[nrwork], &iwork[iwk], info);
		if (*info != 0) {
		    return;
		}
		bxst = bx + st1;
		Clalsa(icmpq2, smlsiz, nsize, nrhs, &B[st + ldb], ldb, &work[bxst], n, &rwork[u + st1], n, &rwork[vt + st1], &iwork[k + st1], &rwork[difl + st1], &rwork[difr + st1]
		       , &rwork[z + st1], &rwork[poles + st1],
		       &iwork[givptr + st1], &iwork[givcol + st1], n,
		       &iwork[perm + st1], &rwork[givnum + st1], &rwork[c + st1], &rwork[s + st1], &rwork[nrwork], &iwork[iwk], info);
		if (*info != 0) {
		    return;
		}
	    }
	    st = i + 1;
	}
    }
//Apply the singular values and treat the tiny ones as zero.
    tol = rcnd * abs(d[iRamax(n, &d[0], 1)]);

    for (i = 0; i < n; i++) {
//Some of the elements in D can be negative because 1-by-1
//subproblems were not solved explicitly.
	if (abs(d[i]) <= tol) {
	    Claset("A", 1, nrhs, Zero, Zero, &work[bx + i - 1], n);
	} else {
	    ++(*rank);
	    Clascl("G", 0, 0, d[i], One, 1, nrhs, &work[bx + i - 1], n, info);
	}
	d[i] = abs(d[i]);
    }
//Now apply back the right singular vectors.
    icmpq2 = 1;
    for (i = 0; i < nsub; i++) {
	st = iwork[i];
	st1 = st - 1;
	nsize = iwork[sizei + i - 1];
	bxst = bx + st1;
	if (nsize == 1) {
	    Ccopy(nrhs, &work[bxst], n, &B[st + ldb], ldb);
	} else if (nsize <= smlsiz) {
//Since B and BX are complex, the following call to DGEMM
//is performed in two steps (real and imaginary parts).
//CALL DGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE,
//            RWORK( VT+ST1 ), N, RWORK( BXST ), N, ZERO,
//            B( ST, 1 ), LDB )
	    j = bxst - n - 1;
	    jreal = irwb - 1;
	    for (jcol = 0; jcol <= nrhs; jcol++) {
		j += n;
		for (jrow = 1; jrow <= nsize; jrow++) {
		    jreal++;
		    rwork[jreal] = work[j + jrow].real();
		}
	    }
	    Rgemm("T", "N", nsize, nrhs, nsize, One, &rwork[vt + st1], n, &rwork[irwb], nsize, Zero, &rwork[irwrb], nsize);
	    j = bxst - n - 1;
	    jimag = irwb - 1;
	    for (jcol = 0; jcol <= nrhs; jcol++) {
		j += n;
		for (jrow = 1; jrow <= nsize; jrow++) {
		    jimag++;
		    rwork[jimag] = work[j + jrow].imag();
		}
	    }
	    Rgemm("T", "N", nsize, nrhs, nsize, One, &rwork[vt + st1], n, &rwork[irwb], nsize, Zero, &rwork[irwib], nsize);
	    jreal = irwrb - 1;
	    jimag = irwib - 1;
	    for (jcol = 0; jcol <= nrhs; jcol++) {
		for (jrow = st; jrow <= st + nsize - 1; jrow++) {
		    jreal++;
		    jimag++;
		    B[jrow + jcol * ldb] = rwork[jreal];
		}
	    }
	} else {
	    Clalsa(icmpq2, smlsiz, nsize, nrhs, &work[bxst], n, &B[st +
								   ldb], ldb, &rwork[u + st1], n, &rwork[vt + st1],
		   &iwork[k + st1], &rwork[difl + st1], &rwork[difr + st1],
		   &rwork[z + st1], &rwork[poles + st1], &iwork[givptr + st1],
		   &iwork[givcol + st1], n, &iwork[perm + st1], &rwork[givnum + st1], &rwork[c + st1], &rwork[s + st1], &rwork[nrwork], &iwork[iwk], info);
	    if (*info != 0) {
		return;
	    }
	}

    }
//Unscale and sort the singular values.
    Rlascl("G", 0, 0, One, orgnrm, n, 1, &d[0], n, info);
    Rlasrt("D", n, &d[0], info);
    Clascl("G", 0, 0, orgnrm, One, n, nrhs, &B[0], ldb, info);
    return;
}
