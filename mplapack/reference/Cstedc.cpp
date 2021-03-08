/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cstedc.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cstedc(const char *compz, INTEGER n, REAL * d, REAL * e, COMPLEX * z, INTEGER ldz, COMPLEX * work, INTEGER lwork, REAL * rwork,
	    INTEGER lrwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER i, j, k, m;
    REAL p;
    INTEGER ii, ll, lgn;
    REAL eps, tiny;
    INTEGER lwmin, start;
    INTEGER finish;
    INTEGER liwmin, icompz;
    REAL orgnrm;
    INTEGER lrwmin;
    INTEGER lquery;
    INTEGER smlsiz;
    REAL Two = 2.0, One = 1.0, Zero = 0.0;
//Test the input parameters.
    *info = 0;
    lquery = lwork == -1 || lrwork == -1 || liwork == -1;
    if (Mlsame(compz, "N")) {
	icompz = 0;
    } else if (Mlsame(compz, "V")) {
	icompz = 1;
    } else if (Mlsame(compz, "I")) {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (icompz < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (ldz < 1 || (icompz > 0 && ldz < max((INTEGER) 1, n))) {
	*info = -6;
    }

    if (*info == 0) {
//Compute the workspace requirements
	smlsiz = iMlaenv(9, "Cstedc", " ", 0, 0, 0, 0);
	if (n <= 1 || icompz == 0) {
	    lwmin = 1;
	    liwmin = 1;
	    lrwmin = 1;
	} else if (n <= smlsiz) {
	    lwmin = 1;
	    liwmin = 1;
	    lrwmin = (n - 1) * 2;
	} else if (icompz == 1) {
	    lgn = (INTEGER) cast2double(log((double) (n)) / log(Two));
	    if ((2 ^ (lgn)) < n) {
		lgn++;
	    }
	    if ((2 ^ (lgn)) < n) {
		lgn++;
	    }
	    lwmin = n * n;
	    lrwmin = n * 3 + 1 + (n * 2) * lgn + n * n * 3;
	    liwmin = n * 6 + 6 + n * 5 * lgn;
	} else if (icompz == 2) {
	    lwmin = 1;
	    lrwmin = (n * 4) + 1 + (n * n * 2);
	    liwmin = n * 5 + 3;
	}
	work[1] = lwmin;
	rwork[1] = (double) lrwmin;
	iwork[1] = liwmin;

	if (lwork < lwmin && !lquery) {
	    *info = -8;
	} else if (lrwork < lrwmin && !lquery) {
	    *info = -10;
	} else if (liwork < liwmin && !lquery) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	Mxerbla("Cstedc", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (n == 1) {
	if (icompz != 0) {
	    z[ldz + 1] = One;
	}
	return;
    }
//If the following conditional clause is removed, then the routine
//will use the Divide and Conquer routine to compute only the
//eigenvalues, which requires (3N + 3N**2) real workspace and
//(2 + 5N + 2N lg(N)) INTEGER workspace.
//Since on many architectures DSTERF is much faster than any other
//algorithm for finding eigenvalues only, it is used here
//as the default. If the conditional clause is removed, then
//information on the size of workspace needs to be changed.
//If COMPZ = 'N', use DSTERF to compute the eigenvalues.
    if (icompz == 0) {
	Rsterf(n, &d[0], &e[0], info);
	goto L70;
    }
//If N is smaller than the minimum divide size (SMLSIZ+1), then
//solve the problem with another solver.
    if (n <= smlsiz) {
	Csteqr(compz, n, &d[0], &e[0], &z[0], ldz, &rwork[1], info);
    } else {
//If COMPZ = 'I', we simply call DSTEDC instead.
	if (icompz == 2) {
	    Rlaset("Full", n, n, Zero, One, &rwork[1], n);
	    ll = n * n + 1;
	    Rstedc("I", n, &d[0], &e[0], &rwork[1], n, &rwork[ll], lrwork - ll + 1, &iwork[1], liwork, info);
	    for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) {
		    z[i + j * ldz] = rwork[(j - 1) * n + i];
		}
	    }
	    goto L70;
	}
//From now on, only option left to be handled is COMPZ = 'V',
//i.e. ICOMPZ = One
//Scale.
	orgnrm = Rlanst("M", n, &d[0], &e[0]);
	if (orgnrm == Zero) {
	    goto L70;
	}
	eps = Rlamch("Epsilon");
	start = 1;
      L30:
	if (start <= n) {
//Let FINISH be the position of the next subdiagonal entry
//such that E( FINISH ) <= TINY or FINISH = N if no such
//subdiagonal exists.  The matrix identified by the elements
//between START and FINISH constitutes an independent
//sub-problem.
	    finish = start;
	  L40:
	    if (finish < n) {
		tiny = eps * sqrt(abs(d[finish])) * sqrt(abs(d[finish + 1]));
		if (abs(e[finish]) > tiny) {
		    ++finish;
		    goto L40;
		}
	    }
//(Sub) Problem determined.  Compute its size and solve it.
	    m = finish - start + 1;
	    if (m > smlsiz) {
//Scale.
		orgnrm = Rlanst("M", m, &d[start], &e[start]);
		Rlascl("G", 0, 0, orgnrm, One, m, 1, &d[start], m, info);
		Rlascl("G", 0, 0, orgnrm, One, m - 1, 1, &e[start], m - 1, info);
		Claed0(n, m, &d[start], &e[start], &z[start * ldz + 1], ldz, &work[0], n, &rwork[1], &iwork[1], info);
		if (*info > 0) {
		    *info = (*info / (m + 1) + start - 1) * (n + 1) + *info % (m + 1) + start - 1;
		    goto L70;
		}
//Scale back.
		Rlascl("G", 0, 0, One, orgnrm, m, 1, &d[start], m, info);
	    } else {
		Rsteqr("I", m, &d[start], &e[start], &rwork[1], m, &rwork[m * m + 1], info);
		Clacrm(n, m, &z[start * ldz + 1], ldz, &rwork[1], m, &work[0], n, &rwork[m * m + 1]);
		Clacpy("A", n, m, &work[0], n, &z[start * ldz + 1], ldz);
		if (*info > 0) {
		    *info = start * (n + 1) + finish;
		    goto L70;
		}
	    }
	    start = finish + 1;
	    goto L30;
	}
//endwhile
//If the problem split any number of times, then the eigenvalues
//will not be properly ordered.  Here we permute the eigenvalues
//(and the associated eigenvectors) INTEGERo ascending order.
	if (m != n) {
//Use Selection Sort to minimize swaps of eigenvectors
	    for (ii = 1; ii <= n; ii++) {
		i = ii - 1;
		k = i;
		p = d[i];
		for (j = ii; j <= n; j++) {
		    if (d[j] < p) {
			k = j;
			p = d[j];
		    }
		}
		if (k != i) {
		    d[k] = d[i];
		    d[i] = p;
		    Cswap(n, &z[i * ldz + 1], 1, &z[k * ldz + 1], 1);
		}
	    }
	}
    }

  L70:
    work[1] = lwmin;
    rwork[1] = (double) lrwmin;
    iwork[1] = liwmin;
    return;
}
