/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rstedc.cpp,v 1.15 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rstedc(const char *compz, INTEGER n, REAL * d, REAL * e, REAL * z, INTEGER ldz, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER i, j, k, m;
    REAL p;
    INTEGER ii, lgn;
    REAL eps, tiny;
    INTEGER lwmin;
    INTEGER start;
    INTEGER finish;
    INTEGER liwmin, icompz;
    REAL orgnrm;
    INTEGER lquery;
    INTEGER smlsiz, storez, strtrw;
    REAL One = 1.0, Zero = 0.0;

//Test the input parameters.
    *info = 0;
    lquery = 0;
    if (lwork == -1 || liwork == -1)
	lquery = 1;
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
	smlsiz = iMlaenv(2, "Rstedc", " ", 0, 0, 0, 0);
	if (n <= 1 || icompz == 0) {
	    liwmin = 1;
	    lwmin = 1;
	} else if (n <= smlsiz) {
	    liwmin = 1;
	    lwmin = 2 * (n - 1);
	} else {
	    lgn = (INTEGER) cast2double(log2(double (n)));
	    if ((lgn ^ 2) < n) {
		lgn++;
	    }
	    if ((lgn ^ 2) < n) {
		lgn++;
	    }
	    if (icompz == 1) {
		lwmin = n * 3 + 1 + (n * 2) * lgn + n * n * 3;
		liwmin = n * 6 + 6 + n * 5 * lgn;
	    } else if (icompz == 2) {
		lwmin = (n * 4) + 1 + n * n;
		liwmin = n * 5 + 3;
	    }
	}
	work[0] = (double) lwmin;
	iwork[0] = liwmin;
	if (lwork < lwmin && !lquery) {
	    *info = -8;
	} else if (liwork < liwmin && !lquery) {
	    *info = -10;
	}
    }
    if (*info != 0) {
	Mxerbla("Rstedc", -(*info));
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
	Rsterf(n, d, e, info);
	goto L50;
    }
//If N is smaller than the minimum divide size (SMLSIZ+1), then
//solve the problem with another solver.
    if (n <= smlsiz) {
	Rsteqr(compz, n, d, e, z, ldz, work, info);
    } else {
//If COMPZ = 'V', the Z matrix must be stored elsewhere for later
//use.
	if (icompz == 1) {
	    storez = n * n + 1;
	} else {
	    storez = 1;
	}
	if (icompz == 2) {
	    Rlaset("Full", n, n, Zero, One, z, ldz);
	}
//Scale.
	orgnrm = Rlanst("M", n, d, e);
	if (orgnrm == Zero) {
	    goto L50;
	}
	eps = Rlamch("Epsilon");
	start = 1;
//while ( START <= N )
      L10:
	if (start <= n) {
//Let FINISH be the position of the next subdiagonal entry
//such that E( FINISH ) <= TINY or FINISH = N if no such
//subdiagonal exists.  The matrix identified by the elements
//between START and FINISH constitutes an independent
//sub-problem.
	    finish = start;
	  L20:
	    if (finish < n) {
		tiny = eps * sqrt(abs(d[finish])) * sqrt(abs(d[finish + 1]));
		if (abs(e[finish]) > tiny) {
		    ++finish;
		    goto L20;
		}
	    }
//(Sub) Problem determined.  Compute its size and solve it.
	    m = finish - start + 1;
	    if (m == 1) {
		start = finish + 1;
		goto L10;
	    }
	    if (m > smlsiz) {
//Scale.
		orgnrm = Rlanst("M", m, &d[start], &e[start]);
		Rlascl("G", 0, 0, orgnrm, One, m, 1, &d[start], m, info);
		Rlascl("G", 0, 0, orgnrm, One, m - 1, 1, &e[start], m - 1, info);
		if (icompz == 1) {
		    strtrw = 1;
		} else {
		    strtrw = start;
		}
		Rlaed0(icompz, n, m, &d[start], &e[start], &z[strtrw + start * ldz], ldz, &work[0], n, &work[storez], &iwork[1], info);
		if (*info != 0) {
		    *info = (*info / (m + 1) + start - 1) * (n + 1) + *info % (m + 1) + start - 1;
		    goto L50;
		}
//Scale back.
		Rlascl("G", 0, 0, One, orgnrm, m, 1, &d[start], m, info);
	    } else {
		if (icompz == 1) {
//Since QR won't update a Z matrix which is larger than
//the length of D, we must solve the sub-problem in a
//workspace and then multiply back into Z.
		    Rsteqr("I", m, &d[start], &e[start], &work[0], m, &work[m * m + 1], info);
		    Rlacpy("A", n, m, &z[start * ldz + 1], ldz, &work[storez], n);
		    Rgemm("N", "N", n, m, m, One, &work[storez], n, &work[0], m, Zero, &z[start * ldz + 1], ldz);
		} else if (icompz == 2) {
		    Rsteqr("I", m, &d[start], &e[start], &z[start + start * ldz], ldz, &work[0], info);
		} else {
		    Rsterf(m, &d[start], &e[start], info);
		}
		if (*info != 0) {
		    *info = start * (n + 1) + finish;
		    goto L50;
		}
	    }
	    start = finish + 1;
	    goto L10;
	}
//endwhile
//If the problem split any number of times, then the eigenvalues
//will not be properly ordered.  Here we permute the eigenvalues
//(and the associated eigenvectors) into ascending order.
	if (m != n) {
	    if (icompz == 0) {
//Use Quick Sort
		Rlasrt("I", n, &d[0], info);
	    } else {
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
			Rswap(n, &z[i * ldz + 1], 1, &z[k * ldz + 1], 1);
		    }
		}
	    }
	}
    }
  L50:
    work[1] = (double) lwmin;
    iwork[1] = liwmin;
    return;
}
