/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgebal.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Cgebal(const char *job, INTEGER n, COMPLEX * A, INTEGER lda, INTEGER * ilo, INTEGER * ihi, REAL * scale, INTEGER * info)
{
    REAL c, f, g;
    INTEGER i, j, k, l, m;
    REAL r, s, ca, ra;
    INTEGER ica, ira, iexc;
    REAL sfmin1, sfmin2, sfmax1, sfmax2;
    INTEGER noconv;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;
    REAL mtemp1, mtemp2;

//Test the input parameters
    *info = 0;
    if (!Mlsame(job, "N") && !Mlsame(job, "P") && !Mlsame(job, "S")
	&& !Mlsame(job, "B")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Cgebal", -(*info));
	return;
    }
    k = 0;
    l = n;
    if (n == 0) {
	goto L210;
    }
    if (Mlsame(job, "N")) {
	for (i = 0; i < n; i++) {
	    scale[i] = One;
	}
	goto L210;
    }
    if (Mlsame(job, "S")) {
	goto L120;
    }
//Permutation to isolate eigenvalues if possible
    goto L50;
//Row and column exchange.
  L20:
    scale[m] = (double) j;
    if (j == m) {
	goto L30;
    }
    Cswap(l, &A[j * lda], 1, &A[m * lda], 1);
    Cswap(n - k + 1, &A[j + k * lda], lda, &A[m + k * lda], lda);
  L30:
    switch (iexc) {
    case 1:
	goto L40;
    case 2:
	goto L80;
    }
//Search for rows isolating an eigenvalue and push them down.
  L40:
    if (l == 1) {
	goto L210;
    }
    l--;
  L50:
    for (j = l; j >= 1; j--) {
	for (i = 0; i < l; i++) {
	    if (i == j) {
		goto L60;
	    }
	    if (A[j + i * lda].real() != Zero || A[j + i * lda].imag() != Zero) {
		goto L70;
	    }
	  L60:
	    ;
	}
	m = l;
	iexc = 1;
	goto L20;
      L70:
	;
    }
    goto L90;
//Search for columns isolating an eigenvalue and push them left.
  L80:
    k++;
  L90:
    for (j = k; j <= l; j++) {
	for (i = k; i <= l; i++) {
	    if (i == j) {
		goto L100;
	    }
	    if (A[i + j * lda].real() != Zero || A[i + j * lda].imag() != Zero) {
		goto L110;
	    }
	  L100:
	    ;
	}
	m = k;
	iexc = 2;
	goto L20;
      L110:
	;
    }
  L120:
    for (i = k; i <= l; i++) {
	scale[i] = One;
    }
    if (Mlsame(job, "P")) {
	goto L210;
    }
//Balance the submatrix in rows K to L.
//Iterative loop for norm reduction
    sfmin1 = Rlamch("S") / Rlamch("P");
    sfmax1 = One / sfmin1;
    sfmin2 = sfmin1 * Two;
    sfmax2 = One / sfmin2;
  L140:
    noconv = MFALSE;

    for (i = k; i <= l; i++) {
	c = Zero;
	r = Zero;
	for (j = k; j <= l; j++) {
	    if (j == i) {
		goto L150;
	    }
	    c = c + abs(A[j + i * lda].real()) + abs(A[j + i * lda].imag());
	    r = r + abs(A[i + j * lda].real()) + abs(A[i + j * lda].imag());
	  L150:
	    ;
	}
	ica = iCamax(l, &A[i * lda], 1);
	ca = abs(A[ica + i * lda]);
	ira = iCamax(n - k + 1, &A[i + k * lda], lda);
	ra = abs(A[i + (ira + k - 1) * lda]);
//Guard against zero C or R due to underflow.
	if (c == Zero || r == Zero) {
	    goto L200;
	}
	g = r / Two;
	f = One;
	s = c + r;
      L160:
	mtemp1 = max(f, c);
	mtemp2 = min(r, g);
	if (c >= g || max(mtemp1, ca) >= sfmax2 || min(mtemp2, ra) <= sfmin2) {
	    goto L170;
	}
	f = f * Two;
	c = c * Two;
	ca = ca * Two;
	r = r / Two;
	g = g / Two;
	ra = ra / Two;
	goto L160;

      L170:
	g = c / Two;
      L180:
	mtemp1 = min(f, c);
	mtemp2 = min(mtemp1, g);
	if (g < r || max(r, ra) >= sfmax2 || min(mtemp1, ca) <= sfmin2) {
	    goto L190;
	}
	f = f / Two;
	c = c / Two;
	g = g / Two;
	ca = ca / Two;
	r = r * Two;
	ra = ra * Two;
	goto L180;
//Now balance.
      L190:
	if (c + r >= s * .95) {
	    goto L200;
	}
	if (f < One && scale[i] < One) {
	    if (f * scale[i] <= sfmin1) {
		goto L200;
	    }
	}
	if (f > One && scale[i] > One) {
	    if (scale[i] >= sfmax1 / f) {
		goto L200;
	    }
	}
	g = One / f;
	scale[i] = scale[i] * f;
	noconv = MTRUE;
	CRscal(n - k + 1, g, &A[i + k * lda], lda);
	CRscal(l, f, &A[i * lda], 1);
      L200:
	;
    }
    if (noconv) {
	goto L140;
    }
  L210:
    *ilo = k;
    *ihi = l;
    return;
}
