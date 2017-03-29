/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtgevc.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Rtgevc(const char *side, const char *howmny, LOGICAL * select,
	    INTEGER n, REAL * s, INTEGER lds, REAL * p, INTEGER ldp, REAL * vl, INTEGER ldvl, REAL * vr, INTEGER ldvr, INTEGER mm, INTEGER * m, REAL * work, INTEGER * info)
{
    INTEGER i, j, ja, jc, je, na, im, jr, jw, nw;
    REAL big;
    LOGICAL lsa, lsb;
    REAL ulp, sum[4];
    INTEGER ibeg, ieig, iend;
    REAL dmin, temp, xmax, sums[4], sump[4];
    REAL cim2a, cim2b, cre2a, cre2b, temp2, bdiag[2] = {0.0, 0.0}, acoef, scale;
    LOGICAL ilall;
    INTEGER iside;
    REAL sbeta;
    LOGICAL il2by2;
    INTEGER iinfo;
    REAL small;
    REAL anorm, bnorm;
    REAL temp2i;
    REAL temp2r;
    LOGICAL compll, compr;
    LOGICAL ilabad, ilbbad;
    REAL acoefa, bcoefa, cimaga, cimagb;
    LOGICAL ilback = 0;
    REAL bcoefi, ascale, bscale, creala, crealb;
    REAL bcoefr, salfar, safmin;
    REAL xscale, bignum;
    LOGICAL ilcomp, ilcplx;
    INTEGER ihwmny;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2, mtemp3, mtemp4;

//Decode and Test the input parameters
    if (Mlsame(howmny, "A")) {
	ihwmny = 1;
	ilall = MTRUE;
	ilback = MFALSE;
    } else if (Mlsame(howmny, "S")) {
	ihwmny = 2;
	ilall = MFALSE;
	ilback = MFALSE;
    } else if (Mlsame(howmny, "B")) {
	ihwmny = 3;
	ilall = MTRUE;
	ilback = MTRUE;
    } else {
	ihwmny = -1;
	ilall = MTRUE;
    }
    if (Mlsame(side, "R")) {
	iside = 1;
	compll = MFALSE;
	compr = MTRUE;
    } else if (Mlsame(side, "L")) {
	iside = 2;
	compll = MTRUE;
	compr = MFALSE;
    } else if (Mlsame(side, "B")) {
	iside = 3;
	compll = MTRUE;
	compr = MTRUE;
    } else {
	iside = -1;
    }
    *info = 0;
    if (iside < 0) {
	*info = -1;
    } else if (ihwmny < 0) {
	*info = -2;
    } else if (n < 0) {
	*info = -4;
    } else if (lds < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldp < max((INTEGER) 1, n)) {
	*info = -8;
    }
    if (*info != 0) {
	Mxerbla("Rtgevc", -(*info));
	return;
    }
//Count the number of eigenvectors to be computed
    if (!ilall) {
	im = 0;
	ilcplx = MFALSE;
	for (j = 0; j < n; j++) {
	    if (ilcplx) {
		ilcplx = MFALSE;
		goto L10;
	    }
	    if (j < n) {
		if (s[j + 1 + j * lds] != Zero) {
		    ilcplx = MTRUE;
		}
	    }
	    if (ilcplx) {
		if (select[j] || select[j + 1]) {
		    im = im + 2;
		}
	    } else {
		if (select[j]) {
		    im++;
		}
	    }
	  L10:
	    ;
	}
    } else {
	im = n;
    }
//Check 2-by-2 diagonal blocks of A, B
    ilabad = MFALSE;
    ilbbad = MFALSE;
    for (j = 0; j < n - 1; j++) {
	if (s[j + 1 + j * lds] != Zero) {
	    if (p[j + j * ldp] == Zero || p[j + 1 + (j + 1) * ldp] == Zero || p[j + (j + 1) * ldp] != Zero) {
		ilbbad = MTRUE;
	    }
	    if (j < n - 1) {
		if (s[j + 2 + (j + 1) * lds] != Zero) {
		    ilabad = MTRUE;
		}
	    }
	}
    }
    if (ilabad) {
	*info = -5;
    } else if (ilbbad) {
	*info = -7;
    } else if (compll && (ldvl < n || ldvl < 1)) {
	*info = -10;
    } else if (compr && (ldvr < n || ldvr < 1)) {
	*info = -12;
    } else if (mm < im) {
	*info = -13;
    }
    if (*info != 0) {
	Mxerbla("Rtgevc", -(*info));
	return;
    }
//Quick return if possible
    *m = im;
    if (n == 0) {
	return;
    }
//Machine Constants
    safmin = Rlamch("Safe minimum");
    big = One / safmin;
    // Rlabad(&safmin, &big);
    ulp = Rlamch("Epsilon") * Rlamch("Base");
    small = safmin * n / ulp;
    big = One / small;
    bignum = One / (safmin * n);
//Compute the 1-norm of each column of the strictly upper triangular
//part (i.e., excluding all elements belonging to the diagonal
//blocks) of A and B to check for possible overflow in the
//triangular solver.
    anorm = abs(s[lds + 1]);
    if (n > 1) {
	anorm = anorm + abs(s[lds + 2]);
    }
    bnorm = abs(p[ldp + 1]);
    work[1] = Zero;
    work[n + 1] = Zero;
    for (j = 2; j <= n; j++) {
	temp = Zero;
	temp2 = Zero;
	if (s[j + (j - 1) * lds] == Zero) {
	    iend = j - 1;
	} else {
	    iend = j - 2;
	}
	for (i = 0; i < iend; i++) {
	    temp = temp + abs(s[i + j * lds]);
	    temp2 = temp2 + abs(p[i + j * ldp]);
	}
	work[j] = temp;
	work[n + j] = temp2;
	for (i = iend + 1; i <= min(j + 1, n); i++) {
	    temp = temp + abs(s[i + j * lds]);
	    temp2 = temp2 + abs(p[i + j * ldp]);
	}
	anorm = max(anorm, temp);
	bnorm = max(bnorm, temp2);
    }
    ascale = One / max(anorm, safmin);
    bscale = One / max(bnorm, safmin);
//Left eigenvectors
    if (compll) {
	ieig = 0;
//Main loop over eigenvalues
	ilcplx = MFALSE;
	for (je = 1; je <= n; je++) {
//Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
//(b) this would be the second of a complex pair.
//Check for complex eigenvalue, so as to be sure of which
//entry(-ies) of SELECT to look at.
	    if (ilcplx) {
		ilcplx = MFALSE;
		goto L220;
	    }
	    nw = 1;
	    if (je < n) {
		if (s[je + 1 + je * lds] != Zero) {
		    ilcplx = MTRUE;
		    nw = 2;
		}
	    }
	    if (ilall) {
		ilcomp = MTRUE;
	    } else if (ilcplx) {
		ilcomp = select[je] || select[je + 1];
	    } else {
		ilcomp = select[je];
	    }
	    if (!ilcomp) {
		goto L220;
	    }
//Decide if (a) singular pencil, (b) real eigenvalue, or
//(c) complex eigenvalue.
	    if (!ilcplx) {
		if (abs(s[je + je * lds]) <= safmin && abs(p[je + je * ldp]) <= safmin) {
//Singular matrix pencil -- return unit eigenvector
		    ieig++;
		    for (jr = 1; jr <= n; jr++) {
			vl[jr + ieig * ldvl] = Zero;
		    }
		    vl[ieig + ieig * ldvl] = One;
		    goto L220;
		}
	    }
//Clear vector
	    for (jr = 1; jr <= nw * n; jr++) {
		work[(n << 1) + jr] = Zero;
	    }
//                                      T
//Compute coefficients in  ( a A - b B )  y = 0
//   a  is  ACOEF
//   b  is  BCOEFR + i*BCOEFI
	    if (!ilcplx) {
//Real eigenvalue
		mtemp1 = abs(s[je + je * lds]) * ascale, mtemp2 = abs(p[je + je * ldp]) * bscale, mtemp3 = max(mtemp1, mtemp2);
		temp = One / max(mtemp1, safmin);
		salfar = temp * s[je + je * lds] * ascale;
		sbeta = temp * p[je + je * ldp] * bscale;
		acoef = sbeta * ascale;
		bcoefr = salfar * bscale;
		bcoefi = Zero;
//Scale to avoid underflow
		scale = One;
		lsa = abs(sbeta) >= safmin && abs(acoef) < small;
		lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
		if (lsa) {
		    scale = small / abs(sbeta) * min(anorm, big);
		}
		if (lsb) {
		    mtemp1 = scale, mtemp2 = small / abs(salfar) * min(bnorm, big);
		    scale = max(mtemp1, mtemp2);
		}
		if (lsa || lsb) {
		    mtemp1 = One, mtemp2 = abs(acoef), mtemp3 = max(mtemp1, mtemp2), mtemp4 = abs(bcoefr);
		    mtemp1 = scale, mtemp2 = One / (safmin * max(mtemp3, mtemp4));
		    scale = min(mtemp1, mtemp2);
		    if (lsa) {
			acoef = ascale * (scale * sbeta);
		    } else {
			acoef = scale * acoef;
		    }
		    if (lsb) {
			bcoefr = bscale * (scale * salfar);
		    } else {
			bcoefr = scale * bcoefr;
		    }
		}
		acoefa = abs(acoef);
		bcoefa = abs(bcoefr);
//First component is 1
		work[(n << 1) + je] = One;
		xmax = One;
	    } else {
//Compllex eigenvalue
		temp = safmin * 100;
		Rlag2(&s[je + je * lds], lds, &p[je + je * ldp], ldp, temp, &acoef, &temp, &bcoefr, &temp2, &bcoefi);
		bcoefi = -bcoefi;
		if (bcoefi == Zero) {
		    *info = je;
		    return;
		}
//Scale to avoid over/underflow
		acoefa = abs(acoef);
		bcoefa = abs(bcoefr) + abs(bcoefi);
		scale = One;
		if (acoefa * ulp < safmin && acoefa >= safmin) {
		    scale = safmin / ulp / acoefa;
		}
		if (bcoefa * ulp < safmin && bcoefa >= safmin) {
		    mtemp1 = scale, mtemp2 = safmin / ulp / bcoefa;
		    scale = max(mtemp1, mtemp2);
		}
		if (safmin * acoefa > ascale) {
		    scale = ascale / (safmin * acoefa);
		}
		if (safmin * bcoefa > bscale) {
		    mtemp1 = scale, mtemp2 = bscale / (safmin * bcoefa);
		    scale = min(mtemp1, mtemp2);
		}
		if (scale != One) {
		    acoef = scale * acoef;
		    acoefa = abs(acoef);
		    bcoefr = scale * bcoefr;
		    bcoefi = scale * bcoefi;
		    bcoefa = abs(bcoefr) + abs(bcoefi);
		}
//Compute first two components of eigenvector
		temp = acoef * s[je + 1 + je * lds];
		temp2r = acoef * s[je + je * lds] - bcoefr * p[je + je * ldp];
		temp2i = -bcoefi * p[je + je * ldp];
		if (abs(temp) > abs(temp2r) + abs(temp2i)) {
		    work[n * 2 + je] = One;
		    work[n * 3 + je] = Zero;
		    work[n * 2 + je + 1] = -temp2r / temp;
		    work[n * 3 + je + 1] = -temp2i / temp;
		} else {
		    work[n * 2 + je + 1] = One;
		    work[n * 3 + je + 1] = Zero;
		    temp = acoef * s[je + (je + 1) * lds];
		    work[(n << 1) + je] = (bcoefr * p[je + 1 + (je + 1) * ldp] - acoef * s[je + 1 + (je + 1) * lds]) / temp;
		    work[n * 3 + je] = bcoefi * p[je + 1 + (je + 1) * ldp] / temp;
		}
		mtemp1 = abs(work[n * 2 + je]) + abs(work[n * 3 + je]);
		mtemp2 = abs(work[n * 2 + je + 1]) + abs(work[n * 3 + je + 1]);
		xmax = max(mtemp1, mtemp2);
	    }
	    mtemp1 = ulp * acoefa * anorm, mtemp2 = ulp * bcoefa * bnorm;
	    mtemp3 = max(mtemp1, mtemp2);
	    dmin = max(mtemp3, safmin);
//                                T
//Triangular solve of  (a A - b B)  y = 0
//                        T
//(rowwise in  (a A - b B) , or columnwise in (a A - b B) )
	    il2by2 = MFALSE;
	    for (j = je + nw; j <= n; j++) {
		if (il2by2) {
		    il2by2 = MFALSE;
		    goto L160;
		}
		na = 1;
		bdiag[0] = p[j + j * ldp];
		if (j < n) {
		    if (s[j + 1 + j * lds] != Zero) {
			il2by2 = MTRUE;
			bdiag[1] = p[j + 1 + (j + 1) * ldp];
			na = 2;
		    }
		}
//Check whether scaling is necessary for dot products
		xscale = One / max(One, xmax);
		mtemp1 = work[j], mtemp2 = work[n + j];
		mtemp3 = max(mtemp1, mtemp2), mtemp4 = acoefa * work[j] + bcoefa * work[n + j];
		temp = max(mtemp3, mtemp4);
		if (il2by2) {
		    mtemp1 = temp, mtemp2 = work[j + 1];
		    mtemp3 = max(mtemp1, mtemp2), mtemp4 = work[n + j + 1], mtemp1 = max(mtemp3, mtemp4), mtemp2 = acoefa * work[j + 1] + bcoefa * work[n + j + 1];
		    temp = max(mtemp1, mtemp2);
		}
		if (temp > bignum * xscale) {
		    for (jw = 0; jw <= nw - 1; jw++) {
			for (jr = je; jr <= j - 1; jr++) {
			    work[(jw + 2) * n + jr] = xscale * work[(jw + 2) * n + jr];
			}
		    }
		    xmax = xmax * xscale;
		}
//Compute dot products
//      j-1
//SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
//      k=je
//To reduce the op count, this is done as
//_        j-1                  _        j-1
//a*conjg( sum  S(k,j)*x(k) ) - b*conjg( sum  P(k,j)*x(k) )
//         k=je                          k=je
//which may cause underflow problems if A or B are close
//to underflow.  (E.g., less than SMALL.)
//A series of compiler directives to defeat vectorization
//for the next loop
		for (jw = 1; jw <= nw; jw++) {
		    for (ja = 1; ja <= na; ja++) {
			sums[ja + (jw << 1) - 3] = Zero;
			sump[ja + (jw << 1) - 3] = Zero;
			for (jr = je; jr <= j - 1; jr++) {
			    sums[ja + (jw << 1) - 3] = sums[ja + (jw << 1) - 3] + s[jr + (j + ja - 1) * lds] * work[(jw + 1) * n + jr];
			    sump[ja + (jw << 1) - 3] = sump[ja + (jw << 1) - 3] + p[jr + (j + ja - 1) * ldp] * work[(jw + 1) * n + jr];
			}
		    }
		}
		for (ja = 1; ja <= na; ja++) {
		    if (ilcplx) {
			sum[ja - 1] = -acoef * sums[ja - 1] + bcoefr * sump[ja - 1] - bcoefi * sump[ja + 1];
			sum[ja + 1] = -acoef * sums[ja + 1] + bcoefr * sump[ja + 1] + bcoefi * sump[ja - 1];
		    } else {
			sum[ja - 1] = -acoef * sums[ja - 1] + bcoefr * sump[ja - 1];
		    }
		}
//                    T
//Solve  ( a A - b B )  y = SUM(,)
//with scaling and perturbation of the denominator
		Rlaln2(MTRUE, na, nw, dmin, acoef, &s[j + j * lds]
		       , lds, bdiag[0], bdiag[1], sum, 2, bcoefr, bcoefi, &work[(n << 1) + j], n, &scale, &temp, &iinfo);
		if (scale < One) {
		    for (jw = 0; jw <= nw - 1; jw++) {
			for (jr = je; jr <= j - 1; jr++) {
			    work[(jw + 2) * n + jr] = scale * work[(jw + 2) * n + jr];
			}
		    }
		    xmax = scale * xmax;
		}
		xmax = max(xmax, temp);
	      L160:
		;
	    }
//Copy eigenvector to VL, back transforming if
//HOWMNY='B'.
	    ieig++;
	    if (ilback) {
		for (jw = 0; jw <= nw - 1; jw++) {
		    Rgemv("N", n, n + 1 - je, One, &vl[je * ldvl + 1], ldvl, &work[(jw + 2) * n + je], 1, Zero, &work[(jw + 4) * n + 1], 1);
		}
		Rlacpy(" ", n, nw, &work[(n << 2) + 1], n, &vl[je * ldvl + 1], ldvl);
		ibeg = 1;
	    } else {
		Rlacpy(" ", n, nw, &work[(n << 1) + 1], n, &vl[ieig * ldvl + 1], ldvl);
		ibeg = je;
	    }
//Scale eigenvector
	    xmax = Zero;
	    if (ilcplx) {
		for (j = ibeg; j <= n; j++) {
		    mtemp1 = xmax, mtemp2 = abs(vl[j + ieig * ldvl]) + abs(vl[j + (ieig + 1) * ldvl]);
		    xmax = max(mtemp1, mtemp2);
		}
	    } else {
		for (j = ibeg; j <= n; j++) {
		    mtemp1 = xmax, mtemp2 = abs(vl[j + ieig * ldvl]);
		    xmax = max(mtemp1, mtemp2);
		}
	    }
	    if (xmax > safmin) {
		xscale = One / xmax;
		for (jw = 0; jw <= nw - 1; jw++) {
		    for (jr = ibeg; jr <= n; jr++) {
			vl[jr + (ieig + jw) * ldvl] = xscale * vl[jr + (ieig + jw) * ldvl];
		    }
		}
	    }
	    ieig = ieig + nw - 1;
	  L220:
	    ;
	}
    }
//Right eigenvectors
    if (compr) {
	ieig = im + 1;
//Main loop over eigenvalues
	ilcplx = MFALSE;
	for (je = n; je >= 1; je--) {
//Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
//(b) this would be the second of a compllex pair.
//Check for compllex eigenvalue, so as to be sure of which
//entry(-ies) of SELECT to look at -- if compllex, SELECT(JE)
//or SELECT(JE-1).
//If this is a compllex pair, the 2-by-2 diagonal block
//corresponding to the eigenvalue is in rows/columns JE-1:JE
	    if (ilcplx) {
		ilcplx = MFALSE;
		goto L500;
	    }
	    nw = 1;
	    if (je > 1) {
		if (s[je + (je - 1) * lds] != Zero) {
		    ilcplx = MTRUE;
		    nw = 2;
		}
	    }
	    if (ilall) {
		ilcomp = MTRUE;
	    } else if (ilcplx) {
		ilcomp = select[je] || select[je - 1];
	    } else {
		ilcomp = select[je];
	    }
	    if (!ilcomp) {
		goto L500;
	    }
//Decide if (a) singular pencil, (b) real eigenvalue, or
//(c) compllex eigenvalue.
	    if (!ilcplx) {
		if (abs(s[je + je * lds]) <= safmin && abs(p[je + je * ldp]) <= safmin) {
//Singular matrix pencil -- unit eigenvector
		    ieig--;
		    for (jr = 1; jr <= n; jr++) {
			vr[jr + ieig * ldvr] = Zero;
		    }
		    vr[ieig + ieig * ldvr] = One;
		    goto L500;
		}
	    }
//Clear vector
	    for (jw = 0; jw <= nw - 1; jw++) {
		for (jr = 1; jr <= n; jr++) {
		    work[(jw + 2) * n + jr] = Zero;
		}
	    }
//Compute coefficients in  ( a A - b B ) x = 0
//   a  is  ACOEF
//   b  is  BCOEFR + i*BCOEFI
	    if (!ilcplx) {
//Real eigenvalue
		mtemp1 = abs(s[je + je * lds]) * ascale, mtemp2 = abs(p[je + je * ldp]) * bscale;
		mtemp3 = max(mtemp1, mtemp2);
		temp = One / max(mtemp3, safmin);
		salfar = temp * s[je + je * lds] * ascale;
		sbeta = temp * p[je + je * ldp] * bscale;
		acoef = sbeta * ascale;
		bcoefr = salfar * bscale;
		bcoefi = Zero;
//Scale to avoid underflow
		scale = One;
		lsa = abs(sbeta) >= safmin && abs(acoef) < small;
		lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
		if (lsa) {
		    scale = small / abs(sbeta) * min(anorm, big);
		}
		if (lsb) {
		    mtemp1 = scale, mtemp2 = small / abs(salfar) * min(bnorm, big);
		    scale = max(mtemp1, mtemp2);
		}
		if (lsa || lsb) {
		    mtemp1 = One, mtemp2 = abs(acoef), mtemp3 = max(mtemp1, mtemp2), mtemp4 = abs(bcoefr);
		    mtemp1 = scale, mtemp2 = One / (safmin * max(mtemp3, mtemp4));
		    scale = min(mtemp1, mtemp2);
		    if (lsa) {
			acoef = ascale * (scale * sbeta);
		    } else {
			acoef = scale * acoef;
		    }
		    if (lsb) {
			bcoefr = bscale * (scale * salfar);
		    } else {
			bcoefr = scale * bcoefr;
		    }
		}
		acoefa = abs(acoef);
		bcoefa = abs(bcoefr);
//First component is 1
		work[(n << 1) + je] = One;
		xmax = One;
//Compute contribution from column JE of A and B to sum
//(See "Further Details", above.)
		for (jr = 1; jr <= je - 1; jr++) {
		    work[(n << 1) + jr] = bcoefr * p[jr + je * ldp] - acoef * s[jr + je * lds];
		}
	    } else {
//Complex eigenvalue
		mtemp1 = safmin * 100;
		Rlag2(&s[je - 1 + (je - 1) * lds], lds, &p[je - 1 + (je - 1) * ldp], ldp, mtemp1, &acoef, &temp, &bcoefr, &temp2, &bcoefi);
		if (bcoefi == Zero) {
		    *info = je - 1;
		    return;
		}
//Scale to avoid over/underflow
		acoefa = abs(acoef);
		bcoefa = abs(bcoefr) + abs(bcoefi);
		scale = One;
		if (acoefa * ulp < safmin && acoefa >= safmin) {
		    scale = safmin / ulp / acoefa;
		}
		if (bcoefa * ulp < safmin && bcoefa >= safmin) {
		    mtemp1 = scale, mtemp2 = safmin / ulp / bcoefa;
		    scale = max(mtemp1, mtemp2);
		}
		if (safmin * acoefa > ascale) {
		    scale = ascale / (safmin * acoefa);
		}
		if (safmin * bcoefa > bscale) {
		    mtemp1 = scale, mtemp2 = bscale / (safmin * bcoefa);
		    scale = min(mtemp1, mtemp2);
		}
		if (scale != One) {
		    acoef = scale * acoef;
		    acoefa = abs(acoef);
		    bcoefr = scale * bcoefr;
		    bcoefi = scale * bcoefi;
		    bcoefa = abs(bcoefr) + abs(bcoefi);
		}
//Compute first two components of eigenvector
//and contribution to sums
		temp = acoef * s[je + (je - 1) * lds];
		temp2r = acoef * s[je + je * lds] - bcoefr * p[je + je * ldp];
		temp2i = -bcoefi * p[je + je * ldp];
		if (abs(temp) >= abs(temp2r) + abs(temp2i)) {
		    work[n * 2 + je] = One;
		    work[n * 3 + je] = Zero;
		    work[n * 4 + je - 1] = -temp2r / temp;
		    work[n * 3 + je - 1] = -temp2i / temp;
		} else {
		    work[n * 2 + je - 1] = One;
		    work[n * 3 + je - 1] = Zero;
		    temp = acoef * s[je - 1 + je * lds];
		    work[(n << 1) + je] = (bcoefr * p[je - 1 + (je - 1) * ldp] - acoef * s[je - 1 + (je - 1) * lds]) / temp;
		    work[n * 3 + je] = bcoefi * p[je - 1 + (je - 1) * ldp] / temp;
		}

		mtemp1 = abs(work[(n << 1) + je]) + abs(work[n * 3 + je]);
		mtemp2 = abs(work[(n << 1) + je - 1]) + abs(work[n * 3 + je - 1]);
		xmax = max(mtemp1, mtemp2);
//Compute contribution from columns JE and JE-1
//of A and B to the sums.
		creala = acoef * work[(n << 1) + je - 1];
		cimaga = acoef * work[n * 3 + je - 1];
		crealb = bcoefr * work[(n << 1) + je - 1] - bcoefi * work[n * 3 + je - 1];
		cimagb = bcoefi * work[(n << 1) + je - 1] + bcoefr * work[n * 3 + je - 1];
		cre2a = acoef * work[(n << 1) + je];
		cim2a = acoef * work[n * 3 + je];
		cre2b = bcoefr * work[(n << 1) + je] - bcoefi * work[n * 3 + je];
		cim2b = bcoefi * work[(n << 1) + je] + bcoefr * work[n * 3 + je];
		for (jr = 1; jr <= je - 2; jr++) {
		    work[(n << 1) + jr] = -creala * s[jr + (je - 1) * lds]
			+ crealb * p[jr + (je - 1) * ldp] - cre2a * s[jr + je * lds] + cre2b * p[jr + je * ldp];
		    work[n * 3 + jr] = -cimaga * s[jr + (je - 1) * lds] + cimagb * p[jr + (je - 1) * ldp] - cim2a * s[jr + je * lds] + cim2b * p[jr + je * ldp];
		}
	    }
	    mtemp1 = ulp * acoefa * anorm, mtemp2 = ulp * bcoefa * bnorm;
	    mtemp3 = max(mtemp1, mtemp2);
	    dmin = max(mtemp3, safmin);
//Columnwise triangular solve of  (a A - b B)  x = 0
	    il2by2 = MFALSE;
	    for (j = je - nw; j >= 1; j--) {
//If a 2-by-2 block, is in position j-1:j, wait until
//next iteration to process it (when it will be j:j+1)
		if (!il2by2 && j > 1) {
		    if (s[j + (j - 1) * lds] != Zero) {
			il2by2 = MTRUE;
			goto L370;
		    }
		}
		bdiag[0] = p[j + j * ldp];
		if (il2by2) {
		    na = 2;
		    bdiag[1] = p[j + 1 + (j + 1) * ldp];
		} else {
		    na = 1;
		}
//Compute x(j) (and x(j+1), if 2-by-2 block)
		Rlaln2(MFALSE, na, nw, dmin, acoef, &s[j + j * lds], lds, bdiag[0], bdiag[1], &work[(n << 1) + j], n, bcoefr, bcoefi, sum, 2, &scale, &temp, &iinfo);
		if (scale < One) {
		    for (jw = 0; jw <= nw - 1; jw++) {
			for (jr = 1; jr <= je; jr++) {
			    work[(jw + 2) * n + jr] = scale * work[(jw + 2) * n + jr];
			}
		    }
		}
		mtemp1 = scale * xmax;
		xmax = max(mtemp1, temp);
		for (jw = 1; jw <= nw; jw++) {
		    for (ja = 1; ja <= na; ja++) {
			work[(jw + 1) * n + j + ja - 1] = sum[ja + (jw << 1) - 3];
		    }
		}
//w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
		if (j > 1) {
//Check whether scaling is necessary for sum.
		    xscale = One / max(One, xmax);
		    temp = acoefa * work[j] + bcoefa * work[n + j];
		    if (il2by2) {
			mtemp1 = temp, mtemp2 = acoefa * work[j + 1] + bcoefa * work[n + j + 1];
			temp = max(mtemp1, mtemp2);
		    }
		    mtemp1 = max(temp, acoefa);
		    temp = max(mtemp1, bcoefa);
		    if (temp > bignum * xscale) {
			for (jw = 0; jw <= nw - 1; jw++) {
			    for (jr = 1; jr <= je; jr++) {
				work[(jw + 2) * n + jr] = xscale * work[(jw + 2) * n + jr];
			    }
			}
			xmax = xmax * xscale;
		    }
//Compute the contributions of the off-diagonals of
//column j (and j+1, if 2-by-2 block) of A and B to the
//sums.
		    for (ja = 1; ja <= na; ja++) {
			if (ilcplx) {
			    creala = acoef * work[(n << 1) + j + ja - 1];
			    cimaga = acoef * work[n * 3 + j + ja - 1];
			    crealb = bcoefr * work[(n << 1) + j + ja - 1] - bcoefi * work[n * 3 + j + ja - 1];
			    cimagb = bcoefi * work[(n << 1) + j + ja - 1] + bcoefr * work[n * 3 + j + ja - 1];
			    for (jr = 1; jr <= j - 1; jr++) {
				work[(n << 1) + jr] = work[(n << 1) + jr] - creala * s[jr + (j + ja - 1) * lds]
				    + crealb * p[jr + (j + ja - 1) * ldp];
				work[n * 3 + jr] = work[n * 3 + jr] - cimaga * s[jr + (j + ja - 1) * lds]
				    + cimagb * p[jr + (j + ja - 1) * ldp];
			    }
			} else {
			    creala = acoef * work[(n << 1) + j + ja - 1];
			    crealb = bcoefr * work[(n << 1) + j + ja - 1];
			    for (jr = 1; jr <= j - 1; jr++) {
				work[(n << 1) + jr] = work[(n << 1) + jr] - creala * s[jr + (j + ja - 1) * lds]
				    + crealb * p[jr + (j + ja - 1) * ldp];
			    }
			}
		    }
		}
		il2by2 = MFALSE;
	      L370:
		;
	    }
//Copy eigenvector to VR, back transforming if
//HOWMNY='B'.
	    ieig = ieig - nw;
	    if (ilback) {
		for (jw = 0; jw <= nw - 1; jw++) {
		    for (jr = 1; jr <= n; jr++) {
			work[(jw + 4) * n + jr] = work[(jw + 2) * n + 1] * vr[jr + ldvr];
		    }
//A series of compiler directives to defeat
//vectorization for the next loop
		    for (jc = 2; jc <= je; jc++) {
			for (jr = 1; jr <= n; jr++) {
			    work[(jw + 4) * n + jr] = work[(jw + 4) * n + jr] + work[(jw + 2) * n + jc] * vr[jr + jc * ldvr];
			}
		    }
		}
		for (jw = 0; jw <= nw - 1; jw++) {
		    for (jr = 1; jr <= n; jr++) {
			vr[jr + (ieig + jw) * ldvr] = work[(jw + 4) * n + jr];
		    }
		}
		iend = n;
	    } else {
		for (jw = 0; jw <= nw - 1; jw++) {
		    for (jr = 1; jr <= n; jr++) {
			vr[jr + (ieig + jw) * ldvr] = work[(jw + 2) * n + jr];
		    }
		}
		iend = je;
	    }
//Scale eigenvector
	    xmax = Zero;
	    if (ilcplx) {
		for (j = 0; j < iend; j++) {
		    mtemp1 = xmax, mtemp2 = abs(vr[j + ieig * ldvr]) + abs(vr[j + (ieig + 1) * ldvr]);
		    xmax = max(mtemp1, mtemp2);
		}
	    } else {
		for (j = 0; j < iend; j++) {
		    mtemp1 = xmax, mtemp2 = abs(vr[j + ieig * ldvr]);
		    xmax = max(mtemp1, mtemp2);
		}
	    }
	    if (xmax > safmin) {
		xscale = One / xmax;
		for (jw = 0; jw <= nw - 1; jw++) {
		    for (jr = 1; jr <= iend; jr++) {
			vr[jr + (ieig + jw) * ldvr] = xscale * vr[jr + (ieig + jw) * ldvr];
		    }
		}
	    }
	  L500:
	    ;
	}
    }
    return;
}
