/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctgevc.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctgevc(const char *side, const char *howmny, LOGICAL * select,
	    INTEGER n, COMPLEX * s, INTEGER lds, COMPLEX * p, INTEGER ldp, COMPLEX * vl, INTEGER ldvl, COMPLEX * vr, INTEGER ldvr,
	    INTEGER mm, INTEGER * m, COMPLEX * work, REAL * rwork, INTEGER * info)
{
    COMPLEX d;
    INTEGER i, j;
    COMPLEX ca, cb;
    INTEGER je, im, jr;
    REAL big;
    LOGICAL lsa, lsb;
    REAL ulp;
    COMPLEX sum;
    INTEGER ibeg, ieig, iend;
    REAL dmin;
    INTEGER isrc;
    REAL temp;
    COMPLEX suma, sumb;
    REAL xmax, scale;
    LOGICAL ilall = 0;
    INTEGER iside;
    REAL sbeta;
    REAL small;
    LOGICAL compll;
    REAL anorm, bnorm;
    LOGICAL compr;
    LOGICAL ilbbad;
    REAL acoefa, bcoefa, acoeff;
    COMPLEX bcoeff;
    LOGICAL ilback = 0;
    REAL ascale, bscale;
    COMPLEX salpha;
    REAL safmin;
    REAL bignum;
    LOGICAL ilcomp;
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
	Mxerbla("Ctgevc", -(*info));
	return;
    }
//Count the number of eigenvectors
    if (!ilall) {
	im = 0;
	for (j = 0; j < n; j++) {
	    if (select[j]) {
		im++;
	    }
	}
    } else {
	im = n;
    }
//Check diagonal of B
    ilbbad = MFALSE;
    for (j = 0; j < n; j++) {
	if (p[j + j * ldp].imag() != Zero) {
	    ilbbad = MTRUE;
	}
    }
    if (ilbbad) {
	*info = -7;
    } else if ((compll && ldvl < n) || ldvl < 1) {
	*info = -10;
    } else if ((compr && ldvr < n) || ldvr < 1) {
	*info = -12;
    } else if (mm < im) {
	*info = -13;
    }
    if (*info != 0) {
	Mxerbla("Ctgevc", -(*info));
	return;
    }
//Quick return if possible
    (*m) = im;
    if (n == 0) {
	return;
    }
//Machine Constants
    safmin = Rlamch("Safe minimum");
    big = One / safmin;
    //Rlabad(&safmin, &big);
    ulp = Rlamch("Epsilon") * Rlamch("Base");
    small = safmin * n / ulp;
    big = One / small;
    bignum = One / (safmin * n);
//Compute the 1-norm of each column of the strictly upper triangular
//part of A and B to check for possible overflow in the triangular
//solver.
    anorm = Cabs1(s[lds + 1]);
    bnorm = Cabs1(p[ldp + 1]);
    rwork[1] = Zero;
    rwork[n + 1] = Zero;
    for (j = 2; j <= n; j++) {
	rwork[j] = Zero;
	rwork[n + j] = Zero;
	for (i = 0; i < j - 1; i++) {
	    rwork[j] = rwork[j] + Cabs1(s[i + j * lds]);
	    rwork[n + j] = rwork[n + j] + Cabs1(p[i + j * ldp]);
	}
	mtemp1 = anorm, mtemp2 = rwork[j] + Cabs1(s[j + j * lds]);
	anorm = max(mtemp1, mtemp2);
	mtemp1 = bnorm, mtemp2 = rwork[n + j] + Cabs1(p[j + j * ldp]);
	bnorm = max(mtemp1, mtemp2);
    }
    ascale = One / max(anorm, safmin);
    bscale = One / max(bnorm, safmin);
//Left eigenvectors
    if (compll) {
	ieig = 0;
//Main loop over eigenvalues
	for (je = 1; je <= n; je++) {
	    if (ilall) {
		ilcomp = MTRUE;
	    } else {
		ilcomp = select[je];
	    }
	    if (ilcomp) {
		ieig++;
		if (Cabs1(s[je + je * lds]) <= safmin && abs(p[je + je * ldp].real()) <= safmin) {
//Singular matrix pencil -- return unit eigenvector
		    for (jr = 1; jr <= n; jr++) {
			vl[jr + ieig * ldvl] = Zero;
		    }
		    vl[ieig + ieig * ldvl] = One;
		    goto L140;
		}
//Non-singular eigenvalue:
//Compute coefficients  a  and  b  in
//     H
//   y  ( a A - b B ) = 0
		mtemp1 = Cabs1(s[je + je * lds]) * ascale, mtemp2 = Cabs1(p[je + je * ldp]) * bscale;
		mtemp3 = max(mtemp1, mtemp2);
		temp = One / max(mtemp3, safmin);
		salpha = temp * s[je + je * lds] * ascale;
		sbeta = temp * p[je + je * ldp].real() * bscale;
		acoeff = sbeta * ascale;
		bcoeff = salpha * bscale;
//Scale to avoid underflow
		lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
		lsb = Cabs1(salpha) >= safmin && Cabs1(bcoeff) < small;
		scale = One;
		if (lsa) {
		    scale = small / abs(sbeta) * min(anorm, big);
		}
		if (lsb) {
		    mtemp1 = scale, mtemp2 = small / Cabs1(salpha) * min(bnorm, big);
		    scale = max(mtemp1, mtemp2);
		}
		if (lsa || lsb) {
		    mtemp3 = One, mtemp4 = abs(acoeff), mtemp1 = max(mtemp3, mtemp4), mtemp2 = Cabs1(bcoeff);
		    mtemp3 = scale, mtemp4 = One / (safmin * max(mtemp1, mtemp2));
		    scale = min(mtemp3, mtemp4);
		    if (lsa) {
			acoeff = ascale * (scale * sbeta);
		    } else {
			acoeff = scale * acoeff;
		    }
		    if (lsb) {
			bcoeff = bscale * (scale * salpha);
		    } else {
			bcoeff = scale * bcoeff;
		    }
		}
		acoefa = abs(acoeff);
		bcoefa = Cabs1(bcoeff);
		xmax = One;
		for (jr = 1; jr <= n; jr++) {
		    work[jr] = Zero;
		}
		work[je] = One;
		mtemp1 = ulp * acoefa * anorm, mtemp2 = ulp * bcoefa * bnorm;
		mtemp3 = max(mtemp1, mtemp2);
		dmin = max(mtemp3, safmin);
//                                H
//Triangular solve of  (a A - b B)  y = 0
//                        H
//(rowwise in  (a A - b B) , or columnwise in a A - b B)
		for (j = je + 1; j <= n; j++) {
//   Compute
//         j-1
//   SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
//         k=je
//   (Scale if necessary)
		    temp = One / xmax;
		    if (acoefa * rwork[j] + bcoefa * rwork[n + j] > bignum * temp) {
			for (jr = je; jr <= j - 1; jr++) {
			    work[jr] = temp * work[jr];
			}
			xmax = One;
		    }
		    suma = Zero;
		    sumb = Zero;
		    for (jr = je; jr <= j - 1; jr++) {
			suma = suma + conj(s[jr + j * lds]) * work[jr];
			sumb = sumb + conj(p[jr + j * ldp]) * work[jr];
		    }
		    sum = acoeff * suma - conj(bcoeff) * sumb;
//Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) )
//with scaling and perturbation of the denominator
		    d = conj(acoeff * s[j + j * lds] - bcoeff * p[j + j * ldp]);
		    if (Cabs1(d) <= dmin) {
			d = dmin;
		    }
		    if (Cabs1(d) < One) {
			if (Cabs1(sum) >= bignum * Cabs1(d)) {
			    temp = One / Cabs1(sum);
			    for (jr = je; jr <= j - 1; jr++) {
				work[jr] = temp * work[jr];
			    }
			    xmax = temp * xmax;
			    sum = temp * sum;
			}
		    }
		    work[j] = Cladiv(-sum, d);
		    mtemp1 = xmax, mtemp2 = Cabs1(work[j]);
		    xmax = max(mtemp1, mtemp2);
		}
//Back transform eigenvector if HOWMNY='B'.
		if (ilback) {
		    Cgemv("N", n, n + 1 - je, One, &vl[je * ldvl + 1], ldvl, &work[je], 1, Zero, &work[n + 1], 1);
		    isrc = 2;
		    ibeg = 1;
		} else {
		    isrc = 1;
		    ibeg = je;
		}
//Copy and scale eigenvector into column of VL
		xmax = Zero;
		for (jr = ibeg; jr <= n; jr++) {
		    mtemp1 = xmax, mtemp2 = Cabs1(work[(isrc - 1) * n + jr]);
		    xmax = max(mtemp1, mtemp2);
		}
		if (xmax > safmin) {
		    temp = One / xmax;
		    for (jr = ibeg; jr <= n; jr++) {
			vl[jr + ieig * ldvl] = temp * work[(isrc - 1) * n + jr];
		    }
		} else {
		    ibeg = n + 1;
		}
		for (jr = 1; jr <= ibeg - 1; jr++) {
		    vl[jr + ieig * ldvl] = Zero;
		}
	    }
	  L140:
	    ;
	}
    }
//Right eigenvectors
    if (compr) {
	ieig = im + 1;
//Main loop over eigenvalues
	for (je = n; je >= 1; je--) {
	    if (ilall) {
		ilcomp = MTRUE;
	    } else {
		ilcomp = select[je];
	    }
	    if (ilcomp) {
		ieig--;
		if (Cabs1(s[je + je * lds]) <= safmin && abs(p[je + je * ldp].real()) <= safmin) {
//Singular matrix pencil -- return unit eigenvector
		    for (jr = 1; jr <= n; jr++) {
			vr[jr + ieig * ldvr] = Zero;
		    }
		    vr[ieig + ieig * ldvr] = One;
		    goto L250;
		}
//Non-singular eigenvalue:
//Compute coefficients  a  and  b  in
//( a A - b B ) x  = 0
		mtemp1 = Cabs1(s[je + je * lds]), mtemp2 = abs(p[je + je * ldp].real()) * bscale;
		mtemp3 = max(mtemp1, mtemp2);
		temp = One / max(mtemp3, safmin);
		salpha = temp * s[je + je * lds] * ascale;
		sbeta = temp * p[je + je * ldp].real() * bscale;
		acoeff = sbeta * ascale;
		bcoeff = salpha * bscale;
//Scale to avoid underflow
		lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
		lsb = Cabs1(salpha) >= safmin && Cabs1(bcoeff) < small;
		scale = One;
		if (lsa) {
		    scale = small / abs(sbeta) * min(anorm, big);
		}
		if (lsb) {
		    mtemp1 = scale, mtemp2 = small / Cabs1(salpha) * min(bnorm, big);
		    scale = max(mtemp1, mtemp2);
		}
		if (lsa || lsb) {
		    mtemp1 = One, mtemp2 = abs(acoeff);
		    mtemp3 = max(mtemp1, mtemp2), mtemp4 = Cabs1(bcoeff);
		    mtemp1 = scale, mtemp2 = One / (safmin * max(mtemp3, mtemp4));
		    scale = min(mtemp1, mtemp2);
		    if (lsa) {
			acoeff = ascale * (scale * sbeta);
		    } else {
			acoeff = scale * acoeff;
		    }
		    if (lsb) {
			bcoeff = bscale * (scale * salpha);
		    } else {
			bcoeff = scale * bcoeff;
		    }
		}
		acoefa = abs(acoeff);
		bcoefa = Cabs1(bcoeff);
		xmax = One;
		for (jr = 1; jr <= n; jr++) {
		    work[jr] = Zero;
		}
		work[je] = One;
		mtemp1 = ulp * acoefa * anorm, mtemp2 = ulp * bcoefa * bnorm;
		mtemp3 = max(mtemp1, mtemp2);
		dmin = max(mtemp3, safmin);
//Triangular solve of  (a A - b B) x = 0  (columnwise)
//WORK(1:j-1) contains sums w,
//WORK(j+1:JE) contains x
		for (jr = 1; jr <= je - 1; jr++) {
		    work[jr] = acoeff * s[jr + je * lds] - bcoeff * p[jr + je * ldp];
		}
		work[je] = One;
		for (j = je - 1; j >= 1; j--) {
//Form x(j) := - w(j) / d
//with scaling and perturbation of the denominator
		    d = acoeff * s[j + j * lds] - bcoeff * p[j + j * ldp];
		    if (Cabs1(d) <= dmin) {
			d = dmin;
		    }
		    if (Cabs1(d) < One) {
			if (Cabs1(work[j]) >= bignum * Cabs1(d)) {
			    temp = One / Cabs1(work[j]);
			    for (jr = 1; jr <= je; jr++) {
				work[jr] = temp * work[jr];
			    }
			}
		    }
		    work[j] = Cladiv(-work[j], d);
		    if (j > 1) {
//w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
			if (Cabs1(work[j]) > One) {
			    temp = One / Cabs1(work[j]);
			    if (acoefa * rwork[j] + bcoefa * rwork[n + j] >= bignum * temp) {
				for (jr = 1; jr <= je; jr++) {
				    work[jr] = temp * work[jr];
				}
			    }
			}
			ca = acoeff * work[j];
			cb = bcoeff * work[j];
			for (jr = 1; jr <= j - 1; jr++) {
			    work[jr] = work[jr] + ca * s[jr + j * lds] - cb * p[jr + j * ldp];
			}
		    }
		}
//Back transform eigenvector if HOWMNY='B'.
		if (ilback) {
		    Cgemv("N", n, je, One, &vr[0], ldvr, &work[0], 1, Zero, &work[n + 1], 1);
		    isrc = 2;
		    iend = n;
		} else {
		    isrc = 1;
		    iend = je;
		}
//Copy and scale eigenvector into column of VR
		xmax = Zero;
		for (jr = 1; jr <= iend; jr++) {
		    mtemp1 = xmax, mtemp2 = Cabs1(work[(isrc - 1) * n + jr]);
		    xmax = max(mtemp1, mtemp2);
		}
		if (xmax > safmin) {
		    temp = One / xmax;
		    for (jr = 1; jr <= iend; jr++) {
			vr[jr + ieig * ldvr] = temp * work[(isrc - 1) * n + jr];
		    }
		} else {
		    iend = 0;
		}
		for (jr = iend + 1; jr <= n; jr++) {
		    vr[jr + ieig * ldvr] = Zero;
		}
	    }
	  L250:
	    ;
	}
    }
    return;
}
