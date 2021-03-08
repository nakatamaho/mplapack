/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chgeqz.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Chgeqz(const char *job, const char *compq, const char *compz, INTEGER n,
	    INTEGER ilo, INTEGER ihi, COMPLEX * h, INTEGER ldh,
	    COMPLEX * t, INTEGER ldt, COMPLEX * alpha, COMPLEX * beta, COMPLEX * q, INTEGER ldq, COMPLEX * z, INTEGER ldz,
	    COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER * info)
{
    REAL c;
    INTEGER j;
    COMPLEX s, t1;
    INTEGER jc, in;
    COMPLEX u12;
    INTEGER jr;
    COMPLEX ad11, ad12, ad21, ad22;
    INTEGER jch;
    INTEGER ilq = 0, ilz = 0;
    REAL ulp;
    COMPLEX abi22;
    REAL absb, atol, btol, temp;
    REAL temp2;
    COMPLEX ctemp;
    INTEGER iiter, ilast, jiter;
    REAL anorm, bnorm;
    INTEGER maxit;
    COMPLEX shift;
    REAL tempr;
    COMPLEX ctemp2, ctemp3;
    INTEGER ilazr2;
    REAL ascale, bscale;
    COMPLEX signbc;
    REAL safmin;
    COMPLEX eshift;
    INTEGER ilschr = 0;
    INTEGER icompq, ilastm;
    COMPLEX rtdisc;
    INTEGER ischur;
    INTEGER ilazro;
    INTEGER icompz, ifirst;
    INTEGER ifrstm;
    INTEGER istart;
    INTEGER lquery;
    REAL Zero = 0.0, Half = .5, One = 1.0;
    REAL mtemp1, mtemp2;

//Decode JOB, COMPQ, COMPZ
    if (Mlsame(job, "E")) {
	ilschr = MFALSE;
	ischur = 1;
    } else if (Mlsame(job, "S")) {
	ilschr = MTRUE;
	ischur = 2;
    } else {
	ischur = 0;
    }
    if (Mlsame(compq, "N")) {
	ilq = MFALSE;
	icompq = 1;
    } else if (Mlsame(compq, "V")) {
	ilq = MTRUE;
	icompq = 2;
    } else if (Mlsame(compq, "I")) {
	ilq = MTRUE;
	icompq = 3;
    } else {
	icompq = 0;
    }
    if (Mlsame(compz, "N")) {
	ilz = MFALSE;
	icompz = 1;
    } else if (Mlsame(compz, "V")) {
	ilz = MTRUE;
	icompz = 2;
    } else if (Mlsame(compz, "I")) {
	ilz = MTRUE;
	icompz = 3;
    } else {
	icompz = 0;
    }
//Check Argument Values
    *info = 0;
    work[1] = max((INTEGER) 1, n);
    lquery = lwork == -1;
    if (ischur == 0) {
	*info = -1;
    } else if (icompq == 0) {
	*info = -2;
    } else if (icompz == 0) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (ilo < 1) {
	*info = -5;
    } else if (ihi > n || ihi < ilo - 1) {
	*info = -6;
    } else if (ldh < n) {
	*info = -8;
    } else if (ldt < n) {
	*info = -10;
    } else if (ldq < 1 || (ilq && ldq < n)) {
	*info = -14;
    } else if (ldz < 1 || (ilz && ldz < n)) {
	*info = -16;
    } else if (lwork < max((INTEGER) 1, n) && !lquery) {
	*info = -18;
    }
    if (*info != 0) {
	Mxerbla("Chgeqz", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
//WORK( 1 ) = CMPLX( 1 )
    if (n <= 0) {
	work[1] = One;
	return;
    }
//Initialize Q and Z
    if (icompq == 3) {
	Claset("Full", n, n, Zero, One, &q[0], ldq);
    }
    if (icompz == 3) {
	Claset("Full", n, n, Zero, One, &z[0], ldz);
    }
//Machine Constants
    in = ihi + 1 - ilo;
    safmin = Rlamch("S");
    ulp = Rlamch("E") * Rlamch("B");
    anorm = Clanhs("F", in, &h[ilo + ilo * ldh], ldh, &rwork[1]);
    bnorm = Clanhs("F", in, &t[ilo + ilo * ldt], ldt, &rwork[1]);
    mtemp1 = safmin, mtemp2 = ulp * anorm;
    atol = max(mtemp1, mtemp2);
    mtemp1 = safmin, mtemp2 = ulp * bnorm;
    btol = max(mtemp1, mtemp2);
    ascale = One / max(safmin, anorm);
    bscale = One / max(safmin, bnorm);
//Set Eigenvalues IHI+1:N
    for (j = ihi + 1; j <= n; j++) {
	absb = abs(t[j + j * ldt]);
	if (absb > safmin) {
	    signbc = conj(t[j + j * ldt] / absb);
	    t[j + j * ldt] = absb;
	    if (ilschr) {
		Cscal(j - 1, signbc, &t[j * ldt + 1], 1);
		Cscal(j, signbc, &h[j * ldh + 1], 1);
	    } else {
		h[j + j * ldh] = h[j + j * ldh] * signbc;
	    }
	    if (ilz) {
		Cscal(n, signbc, &z[j * ldz + 1], 1);
	    }
	} else {
	    t[j + j * ldt] = Zero;
	}
	alpha[j] = h[j + j * ldh];
	beta[j] = t[j + j * ldt];
    }
//If IHI < ILO, skip QZ steps
    if (ihi < ilo) {
	goto L190;
    }
//MAIN QZ ITERATION LOOP
//Initialize dynamic indices
//Eigenvalues ILAST+1:N have been found.
//   Column operations modify rows IFRSTM:whatever
//   Row operations modify columns whatever:ILASTM
//If only eigenvalues are being computed, then
//   IFRSTM is the row of the last splitting row above row ILAST;
//   this is always at least ILO.
//IITER counts iterations since the last eigenvalue was found,
//   to tell when to use an extraordinary shift.
//MAXIT is the maximum number of QZ sweeps allowed.
    ilast = ihi;
    if (ilschr) {
	ifrstm = 1;
	ilastm = n;
    } else {
	ifrstm = ilo;
	ilastm = ihi;
    }
    iiter = 0;
    eshift = Zero;
    maxit = (ihi - ilo + 1) * 30;
    for (jiter = 1; jiter <= maxit; jiter++) {
//Check for too many iterations.
	if (jiter > maxit) {
	    goto L180;
	}
//Split the matrix if possible.
//Two tests:
//   1: H(j,j-1)=0  or  j=ILO
//   2: T(j,j)=0
//Special case: j=ILAST
	if (ilast == ilo) {
	    goto L60;
	} else {
	    if (Cabs1(h[ilast + (ilast - 1) * ldh]) <= atol) {
		h[ilast + (ilast - 1) * ldh] = Zero;
		goto L60;
	    }
	}
	if (abs(t[ilast + ilast * ldt]) <= btol) {
	    t[ilast + ilast * ldt] = Zero;
	    goto L50;
	}
//General case: j<ILAST
	for (j = ilast - 1; j >= ilo; j--) {
//Test 1: for H(j,j-1)=0 or j=ILO
	    if (j == ilo) {
		ilazro = MTRUE;
	    } else {
		if (Cabs1(h[j + (j - 1) * ldh]) <= atol) {
		    h[j + (j - 1) * ldh] = Zero;
		    ilazro = MTRUE;
		} else {
		    ilazro = MFALSE;
		}
	    }
//Test 2: for T(j,j)=0
	    if (abs(t[j + j * ldt]) < btol) {
		t[j + j * ldt] = Zero;
//Test 1a: Check for 2 consecutive small subdiagonals in A
		ilazr2 = MFALSE;
		if (!ilazro) {
		    if ((Cabs1(h[j + (j - 1) * ldh]) * ascale * Cabs1(h[j + 1 + j * ldh])) <= (Cabs1(h[j + j * ldh]) * (ascale * atol))) {
			ilazr2 = MTRUE;
		    }
		}
//If both tests pass (1 & 2), i.e., the leading diagonal
//element of B in the block is zero, split a 1x1 block off
//at the top. (I.e., at the J-th row/column) The leading
//diagonal element of the remainder can also be zero, so
//this may have to be done repeatedly.
		if (ilazro || ilazr2) {
		    for (jch = j; jch <= ilast - 1; jch++) {
			ctemp = h[jch + jch * ldh];
			Clartg(ctemp, h[jch + 1 + jch * ldh], &c, &s, &h[jch + jch * ldh]);
			h[jch + 1 + jch * ldh] = Zero;
			Crot(ilastm - jch, &h[jch + (jch + 1) * ldh], ldh, &h[jch + 1 + (jch + 1) * ldh], ldh, c, s);
			Crot(ilastm - jch, &t[jch + (jch + 1) * ldt], ldt, &t[jch + 1 + (jch + 1) * ldt], ldt, c, s);
			if (ilq) {
			    Crot(n, &q[jch * ldq + 1], 1, &q[(jch + 1) * ldq + 1], 1, c, conj(s));
			}
			if (ilazr2) {
			    h[jch + (jch - 1) * ldh] = h[jch + (jch - 1) * ldh] * c;
			}
			ilazr2 = MFALSE;
			if (Cabs1(t[jch + 1 + (jch + 1) * ldt]) >= btol) {
			    if (jch + 1 >= ilast) {
				goto L60;
			    } else {
				ifirst = jch + 1;
				goto L70;
			    }
			}
			t[jch + 1 + (jch + 1) * ldt] = Zero;
		    }
		    goto L50;
		} else {
//Only test 2 passed -- chase the zero to T(ILAST,ILAST)
//Then process as in the case T(ILAST,ILAST)=0
		    for (jch = j; jch <= ilast - 1; jch++) {
			ctemp = t[jch + (jch + 1) * ldt];
			Clartg(ctemp, t[jch + 1 + (jch + 1) * ldt], &c, &s, &t[jch + (jch + 1) * ldt]);
			t[jch + 1 + (jch + 1) * ldt] = Zero;
			if (jch < ilastm - 1) {
			    Crot(ilastm - jch - 1, &t[jch + (jch + 2) * ldt], ldt, &t[jch + 1 + (jch + 2) * ldt], ldt, c, s);
			}
			Crot(ilastm - jch + 2, &h[jch + (jch - 1) * ldh], ldh, &h[jch + 1 + (jch - 1) * ldh], ldh, c, s);
			if (ilq) {
			    Crot(n, &q[jch * ldq + 1], 1, &q[(jch + 1) * ldq + 1], 1, c, conj(s));
			}
			ctemp = h[jch + 1 + jch * ldh];
			Clartg(ctemp, h[jch + 1 + (jch - 1) * ldh], &c, &s, &h[jch + 1 + jch * ldh]);
			h[jch + 1 + (jch - 1) * ldh] = Zero;
			Crot(jch + 1 - ifrstm, &h[ifrstm + jch * ldh], 1, &h[ifrstm + (jch - 1) * ldh], 1, c, s);
			Crot(jch - ifrstm, &t[ifrstm + jch * ldt], 1, &t[ifrstm + (jch - 1) * ldt], 1, c, s);
			if (ilz) {
			    Crot(n, &z[jch * ldz + 1], 1, &z[(jch - 1) * ldz + 1], 1, c, s);
			}
		    }
		    goto L50;
		}
	    } else if (ilazro) {
//Only test 1 passed -- work on J:ILAST
		ifirst = j;
		goto L70;
	    }
//Neither test passed -- try next J
	}
//(Drop-through is "impossible")
	*info = (n << 1) + 1;
	goto L210;
//T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
//1x1 block.
      L50:
	ctemp = h[ilast + ilast * ldh];
	Clartg(ctemp, h[ilast + (ilast - 1) * ldh], &c, &s, &h[ilast + ilast * ldh]);
	h[ilast + (ilast - 1) * ldh] = Zero;
	Crot(ilast - ifrstm, &h[ifrstm + ilast * ldh], 1, &h[ifrstm + (ilast - 1) * ldh], 1, c, s);
	Crot(ilast - ifrstm, &t[ifrstm + ilast * ldt], 1, &t[ifrstm + (ilast - 1) * ldt], 1, c, s);
	if (ilz) {
	    Crot(n, &z[ilast * ldz + 1], 1, &z[(ilast - 1) * ldz + 1], 1, c, s);
	}
//H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA
      L60:
	absb = abs(t[ilast + ilast * ldt]);
	if (absb > safmin) {
	    signbc = conj(t[ilast + ilast * ldt] / absb);
	    t[ilast + ilast * ldt] = absb;
	    if (ilschr) {
		Cscal(ilast - ifrstm, signbc, &t[ifrstm + ilast * ldt], 1);
		Cscal(ilast + 1 - ifrstm, signbc, &h[ifrstm + ilast * ldh], 1);
	    } else {
		h[ilast + ilast * ldh] = h[ilast + ilast * ldh] * signbc;
	    }
	    if (ilz) {
		Cscal(n, signbc, &z[ilast * ldz + 1], 1);
	    }
	} else {
	    t[ilast + ilast * ldt] = Zero;
	}
	alpha[ilast] = h[ilast + ilast * ldh];
	beta[ilast] = t[ilast + ilast * ldt];
//Go to next block -- exit if finished.
	ilast--;
	if (ilast < ilo) {
	    goto L190;
	}
//Reset counters
	iiter = 0;
	eshift = Zero;
	if (!ilschr) {
	    ilastm = ilast;
	    if (ifrstm > ilast) {
		ifrstm = ilo;
	    }
	}
	goto L160;
//QZ step
//This iteration only involves rows/columns IFIRST:ILAST.  We
//assume IFIRST < ILAST, and that the diagonal of B is non-zero.
      L70:
	iiter++;
	if (!ilschr) {
	    ifrstm = ifirst;
	}
//Compute the Shift.
//At this point, IFIRST < ILAST, and the diagonal elements of
//T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
//magnitude)
	if ((iiter / 10) * 10 != iiter) {
//The Wilkinson shift (AEP p.512), i.e., the eigenvalue of
//the bottom-right 2x2 block of A inv(B) which is nearest to
//the bottom-right element.
//We factor B as U*D, where U has unit diagonals, and
//compute (A*inv(D))*inv(U).
	    u12 = bscale * t[ilast - 1 + ilast * ldt] / (bscale * t[ilast + ilast * ldt]);
	    ad11 = ascale * h[ilast - 1 + (ilast - 1) * ldh] / (bscale * t[ilast - 1 + (ilast - 1) * ldt]);
	    ad21 = ascale * h[ilast + (ilast - 1) * ldh] / (bscale * t[ilast - 1 + (ilast - 1) * ldt]);
	    ad12 = ascale * h[ilast - 1 + ilast * ldh] / (bscale * t[ilast + ilast * ldt]);
	    ad22 = ascale * h[ilast + ilast * ldh] / (bscale * t[ilast + ilast * ldt]);
	    abi22 = ad22 - u12 * ad21;
	    t1 = (ad11 + abi22) * Half;
	    rtdisc = sqrt(t1 * t1 + ad12 * ad21 - ad11 * ad22);
	    temp = (t1 - abi22).real() * (rtdisc).real() + (t1 - abi22).imag() * (rtdisc).imag();
	    if (temp <= Zero) {
		shift = t1 + rtdisc;
	    } else {
		shift = t1 - rtdisc;
	    }
	} else {
//Exceptional shift.  Chosen for no particularly good reason.
	    eshift = eshift + conj(ascale * h[ilast - 1 + ilast * ldh] / (bscale * t[ilast - 1 + (ilast - 1) * ldt]));
	    shift = eshift;
	}
//Now check for two consecutive small subdiagonals.
	for (j = ilast - 1; j >= ifirst + 1; j--) {
	    istart = j;
	    ctemp = ascale * h[j + j * ldh] - shift * (bscale * t[j + j * ldt]);
	    temp = Cabs1(ctemp);
	    temp2 = ascale * Cabs1(h[j + 1 + j * ldh]);
	    tempr = max(temp, temp2);
	    if (tempr < One && tempr != Zero) {
		temp = temp / tempr;
		temp2 = temp2 / tempr;
	    }
	    if (Cabs1(h[j + (j - 1) * ldh]) * temp2 <= temp * atol) {
		goto L90;
	    }
	}
	istart = ifirst;
	ctemp = ascale * h[ifirst + ifirst * ldh] - shift * (bscale * t[ifirst + ifirst * ldt]);
      L90:
//Do an implicit-shift QZ sweep.
//Initial Q
	ctemp2 = ascale * h[istart + 1 + istart * ldh];
	Clartg(ctemp, ctemp2, &c, &s, &ctemp3);
//Sweep
	for (j = istart; j <= ilast - 1; j++) {
	    if (j > istart) {
		ctemp = h[j + (j - 1) * ldh];
		Clartg(ctemp, h[j + 1 + (j - 1) * ldh], &c, &s, &h[j + (j - 1) * ldh]);
		h[j + 1 + (j - 1) * ldh] = Zero;
	    }
	    for (jc = j; jc <= ilastm; jc++) {
		ctemp = c * h[j + jc * ldh] + s * h[j + 1 + jc * ldh];
		h[j + 1 + jc * ldh] = -conj(s) * h[j + jc * ldh] + c * h[j + 1 + jc * ldh];
		h[j + jc * ldh] = ctemp;
		ctemp2 = c * t[j + jc * ldt] + s * t[j + 1 + jc * ldt];
		t[j + 1 + jc * ldt] = -conj(s) * t[j + jc * ldt] + c * t[j + 1 + jc * ldt];
		t[j + jc * ldt] = ctemp2;
	    }
	    if (ilq) {
		for (jr = 1; jr <= n; jr++) {
		    ctemp = c * q[jr + j * ldq] + conj(s) * q[jr + (j + 1) * ldq];
		    q[jr + (j + 1) * ldq] = -s * q[jr + j * ldq] + c * q[jr + (j + 1) * ldq];
		    q[jr + j * ldq] = ctemp;
		}
	    }
	    ctemp = t[j + 1 + (j + 1) * ldt];
	    Clartg(ctemp, t[j + 1 + j * ldt], &c, &s, &t[j + 1 + (j + 1) * ldt]);
	    t[j + 1 + j * ldt] = Zero;
	    for (jr = ifrstm; jr <= min(j + 2, ilast); jr++) {
		ctemp = c * h[jr + (j + 1) * ldh] + s * h[jr + j * ldh];
		h[jr + j * ldh] = -conj(s) * h[jr + (j + 1) * ldh] + c * h[jr + j * ldh];
		h[jr + (j + 1) * ldh] = ctemp;
	    }
	    for (jr = ifrstm; jr <= j; jr++) {
		ctemp = c * t[jr + (j + 1) * ldt] + s * t[jr + j * ldt];
		t[jr + j * ldt] = -conj(s) * t[jr + (j + 1) * ldt] + c * t[jr + j * ldt];
		t[jr + (j + 1) * ldt] = ctemp;
	    }
	    if (ilz) {
		for (jr = 1; jr <= n; jr++) {
		    ctemp = c * z[jr + (j + 1) * ldz] + s * z[jr + j * ldz];
		    z[jr + j * ldz] = -conj(s) * z[jr + (j + 1) * ldz] + c * z[jr + j * ldz];
		    z[jr + (j + 1) * ldz] = ctemp;
		}
	    }
	}
      L160:
	;
    }
//Drop-through = non-convergence
  L180:
    *info = ilast;
    goto L210;
//Successful completion of all QZ steps
  L190:
//Set Eigenvalues 1:ILO-1
    for (j = 0; j < ilo - 1; j++) {
	absb = abs(t[j + j * ldt]);
	if (absb > safmin) {
	    signbc = conj(t[j + j * ldt] / absb);
	    t[j + j * ldt] = absb;
	    if (ilschr) {
		Cscal(j - 1, signbc, &t[j * ldt + 1], 1);
		Cscal(j, signbc, &h[j * ldh + 1], 1);
	    } else {
		h[j + j * ldh] = h[j + j * ldh] * signbc;
	    }
	    if (ilz) {
		Cscal(n, signbc, &z[j * ldz + 1], 1);
	    }
	} else {
	    t[j + j * ldt] = Zero;
	}
	alpha[j] = h[j + j * ldh];
	beta[j] = t[j + j * ldt];
    }
//Normal Termination
    *info = 0;
//Exit (other than argument error) -- return optimal workspace size
  L210:
    work[1] = n;
    return;
}
