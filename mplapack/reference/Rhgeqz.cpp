/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rhgeqz.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rhgeqz(const char *job, const char *compq, const char *compz, INTEGER n,
	    INTEGER ilo, INTEGER ihi, REAL * h, INTEGER ldh, REAL
	    * t, INTEGER ldt, REAL * alphar, REAL * alphai, REAL * beta, REAL * q, INTEGER ldq, REAL * z, INTEGER ldz, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER j;
    REAL s, c, v[3], s1, s2, t1, u1, u2, a11, a12, a21, a22, b11, b22, c12, c21;
    INTEGER jc;
    REAL an, bn, cl, cq, cr;
    INTEGER in;
    REAL u12, w11, w12, w21;
    INTEGER jr;
    REAL cz, w22, sl, wi, sr, vs, wr, b1a, b2a, a1i, a2i, b1i, b2i, a1r, a2r, b1r, b2r, wr2, ad11, ad12, ad21, ad22, c11i, c22i;
    INTEGER jch;
    REAL c11r, c22r;
    INTEGER ilq = 0;
    REAL u12l, tau, sqi;
    INTEGER ilz = 0;
    REAL ulp, sqr, szi, szr, ad11l, ad12l, ad21l, ad22l, ad32l, wabs, atol, btol, temp;
    INTEGER iiter, ilast, jiter;
    REAL anorm, bnorm;
    INTEGER maxit;
    REAL temp2, tempi, tempr;
    INTEGER ilazr2;
    REAL ascale, bscale, scale;
    REAL safmin;
    REAL safmax;
    REAL eshift, s1inv;
    INTEGER ilschr = 0;
    INTEGER icompq, ilastm, ischur;
    INTEGER ilazro;
    INTEGER icompz, ifirst, ifrstm, istart;
    INTEGER ilpivt, lquery;
    REAL Zero = 0.0, Half = .5, One = 1.0, safety = 100.0;
    REAL mtemp1, mtemp2, mtemp3, mtemp4;

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
	*info = -15;
    } else if (ldz < 1 || (ilz && ldz < n)) {
	*info = -17;
    } else if (lwork < max((INTEGER) 1, n) && !lquery) {
	*info = -19;
    }
    if (*info != 0) {
	Mxerbla("Rhgeqz", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n <= 0) {
	work[1] = One;
	return;
    }
//Initialize Q and Z
    if (icompq == 3) {
	Rlaset("Full", n, n, Zero, One, &q[0], ldq);
    }
    if (icompz == 3) {
	Rlaset("Full", n, n, Zero, One, &z[0], ldz);
    }
//Machine Constants
    in = ihi + 1 - ilo;
    safmin = Rlamch("S");
    safmax = One / safmin;
    ulp = Rlamch("E") * Rlamch("B");
    anorm = Rlanhs("F", in, &h[ilo + ilo * ldh], ldh, &work[0]);
    bnorm = Rlanhs("F", in, &t[ilo + ilo * ldt], ldt, &work[0]);
    mtemp1 = safmin, mtemp2 = ulp * anorm;
    atol = max(mtemp1, mtemp2);
    mtemp1 = safmin, mtemp2 = ulp * bnorm;
    btol = max(mtemp1, mtemp2);
    ascale = One / max(safmin, anorm);
    bscale = One / max(safmin, bnorm);
//Set Eigenvalues IHI+1:N
    for (j = ihi + 1; j <= n; j++) {
	if (t[j + j * ldt] < Zero) {
	    if (ilschr) {
		for (jr = 1; jr <= j; jr++) {
		    h[jr + j * ldh] = -h[jr + j * ldh];
		    t[jr + j * ldt] = -t[jr + j * ldt];
		}
	    } else {
		h[j + j * ldh] = -h[j + j * ldh];
		t[j + j * ldt] = -t[j + j * ldt];
	    }
	    if (ilz) {
		for (jr = 1; jr <= n; jr++) {
		    z[jr + j * ldz] = -z[jr + j * ldz];
		}
	    }
	}
	alphar[j] = h[j + j * ldh];
	alphai[j] = Zero;
	beta[j] = t[j + j * ldt];
    }
//If IHI < ILO, skip QZ steps
    if (ihi < ilo) {
	goto L380;
    }
//MAIN QZ ITERATION LOOP
//Initialize dynamic indices
//Eigenvalues ILAST+1:N have been found.
//   Column operations modify rows IFRSTM:whatever.
//   Row operations modify columns whatever:ILASTM.
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
//Split the matrix if possible.
//Two tests:
//   1: H(j,j-1)=0  or  j=ILO
//   2: T(j,j)=0
	if (ilast == ilo) {
//Special case: j=ILAST
	    goto L80;
	} else {
	    if (abs(h[ilast + (ilast - 1) * ldh]) <= atol) {
		h[ilast + (ilast - 1) * ldh] = Zero;
		goto L80;
	    }
	}
	if (abs(t[ilast + ilast * ldt]) <= btol) {
	    t[ilast + ilast * ldt] = Zero;
	    goto L70;
	}
//General case: j<ILAST
	for (j = ilast - 1; j >= ilo; j--) {
//Test 1: for H(j,j-1)=0 or j=ILO
	    if (j == ilo) {
		ilazro = MTRUE;
	    } else {
		if (abs(h[j + (j - 1) * ldh]) <= atol) {
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
		    temp = abs(h[j + (j - 1) * ldh]);
		    temp2 = abs(h[j + j * ldh]);
		    tempr = max(temp, temp2);
		    if (tempr < One && tempr != Zero) {
			temp /= tempr;
			temp2 /= tempr;
		    }
		    if (temp * (ascale * abs(h[j + 1 + j * ldh])) <= temp2 * (ascale * atol)) {
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
			temp = h[jch + jch * ldh];
			Rlartg(temp, h[jch + 1 + jch * ldh], &c, &s, &h[jch + jch * ldh]);
			h[jch + 1 + jch * ldh] = Zero;
			Rrot(ilastm - jch, &h[jch + (jch + 1) * ldh], ldh, &h[jch + 1 + (jch + 1) * ldh], ldh, c, s);
			Rrot(ilastm - jch, &t[jch + (jch + 1) * ldt], ldt, &t[jch + 1 + (jch + 1) * ldt], ldt, c, s);
			if (ilq) {
			    Rrot(n, &q[jch * ldq + 1], 1, &q[(jch + 1) * ldq + 1], 1, c, s);
			}
			if (ilazr2) {
			    h[jch + (jch - 1) * ldh] = h[jch + (jch - 1) * ldh] * c;
			}
			ilazr2 = MFALSE;
			if (abs(t[jch + 1 + (jch + 1) * ldt]) >= btol) {
			    if (jch + 1 >= ilast) {
				goto L80;
			    } else {
				ifirst = jch + 1;
				goto L110;
			    }
			}
			t[jch + 1 + (jch + 1) * ldt] = Zero;
		    }
		    goto L70;
		} else {
//Only test 2 passed -- chase the zero to T(ILAST,ILAST)
//Then process as in the case T(ILAST,ILAST)=0
		    for (jch = j; jch <= ilast - 1; jch++) {
			temp = t[jch + (jch + 1) * ldt];
			Rlartg(temp, t[jch + 1 + (jch + 1) * ldt], &c, &s, &t[jch + (jch + 1) * ldt]);
			t[jch + 1 + (jch + 1) * ldt] = Zero;
			if (jch < ilastm - 1) {
			    Rrot(ilastm - jch - 1, &t[jch + (jch + 2) * ldt], ldt, &t[jch + 1 + (jch + 2) * ldt], ldt, c, s);
			}
			Rrot(ilastm - jch + 2, &h[jch + (jch - 1) * ldh], ldh, &h[jch + 1 + (jch - 1) * ldh], ldh, c, s);
			if (ilq) {
			    Rrot(n, &q[jch * ldq + 1], 1, &q[(jch + 1) * ldq + 1], 1, c, s);
			}
			temp = h[jch + 1 + jch * ldh];
			Rlartg(temp, h[jch + 1 + (jch - 1) * ldh], &c, &s, &h[jch + 1 + jch * ldh]);
			h[jch + 1 + (jch - 1) * ldh] = Zero;
			Rrot(jch + 1 - ifrstm, &h[ifrstm + jch * ldh], 1, &h[ifrstm + (jch - 1) * ldh], 1, c, s);
			Rrot(jch - ifrstm, &t[ifrstm + jch * ldt], 1, &t[ifrstm + (jch - 1) * ldt], 1, c, s);
			if (ilz) {
			    Rrot(n, &z[jch * ldz + 1], 1, &z[(jch - 1) * ldz + 1], 1, c, s);
			}
		    }
		    goto L70;
		}
	    } else if (ilazro) {
//Only test 1 passed -- work on J:ILAST
		ifirst = j;
		goto L110;
	    }
//Neither test passed -- try next J
	}
//(Drop-through is "impossible")
	*info = n + 1;
	goto L420;
//T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
//1x1 block.
      L70:
	temp = h[ilast + ilast * ldh];
	Rlartg(temp, h[ilast + (ilast - 1) * ldh], &c, &s, &h[ilast + ilast * ldh]);
	h[ilast + (ilast - 1) * ldh] = Zero;
	Rrot(ilast - ifrstm, &h[ifrstm + ilast * ldh], 1, &h[ifrstm + (ilast - 1) * ldh], 1, c, s);
	Rrot(ilast - ifrstm, &t[ifrstm + ilast * ldt], 1, &t[ifrstm + (ilast - 1) * ldt], 1, c, s);
	if (ilz) {
	    Rrot(n, &z[ilast * ldz + 1], 1, &z[(ilast - 1) * ldz + 1], 1, c, s);
	}
//H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
//                      and BETA
      L80:
	if (t[ilast + ilast * ldt] < Zero) {
	    if (ilschr) {
		for (j = ifrstm; j <= ilast; j++) {
		    h[j + ilast * ldh] = -h[j + ilast * ldh];
		    t[j + ilast * ldt] = -t[j + ilast * ldt];
		}
	    } else {
		h[ilast + ilast * ldh] = -h[ilast + ilast * ldh];
		t[ilast + ilast * ldt] = -t[ilast + ilast * ldt];
	    }
	    if (ilz) {
		for (j = 0; j < n; j++) {
		    z[j + ilast * ldz] = -z[j + ilast * ldz];
		}
	    }
	}
	alphar[ilast] = h[ilast + ilast * ldh];
	alphai[ilast] = Zero;
	beta[ilast] = t[ilast + ilast * ldt];
//Go to next block -- exit if finished.
	ilast--;
	if (ilast < ilo) {
	    goto L380;
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
	goto L350;
//QZ step
//This iteration only involves rows/columns IFIRST:ILAST. We
//assume IFIRST < ILAST, and that the diagonal of B is non-zero.
      L110:
	iiter++;
	if (!ilschr) {
	    ifrstm = ifirst;
	}
//Compute single shifts.
//At this point, IFIRST < ILAST, and the diagonal elements of
//T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
//magnitude)
	if ((iiter / 10) * 10 == iiter) {
//Exceptional shift.  Chosen for no particularly good reason.
//(Single shift only.)
	    if (maxit * safmin * abs(h[ilast - 1 + ilast * ldh]) < abs(t[ilast - 1 + (ilast - 1) * ldt])) {
		eshift = eshift + h[ilast - 1 + ilast * ldh] / t[ilast - 1 + (ilast - 1) * ldt];
	    } else {
		eshift = eshift + One / (safmin * maxit);
	    }
	    s1 = One;
	    wr = eshift;
	} else {
//Shifts based on the generalized eigenvalues of the
//bottom-right 2x2 block of A and B. The first eigenvalue
//returned by DLAG2 is the Wilkinson shift (AEP p.512),
	    Rlag2(&h[ilast - 1 + (ilast - 1) * ldh], ldh, &t[ilast - 1 + (ilast - 1) * ldt], ldt, safmin * safety, &s1, &s2, &wr, &wr2, &wi);
	    mtemp1 = One, mtemp2 = abs(wr), mtemp3 = max(mtemp1, mtemp2), mtemp4 = abs(wi);
	    mtemp1 = s1, mtemp2 = safmin * max(mtemp3, mtemp4);
	    temp = max(mtemp1, mtemp2);
	    if (wi != Zero) {
		goto L200;
	    }
	}
//Fiddle with shift to avoid overflow
	temp = min(ascale, One) * (safmax * Half);
	if (s1 > temp) {
	    scale = temp / s1;
	} else {
	    scale = One;
	}
	temp = min(bscale, One) * (safmax * Half);
	if (abs(wr) > temp) {
	    mtemp1 = scale, mtemp2 = temp / abs(wr);
	    scale = min(mtemp1, mtemp2);
	}
	s1 = scale * s1;
	wr = scale * wr;
//Now check for two consecutive small subdiagonals.
	for (j = ilast - 1; j >= ifirst + 1; j--) {
	    istart = j;
	    temp = abs(s1 * h[j + (j - 1) * ldh]);
	    temp2 = abs(s1 * h[j + j * ldh] - wr * t[j + j * ldt]);
	    tempr = max(temp, temp2);
	    if (tempr < One && tempr != Zero) {
		temp = temp / tempr;
		temp2 = temp2 / tempr;
	    }
	    if (abs(ascale * h[j + 1 + j * ldh] * temp) <= ascale * atol * temp2) {
		goto L130;
	    }
	}
	istart = ifirst;
      L130:
//Do an implicit single-shift QZ sweep.
//Initial Q
	temp = s1 * h[istart + istart * ldh] - wr * t[istart + istart * ldt];
	temp2 = s1 * h[istart + 1 + istart * ldh];
	Rlartg(temp, temp2, &c, &s, &tempr);
//Sweep
	for (j = istart; j <= ilast - 1; j++) {
	    if (j > istart) {
		temp = h[j + (j - 1) * ldh];
		Rlartg(temp, h[j + 1 + (j - 1) * ldh], &c, &s, &h[j + (j - 1) * ldh]);
		h[j + 1 + (j - 1) * ldh] = Zero;
	    }
	    for (jc = j; jc <= ilastm; jc++) {
		temp = c * h[j + jc * ldh] + s * h[j + 1 + jc * ldh];
		h[j + 1 + jc * ldh] = -s * h[j + jc * ldh] + c * h[j + 1 + jc * ldh];
		h[j + jc * ldh] = temp;
		temp2 = c * t[j + jc * ldt] + s * t[j + 1 + jc * ldt];
		t[j + 1 + jc * ldt] = -s * t[j + jc * ldt] + c * t[j + 1 + jc * ldt];
		t[j + jc * ldt] = temp2;
	    }
	    if (ilq) {
		for (jr = 1; jr <= n; jr++) {
		    temp = c * q[jr + j * ldq] + s * q[jr + (j + 1) * ldq];
		    q[jr + (j + 1) * ldq] = -s * q[jr + j * ldq] + c * q[jr + (j + 1) * ldq];
		    q[jr + j * ldq] = temp;
		}
	    }
	    temp = t[j + 1 + (j + 1) * ldt];
	    Rlartg(temp, t[j + 1 + j * ldt], &c, &s, &t[j + 1 + (j + 1) * ldt]);
	    t[j + 1 + j * ldt] = Zero;
	    for (jr = ifrstm; jr <= min(j + 2, ilast); jr++) {
		temp = c * h[jr + (j + 1) * ldh] + s * h[jr + j * ldh];
		h[jr + j * ldh] = -s * h[jr + (j + 1) * ldh] + c * h[jr + j * ldh];
		h[jr + (j + 1) * ldh] = temp;
	    }
	    for (jr = ifrstm; jr <= j; jr++) {
		temp = c * t[jr + (j + 1) * ldt] + s * t[jr + j * ldt];
		t[jr + j * ldt] = -s * t[jr + (j + 1) * ldt] + c * t[jr + j * ldt];
		t[jr + (j + 1) * ldt] = temp;
	    }
	    if (ilz) {
		for (jr = 1; jr <= n; jr++) {
		    temp = c * z[jr + (j + 1) * ldz] + s * z[jr + j * ldz];
		    z[jr + j * ldz] = -s * z[jr + (j + 1) * ldz] + c * z[jr + j * ldz];
		    z[jr + (j + 1) * ldz] = temp;
		}
	    }
	}
	goto L350;
//Use Francis double-shift
//Note: the Francis double-shift should work with real shifts,
//      but only if the block is at least 3x3.
//      This code may break if this point is reached with
//      a 2x2 block with real eigenvalues.
      L200:
	if (ifirst + 1 == ilast) {
//Special case -- 2x2 block with complex eigenvectors
//Step 1: Standardize, that is, rotate so that
//            ( B11  0  )
//        B = (         )  with B11 non-negative.
//            (  0  B22 )
	    Rlasv2(t[ilast - 1 + (ilast - 1) * ldt], t[ilast - 1 + ilast * ldt], t[ilast + ilast * ldt], &b22, &b11, &sr, &cr, &sl, &cl);
	    if (b11 < Zero) {
		cr = -cr;
		sr = -sr;
		b11 = -b11;
		b22 = -b22;
	    }
	    Rrot(ilastm + 1 - ifirst, &h[ilast - 1 + (ilast - 1) * ldh], ldh, &h[ilast + (ilast - 1) * ldh], ldh, cl, sl);
	    Rrot(ilast + 1 - ifrstm, &h[ifrstm + (ilast - 1) * ldh], 1, &h[ifrstm + ilast * ldh], 1, cr, sr);
	    if (ilast < ilastm) {
		Rrot(ilastm - ilast, &t[ilast - 1 + (ilast + 1) * ldt], ldt, &t[ilast + (ilast + 1) * ldt], ldh, cl, sl);
	    }
	    if (ifrstm < ilast - 1) {
		Rrot(ifirst - ifrstm, &t[ifrstm + (ilast - 1) * ldt], 1, &t[ifrstm + ilast * ldt], 1, cr, sr);
	    }
	    if (ilq) {
		Rrot(n, &q[(ilast - 1) * ldq + 1], 1, &q[ilast * ldq + 1], 1, cl, sl);
	    }
	    if (ilz) {
		Rrot(n, &z[(ilast - 1) * ldz + 1], 1, &z[ilast * ldz + 1], 1, cr, sr);
	    }
	    t[ilast - 1 + (ilast - 1) * ldt] = b11;
	    t[ilast - 1 + ilast * ldt] = Zero;
	    t[ilast + (ilast - 1) * ldt] = Zero;
	    t[ilast + ilast * ldt] = b22;
//If B22 is negative, negate column ILAST
	    if (b22 < Zero) {
		for (j = ifrstm; j <= ilast; j++) {
		    h[j + ilast * ldh] = -h[j + ilast * ldh];
		    t[j + ilast * ldt] = -t[j + ilast * ldt];
		}
		if (ilz) {
		    for (j = 0; j < n; j++) {
			z[j + ilast * ldz] = -z[j + ilast * ldz];
		    }
		}
	    }
//Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)
//Recompute shift
	    Rlag2(&h[ilast - 1 + (ilast - 1) * ldh], ldh, &t[ilast - 1 + (ilast - 1) * ldt], ldt, safmin * safety, &s1, &temp, &wr, &temp2, &wi);
//If standardization has perturbed the shift onto real line,
//do another (real single-shift) QR step.
	    if (wi == Zero) {
		goto L350;
	    }
	    s1inv = One / s1;
//Do EISPACK (QZVAL) computation of alpha and beta
	    a11 = h[ilast - 1 + (ilast - 1) * ldh];
	    a21 = h[ilast + (ilast - 1) * ldh];
	    a12 = h[ilast - 1 + ilast * ldh];
	    a22 = h[ilast + ilast * ldh];
//Compute complex Givens rotation on right
//(Assume some element of C = (sA - wB) > unfl )
//                 __
//(sA - wB) ( CZ   -SZ )
//          ( SZ    CZ )
	    c11r = s1 * a11 - wr * b11;
	    c11i = -wi * b11;
	    c12 = s1 * a12;
	    c21 = s1 * a21;
	    c22r = s1 * a22 - wr * b22;
	    c22i = -wi * b22;
	    if (abs(c11r) + abs(c11i) + abs(c12) > abs(c21) + abs(c22r) + abs(c22i)) {
		t1 = Rlapy3(c12, c11r, c11i);
		cz = c12 / t1;
		szr = -c11r / t1;
		szi = -c11i / t1;
	    } else {
		cz = Rlapy2(c22r, c22i);
		if (cz <= safmin) {
		    cz = Zero;
		    szr = One;
		    szi = Zero;
		} else {
		    tempr = c22r / cz;
		    tempi = c22i / cz;
		    t1 = Rlapy2(cz, c21);
		    cz = cz / t1;
		    szr = -c21 * tempr / t1;
		    szi = c21 * tempi / t1;
		}
	    }
//Compute Givens rotation on left
//(  CQ   SQ )
//(  __      )  A or B
//( -SQ   CQ )
	    an = abs(a11) + abs(a12) + abs(a21) + abs(a22);
	    bn = abs(b11) + abs(b22);
	    wabs = abs(wr) + abs(wi);
	    if (s1 * an > wabs * bn) {
		cq = cz * b11;
		sqr = szr * b22;
		sqi = -szi * b22;
	    } else {
		a1r = cz * a11 + szr * a12;
		a1i = szi * a12;
		a2r = cz * a21 + szr * a22;
		a2i = szi * a22;
		cq = Rlapy2(a1r, a1i);
		if (cq <= safmin) {
		    cq = Zero;
		    sqr = One;
		    sqi = Zero;
		} else {
		    tempr = a1r / cq;
		    tempi = a1i / cq;
		    sqr = tempr * a2r + tempi * a2i;
		    sqi = tempi * a2r - tempr * a2i;
		}
	    }
	    t1 = Rlapy3(cq, sqr, sqi);
	    cq = cq / t1;
	    sqr = sqr / t1;
	    sqi = sqi / t1;
//Compute diagonal elements of QBZ
	    tempr = sqr * szr - sqi * szi;
	    tempi = sqr * szi + sqi * szr;
	    b1r = cq * cz * b11 + tempr * b22;
	    b1i = tempi * b22;
	    b1a = Rlapy2(b1r, b1i);
	    b2r = cq * cz * b22 + tempr * b11;
	    b2i = -tempi * b11;
	    b2a = Rlapy2(b2r, b2i);
//Normalize so beta > 0, and Im( alpha1 ) > 0
	    beta[ilast - 1] = b1a;
	    beta[ilast] = b2a;
	    alphar[ilast - 1] = wr * b1a * s1inv;
	    alphai[ilast - 1] = wi * b1a * s1inv;
	    alphar[ilast] = wr * b2a * s1inv;
	    alphai[ilast] = -(wi * b2a) * s1inv;
//Step 3: Go to next block -- exit if finished.
	    ilast = ifirst - 1;
	    if (ilast < ilo) {
		goto L380;
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
	    goto L350;
	} else {
//Usual case: 3x3 or larger block, using Francis implicit
//            double-shift
//                         2
//Eigenvalue equation is  w  - c w + d = 0,
//                              -1 2        -1
//so compute 1st column of  (A B  )  - c A B   + d
//using the formula in QZIT (from EISPACK)
//We assume that the block is at least 3x3
	    ad11 = ascale * h[ilast - 1 + (ilast - 1) * ldh] / (bscale * t[ilast - 1 + (ilast - 1) * ldt]);
	    ad21 = ascale * h[ilast + (ilast - 1) * ldh] / (bscale * t[ilast - 1 + (ilast - 1) * ldt]);
	    ad12 = ascale * h[ilast - 1 + ilast * ldh] / (bscale * t[ilast + ilast * ldt]);
	    ad22 = ascale * h[ilast + ilast * ldh] / (bscale * t[ilast + ilast * ldt]);
	    u12 = t[ilast - 1 + ilast * ldt] / t[ilast + ilast * ldt];
	    ad11l = ascale * h[ifirst + ifirst * ldh] / (bscale * t[ifirst + ifirst * ldt]);
	    ad21l = ascale * h[ifirst + 1 + ifirst * ldh] / (bscale * t[ifirst + ifirst * ldt]);
	    ad12l = ascale * h[ifirst + (ifirst + 1) * ldh] / (bscale * t[ifirst + 1 + (ifirst + 1) * ldt]);
	    ad22l = ascale * h[ifirst + 1 + (ifirst + 1) * ldh] / (bscale * t[ifirst + 1 + (ifirst + 1) * ldt]);
	    ad32l = ascale * h[ifirst + 2 + (ifirst + 1) * ldh] / (bscale * t[ifirst + 1 + (ifirst + 1) * ldt]);
	    u12l = t[ifirst + (ifirst + 1) * ldt] / t[ifirst + 1 + (ifirst + 1) * ldt];
	    v[0] = (ad11 - ad11l) * (ad22 - ad11l) - ad12 * ad21 + ad21 * u12 * ad11l + (ad12l - ad11l * u12l) * ad21l;
	    v[1] = (ad22l - ad11l - ad21l * u12l - (ad11 - ad11l) - (ad22 - ad11l) + ad21 * u12) * ad21l;
	    v[2] = ad32l * ad21l;
	    istart = ifirst;
	    Rlarfg(3, v, &v[1], 1, &tau);
	    v[0] = One;
//Sweep
	    for (j = istart; j <= ilast - 2; j++) {
//All but last elements: use 3x3 Householder transforms.
//Zero (j-1)st column of A
		if (j > istart) {
		    v[0] = h[j + (j - 1) * ldh];
		    v[1] = h[j + 1 + (j - 1) * ldh];
		    v[2] = h[j + 2 + (j - 1) * ldh];
		    Rlarfg(3, &h[j + (j - 1) * ldh], &v[1], 1, &tau);
		    v[0] = One;
		    h[j + 1 + (j - 1) * ldh] = Zero;
		    h[j + 2 + (j - 1) * ldh] = Zero;
		}
		for (jc = j; jc <= ilastm; jc++) {
		    temp = tau * (h[j + jc * ldh] + v[1] * h[j + 1 + jc * ldh] + v[2] * h[j + 2 + jc * ldh]);
		    h[j + jc * ldh] = h[j + jc * ldh] - temp;
		    h[j + 1 + jc * ldh] = h[j + 1 + jc * ldh] - temp * v[1];
		    h[j + 2 + jc * ldh] = h[j + 2 + jc * ldh] - temp * v[2];
		    temp2 = tau * (t[j + jc * ldt] + v[1] * t[j + 1 + jc * ldt] + v[2] * t[j + 2 + jc * ldt]);
		    t[j + jc * ldt] = t[j + jc * ldt] - temp2;
		    t[j + 1 + jc * ldt] = t[j + 1 + jc * ldt] - temp2 * v[1];
		    t[j + 2 + jc * ldt] = t[j + 2 + jc * ldt] - temp2 * v[2];
		}
		if (ilq) {
		    for (jr = 1; jr <= n; jr++) {
			temp = tau * (q[jr + j * ldq] + v[1] * q[jr + (j + 1) * ldq] + v[2] * q[jr + (j + 2) * ldq]);
			q[jr + j * ldq] = q[jr + j * ldq] - temp;
			q[jr + (j + 1) * ldq] = q[jr + (j + 1) * ldq] - temp * v[1];
			q[jr + (j + 2) * ldq] = q[jr + (j + 2) * ldq] - temp * v[2];
		    }
		}
//Zero j-th column of B (see DLAGBC for details)
//Swap rows to pivot
		ilpivt = MFALSE;
		mtemp1 = abs(t[j + 1 + (j + 1) * ldt]), mtemp2 = abs(t[j + 1 + (j + 2) * ldt]);
		temp = max(mtemp1, mtemp2);
		mtemp1 = abs(t[j + 2 + (j + 1) * ldt]), mtemp2 = abs(t[j + 2 + (j + 2) * ldt]);
		temp2 = max(mtemp1, mtemp2);
		if (max(temp, temp2) < safmin) {
		    scale = Zero;
		    u1 = One;
		    u2 = Zero;
		    goto L250;
		} else if (temp >= temp2) {
		    w11 = t[j + 1 + (j + 1) * ldt];
		    w21 = t[j + 2 + (j + 1) * ldt];
		    w12 = t[j + 1 + (j + 2) * ldt];
		    w22 = t[j + 2 + (j + 2) * ldt];
		    u1 = t[j + 1 + j * ldt];
		    u2 = t[j + 2 + j * ldt];
		} else {
		    w21 = t[j + 1 + (j + 1) * ldt];
		    w11 = t[j + 2 + (j + 1) * ldt];
		    w22 = t[j + 1 + (j + 2) * ldt];
		    w12 = t[j + 2 + (j + 2) * ldt];
		    u2 = t[j + 1 + j * ldt];
		    u1 = t[j + 2 + j * ldt];
		}
//Swap columns if nec.
		if (abs(w12) > abs(w11)) {
		    ilpivt = MTRUE;
		    temp = w12;
		    temp2 = w22;
		    w12 = w11;
		    w22 = w21;
		    w11 = temp;
		    w21 = temp2;
		}
//LU-factor
		temp = w21 / w11;
		u2 = u2 - temp * u1;
		w22 = w22 - temp * w12;
		w21 = Zero;
//Compute SCALE
		scale = One;
		if (abs(w22) < safmin) {
		    scale = Zero;
		    u2 = One;
		    u1 = -w12 / w11;
		    goto L250;
		}
		if (abs(w22) < abs(u2)) {
		    scale = abs(w22 / u2);
		}
		if (abs(w11) < abs(u1)) {
		    mtemp1 = scale, mtemp2 = abs(w11 / u1);
		    scale = min(mtemp1, mtemp2);
		}
//Solve
		u2 = scale * u2 / w22;
		u1 = (scale * u1 - w12 * u2) / w11;
	      L250:
		if (ilpivt) {
		    temp = u2;
		    u2 = u1;
		    u1 = temp;
		}
//Compute Householder Vector
		t1 = sqrt(scale * scale + u1 * u1 + u2 * u2);
		tau = scale / t1 + One;
		vs = -One / (scale + t1);
		v[0] = One;
		v[1] = vs * u1;
		v[2] = vs * u2;
//Apply transformations from the right.
		for (jr = ifrstm; jr <= min(j + 3, ilast); jr++) {
		    temp = tau * (h[jr + j * ldh] + v[1] * h[jr + (j + 1) * ldh] + v[2] * h[jr + (j + 2) * ldh]);
		    h[jr + j * ldh] = h[jr + j * ldh] - temp;
		    h[jr + (j + 1) * ldh] = h[jr + (j + 1) * ldh] - temp * v[1];
		    h[jr + (j + 2) * ldh] = h[jr + (j + 2) * ldh] - temp * v[2];
		}
		for (jr = ifrstm; jr <= j + 2; jr++) {
		    temp = tau * (t[jr + j * ldt] + v[1] * t[jr + (j + 1) * ldt] + v[2] * t[jr + (j + 2) * ldt]);
		    t[jr + j * ldt] = t[jr + j * ldt] - temp;
		    t[jr + (j + 1) * ldt] = t[jr + (j + 1) * ldt] - temp * v[1];
		    t[jr + (j + 2) * ldt] = t[jr + (j + 2) * ldt] - temp * v[2];
		}
		if (ilz) {
		    for (jr = 1; jr <= n; jr++) {
			temp = tau * (z[jr + j * ldz] + v[1] * z[jr + (j + 1) * ldz] + v[2] * z[jr + (j + 2) * ldz]);
			z[jr + j * ldz] = z[jr + j * ldz] - temp;
			z[jr + (j + 1) * ldz] = z[jr + (j + 1) * ldz] - temp * v[1];
			z[jr + (j + 2) * ldz] = z[jr + (j + 2) * ldz] - temp * v[2];
		    }
		}
		t[j + 1 + j * ldt] = Zero;
		t[j + 2 + j * ldt] = Zero;
	    }
//Last elements: Use Givens rotations
//Rotations from the left
	    j = ilast - 1;
	    temp = h[j + (j - 1) * ldh];
	    Rlartg(temp, h[j + 1 + (j - 1) * ldh], &c, &s, &h[j + (j - 1) * ldh]);
	    h[j + 1 + (j - 1) * ldh] = Zero;
	    for (jc = j; jc <= ilastm; jc++) {
		temp = c * h[j + jc * ldh] + s * h[j + 1 + jc * ldh];
		h[j + 1 + jc * ldh] = -s * h[j + jc * ldh] + c * h[j + 1 + jc * ldh];
		h[j + jc * ldh] = temp;
		temp2 = c * t[j + jc * ldt] + s * t[j + 1 + jc * ldt];
		t[j + 1 + jc * ldt] = -s * t[j + jc * ldt] + c * t[j + 1 + jc * ldt];
		t[j + jc * ldt] = temp2;
	    }
	    if (ilq) {
		for (jr = 1; jr <= n; jr++) {
		    temp = c * q[jr + j * ldq] + s * q[jr + (j + 1) * ldq];
		    q[jr + (j + 1) * ldq] = -s * q[jr + j * ldq] + c * q[jr + (j + 1) * ldq];
		    q[jr + j * ldq] = temp;
		}
	    }
//Rotations from the right.
	    temp = t[j + 1 + (j + 1) * ldt];
	    Rlartg(temp, t[j + 1 + j * ldt], &c, &s, &t[j + 1 + (j + 1) * ldt]);
	    t[j + 1 + j * ldt] = Zero;
	    for (jr = ifrstm; jr <= ilast; jr++) {
		temp = c * h[jr + (j + 1) * ldh] + s * h[jr + j * ldh];
		h[jr + j * ldh] = -s * h[jr + (j + 1) * ldh] + c * h[jr + j * ldh];
		h[jr + (j + 1) * ldh] = temp;
	    }
	    for (jr = ifrstm; jr <= ilast - 1; jr++) {
		temp = c * t[jr + (j + 1) * ldt] + s * t[jr + j * ldt];
		t[jr + j * ldt] = -s * t[jr + (j + 1) * ldt] + c * t[jr + j * ldt];
		t[jr + (j + 1) * ldt] = temp;
	    }
	    if (ilz) {
		for (jr = 1; jr <= n; jr++) {
		    temp = c * z[jr + (j + 1) * ldz] + s * z[jr + j * ldz];
		    z[jr + j * ldz] = -s * z[jr + (j + 1) * ldz] + c * z[jr + j * ldz];
		    z[jr + (j + 1) * ldz] = temp;
		}
	    }
//End of Double-Shift code
	}
	goto L350;
//End of iteration loop
      L350:
	;
    }
//Drop-through = non-convergence
    *info = ilast;
    goto L420;
//Successful completion of all QZ steps
  L380:
//Set Eigenvalues 1:ILO-1
    for (j = 0; j < ilo - 1; j++) {
	if (t[j + j * ldt] < Zero) {
	    if (ilschr) {
		for (jr = 1; jr <= j; jr++) {
		    h[jr + j * ldh] = -h[jr + j * ldh];
		    t[jr + j * ldt] = -t[jr + j * ldt];
		}
	    } else {
		h[j + j * ldh] = -h[j + j * ldh];
		t[j + j * ldt] = -t[j + j * ldt];
	    }
	    if (ilz) {
		for (jr = 1; jr <= n; jr++) {
		    z[jr + j * ldz] = -z[jr + j * ldz];
		}
	    }
	}
	alphar[j] = h[j + j * ldh];
	alphai[j] = Zero;
	beta[j] = t[j + j * ldt];
    }
//Normal Termination
    *info = 0;
//Exit (other than argument error) -- return optimal workspace size
  L420:
    work[1] = n;
    return;
}
