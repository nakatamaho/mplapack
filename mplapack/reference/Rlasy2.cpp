/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasy2.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void
Rlasy2(INTEGER ltranl, INTEGER ltranr, INTEGER isgn, INTEGER n1, INTEGER n2, REAL * tl,
       INTEGER ldtl, REAL * tr, INTEGER ldtr, REAL * B, INTEGER ldb, REAL * scale, REAL * x, INTEGER ldx, REAL * xnorm, INTEGER * info)
{
#if defined done 
    INTEGER locu12[4] = { 3, 4, 1, 2 };
    INTEGER locl21[4] = { 2, 1, 4, 3 };
    INTEGER locu22[4] = { 4, 3, 2, 1 };
    INTEGER xswpiv[4] = { MFALSE, MFALSE, MTRUE, MTRUE };
    INTEGER bswpiv[4] = { MFALSE, MTRUE, MFALSE, MTRUE };
    INTEGER i, j, k;
    REAL x2[2], l21, u11, u12;
    INTEGER ip, jp;
    REAL u22, t16[16] /* was [4][4] */ , gam, bet, eps, sgn, tmp[4], tau1, btmp[4], smin;
    INTEGER ipiv;
    REAL temp;
    INTEGER jpiv[4];
    REAL xmax;
    INTEGER ipsv = 0, jpsv = 0;
    INTEGER bswap;
    INTEGER xswap;
    REAL smlnum;
    REAL One = 1.0, Half = 0.5, Zero = 0.0, Two = 2.0, Eight = 8.0, Fourth = 0.125;
    REAL mtemp1, mtemp2, mtemp3;

//Do not check the input parameters for errors
    *info = 0;
//Quick return if possible
    if (n1 == 0 || n2 == 0) {
	return;
    }
//Set constants to control overflow
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    sgn = isgn;
    k = n1 + n1 + n2 - 2;
    switch (k) {
    case 1:
	goto L10;
    case 2:
	goto L20;
    case 3:
	goto L30;
    case 4:
	goto L50;
    }
//1 by 1: TL11*X + SGN*X*TR11 = B11
  L10:
    tau1 = tl[ldtl + 1] + sgn * tr[ldtr + 1];
    bet = abs(tau1);
    if (bet <= smlnum) {
	tau1 = smlnum;
	bet = smlnum;
	*info = 1;
    }
    *scale = One;
    gam = abs(B[ldb + 1]);
    if (smlnum * gam > bet) {
	*scale = One / gam;
    }
    x[ldx + 1] = B[ldb + 1] * *scale / tau1;
    *xnorm = abs(x[ldx + 1]);
    return;
//1 by 2:
//TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
//                                  [TR21 TR22]
  L20:
    mtemp1 = abs(tl[ldtl + 1]), mtemp2 = abs(tr[ldtr + 1]), mtemp1 = max(mtemp1, mtemp2);
    mtemp2 = abs(tr[(ldtr << 1) + 1]), mtemp1 = max(mtemp1, mtemp2);
    mtemp2 = abs(tr[ldtr + 2]), mtemp1 = max(mtemp1, mtemp2);
    mtemp2 = abs(tr[(ldtr << 1) + 2]), mtemp1 = eps * max(mtemp1, mtemp2);
    smin = max(mtemp1, smlnum);
    tmp[0] = tl[ldtl + 1] + sgn * tr[ldtr + 1];
    tmp[3] = tl[ldtl + 1] + sgn * tr[(ldtr << 1) + 2];
    if (ltranr) {
	tmp[1] = sgn * tr[ldtr + 2];
	tmp[2] = sgn * tr[(ldtr << 1) + 1];
    } else {
	tmp[1] = sgn * tr[(ldtr << 1) + 1];
	tmp[2] = sgn * tr[ldtr + 2];
    }
    btmp[0] = B[ldb + 1];
    btmp[1] = B[(ldb << 1) + 1];
    goto L40;
//2 by 1:
//     op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
//       [TL21 TL22] [X21]         [X21]         [B21]
  L30:
    mtemp1 = abs(tr[ldtr + 1]), mtemp2 = abs(tl[ldtl + 1]), mtemp3 = max(mtemp1, mtemp2);
    mtemp1 = abs(tl[(ldtl << 1) + 1]), mtemp2 = max(mtemp1, mtemp3);
    mtemp1 = abs(tl[ldtl + 2]), mtemp3 = max(mtemp1, mtemp2);
    mtemp1 = abs(tl[(ldtl << 1) + 2]);
    mtemp2 = eps * max(mtemp1, mtemp3);
    smin = max(mtemp2, smlnum);
    tmp[0] = tl[ldtl + 1] + sgn * tr[ldtr + 1];
    tmp[3] = tl[(ldtl << 1) + 2] + sgn * tr[ldtr + 1];
    if (ltranl) {
	tmp[1] = tl[(ldtl << 1) + 1];
	tmp[2] = tl[ldtl + 2];
    } else {
	tmp[1] = tl[ldtl + 2];
	tmp[2] = tl[(ldtl << 1) + 1];
    }
    btmp[0] = B[ldb + 1];
    btmp[1] = B[ldb + 2];
  L40:
//Solve 2 by 2 system using complete pivoting.
//Set pivots less than SMIN to SMIN.
    ipiv = iRamax(4, tmp, 1);
    u11 = tmp[ipiv - 1];
    if (abs(u11) <= smin) {
	*info = 1;
	u11 = smin;
    }
    u12 = tmp[locu12[ipiv - 1] - 1];
    l21 = tmp[locl21[ipiv - 1] - 1] / u11;
    u22 = tmp[locu22[ipiv - 1] - 1] - u12 * l21;
    xswap = xswpiv[ipiv - 1];
    bswap = bswpiv[ipiv - 1];
    if (abs(u22) <= smin) {
	*info = 1;
	u22 = smin;
    }
    if (bswap) {
	temp = btmp[1];
	btmp[1] = btmp[0] - l21 * temp;
	btmp[0] = temp;
    } else {
	btmp[1] -= l21 * btmp[0];
    }
    *scale = One;
    if (smlnum * Two * abs(btmp[1]) > abs(u22)
	|| smlnum * Two * abs(btmp[0]) > abs(u11)) {
	mtemp1 = abs(btmp[0]), mtemp2 = abs(btmp[1]);
	*scale = Half / max(mtemp1, mtemp2);
	btmp[0] *= *scale;
	btmp[1] *= *scale;
    }
    x2[1] = btmp[1] / u22;
    x2[0] = btmp[0] / u11 - u12 / u11 * x2[1];
    if (xswap) {
	temp = x2[1];
	x2[1] = x2[0];
	x2[0] = temp;
    }
    x[ldx + 1] = x2[0];
    if (n1 == 1) {
	x[(ldx << 1) + 1] = x2[1];
	*xnorm = abs(x[ldx + 1]) + abs(x[(ldx << 1) + 1]);
    } else {
	x[ldx + 2] = x2[1];
	mtemp1 = abs(x[ldx + 1]), mtemp2 = abs(x[ldx + 2]);
	*xnorm = max(mtemp1, mtemp2);
    }
    return;
//2 by 2:
//op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
//  [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
//Solve equivalent 4 by 4 system using complete pivoting.
//Set pivots less than SMIN to SMIN.
  L50:
    mtemp1 = abs(tr[ldtr + 1]);
    mtemp2 = abs(tr[(ldtr << 1) + 1]);
    mtemp1 = max(mtemp1, mtemp2);
    mtemp2 = abs(tr[ldtr + 2]);
    mtemp1 = max(mtemp1, mtemp2);
    mtemp2 = abs(tr[(ldtr << 1) + 2]);
    smin = max(mtemp1, mtemp2);

    mtemp1 = smin;
    mtemp2 = abs(tl[ldtl + 1]), mtemp1 = max(mtemp1, mtemp2);
    mtemp2 = abs(tl[(ldtl << 1) + 1]), mtemp1 = max(mtemp1, mtemp2);
    mtemp2 = abs(tl[ldtl + 2]), mtemp1 = max(mtemp1, mtemp2);
    mtemp2 = abs(tl[(ldtl << 1) + 2]);
    smin = max(mtemp1, mtemp2);

    mtemp1 = eps * smin;
    smin = max(mtemp1, smlnum);
    btmp[0] = Zero;
    Rcopy(16, btmp, 0, t16, 1);
    t16[0] = tl[ldtl + 1] + sgn * tr[ldtr + 1];
    t16[5] = tl[(ldtl << 1) + 2] + sgn * tr[ldtr + 1];
    t16[10] = tl[ldtl + 1] + sgn * tr[(ldtr << 1) + 2];
    t16[15] = tl[(ldtl << 1) + 2] + sgn * tr[(ldtr << 1) + 2];
    if (ltranl) {
	t16[4] = tl[ldtl + 2];
	t16[1] = tl[(ldtl << 1) + 1];
	t16[14] = tl[ldtl + 2];
	t16[11] = tl[(ldtl << 1) + 1];
    } else {
	t16[4] = tl[(ldtl << 1) + 1];
	t16[1] = tl[ldtl + 2];
	t16[14] = tl[(ldtl << 1) + 1];
	t16[11] = tl[ldtl + 2];
    }
    if (ltranr) {
	t16[8] = sgn * tr[(ldtr << 1) + 1];
	t16[13] = sgn * tr[(ldtr << 1) + 1];
	t16[2] = sgn * tr[ldtr + 2];
	t16[7] = sgn * tr[ldtr + 2];
    } else {
	t16[8] = sgn * tr[ldtr + 2];
	t16[13] = sgn * tr[ldtr + 2];
	t16[2] = sgn * tr[(ldtr << 1) + 1];
	t16[7] = sgn * tr[(ldtr << 1) + 1];
    }
    btmp[0] = B[ldb + 1];
    btmp[1] = B[ldb + 2];
    btmp[2] = B[(ldb << 1) + 1];
    btmp[3] = B[(ldb << 1) + 2];
//Perform elimination
    for (i = 0; i < 3; i++) {
	xmax = Zero;
	for (ip = i; ip <= 4; ip++) {
	    for (jp = i; jp <= 4; jp++) {
		if (abs(t16[ip + (jp << 2) - 5]) >= xmax) {
		    xmax = abs(t16[ip + (jp << 2) - 5]);
		    ipsv = ip;
		    jpsv = jp;
		}
	    }
	}
	if (ipsv != i) {
	    Rswap(4, &t16[ipsv - 1], 4, &t16[i - 1], 4);
	    temp = btmp[i - 1];
	    btmp[i - 1] = btmp[ipsv - 1];
	    btmp[ipsv - 1] = temp;
	}
	if (jpsv != i) {
	    Rswap(4, &t16[(jpsv << 2) - 4], 1, &t16[(i << 2) - 4], 1);
	}
	jpiv[i - 1] = jpsv;
	if (abs(t16[i + (i << 2) - 5]) < smin) {
	    *info = 1;
	    t16[i + (i << 2) - 5] = smin;
	}
	for (j = i + 1; j <= 4; j++) {
	    t16[j + (i << 2) - 5] /= t16[i + (i << 2) - 5];
	    btmp[j - 1] -= t16[j + (i << 2) - 5] * btmp[i - 1];
	    for (k = i + 1; k <= 4; k++) {
		t16[j + (k << 2) - 5] -= t16[j + (i << 2) - 5] * t16[i + (k << 2) - 5];
	    }
	}
    }
    if (abs(t16[15]) < smin) {
	t16[15] = smin;
    }
    *scale = One;
    if (smlnum * Eight * abs(btmp[0]) > abs(t16[0])
	|| smlnum * Eight * abs(btmp[1]) > abs(t16[5])
	|| smlnum * Eight * abs(btmp[2]) > abs(t16[10])
	|| smlnum * Eight * abs(btmp[3]) > abs(t16[15])) {
	mtemp1 = abs(btmp[0]), mtemp2 = abs(btmp[1]), mtemp1 = max(mtemp1, mtemp2);
	mtemp2 = abs(btmp[2]), mtemp1 = max(mtemp1, mtemp2);
	mtemp2 = abs(btmp[3]);
	*scale = Fourth / max(mtemp1, mtemp2);
	btmp[0] = btmp[0] * (*scale);
	btmp[1] = btmp[1] * (*scale);
	btmp[2] = btmp[2] * (*scale);
	btmp[3] = btmp[3] * (*scale);
    }
    for (i = 0; i < 4; i++) {
	k = 5 - i;
	temp = One / t16[k + (k << 2) - 5];
	tmp[k - 1] = btmp[k - 1] * temp;
	for (j = k + 1; j <= 4; j++) {
	    tmp[k - 1] -= temp * t16[k + (j << 2) - 5] * tmp[j - 1];
	}
    }
    for (i = 0; i < 3; i++) {
	if (jpiv[4 - i - 1] != 4 - i) {
	    temp = tmp[4 - i - 1];
	    tmp[4 - i - 1] = tmp[jpiv[4 - i - 1] - 1];
	    tmp[jpiv[4 - i - 1] - 1] = temp;
	}
    }
    x[ldx + 1] = tmp[0];
    x[ldx + 2] = tmp[1];
    x[(ldx << 1) + 1] = tmp[2];
    x[(ldx << 1) + 2] = tmp[3];
    mtemp1 = abs(tmp[0]) + abs(tmp[2]), mtemp2 = abs(tmp[1]) + abs(tmp[3]);
    *xnorm = max(mtemp1, mtemp2);
#endif
    return;
}
