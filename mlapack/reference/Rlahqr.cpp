/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlahqr.cpp,v 1.12 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rlahqr(INTEGER wantt, INTEGER wantz, INTEGER n, INTEGER ilo, INTEGER ihi, REAL * h, INTEGER ldh,
       REAL * wr, REAL * wi, INTEGER iloz, INTEGER ihiz, REAL * z, INTEGER ldz, INTEGER * info)
{
    INTEGER i, j, k, l, m = 0;
    REAL s = 0.0, v[3];
    INTEGER i1, i2 = 0;
    REAL t1, t2, t3, v2, v3, aa, ab, ba, bb, h11, h12, h21, h22, cs;
    INTEGER nh;
    REAL sn;
    INTEGER nr;
    REAL tr;
    INTEGER nz;
    REAL det, h21s;
    INTEGER its;
    REAL ulp, sum, tst, rt1i, rt2i, rt1r, rt2r;
    REAL safmin, safmax, rtdisc, smlnum;
    REAL Zero = 0.0, One = 1.0, Two = 2.0, Dat1 = 3.0 / 4.0, Dat2 = -0.4375;
    REAL mtemp1, mtemp2;

    *info = 0;
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (ilo == ihi) {
	wr[ilo] = h[ilo + ilo * ldh];
	wi[ilo] = Zero;
	return;
    }
//==== clear out the trash ====
    for (j = ilo; j <= ihi - 3; j++) {
	h[j + 2 + j * ldh] = Zero;
	h[j + 3 + j * ldh] = Zero;

    }
    if (ilo <= ihi - 2) {
	h[ihi + (ihi - 2) * ldh] = Zero;
    }
    nh = ihi - ilo + 1;
    nz = ihiz - iloz + 1;
//Set machine-dependent constants for the stopping criterion.
    safmin = Rlamch("SAFE MINIMUM");
    safmax = One / safmin;
    ulp = Rlamch("PRECISION");
    smlnum = safmin * ((REAL) double (nh) / ulp);
//I1 and I2 are the indices of the first row and last column of H
//to which transformations must be applied. If eigenvalues only are
//being computed, I1 and I2 are set inside the main loop.
    if (wantt) {
	i1 = 1;
	i2 = n;
    }
//The main loop begins here. I is the loop index and decreases from
//IHI to ILO in steps of 1 or Two Each iteration of the loop works
//with the active submatrix in rows and columns L to I.
//Eigenvalues I+1 to IHI have already converged. Either L = ILO or
//H(L,L-1) is negligible so that the matrix splits.
    i = ihi;
  L20:
    l = ilo;
    if (i < ilo) {
	goto L160;
    }
//Perform QR iterations on rows and columns ILO to I until a
//submatrix of order 1 or 2 splits off at the bottom because a
//subdiagonal element has become negligible.
    for (its = 0; its <= 30; its++) {
//Look for a single small subdiagonal element.
	for (k = i; k >= 0; k--) {
	    if (abs(h[k + (k - 1) * ldh]) <= smlnum) {
		goto L40;
	    }
	    tst = abs(h[k - 1 + (k - 1) * ldh]) + abs(h[k + k * ldh]);
	    if (tst == Zero) {
		if (k - 2 >= ilo) {
		    tst += abs(h[k - 1 + (k - 2) * ldh]);
		}
		if (k + 1 <= ihi) {
		    tst += abs(h[k + 1 + k * ldh]);
		}
//==== The following is a conservative small subdiagonal
//.    deflation  criterion due to Ahues & Tisseur (LAWN 122,
//.    1997). It has better mathematical foundation and
//.    improves accuracy in some cases.  ====
		if (abs(h[k + (k - 1) * ldh]) <= ulp * tst) {
		    mtemp1 = abs(h[k + (k - 1) * ldh]), mtemp2 = abs(h[k - 1 + k * ldh]);
		    ab = max(mtemp1, mtemp2);
		    mtemp1 = abs(h[k + (k - 1) * ldh]), mtemp2 = abs(h[k - 1 + k * ldh]);
		    ba = min(mtemp1, mtemp2);
		    mtemp1 = abs(h[k + k * ldh]), mtemp2 = abs(h[k - 1 + (k - 1) * ldh] - h[k + k * ldh]);
		    aa = max(mtemp1, mtemp2);
		    mtemp1 = abs(h[k + k * ldh]), mtemp2 = abs(h[k - 1 + (k - 1) * ldh] - h[k + k * ldh]);
		    bb = min(mtemp1, mtemp2);
		    s = aa + ab;
		    mtemp1 = smlnum, mtemp2 = ulp * (bb * (aa / s));
		    if (ba * (ab / s) <= max(mtemp1, mtemp2)) {
			goto L40;
		    }
		}
	    }
	  L40:
	    l = k;
	    if (l > ilo) {
//H(L,L-1) is negligible
		h[l + (l - 1) * ldh] = Zero;
	    }
//Exit from loop if a submatrix of order 1 or 2 has split off.
	    if (l >= i - 1) {
		goto L150;
	    }
//Now the active submatrix is in rows and columns L to I. If
//eigenvalues only are being computed, only the active submatrix
//need be transformed.
	    if (!(wantt)) {
		i1 = l;
		i2 = i;
	    }
	    if (its == 10 || its == 20) {
//Exceptional shift.
		h11 = s * Dat1 + h[i + i * ldh];
		h12 = s * Dat2;
		h21 = s;
		h22 = h11;
	    } else {
//Prepare to use Francis' double shift
//(i.e. 2nd degree generalized Rayleigh quotient)
		h11 = h[i - 1 + (i - 1) * ldh];
		h21 = h[i + (i - 1) * ldh];
		h12 = h[i - 1 + i * ldh];
		h22 = h[i + i * ldh];
	    }
	    s = abs(h11) + abs(h12) + abs(h21) + abs(h22);
	    if (s == Zero) {
		rt1r = Zero;
		rt1i = Zero;
		rt2r = Zero;
		rt2i = Zero;
	    } else {
		h11 = h11 / s;
		h21 = h21 / s;
		h12 = h12 / s;
		h22 = h22 / s;
		tr = (h11 + h22) / Two;
		det = (h11 - tr) * (h22 - tr) - h12 * h21;
		rtdisc = sqrt((abs(det)));
		if (det >= Zero) {
//==== complex conjugate shifts ====
		    rt1r = tr * s;
		    rt2r = rt1r;
		    rt1i = rtdisc * s;
		    rt2i = -rt1i;
		} else {
//==== real shifts (use only one of them)  ====
		    rt1r = tr + rtdisc;
		    rt2r = tr - rtdisc;
		    if (abs(rt1r - h22) <= abs(rt2r - h22)) {
			rt1r = rt1r * s;
			rt2r = rt1r;
		    } else {
			rt2r = rt2r * s;
			rt1r = rt2r;
		    }
		    rt1i = Zero;
		    rt2i = Zero;
		}
	    }
//Look for two consecutive small subdiagonal elements.
	    for (m = i - 2; m >= l; m--) {
//Determine the effect of starting the double-shift QR
//iteration at row M, and see if this would make H(M,M-1)
//negligible.  (The following uses scaling to avoid
//overflows and most underflows.)

		h21s = h[m + 1 + m * ldh];
		s = abs(h[m + m * ldh] - rt2r) + abs(rt2i) + abs(h21s);
		h21s = h[m + 1 + m * ldh] / s;
		v[0] = h21s * h[m + (m + 1) * ldh] + (h[m + m * ldh] - rt1r) * ((h[m + m * ldh] - rt2r) / s) - rt1i * (rt2i / s);
		v[1] = h21s * (h[m + m * ldh] + h[m + 1 + (m + 1) * ldh]
			       - rt1r - rt2r);
		v[2] = h21s * h[m + 2 + (m + 1) * ldh];
		s = abs(v[0]) + abs(v[1]) + abs(v[2]);
		v[0] = v[0] / s;
		v[1] = v[1] / s;
		v[2] = v[2] / s;
		if (m == l) {
		    goto L60;
		}
		if (abs(h[m + (m - 1) * ldh]) * (abs(v[1]) + abs(v[2])) <= ulp * abs(v[0]) * abs(h[m - 1 + (m - 1)]) + abs(h[m + m * ldh]) + abs(h[m + 1 + (m + 1) * ldh]))
		    goto L60;
	    }

	}
      L60:
//Double-shift QR step
	for (k = m; k <= i - 1; k++) {
//The first iteration of this loop determines a reflection G
//from the vector V and applies it from left and right to H,
//thus creating a nonzero bulge below the subdiagonal.
//Each subsequent iteration determines a reflection G to
//restore the Hessenberg form in the (K-1)th column, and thus
//chases the bulge one step toward the bottom of the active
//submatrix. NR is the order of G.
	    nr = min((INTEGER) 3, i - k + 1);
	    if (k > m) {
		Rcopy(nr, &h[k + (k - 1) * ldh], 1, v, 1);
	    }
	    Rlarfg(nr, v, &v[1], 1, &t1);
	    if (k > m) {
		h[k + (k - 1) * ldh] = v[0];
		h[k + 1 + (k - 1) * ldh] = Zero;
		if (k < i - 1) {
		    h[k + 2 + (k - 1) * ldh] = Zero;
		}
	    } else if (m > l) {
		h[k + (k - 1) * ldh] = -h[k + (k - 1) * ldh];
	    }
	    v2 = v[1];
	    t2 = t1 * v2;
	    if (nr == 3) {
		v3 = v[2];
		t3 = t1 * v3;
//Apply G from the left to transform the rows of the matrix
//in columns K to I2
		for (j = k; j <= i2; j++) {
		    sum = h[k + j * ldh] + v2 * h[k + 1 + j * ldh]
			+ v3 * h[k + 2 + j * ldh];
		    h[k + j * ldh] = h[k + j * ldh] - sum * t1;
		    h[k + 1 + j * ldh] = h[k + 1 + j * ldh] - sum * t2;
		    h[k + 2 + j * ldh] = h[k + 2 + j * ldh] - sum * t3;
		}
//Apply G from the right to transform the columns of the
//matrix in rows I1 to min(K+3,I).
		i2 = min(k + 3, i);
		for (j = i1; j <= i2; j++) {
		    sum = h[j + k * ldh] + v2 * h[j + (k + 1) * ldh]
			+ v3 * h[j + (k + 2) * ldh];
		    h[j + k * ldh] = h[j + k * ldh] - sum * t1;
		    h[j + (k + 1) * ldh] = h[j + (k + 1) * ldh] - sum * t2;
		    h[j + (k + 2) * ldh] = h[j + (k + 2) * ldh] - sum * t3;
		}
		if (wantz) {
//Accumulate transformations in the matrix Z
		    for (j = iloz; j <= ihiz; j++) {
			sum = z[j + k * ldz] + v2 * z[j + (k + 1) * ldz] + v3 * z[j + (k + 2) * ldz];
			z[j + k * ldz] = z[j + k * ldz] - sum * t1;
			z[j + (k + 1) * ldz] = z[j + (k + 1) * ldz] - sum * t2;
			z[j + (k + 2) * ldz] = z[j + (k + 2) * ldz] - sum * t3;
		    }
		}
	    } else if (nr == 2) {
//Apply G from the left to transform the rows of the matrix
//in columns K to I2
		for (j = k; j <= i2; j++) {
		    sum = h[k + j * ldh] + v2 * h[k + 1 + j * ldh];
		    h[k + j * ldh] = h[k + j * ldh] - sum * t1;
		    h[k + 1 + j * ldh] = h[k + 1 + j * ldh] - sum * t2;
		}
//Apply G from the right to transform the columns of the
//matrix in rows I1 to min(K+3,I).
		for (j = i1; j <= i; j++) {
		    sum = h[j + k * ldh] + v2 * h[j + (k + 1) * ldh];
		    h[j + k * ldh] = h[j + k * ldh] - sum * t1;
		    h[j + (k + 1) * ldh] = h[j + (k + 1) * ldh] - sum * t2;
		}
		if (wantz) {
//Accumulate transformations in the matrix Z
		    for (j = iloz; j <= ihiz; j++) {
			sum = z[j + k * ldz] + v2 * z[j + (k + 1) * ldz];
			z[j + k * ldz] = z[j + k * ldz] - sum * t1;
			z[j + (k + 1) * ldz] = z[j + (k + 1) * ldz] - sum * t2;
		    }
		}
	    }
	}
    }
//Failure to converge in remaining number of iterations
    *info = i;
    return;
  L150:
    if (l == i) {
//H(I,I-1) is negligible: one eigenvalue has converged.
	wr[i] = h[i + i * ldh];
	wi[i] = Zero;
    } else if (l == i - 1) {

//H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
//Transform the 2-by-2 submatrix to standard Schur form, 
//and compute and store the eigenvalues.
	Rlanv2(&h[i - 1 + (i - 1) * ldh], &h[i - 1 + i * ldh], &h[i + (i - 1) * ldh], &h[i + i * ldh], &wr[i - 1], &wi[i - 1], &wr[i], &wi[i], &cs, &sn);
	if (wantt) {
//Apply the transformation to the rest of H.
	    if (i2 > i) {
		Rrot(i2 - i, &h[i - 1 + (i + 1) * ldh], ldh, &h[i + (i + 1) * ldh], ldh, cs, sn);
	    }
	    Rrot(i - i1 - 1, &h[i1 + (i - 1) * ldh], 1, &h[i1 + i * ldh], 1, cs, sn);
	}
	if (wantz) {
//Apply the transformation to Z.
	    Rrot(nz, &z[iloz + (i - 1) * ldz], 1, &z[iloz + i * ldz], 1, cs, sn);
	}
    }
//return to start of the main loop with new value of I.
    i = l - 1;
    goto L20;

  L160:
    return;
}
