/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clahqr.cpp,v 1.10 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Clahqr(INTEGER wantt, INTEGER wantz, INTEGER n, INTEGER ilo, INTEGER ihi, COMPLEX * h, INTEGER ldh,
       COMPLEX * w, INTEGER iloz, INTEGER ihiz, COMPLEX * z, INTEGER ldz, INTEGER * info)
{
    INTEGER i, j, k, l, m;
    REAL s;
    COMPLEX t, u, v[2], x, y;
    INTEGER i1, i2 = 0;
    COMPLEX t1;
    REAL t2;
    COMPLEX v2;
    REAL aa, ab, ba, bb, h10;
    COMPLEX h11;
    COMPLEX h21;
    COMPLEX h22, sc;
    INTEGER nh, nz;
    REAL sx;
    INTEGER jhi;
    COMPLEX h11s;
    INTEGER jlo, its;
    REAL ulp;
    COMPLEX sum;
    REAL tst;
    COMPLEX temp;
    REAL rtemp;
    REAL safmin, safmax;
    REAL smlnum;
    REAL Zero = 0.0, Dat1 = 0.75, One = 1.0, Half = 0.5;
    REAL mtemp1, mtemp2;

    *info = 0;
//Quick return if possible
    if (n == 0) {
	return;
    }
    if (ilo == ihi) {
	i1 = ilo;
	i2 = ilo + ilo * ldh;
	w[i1] = h[i2];
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
//==== ensure that subdiagonal entries are real ====
    for (i = ilo + 1; i <= ihi; i++) {
	if (h[i + (i - 1) * ldh].imag() != Zero) {
//==== The following redundant normalization
//.    avoids problems with both gradual and
//.    sudden underflow in ABS(H(I,I-1)) ====
	    sc = h[i + (i - 1) * ldh] / RCabs1(h[i + (i - 1) * ldh]);
	    sc = conj(sc) / abs(sc);
	    h[i + (i - 1) * ldh] = abs(h[i + (i - 1) * ldh]);
	    if (wantt) {
		jlo = 1;
		jhi = n;
	    } else {
		jlo = ilo;
		jhi = ihi;
	    }
	    Cscal(jhi - i + 1, sc, &h[i + i * ldh], ldh);
	    Cscal(min(jhi, i + 1) - jlo + 1, conj(sc), &h[jlo + i * ldh], 1);
	    if (wantz) {
		Cscal(ihiz - iloz + 1, conj(sc), &z[iloz + i * ldz], 1);
	    }
	}
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
//IHI to ILO in steps of One Each iteration of the loop works
//with the active submatrix in rows and columns L to I.
//Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
//H(L,L-1) is negligible so that the matrix splits.

    i = ihi;
  L30:
    if (i < ilo) {
	goto L150;
    }
//Perform QR iterations on rows and columns ILO to I until a
//submatrix of order 1 splits off at the bottom because a
//subdiagonal element has become negligible.
    l = ilo;
    for (its = 0; its <= 30; its++) {
//Look for a single small subdiagonal element.
	for (k = i; k >= l + 1; k--) {
	    i2 = k + (k - 1) * ldh;
	    if (RCabs1(h[k + (k - 1) * ldh]) <= smlnum) {
		goto L50;
	    }
	    tst = RCabs1(h[k - 1 + (k - 1) * ldh]) + RCabs1(h[k + k * ldh]);
	    if (tst == Zero) {
		if (k - 2 >= ilo) {
		    tst = tst + abs(h[k - 1 + (k - 2) * ldh].real());
		}
		if (k + 1 <= ihi) {
		    tst = tst + abs(h[k + 1 + k * ldh].real());
		}
	    }
//==== The following is a conservative small subdiagonal
//.    deflation criterion due to Ahues & Tisseur (LAWN 122,
//.    1997). It has better mathematical foundation and
//.    improves accuracy in some examples.  ====
	    if (abs(h[k + (k - 1) * ldh].real()) <= ulp * tst) {
		mtemp1 = RCabs1(h[k + (k - 1) * ldh]), mtemp2 = RCabs1(h[k - 1 + k * ldh]);
		ab = max(mtemp1, mtemp2);
		ba = min(mtemp1, mtemp2);
		mtemp1 = RCabs1(h[k + k * ldh]);
		mtemp2 = RCabs1(h[(k - 1) + (k - 1) * ldh] - h[k + k * ldh]);
		aa = max(mtemp1, mtemp2);
		bb = min(mtemp1, mtemp2);
		s = aa + ab;
		mtemp1 = smlnum, mtemp2 = ulp * (bb * (aa / s));
		if (ba * (ab / s) <= max(mtemp1, mtemp2)) {
		    goto L50;
		}
	    }

	}
      L50:
	l = k;
	if (l > ilo) {
//H(L,L-1) is negligible
	    h[l + (l - 1) * ldh] = Zero;
	}
//Exit from loop if a submatrix of order 1 has split off.
	if (l >= i) {
	    goto L140;
	}
//Now the active submatrix is in rows and columns L to I. If
//eigenvalues only are being computed, only the active submatrix
//need be transformed.
	if (!(wantt)) {
	    i1 = l;
	    i2 = i;
	}
	if (its == 10 || its == 20) {
	    s = Dat1 * abs(h[i + (i - 1) * ldh].real());
	    t = s + h[i + i * ldh];
	} else {
//Wilkinson's shift.
	    t = h[i + i * ldh];
	    u = sqrt(h[i - 1 + i * ldh]) + sqrt(h[i + (i - 1) * ldh]);
	    s = RCabs1(u);
	    if (s != Zero) {
		x = Half * (h[i - 1 + (i - 1) * ldh] - t);
		sx = RCabs1(abs(x));
		mtemp1 = s, mtemp2 = RCabs1(x);
		s = max(mtemp1, mtemp2);
		y = s * sqrt((x / s) * (x / s) + (u / s) * (u / s));
		if (sx > Zero) {
		    if ((x / sx).real() * (y).real() + (x / sx).imag() * (y).imag() < Zero)
			y = -y;
		}
	    }
	    t = t - u * Cladiv(u, (x + y));
	}
    }
//Look for two consecutive small subdiagonal elements.
    for (m = i - 1; m >= l + 1; m--) {
//Determine the effect of starting the single-shift QR
//iteration at row M, and see if this would make H(M,M-1)
//negligible.
	h11 = h[m + m * ldh];
	h22 = h[m + 1 + (m + 1) * ldh];
	h11s = h11 - t;
	h21 = h[m + 1 + m * ldh].real();
	s = RCabs1(h11s) + abs(h21);
	h11s = h11s / s;
	h21 = h21 / s;
	v[1] = h11s;
	v[2] = h21;
	h10 = h[m + (m - 1) * ldh].real();
	if (abs(h10) * abs(h21) <= ulp * (RCabs1(h11s) * (RCabs1(h11) + RCabs1(h22))))
	    goto L70;
    }

    h11 = h[l + l * ldh];
    h22 = h[l + 1 + (l + 1) * ldh];
    h11s = h11 - t;
    h21 = h[l + 1 + l * ldh];
    s = RCabs1(h11s) + abs(h21);
    h11s = h11s / s;
    h21 = h21 / s;
    v[1] = h11s;
    v[2] = h21;

  L70:

//Single-shift QR step
    for (k = m; k <= i - 1; k++) {
//The first iteration of this loop determines a reflection G
//from the vector V and applies it from left and right to H,
//thus creating a nonzero bulge below the subdiagonal.
//Each subsequent iteration determines a reflection G to
//restore the Hessenberg form in the (K-1)th column, and thus
//chases the bulge one step toward the bottom of the active
//submatrix.
//V(2) is always real before the call to ZLARFG, and hence
//after the call T2 ( = T1*V(2) ) is also real.
	if (k > m) {
	    Ccopy(2, &h[k + (k - 1) * ldh], 1, v, 1);
	}
	Clarfg(2, v, &v[1], 1, &t1);
	if (k > m) {
	    h[k + (k - 1) * ldh] = v[1];
	    h[k + 1 + (k - 1) * ldh] = Zero;
	}
	v2 = v[1];
	t2 = (t1 * v2).real();
//Apply G from the left to transform the rows of the matrix
//in columns K to I2
	for (j = k; j <= i2; j++) {
	    sum = conj(t1) * h[k + j * ldh] + t2 * h[k + 1 + j * ldh];
	    h[k + j * ldh] = h[k + j * ldh] - sum;
	    h[k + 1 + j * ldh] = h[k + 1 + j * ldh] - sum * v2;
	}
//Apply G from the right to transform the columns of the
//matrix in rows I1 to min(K+2,I).
	for (j = i1; j <= min(k + 2, i); j++) {
	    sum = t1 * h[j + k * ldh] + t2 * h[j + (k + 1) * ldh];
	    h[j + k * ldh] = h[j + k * ldh] - sum;
	    h[j + (k + 1) * ldh] = h[j + (k + 1) * ldh] - sum * conj(v2);
	}
	if (wantz) {
//Accumulate transformations in the matrix Z
	    for (j = iloz; j <= ihiz; j++) {
		sum = t1 * z[j + k * ldz] + t2 * z[j + (k + 1) * ldz];
		z[j + k * ldz] = z[j + k * ldz] - sum;
		z[j + (k + 1) * ldz] = z[j + (k + 1) * ldz] - sum * conj(v2);
	    }
	}
	if (k == m && m > l) {
//If the QR step was started at row M > L because two
//consecutive small subdiagonals were found, then extra
//scaling must be performed to ensure that H(M,M-1) remains
	    temp = One - t1;
	    temp = temp / abs(temp);
	    h[m + 1 + m * ldh] = h[m + 1 + m * ldh] * conj(temp);
	    if (m + 2 <= i) {
		h[m + 2 + (m + 1) * ldh] = h[m + 2 + (m + 1) * ldh] * temp;
	    }
	    for (j = m; j <= i; j++) {
		if (j != m + 1) {
		    if (i2 > j) {
			Cscal(i2 - j, temp, &h[j + (j + 1) * ldh], ldh);
			Cscal(j - i1, conj(temp), &h[i1 + j * ldh], 1);
		    }
		    if (wantz) {
			Cscal(nz, conj(temp), &z[iloz + j * ldz], 1);
		    }
		}
	    }
	}
    }
//Ensure that H(I,I-1) is real.
    temp = h[i + (i - 1) * ldh];
    if (temp.imag() != Zero) {
	rtemp = abs(temp);
	h[i + (i - 1) * ldh] = rtemp;
	temp = temp / rtemp;
	if (i2 > i) {
	    Cscal(i2 - i, conj(temp), &h[i + (i + 1) * ldh], ldh);
	}
	Cscal(i - i1, temp, &h[i1 + i * ldh], 1);
	if (wantz) {
	    Cscal(nz, temp, &z[iloz + i * ldz], 1);
	}
    }
//Failure to converge in remaining number of iterations
    *info = i;
    return;
  L140:

//H(I,I-1) is negligible: one eigenvalue has converged.
    w[i] = h[i + i * ldh];
//return to start of the main loop with new value of I.
    i = l - 1;
    goto L30;
  L150:
    return;
}
