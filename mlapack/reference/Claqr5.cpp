/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Claqr5.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Claqr5(LOGICAL wantt, LOGICAL wantz, INTEGER kacc22,
	    INTEGER n, INTEGER ktop, INTEGER kbot, INTEGER nshfts,
	    COMPLEX * s, COMPLEX * h, INTEGER ldh, INTEGER iloz, INTEGER ihiz, COMPLEX * z, INTEGER ldz, COMPLEX * v, INTEGER ldv,
	    COMPLEX * u, INTEGER ldu, INTEGER nv, COMPLEX * wv, INTEGER ldwv, INTEGER nh, COMPLEX * wh, INTEGER ldwh)
{
    INTEGER j, k, m, i2, j2, i4, j4, k1;
    REAL h11, h12, h21, h22;
    INTEGER m22, ns, nu;
    COMPLEX vt[3];
    REAL scl;
    INTEGER kdu, kms;
    REAL ulp;
    INTEGER knz, kzs;
    REAL tst1, tst2;
    COMPLEX beta;
    INTEGER blk22, bmp22;
    INTEGER mend, jcol, jlen, jbot, mbot, jtop, jrow, mtop;
    COMPLEX alpha;
    INTEGER accum;
    INTEGER ndcol, incol, krcol, nbmps;
    REAL safmin, safmax;
    COMPLEX refsum;
    INTEGER mstart;
    REAL smlnum;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;

//==== If there are no shifts, then there is nothing to do. ====
    if (nshfts < 2) {
	return;
    }
//==== If the active block is empty or 1-by-1, then there
//.    is nothing to do. ====
    if (ktop >= kbot) {
	return;
    }
//==== NSHFTS is supposed to be even, but if is odd,
//.    then simply reduce it by one.  ====
    ns = nshfts - nshfts % 2;
//==== Machine constants for deflation ====
    safmin = Rlamch("SAFE MINIMUM");
    safmax = One / safmin;
    ulp = Rlamch("PRECISION");
    smlnum = safmin * ((double) (n) / ulp);
//==== Use accumulated reflections to update far-from-diagonal
//.    entries ? ====
    accum = kacc22 == 1 || kacc22 == 2;
//==== If so, exploit the 2-by-2 block structure? ====
    blk22 = ns > 2 && kacc22 == 2;
//==== clear trash ====
    if (ktop + 2 <= kbot) {
	h[ktop + 2 + ktop * ldh] = Zero;
    }
//==== NBMPS = number of 2-shift bulges in the chain ====
    nbmps = ns / 2;
//==== KDU = width of slab ====
    kdu = nbmps * 6 - 3;
//==== Create and chase chains of NBMPS bulges ====
    for (incol = (1 - nbmps) * 3 + ktop - 1; incol <= kbot - 2; incol = incol + nbmps * 3 - 2) {
	ndcol = incol + kdu;
	if (accum) {
	    Claset("ALL", kdu, kdu, Zero, One, &u[0], ldu);
	}
//==== Near-the-diagonal bulge chase.  The following loop
//.    performs the near-the-diagonal part of a small bulge
//.    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal
//.    chunk extends from column INCOL to column NDCOL
//.    (including both column INCOL and column NDCOL). The
//.    following loop chases a 3*NBMPS column long chain of
//.    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL
//.    may be less than KTOP and and NDCOL may be greater than
//.    KBOT indicating phantom columns from which to chase
//.    bulges before they are actually introduced or to which
//.    to chase bulges beyond column KBOT.)  ====
	for (krcol = incol; krcol <= min(incol + nbmps * 3 - 3, kbot - 2); krcol++) {
//==== Bulges number MTOP to MBOT are active double implicit
//.    shift bulges.  There may or may not also be small
//.    2-by-2 bulge, if there is room.  The inactive bulges
//.    (if any) must wait until the active bulges have moved
//.    down the diagonal to make room.  The phantom matrix
//.    paradigm described above helps keep track.  ====
	    mtop = max((INTEGER) 1, (ktop - 1 - krcol + 2) / 3 + 1);
	    mbot = min(nbmps, (kbot - krcol) / 3);
	    m22 = mbot + 1;
	    bmp22 = mbot < nbmps && krcol + (m22 - 1) * 3 == kbot - 2;
//==== Generate reflections to chase the chain right
//.    one column.  (The minimum value of K is KTOP-One) ====
	    for (m = mtop; m <= mbot; m++) {
		k = krcol + (m - 1) * 3;
		if (k == ktop - 1) {
		    Claqr1(3, &h[ktop + ktop * ldh], ldh, s[(m * 2) - 1], s[m * 2], &v[m * ldv + 1]);
		    alpha = v[m * ldv + 1];
		    Clarfg(3, &alpha, &v[m * ldv + 2], 1, &v[m * ldv + 1]);
		} else {
		    beta = h[k + 1 + k * ldh];
		    v[m * ldv + 2] = h[k + 2 + k * ldh];
		    v[m * ldv + 3] = h[k + 3 + k * ldh];
		    Clarfg(3, &beta, &v[m * ldv + 2], 1, &v[m * ldv + 1]);
//==== A Bulge may collapse because of vigilant
//.    deflation or destructive underflow.  (The
//.    initial bulge is always collapsed.) Use
//.    the two-small-subdiagonals trick to try
//.    to get it started again. If V(2,M).NE.0 and
//.    V(3,M) = H(K+3,K+1) = H(K+3,K+2) = 0, then
//.    this bulge is collapsing into a zero
//.    subdiagonal.  It will be restarted next
//.    trip through the loop.)
		    if ((v[m * ldv + 1] != Zero && (v[m * ldv + 3] != Zero))
			|| (h[k + 3 + (k + 1) * ldh] == Zero && h[k + 3 + (k + 2) * ldh] == Zero)) {
//==== Typical case: not collapsed (yet). ====
			h[k + 1 + k * ldh] = beta;
			h[k + 2 + k * ldh] = Zero;
			h[k + 3 + k * ldh] = Zero;
		    } else {
//==== Atypical case: collapsed.  Attempt to
//.    reintroduce ignoring H(K+1,K).  If the
//.    fill resulting from the new reflector
//.    is too large, then abandon it.
//.    Otherwise, use the new one. ====
			Claqr1(3, &h[k + 1 + (k + 1) * ldh], ldh, s[(m * 2) - 1], s[m * 2], vt);
			scl = Cabs1(vt[0]) + Cabs1(vt[1]) + Cabs1(vt[2]);
			if (scl != Zero) {
			    vt[0] = vt[0] / scl;
			    vt[1] = vt[1] / scl;
			    vt[2] = vt[2] / scl;
			}
//==== The following is the traditional and
//.    conservative two-small-subdiagonals
//.    test.  ====
//.
			if ((Cabs1(h[k + 1 + k * ldh]) * Cabs1(vt[1])) + Cabs1(vt[2]) >
			    ulp * Cabs1(vt[0]) * Cabs1(h[k + k * ldh]) + Cabs1(h[k + 1 + (k + 1) * ldh]) + Cabs1(h[k + 2 + (k + 2) * ldh])) {
//==== Starting a new bulge here would
//.    create non-negligible fill.   If
//.    the old reflector is diagonal (only
//.    possible with underflows), then
//.    change it to I.  Otherwise, use
//.    it with trepidation. ====
			    if (v[m * ldv + 2] == Zero && v[m * ldv + 3] == Zero) {
				v[m * ldv + 1] = Zero;
			    } else {
				h[k + 1 + k * ldh] = beta;
				h[k + 2 + k * ldh] = Zero;
				h[k + 3 + k * ldh] = Zero;
			    }
			} else {
//==== Stating a new bulge here would
//.    create only negligible fill.
//.    Replace the old reflector with
//.    the new one. ====
			    alpha = vt[0];
			    Clarfg(3, &alpha, &vt[1], 1, vt);
			    refsum = h[k + 1 + k * ldh] + h[k + 2 + k * ldh] * conj(vt[1]) + h[k + 3 + k * ldh] * conj(vt[2]);
			    h[k + 1 + k * ldh] = h[k + 1 + k * ldh] - conj(*vt) * refsum;
			    h[k + 2 + k * ldh] = Zero;
			    h[k + 3 + k * ldh] = Zero;
			    v[m * ldv + 1] = vt[0];
			    v[m * ldv + 2] = vt[1];
			    v[m * ldv + 3] = vt[2];
			}
		    }
		}
	    }
//==== Generate a 2-by-2 reflection, if needed. ====
	    k = krcol + (m22 - 1) * 3;
	    if (bmp22) {
		if (k == ktop - 1) {
		    Claqr1(2, &h[k + 1 + (k + 1) * ldh], ldh, s[(m22 * 2) - 1], s[m22 * 2], &v[m22 * ldv + 1]);
		    beta = v[m22 * ldv + 1];
		    Clarfg(2, &beta, &v[m22 * ldv + 2], 1, &v[m22 * ldv + 1]);
		} else {
		    beta = h[k + 1 + k * ldh];
		    v[m22 * ldv + 2] = h[k + 2 + k * ldh];
		    Clarfg(2, &beta, &v[m22 * ldv + 2], 1, &v[m22 * ldv + 1]);
		    h[k + 1 + k * ldh] = beta;
		    h[k + 2 + k * ldh] = Zero;
		}
	    } else {
//==== Initialize V(1,M22) here to avoid possible undefined
//.    variable problems later. ====
		v[m22 * ldv + 1] = Zero;
	    }
//==== Multiply H by reflections from the left ====
	    if (accum) {
		jbot = min(ndcol, kbot);
	    } else if (wantt) {
		jbot = n;
	    } else {
		jbot = kbot;
	    }
	    for (j = max(ktop, krcol); j <= jbot; j++) {
		mend = min(mbot, (j - krcol + 2) / 3);
		for (m = mtop; m <= mend; m++) {
		    k = krcol + (m - 1) * 3;
		    refsum = conj(v[m * ldv + 1]) * (h[k + 1 + j * ldh] + conj(v[m * ldv + 2]) * h[k + 2 + j * ldh] + conj(v[m * ldv + 3]) * h[k + 3 + j * ldh]);
		    h[k + 1 + j * ldh] = h[k + 1 + j * ldh] - refsum;
		    h[k + 2 + j * ldh] = h[k + 2 + j * ldh] - refsum * v[m * ldv + 2];
		    h[k + 3 + j * ldh] = h[k + 3 + j * ldh] - refsum * v[m * ldv + 3];
		}
	    }
	    if (bmp22) {
		k = krcol + (m22 - 1) * 3;
		for (j = max(k + 1, ktop); j <= jbot; j++) {
		    refsum = conj(v[m22 * ldv + 1]) * (h[k + 1 + j * ldh] + conj(v[m22 * ldv + 2]) * h[k + 2 + j * ldh]);
		    h[k + 1 + j * ldh] = h[k + 1 + j * ldh] - refsum;
		    h[k + 2 + j * ldh] = h[k + 2 + j * ldh] - refsum * v[m22 * ldv + 2];
		}
	    }
//==== Multiply H by reflections from the right.
//.    Delay filling in the last row until the
//.    vigilant deflation check is complete. ====
	    if (accum) {
		jtop = max(ktop, incol);
	    } else if (wantt) {
		jtop = 1;
	    } else {
		jtop = ktop;
	    }
	    for (m = mtop; m <= mbot; m++) {
		if (v[m * ldv + 1] != Zero) {
		    k = krcol + (m - 1) * 3;
		    for (j = jtop; j <= min(kbot, k + 3); j++) {
			refsum = v[m * ldv + 1] * (h[j + (k + 1) * ldh] + v[m * ldv + 2] * h[j + (k + 2)
											     * ldh] + v[m * ldv + 3] * h[j + (k + 3) * ldh]);
			h[j + (k + 1) * ldh] = h[j + (k + 1) * ldh] - refsum;
			h[j + (k + 2) * ldh] = h[j + (k + 2) * ldh] - refsum * conj(v[m * ldv + 2]);
			h[j + (k + 3) * ldh] = h[j + (k + 3) * ldh] - refsum * conj(v[m * ldv + 3]);
		    }
		    if (accum) {
//==== Accumulate U. (If necessary, update Z later
//.    with with an efficient matrix-matrix
//.    multiply.) ====
			kms = k - incol;
			for (j = max((INTEGER) 1, ktop - incol); j <= kdu; j++) {
			    refsum = v[m * ldv + 1] * (u[j + (kms + 1) * ldu] + v[m * ldv + 2] * u[j + (kms + 2) * ldu] + v[m * ldv + 3] * u[j + (kms + 3) * ldu]);
			    u[j + (kms + 1) * ldu] = u[j + (kms + 1) * ldu] - refsum;
			    u[j + (kms + 2) * ldu] = u[j + (kms + 2) * ldu] - refsum * conj(v[m * ldv + 2]);
			    u[j + (kms + 3) * ldu] = u[j + (kms + 3) * ldu] - refsum * conj(v[m * ldv + 3]);
			}
		    } else if (wantz) {
//==== U is not accumulated, so update Z
//.    now by multiplying by reflections
//.    from the right. ====
			for (j = iloz; j <= ihiz; j++) {
			    refsum = v[m * ldv + 1] * (z[j + (k + 1) * ldz] + v[m * ldv + 2] * z[j + (k + 2) * ldz] + v[m * ldv + 3] * z[j + (k + 3) * ldz]);
			    z[j + (k + 1) * ldz] = z[j + (k + 1) * ldz] - refsum;
			    z[j + (k + 2) * ldz] = z[j + (k + 2) * ldz] - refsum * conj(v[m * ldv + 2]);
			    z[j + (k + 3) * ldz] = z[j + (k + 3) * ldz] - refsum * conj(v[m * ldv + 3]);
			}
		    }
		}
	    }
//==== Special case: 2-by-2 reflection (if needed) ====
	    k = krcol + (m22 - 1) * 3;
	    if (bmp22 && v[m22 * ldv + 1] != Zero) {
		for (j = jtop; j <= min(kbot, k + 3); j++) {
		    refsum = v[m22 * ldv + 1] * (h[j + (k + 1) * ldh]
						 + v[m22 * ldv + 2] * h[j + (k + 2) * ldh]);
		    h[j + (k + 1) * ldh] = h[j + (k + 1) * ldh] - refsum;
		    h[j + (k + 2) * ldh] = h[j + (k + 2) * ldh] - refsum * conj(v[m22 * ldv + 2]);
		}
		if (accum) {
		    kms = k - incol;
		    for (j = max((INTEGER) 1, ktop - incol); j <= kdu; j++) {
			refsum = v[m22 * ldv + 1] * (u[j + (kms + 1) * ldu] + v[m22 * ldv + 2] * u[j + (kms + 2) * ldu]);
			u[j + (kms + 1) * ldu] = u[j + (kms + 1) * ldu] - refsum;
			u[j + (kms + 2) * ldu] = u[j + (kms + 2) * ldu] - refsum * conj(v[m22 * ldv + 2]);
		    }
		} else if (wantz) {
		    for (j = iloz; j <= ihiz; j++) {
			refsum = v[m22 * ldv + 1] * (z[j + (k + 1) * ldz] + v[m22 * ldv + 2] * z[j + (k + 2) * ldz]);
			z[j + (k + 1) * ldz] = z[j + (k + 1) * ldz] - refsum;
			z[j + (k + 2) * ldz] = z[j + (k + 2) * ldz] - refsum * conj(v[m22 * ldv + 2]);
		    }
		}
	    }
//==== Vigilant deflation check ====
	    mstart = mtop;
	    if (krcol + (mstart - 1) * 3 < ktop) {
		mstart++;
	    }
	    mend = mbot;
	    if (bmp22) {
		mend++;
	    }
	    if (krcol == kbot - 2) {
		mend++;
	    }
	    for (m = mstart; m <= mend; m++) {
		k = min(kbot - 1, krcol + (m - 1) * 3);
//==== The following convergence test requires that
//.    the tradition small-compared-to-nearby-diagonals
//.    criterion and the Ahues & Tisseur (LAWN 122, 1997)
//.    criteria both be satisfied.  The latter improves
//.    accuracy in some examples. Falling back on an
//.    alternate convergence criterion when TST1 or TST2
//.    is zero (as done here) is traditional but probably
//.    unnecessary. ====
		if (h[k + 1 + k * ldh] != Zero) {
		    tst1 = Cabs1(h[k + k * ldh]) + Cabs1(h[k + 1 + (k + 1) * ldh]);
		    if (tst1 == Zero) {
			if (k >= ktop + 1) {
			    tst1 = tst1 + Cabs1(h[k + (k - 1) * ldh]);
			}
			if (k >= ktop + 2) {
			    tst1 = tst1 + Cabs1(h[k + (k - 2) * ldh]);
			}
			if (k >= ktop + 3) {
			    tst1 = tst1 + Cabs1(h[k + (k - 3) * ldh]);
			}
			if (k <= kbot - 2) {
			    tst1 = tst1 + Cabs1(h[k + 2 + (k + 1) * ldh]);
			}
			if (k <= kbot - 3) {
			    tst1 = tst1 + Cabs1(h[k + 3 + (k + 1) * ldh]);
			}
			if (k <= kbot - 4) {
			    tst1 = tst1 + Cabs1(h[k + 4 + (k + 1) * ldh]);
			}
		    }
		    mtemp1 = smlnum, mtemp2 = ulp * tst1;
		    if (Cabs1(h[k + 1 + k * ldh]) <= max(mtemp1, mtemp2)) {
			mtemp1 = Cabs1(h[k + 1 + k * ldh]);
			mtemp2 = Cabs1(h[k + (k + 1) * ldh]);
			h12 = max(mtemp1, mtemp2);
			mtemp1 = Cabs1(h[k + 1 + k * ldh]);
			mtemp2 = Cabs1(h[k + (k + 1) * ldh]);
			h21 = min(mtemp1, mtemp2);
			mtemp1 = Cabs1(h[k + 1 + (k + 1) * ldh]);
			mtemp2 = Cabs1(h[k + k * ldh] - h[k + 1 + (k + 1) * ldh]);
			h11 = max(mtemp1, mtemp2);
			mtemp1 = Cabs1(h[k + 1 + (k + 1) * ldh]);
			mtemp2 = Cabs1(h[k + k * ldh] - h[k + 1 + (k + 1) * ldh]);
			h22 = min(mtemp1, mtemp2);
			scl = h11 + h12;
			tst2 = h22 * (h11 / scl);
			mtemp1 = smlnum, mtemp2 = ulp * tst2;
			if (tst2 == Zero || h21 * (h12 / scl) <= max(mtemp1, mtemp2)) {
			    h[k + 1 + k * ldh] = Zero;
			}
		    }
		}
	    }
//==== Fill in the last row of each bulge. ====
	    mend = min(nbmps, (kbot - krcol - 1) / 3);
	    for (m = mtop; m <= mend; m++) {
		k = krcol + (m - 1) * 3;
		refsum = v[m * ldv + 1] * v[m * ldv + 3] * h[k + 4 + (k + 3) * ldh];
		h[k + 4 + (k + 1) * ldh] = -refsum;
		h[k + 4 + (k + 2) * ldh] = -refsum * conj(v[m * ldv + 2]);
		h[k + 4 + (k + 3) * ldh] = h[k + 4 + (k + 3) * ldh] - refsum * conj(v[m * ldv + 3]);
	    }
//==== End of near-the-diagonal bulge chase. ====
	}
//==== Use U (if accumulated) to update far-from-diagonal
//.    entries in H.  If required, use U to update Z as
//.    well. ====
	if (accum) {
	    if (wantt) {
		jtop = 1;
		jbot = n;
	    } else {
		jtop = ktop;
		jbot = kbot;
	    }
	    if (!blk22 || incol < ktop || ndcol > kbot || ns <= 2) {
//==== Updates not exploiting the 2-by-2 block
//.    structure of U.  K1 and NU keep track of
//.    the location and size of U in the special
//.    cases of introducing bulges and chasing
//.    bulges off the bottom.  In these special
//.    cases and in case the number of shifts
//.    is NS = 2, there is no 2-by-2 block
//.    structure to exploit.  ====
		k1 = max((INTEGER) 1, ktop - incol);
		nu = kdu - max((INTEGER) 0, ndcol - kbot) - k1 + 1;
//==== Horizontal Multiply ====
		for (jcol = min(ndcol, kbot) + 1; jcol <= jbot; jcol += nh) {
		    jlen = min(nh, jbot - jcol + 1);
		    Cgemm("C", "N", nu, jlen, nu, (COMPLEX) One, &u[k1 + k1 * ldu], ldu, &h[incol + k1 + jcol * ldh], ldh, (COMPLEX) Zero, &wh[0], ldwh);
		    Clacpy("ALL", nu, jlen, &wh[0], ldwh, &h[incol + k1 + jcol * ldh], ldh);
		}
//==== Vertical multiply ====
		for (jrow = jtop; jrow <= max(ktop, incol) - 1; jrow = jrow + nv) {
		    jlen = min(nv, max(ktop, incol) - jrow);
		    Cgemm("N", "N", jlen, nu, nu, (COMPLEX) One, &h[jrow + (incol + k1) * ldh], ldh, &u[k1 + k1 * ldu], ldu, (COMPLEX) Zero, &wv[0], ldwv);
		    Clacpy("ALL", jlen, nu, &wv[0], ldwv, &h[jrow + (incol + k1) * ldh], ldh);
		}
//==== Z multiply (also vertical) ====
		if (wantz) {
		    for (jrow = iloz; jrow <= ihiz; jrow += nv) {
			jlen = min(nv, ihiz - jrow + 1);
			Cgemm("N", "N", jlen, nu, nu, (COMPLEX) One, &z[jrow + (incol + k1) * ldz], ldz, &u[k1 + k1 * ldu], ldu, (COMPLEX) Zero, &wv[0], ldwv);
			Clacpy("ALL", jlen, nu, &wv[0], ldwv, &z[jrow + (incol + k1) * ldz], ldz);
		    }
		}
	    } else {
//==== Updates exploiting U's 2-by-2 block structure.
//.    (I2, I4, J2, J4 are the last rows and columns
//.    of the blocks.) ====
		i2 = (kdu + 1) / 2;
		i4 = kdu;
		j2 = i4 - i2;
		j4 = kdu;
//==== KZS and KNZ deal with the band of zeros
//.    along the diagonal of one of the triangular
//.    blocks. ====
		kzs = j4 - j2 - (ns + 1);
		knz = ns + 1;
//==== Horizontal multiply ====
		for (jcol = min(ndcol, kbot) + 1; jcol <= jbot; jcol += nh) {
		    jlen = min(nh, jbot - jcol + 1);
//==== Copy bottom of H to top+KZS of scratch ====
// (The first KZS rows get multiplied by zero.) ====
		    Clacpy("ALL", knz, jlen, &h[incol + 1 + j2 + jcol * ldh], ldh, &wh[kzs + 1 + ldwh], ldwh);
//==== Multiply by U21' ====
		    Claset("ALL", kzs, jlen, (COMPLEX) Zero, (COMPLEX) Zero, &wh[0], ldwh);
		    Ctrmm("L", "U", "C", "N", knz, jlen, (COMPLEX) One, &u[j2 + 1 + (kzs + 1) * ldu], ldu, &wh[kzs + 1 + ldwh], ldwh);
//==== Multiply top of H by U11' ====
		    Cgemm("C", "N", i2, jlen, j2, (COMPLEX) One, &u[0], ldu, &h[incol + 1 + jcol * ldh], ldh, (COMPLEX) One, &wh[0], ldwh);
//==== Copy top of H bottom of WH ====
		    Clacpy("ALL", j2, jlen, &h[incol + 1 + jcol * ldh], ldh, &wh[i2 + 1 + ldwh], ldwh);
//==== Multiply by U21' ====
		    Ctrmm("L", "L", "C", "N", j2, jlen, (COMPLEX) One, &u[(i2 + 1) * ldu + 1], ldu, &wh[i2 + 1 + ldwh], ldwh);
//==== Multiply by U22 ====
		    Cgemm("C", "N", i4 - i2, jlen, j4 - j2, (COMPLEX) One, &u[j2 + 1 + (i2 + 1) * ldu], ldu,
			  &h[incol + 1 + j2 + jcol * ldh], ldh, (COMPLEX) One, &wh[i2 + 1 + ldwh], ldwh);
//==== Copy it back ====
		    Clacpy("ALL", kdu, jlen, &wh[0], ldwh, &h[incol + 1 + jcol * ldh], ldh);
		}
//==== Vertical multiply ====
		for (jrow = jtop; jrow <= max(incol, ktop) - 1; jrow += nv) {
		    jlen = min(nv, max(incol, ktop) - jrow);
//==== Copy right of H to scratch (the first KZS
//.    columns get multiplied by zero) ====
		    Clacpy("ALL", jlen, knz, &h[jrow + (incol + 1 + j2) * ldh], ldh, &wv[(kzs + 1) * ldwv + 1], ldwv);
//==== Multiply by U21 ====
		    Claset("ALL", jlen, kzs, Zero, Zero, &wv[0], ldwv);
		    Ctrmm("R", "U", "N", "N", jlen, knz, One, &u[j2 + 1 + (kzs + 1) * ldu], ldu, &wv[(kzs + 1) * ldwv + 1], ldwv);
//==== Multiply by U11 ====
		    Cgemm("N", "N", jlen, i2, j2, One, &h[jrow + (incol + 1) * ldh], ldh, &u[0], ldu, (COMPLEX) One, &wv[0], ldwv);
//==== Copy left of H to right of scratch ====
		    Clacpy("ALL", jlen, j2, &h[jrow + (incol + 1) * ldh], ldh, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Multiply by U21 ====
		    Ctrmm("R", "L", "N", "N", jlen, i4 - i2, One, &u[(i2 + 1) * ldu + 1], ldu, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Multiply by U22 ====
		    Cgemm("N", "N", jlen, i4 - i2, j4 - j2, One, &h[jrow + (incol + 1 + j2) * ldh], ldh, &u[j2 + 1 + (i2 + 1) * ldu], ldu, One, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Copy it back ====
		    Clacpy("ALL", jlen, kdu, &wv[0], ldwv, &h[jrow + (incol + 1) * ldh], ldh);
		}
//==== Multiply Z (also vertical) ====
		if (wantz) {
		    for (jrow = iloz; jrow <= ihiz; jrow += nv) {
			jlen = min(nv, ihiz - jrow + 1);
//==== Copy right of Z to left of scratch (first
//.     KZS columns get multiplied by zero) ====
			Clacpy("ALL", jlen, knz, &z[jrow + (incol + 1 + j2) * ldz], ldz, &wv[(kzs + 1) * ldwv + 1], ldwv);
//==== Multiply by U12 ====
			Claset("ALL", jlen, kzs, Zero, Zero, &wv[0], ldwv);
			Ctrmm("R", "U", "N", "N", jlen, knz, One, &u[j2 + 1 + (kzs + 1) * ldu], ldu, &wv[(kzs + 1) * ldwv + 1], ldwv);
//==== Multiply by U11 ====
			Cgemm("N", "N", jlen, i2, j2, One, &z[jrow + (incol + 1) * ldz], ldz, &u[0], ldu, (COMPLEX) One, &wv[0], ldwv);
//==== Copy left of Z to right of scratch ====
			Clacpy("ALL", jlen, j2, &z[jrow + (incol + 1) * ldz], ldz, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Multiply by U21 ====
			Ctrmm("R", "L", "N", "N", jlen, i4 - i2, One, &u[(i2 + 1) * ldu + 1], ldu, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Multiply by U22 ====
			Cgemm("N", "N", jlen, i4 - i2, j4 - j2, One, &z[jrow + (incol + 1 + j2) * ldz], ldz, &u[j2 + 1 + (i2 + 1) * ldu], ldu, One, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Copy the result back to Z ====
			Clacpy("ALL", jlen, kdu, &wv[0], ldwv, &z[jrow + (incol + 1) * ldz], ldz);
		    }
		}
	    }
	}
    }
//==== End of ZLAQR5 ====
    return;
}
