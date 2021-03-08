/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaqr5.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rlaqr5(LOGICAL wantt, LOGICAL wantz, INTEGER kac22, INTEGER n, INTEGER ktop, INTEGER kbot,
       INTEGER nshfts, REAL * sr, REAL * si, REAL * h, INTEGER ldh,
       INTEGER iloz, INTEGER ihiz, REAL * z, INTEGER ldz, REAL * v, INTEGER ldv, REAL * u, INTEGER ldu, INTEGER nv, REAL * wv, INTEGER ldwv, INTEGER nh, REAL * wh, INTEGER ldwh)
{
    INTEGER i, j, k, m, i2, j2, i4, j4, k1;
    REAL h11, h12, h21, h22;
    INTEGER m22, ns, nu;
    REAL vt[3], scl;
    INTEGER kdu, kms;
    REAL ulp;
    INTEGER knz, kzs;
    REAL tst1, tst2, beta;
    INTEGER blk22, bmp22;
    INTEGER mend, jcol, jlen, jbot, mbot;
    REAL swap;
    INTEGER jtop, jrow, mtop;
    REAL alpha;
    INTEGER accum;
    INTEGER ndcol, incol, krcol, nbmps;
    REAL safmin;
    REAL safmax, refsum;
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
//==== Shuffle shifts into pairs of real shifts and pairs
//.    of complex conjugate shifts assuming complex
//.    conjugate shifts are already adjacent to one
//.    another. ====
    for (i = 0; i < nshfts - 2; i = i + 2) {
	if (si[i] != -si[i + 1]) {
	    swap = sr[i];
	    sr[i] = sr[i + 1];
	    sr[i + 1] = sr[i + 2];
	    sr[i + 2] = swap;
	    swap = si[i];
	    si[i] = si[i + 1];
	    si[i + 1] = si[i + 2];
	    si[i + 2] = swap;
	}
    }
//==== NSHFTS is supposed to be even, but if is odd,
//.   then simply reduce it by one.  The shuffle above
//.    ensures that the dropped shift is real and that
//.    the remaining shifts are paired. ====
    ns = nshfts - nshfts % 2;
//==== Machine constants for deflation ====
    safmin = Rlamch("SAFE MINIMUM");
    safmax = One / safmin;
    ulp = Rlamch("PRECISION");
    smlnum = safmin * ((REAL) double (n) / ulp);
//==== Use accumulated reflections to update far-from-diagonal
//.    entries ? ====
    accum = kac22 == 1 || kac22 == 2;
//==== If so, exploit the 2-by-2 block structure? ====
    blk22 = ns > 2 && kac22 == 2;
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
	    Rlaset("ALL", kdu, kdu, Zero, One, u, ldu);
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
		    Rlaqr1(3, &h[ktop + ktop * ldh], ldh, sr[(m * 2) - 1], si[(m * 2) - 1], sr[m * 2], si[m * 2], &v[m * ldv + 1]);
		    alpha = v[m * ldv + 1];
		    Rlarfg(3, &alpha, &v[m * ldv + 2], 1, &v[m * ldv + 1]);
		} else {
		    beta = h[k + 1 + k * ldh];
		    v[m * ldv + 2] = h[k + 2 + k * ldh];
		    v[m * ldv + 3] = h[k + 3 + k * ldh];
		    Rlarfg(3, &beta, &v[m * ldv + 2], 1, &v[m * ldv + 1]);
//==== A Bulge may collapse because of vigilant
//.    deflation or destructive underflow.  (The
//.    initial bulge is always collapsed.) Use
//.    the two-small-subdiagonals trick to try
//.    to get it started again. If V(2,M).NE.0 and
//.    V(3,M) = H(K+3,K+1) = H(K+3,K+2) = 0, then
//.    this bulge is collapsing into a zero
//.    subdiagonal.  It will be restarted next
//.    trip through the loop.)
		    if (v[m * ldv + 1] != Zero && (v[m * ldv + 3] != Zero || (h[k + 3 + (k + 1) * ldh] == Zero && h[k + 3 + (k + 2) * ldh] == Zero))) {
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
			Rlaqr1(3, &h[k + 1 + (k + 1) * ldh], ldh, sr[(m * 2) - 1], si[(m * 2) - 1], sr[m * 2], si[m * 2], vt);
			scl = abs(vt[0]) + abs(vt[1]) + abs(vt[2]);
			if (scl != Zero) {
			    vt[0] = vt[0] / scl;
			    vt[1] = vt[1] / scl;
			    vt[2] = vt[2] / scl;
			}
//==== The following is the traditional and
//.    conservative two-small-subdiagonals
//.    test.  ====
//.
			if (abs(h[k + 1 + k * ldh]) * (abs(vt[1]) +
						       abs(vt[2])) > ulp * abs(vt[0]) * (abs(h[k + k * ldh]) + abs(h[k + 1 + (k + 1) * ldh]) + abs(h[k + 2 + (k + 2) * ldh]))) {
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
			    Rlarfg(3, &alpha, &vt[1], 1, vt);
			    refsum = h[k + 1 + k * ldh] + h[k + 2 + k * ldh] * vt[1] + h[k + 3 + k * ldh]
				* vt[2];
			    h[k + 1 + k * ldh] = h[k + 1 + k * ldh] - vt[0] * refsum;
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
		    Rlaqr1(2, &h[k + 1 + (k + 1) * ldh], ldh, sr[(m22 * 2) - 1], si[(m22 * 2) - 1], sr[m22 * 2], si[m22 * 2], &v[m22 * ldv + 1]);
		    beta = v[m22 * ldv + 1];
		    Rlarfg(2, &beta, &v[m22 * ldv + 2], 1, &v[m22 * ldv + 1]);
		} else {
		    beta = h[k + 1 + k * ldh];
		    v[m22 * ldv + 2] = h[k + 2 + k * ldh];
		    Rlarfg(2, &beta, &v[m22 * ldv + 2], 1, &v[m22 * ldv + 1]);
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
		    refsum = v[m * ldv + 1] * (h[k + 1 + j * ldh] + v[m * ldv + 2] * h[k + 2 + j * ldh] + v[m * ldv + 3] * h[k + 3 + j * ldh]);
		    h[k + 1 + j * ldh] = h[k + 1 + j * ldh] - refsum;
		    h[k + 2 + j * ldh] = h[k + 2 + j * ldh] - refsum * v[m * ldv + 2];
		    h[k + 3 + j * ldh] = h[k + 3 + j * ldh] - refsum * v[m * ldv + 3];

		}

	    }
	    if (bmp22) {
		k = krcol + (m22 - 1) * 3;
		for (j = max(k + 1, ktop); j <= jbot; j++) {
		    refsum = v[m22 * ldv + 1] * (h[k + 1 + j * ldh] + v[m22 * ldv + 2] * h[k + 2 + j * ldh]);
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
			h[j + (k + 2) * ldh] = h[j + (k + 2) * ldh] - refsum * v[m * ldv + 2];
			h[j + (k + 3) * ldh] = h[j + (k + 3) * ldh] - refsum * v[m * ldv + 3];
		    }
		    if (accum) {
//==== Accumulate U. (If necessary, update Z later
//.    with with an efficient matrix-matrix
//.    multiply.) ====
			kms = k - incol;
			for (j = max((INTEGER) 1, ktop - incol); j <= kdu; j++) {
			    refsum = v[m * ldv + 1] * (u[j + (kms + 1) * ldu] + v[m * ldv + 2] * u[j + (kms + 2) * ldu] + v[m * ldv + 3] * u[j + (kms + 3) * ldu]);
			    u[j + (kms + 1) * ldu] = u[j + (kms + 1) * ldu] - refsum;
			    u[j + (kms + 2) * ldu] = u[j + (kms + 2) * ldu] - refsum * v[m * ldv + 2];
			    u[j + (kms + 3) * ldu] = u[j + (kms + 3) * ldu] - refsum * v[m * ldv + 3];
			}
		    } else if (wantz) {
//==== U is not accumulated, so update Z
//.    now by multiplying by reflections
//.    from the right. ====
			for (j = iloz; j <= ihiz; j++) {
			    refsum = v[m * ldv + 1] * (z[j + (k + 1) * ldz] + v[m * ldv + 2] * z[j + (k + 2) * ldz] + v[m * ldv + 3] * z[j + (k + 3) * ldz]);
			    z[j + (k + 1) * ldz] = z[j + (k + 1) * ldz] - refsum;
			    z[j + (k + 2) * ldz] = z[j + (k + 2) * ldz] - refsum * v[m * ldv + 2];
			    z[j + (k + 3) * ldz] = z[j + (k + 3) * ldz] - refsum * v[m * ldv + 3];
			}
		    }
		}
	    }
//== Special case: 2-by-2 reflection (if needed) ====
	    k = krcol + (m22 - 1) * 3;
	    if (bmp22 && v[m22 * ldv + 1] != Zero) {
		for (j = jtop; j <= min(kbot, k + 3); j++) {
		    refsum = v[m22 * ldv + 1] * (h[j + (k + 1) * ldh]
						 + v[m22 * ldv + 2] * h[j + (k + 2) * ldh]);
		    h[j + (k + 1) * ldh] = h[j + (k + 1) * ldh] - refsum;
		    h[j + (k + 2) * ldh] = h[j + (k + 2) * ldh] - refsum * v[m22 * ldv + 2];
		}
		if (accum) {
		    kms = k - incol;
		    for (j = max((INTEGER) 1, ktop - incol); j <= kdu; j++) {
			refsum = v[m22 * ldv + 1] * (u[j + (kms + 1) * ldu] + v[m22 * ldv + 2] * u[j + (kms + 2) * ldu]);
			u[j + (kms + 1) * ldu] = u[j + (kms + 1) * ldu] - refsum;
			u[j + (kms + 2) * ldu] = u[j + (kms + 2) * ldu] - refsum * v[m22 * ldv + 2];

		    }
		} else if (wantz) {
		    for (j = iloz; j <= ihiz; j++) {
			refsum = v[m22 * ldv + 1] * (z[j + (k + 1) * ldz] + v[m22 * ldv + 2] * z[j + (k + 2) * ldz]);
			z[j + (k + 1) * ldz] = z[j + (k + 1) * ldz] - refsum;
			z[j + (k + 2) * ldz] = z[j + (k + 2) * ldz] - refsum * v[m22 * ldv + 2];
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
		    tst1 = abs(h[k + k * ldh]) + abs(h[k + 1 + (k + 1) * ldh]);
		    if (tst1 == Zero) {
			if (k >= ktop + 1) {
			    tst1 += abs(h[k + (k - 1) * ldh]);
			}
			if (k >= ktop + 2) {
			    tst1 += abs(h[k + (k - 2) * ldh]);
			}
			if (k >= ktop + 3) {
			    tst1 += abs(h[k + (k - 3) * ldh]);
			}
			if (k <= kbot - 2) {
			    tst1 += abs(h[k + 2 + (k + 1) * ldh]);
			}
			if (k <= kbot - 3) {
			    tst1 += abs(h[k + 3 + (k + 1) * ldh]);
			}
			if (k <= kbot - 4) {
			    tst1 += abs(h[k + 4 + (k + 1) * ldh]);
			}
		    }
		    mtemp1 = smlnum, mtemp2 = ulp * tst1;
		    if (abs(h[k + 1 + k * ldh]) <= max(mtemp1, mtemp2)) {
			mtemp1 = abs(h[k + 1 + k * ldh]), mtemp2 = abs(h[k + (k + 1) * ldh]);
			h12 = max(mtemp1, mtemp2);
			mtemp1 = abs(h[k + 1 + k * ldh]), mtemp2 = abs(h[k + (k + 1) * ldh]);
			h21 = min(mtemp1, mtemp2);
			mtemp1 = abs(h[k + 1 + (k + 1) * ldh]), mtemp2 = abs(h[k + k * ldh] - h[k + 1 + (k + 1) * ldh]);
			h11 = max(mtemp1, mtemp2);
			mtemp1 = abs(h[k + 1 + (k + 1) * ldh]), mtemp2 = abs(h[k + k * ldh] - h[k + 1 + (k + 1) * ldh]);
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
		h[k + 4 + (k + 2) * ldh] = -refsum * v[m * ldv + 2];
		h[k + 4 + (k + 3) * ldh] = h[k + 4 + (k + 3) * ldh] - refsum * v[m * ldv + 3];
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
		for (jcol = min(ndcol, kbot) + 1; jcol <= jbot; jcol = jcol + nh) {
		    jlen = min(nh, jbot - jcol + 1);
		    Rgemm("C", "N", nu, jlen, nu, One, &u[k1 + k1 * ldu], ldu, &h[incol + k1 + jcol * ldh], ldh, Zero, wh, ldwh);
		    Rlacpy("ALL", nu, jlen, wh, ldwh, &h[incol + k1 + jcol * ldh], ldh);
		}
//==== Vertical multiply ====
		for (jrow = jtop; jrow <= max(ktop, incol) - 1; jrow = jrow + nv) {
		    jlen = min(nv, max(ktop, incol) - jrow);
		    Rgemm("N", "N", jlen, nu, nu, One, &h[jrow + (incol + k1) * ldh], ldh, &u[k1 + k1 * ldu], ldu, Zero, wv, ldwv);
		    Rlacpy("ALL", jlen, nu, wv, ldwv, &h[jrow + (incol + k1) * ldh], ldh);
		}
//==== Z multiply (also vertical) ====
		if (wantz) {
		    for (jrow = iloz; jrow <= ihiz; jrow = jrow + nv) {
			jlen = min(nv, ihiz - jrow + 1);
			Rgemm("N", "N", jlen, nu, nu, One, &z[jrow + (incol + k1) * ldz], ldz, &u[k1 + k1 * ldu], ldu, Zero, wv, ldwv);
			Rlacpy("ALL", jlen, nu, wv, ldwv, &z[jrow + (incol + k1) * ldz], ldz);
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
		for (jcol = min(ndcol, kbot) + 1; jcol <= jbot; jcol = jcol + nh) {
		    jlen = min(nh, jbot - jcol + 1);
//==== Copy bottom of H to top+KZS of scratch ====
// (The first KZS rows get multiplied by zero.) ====
		    Rlacpy("ALL", knz, jlen, &h[incol + 1 + j2 + jcol * ldh], ldh, &wh[kzs + 1 + ldwh], ldwh);
//==== Multiply by U21' ====
		    Rlaset("ALL", kzs, jlen, Zero, Zero, wh, ldwh);
		    Rtrmm("L", "U", "C", "N", knz, jlen, One, &u[j2 + 1 + (kzs + 1) * ldu], ldu, &wh[kzs + 1 + ldwh], ldwh);
//==== Multiply top of H by U11' ====
		    Rgemm("C", "N", i2, jlen, j2, One, u, ldu, &h[incol + 1 + jcol * ldh], ldh, One, wh, ldwh);
//==== Copy top of H bottom of WH ====
		    Rlacpy("ALL", j2, jlen, &h[incol + 1 + jcol * ldh]
			   , ldh, &wh[i2 + 1 + ldwh], ldwh);
//==== Multiply by U21' ====
		    Rtrmm("L", "L", "C", "N", j2, jlen, One, &u[(i2 + 1)
								* ldu + 1], ldu, &wh[i2 + 1 + ldwh], ldwh);
//==== Multiply by U22 ====
		    Rgemm("C", "N", i4 - i2, jlen, j4 - j2, One, &u[j2 + 1 + (i2 + 1) * ldu], ldu, &h[incol + 1 + j2 + jcol * ldh], ldh, One, &wh[i2 + 1 + ldwh], ldwh);
//==== Copy it back ====
		    Rlacpy("ALL", kdu, jlen, wh, ldwh, &h[incol + 1 + jcol * ldh], ldh);
		}
//==== Vertical multiply ====
		for (jrow = jtop; jrow <= max(incol, ktop) - 1; jrow = jrow + nv) {
		    jlen = min(nv, max(incol, ktop) - jrow);
//==== Copy right of H to scratch (the first KZS
//.    columns get multiplied by zero) ====
		    Rlacpy("ALL", jlen, knz, &h[jrow + (incol + 1 + j2) * ldh], ldh, &wv[(kzs + 1) * ldwv + 1], ldwv);
//==== Multiply by U21 ====
		    Rlaset("ALL", jlen, kzs, Zero, Zero, wv, ldwv);
		    Rtrmm("R", "U", "N", "N", jlen, knz, One, &u[j2 + 1 + (kzs + 1) * ldu], ldu, &wv[(kzs + 1) * ldwv + 1], ldwv);
//==== Multiply by U11 ====
		    Rgemm("N", "N", jlen, i2, j2, One, &h[jrow + (incol + 1) * ldh], ldh, u, ldu, One, wv, ldwv);
//=== Copy left of H to right of scratch ====
		    Rlacpy("ALL", jlen, j2, &h[jrow + (incol + 1) * ldh], ldh, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Multiply by U21 ====
		    Rtrmm("R", "L", "N", "N", jlen, i4 - i2, One, &u[(i2 + 1) * ldu + 1], ldu, &wv[(i2 + 1) * ldwv + 1]
			  , ldwv);
//==== Multiply by U22 ====
		    Rgemm("N", "N", jlen, i4 - i2, j4 - j2, One, &h[jrow + (incol + 1 + j2) * ldh], ldh, &u[j2 + 1 + (i2 + 1) * ldu], ldu, One, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Copy it back ====
		    Rlacpy("ALL", jlen, kdu, wv, ldwv, &h[jrow + (incol + 1) * ldh], ldh);
		}
//==== Multiply Z (also vertical) ====
		if (wantz) {
		    for (jrow = iloz; jrow <= ihiz; jrow = jrow + nv) {
			jlen = min(nv, ihiz - jrow + 1);
//==== Copy right of Z to left of scratch (first
//.     KZS columns get multiplied by zero) ====
			Rlacpy("ALL", jlen, knz, &z[jrow + (incol + 1 + j2) * ldz], ldz, &wv[(kzs + 1) * ldwv + 1], ldwv);
//==== Multiply by U12 ====
			Rlaset("ALL", jlen, kzs, Zero, Zero, wv, ldwv);
			Rtrmm("R", "U", "N", "N", jlen, knz, One, &u[j2 + 1 + (kzs + 1) * ldu], ldu, &wv[(kzs + 1) * ldwv + 1], ldwv);
//==== Multiply by U11 ====
			Rgemm("N", "N", jlen, i2, j2, One, &z[jrow + (incol + 1) * ldz], ldz, u, ldu, One, wv, ldwv);
//==== Copy left of Z to right of scratch ====
			Rlacpy("ALL", jlen, j2, &z[jrow + (incol + 1) * ldz], ldz, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Multiply by U21 ====
			Rtrmm("R", "L", "N", "N", jlen, i4 - i2, One, &u[(i2 + 1) * ldu + 1], ldu, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Multiply by U22 ====
			Rgemm("N", "N", jlen, i4 - i2, j4 - j2, One, &z[jrow + (incol + 1 + j2) * ldz], ldz, &u[j2 + 1 + (i2 + 1) * ldu], ldu, One, &wv[(i2 + 1) * ldwv + 1], ldwv);
//==== Copy the result back to Z ====
			Rlacpy("ALL", jlen, kdu, wv, ldwv, &z[jrow + (incol + 1) * ldz], ldz);
		    }
		}
	    }
	}
    }
//==== End of DLAQR5 ====
    return;
}
