/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Claqr0.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Claqr0(LOGICAL wantt, LOGICAL wantz, INTEGER n, INTEGER ilo, INTEGER ihi, COMPLEX * h, INTEGER ldh, COMPLEX * w, INTEGER * iloz,
       INTEGER * ihiz, COMPLEX * z, INTEGER ldz, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, k;
    REAL s;
    COMPLEX aa, bb, cc, dd;
    INTEGER ld, nh, it, ks, kt, ku, kv, ls, ns, nw;
    COMPLEX tr2, det;
    INTEGER inf, kdu, nho, nve, kwh, nsr, nwr, kwv, ndfl, kbot, nmin;
    COMPLEX swap;
    INTEGER ktop;
    COMPLEX zdum[1];
    INTEGER kacc22;
    INTEGER nwinc = 0;
    INTEGER itmax, nsmax, nwmax, kwtop;
    INTEGER nibble;
    char jbcmpz[2];
    COMPLEX rtdisc;
    INTEGER sorted;
    INTEGER lwkopt;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;

    *info = 0;
//==== Quick return for N = 0: nothing to do. ====
    if (n == 0) {
	work[1] = One;
	return;
    }
//==== Set up job flags for iMlaenv. ====
    if (wantt) {
	jbcmpz[0] = 'S';
    } else {
	jbcmpz[0] = 'E';
    }
    if (wantz) {
	jbcmpz[1] = 'V';
    } else {
	jbcmpz[1] = 'N';
    }
//==== Tiny matrices must use ZLAHQR. ====
    if (n <= 11) {
//==== Estimate optimal workspace. ====
	lwkopt = 1;
	if (lwork != -1) {
	    Clahqr(wantt, wantz, n, ilo, ihi, &h[0], ldh, &w[1], (*iloz), (*ihiz), &z[0], ldz, info);
	}
    } else {
//==== Use small bulge multi-shift QR with aggressive early
//.    deflation on larger-than-tiny matrices. ====
//==== Hope for the best. ====
	*info = 0;
//==== NWR = recommended deflation window size.  At this
//.    point,  N .GT. NTINY = 11, so there is enough
//.    subdiagonal workspace for NWR.GE.2 as required.
//.    (In fact, there is enough subdiagonal space for
//.    NWR.GE.3.) ====
	nwr = iMlaenv(13, "Claqr0", jbcmpz, n, ilo, ihi, lwork);
	nwr = max((INTEGER) 2, nwr);
	nwr = min(min(ihi - ilo + 1, (n - 1) / 3), nwr);
	nw = nwr;
//==== NSR = recommended number of simultaneous shifts.
//.    At this point N .GT. NTINY = 11, so there is at
//.    enough subdiagonal workspace for NSR to be even
//.    and greater than or equal to two as required. ====
	nsr = iMlaenv(15, "Claqr0", jbcmpz, n, ilo, ihi, lwork);
	nsr = min(min(nsr, (n + 6) / 9), ihi - ilo);
	nsr = max((INTEGER) 2, nsr - nsr % 2);
//==== Estimate optimal workspace ====
//==== Workspace query call to ZLAQR3 ====
	Claqr3(wantt, wantz, n, ilo, ihi, nwr + 1, &h[0], ldh, *iloz, *ihiz, &z[0], ldz, &ls, &ld, &w[1], &h[0], ldh, n, &h[0], ldh, n, &h[0], ldh, &work[0], -1);
//==== Optimal workspace = MAX(ZLAQR5, ZLAQR3) ====
	lwkopt = max(nsr * 3 / 2, (INTEGER) cast2double(work[1].real()));
//==== Quick return in case of workspace query. ====
	if (lwork == -1) {
	    work[1] = lwkopt;
	    return;
	}
//==== ZLAHQR/ZLAQR0 crossover point ====
	nmin = iMlaenv(12, "Claqr0", jbcmpz, n, ilo, ihi, lwork);
	nmin = max((INTEGER) 11, nmin);
//==== Nibble crossover point ====
	nibble = iMlaenv(14, "Claqr0", jbcmpz, n, ilo, ihi, lwork);
	nibble = max((INTEGER) 0, nibble);
//==== Accumulate reflections during ttswp?  Use block
//.    2-by-2 structure during matrix-matrix multiply? ====
	kacc22 = iMlaenv(16, "Claqr0", jbcmpz, n, ilo, ihi, lwork);
	kacc22 = max((INTEGER) 0, kacc22);
	kacc22 = min((INTEGER) 2, kacc22);
//==== NWMAX = the largest possible deflation window for
//.    which there is sufficient workspace. ====
	nwmax = min((n - 1) / 3, lwork / 2);
//==== NSMAX = the Largest number of simultaneous shifts
//.    for which there is sufficient workspace. ====
	nsmax = min((n + 6) / 9, (lwork << 1) / 3);
	nsmax = nsmax - nsmax % 2;
//==== NDFL: an iteration count restarted at deflation. ====
	ndfl = 0;
//==== ITMAX = iteration limit ====
	itmax = max((INTEGER) 10, ihi - ilo + 1) * 30;
//==== Last row and column in the active block ====
	kbot = ihi;
//==== Main Loop ====
	for (it = 1; it <= itmax; it++) {
//==== Done when KBOT falls below ILO ====
	    if (kbot < ilo) {
		goto L80;
	    }
//==== Locate active block ====
	    for (k = kbot; k >= ilo + 1; k--) {
		if (h[k + (k - 1) * ldh] == Zero) {
		    goto L20;
		}
	    }
	    k = ilo;
	  L20:
	    ktop = k;
//==== Select deflation window size ====
	    nh = kbot - ktop + 1;
	    if (ndfl < 5 || nh < nw) {
/*              ==== Typical deflation window.  If possible and */
/*              .    advisable, nibble the entire active block. */
/*              .    If not, use size NWR or NWR+1 depending upon */
/*              .    which has the smaller corresponding subdiagonal */
/*              .    entry (a heuristic). ==== */
		nwinc = MTRUE;
		if (nh <= min(nmin, nwmax)) {
		    nw = nh;
		} else {
		    nw = min(min(nwr, nh), nwmax);
		    if (nw < nwmax) {
			if (nw >= nh - 1) {
			    nw = nh;
			} else {
			    kwtop = kbot - nw + 1;
			    if (Cabs1(h[kwtop + (kwtop - 1) * ldh]) > Cabs1(h[kwtop - 1 + (kwtop - 2) * ldh])) {
				++nw;
			    }
			}
		    }
		}
	    } else {
/*              ==== Exceptional deflation window.  If there have */
/*              .    been no deflations in KEXNW or more iterations, */
/*              .    then vary the deflation window size.   At first, */
/*              .    because, larger windows are, in general, more */
/*              .    powerful than smaller ones, rapidly increase the */
/*              .    window up to the maximum reasonable and possible. */
/*              .    Then maybe try a slightly smaller window.  ==== */
		if (nwinc && nw < min(nwmax, nh)) {
		    nw = min(min(nwmax, nh), nw << 1);
		} else {
		    nwinc = MFALSE;
		    if (nw == nh && nh > 2) {
			nw = nh - 1;
		    }
		}
	    }
/*           ==== Aggressive early deflation: */
/*           .    split workspace under the subdiagonal into */
/*           .      - an nw-by-nw work array V in the lower */
/*           .        left-hand-corner, */
/*           .      - an NW-by-at-least-NW-but-more-is-better */
/*           .        (NW-by-NHO) horizontal work array along */
/*           .        the bottom edge, */
/*           .      - an at-least-NW-but-more-is-better (NHV-by-NW) */
/*           .        vertical work array along the left-hand-edge. */
/*           .        ==== */
	    kv = n - nw + 1;
	    kt = nw + 1;
	    nho = n - nw - 1 - kt + 1;
	    kwv = nw + 2;
	    nve = n - nw - kwv + 1;
//==== Aggressive early deflation ====
	    Claqr3(wantt, wantz, n, ktop, kbot, nw, &h[0], ldh, *iloz, *ihiz, &z[0], ldz, &ls, &ld, &w[1], &h[kv + ldh], ldh, nho,
		   &h[kv + kt * ldh], ldh, nve, &h[kwv + ldh], ldh, &work[0], lwork);
/*           ==== Adjust KBOT accounting for new deflations. ==== */
	    kbot = kbot - ld;
//==== KS points to the shifts. ====
	    ks = kbot - ls + 1;
//==== Skip an expensive QR sweep if there is a (partly
//.    heuristic) reason to expect that many eigenvalues
//.    will deflate without it.  Here, the QR sweep is
//.    skipped if many eigenvalues have just been deflated
//.    or if the remaining active block is small.
	    if (ld == 0 || (ld * 100 <= nw * nibble && kbot - ktop + 1) > min(nmin, nwmax)) {
//==== NS = nominal number of simultaneous shifts.
//.    This may be lowered (slightly) if ZLAQR3
//.    did not provide that many shifts. ====
		ns = min(min(nsmax, nsr), max((INTEGER) 2, kbot - ktop));
		ns = ns - ns % 2;
//==== If there have been no deflations
//.    in a multiple of KEXSH iterations,
//.    then try exceptional shifts.
//.    Otherwise use shifts provided by
//.    ZLAQR3 above or from the eigenvalues
//.    of a trailing principal submatrix. ====
		if (ndfl % 6 == 0) {
		    ks = kbot - ns + 1;
		    for (i = kbot; i >= ks + 1; i = i - 2) {
			w[i] = h[i + i * ldh] + Cabs1(h[i + (i - 1) * ldh]) * .75;
			w[i - 1] = w[i];
		    }
		} else {

//==== Got NS/2 or fewer shifts? Use ZLAQR4 or
//.    ZLAHQR on a trailing principal submatrix to
//.    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
//.    there is enough space below the subdiagonal
//.    to fit an NS-by-NS scratch array.) ====
		    if (kbot - ks + 1 <= ns / 2) {
			ks = kbot - ns + 1;
			kt = n - ns + 1;
			Clacpy("A", ns, ns, &h[ks + ks * ldh], ldh, &h[kt + ldh], ldh);
			if (ns > nmin) {
			    Claqr4(MFALSE, MFALSE, ns, 1, ns, &h[kt + ldh], ldh, &w[ks], 1, 1, zdum, 1, &work[0], lwork, &inf);
			} else {
			    Clahqr(MFALSE, MFALSE, ns, 1, ns, &h[kt + ldh], ldh, &w[ks], 1, 1, zdum, 1, &inf);
			}
			ks = ks + inf;
//==== In case of a rare QR failure use
//.    eigenvalues of the trailing 2-by-2
//.    principal submatrix.  Scale to avoid
//.    overflows, underflows and subnormals.
//.    (The scale factor S can not be zero,
//.    because H(KBOT,KBOT-1) is nonzero.) ====
			if (ks >= kbot) {
			    s = Cabs1(h[kbot - 1 + (kbot - 1) * ldh]) + Cabs1(h[kbot + (kbot - 1) * ldh]) + Cabs1(h[kbot - 1 + kbot * ldh]) + Cabs1(h[kbot + kbot * ldh]);
			    aa = h[kbot - 1 + (kbot - 1) * ldh] / s;
			    cc = h[kbot + (kbot - 1) * ldh] / s;
			    bb = h[kbot - 1 + kbot * ldh] / s;
			    dd = h[kbot + kbot * ldh] / s;
			    tr2 = (aa + dd) / Two;
			    det = (aa - tr2) * (dd - tr2) - bb * cc;
			    rtdisc = sqrt(-det);
			    w[kbot - 1] = (tr2 + rtdisc) * s;
			    w[kbot] = (tr2 - rtdisc) * s;
			    ks = kbot - 1;
			}
		    }
		    if (kbot - ks + 1 > ns) {
//==== Sort the shifts (Helps a little) ====
			sorted = MFALSE;
			for (k = kbot; k >= ks + 1; k--) {
			    if (sorted) {
				goto L60;
			    }
			    sorted = MTRUE;
			    for (i = ks; i <= k - 1; i++) {
				if (Cabs1(w[i]) < Cabs1(w[i + 1])) {
				    sorted = MFALSE;
				    swap = w[i];
				    w[i] = w[i + 1];
				    w[i + 1] = swap;
				}
			    }
			}
		      L60:
			;
		    }
		}
//==== If there are only two shifts, then use
//.    only one.  ====
		if (kbot - ks + 1 == 2) {
		    if (Cabs1(w[kbot] - h[kbot + kbot * ldh]) < Cabs1(w[kbot - 1] - h[kbot + kbot * ldh])) {
			w[kbot - 1] = w[kbot];
		    } else {
			w[kbot] = w[kbot - 1];
		    }
		}
//==== Use up to NS of the the smallest magnatiude
//.    shifts.  If there aren't NS shifts available,
//.    then use them all, possibly dropping one to
//.    make the number of shifts even. ====
		ns = min(ns, kbot - ks + 1);
		ns = ns - ns % 2;
		ks = kbot - ns + 1;
//==== Small-bulge multi-shift QR sweep:
//.    split workspace under the subdiagonal into
//.    - a KDU-by-KDU work array U in the lower
//.      left-hand-corner,
//.    - a KDU-by-at-least-KDU-but-more-is-better
//.      (KDU-by-NHo) horizontal work array WH along
//.      the bottom edge,
//.    - and an at-least-KDU-but-more-is-better-by-KDU
//.      (NVE-by-KDU) vertical work WV arrow along
//.      the left-hand-edge. ====
		kdu = ns * 3 - 3;
		ku = n - kdu + 1;
		kwh = kdu + 1;
		nho = n - kdu - 3 - (kdu + 1) + 1;
		kwv = kdu + 4;
		nve = n - kdu - kwv + 1;
//==== Small-bulge multi-shift QR sweep ====
		Claqr5(wantt, wantz, kacc22, n, ktop, kbot, ns, &w[ks], &h[0], ldh, *iloz, *ihiz, &z[0], ldz, &work[0], 3,
		       &h[ku + ldh], ldh, nve, &h[kwv + ldh], ldh, nho, &h[ku + kwh * ldh], ldh);
	    }
//==== Note progress (or the lack of it). ====
	    if (ld > 0) {
		ndfl = 0;
	    } else {
		++ndfl;
	    }
//==== End of main loop ====
	}
//==== Iteration limit exceeded.  Set INFO to show where
//.    the problem occurred and exit. ====
	*info = kbot;
      L80:
	;
    }
//==== Return the optimal value of LWORK. ====
    work[1] = lwkopt;
//==== End of ZLAQR0 ====
    return;
}
