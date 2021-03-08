/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaqr3.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rlaqr3(LOGICAL wantt, LOGICAL wantz, INTEGER n, INTEGER ktop, INTEGER kbot, INTEGER nw,
       REAL * h, INTEGER ldh, INTEGER iloz, INTEGER ihiz, REAL * z, INTEGER ldz, INTEGER ns,
       INTEGER nd, REAL * sr, REAL * si, REAL * v, INTEGER ldv, INTEGER nh, REAL * t, INTEGER ldt, INTEGER nv, REAL * wv, INTEGER ldwv, REAL * work, INTEGER lwork)
{
    INTEGER i, j, k;
    REAL s, aa, bb, cc, dd, cs, sn;
    INTEGER jw;
    REAL evi, evk, foo;
    INTEGER kln;
    REAL tau, ulp;
    INTEGER lwk1, lwk2, lwk3;
    REAL beta;
    INTEGER kend, kcol, info, nmin, ifst, ilst, ltop, krow;
    INTEGER bulge;
    INTEGER infqr, kwtop;
    REAL safmin;
    REAL safmax;
    INTEGER sorted;
    REAL smlnum;
    INTEGER lwkopt;
    REAL mtemp1, mtemp2;
    REAL mtemp3, mtemp4;
    REAL Zero = 0.0, One = 1.0;

//==== Estimate optimal workspace. ====
    jw = min(nw, kbot - ktop + 1);
    if (jw <= 2) {
	lwkopt = 1;
    } else {
//==== Workspace query call to DGEHRD ====
	Rgehrd(jw, 1, jw - 1, t, ldt, &work[0], &work[0], -1, &info);
	lwk1 = (INTEGER) cast2double(work[1]);
//==== Workspace query call to DORGHR ====
	Rorghr(jw - 1, 1, jw - 1, t, ldt, &work[0], &work[0], -1, &info);
	lwk2 = (INTEGER) cast2double(work[1]);
//==== Workspace query call to DLAQR4 ====
	Rlaqr4(MTRUE, MTRUE, jw, 1, jw, t, ldt, &sr[1], &si[1], 1, jw, &v[0], ldv, &work[0], -1, &infqr);
	lwk3 = (INTEGER) cast2double(work[1]);
//==== Optimal workspace ====
	lwkopt = max(jw + max(lwk1, lwk2), lwk3);
    }
//==== Quick return in case of workspace query. ====
    if (lwork == -1) {
	work[1] = lwkopt;
	return;
    }
//==== Nothing to do ...
//... for an empty active block ... ====
    ns = 0;
    nd = 0;
    if (ktop > kbot) {
	return;
    }
//... nor for an empty deflation window. ====
    if (nw < 1) {
	return;
    }
//==== Machine constants ====
    safmin = Rlamch("SAFE MINIMUM");
    safmax = One / safmin;
    ulp = Rlamch("PRECISION");
    smlnum = safmin * ((REAL) double (n) / ulp);
//==== Setup deflation window ====
    jw = min(nw, kbot - ktop + 1);
    kwtop = kbot - jw + 1;
    if (kwtop == ktop) {
	s = Zero;
    } else {
	s = h[kwtop + (kwtop - 1) * ldh];
    }
    if (kbot == kwtop) {
//==== 1-by-1 deflation window: not much to do ====
	sr[kwtop] = h[kwtop + kwtop * ldh];
	si[kwtop] = Zero;
	ns = 1;
	nd = 0;
	mtemp1 = smlnum, mtemp2 = ulp * abs(h[kwtop + kwtop * ldh]);
	if (abs(s) <= max(mtemp1, mtemp2)) {
	    ns = 0;
	    nd = 1;
	    if (kwtop > ktop) {
		h[kwtop + (kwtop - 1) * ldh] = Zero;
	    }
	}
	return;
    }
//==== Convert to spike-triangular form.  (In case of a
//.    rare QR failure, this routine continues to do
//.    aggressive early deflation using that part of
//.    the deflation window that converged using INFQR
//.    here and there to keep track.) ====
    Rlacpy("U", jw, jw, &h[kwtop + kwtop * ldh], ldh, t, ldt);
    Rcopy(jw - 1, &h[kwtop + 1 + kwtop * ldh], ldh + 1, &t[ldt + 2], ldt + 1);
    Rlaset("A", jw, jw, Zero, One, v, ldv);
    nmin = iMlaenv(12, "Rlaqr3", "SV", jw, 1, jw, lwork);
    if (jw > nmin) {
	Rlaqr4(MTRUE, MTRUE, jw, 1, jw, t, ldt, &sr[kwtop], &si[kwtop], 1, jw, v, ldv, work, lwork, &infqr);
    } else {
	Rlahqr(MTRUE, MTRUE, jw, 1, jw, t, ldt, &sr[kwtop], &si[kwtop], 1, jw, v, ldv, &infqr);
    }
//==== DTREXC needs a clean margin near the diagonal ====
    for (j = 0; j < jw - 3; j++) {
	t[j + 2 + j * ldt] = Zero;
	t[j + 3 + j * ldt] = Zero;
    }
    if (jw > 2) {
	t[jw + (jw - 2) * ldt] = Zero;
    }
//==== Deflation detection loop ====
    ns = jw;
    ilst = infqr + 1;
  L20:
    if (ilst <= ns) {
	if (ns == 1) {
	    bulge = MFALSE;
	} else {
	    bulge = t[ns + (ns - 1) * ldt] != Zero;
	}
//==== Small spike tip test for deflation ====
	if (!bulge) {
//==== Real eigenvalue ====
	    foo = abs(t[ns + ns * ldt]);
	    if (foo == Zero) {
		foo = abs(s);
	    }
	    if (abs(s * v[ns * ldv + 1]) <= max(mtemp1 = smlnum, mtemp2 = ulp * foo)) {
//==== Deflatable ====
		--(ns);
	    } else {
//==== Undeflatable.   Move it up out of the way.
//.    (DTREXC can not fail in this case.) ====
		ifst = ns;
		Rtrexc("V", jw, t, ldt, v, ldv, &ifst, &ilst, &work[0], &info);
		ilst++;
	    }
	} else {
//==== Complex conjugate pair ====
	    foo = abs(t[ns + ns * ldt]) + sqrt(abs(t[ns + (ns - 1) * ldt])) * sqrt(abs(t[ns - 1 + ns * ldt]));
	    if (foo == Zero) {
		foo = abs(s);
	    }
	    mtemp1 = abs(s * v[ns * ldv + 1]), mtemp2 = abs(s * v[(ns - 1) * ldv + 1]);
	    mtemp3 = smlnum, mtemp4 = ulp * foo;
	    if (max(mtemp1, mtemp2) <= max(mtemp3, mtemp4)) {
//==== Deflatable ====
		ns = ns - 2;
	    } else {
//==== Undflatable. Move them up out of the way.
//.    Fortunately, DTREXC does the right thing with
//.    ILST in case of a rare exchange failure. ====
		ifst = ns;
		Rtrexc("V", jw, t, ldt, &v[0], ldv, &ifst, &ilst, &work[0], &info);
		ilst = ilst + 2;
	    }
	}
//==== End deflation detection loop ====
	goto L20;
    }
//==== Return to Hessenberg form ====
    if (ns == 0) {
	s = Zero;
    }
    if (ns < jw) {
//==== sorting diagonal blocks of T improves accuracy for
//.    graded matrices.  Bubble sort deals well with
//.    exchange failures. ====
	sorted = MFALSE;
	i = ns + 1;
      L30:
	if (sorted) {
	    goto L50;
	}
	sorted = MTRUE;
	kend = i - 1;
	i = infqr + 1;
	if (i == ns) {
	    k = i + 1;
	} else if (t[i + 1 + i * ldt] == Zero) {
	    k = i + 1;
	} else {
	    k = i + 2;
	}
      L40:
	if (k <= kend) {
	    if (k == i + 1) {
		evi = abs(t[i + i * ldt]);
	    } else {
		evi = abs(t[i + i * ldt]) + sqrt(abs(t[i + 1 + i * ldt])) * sqrt(abs(t[i + (i + 1) * ldt]));
	    }
	    if (k == kend) {
		evk = abs(t[k + k * ldt]);
	    } else if (t[k + 1 + k * ldt] == Zero) {
		evk = abs(t[k + k * ldt]);
	    } else {
		evk = abs(t[k + k * ldt]) + sqrt(abs(t[k + 1 + k * ldt])) * sqrt(abs(t[k + (k + 1) * ldt]));
	    }
	    if (evi >= evk) {
		i = k;
	    } else {
		sorted = MFALSE;
		ifst = i;
		ilst = k;
		Rtrexc("V", jw, t, ldt, &v[0], ldv, &ifst, &ilst, &work[0], &info);
		if (info == 0) {
		    i = ilst;
		} else {
		    i = k;
		}
	    }
	    if (i == kend) {
		k = i + 1;
	    } else if (t[i + 1 + i * ldt] == Zero) {
		k = i + 1;
	    } else {
		k = i + 2;
	    }
	    goto L40;
	}
	goto L30;
      L50:
	;
    }
//==== Restore shift/eigenvalue array from T ====
    i = jw;
  L60:
    if (i >= infqr + 1) {
	if (i == infqr + 1) {
	    sr[kwtop + i - 1] = t[i + i * ldt];
	    si[kwtop + i - 1] = Zero;
	    i--;
	} else if (t[i + (i - 1) * ldt] == Zero) {
	    sr[kwtop + i - 1] = t[i + i * ldt];
	    si[kwtop + i - 1] = Zero;
	    i--;
	} else {
	    aa = t[i - 1 + (i - 1) * ldt];
	    cc = t[i + (i - 1) * ldt];
	    bb = t[i - 1 + i * ldt];
	    dd = t[i + i * ldt];
	    Rlanv2(&aa, &bb, &cc, &dd, &sr[kwtop + i - 2], &si[kwtop + i - 2], &sr[kwtop + i - 1], &si[kwtop + i - 1], &cs, &sn);
	    i += -2;
	}
	goto L60;
    }
    if (ns < jw || s == Zero) {
	if (ns > 1 && s != Zero) {
//==== Reflect spike back into lower triangle ====
	    Rcopy(ns, &v[0], ldv, &work[0], 1);
	    beta = work[1];
	    Rlarfg(ns, &beta, &work[2], 1, &tau);
	    work[1] = One;

	    Rlaset("L", jw - 2, jw - 2, Zero, Zero, &t[ldt + 3], ldt);
	    Rlarf("L", ns, jw, work, 1, tau, t, ldt, &work[jw + 1]);
	    Rlarf("R", ns, ns, work, 1, tau, t, ldt, &work[jw + 1]);
	    Rlarf("R", jw, ns, work, 1, tau, v, ldv, &work[jw + 1]);
	    Rgehrd(jw, 1, ns, t, ldt, work, &work[jw + 1], lwork - jw, &info);
	}
//==== Copy updated reduced window into place ====
	if (kwtop > 1) {
	    h[kwtop + (kwtop - 1) * ldh] = s * v[ldv + 1];
	}
	Rlacpy("U", jw, jw, t, ldt, &h[kwtop + kwtop * ldh], ldh);
	Rcopy(jw - 1, &t[ldt + 2], ldt + 1, &h[kwtop + 1 + kwtop * ldh], ldh + 1);
//==== Accumulate orthogonal matrix in order update
//.    H and Z, if requested.  (A modified version
//.    of  DORGHR that accumulates block Householder
//.    transformations into V directly might be
//.    marginally more efficient than the following.) ====
	if (ns > 1 && s != Zero) {
	    Rorghr(jw, 1, ns, t, ldt, work, &work[jw + 1], lwork - jw, &info);
	    Rgemm("N", "N", jw, ns, ns, One, v, ldv, t, ldt, Zero, wv, ldwv);
	    Rlacpy("A", jw, ns, wv, ldwv, v, ldv);
	}
//==== Update vertical slab in H ====
	if (wantt) {
	    ltop = 1;
	} else {
	    ltop = ktop;
	}
	for (krow = ltop; krow <= kwtop - 1; krow = krow + nv) {
	    kln = min(nv, kwtop - krow);
	    Rgemm("N", "N", kln, jw, jw, One, &h[krow + kwtop * ldh], ldh, v, ldv, Zero, wv, ldwv);
	    Rlacpy("A", kln, jw, wv, ldwv, &h[krow + kwtop * ldh], ldh);
	}
//==== Update horizontal slab in H ====
	if (wantt) {
	    for (kcol = kbot + 1; kcol <= n; kcol = kcol + nh) {
		kln = min(nh, n - kcol + 1);
		Rgemm("C", "N", jw, kln, jw, One, v, ldv, &h[kwtop + kcol * ldh], ldh, Zero, t, ldt);
		Rlacpy("A", jw, kln, t, ldt, &h[kwtop + kcol * ldh], ldh);
	    }
	}
//==== Update vertical slab in Z ====
	if (wantz) {
	    for (krow = iloz; krow <= ihiz; krow = krow + nv) {
		kln = min(nv, ihiz - krow + 1);
		Rgemm("N", "N", kln, jw, jw, One, &z[krow + kwtop * ldz], ldz, v, ldv, Zero, wv, ldwv);
		Rlacpy("A", kln, jw, wv, ldwv, &z[krow + kwtop * ldz], ldz);
	    }
	}
    }
//==== Return the number of deflations ... ====
    nd = jw - ns;
//==== ... and the number of shifts. (Subtracting
//.    INFQR from the spike length takes care
//.    of the case of a rare QR failure while
//.    calculating eigenvalues of the deflation
//.    window.)  ====

    ns = ns - infqr;
//==== Return optimal workspace. ====
    work[1] = (REAL) double (lwkopt);
//==== End of DLAQR3 ====
    return;
}
