/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Claqr3.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Claqr3(LOGICAL wantt, LOGICAL wantz, INTEGER n,
	    INTEGER ktop, INTEGER kbot, INTEGER nw, COMPLEX * h,
	    INTEGER ldh, INTEGER iloz, INTEGER ihiz, COMPLEX * z,
	    INTEGER ldz, INTEGER * ns, INTEGER * nd, COMPLEX * sh, COMPLEX * v, INTEGER ldv, INTEGER nh, COMPLEX * t, INTEGER ldt,
	    INTEGER nv, COMPLEX * wv, INTEGER ldwv, COMPLEX * work, INTEGER lwork)
{
    INTEGER i, j;
    COMPLEX s;
    INTEGER jw;
    REAL foo;
    INTEGER kln;
    COMPLEX tau;
    INTEGER knt;
    REAL ulp;
    INTEGER lwk1, lwk2, lwk3;
    COMPLEX beta;
    INTEGER kcol, info, nmin, ifst, ilst, ltop, krow;
    INTEGER infqr;
    INTEGER kwtop;
    REAL safmin;
    REAL safmax, smlnum;
    INTEGER lwkopt;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;

//==== Estimate optimal workspace. ====
    jw = min(nw, kbot - ktop + 1);
    if (jw <= 2) {
	lwkopt = 1;
    } else {
//==== Workspace query call to ZGEHRD ====
	Cgehrd(jw, (INTEGER) 1, jw - 1, &t[0], ldt, &work[0], &work[0], (INTEGER) - 1, &info);
	lwk1 = (INTEGER) cast2double(work[1].real());
//==== Workspace query call to ZUNGHR ====
	Cunghr(jw, (INTEGER) 1, jw - 1, &t[0], ldt, &work[0], &work[0], (INTEGER) - 1, &info);
	lwk2 = (INTEGER) cast2double(work[1].real());
//==== Workspace query call to ZLAQR4 ====
	Claqr4(MTRUE, MTRUE, jw, 1, jw, &t[0], ldt, &sh[1], 1, jw, &v[0], ldv, &work[0], -1, &infqr);
	lwk3 = (INTEGER) cast2double(work[1].real());
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
    //Rlabad(&safmin, &safmax);
    ulp = Rlamch("PRECISION");
    smlnum = safmin * (n / ulp);
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
	sh[kwtop] = h[kwtop + kwtop * ldh];
	*ns = 1;
	*nd = 0;
	mtemp1 = smlnum, mtemp2 = ulp * Cabs1(h[kwtop + kwtop * ldh]);
	if (Cabs1(s) <= max(mtemp1, mtemp2)) {
	    *ns = 0;
	    *nd = 1;
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
    Clacpy("U", jw, jw, &h[kwtop + kwtop * ldh], ldh, &t[0], ldt);
    Ccopy(jw - 1, &h[kwtop + 1 + kwtop * ldh], ldh + 1, &t[ldt + 2], ldt + 1);
    Claset("A", jw, jw, Zero, One, &v[0], ldv);
    nmin = iMlaenv(12, "Claqr3", "SV", jw, 1, jw, lwork);
    if (jw > nmin) {
	Claqr4(MTRUE, MTRUE, jw, 1, jw, &t[0], ldt, &sh[kwtop], 1, jw, &v[0], ldv, &work[0], lwork, &infqr);
    } else {
	Clahqr(MTRUE, MTRUE, jw, 1, jw, &t[0], ldt, &sh[kwtop], 1, jw, &v[0], ldv, &infqr);
    }
//==== Deflation detection loop ====
    *ns = jw;
    ilst = infqr + 1;
    for (knt = infqr + 1; knt <= jw; knt++) {
//==== Small spike tip deflation test ====
	foo = Cabs1(t[(*ns) + (*ns) * ldt]);
	if (foo == Zero) {
	    foo = Cabs1(s);
	}
	mtemp1 = smlnum, mtemp2 = ulp * foo;
	if (Cabs1(s) * Cabs1(v[(*ns) * ldv + 1]) <= max(mtemp1, mtemp2)) {
//==== One more converged eigenvalue ====
	    --(ns);
	} else {
//==== One undflatable eigenvalue.  Move it up out of the
//.    way.   (ZTREXC can not fail in this case.) ====
	    ifst = (*ns);
	    Ctrexc("V", jw, &t[0], ldt, &v[0], ldv, ifst, ilst, &info);
	    ilst++;
	}
    }
//==== Return to Hessenberg form ====
    if (ns == 0) {
	s = Zero;
    }
    if ((*ns) < jw) {
//==== sorting the diagonal of T improves accuracy for
//.    graded matrices.  ====
	for (i = infqr + 1; i <= (*ns); i++) {
	    ifst = i;
	    for (j = i + 1; j <= (*ns); j++) {
		if (Cabs1(t[j + j * ldt]) > Cabs1(t[ifst + ifst * ldt])) {
		    ifst = j;
		}
	    }
	    ilst = i;
	    if (ifst != ilst) {
		Ctrexc("V", jw, &t[0], ldt, &v[0], ldv, ifst, ilst, &info);
	    }
	}
    }
//==== Restore shift/eigenvalue array from T ====
    for (i = infqr + 1; i <= jw; i++) {
	sh[kwtop + i - 1] = t[i + i * ldt];
    }
    if ((*ns) < jw || s == Zero) {
	if ((*ns) > 1 && s != Zero) {
//==== Reflect spike back into lower triangle ====
	    Ccopy((*ns), &v[0], ldv, &work[0], 1);
	    for (i = 0; i < (*ns); i++) {
		work[i] = conj(work[i]);
	    }
	    beta = work[1];
	    Clarfg((*ns), &beta, &work[2], 1, &tau);
	    work[1] = One;
	    Claset("L", jw - 2, jw - 2, Zero, Zero, &t[ldt + 3], ldt);
	    Clarf("L", (*ns), jw, &work[0], 1, conj(tau), &t[0], ldt, &work[jw + 1]);
	    Clarf("R", (*ns), (*ns), &work[0], 1, tau, &t[0], ldt, &work[jw + 1]);
	    Clarf("R", jw, (*ns), &work[0], 1, tau, &v[0], ldv, &work[jw + 1]);
	    Cgehrd(jw, 1, (*ns), &t[0], ldt, &work[0], &work[jw + 1], lwork - jw, &info);
	}
//==== Copy updated reduced window into place ====
	if (kwtop > 1) {
	    h[kwtop + (kwtop - 1) * ldh] = s * conj(v[ldv + 1]);
	}
	Clacpy("U", jw, jw, &t[0], ldt, &h[kwtop + kwtop * ldh], ldh);
	Ccopy(jw - 1, &t[ldt + 2], ldt + 1, &h[kwtop + 1 + kwtop * ldh], ldh + 1);
//==== Accumulate orthogonal matrix in order update
//.    H and Z, if requested.  (A modified version
//.    of  ZUNGHR that accumulates block Householder
//.    transformations into V directly might be
//.    marginally more efficient than the following.) ====
	if ((*ns) > 1 && s != Zero) {
	    Cunghr(jw, 1, (*ns), &t[0], ldt, &work[0], &work[jw + 1], lwork - jw, &info);
	    Cgemm("N", "N", jw, (*ns), (*ns), (COMPLEX) One, &v[0], ldv, &t[0], ldt, (COMPLEX) Zero, &wv[0], ldwv);
	    Clacpy("A", jw, (*ns), &wv[0], ldwv, &v[0], ldv);
	}
//==== Update vertical slab in H ====
	if (wantt) {
	    ltop = 1;
	} else {
	    ltop = ktop;
	}
	for (krow = ltop; krow <= kwtop - 1; krow = krow + nv) {
	    kln = min(nv, kwtop - krow);
	    Cgemm("N", "N", kln, jw, jw, One, &h[krow + kwtop * ldh], ldh, &v[0], ldv, Zero, &wv[0], ldwv);
	    Clacpy("A", kln, jw, &wv[0], ldwv, &h[krow + kwtop * ldh], ldh);
	}
//==== Update horizontal slab in H ====
	if (wantt) {
	    for (kcol = kbot + 1; kcol <= n; kcol = kcol + nh) {
		kln = min(nh, n - kcol + 1);
		Cgemm("C", "N", jw, kln, jw, (COMPLEX) One, &v[0], ldv, &h[kwtop + kcol * ldh], ldh, (COMPLEX) Zero, &t[0], ldt);
		Clacpy("A", jw, kln, &t[0], ldt, &h[kwtop + kcol * ldh], ldh);
	    }
	}
//==== Update vertical slab in Z ====
	if (wantz) {
	    for (krow = iloz; krow <= ihiz; krow = krow + nv) {
		kln = min(nv, ihiz - krow + 1);
		Cgemm("N", "N", kln, jw, jw, One, &z[krow + kwtop * ldz], ldz, &v[0], ldv, Zero, &wv[0]
		      , ldwv);
		Clacpy("A", kln, jw, &wv[0], ldwv, &z[krow + kwtop * ldz], ldz);
	    }
	}
    }
//==== Return the number of deflations ... ====
    (*nd) = jw - (*ns);
//==== ... and the number of shifts. (Subtracting */
//.    INFQR from the spike length takes care */
//.    of the case of a rare QR failure while */
//.    calculating eigenvalues of the deflation */
//.    window.)  ==== */
    (*ns) = (*ns) - infqr;
//==== Return optimal workspace. ====
    work[1] = lwkopt;
//==== End of ZLAQR3 ====
    return;
}
