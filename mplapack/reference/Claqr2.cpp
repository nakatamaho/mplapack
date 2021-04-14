/*
 * Copyright (c) 2008-2021
 *      Nakata, Maho
 *      All rights reserved.
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

#include <mpblas.h>
#include <mplapack.h>

void Claqr2(bool const wantt, bool const wantz, INTEGER const n, INTEGER const ktop, INTEGER const kbot, INTEGER const nw, COMPLEX *h, INTEGER const ldh, INTEGER const iloz, INTEGER const ihiz, COMPLEX *z, INTEGER const ldz, INTEGER &ns, INTEGER &nd, COMPLEX *sh, COMPLEX *v, INTEGER const ldv, INTEGER const nh, COMPLEX *t, INTEGER const ldt, INTEGER const nv, COMPLEX *wv, INTEGER const ldwv, COMPLEX *work, INTEGER const lwork) {
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  ================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    COMPLEX cdum = 0.0;
    abs1[cdum - 1] = abs(cdum.real()) + abs(cdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     ==== Estimate optimal workspace. ====
    //
    INTEGER jw = min(nw, kbot - ktop + 1);
    INTEGER lwkopt = 0;
    INTEGER info = 0;
    INTEGER lwk1 = 0;
    INTEGER lwk2 = 0;
    if (jw <= 2) {
        lwkopt = 1;
    } else {
        //
        //        ==== Workspace query call to Cgehrd ====
        //
        Cgehrd(jw, 1, jw - 1, t, ldt, work, work, -1, info);
        lwk1 = int(work[1 - 1]);
        //
        //        ==== Workspace query call to Cunmhr ====
        //
        Cunmhr("R", "N", jw, jw, 1, jw - 1, t, ldt, work, v, ldv, work, -1, info);
        lwk2 = int(work[1 - 1]);
        //
        //        ==== Optimal workspace ====
        //
        lwkopt = jw + max(lwk1, lwk2);
    }
    //
    //     ==== Quick return in case of workspace query. ====
    //
    if (lwork == -1) {
        work[1 - 1] = COMPLEX(lwkopt, 0);
        return;
    }
    //
    //     ==== Nothing to do ...
    //     ... for an empty active block ... ====
    ns = 0;
    nd = 0;
    const COMPLEX one = (1.0, 0.0);
    work[1 - 1] = one;
    if (ktop > kbot) {
        return;
    }
    //     ... nor for an empty deflation window. ====
    if (nw < 1) {
        return;
    }
    //
    //     ==== Machine constants ====
    //
    REAL safmin = Rlamch("SAFE MINIMUM");
    const REAL rone = 1.0;
    REAL safmax = rone / safmin;
    Rlabad(safmin, safmax);
    REAL ulp = Rlamch("PRECISION");
    REAL smlnum = safmin * (n.real() / ulp);
    //
    //     ==== Setup deflation window ====
    //
    jw = min(nw, kbot - ktop + 1);
    INTEGER kwtop = kbot - jw + 1;
    const COMPLEX zero = (0.0, 0.0);
    COMPLEX s = 0.0;
    if (kwtop == ktop) {
        s = zero;
    } else {
        s = h[(kwtop - 1) + ((kwtop - 1) - 1) * ldh];
    }
    //
    if (kbot == kwtop) {
        //
        //        ==== 1-by-1 deflation window: not much to do ====
        //
        sh[kwtop - 1] = h[(kwtop - 1) + (kwtop - 1) * ldh];
        ns = 1;
        nd = 0;
        if (abs1[s - 1] <= max(smlnum, ulp * abs1[h[(kwtop - 1) + (kwtop - 1) * ldh] - 1])) {
            ns = 0;
            nd = 1;
            if (kwtop > ktop) {
                h[(kwtop - 1) + ((kwtop - 1) - 1) * ldh] = zero;
            }
        }
        work[1 - 1] = one;
        return;
    }
    //
    //     ==== Convert to spike-triangular form.  (In case of a
    //     .    rare QR failure, this routine continues to do
    //     .    aggressive early deflation using that part of
    //     .    the deflation window that converged using INFQR
    //     .    here and there to keep track.) ====
    //
    Clacpy("U", jw, jw, &h[(kwtop - 1) + (kwtop - 1) * ldh], ldh, t, ldt);
    Ccopy(jw - 1, &h[((kwtop + 1) - 1) + (kwtop - 1) * ldh], ldh + 1, &t[(2 - 1)], ldt + 1);
    //
    Claset("A", jw, jw, zero, one, v, ldv);
    INTEGER infqr = 0;
    Clahqr(true, true, jw, 1, jw, t, ldt, sh[kwtop - 1], 1, jw, v, ldv, infqr);
    //
    //     ==== Deflation detection loop ====
    //
    ns = jw;
    INTEGER ilst = infqr + 1;
    INTEGER knt = 0;
    REAL foo = 0.0;
    const REAL rzero = 0.0;
    INTEGER ifst = 0;
    for (knt = infqr + 1; knt <= jw; knt = knt + 1) {
        //
        //        ==== Small spike tip deflation test ====
        //
        foo = abs1[t[(ns - 1) + (ns - 1) * ldt] - 1];
        if (foo == rzero) {
            foo = abs1[s - 1];
        }
        if (abs1[s - 1] * abs1[v[(ns - 1) * ldv] - 1] <= max(smlnum, ulp * foo)) {
            //
            //           ==== One more converged eigenvalue ====
            //
            ns = ns - 1;
        } else {
            //
            //           ==== One undeflatable eigenvalue.  Move it up out of the
            //           .    way.   (Ctrexc can not fail in this case.) ====
            //
            ifst = ns;
            Ctrexc("V", jw, t, ldt, v, ldv, ifst, ilst, info);
            ilst++;
        }
    }
    //
    //        ==== Return to Hessenberg form ====
    //
    if (ns == 0) {
        s = zero;
    }
    //
    INTEGER i = 0;
    INTEGER j = 0;
    if (ns < jw) {
        //
        //        ==== sorting the diagonal of T improves accuracy for
        //        .    graded matrices.  ====
        //
        for (i = infqr + 1; i <= ns; i = i + 1) {
            ifst = i;
            for (j = i + 1; j <= ns; j = j + 1) {
                if (abs1[t[(j - 1) + (j - 1) * ldt] - 1] > abs1[t[(ifst - 1) + (ifst - 1) * ldt] - 1]) {
                    ifst = j;
                }
            }
            ilst = i;
            if (ifst != ilst) {
                Ctrexc("V", jw, t, ldt, v, ldv, ifst, ilst, info);
            }
        }
    }
    //
    //     ==== Restore shift/eigenvalue array from T ====
    //
    for (i = infqr + 1; i <= jw; i = i + 1) {
        sh[(kwtop + i - 1) - 1] = t[(i - 1) + (i - 1) * ldt];
    }
    //
    COMPLEX beta = 0.0;
    COMPLEX tau = 0.0;
    INTEGER ltop = 0;
    INTEGER krow = 0;
    INTEGER kln = 0;
    INTEGER kcol = 0;
    if (ns < jw || s == zero) {
        if (ns > 1 && s != zero) {
            //
            //           ==== Reflect spike back into lower triangle ====
            //
            Ccopy(ns, v, ldv, work, 1);
            for (i = 1; i <= ns; i = i + 1) {
                work[i - 1] = conj(work[i - 1]);
            }
            beta = work[1 - 1];
            Clarfg(ns, beta, &work[2 - 1], 1, tau);
            work[1 - 1] = one;
            //
            Claset("L", jw - 2, jw - 2, zero, zero, &t[(3 - 1)], ldt);
            //
            Clarf("L", ns, jw, work, 1, conj(tau), t, ldt, &work[(jw + 1) - 1]);
            Clarf("R", ns, ns, work, 1, tau, t, ldt, &work[(jw + 1) - 1]);
            Clarf("R", jw, ns, work, 1, tau, v, ldv, &work[(jw + 1) - 1]);
            //
            Cgehrd(jw, 1, ns, t, ldt, work, &work[(jw + 1) - 1], lwork - jw, info);
        }
        //
        //        ==== Copy updated reduced window into place ====
        //
        if (kwtop > 1) {
            h[(kwtop - 1) + ((kwtop - 1) - 1) * ldh] = s * conj(v[(1 - 1)]);
        }
        Clacpy("U", jw, jw, t, ldt, &h[(kwtop - 1) + (kwtop - 1) * ldh], ldh);
        Ccopy(jw - 1, &t[(2 - 1)], ldt + 1, &h[((kwtop + 1) - 1) + (kwtop - 1) * ldh], ldh + 1);
        //
        //        ==== Accumulate orthogonal matrix in order update
        //        .    H and Z, if requested.  ====
        //
        if (ns > 1 && s != zero) {
            Cunmhr("R", "N", jw, ns, 1, ns, t, ldt, work, v, ldv, &work[(jw + 1) - 1], lwork - jw, info);
        }
        //
        //        ==== Update vertical slab in H ====
        //
        if (wantt) {
            ltop = 1;
        } else {
            ltop = ktop;
        }
        for (krow = ltop; krow <= kwtop - 1; krow = krow + nv) {
            kln = min(nv, kwtop - krow);
            Cgemm("N", "N", kln, jw, jw, one, &h[(krow - 1) + (kwtop - 1) * ldh], ldh, v, ldv, zero, wv, ldwv);
            Clacpy("A", kln, jw, wv, ldwv, &h[(krow - 1) + (kwtop - 1) * ldh], ldh);
        }
        //
        //        ==== Update horizontal slab in H ====
        //
        if (wantt) {
            for (kcol = kbot + 1; kcol <= n; kcol = kcol + nh) {
                kln = min(nh, n - kcol + 1);
                Cgemm("C", "N", jw, kln, jw, one, v, ldv, &h[(kwtop - 1) + (kcol - 1) * ldh], ldh, zero, t, ldt);
                Clacpy("A", jw, kln, t, ldt, &h[(kwtop - 1) + (kcol - 1) * ldh], ldh);
            }
        }
        //
        //        ==== Update vertical slab in Z ====
        //
        if (wantz) {
            for (krow = iloz; krow <= ihiz; krow = krow + nv) {
                kln = min(nv, ihiz - krow + 1);
                Cgemm("N", "N", kln, jw, jw, one, &z[(krow - 1) + (kwtop - 1) * ldz], ldz, v, ldv, zero, wv, ldwv);
                Clacpy("A", kln, jw, wv, ldwv, &z[(krow - 1) + (kwtop - 1) * ldz], ldz);
            }
        }
    }
    //
    //     ==== Return the number of deflations ... ====
    //
    nd = jw - ns;
    //
    //     ==== ... and the number of shifts. (Subtracting
    //     .    INFQR from the spike length takes care
    //     .    of the case of a rare QR failure while
    //     .    calculating eigenvalues of the deflation
    //     .    window.)  ====
    //
    ns = ns - infqr;
    //
    //      ==== Return optimal workspace. ====
    //
    work[1 - 1] = COMPLEX(lwkopt, 0);
    //
    //     ==== End of Claqr2 ====
    //
}
