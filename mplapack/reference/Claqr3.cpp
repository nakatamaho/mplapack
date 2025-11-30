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

inline REAL cabs1(COMPLEX cdum) { return (abs(cdum.real()) + abs(cdum.imag())); }

void Claqr3(bool const wantt, bool const wantz, INTEGER const n, INTEGER const ktop, INTEGER const kbot, INTEGER const nw, COMPLEX *h, INTEGER const ldh, INTEGER const iloz, INTEGER const ihiz, COMPLEX *z, INTEGER const ldz, INTEGER &ns, INTEGER &nd, COMPLEX *sh, COMPLEX *v, INTEGER const ldv, INTEGER const nh, COMPLEX *t, INTEGER const ldt, INTEGER const nv, COMPLEX *wv, INTEGER const ldwv, COMPLEX *work, INTEGER const lwork) {
    //
    COMPLEX cdum = 0.0;
    //
    //
    INTEGER jw = min(nw, kbot - ktop + 1);
    INTEGER lwkopt = 0;
    INTEGER info = 0;
    INTEGER lwk1 = 0;
    INTEGER lwk2 = 0;
    INTEGER infqr = 0;
    INTEGER lwk3 = 0;
    if (jw <= 2) {
        lwkopt = 1;
    } else {
        //
        //
        Cgehrd(jw, 1, jw - 1, t, ldt, work, work, -1, info);
        lwk1 = castINTEGER(work[1 - 1].real());
        //
        //
        Cunmhr("R", "N", jw, jw, 1, jw - 1, t, ldt, work, v, ldv, work, -1, info);
        lwk2 = castINTEGER(work[1 - 1].real());
        //
        //
        Claqr4(true, true, jw, 1, jw, t, ldt, sh, 1, jw, v, ldv, work, -1, infqr);
        lwk3 = castINTEGER(work[1 - 1].real());
        //
        //
        lwkopt = max(jw + max(lwk1, lwk2), lwk3);
    }
    //
    //
    if (lwork == -1) {
        work[1 - 1] = COMPLEX(lwkopt, 0.0);
        return;
    }
    //
    ns = 0;
    nd = 0;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    work[1 - 1] = one;
    if (ktop > kbot) {
        return;
    }
    if (nw < 1) {
        return;
    }
    //
    //
    REAL safmin = Rlamch("SAFE MINIMUM");
    const REAL rone = 1.0;
    REAL safmax = rone / safmin;
    REAL ulp = Rlamch("PRECISION");
    REAL smlnum = safmin * (castREAL(n) / ulp);
    //
    //
    jw = min(nw, kbot - ktop + 1);
    INTEGER kwtop = kbot - jw + 1;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    COMPLEX s = 0.0;
    if (kwtop == ktop) {
        s = zero;
    } else {
        s = h[(kwtop - 1) + ((kwtop - 1) - 1) * ldh];
    }
    //
    if (kbot == kwtop) {
        //
        //
        sh[kwtop - 1] = h[(kwtop - 1) + (kwtop - 1) * ldh];
        ns = 1;
        nd = 0;
        if (cabs1(s) <= max(smlnum, REAL(ulp * cabs1(h[(kwtop - 1) + (kwtop - 1) * ldh])))) {
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
    // .    rare QR failure, this routine continues to do
    // .    aggressive early deflation using that part of
    // .    the deflation window that converged using INFQR
    //
    Clacpy("U", jw, jw, &h[(kwtop - 1) + (kwtop - 1) * ldh], ldh, t, ldt);
    Ccopy(jw - 1, &h[((kwtop + 1) - 1) + (kwtop - 1) * ldh], ldh + 1, &t[(2 - 1)], ldt + 1);
    //
    Claset("A", jw, jw, zero, one, v, ldv);
    INTEGER nmin = iMlaenv(12, "Claqr3", "SV", jw, 1, jw, lwork);
    if (jw > nmin) {
        Claqr4(true, true, jw, 1, jw, t, ldt, &sh[kwtop - 1], 1, jw, v, ldv, work, lwork, infqr);
    } else {
        Clahqr(true, true, jw, 1, jw, t, ldt, &sh[kwtop - 1], 1, jw, v, ldv, infqr);
    }
    //
    //
    ns = jw;
    INTEGER ilst = infqr + 1;
    INTEGER knt = 0;
    REAL foo = 0.0;
    const REAL rzero = 0.0;
    INTEGER ifst = 0;
    for (knt = infqr + 1; knt <= jw; knt = knt + 1) {
        //
        //
        foo = cabs1(t[(ns - 1) + (ns - 1) * ldt]);
        if (foo == rzero) {
            foo = cabs1(s);
        }
        if (cabs1(s) * cabs1(v[(ns - 1) * ldv]) <= max(smlnum, REAL(ulp * foo))) {
            //
            //
            ns = ns - 1;
        } else {
            //
            //
            ifst = ns;
            Ctrexc("V", jw, t, ldt, v, ldv, ifst, ilst, info);
            ilst++;
        }
    }
    //
    //
    if (ns == 0) {
        s = zero;
    }
    //
    INTEGER i = 0;
    INTEGER j = 0;
    if (ns < jw) {
        //
        //
        for (i = infqr + 1; i <= ns; i = i + 1) {
            ifst = i;
            for (j = i + 1; j <= ns; j = j + 1) {
                if (cabs1(t[(j - 1) + (j - 1) * ldt]) > cabs1(t[(ifst - 1) + (ifst - 1) * ldt])) {
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
        //
        if (kwtop > 1) {
            h[(kwtop - 1) + ((kwtop - 1) - 1) * ldh] = s * conj(v[0]);
        }
        Clacpy("U", jw, jw, t, ldt, &h[(kwtop - 1) + (kwtop - 1) * ldh], ldh);
        Ccopy(jw - 1, &t[(2 - 1)], ldt + 1, &h[((kwtop + 1) - 1) + (kwtop - 1) * ldh], ldh + 1);
        //
        //
        if (ns > 1 && s != zero) {
            Cunmhr("R", "N", jw, ns, 1, ns, t, ldt, work, v, ldv, &work[(jw + 1) - 1], lwork - jw, info);
        }
        //
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
        //
        if (wantt) {
            for (kcol = kbot + 1; kcol <= n; kcol = kcol + nh) {
                kln = min(nh, n - kcol + 1);
                Cgemm("C", "N", jw, kln, jw, one, v, ldv, &h[(kwtop - 1) + (kcol - 1) * ldh], ldh, zero, t, ldt);
                Clacpy("A", jw, kln, t, ldt, &h[(kwtop - 1) + (kcol - 1) * ldh], ldh);
            }
        }
        //
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
    //
    nd = jw - ns;
    //
    // .    INFQR from the spike length takes care
    // .    of the case of a rare QR failure while
    // .    calculating eigenvalues of the deflation
    //
    ns = ns - infqr;
    //
    //
    work[1 - 1] = COMPLEX(lwkopt, 0.0);
    //
    //
}
