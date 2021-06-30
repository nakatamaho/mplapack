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

void Rlaqr3(bool const wantt, bool const wantz, INTEGER const n, INTEGER const ktop, INTEGER const kbot, INTEGER const nw, REAL *h, INTEGER const ldh, INTEGER const iloz, INTEGER const ihiz, REAL *z, INTEGER const ldz, INTEGER &ns, INTEGER &nd, REAL *sr, REAL *si, REAL *v, INTEGER const ldv, INTEGER const nh, REAL *t, INTEGER const ldt, INTEGER const nv, REAL *wv, INTEGER const ldwv, REAL *work, INTEGER const lwork) {
    INTEGER jw = 0;
    INTEGER lwkopt = 0;
    INTEGER info = 0;
    INTEGER lwk1 = 0;
    INTEGER lwk2 = 0;
    INTEGER infqr = 0;
    INTEGER lwk3 = 0;
    const REAL one = 1.0;
    REAL safmin = 0.0;
    REAL safmax = 0.0;
    REAL ulp = 0.0;
    REAL smlnum = 0.0;
    INTEGER kwtop = 0;
    const REAL zero = 0.0;
    REAL s = 0.0;
    INTEGER nmin = 0;
    INTEGER j = 0;
    INTEGER ilst = 0;
    bool bulge = false;
    REAL foo = 0.0;
    INTEGER ifst = 0;
    bool sorted = false;
    INTEGER i = 0;
    INTEGER kend = 0;
    INTEGER k = 0;
    REAL evi = 0.0;
    REAL evk = 0.0;
    REAL aa = 0.0;
    REAL cc = 0.0;
    REAL bb = 0.0;
    REAL dd = 0.0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    REAL beta = 0.0;
    REAL tau = 0.0;
    INTEGER ltop = 0;
    INTEGER krow = 0;
    INTEGER kln = 0;
    INTEGER kcol = 0;
    //
    //     ==== Estimate optimal workspace. ====
    //
    jw = min(nw, kbot - ktop + 1);
    if (jw <= 2) {
        lwkopt = 1;
    } else {
        //
        //        ==== Workspace query call to Rgehrd ====
        //
        Rgehrd(jw, 1, jw - 1, t, ldt, work, work, -1, info);
        lwk1 = castINTEGER(work[1 - 1]);
        //
        //        ==== Workspace query call to Rormhr ====
        //
        Rormhr("R", "N", jw, jw, 1, jw - 1, t, ldt, work, v, ldv, work, -1, info);
        lwk2 = castINTEGER(work[1 - 1]);
        //
        //        ==== Workspace query call to Rlaqr4 ====
        //
        Rlaqr4(true, true, jw, 1, jw, t, ldt, sr, si, 1, jw, v, ldv, work, -1, infqr);
        lwk3 = castINTEGER(work[1 - 1]);
        //
        //        ==== Optimal workspace ====
        //
        lwkopt = max(jw + max(lwk1, lwk2), lwk3);
    }
    //
    //     ==== Quick return in case of workspace query. ====
    //
    if (lwork == -1) {
        work[1 - 1] = castREAL(lwkopt);
        return;
    }
    //
    //     ==== Nothing to do ...
    //     ... for an empty active block ... ====
    ns = 0;
    nd = 0;
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
    safmin = Rlamch("SAFE MINIMUM");
    safmax = one / safmin;
    ulp = Rlamch("PRECISION");
    smlnum = safmin * (castREAL(n) / ulp);
    //
    //     ==== Setup deflation window ====
    //
    jw = min(nw, kbot - ktop + 1);
    kwtop = kbot - jw + 1;
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
        sr[kwtop - 1] = h[(kwtop - 1) + (kwtop - 1) * ldh];
        si[kwtop - 1] = zero;
        ns = 1;
        nd = 0;
        if (abs(s) <= max(smlnum, REAL(ulp * abs(h[(kwtop - 1) + (kwtop - 1) * ldh])))) {
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
    Rlacpy("U", jw, jw, &h[(kwtop - 1) + (kwtop - 1) * ldh], ldh, t, ldt);
    Rcopy(jw - 1, &h[((kwtop + 1) - 1) + (kwtop - 1) * ldh], ldh + 1, &t[(2 - 1)], ldt + 1);
    //
    Rlaset("A", jw, jw, zero, one, v, ldv);
    nmin = iMlaenv(12, "Rlaqr3", "SV", jw, 1, jw, lwork);
    if (jw > nmin) {
        Rlaqr4(true, true, jw, 1, jw, t, ldt, &sr[kwtop - 1], &si[kwtop - 1], 1, jw, v, ldv, work, lwork, infqr);
    } else {
        Rlahqr(true, true, jw, 1, jw, t, ldt, &sr[kwtop - 1], &si[kwtop - 1], 1, jw, v, ldv, infqr);
    }
    //
    //     ==== Rtrexc needs a clean margin near the diagonal ====
    //
    for (j = 1; j <= jw - 3; j = j + 1) {
        t[((j + 2) - 1) + (j - 1) * ldt] = zero;
        t[((j + 3) - 1) + (j - 1) * ldt] = zero;
    }
    if (jw > 2) {
        t[(jw - 1) + ((jw - 2) - 1) * ldt] = zero;
    }
    //
    //     ==== Deflation detection loop ====
    //
    ns = jw;
    ilst = infqr + 1;
statement_20:
    if (ilst <= ns) {
        if (ns == 1) {
            bulge = false;
        } else {
            bulge = t[(ns - 1) + ((ns - 1) - 1) * ldt] != zero;
        }
        //
        //        ==== Small spike tip test for deflation ====
        //
        if (!bulge) {
            //
            //           ==== Real eigenvalue ====
            //
            foo = abs(t[(ns - 1) + (ns - 1) * ldt]);
            if (foo == zero) {
                foo = abs(s);
            }
            if (abs(s * v[(ns - 1) * ldv]) <= max(smlnum, REAL(ulp * foo))) {
                //
                //              ==== Deflatable ====
                //
                ns = ns - 1;
            } else {
                //
                //              ==== Undeflatable.   Move it up out of the way.
                //              .    (Rtrexc can not fail in this case.) ====
                //
                ifst = ns;
                Rtrexc("V", jw, t, ldt, v, ldv, ifst, ilst, work, info);
                ilst++;
            }
        } else {
            //
            //           ==== Complex conjugate pair ====
            //
            foo = abs(t[(ns - 1) + (ns - 1) * ldt]) + sqrt(abs(t[(ns - 1) + ((ns - 1) - 1) * ldt])) * sqrt(abs(t[((ns - 1) - 1) + (ns - 1) * ldt]));
            if (foo == zero) {
                foo = abs(s);
            }
            if (max(abs(s * v[(ns - 1) * ldv]), abs(s * v[((ns - 1) - 1) * ldv])) <= max(smlnum, REAL(ulp * foo))) {
                //
                //              ==== Deflatable ====
                //
                ns = ns - 2;
            } else {
                //
                //              ==== Undeflatable. Move them up out of the way.
                //              .    Fortunately, Rtrexc does the right thing with
                //              .    ILST in case of a rare exchange failure. ====
                //
                ifst = ns;
                Rtrexc("V", jw, t, ldt, v, ldv, ifst, ilst, work, info);
                ilst += 2;
            }
        }
        //
        //        ==== End deflation detection loop ====
        //
        goto statement_20;
    }
    //
    //        ==== Return to Hessenberg form ====
    //
    if (ns == 0) {
        s = zero;
    }
    //
    if (ns < jw) {
        //
        //        ==== sorting diagonal blocks of T improves accuracy for
        //        .    graded matrices.  Bubble sort deals well with
        //        .    exchange failures. ====
        //
        sorted = false;
        i = ns + 1;
    statement_30:
        if (sorted) {
            goto statement_50;
        }
        sorted = true;
        //
        kend = i - 1;
        i = infqr + 1;
        if (i == ns) {
            k = i + 1;
        } else if (t[((i + 1) - 1) + (i - 1) * ldt] == zero) {
            k = i + 1;
        } else {
            k = i + 2;
        }
    statement_40:
        if (k <= kend) {
            if (k == i + 1) {
                evi = abs(t[(i - 1) + (i - 1) * ldt]);
            } else {
                evi = abs(t[(i - 1) + (i - 1) * ldt]) + sqrt(abs(t[((i + 1) - 1) + (i - 1) * ldt])) * sqrt(abs(t[(i - 1) + ((i + 1) - 1) * ldt]));
            }
            //
            if (k == kend) {
                evk = abs(t[(k - 1) + (k - 1) * ldt]);
            } else if (t[((k + 1) - 1) + (k - 1) * ldt] == zero) {
                evk = abs(t[(k - 1) + (k - 1) * ldt]);
            } else {
                evk = abs(t[(k - 1) + (k - 1) * ldt]) + sqrt(abs(t[((k + 1) - 1) + (k - 1) * ldt])) * sqrt(abs(t[(k - 1) + ((k + 1) - 1) * ldt]));
            }
            //
            if (evi >= evk) {
                i = k;
            } else {
                sorted = false;
                ifst = i;
                ilst = k;
                Rtrexc("V", jw, t, ldt, v, ldv, ifst, ilst, work, info);
                if (info == 0) {
                    i = ilst;
                } else {
                    i = k;
                }
            }
            if (i == kend) {
                k = i + 1;
            } else if (t[((i + 1) - 1) + (i - 1) * ldt] == zero) {
                k = i + 1;
            } else {
                k = i + 2;
            }
            goto statement_40;
        }
        goto statement_30;
    statement_50:;
    }
    //
    //     ==== Restore shift/eigenvalue array from T ====
    //
    i = jw;
statement_60:
    if (i >= infqr + 1) {
        if (i == infqr + 1) {
            sr[(kwtop + i - 1) - 1] = t[(i - 1) + (i - 1) * ldt];
            si[(kwtop + i - 1) - 1] = zero;
            i = i - 1;
        } else if (t[(i - 1) + ((i - 1) - 1) * ldt] == zero) {
            sr[(kwtop + i - 1) - 1] = t[(i - 1) + (i - 1) * ldt];
            si[(kwtop + i - 1) - 1] = zero;
            i = i - 1;
        } else {
            aa = t[((i - 1) - 1) + ((i - 1) - 1) * ldt];
            cc = t[(i - 1) + ((i - 1) - 1) * ldt];
            bb = t[((i - 1) - 1) + (i - 1) * ldt];
            dd = t[(i - 1) + (i - 1) * ldt];
            Rlanv2(aa, bb, cc, dd, sr[(kwtop + i - 2) - 1], si[(kwtop + i - 2) - 1], sr[(kwtop + i - 1) - 1], si[(kwtop + i - 1) - 1], cs, sn);
            i = i - 2;
        }
        goto statement_60;
    }
    //
    if (ns < jw || s == zero) {
        if (ns > 1 && s != zero) {
            //
            //           ==== Reflect spike back into lower triangle ====
            //
            Rcopy(ns, v, ldv, work, 1);
            beta = work[1 - 1];
            Rlarfg(ns, beta, &work[2 - 1], 1, tau);
            work[1 - 1] = one;
            //
            Rlaset("L", jw - 2, jw - 2, zero, zero, &t[(3 - 1)], ldt);
            //
            Rlarf("L", ns, jw, work, 1, tau, t, ldt, &work[(jw + 1) - 1]);
            Rlarf("R", ns, ns, work, 1, tau, t, ldt, &work[(jw + 1) - 1]);
            Rlarf("R", jw, ns, work, 1, tau, v, ldv, &work[(jw + 1) - 1]);
            //
            Rgehrd(jw, 1, ns, t, ldt, work, &work[(jw + 1) - 1], lwork - jw, info);
        }
        //
        //        ==== Copy updated reduced window into place ====
        //
        if (kwtop > 1) {
            h[(kwtop - 1) + ((kwtop - 1) - 1) * ldh] = s * v[(1 - 1)];
        }
        Rlacpy("U", jw, jw, t, ldt, &h[(kwtop - 1) + (kwtop - 1) * ldh], ldh);
        Rcopy(jw - 1, &t[(2 - 1)], ldt + 1, &h[((kwtop + 1) - 1) + (kwtop - 1) * ldh], ldh + 1);
        //
        //        ==== Accumulate orthogonal matrix in order update
        //        .    H and Z, if requested.  ====
        //
        if (ns > 1 && s != zero) {
            Rormhr("R", "N", jw, ns, 1, ns, t, ldt, work, v, ldv, &work[(jw + 1) - 1], lwork - jw, info);
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
            Rgemm("N", "N", kln, jw, jw, one, &h[(krow - 1) + (kwtop - 1) * ldh], ldh, v, ldv, zero, wv, ldwv);
            Rlacpy("A", kln, jw, wv, ldwv, &h[(krow - 1) + (kwtop - 1) * ldh], ldh);
        }
        //
        //        ==== Update horizontal slab in H ====
        //
        if (wantt) {
            for (kcol = kbot + 1; kcol <= n; kcol = kcol + nh) {
                kln = min(nh, n - kcol + 1);
                Rgemm("C", "N", jw, kln, jw, one, v, ldv, &h[(kwtop - 1) + (kcol - 1) * ldh], ldh, zero, t, ldt);
                Rlacpy("A", jw, kln, t, ldt, &h[(kwtop - 1) + (kcol - 1) * ldh], ldh);
            }
        }
        //
        //        ==== Update vertical slab in Z ====
        //
        if (wantz) {
            for (krow = iloz; krow <= ihiz; krow = krow + nv) {
                kln = min(nv, ihiz - krow + 1);
                Rgemm("N", "N", kln, jw, jw, one, &z[(krow - 1) + (kwtop - 1) * ldz], ldz, v, ldv, zero, wv, ldwv);
                Rlacpy("A", kln, jw, wv, ldwv, &z[(krow - 1) + (kwtop - 1) * ldz], ldz);
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
    work[1 - 1] = castREAL(lwkopt);
    //
    //     ==== End of Rlaqr3 ====
    //
}
