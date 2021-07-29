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

inline REAL abs1(COMPLEX x) { return abs(x.real()) + abs(x.imag()); }

void Chgeqz(const char *job, const char *compq, const char *compz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, COMPLEX *h, INTEGER const ldh, COMPLEX *t, INTEGER const ldt, COMPLEX *alpha, COMPLEX *beta, COMPLEX *q, INTEGER const ldq, COMPLEX *z, INTEGER const ldz, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER &info) {
    COMPLEX x = 0.0;
    bool ilschr = false;
    INTEGER ischur = 0;
    bool ilq = false;
    INTEGER icompq = 0;
    bool ilz = false;
    INTEGER icompz = 0;
    bool lquery = false;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER in = 0;
    REAL safmin = 0.0;
    REAL ulp = 0.0;
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    REAL atol = 0.0;
    REAL btol = 0.0;
    const REAL one = 1.0;
    REAL ascale = 0.0;
    REAL bscale = 0.0;
    INTEGER j = 0;
    REAL absb = 0.0;
    COMPLEX signbc = 0.0;
    INTEGER ilast = 0;
    INTEGER ifrstm = 0;
    INTEGER ilastm = 0;
    INTEGER iiter = 0;
    COMPLEX eshift = 0.0;
    INTEGER maxit = 0;
    INTEGER jiter = 0;
    bool ilazro = false;
    REAL temp = 0.0;
    bool ilazr2 = false;
    INTEGER jch = 0;
    COMPLEX ctemp = 0.0;
    REAL c = 0.0;
    COMPLEX s = 0.0;
    INTEGER ifirst = 0;
    COMPLEX u12 = 0.0;
    COMPLEX ad11 = 0.0;
    COMPLEX ad21 = 0.0;
    COMPLEX ad12 = 0.0;
    COMPLEX ad22 = 0.0;
    COMPLEX abi22 = 0.0;
    COMPLEX abi12 = 0.0;
    COMPLEX shift = 0.0;
    const REAL zero = 0.0;
    const REAL half = 0.5e+0;
    REAL temp2 = 0.0;
    COMPLEX y = 0.0;
    INTEGER istart = 0;
    REAL tempr = 0.0;
    COMPLEX ctemp2 = 0.0;
    COMPLEX ctemp3 = 0.0;
    INTEGER jc = 0;
    INTEGER jr = 0;
    //
    //     Decode JOB, COMPQ, COMPZ
    //
    if (Mlsame(job, "E")) {
        ilschr = false;
        ischur = 1;
    } else if (Mlsame(job, "S")) {
        ilschr = true;
        ischur = 2;
    } else {
        ilschr = true;
        ischur = 0;
    }
    //
    if (Mlsame(compq, "N")) {
        ilq = false;
        icompq = 1;
    } else if (Mlsame(compq, "V")) {
        ilq = true;
        icompq = 2;
    } else if (Mlsame(compq, "I")) {
        ilq = true;
        icompq = 3;
    } else {
        ilq = true;
        icompq = 0;
    }
    //
    if (Mlsame(compz, "N")) {
        ilz = false;
        icompz = 1;
    } else if (Mlsame(compz, "V")) {
        ilz = true;
        icompz = 2;
    } else if (Mlsame(compz, "I")) {
        ilz = true;
        icompz = 3;
    } else {
        ilz = true;
        icompz = 0;
    }
    //
    //     Check Argument Values
    //
    info = 0;
    work[1 - 1] = max((INTEGER)1, n);
    lquery = (lwork == -1);
    if (ischur == 0) {
        info = -1;
    } else if (icompq == 0) {
        info = -2;
    } else if (icompz == 0) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (ilo < 1) {
        info = -5;
    } else if (ihi > n || ihi < ilo - 1) {
        info = -6;
    } else if (ldh < n) {
        info = -8;
    } else if (ldt < n) {
        info = -10;
    } else if (ldq < 1 || (ilq && ldq < n)) {
        info = -14;
    } else if (ldz < 1 || (ilz && ldz < n)) {
        info = -16;
    } else if (lwork < max((INTEGER)1, n) && !lquery) {
        info = -18;
    }
    if (info != 0) {
        Mxerbla("Chgeqz", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    //     WORK( 1 ) = CMPLX( 1 )
    if (n <= 0) {
        work[1 - 1] = COMPLEX(1);
        return;
    }
    //
    //     Initialize Q and Z
    //
    if (icompq == 3) {
        Claset("Full", n, n, czero, cone, q, ldq);
    }
    if (icompz == 3) {
        Claset("Full", n, n, czero, cone, z, ldz);
    }
    //
    //     Machine Constants
    //
    in = ihi + 1 - ilo;
    safmin = Rlamch("S");
    ulp = Rlamch("E") * Rlamch("B");
    anorm = Clanhs("F", in, &h[(ilo - 1) + (ilo - 1) * ldh], ldh, rwork);
    bnorm = Clanhs("F", in, &t[(ilo - 1) + (ilo - 1) * ldt], ldt, rwork);
    atol = max(safmin, REAL(ulp * anorm));
    btol = max(safmin, REAL(ulp * bnorm));
    ascale = one / max(safmin, anorm);
    bscale = one / max(safmin, bnorm);
    //
    //     Set Eigenvalues IHI+1:N
    //
    for (j = ihi + 1; j <= n; j = j + 1) {
        absb = abs(t[(j - 1) + (j - 1) * ldt]);
        if (absb > safmin) {
            signbc = conj(t[(j - 1) + (j - 1) * ldt] / absb);
            t[(j - 1) + (j - 1) * ldt] = absb;
            if (ilschr) {
                Cscal(j - 1, signbc, &t[(j - 1) * ldt], 1);
                Cscal(j, signbc, &h[(j - 1) * ldh], 1);
            } else {
                Cscal(1, signbc, &h[(j - 1) + (j - 1) * ldh], 1);
            }
            if (ilz) {
                Cscal(n, signbc, &z[(j - 1) * ldz], 1);
            }
        } else {
            t[(j - 1) + (j - 1) * ldt] = czero;
        }
        alpha[j - 1] = h[(j - 1) + (j - 1) * ldh];
        beta[j - 1] = t[(j - 1) + (j - 1) * ldt];
    }
    //
    //     If IHI < ILO, skip QZ steps
    //
    if (ihi < ilo) {
        goto statement_190;
    }
    //
    //     MAIN QZ ITERATION LOOP
    //
    //     Initialize dynamic indices
    //
    //     Eigenvalues ILAST+1:N have been found.
    //        Column operations modify rows IFRSTM:whatever
    //        Row operations modify columns whatever:ILASTM
    //
    //     If only eigenvalues are being computed, then
    //        IFRSTM is the row of the last splitting row above row ILAST;
    //        this is always at least ILO.
    //     IITER counts iterations since the last eigenvalue was found,
    //        to tell when to use an extraordinary shift.
    //     MAXIT is the maximum number of QZ sweeps allowed.
    //
    ilast = ihi;
    if (ilschr) {
        ifrstm = 1;
        ilastm = n;
    } else {
        ifrstm = ilo;
        ilastm = ihi;
    }
    iiter = 0;
    eshift = czero;
    maxit = 30 * (ihi - ilo + 1);
    //
    for (jiter = 1; jiter <= maxit; jiter = jiter + 1) {
        //
        //        Check for too many iterations.
        //
        if (jiter > maxit) {
            goto statement_180;
        }
        //
        //        Split the matrix if possible.
        //
        //        Two tests:
        //           1: H(j,j-1)=0  or  j=ILO
        //           2: T(j,j)=0
        //
        //        Special case: j=ILAST
        //
        if (ilast == ilo) {
            goto statement_60;
        } else {
            if (abs1(h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) <= max(safmin, REAL(ulp * (abs1(h[(ilast - 1) + (ilast - 1) * ldh]) + abs1(h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh]))))) {
                h[(ilast - 1) + ((ilast - 1) - 1) * ldh] = czero;
                goto statement_60;
            }
        }
        //
        if (abs(t[(ilast - 1) + (ilast - 1) * ldt]) <= max(safmin, REAL(ulp * (abs(t[((ilast - 1) - 1) + (ilast - 1) * ldt]) + abs(t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt]))))) {
            t[(ilast - 1) + (ilast - 1) * ldt] = czero;
            goto statement_50;
        }
        //
        //        General case: j<ILAST
        //
        for (j = ilast - 1; j >= ilo; j = j - 1) {
            //
            //           Test 1: for H(j,j-1)=0 or j=ILO
            //
            if (j == ilo) {
                ilazro = true;
            } else {
                if (abs1(h[(j - 1) + ((j - 1) - 1) * ldh]) <= max(safmin, REAL(ulp * (abs1(h[(j - 1) + (j - 1) * ldh]) + abs1(h[((j - 1) - 1) + ((j - 1) - 1) * ldh]))))) {
                    h[(j - 1) + ((j - 1) - 1) * ldh] = czero;
                    ilazro = true;
                } else {
                    ilazro = false;
                }
            }
            //
            //           Test 2: for T(j,j)=0
            //
            temp = abs(t[(j - 1) + ((j + 1) - 1) * ldt]);
            if (j > ilo) {
                temp += abs(t[((j - 1) - 1) + (j - 1) * ldt]);
            }
            if (abs(t[(j - 1) + (j - 1) * ldt]) < max(safmin, REAL(ulp * temp))) {
                t[(j - 1) + (j - 1) * ldt] = czero;
                //
                //              Test 1a: Check for 2 consecutive small subdiagonals in A
                //
                ilazr2 = false;
                if (!ilazro) {
                    if (abs1(h[(j - 1) + ((j - 1) - 1) * ldh]) * (ascale * abs1(h[((j + 1) - 1) + (j - 1) * ldh])) <= abs1(h[(j - 1) + (j - 1) * ldh]) * (ascale * atol)) {
                        ilazr2 = true;
                    }
                }
                //
                //              If both tests pass (1 & 2), i.e., the leading diagonal
                //              element of B in the block is zero, split a 1x1 block off
                //              at the top. (I.e., at the J-th row/column) The leading
                //              diagonal element of the remainder can also be zero, so
                //              this may have to be done repeatedly.
                //
                if (ilazro || ilazr2) {
                    for (jch = j; jch <= ilast - 1; jch = jch + 1) {
                        ctemp = h[(jch - 1) + (jch - 1) * ldh];
                        Clartg(ctemp, h[((jch + 1) - 1) + (jch - 1) * ldh], c, s, h[(jch - 1) + (jch - 1) * ldh]);
                        h[((jch + 1) - 1) + (jch - 1) * ldh] = czero;
                        Crot(ilastm - jch, &h[(jch - 1) + ((jch + 1) - 1) * ldh], ldh, &h[((jch + 1) - 1) + ((jch + 1) - 1) * ldh], ldh, c, s);
                        Crot(ilastm - jch, &t[(jch - 1) + ((jch + 1) - 1) * ldt], ldt, &t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt], ldt, c, s);
                        if (ilq) {
                            Crot(n, &q[(jch - 1) * ldq], 1, &q[((jch + 1) - 1) * ldq], 1, c, conj(s));
                        }
                        if (ilazr2) {
                            h[(jch - 1) + ((jch - 1) - 1) * ldh] = h[(jch - 1) + ((jch - 1) - 1) * ldh] * c;
                        }
                        ilazr2 = false;
                        if (abs1(t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt]) >= btol) {
                            if (jch + 1 >= ilast) {
                                goto statement_60;
                            } else {
                                ifirst = jch + 1;
                                goto statement_70;
                            }
                        }
                        t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt] = czero;
                    }
                    goto statement_50;
                } else {
                    //
                    //                 Only test 2 passed -- chase the zero to T(ILAST,ILAST)
                    //                 Then process as in the case T(ILAST,ILAST)=0
                    //
                    for (jch = j; jch <= ilast - 1; jch = jch + 1) {
                        ctemp = t[(jch - 1) + ((jch + 1) - 1) * ldt];
                        Clartg(ctemp, t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt], c, s, t[(jch - 1) + ((jch + 1) - 1) * ldt]);
                        t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt] = czero;
                        if (jch < ilastm - 1) {
                            Crot(ilastm - jch - 1, &t[(jch - 1) + ((jch + 2) - 1) * ldt], ldt, &t[((jch + 1) - 1) + ((jch + 2) - 1) * ldt], ldt, c, s);
                        }
                        Crot(ilastm - jch + 2, &h[(jch - 1) + ((jch - 1) - 1) * ldh], ldh, &h[((jch + 1) - 1) + ((jch - 1) - 1) * ldh], ldh, c, s);
                        if (ilq) {
                            Crot(n, &q[(jch - 1) * ldq], 1, &q[((jch + 1) - 1) * ldq], 1, c, conj(s));
                        }
                        ctemp = h[((jch + 1) - 1) + (jch - 1) * ldh];
                        Clartg(ctemp, h[((jch + 1) - 1) + ((jch - 1) - 1) * ldh], c, s, h[((jch + 1) - 1) + (jch - 1) * ldh]);
                        h[((jch + 1) - 1) + ((jch - 1) - 1) * ldh] = czero;
                        Crot(jch + 1 - ifrstm, &h[(ifrstm - 1) + (jch - 1) * ldh], 1, &h[(ifrstm - 1) + ((jch - 1) - 1) * ldh], 1, c, s);
                        Crot(jch - ifrstm, &t[(ifrstm - 1) + (jch - 1) * ldt], 1, &t[(ifrstm - 1) + ((jch - 1) - 1) * ldt], 1, c, s);
                        if (ilz) {
                            Crot(n, &z[(jch - 1) * ldz], 1, &z[((jch - 1) - 1) * ldz], 1, c, s);
                        }
                    }
                    goto statement_50;
                }
            } else if (ilazro) {
                //
                //              Only test 1 passed -- work on J:ILAST
                //
                ifirst = j;
                goto statement_70;
            }
            //
            //           Neither test passed -- try next J
            //
        }
        //
        //        (Drop-through is "impossible")
        //
        info = 2 * n + 1;
        goto statement_210;
    //
    //        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
    //        1x1 block.
    //
    statement_50:
        ctemp = h[(ilast - 1) + (ilast - 1) * ldh];
        Clartg(ctemp, h[(ilast - 1) + ((ilast - 1) - 1) * ldh], c, s, h[(ilast - 1) + (ilast - 1) * ldh]);
        h[(ilast - 1) + ((ilast - 1) - 1) * ldh] = czero;
        Crot(ilast - ifrstm, &h[(ifrstm - 1) + (ilast - 1) * ldh], 1, &h[(ifrstm - 1) + ((ilast - 1) - 1) * ldh], 1, c, s);
        Crot(ilast - ifrstm, &t[(ifrstm - 1) + (ilast - 1) * ldt], 1, &t[(ifrstm - 1) + ((ilast - 1) - 1) * ldt], 1, c, s);
        if (ilz) {
            Crot(n, &z[(ilast - 1) * ldz], 1, &z[((ilast - 1) - 1) * ldz], 1, c, s);
        }
    //
    //        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA
    //
    statement_60:
        absb = abs(t[(ilast - 1) + (ilast - 1) * ldt]);
        if (absb > safmin) {
            signbc = conj(t[(ilast - 1) + (ilast - 1) * ldt] / absb);
            t[(ilast - 1) + (ilast - 1) * ldt] = absb;
            if (ilschr) {
                Cscal(ilast - ifrstm, signbc, &t[(ifrstm - 1) + (ilast - 1) * ldt], 1);
                Cscal(ilast + 1 - ifrstm, signbc, &h[(ifrstm - 1) + (ilast - 1) * ldh], 1);
            } else {
                Cscal(1, signbc, &h[(ilast - 1) + (ilast - 1) * ldh], 1);
            }
            if (ilz) {
                Cscal(n, signbc, &z[(ilast - 1) * ldz], 1);
            }
        } else {
            t[(ilast - 1) + (ilast - 1) * ldt] = czero;
        }
        alpha[ilast - 1] = h[(ilast - 1) + (ilast - 1) * ldh];
        beta[ilast - 1] = t[(ilast - 1) + (ilast - 1) * ldt];
        //
        //        Go to next block -- exit if finished.
        //
        ilast = ilast - 1;
        if (ilast < ilo) {
            goto statement_190;
        }
        //
        //        Reset counters
        //
        iiter = 0;
        eshift = czero;
        if (!ilschr) {
            ilastm = ilast;
            if (ifrstm > ilast) {
                ifrstm = ilo;
            }
        }
        goto statement_160;
    //
    //        QZ step
    //
    //        This iteration only involves rows/columns IFIRST:ILAST.  We
    //        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
    //
    statement_70:
        iiter++;
        if (!ilschr) {
            ifrstm = ifirst;
        }
        //
        //        Compute the Shift.
        //
        //        At this point, IFIRST < ILAST, and the diagonal elements of
        //        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
        //        magnitude)
        //
        if ((iiter / 10) * 10 != iiter) {
            //
            //           The Wilkinson shift (AEP p.512), i.e., the eigenvalue of
            //           the bottom-right 2x2 block of A inv(B) which is nearest to
            //           the bottom-right element.
            //
            //           We factor B as U*D, where U has unit diagonals, and
            //           compute (A*inv(D))*inv(U).
            //
            u12 = (bscale * t[((ilast - 1) - 1) + (ilast - 1) * ldt]) / (bscale * t[(ilast - 1) + (ilast - 1) * ldt]);
            ad11 = (ascale * h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh]) / (bscale * t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt]);
            ad21 = (ascale * h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) / (bscale * t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt]);
            ad12 = (ascale * h[((ilast - 1) - 1) + (ilast - 1) * ldh]) / (bscale * t[(ilast - 1) + (ilast - 1) * ldt]);
            ad22 = (ascale * h[(ilast - 1) + (ilast - 1) * ldh]) / (bscale * t[(ilast - 1) + (ilast - 1) * ldt]);
            abi22 = ad22 - u12 * ad21;
            abi12 = ad12 - u12 * ad11;
            //
            shift = abi22;
            ctemp = sqrt(abi12) * sqrt(ad21);
            temp = abs1(ctemp);
            if (ctemp != zero) {
                x = half * (ad11 - shift);
                temp2 = abs1(x);
                temp = max(temp, abs1(x));
                y = temp * sqrt((x / temp) * (x / temp) + (ctemp / temp) * (ctemp / temp));
                if (temp2 > zero) {
                    if ((x / temp2).real() * y.real() + (x / temp2).imag() * y.imag() < zero) {
                        y = -y;
                    }
                }
                shift = shift - ctemp * Cladiv(ctemp, (x + y));
            }
        } else {
            //
            //           Exceptional shift.  Chosen for no particularly good reason.
            //
            if ((iiter / 20) * 20 == iiter && bscale * abs1(t[(ilast - 1) + (ilast - 1) * ldt]) > safmin) {
                eshift += (ascale * h[(ilast - 1) + (ilast - 1) * ldh]) / (bscale * t[(ilast - 1) + (ilast - 1) * ldt]);
            } else {
                eshift += (ascale * h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) / (bscale * t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt]);
            }
            shift = eshift;
        }
        //
        //        Now check for two consecutive small subdiagonals.
        //
        for (j = ilast - 1; j >= ifirst + 1; j = j - 1) {
            istart = j;
            ctemp = ascale * h[(j - 1) + (j - 1) * ldh] - shift * (bscale * t[(j - 1) + (j - 1) * ldt]);
            temp = abs1(ctemp);
            temp2 = ascale * abs1(h[((j + 1) - 1) + (j - 1) * ldh]);
            tempr = max(temp, temp2);
            if (tempr < one && tempr != zero) {
                temp = temp / tempr;
                temp2 = temp2 / tempr;
            }
            if (abs1(h[(j - 1) + ((j - 1) - 1) * ldh]) * temp2 <= temp * atol) {
                goto statement_90;
            }
        }
        //
        istart = ifirst;
        ctemp = ascale * h[(ifirst - 1) + (ifirst - 1) * ldh] - shift * (bscale * t[(ifirst - 1) + (ifirst - 1) * ldt]);
    statement_90:
        //
        //        Do an implicit-shift QZ sweep.
        //
        //        Initial Q
        //
        ctemp2 = ascale * h[((istart + 1) - 1) + (istart - 1) * ldh];
        Clartg(ctemp, ctemp2, c, s, ctemp3);
        //
        //        Sweep
        //
        for (j = istart; j <= ilast - 1; j = j + 1) {
            if (j > istart) {
                ctemp = h[(j - 1) + ((j - 1) - 1) * ldh];
                Clartg(ctemp, h[((j + 1) - 1) + ((j - 1) - 1) * ldh], c, s, h[(j - 1) + ((j - 1) - 1) * ldh]);
                h[((j + 1) - 1) + ((j - 1) - 1) * ldh] = czero;
            }
            //
            for (jc = j; jc <= ilastm; jc = jc + 1) {
                ctemp = c * h[(j - 1) + (jc - 1) * ldh] + s * h[((j + 1) - 1) + (jc - 1) * ldh];
                h[((j + 1) - 1) + (jc - 1) * ldh] = -conj(s) * h[(j - 1) + (jc - 1) * ldh] + c * h[((j + 1) - 1) + (jc - 1) * ldh];
                h[(j - 1) + (jc - 1) * ldh] = ctemp;
                ctemp2 = c * t[(j - 1) + (jc - 1) * ldt] + s * t[((j + 1) - 1) + (jc - 1) * ldt];
                t[((j + 1) - 1) + (jc - 1) * ldt] = -conj(s) * t[(j - 1) + (jc - 1) * ldt] + c * t[((j + 1) - 1) + (jc - 1) * ldt];
                t[(j - 1) + (jc - 1) * ldt] = ctemp2;
            }
            if (ilq) {
                for (jr = 1; jr <= n; jr = jr + 1) {
                    ctemp = c * q[(jr - 1) + (j - 1) * ldq] + conj(s) * q[(jr - 1) + ((j + 1) - 1) * ldq];
                    q[(jr - 1) + ((j + 1) - 1) * ldq] = -s * q[(jr - 1) + (j - 1) * ldq] + c * q[(jr - 1) + ((j + 1) - 1) * ldq];
                    q[(jr - 1) + (j - 1) * ldq] = ctemp;
                }
            }
            //
            ctemp = t[((j + 1) - 1) + ((j + 1) - 1) * ldt];
            Clartg(ctemp, t[((j + 1) - 1) + (j - 1) * ldt], c, s, t[((j + 1) - 1) + ((j + 1) - 1) * ldt]);
            t[((j + 1) - 1) + (j - 1) * ldt] = czero;
            //
            for (jr = ifrstm; jr <= min(j + 2, ilast); jr = jr + 1) {
                ctemp = c * h[(jr - 1) + ((j + 1) - 1) * ldh] + s * h[(jr - 1) + (j - 1) * ldh];
                h[(jr - 1) + (j - 1) * ldh] = -conj(s) * h[(jr - 1) + ((j + 1) - 1) * ldh] + c * h[(jr - 1) + (j - 1) * ldh];
                h[(jr - 1) + ((j + 1) - 1) * ldh] = ctemp;
            }
            for (jr = ifrstm; jr <= j; jr = jr + 1) {
                ctemp = c * t[(jr - 1) + ((j + 1) - 1) * ldt] + s * t[(jr - 1) + (j - 1) * ldt];
                t[(jr - 1) + (j - 1) * ldt] = -conj(s) * t[(jr - 1) + ((j + 1) - 1) * ldt] + c * t[(jr - 1) + (j - 1) * ldt];
                t[(jr - 1) + ((j + 1) - 1) * ldt] = ctemp;
            }
            if (ilz) {
                for (jr = 1; jr <= n; jr = jr + 1) {
                    ctemp = c * z[(jr - 1) + ((j + 1) - 1) * ldz] + s * z[(jr - 1) + (j - 1) * ldz];
                    z[(jr - 1) + (j - 1) * ldz] = -conj(s) * z[(jr - 1) + ((j + 1) - 1) * ldz] + c * z[(jr - 1) + (j - 1) * ldz];
                    z[(jr - 1) + ((j + 1) - 1) * ldz] = ctemp;
                }
            }
        }
    //
    statement_160:;
        //
    }
//
//     Drop-through = non-convergence
//
statement_180:
    info = ilast;
    goto statement_210;
//
//     Successful completion of all QZ steps
//
statement_190:
    //
    //     Set Eigenvalues 1:ILO-1
    //
    for (j = 1; j <= ilo - 1; j = j + 1) {
        absb = abs(t[(j - 1) + (j - 1) * ldt]);
        if (absb > safmin) {
            signbc = conj(t[(j - 1) + (j - 1) * ldt] / absb);
            t[(j - 1) + (j - 1) * ldt] = absb;
            if (ilschr) {
                Cscal(j - 1, signbc, &t[(j - 1) * ldt], 1);
                Cscal(j, signbc, &h[(j - 1) * ldh], 1);
            } else {
                Cscal(1, signbc, &h[(j - 1) + (j - 1) * ldh], 1);
            }
            if (ilz) {
                Cscal(n, signbc, &z[(j - 1) * ldz], 1);
            }
        } else {
            t[(j - 1) + (j - 1) * ldt] = czero;
        }
        alpha[j - 1] = h[(j - 1) + (j - 1) * ldh];
        beta[j - 1] = t[(j - 1) + (j - 1) * ldt];
    }
    //
    //     Normal Termination
    //
    info = 0;
//
//     Exit (other than argument error) -- return optimal workspace size
//
statement_210:
    work[1 - 1] = COMPLEX(n);
    //
    //     End of Chgeqz
    //
}
