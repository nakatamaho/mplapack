/*
 * Copyright (c) 2021
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

void Rhgeqz(const char *job, const char *compq, const char *compz, INTEGER const &n, INTEGER const &ilo, INTEGER const &ihi, REAL *h, INTEGER const &ldh, REAL *t, INTEGER const &ldt, REAL *alphar, REAL *alphai, REAL *beta, REAL *q, INTEGER const &ldq, REAL *z, INTEGER const &ldz, REAL *work, INTEGER const &lwork, INTEGER &info) {
    bool ilschr = false;
    INTEGER ischur = 0;
    bool ilq = false;
    INTEGER icompq = 0;
    bool ilz = false;
    INTEGER icompz = 0;
    bool lquery = false;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    INTEGER in = 0;
    REAL safmin = 0.0;
    REAL safmax = 0.0;
    REAL ulp = 0.0;
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    REAL atol = 0.0;
    REAL btol = 0.0;
    REAL ascale = 0.0;
    REAL bscale = 0.0;
    INTEGER j = 0;
    INTEGER jr = 0;
    INTEGER ilast = 0;
    INTEGER ifrstm = 0;
    INTEGER ilastm = 0;
    INTEGER iiter = 0;
    REAL eshift = 0.0;
    INTEGER maxit = 0;
    INTEGER jiter = 0;
    bool ilazro = false;
    REAL temp = 0.0;
    bool ilazr2 = false;
    REAL temp2 = 0.0;
    REAL tempr = 0.0;
    INTEGER jch = 0;
    REAL c = 0.0;
    REAL s = 0.0;
    INTEGER ifirst = 0;
    REAL s1 = 0.0;
    REAL wr = 0.0;
    const REAL safety = 1.0e+2;
    REAL s2 = 0.0;
    REAL wr2 = 0.0;
    REAL wi = 0.0;
    const REAL half = 0.5e+0;
    REAL scale = 0.0;
    INTEGER istart = 0;
    INTEGER jc = 0;
    REAL b22 = 0.0;
    REAL b11 = 0.0;
    REAL sr = 0.0;
    REAL cr = 0.0;
    REAL sl = 0.0;
    REAL cl = 0.0;
    REAL s1inv = 0.0;
    REAL a11 = 0.0;
    REAL a21 = 0.0;
    REAL a12 = 0.0;
    REAL a22 = 0.0;
    REAL c11r = 0.0;
    REAL c11i = 0.0;
    REAL c12 = 0.0;
    REAL c21 = 0.0;
    REAL c22r = 0.0;
    REAL c22i = 0.0;
    REAL t1 = 0.0;
    REAL cz = 0.0;
    REAL szr = 0.0;
    REAL szi = 0.0;
    REAL tempi = 0.0;
    REAL an = 0.0;
    REAL bn = 0.0;
    REAL wabs = 0.0;
    REAL cq = 0.0;
    REAL sqr = 0.0;
    REAL sqi = 0.0;
    REAL a1r = 0.0;
    REAL a1i = 0.0;
    REAL a2r = 0.0;
    REAL a2i = 0.0;
    REAL b1r = 0.0;
    REAL b1i = 0.0;
    REAL b1a = 0.0;
    REAL b2r = 0.0;
    REAL b2i = 0.0;
    REAL b2a = 0.0;
    REAL ad11 = 0.0;
    REAL ad21 = 0.0;
    REAL ad12 = 0.0;
    REAL ad22 = 0.0;
    REAL u12 = 0.0;
    REAL ad11l = 0.0;
    REAL ad21l = 0.0;
    REAL ad12l = 0.0;
    REAL ad22l = 0.0;
    REAL ad32l = 0.0;
    REAL u12l = 0.0;
    arr_1d<3, REAL> v(fill0);
    REAL tau = 0.0;
    bool ilpivt = false;
    REAL u1 = 0.0;
    REAL u2 = 0.0;
    REAL w11 = 0.0;
    REAL w21 = 0.0;
    REAL w12 = 0.0;
    REAL w22 = 0.0;
    REAL vs = 0.0;
    //
    //  -- LAPACK computational routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //    $                     SAFETY = 1.0E+0 )
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
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
        info = -15;
    } else if (ldz < 1 || (ilz && ldz < n)) {
        info = -17;
    } else if (lwork < max((INTEGER)1, n) && !lquery) {
        info = -19;
    }
    if (info != 0) {
        Mxerbla("Rhgeqz", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        work[1 - 1] = 1.real();
        return;
    }
    //
    //     Initialize Q and Z
    //
    if (icompq == 3) {
        Rlaset("Full", n, n, zero, one, q, ldq);
    }
    if (icompz == 3) {
        Rlaset("Full", n, n, zero, one, z, ldz);
    }
    //
    //     Machine Constants
    //
    in = ihi + 1 - ilo;
    safmin = dlamch("S");
    safmax = one / safmin;
    ulp = dlamch("E") * dlamch("B");
    anorm = Rlanhs[("F" - 1) + (in - 1) * ldRlanhs];
    bnorm = Rlanhs[("F" - 1) + (in - 1) * ldRlanhs];
    atol = max(safmin, ulp * anorm);
    btol = max(safmin, ulp * bnorm);
    ascale = one / max(safmin, anorm);
    bscale = one / max(safmin, bnorm);
    //
    //     Set Eigenvalues IHI+1:N
    //
    for (j = ihi + 1; j <= n; j = j + 1) {
        if (t[(j - 1) + (j - 1) * ldt] < zero) {
            if (ilschr) {
                for (jr = 1; jr <= j; jr = jr + 1) {
                    h[(jr - 1) + (j - 1) * ldh] = -h[(jr - 1) + (j - 1) * ldh];
                    t[(jr - 1) + (j - 1) * ldt] = -t[(jr - 1) + (j - 1) * ldt];
                }
            } else {
                h[(j - 1) + (j - 1) * ldh] = -h[(j - 1) + (j - 1) * ldh];
                t[(j - 1) + (j - 1) * ldt] = -t[(j - 1) + (j - 1) * ldt];
            }
            if (ilz) {
                for (jr = 1; jr <= n; jr = jr + 1) {
                    z[(jr - 1) + (j - 1) * ldz] = -z[(jr - 1) + (j - 1) * ldz];
                }
            }
        }
        alphar[j - 1] = h[(j - 1) + (j - 1) * ldh];
        alphai[j - 1] = zero;
        beta[j - 1] = t[(j - 1) + (j - 1) * ldt];
    }
    //
    //     If IHI < ILO, skip QZ steps
    //
    if (ihi < ilo) {
        goto statement_380;
    }
    //
    //     MAIN QZ ITERATION LOOP
    //
    //     Initialize dynamic indices
    //
    //     Eigenvalues ILAST+1:N have been found.
    //        Column operations modify rows IFRSTM:whatever.
    //        Row operations modify columns whatever:ILASTM.
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
    eshift = zero;
    maxit = 30 * (ihi - ilo + 1);
    //
    for (jiter = 1; jiter <= maxit; jiter = jiter + 1) {
        //
        //        Split the matrix if possible.
        //
        //        Two tests:
        //           1: H(j,j-1)=0  or  j=ILO
        //           2: T(j,j)=0
        //
        if (ilast == ilo) {
            //
            //           Special case: j=ILAST
            //
            goto statement_80;
        } else {
            if (abs(h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) <= max(safmin, ulp * (abs(h[(ilast - 1) + (ilast - 1) * ldh]) + abs(h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh])))) {
                h[(ilast - 1) + ((ilast - 1) - 1) * ldh] = zero;
                goto statement_80;
            }
        }
        //
        if (abs(t[(ilast - 1) + (ilast - 1) * ldt]) <= max(safmin, ulp * (abs(t[((ilast - 1) - 1) + (ilast - 1) * ldt]) + abs(t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt])))) {
            t[(ilast - 1) + (ilast - 1) * ldt] = zero;
            goto statement_70;
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
                if (abs(h[(j - 1) + ((j - 1) - 1) * ldh]) <= max(safmin, ulp * (abs(h[(j - 1) + (j - 1) * ldh]) + abs(h[((j - 1) - 1) + ((j - 1) - 1) * ldh])))) {
                    h[(j - 1) + ((j - 1) - 1) * ldh] = zero;
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
            if (abs(t[(j - 1) + (j - 1) * ldt]) < max(safmin, ulp * temp)) {
                t[(j - 1) + (j - 1) * ldt] = zero;
                //
                //              Test 1a: Check for 2 consecutive small subdiagonals in A
                //
                ilazr2 = false;
                if (!ilazro) {
                    temp = abs(h[(j - 1) + ((j - 1) - 1) * ldh]);
                    temp2 = abs(h[(j - 1) + (j - 1) * ldh]);
                    tempr = max(temp, temp2);
                    if (tempr < one && tempr != zero) {
                        temp = temp / tempr;
                        temp2 = temp2 / tempr;
                    }
                    if (temp * (ascale * abs(h[((j + 1) - 1) + (j - 1) * ldh])) <= temp2 * (ascale * atol)) {
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
                        temp = h[(jch - 1) + (jch - 1) * ldh];
                        Rlartg(temp, h[((jch + 1) - 1) + (jch - 1) * ldh], c, s, h[(jch - 1) + (jch - 1) * ldh]);
                        h[((jch + 1) - 1) + (jch - 1) * ldh] = zero;
                        Rrot(ilastm - jch, h[(jch - 1) + ((jch + 1) - 1) * ldh], ldh, h[((jch + 1) - 1) + ((jch + 1) - 1) * ldh], ldh, c, s);
                        Rrot(ilastm - jch, t[(jch - 1) + ((jch + 1) - 1) * ldt], ldt, t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt], ldt, c, s);
                        if (ilq) {
                            Rrot(n, q[(jch - 1) * ldq], 1, q[((jch + 1) - 1) * ldq], 1, c, s);
                        }
                        if (ilazr2) {
                            h[(jch - 1) + ((jch - 1) - 1) * ldh] = h[(jch - 1) + ((jch - 1) - 1) * ldh] * c;
                        }
                        ilazr2 = false;
                        if (abs(t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt]) >= btol) {
                            if (jch + 1 >= ilast) {
                                goto statement_80;
                            } else {
                                ifirst = jch + 1;
                                goto statement_110;
                            }
                        }
                        t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt] = zero;
                    }
                    goto statement_70;
                } else {
                    //
                    //                 Only test 2 passed -- chase the zero to T(ILAST,ILAST)
                    //                 Then process as in the case T(ILAST,ILAST)=0
                    //
                    for (jch = j; jch <= ilast - 1; jch = jch + 1) {
                        temp = t[(jch - 1) + ((jch + 1) - 1) * ldt];
                        Rlartg(temp, t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt], c, s, t[(jch - 1) + ((jch + 1) - 1) * ldt]);
                        t[((jch + 1) - 1) + ((jch + 1) - 1) * ldt] = zero;
                        if (jch < ilastm - 1) {
                            Rrot(ilastm - jch - 1, t[(jch - 1) + ((jch + 2) - 1) * ldt], ldt, t[((jch + 1) - 1) + ((jch + 2) - 1) * ldt], ldt, c, s);
                        }
                        Rrot(ilastm - jch + 2, h[(jch - 1) + ((jch - 1) - 1) * ldh], ldh, h[((jch + 1) - 1) + ((jch - 1) - 1) * ldh], ldh, c, s);
                        if (ilq) {
                            Rrot(n, q[(jch - 1) * ldq], 1, q[((jch + 1) - 1) * ldq], 1, c, s);
                        }
                        temp = h[((jch + 1) - 1) + (jch - 1) * ldh];
                        Rlartg(temp, h[((jch + 1) - 1) + ((jch - 1) - 1) * ldh], c, s, h[((jch + 1) - 1) + (jch - 1) * ldh]);
                        h[((jch + 1) - 1) + ((jch - 1) - 1) * ldh] = zero;
                        Rrot(jch + 1 - ifrstm, h[(ifrstm - 1) + (jch - 1) * ldh], 1, h[(ifrstm - 1) + ((jch - 1) - 1) * ldh], 1, c, s);
                        Rrot(jch - ifrstm, t[(ifrstm - 1) + (jch - 1) * ldt], 1, t[(ifrstm - 1) + ((jch - 1) - 1) * ldt], 1, c, s);
                        if (ilz) {
                            Rrot(n, z[(jch - 1) * ldz], 1, z[((jch - 1) - 1) * ldz], 1, c, s);
                        }
                    }
                    goto statement_70;
                }
            } else if (ilazro) {
                //
                //              Only test 1 passed -- work on J:ILAST
                //
                ifirst = j;
                goto statement_110;
            }
            //
            //           Neither test passed -- try next J
            //
        }
        //
        //        (Drop-through is "impossible")
        //
        info = n + 1;
        goto statement_420;
    //
    //        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
    //        1x1 block.
    //
    statement_70:
        temp = h[(ilast - 1) + (ilast - 1) * ldh];
        Rlartg(temp, h[(ilast - 1) + ((ilast - 1) - 1) * ldh], c, s, h[(ilast - 1) + (ilast - 1) * ldh]);
        h[(ilast - 1) + ((ilast - 1) - 1) * ldh] = zero;
        Rrot(ilast - ifrstm, h[(ifrstm - 1) + (ilast - 1) * ldh], 1, h[(ifrstm - 1) + ((ilast - 1) - 1) * ldh], 1, c, s);
        Rrot(ilast - ifrstm, t[(ifrstm - 1) + (ilast - 1) * ldt], 1, t[(ifrstm - 1) + ((ilast - 1) - 1) * ldt], 1, c, s);
        if (ilz) {
            Rrot(n, z[(ilast - 1) * ldz], 1, z[((ilast - 1) - 1) * ldz], 1, c, s);
        }
    //
    //        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
    //                              and BETA
    //
    statement_80:
        if (t[(ilast - 1) + (ilast - 1) * ldt] < zero) {
            if (ilschr) {
                for (j = ifrstm; j <= ilast; j = j + 1) {
                    h[(j - 1) + (ilast - 1) * ldh] = -h[(j - 1) + (ilast - 1) * ldh];
                    t[(j - 1) + (ilast - 1) * ldt] = -t[(j - 1) + (ilast - 1) * ldt];
                }
            } else {
                h[(ilast - 1) + (ilast - 1) * ldh] = -h[(ilast - 1) + (ilast - 1) * ldh];
                t[(ilast - 1) + (ilast - 1) * ldt] = -t[(ilast - 1) + (ilast - 1) * ldt];
            }
            if (ilz) {
                for (j = 1; j <= n; j = j + 1) {
                    z[(j - 1) + (ilast - 1) * ldz] = -z[(j - 1) + (ilast - 1) * ldz];
                }
            }
        }
        alphar[ilast - 1] = h[(ilast - 1) + (ilast - 1) * ldh];
        alphai[ilast - 1] = zero;
        beta[ilast - 1] = t[(ilast - 1) + (ilast - 1) * ldt];
        //
        //        Go to next block -- exit if finished.
        //
        ilast = ilast - 1;
        if (ilast < ilo) {
            goto statement_380;
        }
        //
        //        Reset counters
        //
        iiter = 0;
        eshift = zero;
        if (!ilschr) {
            ilastm = ilast;
            if (ifrstm > ilast) {
                ifrstm = ilo;
            }
        }
        goto statement_350;
    //
    //        QZ step
    //
    //        This iteration only involves rows/columns IFIRST:ILAST. We
    //        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
    //
    statement_110:
        iiter++;
        if (!ilschr) {
            ifrstm = ifirst;
        }
        //
        //        Compute single shifts.
        //
        //        At this poINTEGER, IFIRST < ILAST, and the diagonal elements of
        //        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
        //        magnitude)
        //
        if ((iiter / 10) * 10 == iiter) {
            //
            //           Exceptional shift.  Chosen for no particularly good reason.
            //           (Single shift only.)
            //
            if ((maxit.real() * safmin) * abs(h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) < abs(t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt])) {
                eshift = h[(ilast - 1) + ((ilast - 1) - 1) * ldh] / t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt];
            } else {
                eshift += one / (safmin * maxit.real());
            }
            s1 = one;
            wr = eshift;
            //
        } else {
            //
            //           Shifts based on the generalized eigenvalues of the
            //           bottom-right 2x2 block of A and B. The first eigenvalue
            //           returned by Rlag2 is the Wilkinson shift (AEP p.512),
            //
            Rlag2(h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh], ldh, t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt], ldt, safmin * safety, s1, s2, wr, wr2, wi);
            //
            if (abs((wr / s1) * t[(ilast - 1) + (ilast - 1) * ldt] - h[(ilast - 1) + (ilast - 1) * ldh]) > abs((wr2 / s2) * t[(ilast - 1) + (ilast - 1) * ldt] - h[(ilast - 1) + (ilast - 1) * ldh])) {
                temp = wr;
                wr = wr2;
                wr2 = temp;
                temp = s1;
                s1 = s2;
                s2 = temp;
            }
            temp = max(s1, safmin * max(one, abs(wr), abs(wi)));
            if (wi != zero) {
                goto statement_200;
            }
        }
        //
        //        Fiddle with shift to avoid overflow
        //
        temp = min(ascale, one) * (half * safmax);
        if (s1 > temp) {
            scale = temp / s1;
        } else {
            scale = one;
        }
        //
        temp = min(bscale, one) * (half * safmax);
        if (abs(wr) > temp) {
            scale = min(scale, temp / abs(wr));
        }
        s1 = scale * s1;
        wr = scale * wr;
        //
        //        Now check for two consecutive small subdiagonals.
        //
        for (j = ilast - 1; j >= ifirst + 1; j = j - 1) {
            istart = j;
            temp = abs(s1 * h[(j - 1) + ((j - 1) - 1) * ldh]);
            temp2 = abs(s1 * h[(j - 1) + (j - 1) * ldh] - wr * t[(j - 1) + (j - 1) * ldt]);
            tempr = max(temp, temp2);
            if (tempr < one && tempr != zero) {
                temp = temp / tempr;
                temp2 = temp2 / tempr;
            }
            if (abs((ascale * h[((j + 1) - 1) + (j - 1) * ldh]) * temp) <= (ascale * atol) * temp2) {
                goto statement_130;
            }
        }
        //
        istart = ifirst;
    statement_130:
        //
        //        Do an implicit single-shift QZ sweep.
        //
        //        Initial Q
        //
        temp = s1 * h[(istart - 1) + (istart - 1) * ldh] - wr * t[(istart - 1) + (istart - 1) * ldt];
        temp2 = s1 * h[((istart + 1) - 1) + (istart - 1) * ldh];
        Rlartg(temp, temp2, c, s, tempr);
        //
        //        Sweep
        //
        for (j = istart; j <= ilast - 1; j = j + 1) {
            if (j > istart) {
                temp = h[(j - 1) + ((j - 1) - 1) * ldh];
                Rlartg(temp, h[((j + 1) - 1) + ((j - 1) - 1) * ldh], c, s, h[(j - 1) + ((j - 1) - 1) * ldh]);
                h[((j + 1) - 1) + ((j - 1) - 1) * ldh] = zero;
            }
            //
            for (jc = j; jc <= ilastm; jc = jc + 1) {
                temp = c * h[(j - 1) + (jc - 1) * ldh] + s * h[((j + 1) - 1) + (jc - 1) * ldh];
                h[((j + 1) - 1) + (jc - 1) * ldh] = -s * h[(j - 1) + (jc - 1) * ldh] + c * h[((j + 1) - 1) + (jc - 1) * ldh];
                h[(j - 1) + (jc - 1) * ldh] = temp;
                temp2 = c * t[(j - 1) + (jc - 1) * ldt] + s * t[((j + 1) - 1) + (jc - 1) * ldt];
                t[((j + 1) - 1) + (jc - 1) * ldt] = -s * t[(j - 1) + (jc - 1) * ldt] + c * t[((j + 1) - 1) + (jc - 1) * ldt];
                t[(j - 1) + (jc - 1) * ldt] = temp2;
            }
            if (ilq) {
                for (jr = 1; jr <= n; jr = jr + 1) {
                    temp = c * q[(jr - 1) + (j - 1) * ldq] + s * q[(jr - 1) + ((j + 1) - 1) * ldq];
                    q[(jr - 1) + ((j + 1) - 1) * ldq] = -s * q[(jr - 1) + (j - 1) * ldq] + c * q[(jr - 1) + ((j + 1) - 1) * ldq];
                    q[(jr - 1) + (j - 1) * ldq] = temp;
                }
            }
            //
            temp = t[((j + 1) - 1) + ((j + 1) - 1) * ldt];
            Rlartg(temp, t[((j + 1) - 1) + (j - 1) * ldt], c, s, t[((j + 1) - 1) + ((j + 1) - 1) * ldt]);
            t[((j + 1) - 1) + (j - 1) * ldt] = zero;
            //
            for (jr = ifrstm; jr <= min(j + 2, ilast); jr = jr + 1) {
                temp = c * h[(jr - 1) + ((j + 1) - 1) * ldh] + s * h[(jr - 1) + (j - 1) * ldh];
                h[(jr - 1) + (j - 1) * ldh] = -s * h[(jr - 1) + ((j + 1) - 1) * ldh] + c * h[(jr - 1) + (j - 1) * ldh];
                h[(jr - 1) + ((j + 1) - 1) * ldh] = temp;
            }
            for (jr = ifrstm; jr <= j; jr = jr + 1) {
                temp = c * t[(jr - 1) + ((j + 1) - 1) * ldt] + s * t[(jr - 1) + (j - 1) * ldt];
                t[(jr - 1) + (j - 1) * ldt] = -s * t[(jr - 1) + ((j + 1) - 1) * ldt] + c * t[(jr - 1) + (j - 1) * ldt];
                t[(jr - 1) + ((j + 1) - 1) * ldt] = temp;
            }
            if (ilz) {
                for (jr = 1; jr <= n; jr = jr + 1) {
                    temp = c * z[(jr - 1) + ((j + 1) - 1) * ldz] + s * z[(jr - 1) + (j - 1) * ldz];
                    z[(jr - 1) + (j - 1) * ldz] = -s * z[(jr - 1) + ((j + 1) - 1) * ldz] + c * z[(jr - 1) + (j - 1) * ldz];
                    z[(jr - 1) + ((j + 1) - 1) * ldz] = temp;
                }
            }
        }
        //
        goto statement_350;
    //
    //        Use Francis REAL-shift
    //
    //        Note: the Francis REAL-shift should work with real shifts,
    //              but only if the block is at least 3x3.
    //              This code may break if this poINTEGER is reached with
    //              a 2x2 block with real eigenvalues.
    //
    statement_200:
        if (ifirst + 1 == ilast) {
            //
            //           Special case -- 2x2 block with complex eigenvectors
            //
            //           Step 1: Standardize, that is, rotate so that
            //
            //                       ( B11  0  )
            //                   B = (         )  with B11 non-negative.
            //                       (  0  B22 )
            //
            Rlasv2(t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt], t[((ilast - 1) - 1) + (ilast - 1) * ldt], t[(ilast - 1) + (ilast - 1) * ldt], b22, b11, sr, cr, sl, cl);
            //
            if (b11 < zero) {
                cr = -cr;
                sr = -sr;
                b11 = -b11;
                b22 = -b22;
            }
            //
            Rrot(ilastm + 1 - ifirst, h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh], ldh, h[(ilast - 1) + ((ilast - 1) - 1) * ldh], ldh, cl, sl);
            Rrot(ilast + 1 - ifrstm, h[(ifrstm - 1) + ((ilast - 1) - 1) * ldh], 1, h[(ifrstm - 1) + (ilast - 1) * ldh], 1, cr, sr);
            //
            if (ilast < ilastm) {
                Rrot(ilastm - ilast, t[((ilast - 1) - 1) + ((ilast + 1) - 1) * ldt], ldt, t[(ilast - 1) + ((ilast + 1) - 1) * ldt], ldt, cl, sl);
            }
            if (ifrstm < ilast - 1) {
                Rrot(ifirst - ifrstm, t[(ifrstm - 1) + ((ilast - 1) - 1) * ldt], 1, t[(ifrstm - 1) + (ilast - 1) * ldt], 1, cr, sr);
            }
            //
            if (ilq) {
                Rrot(n, q[((ilast - 1) - 1) * ldq], 1, q[(ilast - 1) * ldq], 1, cl, sl);
            }
            if (ilz) {
                Rrot(n, z[((ilast - 1) - 1) * ldz], 1, z[(ilast - 1) * ldz], 1, cr, sr);
            }
            //
            t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt] = b11;
            t[((ilast - 1) - 1) + (ilast - 1) * ldt] = zero;
            t[(ilast - 1) + ((ilast - 1) - 1) * ldt] = zero;
            t[(ilast - 1) + (ilast - 1) * ldt] = b22;
            //
            //           If B22 is negative, negate column ILAST
            //
            if (b22 < zero) {
                for (j = ifrstm; j <= ilast; j = j + 1) {
                    h[(j - 1) + (ilast - 1) * ldh] = -h[(j - 1) + (ilast - 1) * ldh];
                    t[(j - 1) + (ilast - 1) * ldt] = -t[(j - 1) + (ilast - 1) * ldt];
                }
                //
                if (ilz) {
                    for (j = 1; j <= n; j = j + 1) {
                        z[(j - 1) + (ilast - 1) * ldz] = -z[(j - 1) + (ilast - 1) * ldz];
                    }
                }
                b22 = -b22;
            }
            //
            //           Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)
            //
            //           Recompute shift
            //
            Rlag2(h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh], ldh, t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt], ldt, safmin * safety, s1, temp, wr, temp2, wi);
            //
            //           If standardization has perturbed the shift onto real line,
            //           do another (real single-shift) QR step.
            //
            if (wi == zero) {
                goto statement_350;
            }
            s1inv = one / s1;
            //
            //           Do EISPACK (QZVAL) computation of alpha and beta
            //
            a11 = h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh];
            a21 = h[(ilast - 1) + ((ilast - 1) - 1) * ldh];
            a12 = h[((ilast - 1) - 1) + (ilast - 1) * ldh];
            a22 = h[(ilast - 1) + (ilast - 1) * ldh];
            //
            //           Compute complex Givens rotation on right
            //           (Assume some element of C = (sA - wB) > unfl )
            //                            __
            //           (sA - wB) ( CZ   -SZ )
            //                     ( SZ    CZ )
            //
            c11r = s1 * a11 - wr * b11;
            c11i = -wi * b11;
            c12 = s1 * a12;
            c21 = s1 * a21;
            c22r = s1 * a22 - wr * b22;
            c22i = -wi * b22;
            //
            if (abs(c11r) + abs(c11i) + abs(c12) > abs(c21) + abs(c22r) + abs(c22i)) {
                t1 = Rlapy3[(c12 - 1) + (c11r - 1) * ldRlapy3];
                cz = c12 / t1;
                szr = -c11r / t1;
                szi = -c11i / t1;
            } else {
                cz = Rlapy2[(c22r - 1) + (c22i - 1) * ldRlapy2];
                if (cz <= safmin) {
                    cz = zero;
                    szr = one;
                    szi = zero;
                } else {
                    tempr = c22r / cz;
                    tempi = c22i / cz;
                    t1 = Rlapy2[(cz - 1) + (c21 - 1) * ldRlapy2];
                    cz = cz / t1;
                    szr = -c21 * tempr / t1;
                    szi = c21 * tempi / t1;
                }
            }
            //
            //           Compute Givens rotation on left
            //
            //           (  CQ   SQ )
            //           (  __      )  A or B
            //           ( -SQ   CQ )
            //
            an = abs(a11) + abs(a12) + abs(a21) + abs(a22);
            bn = abs(b11) + abs(b22);
            wabs = abs(wr) + abs(wi);
            if (s1 * an > wabs * bn) {
                cq = cz * b11;
                sqr = szr * b22;
                sqi = -szi * b22;
            } else {
                a1r = cz * a11 + szr * a12;
                a1i = szi * a12;
                a2r = cz * a21 + szr * a22;
                a2i = szi * a22;
                cq = Rlapy2[(a1r - 1) + (a1i - 1) * ldRlapy2];
                if (cq <= safmin) {
                    cq = zero;
                    sqr = one;
                    sqi = zero;
                } else {
                    tempr = a1r / cq;
                    tempi = a1i / cq;
                    sqr = tempr * a2r + tempi * a2i;
                    sqi = tempi * a2r - tempr * a2i;
                }
            }
            t1 = Rlapy3[(cq - 1) + (sqr - 1) * ldRlapy3];
            cq = cq / t1;
            sqr = sqr / t1;
            sqi = sqi / t1;
            //
            //           Compute diagonal elements of QBZ
            //
            tempr = sqr * szr - sqi * szi;
            tempi = sqr * szi + sqi * szr;
            b1r = cq * cz * b11 + tempr * b22;
            b1i = tempi * b22;
            b1a = Rlapy2[(b1r - 1) + (b1i - 1) * ldRlapy2];
            b2r = cq * cz * b22 + tempr * b11;
            b2i = -tempi * b11;
            b2a = Rlapy2[(b2r - 1) + (b2i - 1) * ldRlapy2];
            //
            //           Normalize so beta > 0, and Im( alpha1 ) > 0
            //
            beta[(ilast - 1) - 1] = b1a;
            beta[ilast - 1] = b2a;
            alphar[(ilast - 1) - 1] = (wr * b1a) * s1inv;
            alphai[(ilast - 1) - 1] = (wi * b1a) * s1inv;
            alphar[ilast - 1] = (wr * b2a) * s1inv;
            alphai[ilast - 1] = -(wi * b2a) * s1inv;
            //
            //           Step 3: Go to next block -- exit if finished.
            //
            ilast = ifirst - 1;
            if (ilast < ilo) {
                goto statement_380;
            }
            //
            //           Reset counters
            //
            iiter = 0;
            eshift = zero;
            if (!ilschr) {
                ilastm = ilast;
                if (ifrstm > ilast) {
                    ifrstm = ilo;
                }
            }
            goto statement_350;
        } else {
            //
            //           Usual case: 3x3 or larger block, using Francis implicit
            //                       REAL-shift
            //
            //                                    2
            //           Eigenvalue equation is  w  - c w + d = 0,
            //
            //                                         -1 2        -1
            //           so compute 1st column of  (A B  )  - c A B   + d
            //           using the formula in QZIT (from EISPACK)
            //
            //           We assume that the block is at least 3x3
            //
            ad11 = (ascale * h[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldh]) / (bscale * t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt]);
            ad21 = (ascale * h[(ilast - 1) + ((ilast - 1) - 1) * ldh]) / (bscale * t[((ilast - 1) - 1) + ((ilast - 1) - 1) * ldt]);
            ad12 = (ascale * h[((ilast - 1) - 1) + (ilast - 1) * ldh]) / (bscale * t[(ilast - 1) + (ilast - 1) * ldt]);
            ad22 = (ascale * h[(ilast - 1) + (ilast - 1) * ldh]) / (bscale * t[(ilast - 1) + (ilast - 1) * ldt]);
            u12 = t[((ilast - 1) - 1) + (ilast - 1) * ldt] / t[(ilast - 1) + (ilast - 1) * ldt];
            ad11l = (ascale * h[(ifirst - 1) + (ifirst - 1) * ldh]) / (bscale * t[(ifirst - 1) + (ifirst - 1) * ldt]);
            ad21l = (ascale * h[((ifirst + 1) - 1) + (ifirst - 1) * ldh]) / (bscale * t[(ifirst - 1) + (ifirst - 1) * ldt]);
            ad12l = (ascale * h[(ifirst - 1) + ((ifirst + 1) - 1) * ldh]) / (bscale * t[((ifirst + 1) - 1) + ((ifirst + 1) - 1) * ldt]);
            ad22l = (ascale * h[((ifirst + 1) - 1) + ((ifirst + 1) - 1) * ldh]) / (bscale * t[((ifirst + 1) - 1) + ((ifirst + 1) - 1) * ldt]);
            ad32l = (ascale * h[((ifirst + 2) - 1) + ((ifirst + 1) - 1) * ldh]) / (bscale * t[((ifirst + 1) - 1) + ((ifirst + 1) - 1) * ldt]);
            u12l = t[(ifirst - 1) + ((ifirst + 1) - 1) * ldt] / t[((ifirst + 1) - 1) + ((ifirst + 1) - 1) * ldt];
            //
            v[1 - 1] = (ad11 - ad11l) * (ad22 - ad11l) - ad12 * ad21 + ad21 * u12 * ad11l + (ad12l - ad11l * u12l) * ad21l;
            v[2 - 1] = ((ad22l - ad11l) - ad21l * u12l - (ad11 - ad11l) - (ad22 - ad11l) + ad21 * u12) * ad21l;
            v[3 - 1] = ad32l * ad21l;
            //
            istart = ifirst;
            //
            Rlarfg(3, v[1 - 1], v[2 - 1], 1, tau);
            v[1 - 1] = one;
            //
            //           Sweep
            //
            for (j = istart; j <= ilast - 2; j = j + 1) {
                //
                //              All but last elements: use 3x3 Householder transforms.
                //
                //              Zero (j-1)st column of A
                //
                if (j > istart) {
                    v[1 - 1] = h[(j - 1) + ((j - 1) - 1) * ldh];
                    v[2 - 1] = h[((j + 1) - 1) + ((j - 1) - 1) * ldh];
                    v[3 - 1] = h[((j + 2) - 1) + ((j - 1) - 1) * ldh];
                    //
                    Rlarfg(3, h[(j - 1) + ((j - 1) - 1) * ldh], v[2 - 1], 1, tau);
                    v[1 - 1] = one;
                    h[((j + 1) - 1) + ((j - 1) - 1) * ldh] = zero;
                    h[((j + 2) - 1) + ((j - 1) - 1) * ldh] = zero;
                }
                //
                for (jc = j; jc <= ilastm; jc = jc + 1) {
                    temp = tau * (h[(j - 1) + (jc - 1) * ldh] + v[2 - 1] * h[((j + 1) - 1) + (jc - 1) * ldh] + v[3 - 1] * h[((j + 2) - 1) + (jc - 1) * ldh]);
                    h[(j - 1) + (jc - 1) * ldh] = h[(j - 1) + (jc - 1) * ldh] - temp;
                    h[((j + 1) - 1) + (jc - 1) * ldh] = h[((j + 1) - 1) + (jc - 1) * ldh] - temp * v[2 - 1];
                    h[((j + 2) - 1) + (jc - 1) * ldh] = h[((j + 2) - 1) + (jc - 1) * ldh] - temp * v[3 - 1];
                    temp2 = tau * (t[(j - 1) + (jc - 1) * ldt] + v[2 - 1] * t[((j + 1) - 1) + (jc - 1) * ldt] + v[3 - 1] * t[((j + 2) - 1) + (jc - 1) * ldt]);
                    t[(j - 1) + (jc - 1) * ldt] = t[(j - 1) + (jc - 1) * ldt] - temp2;
                    t[((j + 1) - 1) + (jc - 1) * ldt] = t[((j + 1) - 1) + (jc - 1) * ldt] - temp2 * v[2 - 1];
                    t[((j + 2) - 1) + (jc - 1) * ldt] = t[((j + 2) - 1) + (jc - 1) * ldt] - temp2 * v[3 - 1];
                }
                if (ilq) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        temp = tau * (q[(jr - 1) + (j - 1) * ldq] + v[2 - 1] * q[(jr - 1) + ((j + 1) - 1) * ldq] + v[3 - 1] * q[(jr - 1) + ((j + 2) - 1) * ldq]);
                        q[(jr - 1) + (j - 1) * ldq] = q[(jr - 1) + (j - 1) * ldq] - temp;
                        q[(jr - 1) + ((j + 1) - 1) * ldq] = q[(jr - 1) + ((j + 1) - 1) * ldq] - temp * v[2 - 1];
                        q[(jr - 1) + ((j + 2) - 1) * ldq] = q[(jr - 1) + ((j + 2) - 1) * ldq] - temp * v[3 - 1];
                    }
                }
                //
                //              Zero j-th column of B (see DLAGBC for details)
                //
                //              Swap rows to pivot
                //
                ilpivt = false;
                temp = max(abs(t[((j + 1) - 1) + ((j + 1) - 1) * ldt]), abs(t[((j + 1) - 1) + ((j + 2) - 1) * ldt]));
                temp2 = max(abs(t[((j + 2) - 1) + ((j + 1) - 1) * ldt]), abs(t[((j + 2) - 1) + ((j + 2) - 1) * ldt]));
                if (max(temp, temp2) < safmin) {
                    scale = zero;
                    u1 = one;
                    u2 = zero;
                    goto statement_250;
                } else if (temp >= temp2) {
                    w11 = t[((j + 1) - 1) + ((j + 1) - 1) * ldt];
                    w21 = t[((j + 2) - 1) + ((j + 1) - 1) * ldt];
                    w12 = t[((j + 1) - 1) + ((j + 2) - 1) * ldt];
                    w22 = t[((j + 2) - 1) + ((j + 2) - 1) * ldt];
                    u1 = t[((j + 1) - 1) + (j - 1) * ldt];
                    u2 = t[((j + 2) - 1) + (j - 1) * ldt];
                } else {
                    w21 = t[((j + 1) - 1) + ((j + 1) - 1) * ldt];
                    w11 = t[((j + 2) - 1) + ((j + 1) - 1) * ldt];
                    w22 = t[((j + 1) - 1) + ((j + 2) - 1) * ldt];
                    w12 = t[((j + 2) - 1) + ((j + 2) - 1) * ldt];
                    u2 = t[((j + 1) - 1) + (j - 1) * ldt];
                    u1 = t[((j + 2) - 1) + (j - 1) * ldt];
                }
                //
                //              Swap columns if nec.
                //
                if (abs(w12) > abs(w11)) {
                    ilpivt = true;
                    temp = w12;
                    temp2 = w22;
                    w12 = w11;
                    w22 = w21;
                    w11 = temp;
                    w21 = temp2;
                }
                //
                //              LU-factor
                //
                temp = w21 / w11;
                u2 = u2 - temp * u1;
                w22 = w22 - temp * w12;
                w21 = zero;
                //
                //              Compute SCALE
                //
                scale = one;
                if (abs(w22) < safmin) {
                    scale = zero;
                    u2 = one;
                    u1 = -w12 / w11;
                    goto statement_250;
                }
                if (abs(w22) < abs(u2)) {
                    scale = abs(w22 / u2);
                }
                if (abs(w11) < abs(u1)) {
                    scale = min(scale, abs(w11 / u1));
                }
                //
                //              Solve
                //
                u2 = (scale * u2) / w22;
                u1 = (scale * u1 - w12 * u2) / w11;
            //
            statement_250:
                if (ilpivt) {
                    temp = u2;
                    u2 = u1;
                    u1 = temp;
                }
                //
                //              Compute Householder Vector
                //
                t1 = sqrt(pow2(scale) + pow2(u1) + pow2(u2));
                tau = one + scale / t1;
                vs = -one / (scale + t1);
                v[1 - 1] = one;
                v[2 - 1] = vs * u1;
                v[3 - 1] = vs * u2;
                //
                //              Apply transformations from the right.
                //
                for (jr = ifrstm; jr <= min(j + 3, ilast); jr = jr + 1) {
                    temp = tau * (h[(jr - 1) + (j - 1) * ldh] + v[2 - 1] * h[(jr - 1) + ((j + 1) - 1) * ldh] + v[3 - 1] * h[(jr - 1) + ((j + 2) - 1) * ldh]);
                    h[(jr - 1) + (j - 1) * ldh] = h[(jr - 1) + (j - 1) * ldh] - temp;
                    h[(jr - 1) + ((j + 1) - 1) * ldh] = h[(jr - 1) + ((j + 1) - 1) * ldh] - temp * v[2 - 1];
                    h[(jr - 1) + ((j + 2) - 1) * ldh] = h[(jr - 1) + ((j + 2) - 1) * ldh] - temp * v[3 - 1];
                }
                for (jr = ifrstm; jr <= j + 2; jr = jr + 1) {
                    temp = tau * (t[(jr - 1) + (j - 1) * ldt] + v[2 - 1] * t[(jr - 1) + ((j + 1) - 1) * ldt] + v[3 - 1] * t[(jr - 1) + ((j + 2) - 1) * ldt]);
                    t[(jr - 1) + (j - 1) * ldt] = t[(jr - 1) + (j - 1) * ldt] - temp;
                    t[(jr - 1) + ((j + 1) - 1) * ldt] = t[(jr - 1) + ((j + 1) - 1) * ldt] - temp * v[2 - 1];
                    t[(jr - 1) + ((j + 2) - 1) * ldt] = t[(jr - 1) + ((j + 2) - 1) * ldt] - temp * v[3 - 1];
                }
                if (ilz) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        temp = tau * (z[(jr - 1) + (j - 1) * ldz] + v[2 - 1] * z[(jr - 1) + ((j + 1) - 1) * ldz] + v[3 - 1] * z[(jr - 1) + ((j + 2) - 1) * ldz]);
                        z[(jr - 1) + (j - 1) * ldz] = z[(jr - 1) + (j - 1) * ldz] - temp;
                        z[(jr - 1) + ((j + 1) - 1) * ldz] = z[(jr - 1) + ((j + 1) - 1) * ldz] - temp * v[2 - 1];
                        z[(jr - 1) + ((j + 2) - 1) * ldz] = z[(jr - 1) + ((j + 2) - 1) * ldz] - temp * v[3 - 1];
                    }
                }
                t[((j + 1) - 1) + (j - 1) * ldt] = zero;
                t[((j + 2) - 1) + (j - 1) * ldt] = zero;
            }
            //
            //           Last elements: Use Givens rotations
            //
            //           Rotations from the left
            //
            j = ilast - 1;
            temp = h[(j - 1) + ((j - 1) - 1) * ldh];
            Rlartg(temp, h[((j + 1) - 1) + ((j - 1) - 1) * ldh], c, s, h[(j - 1) + ((j - 1) - 1) * ldh]);
            h[((j + 1) - 1) + ((j - 1) - 1) * ldh] = zero;
            //
            for (jc = j; jc <= ilastm; jc = jc + 1) {
                temp = c * h[(j - 1) + (jc - 1) * ldh] + s * h[((j + 1) - 1) + (jc - 1) * ldh];
                h[((j + 1) - 1) + (jc - 1) * ldh] = -s * h[(j - 1) + (jc - 1) * ldh] + c * h[((j + 1) - 1) + (jc - 1) * ldh];
                h[(j - 1) + (jc - 1) * ldh] = temp;
                temp2 = c * t[(j - 1) + (jc - 1) * ldt] + s * t[((j + 1) - 1) + (jc - 1) * ldt];
                t[((j + 1) - 1) + (jc - 1) * ldt] = -s * t[(j - 1) + (jc - 1) * ldt] + c * t[((j + 1) - 1) + (jc - 1) * ldt];
                t[(j - 1) + (jc - 1) * ldt] = temp2;
            }
            if (ilq) {
                for (jr = 1; jr <= n; jr = jr + 1) {
                    temp = c * q[(jr - 1) + (j - 1) * ldq] + s * q[(jr - 1) + ((j + 1) - 1) * ldq];
                    q[(jr - 1) + ((j + 1) - 1) * ldq] = -s * q[(jr - 1) + (j - 1) * ldq] + c * q[(jr - 1) + ((j + 1) - 1) * ldq];
                    q[(jr - 1) + (j - 1) * ldq] = temp;
                }
            }
            //
            //           Rotations from the right.
            //
            temp = t[((j + 1) - 1) + ((j + 1) - 1) * ldt];
            Rlartg(temp, t[((j + 1) - 1) + (j - 1) * ldt], c, s, t[((j + 1) - 1) + ((j + 1) - 1) * ldt]);
            t[((j + 1) - 1) + (j - 1) * ldt] = zero;
            //
            for (jr = ifrstm; jr <= ilast; jr = jr + 1) {
                temp = c * h[(jr - 1) + ((j + 1) - 1) * ldh] + s * h[(jr - 1) + (j - 1) * ldh];
                h[(jr - 1) + (j - 1) * ldh] = -s * h[(jr - 1) + ((j + 1) - 1) * ldh] + c * h[(jr - 1) + (j - 1) * ldh];
                h[(jr - 1) + ((j + 1) - 1) * ldh] = temp;
            }
            for (jr = ifrstm; jr <= ilast - 1; jr = jr + 1) {
                temp = c * t[(jr - 1) + ((j + 1) - 1) * ldt] + s * t[(jr - 1) + (j - 1) * ldt];
                t[(jr - 1) + (j - 1) * ldt] = -s * t[(jr - 1) + ((j + 1) - 1) * ldt] + c * t[(jr - 1) + (j - 1) * ldt];
                t[(jr - 1) + ((j + 1) - 1) * ldt] = temp;
            }
            if (ilz) {
                for (jr = 1; jr <= n; jr = jr + 1) {
                    temp = c * z[(jr - 1) + ((j + 1) - 1) * ldz] + s * z[(jr - 1) + (j - 1) * ldz];
                    z[(jr - 1) + (j - 1) * ldz] = -s * z[(jr - 1) + ((j + 1) - 1) * ldz] + c * z[(jr - 1) + (j - 1) * ldz];
                    z[(jr - 1) + ((j + 1) - 1) * ldz] = temp;
                }
            }
            //
            //           End of Double-Shift code
            //
        }
        //
        goto statement_350;
    //
    //        End of iteration loop
    //
    statement_350:;
    }
    //
    //     Drop-through = non-convergence
    //
    info = ilast;
    goto statement_420;
//
//     Successful completion of all QZ steps
//
statement_380:
    //
    //     Set Eigenvalues 1:ILO-1
    //
    for (j = 1; j <= ilo - 1; j = j + 1) {
        if (t[(j - 1) + (j - 1) * ldt] < zero) {
            if (ilschr) {
                for (jr = 1; jr <= j; jr = jr + 1) {
                    h[(jr - 1) + (j - 1) * ldh] = -h[(jr - 1) + (j - 1) * ldh];
                    t[(jr - 1) + (j - 1) * ldt] = -t[(jr - 1) + (j - 1) * ldt];
                }
            } else {
                h[(j - 1) + (j - 1) * ldh] = -h[(j - 1) + (j - 1) * ldh];
                t[(j - 1) + (j - 1) * ldt] = -t[(j - 1) + (j - 1) * ldt];
            }
            if (ilz) {
                for (jr = 1; jr <= n; jr = jr + 1) {
                    z[(jr - 1) + (j - 1) * ldz] = -z[(jr - 1) + (j - 1) * ldz];
                }
            }
        }
        alphar[j - 1] = h[(j - 1) + (j - 1) * ldh];
        alphai[j - 1] = zero;
        beta[j - 1] = t[(j - 1) + (j - 1) * ldt];
    }
    //
    //     Normal Termination
    //
    info = 0;
//
//     Exit (other than argument error) -- return optimal workspace size
//
statement_420:
    work[1 - 1] = n.real();
    //
    //     End of Rhgeqz
    //
}
