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

void Claqr5(bool const wantt, bool const wantz, INTEGER const kacc22, INTEGER const n, INTEGER const ktop, INTEGER const kbot, INTEGER const nshfts, COMPLEX *s, COMPLEX *h, INTEGER const ldh, INTEGER const iloz, INTEGER const ihiz, COMPLEX *z, INTEGER const ldz, COMPLEX *v, INTEGER const ldv, COMPLEX *u, INTEGER const ldu, INTEGER const nv, COMPLEX *wv, INTEGER const ldwv, INTEGER const nh, COMPLEX *wh, INTEGER const ldwh) {
    COMPLEX cdum = 0.0;
    //
    //
    if (nshfts < 2) {
        return;
    }
    //
    //
    if (ktop >= kbot) {
        return;
    }
    //
    //
    INTEGER ns = nshfts - mod(nshfts, 2);
    //
    //
    REAL safmin = Rlamch("SAFE MINIMUM");
    const REAL rone = 1.0;
    REAL safmax = rone / safmin;
    REAL ulp = Rlamch("PRECISION");
    REAL smlnum = safmin * (castREAL(n) / ulp);
    //
    //
    bool accum = (kacc22 == 1) || (kacc22 == 2);
    //
    //
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    if (ktop + 2 <= kbot) {
        h[((ktop + 2) - 1) + (ktop - 1) * ldh] = zero;
    }
    //
    //
    INTEGER nbmps = ns / 2;
    //
    //
    INTEGER kdu = 4 * nbmps;
    //
    //
    INTEGER incol = 0;
    INTEGER jtop = 0;
    INTEGER ndcol = 0;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    INTEGER krcol = 0;
    INTEGER mtop = 0;
    INTEGER mbot = 0;
    INTEGER m22 = 0;
    bool bmp22 = false;
    INTEGER k = 0;
    COMPLEX beta = 0.0;
    INTEGER j = 0;
    COMPLEX refsum = 0.0;
    INTEGER jbot = 0;
    REAL tst1 = 0.0;
    const REAL rzero = 0.0;
    REAL h12 = 0.0;
    REAL h21 = 0.0;
    REAL h11 = 0.0;
    REAL h22 = 0.0;
    REAL scl = 0.0;
    REAL tst2 = 0.0;
    INTEGER kms = 0;
    INTEGER m = 0;
    COMPLEX alpha = 0.0;
    COMPLEX vt[3];
    INTEGER i2 = 0;
    INTEGER i4 = 0;
    INTEGER k1 = 0;
    INTEGER nu = 0;
    INTEGER jcol = 0;
    INTEGER jlen = 0;
    INTEGER jrow = 0;
    for (incol = ktop - 2 * nbmps + 1; incol <= kbot - 2; incol = incol + 2 * nbmps) {
        //
        // JTOP = Index from which updates from the right start.
        //
        if (accum) {
            jtop = max(ktop, incol);
        } else if (wantt) {
            jtop = 1;
        } else {
            jtop = ktop;
        }
        //
        ndcol = incol + kdu;
        if (accum) {
            Claset("ALL", kdu, kdu, zero, one, u, ldu);
        }
        //
        // .    performs the near-the-diagonal part of a small bulge
        // .    multi-shift QR sweep.  Each 4*NBMPS column diagonal
        // .    chunk extends from column INCOL to column NDCOL
        // .    (including both column INCOL and column NDCOL). The
        // .    following loop chases a 2*NBMPS+1 column long chain of
        // .    NBMPS bulges 2*NBMPS columns to the right.  (INCOL
        // .    may be less than KTOP and and NDCOL may be greater than
        // .    KBOT indicating phantom columns from which to chase
        // .    bulges before they are actually introduced or to which
        //
        for (krcol = incol; krcol <= min(incol + 2 * nbmps - 1, kbot - 2); krcol = krcol + 1) {
            //
            // .    shift bulges.  There may or may not also be small
            // .    2-by-2 bulge, if there is room.  The inactive bulges
            // .    (if any) must wait until the active bulges have moved
            // .    down the diagonal to make room.  The phantom matrix
            //
            mtop = max((INTEGER)1, (ktop - krcol) / 2 + 1);
            mbot = min(nbmps, (kbot - krcol - 1) / 2);
            m22 = mbot + 1;
            bmp22 = (mbot < nbmps) && (krcol + 2 * (m22 - 1)) == (kbot - 2);
            //
            //
            if (bmp22) {
                //
                //
                k = krcol + 2 * (m22 - 1);
                if (k == ktop - 1) {
                    Claqr1(2, &h[((k + 1) - 1) + ((k + 1) - 1) * ldh], ldh, s[(2 * m22 - 1) - 1], s[(2 * m22) - 1], &v[(m22 - 1) * ldv]);
                    beta = v[(m22 - 1) * ldv];
                    Clarfg(2, beta, &v[(2 - 1) + (m22 - 1) * ldv], 1, v[(m22 - 1) * ldv]);
                } else {
                    beta = h[((k + 1) - 1) + (k - 1) * ldh];
                    v[(2 - 1) + (m22 - 1) * ldv] = h[((k + 2) - 1) + (k - 1) * ldh];
                    Clarfg(2, beta, &v[(2 - 1) + (m22 - 1) * ldv], 1, v[(m22 - 1) * ldv]);
                    h[((k + 1) - 1) + (k - 1) * ldh] = beta;
                    h[((k + 2) - 1) + (k - 1) * ldh] = zero;
                }
                //
                //
                for (j = jtop; j <= min(kbot, k + 3); j = j + 1) {
                    refsum = v[(m22 - 1) * ldv] * (h[(j - 1) + ((k + 1) - 1) * ldh] + v[(2 - 1) + (m22 - 1) * ldv] * h[(j - 1) + ((k + 2) - 1) * ldh]);
                    h[(j - 1) + ((k + 1) - 1) * ldh] = h[(j - 1) + ((k + 1) - 1) * ldh] - refsum;
                    h[(j - 1) + ((k + 2) - 1) * ldh] = h[(j - 1) + ((k + 2) - 1) * ldh] - refsum * conj(v[(2 - 1) + (m22 - 1) * ldv]);
                }
                //
                //
                if (accum) {
                    jbot = min(ndcol, kbot);
                } else if (wantt) {
                    jbot = n;
                } else {
                    jbot = kbot;
                }
                for (j = k + 1; j <= jbot; j = j + 1) {
                    refsum = conj(v[(m22 - 1) * ldv]) * (h[((k + 1) - 1) + (j - 1) * ldh] + conj(v[(2 - 1) + (m22 - 1) * ldv]) * h[((k + 2) - 1) + (j - 1) * ldh]);
                    h[((k + 1) - 1) + (j - 1) * ldh] = h[((k + 1) - 1) + (j - 1) * ldh] - refsum;
                    h[((k + 2) - 1) + (j - 1) * ldh] = h[((k + 2) - 1) + (j - 1) * ldh] - refsum * v[(2 - 1) + (m22 - 1) * ldv];
                }
                //
                // .    the tradition small-compared-to-nearby-diagonals
                // .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
                // .    criteria both be satisfied.  The latter improves
                // .    accuracy in some examples. Falling back on an
                // .    alternate convergence criterion when TST1 or TST2
                // .    is zero (as done here) is traditional but probably
                //
                if (k >= ktop) {
                    if (h[((k + 1) - 1) + (k - 1) * ldh] != zero) {
                        tst1 = cabs1(h[(k - 1) + (k - 1) * ldh]) + cabs1(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]);
                        if (tst1 == rzero) {
                            if (k >= ktop + 1) {
                                tst1 += cabs1(h[(k - 1) + ((k - 1) - 1) * ldh]);
                            }
                            if (k >= ktop + 2) {
                                tst1 += cabs1(h[(k - 1) + ((k - 2) - 1) * ldh]);
                            }
                            if (k >= ktop + 3) {
                                tst1 += cabs1(h[(k - 1) + ((k - 3) - 1) * ldh]);
                            }
                            if (k <= kbot - 2) {
                                tst1 += cabs1(h[((k + 2) - 1) + ((k + 1) - 1) * ldh]);
                            }
                            if (k <= kbot - 3) {
                                tst1 += cabs1(h[((k + 3) - 1) + ((k + 1) - 1) * ldh]);
                            }
                            if (k <= kbot - 4) {
                                tst1 += cabs1(h[((k + 4) - 1) + ((k + 1) - 1) * ldh]);
                            }
                        }
                        if (cabs1(h[((k + 1) - 1) + (k - 1) * ldh]) <= max(smlnum, REAL(ulp * tst1))) {
                            h12 = max(cabs1(h[((k + 1) - 1) + (k - 1) * ldh]), cabs1(h[(k - 1) + ((k + 1) - 1) * ldh]));
                            h21 = min(cabs1(h[((k + 1) - 1) + (k - 1) * ldh]), cabs1(h[(k - 1) + ((k + 1) - 1) * ldh]));
                            h11 = max(cabs1(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]), cabs1(h[(k - 1) + (k - 1) * ldh] - h[((k + 1) - 1) + ((k + 1) - 1) * ldh]));
                            h22 = min(cabs1(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]), cabs1(h[(k - 1) + (k - 1) * ldh] - h[((k + 1) - 1) + ((k + 1) - 1) * ldh]));
                            scl = h11 + h12;
                            tst2 = h22 * (h11 / scl);
                            //
                            if (tst2 == rzero || h21 * (h12 / scl) <= max(smlnum, REAL(ulp * tst2))) {
                                h[((k + 1) - 1) + (k - 1) * ldh] = zero;
                            }
                        }
                    }
                }
                //
                //
                if (accum) {
                    kms = k - incol;
                    for (j = max((INTEGER)1, ktop - incol); j <= kdu; j = j + 1) {
                        refsum = v[(m22 - 1) * ldv] * (u[(j - 1) + ((kms + 1) - 1) * ldu] + v[(2 - 1) + (m22 - 1) * ldv] * u[(j - 1) + ((kms + 2) - 1) * ldu]);
                        u[(j - 1) + ((kms + 1) - 1) * ldu] = u[(j - 1) + ((kms + 1) - 1) * ldu] - refsum;
                        u[(j - 1) + ((kms + 2) - 1) * ldu] = u[(j - 1) + ((kms + 2) - 1) * ldu] - refsum * conj(v[(2 - 1) + (m22 - 1) * ldv]);
                    }
                } else if (wantz) {
                    for (j = iloz; j <= ihiz; j = j + 1) {
                        refsum = v[(m22 - 1) * ldv] * (z[(j - 1) + ((k + 1) - 1) * ldz] + v[(2 - 1) + (m22 - 1) * ldv] * z[(j - 1) + ((k + 2) - 1) * ldz]);
                        z[(j - 1) + ((k + 1) - 1) * ldz] = z[(j - 1) + ((k + 1) - 1) * ldz] - refsum;
                        z[(j - 1) + ((k + 2) - 1) * ldz] = z[(j - 1) + ((k + 2) - 1) * ldz] - refsum * conj(v[(2 - 1) + (m22 - 1) * ldv]);
                    }
                }
            }
            //
            //
            for (m = mbot; m >= mtop; m = m - 1) {
                k = krcol + 2 * (m - 1);
                if (k == ktop - 1) {
                    Claqr1(3, &h[(ktop - 1) + (ktop - 1) * ldh], ldh, s[(2 * m - 1) - 1], s[(2 * m) - 1], &v[(m - 1) * ldv]);
                    alpha = v[(m - 1) * ldv];
                    Clarfg(3, alpha, &v[(2 - 1) + (m - 1) * ldv], 1, v[(m - 1) * ldv]);
                } else {
                    //
                    // .    Mth bulge. Exploit fact that first two elements
                    //
                    refsum = v[(m - 1) * ldv] * v[(3 - 1) + (m - 1) * ldv] * h[((k + 3) - 1) + ((k + 2) - 1) * ldh];
                    h[((k + 3) - 1) + (k - 1) * ldh] = -refsum;
                    h[((k + 3) - 1) + ((k + 1) - 1) * ldh] = -refsum * conj(v[(2 - 1) + (m - 1) * ldv]);
                    h[((k + 3) - 1) + ((k + 2) - 1) * ldh] = h[((k + 3) - 1) + ((k + 2) - 1) * ldh] - refsum * conj(v[(3 - 1) + (m - 1) * ldv]);
                    //
                    //
                    beta = h[((k + 1) - 1) + (k - 1) * ldh];
                    v[(2 - 1) + (m - 1) * ldv] = h[((k + 2) - 1) + (k - 1) * ldh];
                    v[(3 - 1) + (m - 1) * ldv] = h[((k + 3) - 1) + (k - 1) * ldh];
                    Clarfg(3, beta, &v[(2 - 1) + (m - 1) * ldv], 1, v[(m - 1) * ldv]);
                    //
                    // .    deflation or destructive underflow.  In the
                    // .    underflow case, try the two-small-subdiagonals
                    //
                    if (h[((k + 3) - 1) + (k - 1) * ldh] != zero || h[((k + 3) - 1) + ((k + 1) - 1) * ldh] != zero || h[((k + 3) - 1) + ((k + 2) - 1) * ldh] == zero) {
                        //
                        //
                        h[((k + 1) - 1) + (k - 1) * ldh] = beta;
                        h[((k + 2) - 1) + (k - 1) * ldh] = zero;
                        h[((k + 3) - 1) + (k - 1) * ldh] = zero;
                    } else {
                        //
                        // .    reintroduce ignoring H(K+1,K) and H(K+2,K).
                        // .    If the fill resulting from the new
                        // .    reflector is too large, then abandon it.
                        //
                        Claqr1(3, &h[((k + 1) - 1) + ((k + 1) - 1) * ldh], ldh, s[(2 * m - 1) - 1], s[(2 * m) - 1], vt);
                        alpha = vt[0];
                        Clarfg(3, alpha, &vt[2 - 1], 1, vt[0]);
                        refsum = conj(vt[0]) * (h[((k + 1) - 1) + (k - 1) * ldh] + conj(vt[2 - 1]) * h[((k + 2) - 1) + (k - 1) * ldh]);
                        //
                        if (cabs1(h[((k + 2) - 1) + (k - 1) * ldh] - refsum * vt[2 - 1]) + cabs1(refsum * vt[3 - 1]) > ulp * (cabs1(h[(k - 1) + (k - 1) * ldh]) + cabs1(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) + cabs1(h[((k + 2) - 1) + ((k + 2) - 1) * ldh]))) {
                            //
                            // .    create non-negligible fill.  Use
                            //
                            h[((k + 1) - 1) + (k - 1) * ldh] = beta;
                            h[((k + 2) - 1) + (k - 1) * ldh] = zero;
                            h[((k + 3) - 1) + (k - 1) * ldh] = zero;
                        } else {
                            //
                            // .    create only negligible fill.
                            // .    Replace the old reflector with
                            //
                            h[((k + 1) - 1) + (k - 1) * ldh] = h[((k + 1) - 1) + (k - 1) * ldh] - refsum;
                            h[((k + 2) - 1) + (k - 1) * ldh] = zero;
                            h[((k + 3) - 1) + (k - 1) * ldh] = zero;
                            v[(m - 1) * ldv] = vt[0];
                            v[(2 - 1) + (m - 1) * ldv] = vt[2 - 1];
                            v[(3 - 1) + (m - 1) * ldv] = vt[3 - 1];
                        }
                    }
                }
                //
                // .     the first column of update from the left.
                // .     These updates are required for the vigilant
                // .     deflation check. We still delay most of the
                //
                for (j = jtop; j <= min(kbot, k + 3); j = j + 1) {
                    refsum = v[(m - 1) * ldv] * (h[(j - 1) + ((k + 1) - 1) * ldh] + v[(2 - 1) + (m - 1) * ldv] * h[(j - 1) + ((k + 2) - 1) * ldh] + v[(3 - 1) + (m - 1) * ldv] * h[(j - 1) + ((k + 3) - 1) * ldh]);
                    h[(j - 1) + ((k + 1) - 1) * ldh] = h[(j - 1) + ((k + 1) - 1) * ldh] - refsum;
                    h[(j - 1) + ((k + 2) - 1) * ldh] = h[(j - 1) + ((k + 2) - 1) * ldh] - refsum * conj(v[(2 - 1) + (m - 1) * ldv]);
                    h[(j - 1) + ((k + 3) - 1) * ldh] = h[(j - 1) + ((k + 3) - 1) * ldh] - refsum * conj(v[(3 - 1) + (m - 1) * ldv]);
                }
                //
                //
                refsum = conj(v[(m - 1) * ldv]) * (h[((k + 1) - 1) + ((k + 1) - 1) * ldh] + conj(v[(2 - 1) + (m - 1) * ldv]) * h[((k + 2) - 1) + ((k + 1) - 1) * ldh] + conj(v[(3 - 1) + (m - 1) * ldv]) * h[((k + 3) - 1) + ((k + 1) - 1) * ldh]);
                h[((k + 1) - 1) + ((k + 1) - 1) * ldh] = h[((k + 1) - 1) + ((k + 1) - 1) * ldh] - refsum;
                h[((k + 2) - 1) + ((k + 1) - 1) * ldh] = h[((k + 2) - 1) + ((k + 1) - 1) * ldh] - refsum * v[(2 - 1) + (m - 1) * ldv];
                h[((k + 3) - 1) + ((k + 1) - 1) * ldh] = h[((k + 3) - 1) + ((k + 1) - 1) * ldh] - refsum * v[(3 - 1) + (m - 1) * ldv];
                //
                // .    the tradition small-compared-to-nearby-diagonals
                // .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
                // .    criteria both be satisfied.  The latter improves
                // .    accuracy in some examples. Falling back on an
                // .    alternate convergence criterion when TST1 or TST2
                // .    is zero (as done here) is traditional but probably
                //
                if (k < ktop) {
                    continue;
                }
                if (h[((k + 1) - 1) + (k - 1) * ldh] != zero) {
                    tst1 = cabs1(h[(k - 1) + (k - 1) * ldh]) + cabs1(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]);
                    if (tst1 == rzero) {
                        if (k >= ktop + 1) {
                            tst1 += cabs1(h[(k - 1) + ((k - 1) - 1) * ldh]);
                        }
                        if (k >= ktop + 2) {
                            tst1 += cabs1(h[(k - 1) + ((k - 2) - 1) * ldh]);
                        }
                        if (k >= ktop + 3) {
                            tst1 += cabs1(h[(k - 1) + ((k - 3) - 1) * ldh]);
                        }
                        if (k <= kbot - 2) {
                            tst1 += cabs1(h[((k + 2) - 1) + ((k + 1) - 1) * ldh]);
                        }
                        if (k <= kbot - 3) {
                            tst1 += cabs1(h[((k + 3) - 1) + ((k + 1) - 1) * ldh]);
                        }
                        if (k <= kbot - 4) {
                            tst1 += cabs1(h[((k + 4) - 1) + ((k + 1) - 1) * ldh]);
                        }
                    }
                    if (cabs1(h[((k + 1) - 1) + (k - 1) * ldh]) <= max(smlnum, REAL(ulp * tst1))) {
                        h12 = max(cabs1(h[((k + 1) - 1) + (k - 1) * ldh]), cabs1(h[(k - 1) + ((k + 1) - 1) * ldh]));
                        h21 = min(cabs1(h[((k + 1) - 1) + (k - 1) * ldh]), cabs1(h[(k - 1) + ((k + 1) - 1) * ldh]));
                        h11 = max(cabs1(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]), cabs1(h[(k - 1) + (k - 1) * ldh] - h[((k + 1) - 1) + ((k + 1) - 1) * ldh]));
                        h22 = min(cabs1(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]), cabs1(h[(k - 1) + (k - 1) * ldh] - h[((k + 1) - 1) + ((k + 1) - 1) * ldh]));
                        scl = h11 + h12;
                        tst2 = h22 * (h11 / scl);
                        //
                        if (tst2 == rzero || h21 * (h12 / scl) <= max(smlnum, REAL(ulp * tst2))) {
                            h[((k + 1) - 1) + (k - 1) * ldh] = zero;
                        }
                    }
                }
            }
            //
            //
            if (accum) {
                jbot = min(ndcol, kbot);
            } else if (wantt) {
                jbot = n;
            } else {
                jbot = kbot;
            }
            //
            for (m = mbot; m >= mtop; m = m - 1) {
                k = krcol + 2 * (m - 1);
                for (j = max(ktop, krcol + 2 * m); j <= jbot; j = j + 1) {
                    refsum = conj(v[(m - 1) * ldv]) * (h[((k + 1) - 1) + (j - 1) * ldh] + conj(v[(2 - 1) + (m - 1) * ldv]) * h[((k + 2) - 1) + (j - 1) * ldh] + conj(v[(3 - 1) + (m - 1) * ldv]) * h[((k + 3) - 1) + (j - 1) * ldh]);
                    h[((k + 1) - 1) + (j - 1) * ldh] = h[((k + 1) - 1) + (j - 1) * ldh] - refsum;
                    h[((k + 2) - 1) + (j - 1) * ldh] = h[((k + 2) - 1) + (j - 1) * ldh] - refsum * v[(2 - 1) + (m - 1) * ldv];
                    h[((k + 3) - 1) + (j - 1) * ldh] = h[((k + 3) - 1) + (j - 1) * ldh] - refsum * v[(3 - 1) + (m - 1) * ldv];
                }
            }
            //
            //
            if (accum) {
                //
                // .    with an efficient matrix-matrix
                //
                for (m = mbot; m >= mtop; m = m - 1) {
                    k = krcol + 2 * (m - 1);
                    kms = k - incol;
                    i2 = max((INTEGER)1, ktop - incol);
                    i2 = max(i2, kms - (krcol - incol) + 1);
                    i4 = min(kdu, krcol + 2 * (mbot - 1) - incol + 5);
                    for (j = i2; j <= i4; j = j + 1) {
                        refsum = v[(m - 1) * ldv] * (u[(j - 1) + ((kms + 1) - 1) * ldu] + v[(2 - 1) + (m - 1) * ldv] * u[(j - 1) + ((kms + 2) - 1) * ldu] + v[(3 - 1) + (m - 1) * ldv] * u[(j - 1) + ((kms + 3) - 1) * ldu]);
                        u[(j - 1) + ((kms + 1) - 1) * ldu] = u[(j - 1) + ((kms + 1) - 1) * ldu] - refsum;
                        u[(j - 1) + ((kms + 2) - 1) * ldu] = u[(j - 1) + ((kms + 2) - 1) * ldu] - refsum * conj(v[(2 - 1) + (m - 1) * ldv]);
                        u[(j - 1) + ((kms + 3) - 1) * ldu] = u[(j - 1) + ((kms + 3) - 1) * ldu] - refsum * conj(v[(3 - 1) + (m - 1) * ldv]);
                    }
                }
            } else if (wantz) {
                //
                // .    now by multiplying by reflections
                //
                for (m = mbot; m >= mtop; m = m - 1) {
                    k = krcol + 2 * (m - 1);
                    for (j = iloz; j <= ihiz; j = j + 1) {
                        refsum = v[(m - 1) * ldv] * (z[(j - 1) + ((k + 1) - 1) * ldz] + v[(2 - 1) + (m - 1) * ldv] * z[(j - 1) + ((k + 2) - 1) * ldz] + v[(3 - 1) + (m - 1) * ldv] * z[(j - 1) + ((k + 3) - 1) * ldz]);
                        z[(j - 1) + ((k + 1) - 1) * ldz] = z[(j - 1) + ((k + 1) - 1) * ldz] - refsum;
                        z[(j - 1) + ((k + 2) - 1) * ldz] = z[(j - 1) + ((k + 2) - 1) * ldz] - refsum * conj(v[(2 - 1) + (m - 1) * ldv]);
                        z[(j - 1) + ((k + 3) - 1) * ldz] = z[(j - 1) + ((k + 3) - 1) * ldz] - refsum * conj(v[(3 - 1) + (m - 1) * ldv]);
                    }
                }
            }
            //
            //
        }
        //
        // .    entries in H.  If required, use U to update Z as
        //
        if (accum) {
            if (wantt) {
                jtop = 1;
                jbot = n;
            } else {
                jtop = ktop;
                jbot = kbot;
            }
            k1 = max((INTEGER)1, ktop - incol);
            nu = (kdu - max((INTEGER)0, ndcol - kbot)) - k1 + 1;
            //
            //
            for (jcol = min(ndcol, kbot) + 1; jcol <= jbot; jcol = jcol + nh) {
                jlen = min(nh, jbot - jcol + 1);
                Cgemm("C", "N", nu, jlen, nu, one, &u[(k1 - 1) + (k1 - 1) * ldu], ldu, &h[((incol + k1) - 1) + (jcol - 1) * ldh], ldh, zero, wh, ldwh);
                Clacpy("ALL", nu, jlen, wh, ldwh, &h[((incol + k1) - 1) + (jcol - 1) * ldh], ldh);
            }
            //
            //
            for (jrow = jtop; jrow <= max(ktop, incol) - 1; jrow = jrow + nv) {
                jlen = min(nv, max(ktop, incol) - jrow);
                Cgemm("N", "N", jlen, nu, nu, one, &h[(jrow - 1) + ((incol + k1) - 1) * ldh], ldh, &u[(k1 - 1) + (k1 - 1) * ldu], ldu, zero, wv, ldwv);
                Clacpy("ALL", jlen, nu, wv, ldwv, &h[(jrow - 1) + ((incol + k1) - 1) * ldh], ldh);
            }
            //
            //
            if (wantz) {
                for (jrow = iloz; jrow <= ihiz; jrow = jrow + nv) {
                    jlen = min(nv, ihiz - jrow + 1);
                    Cgemm("N", "N", jlen, nu, nu, one, &z[(jrow - 1) + ((incol + k1) - 1) * ldz], ldz, &u[(k1 - 1) + (k1 - 1) * ldu], ldu, zero, wv, ldwv);
                    Clacpy("ALL", jlen, nu, wv, ldwv, &z[(jrow - 1) + ((incol + k1) - 1) * ldz], ldz);
                }
            }
        }
    }
    //
    //
}
