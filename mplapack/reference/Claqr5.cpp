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

void Claqr5(bool const wantt, bool const wantz, INTEGER const kacc22, INTEGER const n, INTEGER const ktop, INTEGER const kbot, INTEGER const nshfts, COMPLEX *s, COMPLEX *h, INTEGER const ldh, INTEGER const iloz, INTEGER const ihiz, COMPLEX *z, INTEGER const ldz, COMPLEX *v, INTEGER const ldv, COMPLEX *u, INTEGER const ldu, INTEGER const nv, COMPLEX *wv, INTEGER const ldwv, INTEGER const nh, COMPLEX *wh, INTEGER const ldwh) {
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
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    COMPLEX cdum = 0.0;
    abs1[cdum - 1] = abs(cdum.real()) + abs(cdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     ==== If there are no shifts, then there is nothing to do. ====
    //
    if (nshfts < 2) {
        return;
    }
    //
    //     ==== If the active block is empty or 1-by-1, then there
    //     .    is nothing to do. ====
    //
    if (ktop >= kbot) {
        return;
    }
    //
    //     ==== NSHFTS is supposed to be even, but if it is odd,
    //     .    then simply reduce it by one.  ====
    //
    INTEGER ns = nshfts - mod(nshfts, 2);
    //
    //     ==== Machine constants for deflation ====
    //
    REAL safmin = Rlamch("SAFE MINIMUM");
    const REAL rone = 1.0;
    REAL safmax = rone / safmin;
    Rlabad(safmin, safmax);
    REAL ulp = Rlamch("PRECISION");
    REAL smlnum = safmin * (n.real() / ulp);
    //
    //     ==== Use accumulated reflections to update far-from-diagonal
    //     .    entries ? ====
    //
    bool accum = (kacc22 == 1) || (kacc22 == 2);
    //
    //     ==== clear trash ====
    //
    const COMPLEX zero = (0.0, 0.0);
    if (ktop + 2 <= kbot) {
        h[((ktop + 2) - 1) + (ktop - 1) * ldh] = zero;
    }
    //
    //     ==== NBMPS = number of 2-shift bulges in the chain ====
    //
    INTEGER nbmps = ns / 2;
    //
    //     ==== KDU = width of slab ====
    //
    INTEGER kdu = 4 * nbmps;
    //
    //     ==== Create and chase chains of NBMPS bulges ====
    //
    INTEGER incol = 0;
    INTEGER jtop = 0;
    INTEGER ndcol = 0;
    const COMPLEX one = (1.0, 0.0);
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
    arr_1d<3, COMPLEX> vt(fill0);
    INTEGER i2 = 0;
    INTEGER i4 = 0;
    INTEGER k1 = 0;
    INTEGER nu = 0;
    INTEGER jcol = 0;
    INTEGER jlen = 0;
    INTEGER jrow = 0;
    for (incol = ktop - 2 * nbmps + 1; incol <= kbot - 2; incol = incol + 2 * nbmps) {
        //
        //        JTOP = Index from which updates from the right start.
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
        //        ==== Near-the-diagonal bulge chase.  The following loop
        //        .    performs the near-the-diagonal part of a small bulge
        //        .    multi-shift QR sweep.  Each 4*NBMPS column diagonal
        //        .    chunk extends from column INCOL to column NDCOL
        //        .    (including both column INCOL and column NDCOL). The
        //        .    following loop chases a 2*NBMPS+1 column long chain of
        //        .    NBMPS bulges 2*NBMPS columns to the right.  (INCOL
        //        .    may be less than KTOP and and NDCOL may be greater than
        //        .    KBOT indicating phantom columns from which to chase
        //        .    bulges before they are actually introduced or to which
        //        .    to chase bulges beyond column KBOT.)  ====
        //
        for (krcol = incol; krcol <= min(incol + 2 * nbmps - 1, kbot - 2); krcol = krcol + 1) {
            //
            //           ==== Bulges number MTOP to MBOT are active REAL implicit
            //           .    shift bulges.  There may or may not also be small
            //           .    2-by-2 bulge, if there is room.  The inactive bulges
            //           .    (if any) must wait until the active bulges have moved
            //           .    down the diagonal to make room.  The phantom matrix
            //           .    paradigm described above helps keep track.  ====
            //
            mtop = max((INTEGER)1, (ktop - krcol) / 2 + 1);
            mbot = min(nbmps, (kbot - krcol - 1) / 2);
            m22 = mbot + 1;
            bmp22 = (mbot < nbmps) && (krcol + 2 * (m22 - 1)) == (kbot - 2);
            //
            //           ==== Generate reflections to chase the chain right
            //           .    one column.  (The minimum value of K is KTOP-1.) ====
            //
            if (bmp22) {
                //
                //              ==== Special case: 2-by-2 reflection at bottom treated
                //              .    separately ====
                //
                k = krcol + 2 * (m22 - 1);
                if (k == ktop - 1) {
                    Claqr1(2, &h[((k + 1) - 1) + ((k + 1) - 1) * ldh], ldh, s[(2 * m22 - 1) - 1], s[(2 * m22) - 1], &v[(m22 - 1) * ldv]);
                    beta = v[(m22 - 1) * ldv];
                    Clarfg(2, beta, &v[(2 - 1) + (m22 - 1) * ldv], 1, &v[(m22 - 1) * ldv]);
                } else {
                    beta = h[((k + 1) - 1) + (k - 1) * ldh];
                    v[(2 - 1) + (m22 - 1) * ldv] = h[((k + 2) - 1) + (k - 1) * ldh];
                    Clarfg(2, beta, &v[(2 - 1) + (m22 - 1) * ldv], 1, &v[(m22 - 1) * ldv]);
                    h[((k + 1) - 1) + (k - 1) * ldh] = beta;
                    h[((k + 2) - 1) + (k - 1) * ldh] = zero;
                }
                //
                //              ==== Perform update from right within
                //              .    computational window. ====
                //
                for (j = jtop; j <= min(kbot, k + 3); j = j + 1) {
                    refsum = v[(m22 - 1) * ldv] * (h[(j - 1) + ((k + 1) - 1) * ldh] + v[(2 - 1) + (m22 - 1) * ldv] * h[(j - 1) + ((k + 2) - 1) * ldh]);
                    h[(j - 1) + ((k + 1) - 1) * ldh] = h[(j - 1) + ((k + 1) - 1) * ldh] - refsum;
                    h[(j - 1) + ((k + 2) - 1) * ldh] = h[(j - 1) + ((k + 2) - 1) * ldh] - refsum * conj(v[(2 - 1) + (m22 - 1) * ldv]);
                }
                //
                //              ==== Perform update from left within
                //              .    computational window. ====
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
                //              ==== The following convergence test requires that
                //              .    the tradition small-compared-to-nearby-diagonals
                //              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
                //              .    criteria both be satisfied.  The latter improves
                //              .    accuracy in some examples. Falling back on an
                //              .    alternate convergence criterion when TST1 or TST2
                //              .    is zero (as done here) is traditional but probably
                //              .    unnecessary. ====
                //
                if (k >= ktop) {
                    if (h[((k + 1) - 1) + (k - 1) * ldh] != zero) {
                        tst1 = abs1[h[(k - 1) + (k - 1) * ldh] - 1] + abs1[(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1];
                        if (tst1 == rzero) {
                            if (k >= ktop + 1) {
                                tst1 += abs1[(h[(k - 1) + ((k - 1) - 1) * ldh]) - 1];
                            }
                            if (k >= ktop + 2) {
                                tst1 += abs1[(h[(k - 1) + ((k - 2) - 1) * ldh]) - 1];
                            }
                            if (k >= ktop + 3) {
                                tst1 += abs1[(h[(k - 1) + ((k - 3) - 1) * ldh]) - 1];
                            }
                            if (k <= kbot - 2) {
                                tst1 += abs1[(h[((k + 2) - 1) + ((k + 1) - 1) * ldh]) - 1];
                            }
                            if (k <= kbot - 3) {
                                tst1 += abs1[(h[((k + 3) - 1) + ((k + 1) - 1) * ldh]) - 1];
                            }
                            if (k <= kbot - 4) {
                                tst1 += abs1[(h[((k + 4) - 1) + ((k + 1) - 1) * ldh]) - 1];
                            }
                        }
                        if (abs1[(h[((k + 1) - 1) + (k - 1) * ldh]) - 1] <= max(smlnum, ulp * tst1)) {
                            h12 = max(abs1[(h[((k + 1) - 1) + (k - 1) * ldh]) - 1], abs1[(h[(k - 1) + ((k + 1) - 1) * ldh]) - 1]);
                            h21 = min(abs1[(h[((k + 1) - 1) + (k - 1) * ldh]) - 1], abs1[(h[(k - 1) + ((k + 1) - 1) * ldh]) - 1]);
                            h11 = max(abs1[(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1], abs1[(h[(k - 1) + (k - 1) * ldh] - h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1]);
                            h22 = min(abs1[(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1], abs1[(h[(k - 1) + (k - 1) * ldh] - h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1]);
                            scl = h11 + h12;
                            tst2 = h22 * (h11 / scl);
                            //
                            if (tst2 == rzero || h21 * (h12 / scl) <= max(smlnum, ulp * tst2)) {
                                h[((k + 1) - 1) + (k - 1) * ldh] = zero;
                            }
                        }
                    }
                }
                //
                //              ==== Accumulate orthogonal transformations. ====
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
            //           ==== Normal case: Chain of 3-by-3 reflections ====
            //
            for (m = mbot; m >= mtop; m = m - 1) {
                k = krcol + 2 * (m - 1);
                if (k == ktop - 1) {
                    Claqr1(3, &h[(ktop - 1) + (ktop - 1) * ldh], ldh, s[(2 * m - 1) - 1], s[(2 * m) - 1], &v[(m - 1) * ldv]);
                    alpha = v[(m - 1) * ldv];
                    Clarfg(3, alpha, &v[(2 - 1) + (m - 1) * ldv], 1, &v[(m - 1) * ldv]);
                } else {
                    //
                    //                 ==== Perform delayed transformation of row below
                    //                 .    Mth bulge. Exploit fact that first two elements
                    //                 .    of row are actually zero. ====
                    //
                    refsum = v[(m - 1) * ldv] * v[(3 - 1) + (m - 1) * ldv] * h[((k + 3) - 1) + ((k + 2) - 1) * ldh];
                    h[((k + 3) - 1) + (k - 1) * ldh] = -refsum;
                    h[((k + 3) - 1) + ((k + 1) - 1) * ldh] = -refsum * conj(v[(2 - 1) + (m - 1) * ldv]);
                    h[((k + 3) - 1) + ((k + 2) - 1) * ldh] = h[((k + 3) - 1) + ((k + 2) - 1) * ldh] - refsum * conj(v[(3 - 1) + (m - 1) * ldv]);
                    //
                    //                 ==== Calculate reflection to move
                    //                 .    Mth bulge one step. ====
                    //
                    beta = h[((k + 1) - 1) + (k - 1) * ldh];
                    v[(2 - 1) + (m - 1) * ldv] = h[((k + 2) - 1) + (k - 1) * ldh];
                    v[(3 - 1) + (m - 1) * ldv] = h[((k + 3) - 1) + (k - 1) * ldh];
                    Clarfg(3, beta, &v[(2 - 1) + (m - 1) * ldv], 1, &v[(m - 1) * ldv]);
                    //
                    //                 ==== A Bulge may collapse because of vigilant
                    //                 .    deflation or destructive underflow.  In the
                    //                 .    underflow case, try the two-small-subdiagonals
                    //                 .    trick to try to reinflate the bulge.  ====
                    //
                    if (h[((k + 3) - 1) + (k - 1) * ldh] != zero || h[((k + 3) - 1) + ((k + 1) - 1) * ldh] != zero || h[((k + 3) - 1) + ((k + 2) - 1) * ldh] == zero) {
                        //
                        //                    ==== Typical case: not collapsed (yet). ====
                        //
                        h[((k + 1) - 1) + (k - 1) * ldh] = beta;
                        h[((k + 2) - 1) + (k - 1) * ldh] = zero;
                        h[((k + 3) - 1) + (k - 1) * ldh] = zero;
                    } else {
                        //
                        //                    ==== Atypical case: collapsed.  Attempt to
                        //                    .    reintroduce ignoring H(K+1,K) and H(K+2,K).
                        //                    .    If the fill resulting from the new
                        //                    .    reflector is too large, then abandon it.
                        //                    .    Otherwise, use the new one. ====
                        //
                        Claqr1(3, &h[((k + 1) - 1) + ((k + 1) - 1) * ldh], ldh, s[(2 * m - 1) - 1], s[(2 * m) - 1], vt);
                        alpha = vt[1 - 1];
                        Clarfg(3, alpha, vt[2 - 1], 1, vt[1 - 1]);
                        refsum = conj(vt[1 - 1]) * (h[((k + 1) - 1) + (k - 1) * ldh] + conj(vt[2 - 1]) * h[((k + 2) - 1) + (k - 1) * ldh]);
                        //
                        if (abs1[(h[((k + 2) - 1) + (k - 1) * ldh] - refsum * vt[2 - 1]) - 1] + abs1[(refsum * vt[3 - 1]) - 1] > ulp * (abs1[h[(k - 1) + (k - 1) * ldh] - 1] + abs1[(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1] + abs1[(h[((k + 2) - 1) + ((k + 2) - 1) * ldh]) - 1])) {
                            //
                            //                       ==== Starting a new bulge here would
                            //                       .    create non-negligible fill.  Use
                            //                       .    the old one with trepidation. ====
                            //
                            h[((k + 1) - 1) + (k - 1) * ldh] = beta;
                            h[((k + 2) - 1) + (k - 1) * ldh] = zero;
                            h[((k + 3) - 1) + (k - 1) * ldh] = zero;
                        } else {
                            //
                            //                       ==== Starting a new bulge here would
                            //                       .    create only negligible fill.
                            //                       .    Replace the old reflector with
                            //                       .    the new one. ====
                            //
                            h[((k + 1) - 1) + (k - 1) * ldh] = h[((k + 1) - 1) + (k - 1) * ldh] - refsum;
                            h[((k + 2) - 1) + (k - 1) * ldh] = zero;
                            h[((k + 3) - 1) + (k - 1) * ldh] = zero;
                            v[(m - 1) * ldv] = vt[1 - 1];
                            v[(2 - 1) + (m - 1) * ldv] = vt[2 - 1];
                            v[(3 - 1) + (m - 1) * ldv] = vt[3 - 1];
                        }
                    }
                }
                //
                //              ====  Apply reflection from the right and
                //              .     the first column of update from the left.
                //              .     These updates are required for the vigilant
                //              .     deflation check. We still delay most of the
                //              .     updates from the left for efficiency. ====
                //
                for (j = jtop; j <= min(kbot, k + 3); j = j + 1) {
                    refsum = v[(m - 1) * ldv] * (h[(j - 1) + ((k + 1) - 1) * ldh] + v[(2 - 1) + (m - 1) * ldv] * h[(j - 1) + ((k + 2) - 1) * ldh] + v[(3 - 1) + (m - 1) * ldv] * h[(j - 1) + ((k + 3) - 1) * ldh]);
                    h[(j - 1) + ((k + 1) - 1) * ldh] = h[(j - 1) + ((k + 1) - 1) * ldh] - refsum;
                    h[(j - 1) + ((k + 2) - 1) * ldh] = h[(j - 1) + ((k + 2) - 1) * ldh] - refsum * conj(v[(2 - 1) + (m - 1) * ldv]);
                    h[(j - 1) + ((k + 3) - 1) * ldh] = h[(j - 1) + ((k + 3) - 1) * ldh] - refsum * conj(v[(3 - 1) + (m - 1) * ldv]);
                }
                //
                //              ==== Perform update from left for subsequent
                //              .    column. ====
                //
                refsum = conj(v[(m - 1) * ldv]) * (h[((k + 1) - 1) + ((k + 1) - 1) * ldh] + conj(v[(2 - 1) + (m - 1) * ldv]) * h[((k + 2) - 1) + ((k + 1) - 1) * ldh] + conj(v[(3 - 1) + (m - 1) * ldv]) * h[((k + 3) - 1) + ((k + 1) - 1) * ldh]);
                h[((k + 1) - 1) + ((k + 1) - 1) * ldh] = h[((k + 1) - 1) + ((k + 1) - 1) * ldh] - refsum;
                h[((k + 2) - 1) + ((k + 1) - 1) * ldh] = h[((k + 2) - 1) + ((k + 1) - 1) * ldh] - refsum * v[(2 - 1) + (m - 1) * ldv];
                h[((k + 3) - 1) + ((k + 1) - 1) * ldh] = h[((k + 3) - 1) + ((k + 1) - 1) * ldh] - refsum * v[(3 - 1) + (m - 1) * ldv];
                //
                //              ==== The following convergence test requires that
                //              .    the tradition small-compared-to-nearby-diagonals
                //              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
                //              .    criteria both be satisfied.  The latter improves
                //              .    accuracy in some examples. Falling back on an
                //              .    alternate convergence criterion when TST1 or TST2
                //              .    is zero (as done here) is traditional but probably
                //              .    unnecessary. ====
                //
                if (k < ktop) {
                    continue;
                }
                if (h[((k + 1) - 1) + (k - 1) * ldh] != zero) {
                    tst1 = abs1[h[(k - 1) + (k - 1) * ldh] - 1] + abs1[(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1];
                    if (tst1 == rzero) {
                        if (k >= ktop + 1) {
                            tst1 += abs1[(h[(k - 1) + ((k - 1) - 1) * ldh]) - 1];
                        }
                        if (k >= ktop + 2) {
                            tst1 += abs1[(h[(k - 1) + ((k - 2) - 1) * ldh]) - 1];
                        }
                        if (k >= ktop + 3) {
                            tst1 += abs1[(h[(k - 1) + ((k - 3) - 1) * ldh]) - 1];
                        }
                        if (k <= kbot - 2) {
                            tst1 += abs1[(h[((k + 2) - 1) + ((k + 1) - 1) * ldh]) - 1];
                        }
                        if (k <= kbot - 3) {
                            tst1 += abs1[(h[((k + 3) - 1) + ((k + 1) - 1) * ldh]) - 1];
                        }
                        if (k <= kbot - 4) {
                            tst1 += abs1[(h[((k + 4) - 1) + ((k + 1) - 1) * ldh]) - 1];
                        }
                    }
                    if (abs1[(h[((k + 1) - 1) + (k - 1) * ldh]) - 1] <= max(smlnum, ulp * tst1)) {
                        h12 = max(abs1[(h[((k + 1) - 1) + (k - 1) * ldh]) - 1], abs1[(h[(k - 1) + ((k + 1) - 1) * ldh]) - 1]);
                        h21 = min(abs1[(h[((k + 1) - 1) + (k - 1) * ldh]) - 1], abs1[(h[(k - 1) + ((k + 1) - 1) * ldh]) - 1]);
                        h11 = max(abs1[(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1], abs1[(h[(k - 1) + (k - 1) * ldh] - h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1]);
                        h22 = min(abs1[(h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1], abs1[(h[(k - 1) + (k - 1) * ldh] - h[((k + 1) - 1) + ((k + 1) - 1) * ldh]) - 1]);
                        scl = h11 + h12;
                        tst2 = h22 * (h11 / scl);
                        //
                        if (tst2 == rzero || h21 * (h12 / scl) <= max(smlnum, ulp * tst2)) {
                            h[((k + 1) - 1) + (k - 1) * ldh] = zero;
                        }
                    }
                }
            }
            //
            //           ==== Multiply H by reflections from the left ====
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
            //           ==== Accumulate orthogonal transformations. ====
            //
            if (accum) {
                //
                //              ==== Accumulate U. (If needed, update Z later
                //              .    with an efficient matrix-matrix
                //              .    multiply.) ====
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
                //              ==== U is not accumulated, so update Z
                //              .    now by multiplying by reflections
                //              .    from the right. ====
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
            //           ==== End of near-the-diagonal bulge chase. ====
            //
        }
        //
        //        ==== Use U (if accumulated) to update far-from-diagonal
        //        .    entries in H.  If required, use U to update Z as
        //        .    well. ====
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
            nu = (kdu - max(0, ndcol - kbot)) - k1 + 1;
            //
            //           ==== Horizontal Multiply ====
            //
            for (jcol = min(ndcol, kbot) + 1; jcol <= jbot; jcol = jcol + nh) {
                jlen = min(nh, jbot - jcol + 1);
                Cgemm("C", "N", nu, jlen, nu, one, u[(k1 - 1) + (k1 - 1) * ldu], ldu, &h[((incol + k1) - 1) + (jcol - 1) * ldh], ldh, zero, wh, ldwh);
                Clacpy("ALL", nu, jlen, wh, ldwh, &h[((incol + k1) - 1) + (jcol - 1) * ldh], ldh);
            }
            //
            //           ==== Vertical multiply ====
            //
            for (jrow = jtop; jrow <= max(ktop, incol) - 1; jrow = jrow + nv) {
                jlen = min(nv, max(ktop, incol) - jrow);
                Cgemm("N", "N", jlen, nu, nu, one, &h[(jrow - 1) + ((incol + k1) - 1) * ldh], ldh, u[(k1 - 1) + (k1 - 1) * ldu], ldu, zero, wv, ldwv);
                Clacpy("ALL", jlen, nu, wv, ldwv, &h[(jrow - 1) + ((incol + k1) - 1) * ldh], ldh);
            }
            //
            //           ==== Z multiply (also vertical) ====
            //
            if (wantz) {
                for (jrow = iloz; jrow <= ihiz; jrow = jrow + nv) {
                    jlen = min(nv, ihiz - jrow + 1);
                    Cgemm("N", "N", jlen, nu, nu, one, &z[(jrow - 1) + ((incol + k1) - 1) * ldz], ldz, u[(k1 - 1) + (k1 - 1) * ldu], ldu, zero, wv, ldwv);
                    Clacpy("ALL", jlen, nu, wv, ldwv, &z[(jrow - 1) + ((incol + k1) - 1) * ldz], ldz);
                }
            }
        }
    }
    //
    //     ==== End of Claqr5 ====
    //
}
