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

void Claqr4(bool const wantt, bool const wantz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, COMPLEX *h, INTEGER const ldh, COMPLEX *w, INTEGER const iloz, INTEGER const ihiz, COMPLEX *z, INTEGER const ldz, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
    COMPLEX cdum = 0.0;
    const COMPLEX one = (1.0, 0.0);
    const INTEGER ntiny = 15;
    INTEGER lwkopt = 0;
    str<2> jbcmpz = char0;
    INTEGER nwr = 0;
    INTEGER nsr = 0;
    INTEGER ls = 0;
    INTEGER ld = 0;
    INTEGER nmin = 0;
    INTEGER nibble = 0;
    INTEGER kacc22 = 0;
    INTEGER nwmax = 0;
    INTEGER nw = 0;
    INTEGER nsmax = 0;
    INTEGER ndfl = 0;
    const INTEGER kexsh = 6;
    INTEGER itmax = 0;
    INTEGER kbot = 0;
    INTEGER it = 0;
    INTEGER k = 0;
    const COMPLEX zero = (0.0, 0.0);
    INTEGER ktop = 0;
    INTEGER nh = 0;
    INTEGER nwupbd = 0;
    const INTEGER kexnw = 5;
    INTEGER kwtop = 0;
    INTEGER ndec = 0;
    INTEGER kv = 0;
    INTEGER kt = 0;
    INTEGER nho = 0;
    INTEGER kwv = 0;
    INTEGER nve = 0;
    INTEGER ks = 0;
    INTEGER ns = 0;
    INTEGER i = 0;
    const REAL wilk1 = 0.75e0;
    arr_2d<1, 1, COMPLEX> zdum(fill0);
    INTEGER inf = 0;
    REAL s = 0.0;
    COMPLEX aa = 0.0;
    COMPLEX cc = 0.0;
    COMPLEX bb = 0.0;
    COMPLEX dd = 0.0;
    const REAL two = 2.0;
    COMPLEX tr2 = 0.0;
    COMPLEX det = 0.0;
    COMPLEX rtdisc = 0.0;
    bool sorted = false;
    COMPLEX swap = 0.0;
    INTEGER kdu = 0;
    INTEGER ku = 0;
    INTEGER kwh = 0;
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
    //
    //     ==== Matrices of order NTINY or smaller must be processed by
    //     .    Clahqr because of insufficient subdiagonal scratch space.
    //     .    (This is a hard limit.) ====
    //
    //     ==== Exceptional deflation windows:  try to cure rare
    //     .    slow convergence by varying the size of the
    //     .    deflation window after KEXNW iterations. ====
    //
    //     ==== Exceptional shifts: try to cure rare slow convergence
    //     .    with ad-hoc exceptional shifts every KEXSH iterations.
    //     .    ====
    //
    //     ==== The constant WILK1 is used to form the exceptional
    //     .    shifts. ====
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    abs1[cdum - 1] = abs(cdum.real()) + abs(cdum.imag());
    //     ..
    //     .. Executable Statements ..
    info = 0;
    //
    //     ==== Quick return for N = 0: nothing to do. ====
    //
    if (n == 0) {
        work[1 - 1] = one;
        return;
    }
    //
    if (n <= ntiny) {
        //
        //        ==== Tiny matrices must use Clahqr. ====
        //
        lwkopt = 1;
        if (lwork != -1) {
            Clahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info);
        }
    } else {
        //
        //        ==== Use small bulge multi-shift QR with aggressive early
        //        .    deflation on larger-than-tiny matrices. ====
        //
        //        ==== Hope for the best. ====
        //
        info = 0;
        //
        //        ==== Set up job flags for iMlaenv. ====
        //
        if (wantt) {
            jbcmpz[(1 - 1)] = "S";
        } else {
            jbcmpz[(1 - 1)] = "E";
        }
        if (wantz) {
            jbcmpz[(2 - 1) + (2 - 1) * ldjbcmpz] = "V";
        } else {
            jbcmpz[(2 - 1) + (2 - 1) * ldjbcmpz] = "N";
        }
        //
        //        ==== NWR = recommended deflation window size.  At this
        //        .    point,  N .GT. NTINY = 15, so there is enough
        //        .    subdiagonal workspace for NWR.GE.2 as required.
        //        .    (In fact, there is enough subdiagonal space for
        //        .    NWR.GE.4.) ====
        //
        nwr = iMlaenv(13, "Claqr4", jbcmpz, n, ilo, ihi, lwork);
        nwr = max(2, nwr);
        nwr = min(ihi - ilo + 1, (n - 1) / 3, nwr);
        //
        //        ==== NSR = recommended number of simultaneous shifts.
        //        .    At this poINTEGER N .GT. NTINY = 15, so there is at
        //        .    enough subdiagonal workspace for NSR to be even
        //        .    and greater than or equal to two as required. ====
        //
        nsr = iMlaenv(15, "Claqr4", jbcmpz, n, ilo, ihi, lwork);
        nsr = min(nsr, (n - 3) / 6, ihi - ilo);
        nsr = max(2, nsr - mod(nsr, 2));
        //
        //        ==== Estimate optimal workspace ====
        //
        //        ==== Workspace query call to Claqr2 ====
        //
        Claqr2(wantt, wantz, n, ilo, ihi, nwr + 1, h, ldh, iloz, ihiz, z, ldz, ls, ld, w, h, ldh, n, h, ldh, n, h, ldh, work, -1);
        //
        //        ==== Optimal workspace = MAX(Claqr5, Claqr2) ====
        //
        lwkopt = max(3 * nsr / 2, int(work[1 - 1]));
        //
        //        ==== Quick return in case of workspace query. ====
        //
        if (lwork == -1) {
            work[1 - 1] = COMPLEX(lwkopt, 0);
            return;
        }
        //
        //        ==== Clahqr/Claqr0 crossover poINTEGER ====
        //
        nmin = iMlaenv(12, "Claqr4", jbcmpz, n, ilo, ihi, lwork);
        nmin = max(ntiny, nmin);
        //
        //        ==== Nibble crossover poINTEGER ====
        //
        nibble = iMlaenv(14, "Claqr4", jbcmpz, n, ilo, ihi, lwork);
        nibble = max(0, nibble);
        //
        //        ==== Accumulate reflections during ttswp?  Use block
        //        .    2-by-2 structure during matrix-matrix multiply? ====
        //
        kacc22 = iMlaenv(16, "Claqr4", jbcmpz, n, ilo, ihi, lwork);
        kacc22 = max(0, kacc22);
        kacc22 = min(2, kacc22);
        //
        //        ==== NWMAX = the largest possible deflation window for
        //        .    which there is sufficient workspace. ====
        //
        nwmax = min((n - 1) / 3, lwork / 2);
        nw = nwmax;
        //
        //        ==== NSMAX = the Largest number of simultaneous shifts
        //        .    for which there is sufficient workspace. ====
        //
        nsmax = min((n - 3) / 6, 2 * lwork / 3);
        nsmax = nsmax - mod(nsmax, 2);
        //
        //        ==== NDFL: an iteration count restarted at deflation. ====
        //
        ndfl = 1;
        //
        //        ==== ITMAX = iteration limit ====
        //
        itmax = max(30, 2 * kexsh) * max((INTEGER)10, (ihi - ilo + 1));
        //
        //        ==== Last row and column in the active block ====
        //
        kbot = ihi;
        //
        //        ==== Main Loop ====
        //
        for (it = 1; it <= itmax; it = it + 1) {
            //
            //           ==== Done when KBOT falls below ILO ====
            //
            if (kbot < ilo) {
                goto statement_80;
            }
            //
            //           ==== Locate active block ====
            //
            for (k = kbot; k >= ilo + 1; k = k - 1) {
                if (h[(k - 1) + ((k - 1) - 1) * ldh] == zero) {
                    goto statement_20;
                }
            }
            k = ilo;
        statement_20:
            ktop = k;
            //
            //           ==== Select deflation window size:
            //           .    Typical Case:
            //           .      If possible and advisable, nibble the entire
            //           .      active block.  If not, use size MIN(NWR,NWMAX)
            //           .      or MIN(NWR+1,NWMAX) depending upon which has
            //           .      the smaller corresponding subdiagonal entry
            //           .      (a heuristic).
            //           .
            //           .    Exceptional Case:
            //           .      If there have been no deflations in KEXNW or
            //           .      more iterations, then vary the deflation window
            //           .      size.   At first, because, larger windows are,
            //           .      in general, more powerful than smaller ones,
            //           .      rapidly increase the window to the maximum possible.
            //           .      Then, gradually reduce the window size. ====
            //
            nh = kbot - ktop + 1;
            nwupbd = min(nh, nwmax);
            if (ndfl < kexnw) {
                nw = min(nwupbd, nwr);
            } else {
                nw = min(nwupbd, 2 * nw);
            }
            if (nw < nwmax) {
                if (nw >= nh - 1) {
                    nw = nh;
                } else {
                    kwtop = kbot - nw + 1;
                    if (abs1[(h[(kwtop - 1) + ((kwtop - 1) - 1) * ldh]) - 1] > abs1[(h[((kwtop - 1) - 1) + ((kwtop - 2) - 1) * ldh]) - 1]) {
                        nw++;
                    }
                }
            }
            if (ndfl < kexnw) {
                ndec = -1;
            } else if (ndec >= 0 || nw >= nwupbd) {
                ndec++;
                if (nw - ndec < 2) {
                    ndec = 0;
                }
                nw = nw - ndec;
            }
            //
            //           ==== Aggressive early deflation:
            //           .    split workspace under the subdiagonal into
            //           .      - an nw-by-nw work array V in the lower
            //           .        left-hand-corner,
            //           .      - an NW-by-at-least-NW-but-more-is-better
            //           .        (NW-by-NHO) horizontal work array along
            //           .        the bottom edge,
            //           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
            //           .        vertical work array along the left-hand-edge.
            //           .        ====
            //
            kv = n - nw + 1;
            kt = nw + 1;
            nho = (n - nw - 1) - kt + 1;
            kwv = nw + 2;
            nve = (n - nw) - kwv + 1;
            //
            //           ==== Aggressive early deflation ====
            //
            Claqr2(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ls, ld, w, &h[(kv - 1)], ldh, nho, &h[(kv - 1) + (kt - 1) * ldh], ldh, nve, &h[(kwv - 1)], ldh, work, lwork);
            //
            //           ==== Adjust KBOT accounting for new deflations. ====
            //
            kbot = kbot - ld;
            //
            //           ==== KS points to the shifts. ====
            //
            ks = kbot - ls + 1;
            //
            //           ==== Skip an expensive QR sweep if there is a (partly
            //           .    heuristic) reason to expect that many eigenvalues
            //           .    will deflate without it.  Here, the QR sweep is
            //           .    skipped if many eigenvalues have just been deflated
            //           .    or if the remaining active block is small.
            //
            if ((ld == 0) || ((100 * ld <= nw * nibble) && (kbot - ktop + 1 > min(nmin, nwmax)))) {
                //
                //              ==== NS = nominal number of simultaneous shifts.
                //              .    This may be lowered (slightly) if Claqr2
                //              .    did not provide that many shifts. ====
                //
                ns = min(nsmax, nsr, max(2, kbot - ktop));
                ns = ns - mod(ns, 2);
                //
                //              ==== If there have been no deflations
                //              .    in a multiple of KEXSH iterations,
                //              .    then try exceptional shifts.
                //              .    Otherwise use shifts provided by
                //              .    Claqr2 above or from the eigenvalues
                //              .    of a trailing principal submatrix. ====
                //
                if (mod(ndfl, kexsh) == 0) {
                    ks = kbot - ns + 1;
                    for (i = kbot; i >= ks + 1; i = i - 2) {
                        w[i - 1] = h[(i - 1) + (i - 1) * ldh] + wilk1 * abs1[(h[(i - 1) + ((i - 1) - 1) * ldh]) - 1];
                        w[(i - 1) - 1] = w[i - 1];
                    }
                } else {
                    //
                    //                 ==== Got NS/2 or fewer shifts? Use Clahqr
                    //                 .    on a trailing principal submatrix to
                    //                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
                    //                 .    there is enough space below the subdiagonal
                    //                 .    to fit an NS-by-NS scratch array.) ====
                    //
                    if (kbot - ks + 1 <= ns / 2) {
                        ks = kbot - ns + 1;
                        kt = n - ns + 1;
                        Clacpy("A", ns, ns, &h[(ks - 1) + (ks - 1) * ldh], ldh, &h[(kt - 1)], ldh);
                        Clahqr(false, false, ns, 1, ns, &h[(kt - 1)], ldh, &w[ks - 1], 1, 1, zdum, 1, inf);
                        ks += inf;
                        //
                        //                    ==== In case of a rare QR failure use
                        //                    .    eigenvalues of the trailing 2-by-2
                        //                    .    principal submatrix.  Scale to avoid
                        //                    .    overflows, underflows and subnormals.
                        //                    .    (The scale factor S can not be zero,
                        //                    .    because H(KBOT,KBOT-1) is nonzero.) ====
                        //
                        if (ks >= kbot) {
                            s = abs1[(h[((kbot - 1) - 1) + ((kbot - 1) - 1) * ldh]) - 1] + abs1[(h[(kbot - 1) + ((kbot - 1) - 1) * ldh]) - 1] + abs1[(h[((kbot - 1) - 1) + (kbot - 1) * ldh]) - 1] + abs1[h[(kbot - 1) + (kbot - 1) * ldh] - 1];
                            aa = h[((kbot - 1) - 1) + ((kbot - 1) - 1) * ldh] / s;
                            cc = h[(kbot - 1) + ((kbot - 1) - 1) * ldh] / s;
                            bb = h[((kbot - 1) - 1) + (kbot - 1) * ldh] / s;
                            dd = h[(kbot - 1) + (kbot - 1) * ldh] / s;
                            tr2 = (aa + dd) / two;
                            det = (aa - tr2) * (dd - tr2) - bb * cc;
                            rtdisc = sqrt(-det);
                            w[(kbot - 1) - 1] = (tr2 + rtdisc) * s;
                            w[kbot - 1] = (tr2 - rtdisc) * s;
                            //
                            ks = kbot - 1;
                        }
                    }
                    //
                    if (kbot - ks + 1 > ns) {
                        //
                        //                    ==== Sort the shifts (Helps a little) ====
                        //
                        sorted = false;
                        for (k = kbot; k >= ks + 1; k = k - 1) {
                            if (sorted) {
                                goto statement_60;
                            }
                            sorted = true;
                            for (i = ks; i <= k - 1; i = i + 1) {
                                if (abs1[w[i - 1] - 1] < abs1[(w[(i + 1) - 1]) - 1]) {
                                    sorted = false;
                                    swap = w[i - 1];
                                    w[i - 1] = w[(i + 1) - 1];
                                    w[(i + 1) - 1] = swap;
                                }
                            }
                        }
                    statement_60:;
                    }
                }
                //
                //              ==== If there are only two shifts, then use
                //              .    only one.  ====
                //
                if (kbot - ks + 1 == 2) {
                    if (abs1[(w[kbot - 1] - h[(kbot - 1) + (kbot - 1) * ldh]) - 1] < abs1[(w[(kbot - 1) - 1] - h[(kbot - 1) + (kbot - 1) * ldh]) - 1]) {
                        w[(kbot - 1) - 1] = w[kbot - 1];
                    } else {
                        w[kbot - 1] = w[(kbot - 1) - 1];
                    }
                }
                //
                //              ==== Use up to NS of the the smallest magnitude
                //              .    shifts.  If there aren't NS shifts available,
                //              .    then use them all, possibly dropping one to
                //              .    make the number of shifts even. ====
                //
                ns = min(ns, kbot - ks + 1);
                ns = ns - mod(ns, 2);
                ks = kbot - ns + 1;
                //
                //              ==== Small-bulge multi-shift QR sweep:
                //              .    split workspace under the subdiagonal into
                //              .    - a KDU-by-KDU work array U in the lower
                //              .      left-hand-corner,
                //              .    - a KDU-by-at-least-KDU-but-more-is-better
                //              .      (KDU-by-NHo) horizontal work array WH along
                //              .      the bottom edge,
                //              .    - and an at-least-KDU-but-more-is-better-by-KDU
                //              .      (NVE-by-KDU) vertical work WV arrow along
                //              .      the left-hand-edge. ====
                //
                kdu = 2 * ns;
                ku = n - kdu + 1;
                kwh = kdu + 1;
                nho = (n - kdu + 1 - 4) - (kdu + 1) + 1;
                kwv = kdu + 4;
                nve = n - kdu - kwv + 1;
                //
                //              ==== Small-bulge multi-shift QR sweep ====
                //
                Claqr5(wantt, wantz, kacc22, n, ktop, kbot, ns, &w[ks - 1], h, ldh, iloz, ihiz, z, ldz, work, 3, &h[(ku - 1)], ldh, nve, &h[(kwv - 1)], ldh, nho, &h[(ku - 1) + (kwh - 1) * ldh], ldh);
            }
            //
            //           ==== Note progress (or the lack of it). ====
            //
            if (ld > 0) {
                ndfl = 1;
            } else {
                ndfl++;
            }
            //
            //           ==== End of main loop ====
        }
        //
        //        ==== Iteration limit exceeded.  Set INFO to show where
        //        .    the problem occurred and exit. ====
        //
        info = kbot;
    statement_80:;
    }
    //
    //     ==== Return the optimal value of LWORK. ====
    //
    work[1 - 1] = COMPLEX(lwkopt, 0);
    //
    //     ==== End of Claqr4 ====
    //
}
