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

void Rlarrd(const char *range, const char *order, INTEGER const n, REAL const vl, REAL const vu, INTEGER const il, INTEGER const iu, REAL *gers, REAL const reltol, REAL *d, REAL *e, REAL *e2, REAL const pivmin, INTEGER const nsplit, INTEGER *isplit, INTEGER &m, REAL *w, REAL *werr, REAL &wl, REAL &wu, INTEGER *iblock, INTEGER *indexw, REAL *work, INTEGER *iwork, INTEGER &info) {
    const INTEGER allrng = 1;
    INTEGER irange = 0;
    const INTEGER valrng = 2;
    const INTEGER indrng = 3;
    bool ncnvrg = false;
    bool toofew = false;
    REAL eps = 0.0;
    REAL uflow = 0.0;
    const REAL zero = 0.0;
    INTEGER nb = 0;
    REAL gl = 0.0;
    REAL gu = 0.0;
    INTEGER i = 0;
    REAL tnorm = 0.0;
    const REAL two = 2.0;
    const REAL fudge = two;
    REAL rtoli = 0.0;
    REAL atoli = 0.0;
    INTEGER itmax = 0;
    INTEGER iout = 0;
    INTEGER iinfo = 0;
    REAL wlu = 0.0;
    INTEGER nwl = 0;
    REAL wul = 0.0;
    INTEGER nwu = 0;
    INTEGER iend = 0;
    INTEGER jblk = 0;
    INTEGER ioff = 0;
    INTEGER ibegin = 0;
    INTEGER in = 0;
    REAL tmp1 = 0.0;
    INTEGER j = 0;
    INTEGER idumma[1];
    INTEGER im = 0;
    INTEGER iwoff = 0;
    const REAL one = 1.0;
    const REAL half = one / two;
    REAL tmp2 = 0.0;
    INTEGER ib = 0;
    INTEGER je = 0;
    INTEGER idiscl = 0;
    INTEGER idiscu = 0;
    INTEGER jee = 0;
    REAL wkill = 0.0;
    INTEGER jdisc = 0;
    INTEGER iw = 0;
    INTEGER ie = 0;
    INTEGER itmp1 = 0;
    INTEGER itmp2 = 0;
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    //     Decode RANGE
    //
    if (Mlsame(range, "A")) {
        irange = allrng;
    } else if (Mlsame(range, "V")) {
        irange = valrng;
    } else if (Mlsame(range, "I")) {
        irange = indrng;
    } else {
        irange = 0;
    }
    //
    //     Check for Errors
    //
    if (irange <= 0) {
        info = -1;
    } else if (!(Mlsame(order, "B") || Mlsame(order, "E"))) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (irange == valrng) {
        if (vl >= vu) {
            info = -5;
        }
    } else if (irange == indrng && (il < 1 || il > max((INTEGER)1, n))) {
        info = -6;
    } else if (irange == indrng && (iu < min(n, il) || iu > n)) {
        info = -7;
    }
    //
    if (info != 0) {
        return;
    }
    //
    //     Initialize error flags
    info = 0;
    ncnvrg = false;
    toofew = false;
    //
    //     Quick return if possible
    m = 0;
    if (n == 0) {
        return;
    }
    //
    //     Simplification:
    if (irange == indrng && il == 1 && iu == n) {
        irange = 1;
    }
    //
    //     Get machine constants
    eps = Rlamch("P");
    uflow = Rlamch("U");
    //
    //     Special Case when N=1
    //     Treat case of 1x1 matrix for quick return
    if (n == 1) {
        if ((irange == allrng) || ((irange == valrng) && (d[1 - 1] > vl) && (d[1 - 1] <= vu)) || ((irange == indrng) && (il == 1) && (iu == 1))) {
            m = 1;
            w[1 - 1] = d[1 - 1];
            //           The computation error of the eigenvalue is zero
            werr[1 - 1] = zero;
            iblock[1 - 1] = 1;
            indexw[1 - 1] = 1;
        }
        return;
    }
    //
    //     NB is the minimum vector length for vector bisection, or 0
    //     if only scalar is to be done.
    nb = iMlaenv(1, "Rstebz", " ", n, -1, -1, -1);
    if (nb <= 1) {
        nb = 0;
    }
    //
    //     Find global spectral radius
    gl = d[1 - 1];
    gu = d[1 - 1];
    for (i = 1; i <= n; i = i + 1) {
        gl = min(gl, gers[(2 * i - 1) - 1]);
        gu = max(gu, gers[(2 * i) - 1]);
    }
    //     Compute global Gerschgorin bounds and spectral diameter
    tnorm = max(abs(gl), abs(gu));
    gl = gl - fudge * tnorm * eps * n - fudge * two * pivmin;
    gu += fudge * tnorm * eps * n + fudge * two * pivmin;
    //     [JAN/28/2009] remove the line below since SPDIAM variable not use
    //     SPDIAM = GU - GL
    //     Input arguments for Rlaebz:
    //     The relative tolerance.  An interval (a,b] lies within
    //     "relative tolerance" if  b-a < RELTOL*max(|a|,|b|),
    rtoli = reltol;
    //     Set the absolute tolerance for interval convergence to zero to force
    //     interval convergence based on relative size of the interval.
    //     This is dangerous because intervals might not converge when RELTOL is
    //     small. But at least a very small number should be selected so that for
    //     strongly graded matrices, the code can get relatively accurate
    //     eigenvalues.
    atoli = fudge * two * uflow + fudge * two * pivmin;
    //
    if (irange == indrng) {
        //
        //        RANGE='I': Compute an interval containing eigenvalues
        //        IL through IU. The initial interval [GL,GU] from the global
        //        Gerschgorin bounds GL and GU is refined by Rlaebz.
        itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
        if (itmax >= 1024)
            itmax = 1024; // XXX itmax can be too large for MPFR (=10^8)
        work[(n + 1) - 1] = gl;
        work[(n + 2) - 1] = gl;
        work[(n + 3) - 1] = gu;
        work[(n + 4) - 1] = gu;
        work[(n + 5) - 1] = gl;
        work[(n + 6) - 1] = gu;
        iwork[1 - 1] = -1;
        iwork[2 - 1] = -1;
        iwork[3 - 1] = n + 1;
        iwork[4 - 1] = n + 1;
        iwork[5 - 1] = il - 1;
        iwork[6 - 1] = iu;
        //
        Rlaebz(3, itmax, n, 2, 2, nb, atoli, rtoli, pivmin, d, e, e2, &iwork[5 - 1], &work[(n + 1) - 1], &work[(n + 5) - 1], iout, iwork, w, iblock, iinfo);
        if (iinfo != 0) {
            info = iinfo;
            return;
        }
        //        On exit, output intervals may not be ordered by ascending negcount
        if (iwork[6 - 1] == iu) {
            wl = work[(n + 1) - 1];
            wlu = work[(n + 3) - 1];
            nwl = iwork[1 - 1];
            wu = work[(n + 4) - 1];
            wul = work[(n + 2) - 1];
            nwu = iwork[4 - 1];
        } else {
            wl = work[(n + 2) - 1];
            wlu = work[(n + 4) - 1];
            nwl = iwork[2 - 1];
            wu = work[(n + 3) - 1];
            wul = work[(n + 1) - 1];
            nwu = iwork[3 - 1];
        }
        //        On exit, the interval [WL, WLU] contains a value with negcount NWL,
        //        and [WUL, WU] contains a value with negcount NWU.
        if (nwl < 0 || nwl >= n || nwu < 1 || nwu > n) {
            info = 4;
            return;
        }
        //
    } else if (irange == valrng) {
        wl = vl;
        wu = vu;
        //
    } else if (irange == allrng) {
        wl = gl;
        wu = gu;
    }
    //
    //     Find Eigenvalues -- Loop Over blocks and recompute NWL and NWU.
    //     NWL accumulates the number of eigenvalues .le. WL,
    //     NWU accumulates the number of eigenvalues .le. WU
    m = 0;
    iend = 0;
    info = 0;
    nwl = 0;
    nwu = 0;
    //
    for (jblk = 1; jblk <= nsplit; jblk = jblk + 1) {
        ioff = iend;
        ibegin = ioff + 1;
        iend = isplit[jblk - 1];
        in = iend - ioff;
        //
        if (in == 1) {
            //           1x1 block
            if (wl >= d[ibegin - 1] - pivmin) {
                nwl++;
            }
            if (wu >= d[ibegin - 1] - pivmin) {
                nwu++;
            }
            if (irange == allrng || (wl < d[ibegin - 1] - pivmin && wu >= d[ibegin - 1] - pivmin)) {
                m++;
                w[m - 1] = d[ibegin - 1];
                werr[m - 1] = zero;
                //              The gap for a single block doesn't matter for the later
                //              algorithm and is assigned an arbitrary large value
                iblock[m - 1] = jblk;
                indexw[m - 1] = 1;
            }
            //
            //        Disabled 2x2 case because of a failure on the following matrix
            //        RANGE = 'I', IL = IU = 4
            //          Original Tridiagonal, d = [
            //           -0.150102010615740E+00
            //           -0.849897989384260E+00
            //           -0.128208148052635E-15
            //            0.128257718286320E-15
            //          ];
            //          e = [
            //           -0.357171383266986E+00
            //           -0.180411241501588E-15
            //           -0.175152352710251E-15
            //          ];
            //
            //         ELSE IF( IN.EQ.2 ) THEN
            //*           2x2 block
            //            DISC = SQRT( (HALF*(D(IBEGIN)-D(IEND)))**2 + E(IBEGIN)**2 )
            //            TMP1 = HALF*(D(IBEGIN)+D(IEND))
            //            L1 = TMP1 - DISC
            //            IF( WL.GE. L1-PIVMIN )
            //     $         NWL = NWL + 1
            //            IF( WU.GE. L1-PIVMIN )
            //     $         NWU = NWU + 1
            //            IF( IRANGE.EQ.ALLRNG .OR. ( WL.LT.L1-PIVMIN .AND. WU.GE.
            //     $          L1-PIVMIN ) ) THEN
            //               M = M + 1
            //               W( M ) = L1
            //*              The uncertainty of eigenvalues of a 2x2 matrix is very small
            //               WERR( M ) = EPS * ABS( W( M ) ) * TWO
            //               IBLOCK( M ) = JBLK
            //               INDEXW( M ) = 1
            //            ENDIF
            //            L2 = TMP1 + DISC
            //            IF( WL.GE. L2-PIVMIN )
            //     $         NWL = NWL + 1
            //            IF( WU.GE. L2-PIVMIN )
            //     $         NWU = NWU + 1
            //            IF( IRANGE.EQ.ALLRNG .OR. ( WL.LT.L2-PIVMIN .AND. WU.GE.
            //     $          L2-PIVMIN ) ) THEN
            //               M = M + 1
            //               W( M ) = L2
            //*              The uncertainty of eigenvalues of a 2x2 matrix is very small
            //               WERR( M ) = EPS * ABS( W( M ) ) * TWO
            //               IBLOCK( M ) = JBLK
            //               INDEXW( M ) = 2
            //            ENDIF
        } else {
            //           General Case - block of size IN >= 2
            //           Compute local Gerschgorin interval and use it as the initial
            //           interval for Rlaebz
            gu = d[ibegin - 1];
            gl = d[ibegin - 1];
            tmp1 = zero;
            //
            for (j = ibegin; j <= iend; j = j + 1) {
                gl = min(gl, gers[(2 * j - 1) - 1]);
                gu = max(gu, gers[(2 * j) - 1]);
            }
            //           [JAN/28/2009]
            //           change SPDIAM by TNORM in lines 2 and 3 thereafter
            //           line 1: remove computation of SPDIAM (not useful anymore)
            //           SPDIAM = GU - GL
            //           GL = GL - FUDGE*SPDIAM*EPS*IN - FUDGE*PIVMIN
            //           GU = GU + FUDGE*SPDIAM*EPS*IN + FUDGE*PIVMIN
            gl = gl - fudge * tnorm * eps * in - fudge * pivmin;
            gu += fudge * tnorm * eps * in + fudge * pivmin;
            //
            if (irange > 1) {
                if (gu < wl) {
                    //                 the local block contains none of the wanted eigenvalues
                    nwl += in;
                    nwu += in;
                    goto statement_70;
                }
                //              refine search interval if possible, only range (WL,WU] matters
                gl = max(gl, wl);
                gu = min(gu, wu);
                if (gl >= gu) {
                    goto statement_70;
                }
            }
            //
            //           Find negcount of initial interval boundaries GL and GU
            work[(n + 1) - 1] = gl;
            work[(n + in + 1) - 1] = gu;
            Rlaebz(1, 0, in, in, 1, nb, atoli, rtoli, pivmin, &d[ibegin - 1], &e[ibegin - 1], &e2[ibegin - 1], idumma, &work[(n + 1) - 1], &work[(n + 2 * in + 1) - 1], im, iwork, &w[(m + 1) - 1], &iblock[(m + 1) - 1], iinfo);
            if (iinfo != 0) {
                info = iinfo;
                return;
            }
            //
            nwl += iwork[1 - 1];
            nwu += iwork[(in + 1) - 1];
            iwoff = m - iwork[1 - 1];
            //
            //           Compute Eigenvalues
            itmax = castINTEGER((log(gu - gl + pivmin) - log(pivmin)) / log(two)) + 2;
            Rlaebz(2, itmax, in, in, 1, nb, atoli, rtoli, pivmin, &d[ibegin - 1], &e[ibegin - 1], &e2[ibegin - 1], idumma, &work[(n + 1) - 1], &work[(n + 2 * in + 1) - 1], iout, iwork, &w[(m + 1) - 1], &iblock[(m + 1) - 1], iinfo);
            if (iinfo != 0) {
                info = iinfo;
                return;
            }
            //
            //           Copy eigenvalues into W and IBLOCK
            //           Use -JBLK for block number for unconverged eigenvalues.
            //           Loop over the number of output intervals from Rlaebz
            for (j = 1; j <= iout; j = j + 1) {
                //              eigenvalue approximation is middle point of interval
                tmp1 = half * (work[(j + n) - 1] + work[(j + in + n) - 1]);
                //              semi length of error interval
                tmp2 = half * abs(work[(j + n) - 1] - work[(j + in + n) - 1]);
                if (j > iout - iinfo) {
                    //                 Flag non-convergence.
                    ncnvrg = true;
                    ib = -jblk;
                } else {
                    ib = jblk;
                }
                for (je = iwork[j - 1] + 1 + iwoff; je <= iwork[(j + in) - 1] + iwoff; je = je + 1) {
                    w[je - 1] = tmp1;
                    werr[je - 1] = tmp2;
                    indexw[je - 1] = je - iwoff;
                    iblock[je - 1] = ib;
                }
            }
            //
            m += im;
        }
    statement_70:;
    }
    //
    //     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
    //     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
    if (irange == indrng) {
        idiscl = il - 1 - nwl;
        idiscu = nwu - iu;
        //
        if (idiscl > 0) {
            im = 0;
            for (je = 1; je <= m; je = je + 1) {
                //              Remove some of the smallest eigenvalues from the left so that
                //              at the end IDISCL =0. Move all eigenvalues up to the left.
                if (w[je - 1] <= wlu && idiscl > 0) {
                    idiscl = idiscl - 1;
                } else {
                    im++;
                    w[im - 1] = w[je - 1];
                    werr[im - 1] = werr[je - 1];
                    indexw[im - 1] = indexw[je - 1];
                    iblock[im - 1] = iblock[je - 1];
                }
            }
            m = im;
        }
        if (idiscu > 0) {
            //           Remove some of the largest eigenvalues from the right so that
            //           at the end IDISCU =0. Move all eigenvalues up to the left.
            im = m + 1;
            for (je = m; je >= 1; je = je - 1) {
                if (w[je - 1] >= wul && idiscu > 0) {
                    idiscu = idiscu - 1;
                } else {
                    im = im - 1;
                    w[im - 1] = w[je - 1];
                    werr[im - 1] = werr[je - 1];
                    indexw[im - 1] = indexw[je - 1];
                    iblock[im - 1] = iblock[je - 1];
                }
            }
            jee = 0;
            for (je = im; je <= m; je = je + 1) {
                jee++;
                w[jee - 1] = w[je - 1];
                werr[jee - 1] = werr[je - 1];
                indexw[jee - 1] = indexw[je - 1];
                iblock[jee - 1] = iblock[je - 1];
            }
            m = m - im + 1;
        }
        //
        if (idiscl > 0 || idiscu > 0) {
            //           Code to deal with effects of bad arithmetic. (If N(w) is
            //           monotone non-decreasing, this should never happen.)
            //           Some low eigenvalues to be discarded are not in (WL,WLU],
            //           or high eigenvalues to be discarded are not in (WUL,WU]
            //           so just kill off the smallest IDISCL/largest IDISCU
            //           eigenvalues, by marking the corresponding IBLOCK = 0
            if (idiscl > 0) {
                wkill = wu;
                for (jdisc = 1; jdisc <= idiscl; jdisc = jdisc + 1) {
                    iw = 0;
                    for (je = 1; je <= m; je = je + 1) {
                        if (iblock[je - 1] != 0 && (w[je - 1] < wkill || iw == 0)) {
                            iw = je;
                            wkill = w[je - 1];
                        }
                    }
                    iblock[iw - 1] = 0;
                }
            }
            if (idiscu > 0) {
                wkill = wl;
                for (jdisc = 1; jdisc <= idiscu; jdisc = jdisc + 1) {
                    iw = 0;
                    for (je = 1; je <= m; je = je + 1) {
                        if (iblock[je - 1] != 0 && (w[je - 1] >= wkill || iw == 0)) {
                            iw = je;
                            wkill = w[je - 1];
                        }
                    }
                    iblock[iw - 1] = 0;
                }
            }
            //           Now erase all eigenvalues with IBLOCK set to zero
            im = 0;
            for (je = 1; je <= m; je = je + 1) {
                if (iblock[je - 1] != 0) {
                    im++;
                    w[im - 1] = w[je - 1];
                    werr[im - 1] = werr[je - 1];
                    indexw[im - 1] = indexw[je - 1];
                    iblock[im - 1] = iblock[je - 1];
                }
            }
            m = im;
        }
        if (idiscl < 0 || idiscu < 0) {
            toofew = true;
        }
    }
    //
    if ((irange == allrng && m != n) || (irange == indrng && m != iu - il + 1)) {
        toofew = true;
    }
    //
    //     If ORDER='B', do nothing the eigenvalues are already sorted by
    //        block.
    //     If ORDER='E', sort the eigenvalues from smallest to largest
    //
    if (Mlsame(order, "E") && nsplit > 1) {
        for (je = 1; je <= m - 1; je = je + 1) {
            ie = 0;
            tmp1 = w[je - 1];
            for (j = je + 1; j <= m; j = j + 1) {
                if (w[j - 1] < tmp1) {
                    ie = j;
                    tmp1 = w[j - 1];
                }
            }
            if (ie != 0) {
                tmp2 = werr[ie - 1];
                itmp1 = iblock[ie - 1];
                itmp2 = indexw[ie - 1];
                w[ie - 1] = w[je - 1];
                werr[ie - 1] = werr[je - 1];
                iblock[ie - 1] = iblock[je - 1];
                indexw[ie - 1] = indexw[je - 1];
                w[je - 1] = tmp1;
                werr[je - 1] = tmp2;
                iblock[je - 1] = itmp1;
                indexw[je - 1] = itmp2;
            }
        }
    }
    //
    info = 0;
    if (ncnvrg) {
        info++;
    }
    if (toofew) {
        info += 2;
    }
    //
    //     End of Rlarrd
    //
}
