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

void Rstebz(const char *range, const char *order, INTEGER const n, REAL const vl, REAL const vu, INTEGER const il, INTEGER const iu, REAL const abstol, REAL *d, REAL *e, INTEGER &m, INTEGER &nsplit, REAL *w, INTEGER *iblock, INTEGER *isplit, REAL *work, INTEGER *iwork, INTEGER &info) {
    INTEGER irange = 0;
    INTEGER iorder = 0;
    bool ncnvrg = false;
    bool toofew = false;
    REAL safemn = 0.0;
    REAL ulp = 0.0;
    const REAL relfac = 2.0;
    REAL rtoli = 0.0;
    INTEGER nb = 0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL pivmin = 0.0;
    INTEGER j = 0;
    REAL tmp1 = 0.0;
    REAL gu = 0.0;
    REAL gl = 0.0;
    REAL tmp2 = 0.0;
    REAL tnorm = 0.0;
    const REAL fudge = 2.1e0;
    const REAL two = 2.0;
    INTEGER itmax = 0;
    REAL atoli = 0.0;
    INTEGER iout = 0;
    INTEGER iinfo = 0;
    REAL wl = 0.0;
    REAL wlu = 0.0;
    INTEGER nwl = 0;
    REAL wu = 0.0;
    REAL wul = 0.0;
    INTEGER nwu = 0;
    INTEGER iend = 0;
    INTEGER jb = 0;
    INTEGER ioff = 0;
    INTEGER ibegin = 0;
    INTEGER in = 0;
    REAL bnorm = 0.0;
    INTEGER idumma[1];
    INTEGER im = 0;
    INTEGER iwoff = 0;
    const REAL half = 1.0 / two;
    INTEGER ib = 0;
    INTEGER je = 0;
    INTEGER idiscl = 0;
    INTEGER idiscu = 0;
    REAL wkill = 0.0;
    INTEGER jdisc = 0;
    INTEGER iw = 0;
    INTEGER ie = 0;
    INTEGER itmp1 = 0;
    //
    info = 0;
    //
    //     Decode RANGE
    //
    if (Mlsame(range, "A")) {
        irange = 1;
    } else if (Mlsame(range, "V")) {
        irange = 2;
    } else if (Mlsame(range, "I")) {
        irange = 3;
    } else {
        irange = 0;
    }
    //
    //     Decode ORDER
    //
    if (Mlsame(order, "B")) {
        iorder = 2;
    } else if (Mlsame(order, "E")) {
        iorder = 1;
    } else {
        iorder = 0;
    }
    //
    //     Check for Errors
    //
    if (irange <= 0) {
        info = -1;
    } else if (iorder <= 0) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (irange == 2) {
        if (vl >= vu) {
            info = -5;
        }
    } else if (irange == 3 && (il < 1 || il > max((INTEGER)1, n))) {
        info = -6;
    } else if (irange == 3 && (iu < min(n, il) || iu > n)) {
        info = -7;
    }
    //
    if (info != 0) {
        Mxerbla("Rstebz", -info);
        return;
    }
    //
    //     Initialize error flags
    //
    info = 0;
    ncnvrg = false;
    toofew = false;
    //
    //     Quick return if possible
    //
    m = 0;
    if (n == 0) {
        return;
    }
    //
    //     Simplifications:
    //
    if (irange == 3 && il == 1 && iu == n) {
        irange = 1;
    }
    //
    //     Get machine constants
    //     NB is the minimum vector length for vector bisection, or 0
    //     if only scalar is to be done.
    //
    safemn = Rlamch("S");
    ulp = Rlamch("P");
    rtoli = ulp * relfac;
    nb = iMlaenv(1, "Rstebz", " ", n, -1, -1, -1);
    if (nb <= 1) {
        nb = 0;
    }
    //
    //     Special Case when N=1
    //
    if (n == 1) {
        nsplit = 1;
        isplit[1 - 1] = 1;
        if (irange == 2 && (vl >= d[1 - 1] || vu < d[1 - 1])) {
            m = 0;
        } else {
            w[1 - 1] = d[1 - 1];
            iblock[1 - 1] = 1;
            m = 1;
        }
        return;
    }
    //
    //     Compute Splitting Points
    //
    nsplit = 1;
    work[n - 1] = zero;
    pivmin = one;
    //
    for (j = 2; j <= n; j = j + 1) {
        tmp1 = pow2(e[(j - 1) - 1]);
        if (abs(d[j - 1] * d[(j - 1) - 1]) * pow2(ulp) + safemn > tmp1) {
            isplit[nsplit - 1] = j - 1;
            nsplit++;
            work[(j - 1) - 1] = zero;
        } else {
            work[(j - 1) - 1] = tmp1;
            pivmin = max(pivmin, tmp1);
        }
    }
    isplit[nsplit - 1] = n;
    pivmin = pivmin * safemn;
    //
    //     Compute Interval and ATOLI
    //
    if (irange == 3) {
        //
        //        RANGE='I': Compute the interval containing eigenvalues
        //                   IL through IU.
        //
        //        Compute Gershgorin interval for entire (split) matrix
        //        and use it as the initial interval
        //
        gu = d[1 - 1];
        gl = d[1 - 1];
        tmp1 = zero;
        //
        for (j = 1; j <= n - 1; j = j + 1) {
            tmp2 = sqrt(work[j - 1]);
            gu = max(gu, d[j - 1] + tmp1 + tmp2);
            gl = min(gl, d[j - 1] - tmp1 - tmp2);
            tmp1 = tmp2;
        }
        //
        gu = max(gu, d[n - 1] + tmp1);
        gl = min(gl, d[n - 1] - tmp1);
        tnorm = max(abs(gl), abs(gu));
        gl = gl - fudge * tnorm * ulp * n - fudge * two * pivmin;
        gu += fudge * tnorm * ulp * n + fudge * pivmin;
        //
        //        Compute Iteration parameters
        //
        itmax = castINTEGER((log(tnorm + pivmin) - log(pivmin)) / log(two)) + 2;
        if (itmax >= 1024)
            itmax = 1024; // XXX itmax can be too large for MPFR (=10^8)

        if (abstol <= zero) {
            atoli = ulp * tnorm;
        } else {
            atoli = abstol;
        }
        //
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
        Rlaebz(3, itmax, n, 2, 2, nb, atoli, rtoli, pivmin, d, e, work, &iwork[5 - 1], &work[(n + 1) - 1], &work[(n + 5) - 1], iout, iwork, w, iblock, iinfo);
        //
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
        //
        if (nwl < 0 || nwl >= n || nwu < 1 || nwu > n) {
            info = 4;
            return;
        }
    } else {
        //
        //        RANGE='A' or 'V' -- Set ATOLI
        //
        tnorm = max(abs(d[1 - 1]) + abs(e[1 - 1]), abs(d[n - 1]) + abs(e[(n - 1) - 1]));
        //
        for (j = 2; j <= n - 1; j = j + 1) {
            tnorm = max(tnorm, abs(d[j - 1]) + abs(e[(j - 1) - 1]) + abs(e[j - 1]));
        }
        //
        if (abstol <= zero) {
            atoli = ulp * tnorm;
        } else {
            atoli = abstol;
        }
        //
        if (irange == 2) {
            wl = vl;
            wu = vu;
        } else {
            wl = zero;
            wu = zero;
        }
    }
    //
    //     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
    //     NWL accumulates the number of eigenvalues .le. WL,
    //     NWU accumulates the number of eigenvalues .le. WU
    //
    m = 0;
    iend = 0;
    info = 0;
    nwl = 0;
    nwu = 0;
    //
    for (jb = 1; jb <= nsplit; jb = jb + 1) {
        ioff = iend;
        ibegin = ioff + 1;
        iend = isplit[jb - 1];
        in = iend - ioff;
        //
        if (in == 1) {
            //
            //           Special Case -- IN=1
            //
            if (irange == 1 || wl >= d[ibegin - 1] - pivmin) {
                nwl++;
            }
            if (irange == 1 || wu >= d[ibegin - 1] - pivmin) {
                nwu++;
            }
            if (irange == 1 || (wl < d[ibegin - 1] - pivmin && wu >= d[ibegin - 1] - pivmin)) {
                m++;
                w[m - 1] = d[ibegin - 1];
                iblock[m - 1] = jb;
            }
        } else {
            //
            //           General Case -- IN > 1
            //
            //           Compute Gershgorin Interval
            //           and use it as the initial interval
            //
            gu = d[ibegin - 1];
            gl = d[ibegin - 1];
            tmp1 = zero;
            //
            for (j = ibegin; j <= iend - 1; j = j + 1) {
                tmp2 = abs(e[j - 1]);
                gu = max(gu, d[j - 1] + tmp1 + tmp2);
                gl = min(gl, d[j - 1] - tmp1 - tmp2);
                tmp1 = tmp2;
            }
            //
            gu = max(gu, d[iend - 1] + tmp1);
            gl = min(gl, d[iend - 1] - tmp1);
            bnorm = max(abs(gl), abs(gu));
            gl = gl - fudge * bnorm * ulp * in - fudge * pivmin;
            gu += fudge * bnorm * ulp * in + fudge * pivmin;
            //
            //           Compute ATOLI for the current submatrix
            //
            if (abstol <= zero) {
                atoli = ulp * max(abs(gl), abs(gu));
            } else {
                atoli = abstol;
            }
            //
            if (irange > 1) {
                if (gu < wl) {
                    nwl += in;
                    nwu += in;
                    goto statement_70;
                }
                gl = max(gl, wl);
                gu = min(gu, wu);
                if (gl >= gu) {
                    goto statement_70;
                }
            }
            //
            //           Set Up Initial Interval
            //
            work[(n + 1) - 1] = gl;
            work[(n + in + 1) - 1] = gu;
            Rlaebz(1, 0, in, in, 1, nb, atoli, rtoli, pivmin, &d[ibegin - 1], &e[ibegin - 1], &work[ibegin - 1], idumma, &work[(n + 1) - 1], &work[(n + 2 * in + 1) - 1], im, iwork, &w[(m + 1) - 1], &iblock[(m + 1) - 1], iinfo);
            //
            nwl += iwork[1 - 1];
            nwu += iwork[(in + 1) - 1];
            iwoff = m - iwork[1 - 1];
            //
            //           Compute Eigenvalues
            //
            itmax = castINTEGER((log(gu - gl + pivmin) - log(pivmin)) / log(two)) + 2;
            if (itmax >= 1024)
                itmax = 1024; // XXX itmax can be too large for MPFR (=10^8)
            Rlaebz(2, itmax, in, in, 1, nb, atoli, rtoli, pivmin, &d[ibegin - 1], &e[ibegin - 1], &work[ibegin - 1], idumma, &work[(n + 1) - 1], &work[(n + 2 * in + 1) - 1], iout, iwork, &w[(m + 1) - 1], &iblock[(m + 1) - 1], iinfo);
            //
            //           Copy Eigenvalues Into W and IBLOCK
            //           Use -JB for block number for unconverged eigenvalues.
            //
            for (j = 1; j <= iout; j = j + 1) {
                tmp1 = half * (work[(j + n) - 1] + work[(j + in + n) - 1]);
                //
                //              Flag non-convergence.
                //
                if (j > iout - iinfo) {
                    ncnvrg = true;
                    ib = -jb;
                } else {
                    ib = jb;
                }
                for (je = iwork[j - 1] + 1 + iwoff; je <= iwork[(j + in) - 1] + iwoff; je = je + 1) {
                    w[je - 1] = tmp1;
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
    //
    if (irange == 3) {
        im = 0;
        idiscl = il - 1 - nwl;
        idiscu = nwu - iu;
        //
        if (idiscl > 0 || idiscu > 0) {
            for (je = 1; je <= m; je = je + 1) {
                if (w[je - 1] <= wlu && idiscl > 0) {
                    idiscl = idiscl - 1;
                } else if (w[je - 1] >= wul && idiscu > 0) {
                    idiscu = idiscu - 1;
                } else {
                    im++;
                    w[im - 1] = w[je - 1];
                    iblock[im - 1] = iblock[je - 1];
                }
            }
            m = im;
        }
        if (idiscl > 0 || idiscu > 0) {
            //
            //           Code to deal with effects of bad arithmetic:
            //           Some low eigenvalues to be discarded are not in (WL,WLU],
            //           or high eigenvalues to be discarded are not in (WUL,WU]
            //           so just kill off the smallest IDISCL/largest IDISCU
            //           eigenvalues, by simply finding the smallest/largest
            //           eigenvalue(s).
            //
            //           (If N(w) is monotone non-decreasing, this should never
            //               happen.)
            //
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
                //
                wkill = wl;
                for (jdisc = 1; jdisc <= idiscu; jdisc = jdisc + 1) {
                    iw = 0;
                    for (je = 1; je <= m; je = je + 1) {
                        if (iblock[je - 1] != 0 && (w[je - 1] > wkill || iw == 0)) {
                            iw = je;
                            wkill = w[je - 1];
                        }
                    }
                    iblock[iw - 1] = 0;
                }
            }
            im = 0;
            for (je = 1; je <= m; je = je + 1) {
                if (iblock[je - 1] != 0) {
                    im++;
                    w[im - 1] = w[je - 1];
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
    //     If ORDER='B', do nothing -- the eigenvalues are already sorted
    //        by block.
    //     If ORDER='E', sort the eigenvalues from smallest to largest
    //
    if (iorder == 1 && nsplit > 1) {
        for (je = 1; je <= m - 1; je = je + 1) {
            ie = 0;
            tmp1 = w[je - 1];
            for (j = je + 1; j <= m; j = j + 1) {
                if (w[j - 1] < tmp1) {
                    ie = j;
                    tmp1 = w[j - 1];
                }
            }
            //
            if (ie != 0) {
                itmp1 = iblock[ie - 1];
                w[ie - 1] = w[je - 1];
                iblock[ie - 1] = iblock[je - 1];
                w[je - 1] = tmp1;
                iblock[je - 1] = itmp1;
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
    //     End of Rstebz
    //
}
