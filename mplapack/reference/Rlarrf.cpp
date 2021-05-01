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

void Rlarrf(INTEGER const n, REAL *d, REAL *l, REAL *ld, INTEGER const clstrt, INTEGER const clend, REAL *w, REAL *wgap, REAL *werr, REAL const spdiam, REAL const clgapl, REAL const clgapr, REAL const pivmin, REAL &sigma, REAL *dplus, REAL *lplus, REAL *work, INTEGER &info) {
    const INTEGER ktrymax = 1;
    REAL fact = 0.0;
    REAL eps = 0.0;
    INTEGER shift = 0;
    bool forcer = false;
    bool nofail = false;
    REAL clwdth = 0.0;
    REAL avgap = 0.0;
    REAL mingap = 0.0;
    REAL lsigma = 0.0;
    REAL rsigma = 0.0;
    const REAL four = 4.0;
    const REAL quart = 0.25e0;
    const REAL two = 2.0;
    REAL ldmax = 0.0;
    REAL rdmax = 0.0;
    REAL ldelta = 0.0;
    REAL rdelta = 0.0;
    REAL s = 0.0;
    const REAL one = 1.0;
    REAL smlgrowth = 0.0;
    REAL fail = 0.0;
    REAL fail2 = 0.0;
    REAL bestshift = 0.0;
    INTEGER ktry = 0;
    const REAL maxgrowth1 = 8.e0;
    REAL growthbound = 0.0;
    bool sawnan1 = false;
    bool sawnan2 = false;
    REAL max1 = 0.0;
    INTEGER i = 0;
    const INTEGER sleft = 1;
    REAL max2 = 0.0;
    const INTEGER sright = 2;
    INTEGER indx = 0;
    bool dorrr1 = false;
    bool tryrrr1 = false;
    REAL tmp = 0.0;
    REAL znm2 = 0.0;
    REAL prod = 0.0;
    REAL oldp = 0.0;
    REAL rrr1 = 0.0;
    const REAL maxgrowth2 = 8.e0;
    REAL rrr2 = 0.0;
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
    //  =====================================================================
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
    //     .. Executable Statements ..
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    fact = castREAL(pow(2, ktrymax));
    eps = Rlamch("Precision");
    shift = 0;
    forcer = false;
    //
    //     Note that we cannot guarantee that for any of the shifts tried,
    //     the factorization has a small or even moderate element growth.
    //     There could be Ritz values at both ends of the cluster and despite
    //     backing off, there are examples where all factorizations tried
    //     (in IEEE mode, allowing zero pivots & infinities) have INFINITE
    //     element growth.
    //     For this reason, we should use PIVMIN in this subroutine so that at
    //     least the L D L^T factorization exists. It can be checked afterwards
    //     whether the element growth caused bad residuals/orthogonality.
    //
    //     Decide whether the code should accept the best among all
    //     representations despite large element growth or signal INFO=1
    //     Setting NOFAIL to .FALSE. for quick fix for bug 113
    nofail = false;
    //
    //     Compute the average gap length of the cluster
    clwdth = abs(w[clend - 1] - w[clstrt - 1]) + werr[clend - 1] + werr[clstrt - 1];
    avgap = clwdth / castREAL(clend - clstrt);
    mingap = min(clgapl, clgapr);
    //     Initial values for shifts to both ends of cluster
    lsigma = min(w[clstrt - 1], w[clend - 1]) - werr[clstrt - 1];
    rsigma = max(w[clstrt - 1], w[clend - 1]) + werr[clend - 1];
    //
    //     Use a small fudge to make sure that we really shift to the outside
    lsigma = lsigma - abs(lsigma) * four * eps;
    rsigma += abs(rsigma) * four * eps;
    //
    //     Compute upper bounds for how much to back off the initial shifts
    ldmax = quart * mingap + two * pivmin;
    rdmax = quart * mingap + two * pivmin;
    //
    ldelta = max(avgap, wgap[clstrt - 1]) / fact;
    rdelta = max(avgap, wgap[(clend - 1) - 1]) / fact;
    //
    //     Initialize the record of the best representation found
    //
    s = Rlamch("S");
    smlgrowth = one / s;
    fail = castREAL(n - 1) * mingap / (spdiam * eps);
    fail2 = castREAL(n - 1) * mingap / (spdiam * sqrt(eps));
    bestshift = lsigma;
    //
    //     while (KTRY <= KTRYMAX)
    ktry = 0;
    growthbound = maxgrowth1 * spdiam;
//
statement_5:
    sawnan1 = false;
    sawnan2 = false;
    //     Ensure that we do not back off too much of the initial shifts
    ldelta = min(ldmax, ldelta);
    rdelta = min(rdmax, rdelta);
    //
    //     Compute the element growth when shifting to both ends of the cluster
    //     accept the shift if there is no element growth at one of the two ends
    //
    //     Left end
    s = -lsigma;
    dplus[1 - 1] = d[1 - 1] + s;
    if (abs(dplus[1 - 1]) < pivmin) {
        dplus[1 - 1] = -pivmin;
        //        Need to set SAWNAN1 because refined RRR test should not be used
        //        in this case
        sawnan1 = true;
    }
    max1 = abs(dplus[1 - 1]);
    for (i = 1; i <= n - 1; i = i + 1) {
        lplus[i - 1] = ld[i - 1] / dplus[i - 1];
        s = s * lplus[i - 1] * l[i - 1] - lsigma;
        dplus[(i + 1) - 1] = d[(i + 1) - 1] + s;
        if (abs(dplus[(i + 1) - 1]) < pivmin) {
            dplus[(i + 1) - 1] = -pivmin;
            //           Need to set SAWNAN1 because refined RRR test should not be used
            //           in this case
            sawnan1 = true;
        }
        max1 = max(max1, abs(dplus[(i + 1) - 1]));
    }
    sawnan1 = sawnan1 || Risnan(max1);
    //
    if (forcer || (max1 <= growthbound && !sawnan1)) {
        sigma = lsigma;
        shift = sleft;
        goto statement_100;
    }
    //
    //     Right end
    s = -rsigma;
    work[1 - 1] = d[1 - 1] + s;
    if (abs(work[1 - 1]) < pivmin) {
        work[1 - 1] = -pivmin;
        //        Need to set SAWNAN2 because refined RRR test should not be used
        //        in this case
        sawnan2 = true;
    }
    max2 = abs(work[1 - 1]);
    for (i = 1; i <= n - 1; i = i + 1) {
        work[(n + i) - 1] = ld[i - 1] / work[i - 1];
        s = s * work[(n + i) - 1] * l[i - 1] - rsigma;
        work[(i + 1) - 1] = d[(i + 1) - 1] + s;
        if (abs(work[(i + 1) - 1]) < pivmin) {
            work[(i + 1) - 1] = -pivmin;
            //           Need to set SAWNAN2 because refined RRR test should not be used
            //           in this case
            sawnan2 = true;
        }
        max2 = max(max2, abs(work[(i + 1) - 1]));
    }
    sawnan2 = sawnan2 || Risnan(max2);
    //
    if (forcer || (max2 <= growthbound && !sawnan2)) {
        sigma = rsigma;
        shift = sright;
        goto statement_100;
    }
    //     If we are at this point, both shifts led to too much element growth
    //
    //     Record the better of the two shifts (provided it didn't lead to NaN)
    if (sawnan1 && sawnan2) {
        //        both MAX1 and MAX2 are NaN
        goto statement_50;
    } else {
        if (!sawnan1) {
            indx = 1;
            if (max1 <= smlgrowth) {
                smlgrowth = max1;
                bestshift = lsigma;
            }
        }
        if (!sawnan2) {
            if (sawnan1 || max2 <= max1) {
                indx = 2;
            }
            if (max2 <= smlgrowth) {
                smlgrowth = max2;
                bestshift = rsigma;
            }
        }
    }
    //
    //     If we are here, both the left and the right shift led to
    //     element growth. If the element growth is moderate, then
    //     we may still accept the representation, if it passes a
    //     refined test for RRR. This test supposes that no NaN occurred.
    //     Moreover, we use the refined RRR test only for isolated clusters.
    if ((clwdth < mingap / 128.0) && (min(max1, max2) < fail2) && (!sawnan1) && (!sawnan2)) {
        dorrr1 = true;
    } else {
        dorrr1 = false;
    }
    tryrrr1 = true;
    if (tryrrr1 && dorrr1) {
        if (indx == 1) {
            tmp = abs(dplus[n - 1]);
            znm2 = one;
            prod = one;
            oldp = one;
            for (i = n - 1; i >= 1; i = i - 1) {
                if (prod <= eps) {
                    prod = ((dplus[(i + 1) - 1] * work[(n + i + 1) - 1]) / (dplus[i - 1] * work[(n + i) - 1])) * oldp;
                } else {
                    prod = prod * abs(work[(n + i) - 1]);
                }
                oldp = prod;
                znm2 += pow2(prod);
                tmp = max(tmp, abs(dplus[i - 1] * prod));
            }
            rrr1 = tmp / (spdiam * sqrt(znm2));
            if (rrr1 <= maxgrowth2) {
                sigma = lsigma;
                shift = sleft;
                goto statement_100;
            }
        } else if (indx == 2) {
            tmp = abs(work[n - 1]);
            znm2 = one;
            prod = one;
            oldp = one;
            for (i = n - 1; i >= 1; i = i - 1) {
                if (prod <= eps) {
                    prod = ((work[(i + 1) - 1] * lplus[(i + 1) - 1]) / (work[i - 1] * lplus[i - 1])) * oldp;
                } else {
                    prod = prod * abs(lplus[i - 1]);
                }
                oldp = prod;
                znm2 += pow2(prod);
                tmp = max(tmp, abs(work[i - 1] * prod));
            }
            rrr2 = tmp / (spdiam * sqrt(znm2));
            if (rrr2 <= maxgrowth2) {
                sigma = rsigma;
                shift = sright;
                goto statement_100;
            }
        }
    }
//
statement_50:
    //
    if (ktry < ktrymax) {
        //        If we are here, both shifts failed also the RRR test.
        //        Back off to the outside
        lsigma = max(lsigma - ldelta, lsigma - ldmax);
        rsigma = min(rsigma + rdelta, rsigma + rdmax);
        ldelta = two * ldelta;
        rdelta = two * rdelta;
        ktry++;
        goto statement_5;
    } else {
        //        None of the representations investigated satisfied our
        //        criteria. Take the best one we found.
        if ((smlgrowth < fail) || nofail) {
            lsigma = bestshift;
            rsigma = bestshift;
            forcer = true;
            goto statement_5;
        } else {
            info = 1;
            return;
        }
    }
//
statement_100:
    if (shift == sleft) {
    } else if (shift == sright) {
        //        store new L and D back into DPLUS, LPLUS
        Rcopy(n, work, 1, dplus, 1);
        Rcopy(n - 1, &work[(n + 1) - 1], 1, lplus, 1);
    }
    //
    //     End of Rlarrf
    //
}
