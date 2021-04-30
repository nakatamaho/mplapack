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

void Rlarrb(INTEGER const n, REAL *d, REAL *lld, INTEGER const ifirst, INTEGER const ilast, REAL const rtol1, REAL const rtol2, INTEGER const offset, REAL *w, REAL *wgap, REAL *werr, REAL *work, INTEGER *iwork, REAL const pivmin, REAL const spdiam, INTEGER const twist, INTEGER &info) {
    const REAL two = 2.0;
    INTEGER maxitr = 0;
    REAL mnwdth = 0.0;
    INTEGER r = 0;
    INTEGER i1 = 0;
    INTEGER nint = 0;
    INTEGER prev = 0;
    REAL rgap = 0.0;
    INTEGER i = 0;
    INTEGER k = 0;
    INTEGER ii = 0;
    REAL left = 0.0;
    REAL right = 0.0;
    REAL lgap = 0.0;
    REAL gap = 0.0;
    REAL back = 0.0;
    INTEGER negcnt = 0;
    const REAL half = 0.5e0;
    REAL width = 0.0;
    REAL tmp = 0.0;
    REAL cvrgd = 0.0;
    INTEGER iter = 0;
    INTEGER olnint = 0;
    INTEGER ip = 0;
    INTEGER next = 0;
    REAL mid = 0.0;
    const REAL zero = 0.0;
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
    //
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
    maxitr = castINTEGER((log(spdiam + pivmin) - log(pivmin)) / log(two)) + (INTEGER)2;
    mnwdth = two * pivmin;
    //
    r = twist;
    if ((r < 1) || (r > n)) {
        r = n;
    }
    //
    //     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
    //     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
    //     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )
    //     for an unconverged interval is set to the index of the next unconverged
    //     interval, and is -1 or 0 for a converged interval. Thus a linked
    //     list of unconverged intervals is set up.
    //
    i1 = ifirst;
    //     The number of unconverged intervals
    nint = 0;
    //     The last unconverged interval found
    prev = 0;
    //
    rgap = wgap[(i1 - offset) - 1];
    for (i = i1; i <= ilast; i = i + 1) {
        k = 2 * i;
        ii = i - offset;
        left = w[ii - 1] - werr[ii - 1];
        right = w[ii - 1] + werr[ii - 1];
        lgap = rgap;
        rgap = wgap[ii - 1];
        gap = min(lgap, rgap);
        //
        //        Make sure that [LEFT,RIGHT] contains the desired eigenvalue
        //        Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT
        //
        //        Do while( NEGCNT(LEFT).GT.I-1 )
        //
        back = werr[ii - 1];
    statement_20:
        negcnt = Rlaneg(n, d, lld, left, pivmin, r);
        if (negcnt > i - 1) {
            left = left - back;
            back = two * back;
            goto statement_20;
        }
        //
        //        Do while( NEGCNT(RIGHT).LT.I )
        //        Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT
        //
        back = werr[ii - 1];
    statement_50:
        //
        negcnt = Rlaneg(n, d, lld, right, pivmin, r);
        if (negcnt < i) {
            right += back;
            back = two * back;
            goto statement_50;
        }
        width = half * abs(left - right);
        tmp = max(abs(left), abs(right));
        cvrgd = max(rtol1 * gap, rtol2 * tmp);
        if (width <= cvrgd || width <= mnwdth) {
            //           This interval has already converged and does not need refinement.
            //           (Note that the gaps might change through refining the
            //            eigenvalues, however, they can only get bigger.)
            //           Remove it from the list.
            iwork[(k - 1) - 1] = -1;
            //           Make sure that I1 always points to the first unconverged interval
            if ((i == i1) && (i < ilast)) {
                i1 = i + 1;
            }
            if ((prev >= i1) && (i <= ilast)) {
                iwork[(2 * prev - 1) - 1] = i + 1;
            }
        } else {
            //           unconverged interval found
            prev = i;
            nint++;
            iwork[(k - 1) - 1] = i + 1;
            iwork[k - 1] = negcnt;
        }
        work[(k - 1) - 1] = left;
        work[k - 1] = right;
    }
    //
    //     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals
    //     and while (ITER.LT.MAXITR)
    //
    iter = 0;
statement_80:
    prev = i1 - 1;
    i = i1;
    olnint = nint;
    //
    for (ip = 1; ip <= olnint; ip = ip + 1) {
        k = 2 * i;
        ii = i - offset;
        rgap = wgap[ii - 1];
        lgap = rgap;
        if (ii > 1) {
            lgap = wgap[(ii - 1) - 1];
        }
        gap = min(lgap, rgap);
        next = iwork[(k - 1) - 1];
        left = work[(k - 1) - 1];
        right = work[k - 1];
        mid = half * (left + right);
        //
        //        semiwidth of interval
        width = right - mid;
        tmp = max(abs(left), abs(right));
        cvrgd = max(rtol1 * gap, rtol2 * tmp);
        if ((width <= cvrgd) || (width <= mnwdth) || (iter == maxitr)) {
            //           reduce number of unconverged intervals
            nint = nint - 1;
            //           Mark interval as converged.
            iwork[(k - 1) - 1] = 0;
            if (i1 == i) {
                i1 = next;
            } else {
                //              Prev holds the last unconverged interval previously examined
                if (prev >= i1) {
                    iwork[(2 * prev - 1) - 1] = next;
                }
            }
            i = next;
            goto statement_100;
        }
        prev = i;
        //
        //        Perform one bisection step
        //
        negcnt = Rlaneg(n, d, lld, mid, pivmin, r);
        if (negcnt <= i - 1) {
            work[(k - 1) - 1] = mid;
        } else {
            work[k - 1] = mid;
        }
        i = next;
    statement_100:;
    }
    iter++;
    //     do another loop if there are still unconverged intervals
    //     However, in the last iteration, all intervals are accepted
    //     since this is the best we can do.
    if ((nint > 0) && (iter <= maxitr)) {
        goto statement_80;
    }
    //
    //     At this point, all the intervals have converged
    for (i = ifirst; i <= ilast; i = i + 1) {
        k = 2 * i;
        ii = i - offset;
        //        All intervals marked by '0' have been refined.
        if (iwork[(k - 1) - 1] == 0) {
            w[ii - 1] = half * (work[(k - 1) - 1] + work[k - 1]);
            werr[ii - 1] = work[k - 1] - w[ii - 1];
        }
    }
    //
    for (i = ifirst + 1; i <= ilast; i = i + 1) {
        k = 2 * i;
        ii = i - offset;
        wgap[(ii - 1) - 1] = max(zero, w[ii - 1] - werr[ii - 1] - w[(ii - 1) - 1] - werr[(ii - 1) - 1]);
    }
    //
    //     End of Rlarrb
    //
}
