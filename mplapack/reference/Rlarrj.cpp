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

void Rlarrj(INTEGER const n, REAL *d, REAL *e2, INTEGER const ifirst, INTEGER const ilast, REAL const rtol, INTEGER const offset, REAL *w, REAL *werr, REAL *work, INTEGER *iwork, REAL const pivmin, REAL const spdiam, INTEGER &info) {
    const REAL two = 2.0;
    INTEGER maxitr = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER nINTEGER = 0;
    INTEGER prev = 0;
    INTEGER i = 0;
    INTEGER k = 0;
    INTEGER ii = 0;
    REAL left = 0.0;
    REAL mid = 0.0;
    REAL right = 0.0;
    REAL width = 0.0;
    REAL tmp = 0.0;
    const REAL one = 1.0;
    REAL fac = 0.0;
    INTEGER cnt = 0;
    REAL s = 0.0;
    REAL dplus = 0.0;
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER savi1 = 0;
    INTEGER iter = 0;
    INTEGER olnINTEGER = 0;
    INTEGER p = 0;
    INTEGER next = 0;
    const REAL half = 0.5e0;
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
    maxitr = INTEGER((log(spdiam + pivmin) - log(pivmin)) / log(two)) + 2;
    //
    //     Initialize unconverged INTEGERervals in [ WORK(2*I-1), WORK(2*I) ].
    //     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
    //     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The INTEGEReger IWORK( 2*I-1 )
    //     for an unconverged INTEGERerval is set to the index of the next unconverged
    //     INTEGERerval, and is -1 or 0 for a converged INTEGERerval. Thus a linked
    //     list of unconverged INTEGERervals is set up.
    //
    i1 = ifirst;
    i2 = ilast;
    //     The number of unconverged INTEGERervals
    nINTEGER = 0;
    //     The last unconverged INTEGERerval found
    prev = 0;
    for (i = i1; i <= i2; i = i + 1) {
        k = 2 * i;
        ii = i - offset;
        left = w[ii - 1] - werr[ii - 1];
        mid = w[ii - 1];
        right = w[ii - 1] + werr[ii - 1];
        width = right - mid;
        tmp = max(abs(left), abs(right));
        //
        //        The following test prevents the test of converged INTEGERervals
        if (width < rtol * tmp) {
            //           This INTEGERerval has already converged and does not need refinement.
            //           (Note that the gaps might change through refining the
            //            eigenvalues, however, they can only get bigger.)
            //           Remove it from the list.
            iwork[(k - 1) - 1] = -1;
            //           Make sure that I1 always poINTEGERs to the first unconverged INTEGERerval
            if ((i == i1) && (i < i2)) {
                i1 = i + 1;
            }
            if ((prev >= i1) && (i <= i2)) {
                iwork[(2 * prev - 1) - 1] = i + 1;
            }
        } else {
            //           unconverged INTEGERerval found
            prev = i;
            //           Make sure that [LEFT,RIGHT] contains the desired eigenvalue
            //
            //           Do while( CNT(LEFT).GT.I-1 )
            //
            fac = one;
        statement_20:
            cnt = 0;
            s = left;
            dplus = d[1 - 1] - s;
            if (dplus < zero) {
                cnt++;
            }
            for (j = 2; j <= n; j = j + 1) {
                dplus = d[j - 1] - s - e2[(j - 1) - 1] / dplus;
                if (dplus < zero) {
                    cnt++;
                }
            }
            if (cnt > i - 1) {
                left = left - werr[ii - 1] * fac;
                fac = two * fac;
                goto statement_20;
            }
            //
            //           Do while( CNT(RIGHT).LT.I )
            //
            fac = one;
        statement_50:
            cnt = 0;
            s = right;
            dplus = d[1 - 1] - s;
            if (dplus < zero) {
                cnt++;
            }
            for (j = 2; j <= n; j = j + 1) {
                dplus = d[j - 1] - s - e2[(j - 1) - 1] / dplus;
                if (dplus < zero) {
                    cnt++;
                }
            }
            if (cnt < i) {
                right += werr[ii - 1] * fac;
                fac = two * fac;
                goto statement_50;
            }
            nINTEGER++;
            iwork[(k - 1) - 1] = i + 1;
            iwork[k - 1] = cnt;
        }
        work[(k - 1) - 1] = left;
        work[k - 1] = right;
    }
    //
    savi1 = i1;
    //
    //     Do while( NINT.GT.0 ), i.e. there are still unconverged INTEGERervals
    //     and while (ITER.LT.MAXITR)
    //
    iter = 0;
statement_80:
    prev = i1 - 1;
    i = i1;
    olnINTEGER = nINTEGER;
    //
    for (p = 1; p <= olnINTEGER; p = p + 1) {
        k = 2 * i;
        ii = i - offset;
        next = iwork[(k - 1) - 1];
        left = work[(k - 1) - 1];
        right = work[k - 1];
        mid = half * (left + right);
        //
        //        semiwidth of INTEGERerval
        width = right - mid;
        tmp = max(abs(left), abs(right));
        //
        if ((width < rtol * tmp) || (iter == maxitr)) {
            //           reduce number of unconverged INTEGERervals
            nINTEGER = nINTEGER - 1;
            //           Mark INTEGERerval as converged.
            iwork[(k - 1) - 1] = 0;
            if (i1 == i) {
                i1 = next;
            } else {
                //              Prev holds the last unconverged INTEGERerval previously examined
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
        cnt = 0;
        s = mid;
        dplus = d[1 - 1] - s;
        if (dplus < zero) {
            cnt++;
        }
        for (j = 2; j <= n; j = j + 1) {
            dplus = d[j - 1] - s - e2[(j - 1) - 1] / dplus;
            if (dplus < zero) {
                cnt++;
            }
        }
        if (cnt <= i - 1) {
            work[(k - 1) - 1] = mid;
        } else {
            work[k - 1] = mid;
        }
        i = next;
    //
    statement_100:;
    }
    iter++;
    //     do another loop if there are still unconverged INTEGERervals
    //     However, in the last iteration, all INTEGERervals are accepted
    //     since this is the best we can do.
    if ((nINTEGER > 0) && (iter <= maxitr)) {
        goto statement_80;
    }
    //
    //     At this poINTEGER, all the INTEGERervals have converged
    for (i = savi1; i <= ilast; i = i + 1) {
        k = 2 * i;
        ii = i - offset;
        //        All INTEGERervals marked by '0' have been refined.
        if (iwork[(k - 1) - 1] == 0) {
            w[ii - 1] = half * (work[(k - 1) - 1] + work[k - 1]);
            werr[ii - 1] = work[k - 1] - w[ii - 1];
        }
    }
    //
    //     End of Rlarrj
    //
}
