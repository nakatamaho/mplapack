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

void Rlarrc(const char *jobt, INTEGER const n, REAL const vl, REAL const vu, REAL *d, REAL *e, REAL const  /* pivmin */, INTEGER &eigcnt, INTEGER &lcnt, INTEGER &rcnt, INTEGER &info) {
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
    //     .. External Functions ..
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
    lcnt = 0;
    rcnt = 0;
    eigcnt = 0;
    bool matt = Mlsame(jobt, "T");
    //
    REAL lpivot = 0.0;
    REAL rpivot = 0.0;
    const REAL zero = 0.0;
    INTEGER i = 0;
    REAL tmp = 0.0;
    REAL sl = 0.0;
    REAL su = 0.0;
    REAL tmp2 = 0.0;
    if (matt) {
        //        Sturm sequence count on T
        lpivot = d[1 - 1] - vl;
        rpivot = d[1 - 1] - vu;
        if (lpivot <= zero) {
            lcnt++;
        }
        if (rpivot <= zero) {
            rcnt++;
        }
        for (i = 1; i <= n - 1; i = i + 1) {
            tmp = pow2(e[i - 1]);
            lpivot = (d[(i + 1) - 1] - vl) - tmp / lpivot;
            rpivot = (d[(i + 1) - 1] - vu) - tmp / rpivot;
            if (lpivot <= zero) {
                lcnt++;
            }
            if (rpivot <= zero) {
                rcnt++;
            }
        }
    } else {
        //        Sturm sequence count on L D L^T
        sl = -vl;
        su = -vu;
        for (i = 1; i <= n - 1; i = i + 1) {
            lpivot = d[i - 1] + sl;
            rpivot = d[i - 1] + su;
            if (lpivot <= zero) {
                lcnt++;
            }
            if (rpivot <= zero) {
                rcnt++;
            }
            tmp = e[i - 1] * d[i - 1] * e[i - 1];
            //
            tmp2 = tmp / lpivot;
            if (tmp2 == zero) {
                sl = tmp - vl;
            } else {
                sl = sl * tmp2 - vl;
            }
            //
            tmp2 = tmp / rpivot;
            if (tmp2 == zero) {
                su = tmp - vu;
            } else {
                su = su * tmp2 - vu;
            }
        }
        lpivot = d[n - 1] + sl;
        rpivot = d[n - 1] + su;
        if (lpivot <= zero) {
            lcnt++;
        }
        if (rpivot <= zero) {
            rcnt++;
        }
    }
    eigcnt = rcnt - lcnt;
    //
    //     end of Rlarrc
    //
}
