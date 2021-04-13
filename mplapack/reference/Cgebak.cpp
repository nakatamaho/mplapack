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

void Cgebak(const char *job, const char *side, INTEGER const n, INTEGER const ilo, INTEGER const ihi, REAL *scale, INTEGER const m, COMPLEX *v, INTEGER const ldv, INTEGER &info) {
    bool rightv = false;
    bool leftv = false;
    INTEGER i = 0;
    REAL s = 0.0;
    const REAL one = 1.0;
    INTEGER ii = 0;
    INTEGER k = 0;
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
    //     Decode and Test the input parameters
    //
    rightv = Mlsame(side, "R");
    leftv = Mlsame(side, "L");
    //
    info = 0;
    if (!Mlsame(job, "N") && !Mlsame(job, "P") && !Mlsame(job, "S") && !Mlsame(job, "B")) {
        info = -1;
    } else if (!rightv && !leftv) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (ilo < 1 || ilo > max((INTEGER)1, n)) {
        info = -4;
    } else if (ihi < min(ilo, n) || ihi > n) {
        info = -5;
    } else if (m < 0) {
        info = -7;
    } else if (ldv < max((INTEGER)1, n)) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Cgebak", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    if (m == 0) {
        return;
    }
    if (Mlsame(job, "N")) {
        return;
    }
    //
    if (ilo == ihi) {
        goto statement_30;
    }
    //
    //     Backward balance
    //
    if (Mlsame(job, "S") || Mlsame(job, "B")) {
        //
        if (rightv) {
            for (i = ilo; i <= ihi; i = i + 1) {
                s = scale[i - 1];
                CRscal(m, s, &v[(i - 1)], ldv);
            }
        }
        //
        if (leftv) {
            for (i = ilo; i <= ihi; i = i + 1) {
                s = one / scale[i - 1];
                CRscal(m, s, &v[(i - 1)], ldv);
            }
        }
        //
    }
//
//     Backward permutation
//
//     For  I = ILO-1 step -1 until 1,
//              IHI+1 step 1 until N do --
//
statement_30:
    if (Mlsame(job, "P") || Mlsame(job, "B")) {
        if (rightv) {
            for (ii = 1; ii <= n; ii = ii + 1) {
                i = ii;
                if (i >= ilo && i <= ihi) {
                    goto statement_40;
                }
                if (i < ilo) {
                    i = ilo - ii;
                }
                k = castINTEGER(scale[i - 1]);
                if (k == i) {
                    goto statement_40;
                }
                Cswap(m, &v[(i - 1)], ldv, &v[(k - 1)], ldv);
            statement_40:;
            }
        }
        //
        if (leftv) {
            for (ii = 1; ii <= n; ii = ii + 1) {
                i = ii;
                if (i >= ilo && i <= ihi) {
                    goto statement_50;
                }
                if (i < ilo) {
                    i = ilo - ii;
                }
                k = castINTEGER(scale[i - 1]);
                if (k == i) {
                    goto statement_50;
                }
                Cswap(m, &v[(i - 1)], ldv, &v[(k - 1)], ldv);
            statement_50:;
            }
        }
    }
    //
    //     End of Cgebak
    //
}
