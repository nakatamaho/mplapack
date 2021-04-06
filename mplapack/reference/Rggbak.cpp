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

void Rggbak(const char *job, const char *side, INTEGER const n, INTEGER const ilo, INTEGER const ihi, REAL *lscale, REAL *rscale, INTEGER const m, REAL *v, INTEGER const ldv, INTEGER &info) {
    bool rightv = false;
    bool leftv = false;
    INTEGER i = 0;
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
    //     Test the input parameters
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
    } else if (ilo < 1) {
        info = -4;
    } else if (n == 0 && ihi == 0 && ilo != 1) {
        info = -4;
    } else if (n > 0 && (ihi < ilo || ihi > max((INTEGER)1, n))) {
        info = -5;
    } else if (n == 0 && ilo == 1 && ihi != 0) {
        info = -5;
    } else if (m < 0) {
        info = -8;
    } else if (ldv < max((INTEGER)1, n)) {
        info = -10;
    }
    if (info != 0) {
        Mxerbla("Rggbak", -info);
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
        //        Backward transformation on right eigenvectors
        //
        if (rightv) {
            for (i = ilo; i <= ihi; i = i + 1) {
                Rscal(m, rscale[i - 1], &v[(i - 1)], ldv);
            }
        }
        //
        //        Backward transformation on left eigenvectors
        //
        if (leftv) {
            for (i = ilo; i <= ihi; i = i + 1) {
                Rscal(m, lscale[i - 1], &v[(i - 1)], ldv);
            }
        }
    }
//
//     Backward permutation
//
statement_30:
    if (Mlsame(job, "P") || Mlsame(job, "B")) {
        //
        //        Backward permutation on right eigenvectors
        //
        if (rightv) {
            if (ilo == 1) {
                goto statement_50;
            }
            //
            for (i = ilo - 1; i >= 1; i = i - 1) {
                k = INTEGER(rscale[i - 1]);
                if (k == i) {
                    goto statement_40;
                }
                Rswap(m, &v[(i - 1)], ldv, &v[(k - 1)], ldv);
            statement_40:;
            }
        //
        statement_50:
            if (ihi == n) {
                goto statement_70;
            }
            for (i = ihi + 1; i <= n; i = i + 1) {
                k = INTEGER(rscale[i - 1]);
                if (k == i) {
                    goto statement_60;
                }
                Rswap(m, &v[(i - 1)], ldv, &v[(k - 1)], ldv);
            statement_60:;
            }
        }
    //
    //        Backward permutation on left eigenvectors
    //
    statement_70:
        if (leftv) {
            if (ilo == 1) {
                goto statement_90;
            }
            for (i = ilo - 1; i >= 1; i = i - 1) {
                k = INTEGER(lscale[i - 1]);
                if (k == i) {
                    goto statement_80;
                }
                Rswap(m, &v[(i - 1)], ldv, &v[(k - 1)], ldv);
            statement_80:;
            }
        //
        statement_90:
            if (ihi == n) {
                goto statement_110;
            }
            for (i = ihi + 1; i <= n; i = i + 1) {
                k = INTEGER(lscale[i - 1]);
                if (k == i) {
                    goto statement_100;
                }
                Rswap(m, &v[(i - 1)], ldv, &v[(k - 1)], ldv);
            statement_100:;
            }
        }
    }
//
statement_110:;
    //
    //     End of Rggbak
    //
}
