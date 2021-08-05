/*
 * Copyright (c) 2021
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

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_eig.h>

#include <mplapack_debug.h>

void Cunt03(const char *rc, INTEGER const mu, INTEGER const mv, INTEGER const n, INTEGER const k, COMPLEX *u, INTEGER const ldu, COMPLEX *v, INTEGER const ldv, COMPLEX *work, INTEGER const lwork, REAL *rwork, REAL &result, INTEGER &info) {
    //
    //  -- LAPACK test routine --
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Check inputs
    //
    info = 0;
    INTEGER irc = 0;
    if (Mlsame(rc, "R")) {
        irc = 0;
    } else if (Mlsame(rc, "C")) {
        irc = 1;
    } else {
        irc = -1;
    }
    if (irc == -1) {
        info = -1;
    } else if (mu < 0) {
        info = -2;
    } else if (mv < 0) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (k < 0 || k > max(mu, mv)) {
        info = -5;
    } else if ((irc == 0 && ldu < max((INTEGER)1, mu)) || (irc == 1 && ldu < max((INTEGER)1, n))) {
        info = -7;
    } else if ((irc == 0 && ldv < max((INTEGER)1, mv)) || (irc == 1 && ldv < max((INTEGER)1, n))) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Cunt03", -info);
        return;
    }
    //
    //     Initialize result
    //
    const REAL zero = 0.0;
    result = zero;
    if (mu == 0 || mv == 0 || n == 0) {
        return;
    }
    //
    //     Machine constants
    //
    REAL ulp = Rlamch("Precision");
    //
    REAL res1 = 0.0;
    INTEGER i = 0;
    INTEGER lmx = 0;
    const REAL one = 1.0;
    COMPLEX sv = 0.0;
    COMPLEX su = 0.0;
    COMPLEX s = 0.0;
    INTEGER j = 0;
    REAL res2 = 0.0;
    if (irc == 0) {
        //
        //        Compare rows
        //
        res1 = zero;
        for (i = 1; i <= k; i = i + 1) {
            lmx = iCamax(n, &u[(i - 1)], ldu);
            if (v[(i - 1) + (lmx - 1) * ldv] == COMPLEX(zero)) {
                sv = one;
            } else {
                sv = abs(v[(i - 1) + (lmx - 1) * ldv]) / v[(i - 1) + (lmx - 1) * ldv];
            }
            if (u[(i - 1) + (lmx - 1) * ldu] == COMPLEX(zero)) {
                su = one;
            } else {
                su = abs(u[(i - 1) + (lmx - 1) * ldu]) / u[(i - 1) + (lmx - 1) * ldu];
            }
            s = sv / su;
            for (j = 1; j <= n; j = j + 1) {
                res1 = max(res1, abs(u[(i - 1) + (j - 1) * ldu] - s * v[(i - 1) + (j - 1) * ldv]));
            }
        }
        res1 = res1 / (castREAL(n) * ulp);
        //
        //        Compute orthogonality of rows of V.
        //
        Cunt01("Rows", mv, n, v, ldv, work, lwork, rwork, res2);
        //
    } else {
        //
        //        Compare columns
        //
        res1 = zero;
        for (i = 1; i <= k; i = i + 1) {
            lmx = iCamax(n, &u[(i - 1) * ldu], 1);
            if (v[(lmx - 1) + (i - 1) * ldv] == COMPLEX(zero)) {
                sv = one;
            } else {
                sv = abs(v[(lmx - 1) + (i - 1) * ldv]) / v[(lmx - 1) + (i - 1) * ldv];
            }
            if (u[(lmx - 1) + (i - 1) * ldu] == COMPLEX(zero)) {
                su = one;
            } else {
                su = abs(u[(lmx - 1) + (i - 1) * ldu]) / u[(lmx - 1) + (i - 1) * ldu];
            }
            s = sv / su;
            for (j = 1; j <= n; j = j + 1) {
                res1 = max(res1, abs(u[(j - 1) + (i - 1) * ldu] - s * v[(j - 1) + (i - 1) * ldv]));
            }
        }
        res1 = res1 / (castREAL(n) * ulp);
        //
        //        Compute orthogonality of columns of V.
        //
        Cunt01("Columns", n, mv, v, ldv, work, lwork, rwork, res2);
    }
    //
    result = min(REAL(max(res1, res2)), REAL(one / ulp));
    //
    //     End of Cunt03
    //
}
