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

void Rsvdch(INTEGER const n, REAL *s, REAL *e, REAL *svd, REAL const tol, INTEGER &info) {
    REAL unfl = 0.0;
    REAL ovfl = 0.0;
    REAL eps = 0.0;
    REAL unflep = 0.0;
    INTEGER tpnt = 0;
    INTEGER bpnt = 0;
    const REAL one = 1.0;
    REAL upper = 0.0;
    REAL lower = 0.0;
    REAL tuppr = 0.0;
    INTEGER numl = 0;
    INTEGER numu = 0;
    INTEGER count = 0;
    const REAL zero = 0.0;
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Get machine constants
    //
    info = 0;
    if (n <= 0) {
        return;
    }
    unfl = Rlamch("Safe minimum");
    ovfl = Rlamch("Overflow");
    eps = Rlamch("Epsilon") * Rlamch("Base");
    //
    //     UNFLEP is chosen so that when an eigenvalue is multiplied by the
    //     scale factor sqrt(OVFL)*sqrt(sqrt(UNFL))/MX in Rsvdct, it exceeds
    //     sqrt(UNFL), which is the lower limit for Rsvdct.
    //
    unflep = (sqrt(sqrt(unfl)) / sqrt(ovfl)) * svd[1 - 1] + unfl / eps;
    //
    //     The value of EPS works best when TOL .GE. 10.
    //
    eps = tol * max(REAL(castREAL(n) / 10.0), REAL(1.0)) * eps;
    //
    //     TPNT points to singular value at right endpoint of interval
    //     BPNT points to singular value at left  endpoint of interval
    //
    tpnt = 1;
    bpnt = 1;
//
//     Begin loop over all intervals
//
statement_10:
    upper = (one + eps) * svd[tpnt - 1] + unflep;
    lower = (one - eps) * svd[bpnt - 1] - unflep;
    if (lower <= unflep) {
        lower = -upper;
    }
//
//     Begin loop merging overlapping intervals
//
statement_20:
    if (bpnt == n) {
        goto statement_30;
    }
    tuppr = (one + eps) * svd[(bpnt + 1) - 1] + unflep;
    if (tuppr < lower) {
        goto statement_30;
    }
    //
    //     Merge
    //
    bpnt++;
    lower = (one - eps) * svd[bpnt - 1] - unflep;
    if (lower <= unflep) {
        lower = -upper;
    }
    goto statement_20;
statement_30:
    //
    //     Count singular values in interval [ LOWER, UPPER ]
    //
    Rsvdct(n, s, e, lower, numl);
    Rsvdct(n, s, e, upper, numu);
    count = numu - numl;
    if (lower < zero) {
        count = count / 2;
    }
    if (count != bpnt - tpnt + 1) {
        //
        //        Wrong number of singular values in interval
        //
        info = tpnt;
        goto statement_40;
    }
    tpnt = bpnt + 1;
    bpnt = tpnt;
    if (tpnt <= n) {
        goto statement_10;
    }
statement_40:;
    //
    //     End of Rsvdch
    //
}
