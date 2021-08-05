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

void Rstech(INTEGER const n, REAL *a, REAL *b, REAL *eig, REAL const tol, REAL *work, INTEGER &info) {
    const REAL zero = 0.0;
    REAL eps = 0.0;
    REAL unflep = 0.0;
    REAL mx = 0.0;
    INTEGER i = 0;
    INTEGER isub = 0;
    REAL emin = 0.0;
    INTEGER j = 0;
    INTEGER tpnt = 0;
    INTEGER bpnt = 0;
    REAL upper = 0.0;
    REAL lower = 0.0;
    REAL tuppr = 0.0;
    INTEGER numl = 0;
    INTEGER numu = 0;
    INTEGER count = 0;
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
    //     Check input parameters
    //
    info = 0;
    if (n == 0) {
        return;
    }
    if (n < 0) {
        info = -1;
        return;
    }
    if (tol < zero) {
        info = -5;
        return;
    }
    //
    //     Get machine constants
    //
    eps = Rlamch("Epsilon") * Rlamch("Base");
    unflep = Rlamch("Safe minimum") / eps;
    eps = tol * eps;
    //
    //     Compute maximum absolute eigenvalue, error tolerance
    //
    mx = abs(eig[1 - 1]);
    for (i = 2; i <= n; i = i + 1) {
        mx = max(mx, REAL(abs(eig[i - 1])));
    }
    eps = max(REAL(eps * mx), unflep);
    //
    //     Sort eigenvalues from EIG into WORK
    //
    for (i = 1; i <= n; i = i + 1) {
        work[i - 1] = eig[i - 1];
    }
    for (i = 1; i <= n - 1; i = i + 1) {
        isub = 1;
        emin = work[1 - 1];
        for (j = 2; j <= n + 1 - i; j = j + 1) {
            if (work[j - 1] < emin) {
                isub = j;
                emin = work[j - 1];
            }
        }
        if (isub != n + 1 - i) {
            work[isub - 1] = work[(n + 1 - i) - 1];
            work[(n + 1 - i) - 1] = emin;
        }
    }
    //
    //     TPNT points to singular value at right endpoint of interval
    //     BPNT points to singular value at left  endpoint of interval
    //
    tpnt = 1;
    bpnt = 1;
//
//     Begin loop over all intervals
//
statement_50:
    upper = work[tpnt - 1] + eps;
    lower = work[bpnt - 1] - eps;
//
//     Begin loop merging overlapping intervals
//
statement_60:
    if (bpnt == n) {
        goto statement_70;
    }
    tuppr = work[(bpnt + 1) - 1] + eps;
    if (tuppr < lower) {
        goto statement_70;
    }
    //
    //     Merge
    //
    bpnt++;
    lower = work[bpnt - 1] - eps;
    goto statement_60;
statement_70:
    //
    //     Count singular values in interval [ LOWER, UPPER ]
    //
    Rstect(n, a, b, lower, numl);
    Rstect(n, a, b, upper, numu);
    count = numu - numl;
    if (count != bpnt - tpnt + 1) {
        //
        //        Wrong number of singular values in interval
        //
        info = tpnt;
        goto statement_80;
    }
    tpnt = bpnt + 1;
    bpnt = tpnt;
    if (tpnt <= n) {
        goto statement_50;
    }
statement_80:;
    //
    //     End of Rstech
    //
}
