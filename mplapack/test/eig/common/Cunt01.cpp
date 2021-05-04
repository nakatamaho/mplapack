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

void Cunt01(const char *rowcol, INTEGER const m, INTEGER const n, COMPLEX *u, INTEGER const ldu, COMPLEX *work, INTEGER const lwork, REAL *rwork, REAL &resid) {
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
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    COMPLEX zdum = 0.0;
    abs1(zdum) = abs(zdum.real()) + abs(zdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    const REAL zero = 0.0;
    resid = zero;
    //
    //     Quick return if possible
    //
    if (m <= 0 || n <= 0) {
        return;
    }
    //
    REAL eps = Rlamch("Precision");
    char transu[1];
    INTEGER k = 0;
    if (m < n || (m == n && Mlsame(rowcol, "R"))) {
        transu = "N";
        k = n;
    } else {
        transu = "C";
        k = m;
    }
    INTEGER mnmin = min(m, n);
    //
    INTEGER ldwork = 0;
    if ((mnmin + 1) * mnmin <= lwork) {
        ldwork = mnmin;
    } else {
        ldwork = 0;
    }
    const REAL one = 1.0;
    INTEGER j = 0;
    INTEGER i = 0;
    COMPLEX tmp = 0.0;
    if (ldwork > 0) {
        //
        //        Compute I - U*U' or I - U'*U.
        //
        Claset("Upper", mnmin, mnmin, COMPLEX(zero), COMPLEX(one), work, ldwork);
        Cherk("Upper", transu, mnmin, k, -one, u, ldu, one, work, ldwork);
        //
        //        Compute norm( I - U*U' ) / ( K * EPS ) .
        //
        resid = Clansy("1", "Upper", mnmin, work, ldwork, rwork);
        resid = (resid / k.real()) / eps;
    } else if (transu == "C") {
        //
        //        Find the maximum element in abs( I - U'*U ) / ( m * EPS )
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= j; i = i + 1) {
                if (i != j) {
                    tmp = zero;
                } else {
                    tmp = one;
                }
                tmp = tmp - Cdotc(m, &u[(i - 1) * ldu], 1, &u[(j - 1) * ldu], 1);
                resid = max(resid, abs1(tmp));
            }
        }
        resid = (resid / m.real()) / eps;
    } else {
        //
        //        Find the maximum element in abs( I - U*U' ) / ( n * EPS )
        //
        for (j = 1; j <= m; j = j + 1) {
            for (i = 1; i <= j; i = i + 1) {
                if (i != j) {
                    tmp = zero;
                } else {
                    tmp = one;
                }
                tmp = tmp - Cdotc(n, &u[(j - 1)], ldu, &u[(i - 1)], ldu);
                resid = max(resid, abs1(tmp));
            }
        }
        resid = (resid / n.real()) / eps;
    }
    //
    //     End of Cunt01
    //
}
