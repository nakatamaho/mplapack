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

void Cla_lin_berr(INTEGER const &n, INTEGER const &nz, INTEGER const &nrhs, COMPLEX *res, REAL *ayb, REAL *berr) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function Definitions ..
    COMPLEX cdum = 0.0;
    abs1[cdum - 1] = abs(cdum.real()) + abs(cdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     Adding SAFE1 to the numerator guards against spuriously zero
    //     residuals.  A similar safeguard is in the CLA_yyAMV routine used
    //     to compute AYB.
    //
    REAL safe1 = dlamch("Safe minimum");
    safe1 = (nz + 1) * safe1;
    //
    INTEGER j = 0;
    INTEGER i = 0;
    REAL tmp = 0.0;
    for (j = 1; j <= nrhs; j = j + 1) {
        berr[j - 1] = 0.0;
        for (i = 1; i <= n; i = i + 1) {
            if (ayb[(i - 1) + (j - 1) * ldayb] != 0.0) {
                tmp = (safe1 + abs1[res[(i - 1) + (j - 1) * ldres] - 1]) / ayb[(i - 1) + (j - 1) * ldayb];
                berr[j - 1] = max(berr[j - 1], tmp);
            }
            //
            //     If AYB is exactly 0.0 (and if computed by CLA_yyAMV), then we know
            //     the true residual also must be exactly 0.0.
            //
        }
    }
}
