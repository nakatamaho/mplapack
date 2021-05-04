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

void Rcsdts(INTEGER const m, INTEGER const p, INTEGER const q, REAL *x, REAL *xf, INTEGER const ldx, REAL *u1, INTEGER const ldu1, REAL *u2, INTEGER const ldu2, REAL *v1t, INTEGER const ldv1t, REAL *v2t, INTEGER const ldv2t, REAL *theta, INTEGER *iwork, REAL *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    REAL ulp = Rlamch("Precision");
    const REAL realone = 1.0;
    REAL ulpinv = realone / ulp;
    //
    //     The first half of the routine checks the 2-by-2 CSD
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    Rlaset("Full", m, m, zero, one, work, ldx);
    Rsyrk("Upper", "Conjugate transpose", m, m, -one, x, ldx, one, work, ldx);
    REAL eps2 = 0.0;
    if (m > 0) {
        eps2 = max({ulp, Rlange("1", m, m, work, ldx, rwork) / m.real()});
    } else {
        eps2 = ulp;
    }
    INTEGER r = min({p, m - p, q, m - q});
    //
    //     Copy the matrix X to the array XF.
    //
    Rlacpy("Full", m, m, x, ldx, xf, ldx);
    //
    //     Compute the CSD
    //
    INTEGER info = 0;
    Rorcsd("Y", "Y", "Y", "Y", "N", "D", m, p, q, xf[(1 - 1)], ldx, xf[((q + 1) - 1) * ldxf], ldx, xf[((p + 1) - 1)], ldx, xf[((p + 1) - 1) + ((q + 1) - 1) * ldxf], ldx, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, iwork, info);
    //
    //     Compute XF := diag(U1,U2)'*X*diag(V1,V2) - [D11 D12; D21 D22]
    //
    Rlacpy("Full", m, m, x, ldx, xf, ldx);
    //
    Rgemm("No transpose", "Conjugate transpose", p, q, q, one, xf, ldx, v1t, ldv1t, zero, work, ldx);
    //
    Rgemm("Conjugate transpose", "No transpose", p, q, p, one, u1, ldu1, work, ldx, zero, xf, ldx);
    //
    INTEGER i = 0;
    for (i = 1; i <= min(p, q) - r; i = i + 1) {
        xf[(i - 1) + (i - 1) * ldxf] = xf[(i - 1) + (i - 1) * ldxf] - one;
    }
    for (i = 1; i <= r; i = i + 1) {
        xf(min(p, q) - r + i, min(p, q) - r + i) = xf(min(p, q) - r + i, min(p, q) - r + i) - cos(theta[i - 1]);
    }
    //
    Rgemm("No transpose", "Conjugate transpose", p, m - q, m - q, one, xf[((q + 1) - 1) * ldxf], ldx, v2t, ldv2t, zero, work, ldx);
    //
    Rgemm("Conjugate transpose", "No transpose", p, m - q, p, one, u1, ldu1, work, ldx, zero, xf[((q + 1) - 1) * ldxf], ldx);
    //
    for (i = 1; i <= min(p, m - q) - r; i = i + 1) {
        xf[((p - i + 1) - 1) + ((m - i + 1) - 1) * ldxf] += one;
    }
    for (i = 1; i <= r; i = i + 1) {
        xf(p - (min(p, m - q) - r) + 1 - i, m - (min(p, m - q) - r) + 1 - i) += sin(theta[(r - i + 1) - 1]);
    }
    //
    Rgemm("No transpose", "Conjugate transpose", m - p, q, q, one, xf[((p + 1) - 1)], ldx, v1t, ldv1t, zero, work, ldx);
    //
    Rgemm("Conjugate transpose", "No transpose", m - p, q, m - p, one, u2, ldu2, work, ldx, zero, xf[((p + 1) - 1)], ldx);
    //
    for (i = 1; i <= min(m - p, q) - r; i = i + 1) {
        xf[((m - i + 1) - 1) + ((q - i + 1) - 1) * ldxf] = xf[((m - i + 1) - 1) + ((q - i + 1) - 1) * ldxf] - one;
    }
    for (i = 1; i <= r; i = i + 1) {
        xf(m - (min(m - p, q) - r) + 1 - i, q - (min(m - p, q) - r) + 1 - i) = xf(m - (min(m - p, q) - r) + 1 - i, q - (min(m - p, q) - r) + 1 - i) - sin(theta[(r - i + 1) - 1]);
    }
    //
    Rgemm("No transpose", "Conjugate transpose", m - p, m - q, m - q, one, xf[((p + 1) - 1) + ((q + 1) - 1) * ldxf], ldx, v2t, ldv2t, zero, work, ldx);
    //
    Rgemm("Conjugate transpose", "No transpose", m - p, m - q, m - p, one, u2, ldu2, work, ldx, zero, xf[((p + 1) - 1) + ((q + 1) - 1) * ldxf], ldx);
    //
    for (i = 1; i <= min(m - p, m - q) - r; i = i + 1) {
        xf[((p + i) - 1) + ((q + i) - 1) * ldxf] = xf[((p + i) - 1) + ((q + i) - 1) * ldxf] - one;
    }
    for (i = 1; i <= r; i = i + 1) {
        xf(p + (min(m - p, m - q) - r) + i, q + (min(m - p, m - q) - r) + i) = xf(p + (min(m - p, m - q) - r) + i, q + (min(m - p, m - q) - r) + i) - cos(theta[i - 1]);
    }
    //
    //     Compute norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 ) .
    //
    REAL resid = Rlange("1", p, q, xf, ldx, rwork);
    result[1 - 1] = (resid / (max({(INTEGER)1, p, q})).real()) / eps2;
    //
    //     Compute norm( U1'*X12*V2 - D12 ) / ( MAX(1,P,M-Q)*EPS2 ) .
    //
    resid = Rlange("1", p, m - q, xf[((q + 1) - 1) * ldxf], ldx, rwork);
    result[2 - 1] = (resid / (max({(INTEGER)1, p, m - q})).real()) / eps2;
    //
    //     Compute norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 ) .
    //
    resid = Rlange("1", m - p, q, xf[((p + 1) - 1)], ldx, rwork);
    result[3 - 1] = (resid / (max({(INTEGER)1, m - p, q})).real()) / eps2;
    //
    //     Compute norm( U2'*X22*V2 - D22 ) / ( MAX(1,M-P,M-Q)*EPS2 ) .
    //
    resid = Rlange("1", m - p, m - q, xf[((p + 1) - 1) + ((q + 1) - 1) * ldxf], ldx, rwork);
    result[4 - 1] = (resid / (max({(INTEGER)1, m - p, m - q})).real()) / eps2;
    //
    //     Compute I - U1'*U1
    //
    Rlaset("Full", p, p, zero, one, work, ldu1);
    Rsyrk("Upper", "Conjugate transpose", p, p, -one, u1, ldu1, one, work, ldu1);
    //
    //     Compute norm( I - U'*U ) / ( MAX(1,P) * ULP ) .
    //
    resid = Rlansy("1", "Upper", p, work, ldu1, rwork);
    result[5 - 1] = (resid / (max((INTEGER)1, p)).real()) / ulp;
    //
    //     Compute I - U2'*U2
    //
    Rlaset("Full", m - p, m - p, zero, one, work, ldu2);
    Rsyrk("Upper", "Conjugate transpose", m - p, m - p, -one, u2, ldu2, one, work, ldu2);
    //
    //     Compute norm( I - U2'*U2 ) / ( MAX(1,M-P) * ULP ) .
    //
    resid = Rlansy("1", "Upper", m - p, work, ldu2, rwork);
    result[6 - 1] = (resid / (max((INTEGER)1, m - p)).real()) / ulp;
    //
    //     Compute I - V1T*V1T'
    //
    Rlaset("Full", q, q, zero, one, work, ldv1t);
    Rsyrk("Upper", "No transpose", q, q, -one, v1t, ldv1t, one, work, ldv1t);
    //
    //     Compute norm( I - V1T*V1T' ) / ( MAX(1,Q) * ULP ) .
    //
    resid = Rlansy("1", "Upper", q, work, ldv1t, rwork);
    result[7 - 1] = (resid / (max((INTEGER)1, q)).real()) / ulp;
    //
    //     Compute I - V2T*V2T'
    //
    Rlaset("Full", m - q, m - q, zero, one, work, ldv2t);
    Rsyrk("Upper", "No transpose", m - q, m - q, -one, v2t, ldv2t, one, work, ldv2t);
    //
    //     Compute norm( I - V2T*V2T' ) / ( MAX(1,M-Q) * ULP ) .
    //
    resid = Rlansy("1", "Upper", m - q, work, ldv2t, rwork);
    result[8 - 1] = (resid / (max((INTEGER)1, m - q)).real()) / ulp;
    //
    //     Check sorting
    //
    const REAL realzero = 0.0;
    result[9 - 1] = realzero;
    const REAL piover2 = 1.57079632679489661923132169163975144210e0;
    for (i = 1; i <= r; i = i + 1) {
        if (theta[i - 1] < realzero || theta[i - 1] > piover2) {
            result[9 - 1] = ulpinv;
        }
        if (i > 1) {
            if (theta[i - 1] < theta[(i - 1) - 1]) {
                result[9 - 1] = ulpinv;
            }
        }
    }
    //
    //     The second half of the routine checks the 2-by-1 CSD
    //
    Rlaset("Full", q, q, zero, one, work, ldx);
    Rsyrk("Upper", "Conjugate transpose", q, m, -one, x, ldx, one, work, ldx);
    if (m > 0) {
        eps2 = max({ulp, Rlange("1", q, q, work, ldx, rwork) / m.real()});
    } else {
        eps2 = ulp;
    }
    r = min({p, m - p, q, m - q});
    //
    //     Copy the matrix [ X11; X21 ] to the array XF.
    //
    Rlacpy("Full", m, q, x, ldx, xf, ldx);
    //
    //     Compute the CSD
    //
    Rorcsd2by1("Y", "Y", "Y", m, p, q, xf[(1 - 1)], ldx, xf[((p + 1) - 1)], ldx, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork, iwork, info);
    //
    //     Compute [X11;X21] := diag(U1,U2)'*[X11;X21]*V1 - [D11;D21]
    //
    Rgemm("No transpose", "Conjugate transpose", p, q, q, one, x, ldx, v1t, ldv1t, zero, work, ldx);
    //
    Rgemm("Conjugate transpose", "No transpose", p, q, p, one, u1, ldu1, work, ldx, zero, x, ldx);
    //
    for (i = 1; i <= min(p, q) - r; i = i + 1) {
        x[(i - 1) + (i - 1) * ldx] = x[(i - 1) + (i - 1) * ldx] - one;
    }
    for (i = 1; i <= r; i = i + 1) {
        x(min(p, q) - r + i, min(p, q) - r + i) = x(min(p, q) - r + i, min(p, q) - r + i) - cos(theta[i - 1]);
    }
    //
    Rgemm("No transpose", "Conjugate transpose", m - p, q, q, one, &x[((p + 1) - 1)], ldx, v1t, ldv1t, zero, work, ldx);
    //
    Rgemm("Conjugate transpose", "No transpose", m - p, q, m - p, one, u2, ldu2, work, ldx, zero, &x[((p + 1) - 1)], ldx);
    //
    for (i = 1; i <= min(m - p, q) - r; i = i + 1) {
        x[((m - i + 1) - 1) + ((q - i + 1) - 1) * ldx] = x[((m - i + 1) - 1) + ((q - i + 1) - 1) * ldx] - one;
    }
    for (i = 1; i <= r; i = i + 1) {
        x(m - (min(m - p, q) - r) + 1 - i, q - (min(m - p, q) - r) + 1 - i) = x(m - (min(m - p, q) - r) + 1 - i, q - (min(m - p, q) - r) + 1 - i) - sin(theta[(r - i + 1) - 1]);
    }
    //
    //     Compute norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 ) .
    //
    resid = Rlange("1", p, q, x, ldx, rwork);
    result[10 - 1] = (resid / (max({(INTEGER)1, p, q})).real()) / eps2;
    //
    //     Compute norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 ) .
    //
    resid = Rlange("1", m - p, q, &x[((p + 1) - 1)], ldx, rwork);
    result[11 - 1] = (resid / (max({(INTEGER)1, m - p, q})).real()) / eps2;
    //
    //     Compute I - U1'*U1
    //
    Rlaset("Full", p, p, zero, one, work, ldu1);
    Rsyrk("Upper", "Conjugate transpose", p, p, -one, u1, ldu1, one, work, ldu1);
    //
    //     Compute norm( I - U1'*U1 ) / ( MAX(1,P) * ULP ) .
    //
    resid = Rlansy("1", "Upper", p, work, ldu1, rwork);
    result[12 - 1] = (resid / (max((INTEGER)1, p)).real()) / ulp;
    //
    //     Compute I - U2'*U2
    //
    Rlaset("Full", m - p, m - p, zero, one, work, ldu2);
    Rsyrk("Upper", "Conjugate transpose", m - p, m - p, -one, u2, ldu2, one, work, ldu2);
    //
    //     Compute norm( I - U2'*U2 ) / ( MAX(1,M-P) * ULP ) .
    //
    resid = Rlansy("1", "Upper", m - p, work, ldu2, rwork);
    result[13 - 1] = (resid / (max((INTEGER)1, m - p)).real()) / ulp;
    //
    //     Compute I - V1T*V1T'
    //
    Rlaset("Full", q, q, zero, one, work, ldv1t);
    Rsyrk("Upper", "No transpose", q, q, -one, v1t, ldv1t, one, work, ldv1t);
    //
    //     Compute norm( I - V1T*V1T' ) / ( MAX(1,Q) * ULP ) .
    //
    resid = Rlansy("1", "Upper", q, work, ldv1t, rwork);
    result[14 - 1] = (resid / (max((INTEGER)1, q)).real()) / ulp;
    //
    //     Check sorting
    //
    result[15 - 1] = realzero;
    for (i = 1; i <= r; i = i + 1) {
        if (theta[i - 1] < realzero || theta[i - 1] > piover2) {
            result[15 - 1] = ulpinv;
        }
        if (i > 1) {
            if (theta[i - 1] < theta[(i - 1) - 1]) {
                result[15 - 1] = ulpinv;
            }
        }
    }
    //
    //     End of Rcsdts
    //
}
