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

void Clarfg(INTEGER const n, COMPLEX &alpha, COMPLEX *x, INTEGER const incx, COMPLEX &tau) {
    const REAL zero = 0.0;
    REAL xnorm = 0.0;
    REAL alphr = 0.0;
    REAL alphi = 0.0;
    REAL beta = 0.0;
    REAL safmin = 0.0;
    const REAL one = 1.0;
    REAL rsafmn = 0.0;
    INTEGER knt = 0;
    INTEGER j = 0;
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
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    if (n <= 0) {
        tau = zero;
        return;
    }
    //
    xnorm = RCnrm2(n - 1, x, incx);
    alphr = alpha.real();
    alphi = alpha.imag();
    //
    if (xnorm == zero && alphi == zero) {
        //
        //        H  =  I
        //
        tau = zero;
    } else {
        //
        //        general case
        //
        beta = -sign(Rlapy3(alphr, alphi, xnorm), alphr);
        safmin = Rlamch("S") / Rlamch("E");
        rsafmn = one / safmin;
        //
        knt = 0;
        if (abs(beta) < safmin) {
        //
        //           XNORM, BETA may be inaccurate; scale X and recompute them
        //
        statement_10:
            knt++;
            CRscal(n - 1, rsafmn, x, incx);
            beta = beta * rsafmn;
            alphi = alphi * rsafmn;
            alphr = alphr * rsafmn;
            if ((abs(beta) < safmin) && (knt < 20)) {
                goto statement_10;
            }
            //
            //           New BETA is at most 1, at least SAFMIN
            //
            xnorm = RCnrm2(n - 1, x, incx);
            alpha = COMPLEX(alphr, alphi);
            beta = -sign(Rlapy3(alphr, alphi, xnorm), alphr);
        }
        tau = COMPLEX((beta - alphr) / beta, -alphi / beta);
        alpha = Cladiv(COMPLEX(one), alpha - beta);
        Cscal(n - 1, alpha, x, incx);
        //
        //        If ALPHA is subnormal, it may lose relative accuracy
        //
        for (j = 1; j <= knt; j = j + 1) {
            beta = beta * safmin;
        }
        alpha = beta;
    }
    //
    //     End of Clarfg
    //
}
