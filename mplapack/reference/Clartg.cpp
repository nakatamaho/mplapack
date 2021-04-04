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

void Clartg(COMPLEX const &f, COMPLEX const &g, REAL &cs, COMPLEX &sn, COMPLEX &r) {
    COMPLEX ff = 0.0;
    REAL safmin = 0.0;
    REAL eps = 0.0;
    const REAL two = 2.0e+0;
    REAL safmn2 = 0.0;
    const REAL one = 1.0;
    REAL safmx2 = 0.0;
    REAL scale = 0.0;
    COMPLEX fs = 0.0;
    COMPLEX gs = 0.0;
    INTEGER count = 0;
    const COMPLEX czero = (0.0, 0.0);
    REAL f2 = 0.0;
    REAL g2 = 0.0;
    const REAL zero = 0.0;
    REAL d = 0.0;
    REAL f2s = 0.0;
    REAL g2s = 0.0;
    REAL dr = 0.0;
    REAL di = 0.0;
    INTEGER i = 0;
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     LOGICAL            FIRST
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    abs1[ff - 1] = max(abs(ff.real()), abs(ff.imag()));
    abssq[ff - 1] = pow2(ff.real()) + pow2(ff.imag());
    //     ..
    //     .. Executable Statements ..
    //
    safmin = dlamch("S");
    eps = dlamch("E");
    safmn2 = pow(dlamch("B"), INTEGER(log[(safmin / eps) - 1] / log[dlamch("B") - 1] / two));
    safmx2 = one / safmn2;
    scale = max(abs1[f - 1], abs1[g - 1]);
    fs = f;
    gs = g;
    count = 0;
    if (scale >= safmx2) {
    statement_10:
        count++;
        fs = fs * safmn2;
        gs = gs * safmn2;
        scale = scale * safmn2;
        if (scale >= safmx2 && count < 20) {
            goto statement_10;
        }
    } else if (scale <= safmn2) {
        if (g == czero || Risnan(abs(g))) {
            cs = one;
            sn = czero;
            r = f;
            return;
        }
    statement_20:
        count = count - 1;
        fs = fs * safmx2;
        gs = gs * safmx2;
        scale = scale * safmx2;
        if (scale <= safmn2) {
            goto statement_20;
        }
    }
    f2 = abssq[fs - 1];
    g2 = abssq[gs - 1];
    if (f2 <= max(g2, one) * safmin) {
        //
        //        This is a rare case: F is very small.
        //
        if (f == czero) {
            cs = zero;
            r = Rlapy2[(g.real() - 1) + (g.imag() - 1) * ldRlapy2];
            //           Do complex/real division explicitly with two real divisions
            d = Rlapy2[(gs.real() - 1) + (gs.imag() - 1) * ldRlapy2];
            sn = COMPLEX(gs.real() / d, -gs.imag() / d);
            return;
        }
        f2s = Rlapy2[(fs.real() - 1) + (fs.imag() - 1) * ldRlapy2];
        //        G2 and G2S are accurate
        //        G2 is at least SAFMIN, and G2S is at least SAFMN2
        g2s = sqrt(g2);
        //        Error in CS from underflow in F2S is at most
        //        UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS
        //        If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN,
        //        and so CS .lt. sqrt(SAFMIN)
        //        If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN
        //        and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
        //        Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
        cs = f2s / g2s;
        //        Make sure abs(FF) = 1
        //        Do complex/real division explicitly with 2 real divisions
        if (abs1[f - 1] > one) {
            d = Rlapy2[(f.real() - 1) + (f.imag() - 1) * ldRlapy2];
            ff = COMPLEX(f.real() / d, f.imag() / d);
        } else {
            dr = safmx2 * f.real();
            di = safmx2 * f.imag();
            d = Rlapy2[(dr - 1) + (di - 1) * ldRlapy2];
            ff = COMPLEX(dr / d, di / d);
        }
        sn = ff * COMPLEX(gs.real() / g2s, -gs.imag() / g2s);
        r = cs * f + sn * g;
    } else {
        //
        //        This is the most common case.
        //        Neither F2 nor F2/G2 are less than SAFMIN
        //        F2S cannot overflow, and it is accurate
        //
        f2s = sqrt(one + g2 / f2);
        //        Do the F2S(real)*FS(complex) multiply with two real multiplies
        r = COMPLEX(f2s * fs.real(), f2s * fs.imag());
        cs = one / f2s;
        d = f2 + g2;
        //        Do complex/real division explicitly with two real divisions
        sn = COMPLEX(r.real() / d, r.imag() / d);
        sn = sn * conj(gs);
        if (count != 0) {
            if (count > 0) {
                for (i = 1; i <= count; i = i + 1) {
                    r = r * safmx2;
                }
            } else {
                for (i = 1; i <= -count; i = i + 1) {
                    r = r * safmn2;
                }
            }
        }
    }
    //
    //     End of Clartg
    //
}
