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

void Rlasv2(REAL const &f, REAL const &g, REAL const &h, REAL &ssmin, REAL &ssmax, REAL &snr, REAL &csr, REAL &snl, REAL &csl) {
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    // =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    REAL ft = f;
    REAL fa = abs(ft);
    REAL ht = h;
    REAL ha = abs(h);
    //
    //     PMAX poINTEGERs to the maximum absolute element of matrix
    //       PMAX = 1 if F largest in absolute values
    //       PMAX = 2 if G largest in absolute values
    //       PMAX = 3 if H largest in absolute values
    //
    INTEGER pmax = 1;
    bool swap = (ha > fa);
    REAL temp = 0.0;
    if (swap) {
        pmax = 3;
        temp = ft;
        ft = ht;
        ht = temp;
        temp = fa;
        fa = ha;
        ha = temp;
        //
        //        Now FA .ge. HA
        //
    }
    REAL gt = g;
    REAL ga = abs(gt);
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL clt = 0.0;
    REAL crt = 0.0;
    REAL slt = 0.0;
    REAL srt = 0.0;
    bool gasmal = false;
    REAL d = 0.0;
    REAL l = 0.0;
    REAL m = 0.0;
    const REAL two = 2.0;
    REAL t = 0.0;
    REAL mm = 0.0;
    REAL tt = 0.0;
    REAL s = 0.0;
    REAL r = 0.0;
    const REAL half = 0.5e0;
    REAL a = 0.0;
    const REAL four = 4.0;
    if (ga == zero) {
        //
        //        Diagonal matrix
        //
        ssmin = ha;
        ssmax = fa;
        clt = one;
        crt = one;
        slt = zero;
        srt = zero;
    } else {
        gasmal = true;
        if (ga > fa) {
            pmax = 2;
            if ((fa / ga) < dlamch("EPS")) {
                //
                //              Case of very large GA
                //
                gasmal = false;
                ssmax = ga;
                if (ha > one) {
                    ssmin = fa / (ga / ha);
                } else {
                    ssmin = (fa / ga) * ha;
                }
                clt = one;
                slt = ht / gt;
                srt = one;
                crt = ft / gt;
            }
        }
        if (gasmal) {
            //
            //           Normal case
            //
            d = fa - ha;
            if (d == fa) {
                //
                //              Copes with infinite F or H
                //
                l = one;
            } else {
                l = d / fa;
            }
            //
            //           Note that 0 .le. L .le. 1
            //
            m = gt / ft;
            //
            //           Note that abs(M) .le. 1/macheps
            //
            t = two - l;
            //
            //           Note that T .ge. 1
            //
            mm = m * m;
            tt = t * t;
            s = sqrt(tt + mm);
            //
            //           Note that 1 .le. S .le. 1 + 1/macheps
            //
            if (l == zero) {
                r = abs(m);
            } else {
                r = sqrt(l * l + mm);
            }
            //
            //           Note that 0 .le. R .le. 1 + 1/macheps
            //
            a = half * (s + r);
            //
            //           Note that 1 .le. A .le. 1 + abs(M)
            //
            ssmin = ha / a;
            ssmax = fa * a;
            if (mm == zero) {
                //
                //              Note that M is very tiny
                //
                if (l == zero) {
                    t = sign[(two - 1) + (ft - 1) * ldsign] * sign[(one - 1) + (gt - 1) * ldsign];
                } else {
                    t = gt / sign[(d - 1) + (ft - 1) * ldsign] + m / t;
                }
            } else {
                t = (m / (s + t) + m / (r + l)) * (one + a);
            }
            l = sqrt(t * t + four);
            crt = two / l;
            srt = t / l;
            clt = (crt + srt * m) / a;
            slt = (ht / ft) * srt / a;
        }
    }
    if (swap) {
        csl = srt;
        snl = crt;
        csr = slt;
        snr = clt;
    } else {
        csl = clt;
        snl = slt;
        csr = crt;
        snr = srt;
    }
    //
    //     Correct signs of SSMAX and SSMIN
    //
    REAL tsign = 0.0;
    if (pmax == 1) {
        tsign = sign[(one - 1) + (csr - 1) * ldsign] * sign[(one - 1) + (csl - 1) * ldsign] * sign[(one - 1) + (f - 1) * ldsign];
    }
    if (pmax == 2) {
        tsign = sign[(one - 1) + (snr - 1) * ldsign] * sign[(one - 1) + (csl - 1) * ldsign] * sign[(one - 1) + (g - 1) * ldsign];
    }
    if (pmax == 3) {
        tsign = sign[(one - 1) + (snr - 1) * ldsign] * sign[(one - 1) + (snl - 1) * ldsign] * sign[(one - 1) + (h - 1) * ldsign];
    }
    ssmax = sign[(ssmax - 1) + (tsign - 1) * ldsign];
    ssmin = sign[(ssmin - 1) + ((tsign * sign[(one - 1) + (f - 1) * ldsign] * sign[(one - 1) + (h - 1) * ldsign]) - 1) * ldsign];
    //
    //     End of Rlasv2
    //
}
