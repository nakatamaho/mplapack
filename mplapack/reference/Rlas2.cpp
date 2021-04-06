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

void Rlas2(REAL const f, REAL const g, REAL const h, REAL &ssmin, REAL &ssmax) {
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  ====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    REAL fa = abs(f);
    REAL ga = abs(g);
    REAL ha = abs(h);
    REAL fhmn = min(fa, ha);
    REAL fhmx = max(fa, ha);
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL as = 0.0;
    REAL at = 0.0;
    REAL au = 0.0;
    const REAL two = 2.0;
    REAL c = 0.0;
    if (fhmn == zero) {
        ssmin = zero;
        if (fhmx == zero) {
            ssmax = ga;
        } else {
            ssmax = max(fhmx, ga) * sqrt(one + pow2((min(fhmx, ga) / max(fhmx, ga))));
        }
    } else {
        if (ga < fhmx) {
            as = one + fhmn / fhmx;
            at = (fhmx - fhmn) / fhmx;
            au = pow2((ga / fhmx));
            c = two / (sqrt(as * as + au) + sqrt(at * at + au));
            ssmin = fhmn * c;
            ssmax = fhmx / c;
        } else {
            au = fhmx / ga;
            if (au == zero) {
                //
                //              Avoid possible harmful underflow if exponent range
                //              asymmetric (true SSMIN may not underflow even if
                //              AU underflows)
                //
                ssmin = (fhmn * fhmx) / ga;
                ssmax = ga;
            } else {
                as = one + fhmn / fhmx;
                at = (fhmx - fhmn) / fhmx;
                c = one / (sqrt(one + pow2((as * au))) + sqrt(one + pow2((at * au))));
                ssmin = (fhmn * c) * au;
                ssmin += ssmin;
                ssmax = ga / (c + c);
            }
        }
    }
    //
    //     End of Rlas2
    //
}
