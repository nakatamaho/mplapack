/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ZRSCL.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Crscl(INTEGER const n, COMPLEX const a, COMPLEX *x, INTEGER const incx) {
    //
    // Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    // Get machine parameters
    //
    REAL safmin = Rlamch("S");
    const REAL one = 1.0;
    REAL safmax = one / safmin;
    REAL ov = Rlamch("O");
    //
    // Initialize constants related to A.
    //
    REAL ar = a.real();
    REAL ai = a.imag();
    REAL absr = abs(ar);
    REAL absi = abs(ai);
    //
    const REAL zero = 0.0;
    REAL ur = 0.0;
    REAL ui = 0.0;
    if (ai == zero) {
        // If alpha is real, then we can use csrscl
        CRrscl(n, ar, x, incx);
        //
    } else if (ar == zero) {
        // If alpha has a zero real part, then we follow the same rules as if
        // alpha were real.
        if (absi > safmax) {
            CRscal(n, safmin, x, incx);
            Cscal(n, COMPLEX(zero, -safmax / ai), x, incx);
        } else if (absi < safmin) {
            Cscal(n, COMPLEX(zero, -safmin / ai), x, incx);
            CRscal(n, safmax, x, incx);
        } else {
            Cscal(n, COMPLEX(zero, -one / ai), x, incx);
        }
        //
    } else {
        // The following numbers can be computed.
        // They are the inverse of the real and imaginary parts of 1/alpha.
        // Note that a and b are always different from zero.
        // NaNs are only possible if either:
        // 1. alphaR or alphaI is NaN.
        // 2. alphaR and alphaI are both infinite, in which case it makes sense
        // to propagate a NaN.
        ur = ar + ai * (ai / ar);
        ui = ai + ar * (ar / ai);
        //
        if ((abs(ur) < safmin) || (abs(ui) < safmin)) {
            // This means that both alphaR and alphaI are very small.
            Cscal(n, COMPLEX(safmin / ur, -safmin / ui), x, incx);
            CRscal(n, safmax, x, incx);
        } else if ((abs(ur) > safmax) || (abs(ui) > safmax)) {
            if ((absr > ov) || (absi > ov)) {
                // This means that a and b are both Inf. No need for scaling.
                Cscal(n, COMPLEX(one / ur, -one / ui), x, incx);
            } else {
                CRscal(n, safmin, x, incx);
                if ((abs(ur) > ov) || (abs(ui) > ov)) {
                    // Infs were generated. We do proper scaling to avoid them.
                    if (absr >= absi) {
                        // ABS( UR ) <= ABS( UI )
                        ur = (safmin * ar) + safmin * (ai * (ai / ar));
                        ui = (safmin * ai) + ar * ((safmin * ar) / ai);
                    } else {
                        // ABS( UR ) > ABS( UI )
                        ur = (safmin * ar) + ai * ((safmin * ai) / ar);
                        ui = (safmin * ai) + safmin * (ar * (ar / ai));
                    }
                    Cscal(n, COMPLEX(one / ur, -one / ui), x, incx);
                } else {
                    Cscal(n, COMPLEX(safmax / ur, -safmax / ui), x, incx);
                }
            }
        } else {
            Cscal(n, COMPLEX(one / ur, -one / ui), x, incx);
        }
    }
    //
    // End of Crscl
    //
}
