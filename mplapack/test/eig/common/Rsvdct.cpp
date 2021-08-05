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

void Rsvdct(INTEGER const n, REAL *s, REAL *e, REAL const shift, INTEGER &num) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Get machine constants
    //
    REAL unfl = 2 * Rlamch("Safe minimum");
    const REAL one = 1.0;
    REAL ovfl = one / unfl;
    //
    //     Find largest entry
    //
    REAL mx = abs(s[1 - 1]);
    INTEGER i = 0;
    for (i = 1; i <= n - 1; i = i + 1) {
        mx = max({mx, REAL(abs(s[(i + 1) - 1])), REAL(abs(e[i - 1]))});
    }
    //
    const REAL zero = 0.0;
    if (mx == zero) {
        if (shift < zero) {
            num = 0;
        } else {
            num = 2 * n;
        }
        return;
    }
    //
    //     Compute scale factors as in Kahan's report
    //
    REAL sun = sqrt(unfl);
    REAL ssun = sqrt(sun);
    REAL sov = sqrt(ovfl);
    REAL tom = ssun * sov;
    REAL m1 = 0.0;
    REAL m2 = 0.0;
    if (mx <= one) {
        m1 = one / mx;
        m2 = tom;
    } else {
        m1 = one;
        m2 = tom / mx;
    }
    //
    //     Begin counting
    //
    REAL u = one;
    num = 0;
    REAL sshift = (shift * m1) * m2;
    u = -sshift;
    if (u <= sun) {
        if (u <= zero) {
            num++;
            if (u > -sun) {
                u = -sun;
            }
        } else {
            u = sun;
        }
    }
    REAL tmp = (s[1 - 1] * m1) * m2;
    u = -tmp * (tmp / u) - sshift;
    if (u <= sun) {
        if (u <= zero) {
            num++;
            if (u > -sun) {
                u = -sun;
            }
        } else {
            u = sun;
        }
    }
    for (i = 1; i <= n - 1; i = i + 1) {
        tmp = (e[i - 1] * m1) * m2;
        u = -tmp * (tmp / u) - sshift;
        if (u <= sun) {
            if (u <= zero) {
                num++;
                if (u > -sun) {
                    u = -sun;
                }
            } else {
                u = sun;
            }
        }
        tmp = (s[(i + 1) - 1] * m1) * m2;
        u = -tmp * (tmp / u) - sshift;
        if (u <= sun) {
            if (u <= zero) {
                num++;
                if (u > -sun) {
                    u = -sun;
                }
            } else {
                u = sun;
            }
        }
    }
    //
    //     End of Rsvdct
    //
}
