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

void Rget53(REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const scale, REAL const wr, REAL const wi, REAL &result, INTEGER &info) {
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
    //     Initialize
    //
    info = 0;
    const REAL zero = 0.0;
    result = zero;
    REAL scales = scale;
    REAL wrs = wr;
    REAL wis = wi;
    //
    //     Machine constants and norms
    //
    REAL safmin = Rlamch("Safe minimum");
    REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
    REAL absw = abs(wrs) + abs(wis);
    REAL anorm = max({REAL(abs(a[(1 - 1) + (1 - 1) * lda]) + abs(a[(2 - 1)])), REAL(abs(a[(1 - 1) + (2 - 1) * lda]) + abs(a[(2 - 1) + (2 - 1) * lda])), safmin});
    REAL bnorm = max({REAL(abs(b[(1 - 1) + (1 - 1) * ldb])), REAL(abs(b[(1 - 1) + (2 - 1) * ldb]) + abs(b[(2 - 1) + (2 - 1) * ldb])), safmin});
    //
    //     Check for possible overflow.
    //
    REAL temp = (safmin * bnorm) * absw + (safmin * anorm) * scales;
    const REAL one = 1.0;
    if (temp >= one) {
        //
        //        Scale down to avoid overflow
        //
        info = 1;
        temp = one / temp;
        scales = scales * temp;
        wrs = wrs * temp;
        wis = wis * temp;
        absw = abs(wrs) + abs(wis);
    }
    REAL s1 = max(REAL(ulp * max(scales * anorm, absw * bnorm)), REAL(safmin * max(scales, absw)));
    //
    //     Check for W and SCALE essentially zero.
    //
    if (s1 < safmin) {
        info = 2;
        if (scales < safmin && absw < safmin) {
            info = 3;
            result = one / ulp;
            return;
        }
        //
        //        Scale up to avoid underflow
        //
        temp = one / max(REAL(scales * anorm + absw * bnorm), safmin);
        scales = scales * temp;
        wrs = wrs * temp;
        wis = wis * temp;
        absw = abs(wrs) + abs(wis);
        s1 = max(REAL(ulp * max(scales * anorm, absw * bnorm)), REAL(safmin * max(scales, absw)));
        if (s1 < safmin) {
            info = 3;
            result = one / ulp;
            return;
        }
    }
    //
    //     Compute C = s A - w B
    //
    REAL cr11 = scales * a[(1 - 1) + (1 - 1) * lda] - wrs * b[(1 - 1) + (1 - 1) * ldb];
    REAL ci11 = -wis * b[(1 - 1) + (1 - 1) * ldb];
    REAL cr21 = scales * a[(2 - 1)];
    REAL cr12 = scales * a[(1 - 1) + (2 - 1) * lda] - wrs * b[(1 - 1) + (2 - 1) * ldb];
    REAL ci12 = -wis * b[(1 - 1) + (2 - 1) * ldb];
    REAL cr22 = scales * a[(2 - 1) + (2 - 1) * lda] - wrs * b[(2 - 1) + (2 - 1) * ldb];
    REAL ci22 = -wis * b[(2 - 1) + (2 - 1) * ldb];
    //
    //     Compute the smallest singular value of s A - w B:
    //
    //                 |det( s A - w B )|
    //     sigma_min = ------------------
    //                 norm( s A - w B )
    //
    REAL cnorm = max({REAL(abs(cr11) + abs(ci11) + abs(cr21)), REAL(abs(cr12) + abs(ci12) + abs(cr22) + abs(ci22)), safmin});
    REAL cscale = one / sqrt(cnorm);
    REAL detr = (cscale * cr11) * (cscale * cr22) - (cscale * ci11) * (cscale * ci22) - (cscale * cr12) * (cscale * cr21);
    REAL deti = (cscale * cr11) * (cscale * ci22) + (cscale * ci11) * (cscale * cr22) - (cscale * ci12) * (cscale * cr21);
    REAL sigmin = abs(detr) + abs(deti);
    result = sigmin / s1;
    //
    //     End of Rget53
    //
}
