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

void Cpteqr(const char *compz, INTEGER const &n, REAL *d, REAL *e, COMPLEX *z, INTEGER const &ldz, REAL *work, INTEGER &info) {
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
    //  ====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    //
    INTEGER icompz = 0;
    if (Mlsame(compz, "N")) {
        icompz = 0;
    } else if (Mlsame(compz, "V")) {
        icompz = 1;
    } else if (Mlsame(compz, "I")) {
        icompz = 2;
    } else {
        icompz = -1;
    }
    if (icompz < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if ((ldz < 1) || (icompz > 0 && ldz < max((INTEGER)1, n))) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Cpteqr", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    const COMPLEX cone = (1.0, 0.0);
    if (n == 1) {
        if (icompz > 0) {
            z[(1 - 1)] = cone;
        }
        return;
    }
    const COMPLEX czero = (0.0, 0.0);
    if (icompz == 2) {
        Claset("Full", n, n, czero, cone, z, ldz);
    }
    //
    //     Call Rpttrf to factor the matrix.
    //
    Rpttrf(n, d, e, info);
    if (info != 0) {
        return;
    }
    INTEGER i = 0;
    for (i = 1; i <= n; i = i + 1) {
        d[i - 1] = sqrt(d[i - 1]);
    }
    for (i = 1; i <= n - 1; i = i + 1) {
        e[i - 1] = e[i - 1] * d[i - 1];
    }
    //
    //     Call Cbdsqr to compute the singular values/vectors of the
    //     bidiagonal factor.
    //
    INTEGER nru = 0;
    if (icompz > 0) {
        nru = n;
    } else {
        nru = 0;
    }
    arr_2d<1, 1, COMPLEX> vt(fill0);
    arr_2d<1, 1, COMPLEX> c(fill0);
    Cbdsqr("Lower", n, 0, nru, 0, d, e, vt, 1, z, ldz, c, 1, work, info);
    //
    //     Square the singular values.
    //
    if (info == 0) {
        for (i = 1; i <= n; i = i + 1) {
            d[i - 1] = d[i - 1] * d[i - 1];
        }
    } else {
        info += n;
    }
    //
    //     End of Cpteqr
    //
}
