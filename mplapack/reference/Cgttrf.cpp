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

inline REAL abs1(COMPLEX zdum) { return abs(zdum.real()) + abs(zdum.imag()); }

void Cgttrf(INTEGER const n, COMPLEX *dl, COMPLEX *d, COMPLEX *du, COMPLEX *du2, INTEGER *ipiv, INTEGER &info) {
    COMPLEX zdum = 0.0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    COMPLEX fact = 0.0;
    COMPLEX temp = 0.0;
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
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    if (n < 0) {
        info = -1;
        Mxerbla("Cgttrf", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Initialize IPIV(i) = i and DU2(i) = 0
    //
    for (i = 1; i <= n; i = i + 1) {
        ipiv[i - 1] = i;
    }
    for (i = 1; i <= n - 2; i = i + 1) {
        du2[i - 1] = zero;
    }
    //
    for (i = 1; i <= n - 2; i = i + 1) {
        if (abs1(d[i - 1]) >= abs1(dl[i - 1])) {
            //
            //           No row interchange required, eliminate DL(I)
            //
            if (abs1(d[i - 1]) != zero) {
                fact = dl[i - 1] / d[i - 1];
                dl[i - 1] = fact;
                d[(i + 1) - 1] = d[(i + 1) - 1] - fact * du[i - 1];
            }
        } else {
            //
            //           Interchange rows I and I+1, eliminate DL(I)
            //
            fact = d[i - 1] / dl[i - 1];
            d[i - 1] = dl[i - 1];
            dl[i - 1] = fact;
            temp = du[i - 1];
            du[i - 1] = d[(i + 1) - 1];
            d[(i + 1) - 1] = temp - fact * d[(i + 1) - 1];
            du2[i - 1] = du[(i + 1) - 1];
            du[(i + 1) - 1] = -fact * du[(i + 1) - 1];
            ipiv[i - 1] = i + 1;
        }
    }
    if (n > 1) {
        i = n - 1;
        if (abs1(d[i - 1]) >= abs1(dl[i - 1])) {
            if (abs1(d[i - 1]) != zero) {
                fact = dl[i - 1] / d[i - 1];
                dl[i - 1] = fact;
                d[(i + 1) - 1] = d[(i + 1) - 1] - fact * du[i - 1];
            }
        } else {
            fact = d[i - 1] / dl[i - 1];
            d[i - 1] = dl[i - 1];
            dl[i - 1] = fact;
            temp = du[i - 1];
            du[i - 1] = d[(i + 1) - 1];
            d[(i + 1) - 1] = temp - fact * d[(i + 1) - 1];
            ipiv[i - 1] = i + 1;
        }
    }
    //
    //     Check for a zero on the diagonal of U.
    //
    for (i = 1; i <= n; i = i + 1) {
        if (abs1(d[i - 1]) == zero) {
            info = i;
            goto statement_50;
        }
    }
statement_50:;
    //
    //     End of Cgttrf
    //
}
