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

void Cgtsv(INTEGER const &n, INTEGER const &nrhs, COMPLEX *dl, COMPLEX *d, COMPLEX *du, COMPLEX *b, INTEGER const &ldb, INTEGER &info) {
    //
    //  -- LAPACK driver routine --
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    COMPLEX zdum = 0.0;
    abs1[zdum - 1] = abs(zdum.real()) + abs(zdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    if (n < 0) {
        info = -1;
    } else if (nrhs < 0) {
        info = -2;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -7;
    }
    if (info != 0) {
        Mxerbla("Cgtsv ", -info);
        return;
    }
    //
    if (n == 0) {
        return;
    }
    //
    INTEGER k = 0;
    const COMPLEX zero = (0.0, 0.0);
    COMPLEX mult = 0.0;
    INTEGER j = 0;
    COMPLEX temp = 0.0;
    for (k = 1; k <= n - 1; k = k + 1) {
        if (dl[k - 1] == zero) {
            //
            //           Subdiagonal is zero, no elimination is required.
            //
            if (d[k - 1] == zero) {
                //
                //              Diagonal is zero: set INFO = K and return; a unique
                //              solution can not be found.
                //
                info = k;
                return;
            }
        } else if (abs1[d[k - 1] - 1] >= abs1[dl[k - 1] - 1]) {
            //
            //           No row INTEGERerchange required
            //
            mult = dl[k - 1] / d[k - 1];
            d[(k + 1) - 1] = d[(k + 1) - 1] - mult * du[k - 1];
            for (j = 1; j <= nrhs; j = j + 1) {
                b[((k + 1) - 1) + (j - 1) * ldb] = b[((k + 1) - 1) + (j - 1) * ldb] - mult * b[(k - 1) + (j - 1) * ldb];
            }
            if (k < (n - 1)) {
                dl[k - 1] = zero;
            }
        } else {
            //
            //           Interchange rows K and K+1
            //
            mult = d[k - 1] / dl[k - 1];
            d[k - 1] = dl[k - 1];
            temp = d[(k + 1) - 1];
            d[(k + 1) - 1] = du[k - 1] - mult * temp;
            if (k < (n - 1)) {
                dl[k - 1] = du[(k + 1) - 1];
                du[(k + 1) - 1] = -mult * dl[k - 1];
            }
            du[k - 1] = temp;
            for (j = 1; j <= nrhs; j = j + 1) {
                temp = b[(k - 1) + (j - 1) * ldb];
                b[(k - 1) + (j - 1) * ldb] = b[((k + 1) - 1) + (j - 1) * ldb];
                b[((k + 1) - 1) + (j - 1) * ldb] = temp - mult * b[((k + 1) - 1) + (j - 1) * ldb];
            }
        }
    }
    if (d[n - 1] == zero) {
        info = n;
        return;
    }
    //
    //     Back solve with the matrix U from the factorization.
    //
    for (j = 1; j <= nrhs; j = j + 1) {
        b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] / d[n - 1];
        if (n > 1) {
            b[((n - 1) - 1) + (j - 1) * ldb] = (b[((n - 1) - 1) + (j - 1) * ldb] - du[(n - 1) - 1] * b[(n - 1) + (j - 1) * ldb]) / d[(n - 1) - 1];
        }
        for (k = n - 2; k >= 1; k = k - 1) {
            b[(k - 1) + (j - 1) * ldb] = (b[(k - 1) + (j - 1) * ldb] - du[k - 1] * b[((k + 1) - 1) + (j - 1) * ldb] - dl[k - 1] * b[((k + 2) - 1) + (j - 1) * ldb]) / d[k - 1];
        }
    }
    //
    //     End of Cgtsv
    //
}
