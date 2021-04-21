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

void Cupgtr(const char *uplo, INTEGER const n, COMPLEX *ap, COMPLEX *tau, COMPLEX *q, INTEGER const ldq, COMPLEX *work, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    bool upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (ldq < max((INTEGER)1, n)) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Cupgtr", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    INTEGER ij = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER iinfo = 0;
    if (upper) {
        //
        //        Q was determined by a call to Chptrd with UPLO = 'U'
        //
        //        Unpack the vectors which define the elementary reflectors and
        //        set the last row and column of Q equal to those of the unit
        //        matrix
        //
        ij = 2;
        for (j = 1; j <= n - 1; j = j + 1) {
            for (i = 1; i <= j - 1; i = i + 1) {
                q[(i - 1) + (j - 1) * ldq] = ap[ij - 1];
                ij++;
            }
            ij += 2;
            q[(n - 1) + (j - 1) * ldq] = czero;
        }
        for (i = 1; i <= n - 1; i = i + 1) {
            q[(i - 1) + (n - 1) * ldq] = czero;
        }
        q[(n - 1) + (n - 1) * ldq] = cone;
        //
        //        Generate Q(1:n-1,1:n-1)
        //
        Cung2l(n - 1, n - 1, n - 1, q, ldq, tau, work, iinfo);
        //
    } else {
        //
        //        Q was determined by a call to Chptrd with UPLO = 'L'.
        //
        //        Unpack the vectors which define the elementary reflectors and
        //        set the first row and column of Q equal to those of the unit
        //        matrix
        //
        q[(1 - 1)] = cone;
        for (i = 2; i <= n; i = i + 1) {
            q[(i - 1)] = czero;
        }
        ij = 3;
        for (j = 2; j <= n; j = j + 1) {
            q[(j - 1) * ldq] = czero;
            for (i = j + 1; i <= n; i = i + 1) {
                q[(i - 1) + (j - 1) * ldq] = ap[ij - 1];
                ij++;
            }
            ij += 2;
        }
        if (n > 1) {
            //
            //           Generate Q(2:n,2:n)
            //
            Cung2r(n - 1, n - 1, n - 1, &q[(2 - 1) + (2 - 1) * ldq], ldq, tau, work, iinfo);
        }
    }
    //
    //     End of Cupgtr
    //
}
