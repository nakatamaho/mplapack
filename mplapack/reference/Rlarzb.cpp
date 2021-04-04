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

void Rlarzb(const char *side, const char *trans, const char *direct, const char *storev, INTEGER const &m, INTEGER const &n, INTEGER const &k, INTEGER const &l, REAL *v, INTEGER const &ldv, REAL *t, INTEGER const &ldt, REAL *c, INTEGER const &ldc, REAL *work, INTEGER const &ldwork) {
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
    //     .. Executable Statements ..
    //
    //     Quick return if possible
    //
    if (m <= 0 || n <= 0) {
        return;
    }
    //
    //     Check for currently supported options
    //
    INTEGER info = 0;
    if (!Mlsame(direct, "B")) {
        info = -3;
    } else if (!Mlsame(storev, "R")) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Rlarzb", -info);
        return;
    }
    //
    str<1> transt = char0;
    if (Mlsame(trans, "N")) {
        transt = "T";
    } else {
        transt = "N";
    }
    //
    INTEGER j = 0;
    const REAL one = 1.0;
    INTEGER i = 0;
    if (Mlsame(side, "L")) {
        //
        //        Form  H * C  or  H**T * C
        //
        //        W( 1:n, 1:k ) = C( 1:k, 1:n )**T
        //
        for (j = 1; j <= k; j = j + 1) {
            Rcopy(n, c[(j - 1)], ldc, work[(j - 1) * ldwork], 1);
        }
        //
        //        W( 1:n, 1:k ) = W( 1:n, 1:k ) + ...
        //                        C( m-l+1:m, 1:n )**T * V( 1:k, 1:l )**T
        //
        if (l > 0) {
            Rgemm("Transpose", "Transpose", n, k, l, one, c[((m - l + 1) - 1)], ldc, v, ldv, one, work, ldwork);
        }
        //
        //        W( 1:n, 1:k ) = W( 1:n, 1:k ) * T**T  or  W( 1:m, 1:k ) * T
        //
        Rtrmm("Right", "Lower", transt, "Non-unit", n, k, one, t, ldt, work, ldwork);
        //
        //        C( 1:k, 1:n ) = C( 1:k, 1:n ) - W( 1:n, 1:k )**T
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - work[(j - 1) + (i - 1) * ldwork];
            }
        }
        //
        //        C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
        //                            V( 1:k, 1:l )**T * W( 1:n, 1:k )**T
        //
        if (l > 0) {
            Rgemm("Transpose", "Transpose", l, n, k, -one, v, ldv, work, ldwork, one, c[((m - l + 1) - 1)], ldc);
        }
        //
    } else if (Mlsame(side, "R")) {
        //
        //        Form  C * H  or  C * H**T
        //
        //        W( 1:m, 1:k ) = C( 1:m, 1:k )
        //
        for (j = 1; j <= k; j = j + 1) {
            Rcopy(m, c[(j - 1) * ldc], 1, work[(j - 1) * ldwork], 1);
        }
        //
        //        W( 1:m, 1:k ) = W( 1:m, 1:k ) + ...
        //                        C( 1:m, n-l+1:n ) * V( 1:k, 1:l )**T
        //
        if (l > 0) {
            Rgemm("No transpose", "Transpose", m, k, l, one, c[((n - l + 1) - 1) * ldc], ldc, v, ldv, one, work, ldwork);
        }
        //
        //        W( 1:m, 1:k ) = W( 1:m, 1:k ) * T  or  W( 1:m, 1:k ) * T**T
        //
        Rtrmm("Right", "Lower", trans, "Non-unit", m, k, one, t, ldt, work, ldwork);
        //
        //        C( 1:m, 1:k ) = C( 1:m, 1:k ) - W( 1:m, 1:k )
        //
        for (j = 1; j <= k; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                c[(i - 1) + (j - 1) * ldc] = c[(i - 1) + (j - 1) * ldc] - work[(i - 1) + (j - 1) * ldwork];
            }
        }
        //
        //        C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
        //                            W( 1:m, 1:k ) * V( 1:k, 1:l )
        //
        if (l > 0) {
            Rgemm("No transpose", "No transpose", m, l, k, -one, work, ldwork, v, ldv, one, c[((n - l + 1) - 1) * ldc], ldc);
        }
        //
    }
    //
    //     End of Rlarzb
    //
}
