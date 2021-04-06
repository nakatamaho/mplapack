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

void Rlarz(const char *side, INTEGER const m, INTEGER const n, INTEGER const l, REAL *v, INTEGER const incv, REAL const tau, REAL *c, INTEGER const ldc, REAL *work) {
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
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if (Mlsame(side, "L")) {
        //
        //        Form  H * C
        //
        if (tau != zero) {
            //
            //           w( 1:n ) = C( 1, 1:n )
            //
            Rcopy(n, c, ldc, work, 1);
            //
            //           w( 1:n ) = w( 1:n ) + C( m-l+1:m, 1:n )**T * v( 1:l )
            //
            Rgemv("Transpose", l, n, one, &c[((m - l + 1) - 1)], ldc, v, incv, one, work, 1);
            //
            //           C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )
            //
            Raxpy(n, -tau, work, 1, c, ldc);
            //
            //           C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
            //                               tau * v( 1:l ) * w( 1:n )**T
            //
            Rger(l, n, -tau, v, incv, work, 1, &c[((m - l + 1) - 1)], ldc);
        }
        //
    } else {
        //
        //        Form  C * H
        //
        if (tau != zero) {
            //
            //           w( 1:m ) = C( 1:m, 1 )
            //
            Rcopy(m, c, 1, work, 1);
            //
            //           w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )
            //
            Rgemv("No transpose", m, l, one, &c[((n - l + 1) - 1) * ldc], ldc, v, incv, one, work, 1);
            //
            //           C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )
            //
            Raxpy(m, -tau, work, 1, c, 1);
            //
            //           C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
            //                               tau * w( 1:m ) * v( 1:l )**T
            //
            Rger(m, l, -tau, work, 1, v, incv, &c[((n - l + 1) - 1) * ldc], ldc);
            //
        }
        //
    }
    //
    //     End of Rlarz
    //
}
