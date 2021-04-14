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

void Chptrd(const char *uplo, INTEGER const n, COMPLEX *ap, REAL *d, REAL *e, COMPLEX *tau, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    info = 0;
    bool upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    }
    if (info != 0) {
        Mxerbla("Chptrd", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    INTEGER i1 = 0;
    INTEGER i = 0;
    COMPLEX alpha = 0.0;
    COMPLEX taui = 0.0;
    const COMPLEX zero = (0.0, 0.0);
    const COMPLEX one = (1.0, 0.0);
    const COMPLEX half = (0.5e+0, 0.0);
    INTEGER ii = 0;
    INTEGER i1i1 = 0;
    if (upper) {
        //
        //        Reduce the upper triangle of A.
        //        I1 is the index in AP of A(1,I+1).
        //
        i1 = n * (n - 1) / 2 + 1;
        ap[(i1 + n - 1) - 1] = ap[(i1 + n - 1) - 1].real();
        for (i = n - 1; i >= 1; i = i - 1) {
            //
            //           Generate elementary reflector H(i) = I - tau * v * v**H
            //           to annihilate A(1:i-1,i+1)
            //
            alpha = ap[(i1 + i - 1) - 1];
            Clarfg(i, alpha, &ap[i1 - 1], 1, taui);
            e[i - 1] = alpha.real();
            //
            if (taui != zero) {
                //
                //              Apply H(i) from both sides to A(1:i,1:i)
                //
                ap[(i1 + i - 1) - 1] = one;
                //
                //              Compute  y := tau * A * v  storing y in TAU(1:i)
                //
                Chpmv(uplo, i, taui, ap, &ap[i1 - 1], 1, zero, tau, 1);
                //
                //              Compute  w := y - 1/2 * tau * (y**H *v) * v
                //
                alpha = -half * taui * Cdotc(i, tau, 1, &ap[i1 - 1], 1);
                Caxpy(i, alpha, &ap[i1 - 1], 1, tau, 1);
                //
                //              Apply the transformation as a rank-2 update:
                //                 A := A - v * w**H - w * v**H
                //
                Chpr2(uplo, i, -one, &ap[i1 - 1], 1, tau, 1, ap);
                //
            }
            ap[(i1 + i - 1) - 1] = e[i - 1];
            d[(i + 1) - 1] = ap[(i1 + i) - 1].real();
            tau[i - 1] = taui;
            i1 = i1 - i;
        }
        d[1 - 1] = ap[1 - 1].real();
    } else {
        //
        //        Reduce the lower triangle of A. II is the index in AP of
        //        A(i,i) and I1I1 is the index of A(i+1,i+1).
        //
        ii = 1;
        ap[1 - 1] = ap[1 - 1].real();
        for (i = 1; i <= n - 1; i = i + 1) {
            i1i1 = ii + n - i + 1;
            //
            //           Generate elementary reflector H(i) = I - tau * v * v**H
            //           to annihilate A(i+2:n,i)
            //
            alpha = ap[(ii + 1) - 1];
            Clarfg(n - i, alpha, &ap[(ii + 2) - 1], 1, taui);
            e[i - 1] = alpha.real();
            //
            if (taui != zero) {
                //
                //              Apply H(i) from both sides to A(i+1:n,i+1:n)
                //
                ap[(ii + 1) - 1] = one;
                //
                //              Compute  y := tau * A * v  storing y in TAU(i:n-1)
                //
                Chpmv(uplo, n - i, taui, &ap[i1i1 - 1], &ap[(ii + 1) - 1], 1, zero, &tau[i - 1], 1);
                //
                //              Compute  w := y - 1/2 * tau * (y**H *v) * v
                //
                alpha = -half * taui * Cdotc(n - i, &tau[i - 1], 1, &ap[(ii + 1) - 1], 1);
                Caxpy(n - i, alpha, &ap[(ii + 1) - 1], 1, &tau[i - 1], 1);
                //
                //              Apply the transformation as a rank-2 update:
                //                 A := A - v * w**H - w * v**H
                //
                Chpr2(uplo, n - i, -one, &ap[(ii + 1) - 1], 1, &tau[i - 1], 1, &ap[i1i1 - 1]);
                //
            }
            ap[(ii + 1) - 1] = e[i - 1];
            d[i - 1] = ap[ii - 1].real();
            tau[i - 1] = taui;
            ii = i1i1;
        }
        d[n - 1] = ap[ii - 1].real();
    }
    //
    //     End of Chptrd
    //
}
