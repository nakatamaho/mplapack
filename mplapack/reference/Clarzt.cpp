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

void Clarzt(const char *direct, const char *storev, INTEGER const n, INTEGER const k, COMPLEX *v, INTEGER const ldv, COMPLEX *tau, COMPLEX *t, INTEGER const ldt) {
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
    //     .. Executable Statements ..
    //
    //     Check for currently supported options
    //
    INTEGER info = 0;
    if (!Mlsame(direct, "B")) {
        info = -1;
    } else if (!Mlsame(storev, "R")) {
        info = -2;
    }
    if (info != 0) {
        Mxerbla("Clarzt", -info);
        return;
    }
    //
    INTEGER i = 0;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    INTEGER j = 0;
    for (i = k; i >= 1; i = i - 1) {
        if (tau[i - 1] == zero) {
            //
            //           H(i)  =  I
            //
            for (j = i; j <= k; j = j + 1) {
                t[(j - 1) + (i - 1) * ldt] = zero;
            }
        } else {
            //
            //           general case
            //
            if (i < k) {
                //
                //              T(i+1:k,i) = - tau(i) * V(i+1:k,1:n) * V(i,1:n)**H
                //
                Clacgv(n, &v[(i - 1)], ldv);
                Cgemv("No transpose", k - i, n, -tau[i - 1], &v[((i + 1) - 1)], ldv, &v[(i - 1)], ldv, zero, &t[((i + 1) - 1) + (i - 1) * ldt], 1);
                Clacgv(n, &v[(i - 1)], ldv);
                //
                //              T(i+1:k,i) = T(i+1:k,i+1:k) * T(i+1:k,i)
                //
                Ctrmv("Lower", "No transpose", "Non-unit", k - i, &t[((i + 1) - 1) + ((i + 1) - 1) * ldt], ldt, &t[((i + 1) - 1) + (i - 1) * ldt], 1);
            }
            t[(i - 1) + (i - 1) * ldt] = tau[i - 1];
        }
    }
    //
    //     End of Clarzt
    //
}
