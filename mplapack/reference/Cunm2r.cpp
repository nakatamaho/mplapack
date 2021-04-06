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

void Cunm2r(const char *side, const char *trans, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *tau, COMPLEX *c, INTEGER const ldc, COMPLEX *work, INTEGER &info) {
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
    bool left = Mlsame(side, "L");
    bool notran = Mlsame(trans, "N");
    //
    //     NQ is the order of Q
    //
    INTEGER nq = 0;
    if (left) {
        nq = m;
    } else {
        nq = n;
    }
    if (!left && !Mlsame(side, "R")) {
        info = -1;
    } else if (!notran && !Mlsame(trans, "C")) {
        info = -2;
    } else if (m < 0) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (k < 0 || k > nq) {
        info = -5;
    } else if (lda < max((INTEGER)1, nq)) {
        info = -7;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -10;
    }
    if (info != 0) {
        Mxerbla("Cunm2r", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0 || k == 0) {
        return;
    }
    //
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER i3 = 0;
    if ((left && !notran || !left && notran)) {
        i1 = 1;
        i2 = k;
        i3 = 1;
    } else {
        i1 = k;
        i2 = 1;
        i3 = -1;
    }
    //
    INTEGER ni = 0;
    INTEGER jc = 0;
    INTEGER mi = 0;
    INTEGER ic = 0;
    if (left) {
        ni = n;
        jc = 1;
    } else {
        mi = m;
        ic = 1;
    }
    //
    INTEGER i = 0;
    COMPLEX taui = 0.0;
    COMPLEX aii = 0.0;
    const COMPLEX one = (1.0, 0.0);
    for (i = i1; i <= i2; i = i + i3) {
        if (left) {
            //
            //           H(i) or H(i)**H is applied to C(i:m,1:n)
            //
            mi = m - i + 1;
            ic = i;
        } else {
            //
            //           H(i) or H(i)**H is applied to C(1:m,i:n)
            //
            ni = n - i + 1;
            jc = i;
        }
        //
        //        Apply H(i) or H(i)**H
        //
        if (notran) {
            taui = tau[i - 1];
        } else {
            taui = conj(tau[i - 1]);
        }
        aii = a[(i - 1) + (i - 1) * lda];
        a[(i - 1) + (i - 1) * lda] = one;
        Clarf(side, mi, ni, &a[(i - 1) + (i - 1) * lda], 1, taui, &c[(ic - 1) + (jc - 1) * ldc], ldc, work);
        a[(i - 1) + (i - 1) * lda] = aii;
    }
    //
    //     End of Cunm2r
    //
}
