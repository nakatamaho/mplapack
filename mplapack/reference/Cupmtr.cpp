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

void Cupmtr(const char *side, const char *uplo, const char *trans, INTEGER const &m, INTEGER const &n, COMPLEX *ap, COMPLEX *tau, COMPLEX *c, INTEGER const &ldc, COMPLEX *work, INTEGER &info) {
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
    bool upper = Mlsame(uplo, "U");
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
    } else if (!upper && !Mlsame(uplo, "L")) {
        info = -2;
    } else if (!notran && !Mlsame(trans, "C")) {
        info = -3;
    } else if (m < 0) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Cupmtr", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    bool forwrd = false;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER i3 = 0;
    INTEGER ii = 0;
    INTEGER ni = 0;
    INTEGER mi = 0;
    INTEGER i = 0;
    COMPLEX taui = 0.0;
    COMPLEX aii = 0.0;
    const COMPLEX one = (1.0, 0.0);
    INTEGER jc = 0;
    INTEGER ic = 0;
    if (upper) {
        //
        //        Q was determined by a call to Chptrd with UPLO = 'U'
        //
        forwrd = (left && notran) || (!left && !notran);
        //
        if (forwrd) {
            i1 = 1;
            i2 = nq - 1;
            i3 = 1;
            ii = 2;
        } else {
            i1 = nq - 1;
            i2 = 1;
            i3 = -1;
            ii = nq * (nq + 1) / 2 - 1;
        }
        //
        if (left) {
            ni = n;
        } else {
            mi = m;
        }
        //
        for (i = i1; i <= i2; i = i + i3) {
            if (left) {
                //
                //              H(i) or H(i)**H is applied to C(1:i,1:n)
                //
                mi = i;
            } else {
                //
                //              H(i) or H(i)**H is applied to C(1:m,1:i)
                //
                ni = i;
            }
            //
            //           Apply H(i) or H(i)**H
            //
            if (notran) {
                taui = tau[i - 1];
            } else {
                taui = conj(tau[i - 1]);
            }
            aii = ap[ii - 1];
            ap[ii - 1] = one;
            Clarf(side, mi, ni, ap[(ii - i + 1) - 1], 1, taui, c, ldc, work);
            ap[ii - 1] = aii;
            //
            if (forwrd) {
                ii += i + 2;
            } else {
                ii = ii - i - 1;
            }
        }
    } else {
        //
        //        Q was determined by a call to Chptrd with UPLO = 'L'.
        //
        forwrd = (left && !notran) || (!left && notran);
        //
        if (forwrd) {
            i1 = 1;
            i2 = nq - 1;
            i3 = 1;
            ii = 2;
        } else {
            i1 = nq - 1;
            i2 = 1;
            i3 = -1;
            ii = nq * (nq + 1) / 2 - 1;
        }
        //
        if (left) {
            ni = n;
            jc = 1;
        } else {
            mi = m;
            ic = 1;
        }
        //
        for (i = i1; i <= i2; i = i + i3) {
            aii = ap[ii - 1];
            ap[ii - 1] = one;
            if (left) {
                //
                //              H(i) or H(i)**H is applied to C(i+1:m,1:n)
                //
                mi = m - i;
                ic = i + 1;
            } else {
                //
                //              H(i) or H(i)**H is applied to C(1:m,i+1:n)
                //
                ni = n - i;
                jc = i + 1;
            }
            //
            //           Apply H(i) or H(i)**H
            //
            if (notran) {
                taui = tau[i - 1];
            } else {
                taui = conj(tau[i - 1]);
            }
            Clarf(side, mi, ni, ap[ii - 1], 1, taui, c[(ic - 1) + (jc - 1) * ldc], ldc, work);
            ap[ii - 1] = aii;
            //
            if (forwrd) {
                ii += nq - i + 1;
            } else {
                ii = ii - nq + i - 2;
            }
        }
    }
    //
    //     End of Cupmtr
    //
}
