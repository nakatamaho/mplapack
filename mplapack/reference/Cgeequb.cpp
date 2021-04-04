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

void Cgeequb(INTEGER const &m, INTEGER const &n, COMPLEX *a, INTEGER const &lda, REAL *r, REAL *c, REAL &rowcnd, REAL &colcnd, REAL &amax, INTEGER &info) {
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
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    COMPLEX zdum = 0.0;
    abs1[zdum - 1] = abs(zdum.real()) + abs(zdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Cgeequb", -info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const REAL one = 1.0;
    const REAL zero = 0.0;
    if (m == 0 || n == 0) {
        rowcnd = one;
        colcnd = one;
        amax = zero;
        return;
    }
    //
    //     Get machine constants.  Assume SMLNUM is a power of the radix.
    //
    REAL smlnum = dlamch("S");
    REAL bignum = one / smlnum;
    REAL radix = dlamch("B");
    REAL logrdx = log[radix - 1];
    //
    //     Compute row scale factors.
    //
    INTEGER i = 0;
    for (i = 1; i <= m; i = i + 1) {
        r[i - 1] = zero;
    }
    //
    //     Find the maximum element in each row.
    //
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            r[i - 1] = max(r[i - 1], abs1[a[(i - 1) + (j - 1) * lda] - 1]);
        }
    }
    for (i = 1; i <= m; i = i + 1) {
        if (r[i - 1] > zero) {
            r[i - 1] = pow(radix, INTEGER(log[r[i - 1] - 1] / logrdx));
        }
    }
    //
    //     Find the maximum and minimum scale factors.
    //
    REAL rcmin = bignum;
    REAL rcmax = zero;
    for (i = 1; i <= m; i = i + 1) {
        rcmax = max(rcmax, r[i - 1]);
        rcmin = min(rcmin, r[i - 1]);
    }
    amax = rcmax;
    //
    if (rcmin == zero) {
        //
        //        Find the first zero scale factor and return an error code.
        //
        for (i = 1; i <= m; i = i + 1) {
            if (r[i - 1] == zero) {
                info = i;
                return;
            }
        }
    } else {
        //
        //        Invert the scale factors.
        //
        for (i = 1; i <= m; i = i + 1) {
            r[i - 1] = one / min(max(r[i - 1], smlnum), bignum);
        }
        //
        //        Compute ROWCND = min(R(I)) / max(R(I)).
        //
        rowcnd = max(rcmin, smlnum) / min(rcmax, bignum);
    }
    //
    //     Compute column scale factors.
    //
    for (j = 1; j <= n; j = j + 1) {
        c[j - 1] = zero;
    }
    //
    //     Find the maximum element in each column,
    //     assuming the row scaling computed above.
    //
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= m; i = i + 1) {
            c[j - 1] = max(c[j - 1], abs1[a[(i - 1) + (j - 1) * lda] - 1] * r[i - 1]);
        }
        if (c[j - 1] > zero) {
            c[j - 1] = pow(radix, INTEGER(log[c[j - 1] - 1] / logrdx));
        }
    }
    //
    //     Find the maximum and minimum scale factors.
    //
    rcmin = bignum;
    rcmax = zero;
    for (j = 1; j <= n; j = j + 1) {
        rcmin = min(rcmin, c[j - 1]);
        rcmax = max(rcmax, c[j - 1]);
    }
    //
    if (rcmin == zero) {
        //
        //        Find the first zero scale factor and return an error code.
        //
        for (j = 1; j <= n; j = j + 1) {
            if (c[j - 1] == zero) {
                info = m + j;
                return;
            }
        }
    } else {
        //
        //        Invert the scale factors.
        //
        for (j = 1; j <= n; j = j + 1) {
            c[j - 1] = one / min(max(c[j - 1], smlnum), bignum);
        }
        //
        //        Compute COLCND = min(C(J)) / max(C(J)).
        //
        colcnd = max(rcmin, smlnum) / min(rcmax, bignum);
    }
    //
    //     End of Cgeequb
    //
}
