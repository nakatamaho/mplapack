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

void Rlaed9(INTEGER const k, INTEGER const kstart, INTEGER const kstop, INTEGER const n, REAL *d, REAL *q, INTEGER const ldq, REAL const rho, REAL *dlamda, REAL *w, REAL *s, INTEGER const lds, INTEGER &info) {
    INTEGER i = 0;
    INTEGER j = 0;
    REAL temp = 0.0;
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
    //     Test the input parameters.
    //
    info = 0;
    //
    if (k < 0) {
        info = -1;
    } else if (kstart < 1 || kstart > max((INTEGER)1, k)) {
        info = -2;
    } else if (max((INTEGER)1, kstop) < kstart || kstop > max((INTEGER)1, k)) {
        info = -3;
    } else if (n < k) {
        info = -4;
    } else if (ldq < max((INTEGER)1, k)) {
        info = -7;
    } else if (lds < max((INTEGER)1, k)) {
        info = -12;
    }
    if (info != 0) {
        Mxerbla("Rlaed9", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (k == 0) {
        return;
    }
    //
    //     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can
    //     be computed with high relative accuracy (barring over/underflow).
    //     This is a problem on machines without a guard digit in
    //     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
    //     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I),
    //     which on any of these machines zeros out the bottommost
    //     bit of DLAMDA(I) if it is 1; this makes the subsequent
    //     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation
    //     occurs. On binary machines with a guard digit (almost all
    //     machines) it does not change DLAMDA(I) at all. On hexadecimal
    //     and decimal machines with a guard digit, it slightly
    //     changes the bottommost bits of DLAMDA(I). It does not account
    //     for hexadecimal or decimal machines without guard digits
    //     (we know of none). We use a subroutine call to compute
    //     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
    //     this code.
    //
    for (i = 1; i <= n; i = i + 1) {
        dlamda[i - 1] = Rlamc3(dlamda[i - 1], (dlamda[i - 1] - 1)) - dlamda[i - 1];
    }
    //
    for (j = kstart; j <= kstop; j = j + 1) {
        Rlaed4(k, j, dlamda, w, &q[(j - 1) * ldq], rho, d[j - 1], info);
        //
        //        If the zero finder fails, the computation is terminated.
        //
        if (info != 0) {
            goto statement_120;
        }
    }
    //
    if (k == 1 || k == 2) {
        for (i = 1; i <= k; i = i + 1) {
            for (j = 1; j <= k; j = j + 1) {
                s[(j - 1) + (i - 1) * lds] = q[(j - 1) + (i - 1) * ldq];
            }
        }
        goto statement_120;
    }
    //
    //     Compute updated W.
    //
    Rcopy(k, w, 1, s, 1);
    //
    //     Initialize W(I) = Q(I,I)
    //
    Rcopy(k, q, ldq + 1, w, 1);
    for (j = 1; j <= k; j = j + 1) {
        for (i = 1; i <= j - 1; i = i + 1) {
            w[i - 1] = w[i - 1] * (q[(i - 1) + (j - 1) * ldq] / (dlamda[i - 1] - dlamda[j - 1]));
        }
        for (i = j + 1; i <= k; i = i + 1) {
            w[i - 1] = w[i - 1] * (q[(i - 1) + (j - 1) * ldq] / (dlamda[i - 1] - dlamda[j - 1]));
        }
    }
    for (i = 1; i <= k; i = i + 1) {
        w[i - 1] = sign(sqrt(-w[i - 1]), s[(i - 1)]);
    }
    //
    //     Compute eigenvectors of the modified rank-1 modification.
    //
    for (j = 1; j <= k; j = j + 1) {
        for (i = 1; i <= k; i = i + 1) {
            q[(i - 1) + (j - 1) * ldq] = w[i - 1] / q[(i - 1) + (j - 1) * ldq];
        }
        temp = Rnrm2(k, &q[(j - 1) * ldq], 1);
        for (i = 1; i <= k; i = i + 1) {
            s[(i - 1) + (j - 1) * lds] = q[(i - 1) + (j - 1) * ldq] / temp;
        }
    }
//
statement_120:;
    //
    //     End of Rlaed9
    //
}
