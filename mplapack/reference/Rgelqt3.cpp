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

void Rgelqt3(INTEGER const &m, INTEGER const &n, REAL *a, INTEGER const &lda, REAL *t, INTEGER const &ldt, INTEGER &info) {
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
    //     .. Executable Statements ..
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < m) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    } else if (ldt < max((INTEGER)1, m)) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Rgelqt3", -info);
        return;
    }
    //
    INTEGER m1 = 0;
    INTEGER m2 = 0;
    INTEGER i1 = 0;
    INTEGER j1 = 0;
    INTEGER iinfo = 0;
    INTEGER i = 0;
    INTEGER j = 0;
    const REAL one = 1.00;
    if (m == 1) {
        //
        //        Compute Householder transform when N=1
        //
    Rlarfg(n, a, a[(min(2-1)*lda], lda, t);
    //
    } else {
        //
        //        Otherwise, split A INTEGERo blocks...
        //
        m1 = m / 2;
        m2 = m - m1;
        i1 = min(m1 + 1, m);
        j1 = min(m + 1, n);
        //
        //        Compute A(1:M1,1:N) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H
        //
        Rgelqt3(m1, n, a, lda, t, ldt, iinfo);
        //
        //        Compute A(J1:M,1:N) = Q1^H A(J1:M,1:N) [workspace: T(1:N1,J1:N)]
        //
        for (i = 1; i <= m2; i = i + 1) {
            for (j = 1; j <= m1; j = j + 1) {
                t[((i + m1) - 1) + (j - 1) * ldt] = a[((i + m1) - 1) + (j - 1) * lda];
            }
        }
        Rtrmm("R", "U", "T", "U", m2, m1, one, a, lda, t[(i1 - 1)], ldt);
        //
        Rgemm("N", "T", m2, m1, n - m1, one, a[(i1 - 1) + (i1 - 1) * lda], lda, a[(i1 - 1) * lda], lda, one, t[(i1 - 1)], ldt);
        //
        Rtrmm("R", "U", "N", "N", m2, m1, one, t, ldt, t[(i1 - 1)], ldt);
        //
        Rgemm("N", "N", m2, n - m1, m1, -one, t[(i1 - 1)], ldt, a[(i1 - 1) * lda], lda, one, a[(i1 - 1) + (i1 - 1) * lda], lda);
        //
        Rtrmm("R", "U", "N", "U", m2, m1, one, a, lda, t[(i1 - 1)], ldt);
        //
        for (i = 1; i <= m2; i = i + 1) {
            for (j = 1; j <= m1; j = j + 1) {
                a[((i + m1) - 1) + (j - 1) * lda] = a[((i + m1) - 1) + (j - 1) * lda] - t[((i + m1) - 1) + (j - 1) * ldt];
                t[((i + m1) - 1) + (j - 1) * ldt] = 0;
            }
        }
        //
        //        Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H
        //
        Rgelqt3(m2, n - m1, a[(i1 - 1) + (i1 - 1) * lda], lda, t[(i1 - 1) + (i1 - 1) * ldt], ldt, iinfo);
        //
        //        Compute T3 = T(J1:N1,1:N) = -T1 Y1^H Y2 T2
        //
        for (i = 1; i <= m2; i = i + 1) {
            for (j = 1; j <= m1; j = j + 1) {
                t[(j - 1) + ((i + m1) - 1) * ldt] = (a[(j - 1) + ((i + m1) - 1) * lda]);
            }
        }
        //
        Rtrmm("R", "U", "T", "U", m1, m2, one, a[(i1 - 1) + (i1 - 1) * lda], lda, t[(i1 - 1) * ldt], ldt);
        //
        Rgemm("N", "T", m1, m2, n - m, one, a[(j1 - 1) * lda], lda, a[(i1 - 1) + (j1 - 1) * lda], lda, one, t[(i1 - 1) * ldt], ldt);
        //
        Rtrmm("L", "U", "N", "N", m1, m2, -one, t, ldt, t[(i1 - 1) * ldt], ldt);
        //
        Rtrmm("R", "U", "N", "N", m1, m2, one, t[(i1 - 1) + (i1 - 1) * ldt], ldt, t[(i1 - 1) * ldt], ldt);
        //
        //        Y = (Y1,Y2); L = [ L1            0  ];  T = [T1 T3]
        //                         [ A(1:N1,J1:N)  L2 ]       [ 0 T2]
        //
    }
    //
    //     End of Rgelqt3
    //
}
