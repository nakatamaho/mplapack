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

void Cgeqrt3(INTEGER const &m, INTEGER const &n, COMPLEX *a, INTEGER const &lda, COMPLEX *t, INTEGER const &ldt, INTEGER &info) {
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
    if (n < 0) {
        info = -2;
    } else if (m < n) {
        info = -1;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    } else if (ldt < max((INTEGER)1, n)) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Cgeqrt3", -info);
        return;
    }
    //
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER j1 = 0;
    INTEGER i1 = 0;
    INTEGER iinfo = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    const COMPLEX one = (1.00, 0.00);
    if (n == 1) {
        //
        //        Compute Householder transform when N=1
        //
        Clarfg(m, a[(1 - 1)], a[(min(2 - 1) + (m)-1) * lda], 1, t[(1 - 1)]);
        //
    } else {
        //
        //        Otherwise, split A INTEGERo blocks...
        //
        n1 = n / 2;
        n2 = n - n1;
        j1 = min(n1 + 1, n);
        i1 = min(n + 1, m);
        //
        //        Compute A(1:M,1:N1) <- (Y1,R1,T1), where Q1 = I - Y1 T1 Y1^H
        //
        Cgeqrt3(m, n1, a, lda, t, ldt, iinfo);
        //
        //        Compute A(1:M,J1:N) = Q1^H A(1:M,J1:N) [workspace: T(1:N1,J1:N)]
        //
        for (j = 1; j <= n2; j = j + 1) {
            for (i = 1; i <= n1; i = i + 1) {
                t[(i - 1) + ((j + n1) - 1) * ldt] = a[(i - 1) + ((j + n1) - 1) * lda];
            }
        }
        Ctrmm("L", "L", "C", "U", n1, n2, one, a, lda, t[(j1 - 1) * ldt], ldt);
        //
        Cgemm("C", "N", n1, n2, m - n1, one, a[(j1 - 1)], lda, a[(j1 - 1) + (j1 - 1) * lda], lda, one, t[(j1 - 1) * ldt], ldt);
        //
        Ctrmm("L", "U", "C", "N", n1, n2, one, t, ldt, t[(j1 - 1) * ldt], ldt);
        //
        Cgemm("N", "N", m - n1, n2, n1, -one, a[(j1 - 1)], lda, t[(j1 - 1) * ldt], ldt, one, a[(j1 - 1) + (j1 - 1) * lda], lda);
        //
        Ctrmm("L", "L", "N", "U", n1, n2, one, a, lda, t[(j1 - 1) * ldt], ldt);
        //
        for (j = 1; j <= n2; j = j + 1) {
            for (i = 1; i <= n1; i = i + 1) {
                a[(i - 1) + ((j + n1) - 1) * lda] = a[(i - 1) + ((j + n1) - 1) * lda] - t[(i - 1) + ((j + n1) - 1) * ldt];
            }
        }
        //
        //        Compute A(J1:M,J1:N) <- (Y2,R2,T2) where Q2 = I - Y2 T2 Y2^H
        //
        Cgeqrt3(m - n1, n2, a[(j1 - 1) + (j1 - 1) * lda], lda, t[(j1 - 1) + (j1 - 1) * ldt], ldt, iinfo);
        //
        //        Compute T3 = T(1:N1,J1:N) = -T1 Y1^H Y2 T2
        //
        for (i = 1; i <= n1; i = i + 1) {
            for (j = 1; j <= n2; j = j + 1) {
                t[(i - 1) + ((j + n1) - 1) * ldt] = conjg[(a[((j + n1) - 1) + (i - 1) * lda]) - 1];
            }
        }
        //
        Ctrmm("R", "L", "N", "U", n1, n2, one, a[(j1 - 1) + (j1 - 1) * lda], lda, t[(j1 - 1) * ldt], ldt);
        //
        Cgemm("C", "N", n1, n2, m - n, one, a[(i1 - 1)], lda, a[(i1 - 1) + (j1 - 1) * lda], lda, one, t[(j1 - 1) * ldt], ldt);
        //
        Ctrmm("L", "U", "N", "N", n1, n2, -one, t, ldt, t[(j1 - 1) * ldt], ldt);
        //
        Ctrmm("R", "U", "N", "N", n1, n2, one, t[(j1 - 1) + (j1 - 1) * ldt], ldt, t[(j1 - 1) * ldt], ldt);
        //
        //        Y = (Y1,Y2); R = [ R1  A(1:N1,J1:N) ];  T = [T1 T3]
        //                         [  0        R2     ]       [ 0 T2]
        //
    }
    //
    //     End of Cgeqrt3
    //
}
