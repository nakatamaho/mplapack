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

void Rtplqt2(INTEGER const &m, INTEGER const &n, INTEGER const &l, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL *t, INTEGER const &ldt, INTEGER &info) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (l < 0 || l > min(m, n)) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, m)) {
        info = -7;
    } else if (ldt < max((INTEGER)1, m)) {
        info = -9;
    }
    if (info != 0) {
        Mxerbla("Rtplqt2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0 || m == 0) {
        return;
    }
    //
    INTEGER i = 0;
    INTEGER p = 0;
    INTEGER j = 0;
    const REAL one = 1.0f;
    REAL alpha = 0.0;
    for (i = 1; i <= m; i = i + 1) {
        //
        //        Generate elementary reflector H(I) to annihilate B(I,:)
        //
        p = n - l + min(l, i);
        Rlarfg(p + 1, a[(i - 1) + (i - 1) * lda], b[(i - 1)], ldb, t[(i - 1) * ldt]);
        if (i < m) {
            //
            //           W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)]
            //
            for (j = 1; j <= m - i; j = j + 1) {
                t[(m - 1) + (j - 1) * ldt] = (a[((i + j) - 1) + (i - 1) * lda]);
            }
            Rgemv("N", m - i, p, one, b[((i + 1) - 1)], ldb, b[(i - 1)], ldb, one, t[(m - 1)], ldt);
            //
            //           C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H
            //
            alpha = -(t[(i - 1) * ldt]);
            for (j = 1; j <= m - i; j = j + 1) {
                a[((i + j) - 1) + (i - 1) * lda] += alpha * (t[(m - 1) + (j - 1) * ldt]);
            }
            Rger(m - i, p, alpha, t[(m - 1)], ldt, b[(i - 1)], ldb, b[((i + 1) - 1)], ldb);
        }
    }
    //
    const REAL zero = 0.0f;
    INTEGER np = 0;
    INTEGER mp = 0;
    for (i = 2; i <= m; i = i + 1) {
        //
        //        T(I,1:I-1) := C(I:I-1,1:N) * (alpha * C(I,I:N)^H)
        //
        alpha = -t[(i - 1) * ldt];
        //
        for (j = 1; j <= i - 1; j = j + 1) {
            t[(i - 1) + (j - 1) * ldt] = zero;
        }
        p = min(i - 1, l);
        np = min(n - l + 1, n);
        mp = min(p + 1, m);
        //
        //        Triangular part of B2
        //
        for (j = 1; j <= p; j = j + 1) {
            t[(i - 1) + (j - 1) * ldt] = alpha * b[(i - 1) + ((n - l + j) - 1) * ldb];
        }
        Rtrmv("L", "N", "N", p, b[(np - 1) * ldb], ldb, t[(i - 1)], ldt);
        //
        //        Rectangular part of B2
        //
        Rgemv("N", i - 1 - p, l, alpha, b[(mp - 1) + (np - 1) * ldb], ldb, b[(i - 1) + (np - 1) * ldb], ldb, zero, t[(i - 1) + (mp - 1) * ldt], ldt);
        //
        //        B1
        //
        Rgemv("N", i - 1, n - l, alpha, b, ldb, b[(i - 1)], ldb, one, t[(i - 1)], ldt);
        //
        //        T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1)
        //
        Rtrmv("L", "T", "N", i - 1, t, ldt, t[(i - 1)], ldt);
        //
        //        T(I,I) = tau(I)
        //
        t[(i - 1) + (i - 1) * ldt] = t[(i - 1) * ldt];
        t[(i - 1) * ldt] = zero;
    }
    for (i = 1; i <= m; i = i + 1) {
        for (j = i + 1; j <= m; j = j + 1) {
            t[(i - 1) + (j - 1) * ldt] = t[(j - 1) + (i - 1) * ldt];
            t[(j - 1) + (i - 1) * ldt] = zero;
        }
    }
    //
    //     End of Rtplqt2
    //
}
