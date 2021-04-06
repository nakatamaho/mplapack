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

void Cungtr(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *tau, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    bool lquery = (lwork == -1);
    bool upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    } else if (lwork < max((INTEGER)1, n - 1) && !lquery) {
        info = -7;
    }
    //
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    if (info == 0) {
        if (upper) {
            nb = iMlaenv(1, "Cungql", " ", n - 1, n - 1, n - 1, -1);
        } else {
            nb = iMlaenv(1, "Cungqr", " ", n - 1, n - 1, n - 1, -1);
        }
        lwkopt = max((INTEGER)1, n - 1) * nb;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Cungtr", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        work[1 - 1] = 1;
        return;
    }
    //
    INTEGER j = 0;
    INTEGER i = 0;
    const COMPLEX zero = (0.0, 0.0);
    const COMPLEX one = (1.0, 0.0);
    INTEGER iinfo = 0;
    if (upper) {
        //
        //        Q was determined by a call to Chetrd with UPLO = 'U'
        //
        //        Shift the vectors which define the elementary reflectors one
        //        column to the left, and set the last row and column of Q to
        //        those of the unit matrix
        //
        for (j = 1; j <= n - 1; j = j + 1) {
            for (i = 1; i <= j - 1; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + ((j + 1) - 1) * lda];
            }
            a[(n - 1) + (j - 1) * lda] = zero;
        }
        for (i = 1; i <= n - 1; i = i + 1) {
            a[(i - 1) + (n - 1) * lda] = zero;
        }
        a[(n - 1) + (n - 1) * lda] = one;
        //
        //        Generate Q(1:n-1,1:n-1)
        //
        Cungql(n - 1, n - 1, n - 1, a, lda, tau, work, lwork, iinfo);
        //
    } else {
        //
        //        Q was determined by a call to Chetrd with UPLO = 'L'.
        //
        //        Shift the vectors which define the elementary reflectors one
        //        column to the right, and set the first row and column of Q to
        //        those of the unit matrix
        //
        for (j = n; j >= 2; j = j - 1) {
            a[(j - 1) * lda] = zero;
            for (i = j + 1; i <= n; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = a[(i - 1) + ((j - 1) - 1) * lda];
            }
        }
        a[(1 - 1)] = one;
        for (i = 2; i <= n; i = i + 1) {
            a[(i - 1)] = zero;
        }
        if (n > 1) {
            //
            //           Generate Q(2:n,2:n)
            //
            Cungqr(n - 1, n - 1, n - 1, &a[(2 - 1) + (2 - 1) * lda], lda, tau, work, lwork, iinfo);
        }
    }
    work[1 - 1] = lwkopt;
    //
    //     End of Cungtr
    //
}
