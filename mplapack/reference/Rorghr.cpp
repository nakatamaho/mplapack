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

void Rorghr(INTEGER const &n, INTEGER const &ilo, INTEGER const &ihi, REAL *a, INTEGER const &lda, REAL *tau, REAL *work, INTEGER const &lwork, INTEGER &info) {
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
    //     Test the input arguments
    //
    info = 0;
    INTEGER nh = ihi - ilo;
    bool lquery = (lwork == -1);
    if (n < 0) {
        info = -1;
    } else if (ilo < 1 || ilo > max((INTEGER)1, n)) {
        info = -2;
    } else if (ihi < min(ilo, n) || ihi > n) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (lwork < max((INTEGER)1, nh) && !lquery) {
        info = -8;
    }
    //
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    if (info == 0) {
        nb = iMlaenv[("Rorgqr" - 1) * ldiMlaenv];
        lwkopt = max((INTEGER)1, nh) * nb;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Rorghr", -info);
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
    //     Shift the vectors which define the elementary reflectors one
    //     column to the right, and set the first ilo and the last n-ihi
    //     rows and columns to those of the unit matrix
    //
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    for (j = ihi; j >= ilo + 1; j = j - 1) {
        for (i = 1; i <= j - 1; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
        }
        for (i = j + 1; i <= ihi; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = a[(i - 1) + ((j - 1) - 1) * lda];
        }
        for (i = ihi + 1; i <= n; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
        }
    }
    const REAL one = 1.0;
    for (j = 1; j <= ilo; j = j + 1) {
        for (i = 1; i <= n; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
        }
        a[(j - 1) + (j - 1) * lda] = one;
    }
    for (j = ihi + 1; j <= n; j = j + 1) {
        for (i = 1; i <= n; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
        }
        a[(j - 1) + (j - 1) * lda] = one;
    }
    //
    INTEGER iinfo = 0;
    if (nh > 0) {
        //
        //        Generate Q(ilo+1:ihi,ilo+1:ihi)
        //
        Rorgqr(nh, nh, nh, a[((ilo + 1) - 1) + ((ilo + 1) - 1) * lda], lda, tau[ilo - 1], work, lwork, iinfo);
    }
    work[1 - 1] = lwkopt;
    //
    //     End of Rorghr
    //
}
