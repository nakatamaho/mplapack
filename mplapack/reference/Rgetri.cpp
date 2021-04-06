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

void Rgetri(INTEGER const n, REAL *a, INTEGER const lda, INTEGER *ipiv, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //     Test the input parameters.
    //
    info = 0;
    INTEGER nb = iMlaenv(1, "Rgetri", " ", n, -1, -1, -1);
    INTEGER lwkopt = n * nb;
    work[1 - 1] = lwkopt;
    bool lquery = (lwork == -1);
    if (n < 0) {
        info = -1;
    } else if (lda < max((INTEGER)1, n)) {
        info = -3;
    } else if (lwork < max((INTEGER)1, n) && !lquery) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Rgetri", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Form inv(U).  If INFO > 0 from Rtrtri, then U is singular,
    //     and the inverse is not computed.
    //
    Rtrtri("Upper", "Non-unit", n, a, lda, info);
    if (info > 0) {
        return;
    }
    //
    INTEGER nbmin = 2;
    INTEGER ldwork = n;
    INTEGER iws = 0;
    if (nb > 1 && nb < n) {
        iws = max(ldwork * nb, (INTEGER)1);
        if (lwork < iws) {
            nb = lwork / ldwork;
            nbmin = max((INTEGER)2, iMlaenv(2, "Rgetri", " ", n, -1, -1, -1));
        }
    } else {
        iws = n;
    }
    //
    //     Solve the equation inv(A)*L = inv(U) for inv(A).
    //
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    INTEGER nn = 0;
    INTEGER jb = 0;
    INTEGER jj = 0;
    if (nb < nbmin || nb >= n) {
        //
        //        Use unblocked code.
        //
        for (j = n; j >= 1; j = j - 1) {
            //
            //           Copy current column of L to WORK and replace with zeros.
            //
            for (i = j + 1; i <= n; i = i + 1) {
                work[i - 1] = a[(i - 1) + (j - 1) * lda];
                a[(i - 1) + (j - 1) * lda] = zero;
            }
            //
            //           Compute current column of inv(A).
            //
            if (j < n) {
                Rgemv("No transpose", n, n - j, -one, &a[((j + 1) - 1) * lda], lda, &work[(j + 1) - 1], 1, one, &a[(j - 1) * lda], 1);
            }
        }
    } else {
        //
        //        Use blocked code.
        //
        nn = ((n - 1) / nb) * nb + 1;
        for (j = nn; j <= 1; j = j + -nb) {
            jb = min(nb, n - j + 1);
            //
            //           Copy current block column of L to WORK and replace with
            //           zeros.
            //
            for (jj = j; jj <= j + jb - 1; jj = jj + 1) {
                for (i = jj + 1; i <= n; i = i + 1) {
                    work[(i + (jj - j) * ldwork) - 1] = a[(i - 1) + (jj - 1) * lda];
                    a[(i - 1) + (jj - 1) * lda] = zero;
                }
            }
            //
            //           Compute current block column of inv(A).
            //
            if (j + jb <= n) {
                Rgemm("No transpose", "No transpose", n, jb, n - j - jb + 1, -one, &a[((j + jb) - 1) * lda], lda, &work[(j + jb) - 1], ldwork, one, &a[(j - 1) * lda], lda);
            }
            Rtrsm("Right", "Lower", "No transpose", "Unit", n, jb, one, &work[j - 1], ldwork, &a[(j - 1) * lda], lda);
        }
    }
    //
    //     Apply column INTEGERerchanges.
    //
    INTEGER jp = 0;
    for (j = n - 1; j >= 1; j = j - 1) {
        jp = ipiv[j - 1];
        if (jp != j) {
            Rswap(n, &a[(j - 1) * lda], 1, &a[(jp - 1) * lda], 1);
        }
    }
    //
    work[1 - 1] = iws;
    //
    //     End of Rgetri
    //
}
