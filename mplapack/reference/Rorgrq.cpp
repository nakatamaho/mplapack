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

void Rorgrq(INTEGER const m, INTEGER const n, INTEGER const k, REAL *a, INTEGER const lda, REAL *tau, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    bool lquery = (lwork == -1);
    if (m < 0) {
        info = -1;
    } else if (n < m) {
        info = -2;
    } else if (k < 0 || k > m) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    }
    //
    INTEGER lwkopt = 0;
    INTEGER nb = 0;
    if (info == 0) {
        if (m <= 0) {
            lwkopt = 1;
        } else {
            nb = iMlaenv(1, "Rorgrq", " ", m, n, k, -1);
            lwkopt = m * nb;
        }
        work[1 - 1] = lwkopt;
        //
        if (lwork < max((INTEGER)1, m) && !lquery) {
            info = -8;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rorgrq", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (m <= 0) {
        return;
    }
    //
    INTEGER nbmin = 2;
    INTEGER nx = 0;
    INTEGER iws = m;
    INTEGER ldwork = 0;
    if (nb > 1 && nb < k) {
        //
        //        Determine when to cross over from blocked to unblocked code.
        //
        nx = max(0, iMlaenv(3, "Rorgrq", " ", m, n, k, -1));
        if (nx < k) {
            //
            //           Determine if workspace is large enough for blocked code.
            //
            ldwork = m;
            iws = ldwork * nb;
            if (lwork < iws) {
                //
                //              Not enough workspace to use optimal NB:  reduce NB and
                //              determine the minimum value of NB.
                //
                nb = lwork / ldwork;
                nbmin = max(2, iMlaenv(2, "Rorgrq", " ", m, n, k, -1));
            }
        }
    }
    //
    INTEGER kk = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    if (nb >= nbmin && nb < k && nx < k) {
        //
        //        Use blocked code after the first block.
        //        The last kk rows are handled by the block method.
        //
        kk = min(k, ((k - nx + nb - 1) / nb) * nb);
        //
        //        Set A(1:m-kk,n-kk+1:n) to zero.
        //
        for (j = n - kk + 1; j <= n; j = j + 1) {
            for (i = 1; i <= m - kk; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = zero;
            }
        }
    } else {
        kk = 0;
    }
    //
    //     Use unblocked code for the first or only block.
    //
    INTEGER iinfo = 0;
    Rorgr2(m - kk, n - kk, k - kk, a, lda, tau, work, iinfo);
    //
    INTEGER ib = 0;
    INTEGER ii = 0;
    INTEGER l = 0;
    if (kk > 0) {
        //
        //        Use blocked code
        //
        for (i = k - kk + 1; i <= k; i = i + nb) {
            ib = min(nb, k - i + 1);
            ii = m - k + i;
            if (ii > 1) {
                //
                //              Form the triangular factor of the block reflector
                //              H = H(i+ib-1) . . . H(i+1) H(i)
                //
                Rlarft("Backward", "Rowwise", n - k + i + ib - 1, ib, &a[(ii - 1)], lda, &tau[i - 1], work, ldwork);
                //
                //              Apply H**T to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
                //
                Rlarfb("Right", "Transpose", "Backward", "Rowwise", ii - 1, n - k + i + ib - 1, ib, &a[(ii - 1)], lda, work, ldwork, a, lda, &work[(ib + 1) - 1], ldwork);
            }
            //
            //           Apply H**T to columns 1:n-k+i+ib-1 of current block
            //
            Rorgr2(ib, n - k + i + ib - 1, ib, &a[(ii - 1)], lda, &tau[i - 1], work, iinfo);
            //
            //           Set columns n-k+i+ib:n of current block to zero
            //
            for (l = n - k + i + ib; l <= n; l = l + 1) {
                for (j = ii; j <= ii + ib - 1; j = j + 1) {
                    a[(j - 1) + (l - 1) * lda] = zero;
                }
            }
        }
    }
    //
    work[1 - 1] = iws;
    //
    //     End of Rorgrq
    //
}
