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

void Cungql(INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *tau, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    } else if (n < 0 || n > m) {
        info = -2;
    } else if (k < 0 || k > n) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    }
    //
    INTEGER lwkopt = 0;
    INTEGER nb = 0;
    if (info == 0) {
        if (n == 0) {
            lwkopt = 1;
        } else {
            nb = iMlaenv(1, "Cungql", " ", m, n, k, -1);
            lwkopt = n * nb;
        }
        work[1 - 1] = lwkopt;
        //
        if (lwork < max((INTEGER)1, n) && !lquery) {
            info = -8;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cungql", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    INTEGER nbmin = 2;
    INTEGER nx = 0;
    INTEGER iws = n;
    INTEGER ldwork = 0;
    if (nb > 1 && nb < k) {
        //
        //        Determine when to cross over from blocked to unblocked code.
        //
        nx = max((INTEGER)0, iMlaenv(3, "Cungql", " ", m, n, k, -1));
        if (nx < k) {
            //
            //           Determine if workspace is large enough for blocked code.
            //
            ldwork = n;
            iws = ldwork * nb;
            if (lwork < iws) {
                //
                //              Not enough workspace to use optimal NB:  reduce NB and
                //              determine the minimum value of NB.
                //
                nb = lwork / ldwork;
                nbmin = max((INTEGER)2, iMlaenv(2, "Cungql", " ", m, n, k, -1));
            }
        }
    }
    //
    INTEGER kk = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    const COMPLEX zero = (0.0, 0.0);
    if (nb >= nbmin && nb < k && nx < k) {
        //
        //        Use blocked code after the first block.
        //        The last kk columns are handled by the block method.
        //
        kk = min(k, ((k - nx + nb - 1) / nb) * nb);
        //
        //        Set A(m-kk+1:m,1:n-kk) to zero.
        //
        for (j = 1; j <= n - kk; j = j + 1) {
            for (i = m - kk + 1; i <= m; i = i + 1) {
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
    Cung2l(m - kk, n - kk, k - kk, a, lda, tau, work, iinfo);
    //
    INTEGER ib = 0;
    INTEGER l = 0;
    if (kk > 0) {
        //
        //        Use blocked code
        //
        for (i = k - kk + 1; i <= k; i = i + nb) {
            ib = min(nb, k - i + 1);
            if (n - k + i > 1) {
                //
                //              Form the triangular factor of the block reflector
                //              H = H(i+ib-1) . . . H(i+1) H(i)
                //
                Clarft("Backward", "Columnwise", m - k + i + ib - 1, ib, &a[((n - k + i) - 1) * lda], lda, &tau[i - 1], work, ldwork);
                //
                //              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
                //
                Clarfb("Left", "No transpose", "Backward", "Columnwise", m - k + i + ib - 1, n - k + i - 1, ib, &a[((n - k + i) - 1) * lda], lda, work, ldwork, a, lda, &work[(ib + 1) - 1], ldwork);
            }
            //
            //           Apply H to rows 1:m-k+i+ib-1 of current block
            //
            Cung2l(m - k + i + ib - 1, ib, ib, &a[((n - k + i) - 1) * lda], lda, &tau[i - 1], work, iinfo);
            //
            //           Set rows m-k+i+ib:m of current block to zero
            //
            for (j = n - k + i; j <= n - k + i + ib - 1; j = j + 1) {
                for (l = m - k + i + ib; l <= m; l = l + 1) {
                    a[(l - 1) + (j - 1) * lda] = zero;
                }
            }
        }
    }
    //
    work[1 - 1] = iws;
    //
    //     End of Cungql
    //
}
