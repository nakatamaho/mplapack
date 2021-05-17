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

void Rtzrzf(INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *tau, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    }
    //
    INTEGER lwkopt = 0;
    INTEGER lwkmin = 0;
    INTEGER nb = 0;
    if (info == 0) {
        if (m == 0 || m == n) {
            lwkopt = 1;
            lwkmin = 1;
        } else {
            //
            //           Determine the block size.
            //
            nb = iMlaenv(1, "Rgerqf", " ", m, n, -1, -1);
            lwkopt = m * nb;
            lwkmin = max((INTEGER)1, m);
        }
        work[1 - 1] = lwkopt;
        //
        if (lwork < lwkmin && !lquery) {
            info = -7;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rtzrzf", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    INTEGER i = 0;
    const REAL zero = 0.0;
    if (m == 0) {
        return;
    } else if (m == n) {
        for (i = 1; i <= n; i = i + 1) {
            tau[i - 1] = zero;
        }
        return;
    }
    //
    INTEGER nbmin = 2;
    INTEGER nx = 1;
    INTEGER iws = m;
    INTEGER ldwork = 0;
    if (nb > 1 && nb < m) {
        //
        //        Determine when to cross over from blocked to unblocked code.
        //
        nx = max((INTEGER)0, iMlaenv(3, "Rgerqf", " ", m, n, -1, -1));
        if (nx < m) {
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
                nbmin = max((INTEGER)2, iMlaenv(2, "Rgerqf", " ", m, n, -1, -1));
            }
        }
    }
    //
    INTEGER m1 = 0;
    INTEGER ki = 0;
    INTEGER kk = 0;
    INTEGER ib = 0;
    INTEGER mu = 0;
    if (nb >= nbmin && nb < m && nx < m) {
        //
        //        Use blocked code initially.
        //        The last kk rows are handled by the block method.
        //
        m1 = min(m + 1, n);
        ki = ((m - nx - 1) / nb) * nb;
        kk = min(m, ki + nb);
        //
        for (i = m - kk + ki + 1; i >= m - kk + 1; i = i - nb) {
            ib = min(m - i + 1, nb);
            //
            //           Compute the TZ factorization of the current block
            //           A(i:i+ib-1,i:n)
            //
            Rlatrz(ib, n - i + 1, n - m, &a[(i - 1) + (i - 1) * lda], lda, &tau[i - 1], work);
            if (i > 1) {
                //
                //              Form the triangular factor of the block reflector
                //              H = H(i+ib-1) . . . H(i+1) H(i)
                //
                Rlarzt("Backward", "Rowwise", n - m, ib, &a[(i - 1) + (m1 - 1) * lda], lda, &tau[i - 1], work, ldwork);
                //
                //              Apply H to A(1:i-1,i:n) from the right
                //
                Rlarzb("Right", "No transpose", "Backward", "Rowwise", i - 1, n - i + 1, ib, n - m, &a[(i - 1) + (m1 - 1) * lda], lda, work, ldwork, &a[(i - 1) * lda], lda, &work[(ib + 1) - 1], ldwork);
            }
        }
        mu = i + nb - 1;
    } else {
        mu = m;
    }
    //
    //     Use unblocked code to factor the last or only block
    //
    if (mu > 0) {
        Rlatrz(mu, n, n - m, a, lda, tau, work);
    }
    //
    work[1 - 1] = lwkopt;
    //
    //     End of Rtzrzf
    //
}
