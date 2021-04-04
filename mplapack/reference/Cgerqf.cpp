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

void Cgerqf(INTEGER const &m, INTEGER const &n, COMPLEX *a, INTEGER const &lda, COMPLEX *tau, COMPLEX *work, INTEGER const &lwork, INTEGER &info) {
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
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    }
    //
    INTEGER k = 0;
    INTEGER lwkopt = 0;
    INTEGER nb = 0;
    if (info == 0) {
        k = min(m, n);
        if (k == 0) {
            lwkopt = 1;
        } else {
            nb = iMlaenv[("Cgerqf" - 1) * ldiMlaenv];
            lwkopt = m * nb;
        }
        work[1 - 1] = lwkopt;
        //
        if (lwork < max((INTEGER)1, m) && !lquery) {
            info = -7;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cgerqf", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (k == 0) {
        return;
    }
    //
    INTEGER nbmin = 2;
    INTEGER nx = 1;
    INTEGER iws = m;
    INTEGER ldwork = 0;
    if (nb > 1 && nb < k) {
        //
        //        Determine when to cross over from blocked to unblocked code.
        //
        nx = max(0, iMlaenv[(3 - 1) + ("Cgerqf" - 1) * ldiMlaenv]);
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
                nbmin = max(2, iMlaenv[(2 - 1) + ("Cgerqf" - 1) * ldiMlaenv]);
            }
        }
    }
    //
    INTEGER ki = 0;
    INTEGER kk = 0;
    INTEGER i = 0;
    INTEGER ib = 0;
    INTEGER iinfo = 0;
    INTEGER mu = 0;
    INTEGER nu = 0;
    if (nb >= nbmin && nb < k && nx < k) {
        //
        //        Use blocked code initially.
        //        The last kk rows are handled by the block method.
        //
        ki = ((k - nx - 1) / nb) * nb;
        kk = min(k, ki + nb);
        //
        for (i = k - kk + ki + 1; i <= k - kk + 1; i = i + -nb) {
            ib = min(k - i + 1, nb);
            //
            //           Compute the RQ factorization of the current block
            //           A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1)
            //
            Cgerq2(ib, n - k + i + ib - 1, a[((m - k + i) - 1)], lda, tau[i - 1], work, iinfo);
            if (m - k + i > 1) {
                //
                //              Form the triangular factor of the block reflector
                //              H = H(i+ib-1) . . . H(i+1) H(i)
                //
                Clarft("Backward", "Rowwise", n - k + i + ib - 1, ib, a[((m - k + i) - 1)], lda, tau[i - 1], work, ldwork);
                //
                //              Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the right
                //
                Clarfb("Right", "No transpose", "Backward", "Rowwise", m - k + i - 1, n - k + i + ib - 1, ib, a[((m - k + i) - 1)], lda, work, ldwork, a, lda, work[(ib + 1) - 1], ldwork);
            }
        }
        mu = m - k + i + nb - 1;
        nu = n - k + i + nb - 1;
    } else {
        mu = m;
        nu = n;
    }
    //
    //     Use unblocked code to factor the last or only block
    //
    if (mu > 0 && nu > 0) {
        Cgerq2(mu, nu, a, lda, tau, work, iinfo);
    }
    //
    work[1 - 1] = iws;
    //
    //     End of Cgerqf
    //
}
