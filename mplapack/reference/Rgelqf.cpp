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

void Rgelqf(INTEGER const &m, INTEGER const &n, REAL *a, INTEGER const &lda, REAL *tau, REAL *work, INTEGER const &lwork, INTEGER &info) {
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
    INTEGER nb = iMlaenv[("Rgelqf" - 1) * ldiMlaenv];
    INTEGER lwkopt = m * nb;
    work[1 - 1] = lwkopt;
    bool lquery = (lwork == -1);
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    } else if (lwork < max((INTEGER)1, m) && !lquery) {
        info = -7;
    }
    if (info != 0) {
        Mxerbla("Rgelqf", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    INTEGER k = min(m, n);
    if (k == 0) {
        work[1 - 1] = 1;
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
        nx = max(0, iMlaenv[(3 - 1) + ("Rgelqf" - 1) * ldiMlaenv]);
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
                nbmin = max(2, iMlaenv[(2 - 1) + ("Rgelqf" - 1) * ldiMlaenv]);
            }
        }
    }
    //
    INTEGER i = 0;
    INTEGER ib = 0;
    INTEGER iinfo = 0;
    if (nb >= nbmin && nb < k && nx < k) {
        //
        //        Use blocked code initially
        //
        for (i = 1; i <= k - nx; i = i + nb) {
            ib = min(k - i + 1, nb);
            //
            //           Compute the LQ factorization of the current block
            //           A(i:i+ib-1,i:n)
            //
            Rgelq2(ib, n - i + 1, a[(i - 1) + (i - 1) * lda], lda, tau[i - 1], work, iinfo);
            if (i + ib <= m) {
                //
                //              Form the triangular factor of the block reflector
                //              H = H(i) H(i+1) . . . H(i+ib-1)
                //
                Rlarft("Forward", "Rowwise", n - i + 1, ib, a[(i - 1) + (i - 1) * lda], lda, tau[i - 1], work, ldwork);
                //
                //              Apply H to A(i+ib:m,i:n) from the right
                //
                Rlarfb("Right", "No transpose", "Forward", "Rowwise", m - i - ib + 1, n - i + 1, ib, a[(i - 1) + (i - 1) * lda], lda, work, ldwork, a[((i + ib) - 1) + (i - 1) * lda], lda, work[(ib + 1) - 1], ldwork);
            }
        }
    } else {
        i = 1;
    }
    //
    //     Use unblocked code to factor the last or only block.
    //
    if (i <= k) {
        Rgelq2(m - i + 1, n - i + 1, a[(i - 1) + (i - 1) * lda], lda, tau[i - 1], work, iinfo);
    }
    //
    work[1 - 1] = iws;
    //
    //     End of Rgelqf
    //
}
