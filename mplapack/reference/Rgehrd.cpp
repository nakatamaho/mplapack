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

void Rgehrd(INTEGER const n, INTEGER const ilo, INTEGER const ihi, REAL *a, INTEGER const lda, REAL *tau, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //     Test the input parameters
    //
    info = 0;
    bool lquery = (lwork == -1);
    if (n < 0) {
        info = -1;
    } else if (ilo < 1 || ilo > max((INTEGER)1, n)) {
        info = -2;
    } else if (ihi < min(ilo, n) || ihi > n) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (lwork < max((INTEGER)1, n) && !lquery) {
        info = -8;
    }
    //
    const INTEGER nbmax = 64;
    INTEGER nb = 0;
    const INTEGER ldt = nbmax + 1;
    const INTEGER tsize = ldt * nbmax;
    INTEGER lwkopt = 0;
    if (info == 0) {
        //
        //        Compute the workspace requirements
        //
        nb = min(nbmax, iMlaenv(1, "Rgehrd", " ", n, ilo, ihi, -1));
        lwkopt = n * nb + tsize;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Rgehrd", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
    //
    INTEGER i = 0;
    const REAL zero = 0.0;
    for (i = 1; i <= ilo - 1; i = i + 1) {
        tau[i - 1] = zero;
    }
    for (i = max((INTEGER)1, ihi); i <= n - 1; i = i + 1) {
        tau[i - 1] = zero;
    }
    //
    //     Quick return if possible
    //
    INTEGER nh = ihi - ilo + 1;
    if (nh <= 1) {
        work[1 - 1] = 1;
        return;
    }
    //
    //     Determine the block size
    //
    nb = min(nbmax, iMlaenv(1, "Rgehrd", " ", n, ilo, ihi, -1));
    INTEGER nbmin = 2;
    INTEGER nx = 0;
    if (nb > 1 && nb < nh) {
        //
        //        Determine when to cross over from blocked to unblocked code
        //        (last block is always handled by unblocked code)
        //
        nx = max(nb, iMlaenv(3, "Rgehrd", " ", n, ilo, ihi, -1));
        if (nx < nh) {
            //
            //           Determine if workspace is large enough for blocked code
            //
            if (lwork < n * nb + tsize) {
                //
                //              Not enough workspace to use optimal NB:  determine the
                //              minimum value of NB, and reduce NB or force use of
                //              unblocked code
                //
                nbmin = max(2, iMlaenv(2, "Rgehrd", " ", n, ilo, ihi, -1));
                if (lwork >= (n * nbmin + tsize)) {
                    nb = (lwork - tsize) / n;
                } else {
                    nb = 1;
                }
            }
        }
    }
    INTEGER ldwork = n;
    //
    INTEGER iwt = 0;
    INTEGER ib = 0;
    REAL ei = 0.0;
    const REAL one = 1.0;
    INTEGER j = 0;
    if (nb < nbmin || nb >= nh) {
        //
        //        Use unblocked code below
        //
        i = ilo;
        //
    } else {
        //
        //        Use blocked code
        //
        iwt = 1 + n * nb;
        for (i = ilo; i <= ihi - 1 - nx; i = i + nb) {
            ib = min(nb, ihi - i);
            //
            //           Reduce columns i:i+ib-1 to Hessenberg form, returning the
            //           matrices V and T of the block reflector H = I - V*T*V**T
            //           which performs the reduction, and also the matrix Y = A*V*T
            //
            Rlahr2(ihi, i, ib, &a[(i - 1) * lda], lda, &tau[i - 1], &work[iwt - 1], ldt, work, ldwork);
            //
            //           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
            //           right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set
            //           to 1
            //
            ei = a[((i + ib) - 1) + ((i + ib - 1) - 1) * lda];
            a[((i + ib) - 1) + ((i + ib - 1) - 1) * lda] = one;
            Rgemm("No transpose", "Transpose", ihi, ihi - i - ib + 1, ib, -one, work, ldwork, &a[((i + ib) - 1) + (i - 1) * lda], lda, one, &a[((i + ib) - 1) * lda], lda);
            a[((i + ib) - 1) + ((i + ib - 1) - 1) * lda] = ei;
            //
            //           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
            //           right
            //
            Rtrmm("Right", "Lower", "Transpose", "Unit", i, ib - 1, one, &a[((i + 1) - 1) + (i - 1) * lda], lda, work, ldwork);
            for (j = 0; j <= ib - 2; j = j + 1) {
                Raxpy(i, -one, &work[(ldwork * j + 1) - 1], 1, &a[((i + j + 1) - 1) * lda], 1);
            }
            //
            //           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
            //           left
            //
            Rlarfb("Left", "Transpose", "Forward", "Columnwise", ihi - i, n - i - ib + 1, ib, &a[((i + 1) - 1) + (i - 1) * lda], lda, &work[iwt - 1], ldt, &a[((i + 1) - 1) + ((i + ib) - 1) * lda], lda, work, ldwork);
        }
    }
    //
    //     Use unblocked code to reduce the rest of the matrix
    //
    INTEGER iinfo = 0;
    Rgehd2(n, i, ihi, a, lda, tau, work, iinfo);
    work[1 - 1] = lwkopt;
    //
    //     End of Rgehrd
    //
}
