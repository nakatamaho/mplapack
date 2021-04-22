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

void Chetrd(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL *d, REAL *e, COMPLEX *tau, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    bool upper = Mlsame(uplo, "U");
    bool lquery = (lwork == -1);
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    } else if (lwork < 1 && !lquery) {
        info = -9;
    }
    //
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    if (info == 0) {
        //
        //        Determine the block size.
        //
        nb = iMlaenv(1, "Chetrd", uplo, n, -1, -1, -1);
        lwkopt = n * nb;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Chetrd", -info);
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
    INTEGER nx = n;
    INTEGER iws = 1;
    INTEGER ldwork = 0;
    INTEGER nbmin = 0;
    if (nb > 1 && nb < n) {
        //
        //        Determine when to cross over from blocked to unblocked code
        //        (last block is always handled by unblocked code).
        //
        nx = max(nb, iMlaenv(3, "Chetrd", uplo, n, -1, -1, -1));
        if (nx < n) {
            //
            //           Determine if workspace is large enough for blocked code.
            //
            ldwork = n;
            iws = ldwork * nb;
            if (lwork < iws) {
                //
                //              Not enough workspace to use optimal NB:  determine the
                //              minimum value of NB, and reduce NB or force use of
                //              unblocked code by setting NX = N.
                //
                nb = max(lwork / ldwork, (INTEGER)1);
                nbmin = iMlaenv(2, "Chetrd", uplo, n, -1, -1, -1);
                if (nb < nbmin) {
                    nx = n;
                }
            }
        } else {
            nx = n;
        }
    } else {
        nb = 1;
    }
    //
    INTEGER kk = 0;
    INTEGER i = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const REAL one = 1.0;
    INTEGER j = 0;
    INTEGER iinfo = 0;
    if (upper) {
        //
        //        Reduce the upper triangle of A.
        //        Columns 1:kk are handled by the unblocked method.
        //
        kk = n - ((n - nx + nb - 1) / nb) * nb;
        for (i = n - nb + 1; i >= kk + 1; i = i - nb) {
            //
            //           Reduce columns i:i+nb-1 to tridiagonal form and form the
            //           matrix W which is needed to update the unreduced part of
            //           the matrix
            //
            Clatrd(uplo, i + nb - 1, nb, a, lda, e, tau, work, ldwork);
            //
            //           Update the unreduced submatrix A(1:i-1,1:i-1), using an
            //           update of the form:  A := A - V*W**H - W*V**H
            //
            Cher2k(uplo, "No transpose", i - 1, nb, -cone, &a[(i - 1) * lda], lda, work, ldwork, one, a, lda);
            //
            //           Copy superdiagonal elements back into A, and diagonal
            //           elements into D
            //
            for (j = i; j <= i + nb - 1; j = j + 1) {
                a[((j - 1) - 1) + (j - 1) * lda] = e[(j - 1) - 1];
                d[j - 1] = a[(j - 1) + (j - 1) * lda].real();
            }
        }
        //
        //        Use unblocked code to reduce the last or only block
        //
        Chetd2(uplo, kk, a, lda, d, e, tau, iinfo);
    } else {
        //
        //        Reduce the lower triangle of A
        //
        for (i = 1; i <= n - nx; i = i + nb) {
            //
            //           Reduce columns i:i+nb-1 to tridiagonal form and form the
            //           matrix W which is needed to update the unreduced part of
            //           the matrix
            //
            Clatrd(uplo, n - i + 1, nb, &a[(i - 1) + (i - 1) * lda], lda, &e[i - 1], &tau[i - 1], work, ldwork);
            //
            //           Update the unreduced submatrix A(i+nb:n,i+nb:n), using
            //           an update of the form:  A := A - V*W**H - W*V**H
            //
            Cher2k(uplo, "No transpose", n - i - nb + 1, nb, -cone, &a[((i + nb) - 1) + (i - 1) * lda], lda, &work[(nb + 1) - 1], ldwork, one, &a[((i + nb) - 1) + ((i + nb) - 1) * lda], lda);
            //
            //           Copy subdiagonal elements back into A, and diagonal
            //           elements into D
            //
            for (j = i; j <= i + nb - 1; j = j + 1) {
                a[((j + 1) - 1) + (j - 1) * lda] = e[j - 1];
                d[j - 1] = a[(j - 1) + (j - 1) * lda].real();
            }
        }
        //
        //        Use unblocked code to reduce the last or only block
        //
        Chetd2(uplo, n - i + 1, &a[(i - 1) + (i - 1) * lda], lda, &d[i - 1], &e[i - 1], &tau[i - 1], iinfo);
    }
    //
    work[1 - 1] = lwkopt;
    //
    //     End of Chetrd
    //
}
