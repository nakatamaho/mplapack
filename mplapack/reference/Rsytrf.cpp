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

void Rsytrf(const char *uplo, INTEGER const &n, REAL *a, INTEGER const &lda, arr_ref<INTEGER> ipiv, REAL *work, INTEGER const &lwork, INTEGER &info) {
    bool upper = false;
    bool lquery = false;
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    INTEGER nbmin = 0;
    INTEGER ldwork = 0;
    INTEGER iws = 0;
    INTEGER k = 0;
    INTEGER kb = 0;
    INTEGER iinfo = 0;
    INTEGER j = 0;
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
    upper = Mlsame(uplo, "U");
    lquery = (lwork == -1);
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    } else if (lwork < 1 && !lquery) {
        info = -7;
    }
    //
    if (info == 0) {
        //
        //        Determine the block size
        //
        nb = iMlaenv[("Rsytrf" - 1) * ldiMlaenv];
        lwkopt = n * nb;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Rsytrf", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    nbmin = 2;
    ldwork = n;
    if (nb > 1 && nb < n) {
        iws = ldwork * nb;
        if (lwork < iws) {
            nb = max(lwork / ldwork, 1);
            nbmin = max(2, iMlaenv[(2 - 1) + ("Rsytrf" - 1) * ldiMlaenv]);
        }
    } else {
        iws = 1;
    }
    if (nb < nbmin) {
        nb = n;
    }
    //
    if (upper) {
        //
        //        Factorize A as U**T*D*U using the upper triangle of A
        //
        //        K is the main loop index, decreasing from N to 1 in steps of
        //        KB, where KB is the number of columns factorized by Rlasyf;
        //        KB is either NB or NB-1, or K for the last block
        //
        k = n;
    statement_10:
        //
        //        If K < 1, exit from loop
        //
        if (k < 1) {
            goto statement_40;
        }
        //
        if (k > nb) {
            //
            //           Factorize columns k-kb+1:k of A and use blocked code to
            //           update columns 1:k-kb
            //
            Rlasyf(uplo, k, nb, kb, a, lda, ipiv, work, ldwork, iinfo);
        } else {
            //
            //           Use unblocked code to factorize columns 1:k of A
            //
            Rsytf2(uplo, k, a, lda, ipiv, iinfo);
            kb = k;
        }
        //
        //        Set INFO on the first occurrence of a zero pivot
        //
        if (info == 0 && iinfo > 0) {
            info = iinfo;
        }
        //
        //        Decrease K and return to the start of the main loop
        //
        k = k - kb;
        goto statement_10;
        //
    } else {
        //
        //        Factorize A as L*D*L**T using the lower triangle of A
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        KB, where KB is the number of columns factorized by Rlasyf;
        //        KB is either NB or NB-1, or N-K+1 for the last block
        //
        k = 1;
    statement_20:
        //
        //        If K > N, exit from loop
        //
        if (k > n) {
            goto statement_40;
        }
        //
        if (k <= n - nb) {
            //
            //           Factorize columns k:k+kb-1 of A and use blocked code to
            //           update columns k+kb:n
            //
            Rlasyf(uplo, n - k + 1, nb, kb, a[(k - 1) + (k - 1) * lda], lda, ipiv[k - 1], work, ldwork, iinfo);
        } else {
            //
            //           Use unblocked code to factorize columns k:n of A
            //
            Rsytf2(uplo, n - k + 1, a[(k - 1) + (k - 1) * lda], lda, ipiv[k - 1], iinfo);
            kb = n - k + 1;
        }
        //
        //        Set INFO on the first occurrence of a zero pivot
        //
        if (info == 0 && iinfo > 0) {
            info = iinfo + k - 1;
        }
        //
        //        Adjust IPIV
        //
        for (j = k; j <= k + kb - 1; j = j + 1) {
            if (ipiv[j - 1] > 0) {
                ipiv[j - 1] += k - 1;
            } else {
                ipiv[j - 1] = ipiv[j - 1] - k + 1;
            }
        }
        //
        //        Increase K and return to the start of the main loop
        //
        k += kb;
        goto statement_20;
        //
    }
//
statement_40:
    work[1 - 1] = lwkopt;
    //
    //     End of Rsytrf
    //
}
