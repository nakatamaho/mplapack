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

void Rsytrf_rk(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, REAL *e, INTEGER *ipiv, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    INTEGER i = 0;
    INTEGER ip = 0;
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
        info = -8;
    }
    //
    if (info == 0) {
        //
        //        Determine the block size
        //
        nb = iMlaenv(1, "Rsytrf_rk", uplo, n, -1, -1, -1);
        lwkopt = n * nb;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Rsytrf_rk", -info);
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
            nbmin = max(2, iMlaenv(2, "Rsytrf_rk", uplo, n, -1, -1, -1));
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
        //        Factorize A as U*D*U**T using the upper triangle of A
        //
        //        K is the main loop index, decreasing from N to 1 in steps of
        //        KB, where KB is the number of columns factorized by Rlasyf_rk;
        //        KB is either NB or NB-1, or K for the last block
        //
        k = n;
    statement_10:
        //
        //        If K < 1, exit from loop
        //
        if (k < 1) {
            goto statement_15;
        }
        //
        if (k > nb) {
            //
            //           Factorize columns k-kb+1:k of A and use blocked code to
            //           update columns 1:k-kb
            //
            Rlasyf_rk(uplo, k, nb, kb, a, lda, e, ipiv, work, ldwork, iinfo);
        } else {
            //
            //           Use unblocked code to factorize columns 1:k of A
            //
            Rsytf2_rk(uplo, k, a, lda, e, ipiv, iinfo);
            kb = k;
        }
        //
        //        Set INFO on the first occurrence of a zero pivot
        //
        if (info == 0 && iinfo > 0) {
            info = iinfo;
        }
        //
        //        No need to adjust IPIV
        //
        //        Apply permutations to the leading panel 1:k-1
        //
        //        Read IPIV from the last block factored, i.e.
        //        indices  k-kb+1:k and apply row permutations to the
        //        last k+1 colunms k+1:N after that block
        //        (We can do the simple loop over IPIV with decrement -1,
        //        since the ABS value of IPIV( I ) represents the row index
        //        of the interchange with row i in both 1x1 and 2x2 pivot cases)
        //
        if (k < n) {
            for (i = k; i >= (k - kb + 1); i = i - 1) {
                ip = abs(ipiv[i - 1]);
                if (ip != i) {
                    Rswap(n - k, &a[(i - 1) + ((k + 1) - 1) * lda], lda, &a[(ip - 1) + ((k + 1) - 1) * lda], lda);
                }
            }
        }
        //
        //        Decrease K and return to the start of the main loop
        //
        k = k - kb;
        goto statement_10;
    //
    //        This label is the exit from main loop over K decreasing
    //        from N to 1 in steps of KB
    //
    statement_15:;
        //
    } else {
        //
        //        Factorize A as L*D*L**T using the lower triangle of A
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        KB, where KB is the number of columns factorized by Rlasyf_rk;
        //        KB is either NB or NB-1, or N-K+1 for the last block
        //
        k = 1;
    statement_20:
        //
        //        If K > N, exit from loop
        //
        if (k > n) {
            goto statement_35;
        }
        //
        if (k <= n - nb) {
            //
            //           Factorize columns k:k+kb-1 of A and use blocked code to
            //           update columns k+kb:n
            //
            Rlasyf_rk(uplo, n - k + 1, nb, kb, &a[(k - 1) + (k - 1) * lda], lda, &e[k - 1], &ipiv[k - 1], work, ldwork, iinfo);
            //
        } else {
            //
            //           Use unblocked code to factorize columns k:n of A
            //
            Rsytf2_rk(uplo, n - k + 1, &a[(k - 1) + (k - 1) * lda], lda, &e[k - 1], &ipiv[k - 1], iinfo);
            kb = n - k + 1;
            //
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
        for (i = k; i <= k + kb - 1; i = i + 1) {
            if (ipiv[i - 1] > 0) {
                ipiv[i - 1] += k - 1;
            } else {
                ipiv[i - 1] = ipiv[i - 1] - k + 1;
            }
        }
        //
        //        Apply permutations to the leading panel 1:k-1
        //
        //        Read IPIV from the last block factored, i.e.
        //        indices  k:k+kb-1 and apply row permutations to the
        //        first k-1 colunms 1:k-1 before that block
        //        (We can do the simple loop over IPIV with increment 1,
        //        since the ABS value of IPIV( I ) represents the row index
        //        of the interchange with row i in both 1x1 and 2x2 pivot cases)
        //
        if (k > 1) {
            for (i = k; i <= (k + kb - 1); i = i + 1) {
                ip = abs(ipiv[i - 1]);
                if (ip != i) {
                    Rswap(k - 1, &a[(i - 1)], lda, &a[(ip - 1)], lda);
                }
            }
        }
        //
        //        Increase K and return to the start of the main loop
        //
        k += kb;
        goto statement_20;
    //
    //        This label is the exit from main loop over K increasing
    //        from 1 to N in steps of KB
    //
    statement_35:;
        //
        //     End Lower
        //
    }
    //
    work[1 - 1] = lwkopt;
    //
    //     End of Rsytrf_rk
    //
}
