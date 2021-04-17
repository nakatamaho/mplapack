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

void Rsytrf_aa(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, INTEGER *ipiv, REAL *work, INTEGER const lwork, INTEGER &info) {
    INTEGER nb = 0;
    bool upper = false;
    bool lquery = false;
    INTEGER lwkopt = 0;
    INTEGER j = 0;
    INTEGER j1 = 0;
    INTEGER jb = 0;
    INTEGER k1 = 0;
    INTEGER j2 = 0;
    REAL alpha = 0.0;
    const REAL one = 1.0;
    INTEGER k2 = 0;
    INTEGER nj = 0;
    INTEGER j3 = 0;
    INTEGER mj = 0;
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
    //     .. Parameters ..
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
    //     Determine the block size
    //
    nb = iMlaenv(1, "Rsytrf_aa", uplo, n, -1, -1, -1);
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
    } else if (lwork < max((INTEGER)1, 2 * n) && !lquery) {
        info = -7;
    }
    //
    if (info == 0) {
        lwkopt = (nb + 1) * n;
        work[1 - 1] = lwkopt;
    }
    //
    if (info != 0) {
        Mxerbla("Rsytrf_aa", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return
    //
    if (n == 0) {
        return;
    }
    ipiv[1 - 1] = 1;
    if (n == 1) {
        return;
    }
    //
    //     Adjust block size based on the workspace size
    //
    if (lwork < ((1 + nb) * n)) {
        nb = (lwork - n) / n;
    }
    //
    if (upper) {
        //
        //        .....................................................
        //        Factorize A as U**T*D*U using the upper triangle of A
        //        .....................................................
        //
        //        Copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N))
        //
        Rcopy(n, &a[(1 - 1)], lda, &work[1 - 1], 1);
        //
        //        J is the main loop index, increasing from 1 to N in steps of
        //        JB, where JB is the number of columns factorized by Rlasyf;
        //        JB is either NB, or N-J+1 for the last block
        //
        j = 0;
    statement_10:
        if (j >= n) {
            goto statement_20;
        }
        //
        //        each step of the main loop
        //         J is the last column of the previous panel
        //         J1 is the first column of the current panel
        //         K1 identifies if the previous column of the panel has been
        //          explicitly stored, e.g., K1=1 for the first panel, and
        //          K1=0 for the rest
        //
        j1 = j + 1;
        jb = min(n - j1 + 1, nb);
        k1 = max((INTEGER)1, j) - j;
        //
        //        Panel factorization
        //
        Rlasyf_aa(uplo, 2 - k1, n - j, jb, &a[(max((INTEGER)1, j) - 1) + ((j + 1) - 1) * lda], lda, &ipiv[(j + 1) - 1], work, n, &work[(n * nb + 1) - 1]);
        //
        //        Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)
        //
        for (j2 = j + 2; j2 <= min(n, j + jb + 1); j2 = j2 + 1) {
            ipiv[j2 - 1] += j;
            if ((j2 != ipiv[j2 - 1]) && ((j1 - k1) > 2)) {
                Rswap(j1 - k1 - 2, &a[(j2 - 1) * lda], 1, &a[(ipiv[j2 - 1] - 1) * lda], 1);
            }
        }
        j += jb;
        //
        //        Trailing submatrix update, where
        //         the row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and
        //         WORK stores the current block of the auxiriarly matrix H
        //
        if (j < n) {
            //
            //           If first panel and JB=1 (NB=1), then nothing to do
            //
            if (j1 > 1 || jb > 1) {
                //
                //              Merge rank-1 update with BLAS-3 update
                //
                alpha = a[(j - 1) + ((j + 1) - 1) * lda];
                a[(j - 1) + ((j + 1) - 1) * lda] = one;
                Rcopy(n - j, &a[((j - 1) - 1) + ((j + 1) - 1) * lda], lda, &work[((j + 1 - j1 + 1) + jb * n) - 1], 1);
                Rscal(n - j, alpha, &work[((j + 1 - j1 + 1) + jb * n) - 1], 1);
                //
                //              K1 identifies if the previous column of the panel has been
                //               explicitly stored, e.g., K1=1 and K2= 0 for the first panel,
                //               while K1=0 and K2=1 for the rest
                //
                if (j1 > 1) {
                    //
                    //                 Not first panel
                    //
                    k2 = 1;
                } else {
                    //
                    //                 First panel
                    //
                    k2 = 0;
                    //
                    //                 First update skips the first column
                    //
                    jb = jb - 1;
                }
                //
                for (j2 = j + 1; j2 <= n; j2 = j2 + nb) {
                    nj = min(nb, n - j2 + 1);
                    //
                    //                 Update (J2, J2) diagonal block with Rgemv
                    //
                    j3 = j2;
                    for (mj = nj - 1; mj >= 1; mj = mj - 1) {
                        Rgemv("No transpose", mj, jb + 1, -one, &work[(j3 - j1 + 1 + k1 * n) - 1], n, &a[((j1 - k2) - 1) + (j3 - 1) * lda], 1, one, &a[(j3 - 1) + (j3 - 1) * lda], lda);
                        j3++;
                    }
                    //
                    //                 Update off-diagonal block of J2-th block row with Rgemm
                    //
                    Rgemm("Transpose", "Transpose", nj, n - j3 + 1, jb + 1, -one, &a[((j1 - k2) - 1) + (j2 - 1) * lda], lda, &work[(j3 - j1 + 1 + k1 * n) - 1], n, one, &a[(j2 - 1) + (j3 - 1) * lda], lda);
                }
                //
                //              Recover T( J, J+1 )
                //
                a[(j - 1) + ((j + 1) - 1) * lda] = alpha;
            }
            //
            //           WORK(J+1, 1) stores H(J+1, 1)
            //
            Rcopy(n - j, &a[((j + 1) - 1) + ((j + 1) - 1) * lda], lda, &work[1 - 1], 1);
        }
        goto statement_10;
    } else {
        //
        //        .....................................................
        //        Factorize A as L*D*L**T using the lower triangle of A
        //        .....................................................
        //
        //        copy first column A(1:N, 1) into H(1:N, 1)
        //         (stored in WORK(1:N))
        //
        Rcopy(n, &a[(1 - 1)], 1, &work[1 - 1], 1);
        //
        //        J is the main loop index, increasing from 1 to N in steps of
        //        JB, where JB is the number of columns factorized by Rlasyf;
        //        JB is either NB, or N-J+1 for the last block
        //
        j = 0;
    statement_11:
        if (j >= n) {
            goto statement_20;
        }
        //
        //        each step of the main loop
        //         J is the last column of the previous panel
        //         J1 is the first column of the current panel
        //         K1 identifies if the previous column of the panel has been
        //          explicitly stored, e.g., K1=1 for the first panel, and
        //          K1=0 for the rest
        //
        j1 = j + 1;
        jb = min(n - j1 + 1, nb);
        k1 = max((INTEGER)1, j) - j;
        //
        //        Panel factorization
        //
        Rlasyf_aa(uplo, 2 - k1, n - j, jb, &a[(max((INTEGER)1, j) - 1) + ((j + 1) - 1) * lda], lda, &ipiv[(j + 1) - 1], work, n, &work[(n * nb + 1) - 1]);
        //
        //        Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)
        //
        for (j2 = j + 2; j2 <= min(n, j + jb + 1); j2 = j2 + 1) {
            ipiv[j2 - 1] += j;
            if ((j2 != ipiv[j2 - 1]) && ((j1 - k1) > 2)) {
                Rswap(j1 - k1 - 2, &a[(j2 - 1)], lda, &a[(ipiv[j2 - 1] - 1)], lda);
            }
        }
        j += jb;
        //
        //        Trailing submatrix update, where
        //          A(J2+1, J1-1) stores L(J2+1, J1) and
        //          WORK(J2+1, 1) stores H(J2+1, 1)
        //
        if (j < n) {
            //
            //           if first panel and JB=1 (NB=1), then nothing to do
            //
            if (j1 > 1 || jb > 1) {
                //
                //              Merge rank-1 update with BLAS-3 update
                //
                alpha = a[((j + 1) - 1) + (j - 1) * lda];
                a[((j + 1) - 1) + (j - 1) * lda] = one;
                Rcopy(n - j, &a[((j + 1) - 1) + ((j - 1) - 1) * lda], 1, &work[((j + 1 - j1 + 1) + jb * n) - 1], 1);
                Rscal(n - j, alpha, &work[((j + 1 - j1 + 1) + jb * n) - 1], 1);
                //
                //              K1 identifies if the previous column of the panel has been
                //               explicitly stored, e.g., K1=1 and K2= 0 for the first panel,
                //               while K1=0 and K2=1 for the rest
                //
                if (j1 > 1) {
                    //
                    //                 Not first panel
                    //
                    k2 = 1;
                } else {
                    //
                    //                 First panel
                    //
                    k2 = 0;
                    //
                    //                 First update skips the first column
                    //
                    jb = jb - 1;
                }
                //
                for (j2 = j + 1; j2 <= n; j2 = j2 + nb) {
                    nj = min(nb, n - j2 + 1);
                    //
                    //                 Update (J2, J2) diagonal block with Rgemv
                    //
                    j3 = j2;
                    for (mj = nj - 1; mj >= 1; mj = mj - 1) {
                        Rgemv("No transpose", mj, jb + 1, -one, &work[(j3 - j1 + 1 + k1 * n) - 1], n, &a[(j3 - 1) + ((j1 - k2) - 1) * lda], lda, one, &a[(j3 - 1) + (j3 - 1) * lda], 1);
                        j3++;
                    }
                    //
                    //                 Update off-diagonal block in J2-th block column with Rgemm
                    //
                    Rgemm("No transpose", "Transpose", n - j3 + 1, nj, jb + 1, -one, &work[(j3 - j1 + 1 + k1 * n) - 1], n, &a[(j2 - 1) + ((j1 - k2) - 1) * lda], lda, one, &a[(j3 - 1) + (j2 - 1) * lda], lda);
                }
                //
                //              Recover T( J+1, J )
                //
                a[((j + 1) - 1) + (j - 1) * lda] = alpha;
            }
            //
            //           WORK(J+1, 1) stores H(J+1, 1)
            //
            Rcopy(n - j, &a[((j + 1) - 1) + ((j + 1) - 1) * lda], 1, &work[1 - 1], 1);
        }
        goto statement_11;
    }
//
statement_20:;
    //
    //     End of Rsytrf_aa
    //
}
