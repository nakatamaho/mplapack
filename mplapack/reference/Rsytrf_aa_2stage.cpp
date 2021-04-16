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

void Rsytrf_aa_2stage(const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, REAL *tb, INTEGER const ltb, INTEGER *ipiv, INTEGER *ipiv2, REAL *work, INTEGER const lwork, INTEGER &info) {
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
    //     Test the input parameters.
    //
    info = 0;
    bool upper = Mlsame(uplo, "U");
    bool wquery = (lwork == -1);
    bool tquery = (ltb == -1);
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    } else if (ltb < 4 * n && !tquery) {
        info = -6;
    } else if (lwork < n && !wquery) {
        info = -10;
    }
    //
    if (info != 0) {
        Mxerbla("Rsytrf_aa_2stage", -info);
        return;
    }
    //
    //     Answer the query
    //
    INTEGER nb = iMlaenv(1, "Rsytrf_aa_2stage", uplo, n, -1, -1, -1);
    if (info == 0) {
        if (tquery) {
            tb[1 - 1] = (3 * nb + 1) * n;
        }
        if (wquery) {
            work[1 - 1] = n * nb;
        }
    }
    if (tquery || wquery) {
        return;
    }
    //
    //     Quick return
    //
    if (n == 0) {
        return;
    }
    //
    //     Determine the number of the block size
    //
    INTEGER ldtb = ltb / n;
    if (ldtb < 3 * nb + 1) {
        nb = (ldtb - 1) / 3;
    }
    if (lwork < nb * n) {
        nb = lwork / n;
    }
    //
    //     Determine the number of the block columns
    //
    INTEGER nt = (n + nb - 1) / nb;
    INTEGER td = 2 * nb;
    INTEGER kb = min(nb, n);
    //
    //     Initialize vectors/matrices
    //
    INTEGER j = 0;
    for (j = 1; j <= kb; j = j + 1) {
        ipiv[j - 1] = j;
    }
    //
    //     Save NB
    //
    tb[1 - 1] = nb;
    //
    INTEGER i = 0;
    INTEGER jb = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    INTEGER iinfo = 0;
    INTEGER k = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    REAL piv = 0.0;
    if (upper) {
        //
        //        .....................................................
        //        Factorize A as U**T*D*U using the upper triangle of A
        //        .....................................................
        //
        for (j = 0; j <= nt - 1; j = j + 1) {
            //
            //           Generate Jth column of W and H
            //
            kb = min(nb, n - j * nb);
            for (i = 1; i <= j - 1; i = i + 1) {
                if (i == 1) {
                    //                 H(I,J) = T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J)
                    if (i == (j - 1)) {
                        jb = nb + kb;
                    } else {
                        jb = 2 * nb;
                    }
                    Rgemm("NoTranspose", "NoTranspose", nb, kb, jb, one, &tb[(td + 1 + (i * nb) * ldtb) - 1], ldtb - 1, &a[(((i - 1) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, zero, &work[(i * nb + 1) - 1], n);
                } else {
                    //                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J)
                    if (i == j - 1) {
                        jb = 2 * nb + kb;
                    } else {
                        jb = 3 * nb;
                    }
                    Rgemm("NoTranspose", "NoTranspose", nb, kb, jb, one, &tb[(td + nb + 1 + ((i - 1) * nb) * ldtb) - 1], ldtb - 1, &a[(((i - 2) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, zero, &work[(i * nb + 1) - 1], n);
                }
            }
            //
            //           Compute T(J,J)
            //
            Rlacpy("Upper", kb, kb, &a[((j * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, &tb[(td + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
            if (j > 1) {
                //              T(J,J) = U(1:J,J)'*H(1:J)
                Rgemm("Transpose", "NoTranspose", kb, kb, (j - 1) * nb, -one, &a[((j * nb + 1) - 1) * lda], lda, &work[(nb + 1) - 1], n, one, &tb[(td + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
                //              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J)
                Rgemm("Transpose", "NoTranspose", kb, nb, kb, one, &a[(((j - 1) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, &tb[(td + nb + 1 + ((j - 1) * nb) * ldtb) - 1], ldtb - 1, zero, &work[1 - 1], n);
                Rgemm("NoTranspose", "NoTranspose", kb, kb, nb, -one, &work[1 - 1], n, &a[(((j - 2) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, one, &tb[(td + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
            }
            if (j > 0) {
                Rsygst(1, "Upper", kb, &tb[(td + 1 + (j * nb) * ldtb) - 1], ldtb - 1, &a[(((j - 1) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, iinfo);
            }
            //
            //           Expand T(J,J) into full format
            //
            for (i = 1; i <= kb; i = i + 1) {
                for (k = i + 1; k <= kb; k = k + 1) {
                    tb[(td + (k - i) + 1 + (j * nb + i - 1) * ldtb) - 1] = tb[(td - (k - (i + 1)) + (j * nb + k - 1) * ldtb) - 1];
                }
            }
            //
            if (j < nt - 1) {
                if (j > 0) {
                    //
                    //                 Compute H(J,J)
                    //
                    if (j == 1) {
                        Rgemm("NoTranspose", "NoTranspose", kb, kb, kb, one, &tb[(td + 1 + (j * nb) * ldtb) - 1], ldtb - 1, &a[(((j - 1) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, zero, &work[(j * nb + 1) - 1], n);
                    } else {
                        Rgemm("NoTranspose", "NoTranspose", kb, kb, nb + kb, one, &tb[(td + nb + 1 + ((j - 1) * nb) * ldtb) - 1], ldtb - 1, &a[(((j - 2) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, zero, &work[(j * nb + 1) - 1], n);
                    }
                    //
                    //                 Update with the previous column
                    //
                    Rgemm("Transpose", "NoTranspose", nb, n - (j + 1) * nb, j * nb, -one, &work[(nb + 1) - 1], n, &a[(((j + 1) * nb + 1) - 1) * lda], lda, one, &a[((j * nb + 1) - 1) + (((j + 1) * nb + 1) - 1) * lda], lda);
                }
                //
                //              Copy panel to workspace to call Rgetrf
                //
                for (k = 1; k <= nb; k = k + 1) {
                    Rcopy(n - (j + 1) * nb, &a[((j * nb + k) - 1) + (((j + 1) * nb + 1) - 1) * lda], lda, &work[(1 + (k - 1) * n) - 1], 1);
                }
                //
                //              Factorize panel
                //
                Rgetrf(n - (j + 1) * nb, nb, work, n, &ipiv[((j + 1) * nb + 1) - 1], iinfo);
                //               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN
                //                  INFO = IINFO+(J+1)*NB
                //               END IF
                //
                //              Copy panel back
                //
                for (k = 1; k <= nb; k = k + 1) {
                    Rcopy(n - (j + 1) * nb, &work[(1 + (k - 1) * n) - 1], 1, &a[((j * nb + k) - 1) + (((j + 1) * nb + 1) - 1) * lda], lda);
                }
                //
                //              Compute T(J+1, J), zero out for GEMM update
                //
                kb = min(nb, n - (j + 1) * nb);
                Rlaset("Full", kb, nb, zero, zero, &tb[(td + nb + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
                Rlacpy("Upper", kb, nb, work, n, &tb[(td + nb + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
                if (j > 0) {
                    Rtrsm("R", "U", "N", "U", kb, nb, one, &a[(((j - 1) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, &tb[(td + nb + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
                }
                //
                //              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM
                //              updates
                //
                for (k = 1; k <= nb; k = k + 1) {
                    for (i = 1; i <= kb; i = i + 1) {
                        tb[(td - nb + k - i + 1 + (j * nb + nb + i - 1) * ldtb) - 1] = tb[(td + nb + i - k + 1 + (j * nb + k - 1) * ldtb) - 1];
                    }
                }
                Rlaset("Lower", kb, nb, zero, one, &a[((j * nb + 1) - 1) + (((j + 1) * nb + 1) - 1) * lda], lda);
                //
                //              Apply pivots to trailing submatrix of A
                //
                for (k = 1; k <= kb; k = k + 1) {
                    //                 > Adjust ipiv
                    ipiv[((j + 1) * nb + k) - 1] += (j + 1) * nb;
                    //
                    i1 = (j + 1) * nb + k;
                    i2 = ipiv[((j + 1) * nb + k) - 1];
                    if (i1 != i2) {
                        //                    > Apply pivots to previous columns of L
                        Rswap(k - 1, &a[(((j + 1) * nb + 1) - 1) + (i1 - 1) * lda], 1, &a[(((j + 1) * nb + 1) - 1) + (i2 - 1) * lda], 1);
                        //                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                        if (i2 > (i1 + 1)) {
                            Rswap(i2 - i1 - 1, &a[(i1 - 1) + ((i1 + 1) - 1) * lda], lda, &a[((i1 + 1) - 1) + (i2 - 1) * lda], 1);
                        }
                        //                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                        if (i2 < n) {
                            Rswap(n - i2, &a[(i1 - 1) + ((i2 + 1) - 1) * lda], lda, &a[(i2 - 1) + ((i2 + 1) - 1) * lda], lda);
                        }
                        //                    > Swap A(I1, I1) with A(I2, I2)
                        piv = a[(i1 - 1) + (i1 - 1) * lda];
                        a[(i1 - 1) + (i1 - 1) * lda] = a[(i2 - 1) + (i2 - 1) * lda];
                        a[(i2 - 1) + (i2 - 1) * lda] = piv;
                        //                    > Apply pivots to previous columns of L
                        if (j > 0) {
                            Rswap(j * nb, &a[(i1 - 1) * lda], 1, &a[(i2 - 1) * lda], 1);
                        }
                    }
                }
            }
        }
    } else {
        //
        //        .....................................................
        //        Factorize A as L*D*L**T using the lower triangle of A
        //        .....................................................
        //
        for (j = 0; j <= nt - 1; j = j + 1) {
            //
            //           Generate Jth column of W and H
            //
            kb = min(nb, n - j * nb);
            for (i = 1; i <= j - 1; i = i + 1) {
                if (i == 1) {
                    //                  H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)'
                    if (i == j - 1) {
                        jb = nb + kb;
                    } else {
                        jb = 2 * nb;
                    }
                    Rgemm("NoTranspose", "Transpose", nb, kb, jb, one, &tb[(td + 1 + (i * nb) * ldtb) - 1], ldtb - 1, &a[((j * nb + 1) - 1) + (((i - 1) * nb + 1) - 1) * lda], lda, zero, &work[(i * nb + 1) - 1], n);
                } else {
                    //                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)'
                    if (i == j - 1) {
                        jb = 2 * nb + kb;
                    } else {
                        jb = 3 * nb;
                    }
                    Rgemm("NoTranspose", "Transpose", nb, kb, jb, one, &tb[(td + nb + 1 + ((i - 1) * nb) * ldtb) - 1], ldtb - 1, &a[((j * nb + 1) - 1) + (((i - 2) * nb + 1) - 1) * lda], lda, zero, &work[(i * nb + 1) - 1], n);
                }
            }
            //
            //           Compute T(J,J)
            //
            Rlacpy("Lower", kb, kb, &a[((j * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, &
tb[(td + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
            if (j > 1) {
                //              T(J,J) = L(J,1:J)*H(1:J)
                Rgemm("NoTranspose", "NoTranspose", kb, kb, (j - 1) * nb, -one, &a[((j * nb + 1) - 1)], lda, &work[(nb + 1) - 1], n, one, &tb[(td + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
                //              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)'
                Rgemm("NoTranspose", "NoTranspose", kb, nb, kb, one, &a[((j * nb + 1) - 1) + (((j - 1) * nb + 1) - 1) * lda], lda, &tb[(td + nb + 1 + ((j - 1) * nb) * ldtb) - 1], ldtb - 1, zero, &work[1 - 1], n);
                Rgemm("NoTranspose", "Transpose", kb, kb, nb, -one, &work[1 - 1], n, &a[((j * nb + 1) - 1) + (((j - 2) * nb + 1) - 1) * lda], lda, one, &tb[(td + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
            }
            if (j > 0) {
                Rsygst(1, "Lower", kb, &tb[(td + 1 + (j * nb) * ldtb) - 1], ldtb - 1, &a[((j * nb + 1) - 1) + (((j - 1) * nb + 1) - 1) * lda], lda, iinfo);
            }
            //
            //           Expand T(J,J) into full format
            //
            for (i = 1; i <= kb; i = i + 1) {
                for (k = i + 1; k <= kb; k = k + 1) {
                    tb[(td - (k - (i + 1)) + (j * nb + k - 1) * ldtb) - 1] = tb[(td + (k - i) + 1 + (j * nb + i - 1) * ldtb) - 1];
                }
            }
            //
            if (j < nt - 1) {
                if (j > 0) {
                    //
                    //                 Compute H(J,J)
                    //
                    if (j == 1) {
                        Rgemm("NoTranspose", "Transpose", kb, kb, kb, one, &tb[(td + 1 + (j * nb) * ldtb) - 1], ldtb - 1, &a[((j * nb + 1) - 1) + (((j - 1) * nb + 1) - 1) * lda], lda, zero, &work[(j * nb + 1) - 1], n);
                    } else {
                        Rgemm("NoTranspose", "Transpose", kb, kb, nb + kb, one, &tb[(td + nb + 1 + ((j - 1) * nb) * ldtb) - 1], ldtb - 1, &a[((j * nb + 1) - 1) + (((j - 2) * nb + 1) - 1) * lda], lda, zero, &work[(j * nb + 1) - 1], n);
                    }
                    //
                    //                 Update with the previous column
                    //
                    Rgemm("NoTranspose", "NoTranspose", n - (j + 1) * nb, nb, j * nb, -one, &a[(((j + 1) * nb + 1) - 1)], lda, &work[(nb + 1) - 1], n, one, &a[(((j + 1) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda);
                }
                //
                //              Factorize panel
                //
                Rgetrf(n - (j + 1) * nb, nb, &a[(((j + 1) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, &ipiv[((j + 1) * nb + 1) - 1], iinfo);
                //               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN
                //                  INFO = IINFO+(J+1)*NB
                //               END IF
                //
                //              Compute T(J+1, J), zero out for GEMM update
                //
                kb = min(nb, n - (j + 1) * nb);
                Rlaset("Full", kb, nb, zero, zero, &tb[(td + nb + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
                Rlacpy("Upper", kb, nb, &a[(((j + 1) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda, &tb[(td + nb + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
                if (j > 0) {
                    Rtrsm("R", "L", "T", "U", kb, nb, one, &a[((j * nb + 1) - 1) + (((j - 1) * nb + 1) - 1) * lda], lda, &tb[(td + nb + 1 + (j * nb) * ldtb) - 1], ldtb - 1);
                }
                //
                //              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM
                //              updates
                //
                for (k = 1; k <= nb; k = k + 1) {
                    for (i = 1; i <= kb; i = i + 1) {
                        tb[(td - nb + k - i + 1 + (j * nb + nb + i - 1) * ldtb) - 1] = tb[(td + nb + i - k + 1 + (j * nb + k - 1) * ldtb) - 1];
                    }
                }
                Rlaset("Upper", kb, nb, zero, one, &a[(((j + 1) * nb + 1) - 1) + ((j * nb + 1) - 1) * lda], lda);
                //
                //              Apply pivots to trailing submatrix of A
                //
                for (k = 1; k <= kb; k = k + 1) {
                    //                 > Adjust ipiv
                    ipiv[((j + 1) * nb + k) - 1] += (j + 1) * nb;
                    //
                    i1 = (j + 1) * nb + k;
                    i2 = ipiv[((j + 1) * nb + k) - 1];
                    if (i1 != i2) {
                        //                    > Apply pivots to previous columns of L
                        Rswap(k - 1, &a[(i1 - 1) + (((j + 1) * nb + 1) - 1) * lda], lda, &a[(i2 - 1) + (((j + 1) * nb + 1) - 1) * lda], lda);
                        //                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                        if (i2 > (i1 + 1)) {
                            Rswap(i2 - i1 - 1, &a[((i1 + 1) - 1) + (i1 - 1) * lda], 1, &a[(i2 - 1) + ((i1 + 1) - 1) * lda], lda);
                        }
                        //                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                        if (i2 < n) {
                            Rswap(n - i2, &a[((i2 + 1) - 1) + (i1 - 1) * lda], 1, &a[((i2 + 1) - 1) + (i2 - 1) * lda], 1);
                        }
                        //                    > Swap A(I1, I1) with A(I2, I2)
                        piv = a[(i1 - 1) + (i1 - 1) * lda];
                        a[(i1 - 1) + (i1 - 1) * lda] = a[(i2 - 1) + (i2 - 1) * lda];
                        a[(i2 - 1) + (i2 - 1) * lda] = piv;
                        //                    > Apply pivots to previous columns of L
                        if (j > 0) {
                            Rswap(j * nb, &a[(i1 - 1)], lda, &a[(i2 - 1)], lda);
                        }
                    }
                }
                //
                //              Apply pivots to previous columns of L
                //
                //               CALL Rlaswp( J*NB, A( 1, 1 ), LDA,
                //     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 )
            }
        }
    }
    //
    //     Factor the band matrix
    Rgbtrf(n, n, nb, nb, tb, ldtb, ipiv2, info);
    //
    //     End of Rsytrf_aa_2stage
    //
}
