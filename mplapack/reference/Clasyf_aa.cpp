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

void Clasyf_aa(const char *uplo, INTEGER const &j1, INTEGER const &m, INTEGER const &nb, COMPLEX *a, INTEGER const &lda, arr_ref<INTEGER> ipiv, COMPLEX *h, INTEGER const &ldh, COMPLEX *work) {
    INTEGER j = 0;
    INTEGER k1 = 0;
    INTEGER k = 0;
    INTEGER mj = 0;
    const COMPLEX one = 1.0;
    COMPLEX alpha = 0.0;
    INTEGER i2 = 0;
    COMPLEX piv = 0.0;
    INTEGER i1 = 0;
    const COMPLEX zero = 0.0;
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
    j = 1;
    //
    //     K1 is the first column of the panel to be factorized
    //     i.e.,  K1 is 2 for the first block column, and 1 for the rest of the blocks
    //
    k1 = (2 - j1) + 1;
    //
    if (Mlsame(uplo, "U")) {
    //
    //        .....................................................
    //        Factorize A as U**T*D*U using the upper triangle of A
    //        .....................................................
    //
    statement_10:
        if (j > min(m, nb)) {
            goto statement_20;
        }
        //
        //        K is the column to be factorized
        //         when being called from Csytrf_aa,
        //         > for the first block column, J1 is 1, hence J1+J-1 is J,
        //         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1,
        //
        k = j1 + j - 1;
        if (j == m) {
            //
            //            Only need to compute T(J, J)
            //
            mj = 1;
        } else {
            mj = m - j + 1;
        }
        //
        //        H(J:M, J) := A(J, J:M) - H(J:M, 1:(J-1)) * L(J1:(J-1), J),
        //         where H(J:M, J) has been initialized to be A(J, J:M)
        //
        if (k > 2) {
            //
            //        K is the column to be factorized
            //         > for the first block column, K is J, skipping the first two
            //           columns
            //         > for the rest of the columns, K is J+1, skipping only the
            //           first column
            //
            Cgemv("No transpose", mj, j - k1, -one, h[(j - 1) + (k1 - 1) * ldh], ldh, a[(j - 1) * lda], 1, one, h[(j - 1) + (j - 1) * ldh], 1);
        }
        //
        //        Copy H(i:M, i) INTEGERo WORK
        //
        Ccopy(mj, h[(j - 1) + (j - 1) * ldh], 1, work[1 - 1], 1);
        //
        if (j > k1) {
            //
            //           Compute WORK := WORK - L(J-1, J:M) * T(J-1,J),
            //            where A(J-1, J) stores T(J-1, J) and A(J-2, J:M) stores U(J-1, J:M)
            //
            alpha = -a[((k - 1) - 1) + (j - 1) * lda];
            Caxpy(mj, alpha, a[((k - 2) - 1) + (j - 1) * lda], lda, work[1 - 1], 1);
        }
        //
        //        Set A(J, J) = T(J, J)
        //
        a[(k - 1) + (j - 1) * lda] = work[1 - 1];
        //
        if (j < m) {
            //
            //           Compute WORK(2:M) = T(J, J) L(J, (J+1):M)
            //            where A(J, J) stores T(J, J) and A(J-1, (J+1):M) stores U(J, (J+1):M)
            //
            if (k > 1) {
                alpha = -a[(k - 1) + (j - 1) * lda];
                Caxpy(m - j, alpha, a[((k - 1) - 1) + ((j + 1) - 1) * lda], lda, work[2 - 1], 1);
            }
            //
            //           Find max(|WORK(2:M)|)
            //
            i2 = iCamax[((m - j) - 1) + (work[2 - 1] - 1) * ldiCamax] + 1;
            piv = work[i2 - 1];
            //
            //           Apply symmetric pivot
            //
            if ((i2 != 2) && (piv != 0)) {
                //
                //              Swap WORK(I1) and WORK(I2)
                //
                i1 = 2;
                work[i2 - 1] = work[i1 - 1];
                work[i1 - 1] = piv;
                //
                //              Swap A(I1, I1+1:M) with A(I1+1:M, I2)
                //
                i1 += j - 1;
                i2 += j - 1;
                Cswap(i2 - i1 - 1, a[((j1 + i1 - 1) - 1) + ((i1 + 1) - 1) * lda], lda, a[((j1 + i1) - 1) + (i2 - 1) * lda], 1);
                //
                //              Swap A(I1, I2+1:M) with A(I2, I2+1:M)
                //
                if (i2 < m) {
                    Cswap(m - i2, a[((j1 + i1 - 1) - 1) + ((i2 + 1) - 1) * lda], lda, a[((j1 + i2 - 1) - 1) + ((i2 + 1) - 1) * lda], lda);
                }
                //
                //              Swap A(I1, I1) with A(I2,I2)
                //
                piv = a[((i1 + j1 - 1) - 1) + (i1 - 1) * lda];
                a[((j1 + i1 - 1) - 1) + (i1 - 1) * lda] = a[((j1 + i2 - 1) - 1) + (i2 - 1) * lda];
                a[((j1 + i2 - 1) - 1) + (i2 - 1) * lda] = piv;
                //
                //              Swap H(I1, 1:J1) with H(I2, 1:J1)
                //
                Cswap(i1 - 1, h[(i1 - 1)], ldh, h[(i2 - 1)], ldh);
                ipiv[i1 - 1] = i2;
                //
                if (i1 > (k1 - 1)) {
                    //
                    //                 Swap L(1:I1-1, I1) with L(1:I1-1, I2),
                    //                  skipping the first column
                    //
                    Cswap(i1 - k1 + 1, a[(i1 - 1) * lda], 1, a[(i2 - 1) * lda], 1);
                }
            } else {
                ipiv[(j + 1) - 1] = j + 1;
            }
            //
            //           Set A(J, J+1) = T(J, J+1)
            //
            a[(k - 1) + ((j + 1) - 1) * lda] = work[2 - 1];
            //
            if (j < nb) {
                //
                //              Copy A(J+1:M, J+1) INTEGERo H(J:M, J),
                //
                Ccopy(m - j, a[((k + 1) - 1) + ((j + 1) - 1) * lda], lda, h[((j + 1) - 1) + ((j + 1) - 1) * ldh], 1);
            }
            //
            //           Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1),
            //            where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1)
            //
            if (j < (m - 1)) {
                if (a[(k - 1) + ((j + 1) - 1) * lda] != zero) {
                    alpha = one / a[(k - 1) + ((j + 1) - 1) * lda];
                    Ccopy(m - j - 1, work[3 - 1], 1, a[(k - 1) + ((j + 2) - 1) * lda], lda);
                    Cscal(m - j - 1, alpha, a[(k - 1) + ((j + 2) - 1) * lda], lda);
                } else {
                    Claset("Full", 1, m - j - 1, zero, zero, a[(k - 1) + ((j + 2) - 1) * lda], lda);
                }
            }
        }
        j++;
        goto statement_10;
    statement_20:;
        //
    } else {
    //
    //        .....................................................
    //        Factorize A as L*D*L**T using the lower triangle of A
    //        .....................................................
    //
    statement_30:
        if (j > min(m, nb)) {
            goto statement_40;
        }
        //
        //        K is the column to be factorized
        //         when being called from Csytrf_aa,
        //         > for the first block column, J1 is 1, hence J1+J-1 is J,
        //         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1,
        //
        k = j1 + j - 1;
        if (j == m) {
            //
            //            Only need to compute T(J, J)
            //
            mj = 1;
        } else {
            mj = m - j + 1;
        }
        //
        //        H(J:M, J) := A(J:M, J) - H(J:M, 1:(J-1)) * L(J, J1:(J-1))^T,
        //         where H(J:M, J) has been initialized to be A(J:M, J)
        //
        if (k > 2) {
            //
            //        K is the column to be factorized
            //         > for the first block column, K is J, skipping the first two
            //           columns
            //         > for the rest of the columns, K is J+1, skipping only the
            //           first column
            //
            Cgemv("No transpose", mj, j - k1, -one, h[(j - 1) + (k1 - 1) * ldh], ldh, a[(j - 1)], lda, one, h[(j - 1) + (j - 1) * ldh], 1);
        }
        //
        //        Copy H(J:M, J) INTEGERo WORK
        //
        Ccopy(mj, h[(j - 1) + (j - 1) * ldh], 1, work[1 - 1], 1);
        //
        if (j > k1) {
            //
            //           Compute WORK := WORK - L(J:M, J-1) * T(J-1,J),
            //            where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1)
            //
            alpha = -a[(j - 1) + ((k - 1) - 1) * lda];
            Caxpy(mj, alpha, a[(j - 1) + ((k - 2) - 1) * lda], 1, work[1 - 1], 1);
        }
        //
        //        Set A(J, J) = T(J, J)
        //
        a[(j - 1) + (k - 1) * lda] = work[1 - 1];
        //
        if (j < m) {
            //
            //           Compute WORK(2:M) = T(J, J) L((J+1):M, J)
            //            where A(J, J) = T(J, J) and A((J+1):M, J-1) = L((J+1):M, J)
            //
            if (k > 1) {
                alpha = -a[(j - 1) + (k - 1) * lda];
                Caxpy(m - j, alpha, a[((j + 1) - 1) + ((k - 1) - 1) * lda], 1, work[2 - 1], 1);
            }
            //
            //           Find max(|WORK(2:M)|)
            //
            i2 = iCamax[((m - j) - 1) + (work[2 - 1] - 1) * ldiCamax] + 1;
            piv = work[i2 - 1];
            //
            //           Apply symmetric pivot
            //
            if ((i2 != 2) && (piv != 0)) {
                //
                //              Swap WORK(I1) and WORK(I2)
                //
                i1 = 2;
                work[i2 - 1] = work[i1 - 1];
                work[i1 - 1] = piv;
                //
                //              Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                //
                i1 += j - 1;
                i2 += j - 1;
                Cswap(i2 - i1 - 1, a[((i1 + 1) - 1) + ((j1 + i1 - 1) - 1) * lda], 1, a[(i2 - 1) + ((j1 + i1) - 1) * lda], lda);
                //
                //              Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                //
                if (i2 < m) {
                    Cswap(m - i2, a[((i2 + 1) - 1) + ((j1 + i1 - 1) - 1) * lda], 1, a[((i2 + 1) - 1) + ((j1 + i2 - 1) - 1) * lda], 1);
                }
                //
                //              Swap A(I1, I1) with A(I2, I2)
                //
                piv = a[(i1 - 1) + ((j1 + i1 - 1) - 1) * lda];
                a[(i1 - 1) + ((j1 + i1 - 1) - 1) * lda] = a[(i2 - 1) + ((j1 + i2 - 1) - 1) * lda];
                a[(i2 - 1) + ((j1 + i2 - 1) - 1) * lda] = piv;
                //
                //              Swap H(I1, I1:J1) with H(I2, I2:J1)
                //
                Cswap(i1 - 1, h[(i1 - 1)], ldh, h[(i2 - 1)], ldh);
                ipiv[i1 - 1] = i2;
                //
                if (i1 > (k1 - 1)) {
                    //
                    //                 Swap L(1:I1-1, I1) with L(1:I1-1, I2),
                    //                  skipping the first column
                    //
                    Cswap(i1 - k1 + 1, a[(i1 - 1)], lda, a[(i2 - 1)], lda);
                }
            } else {
                ipiv[(j + 1) - 1] = j + 1;
            }
            //
            //           Set A(J+1, J) = T(J+1, J)
            //
            a[((j + 1) - 1) + (k - 1) * lda] = work[2 - 1];
            //
            if (j < nb) {
                //
                //              Copy A(J+1:M, J+1) INTEGERo H(J+1:M, J),
                //
                Ccopy(m - j, a[((j + 1) - 1) + ((k + 1) - 1) * lda], 1, h[((j + 1) - 1) + ((j + 1) - 1) * ldh], 1);
            }
            //
            //           Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1),
            //            where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1)
            //
            if (j < (m - 1)) {
                if (a[((j + 1) - 1) + (k - 1) * lda] != zero) {
                    alpha = one / a[((j + 1) - 1) + (k - 1) * lda];
                    Ccopy(m - j - 1, work[3 - 1], 1, a[((j + 2) - 1) + (k - 1) * lda], 1);
                    Cscal(m - j - 1, alpha, a[((j + 2) - 1) + (k - 1) * lda], 1);
                } else {
                    Claset("Full", m - j - 1, 1, zero, zero, a[((j + 2) - 1) + (k - 1) * lda], lda);
                }
            }
        }
        j++;
        goto statement_30;
    statement_40:;
    }
    //
    //     End of Clasyf_aa
    //
}
