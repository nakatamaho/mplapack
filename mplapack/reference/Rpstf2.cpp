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

void Rpstf2(const char *uplo, INTEGER const &n, REAL *a, INTEGER const &lda, arr_ref<INTEGER> piv, INTEGER &rank, REAL const &tol, REAL *work, INTEGER &info) {
    bool upper = false;
    INTEGER i = 0;
    INTEGER pvt = 0;
    REAL ajj = 0.0;
    const REAL zero = 0.0;
    REAL dstop = 0.0;
    INTEGER j = 0;
    INTEGER itemp = 0;
    REAL dtemp = 0.0;
    const REAL one = 1.0;
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Rpstf2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Initialize PIV
    //
    for (i = 1; i <= n; i = i + 1) {
        piv[i - 1] = i;
    }
    //
    //     Compute stopping value
    //
    pvt = 1;
    ajj = a[(pvt - 1) + (pvt - 1) * lda];
    for (i = 2; i <= n; i = i + 1) {
        if (a[(i - 1) + (i - 1) * lda] > ajj) {
            pvt = i;
            ajj = a[(pvt - 1) + (pvt - 1) * lda];
        }
    }
    if (ajj <= zero || Risnan(ajj)) {
        rank = 0;
        info = 1;
        goto statement_170;
    }
    //
    //     Compute stopping value if not supplied
    //
    if (tol < zero) {
        dstop = n * dlamch("Epsilon") * ajj;
    } else {
        dstop = tol;
    }
    //
    //     Set first half of WORK to zero, holds dot products
    //
    for (i = 1; i <= n; i = i + 1) {
        work[i - 1] = 0;
    }
    //
    if (upper) {
        //
        //        Compute the Cholesky factorization P**T * A * P = U**T * U
        //
        for (j = 1; j <= n; j = j + 1) {
            //
            //        Find pivot, test for exit, else swap rows and columns
            //        Update dot products, compute possible pivots which are
            //        stored in the second half of WORK
            //
            for (i = j; i <= n; i = i + 1) {
                //
                if (j > 1) {
                    work[i - 1] += pow2(a[((j - 1) - 1) + (i - 1) * lda]);
                }
                work[(n + i) - 1] = a[(i - 1) + (i - 1) * lda] - work[i - 1];
                //
            }
            //
            if (j > 1) {
                itemp = maxloc[((work[(((n + j)) - 1) + (((2 * n)) - 1) * ldwork]) - 1)];
                pvt = itemp + j - 1;
                ajj = work[(n + pvt) - 1];
                if (ajj <= dstop || Risnan(ajj)) {
                    a[(j - 1) + (j - 1) * lda] = ajj;
                    goto statement_160;
                }
            }
            //
            if (j != pvt) {
                //
                //              Pivot OK, so can now swap pivot rows and columns
                //
                a[(pvt - 1) + (pvt - 1) * lda] = a[(j - 1) + (j - 1) * lda];
                Rswap(j - 1, a[(j - 1) * lda], 1, a[(pvt - 1) * lda], 1);
                if (pvt < n) {
                    Rswap(n - pvt, a[(j - 1) + ((pvt + 1) - 1) * lda], lda, a[(pvt - 1) + ((pvt + 1) - 1) * lda], lda);
                }
                Rswap(pvt - j - 1, a[(j - 1) + ((j + 1) - 1) * lda], lda, a[((j + 1) - 1) + (pvt - 1) * lda], 1);
                //
                //              Swap dot products and PIV
                //
                dtemp = work[j - 1];
                work[j - 1] = work[pvt - 1];
                work[pvt - 1] = dtemp;
                itemp = piv[pvt - 1];
                piv[pvt - 1] = piv[j - 1];
                piv[j - 1] = itemp;
            }
            //
            ajj = sqrt(ajj);
            a[(j - 1) + (j - 1) * lda] = ajj;
            //
            //           Compute elements J+1:N of row J
            //
            if (j < n) {
                Rgemv("Trans", j - 1, n - j, -one, a[((j + 1) - 1) * lda], lda, a[(j - 1) * lda], 1, one, a[(j - 1) + ((j + 1) - 1) * lda], lda);
                Rscal(n - j, one / ajj, a[(j - 1) + ((j + 1) - 1) * lda], lda);
            }
            //
        }
        //
    } else {
        //
        //        Compute the Cholesky factorization P**T * A * P = L * L**T
        //
        for (j = 1; j <= n; j = j + 1) {
            //
            //        Find pivot, test for exit, else swap rows and columns
            //        Update dot products, compute possible pivots which are
            //        stored in the second half of WORK
            //
            for (i = j; i <= n; i = i + 1) {
                //
                if (j > 1) {
                    work[i - 1] += pow2(a[(i - 1) + ((j - 1) - 1) * lda]);
                }
                work[(n + i) - 1] = a[(i - 1) + (i - 1) * lda] - work[i - 1];
                //
            }
            //
            if (j > 1) {
                itemp = maxloc[((work[(((n + j)) - 1) + (((2 * n)) - 1) * ldwork]) - 1)];
                pvt = itemp + j - 1;
                ajj = work[(n + pvt) - 1];
                if (ajj <= dstop || Risnan(ajj)) {
                    a[(j - 1) + (j - 1) * lda] = ajj;
                    goto statement_160;
                }
            }
            //
            if (j != pvt) {
                //
                //              Pivot OK, so can now swap pivot rows and columns
                //
                a[(pvt - 1) + (pvt - 1) * lda] = a[(j - 1) + (j - 1) * lda];
                Rswap(j - 1, a[(j - 1)], lda, a[(pvt - 1)], lda);
                if (pvt < n) {
                    Rswap(n - pvt, a[((pvt + 1) - 1) + (j - 1) * lda], 1, a[((pvt + 1) - 1) + (pvt - 1) * lda], 1);
                }
                Rswap(pvt - j - 1, a[((j + 1) - 1) + (j - 1) * lda], 1, a[(pvt - 1) + ((j + 1) - 1) * lda], lda);
                //
                //              Swap dot products and PIV
                //
                dtemp = work[j - 1];
                work[j - 1] = work[pvt - 1];
                work[pvt - 1] = dtemp;
                itemp = piv[pvt - 1];
                piv[pvt - 1] = piv[j - 1];
                piv[j - 1] = itemp;
            }
            //
            ajj = sqrt(ajj);
            a[(j - 1) + (j - 1) * lda] = ajj;
            //
            //           Compute elements J+1:N of column J
            //
            if (j < n) {
                Rgemv("No Trans", n - j, j - 1, -one, a[((j + 1) - 1)], lda, a[(j - 1)], lda, one, a[((j + 1) - 1) + (j - 1) * lda], 1);
                Rscal(n - j, one / ajj, a[((j + 1) - 1) + (j - 1) * lda], 1);
            }
            //
        }
        //
    }
    //
    //     Ran to completion, A has full rank
    //
    rank = n;
    //
    goto statement_170;
statement_160:
    //
    //     Rank is number of steps completed.  Set INFO = 1 to signal
    //     that the factorization cannot be used to solve a system.
    //
    rank = j - 1;
    info = 1;
//
statement_170:;
    //
    //     End of Rpstf2
    //
}
