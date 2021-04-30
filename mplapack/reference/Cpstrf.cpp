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

void Cpstrf(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *piv, INTEGER &rank, REAL const tol, REAL *work, INTEGER &info) {
    bool upper = false;
    INTEGER nb = 0;
    INTEGER i = 0;
    INTEGER pvt = 0;
    REAL ajj = 0.0;
    const REAL zero = 0.0;
    REAL dstop = 0.0;
    INTEGER k = 0;
    INTEGER jb = 0;
    INTEGER j = 0;
    INTEGER itemp = 0;
    COMPLEX ztemp = 0.0;
    REAL dtemp = 0.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
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
    //     Test the input parameters.
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
        Mxerbla("Cpstrf", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Get block size
    //
    nb = iMlaenv(1, "Cpotrf", uplo, n, -1, -1, -1);
    if (nb <= 1 || nb >= n) {
        //
        //        Use unblocked code
        //
        Cpstf2(uplo, n, &a[(1 - 1)], lda, piv, rank, tol, work, info);
        goto statement_230;
        //
    } else {
        //
        //     Initialize PIV
        //
        for (i = 1; i <= n; i = i + 1) {
            piv[i - 1] = i;
        }
        //
        //     Compute stopping value
        //
        for (i = 1; i <= n; i = i + 1) {
            work[i - 1] = a[(i - 1) + (i - 1) * lda].real();
        }
        pvt = Mmaxloc(work, 1, n, 1);
        ajj = a[(pvt - 1) + (pvt - 1) * lda].real();
        if (ajj <= zero || Risnan(ajj)) {
            rank = 0;
            info = 1;
            goto statement_230;
        }
        //
        //     Compute stopping value if not supplied
        //
        if (tol < zero) {
            dstop = n * Rlamch("Epsilon") * ajj;
        } else {
            dstop = tol;
        }
        //
        if (upper) {
            //
            //           Compute the Cholesky factorization P**T * A * P = U**H * U
            //
            for (k = 1; k <= n; k = k + nb) {
                //
                //              Account for last block not being NB wide
                //
                jb = min(nb, n - k + 1);
                //
                //              Set relevant part of first half of WORK to zero,
                //              holds dot products
                //
                for (i = k; i <= n; i = i + 1) {
                    work[i - 1] = 0;
                }
                //
                for (j = k; j <= k + jb - 1; j = j + 1) {
                    //
                    //              Find pivot, test for exit, else swap rows and columns
                    //              Update dot products, compute possible pivots which are
                    //              stored in the second half of WORK
                    //
                    for (i = j; i <= n; i = i + 1) {
                        //
                        if (j > k) {
                            work[i - 1] += (conj(a[((j - 1) - 1) + (i - 1) * lda]) * a[((j - 1) - 1) + (i - 1) * lda]).real();
                        }
                        work[(n + i) - 1] = a[(i - 1) + (i - 1) * lda].real() - work[i - 1];
                        //
                    }
                    //
                    if (j > 1) {
                        itemp = Mmaxloc(work, n + j, 2 * n, 1);
                        pvt = itemp + j - 1;
                        ajj = work[(n + pvt) - 1];
                        if (ajj <= dstop || Risnan(ajj)) {
                            a[(j - 1) + (j - 1) * lda] = ajj;
                            goto statement_220;
                        }
                    }
                    //
                    if (j != pvt) {
                        //
                        //                    Pivot OK, so can now swap pivot rows and columns
                        //
                        a[(pvt - 1) + (pvt - 1) * lda] = a[(j - 1) + (j - 1) * lda];
                        Cswap(j - 1, &a[(j - 1) * lda], 1, &a[(pvt - 1) * lda], 1);
                        if (pvt < n) {
                            Cswap(n - pvt, &a[(j - 1) + ((pvt + 1) - 1) * lda], lda, &a[(pvt - 1) + ((pvt + 1) - 1) * lda], lda);
                        }
                        for (i = j + 1; i <= pvt - 1; i = i + 1) {
                            ztemp = conj(a[(j - 1) + (i - 1) * lda]);
                            a[(j - 1) + (i - 1) * lda] = conj(a[(i - 1) + (pvt - 1) * lda]);
                            a[(i - 1) + (pvt - 1) * lda] = ztemp;
                        }
                        a[(j - 1) + (pvt - 1) * lda] = conj(a[(j - 1) + (pvt - 1) * lda]);
                        //
                        //                    Swap dot products and PIV
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
                    //                 Compute elements J+1:N of row J.
                    //
                    if (j < n) {
                        Clacgv(j - 1, &a[(j - 1) * lda], 1);
                        Cgemv("Trans", j - k, n - j, -cone, &a[(k - 1) + ((j + 1) - 1) * lda], lda, &a[(k - 1) + (j - 1) * lda], 1, cone, &a[(j - 1) + ((j + 1) - 1) * lda], lda);
                        Clacgv(j - 1, &a[(j - 1) * lda], 1);
                        CRscal(n - j, one / ajj, &a[(j - 1) + ((j + 1) - 1) * lda], lda);
                    }
                    //
                }
                //
                //              Update trailing matrix, J already incremented
                //
                if (k + jb <= n) {
                    Cherk("Upper", "Conj Trans", n - j + 1, jb, -one, &a[(k - 1) + (j - 1) * lda], lda, one, &a[(j - 1) + (j - 1) * lda], lda);
                }
                //
            }
            //
        } else {
            //
            //        Compute the Cholesky factorization P**T * A * P = L * L**H
            //
            for (k = 1; k <= n; k = k + nb) {
                //
                //              Account for last block not being NB wide
                //
                jb = min(nb, n - k + 1);
                //
                //              Set relevant part of first half of WORK to zero,
                //              holds dot products
                //
                for (i = k; i <= n; i = i + 1) {
                    work[i - 1] = 0;
                }
                //
                for (j = k; j <= k + jb - 1; j = j + 1) {
                    //
                    //              Find pivot, test for exit, else swap rows and columns
                    //              Update dot products, compute possible pivots which are
                    //              stored in the second half of WORK
                    //
                    for (i = j; i <= n; i = i + 1) {
                        //
                        if (j > k) {
                            work[i - 1] += (conj(a[(i - 1) + ((j - 1) - 1) * lda]) * a[(i - 1) + ((j - 1) - 1) * lda]).real();
                        }
                        work[(n + i) - 1] = a[(i - 1) + (i - 1) * lda].real() - work[i - 1];
                        //
                    }
                    //
                    if (j > 1) {
                        itemp = Mmaxloc(work, n + j, 2 * n, 1);
                        pvt = itemp + j - 1;
                        ajj = work[(n + pvt) - 1];
                        if (ajj <= dstop || Risnan(ajj)) {
                            a[(j - 1) + (j - 1) * lda] = ajj;
                            goto statement_220;
                        }
                    }
                    //
                    if (j != pvt) {
                        //
                        //                    Pivot OK, so can now swap pivot rows and columns
                        //
                        a[(pvt - 1) + (pvt - 1) * lda] = a[(j - 1) + (j - 1) * lda];
                        Cswap(j - 1, &a[(j - 1)], lda, &a[(pvt - 1)], lda);
                        if (pvt < n) {
                            Cswap(n - pvt, &a[((pvt + 1) - 1) + (j - 1) * lda], 1, &a[((pvt + 1) - 1) + (pvt - 1) * lda], 1);
                        }
                        for (i = j + 1; i <= pvt - 1; i = i + 1) {
                            ztemp = conj(a[(i - 1) + (j - 1) * lda]);
                            a[(i - 1) + (j - 1) * lda] = conj(a[(pvt - 1) + (i - 1) * lda]);
                            a[(pvt - 1) + (i - 1) * lda] = ztemp;
                        }
                        a[(pvt - 1) + (j - 1) * lda] = conj(a[(pvt - 1) + (j - 1) * lda]);
                        //
                        //                    Swap dot products and PIV
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
                    //                 Compute elements J+1:N of column J.
                    //
                    if (j < n) {
                        Clacgv(j - 1, &a[(j - 1)], lda);
                        Cgemv("No Trans", n - j, j - k, -cone, &a[((j + 1) - 1) + (k - 1) * lda], lda, &a[(j - 1) + (k - 1) * lda], lda, cone, &a[((j + 1) - 1) + (j - 1) * lda], 1);
                        Clacgv(j - 1, &a[(j - 1)], lda);
                        CRscal(n - j, one / ajj, &a[((j + 1) - 1) + (j - 1) * lda], 1);
                    }
                    //
                }
                //
                //              Update trailing matrix, J already incremented
                //
                if (k + jb <= n) {
                    Cherk("Lower", "No Trans", n - j + 1, jb, -one, &a[(j - 1) + (k - 1) * lda], lda, one, &a[(j - 1) + (j - 1) * lda], lda);
                }
                //
            }
            //
        }
    }
    //
    //     Ran to completion, A has full rank
    //
    rank = n;
    //
    goto statement_230;
statement_220:
    //
    //     Rank is the number of steps completed.  Set INFO = 1 to signal
    //     that the factorization cannot be used to solve a system.
    //
    rank = j - 1;
    info = 1;
//
statement_230:;
    //
    //     End of Cpstrf
    //
}
