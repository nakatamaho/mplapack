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

void Rgbtrf(INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, REAL *ab, INTEGER const ldab, INTEGER *ipiv, INTEGER &info) {
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     KV is the number of superdiagonals in the factor U, allowing for
    //     fill-in
    //
    INTEGER kv = ku + kl;
    //
    //     Test the input parameters.
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (kl < 0) {
        info = -3;
    } else if (ku < 0) {
        info = -4;
    } else if (ldab < kl + kv + 1) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Rgbtrf", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     Determine the block size for this environment
    //
    INTEGER nb = iMlaenv(1, "Rgbtrf", " ", m, n, kl, ku);
    //
    //     The block size must not exceed the limit set by the size of the
    //     local arrays WORK13 and WORK31.
    //
    const INTEGER nbmax = 64;
    nb = min(nb, nbmax);
    //
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    const INTEGER ldwork = nbmax + 1;
    REAL work13[ldwork * nbmax];
    INTEGER ldwork13 = ldwork;
    REAL work31[ldwork * nbmax];
    INTEGER ldwork31 = ldwork;
    INTEGER ju = 0;
    INTEGER jb = 0;
    INTEGER i2 = 0;
    INTEGER i3 = 0;
    INTEGER jj = 0;
    INTEGER km = 0;
    INTEGER jp = 0;
    const REAL one = 1.0;
    INTEGER jm = 0;
    INTEGER nw = 0;
    INTEGER j2 = 0;
    INTEGER j3 = 0;
    INTEGER k2 = 0;
    INTEGER ii = 0;
    INTEGER ip = 0;
    REAL temp = 0.0;
    if (nb <= 1 || nb > kl) {
        //
        //        Use unblocked code
        //
        Rgbtf2(m, n, kl, ku, ab, ldab, ipiv, info);
    } else {
        //
        //        Use blocked code
        //
        //        Zero the superdiagonal elements of the work array WORK13
        //
        for (j = 1; j <= nb; j = j + 1) {
            for (i = 1; i <= j - 1; i = i + 1) {
                work13[(i - 1) + (j - 1) * ldwork13] = zero;
            }
        }
        //
        //        Zero the subdiagonal elements of the work array WORK31
        //
        for (j = 1; j <= nb; j = j + 1) {
            for (i = j + 1; i <= nb; i = i + 1) {
                work31[(i - 1) + (j - 1) * ldwork31] = zero;
            }
        }
        //
        //        Gaussian elimination with partial pivoting
        //
        //        Set fill-in elements in columns KU+2 to KV to zero
        //
        for (j = ku + 2; j <= min(kv, n); j = j + 1) {
            for (i = kv - j + 2; i <= kl; i = i + 1) {
                ab[(i - 1) + (j - 1) * ldab] = zero;
            }
        }
        //
        //        JU is the index of the last column affected by the current
        //        stage of the factorization
        //
        ju = 1;
        //
        for (j = 1; j <= min(m, n); j = j + nb) {
            jb = min(nb, min(m, n) - j + 1);
            //
            //           The active part of the matrix is partitioned
            //
            //              A11   A12   A13
            //              A21   A22   A23
            //              A31   A32   A33
            //
            //           Here A11, A21 and A31 denote the current block of JB columns
            //           which is about to be factorized. The number of rows in the
            //           partitioning are JB, I2, I3 respectively, and the numbers
            //           of columns are JB, J2, J3. The superdiagonal elements of A13
            //           and the subdiagonal elements of A31 lie outside the band.
            //
            i2 = min(kl - jb, m - j - jb + 1);
            i3 = min(jb, m - j - kl + 1);
            //
            //           J2 and J3 are computed after JU has been updated.
            //
            //           Factorize the current block of JB columns
            //
            for (jj = j; jj <= j + jb - 1; jj = jj + 1) {
                //
                //              Set fill-in elements in column JJ+KV to zero
                //
                if (jj + kv <= n) {
                    for (i = 1; i <= kl; i = i + 1) {
                        ab[(i - 1) + ((jj + kv) - 1) * ldab] = zero;
                    }
                }
                //
                //              Find pivot and test for singularity. KM is the number of
                //              subdiagonal elements in the current column.
                //
                km = min(kl, m - jj);
                jp = iRamax(km + 1, &ab[((kv + 1) - 1) + (jj - 1) * ldab], 1);
                ipiv[jj - 1] = jp + jj - j;
                if (ab[((kv + jp) - 1) + (jj - 1) * ldab] != zero) {
                    ju = max(ju, min(jj + ku + jp - 1, n));
                    if (jp != 1) {
                        //
                        //                    Apply interchange to columns J to J+JB-1
                        //
                        if (jp + jj - 1 < j + kl) {
                            //
                            Rswap(jb, &ab[((kv + 1 + jj - j) - 1) + (j - 1) * ldab], ldab - 1, &ab[((kv + jp + jj - j) - 1) + (j - 1) * ldab], ldab - 1);
                        } else {
                            //
                            //                       The interchange affects columns J to JJ-1 of A31
                            //                       which are stored in the work array WORK31
                            //
                            Rswap(jj - j, &ab[((kv + 1 + jj - j) - 1) + (j - 1) * ldab], ldab - 1, &work31[((jp + jj - j - kl) - 1)], ldwork);
                            Rswap(j + jb - jj, &ab[((kv + 1) - 1) + (jj - 1) * ldab], ldab - 1, &ab[((kv + jp) - 1) + (jj - 1) * ldab], ldab - 1);
                        }
                    }
                    //
                    //                 Compute multipliers
                    //
                    Rscal(km, one / ab[((kv + 1) - 1) + (jj - 1) * ldab], &ab[((kv + 2) - 1) + (jj - 1) * ldab], 1);
                    //
                    //                 Update trailing submatrix within the band and within
                    //                 the current block. JM is the index of the last column
                    //                 which needs to be updated.
                    //
                    jm = min(ju, j + jb - 1);
                    if (jm > jj) {
                        Rger(km, jm - jj, -one, &ab[((kv + 2) - 1) + (jj - 1) * ldab], 1, &ab[(kv - 1) + ((jj + 1) - 1) * ldab], ldab - 1, &ab[((kv + 1) - 1) + ((jj + 1) - 1) * ldab], ldab - 1);
                    }
                } else {
                    //
                    //                 If pivot is zero, set INFO to the index of the pivot
                    //                 unless a zero pivot has already been found.
                    //
                    if (info == 0) {
                        info = jj;
                    }
                }
                //
                //              Copy current column of A31 into the work array WORK31
                //
                nw = min(jj - j + 1, i3);
                if (nw > 0) {
                    Rcopy(nw, &ab[((kv + kl + 1 - jj + j) - 1) + (jj - 1) * ldab], 1, &work31[((jj - j + 1) - 1) * ldwork31], 1);
                }
            }
            if (j + jb <= n) {
                //
                //              Apply the row interchanges to the other blocks.
                //
                j2 = min(ju - j + 1, kv) - jb;
                j3 = max(0, ju - j - kv + 1);
                //
                //              Use Rlaswp to apply the row interchanges to A12, A22, and
                //              A32.
                //
                Rlaswp(j2, &ab[((kv + 1 - jb) - 1) + ((j + jb) - 1) * ldab], ldab - 1, 1, jb, &ipiv[j - 1], 1);
                //
                //              Adjust the pivot indices.
                //
                for (i = j; i <= j + jb - 1; i = i + 1) {
                    ipiv[i - 1] += j - 1;
                }
                //
                //              Apply the row interchanges to A13, A23, and A33
                //              columnwise.
                //
                k2 = j - 1 + jb + j2;
                for (i = 1; i <= j3; i = i + 1) {
                    jj = k2 + i;
                    for (ii = j + i - 1; ii <= j + jb - 1; ii = ii + 1) {
                        ip = ipiv[ii - 1];
                        if (ip != ii) {
                            temp = ab[((kv + 1 + ii - jj) - 1) + (jj - 1) * ldab];
                            ab[((kv + 1 + ii - jj) - 1) + (jj - 1) * ldab] = ab[((kv + 1 + ip - jj) - 1) + (jj - 1) * ldab];
                            ab[((kv + 1 + ip - jj) - 1) + (jj - 1) * ldab] = temp;
                        }
                    }
                }
                //
                //              Update the relevant part of the trailing submatrix
                //
                if (j2 > 0) {
                    //
                    //                 Update A12
                    //
                    Rtrsm("Left", "Lower", "No transpose", "Unit", jb, j2, one, &ab[((kv + 1) - 1) + (j - 1) * ldab], ldab - 1, &ab[((kv + 1 - jb) - 1) + ((j + jb) - 1) * ldab], ldab - 1);
                    //
                    if (i2 > 0) {
                        //
                        //                    Update A22
                        //
                        Rgemm("No transpose", "No transpose", i2, j2, jb, -one, &ab[((kv + 1 + jb) - 1) + (j - 1) * ldab], ldab - 1, &ab[((kv + 1 - jb) - 1) + ((j + jb) - 1) * ldab], ldab - 1, one, &ab[((kv + 1) - 1) + ((j + jb) - 1) * ldab], ldab - 1);
                    }
                    //
                    if (i3 > 0) {
                        //
                        //                    Update A32
                        //
                        Rgemm("No transpose", "No transpose", i3, j2, jb, -one, work31, ldwork, &ab[((kv + 1 - jb) - 1) + ((j + jb) - 1) * ldab], ldab - 1, one, &ab[((kv + kl + 1 - jb) - 1) + ((j + jb) - 1) * ldab], ldab - 1);
                    }
                }
                //
                if (j3 > 0) {
                    //
                    //                 Copy the lower triangle of A13 into the work array
                    //                 WORK13
                    //
                    for (jj = 1; jj <= j3; jj = jj + 1) {
                        for (ii = jj; ii <= jb; ii = ii + 1) {
                            work13[(ii - 1) + (jj - 1) * ldwork13] = ab[((ii - jj + 1) - 1) + ((jj + j + kv - 1) - 1) * ldab];
                        }
                    }
                    //
                    //                 Update A13 in the work array
                    //
                    Rtrsm("Left", "Lower", "No transpose", "Unit", jb, j3, one, &ab[((kv + 1) - 1) + (j - 1) * ldab], ldab - 1, work13, ldwork);
                    //
                    if (i2 > 0) {
                        //
                        //                    Update A23
                        //
                        Rgemm("No transpose", "No transpose", i2, j3, jb, -one, &ab[((kv + 1 + jb) - 1) + (j - 1) * ldab], ldab - 1, work13, ldwork, one, &ab[((1 + jb) - 1) + ((j + kv) - 1) * ldab], ldab - 1);
                    }
                    //
                    if (i3 > 0) {
                        //
                        //                    Update A33
                        //
                        Rgemm("No transpose", "No transpose", i3, j3, jb, -one, work31, ldwork, work13, ldwork, one, &ab[((1 + kl) - 1) + ((j + kv) - 1) * ldab], ldab - 1);
                    }
                    //
                    //                 Copy the lower triangle of A13 back into place
                    //
                    for (jj = 1; jj <= j3; jj = jj + 1) {
                        for (ii = jj; ii <= jb; ii = ii + 1) {
                            ab[((ii - jj + 1) - 1) + ((jj + j + kv - 1) - 1) * ldab] = work13[(ii - 1) + (jj - 1) * ldwork13];
                        }
                    }
                }
            } else {
                //
                //              Adjust the pivot indices.
                //
                for (i = j; i <= j + jb - 1; i = i + 1) {
                    ipiv[i - 1] += j - 1;
                }
            }
            //
            //           Partially undo the interchanges in the current block to
            //           restore the upper triangular form of A31 and copy the upper
            //           triangle of A31 back into place
            //
            for (jj = j + jb - 1; jj >= j; jj = jj - 1) {
                jp = ipiv[jj - 1] - jj + 1;
                if (jp != 1) {
                    //
                    //                 Apply interchange to columns J to JJ-1
                    //
                    if (jp + jj - 1 < j + kl) {
                        //
                        //                    The interchange does not affect A31
                        //
                        Rswap(jj - j, &ab[((kv + 1 + jj - j) - 1) + (j - 1) * ldab], ldab - 1, &ab[((kv + jp + jj - j) - 1) + (j - 1) * ldab], ldab - 1);
                    } else {
                        //
                        //                    The interchange does affect A31
                        //
                        Rswap(jj - j, &ab[((kv + 1 + jj - j) - 1) + (j - 1) * ldab], ldab - 1, &work31[((jp + jj - j - kl) - 1)], ldwork);
                    }
                }
                //
                //              Copy the current column of A31 back into place
                //
                nw = min(i3, jj - j + 1);
                if (nw > 0) {
                    Rcopy(nw, &work31[((jj - j + 1) - 1) * ldwork31], 1, &ab[((kv + kl + 1 - jj + j) - 1) + (jj - 1) * ldab], 1);
                }
            }
        }
    }
    //
    //     End of Rgbtrf
    //
}
