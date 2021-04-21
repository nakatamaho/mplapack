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

void Cpbtrf(const char *uplo, INTEGER const n, INTEGER const kd, COMPLEX *ab, INTEGER const ldab, INTEGER &info) {
    INTEGER nb = 0;
    const INTEGER nbmax = 32;
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    const INTEGER ldwork = nbmax + 1;
    INTEGER ib = 0;
    INTEGER ii = 0;
    INTEGER i2 = 0;
    INTEGER i3 = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const REAL one = 1.0;
    INTEGER jj = 0;
    COMPLEX work[ldwork * nbmax];
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
    //     Test the input parameters.
    //
    info = 0;
    if ((!Mlsame(uplo, "U")) && (!Mlsame(uplo, "L"))) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (kd < 0) {
        info = -3;
    } else if (ldab < kd + 1) {
        info = -5;
    }
    if (info != 0) {
        Mxerbla("Cpbtrf", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Determine the block size for this environment
    //
    nb = iMlaenv(1, "Cpbtrf", uplo, n, kd, -1, -1);
    //
    //     The block size must not exceed the semi-bandwidth KD, and must not
    //     exceed the limit set by the size of the local array WORK.
    //
    nb = min(nb, nbmax);
    //
    if (nb <= 1 || nb > kd) {
        //
        //        Use unblocked code
        //
        Cpbtf2(uplo, n, kd, ab, ldab, info);
    } else {
        //
        //        Use blocked code
        //
        if (Mlsame(uplo, "U")) {
            //
            //           Compute the Cholesky factorization of a Hermitian band
            //           matrix, given the upper triangle of the matrix in band
            //           storage.
            //
            //           Zero the upper triangle of the work array.
            //
            for (j = 1; j <= nb; j = j + 1) {
                for (i = 1; i <= j - 1; i = i + 1) {
                    work[(i - 1) + (j - 1) * ldwork] = zero;
                }
            }
            //
            //           Process the band matrix one diagonal block at a time.
            //
            for (i = 1; i <= n; i = i + nb) {
                ib = min(nb, n - i + 1);
                //
                //              Factorize the diagonal block
                //
                Cpotf2(uplo, ib, &ab[((kd + 1) - 1) + (i - 1) * ldab], ldab - 1, ii);
                if (ii != 0) {
                    info = i + ii - 1;
                    goto statement_150;
                }
                if (i + ib <= n) {
                    //
                    //                 Update the relevant part of the trailing submatrix.
                    //                 If A11 denotes the diagonal block which has just been
                    //                 factorized, then we need to update the remaining
                    //                 blocks in the diagram:
                    //
                    //                    A11   A12   A13
                    //                          A22   A23
                    //                                A33
                    //
                    //                 The numbers of rows and columns in the partitioning
                    //                 are IB, I2, I3 respectively. The blocks A12, A22 and
                    //                 A23 are empty if IB = KD. The upper triangle of A13
                    //                 lies outside the band.
                    //
                    i2 = min(kd - ib, n - i - ib + 1);
                    i3 = min(ib, n - i - kd + 1);
                    //
                    if (i2 > 0) {
                        //
                        //                    Update A12
                        //
                        Ctrsm("Left", "Upper", "Conjugate transpose", "Non-unit", ib, i2, cone, &ab[((kd + 1) - 1) + (i - 1) * ldab], ldab - 1, &ab[((kd + 1 - ib) - 1) + ((i + ib) - 1) * ldab], ldab - 1);
                        //
                        //                    Update A22
                        //
                        Cherk("Upper", "Conjugate transpose", i2, ib, -one, &ab[((kd + 1 - ib) - 1) + ((i + ib) - 1) * ldab], ldab - 1, one, &ab[((kd + 1) - 1) + ((i + ib) - 1) * ldab], ldab - 1);
                    }
                    //
                    if (i3 > 0) {
                        //
                        //                    Copy the lower triangle of A13 into the work array.
                        //
                        for (jj = 1; jj <= i3; jj = jj + 1) {
                            for (ii = jj; ii <= ib; ii = ii + 1) {
                                work[(ii - 1) + (jj - 1) * ldwork] = ab[((ii - jj + 1) - 1) + ((jj + i + kd - 1) - 1) * ldab];
                            }
                        }
                        //
                        //                    Update A13 (in the work array).
                        //
                        Ctrsm("Left", "Upper", "Conjugate transpose", "Non-unit", ib, i3, cone, &ab[((kd + 1) - 1) + (i - 1) * ldab], ldab - 1, work, ldwork);
                        //
                        //                    Update A23
                        //
                        if (i2 > 0) {
                            Cgemm("Conjugate transpose", "No transpose", i2, i3, ib, -cone, &ab[((kd + 1 - ib) - 1) + ((i + ib) - 1) * ldab], ldab - 1, work, ldwork, cone, &ab[((1 + ib) - 1) + ((i + kd) - 1) * ldab], ldab - 1);
                        }
                        //
                        //                    Update A33
                        //
                        Cherk("Upper", "Conjugate transpose", i3, ib, -one, work, ldwork, one, &ab[((kd + 1) - 1) + ((i + kd) - 1) * ldab], ldab - 1);
                        //
                        //                    Copy the lower triangle of A13 back into place.
                        //
                        for (jj = 1; jj <= i3; jj = jj + 1) {
                            for (ii = jj; ii <= ib; ii = ii + 1) {
                                ab[((ii - jj + 1) - 1) + ((jj + i + kd - 1) - 1) * ldab] = work[(ii - 1) + (jj - 1) * ldwork];
                            }
                        }
                    }
                }
            }
        } else {
            //
            //           Compute the Cholesky factorization of a Hermitian band
            //           matrix, given the lower triangle of the matrix in band
            //           storage.
            //
            //           Zero the lower triangle of the work array.
            //
            for (j = 1; j <= nb; j = j + 1) {
                for (i = j + 1; i <= nb; i = i + 1) {
                    work[(i - 1) + (j - 1) * ldwork] = zero;
                }
            }
            //
            //           Process the band matrix one diagonal block at a time.
            //
            for (i = 1; i <= n; i = i + nb) {
                ib = min(nb, n - i + 1);
                //
                //              Factorize the diagonal block
                //
                Cpotf2(uplo, ib, &ab[(i - 1) * ldab], ldab - 1, ii);
                if (ii != 0) {
                    info = i + ii - 1;
                    goto statement_150;
                }
                if (i + ib <= n) {
                    //
                    //                 Update the relevant part of the trailing submatrix.
                    //                 If A11 denotes the diagonal block which has just been
                    //                 factorized, then we need to update the remaining
                    //                 blocks in the diagram:
                    //
                    //                    A11
                    //                    A21   A22
                    //                    A31   A32   A33
                    //
                    //                 The numbers of rows and columns in the partitioning
                    //                 are IB, I2, I3 respectively. The blocks A21, A22 and
                    //                 A32 are empty if IB = KD. The lower triangle of A31
                    //                 lies outside the band.
                    //
                    i2 = min(kd - ib, n - i - ib + 1);
                    i3 = min(ib, n - i - kd + 1);
                    //
                    if (i2 > 0) {
                        //
                        //                    Update A21
                        //
                        Ctrsm("Right", "Lower", "Conjugate transpose", "Non-unit", i2, ib, cone, &ab[(i - 1) * ldab], ldab - 1, &ab[((1 + ib) - 1) + (i - 1) * ldab], ldab - 1);
                        //
                        //                    Update A22
                        //
                        Cherk("Lower", "No transpose", i2, ib, -one, &ab[((1 + ib) - 1) + (i - 1) * ldab], ldab - 1, one, &ab[((i + ib) - 1) * ldab], ldab - 1);
                    }
                    //
                    if (i3 > 0) {
                        //
                        //                    Copy the upper triangle of A31 into the work array.
                        //
                        for (jj = 1; jj <= ib; jj = jj + 1) {
                            for (ii = 1; ii <= min(jj, i3); ii = ii + 1) {
                                work[(ii - 1) + (jj - 1) * ldwork] = ab[((kd + 1 - jj + ii) - 1) + ((jj + i - 1) - 1) * ldab];
                            }
                        }
                        //
                        //                    Update A31 (in the work array).
                        //
                        Ctrsm("Right", "Lower", "Conjugate transpose", "Non-unit", i3, ib, cone, &ab[(i - 1) * ldab], ldab - 1, work, ldwork);
                        //
                        //                    Update A32
                        //
                        if (i2 > 0) {
                            Cgemm("No transpose", "Conjugate transpose", i3, i2, ib, -cone, work, ldwork, &ab[((1 + ib) - 1) + (i - 1) * ldab], ldab - 1, cone, &ab[((1 + kd - ib) - 1) + ((i + ib) - 1) * ldab], ldab - 1);
                        }
                        //
                        //                    Update A33
                        //
                        Cherk("Lower", "No transpose", i3, ib, -one, work, ldwork, one, &ab[((i + kd) - 1) * ldab], ldab - 1);
                        //
                        //                    Copy the upper triangle of A31 back into place.
                        //
                        for (jj = 1; jj <= ib; jj = jj + 1) {
                            for (ii = 1; ii <= min(jj, i3); ii = ii + 1) {
                                ab[((kd + 1 - jj + ii) - 1) + ((jj + i - 1) - 1) * ldab] = work[(ii - 1) + (jj - 1) * ldwork];
                            }
                        }
                    }
                }
            }
        }
    }
    return;
//
statement_150:;
    //
    //     End of Cpbtrf
    //
}
