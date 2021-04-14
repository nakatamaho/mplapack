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

void Cgbtrs(const char *trans, INTEGER const n, INTEGER const kl, INTEGER const ku, INTEGER const nrhs, COMPLEX *ab, INTEGER const ldab, INTEGER *ipiv, COMPLEX *b, INTEGER const ldb, INTEGER &info) {
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
    bool notran = Mlsame(trans, "N");
    if (!notran && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (kl < 0) {
        info = -3;
    } else if (ku < 0) {
        info = -4;
    } else if (nrhs < 0) {
        info = -5;
    } else if (ldab < (2 * kl + ku + 1)) {
        info = -7;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -10;
    }
    if (info != 0) {
        Mxerbla("Cgbtrs", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0 || nrhs == 0) {
        return;
    }
    //
    INTEGER kd = ku + kl + 1;
    bool lnoti = kl > 0;
    //
    INTEGER j = 0;
    INTEGER lm = 0;
    INTEGER l = 0;
    const COMPLEX one = (1.0, 0.0);
    INTEGER i = 0;
    if (notran) {
        //
        //        Solve  A*X = B.
        //
        //        Solve L*X = B, overwriting B with X.
        //
        //        L is represented as a product of permutations and unit lower
        //        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
        //        where each transformation L(i) is a rank-one modification of
        //        the identity matrix.
        //
        if (lnoti) {
            for (j = 1; j <= n - 1; j = j + 1) {
                lm = min(kl, n - j);
                l = ipiv[j - 1];
                if (l != j) {
                    Cswap(nrhs, &b[(l - 1)], ldb, &b[(j - 1)], ldb);
                }
                Cgeru(lm, nrhs, -one, ab[((kd + 1) - 1) + (j - 1) * ldab], 1, &b[(j - 1)], ldb, &b[((j + 1) - 1)], ldb);
            }
        }
        //
        for (i = 1; i <= nrhs; i = i + 1) {
            //
            //           Solve U*X = B, overwriting B with X.
            //
            Ctbsv("Upper", "No transpose", "Non-unit", n, kl + ku, ab, ldab, &b[(i - 1) * ldb], 1);
        }
        //
    } else if (Mlsame(trans, "T")) {
        //
        //        Solve A**T * X = B.
        //
        for (i = 1; i <= nrhs; i = i + 1) {
            //
            //           Solve U**T * X = B, overwriting B with X.
            //
            Ctbsv("Upper", "Transpose", "Non-unit", n, kl + ku, ab, ldab, &b[(i - 1) * ldb], 1);
        }
        //
        //        Solve L**T * X = B, overwriting B with X.
        //
        if (lnoti) {
            for (j = n - 1; j >= 1; j = j - 1) {
                lm = min(kl, n - j);
                Cgemv("Transpose", lm, nrhs, -one, &b[((j + 1) - 1)], ldb, ab[((kd + 1) - 1) + (j - 1) * ldab], 1, one, &b[(j - 1)], ldb);
                l = ipiv[j - 1];
                if (l != j) {
                    Cswap(nrhs, &b[(l - 1)], ldb, &b[(j - 1)], ldb);
                }
            }
        }
        //
    } else {
        //
        //        Solve A**H * X = B.
        //
        for (i = 1; i <= nrhs; i = i + 1) {
            //
            //           Solve U**H * X = B, overwriting B with X.
            //
            Ctbsv("Upper", "Conjugate transpose", "Non-unit", n, kl + ku, ab, ldab, &b[(i - 1) * ldb], 1);
        }
        //
        //        Solve L**H * X = B, overwriting B with X.
        //
        if (lnoti) {
            for (j = n - 1; j >= 1; j = j - 1) {
                lm = min(kl, n - j);
                Clacgv(nrhs, &b[(j - 1)], ldb);
                Cgemv("Conjugate transpose", lm, nrhs, -one, &b[((j + 1) - 1)], ldb, ab[((kd + 1) - 1) + (j - 1) * ldab], 1, one, &b[(j - 1)], ldb);
                Clacgv(nrhs, &b[(j - 1)], ldb);
                l = ipiv[j - 1];
                if (l != j) {
                    Cswap(nrhs, &b[(l - 1)], ldb, &b[(j - 1)], ldb);
                }
            }
        }
    }
    //
    //     End of Cgbtrs
    //
}
