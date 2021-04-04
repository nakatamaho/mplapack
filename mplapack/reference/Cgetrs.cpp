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

void Cgetrs(const char *trans, INTEGER const &n, INTEGER const &nrhs, COMPLEX *a, INTEGER const &lda, INTEGER *ipiv, COMPLEX *b, INTEGER const &ldb, INTEGER &info) {
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
    } else if (nrhs < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Cgetrs", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0 || nrhs == 0) {
        return;
    }
    //
    const COMPLEX one = (1.0, 0.0);
    if (notran) {
        //
        //        Solve A * X = B.
        //
        //        Apply row INTEGERerchanges to the right hand sides.
        //
        Claswp(nrhs, b, ldb, 1, n, ipiv, 1);
        //
        //        Solve L*X = B, overwriting B with X.
        //
        Ctrsm("Left", "Lower", "No transpose", "Unit", n, nrhs, one, a, lda, b, ldb);
        //
        //        Solve U*X = B, overwriting B with X.
        //
        Ctrsm("Left", "Upper", "No transpose", "Non-unit", n, nrhs, one, a, lda, b, ldb);
    } else {
        //
        //        Solve A**T * X = B  or A**H * X = B.
        //
        //        Solve U**T *X = B or U**H *X = B, overwriting B with X.
        //
        Ctrsm("Left", "Upper", trans, "Non-unit", n, nrhs, one, a, lda, b, ldb);
        //
        //        Solve L**T *X = B, or L**H *X = B overwriting B with X.
        //
        Ctrsm("Left", "Lower", trans, "Unit", n, nrhs, one, a, lda, b, ldb);
        //
        //        Apply row INTEGERerchanges to the solution vectors.
        //
        Claswp(nrhs, b, ldb, 1, n, ipiv, -1);
    }
    //
    //     End of Cgetrs
    //
}
