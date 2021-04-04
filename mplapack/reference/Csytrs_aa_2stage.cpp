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

void Csytrs_aa_2stage(const char *uplo, INTEGER const &n, INTEGER const &nrhs, COMPLEX *a, INTEGER const &lda, COMPLEX *tb, INTEGER const &ltb, INTEGER *ipiv, INTEGER *ipiv2, COMPLEX *b, INTEGER const &ldb, INTEGER &info) {
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
    info = 0;
    bool upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ltb < (4 * n)) {
        info = -7;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -11;
    }
    if (info != 0) {
        Mxerbla("Csytrs_aa_2stage", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0 || nrhs == 0) {
        return;
    }
    //
    //     Read NB and compute LDTB
    //
    INTEGER nb = INTEGER(tb[1 - 1]);
    INTEGER ldtb = ltb / n;
    //
    const COMPLEX one = (1.0f, 0.0f);
    if (upper) {
        //
        //        Solve A*X = B, where A = U**T*T*U.
        //
        if (n > nb) {
            //
            //           Pivot, P**T * B -> B
            //
            Claswp(nrhs, b, ldb, nb + 1, n, ipiv, 1);
            //
            //           Compute (U**T \ B) -> B    [ (U**T \P**T * B) ]
            //
            Ctrsm("L", "U", "T", "U", n - nb, nrhs, one, a[((nb + 1) - 1) * lda], lda, b[((nb + 1) - 1)], ldb);
            //
        }
        //
        //        Compute T \ B -> B   [ T \ (U**T \P**T * B) ]
        //
        Cgbtrs("N", n, nb, nb, nrhs, tb, ldtb, ipiv2, b, ldb, info);
        if (n > nb) {
            //
            //           Compute (U \ B) -> B   [ U \ (T \ (U**T \P**T * B) ) ]
            //
            Ctrsm("L", "U", "N", "U", n - nb, nrhs, one, a[((nb + 1) - 1) * lda], lda, b[((nb + 1) - 1)], ldb);
            //
            //           Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]
            //
            Claswp(nrhs, b, ldb, nb + 1, n, ipiv, -1);
            //
        }
        //
    } else {
        //
        //        Solve A*X = B, where A = L*T*L**T.
        //
        if (n > nb) {
            //
            //           Pivot, P**T * B -> B
            //
            Claswp(nrhs, b, ldb, nb + 1, n, ipiv, 1);
            //
            //           Compute (L \ B) -> B    [ (L \P**T * B) ]
            //
            Ctrsm("L", "L", "N", "U", n - nb, nrhs, one, a[((nb + 1) - 1)], lda, b[((nb + 1) - 1)], ldb);
            //
        }
        //
        //        Compute T \ B -> B   [ T \ (L \P**T * B) ]
        //
        Cgbtrs("N", n, nb, nb, nrhs, tb, ldtb, ipiv2, b, ldb, info);
        if (n > nb) {
            //
            //           Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]
            //
            Ctrsm("L", "L", "T", "U", n - nb, nrhs, one, a[((nb + 1) - 1)], lda, b[((nb + 1) - 1)], ldb);
            //
            //           Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]
            //
            Claswp(nrhs, b, ldb, nb + 1, n, ipiv, -1);
            //
        }
    }
    //
    //     End of Csytrs_aa_2stage
    //
}
