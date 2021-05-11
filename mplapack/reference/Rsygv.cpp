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

void Rsygv(INTEGER const itype, const char *jobz, const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *w, REAL *work, INTEGER const lwork, INTEGER &info) {
    //
    //  -- LAPACK driver routine --
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
    bool wantz = Mlsame(jobz, "V");
    bool upper = Mlsame(uplo, "U");
    bool lquery = (lwork == -1);
    //
    info = 0;
    if (itype < 1 || itype > 3) {
        info = -1;
    } else if (!(wantz || Mlsame(jobz, "N"))) {
        info = -2;
    } else if (!(upper || Mlsame(uplo, "L"))) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -8;
    }
    //
    INTEGER lwkmin = 0;
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    if (info == 0) {
        lwkmin = max((INTEGER)1, 3 * n - 1);
        nb = iMlaenv(1, "Rsytrd", uplo, n, -1, -1, -1);
        lwkopt = max(lwkmin, (nb + 2) * n);
        work[1 - 1] = lwkopt;
        //
        if (lwork < lwkmin && !lquery) {
            info = -11;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rsygv ", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Form a Cholesky factorization of B.
    //
    Rpotrf(uplo, n, b, ldb, info);
    if (info != 0) {
        info += n;
        return;
    }
    //
    //     Transform problem to standard eigenvalue problem and solve.
    //
    Rsygst(itype, uplo, n, a, lda, b, ldb, info);
    Rsyev(jobz, uplo, n, a, lda, w, work, lwork, info);
    //
    INTEGER neig = 0;
    char trans;
    const REAL one = 1.0;
    if (wantz) {
        //
        //        Backtransform eigenvectors to the original problem.
        //
        neig = n;
        if (info > 0) {
            neig = info - 1;
        }
        if (itype == 1 || itype == 2) {
            //
            //           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            //           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
            //
            if (upper) {
                trans = 'N';
            } else {
                trans = 'T';
            }
            //
            Rtrsm("Left", uplo, &trans, "Non-unit", n, neig, one, b, ldb, a, lda);
            //
        } else if (itype == 3) {
            //
            //           For B*A*x=(lambda)*x;
            //           backtransform eigenvectors: x = L*y or U**T*y
            //
            if (upper) {
                trans = 'T';
            } else {
                trans = 'N';
            }
            //
            Rtrmm("Left", uplo, &trans, "Non-unit", n, neig, one, b, ldb, a, lda);
        }
    }
    //
    work[1 - 1] = lwkopt;
    //
    //     End of Rsygv
    //
}
