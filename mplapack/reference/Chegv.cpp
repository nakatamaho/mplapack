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

void Chegv(INTEGER const itype, const char *jobz, const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL *w, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER &info) {
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
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    if (info == 0) {
        nb = iMlaenv(1, "Chetrd", uplo, n, -1, -1, -1);
        lwkopt = max((INTEGER)1, (nb + 1) * n);
        work[1 - 1] = lwkopt;
        //
        if (lwork < max((INTEGER)1, 2 * n - 1) && !lquery) {
            info = -11;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Chegv ", -info);
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
    Cpotrf(uplo, n, b, ldb, info);
    if (info != 0) {
        info += n;
        return;
    }
    //
    //     Transform problem to standard eigenvalue problem and solve.
    //
    Chegst(itype, uplo, n, a, lda, b, ldb, info);
    Cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
    //
    INTEGER neig = 0;
    char trans;
    const COMPLEX one = COMPLEX(1.0, 0.0);
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
            //           backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y
            //
            if (upper) {
                trans = 'N';
            } else {
                trans = 'C';
            }
            //
            Ctrsm("Left", uplo, &trans, "Non-unit", n, neig, one, b, ldb, a, lda);
            //
        } else if (itype == 3) {
            //
            //           For B*A*x=(lambda)*x;
            //           backtransform eigenvectors: x = L*y or U**H *y
            //
            if (upper) {
                trans = 'C';
            } else {
                trans = 'N';
            }
            //
            Ctrmm("Left", uplo, &trans, "Non-unit", n, neig, one, b, ldb, a, lda);
        }
    }
    //
    work[1 - 1] = lwkopt;
    //
    //     End of Chegv
    //
}
