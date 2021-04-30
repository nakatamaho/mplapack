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

void Chegvd(INTEGER const itype, const char *jobz, const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL *w, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER const lrwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
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
    bool lquery = (lwork == -1 || lrwork == -1 || liwork == -1);
    //
    info = 0;
    INTEGER lwmin = 0;
    INTEGER lrwmin = 0;
    INTEGER liwmin = 0;
    if (n <= 1) {
        lwmin = 1;
        lrwmin = 1;
        liwmin = 1;
    } else if (wantz) {
        lwmin = 2 * n + n * n;
        lrwmin = 1 + 5 * n + 2 * n * n;
        liwmin = 3 + 5 * n;
    } else {
        lwmin = n + 1;
        lrwmin = n;
        liwmin = 1;
    }
    INTEGER lopt = lwmin;
    INTEGER lropt = lrwmin;
    INTEGER liopt = liwmin;
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
    if (info == 0) {
        work[1 - 1] = lopt;
        rwork[1 - 1] = lropt;
        iwork[1 - 1] = liopt;
        //
        if (lwork < lwmin && !lquery) {
            info = -11;
        } else if (lrwork < lrwmin && !lquery) {
            info = -13;
        } else if (liwork < liwmin && !lquery) {
            info = -15;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Chegvd", -info);
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
    Cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info);
    lopt = max(lopt.real(), &work[1 - 1].real());
    lropt = max(lropt.real(), &rwork[1 - 1].real());
    liopt = max(liopt.real(), &iwork[1 - 1].real());
    //
    char trans;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    if (wantz && info == 0) {
        //
        //        Backtransform eigenvectors to the original problem.
        //
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
            Ctrsm("Left", uplo, trans, "Non-unit", n, n, cone, b, ldb, a, lda);
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
            Ctrmm("Left", uplo, trans, "Non-unit", n, n, cone, b, ldb, a, lda);
        }
    }
    //
    work[1 - 1] = lopt;
    rwork[1 - 1] = lropt;
    iwork[1 - 1] = liopt;
    //
    //     End of Chegvd
    //
}
