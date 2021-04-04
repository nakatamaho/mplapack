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

void Rsygvd(INTEGER const &itype, const char *jobz, const char *uplo, INTEGER const &n, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL *w, REAL *work, INTEGER const &lwork, arr_ref<INTEGER> iwork, INTEGER const &liwork, INTEGER &info) {
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
    bool lquery = (lwork == -1 || liwork == -1);
    //
    info = 0;
    INTEGER liwmin = 0;
    INTEGER lwmin = 0;
    if (n <= 1) {
        liwmin = 1;
        lwmin = 1;
    } else if (wantz) {
        liwmin = 3 + 5 * n;
        lwmin = 1 + 6 * n + 2 * pow2(n);
    } else {
        liwmin = 1;
        lwmin = 2 * n + 1;
    }
    INTEGER lopt = lwmin;
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
        iwork[1 - 1] = liopt;
        //
        if (lwork < lwmin && !lquery) {
            info = -11;
        } else if (liwork < liwmin && !lquery) {
            info = -13;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rsygvd", -info);
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
    Rsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
    lopt = max(lopt.real(), work[1 - 1].real());
    liopt = max(liopt.real(), iwork[1 - 1].real());
    //
    str<1> trans = char0;
    const REAL one = 1.0;
    if (wantz && info == 0) {
        //
        //        Backtransform eigenvectors to the original problem.
        //
        if (itype == 1 || itype == 2) {
            //
            //           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            //           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
            //
            if (upper) {
                trans = "N";
            } else {
                trans = "T";
            }
            //
            Rtrsm("Left", uplo, trans, "Non-unit", n, n, one, b, ldb, a, lda);
            //
        } else if (itype == 3) {
            //
            //           For B*A*x=(lambda)*x;
            //           backtransform eigenvectors: x = L*y or U**T*y
            //
            if (upper) {
                trans = "T";
            } else {
                trans = "N";
            }
            //
            Rtrmm("Left", uplo, trans, "Non-unit", n, n, one, b, ldb, a, lda);
        }
    }
    //
    work[1 - 1] = lopt;
    iwork[1 - 1] = liopt;
    //
    //     End of Rsygvd
    //
}
