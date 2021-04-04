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

void Chpgvd(INTEGER const &itype, const char *jobz, const char *uplo, INTEGER const &n, COMPLEX *ap, COMPLEX *bp, REAL *w, COMPLEX *z, INTEGER const &ldz, COMPLEX *work, INTEGER const &lwork, REAL *rwork, INTEGER const &lrwork, arr_ref<INTEGER> iwork, INTEGER const &liwork, INTEGER &info) {
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
    if (itype < 1 || itype > 3) {
        info = -1;
    } else if (!(wantz || Mlsame(jobz, "N"))) {
        info = -2;
    } else if (!(upper || Mlsame(uplo, "L"))) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (ldz < 1 || (wantz && ldz < n)) {
        info = -9;
    }
    //
    INTEGER lwmin = 0;
    INTEGER liwmin = 0;
    INTEGER lrwmin = 0;
    if (info == 0) {
        if (n <= 1) {
            lwmin = 1;
            liwmin = 1;
            lrwmin = 1;
        } else {
            if (wantz) {
                lwmin = 2 * n;
                lrwmin = 1 + 5 * n + 2 * pow2(n);
                liwmin = 3 + 5 * n;
            } else {
                lwmin = n;
                lrwmin = n;
                liwmin = 1;
            }
        }
        //
        work[1 - 1] = lwmin;
        rwork[1 - 1] = lrwmin;
        iwork[1 - 1] = liwmin;
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
        Mxerbla("Chpgvd", -info);
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
    Cpptrf(uplo, n, bp, info);
    if (info != 0) {
        info += n;
        return;
    }
    //
    //     Transform problem to standard eigenvalue problem and solve.
    //
    Chpgst(itype, uplo, n, ap, bp, info);
    Chpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info);
    lwmin = max(lwmin.real(), work[1 - 1].real());
    lrwmin = max(lrwmin.real(), rwork[1 - 1].real());
    liwmin = max(liwmin.real(), iwork[1 - 1].real());
    //
    INTEGER neig = 0;
    str<1> trans = char0;
    INTEGER j = 0;
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
                trans = "N";
            } else {
                trans = "C";
            }
            //
            for (j = 1; j <= neig; j = j + 1) {
                Ctpsv(uplo, trans, "Non-unit", n, bp, z[(j - 1) * ldz], 1);
            }
            //
        } else if (itype == 3) {
            //
            //           For B*A*x=(lambda)*x;
            //           backtransform eigenvectors: x = L*y or U**H *y
            //
            if (upper) {
                trans = "C";
            } else {
                trans = "N";
            }
            //
            for (j = 1; j <= neig; j = j + 1) {
                Ctpmv(uplo, trans, "Non-unit", n, bp, z[(j - 1) * ldz], 1);
            }
        }
    }
    //
    work[1 - 1] = lwmin;
    rwork[1 - 1] = lrwmin;
    iwork[1 - 1] = liwmin;
    //
    //     End of Chpgvd
    //
}
