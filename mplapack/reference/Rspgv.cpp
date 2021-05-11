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

void Rspgv(INTEGER const itype, const char *jobz, const char *uplo, INTEGER const n, REAL *ap, REAL *bp, REAL *w, REAL *z, INTEGER const ldz, REAL *work, INTEGER &info) {
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
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    bool wantz = Mlsame(jobz, "V");
    bool upper = Mlsame(uplo, "U");
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
    if (info != 0) {
        Mxerbla("Rspgv ", -info);
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
    Rpptrf(uplo, n, bp, info);
    if (info != 0) {
        info += n;
        return;
    }
    //
    //     Transform problem to standard eigenvalue problem and solve.
    //
    Rspgst(itype, uplo, n, ap, bp, info);
    Rspev(jobz, uplo, n, ap, w, z, ldz, work, info);
    //
    INTEGER neig = 0;
    char trans;
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
            //           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
            //
            if (upper) {
                trans = 'N';
            } else {
                trans = 'T';
            }
            //
            for (j = 1; j <= neig; j = j + 1) {
                Rtpsv(uplo, &trans, "Non-unit", n, bp, &z[(j - 1) * ldz], 1);
            }
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
            for (j = 1; j <= neig; j = j + 1) {
                Rtpmv(uplo, &trans, "Non-unit", n, bp, &z[(j - 1) * ldz], 1);
            }
        }
    }
    //
    //     End of Rspgv
    //
}
