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

void Rspgvx(INTEGER const &itype, const char *jobz, const char *range, const char *uplo, INTEGER const &n, REAL *ap, REAL *bp, REAL const &vl, REAL const &vu, INTEGER const &il, INTEGER const &iu, REAL const &abstol, INTEGER &m, REAL *w, REAL *z, INTEGER const &ldz, REAL *work, INTEGER *iwork, INTEGER *ifail, INTEGER &info) {
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
    // =====================================================================
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
    bool upper = Mlsame(uplo, "U");
    bool wantz = Mlsame(jobz, "V");
    bool alleig = Mlsame(range, "A");
    bool valeig = Mlsame(range, "V");
    bool indeig = Mlsame(range, "I");
    //
    info = 0;
    if (itype < 1 || itype > 3) {
        info = -1;
    } else if (!(wantz || Mlsame(jobz, "N"))) {
        info = -2;
    } else if (!(alleig || valeig || indeig)) {
        info = -3;
    } else if (!(upper || Mlsame(uplo, "L"))) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else {
        if (valeig) {
            if (n > 0 && vu <= vl) {
                info = -9;
            }
        } else if (indeig) {
            if (il < 1) {
                info = -10;
            } else if (iu < min(n, il) || iu > n) {
                info = -11;
            }
        }
    }
    if (info == 0) {
        if (ldz < 1 || (wantz && ldz < n)) {
            info = -16;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rspgvx", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    m = 0;
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
    Rspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info);
    //
    str<1> trans = char0;
    INTEGER j = 0;
    if (wantz) {
        //
        //        Backtransform eigenvectors to the original problem.
        //
        if (info > 0) {
            m = info - 1;
        }
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
            for (j = 1; j <= m; j = j + 1) {
                Rtpsv(uplo, trans, "Non-unit", n, bp, z[(j - 1) * ldz], 1);
            }
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
            for (j = 1; j <= m; j = j + 1) {
                Rtpmv(uplo, trans, "Non-unit", n, bp, z[(j - 1) * ldz], 1);
            }
        }
    }
    //
    //     End of Rspgvx
    //
}
