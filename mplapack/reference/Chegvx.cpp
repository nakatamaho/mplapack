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

void Chegvx(INTEGER const itype, const char *jobz, const char *range, const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL const vl, REAL const vu, INTEGER const il, INTEGER const iu, REAL const abstol, INTEGER &m, REAL *w, COMPLEX *z, INTEGER const ldz, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER *iwork, INTEGER *ifail, INTEGER &info) {
    //
    //     Test the input parameters.
    //
    bool wantz = Mlsame(jobz, "V");
    bool upper = Mlsame(uplo, "U");
    bool alleig = Mlsame(range, "A");
    bool valeig = Mlsame(range, "V");
    bool indeig = Mlsame(range, "I");
    bool lquery = (lwork == -1);
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
    } else if (lda < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -9;
    } else {
        if (valeig) {
            if (n > 0 && vu <= vl) {
                info = -11;
            }
        } else if (indeig) {
            if (il < 1 || il > max((INTEGER)1, n)) {
                info = -12;
            } else if (iu < min(n, il) || iu > n) {
                info = -13;
            }
        }
    }
    if (info == 0) {
        if (ldz < 1 || (wantz && ldz < n)) {
            info = -18;
        }
    }
    //
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    if (info == 0) {
        nb = iMlaenv(1, "Chetrd", uplo, n, -1, -1, -1);
        lwkopt = max((INTEGER)1, (nb + 1) * n);
        work[1 - 1] = lwkopt;
        //
        if (lwork < max((INTEGER)1, 2 * n) && !lquery) {
            info = -20;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Chegvx", -info);
        return;
    } else if (lquery) {
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
    Cpotrf(uplo, n, b, ldb, info);
    if (info != 0) {
        info += n;
        return;
    }
    //
    //     Transform problem to standard eigenvalue problem and solve.
    //
    Chegst(itype, uplo, n, a, lda, b, ldb, info);
    Cheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, rwork, iwork, ifail, info);
    //
    char trans;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
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
            //           backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y
            //
            if (upper) {
                trans = 'N';
            } else {
                trans = 'C';
            }
            //
            Ctrsm("Left", uplo, &trans, "Non-unit", n, m, cone, b, ldb, z, ldz);
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
            Ctrmm("Left", uplo, &trans, "Non-unit", n, m, cone, b, ldb, z, ldz);
        }
    }
    //
    //     Set WORK(1) to optimal complex workspace size.
    //
    work[1 - 1] = lwkopt;
    //
    //     End of Chegvx
    //
}
