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

void Chbgvd(const char *jobz, const char *uplo, INTEGER const n, INTEGER const ka, INTEGER const kb, COMPLEX *ab, INTEGER const ldab, COMPLEX *bb, INTEGER const ldbb, REAL *w, COMPLEX *z, INTEGER const ldz, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER const lrwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
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
        lwmin = 1 + n;
        lrwmin = 1 + n;
        liwmin = 1;
    } else if (wantz) {
        lwmin = 2 * pow2(n);
        lrwmin = 1 + 5 * n + 2 * pow2(n);
        liwmin = 3 + 5 * n;
    } else {
        lwmin = n;
        lrwmin = n;
        liwmin = 1;
    }
    if (!(wantz || Mlsame(jobz, "N"))) {
        info = -1;
    } else if (!(upper || Mlsame(uplo, "L"))) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (ka < 0) {
        info = -4;
    } else if (kb < 0 || kb > ka) {
        info = -5;
    } else if (ldab < ka + 1) {
        info = -7;
    } else if (ldbb < kb + 1) {
        info = -9;
    } else if (ldz < 1 || (wantz && ldz < n)) {
        info = -12;
    }
    //
    if (info == 0) {
        work[1 - 1] = lwmin;
        rwork[1 - 1] = lrwmin;
        iwork[1 - 1] = liwmin;
        //
        if (lwork < lwmin && !lquery) {
            info = -14;
        } else if (lrwork < lrwmin && !lquery) {
            info = -16;
        } else if (liwork < liwmin && !lquery) {
            info = -18;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Chbgvd", -info);
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
    //     Form a split Cholesky factorization of B.
    //
    Cpbstf(uplo, n, kb, bb, ldbb, info);
    if (info != 0) {
        info += n;
        return;
    }
    //
    //     Transform problem to standard eigenvalue problem.
    //
    INTEGER inde = 1;
    INTEGER indwrk = inde + n;
    INTEGER indwk2 = 1 + n * n;
    INTEGER llwk2 = lwork - indwk2 + 2;
    INTEGER llrwk = lrwork - indwrk + 2;
    INTEGER iinfo = 0;
    Chbgst(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, z, ldz, work, rwork, iinfo);
    //
    //     Reduce Hermitian band matrix to tridiagonal form.
    //
    char vect;
    if (wantz) {
        vect = 'U';
    } else {
        vect = 'N';
    }
    Chbtrd(&vect, uplo, n, ka, ab, ldab, w, &rwork[inde - 1], z, ldz, work, iinfo);
    //
    //     For eigenvalues only, call Rsterf.  For eigenvectors, call Cstedc.
    //
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    if (!wantz) {
        Rsterf(n, w, &rwork[inde - 1], info);
    } else {
        Cstedc("I", n, w, &rwork[inde - 1], work, n, &work[indwk2 - 1], llwk2, &rwork[indwrk - 1], llrwk, iwork, liwork, info);
        Cgemm("N", "N", n, n, n, cone, z, ldz, work, n, czero, &work[indwk2 - 1], n);
        Clacpy("A", n, n, &work[indwk2 - 1], n, z, ldz);
    }
    //
    work[1 - 1] = lwmin;
    rwork[1 - 1] = lrwmin;
    iwork[1 - 1] = liwmin;
    //
    //     End of Chbgvd
    //
}
