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

void Cheevd_2stage(const char *jobz, const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL *w, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER const lrwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
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
    //
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
    bool lower = Mlsame(uplo, "L");
    bool lquery = (lwork == -1 || lrwork == -1 || liwork == -1);
    //
    info = 0;
    if (!(Mlsame(jobz, "N"))) {
        info = -1;
    } else if (!(lower || Mlsame(uplo, "U"))) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    }
    //
    INTEGER lwmin = 0;
    INTEGER lrwmin = 0;
    INTEGER liwmin = 0;
    INTEGER kd = 0;
    INTEGER ib = 0;
    INTEGER lhtrd = 0;
    INTEGER lwtrd = 0;
    if (info == 0) {
        if (n <= 1) {
            lwmin = 1;
            lrwmin = 1;
            liwmin = 1;
        } else {
            kd = iMlaenv2stage(1, "Chetrd_2stage", jobz, n, -1, -1, -1);
            ib = iMlaenv2stage(2, "Chetrd_2stage", jobz, n, kd, -1, -1);
            lhtrd = iMlaenv2stage(3, "Chetrd_2stage", jobz, n, kd, ib, -1);
            lwtrd = iMlaenv2stage(4, "Chetrd_2stage", jobz, n, kd, ib, -1);
            if (wantz) {
                lwmin = 2 * n + n * n;
                lrwmin = 1 + 5 * n + 2 * pow2(n);
                liwmin = 3 + 5 * n;
            } else {
                lwmin = n + 1 + lhtrd + lwtrd;
                lrwmin = n;
                liwmin = 1;
            }
        }
        work[1 - 1] = lwmin;
        rwork[1 - 1] = lrwmin;
        iwork[1 - 1] = liwmin;
        //
        if (lwork < lwmin && !lquery) {
            info = -8;
        } else if (lrwork < lrwmin && !lquery) {
            info = -10;
        } else if (liwork < liwmin && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cheevd_2stage", -info);
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
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    if (n == 1) {
        w[1 - 1] = a[(1 - 1)].real();
        if (wantz) {
            a[(1 - 1)] = cone;
        }
        return;
    }
    //
    //     Get machine constants.
    //
    REAL safmin = Rlamch("Safe minimum");
    REAL eps = Rlamch("Precision");
    REAL smlnum = safmin / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    REAL rmin = sqrt(smlnum);
    REAL rmax = sqrt(bignum);
    //
    //     Scale matrix to allowable range, if necessary.
    //
    REAL anrm = Clanhe("M", uplo, n, a, lda, rwork);
    INTEGER iscale = 0;
    const REAL zero = 0.0;
    REAL sigma = 0.0;
    if (anrm > zero && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        Clascl(uplo, 0, 0, one, sigma, n, n, a, lda, info);
    }
    //
    //     Call Chetrd_2stage to reduce Hermitian matrix to tridiagonal form.
    //
    INTEGER inde = 1;
    INTEGER indrwk = inde + n;
    INTEGER llrwk = lrwork - indrwk + 1;
    INTEGER indtau = 1;
    INTEGER indhous = indtau + n;
    INTEGER indwrk = indhous + lhtrd;
    INTEGER llwork = lwork - indwrk + 1;
    INTEGER indwk2 = indwrk + n * n;
    INTEGER llwrk2 = lwork - indwk2 + 1;
    //
    INTEGER iinfo = 0;
    Chetrd_2stage(jobz, uplo, n, a, lda, w, &rwork[inde - 1], &work[indtau - 1], &work[indhous - 1], lhtrd, &work[indwrk - 1], llwork, iinfo);
    //
    //     For eigenvalues only, call Rsterf.  For eigenvectors, first call
    //     Cstedc to generate the eigenvector matrix, WORK(INDWRK), of the
    //     tridiagonal matrix, then call Cunmtr to multiply it to the
    //     Householder transformations represented as Householder vectors in
    //     A.
    //
    if (!wantz) {
        Rsterf(n, w, &rwork[inde - 1], info);
    } else {
        Cstedc("I", n, w, &rwork[inde - 1], &work[indwrk - 1], n, &work[indwk2 - 1], llwrk2, &rwork[indrwk - 1], llrwk, iwork, liwork, info);
        Cunmtr("L", uplo, "N", n, n, a, lda, &work[indtau - 1], &work[indwrk - 1], n, &work[indwk2 - 1], llwrk2, iinfo);
        Clacpy("A", n, n, &work[indwrk - 1], n, a, lda);
    }
    //
    //     If matrix was scaled, then rescale eigenvalues appropriately.
    //
    INTEGER imax = 0;
    if (iscale == 1) {
        if (info == 0) {
            imax = n;
        } else {
            imax = info - 1;
        }
        Rscal(imax, one / sigma, w, 1);
    }
    //
    work[1 - 1] = lwmin;
    rwork[1 - 1] = lrwmin;
    iwork[1 - 1] = liwmin;
    //
    //     End of Cheevd_2stage
    //
}
