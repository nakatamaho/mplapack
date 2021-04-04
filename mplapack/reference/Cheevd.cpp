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

void Cheevd(const char *jobz, const char *uplo, INTEGER const &n, COMPLEX *a, INTEGER const &lda, REAL *w, COMPLEX *work, INTEGER const &lwork, REAL *rwork, INTEGER const &lrwork, arr_ref<INTEGER> iwork, INTEGER const &liwork, INTEGER &info) {
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
    bool lower = Mlsame(uplo, "L");
    bool lquery = (lwork == -1 || lrwork == -1 || liwork == -1);
    //
    info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
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
    INTEGER lopt = 0;
    INTEGER lropt = 0;
    INTEGER liopt = 0;
    if (info == 0) {
        if (n <= 1) {
            lwmin = 1;
            lrwmin = 1;
            liwmin = 1;
            lopt = lwmin;
            lropt = lrwmin;
            liopt = liwmin;
        } else {
            if (wantz) {
                lwmin = 2 * n + n * n;
                lrwmin = 1 + 5 * n + 2 * pow2(n);
                liwmin = 3 + 5 * n;
            } else {
                lwmin = n + 1;
                lrwmin = n;
                liwmin = 1;
            }
            lopt = max(lwmin, n + iMlaenv[("Chetrd" - 1) * ldiMlaenv]);
            lropt = lrwmin;
            liopt = liwmin;
        }
        work[1 - 1] = lopt;
        rwork[1 - 1] = lropt;
        iwork[1 - 1] = liopt;
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
        Mxerbla("Cheevd", -info);
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
    const COMPLEX cone = (1.0, 0.0);
    if (n == 1) {
        w[1 - 1] = a[(1 - 1)];
        if (wantz) {
            a[(1 - 1)] = cone;
        }
        return;
    }
    //
    //     Get machine constants.
    //
    REAL safmin = dlamch("Safe minimum");
    REAL eps = dlamch("Precision");
    REAL smlnum = safmin / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    REAL rmin = sqrt(smlnum);
    REAL rmax = sqrt(bignum);
    //
    //     Scale matrix to allowable range, if necessary.
    //
    REAL anrm = Clanhe[("M" - 1) + (uplo - 1) * ldClanhe];
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
    //     Call Chetrd to reduce Hermitian matrix to tridiagonal form.
    //
    INTEGER inde = 1;
    INTEGER indtau = 1;
    INTEGER indwrk = indtau + n;
    INTEGER indrwk = inde + n;
    INTEGER indwk2 = indwrk + n * n;
    INTEGER llwork = lwork - indwrk + 1;
    INTEGER llwrk2 = lwork - indwk2 + 1;
    INTEGER llrwk = lrwork - indrwk + 1;
    INTEGER iinfo = 0;
    Chetrd(uplo, n, a, lda, w, rwork[inde - 1], work[indtau - 1], work[indwrk - 1], llwork, iinfo);
    //
    //     For eigenvalues only, call Rsterf.  For eigenvectors, first call
    //     Cstedc to generate the eigenvector matrix, WORK(INDWRK), of the
    //     tridiagonal matrix, then call Cunmtr to multiply it to the
    //     Householder transformations represented as Householder vectors in
    //     A.
    //
    if (!wantz) {
        Rsterf(n, w, rwork[inde - 1], info);
    } else {
        Cstedc("I", n, w, rwork[inde - 1], work[indwrk - 1], n, work[indwk2 - 1], llwrk2, rwork[indrwk - 1], llrwk, iwork, liwork, info);
        Cunmtr("L", uplo, "N", n, n, a, lda, work[indtau - 1], work[indwrk - 1], n, work[indwk2 - 1], llwrk2, iinfo);
        Clacpy("A", n, n, work[indwrk - 1], n, a, lda);
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
    work[1 - 1] = lopt;
    rwork[1 - 1] = lropt;
    iwork[1 - 1] = liopt;
    //
    //     End of Cheevd
    //
}
