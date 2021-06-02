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

void Rsyev_2stage(const char *jobz, const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, REAL *w, REAL *work, INTEGER const lwork, INTEGER &info) {
    //
    //     Test the input parameters.
    //
    bool wantz = Mlsame(jobz, "V");
    bool lower = Mlsame(uplo, "L");
    bool lquery = (lwork == -1);
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
    INTEGER kd = 0;
    INTEGER ib = 0;
    INTEGER lhtrd = 0;
    INTEGER lwtrd = 0;
    INTEGER lwmin = 0;
    if (info == 0) {
        kd = iMlaenv2stage(1, "Rsytrd_2stage", jobz, n, -1, -1, -1);
        ib = iMlaenv2stage(2, "Rsytrd_2stage", jobz, n, kd, -1, -1);
        lhtrd = iMlaenv2stage(3, "Rsytrd_2stage", jobz, n, kd, ib, -1);
        lwtrd = iMlaenv2stage(4, "Rsytrd_2stage", jobz, n, kd, ib, -1);
        lwmin = 2 * n + lhtrd + lwtrd;
        work[1 - 1] = lwmin;
        //
        if (lwork < lwmin && !lquery) {
            info = -8;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rsyev_2stage", -info);
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
    const REAL one = 1.0;
    if (n == 1) {
        w[1 - 1] = a[(1 - 1)];
        work[1 - 1] = 2;
        if (wantz) {
            a[(1 - 1)] = one;
        }
        return;
    }
    //
    //     Get machine constants.
    //
    REAL safmin = Rlamch("Safe minimum");
    REAL eps = Rlamch("Precision");
    REAL smlnum = safmin / eps;
    REAL bignum = one / smlnum;
    REAL rmin = sqrt(smlnum);
    REAL rmax = sqrt(bignum);
    //
    //     Scale matrix to allowable range, if necessary.
    //
    REAL anrm = Rlansy("M", uplo, n, a, lda, work);
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
        Rlascl(uplo, 0, 0, one, sigma, n, n, a, lda, info);
    }
    //
    //     Call Rsytrd_2stage to reduce symmetric matrix to tridiagonal form.
    //
    INTEGER inde = 1;
    INTEGER indtau = inde + n;
    INTEGER indhous = indtau + n;
    INTEGER indwrk = indhous + lhtrd;
    INTEGER llwork = lwork - indwrk + 1;
    //
    INTEGER iinfo = 0;
    Rsytrd_2stage(jobz, uplo, n, a, lda, w, &work[inde - 1], &work[indtau - 1], &work[indhous - 1], lhtrd, &work[indwrk - 1], llwork, iinfo);
    //
    //     For eigenvalues only, call Rsterf.  For eigenvectors, first call
    //     Rorgtr to generate the orthogonal matrix, then call Rsteqr.
    //
    if (!wantz) {
        Rsterf(n, w, &work[inde - 1], info);
    } else {
        //        Not available in this release, and argument checking should not
        //        let it getting here
        return;
        Rorgtr(uplo, n, a, lda, &work[indtau - 1], &work[indwrk - 1], llwork, iinfo);
        Rsteqr(jobz, n, w, &work[inde - 1], a, lda, &work[indtau - 1], info);
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
    //     Set WORK(1) to optimal workspace size.
    //
    work[1 - 1] = lwmin;
    //
    //     End of Rsyev_2stage
    //
}
