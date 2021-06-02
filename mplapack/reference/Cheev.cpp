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

void Cheev(const char *jobz, const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL *w, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER &info) {
    //
    //     Test the input parameters.
    //
    bool wantz = Mlsame(jobz, "V");
    bool lower = Mlsame(uplo, "L");
    bool lquery = (lwork == -1);
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
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    if (info == 0) {
        nb = iMlaenv(1, "Chetrd", uplo, n, -1, -1, -1);
        lwkopt = max((INTEGER)1, (nb + 1) * n);
        work[1 - 1] = lwkopt;
        //
        if (lwork < max((INTEGER)1, 2 * n - 1) && !lquery) {
            info = -8;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cheev", -info);
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
        work[1 - 1] = 1;
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
    //     Call Chetrd to reduce Hermitian matrix to tridiagonal form.
    //
    INTEGER inde = 1;
    INTEGER indtau = 1;
    INTEGER indwrk = indtau + n;
    INTEGER llwork = lwork - indwrk + 1;
    INTEGER iinfo = 0;
    Chetrd(uplo, n, a, lda, w, &rwork[inde - 1], &work[indtau - 1], &work[indwrk - 1], llwork, iinfo);
    //
    //     For eigenvalues only, call Rsterf.  For eigenvectors, first call
    //     Cungtr to generate the unitary matrix, then call Csteqr.
    //
    if (!wantz) {
        Rsterf(n, w, &rwork[inde - 1], info);
    } else {
        Cungtr(uplo, n, a, lda, &work[indtau - 1], &work[indwrk - 1], llwork, iinfo);
        indwrk = inde + n;
        Csteqr(jobz, n, w, &rwork[inde - 1], a, lda, &rwork[indwrk - 1], info);
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
    //     Set WORK(1) to optimal complex workspace size.
    //
    work[1 - 1] = lwkopt;
    //
    //     End of Cheev
    //
}
