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

void Chpev(const char *jobz, const char *uplo, INTEGER const n, COMPLEX *ap, REAL *w, COMPLEX *z, INTEGER const ldz, COMPLEX *work, REAL *rwork, INTEGER &info) {
    //
    //     Test the input parameters.
    //
    bool wantz = Mlsame(jobz, "V");
    //
    info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
        info = -1;
    } else if (!(Mlsame(uplo, "L") || Mlsame(uplo, "U"))) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (ldz < 1 || (wantz && ldz < n)) {
        info = -7;
    }
    //
    if (info != 0) {
        Mxerbla("Chpev", -info);
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
        w[1 - 1] = ap[1 - 1].real();
        rwork[1 - 1] = 1;
        if (wantz) {
            z[(1 - 1)] = one;
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
    REAL anrm = Clanhp("M", uplo, n, ap, rwork);
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
        CRscal((n * (n + 1)) / 2, sigma, ap, 1);
    }
    //
    //     Call Chptrd to reduce Hermitian packed matrix to tridiagonal form.
    //
    INTEGER inde = 1;
    INTEGER indtau = 1;
    INTEGER iinfo = 0;
    Chptrd(uplo, n, ap, w, &rwork[inde - 1], &work[indtau - 1], iinfo);
    //
    //     For eigenvalues only, call Rsterf.  For eigenvectors, first call
    //     Cupgtr to generate the orthogonal matrix, then call Csteqr.
    //
    INTEGER indwrk = 0;
    INTEGER indrwk = 0;
    if (!wantz) {
        Rsterf(n, w, &rwork[inde - 1], info);
    } else {
        indwrk = indtau + n;
        Cupgtr(uplo, n, ap, &work[indtau - 1], z, ldz, &work[indwrk - 1], iinfo);
        indrwk = inde + n;
        Csteqr(jobz, n, w, &rwork[inde - 1], z, ldz, &rwork[indrwk - 1], info);
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
    //     End of Chpev
    //
}
