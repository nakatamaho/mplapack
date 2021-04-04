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

void Rsyevd(const char *jobz, const char *uplo, INTEGER const &n, REAL *a, INTEGER const &lda, REAL *w, REAL *work, INTEGER const &lwork, arr_ref<INTEGER> iwork, INTEGER const &liwork, INTEGER &info) {
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
    bool lquery = (lwork == -1 || liwork == -1);
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
    INTEGER liwmin = 0;
    INTEGER lwmin = 0;
    INTEGER lopt = 0;
    INTEGER liopt = 0;
    if (info == 0) {
        if (n <= 1) {
            liwmin = 1;
            lwmin = 1;
            lopt = lwmin;
            liopt = liwmin;
        } else {
            if (wantz) {
                liwmin = 3 + 5 * n;
                lwmin = 1 + 6 * n + 2 * pow2(n);
            } else {
                liwmin = 1;
                lwmin = 2 * n + 1;
            }
            lopt = max(lwmin, 2 * n + iMlaenv[("Rsytrd" - 1) * ldiMlaenv]);
            liopt = liwmin;
        }
        work[1 - 1] = lopt;
        iwork[1 - 1] = liopt;
        //
        if (lwork < lwmin && !lquery) {
            info = -8;
        } else if (liwork < liwmin && !lquery) {
            info = -10;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rsyevd", -info);
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
        if (wantz) {
            a[(1 - 1)] = one;
        }
        return;
    }
    //
    //     Get machine constants.
    //
    REAL safmin = dlamch("Safe minimum");
    REAL eps = dlamch("Precision");
    REAL smlnum = safmin / eps;
    REAL bignum = one / smlnum;
    REAL rmin = sqrt(smlnum);
    REAL rmax = sqrt(bignum);
    //
    //     Scale matrix to allowable range, if necessary.
    //
    REAL anrm = Rlansy[("M" - 1) + (uplo - 1) * ldRlansy];
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
    //     Call Rsytrd to reduce symmetric matrix to tridiagonal form.
    //
    INTEGER inde = 1;
    INTEGER indtau = inde + n;
    INTEGER indwrk = indtau + n;
    INTEGER llwork = lwork - indwrk + 1;
    INTEGER indwk2 = indwrk + n * n;
    INTEGER llwrk2 = lwork - indwk2 + 1;
    //
    INTEGER iinfo = 0;
    Rsytrd(uplo, n, a, lda, w, work[inde - 1], work[indtau - 1], work[indwrk - 1], llwork, iinfo);
    //
    //     For eigenvalues only, call Rsterf.  For eigenvectors, first call
    //     Rstedc to generate the eigenvector matrix, WORK(INDWRK), of the
    //     tridiagonal matrix, then call Rormtr to multiply it by the
    //     Householder transformations stored in A.
    //
    if (!wantz) {
        Rsterf(n, w, work[inde - 1], info);
    } else {
        Rstedc("I", n, w, work[inde - 1], work[indwrk - 1], n, work[indwk2 - 1], llwrk2, iwork, liwork, info);
        Rormtr("L", uplo, "N", n, n, a, lda, work[indtau - 1], work[indwrk - 1], n, work[indwk2 - 1], llwrk2, iinfo);
        Rlacpy("A", n, n, work[indwrk - 1], n, a, lda);
    }
    //
    //     If matrix was scaled, then rescale eigenvalues appropriately.
    //
    if (iscale == 1) {
        Rscal(n, one / sigma, w, 1);
    }
    //
    work[1 - 1] = lopt;
    iwork[1 - 1] = liopt;
    //
    //     End of Rsyevd
    //
}
