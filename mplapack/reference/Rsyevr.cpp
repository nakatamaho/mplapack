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

void Rsyevr(const char *jobz, const char *range, const char *uplo, INTEGER const &n, REAL *a, INTEGER const &lda, REAL const &vl, REAL const &vu, INTEGER const &il, INTEGER const &iu, REAL const &abstol, INTEGER &m, REAL *w, REAL *z, INTEGER const &ldz, arr_ref<INTEGER> isuppz, REAL *work, INTEGER const &lwork, arr_ref<INTEGER> iwork, INTEGER const &liwork, INTEGER &info) {
    INTEGER ieeeok = 0;
    bool lower = false;
    bool wantz = false;
    bool alleig = false;
    bool valeig = false;
    bool indeig = false;
    bool lquery = false;
    INTEGER lwmin = 0;
    INTEGER liwmin = 0;
    INTEGER nb = 0;
    INTEGER lwkopt = 0;
    const REAL one = 1.0;
    REAL safmin = 0.0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    REAL bignum = 0.0;
    REAL rmin = 0.0;
    REAL rmax = 0.0;
    INTEGER iscale = 0;
    REAL abstll = 0.0;
    REAL vll = 0.0;
    REAL vuu = 0.0;
    REAL anrm = 0.0;
    const REAL zero = 0.0;
    REAL sigma = 0.0;
    INTEGER j = 0;
    INTEGER indtau = 0;
    INTEGER indd = 0;
    INTEGER inde = 0;
    INTEGER inddd = 0;
    INTEGER indee = 0;
    INTEGER indwk = 0;
    INTEGER llwork = 0;
    INTEGER indibl = 0;
    INTEGER indisp = 0;
    INTEGER indifl = 0;
    INTEGER indiwo = 0;
    INTEGER iinfo = 0;
    const REAL two = 2.0e+0;
    bool tryrac = false;
    INTEGER indwkn = 0;
    INTEGER llwrkn = 0;
    str<1> order = char0;
    INTEGER nsplit = 0;
    INTEGER imax = 0;
    INTEGER i = 0;
    REAL tmp1 = 0.0;
    INTEGER jj = 0;
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
    ieeeok = iMlaenv[(10 - 1) + ("Rsyevr" - 1) * ldiMlaenv];
    //
    lower = Mlsame(uplo, "L");
    wantz = Mlsame(jobz, "V");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");
    //
    lquery = ((lwork == -1) || (liwork == -1));
    //
    lwmin = max((INTEGER)1, 26 * n);
    liwmin = max((INTEGER)1, 10 * n);
    //
    info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
        info = -1;
    } else if (!(alleig || valeig || indeig)) {
        info = -2;
    } else if (!(lower || Mlsame(uplo, "U"))) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, n)) {
        info = -6;
    } else {
        if (valeig) {
            if (n > 0 && vu <= vl) {
                info = -8;
            }
        } else if (indeig) {
            if (il < 1 || il > max((INTEGER)1, n)) {
                info = -9;
            } else if (iu < min(n, il) || iu > n) {
                info = -10;
            }
        }
    }
    if (info == 0) {
        if (ldz < 1 || (wantz && ldz < n)) {
            info = -15;
        } else if (lwork < lwmin && !lquery) {
            info = -18;
        } else if (liwork < liwmin && !lquery) {
            info = -20;
        }
    }
    //
    if (info == 0) {
        nb = iMlaenv[("Rsytrd" - 1) * ldiMlaenv];
        nb = max(nb, iMlaenv[("Rormtr" - 1) * ldiMlaenv]);
        lwkopt = max((nb + 1) * n, lwmin);
        work[1 - 1] = lwkopt;
        iwork[1 - 1] = liwmin;
    }
    //
    if (info != 0) {
        Mxerbla("Rsyevr", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    m = 0;
    if (n == 0) {
        work[1 - 1] = 1;
        return;
    }
    //
    if (n == 1) {
        work[1 - 1] = 7;
        if (alleig || indeig) {
            m = 1;
            w[1 - 1] = a[(1 - 1)];
        } else {
            if (vl < a[(1 - 1)] && vu >= a[(1 - 1)]) {
                m = 1;
                w[1 - 1] = a[(1 - 1)];
            }
        }
        if (wantz) {
            z[(1 - 1)] = one;
            isuppz[1 - 1] = 1;
            isuppz[2 - 1] = 1;
        }
        return;
    }
    //
    //     Get machine constants.
    //
    safmin = dlamch("Safe minimum");
    eps = dlamch("Precision");
    smlnum = safmin / eps;
    bignum = one / smlnum;
    rmin = sqrt(smlnum);
    rmax = min(sqrt(bignum), one / sqrt(sqrt(safmin)));
    //
    //     Scale matrix to allowable range, if necessary.
    //
    iscale = 0;
    abstll = abstol;
    if (valeig) {
        vll = vl;
        vuu = vu;
    }
    anrm = Rlansy[("M" - 1) + (uplo - 1) * ldRlansy];
    if (anrm > zero && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        if (lower) {
            for (j = 1; j <= n; j = j + 1) {
                Rscal(n - j + 1, sigma, a[(j - 1) + (j - 1) * lda], 1);
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                Rscal(j, sigma, a[(j - 1) * lda], 1);
            }
        }
        if (abstol > 0) {
            abstll = abstol * sigma;
        }
        if (valeig) {
            vll = vl * sigma;
            vuu = vu * sigma;
        }
    }
    //
    //     Initialize indices INTEGERo workspaces.  Note: The IWORK indices are
    //     used only if Rsterf or Rstemr fail.
    //
    //     WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the
    //     elementary reflectors used in Rsytrd.
    indtau = 1;
    //     WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries.
    indd = indtau + n;
    //     WORK(INDE:INDE+N-1) stores the off-diagonal entries of the
    //     tridiagonal matrix from Rsytrd.
    inde = indd + n;
    //     WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over
    //     -written by Rstemr (the Rsterf path copies the diagonal to W).
    inddd = inde + n;
    //     WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over
    //     -written while computing the eigenvalues in Rsterf and Rstemr.
    indee = inddd + n;
    //     INDWK is the starting offset of the left-over workspace, and
    //     LLWORK is the remaining workspace size.
    indwk = indee + n;
    llwork = lwork - indwk + 1;
    //
    //     IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in Rstebz and
    //     stores the block indices of each of the M<=N eigenvalues.
    indibl = 1;
    //     IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in Rstebz and
    //     stores the starting and finishing indices of each block.
    indisp = indibl + n;
    //     IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
    //     that corresponding to eigenvectors that fail to converge in
    //     Rstein.  This information is discarded; if any fail, the driver
    //     returns INFO > 0.
    indifl = indisp + n;
    //     INDIWO is the offset of the remaining INTEGEReger workspace.
    indiwo = indifl + n;
    //
    //     Call Rsytrd to reduce symmetric matrix to tridiagonal form.
    //
    Rsytrd(uplo, n, a, lda, work[indd - 1], work[inde - 1], work[indtau - 1], work[indwk - 1], llwork, iinfo);
    //
    //     If all eigenvalues are desired
    //     then call Rsterf or Rstemr and Rormtr.
    //
    if ((alleig || (indeig && il == 1 && iu == n)) && ieeeok == 1) {
        if (!wantz) {
            Rcopy(n, work[indd - 1], 1, w, 1);
            Rcopy(n - 1, work[inde - 1], 1, work[indee - 1], 1);
            Rsterf(n, w, work[indee - 1], info);
        } else {
            Rcopy(n - 1, work[inde - 1], 1, work[indee - 1], 1);
            Rcopy(n, work[indd - 1], 1, work[inddd - 1], 1);
            //
            if (abstol <= two * n * eps) {
                tryrac = true;
            } else {
                tryrac = false;
            }
            Rstemr(jobz, "A", n, work[inddd - 1], work[indee - 1], vl, vu, il, iu, m, w, z, ldz, n, isuppz, tryrac, work[indwk - 1], lwork, iwork, liwork, info);
            //
            //        Apply orthogonal matrix used in reduction to tridiagonal
            //        form to eigenvectors returned by Rstemr.
            //
            if (wantz && info == 0) {
                indwkn = inde;
                llwrkn = lwork - indwkn + 1;
                Rormtr("L", uplo, "N", n, m, a, lda, work[indtau - 1], z, ldz, work[indwkn - 1], llwrkn, iinfo);
            }
        }
        //
        if (info == 0) {
            //           Everything worked.  Skip Rstebz/Rstein.  IWORK(:) are
            //           undefined.
            m = n;
            goto statement_30;
        }
        info = 0;
    }
    //
    //     Otherwise, call Rstebz and, if eigenvectors are desired, Rstein.
    //     Also call Rstebz and Rstein if Rstemr fails.
    //
    if (wantz) {
        order = "B";
    } else {
        order = "E";
    }
    //
    Rstebz(range, order, n, vll, vuu, il, iu, abstll, work[indd - 1], work[inde - 1], m, nsplit, w, iwork[indibl - 1], iwork[indisp - 1], work[indwk - 1], iwork[indiwo - 1], info);
    //
    if (wantz) {
        Rstein(n, work[indd - 1], work[inde - 1], m, w, iwork[indibl - 1], iwork[indisp - 1], z, ldz, work[indwk - 1], iwork[indiwo - 1], iwork[indifl - 1], info);
        //
        //        Apply orthogonal matrix used in reduction to tridiagonal
        //        form to eigenvectors returned by Rstein.
        //
        indwkn = inde;
        llwrkn = lwork - indwkn + 1;
        Rormtr("L", uplo, "N", n, m, a, lda, work[indtau - 1], z, ldz, work[indwkn - 1], llwrkn, iinfo);
    }
//
//     If matrix was scaled, then rescale eigenvalues appropriately.
//
//  Jump here if Rstemr/Rstein succeeded.
statement_30:
    if (iscale == 1) {
        if (info == 0) {
            imax = m;
        } else {
            imax = info - 1;
        }
        Rscal(imax, one / sigma, w, 1);
    }
    //
    //     If eigenvalues are not in order, then sort them, along with
    //     eigenvectors.  Note: We do not sort the IFAIL portion of IWORK.
    //     It may not be initialized (if Rstemr/Rstein succeeded), and we do
    //     not return this detailed information to the user.
    //
    if (wantz) {
        for (j = 1; j <= m - 1; j = j + 1) {
            i = 0;
            tmp1 = w[j - 1];
            for (jj = j + 1; jj <= m; jj = jj + 1) {
                if (w[jj - 1] < tmp1) {
                    i = jj;
                    tmp1 = w[jj - 1];
                }
            }
            //
            if (i != 0) {
                w[i - 1] = w[j - 1];
                w[j - 1] = tmp1;
                Rswap(n, z[(i - 1) * ldz], 1, z[(j - 1) * ldz], 1);
            }
        }
    }
    //
    //     Set WORK(1) to optimal workspace size.
    //
    work[1 - 1] = lwkopt;
    iwork[1 - 1] = liwmin;
    //
    //     End of Rsyevr
    //
}
