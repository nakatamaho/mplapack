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

void Rstevr(const char *jobz, const char *range, INTEGER const &n, REAL *d, REAL *e, REAL const &vl, REAL const &vu, INTEGER const &il, INTEGER const &iu, REAL const &abstol, INTEGER &m, REAL *w, REAL *z, INTEGER const &ldz, INTEGER *isuppz, REAL *work, INTEGER const &lwork, arr_ref<INTEGER> iwork, INTEGER const &liwork, INTEGER &info) {
    INTEGER ieeeok = 0;
    bool wantz = false;
    bool alleig = false;
    bool valeig = false;
    bool indeig = false;
    bool lquery = false;
    INTEGER lwmin = 0;
    INTEGER liwmin = 0;
    const REAL one = 1.0;
    REAL safmin = 0.0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    REAL bignum = 0.0;
    REAL rmin = 0.0;
    REAL rmax = 0.0;
    INTEGER iscale = 0;
    REAL vll = 0.0;
    REAL vuu = 0.0;
    REAL tnrm = 0.0;
    const REAL zero = 0.0;
    REAL sigma = 0.0;
    INTEGER indibl = 0;
    INTEGER indisp = 0;
    INTEGER indifl = 0;
    INTEGER indiwo = 0;
    bool test = false;
    const REAL two = 2.0e+0;
    bool tryrac = false;
    str<1> order = char0;
    INTEGER nsplit = 0;
    INTEGER imax = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    REAL tmp1 = 0.0;
    INTEGER jj = 0;
    INTEGER itmp1 = 0;
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
    ieeeok = iMlaenv[(10 - 1) + ("Rstevr" - 1) * ldiMlaenv];
    //
    wantz = Mlsame(jobz, "V");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");
    //
    lquery = ((lwork == -1) || (liwork == -1));
    lwmin = max((INTEGER)1, 20 * n);
    liwmin = max((INTEGER)1, 10 * n);
    //
    info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
        info = -1;
    } else if (!(alleig || valeig || indeig)) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else {
        if (valeig) {
            if (n > 0 && vu <= vl) {
                info = -7;
            }
        } else if (indeig) {
            if (il < 1 || il > max((INTEGER)1, n)) {
                info = -8;
            } else if (iu < min(n, il) || iu > n) {
                info = -9;
            }
        }
    }
    if (info == 0) {
        if (ldz < 1 || (wantz && ldz < n)) {
            info = -14;
        }
    }
    //
    if (info == 0) {
        work[1 - 1] = lwmin;
        iwork[1 - 1] = liwmin;
        //
        if (lwork < lwmin && !lquery) {
            info = -17;
        } else if (liwork < liwmin && !lquery) {
            info = -19;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rstevr", -info);
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
    if (n == 1) {
        if (alleig || indeig) {
            m = 1;
            w[1 - 1] = d[1 - 1];
        } else {
            if (vl < d[1 - 1] && vu >= d[1 - 1]) {
                m = 1;
                w[1 - 1] = d[1 - 1];
            }
        }
        if (wantz) {
            z[(1 - 1)] = one;
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
    if (valeig) {
        vll = vl;
        vuu = vu;
    }
    //
    tnrm = Rlanst[("M" - 1) + (n - 1) * ldRlanst];
    if (tnrm > zero && tnrm < rmin) {
        iscale = 1;
        sigma = rmin / tnrm;
    } else if (tnrm > rmax) {
        iscale = 1;
        sigma = rmax / tnrm;
    }
    if (iscale == 1) {
        Rscal(n, sigma, d, 1);
        Rscal(n - 1, sigma, e[1 - 1], 1);
        if (valeig) {
            vll = vl * sigma;
            vuu = vu * sigma;
        }
    }
    //
    //     Initialize indices INTEGERo workspaces.  Note: These indices are used only
    //     if Rsterf or Rstemr fail.
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
    indiwo = indisp + n;
    //
    //     If all eigenvalues are desired, then
    //     call Rsterf or Rstemr.  If this fails for some eigenvalue, then
    //     try Rstebz.
    //
    test = false;
    if (indeig) {
        if (il == 1 && iu == n) {
            test = true;
        }
    }
    if ((alleig || test) && ieeeok == 1) {
        Rcopy(n - 1, e[1 - 1], 1, work[1 - 1], 1);
        if (!wantz) {
            Rcopy(n, d, 1, w, 1);
            Rsterf(n, w, work, info);
        } else {
            Rcopy(n, d, 1, work[(n + 1) - 1], 1);
            if (abstol <= two * n * eps) {
                tryrac = true;
            } else {
                tryrac = false;
            }
            Rstemr(jobz, "A", n, work[(n + 1) - 1], work, vl, vu, il, iu, m, w, z, ldz, n, isuppz, tryrac, work[(2 * n + 1) - 1], lwork - 2 * n, iwork, liwork, info);
            //
        }
        if (info == 0) {
            m = n;
            goto statement_10;
        }
        info = 0;
    }
    //
    //     Otherwise, call Rstebz and, if eigenvectors are desired, Rstein.
    //
    if (wantz) {
        order = "B";
    } else {
        order = "E";
    }
    //
    Rstebz(range, order, n, vll, vuu, il, iu, abstol, d, e, m, nsplit, w, iwork[indibl - 1], iwork[indisp - 1], work, iwork[indiwo - 1], info);
    //
    if (wantz) {
        Rstein(n, d, e, m, w, iwork[indibl - 1], iwork[indisp - 1], z, ldz, work, iwork[indiwo - 1], iwork[indifl - 1], info);
    }
//
//     If matrix was scaled, then rescale eigenvalues appropriately.
//
statement_10:
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
    //     eigenvectors.
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
                itmp1 = iwork[i - 1];
                w[i - 1] = w[j - 1];
                iwork[i - 1] = iwork[j - 1];
                w[j - 1] = tmp1;
                iwork[j - 1] = itmp1;
                Rswap(n, z[(i - 1) * ldz], 1, z[(j - 1) * ldz], 1);
            }
        }
    }
    //
    //      Causes problems with tests 19 & 20:
    //      IF (wantz .and. INDEIG ) Z( 1,1) = Z(1,1) / 1.002 + .002
    //
    work[1 - 1] = lwmin;
    iwork[1 - 1] = liwmin;
    //
    //     End of Rstevr
    //
}
