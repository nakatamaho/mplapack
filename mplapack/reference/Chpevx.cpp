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

void Chpevx(const char *jobz, const char *range, const char *uplo, INTEGER const n, COMPLEX *ap, REAL const vl, REAL const vu, INTEGER const il, INTEGER const iu, REAL const abstol, INTEGER &m, REAL *w, COMPLEX *z, INTEGER const ldz, COMPLEX *work, REAL *rwork, INTEGER *iwork, INTEGER *ifail, INTEGER &info) {
    bool wantz = false;
    bool alleig = false;
    bool valeig = false;
    bool indeig = false;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    REAL safmin = 0.0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL rmin = 0.0;
    REAL rmax = 0.0;
    INTEGER iscale = 0;
    REAL abstll = 0.0;
    REAL vll = 0.0;
    REAL vuu = 0.0;
    const REAL zero = 0.0;
    REAL anrm = 0.0;
    REAL sigma = 0.0;
    INTEGER indd = 0;
    INTEGER inde = 0;
    INTEGER indrwk = 0;
    INTEGER indtau = 0;
    INTEGER indwrk = 0;
    INTEGER iinfo = 0;
    bool test = false;
    INTEGER indee = 0;
    INTEGER i = 0;
    char order;
    INTEGER indibl = 0;
    INTEGER indisp = 0;
    INTEGER indiwk = 0;
    INTEGER nsplit = 0;
    INTEGER imax = 0;
    INTEGER j = 0;
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
    wantz = Mlsame(jobz, "V");
    alleig = Mlsame(range, "A");
    valeig = Mlsame(range, "V");
    indeig = Mlsame(range, "I");
    //
    info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
        info = -1;
    } else if (!(alleig || valeig || indeig)) {
        info = -2;
    } else if (!(Mlsame(uplo, "L") || Mlsame(uplo, "U"))) {
        info = -3;
    } else if (n < 0) {
        info = -4;
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
    if (info != 0) {
        Mxerbla("Chpevx", -info);
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
            w[1 - 1] = ap[1 - 1].real();
        } else {
            if (vl < ap[1 - 1].real() && vu >= ap[1 - 1].real()) {
                m = 1;
                w[1 - 1] = ap[1 - 1].real();
            }
        }
        if (wantz) {
            z[(1 - 1)] = cone;
        }
        return;
    }
    //
    //     Get machine constants.
    //
    safmin = Rlamch("Safe minimum");
    eps = Rlamch("Precision");
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
    } else {
        vll = zero;
        vuu = zero;
    }
    anrm = Clanhp("M", uplo, n, ap, rwork);
    if (anrm > zero && anrm < rmin) {
        iscale = 1;
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        iscale = 1;
        sigma = rmax / anrm;
    }
    if (iscale == 1) {
        CRscal((n * (n + 1)) / 2, sigma, ap, 1);
        if (abstol > 0) {
            abstll = abstol * sigma;
        }
        if (valeig) {
            vll = vl * sigma;
            vuu = vu * sigma;
        }
    }
    //
    //     Call Chptrd to reduce Hermitian packed matrix to tridiagonal form.
    //
    indd = 1;
    inde = indd + n;
    indrwk = inde + n;
    indtau = 1;
    indwrk = indtau + n;
    Chptrd(uplo, n, ap, &rwork[indd - 1], &rwork[inde - 1], &work[indtau - 1], iinfo);
    //
    //     If all eigenvalues are desired and ABSTOL is less than or equal
    //     to zero, then call Rsterf or Cupgtr and Csteqr.  If this fails
    //     for some eigenvalue, then try Rstebz.
    //
    test = false;
    if (indeig) {
        if (il == 1 && iu == n) {
            test = true;
        }
    }
    if ((alleig || test) && (abstol <= zero)) {
        Rcopy(n, &rwork[indd - 1], 1, w, 1);
        indee = indrwk + 2 * n;
        if (!wantz) {
            Rcopy(n - 1, &rwork[inde - 1], 1, &rwork[indee - 1], 1);
            Rsterf(n, w, &rwork[indee - 1], info);
        } else {
            Cupgtr(uplo, n, ap, &work[indtau - 1], z, ldz, &work[indwrk - 1], iinfo);
            Rcopy(n - 1, &rwork[inde - 1], 1, &rwork[indee - 1], 1);
            Csteqr(jobz, n, w, &rwork[indee - 1], z, ldz, &rwork[indrwk - 1], info);
            if (info == 0) {
                for (i = 1; i <= n; i = i + 1) {
                    ifail[i - 1] = 0;
                }
            }
        }
        if (info == 0) {
            m = n;
            goto statement_20;
        }
        info = 0;
    }
    //
    //     Otherwise, call Rstebz and, if eigenvectors are desired, Cstein.
    //
    if (wantz) {
        order = 'B';
    } else {
        order = 'E';
    }
    indibl = 1;
    indisp = indibl + n;
    indiwk = indisp + n;
    Rstebz(range, &order, n, vll, vuu, il, iu, abstll, &rwork[indd - 1], &rwork[inde - 1], m, nsplit, w, &iwork[indibl - 1], &iwork[indisp - 1], &rwork[indrwk - 1], &iwork[indiwk - 1], info);
    //
    if (wantz) {
        Cstein(n, &rwork[indd - 1], &rwork[inde - 1], m, w, &iwork[indibl - 1], &iwork[indisp - 1], z, ldz, &rwork[indrwk - 1], &iwork[indiwk - 1], ifail, info);
        //
        //        Apply unitary matrix used in reduction to tridiagonal
        //        form to eigenvectors returned by Cstein.
        //
        indwrk = indtau + n;
        Cupmtr("L", uplo, "N", n, m, ap, &work[indtau - 1], z, ldz, &work[indwrk - 1], iinfo);
    }
//
//     If matrix was scaled, then rescale eigenvalues appropriately.
//
statement_20:
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
                itmp1 = iwork[(indibl + i - 1) - 1];
                w[i - 1] = w[j - 1];
                iwork[(indibl + i - 1) - 1] = iwork[(indibl + j - 1) - 1];
                w[j - 1] = tmp1;
                iwork[(indibl + j - 1) - 1] = itmp1;
                Cswap(n, &z[(i - 1) * ldz], 1, &z[(j - 1) * ldz], 1);
                if (info != 0) {
                    itmp1 = ifail[i - 1];
                    ifail[i - 1] = ifail[j - 1];
                    ifail[j - 1] = itmp1;
                }
            }
        }
    }
    //
    //     End of Chpevx
    //
}
