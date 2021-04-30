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

void Cstemr(const char *jobz, const char *range, INTEGER const n, REAL *d, REAL *e, REAL const vl, REAL const vu, INTEGER const il, INTEGER const iu, INTEGER &m, REAL *w, COMPLEX *z, INTEGER const ldz, INTEGER const nzc, INTEGER *isuppz, bool &tryrac, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
    bool wantz = false;
    bool alleig = false;
    bool valeig = false;
    bool indeig = false;
    bool lquery = false;
    bool zquery = false;
    INTEGER lwmin = 0;
    INTEGER liwmin = 0;
    const REAL zero = 0.0;
    REAL wl = 0.0;
    REAL wu = 0.0;
    INTEGER iil = 0;
    INTEGER iiu = 0;
    INTEGER nsplit = 0;
    REAL safmin = 0.0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL rmin = 0.0;
    REAL rmax = 0.0;
    INTEGER nzcmin = 0;
    INTEGER itmp = 0;
    INTEGER itmp2 = 0;
    REAL r1 = 0.0;
    REAL r2 = 0.0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    INTEGER indgrs = 0;
    INTEGER inderr = 0;
    INTEGER indgp = 0;
    INTEGER indd = 0;
    INTEGER inde2 = 0;
    INTEGER indwrk = 0;
    INTEGER iinspl = 0;
    INTEGER iindbl = 0;
    INTEGER iindw = 0;
    INTEGER iindwk = 0;
    REAL scale = 0.0;
    REAL tnrm = 0.0;
    INTEGER iinfo = 0;
    REAL thresh = 0.0;
    INTEGER j = 0;
    const REAL four = 4.0;
    REAL rtol1 = 0.0;
    REAL rtol2 = 0.0;
    REAL pivmin = 0.0;
    const REAL minrgp = 1.0e-3;
    INTEGER ibegin = 0;
    INTEGER wbegin = 0;
    INTEGER jblk = 0;
    INTEGER iend = 0;
    INTEGER in = 0;
    INTEGER wend = 0;
    INTEGER offset = 0;
    INTEGER ifirst = 0;
    INTEGER ilast = 0;
    INTEGER i = 0;
    REAL tmp = 0.0;
    INTEGER jj = 0;
    //
    //  -- LAPACK computational routine --
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
    //
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
    lquery = ((lwork == -1) || (liwork == -1));
    zquery = (nzc == -1);
    //
    //     Rstemr needs WORK of size 6*N, IWORK of size 3*N.
    //     In addition, Rlarre needs WORK of size 6*N, IWORK of size 5*N.
    //     Furthermore, Clarrv needs WORK of size 12*N, IWORK of size 7*N.
    if (wantz) {
        lwmin = 18 * n;
        liwmin = 10 * n;
    } else {
        //        need less workspace if only the eigenvalues are wanted
        lwmin = 12 * n;
        liwmin = 8 * n;
    }
    //
    wl = zero;
    wu = zero;
    iil = 0;
    iiu = 0;
    nsplit = 0;
    //
    if (valeig) {
        //        We do not reference VL, VU in the cases RANGE = 'I','A'
        //        The interval (WL, WU] contains all the wanted eigenvalues.
        //        It is either given by the user or computed in Rlarre.
        wl = vl;
        wu = vu;
    } else if (indeig) {
        //        We do not reference IL, IU in the cases RANGE = 'V','A'
        iil = il;
        iiu = iu;
    }
    //
    info = 0;
    if (!(wantz || Mlsame(jobz, "N"))) {
        info = -1;
    } else if (!(alleig || valeig || indeig)) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (valeig && n > 0 && wu <= wl) {
        info = -7;
    } else if (indeig && (iil < 1 || iil > n)) {
        info = -8;
    } else if (indeig && (iiu < iil || iiu > n)) {
        info = -9;
    } else if (ldz < 1 || (wantz && ldz < n)) {
        info = -13;
    } else if (lwork < lwmin && !lquery) {
        info = -17;
    } else if (liwork < liwmin && !lquery) {
        info = -19;
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
    if (info == 0) {
        work[1 - 1] = lwmin;
        iwork[1 - 1] = liwmin;
        //
        if (wantz && alleig) {
            nzcmin = n;
        } else if (wantz && valeig) {
            Rlarrc("T", n, vl, vu, d, e, safmin, nzcmin, itmp, itmp2, info);
        } else if (wantz && indeig) {
            nzcmin = iiu - iil + 1;
        } else {
            //           WANTZ .EQ. FALSE.
            nzcmin = 0;
        }
        if (zquery && info == 0) {
            z[(1 - 1)] = nzcmin;
        } else if (nzc < nzcmin && !zquery) {
            info = -14;
        }
    }
    //
    if (info != 0) {
        //
        Mxerbla("Cstemr", -info);
        //
        return;
    } else if (lquery || zquery) {
        return;
    }
    //
    //     Handle N = 0, 1, and 2 cases immediately
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
            if (wl < d[1 - 1] && wu >= d[1 - 1]) {
                m = 1;
                w[1 - 1] = d[1 - 1];
            }
        }
        if (wantz && (!zquery)) {
            z[(1 - 1)] = one;
            isuppz[1 - 1] = 1;
            isuppz[2 - 1] = 1;
        }
        return;
    }
    //
    if (n == 2) {
        if (!wantz) {
            Rlae2(d[1 - 1], &e[1 - 1], &d[2 - 1], r1, r2);
        } else if (wantz && (!zquery)) {
            Rlaev2(d[1 - 1], &e[1 - 1], &d[2 - 1], r1, r2, cs, sn);
        }
        if (alleig || (valeig && (r2 > wl) && (r2 <= wu)) || (indeig && (iil == 1))) {
            m++;
            w[m - 1] = r2;
            if (wantz && (!zquery)) {
                z[(m - 1) * ldz] = -sn;
                z[(2 - 1) + (m - 1) * ldz] = cs;
                //              Note: At most one of SN and CS can be zero.
                if (sn != zero) {
                    if (cs != zero) {
                        isuppz[(2 * m - 1) - 1] = 1;
                        isuppz[(2 * m) - 1] = 2;
                    } else {
                        isuppz[(2 * m - 1) - 1] = 1;
                        isuppz[(2 * m) - 1] = 1;
                    }
                } else {
                    isuppz[(2 * m - 1) - 1] = 2;
                    isuppz[(2 * m) - 1] = 2;
                }
            }
        }
        if (alleig || (valeig && (r1 > wl) && (r1 <= wu)) || (indeig && (iiu == 2))) {
            m++;
            w[m - 1] = r1;
            if (wantz && (!zquery)) {
                z[(m - 1) * ldz] = cs;
                z[(2 - 1) + (m - 1) * ldz] = sn;
                //              Note: At most one of SN and CS can be zero.
                if (sn != zero) {
                    if (cs != zero) {
                        isuppz[(2 * m - 1) - 1] = 1;
                        isuppz[(2 * m) - 1] = 2;
                    } else {
                        isuppz[(2 * m - 1) - 1] = 1;
                        isuppz[(2 * m) - 1] = 1;
                    }
                } else {
                    isuppz[(2 * m - 1) - 1] = 2;
                    isuppz[(2 * m) - 1] = 2;
                }
            }
        }
    } else {
        //
        //        Continue with general N
        //
        indgrs = 1;
        inderr = 2 * n + 1;
        indgp = 3 * n + 1;
        indd = 4 * n + 1;
        inde2 = 5 * n + 1;
        indwrk = 6 * n + 1;
        //
        iinspl = 1;
        iindbl = n + 1;
        iindw = 2 * n + 1;
        iindwk = 3 * n + 1;
        //
        //        Scale matrix to allowable range, if necessary.
        //        The allowable range is related to the PIVMIN parameter; see the
        //        comments in Rlarrd.  The preference for scaling small values
        //        up is heuristic; we expect users' matrices not to be close to the
        //        RMAX threshold.
        //
        scale = one;
        tnrm = Rlanst("M", n, d, e);
        if (tnrm > zero && tnrm < rmin) {
            scale = rmin / tnrm;
        } else if (tnrm > rmax) {
            scale = rmax / tnrm;
        }
        if (scale != one) {
            Rscal(n, scale, d, 1);
            Rscal(n - 1, scale, e, 1);
            tnrm = tnrm * scale;
            if (valeig) {
                //              If eigenvalues in interval have to be found,
                //              scale (WL, WU] accordingly
                wl = wl * scale;
                wu = wu * scale;
            }
        }
        //
        //        Compute the desired eigenvalues of the tridiagonal after splitting
        //        into smaller subblocks if the corresponding off-diagonal elements
        //        are small
        //        THRESH is the splitting parameter for Rlarre
        //        A negative THRESH forces the old splitting criterion based on the
        //        size of the off-diagonal. A positive THRESH switches to splitting
        //        which preserves relative accuracy.
        //
        if (tryrac) {
            //           Test whether the matrix warrants the more expensive relative approach.
            Rlarrr(n, d, e, iinfo);
        } else {
            //           The user does not care about relative accurately eigenvalues
            iinfo = -1;
        }
        //        Set the splitting criterion
        if (iinfo == 0) {
            thresh = eps;
        } else {
            thresh = -eps;
            //           relative accuracy is desired but T does not guarantee it
            tryrac = false;
        }
        //
        if (tryrac) {
            //           Copy original diagonal, needed to guarantee relative accuracy
            Rcopy(n, d, 1, &work[indd - 1], 1);
        }
        //        Store the squares of the offdiagonal values of T
        for (j = 1; j <= n - 1; j = j + 1) {
            work[(inde2 + j - 1) - 1] = pow2(e[j - 1]);
        }
        //
        //        Set the tolerance parameters for bisection
        if (!wantz) {
            //           Rlarre computes the eigenvalues to full precision.
            rtol1 = four * eps;
            rtol2 = four * eps;
        } else {
            //           Rlarre computes the eigenvalues to less than full precision.
            //           Clarrv will refine the eigenvalue approximations, and we only
            //           need less accurate initial bisection in Rlarre.
            //           Note: these settings do only affect the subset case and Rlarre
            rtol1 = sqrt(eps);
            rtol2 = max(sqrt(eps) * 5.0e-3, four * eps);
        }
        Rlarre(range, n, wl, wu, iil, iiu, d, e, &work[inde2 - 1], rtol1, rtol2, thresh, nsplit, &iwork[iinspl - 1], m, w, &work[inderr - 1], &work[indgp - 1], &iwork[iindbl - 1], &iwork[iindw - 1], &work[indgrs - 1], pivmin, &work[indwrk - 1], &iwork[iindwk - 1], iinfo);
        if (iinfo != 0) {
            info = 10 + abs(iinfo);
            return;
        }
        //        Note that if RANGE .NE. 'V', Rlarre computes bounds on the desired
        //        part of the spectrum. All desired eigenvalues are contained in
        //        (WL,WU]
        //
        if (wantz) {
            //
            //           Compute the desired eigenvectors corresponding to the computed
            //           eigenvalues
            //
            Clarrv(n, wl, wu, d, e, pivmin, &iwork[iinspl - 1], m, 1, m, minrgp, rtol1, rtol2, w, &work[inderr - 1], &work[indgp - 1], &iwork[iindbl - 1], &iwork[iindw - 1], &work[indgrs - 1], z, ldz, isuppz, &work[indwrk - 1], &iwork[iindwk - 1], iinfo);
            if (iinfo != 0) {
                info = 20 + abs(iinfo);
                return;
            }
        } else {
            //           Rlarre computes eigenvalues of the (shifted) root representation
            //           Clarrv returns the eigenvalues of the unshifted matrix.
            //           However, if the eigenvectors are not desired by the user, we need
            //           to apply the corresponding shifts from Rlarre to obtain the
            //           eigenvalues of the original matrix.
            for (j = 1; j <= m; j = j + 1) {
                itmp = iwork[(iindbl + j - 1) - 1];
                w[j - 1] += e[(iwork[(iinspl + itmp - 1) - 1]) - 1];
            }
        }
        //
        if (tryrac) {
            //           Refine computed eigenvalues so that they are relatively accurate
            //           with respect to the original matrix T.
            ibegin = 1;
            wbegin = 1;
            for (jblk = 1; jblk <= iwork[(iindbl + m - 1) - 1]; jblk = jblk + 1) {
                iend = iwork[(iinspl + jblk - 1) - 1];
                in = iend - ibegin + 1;
                wend = wbegin - 1;
            //              check if any eigenvalues have to be refined in this block
            statement_36:
                if (wend < m) {
                    if (iwork[(iindbl + wend) - 1] == jblk) {
                        wend++;
                        goto statement_36;
                    }
                }
                if (wend < wbegin) {
                    ibegin = iend + 1;
                    goto statement_39;
                }
                //
                offset = iwork[(iindw + wbegin - 1) - 1] - 1;
                ifirst = iwork[(iindw + wbegin - 1) - 1];
                ilast = iwork[(iindw + wend - 1) - 1];
                rtol2 = four * eps;
                Rlarrj(in, &work[(indd + ibegin - 1) - 1], &work[(inde2 + ibegin - 1) - 1], ifirst, ilast, rtol2, offset, &w[wbegin - 1], &work[(inderr + wbegin - 1) - 1], &work[indwrk - 1], &iwork[iindwk - 1], pivmin, tnrm, iinfo);
                ibegin = iend + 1;
                wbegin = wend + 1;
            statement_39:;
            }
        }
        //
        //        If matrix was scaled, then rescale eigenvalues appropriately.
        //
        if (scale != one) {
            Rscal(m, one / scale, w, 1);
        }
    }
    //
    //     If eigenvalues are not in increasing order, then sort them,
    //     possibly along with eigenvectors.
    //
    if (nsplit > 1 || n == 2) {
        if (!wantz) {
            Rlasrt("I", m, w, iinfo);
            if (iinfo != 0) {
                info = 3;
                return;
            }
        } else {
            for (j = 1; j <= m - 1; j = j + 1) {
                i = 0;
                tmp = w[j - 1];
                for (jj = j + 1; jj <= m; jj = jj + 1) {
                    if (w[jj - 1] < tmp) {
                        i = jj;
                        tmp = w[jj - 1];
                    }
                }
                if (i != 0) {
                    w[i - 1] = w[j - 1];
                    w[j - 1] = tmp;
                    if (wantz) {
                        Cswap(n, &z[(i - 1) * ldz], 1, &z[(j - 1) * ldz], 1);
                        itmp = isuppz[(2 * i - 1) - 1];
                        isuppz[(2 * i - 1) - 1] = isuppz[(2 * j - 1) - 1];
                        isuppz[(2 * j - 1) - 1] = itmp;
                        itmp = isuppz[(2 * i) - 1];
                        isuppz[(2 * i) - 1] = isuppz[(2 * j) - 1];
                        isuppz[(2 * j) - 1] = itmp;
                    }
                }
            }
        }
    }
    //
    work[1 - 1] = lwmin;
    iwork[1 - 1] = liwmin;
    //
    //     End of Cstemr
    //
}
