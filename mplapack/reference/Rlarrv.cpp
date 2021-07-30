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

void Rlarrv(INTEGER const n, REAL const vl, REAL const /* vu */, REAL *d, REAL *l, REAL const pivmin, INTEGER *isplit, INTEGER const m, INTEGER const dol, INTEGER const dou, REAL const minrgp, REAL &rtol1, REAL &rtol2, REAL *w, REAL *werr, REAL *wgap, INTEGER *iblock, INTEGER *indexw, REAL *gers, REAL *z, INTEGER const ldz, INTEGER *isuppz, REAL *work, INTEGER *iwork, INTEGER &info) {
    INTEGER indld = 0;
    INTEGER indlld = 0;
    INTEGER indwrk = 0;
    INTEGER minwsize = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    INTEGER iindr = 0;
    INTEGER iindc1 = 0;
    INTEGER iindc2 = 0;
    INTEGER iindwk = 0;
    INTEGER miniwsize = 0;
    INTEGER zusedl = 0;
    INTEGER zusedu = 0;
    INTEGER zusedw = 0;
    REAL eps = 0.0;
    const REAL two = 2.0;
    REAL rqtol = 0.0;
    bool tryrqc = false;
    const REAL four = 4.0;
    INTEGER done = 0;
    INTEGER ibegin = 0;
    INTEGER wbegin = 0;
    INTEGER jblk = 0;
    INTEGER iend = 0;
    REAL sigma = 0.0;
    INTEGER wend = 0;
    REAL gl = 0.0;
    REAL gu = 0.0;
    REAL spdiam = 0.0;
    INTEGER oldien = 0;
    INTEGER in = 0;
    INTEGER im = 0;
    const REAL one = 1.0;
    INTEGER ndepth = 0;
    INTEGER parity = 0;
    INTEGER nclus = 0;
    INTEGER idone = 0;
    INTEGER oldncl = 0;
    INTEGER oldcls = 0;
    INTEGER newcls = 0;
    INTEGER j = 0;
    INTEGER oldfst = 0;
    INTEGER oldlst = 0;
    REAL tmp = 0.0;
    INTEGER p = 0;
    INTEGER q = 0;
    INTEGER offset = 0;
    INTEGER iinfo = 0;
    INTEGER newfst = 0;
    INTEGER newlst = 0;
    INTEGER newsiz = 0;
    INTEGER newftt = 0;
    REAL lgap = 0.0;
    REAL rgap = 0.0;
    INTEGER k = 0;
    REAL tau = 0.0;
    REAL ssigma = 0.0;
    const REAL three = 3.0;
    REAL fudge = 0.0;
    INTEGER iter = 0;
    REAL tol = 0.0;
    INTEGER windex = 0;
    INTEGER windmn = 0;
    INTEGER windpl = 0;
    REAL lambda = 0.0;
    bool eskip = false;
    REAL left = 0.0;
    REAL right = 0.0;
    INTEGER indeig = 0;
    REAL gap = 0.0;
    REAL gaptol = 0.0;
    INTEGER isupmn = 0;
    INTEGER isupmx = 0;
    REAL savgap = 0.0;
    bool usedbs = false;
    bool usedrq = false;
    bool needbs = false;
    INTEGER itmp1 = 0;
    INTEGER negcnt = 0;
    REAL ztz = 0.0;
    REAL mingma = 0.0;
    REAL nrminv = 0.0;
    REAL resid = 0.0;
    REAL rqcorr = 0.0;
    REAL bstres = 0.0;
    REAL bstw = 0.0;
    REAL sgndef = 0.0;
    const REAL half = 0.5e0;
    const INTEGER maxitr = 10;
    bool stp2ii = false;
    INTEGER zfrom = 0;
    INTEGER zto = 0;
    INTEGER ii = 0;
    //
    //  -- LAPACK auxiliary routine --
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
    //     ..
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    //     The first N entries of WORK are reserved for the eigenvalues
    indld = n + 1;
    indlld = 2 * n + 1;
    indwrk = 3 * n + 1;
    minwsize = 12 * n;
    //
    for (i = 1; i <= minwsize; i = i + 1) {
        work[i - 1] = zero;
    }
    //
    //     IWORK(IINDR+1:IINDR+N) hold the twist indices R for the
    //     factorization used to compute the FP vector
    iindr = 0;
    //     IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current
    //     layer and the one above.
    iindc1 = n;
    iindc2 = 2 * n;
    iindwk = 3 * n + 1;
    //
    miniwsize = 7 * n;
    for (i = 1; i <= miniwsize; i = i + 1) {
        iwork[i - 1] = 0;
    }
    //
    zusedl = 1;
    if (dol > 1) {
        //        Set lower bound for use of Z
        zusedl = dol - 1;
    }
    zusedu = m;
    if (dou < m) {
        //        Set lower bound for use of Z
        zusedu = dou + 1;
    }
    //     The width of the part of Z that is used
    zusedw = zusedu - zusedl + 1;
    //
    Rlaset("Full", n, zusedw, zero, zero, &z[(zusedl - 1) * ldz], ldz);
    //
    eps = Rlamch("Precision");
    rqtol = two * eps;
    //
    //     Set expert flags for standard code.
    tryrqc = true;
    //
    if ((dol == 1) && (dou == m)) {
    } else {
        //        Only selected eigenpairs are computed. Since the other evalues
        //        are not refined by RQ iteration, bisection has to compute to full
        //        accuracy.
        rtol1 = four * eps;
        rtol2 = four * eps;
    }
    //
    //     The entries WBEGIN:WEND in W, WERR, WGAP correspond to the
    //     desired eigenvalues. The support of the nonzero eigenvector
    //     entries is contained in the interval IBEGIN:IEND.
    //     Remark that if k eigenpairs are desired, then the eigenvectors
    //     are stored in k contiguous columns of Z.
    //
    //     DONE is the number of eigenvectors already computed
    done = 0;
    ibegin = 1;
    wbegin = 1;
    for (jblk = 1; jblk <= iblock[m - 1]; jblk = jblk + 1) {
        iend = isplit[jblk - 1];
        sigma = l[iend - 1];
        //        Find the eigenvectors of the submatrix indexed IBEGIN
        //        through IEND.
        wend = wbegin - 1;
    statement_15:
        if (wend < m) {
            if (iblock[(wend + 1) - 1] == jblk) {
                wend++;
                goto statement_15;
            }
        }
        if (wend < wbegin) {
            ibegin = iend + 1;
            goto statement_170;
        } else if ((wend < dol) || (wbegin > dou)) {
            ibegin = iend + 1;
            wbegin = wend + 1;
            goto statement_170;
        }
        //
        //        Find local spectral diameter of the block
        gl = gers[(2 * ibegin - 1) - 1];
        gu = gers[(2 * ibegin) - 1];
        for (i = ibegin + 1; i <= iend; i = i + 1) {
            gl = min(gers[(2 * i - 1) - 1], gl);
            gu = max(gers[(2 * i) - 1], gu);
        }
        spdiam = gu - gl;
        //
        //        OLDIEN is the last index of the previous block
        oldien = ibegin - 1;
        //        Calculate the size of the current block
        in = iend - ibegin + 1;
        //        The number of eigenvalues in the current block
        im = wend - wbegin + 1;
        //
        //        This is for a 1x1 block
        if (ibegin == iend) {
            done++;
            z[(ibegin - 1) + (wbegin - 1) * ldz] = one;
            isuppz[(2 * wbegin - 1) - 1] = ibegin;
            isuppz[(2 * wbegin) - 1] = ibegin;
            w[wbegin - 1] += sigma;
            work[wbegin - 1] = w[wbegin - 1];
            ibegin = iend + 1;
            wbegin++;
            goto statement_170;
        }
        //
        //        The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND)
        //        Note that these can be approximations, in this case, the corresp.
        //        entries of WERR give the size of the uncertainty interval.
        //        The eigenvalue approximations will be refined when necessary as
        //        high relative accuracy is required for the computation of the
        //        corresponding eigenvectors.
        Rcopy(im, &w[wbegin - 1], 1, &work[wbegin - 1], 1);
        //
        //        We store in W the eigenvalue approximations w.r.t. the original
        //        matrix T.
        for (i = 1; i <= im; i = i + 1) {
            w[(wbegin + i - 1) - 1] += sigma;
        }
        //
        //        NDEPTH is the current depth of the representation tree
        ndepth = 0;
        //        PARITY is either 1 or 0
        parity = 1;
        //        NCLUS is the number of clusters for the next level of the
        //        representation tree, we start with NCLUS = 1 for the root
        nclus = 1;
        iwork[(iindc1 + 1) - 1] = 1;
        iwork[(iindc1 + 2) - 1] = im;
        //
        //        IDONE is the number of eigenvectors already computed in the current
        //        block
        idone = 0;
    //        loop while( IDONE.LT.IM )
    //        generate the representation tree for the current block and
    //        compute the eigenvectors
    statement_40:
        if (idone < im) {
            //           This is a crude protection against infinitely deep trees
            if (ndepth > m) {
                info = -2;
                return;
            }
            //           breadth first processing of the current level of the representation
            //           tree: OLDNCL = number of clusters on current level
            oldncl = nclus;
            //           reset NCLUS to count the number of child clusters
            nclus = 0;
            //
            parity = 1 - parity;
            if (parity == 0) {
                oldcls = iindc1;
                newcls = iindc2;
            } else {
                oldcls = iindc2;
                newcls = iindc1;
            }
            //           Process the clusters on the current level
            for (i = 1; i <= oldncl; i = i + 1) {
                j = oldcls + 2 * i;
                //              OLDFST, OLDLST = first, last index of current cluster.
                //                               cluster indices start with 1 and are relative
                //                               to WBEGIN when accessing W, WGAP, WERR, Z
                oldfst = iwork[(j - 1) - 1];
                oldlst = iwork[j - 1];
                if (ndepth > 0) {
                    //                 Retrieve relatively robust representation (RRR) of cluster
                    //                 that has been computed at the previous level
                    //                 The RRR is stored in Z and overwritten once the eigenvectors
                    //                 have been computed or when the cluster is refined
                    //
                    if ((dol == 1) && (dou == m)) {
                        //                    Get representation from location of the leftmost evalue
                        //                    of the cluster
                        j = wbegin + oldfst - 1;
                    } else {
                        if (wbegin + oldfst - 1 < dol) {
                            //                       Get representation from the left end of Z array
                            j = dol - 1;
                        } else if (wbegin + oldfst - 1 > dou) {
                            //                       Get representation from the right end of Z array
                            j = dou;
                        } else {
                            j = wbegin + oldfst - 1;
                        }
                    }
                    Rcopy(in, &z[(ibegin - 1) + (j - 1) * ldz], 1, &d[ibegin - 1], 1);
                    Rcopy(in - 1, &z[(ibegin - 1) + ((j + 1) - 1) * ldz], 1, &l[ibegin - 1], 1);
                    sigma = z[(iend - 1) + ((j + 1) - 1) * ldz];
                    //
                    //                 Set the corresponding entries in Z to zero
                    Rlaset("Full", in, 2, zero, zero, &z[(ibegin - 1) + (j - 1) * ldz], ldz);
                }
                //
                //              Compute DL and DLL of current RRR
                for (j = ibegin; j <= iend - 1; j = j + 1) {
                    tmp = d[j - 1] * l[j - 1];
                    work[(indld - 1 + j) - 1] = tmp;
                    work[(indlld - 1 + j) - 1] = tmp * l[j - 1];
                }
                //
                if (ndepth > 0) {
                    //                 P and Q are index of the first and last eigenvalue to compute
                    //                 within the current block
                    p = indexw[(wbegin - 1 + oldfst) - 1];
                    q = indexw[(wbegin - 1 + oldlst) - 1];
                    //                 Offset for the arrays WORK, WGAP and WERR, i.e., the P-OFFSET
                    //                 through the Q-OFFSET elements of these arrays are to be used.
                    //                  OFFSET = P-OLDFST
                    offset = indexw[wbegin - 1] - 1;
                    //                 perform limited bisection (if necessary) to get approximate
                    //                 eigenvalues to the precision needed.
                    Rlarrb(in, &d[ibegin - 1], &work[(indlld + ibegin - 1) - 1], p, q, rtol1, rtol2, offset, &work[wbegin - 1], &wgap[wbegin - 1], &werr[wbegin - 1], &work[indwrk - 1], &iwork[iindwk - 1], pivmin, spdiam, in, iinfo);
                    if (iinfo != 0) {
                        info = -1;
                        return;
                    }
                    //                 We also recompute the extremal gaps. W holds all eigenvalues
                    //                 of the unshifted matrix and must be used for computation
                    //                 of WGAP, the entries of WORK might stem from RRRs with
                    //                 different shifts. The gaps from WBEGIN-1+OLDFST to
                    //                 WBEGIN-1+OLDLST are correctly computed in Rlarrb.
                    //                 However, we only allow the gaps to become greater since
                    //                 this is what should happen when we decrease WERR
                    if (oldfst > 1) {
                        wgap[(wbegin + oldfst - 2) - 1] = max(wgap[(wbegin + oldfst - 2) - 1], REAL(w[(wbegin + oldfst - 1) - 1] - werr[(wbegin + oldfst - 1) - 1] - w[(wbegin + oldfst - 2) - 1] - werr[(wbegin + oldfst - 2) - 1]));
                    }
                    if (wbegin + oldlst - 1 < wend) {
                        wgap[(wbegin + oldlst - 1) - 1] = max(wgap[(wbegin + oldlst - 1) - 1], REAL(w[(wbegin + oldlst) - 1] - werr[(wbegin + oldlst) - 1] - w[(wbegin + oldlst - 1) - 1] - werr[(wbegin + oldlst - 1) - 1]));
                    }
                    //                 Each time the eigenvalues in WORK get refined, we store
                    //                 the newly found approximation with all shifts applied in W
                    for (j = oldfst; j <= oldlst; j = j + 1) {
                        w[(wbegin + j - 1) - 1] = work[(wbegin + j - 1) - 1] + sigma;
                    }
                }
                //
                //              Process the current node.
                newfst = oldfst;
                for (j = oldfst; j <= oldlst; j = j + 1) {
                    if (j == oldlst) {
                        //                    we are at the right end of the cluster, this is also the
                        //                    boundary of the child cluster
                        newlst = j;
                    } else if (wgap[(wbegin + j - 1) - 1] >= minrgp * abs(work[(wbegin + j - 1) - 1])) {
                        //                    the right relative gap is big enough, the child cluster
                        //                    (NEWFST,..,NEWLST) is well separated from the following
                        newlst = j;
                    } else {
                        //                    inside a child cluster, the relative gap is not
                        //                    big enough.
                        goto statement_140;
                    }
                    //
                    //                 Compute size of child cluster found
                    newsiz = newlst - newfst + 1;
                    //
                    //                 NEWFTT is the place in Z where the new RRR or the computed
                    //                 eigenvector is to be stored
                    if ((dol == 1) && (dou == m)) {
                        //                    Store representation at location of the leftmost evalue
                        //                    of the cluster
                        newftt = wbegin + newfst - 1;
                    } else {
                        if (wbegin + newfst - 1 < dol) {
                            //                       Store representation at the left end of Z array
                            newftt = dol - 1;
                        } else if (wbegin + newfst - 1 > dou) {
                            //                       Store representation at the right end of Z array
                            newftt = dou;
                        } else {
                            newftt = wbegin + newfst - 1;
                        }
                    }
                    //
                    if (newsiz > 1) {
                        //
                        //                    Current child is not a singleton but a cluster.
                        //                    Compute and store new representation of child.
                        //
                        //                    Compute left and right cluster gap.
                        //
                        //                    LGAP and RGAP are not computed from WORK because
                        //                    the eigenvalue approximations may stem from RRRs
                        //                    different shifts. However, W hold all eigenvalues
                        //                    of the unshifted matrix. Still, the entries in WGAP
                        //                    have to be computed from WORK since the entries
                        //                    in W might be of the same order so that gaps are not
                        //                    exhibited correctly for very close eigenvalues.
                        if (newfst == 1) {
                            lgap = max(zero, REAL(w[wbegin - 1] - werr[wbegin - 1] - vl));
                        } else {
                            lgap = wgap[(wbegin + newfst - 2) - 1];
                        }
                        rgap = wgap[(wbegin + newlst - 1) - 1];
                        //
                        //                    Compute left- and rightmost eigenvalue of child
                        //                    to high precision in order to shift as close
                        //                    as possible and obtain as large relative gaps
                        //                    as possible
                        //
                        for (k = 1; k <= 2; k = k + 1) {
                            if (k == 1) {
                                p = indexw[(wbegin - 1 + newfst) - 1];
                            } else {
                                p = indexw[(wbegin - 1 + newlst) - 1];
                            }
                            offset = indexw[wbegin - 1] - 1;
                            Rlarrb(in, &d[ibegin - 1], &work[(indlld + ibegin - 1) - 1], p, p, rqtol, rqtol, offset, &work[wbegin - 1], &wgap[wbegin - 1], &werr[wbegin - 1], &work[indwrk - 1], &iwork[iindwk - 1], pivmin, spdiam, in, iinfo);
                        }
                        //
                        if ((wbegin + newlst - 1 < dol) || (wbegin + newfst - 1 > dou)) {
                            //                       if the cluster contains no desired eigenvalues
                            //                       skip the computation of that branch of the rep. tree
                            //
                            //                       We could skip before the refinement of the extremal
                            //                       eigenvalues of the child, but then the representation
                            //                       tree could be different from the one when nothing is
                            //                       skipped. For this reason we skip at this place.
                            idone += newlst - newfst + 1;
                            goto statement_139;
                        }
                        //
                        //                    Compute RRR of child cluster.
                        //                    Note that the new RRR is stored in Z
                        //
                        //                    Rlarrf needs LWORK = 2*N
                        Rlarrf(in, &d[ibegin - 1], &l[ibegin - 1], &work[(indld + ibegin - 1) - 1], newfst, newlst, &work[wbegin - 1], &wgap[wbegin - 1], &werr[wbegin - 1], spdiam, lgap, rgap, pivmin, tau, &z[(ibegin - 1) + (newftt - 1) * ldz], &z[(ibegin - 1) + ((newftt + 1) - 1) * ldz], &work[indwrk - 1], iinfo);
                        if (iinfo == 0) {
                            //                       a new RRR for the cluster was found by Rlarrf
                            //                       update shift and store it
                            ssigma = sigma + tau;
                            z[(iend - 1) + ((newftt + 1) - 1) * ldz] = ssigma;
                            //                       WORK() are the midpoints and WERR() the semi-width
                            //                       Note that the entries in W are unchanged.
                            for (k = newfst; k <= newlst; k = k + 1) {
                                fudge = three * eps * abs(work[(wbegin + k - 1) - 1]);
                                work[(wbegin + k - 1) - 1] = work[(wbegin + k - 1) - 1] - tau;
                                fudge += four * eps * abs(work[(wbegin + k - 1) - 1]);
                                //                          Fudge errors
                                werr[(wbegin + k - 1) - 1] += fudge;
                                //                          Gaps are not fudged. Provided that WERR is small
                                //                          when eigenvalues are close, a zero gap indicates
                                //                          that a new representation is needed for resolving
                                //                          the cluster. A fudge could lead to a wrong decision
                                //                          of judging eigenvalues 'separated' which in
                                //                          reality are not. This could have a negative impact
                                //                          on the orthogonality of the computed eigenvectors.
                            }
                            //
                            nclus++;
                            k = newcls + 2 * nclus;
                            iwork[(k - 1) - 1] = newfst;
                            iwork[k - 1] = newlst;
                        } else {
                            info = -2;
                            return;
                        }
                    } else {
                        //
                        //                    Compute eigenvector of singleton
                        //
                        iter = 0;
                        //
                        tol = four * log(castREAL(in)) * eps;
                        //
                        k = newfst;
                        windex = wbegin + k - 1;
                        windmn = max(windex - 1, (INTEGER)1);
                        windpl = min(windex + 1, m);
                        lambda = work[windex - 1];
                        done++;
                        //                    Check if eigenvector computation is to be skipped
                        if ((windex < dol) || (windex > dou)) {
                            eskip = true;
                            goto statement_125;
                        } else {
                            eskip = false;
                        }
                        left = work[windex - 1] - werr[windex - 1];
                        right = work[windex - 1] + werr[windex - 1];
                        indeig = indexw[windex - 1];
                        //                    Note that since we compute the eigenpairs for a child,
                        //                    all eigenvalue approximations are w.r.t the same shift.
                        //                    In this case, the entries in WORK should be used for
                        //                    computing the gaps since they exhibit even very small
                        //                    differences in the eigenvalues, as opposed to the
                        //                    entries in W which might "look" the same.
                        //
                        if (k == 1) {
                            //                       In the case RANGE='I' and with not much initial
                            //                       accuracy in LAMBDA and VL, the formula
                            //                       LGAP = MAX( ZERO, (SIGMA - VL) + LAMBDA )
                            //                       can lead to an overestimation of the left gap and
                            //                       thus to inadequately early RQI 'convergence'.
                            //                       Prevent this by forcing a small left gap.
                            lgap = eps * max(abs(left), abs(right));
                        } else {
                            lgap = wgap[windmn - 1];
                        }
                        if (k == im) {
                            //                       In the case RANGE='I' and with not much initial
                            //                       accuracy in LAMBDA and VU, the formula
                            //                       can lead to an overestimation of the right gap and
                            //                       thus to inadequately early RQI 'convergence'.
                            //                       Prevent this by forcing a small right gap.
                            rgap = eps * max(abs(left), abs(right));
                        } else {
                            rgap = wgap[windex - 1];
                        }
                        gap = min(lgap, rgap);
                        if ((k == 1) || (k == im)) {
                            //                       The eigenvector support can become wrong
                            //                       because significant entries could be cut off due to a
                            //                       large GAPTOL parameter in LAR1V. Prevent this.
                            gaptol = zero;
                        } else {
                            gaptol = gap * eps;
                        }
                        isupmn = in;
                        isupmx = 1;
                        //                    Update WGAP so that it holds the minimum gap
                        //                    to the left or the right. This is crucial in the
                        //                    case where bisection is used to ensure that the
                        //                    eigenvalue is refined up to the required precision.
                        //                    The correct value is restored afterwards.
                        savgap = wgap[windex - 1];
                        wgap[windex - 1] = gap;
                        //                    We want to use the Rayleigh Quotient Correction
                        //                    as often as possible since it converges quadratically
                        //                    when we are close enough to the desired eigenvalue.
                        //                    However, the Rayleigh Quotient can have the wrong sign
                        //                    and lead us away from the desired eigenvalue. In this
                        //                    case, the best we can do is to use bisection.
                        usedbs = false;
                        usedrq = false;
                        //                    Bisection is initially turned off unless it is forced
                        needbs = !tryrqc;
                    statement_120:
                        //                    Check if bisection should be used to refine eigenvalue
                        if (needbs) {
                            //                       Take the bisection as new iterate
                            usedbs = true;
                            itmp1 = iwork[(iindr + windex) - 1];
                            offset = indexw[wbegin - 1] - 1;
                            Rlarrb(in, &d[ibegin - 1], &work[(indlld + ibegin - 1) - 1], indeig, indeig, zero, two * eps, offset, &work[wbegin - 1], &wgap[wbegin - 1], &werr[wbegin - 1], &work[indwrk - 1], &iwork[iindwk - 1], pivmin, spdiam, itmp1, iinfo);
                            if (iinfo != 0) {
                                info = -3;
                                return;
                            }
                            lambda = work[windex - 1];
                            //                       Reset twist index from inaccurate LAMBDA to
                            //                       force computation of true MINGMA
                            iwork[(iindr + windex) - 1] = 0;
                        }
                        //                    Given LAMBDA, compute the eigenvector.
                        Rlar1v(in, 1, in, lambda, &d[ibegin - 1], &l[ibegin - 1], &work[(indld + ibegin - 1) - 1], &work[(indlld + ibegin - 1) - 1], pivmin, gaptol, &z[(ibegin - 1) + (windex - 1) * ldz], !usedbs, negcnt, ztz, mingma, iwork[(iindr + windex) - 1], &isuppz[(2 * windex - 1) - 1], nrminv, resid, rqcorr, &work[indwrk - 1]);
                        if (iter == 0) {
                            bstres = resid;
                            bstw = lambda;
                        } else if (resid < bstres) {
                            bstres = resid;
                            bstw = lambda;
                        }
                        isupmn = min(isupmn, isuppz[(2 * windex - 1) - 1]);
                        isupmx = max(isupmx, isuppz[(2 * windex) - 1]);
                        iter++;
                        //
                        //                    sin alpha <= |resid|/gap
                        //                    Note that both the residual and the gap are
                        //                    proportional to the matrix, so ||T|| doesn't play
                        //                    a role in the quotient
                        //
                        //                    Convergence test for Rayleigh-Quotient iteration
                        //                    (omitted when Bisection has been used)
                        //
                        if (resid > tol * gap && abs(rqcorr) > rqtol * abs(lambda) && !usedbs) {
                            //                       We need to check that the RQCORR update doesn't
                            //                       move the eigenvalue away from the desired one and
                            //                       towards a neighbor. -> protection with bisection
                            if (indeig <= negcnt) {
                                //                          The wanted eigenvalue lies to the left
                                sgndef = -one;
                            } else {
                                //                          The wanted eigenvalue lies to the right
                                sgndef = one;
                            }
                            //                       We only use the RQCORR if it improves the
                            //                       the iterate reasonably.
                            if ((rqcorr * sgndef >= zero) && (lambda + rqcorr <= right) && (lambda + rqcorr >= left)) {
                                usedrq = true;
                                //                          Store new midpoint of bisection interval in WORK
                                if (sgndef == one) {
                                    //                             The current LAMBDA is on the left of the true
                                    //                             eigenvalue
                                    left = lambda;
                                    //                             We prefer to assume that the error estimate
                                    //                             is correct. We could make the interval not
                                    //                             as a bracket but to be modified if the RQCORR
                                    //                             chooses to. In this case, the RIGHT side should
                                    //                             be modified as follows:
                                    //                              RIGHT = MAX(RIGHT, LAMBDA + RQCORR)
                                } else {
                                    //                             The current LAMBDA is on the right of the true
                                    //                             eigenvalue
                                    right = lambda;
                                    //                             See comment about assuming the error estimate is
                                    //                             correct above.
                                    //                              LEFT = MIN(LEFT, LAMBDA + RQCORR)
                                }
                                work[windex - 1] = half * (right + left);
                                //                          Take RQCORR since it has the correct sign and
                                //                          improves the iterate reasonably
                                lambda += rqcorr;
                                //                          Update width of error interval
                                werr[windex - 1] = half * (right - left);
                            } else {
                                needbs = true;
                            }
                            if (right - left < rqtol * abs(lambda)) {
                                //                             The eigenvalue is computed to bisection accuracy
                                //                             compute eigenvector and stop
                                usedbs = true;
                                goto statement_120;
                            } else if (iter < maxitr) {
                                goto statement_120;
                            } else if (iter == maxitr) {
                                needbs = true;
                                goto statement_120;
                            } else {
                                info = 5;
                                return;
                            }
                        } else {
                            stp2ii = false;
                            if (usedrq && usedbs && bstres <= resid) {
                                lambda = bstw;
                                stp2ii = true;
                            }
                            if (stp2ii) {
                                //                          improve error angle by second step
                                Rlar1v(in, 1, in, lambda, &d[ibegin - 1], &l[ibegin - 1], &work[(indld + ibegin - 1) - 1], &work[(indlld + ibegin - 1) - 1], pivmin, gaptol, &z[(ibegin - 1) + (windex - 1) * ldz], !usedbs, negcnt, ztz, mingma, iwork[(iindr + windex) - 1], &isuppz[(2 * windex - 1) - 1], nrminv, resid, rqcorr, &work[indwrk - 1]);
                            }
                            work[windex - 1] = lambda;
                        }
                        //
                        //                    Compute FP-vector support w.r.t. whole matrix
                        //
                        isuppz[(2 * windex - 1) - 1] += oldien;
                        isuppz[(2 * windex) - 1] += oldien;
                        zfrom = isuppz[(2 * windex - 1) - 1];
                        zto = isuppz[(2 * windex) - 1];
                        isupmn += oldien;
                        isupmx += oldien;
                        //                    Ensure vector is ok if support in the RQI has changed
                        if (isupmn < zfrom) {
                            for (ii = isupmn; ii <= zfrom - 1; ii = ii + 1) {
                                z[(ii - 1) + (windex - 1) * ldz] = zero;
                            }
                        }
                        if (isupmx > zto) {
                            for (ii = zto + 1; ii <= isupmx; ii = ii + 1) {
                                z[(ii - 1) + (windex - 1) * ldz] = zero;
                            }
                        }
                        Rscal(zto - zfrom + 1, nrminv, &z[(zfrom - 1) + (windex - 1) * ldz], 1);
                    statement_125:
                        //                    Update W
                        w[windex - 1] = lambda + sigma;
                        //                    Recompute the gaps on the left and right
                        //                    But only allow them to become larger and not
                        //                    smaller (which can only happen through "bad"
                        //                    cancellation and doesn't reflect the theory
                        //                    where the initial gaps are underestimated due
                        //                    to WERR being too crude.)
                        if (!eskip) {
                            if (k > 1) {
                                wgap[windmn - 1] = max(wgap[windmn - 1], REAL(w[windex - 1] - werr[windex - 1] - w[windmn - 1] - werr[windmn - 1]));
                            }
                            if (windex < wend) {
                                wgap[windex - 1] = max(savgap, REAL(w[windpl - 1] - werr[windpl - 1] - w[windex - 1] - werr[windex - 1]));
                            }
                        }
                        idone++;
                    }
                //                 here ends the code for the current child
                //
                statement_139:
                    //                 Proceed to any remaining child nodes
                    newfst = j + 1;
                statement_140:;
                }
            }
            ndepth++;
            goto statement_40;
        }
        ibegin = iend + 1;
        wbegin = wend + 1;
    statement_170:;
    }
    //
    //     End of Rlarrv
    //
}
