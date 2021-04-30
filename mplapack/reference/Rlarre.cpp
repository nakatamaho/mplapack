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

void Rlarre(const char *range, INTEGER const n, REAL &vl, REAL &vu, INTEGER const il, INTEGER const iu, REAL *d, REAL *e, REAL *e2, REAL const rtol1, REAL const rtol2, REAL const spltol, INTEGER nsplit, INTEGER *isplit, INTEGER &m, REAL *w, REAL *werr, REAL *wgap, INTEGER *iblock, INTEGER *indexw, REAL *gers, REAL &pivmin, REAL *work, INTEGER *iwork, INTEGER &info) {
    const INTEGER allrng = 1;
    INTEGER irange = 0;
    const INTEGER valrng = 3;
    const INTEGER indrng = 2;
    REAL safmin = 0.0;
    REAL eps = 0.0;
    REAL rtl = 0.0;
    REAL bsrtol = 0.0;
    const REAL zero = 0.0;
    REAL gl = 0.0;
    REAL gu = 0.0;
    REAL eold = 0.0;
    REAL emax = 0.0;
    INTEGER i = 0;
    REAL eabs = 0.0;
    REAL tmp1 = 0.0;
    const REAL one = 1.0;
    REAL spdiam = 0.0;
    INTEGER iinfo = 0;
    bool forceb = false;
    bool usedqd = false;
    INTEGER mm = 0;
    INTEGER ibegin = 0;
    INTEGER wbegin = 0;
    INTEGER jblk = 0;
    INTEGER iend = 0;
    INTEGER in = 0;
    INTEGER mb = 0;
    const REAL two = 2.0;
    const REAL half = one / two;
    const REAL fac = half;
    INTEGER wend = 0;
    REAL sigma = 0.0;
    INTEGER indl = 0;
    INTEGER indu = 0;
    REAL tmp = 0.0;
    const REAL hndrd = 100.0;
    REAL isleft = 0.0;
    REAL isrght = 0.0;
    const REAL four = 4.0;
    const REAL fourth = one / four;
    REAL s1 = 0.0;
    REAL s2 = 0.0;
    INTEGER cnt = 0;
    INTEGER cnt1 = 0;
    INTEGER cnt2 = 0;
    REAL sgndef = 0.0;
    REAL tau = 0.0;
    REAL clwdth = 0.0;
    REAL avgap = 0.0;
    INTEGER idum = 0;
    const INTEGER maxtry = 6;
    REAL dpivot = 0.0;
    REAL dmax = 0.0;
    INTEGER j = 0;
    const REAL maxgrowth = 64.0;
    bool norep = false;
    const REAL fudge = 2.0;
    INTEGER iseed[4];
    const REAL pert = 8.0;
    REAL rtol = 0.0;
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
    //
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n <= 0) {
        return;
    }
    //
    //     Decode RANGE
    //
    if (Mlsame(range, "A")) {
        irange = allrng;
    } else if (Mlsame(range, "V")) {
        irange = valrng;
    } else if (Mlsame(range, "I")) {
        irange = indrng;
    }
    //
    m = 0;
    //
    //     Get machine constants
    safmin = Rlamch("S");
    eps = Rlamch("P");
    //
    //     Set parameters
    rtl = sqrt(eps);
    bsrtol = sqrt(eps);
    //
    //     Treat case of 1x1 matrix for quick return
    if (n == 1) {
        if ((irange == allrng) || ((irange == valrng) && (d[1 - 1] > vl) && (d[1 - 1] <= vu)) || ((irange == indrng) && (il == 1) && (iu == 1))) {
            m = 1;
            w[1 - 1] = d[1 - 1];
            //           The computation error of the eigenvalue is zero
            werr[1 - 1] = zero;
            wgap[1 - 1] = zero;
            iblock[1 - 1] = 1;
            indexw[1 - 1] = 1;
            gers[1 - 1] = d[1 - 1];
            gers[2 - 1] = d[1 - 1];
        }
        //        store the shift for the initial RRR, which is zero in this case
        e[1 - 1] = zero;
        return;
    }
    //
    //     General case: tridiagonal matrix of order > 1
    //
    //     Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter.
    //     Compute maximum off-diagonal entry and pivmin.
    gl = d[1 - 1];
    gu = d[1 - 1];
    eold = zero;
    emax = zero;
    e[n - 1] = zero;
    for (i = 1; i <= n; i = i + 1) {
        werr[i - 1] = zero;
        wgap[i - 1] = zero;
        eabs = abs(e[i - 1]);
        if (eabs >= emax) {
            emax = eabs;
        }
        tmp1 = eabs + eold;
        gers[(2 * i - 1) - 1] = d[i - 1] - tmp1;
        gl = min(gl, gers[(2 * i - 1) - 1]);
        gers[(2 * i) - 1] = d[i - 1] + tmp1;
        gu = max(gu, gers[(2 * i) - 1]);
        eold = eabs;
    }
    //     The minimum pivot allowed in the Sturm sequence for T
    pivmin = safmin * max(one, pow2(emax));
    //     Compute spectral diameter. The Gerschgorin bounds give an
    //     estimate that is wrong by at most a factor of SQRT(2)
    spdiam = gu - gl;
    //
    //     Compute splitting points
    Rlarra(n, d, e, e2, spltol, spdiam, nsplit, isplit, iinfo);
    //
    //     Can force use of bisection instead of faster DQDS.
    //     Option left in the code for future multisection work.
    forceb = false;
    //
    //     Initialize USEDQD, DQDS should be used for ALLRNG unless someone
    //     explicitly wants bisection.
    usedqd = ((irange == allrng) && (!forceb));
    //
    if ((irange == allrng) && (!forceb)) {
        //        Set interval [VL,VU] that contains all eigenvalues
        vl = gl;
        vu = gu;
    } else {
        //        We call Rlarrd to find crude approximations to the eigenvalues
        //        in the desired range. In case IRANGE = INDRNG, we also obtain the
        //        interval (VL,VU] that contains all the wanted eigenvalues.
        //        An interval [LEFT,RIGHT] has converged if
        //        RIGHT-LEFT.LT.RTOL*MAX(ABS(LEFT),ABS(RIGHT))
        //        Rlarrd needs a WORK of size 4*N, IWORK of size 3*N
        Rlarrd(range, "B", n, vl, vu, il, iu, gers, bsrtol, d, e, e2, pivmin, nsplit, isplit, mm, w, werr, vl, vu, iblock, indexw, work, iwork, iinfo);
        if (iinfo != 0) {
            info = -1;
            return;
        }
        //        Make sure that the entries M+1 to N in W, WERR, IBLOCK, INDEXW are 0
        for (i = mm + 1; i <= n; i = i + 1) {
            w[i - 1] = zero;
            werr[i - 1] = zero;
            iblock[i - 1] = 0;
            indexw[i - 1] = 0;
        }
    }
    //
    //**
    //     Loop over unreduced blocks
    ibegin = 1;
    wbegin = 1;
    for (jblk = 1; jblk <= nsplit; jblk = jblk + 1) {
        iend = isplit[jblk - 1];
        in = iend - ibegin + 1;
        //
        //        1 X 1 block
        if (in == 1) {
            if ((irange == allrng) || ((irange == valrng) && (d[ibegin - 1] > vl) && (d[ibegin - 1] <= vu)) || ((irange == indrng) && (iblock[wbegin - 1] == jblk))) {
                m++;
                w[m - 1] = d[ibegin - 1];
                werr[m - 1] = zero;
                //              The gap for a single block doesn't matter for the later
                //              algorithm and is assigned an arbitrary large value
                wgap[m - 1] = zero;
                iblock[m - 1] = jblk;
                indexw[m - 1] = 1;
                wbegin++;
            }
            //           E( IEND ) holds the shift for the initial RRR
            e[iend - 1] = zero;
            ibegin = iend + 1;
            goto statement_170;
        }
        //
        //        Blocks of size larger than 1x1
        //
        //        E( IEND ) will hold the shift for the initial RRR, for now set it =0
        e[iend - 1] = zero;
        //
        //        Find local outer bounds GL,GU for the block
        gl = d[ibegin - 1];
        gu = d[ibegin - 1];
        for (i = ibegin; i <= iend; i = i + 1) {
            gl = min(gers[(2 * i - 1) - 1], gl);
            gu = max(gers[(2 * i) - 1], gu);
        }
        spdiam = gu - gl;
        //
        if (!((irange == allrng) && (!forceb))) {
            //           Count the number of eigenvalues in the current block.
            mb = 0;
            for (i = wbegin; i <= mm; i = i + 1) {
                if (iblock[i - 1] == jblk) {
                    mb++;
                } else {
                    goto statement_21;
                }
            }
        statement_21:
            //
            if (mb == 0) {
                //              No eigenvalue in the current block lies in the desired range
                //              E( IEND ) holds the shift for the initial RRR
                e[iend - 1] = zero;
                ibegin = iend + 1;
                goto statement_170;
            } else {
                //
                //              Decide whether dqds or bisection is more efficient
                usedqd = ((mb > fac * in) && (!forceb));
                wend = wbegin + mb - 1;
                //              Calculate gaps for the current block
                //              In later stages, when representations for individual
                //              eigenvalues are different, we use SIGMA = E( IEND ).
                sigma = zero;
                for (i = wbegin; i <= wend - 1; i = i + 1) {
                    wgap[i - 1] = max(zero, w[(i + 1) - 1] - werr[(i + 1) - 1] - (w[i - 1] + werr[i - 1]));
                }
                wgap[wend - 1] = max(zero, vu - sigma - (w[wend - 1] + werr[wend - 1]));
                //              Find local index of the first and last desired evalue.
                indl = indexw[wbegin - 1];
                indu = indexw[wend - 1];
            }
        }
        if (((irange == allrng) && (!forceb)) || usedqd) {
            //           Case of DQDS
            //           Find approximations to the extremal eigenvalues of the block
            Rlarrk(in, 1, gl, gu, &d[ibegin - 1], &e2[ibegin - 1], pivmin, rtl, tmp, tmp1, iinfo);
            if (iinfo != 0) {
                info = -1;
                return;
            }
            isleft = max(gl, tmp - tmp1 - hndrd * eps * abs(tmp - tmp1));
            //
            Rlarrk(in, in, gl, gu, &d[ibegin - 1], &e2[ibegin - 1], pivmin, rtl, tmp, tmp1, iinfo);
            if (iinfo != 0) {
                info = -1;
                return;
            }
            isrght = min(gu, tmp + tmp1 + hndrd * eps * abs(tmp + tmp1));
            //           Improve the estimate of the spectral diameter
            spdiam = isrght - isleft;
        } else {
            //           Case of bisection
            //           Find approximations to the wanted extremal eigenvalues
            isleft = max(gl, w[wbegin - 1] - werr[wbegin - 1] - hndrd * eps * abs(w[wbegin - 1] - werr[wbegin - 1]));
            isrght = min(gu, w[wend - 1] + werr[wend - 1] + hndrd * eps * abs(w[wend - 1] + werr[wend - 1]));
        }
        //
        //        Decide whether the base representation for the current block
        //        L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I
        //        should be on the left or the right end of the current block.
        //        The strategy is to shift to the end which is "more populated"
        //        Furthermore, decide whether to use DQDS for the computation of
        //        the eigenvalue approximations at the end of Rlarre or bisection.
        //        dqds is chosen if all eigenvalues are desired or the number of
        //        eigenvalues to be computed is large compared to the blocksize.
        if ((irange == allrng) && (!forceb)) {
            //           If all the eigenvalues have to be computed, we use dqd
            usedqd = true;
            //           INDL is the local index of the first eigenvalue to compute
            indl = 1;
            indu = in;
            //           MB =  number of eigenvalues to compute
            mb = in;
            wend = wbegin + mb - 1;
            //           Define 1/4 and 3/4 points of the spectrum
            s1 = isleft + fourth * spdiam;
            s2 = isrght - fourth * spdiam;
        } else {
            //           Rlarrd has computed IBLOCK and INDEXW for each eigenvalue
            //           approximation.
            //           choose sigma
            if (usedqd) {
                s1 = isleft + fourth * spdiam;
                s2 = isrght - fourth * spdiam;
            } else {
                tmp = min(isrght, vu) - max(isleft, vl);
                s1 = max(isleft, vl) + fourth * tmp;
                s2 = min(isrght, vu) - fourth * tmp;
            }
        }
        //
        //        Compute the negcount at the 1/4 and 3/4 points
        if (mb > 1) {
            Rlarrc("T", in, s1, s2, &d[ibegin - 1], &e[ibegin - 1], pivmin, cnt, cnt1, cnt2, iinfo);
        }
        //
        if (mb == 1) {
            sigma = gl;
            sgndef = one;
        } else if (cnt1 - indl >= indu - cnt2) {
            if ((irange == allrng) && (!forceb)) {
                sigma = max(isleft, gl);
            } else if (usedqd) {
                //              use Gerschgorin bound as shift to get pos def matrix
                //              for dqds
                sigma = isleft;
            } else {
                //              use approximation of the first desired eigenvalue of the
                //              block as shift
                sigma = max(isleft, vl);
            }
            sgndef = one;
        } else {
            if ((irange == allrng) && (!forceb)) {
                sigma = min(isrght, gu);
            } else if (usedqd) {
                //              use Gerschgorin bound as shift to get neg def matrix
                //              for dqds
                sigma = isrght;
            } else {
                //              use approximation of the first desired eigenvalue of the
                //              block as shift
                sigma = min(isrght, vu);
            }
            sgndef = -one;
        }
        //
        //        An initial SIGMA has been chosen that will be used for computing
        //        T - SIGMA I = L D L^T
        //        Define the increment TAU of the shift in case the initial shift
        //        needs to be refined to obtain a factorization with not too much
        //        element growth.
        if (usedqd) {
            //           The initial SIGMA was to the outer end of the spectrum
            //           the matrix is definite and we need not retreat.
            tau = spdiam * eps * n + two * pivmin;
            tau = max(tau, two * eps * abs(sigma));
        } else {
            if (mb > 1) {
                clwdth = w[wend - 1] + werr[wend - 1] - w[wbegin - 1] - werr[wbegin - 1];
                avgap = abs(clwdth / castREAL(wend - wbegin));
                if (sgndef == one) {
                    tau = half * max(wgap[wbegin - 1], avgap);
                    tau = max(tau, werr[wbegin - 1]);
                } else {
                    tau = half * max(wgap[(wend - 1) - 1], avgap);
                    tau = max(tau, werr[wend - 1]);
                }
            } else {
                tau = werr[wbegin - 1];
            }
        }
        //
        for (idum = 1; idum <= maxtry; idum = idum + 1) {
            //           Compute L D L^T factorization of tridiagonal matrix T - sigma I.
            //           Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of
            //           pivots in WORK(2*IN+1:3*IN)
            dpivot = d[ibegin - 1] - sigma;
            work[1 - 1] = dpivot;
            dmax = abs(work[1 - 1]);
            j = ibegin;
            for (i = 1; i <= in - 1; i = i + 1) {
                work[(2 * in + i) - 1] = one / work[i - 1];
                tmp = e[j - 1] * work[(2 * in + i) - 1];
                work[(in + i) - 1] = tmp;
                dpivot = (d[(j + 1) - 1] - sigma) - tmp * e[j - 1];
                work[(i + 1) - 1] = dpivot;
                dmax = max(dmax, abs(dpivot));
                j++;
            }
            //           check for element growth
            if (dmax > maxgrowth * spdiam) {
                norep = true;
            } else {
                norep = false;
            }
            if (usedqd && !norep) {
                //              Ensure the definiteness of the representation
                //              All entries of D (of L D L^T) must have the same sign
                for (i = 1; i <= in; i = i + 1) {
                    tmp = sgndef * work[i - 1];
                    if (tmp < zero) {
                        norep = true;
                    }
                }
            }
            if (norep) {
                //              Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin
                //              shift which makes the matrix definite. So we should end up
                //              here really only in the case of IRANGE = VALRNG or INDRNG.
                if (idum == maxtry - 1) {
                    if (sgndef == one) {
                        //                    The fudged Gerschgorin shift should succeed
                        sigma = gl - fudge * spdiam * eps * n - fudge * two * pivmin;
                    } else {
                        sigma = gu + fudge * spdiam * eps * n + fudge * two * pivmin;
                    }
                } else {
                    sigma = sigma - sgndef * tau;
                    tau = two * tau;
                }
            } else {
                //              an initial RRR is found
                goto statement_83;
            }
        }
        //        if the program reaches this point, no base representation could be
        //        found in MAXTRY iterations.
        info = 2;
        return;
    //
    statement_83:
        //        At this point, we have found an initial base representation
        //        T - SIGMA I = L D L^T with not too much element growth.
        //        Store the shift.
        e[iend - 1] = sigma;
        //        Store D and L.
        Rcopy(in, work, 1, &d[ibegin - 1], 1);
        Rcopy(in - 1, &work[(in + 1) - 1], 1, &e[ibegin - 1], 1);
        //
        if (mb > 1) {
            //
            //           Perturb each entry of the base representation by a small
            //           (but random) relative amount to overcome difficulties with
            //           glued matrices.
            //
            for (i = 1; i <= 4; i = i + 1) {
                iseed[i - 1] = 1;
            }
            //
            Rlarnv(2, iseed, 2 * in - 1, &work[1 - 1]);
            for (i = 1; i <= in - 1; i = i + 1) {
                d[(ibegin + i - 1) - 1] = d[(ibegin + i - 1) - 1] * (one + eps * pert * work[i - 1]);
                e[(ibegin + i - 1) - 1] = e[(ibegin + i - 1) - 1] * (one + eps * pert * work[(in + i) - 1]);
            }
            d[iend - 1] = d[iend - 1] * (one + eps * four * work[in - 1]);
            //
        }
        //
        //        Don't update the Gerschgorin intervals because keeping track
        //        of the updates would be too much work in Rlarrv.
        //        We update W instead and use it to locate the proper Gerschgorin
        //        intervals.
        //
        //        Compute the required eigenvalues of L D L' by bisection or dqds
        if (!usedqd) {
            //           If Rlarrd has been used, shift the eigenvalue approximations
            //           according to their representation. This is necessary for
            //           a uniform Rlarrv since dqds computes eigenvalues of the
            //           shifted representation. In Rlarrv, W will always hold the
            //           UNshifted eigenvalue approximation.
            for (j = wbegin; j <= wend; j = j + 1) {
                w[j - 1] = w[j - 1] - sigma;
                werr[j - 1] += abs(w[j - 1]) * eps;
            }
            //           call Rlarrb to reduce eigenvalue error of the approximations
            //           from Rlarrd
            for (i = ibegin; i <= iend - 1; i = i + 1) {
                work[i - 1] = d[i - 1] * pow2(e[i - 1]);
            }
            //           use bisection to find EV from INDL to INDU
            Rlarrb(in, &d[ibegin - 1], &work[ibegin - 1], indl, indu, rtol1, rtol2, indl - 1, &w[wbegin - 1], &wgap[wbegin - 1], &werr[wbegin - 1], &work[(2 * n + 1) - 1], iwork, pivmin, spdiam, in, iinfo);
            if (iinfo != 0) {
                info = -4;
                return;
            }
            //           Rlarrb computes all gaps correctly except for the last one
            //           Record distance to VU/GU
            wgap[wend - 1] = max(zero, (vu - sigma) - (w[wend - 1] + werr[wend - 1]));
            for (i = indl; i <= indu; i = i + 1) {
                m++;
                iblock[m - 1] = jblk;
                indexw[m - 1] = i;
            }
        } else {
            //           Call dqds to get all eigs (and then possibly delete unwanted
            //           eigenvalues).
            //           Note that dqds finds the eigenvalues of the L D L^T representation
            //           of T to high relative accuracy. High relative accuracy
            //           might be lost when the shift of the RRR is subtracted to obtain
            //           the eigenvalues of T. However, T is not guaranteed to define its
            //           eigenvalues to high relative accuracy anyway.
            //           Set RTOL to the order of the tolerance used in Rlasq2
            //           This is an ESTIMATED error, the worst case bound is 4*N*EPS
            //           which is usually too large and requires unnecessary work to be
            //           done by bisection when computing the eigenvectors
            rtol = log(castREAL(in)) * four * eps;
            j = ibegin;
            for (i = 1; i <= in - 1; i = i + 1) {
                work[(2 * i - 1) - 1] = abs(d[j - 1]);
                work[(2 * i) - 1] = e[j - 1] * e[j - 1] * work[(2 * i - 1) - 1];
                j++;
            }
            work[(2 * in - 1) - 1] = abs(d[iend - 1]);
            work[(2 * in) - 1] = zero;
            Rlasq2(in, work, iinfo);
            if (iinfo != 0) {
                //              If IINFO = -5 then an index is part of a tight cluster
                //              and should be changed. The index is in IWORK(1) and the
                //              gap is in WORK(N+1)
                info = -5;
                return;
            } else {
                //              Test that all eigenvalues are positive as expected
                for (i = 1; i <= in; i = i + 1) {
                    if (work[i - 1] < zero) {
                        info = -6;
                        return;
                    }
                }
            }
            if (sgndef > zero) {
                for (i = indl; i <= indu; i = i + 1) {
                    m++;
                    w[m - 1] = work[(in - i + 1) - 1];
                    iblock[m - 1] = jblk;
                    indexw[m - 1] = i;
                }
            } else {
                for (i = indl; i <= indu; i = i + 1) {
                    m++;
                    w[m - 1] = -work[i - 1];
                    iblock[m - 1] = jblk;
                    indexw[m - 1] = i;
                }
            }
            //
            for (i = m - mb + 1; i <= m; i = i + 1) {
                //              the value of RTOL below should be the tolerance in Rlasq2
                werr[i - 1] = rtol * abs(w[i - 1]);
            }
            for (i = m - mb + 1; i <= m - 1; i = i + 1) {
                //              compute the right gap between the intervals
                wgap[i - 1] = max(zero, w[(i + 1) - 1] - werr[(i + 1) - 1] - (w[i - 1] + werr[i - 1]));
            }
            wgap[m - 1] = max(zero, (vu - sigma) - (w[m - 1] + werr[m - 1]));
        }
        //        proceed with next block
        ibegin = iend + 1;
        wbegin = wend + 1;
    statement_170:;
    }
    //
    //     end of Rlarre
    //
}
