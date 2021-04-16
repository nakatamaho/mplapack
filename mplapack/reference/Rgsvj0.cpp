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

void Rgsvj0(const char *jobv, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *d, REAL *sva, INTEGER const mv, REAL *v, INTEGER const ldv, REAL const eps, REAL const sfmin, REAL const tol, INTEGER const nsweep, REAL *work, INTEGER const lwork, INTEGER &info) {
    bool applv = false;
    bool rsvec = false;
    INTEGER mvl = 0;
    REAL rooteps = 0.0;
    REAL rootsfmin = 0.0;
    REAL small = 0.0;
    const REAL one = 1.0;
    REAL big = 0.0;
    REAL rootbig = 0.0;
    REAL bigtheta = 0.0;
    REAL roottol = 0.0;
    INTEGER emptsw = 0;
    INTEGER notrot = 0;
    const REAL zero = 0.0;
    REAL fastr[5];
    INTEGER swband = 0;
    INTEGER kbl = 0;
    INTEGER nbl = 0;
    INTEGER blskip = 0;
    INTEGER rowskip = 0;
    INTEGER lkahead = 0;
    INTEGER pskipped = 0;
    INTEGER i = 0;
    REAL mxaapq = 0.0;
    REAL mxsinj = 0.0;
    INTEGER iswrot = 0;
    INTEGER ibr = 0;
    INTEGER igl = 0;
    INTEGER ir1 = 0;
    INTEGER p = 0;
    INTEGER q = 0;
    REAL temp1 = 0.0;
    REAL aapp = 0.0;
    REAL aaqq = 0.0;
    REAL aapp0 = 0.0;
    bool rotok = false;
    REAL aapq = 0.0;
    INTEGER ierr = 0;
    REAL aqoap = 0.0;
    REAL apoaq = 0.0;
    const REAL half = 0.5e0;
    REAL theta = 0.0;
    REAL t = 0.0;
    REAL thsign = 0.0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    INTEGER jbc = 0;
    INTEGER jgl = 0;
    INTEGER ijblsk = 0;
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
    //     .. Local Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    applv = Mlsame(jobv, "A");
    rsvec = Mlsame(jobv, "V");
    if (!(rsvec || applv || Mlsame(jobv, "N"))) {
        info = -1;
    } else if (m < 0) {
        info = -2;
    } else if ((n < 0) || (n > m)) {
        info = -3;
    } else if (lda < m) {
        info = -5;
    } else if ((rsvec || applv) && (mv < 0)) {
        info = -8;
    } else if ((rsvec && (ldv < n)) || (applv && (ldv < mv))) {
        info = -10;
    } else if (tol <= eps) {
        info = -13;
    } else if (nsweep < 0) {
        info = -14;
    } else if (lwork < m) {
        info = -16;
    } else {
        info = 0;
    }
    //
    //     #:(
    if (info != 0) {
        Mxerbla("Rgsvj0", -info);
        return;
    }
    //
    if (rsvec) {
        mvl = n;
    } else if (applv) {
        mvl = mv;
    }
    rsvec = rsvec || applv;
    //
    rooteps = sqrt(eps);
    rootsfmin = sqrt(sfmin);
    small = sfmin / eps;
    big = one / sfmin;
    rootbig = one / rootsfmin;
    bigtheta = one / rooteps;
    roottol = sqrt(tol);
    //
    //     -#- Row-cyclic Jacobi SVD algorithm with column pivoting -#-
    //
    emptsw = (n * (n - 1)) / 2;
    notrot = 0;
    fastr[1 - 1] = zero;
    //
    //     -#- Row-cyclic pivot strategy with de Rijk's pivoting -#-
    //
    swband = 0;
    //[TP] SWBAND is a tuning parameter. It is meaningful and effective
    //     if SGESVJ is used as a computational routine in the preconditioned
    //     Jacobi SVD algorithm SGESVJ. For sweeps i=1:SWBAND the procedure
    //     ......
    //
    kbl = min(8, n);
    //[TP] KBL is a tuning parameter that defines the tile size in the
    //     tiling of the p-q loops of pivot pairs. In general, an optimal
    //     parameters of the computer's memory.
    //
    nbl = n / kbl;
    if ((nbl * kbl) != n) {
        nbl++;
    }
    //
    blskip = (pow2(kbl)) + 1;
    //[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.
    //
    rowskip = min(5, kbl);
    //[TP] ROWSKIP is a tuning parameter.
    //
    lkahead = 1;
    //[TP] LKAHEAD is a tuning parameter.
    swband = 0;
    pskipped = 0;
    //
    for (i = 1; i <= nsweep; i = i + 1) {
        //     .. go go go ...
        //
        mxaapq = zero;
        mxsinj = zero;
        iswrot = 0;
        //
        notrot = 0;
        pskipped = 0;
        //
        for (ibr = 1; ibr <= nbl; ibr = ibr + 1) {
            //
            igl = (ibr - 1) * kbl + 1;
            //
            for (ir1 = 0; ir1 <= min(lkahead, nbl - ibr); ir1 = ir1 + 1) {
                //
                igl += ir1 * kbl;
                //
                for (p = igl; p <= min(igl + kbl - 1, n - 1); p = p + 1) {
                    //
                    //     .. de Rijk's pivoting
                    q = iRamax(n - p + 1, &sva[p - 1], 1) + p - 1;
                    if (p != q) {
                        Rswap(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
                        if (rsvec) {
                            Rswap(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
                        }
                        temp1 = sva[p - 1];
                        sva[p - 1] = sva[q - 1];
                        sva[q - 1] = temp1;
                        temp1 = d[p - 1];
                        d[p - 1] = d[q - 1];
                        d[q - 1] = temp1;
                    }
                    //
                    if (ir1 == 0) {
                        //
                        //        Column norms are periodically updated by explicit
                        //        norm computation.
                        //        Caveat:
                        //        Some BLAS implementations compute Rnrm2(M,A(1,p),1)
                        //        as DSQRT(Rdot(M,A(1,p),1,A(1,p),1)), which may result in
                        //        overflow for ||A(:,p)||_2 > DSQRT(overflow_threshold), and
                        //        underflow for ||A(:,p)||_2 < DSQRT(underflow_threshold).
                        //        Hence, Rnrm2 cannot be trusted, not even in the case when
                        //        the true norm is far from the under(over)flow boundaries.
                        //        If properly implemented Rnrm2 is available, the IF-THEN-ELSE
                        //        below should read "AAPP = Rnrm2( M, A(1,p), 1 ) * D(p)".
                        //
                        if ((sva[p - 1] < rootbig) && (sva[p - 1] > rootsfmin)) {
                            sva[p - 1] = Rnrm2(m, &a[(p - 1) * lda], 1) * d[p - 1];
                        } else {
                            temp1 = zero;
                            aapp = one;
                            Rlassq(m, &a[(p - 1) * lda], 1, temp1, aapp);
                            sva[p - 1] = temp1 * sqrt(aapp) * d[p - 1];
                        }
                        aapp = sva[p - 1];
                    } else {
                        aapp = sva[p - 1];
                    }
                    //
                    if (aapp > zero) {
                        //
                        pskipped = 0;
                        //
                        for (q = p + 1; q <= min(igl + kbl - 1, n); q = q + 1) {
                            //
                            aaqq = sva[q - 1];
                            //
                            if (aaqq > zero) {
                                //
                                aapp0 = aapp;
                                if (aaqq >= one) {
                                    rotok = (small * aapp) <= aaqq;
                                    if (aapp < (big / aaqq)) {
                                        aapq = (Rdot(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1) * d[p - 1] * d[q - 1] / aaqq) / aapp;
                                    } else {
                                        Rcopy(m, &a[(p - 1) * lda], 1, work, 1);
                                        Rlascl("G", 0, 0, aapp, d[p - 1], m, 1, work, lda, ierr);
                                        aapq = Rdot(m, work, 1, &a[(q - 1) * lda], 1) * d[q - 1] / aaqq;
                                    }
                                } else {
                                    rotok = aapp <= (aaqq / small);
                                    if (aapp > (small / aaqq)) {
                                        aapq = (Rdot(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1) * d[p - 1] * d[q - 1] / aaqq) / aapp;
                                    } else {
                                        Rcopy(m, &a[(q - 1) * lda], 1, work, 1);
                                        Rlascl("G", 0, 0, aaqq, d[q - 1], m, 1, work, lda, ierr);
                                        aapq = Rdot(m, work, 1, &a[(p - 1) * lda], 1) * d[p - 1] / aapp;
                                    }
                                }
                                //
                                mxaapq = max(mxaapq, abs(aapq));
                                //
                                //        TO rotate or NOT to rotate, THAT is the question ...
                                //
                                if (abs(aapq) > tol) {
                                    //
                                    //           .. rotate
                                    //           ROTATED = ROTATED + ONE
                                    //
                                    if (ir1 == 0) {
                                        notrot = 0;
                                        pskipped = 0;
                                        iswrot++;
                                    }
                                    //
                                    if (rotok) {
                                        //
                                        aqoap = aaqq / aapp;
                                        apoaq = aapp / aaqq;
                                        theta = -half * abs(aqoap - apoaq) / aapq;
                                        //
                                        if (abs(theta) > bigtheta) {
                                            //
                                            t = half / theta;
                                            fastr[3 - 1] = t * d[p - 1] / d[q - 1];
                                            fastr[4 - 1] = -t * d[q - 1] / d[p - 1];
                                            Rrotm(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1, fastr);
                                            if (rsvec) {
                                                Rrotm(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1, fastr);
                                            }
                                            sva[q - 1] = aaqq * sqrt(max(zero, one + t * apoaq * aapq));
                                            aapp = aapp * sqrt(max(zero, one - t * aqoap * aapq));
                                            mxsinj = max(mxsinj, abs(t));
                                            //
                                        } else {
                                            //
                                            //                 .. choose correct signum for THETA and rotate
                                            //
                                            thsign = -sign(one, aapq);
                                            t = one / (theta + thsign * sqrt(one + theta * theta));
                                            cs = sqrt(one / (one + t * t));
                                            sn = t * cs;
                                            //
                                            mxsinj = max(mxsinj, abs(sn));
                                            sva[q - 1] = aaqq * sqrt(max(zero, one + t * apoaq * aapq));
                                            aapp = aapp * sqrt(max(zero, one - t * aqoap * aapq));
                                            //
                                            apoaq = d[p - 1] / d[q - 1];
                                            aqoap = d[q - 1] / d[p - 1];
                                            if (d[p - 1] >= one) {
                                                if (d[q - 1] >= one) {
                                                    fastr[3 - 1] = t * apoaq;
                                                    fastr[4 - 1] = -t * aqoap;
                                                    d[p - 1] = d[p - 1] * cs;
                                                    d[q - 1] = d[q - 1] * cs;
                                                    Rrotm(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1, fastr);
                                                    if (rsvec) {
                                                        Rrotm(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1, fastr);
                                                    }
                                                } else {
                                                    Raxpy(m, -t * aqoap, &a[(q - 1) * lda], 1, &a[(p - 1) * lda], 1);
                                                    Raxpy(m, cs * sn * apoaq, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
                                                    d[p - 1] = d[p - 1] * cs;
                                                    d[q - 1] = d[q - 1] / cs;
                                                    if (rsvec) {
                                                        Raxpy(mvl, -t * aqoap, &v[(q - 1) * ldv], 1, &v[(p - 1) * ldv], 1);
                                                        Raxpy(mvl, cs * sn * apoaq, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
                                                    }
                                                }
                                            } else {
                                                if (d[q - 1] >= one) {
                                                    Raxpy(m, t * apoaq, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
                                                    Raxpy(m, -cs * sn * aqoap, &a[(q - 1) * lda], 1, &a[(p - 1) * lda], 1);
                                                    d[p - 1] = d[p - 1] / cs;
                                                    d[q - 1] = d[q - 1] * cs;
                                                    if (rsvec) {
                                                        Raxpy(mvl, t * apoaq, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
                                                        Raxpy(mvl, -cs * sn * aqoap, &v[(q - 1) * ldv], 1, &v[(p - 1) * ldv], 1);
                                                    }
                                                } else {
                                                    if (d[p - 1] >= d[q - 1]) {
                                                        Raxpy(m, -t * aqoap, &a[(q - 1) * lda], 1, &a[(p - 1) * lda], 1);
                                                        Raxpy(m, cs * sn * apoaq, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
                                                        d[p - 1] = d[p - 1] * cs;
                                                        d[q - 1] = d[q - 1] / cs;
                                                        if (rsvec) {
                                                            Raxpy(mvl, -t * aqoap, &v[(q - 1) * ldv], 1, &v[(p - 1) * ldv], 1);
                                                            Raxpy(mvl, cs * sn * apoaq, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
                                                        }
                                                    } else {
                                                        Raxpy(m, t * apoaq, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
                                                        Raxpy(m, -cs * sn * aqoap, &a[(q - 1) * lda], 1, &a[(p - 1) * lda], 1);
                                                        d[p - 1] = d[p - 1] / cs;
                                                        d[q - 1] = d[q - 1] * cs;
                                                        if (rsvec) {
                                                            Raxpy(mvl, t * apoaq, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
                                                            Raxpy(mvl, -cs * sn * aqoap, &v[(q - 1) * ldv], 1, &v[(p - 1) * ldv], 1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        //
                                    } else {
                                        //              .. have to use modified Gram-Schmidt like transformation
                                        Rcopy(m, &a[(p - 1) * lda], 1, work, 1);
                                        Rlascl("G", 0, 0, aapp, one, m, 1, work, lda, ierr);
                                        Rlascl("G", 0, 0, aaqq, one, m, 1, &a[(q - 1) * lda], lda, ierr);
                                        temp1 = -aapq * d[p - 1] / d[q - 1];
                                        Raxpy(m, temp1, work, 1, &a[(q - 1) * lda], 1);
                                        Rlascl("G", 0, 0, one, aaqq, m, 1, &a[(q - 1) * lda], lda, ierr);
                                        sva[q - 1] = aaqq * sqrt(max(zero, one - aapq * aapq));
                                        mxsinj = max(mxsinj, sfmin);
                                    }
                                    //           END IF ROTOK THEN ... ELSE
                                    //
                                    //           In the case of cancellation in updating SVA(q), SVA(p)
                                    //           recompute SVA(q), SVA(p).
                                    if (pow2((sva[q - 1] / aaqq)) <= rooteps) {
                                        if ((aaqq < rootbig) && (aaqq > rootsfmin)) {
                                            sva[q - 1] = Rnrm2(m, &a[(q - 1) * lda], 1) * d[q - 1];
                                        } else {
                                            t = zero;
                                            aaqq = one;
                                            Rlassq(m, &a[(q - 1) * lda], 1, t, aaqq);
                                            sva[q - 1] = t * sqrt(aaqq) * d[q - 1];
                                        }
                                    }
                                    if ((aapp / aapp0) <= rooteps) {
                                        if ((aapp < rootbig) && (aapp > rootsfmin)) {
                                            aapp = Rnrm2(m, &a[(p - 1) * lda], 1) * d[p - 1];
                                        } else {
                                            t = zero;
                                            aapp = one;
                                            Rlassq(m, &a[(p - 1) * lda], 1, t, aapp);
                                            aapp = t * sqrt(aapp) * d[p - 1];
                                        }
                                        sva[p - 1] = aapp;
                                    }
                                    //
                                } else {
                                    //        A(:,p) and A(:,q) already numerically orthogonal
                                    if (ir1 == 0) {
                                        notrot++;
                                    }
                                    pskipped++;
                                }
                            } else {
                                //        A(:,q) is zero column
                                if (ir1 == 0) {
                                    notrot++;
                                }
                                pskipped++;
                            }
                            //
                            if ((i <= swband) && (pskipped > rowskip)) {
                                if (ir1 == 0) {
                                    aapp = -aapp;
                                }
                                notrot = 0;
                                goto statement_2103;
                            }
                            //
                        }
                    //     END q-LOOP
                    //
                    statement_2103:
                        //     bailed out of q-loop
                        //
                        sva[p - 1] = aapp;
                        //
                    } else {
                        sva[p - 1] = aapp;
                        if ((ir1 == 0) && (aapp == zero)) {
                            notrot += min(igl + kbl - 1, n) - p;
                        }
                    }
                    //
                }
                //     end of the p-loop
                //     end of doing the block ( ibr, ibr )
            }
            //     end of ir1-loop
            //
            //........................................................
            // ... go to the off diagonal blocks
            //
            igl = (ibr - 1) * kbl + 1;
            //
            for (jbc = ibr + 1; jbc <= nbl; jbc = jbc + 1) {
                //
                jgl = (jbc - 1) * kbl + 1;
                //
                //        doing the block at ( ibr, jbc )
                //
                ijblsk = 0;
                for (p = igl; p <= min(igl + kbl - 1, n); p = p + 1) {
                    //
                    aapp = sva[p - 1];
                    //
                    if (aapp > zero) {
                        //
                        pskipped = 0;
                        //
                        for (q = jgl; q <= min(jgl + kbl - 1, n); q = q + 1) {
                            //
                            aaqq = sva[q - 1];
                            //
                            if (aaqq > zero) {
                                aapp0 = aapp;
                                //
                                //     -#- M x 2 Jacobi SVD -#-
                                //
                                //        -#- Safe Gram matrix computation -#-
                                //
                                if (aaqq >= one) {
                                    if (aapp >= aaqq) {
                                        rotok = (small * aapp) <= aaqq;
                                    } else {
                                        rotok = (small * aaqq) <= aapp;
                                    }
                                    if (aapp < (big / aaqq)) {
                                        aapq = (Rdot(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1) * d[p - 1] * d[q - 1] / aaqq) / aapp;
                                    } else {
                                        Rcopy(m, &a[(p - 1) * lda], 1, work, 1);
                                        Rlascl("G", 0, 0, aapp, d[p - 1], m, 1, work, lda, ierr);
                                        aapq = Rdot(m, work, 1, &a[(q - 1) * lda], 1) * d[q - 1] / aaqq;
                                    }
                                } else {
                                    if (aapp >= aaqq) {
                                        rotok = aapp <= (aaqq / small);
                                    } else {
                                        rotok = aaqq <= (aapp / small);
                                    }
                                    if (aapp > (small / aaqq)) {
                                        aapq = (Rdot(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1) * d[p - 1] * d[q - 1] / aaqq) / aapp;
                                    } else {
                                        Rcopy(m, &a[(q - 1) * lda], 1, work, 1);
                                        Rlascl("G", 0, 0, aaqq, d[q - 1], m, 1, work, lda, ierr);
                                        aapq = Rdot(m, work, 1, &a[(p - 1) * lda], 1) * d[p - 1] / aapp;
                                    }
                                }
                                //
                                mxaapq = max(mxaapq, abs(aapq));
                                //
                                //        TO rotate or NOT to rotate, THAT is the question ...
                                //
                                if (abs(aapq) > tol) {
                                    notrot = 0;
                                    //           ROTATED  = ROTATED + 1
                                    pskipped = 0;
                                    iswrot++;
                                    //
                                    if (rotok) {
                                        //
                                        aqoap = aaqq / aapp;
                                        apoaq = aapp / aaqq;
                                        theta = -half * abs(aqoap - apoaq) / aapq;
                                        if (aaqq > aapp0) {
                                            theta = -theta;
                                        }
                                        //
                                        if (abs(theta) > bigtheta) {
                                            t = half / theta;
                                            fastr[3 - 1] = t * d[p - 1] / d[q - 1];
                                            fastr[4 - 1] = -t * d[q - 1] / d[p - 1];
                                            Rrotm(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1, fastr);
                                            if (rsvec) {
                                                Rrotm(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1, fastr);
                                            }
                                            sva[q - 1] = aaqq * sqrt(max(zero, one + t * apoaq * aapq));
                                            aapp = aapp * sqrt(max(zero, one - t * aqoap * aapq));
                                            mxsinj = max(mxsinj, abs(t));
                                        } else {
                                            //
                                            //                 .. choose correct signum for THETA and rotate
                                            //
                                            thsign = -sign(one, aapq);
                                            if (aaqq > aapp0) {
                                                thsign = -thsign;
                                            }
                                            t = one / (theta + thsign * sqrt(one + theta * theta));
                                            cs = sqrt(one / (one + t * t));
                                            sn = t * cs;
                                            mxsinj = max(mxsinj, abs(sn));
                                            sva[q - 1] = aaqq * sqrt(max(zero, one + t * apoaq * aapq));
                                            aapp = aapp * sqrt(max(zero, one - t * aqoap * aapq));
                                            //
                                            apoaq = d[p - 1] / d[q - 1];
                                            aqoap = d[q - 1] / d[p - 1];
                                            if (d[p - 1] >= one) {
                                                //
                                                if (d[q - 1] >= one) {
                                                    fastr[3 - 1] = t * apoaq;
                                                    fastr[4 - 1] = -t * aqoap;
                                                    d[p - 1] = d[p - 1] * cs;
                                                    d[q - 1] = d[q - 1] * cs;
                                                    Rrotm(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1, fastr);
                                                    if (rsvec) {
                                                        Rrotm(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1, fastr);
                                                    }
                                                } else {
                                                    Raxpy(m, -t * aqoap, &a[(q - 1) * lda], 1, &a[(p - 1) * lda], 1);
                                                    Raxpy(m, cs * sn * apoaq, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
                                                    if (rsvec) {
                                                        Raxpy(mvl, -t * aqoap, &v[(q - 1) * ldv], 1, &v[(p - 1) * ldv], 1);
                                                        Raxpy(mvl, cs * sn * apoaq, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
                                                    }
                                                    d[p - 1] = d[p - 1] * cs;
                                                    d[q - 1] = d[q - 1] / cs;
                                                }
                                            } else {
                                                if (d[q - 1] >= one) {
                                                    Raxpy(m, t * apoaq, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
                                                    Raxpy(m, -cs * sn * aqoap, &a[(q - 1) * lda], 1, &a[(p - 1) * lda], 1);
                                                    if (rsvec) {
                                                        Raxpy(mvl, t * apoaq, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
                                                        Raxpy(mvl, -cs * sn * aqoap, &v[(q - 1) * ldv], 1, &v[(p - 1) * ldv], 1);
                                                    }
                                                    d[p - 1] = d[p - 1] / cs;
                                                    d[q - 1] = d[q - 1] * cs;
                                                } else {
                                                    if (d[p - 1] >= d[q - 1]) {
                                                        Raxpy(m, -t * aqoap, &a[(q - 1) * lda], 1, &a[(p - 1) * lda], 1);
                                                        Raxpy(m, cs * sn * apoaq, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
                                                        d[p - 1] = d[p - 1] * cs;
                                                        d[q - 1] = d[q - 1] / cs;
                                                        if (rsvec) {
                                                            Raxpy(mvl, -t * aqoap, &v[(q - 1) * ldv], 1, &v[(p - 1) * ldv], 1);
                                                            Raxpy(mvl, cs * sn * apoaq, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
                                                        }
                                                    } else {
                                                        Raxpy(m, t * apoaq, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
                                                        Raxpy(m, -cs * sn * aqoap, &a[(q - 1) * lda], 1, &a[(p - 1) * lda], 1);
                                                        d[p - 1] = d[p - 1] / cs;
                                                        d[q - 1] = d[q - 1] * cs;
                                                        if (rsvec) {
                                                            Raxpy(mvl, t * apoaq, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
                                                            Raxpy(mvl, -cs * sn * aqoap, &v[(q - 1) * ldv], 1, &v[(p - 1) * ldv], 1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        //
                                    } else {
                                        if (aapp > aaqq) {
                                            Rcopy(m, &a[(p - 1) * lda], 1, work, 1);
                                            Rlascl("G", 0, 0, aapp, one, m, 1, work, lda, ierr);
                                            Rlascl("G", 0, 0, aaqq, one, m, 1, &a[(q - 1) * lda], lda, ierr);
                                            temp1 = -aapq * d[p - 1] / d[q - 1];
                                            Raxpy(m, temp1, work, 1, &a[(q - 1) * lda], 1);
                                            Rlascl("G", 0, 0, one, aaqq, m, 1, &a[(q - 1) * lda], lda, ierr);
                                            sva[q - 1] = aaqq * sqrt(max(zero, one - aapq * aapq));
                                            mxsinj = max(mxsinj, sfmin);
                                        } else {
                                            Rcopy(m, &a[(q - 1) * lda], 1, work, 1);
                                            Rlascl("G", 0, 0, aaqq, one, m, 1, work, lda, ierr);
                                            Rlascl("G", 0, 0, aapp, one, m, 1, &a[(p - 1) * lda], lda, ierr);
                                            temp1 = -aapq * d[q - 1] / d[p - 1];
                                            Raxpy(m, temp1, work, 1, &a[(p - 1) * lda], 1);
                                            Rlascl("G", 0, 0, one, aapp, m, 1, &a[(p - 1) * lda], lda, ierr);
                                            sva[p - 1] = aapp * sqrt(max(zero, one - aapq * aapq));
                                            mxsinj = max(mxsinj, sfmin);
                                        }
                                    }
                                    //           END IF ROTOK THEN ... ELSE
                                    //
                                    //           In the case of cancellation in updating SVA(q)
                                    //           .. recompute SVA(q)
                                    if (pow2((sva[q - 1] / aaqq)) <= rooteps) {
                                        if ((aaqq < rootbig) && (aaqq > rootsfmin)) {
                                            sva[q - 1] = Rnrm2(m, &a[(q - 1) * lda], 1) * d[q - 1];
                                        } else {
                                            t = zero;
                                            aaqq = one;
                                            Rlassq(m, &a[(q - 1) * lda], 1, t, aaqq);
                                            sva[q - 1] = t * sqrt(aaqq) * d[q - 1];
                                        }
                                    }
                                    if (pow2((aapp / aapp0)) <= rooteps) {
                                        if ((aapp < rootbig) && (aapp > rootsfmin)) {
                                            aapp = Rnrm2(m, &a[(p - 1) * lda], 1) * d[p - 1];
                                        } else {
                                            t = zero;
                                            aapp = one;
                                            Rlassq(m, &a[(p - 1) * lda], 1, t, aapp);
                                            aapp = t * sqrt(aapp) * d[p - 1];
                                        }
                                        sva[p - 1] = aapp;
                                    }
                                    //              end of OK rotation
                                } else {
                                    notrot++;
                                    pskipped++;
                                    ijblsk++;
                                }
                            } else {
                                notrot++;
                                pskipped++;
                                ijblsk++;
                            }
                            //
                            if ((i <= swband) && (ijblsk >= blskip)) {
                                sva[p - 1] = aapp;
                                notrot = 0;
                                goto statement_2011;
                            }
                            if ((i <= swband) && (pskipped > rowskip)) {
                                aapp = -aapp;
                                notrot = 0;
                                goto statement_2203;
                            }
                            //
                        }
                    //        end of the q-loop
                    statement_2203:
                        //
                        sva[p - 1] = aapp;
                        //
                    } else {
                        if (aapp == zero) {
                            notrot += min(jgl + kbl - 1, n) - jgl + 1;
                        }
                        if (aapp < zero) {
                            notrot = 0;
                        }
                    }
                    //
                }
                //     end of the p-loop
            }
        //     end of the jbc-loop
        statement_2011:
            // 2011 bailed out of the jbc-loop
            for (p = igl; p <= min(igl + kbl - 1, n); p = p + 1) {
                sva[p - 1] = abs(sva[p - 1]);
            }
            //
        }
        // 2000 :: end of the ibr-loop
        //
        //     .. update SVA(N)
        if ((sva[n - 1] < rootbig) && (sva[n - 1] > rootsfmin)) {
            sva[n - 1] = Rnrm2(m, &a[(n - 1) * lda], 1) * d[n - 1];
        } else {
            t = zero;
            aapp = one;
            Rlassq(m, &a[(n - 1) * lda], 1, t, aapp);
            sva[n - 1] = t * sqrt(aapp) * d[n - 1];
        }
        //
        //     Additional steering devices
        //
        if ((i < swband) && ((mxaapq <= roottol) || (iswrot <= n))) {
            swband = i;
        }
        //
        if ((i > swband + 1) && (mxaapq < castREAL(n) * tol) && (castREAL(n) * mxaapq * mxsinj < tol)) {
            goto statement_1994;
        }
        //
        if (notrot >= emptsw) {
            goto statement_1994;
        }
        //
    }
    //     end i=1:NSWEEP loop
    // #:) Reaching this poINTEGER means that the procedure has completed the given
    //     number of iterations.
    info = nsweep - 1;
    goto statement_1995;
statement_1994:
    // #:) Reaching this poINTEGER means that during the i-th sweep all pivots were
    //     below the given tolerance, causing early exit.
    //
    info = 0;
// #:) INFO = 0 confirms successful iterations.
statement_1995:
    //
    //     Sort the vector D.
    for (p = 1; p <= n - 1; p = p + 1) {
        q = iRamax(n - p + 1, &sva[p - 1], 1) + p - 1;
        if (p != q) {
            temp1 = sva[p - 1];
            sva[p - 1] = sva[q - 1];
            sva[q - 1] = temp1;
            temp1 = d[p - 1];
            d[p - 1] = d[q - 1];
            d[q - 1] = temp1;
            Rswap(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
            if (rsvec) {
                Rswap(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
            }
        }
    }
    //
    //     ..
    //     .. END OF Rgsvj0
    //     ..
}
