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

void Cgsvj1(const char *jobv, INTEGER const m, INTEGER const n, INTEGER const n1, COMPLEX *a, INTEGER const lda, COMPLEX *d, REAL *sva, INTEGER const mv, COMPLEX *v, INTEGER const ldv, REAL const eps, REAL const sfmin, REAL const tol, INTEGER const nsweep, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    INTEGER kbl = 0;
    INTEGER nblr = 0;
    INTEGER nblc = 0;
    INTEGER blskip = 0;
    INTEGER rowskip = 0;
    INTEGER swband = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    REAL mxaapq = 0.0;
    REAL mxsinj = 0.0;
    INTEGER iswrot = 0;
    INTEGER pskipped = 0;
    INTEGER ibr = 0;
    INTEGER igl = 0;
    INTEGER jbc = 0;
    INTEGER jgl = 0;
    INTEGER ijblsk = 0;
    INTEGER p = 0;
    REAL aapp = 0.0;
    INTEGER q = 0;
    REAL aaqq = 0.0;
    REAL aapp0 = 0.0;
    bool rotok = false;
    COMPLEX aapq = 0.0;
    INTEGER ierr = 0;
    REAL aapq1 = 0.0;
    COMPLEX ompq = 0.0;
    REAL aqoap = 0.0;
    REAL apoaq = 0.0;
    const REAL half = 0.5e0;
    REAL theta = 0.0;
    REAL t = 0.0;
    REAL cs = 0.0;
    REAL thsign = 0.0;
    REAL sn = 0.0;
    REAL temp1 = 0.0;
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     .. from BLAS
    //     .. from LAPACK
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
    } else if (n1 < 0) {
        info = -4;
    } else if (lda < m) {
        info = -6;
    } else if ((rsvec || applv) && (mv < 0)) {
        info = -9;
    } else if ((rsvec && (ldv < n)) || (applv && (ldv < mv))) {
        info = -11;
    } else if (tol <= eps) {
        info = -14;
    } else if (nsweep < 0) {
        info = -15;
    } else if (lwork < m) {
        info = -17;
    } else {
        info = 0;
    }
    //
    //     #:(
    if (info != 0) {
        Mxerbla("Cgsvj1", -info);
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
    //     LARGE = BIG / SQRT( DBLE( M*N ) )
    bigtheta = one / rooteps;
    roottol = sqrt(tol);
    //
    //     .. Initialize the right singular vector matrix ..
    //
    //     RSVEC = LSAME( JOBV, 'Y' )
    //
    emptsw = n1 * (n - n1);
    notrot = 0;
    //
    //     .. Row-cyclic pivot strategy with de Rijk's pivoting ..
    //
    kbl = min(8, n);
    nblr = n1 / kbl;
    if ((nblr * kbl) != n1) {
        nblr++;
    }
    //
    //     .. the tiling is nblr-by-nblc [tiles]
    //
    nblc = (n - n1) / kbl;
    if ((nblc * kbl) != (n - n1)) {
        nblc++;
    }
    blskip = (pow2(kbl)) + 1;
    //[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.
    //
    rowskip = min(5, kbl);
    //[TP] ROWSKIP is a tuning parameter.
    swband = 0;
    //[TP] SWBAND is a tuning parameter. It is meaningful and effective
    //     if Cgesvj is used as a computational routine in the preconditioned
    //     Jacobi SVD algorithm Cgejsv.
    //
    //     | *   *   * [x] [x] [x]|
    //     | *   *   * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
    //     | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
    //     |[x] [x] [x] *   *   * |
    //
    for (i = 1; i <= nsweep; i = i + 1) {
        //
        //     .. go go go ...
        //
        mxaapq = zero;
        mxsinj = zero;
        iswrot = 0;
        //
        notrot = 0;
        pskipped = 0;
        //
        //     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs
        //     1 <= p < q <= N. This is the first step toward a blocked implementation
        //     of the rotations. New implementation, based on block transformations,
        //     is under development.
        //
        for (ibr = 1; ibr <= nblr; ibr = ibr + 1) {
            //
            igl = (ibr - 1) * kbl + 1;
            //
            // ... go to the off diagonal blocks
            //
            igl = (ibr - 1) * kbl + 1;
            //
            //            DO 2010 jbc = ibr + 1, NBL
            for (jbc = 1; jbc <= nblc; jbc = jbc + 1) {
                //
                jgl = (jbc - 1) * kbl + n1 + 1;
                //
                //        doing the block at ( ibr, jbc )
                //
                ijblsk = 0;
                for (p = igl; p <= min(igl + kbl - 1, n1); p = p + 1) {
                    //
                    aapp = sva[p - 1];
                    if (aapp > zero) {
                        //
                        pskipped = 0;
                        //
                        for (q = jgl; q <= min(jgl + kbl - 1, n); q = q + 1) {
                            //
                            aaqq = sva[q - 1];
                            if (aaqq > zero) {
                                aapp0 = aapp;
                                //
                                //     .. M x 2 Jacobi SVD ..
                                //
                                //        Safe Gram matrix computation
                                //
                                if (aaqq >= one) {
                                    if (aapp >= aaqq) {
                                        rotok = (small * aapp) <= aaqq;
                                    } else {
                                        rotok = (small * aaqq) <= aapp;
                                    }
                                    if (aapp < (big / aaqq)) {
                                        aapq = (Cdotc(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1) / aaqq) / aapp;
                                    } else {
                                        Ccopy(m, &a[(p - 1) * lda], 1, work, 1);
                                        Clascl("G", 0, 0, aapp, one, m, 1, work, lda, ierr);
                                        aapq = Cdotc(m, work, 1, &a[(q - 1) * lda], 1) / aaqq;
                                    }
                                } else {
                                    if (aapp >= aaqq) {
                                        rotok = aapp <= (aaqq / small);
                                    } else {
                                        rotok = aaqq <= (aapp / small);
                                    }
                                    if (aapp > (small / aaqq)) {
                                        aapq = (Cdotc(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1) / max(aaqq, aapp)) / min(aaqq, aapp);
                                    } else {
                                        Ccopy(m, &a[(q - 1) * lda], 1, work, 1);
                                        Clascl("G", 0, 0, aaqq, one, m, 1, work, lda, ierr);
                                        aapq = Cdotc(m, &a[(p - 1) * lda], 1, work, 1) / aapp;
                                    }
                                }
                                //
                                //                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q)
                                aapq1 = -abs(aapq);
                                mxaapq = max(mxaapq, -aapq1);
                                //
                                //        TO rotate or NOT to rotate, THAT is the question ...
                                //
                                if (abs(aapq1) > tol) {
                                    ompq = aapq / abs(aapq);
                                    notrot = 0;
                                    //[RTD]      ROTATED  = ROTATED + 1
                                    pskipped = 0;
                                    iswrot++;
                                    //
                                    if (rotok) {
                                        //
                                        aqoap = aaqq / aapp;
                                        apoaq = aapp / aaqq;
                                        theta = -half * abs(aqoap - apoaq) / aapq1;
                                        if (aaqq > aapp0) {
                                            theta = -theta;
                                        }
                                        //
                                        if (abs(theta) > bigtheta) {
                                            t = half / theta;
                                            cs = one;
                                            Crot(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1, cs, conj(ompq) * t);
                                            if (rsvec) {
                                                Crot(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1, cs, conj(ompq) * t);
                                            }
                                            sva[q - 1] = aaqq * sqrt(max(zero, one + t * apoaq * aapq1));
                                            aapp = aapp * sqrt(max(zero, one - t * aqoap * aapq1));
                                            mxsinj = max(mxsinj, abs(t));
                                        } else {
                                            //
                                            //                 .. choose correct signum for THETA and rotate
                                            //
                                            thsign = -sign(one, aapq1);
                                            if (aaqq > aapp0) {
                                                thsign = -thsign;
                                            }
                                            t = one / (theta + thsign * sqrt(one + theta * theta));
                                            cs = sqrt(one / (one + t * t));
                                            sn = t * cs;
                                            mxsinj = max(mxsinj, abs(sn));
                                            sva[q - 1] = aaqq * sqrt(max(zero, one + t * apoaq * aapq1));
                                            aapp = aapp * sqrt(max(zero, one - t * aqoap * aapq1));
                                            //
                                            Crot(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1, cs, conj(ompq) * sn);
                                            if (rsvec) {
                                                Crot(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1, cs, conj(ompq) * sn);
                                            }
                                        }
                                        d[p - 1] = -d[q - 1] * ompq;
                                        //
                                    } else {
                                        //              .. have to use modified Gram-Schmidt like transformation
                                        if (aapp > aaqq) {
                                            Ccopy(m, &a[(p - 1) * lda], 1, work, 1);
                                            Clascl("G", 0, 0, aapp, one, m, 1, work, lda, ierr);
                                            Clascl("G", 0, 0, aaqq, one, m, 1, &a[(q - 1) * lda], lda, ierr);
                                            Caxpy(m, -aapq, work, 1, &a[(q - 1) * lda], 1);
                                            Clascl("G", 0, 0, one, aaqq, m, 1, &a[(q - 1) * lda], lda, ierr);
                                            sva[q - 1] = aaqq * sqrt(max(zero, one - aapq1 * aapq1));
                                            mxsinj = max(mxsinj, sfmin);
                                        } else {
                                            Ccopy(m, &a[(q - 1) * lda], 1, work, 1);
                                            Clascl("G", 0, 0, aaqq, one, m, 1, work, lda, ierr);
                                            Clascl("G", 0, 0, aapp, one, m, 1, &a[(p - 1) * lda], lda, ierr);
                                            Caxpy(m, -conj(aapq), work, 1, &a[(p - 1) * lda], 1);
                                            Clascl("G", 0, 0, one, aapp, m, 1, &a[(p - 1) * lda], lda, ierr);
                                            sva[p - 1] = aapp * sqrt(max(zero, one - aapq1 * aapq1));
                                            mxsinj = max(mxsinj, sfmin);
                                        }
                                    }
                                    //           END IF ROTOK THEN ... ELSE
                                    //
                                    //           In the case of cancellation in updating SVA(q), SVA(p)
                                    //           .. recompute SVA(q), SVA(p)
                                    if (pow2((sva[q - 1] / aaqq)) <= rooteps) {
                                        if ((aaqq < rootbig) && (aaqq > rootsfmin)) {
                                            sva[q - 1] = RCnrm2(m, &a[(q - 1) * lda], 1);
                                        } else {
                                            t = zero;
                                            aaqq = one;
                                            Classq(m, &a[(q - 1) * lda], 1, t, aaqq);
                                            sva[q - 1] = t * sqrt(aaqq);
                                        }
                                    }
                                    if (pow2((aapp / aapp0)) <= rooteps) {
                                        if ((aapp < rootbig) && (aapp > rootsfmin)) {
                                            aapp = RCnrm2(m, &a[(p - 1) * lda], 1);
                                        } else {
                                            t = zero;
                                            aapp = one;
                                            Classq(m, &a[(p - 1) * lda], 1, t, aapp);
                                            aapp = t * sqrt(aapp);
                                        }
                                        sva[p - 1] = aapp;
                                    }
                                    //              end of OK rotation
                                } else {
                                    notrot++;
                                    //[RTD]      SKIPPED  = SKIPPED  + 1
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
                        //
                        if (aapp == zero) {
                            notrot += min(jgl + kbl - 1, n) - jgl + 1;
                        }
                        if (aapp < zero) {
                            notrot = 0;
                        }
                        //
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
            //**
        }
        // 2000 :: end of the ibr-loop
        //
        //     .. update SVA(N)
        if ((sva[n - 1] < rootbig) && (sva[n - 1] > rootsfmin)) {
            sva[n - 1] = RCnrm2(m, &a[(n - 1) * lda], 1);
        } else {
            t = zero;
            aapp = one;
            Classq(m, &a[(n - 1) * lda], 1, t, aapp);
            sva[n - 1] = t * sqrt(aapp);
        }
        //
        //     Additional steering devices
        //
        if ((i < swband) && ((mxaapq <= roottol) || (iswrot <= n))) {
            swband = i;
        }
        //
        if ((i > swband + 1) && (mxaapq < sqrt(castREAL(n) * tol) && (castREAL(n) * mxaapq * mxsinj < tol))) {
            goto statement_1994;
        }
        //
        if (notrot >= emptsw) {
            goto statement_1994;
        }
        //
    }
    //     end i=1:NSWEEP loop
    //
    // #:( Reaching this point means that the procedure has not converged.
    info = nsweep - 1;
    goto statement_1995;
//
statement_1994:
    // #:) Reaching this point means numerical convergence after the i-th
    //     sweep.
    //
    info = 0;
// #:) INFO = 0 confirms successful iterations.
statement_1995:
    //
    //     Sort the vector SVA() of column norms.
    for (p = 1; p <= n - 1; p = p + 1) {
        q = iRamax(n - p + 1, &sva[p - 1], 1) + p - 1;
        if (p != q) {
            temp1 = sva[p - 1];
            sva[p - 1] = sva[q - 1];
            sva[q - 1] = temp1;
            aapq = d[p - 1];
            d[p - 1] = d[q - 1];
            d[q - 1] = aapq;
            Cswap(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
            if (rsvec) {
                Cswap(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
            }
        }
    }
    //
    //     ..
    //     .. END OF Cgsvj1
    //     ..
}
