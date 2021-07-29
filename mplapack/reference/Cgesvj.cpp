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

void Cgesvj(const char *joba, const char *jobu, const char *jobv, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL *sva, INTEGER const mv, COMPLEX *v, INTEGER const ldv, COMPLEX *cwork, INTEGER const lwork, REAL *rwork, INTEGER const lrwork, INTEGER &info) {
    bool lsvec = false;
    bool uctol = false;
    bool rsvec = false;
    bool applv = false;
    bool upper = false;
    bool lower = false;
    bool lquery = false;
    const REAL one = 1.0;
    REAL ctol = 0.0;
    REAL epsln = 0.0;
    REAL rooteps = 0.0;
    REAL sfmin = 0.0;
    REAL rootsfmin = 0.0;
    REAL small = 0.0;
    REAL big = 0.0;
    REAL rootbig = 0.0;
    REAL bigtheta = 0.0;
    REAL tol = 0.0;
    REAL roottol = 0.0;
    INTEGER mvl = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    REAL skl = 0.0;
    bool noscale = false;
    bool goscale = false;
    INTEGER p = 0;
    const REAL zero = 0.0;
    REAL aapp = 0.0;
    REAL aaqq = 0.0;
    INTEGER q = 0;
    INTEGER ierr = 0;
    REAL sn = 0.0;
    REAL temp1 = 0.0;
    INTEGER emptsw = 0;
    INTEGER notrot = 0;
    INTEGER swband = 0;
    INTEGER kbl = 0;
    INTEGER nbl = 0;
    INTEGER blskip = 0;
    INTEGER rowskip = 0;
    INTEGER lkahead = 0;
    INTEGER n4 = 0;
    INTEGER n2 = 0;
    INTEGER n34 = 0;
    INTEGER i = 0;
    const INTEGER nsweep = 30;
    REAL mxaapq = 0.0;
    REAL mxsinj = 0.0;
    INTEGER iswrot = 0;
    INTEGER pskipped = 0;
    INTEGER ibr = 0;
    INTEGER igl = 0;
    INTEGER ir1 = 0;
    COMPLEX aapq = 0.0;
    REAL aapp0 = 0.0;
    bool rotok = false;
    REAL aapq1 = 0.0;
    COMPLEX ompq = 0.0;
    REAL aqoap = 0.0;
    REAL apoaq = 0.0;
    const REAL half = 0.5e0;
    REAL theta = 0.0;
    REAL t = 0.0;
    REAL cs = 0.0;
    REAL thsign = 0.0;
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     from BLAS
    //     from LAPACK
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     from BLAS
    //     from LAPACK
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    lsvec = Mlsame(jobu, "U") || Mlsame(jobu, "F");
    uctol = Mlsame(jobu, "C");
    rsvec = Mlsame(jobv, "V") || Mlsame(jobv, "J");
    applv = Mlsame(jobv, "A");
    upper = Mlsame(joba, "U");
    lower = Mlsame(joba, "L");
    //
    lquery = (lwork == -1) || (lrwork == -1);
    if (!(upper || lower || Mlsame(joba, "G"))) {
        info = -1;
    } else if (!(lsvec || uctol || Mlsame(jobu, "N"))) {
        info = -2;
    } else if (!(rsvec || applv || Mlsame(jobv, "N"))) {
        info = -3;
    } else if (m < 0) {
        info = -4;
    } else if ((n < 0) || (n > m)) {
        info = -5;
    } else if (lda < m) {
        info = -7;
    } else if (mv < 0) {
        info = -9;
    } else if ((rsvec && (ldv < n)) || (applv && (ldv < mv))) {
        info = -11;
    } else if (uctol && (rwork[1 - 1] <= one)) {
        info = -12;
    } else if ((lwork < (m + n)) && (!lquery)) {
        info = -13;
    } else if ((lrwork < max(n, (INTEGER)6)) && (!lquery)) {
        info = -15;
    } else {
        info = 0;
    }
    //
    //     #:(
    if (info != 0) {
        Mxerbla("Cgesvj", -info);
        return;
    } else if (lquery) {
        cwork[1 - 1] = m + n;
        rwork[1 - 1] = max(n, (INTEGER)6);
        return;
    }
    //
    // #:) Quick return for void matrix
    //
    if ((m == 0) || (n == 0)) {
        return;
    }
    //
    //     Set numerical parameters
    //     The stopping criterion for Jacobi rotations is
    //
    //     max_{i<>j}|A(:,i)^* * A(:,j)| / (||A(:,i)||*||A(:,j)||) < CTOL*EPS
    //
    //     where EPS is the round-off and CTOL is defined as follows:
    //
    if (uctol) {
        //        ... user controlled
        ctol = rwork[1 - 1];
    } else {
        //        ... default
        if (lsvec || rsvec || applv) {
            ctol = sqrt(castREAL(m));
        } else {
            ctol = castREAL(m);
        }
    }
    //     ... and the machine dependent parameters are
    //[!]  (Make sure that SLAMCH() works properly on the target machine.)
    //
    epsln = Rlamch("Epsilon");
    rooteps = sqrt(epsln);
    sfmin = Rlamch("SafeMinimum");
    rootsfmin = sqrt(sfmin);
    small = sfmin / epsln;
    big = Rlamch("Overflow");
    //     BIG         = ONE    / SFMIN
    rootbig = one / rootsfmin;
    //      LARGE = BIG / SQRT( DBLE( M*N ) )
    bigtheta = one / rooteps;
    //
    tol = ctol * epsln;
    roottol = sqrt(tol);
    //
    if (castREAL(m) * epsln >= one) {
        info = -4;
        Mxerbla("Cgesvj", -info);
        return;
    }
    //
    //     Initialize the right singular vector matrix.
    //
    if (rsvec) {
        mvl = n;
        Claset("A", mvl, n, czero, cone, v, ldv);
    } else if (applv) {
        mvl = mv;
    }
    rsvec = rsvec || applv;
    //
    //     Initialize SVA( 1:N ) = ( ||A e_i||_2, i = 1:N )
    //(!)  If necessary, scale A to protect the largest singular value
    //     from overflow. It is possible that saving the largest singular
    //     value destroys the information about the small ones.
    //     This initial scaling is almost minimal in the sense that the
    //     goal is to make sure that no column norm overflows, and that
    //     SQRT(N)*max_i SVA(i) does not overflow. If INFinite entries
    //     in A are detected, the procedure returns with INFO=-6.
    //
    skl = one / sqrt(castREAL(m) * castREAL(n));
    noscale = true;
    goscale = true;
    //
    if (lower) {
        //        the input matrix is M-by-N lower triangular (trapezoidal)
        for (p = 1; p <= n; p = p + 1) {
            aapp = zero;
            aaqq = one;
            Classq(m - p + 1, &a[(p - 1) + (p - 1) * lda], 1, aapp, aaqq);
            if (aapp > big) {
                info = -6;
                Mxerbla("Cgesvj", -info);
                return;
            }
            aaqq = sqrt(aaqq);
            if ((aapp < (big / aaqq)) && noscale) {
                sva[p - 1] = aapp * aaqq;
            } else {
                noscale = false;
                sva[p - 1] = aapp * (aaqq * skl);
                if (goscale) {
                    goscale = false;
                    for (q = 1; q <= p - 1; q = q + 1) {
                        sva[q - 1] = sva[q - 1] * skl;
                    }
                }
            }
        }
    } else if (upper) {
        //        the input matrix is M-by-N upper triangular (trapezoidal)
        for (p = 1; p <= n; p = p + 1) {
            aapp = zero;
            aaqq = one;
            Classq(p, &a[(p - 1) * lda], 1, aapp, aaqq);
            if (aapp > big) {
                info = -6;
                Mxerbla("Cgesvj", -info);
                return;
            }
            aaqq = sqrt(aaqq);
            if ((aapp < (big / aaqq)) && noscale) {
                sva[p - 1] = aapp * aaqq;
            } else {
                noscale = false;
                sva[p - 1] = aapp * (aaqq * skl);
                if (goscale) {
                    goscale = false;
                    for (q = 1; q <= p - 1; q = q + 1) {
                        sva[q - 1] = sva[q - 1] * skl;
                    }
                }
            }
        }
    } else {
        //        the input matrix is M-by-N general dense
        for (p = 1; p <= n; p = p + 1) {
            aapp = zero;
            aaqq = one;
            Classq(m, &a[(p - 1) * lda], 1, aapp, aaqq);
            if (aapp > big) {
                info = -6;
                Mxerbla("Cgesvj", -info);
                return;
            }
            aaqq = sqrt(aaqq);
            if ((aapp < (big / aaqq)) && noscale) {
                sva[p - 1] = aapp * aaqq;
            } else {
                noscale = false;
                sva[p - 1] = aapp * (aaqq * skl);
                if (goscale) {
                    goscale = false;
                    for (q = 1; q <= p - 1; q = q + 1) {
                        sva[q - 1] = sva[q - 1] * skl;
                    }
                }
            }
        }
    }
    //
    if (noscale) {
        skl = one;
    }
    //
    //     Move the smaller part of the spectrum from the underflow threshold
    //(!)  Start by determining the position of the nonzero entries of the
    //     array SVA() relative to ( SFMIN, BIG ).
    //
    aapp = zero;
    aaqq = big;
    for (p = 1; p <= n; p = p + 1) {
        if (sva[p - 1] != zero) {
            aaqq = min(aaqq, sva[p - 1]);
        }
        aapp = max(aapp, sva[p - 1]);
    }
    //
    // #:) Quick return for zero matrix
    //
    if (aapp == zero) {
        if (lsvec) {
            Claset("G", m, n, czero, cone, a, lda);
        }
        rwork[1 - 1] = one;
        rwork[2 - 1] = zero;
        rwork[3 - 1] = zero;
        rwork[4 - 1] = zero;
        rwork[5 - 1] = zero;
        rwork[6 - 1] = zero;
        return;
    }
    //
    // #:) Quick return for one-column matrix
    //
    if (n == 1) {
        if (lsvec) {
            Clascl("G", 0, 0, sva[1 - 1], skl, m, 1, &a[(1 - 1)], lda, ierr);
        }
        rwork[1 - 1] = one / skl;
        if (sva[1 - 1] >= sfmin) {
            rwork[2 - 1] = one;
        } else {
            rwork[2 - 1] = zero;
        }
        rwork[3 - 1] = zero;
        rwork[4 - 1] = zero;
        rwork[5 - 1] = zero;
        rwork[6 - 1] = zero;
        return;
    }
    //
    //     Protect small singular values from underflow, and try to
    //     avoid underflows/overflows in computing Jacobi rotations.
    //
    sn = sqrt(sfmin / epsln);
    temp1 = sqrt(big / castREAL(n));
    if ((aapp <= sn) || (aaqq >= temp1) || ((sn <= aaqq) && (aapp <= temp1))) {
        temp1 = min(big, REAL(temp1 / aapp));
        //         AAQQ  = AAQQ*TEMP1
        //         AAPP  = AAPP*TEMP1
    } else if ((aaqq <= sn) && (aapp <= temp1)) {
        temp1 = min(REAL(sn / aaqq), REAL(big / (aapp * sqrt(castREAL(n)))));
        //         AAQQ  = AAQQ*TEMP1
        //         AAPP  = AAPP*TEMP1
    } else if ((aaqq >= sn) && (aapp >= temp1)) {
        temp1 = max(sn / aaqq, temp1 / aapp);
        //         AAQQ  = AAQQ*TEMP1
        //         AAPP  = AAPP*TEMP1
    } else if ((aaqq <= sn) && (aapp >= temp1)) {
        temp1 = min(REAL(sn / aaqq), REAL(big / (sqrt(castREAL(n)) * aapp)));
        //         AAQQ  = AAQQ*TEMP1
        //         AAPP  = AAPP*TEMP1
    } else {
        temp1 = one;
    }
    //
    //     Scale, if necessary
    //
    if (temp1 != one) {
        Rlascl("G", 0, 0, one, temp1, n, 1, sva, n, ierr);
    }
    skl = temp1 * skl;
    if (skl != one) {
        Clascl(joba, 0, 0, one, skl, m, n, a, lda, ierr);
        skl = one / skl;
    }
    //
    //     Row-cyclic Jacobi SVD algorithm with column pivoting
    //
    emptsw = (n * (n - 1)) / 2;
    notrot = 0;
    //
    for (q = 1; q <= n; q = q + 1) {
        cwork[q - 1] = cone;
    }
    //
    swband = 3;
    //[TP] SWBAND is a tuning parameter [TP]. It is meaningful and effective
    //     if Cgesvj is used as a computational routine in the preconditioned
    //     Jacobi SVD algorithm Cgejsv. For sweeps i=1:SWBAND the procedure
    //     works on pivots inside a band-like region around the diagonal.
    //     The boundaries are determined dynamically, based on the number of
    //     pivots above a threshold.
    //
    kbl = min((INTEGER)8, n);
    //[TP] KBL is a tuning parameter that defines the tile size in the
    //     tiling of the p-q loops of pivot pairs. In general, an optimal
    //     parameters of the computer's memory.
    //
    nbl = n / kbl;
    if ((nbl * kbl) != n) {
        nbl++;
    }
    //
    blskip = kbl * kbl;
    //[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.
    //
    rowskip = min((INTEGER)5, kbl);
    //[TP] ROWSKIP is a tuning parameter.
    //
    lkahead = 1;
    //[TP] LKAHEAD is a tuning parameter.
    //
    //     Quasi block transformations, using the lower (upper) triangular
    //     structure of the input matrix. The quasi-block-cycling usually
    //     invokes cubic convergence. Big part of this cycle is done inside
    //
    if ((lower || upper) && (n > max((INTEGER)64, 4 * kbl))) {
        //[TP] The number of partition levels and the actual partition are
        //     tuning parameters.
        n4 = n / 4;
        n2 = n / 2;
        n34 = 3 * n4;
        if (applv) {
            q = 0;
        } else {
            q = 1;
        }
        //
        if (lower) {
            //
            //     This works very well on lower triangular matrices, in particular
            //     in the framework of the preconditioned Jacobi SVD (xGEJSV).
            //     The idea is simple:
            //     [+ 0 0 0]   Note that Jacobi transformations of [0 0]
            //     [+ + 0 0]                                       [0 0]
            //     [+ + x 0]   actually work on [x 0]              [x 0]
            //     [+ + x x]                    [x x].             [x x]
            //
            Cgsvj0(jobv, m - n34, n - n34, &a[((n34 + 1) - 1) + ((n34 + 1) - 1) * lda], lda, &cwork[(n34 + 1) - 1], &sva[(n34 + 1) - 1], mvl, &v[((n34 * q + 1) - 1) + ((n34 + 1) - 1) * ldv], ldv, epsln, sfmin, tol, 2, &cwork[(n + 1) - 1], lwork - n, ierr);
            //
            Cgsvj0(jobv, m - n2, n34 - n2, &a[((n2 + 1) - 1) + ((n2 + 1) - 1) * lda], lda, &cwork[(n2 + 1) - 1], &sva[(n2 + 1) - 1], mvl, &v[((n2 * q + 1) - 1) + ((n2 + 1) - 1) * ldv], ldv, epsln, sfmin, tol, 2, &cwork[(n + 1) - 1], lwork - n, ierr);
            //
            Cgsvj1(jobv, m - n2, n - n2, n4, &a[((n2 + 1) - 1) + ((n2 + 1) - 1) * lda], lda, &cwork[(n2 + 1) - 1], &sva[(n2 + 1) - 1], mvl, &v[((n2 * q + 1) - 1) + ((n2 + 1) - 1) * ldv], ldv, epsln, sfmin, tol, 1, &cwork[(n + 1) - 1], lwork - n, ierr);
            //
            Cgsvj0(jobv, m - n4, n2 - n4, &a[((n4 + 1) - 1) + ((n4 + 1) - 1) * lda], lda, &cwork[(n4 + 1) - 1], &sva[(n4 + 1) - 1], mvl, &v[((n4 * q + 1) - 1) + ((n4 + 1) - 1) * ldv], ldv, epsln, sfmin, tol, 1, &cwork[(n + 1) - 1], lwork - n, ierr);
            //
            Cgsvj0(jobv, m, n4, a, lda, cwork, sva, mvl, v, ldv, epsln, sfmin, tol, 1, &cwork[(n + 1) - 1], lwork - n, ierr);
            //
            Cgsvj1(jobv, m, n2, n4, a, lda, cwork, sva, mvl, v, ldv, epsln, sfmin, tol, 1, &cwork[(n + 1) - 1], lwork - n, ierr);
            //
        } else if (upper) {
            //
            Cgsvj0(jobv, n4, n4, a, lda, cwork, sva, mvl, v, ldv, epsln, sfmin, tol, 2, &cwork[(n + 1) - 1], lwork - n, ierr);
            //
            Cgsvj0(jobv, n2, n4, &a[((n4 + 1) - 1) * lda], lda, &cwork[(n4 + 1) - 1], &sva[(n4 + 1) - 1], mvl, &v[((n4 * q + 1) - 1) + ((n4 + 1) - 1) * ldv], ldv, epsln, sfmin, tol, 1, &cwork[(n + 1) - 1], lwork - n, ierr);
            //
            Cgsvj1(jobv, n2, n2, n4, a, lda, cwork, sva, mvl, v, ldv, epsln, sfmin, tol, 1, &cwork[(n + 1) - 1], lwork - n, ierr);
            //
            Cgsvj0(jobv, n2 + n4, n4, &a[((n2 + 1) - 1) * lda], lda, &cwork[(n2 + 1) - 1], &sva[(n2 + 1) - 1], mvl, &v[((n2 * q + 1) - 1) + ((n2 + 1) - 1) * ldv], ldv, epsln, sfmin, tol, 1, &cwork[(n + 1) - 1], lwork - n, ierr);
            //
        }
        //
    }
    //
    //     .. Row-cyclic pivot strategy with de Rijk's pivoting ..
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
                    //
                    q = iRamax(n - p + 1, &sva[p - 1], 1) + p - 1;
                    if (p != q) {
                        Cswap(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
                        if (rsvec) {
                            Cswap(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
                        }
                        temp1 = sva[p - 1];
                        sva[p - 1] = sva[q - 1];
                        sva[q - 1] = temp1;
                        aapq = cwork[p - 1];
                        cwork[p - 1] = cwork[q - 1];
                        cwork[q - 1] = aapq;
                    }
                    //
                    if (ir1 == 0) {
                        //
                        //        Column norms are periodically updated by explicit
                        //        norm computation.
                        //[!]     Caveat:
                        //        Unfortunately, some BLAS implementations compute RCnrm2(M,A(1,p),1)
                        //        as SQRT(S=CDOTC(M,A(1,p),1,A(1,p),1)), which may cause the result to
                        //        overflow for ||A(:,p)||_2 > SQRT(overflow_threshold), and to
                        //        underflow for ||A(:,p)||_2 < SQRT(underflow_threshold).
                        //        Hence, RCnrm2 cannot be trusted, not even in the case when
                        //        the true norm is far from the under(over)flow boundaries.
                        //        If properly implemented SCNRM2 is available, the IF-THEN-ELSE-END IF
                        //        below should be replaced with "AAPP = RCnrm2( M, A(1,p), 1 )".
                        //
                        if ((sva[p - 1] < rootbig) && (sva[p - 1] > rootsfmin)) {
                            sva[p - 1] = RCnrm2(m, &a[(p - 1) * lda], 1);
                        } else {
                            temp1 = zero;
                            aapp = one;
                            Classq(m, &a[(p - 1) * lda], 1, temp1, aapp);
                            sva[p - 1] = temp1 * sqrt(aapp);
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
                                        aapq = (Cdotc(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1) / aaqq) / aapp;
                                    } else {
                                        Ccopy(m, &a[(p - 1) * lda], 1, &cwork[(n + 1) - 1], 1);
                                        Clascl("G", 0, 0, aapp, one, m, 1, &cwork[(n + 1) - 1], lda, ierr);
                                        aapq = Cdotc(m, &cwork[(n + 1) - 1], 1, &a[(q - 1) * lda], 1) / aaqq;
                                    }
                                } else {
                                    rotok = aapp <= (aaqq / small);
                                    if (aapp > (small / aaqq)) {
                                        aapq = (Cdotc(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1) / aapp) / aaqq;
                                    } else {
                                        Ccopy(m, &a[(q - 1) * lda], 1, &cwork[(n + 1) - 1], 1);
                                        Clascl("G", 0, 0, aaqq, one, m, 1, &cwork[(n + 1) - 1], lda, ierr);
                                        aapq = Cdotc(m, &a[(p - 1) * lda], 1, &cwork[(n + 1) - 1], 1) / aapp;
                                    }
                                }
                                //
                                //                           AAPQ = AAPQ * CONJG( CWORK(p) ) * CWORK(q)
                                aapq1 = -abs(aapq);
                                mxaapq = max(mxaapq, REAL(-aapq1));
                                //
                                //        TO rotate or NOT to rotate, THAT is the question ...
                                //
                                if (abs(aapq1) > tol) {
                                    ompq = aapq / abs(aapq);
                                    //
                                    //           .. rotate
                                    //[RTD]      ROTATED = ROTATED + ONE
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
                                        theta = -half * abs(aqoap - apoaq) / aapq1;
                                        //
                                        if (abs(theta) > bigtheta) {
                                            //
                                            t = half / theta;
                                            cs = one;
                                            //
                                            Crot(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1, cs, conj(ompq) * t);
                                            if (rsvec) {
                                                Crot(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1, cs, conj(ompq) * t);
                                            }
                                            //
                                            sva[q - 1] = aaqq * sqrt(max(zero, REAL(one + t * apoaq * aapq1)));
                                            aapp = aapp * sqrt(max(zero, REAL(one - t * aqoap * aapq1)));
                                            mxsinj = max(mxsinj, REAL(abs(t)));
                                            //
                                        } else {
                                            //
                                            //                 .. choose correct signum for THETA and rotate
                                            //
                                            thsign = -sign(one, aapq1);
                                            t = one / (theta + thsign * sqrt(one + theta * theta));
                                            cs = sqrt(one / (one + t * t));
                                            sn = t * cs;
                                            //
                                            mxsinj = max(mxsinj, REAL(abs(sn)));
                                            sva[q - 1] = aaqq * sqrt(max(zero, REAL(one + t * apoaq * aapq1)));
                                            aapp = aapp * sqrt(max(zero, REAL(one - t * aqoap * aapq1)));
                                            //
                                            Crot(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1, cs, conj(ompq) * sn);
                                            if (rsvec) {
                                                Crot(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1, cs, conj(ompq) * sn);
                                            }
                                        }
                                        cwork[p - 1] = -cwork[q - 1] * ompq;
                                        //
                                    } else {
                                        //              .. have to use modified Gram-Schmidt like transformation
                                        Ccopy(m, &a[(p - 1) * lda], 1, &cwork[(n + 1) - 1], 1);
                                        Clascl("G", 0, 0, aapp, one, m, 1, &cwork[(n + 1) - 1], lda, ierr);
                                        Clascl("G", 0, 0, aaqq, one, m, 1, &a[(q - 1) * lda], lda, ierr);
                                        Caxpy(m, -aapq, &cwork[(n + 1) - 1], 1, &a[(q - 1) * lda], 1);
                                        Clascl("G", 0, 0, one, aaqq, m, 1, &a[(q - 1) * lda], lda, ierr);
                                        sva[q - 1] = aaqq * sqrt(max(zero, REAL(one - aapq1 * aapq1)));
                                        mxsinj = max(mxsinj, sfmin);
                                    }
                                    //           END IF ROTOK THEN ... ELSE
                                    //
                                    //           In the case of cancellation in updating SVA(q), SVA(p)
                                    //           recompute SVA(q), SVA(p).
                                    //
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
                                    if ((aapp / aapp0) <= rooteps) {
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
                                    //
                                } else {
                                    //                             A(:,p) and A(:,q) already numerically orthogonal
                                    if (ir1 == 0) {
                                        notrot++;
                                    }
                                    //[RTD]      SKIPPED  = SKIPPED + 1
                                    pskipped++;
                                }
                            } else {
                                //                          A(:,q) is zero column
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
                                        Ccopy(m, &a[(p - 1) * lda], 1, &cwork[(n + 1) - 1], 1);
                                        Clascl("G", 0, 0, aapp, one, m, 1, &cwork[(n + 1) - 1], lda, ierr);
                                        aapq = Cdotc(m, &cwork[(n + 1) - 1], 1, &a[(q - 1) * lda], 1) / aaqq;
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
                                        Ccopy(m, &a[(q - 1) * lda], 1, &cwork[(n + 1) - 1], 1);
                                        Clascl("G", 0, 0, aaqq, one, m, 1, &cwork[(n + 1) - 1], lda, ierr);
                                        aapq = Cdotc(m, &a[(p - 1) * lda], 1, &cwork[(n + 1) - 1], 1) / aapp;
                                    }
                                }
                                //
                                //                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q)
                                aapq1 = -abs(aapq);
                                mxaapq = max(mxaapq, REAL(-aapq1));
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
                                            sva[q - 1] = aaqq * sqrt(max(zero, REAL(one + t * apoaq * aapq1)));
                                            aapp = aapp * sqrt(max(zero, REAL(one - t * aqoap * aapq1)));
                                            mxsinj = max(mxsinj, REAL(abs(t)));
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
                                            mxsinj = max(mxsinj, REAL(abs(sn)));
                                            sva[q - 1] = aaqq * sqrt(max(zero, REAL(one + t * apoaq * aapq1)));
                                            aapp = aapp * sqrt(max(zero, REAL(one - t * aqoap * aapq1)));
                                            //
                                            Crot(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1, cs, conj(ompq) * sn);
                                            if (rsvec) {
                                                Crot(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1, cs, conj(ompq) * sn);
                                            }
                                        }
                                        cwork[p - 1] = -cwork[q - 1] * ompq;
                                        //
                                    } else {
                                        //              .. have to use modified Gram-Schmidt like transformation
                                        if (aapp > aaqq) {
                                            Ccopy(m, &a[(p - 1) * lda], 1, &cwork[(n + 1) - 1], 1);
                                            Clascl("G", 0, 0, aapp, one, m, 1, &cwork[(n + 1) - 1], lda, ierr);
                                            Clascl("G", 0, 0, aaqq, one, m, 1, &a[(q - 1) * lda], lda, ierr);
                                            Caxpy(m, -aapq, &cwork[(n + 1) - 1], 1, &a[(q - 1) * lda], 1);
                                            Clascl("G", 0, 0, one, aaqq, m, 1, &a[(q - 1) * lda], lda, ierr);
                                            sva[q - 1] = aaqq * sqrt(max(zero, REAL(one - aapq1 * aapq1)));
                                            mxsinj = max(mxsinj, sfmin);
                                        } else {
                                            Ccopy(m, &a[(q - 1) * lda], 1, &cwork[(n + 1) - 1], 1);
                                            Clascl("G", 0, 0, aaqq, one, m, 1, &cwork[(n + 1) - 1], lda, ierr);
                                            Clascl("G", 0, 0, aapp, one, m, 1, &a[(p - 1) * lda], lda, ierr);
                                            Caxpy(m, -conj(aapq), &cwork[(n + 1) - 1], 1, &a[(p - 1) * lda], 1);
                                            Clascl("G", 0, 0, one, aapp, m, 1, &a[(p - 1) * lda], lda, ierr);
                                            sva[p - 1] = aapp * sqrt(max(zero, REAL(one - aapq1 * aapq1)));
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
        if ((i > swband + 1) && (mxaapq < sqrt(castREAL(n)) * tol) && (castREAL(n) * mxaapq * mxsinj < tol)) {
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
    //     Sort the singular values and find how many are above
    //     the underflow threshold.
    //
    n2 = 0;
    n4 = 0;
    for (p = 1; p <= n - 1; p = p + 1) {
        q = iRamax(n - p + 1, &sva[p - 1], 1) + p - 1;
        if (p != q) {
            temp1 = sva[p - 1];
            sva[p - 1] = sva[q - 1];
            sva[q - 1] = temp1;
            Cswap(m, &a[(p - 1) * lda], 1, &a[(q - 1) * lda], 1);
            if (rsvec) {
                Cswap(mvl, &v[(p - 1) * ldv], 1, &v[(q - 1) * ldv], 1);
            }
        }
        if (sva[p - 1] != zero) {
            n4++;
            if (sva[p - 1] * skl > sfmin) {
                n2++;
            }
        }
    }
    if (sva[n - 1] != zero) {
        n4++;
        if (sva[n - 1] * skl > sfmin) {
            n2++;
        }
    }
    //
    //     Normalize the left singular vectors.
    //
    if (lsvec || uctol) {
        for (p = 1; p <= n4; p = p + 1) {
            //            CALL CRscal( M, ONE / SVA( p ), A( 1, p ), 1 )
            Clascl("G", 0, 0, sva[p - 1], one, m, 1, &a[(p - 1) * lda], m, ierr);
        }
    }
    //
    //     Scale the product of Jacobi rotations.
    //
    if (rsvec) {
        for (p = 1; p <= n; p = p + 1) {
            temp1 = one / RCnrm2(mvl, &v[(p - 1) * ldv], 1);
            CRscal(mvl, temp1, &v[(p - 1) * ldv], 1);
        }
    }
    //
    //     Undo scaling, if necessary (and possible).
    if (((skl > one) && (sva[1 - 1] < (big / skl))) || ((skl < one) && (sva[(max(n2, (INTEGER)1) - 1)]) > (sfmin / skl))) {
        for (p = 1; p <= n; p = p + 1) {
            sva[p - 1] = skl * sva[p - 1];
        }
        skl = one;
    }
    //
    rwork[1 - 1] = skl;
    //     The singular values of A are SKL*SVA(1:N). If SKL.NE.ONE
    //     then some of the singular values may overflow or underflow and
    //     the spectrum is given in this factored representation.
    //
    rwork[2 - 1] = castREAL(n4);
    //     N4 is the number of computed nonzero singular values of A.
    //
    rwork[3 - 1] = castREAL(n2);
    //     N2 is the number of singular values of A greater than SFMIN.
    //     If N2<N, SVA(N2:N) contains ZEROS and/or denormalized numbers
    //     that may carry some information.
    //
    rwork[4 - 1] = castREAL(i);
    //     i is the index of the last sweep before declaring convergence.
    //
    rwork[5 - 1] = mxaapq;
    //     MXAAPQ is the largest absolute value of scaled pivots in the
    //     last sweep
    //
    rwork[6 - 1] = mxsinj;
    //     MXSINJ is the largest absolute value of the sines of Jacobi angles
    //     in the last sweep
    //
    //     ..
    //     .. END OF Cgesvj
    //     ..
}
