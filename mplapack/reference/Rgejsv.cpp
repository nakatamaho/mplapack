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

void Rgejsv(const char *joba, const char *jobu, const char *jobv, const char *jobr, const char *jobt, const char *jobp, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *sva, REAL *u, INTEGER const ldu, REAL *v, INTEGER const ldv, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER &info) {
    bool lsvec = false;
    bool jracc = false;
    bool rsvec = false;
    bool rowpiv = false;
    bool l2rank = false;
    bool l2aber = false;
    bool errest = false;
    bool l2tran = false;
    bool l2kill = false;
    bool defr = false;
    bool l2pert = false;
    INTEGER n1 = 0;
    REAL epsln = 0.0;
    REAL sfmin = 0.0;
    REAL small = 0.0;
    REAL big = 0.0;
    const REAL one = 1.0;
    REAL scalem = 0.0;
    bool noscal = false;
    bool goscal = false;
    INTEGER p = 0;
    const REAL zero = 0.0;
    REAL aapp = 0.0;
    REAL aaqq = 0.0;
    INTEGER warning = 0;
    INTEGER ierr = 0;
    bool transp = false;
    REAL aatmax = 0.0;
    REAL aatmin = 0.0;
    REAL xsc = 0.0;
    REAL temp1 = 0.0;
    REAL entra = 0.0;
    REAL entrat = 0.0;
    REAL big1 = 0.0;
    INTEGER q = 0;
    bool kill = false;
    REAL uscal1 = 0.0;
    REAL uscal2 = 0.0;
    INTEGER nr = 0;
    bool almort = false;
    REAL maxprj = 0.0;
    REAL sconda = 0.0;
    REAL condr1 = 0.0;
    REAL condr2 = 0.0;
    INTEGER numrank = 0;
    REAL cond_ok = 0.0;
    //
    //     Test the input arguments
    //
    lsvec = Mlsame(jobu, "U") || Mlsame(jobu, "F");
    jracc = Mlsame(jobv, "J");
    rsvec = Mlsame(jobv, "V") || jracc;
    rowpiv = Mlsame(joba, "F") || Mlsame(joba, "G");
    l2rank = Mlsame(joba, "R");
    l2aber = Mlsame(joba, "A");
    errest = Mlsame(joba, "E") || Mlsame(joba, "G");
    l2tran = Mlsame(jobt, "T");
    l2kill = Mlsame(jobr, "R");
    defr = Mlsame(jobr, "N");
    l2pert = Mlsame(jobp, "P");
    //
    if (!(rowpiv || l2rank || l2aber || errest || Mlsame(joba, "C"))) {
        info = -1;
    } else if (!(lsvec || Mlsame(jobu, "N") || Mlsame(jobu, "W"))) {
        info = -2;
    } else if (!(rsvec || Mlsame(jobv, "N") || Mlsame(jobv, "W")) || (jracc && (!lsvec))) {
        info = -3;
    } else if (!(l2kill || defr)) {
        info = -4;
    } else if (!(l2tran || Mlsame(jobt, "N"))) {
        info = -5;
    } else if (!(l2pert || Mlsame(jobp, "N"))) {
        info = -6;
    } else if (m < 0) {
        info = -7;
    } else if ((n < 0) || (n > m)) {
        info = -8;
    } else if (lda < m) {
        info = -10;
    } else if (lsvec && (ldu < m)) {
        info = -13;
    } else if (rsvec && (ldv < n)) {
        info = -15;
    } else if ((!(lsvec || rsvec || errest) && (lwork < max({(INTEGER)7, 4 * n + 1, 2 * m + n}))) || (!(lsvec || rsvec) && errest && (lwork < max({(INTEGER)7, 4 * n + n * n, 2 * m + n}))) || (lsvec && (!rsvec) && (lwork < max({(INTEGER)7, 2 * m + n, 4 * n + 1}))) || (rsvec && (!lsvec) && (lwork < max({(INTEGER)7, 2 * m + n, 4 * n + 1}))) || (lsvec && rsvec && (!jracc) && (lwork < max({2 * m + n, 6 * n + 2 * n * n}))) || (lsvec && rsvec && jracc && lwork < max({2 * m + n, 4 * n + n * n, 2 * n + n * n + 6}))) {
        info = -17;
    } else {
        //        #:)
        info = 0;
    }
    //
    if (info != 0) {
        //       #:(
        Mxerbla("Rgejsv", -info);
        return;
    }
    //
    //     Quick return for void matrix (Y3K safe)
    // #:)
    if ((m == 0) || (n == 0)) {
        iwork[0] = iwork[1] = iwork[2] = 0;
        work[0] = work[1] = work[2] = work[3] = work[4] = work[5] = work[6] = 0.0;
        return;
    }
    //
    //     Determine whether the matrix U should be M x N or M x M
    //
    if (lsvec) {
        n1 = n;
        if (Mlsame(jobu, "F")) {
            n1 = m;
        }
    }
    //
    //     Set numerical parameters
    //
    //!    NOTE: Make sure DLAMCH() does not fail on the target architecture.
    //
    epsln = Rlamch("Epsilon");
    sfmin = Rlamch("SafeMinimum");
    small = sfmin / epsln;
#if defined ___MPLAPACK_BUILD_WITH_DD___
    big = Rlamch_dd("Q");
#elif defined ___MPLAPACK_BUILD_WITH_QD___
    big = Rlamch_qd("Q");
#else
    big = Rlamch("Overflow");
#endif
    //     BIG   = ONE / SFMIN
    //
    //     Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N
    //
    //(!)  If necessary, scale SVA() to protect the largest norm from
    //     overflow. It is possible that this scaling pushes the smallest
    //     column norm left from the underflow threshold (extreme case).
    //
    scalem = one / sqrt(castREAL(m * n));
    noscal = true;
    goscal = true;
    for (p = 1; p <= n; p = p + 1) {
        aapp = zero;
        aaqq = one;
        Rlassq(m, &a[(p - 1) * lda], 1, aapp, aaqq);
        if (aapp > big) {
            info = -9;
            Mxerbla("Rgejsv", -info);
            return;
        }
        aaqq = sqrt(aaqq);
        if ((aapp < (big / aaqq)) && noscal) {
            sva[p - 1] = aapp * aaqq;
        } else {
            noscal = false;
            sva[p - 1] = aapp * (aaqq * scalem);
            if (goscal) {
                goscal = false;
                Rscal(p - 1, scalem, sva, 1);
            }
        }
    }
    //
    if (noscal) {
        scalem = one;
    }
    //
    aapp = zero;
    aaqq = big;
    for (p = 1; p <= n; p = p + 1) {
        aapp = max(aapp, sva[p - 1]);
        if (sva[p - 1] != zero) {
            aaqq = min(aaqq, sva[p - 1]);
        }
    }
    //
    //     Quick return for zero M x N matrix
    // #:)
    if (aapp == zero) {
        if (lsvec) {
            Rlaset("G", m, n1, zero, one, u, ldu);
        }
        if (rsvec) {
            Rlaset("G", n, n, zero, one, v, ldv);
        }
        work[1 - 1] = one;
        work[2 - 1] = one;
        if (errest) {
            work[3 - 1] = one;
        }
        if (lsvec && rsvec) {
            work[4 - 1] = one;
            work[5 - 1] = one;
        }
        if (l2tran) {
            work[6 - 1] = zero;
            work[7 - 1] = zero;
        }
        iwork[1 - 1] = 0;
        iwork[2 - 1] = 0;
        iwork[3 - 1] = 0;
        return;
    }
    //
    //     Issue warning if denormalized column norms detected. Override the
    //     high relative accuracy request. Issue licence to kill columns
    //     (set them to zero) whose norm is less than sigma_max / BIG (roughly).
    // #:(
    warning = 0;
    if (aaqq <= sfmin) {
        l2rank = true;
        l2kill = true;
        warning = 1;
    }
    //
    //     Quick return for one-column matrix
    // #:)
    if (n == 1) {
        //
        if (lsvec) {
            Rlascl("G", 0, 0, sva[1 - 1], scalem, m, 1, &a[(1 - 1)], lda, ierr);
            Rlacpy("A", m, 1, a, lda, u, ldu);
            //           computing all M left singular vectors of the M x 1 matrix
            if (n1 != n) {
                Rgeqrf(m, n, u, ldu, work, &work[(n + 1) - 1], lwork - n, ierr);
                Rorgqr(m, n1, 1, u, ldu, work, &work[(n + 1) - 1], lwork - n, ierr);
                Rcopy(m, &a[(1 - 1)], 1, &u[(1 - 1)], 1);
            }
        }
        if (rsvec) {
            v[(1 - 1)] = one;
        }
        if (sva[1 - 1] < (big * scalem)) {
            sva[1 - 1] = sva[1 - 1] / scalem;
            scalem = one;
        }
        work[1 - 1] = one / scalem;
        work[2 - 1] = one;
        if (sva[1 - 1] != zero) {
            iwork[1 - 1] = 1;
            if ((sva[1 - 1] / scalem) >= sfmin) {
                iwork[2 - 1] = 1;
            } else {
                iwork[2 - 1] = 0;
            }
        } else {
            iwork[1 - 1] = 0;
            iwork[2 - 1] = 0;
        }
        iwork[3 - 1] = 0;
        if (errest) {
            work[3 - 1] = one;
        }
        if (lsvec && rsvec) {
            work[4 - 1] = one;
            work[5 - 1] = one;
        }
        if (l2tran) {
            work[6 - 1] = zero;
            work[7 - 1] = zero;
        }
        return;
        //
    }
    //
    transp = false;
    l2tran = l2tran && (m == n);
    //
    aatmax = -one;
    aatmin = big;
    if (rowpiv || l2tran) {
        //
        //     Compute the row norms, needed to determine row pivoting sequence
        //     (in the case of heavily row weighted A, row pivoting is strongly
        //     advised) and to collect information needed to compare the
        //     structures of A * A^t and A^t * A (in the case L2TRAN.EQ..TRUE.).
        //
        if (l2tran) {
            for (p = 1; p <= m; p = p + 1) {
                xsc = zero;
                temp1 = one;
                Rlassq(n, &a[(p - 1)], lda, xsc, temp1);
                //              Rlassq gets both the ell_2 and the ell_infinity norm
                //              in one pass through the vector
                work[(m + n + p) - 1] = xsc * scalem;
                work[(n + p) - 1] = xsc * (scalem * sqrt(temp1));
                aatmax = max(aatmax, work[(n + p) - 1]);
                if (work[(n + p) - 1] != zero) {
                    aatmin = min(aatmin, work[(n + p) - 1]);
                }
            }
        } else {
            for (p = 1; p <= m; p = p + 1) {
                work[(m + n + p) - 1] = scalem * abs(a[(p - 1) + (iRamax(n, &a[(p - 1)], lda))]);
                aatmax = max(aatmax, work[(m + n + p) - 1]);
                aatmin = min(aatmin, work[(m + n + p) - 1]);
            }
        }
        //
    }
    //
    //     For square matrix A try to determine whether A^t  would be  better
    //     input for the preconditioned Jacobi SVD, with faster convergence.
    //     The decision is based on an O(N) function of the vector of column
    //     and row norms of A, based on the Shannon entropy. This should give
    //     the right choice in most cases when the difference actually matters.
    //     It may fail and pick the slower converging side.
    //
    entra = zero;
    entrat = zero;
    if (l2tran) {
        //
        xsc = zero;
        temp1 = one;
        Rlassq(n, sva, 1, xsc, temp1);
        temp1 = one / temp1;
        //
        entra = zero;
        for (p = 1; p <= n; p = p + 1) {
            big1 = (pow2((sva[p - 1] / xsc))) * temp1;
            if (big1 != zero) {
                entra += big1 * log(big);
            }
        }
        entra = -entra / log(castREAL(n - 1));
        //
        //        Now, SVA().^2/Trace(A^t * A) is a point in the probability simplex.
        //        It is derived from the diagonal of  A^t * A.  Do the same with the
        //        diagonal of A * A^t, compute the entropy of the corresponding
        //        probability distribution. Note that A * A^t and A^t * A have the
        //        same trace.
        //
        entrat = zero;
        for (p = n + 1; p <= n + m; p = p + 1) {
            big1 = (pow2((work[p - 1] / xsc))) * temp1;
            if (big1 != zero) {
                entrat += big1 * log(big1);
            }
        }
        entrat = -entrat / log(castREAL(m - 1));
        //
        //        Analyze the entropies and decide A or A^t. Smaller entropy
        //        usually means better input for the algorithm.
        //
        transp = (entrat < entra);
        //
        //        If A^t is better than A, transpose A.
        //
        if (transp) {
            //           In an optimal implementation, this trivial transpose
            //           should be replaced with faster transpose.
            for (p = 1; p <= n - 1; p = p + 1) {
                for (q = p + 1; q <= n; q = q + 1) {
                    temp1 = a[(q - 1) + (p - 1) * lda];
                    a[(q - 1) + (p - 1) * lda] = a[(p - 1) + (q - 1) * lda];
                    a[(p - 1) + (q - 1) * lda] = temp1;
                }
            }
            for (p = 1; p <= n; p = p + 1) {
                work[(m + n + p) - 1] = sva[p - 1];
                sva[p - 1] = work[(n + p) - 1];
            }
            temp1 = aapp;
            aapp = aatmax;
            aatmax = temp1;
            temp1 = aaqq;
            aaqq = aatmin;
            aatmin = temp1;
            kill = lsvec;
            lsvec = rsvec;
            rsvec = kill;
            if (lsvec) {
                n1 = n;
            }
            //
            rowpiv = true;
        }
        //
    }
    //     END IF L2TRAN
    //
    //     Scale the matrix so that its maximal singular value remains less
    //     than DSQRT(BIG) -- the matrix is scaled so that its maximal column
    //     has Euclidean norm equal to DSQRT(BIG/N). The only reason to keep
    //     DSQRT(BIG) instead of BIG is the fact that Rgejsv uses LAPACK and
    //     BLAS routines that, in some implementations, are not capable of
    //     working in the full interval [SFMIN,BIG] and that they may provoke
    //     overflows in the intermediate results. If the singular values spread
    //     from SFMIN to BIG, then Rgesvj will compute them. So, in that case,
    //     one should use Rgesvj instead of Rgejsv.
    //
    big1 = sqrt(big);
    temp1 = sqrt(big / castREAL(n));
    //
    Rlascl("G", 0, 0, aapp, temp1, n, 1, sva, n, ierr);
    if (aaqq > (aapp * sfmin)) {
        aaqq = (aaqq / aapp) * temp1;
    } else {
        aaqq = (aaqq * temp1) / aapp;
    }
    temp1 = temp1 * scalem;
    Rlascl("G", 0, 0, aapp, temp1, m, n, a, lda, ierr);
    //
    //     To undo scaling at the end of this procedure, multiply the
    //     computed singular values with USCAL2 / USCAL1.
    //
    uscal1 = temp1;
    uscal2 = aapp;
    //
    if (l2kill) {
        //        L2KILL enforces computation of nonzero singular values in
        //        the restricted range of condition number of the initial A,
        //        sigma_max(A) / sigma_min(A) approx. DSQRT(BIG)/DSQRT(SFMIN).
        xsc = sqrt(sfmin);
    } else {
        xsc = small;
        //
        //        Now, if the condition number of A is too big,
        //        sigma_max(A) / sigma_min(A) .GT. DSQRT(BIG/N) * EPSLN / SFMIN,
        //        as a precaution measure, the full SVD is computed using Rgesvj
        //        with accumulated Jacobi rotations. This provides numerically
        //        more robust computation, at the cost of slightly increased run
        //        time. Depending on the concrete implementation of BLAS and LAPACK
        //        (i.e. how they behave in presence of extreme ill-conditioning) the
        //        implementor may decide to remove this switch.
        if ((aaqq < sqrt(sfmin)) && lsvec && rsvec) {
            jracc = true;
        }
        //
    }
    if (aaqq < xsc) {
        for (p = 1; p <= n; p = p + 1) {
            if (sva[p - 1] < xsc) {
                Rlaset("A", m, 1, zero, zero, &a[(p - 1) * lda], lda);
                sva[p - 1] = zero;
            }
        }
    }
    //
    //     Preconditioning using QR factorization with pivoting
    //
    if (rowpiv) {
        //        Optional row permutation (Bjoerck row pivoting):
        //        A result by Cox and Higham shows that the Bjoerck's
        //        row pivoting combined with standard column pivoting
        //        has similar effect as Powell-Reid complete pivoting.
        //        The ell-infinity norms of A are made nonincreasing.
        for (p = 1; p <= m - 1; p = p + 1) {
            q = iRamax(m - p + 1, &work[(m + n + p) - 1], 1) + p - 1;
            iwork[(2 * n + p) - 1] = q;
            if (p != q) {
                temp1 = work[(m + n + p) - 1];
                work[(m + n + p) - 1] = work[(m + n + q) - 1];
                work[(m + n + q) - 1] = temp1;
            }
        }
        Rlaswp(n, a, lda, 1, m - 1, &iwork[(2 * n + 1) - 1], 1);
    }
    //
    //     End of the preparation phase (scaling, optional sorting and
    //     transposing, optional flushing of small columns).
    //
    //     Preconditioning
    //
    //     If the full SVD is needed, the right singular vectors are computed
    //     from a matrix equation, and for that we need theoretical analysis
    //     of the Businger-Golub pivoting. So we use Rgeqp3 as the first RR QRF.
    //     In all other cases the first RR QRF can be chosen by other criteria
    //     (eg speed by replacing global with restricted window pivoting, such
    //     as in SGEQPX from TOMS # 782). Good results will be obtained using
    //     SGEQPX with properly (!) chosen numerical parameters.
    //     Any improvement of Rgeqp3 improves overall performance of Rgejsv.
    //
    //     A * P1 = Q1 * [ R1^t 0]^t:
    for (p = 1; p <= n; p = p + 1) {
        //        .. all columns are free columns
        iwork[p - 1] = 0;
    }
    Rgeqp3(m, n, a, lda, iwork, work, &work[(n + 1) - 1], lwork - n, ierr);
    //
    //     The upper triangular matrix R1 from the first QRF is inspected for
    //     rank deficiency and possibilities for deflation, or possible
    //     ill-conditioning. Depending on the user specified flag L2RANK,
    //     the procedure explores possibilities to reduce the numerical
    //     rank by inspecting the computed upper triangular factor. If
    //     L2RANK or L2ABER are up, then Rgejsv will compute the SVD of
    //     A + dA, where ||dA|| <= f(M,N)*EPSLN.
    //
    nr = 1;
    if (l2aber) {
        //        Standard absolute error bound suffices. All sigma_i with
        //        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
        //        aggressive enforcement of lower numerical rank by introducing a
        //        backward error of the order of N*EPSLN*||A||.
        temp1 = sqrt(castREAL(n)) * epsln;
        for (p = 2; p <= n; p = p + 1) {
            if (abs(a[(p - 1) + (p - 1) * lda]) >= (temp1 * abs(a[(1 - 1)]))) {
                nr++;
            } else {
                goto statement_3002;
            }
        }
    statement_3002:;
    } else if (l2rank) {
        //        .. similarly as above, only slightly more gentle (less aggressive).
        //        Sudden drop on the diagonal of R1 is used as the criterion for
        //        close-to-rank-deficient.
        temp1 = sqrt(sfmin);
        for (p = 2; p <= n; p = p + 1) {
            if ((abs(a[(p - 1) + (p - 1) * lda]) < (epsln * abs(a[((p - 1) - 1) + ((p - 1) - 1) * lda]))) || (abs(a[(p - 1) + (p - 1) * lda]) < small) || (l2kill && (abs(a[(p - 1) + (p - 1) * lda]) < temp1))) {
                goto statement_3402;
            }
            nr++;
        }
    statement_3402:;
        //
    } else {
        //        The goal is high relative accuracy. However, if the matrix
        //        has high scaled condition number the relative accuracy is in
        //        general not feasible. Later on, a condition number estimator
        //        will be deployed to estimate the scaled condition number.
        //        Here we just remove the underflowed part of the triangular
        //        factor. This prevents the situation in which the code is
        //        working hard to get the accuracy not warranted by the data.
        temp1 = sqrt(sfmin);
        for (p = 2; p <= n; p = p + 1) {
            if ((abs(a[(p - 1) + (p - 1) * lda]) < small) || (l2kill && (abs(a[(p - 1) + (p - 1) * lda]) < temp1))) {
                goto statement_3302;
            }
            nr++;
        }
    statement_3302:;
        //
    }
    //
    almort = false;
    if (nr == n) {
        maxprj = one;
        for (p = 2; p <= n; p = p + 1) {
            temp1 = abs(a[(p - 1) + (p - 1) * lda]) / sva[iwork[p - 1] - 1];
            maxprj = min(maxprj, temp1);
        }
        if (pow2(maxprj) >= one - castREAL(n) * epsln) {
            almort = true;
        }
    }
    //
    sconda = -one;
    condr1 = -one;
    condr2 = -one;
    //
    if (errest) {
        if (n == nr) {
            if (rsvec) {
                //              .. V is available as workspace
                Rlacpy("U", n, n, a, lda, v, ldv);
                for (p = 1; p <= n; p = p + 1) {
                    temp1 = sva[iwork[p - 1] - 1];
                    Rscal(p, one / temp1, &v[(p - 1) * ldv], 1);
                }
                Rpocon("U", n, v, ldv, one, temp1, &work[(n + 1) - 1], &iwork[(2 * n + m + 1) - 1], ierr);
            } else if (lsvec) {
                //              .. U is available as workspace
                Rlacpy("U", n, n, a, lda, u, ldu);
                for (p = 1; p <= n; p = p + 1) {
                    temp1 = sva[iwork[p - 1] - 1];
                    Rscal(p, one / temp1, &u[(p - 1) * ldu], 1);
                }
                Rpocon("U", n, u, ldu, one, temp1, &work[(n + 1) - 1], &iwork[(2 * n + m + 1) - 1], ierr);
            } else {
                Rlacpy("U", n, n, a, lda, &work[(n + 1) - 1], n);
                for (p = 1; p <= n; p = p + 1) {
                    temp1 = sva[iwork[p - 1] - 1];
                    Rscal(p, one / temp1, &work[(n + (p - 1) * n + 1) - 1], 1);
                }
                //           .. the columns of R are scaled to have unit Euclidean lengths.
                Rpocon("U", n, &work[(n + 1) - 1], n, one, temp1, &work[(n + n * n + 1) - 1], &iwork[(2 * n + m + 1) - 1], ierr);
            }
            sconda = one / sqrt(temp1);
            //           SCONDA is an estimate of DSQRT(||(R^t * R)^(-1)||_1).
            //           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
        } else {
            sconda = -one;
        }
    }
    //
    l2pert = l2pert && (abs(a[(1 - 1)] / a[(nr - 1) + (nr - 1) * lda]) > sqrt(big1));
    //     If there is no violent scaling, artificial perturbation is not needed.
    //
    //     Phase 3:
    //
    if (!(rsvec || lsvec)) {
        //
        //         Singular Values only
        //
        //         .. transpose A(1:NR,1:N)
        for (p = 1; p <= min(n - 1, nr); p = p + 1) {
            Rcopy(n - p, &a[(p - 1) + ((p + 1) - 1) * lda], lda, &a[((p + 1) - 1) + (p - 1) * lda], 1);
        }
        //
        //        The following two DO-loops introduce small relative perturbation
        //        into the strict upper triangle of the lower triangular matrix.
        //        Small entries below the main diagonal are also changed.
        //        This modification is useful if the computing environment does not
        //        provide/allow FLUSH TO ZERO underflow, for it prevents many
        //        annoying denormalized numbers in case of strongly scaled matrices.
        //        The perturbation is structured so that it does not introduce any
        //        new perturbation of the singular values, and it does not destroy
        //        the job done by the preconditioner.
        //        The licence for this perturbation is in the variable L2PERT, which
        //        should be .FALSE. if FLUSH TO ZERO underflow is active.
        //
        if (!almort) {
            //
            if (l2pert) {
                //              XSC = DSQRT(SMALL)
                xsc = epsln / castREAL(n);
                for (q = 1; q <= nr; q = q + 1) {
                    temp1 = xsc * abs(a[(q - 1) + (q - 1) * lda]);
                    for (p = 1; p <= n; p = p + 1) {
                        if (((p > q) && (abs(a[(p - 1) + (q - 1) * lda]) <= temp1)) || (p < q)) {
                            a[(p - 1) + (q - 1) * lda] = sign(temp1, a[(p - 1) + (q - 1) * lda]);
                        }
                    }
                }
            } else {
                Rlaset("U", nr - 1, nr - 1, zero, zero, &a[(2 - 1) * lda], lda);
            }
            //
            //            .. second preconditioning using the QR factorization
            //
            Rgeqrf(n, nr, a, lda, work, &work[(n + 1) - 1], lwork - n, ierr);
            //
            //           .. and transpose upper to lower triangular
            for (p = 1; p <= nr - 1; p = p + 1) {
                Rcopy(nr - p, &a[(p - 1) + ((p + 1) - 1) * lda], lda, &a[((p + 1) - 1) + (p - 1) * lda], 1);
            }
            //
        }
        //
        //           Row-cyclic Jacobi SVD algorithm with column pivoting
        //
        //           .. again some perturbation (a "background noise") is added
        //           to drown denormals
        if (l2pert) {
            //              XSC = DSQRT(SMALL)
            xsc = epsln / castREAL(n);
            for (q = 1; q <= nr; q = q + 1) {
                temp1 = xsc * abs(a[(q - 1) + (q - 1) * lda]);
                for (p = 1; p <= nr; p = p + 1) {
                    if (((p > q) && (abs(a[(p - 1) + (q - 1) * lda]) <= temp1)) || (p < q)) {
                        a[(p - 1) + (q - 1) * lda] = sign(temp1, a[(p - 1) + (q - 1) * lda]);
                    }
                }
            }
        } else {
            Rlaset("U", nr - 1, nr - 1, zero, zero, &a[(2 - 1) * lda], lda);
        }
        //
        //           .. and one-sided Jacobi rotations are started on a lower
        //           triangular matrix (plus perturbation which is ignored in
        //           the part which destroys triangular form (confusing?!))
        //
        Rgesvj("L", "NoU", "NoV", nr, nr, a, lda, sva, n, v, ldv, work, lwork, info);
        //
        scalem = work[1 - 1];
        numrank = nint(work[2 - 1]);
        //
    } else if (rsvec && (!lsvec)) {
        //
        //        -> Singular Values and Right Singular Vectors <-
        //
        if (almort) {
            //
            //           .. in this case NR equals N
            for (p = 1; p <= nr; p = p + 1) {
                Rcopy(n - p + 1, &a[(p - 1) + (p - 1) * lda], lda, &v[(p - 1) + (p - 1) * ldv], 1);
            }
            Rlaset("Upper", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
            //
            Rgesvj("L", "U", "N", n, nr, v, ldv, sva, nr, a, lda, work, lwork, info);
            scalem = work[1 - 1];
            numrank = castINTEGER(work[2 - 1]);
            //
        } else {
            //
            //        .. two more QR factorizations ( one QRF is not enough, two require
            //        accumulated product of Jacobi rotations, three are perfect )
            //
            Rlaset("Lower", nr - 1, nr - 1, zero, zero, &a[(2 - 1)], lda);
            Rgelqf(nr, n, a, lda, work, &work[(n + 1) - 1], lwork - n, ierr);
            Rlacpy("Lower", nr, nr, a, lda, v, ldv);
            Rlaset("Upper", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
            Rgeqrf(nr, nr, v, ldv, &work[(n + 1) - 1], &work[(2 * n + 1) - 1], lwork - 2 * n, ierr);
            for (p = 1; p <= nr; p = p + 1) {
                Rcopy(nr - p + 1, &v[(p - 1) + (p - 1) * ldv], ldv, &v[(p - 1) + (p - 1) * ldv], 1);
            }
            Rlaset("Upper", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
            //
            Rgesvj("Lower", "U", "N", nr, nr, v, ldv, sva, nr, u, ldu, &work[(n + 1) - 1], lwork, info);
            scalem = work[(n + 1) - 1];
            numrank = nint(work[(n + 2) - 1]);
            if (nr < n) {
                Rlaset("A", n - nr, nr, zero, zero, &v[((nr + 1) - 1)], ldv);
                Rlaset("A", nr, n - nr, zero, zero, &v[((nr + 1) - 1) * ldv], ldv);
                Rlaset("A", n - nr, n - nr, zero, one, &v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
            }
            //
            Rormlq("Left", "Transpose", n, n, nr, a, lda, work, v, ldv, &work[(n + 1) - 1], lwork - n, ierr);
            //
        }
        //
        for (p = 1; p <= n; p = p + 1) {
            Rcopy(n, &v[(p - 1)], ldv, &a[(iwork[p - 1] - 1)], lda);
        }
        Rlacpy("All", n, n, a, lda, v, ldv);
        //
        if (transp) {
            Rlacpy("All", n, n, v, ldv, u, ldu);
        }
        //
    } else if (lsvec && (!rsvec)) {
        //
        //        .. Singular Values and Left Singular Vectors                 ..
        //
        //        .. second preconditioning step to avoid need to accumulate
        //        Jacobi rotations in the Jacobi iterations.
        for (p = 1; p <= nr; p = p + 1) {
            Rcopy(n - p + 1, &a[(p - 1) + (p - 1) * lda], lda, &u[(p - 1) + (p - 1) * ldu], 1);
        }
        Rlaset("Upper", nr - 1, nr - 1, zero, zero, &u[(2 - 1) * ldu], ldu);
        //
        Rgeqrf(n, nr, u, ldu, &work[(n + 1) - 1], &work[(2 * n + 1) - 1], lwork - 2 * n, ierr);
        //
        for (p = 1; p <= nr - 1; p = p + 1) {
            Rcopy(nr - p, &u[(p - 1) + ((p + 1) - 1) * ldu], ldu, &u[((p + 1) - 1) + (p - 1) * ldu], 1);
        }
        Rlaset("Upper", nr - 1, nr - 1, zero, zero, &u[(2 - 1) * ldu], ldu);
        //
        Rgesvj("Lower", "U", "N", nr, nr, u, ldu, sva, nr, a, lda, &work[(n + 1) - 1], lwork - n, info);
        scalem = work[(n + 1) - 1];
        numrank = nint(work[(n + 2) - 1]);
        //
        if (nr < m) {
            Rlaset("A", m - nr, nr, zero, zero, &u[((nr + 1) - 1)], ldu);
            if (nr < n1) {
                Rlaset("A", nr, n1 - nr, zero, zero, &u[((nr + 1) - 1) * ldu], ldu);
                Rlaset("A", m - nr, n1 - nr, zero, one, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
            }
        }
        //
        Rormqr("Left", "No Tr", m, n1, n, a, lda, work, u, ldu, &work[(n + 1) - 1], lwork - n, ierr);
        //
        if (rowpiv) {
            Rlaswp(n1, u, ldu, 1, m - 1, &iwork[(2 * n + 1) - 1], -1);
        }
        //
        for (p = 1; p <= n1; p = p + 1) {
            xsc = one / Rnrm2(m, &u[(p - 1) * ldu], 1);
            Rscal(m, xsc, &u[(p - 1) * ldu], 1);
        }
        //
        if (transp) {
            Rlacpy("All", n, n, u, ldu, v, ldv);
        }
        //
    } else {
        //
        //        .. Full SVD ..
        //
        if (!jracc) {
            //
            if (!almort) {
                //
                //           Second Preconditioning Step (QRF [with pivoting])
                //           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
                //           equivalent to an LQF CALL. Since in many libraries the QRF
                //           seems to be better optimized than the LQF, we do explicit
                //           transpose and use the QRF. This is subject to changes in an
                //           optimized implementation of Rgejsv.
                //
                for (p = 1; p <= nr; p = p + 1) {
                    Rcopy(n - p + 1, &a[(p - 1) + (p - 1) * lda], lda, &v[(p - 1) + (p - 1) * ldv], 1);
                }
                //
                //           .. the following two loops perturb small entries to avoid
                //           denormals in the second QR factorization, where they are
                //           as good as zeros. This is done to avoid painfully slow
                //           computation with denormals. The relative size of the perturbation
                //           is a parameter that can be changed by the implementer.
                //           This perturbation device will be obsolete on machines with
                //           properly implemented arithmetic.
                //           To switch it off, set L2PERT=.FALSE. To remove it from  the
                //           code, remove the action under L2PERT=.TRUE., leave the ELSE part.
                //           The following two loops should be blocked and fused with the
                //           transposed copy above.
                //
                if (l2pert) {
                    xsc = sqrt(small);
                    for (q = 1; q <= nr; q = q + 1) {
                        temp1 = xsc * abs(v[(q - 1) + (q - 1) * ldv]);
                        for (p = 1; p <= n; p = p + 1) {
                            if ((p > q) && (abs(v[(p - 1) + (q - 1) * ldv]) <= temp1) || (p < q)) {
                                v[(p - 1) + (q - 1) * ldv] = sign(temp1, v[(p - 1) + (q - 1) * ldv]);
                            }
                            if (p < q) {
                                v[(p - 1) + (q - 1) * ldv] = -v[(p - 1) + (q - 1) * ldv];
                            }
                        }
                    }
                } else {
                    Rlaset("U", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
                }
                //
                //           Estimate the row scaled condition number of R1
                //           (If R1 is rectangular, N > NR, then the condition number
                //           of the leading NR x NR submatrix is estimated.)
                //
                Rlacpy("L", nr, nr, v, ldv, &work[(2 * n + 1) - 1], nr);
                for (p = 1; p <= nr; p = p + 1) {
                    temp1 = Rnrm2(nr - p + 1, &work[(2 * n + (p - 1) * nr + p) - 1], 1);
                    Rscal(nr - p + 1, one / temp1, &work[(2 * n + (p - 1) * nr + p) - 1], 1);
                }
                Rpocon("Lower", nr, &work[(2 * n + 1) - 1], nr, one, temp1, &work[(2 * n + nr * nr + 1) - 1], &iwork[(m + 2 * n + 1) - 1], ierr);
                condr1 = one / sqrt(temp1);
                //           .. here need a second opinion on the condition number
                //           .. then assume worst case scenario
                //           R1 is OK for inverse <=> CONDR1 .LT. DBLE(N)
                //           more conservative    <=> CONDR1 .LT. DSQRT(DBLE(N))
                //
                cond_ok = sqrt(castREAL(nr));
                //[TP]       COND_OK is a tuning parameter.
                //
                if (condr1 < cond_ok) {
                    //              .. the second QRF without pivoting. Note: in an optimized
                    //              implementation, this QRF should be implemented as the QRF
                    //              of a lower triangular matrix.
                    //              R1^t = Q2 * R2
                    Rgeqrf(n, nr, v, ldv, &work[(n + 1) - 1], &work[(2 * n + 1) - 1], lwork - 2 * n, ierr);
                    //
                    if (l2pert) {
                        xsc = sqrt(small) / epsln;
                        for (p = 2; p <= nr; p = p + 1) {
                            for (q = 1; q <= p - 1; q = q + 1) {
                                temp1 = xsc * min(abs(v[(p - 1) + (p - 1) * ldv]), abs(v[(q - 1) + (q - 1) * ldv]));
                                if (abs(v[(q - 1) + (p - 1) * ldv]) <= temp1) {
                                    v[(q - 1) + (p - 1) * ldv] = sign(temp1, v[(q - 1) + (p - 1) * ldv]);
                                }
                            }
                        }
                    }
                    //
                    if (nr != n) {
                        Rlacpy("A", n, nr, v, ldv, &work[(2 * n + 1) - 1], n);
                    }
                    //              .. save ...
                    //
                    //           .. this transposed copy should be better than naive
                    for (p = 1; p <= nr - 1; p = p + 1) {
                        Rcopy(nr - p, &v[(p - 1) + ((p + 1) - 1) * ldv], ldv, &v[((p + 1) - 1) + (p - 1) * ldv], 1);
                    }
                    //
                    condr2 = condr1;
                    //
                } else {
                    //
                    //              .. ill-conditioned case: second QRF with pivoting
                    //              Note that windowed pivoting would be equally good
                    //              numerically, and more run-time efficient. So, in
                    //              an optimal implementation, the next call to Rgeqp3
                    //              should be replaced with eg. CALL SGEQPX (ACM TOMS #782)
                    //              with properly (carefully) chosen parameters.
                    //
                    //              R1^t * P2 = Q2 * R2
                    for (p = 1; p <= nr; p = p + 1) {
                        iwork[(n + p) - 1] = 0;
                    }
                    Rgeqp3(n, nr, v, ldv, &iwork[(n + 1) - 1], &work[(n + 1) - 1], &work[(2 * n + 1) - 1], lwork - 2 * n, ierr);
                    //*               CALL Rgeqrf( N, NR, V, LDV, WORK(N+1), WORK(2*N+1),
                    //*     $              LWORK-2*N, IERR )
                    if (l2pert) {
                        xsc = sqrt(small);
                        for (p = 2; p <= nr; p = p + 1) {
                            for (q = 1; q <= p - 1; q = q + 1) {
                                temp1 = xsc * min(abs(v[(p - 1) + (p - 1) * ldv]), abs(v[(q - 1) + (q - 1) * ldv]));
                                if (abs(v[(q - 1) + (p - 1) * ldv]) <= temp1) {
                                    v[(q - 1) + (p - 1) * ldv] = sign(temp1, v[(q - 1) + (p - 1) * ldv]);
                                }
                            }
                        }
                    }
                    //
                    Rlacpy("A", n, nr, v, ldv, &work[(2 * n + 1) - 1], n);
                    //
                    if (l2pert) {
                        xsc = sqrt(small);
                        for (p = 2; p <= nr; p = p + 1) {
                            for (q = 1; q <= p - 1; q = q + 1) {
                                temp1 = xsc * min(abs(v[(p - 1) + (p - 1) * ldv]), abs(v[(q - 1) + (q - 1) * ldv]));
                                v[(p - 1) + (q - 1) * ldv] = -sign(temp1, v[(q - 1) + (p - 1) * ldv]);
                            }
                        }
                    } else {
                        Rlaset("L", nr - 1, nr - 1, zero, zero, &v[(2 - 1)], ldv);
                    }
                    //              Now, compute R2 = L3 * Q3, the LQ factorization.
                    Rgelqf(nr, nr, v, ldv, &work[(2 * n + n * nr + 1) - 1], &work[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
                    //              .. and estimate the condition number
                    Rlacpy("L", nr, nr, v, ldv, &work[(2 * n + n * nr + nr + 1) - 1], nr);
                    for (p = 1; p <= nr; p = p + 1) {
                        temp1 = Rnrm2(p, &work[(2 * n + n * nr + nr + p) - 1], nr);
                        Rscal(p, one / temp1, &work[(2 * n + n * nr + nr + p) - 1], nr);
                    }
                    Rpocon("L", nr, &work[(2 * n + n * nr + nr + 1) - 1], nr, one, temp1, &work[(2 * n + n * nr + nr + nr * nr + 1) - 1], &iwork[(m + 2 * n + 1) - 1], ierr);
                    condr2 = one / sqrt(temp1);
                    //
                    if (condr2 >= cond_ok) {
                        //                 .. save the Householder vectors used for Q3
                        //                 (this overwrites the copy of R2, as it will not be
                        //                 needed in this branch, but it does not overwritte the
                        //                 Huseholder vectors of Q2.).
                        Rlacpy("U", nr, nr, v, ldv, &work[(2 * n + 1) - 1], n);
                        //                 .. and the rest of the information on Q3 is in
                        //                 WORK(2*N+N*NR+1:2*N+N*NR+N)
                    }
                    //
                }
                //
                if (l2pert) {
                    xsc = sqrt(small);
                    for (q = 2; q <= nr; q = q + 1) {
                        temp1 = xsc * v[(q - 1) + (q - 1) * ldv];
                        for (p = 1; p <= q - 1; p = p + 1) {
                            //                    V(p,q) = - DSIGN( TEMP1, V(q,p) )
                            v[(p - 1) + (q - 1) * ldv] = -sign(temp1, v[(p - 1) + (q - 1) * ldv]);
                        }
                    }
                } else {
                    Rlaset("U", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
                }
                //
                //        Second preconditioning finished; continue with Jacobi SVD
                //        The input matrix is lower trinagular.
                //
                //        Recover the right singular vectors as solution of a well
                //        conditioned triangular matrix equation.
                //
                if (condr1 < cond_ok) {
                    //
                    Rgesvj("L", "U", "N", nr, nr, v, ldv, sva, nr, u, ldu, &work[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, info);
                    scalem = work[(2 * n + n * nr + nr + 1) - 1];
                    numrank = nint(work[(2 * n + n * nr + nr + 2) - 1]);
                    for (p = 1; p <= nr; p = p + 1) {
                        Rcopy(nr, &v[(p - 1) * ldv], 1, &u[(p - 1) * ldu], 1);
                        Rscal(nr, sva[p - 1], &v[(p - 1) * ldv], 1);
                    }
                    //
                    //        .. pick the right matrix equation and solve it
                    //
                    if (nr == n) {
                        // :))             .. best case, R1 is inverted. The solution of this matrix
                        //                 equation is Q2*V2 = the product of the Jacobi rotations
                        //                 used in Rgesvj, premultiplied with the orthogonal matrix
                        //                 from the second QR factorization.
                        Rtrsm("L", "U", "N", "N", nr, nr, one, a, lda, v, ldv);
                    } else {
                        //                 .. R1 is well conditioned, but non-square. Transpose(R2)
                        //                 is inverted to get the product of the Jacobi rotations
                        //                 used in Rgesvj. The Q-factor from the second QR
                        //                 factorization is then built in explicitly.
                        Rtrsm("L", "U", "T", "N", nr, nr, one, &work[(2 * n + 1) - 1], n, v, ldv);
                        if (nr < n) {
                            Rlaset("A", n - nr, nr, zero, zero, &v[((nr + 1) - 1)], ldv);
                            Rlaset("A", nr, n - nr, zero, zero, &v[((nr + 1) - 1) * ldv], ldv);
                            Rlaset("A", n - nr, n - nr, zero, one, &v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
                        }
                        Rormqr("L", "N", n, n, nr, &work[(2 * n + 1) - 1], n, &work[(n + 1) - 1], v, ldv, &work[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
                    }
                    //
                } else if (condr2 < cond_ok) {
                    //
                    // :)           .. the input matrix A is very likely a relative of
                    //              the Kahan matrix :)
                    //              The matrix R2 is inverted. The solution of the matrix equation
                    //              is Q3^T*V3 = the product of the Jacobi rotations (appplied to
                    //              the lower triangular L3 from the LQ factorization of
                    //              R2=L3*Q3), pre-multiplied with the transposed Q3.
                    Rgesvj("L", "U", "N", nr, nr, v, ldv, sva, nr, u, ldu, &work[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, info);
                    scalem = work[(2 * n + n * nr + nr + 1) - 1];
                    numrank = nint(work[(2 * n + n * nr + nr + 2) - 1]);
                    for (p = 1; p <= nr; p = p + 1) {
                        Rcopy(nr, &v[(p - 1) * ldv], 1, &u[(p - 1) * ldu], 1);
                        Rscal(nr, sva[p - 1], &u[(p - 1) * ldu], 1);
                    }
                    Rtrsm("L", "U", "N", "N", nr, nr, one, &work[(2 * n + 1) - 1], n, u, ldu);
                    //              .. apply the permutation from the second QR factorization
                    for (q = 1; q <= nr; q = q + 1) {
                        for (p = 1; p <= nr; p = p + 1) {
                            work[(2 * n + n * nr + nr + iwork[(n + p) - 1]) - 1] = u[(p - 1) + (q - 1) * ldu];
                        }
                        for (p = 1; p <= nr; p = p + 1) {
                            u[(p - 1) + (q - 1) * ldu] = work[(2 * n + n * nr + nr + p) - 1];
                        }
                    }
                    if (nr < n) {
                        Rlaset("A", n - nr, nr, zero, zero, &v[((nr + 1) - 1)], ldv);
                        Rlaset("A", nr, n - nr, zero, zero, &v[((nr + 1) - 1) * ldv], ldv);
                        Rlaset("A", n - nr, n - nr, zero, one, &v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
                    }
                    Rormqr("L", "N", n, n, nr, &work[(2 * n + 1) - 1], n, &work[(n + 1) - 1], v, ldv, &work[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
                } else {
                    //              Last line of defense.
                    // #:(          This is a rather pathological case: no scaled condition
                    //              improvement after two pivoted QR factorizations. Other
                    //              possibility is that the rank revealing QR factorization
                    //              or the condition estimator has failed, or the COND_OK
                    //              is set very close to ONE (which is unnecessary). Normally,
                    //              this branch should never be executed, but in rare cases of
                    //              failure of the RRQR or condition estimator, the last line of
                    //              defense ensures that Rgejsv completes the task.
                    //              Compute the full SVD of L3 using Rgesvj with explicit
                    //              accumulation of Jacobi rotations.
                    Rgesvj("L", "U", "V", nr, nr, v, ldv, sva, nr, u, ldu, &work[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, info);
                    scalem = work[(2 * n + n * nr + nr + 1) - 1];
                    numrank = nint(work[(2 * n + n * nr + nr + 2) - 1]);
                    if (nr < n) {
                        Rlaset("A", n - nr, nr, zero, zero, &v[((nr + 1) - 1)], ldv);
                        Rlaset("A", nr, n - nr, zero, zero, &v[((nr + 1) - 1) * ldv], ldv);
                        Rlaset("A", n - nr, n - nr, zero, one, &v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
                    }
                    Rormqr("L", "N", n, n, nr, &work[(2 * n + 1) - 1], n, &work[(n + 1) - 1], v, ldv, &work[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
                    //
                    Rormlq("L", "T", nr, nr, nr, &work[(2 * n + 1) - 1], n, &work[(2 * n + n * nr + 1) - 1], u, ldu, &work[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
                    for (q = 1; q <= nr; q = q + 1) {
                        for (p = 1; p <= nr; p = p + 1) {
                            work[(2 * n + n * nr + nr + iwork[(n + p) - 1]) - 1] = u[(p - 1) + (q - 1) * ldu];
                        }
                        for (p = 1; p <= nr; p = p + 1) {
                            u[(p - 1) + (q - 1) * ldu] = work[(2 * n + n * nr + nr + p) - 1];
                        }
                    }
                    //
                }
                //
                //           Permute the rows of V using the (column) permutation from the
                //           first QRF. Also, scale the columns to make them unit in
                //           Euclidean norm. This applies to all cases.
                //
                temp1 = sqrt(castREAL(n)) * epsln;
                for (q = 1; q <= n; q = q + 1) {
                    for (p = 1; p <= n; p = p + 1) {
                        work[(2 * n + n * nr + nr + iwork[p - 1]) - 1] = v[(p - 1) + (q - 1) * ldv];
                    }
                    for (p = 1; p <= n; p = p + 1) {
                        v[(p - 1) + (q - 1) * ldv] = work[(2 * n + n * nr + nr + p) - 1];
                    }
                    xsc = one / Rnrm2(n, &v[(q - 1) * ldv], 1);
                    if ((xsc < (one - temp1)) || (xsc > (one + temp1))) {
                        Rscal(n, xsc, &v[(q - 1) * ldv], 1);
                    }
                }
                //           At this moment, V contains the right singular vectors of A.
                //           Next, assemble the left singular vector matrix U (M x N).
                if (nr < m) {
                    Rlaset("A", m - nr, nr, zero, zero, &u[((nr + 1) - 1)], ldu);
                    if (nr < n1) {
                        Rlaset("A", nr, n1 - nr, zero, zero, &u[((nr + 1) - 1) * ldu], ldu);
                        Rlaset("A", m - nr, n1 - nr, zero, one, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                    }
                }
                //
                //           The Q matrix from the first QRF is built into the left singular
                //           matrix U. This applies to all cases.
                //
                Rormqr("Left", "No_Tr", m, n1, n, a, lda, work, u, ldu, &work[(n + 1) - 1], lwork - n, ierr);
                //
                //           The columns of U are normalized. The cost is O(M*N) flops.
                temp1 = sqrt(castREAL(m)) * epsln;
                for (p = 1; p <= nr; p = p + 1) {
                    xsc = one / Rnrm2(m, &u[(p - 1) * ldu], 1);
                    if ((xsc < (one - temp1)) || (xsc > (one + temp1))) {
                        Rscal(m, xsc, &u[(p - 1) * ldu], 1);
                    }
                }
                //
                //           If the initial QRF is computed with row pivoting, the left
                //           singular vectors must be adjusted.
                //
                if (rowpiv) {
                    Rlaswp(n1, u, ldu, 1, m - 1, &iwork[(2 * n + 1) - 1], -1);
                }
                //
            } else {
                //
                //        .. the initial matrix A has almost orthogonal columns and
                //        the second QRF is not needed
                //
                Rlacpy("Upper", n, n, a, lda, &work[(n + 1) - 1], n);
                if (l2pert) {
                    xsc = sqrt(small);
                    for (p = 2; p <= n; p = p + 1) {
                        temp1 = xsc * work[(n + (p - 1) * n + p) - 1];
                        for (q = 1; q <= p - 1; q = q + 1) {
                            work[(n + (q - 1) * n + p) - 1] = -sign(temp1, work[(n + (p - 1) * n + q) - 1]);
                        }
                    }
                } else {
                    Rlaset("Lower", n - 1, n - 1, zero, zero, &work[(n + 2) - 1], n);
                }
                //
                Rgesvj("Upper", "U", "N", n, n, &work[(n + 1) - 1], n, sva, n, u, ldu, &work[(n + n * n + 1) - 1], lwork - n - n * n, info);
                //
                scalem = work[(n + n * n + 1) - 1];
                numrank = nint(work[(n + n * n + 2) - 1]);
                for (p = 1; p <= n; p = p + 1) {
                    Rcopy(n, &work[(n + (p - 1) * n + 1) - 1], 1, &u[(p - 1) * ldu], 1);
                    Rscal(n, sva[p - 1], &work[(n + (p - 1) * n + 1) - 1], 1);
                }
                //
                Rtrsm("Left", "Upper", "NoTrans", "No UD", n, n, one, a, lda, &work[(n + 1) - 1], n);
                for (p = 1; p <= n; p = p + 1) {
                    Rcopy(n, &work[(n + p) - 1], n, &v[(iwork[p - 1] - 1)], ldv);
                }
                temp1 = sqrt(castREAL(n)) * epsln;
                for (p = 1; p <= n; p = p + 1) {
                    xsc = one / Rnrm2(n, &v[(p - 1) * ldv], 1);
                    if ((xsc < (one - temp1)) || (xsc > (one + temp1))) {
                        Rscal(n, xsc, &v[(p - 1) * ldv], 1);
                    }
                }
                //
                //           Assemble the left singular vector matrix U (M x N).
                //
                if (n < m) {
                    Rlaset("A", m - n, n, zero, zero, &u[((n + 1) - 1)], ldu);
                    if (n < n1) {
                        Rlaset("A", n, n1 - n, zero, zero, &u[((n + 1) - 1) * ldu], ldu);
                        Rlaset("A", m - n, n1 - n, zero, one, &u[((n + 1) - 1) + ((n + 1) - 1) * ldu], ldu);
                    }
                }
                Rormqr("Left", "No Tr", m, n1, n, a, lda, work, u, ldu, &work[(n + 1) - 1], lwork - n, ierr);
                temp1 = sqrt(castREAL(m)) * epsln;
                for (p = 1; p <= n1; p = p + 1) {
                    xsc = one / Rnrm2(m, &u[(p - 1) * ldu], 1);
                    if ((xsc < (one - temp1)) || (xsc > (one + temp1))) {
                        Rscal(m, xsc, &u[(p - 1) * ldu], 1);
                    }
                }
                //
                if (rowpiv) {
                    Rlaswp(n1, u, ldu, 1, m - 1, &iwork[(2 * n + 1) - 1], -1);
                }
                //
            }
            //
            //        end of the  >> almost orthogonal case <<  in the full SVD
            //
        } else {
            //
            //        This branch deploys a preconditioned Jacobi SVD with explicitly
            //        accumulated rotations. It is included as optional, mainly for
            //        experimental purposes. It does perform well, and can also be used.
            //        In this implementation, this branch will be automatically activated
            //        if the  condition number sigma_max(A) / sigma_min(A) is predicted
            //        to be greater than the overflow threshold. This is because the
            //        a posteriori computation of the singular vectors assumes robust
            //        implementation of BLAS and some LAPACK procedures, capable of working
            //        in presence of extreme values. Since that is not always the case, ...
            //
            for (p = 1; p <= nr; p = p + 1) {
                Rcopy(n - p + 1, &a[(p - 1) + (p - 1) * lda], lda, &v[(p - 1) + (p - 1) * ldv], 1);
            }
            //
            if (l2pert) {
                xsc = sqrt(small / epsln);
                for (q = 1; q <= nr; q = q + 1) {
                    temp1 = xsc * abs(v[(q - 1) + (q - 1) * ldv]);
                    for (p = 1; p <= n; p = p + 1) {
                        if ((p > q) && (abs(v[(p - 1) + (q - 1) * ldv]) <= temp1) || (p < q)) {
                            v[(p - 1) + (q - 1) * ldv] = sign(temp1, v[(p - 1) + (q - 1) * ldv]);
                        }
                        if (p < q) {
                            v[(p - 1) + (q - 1) * ldv] = -v[(p - 1) + (q - 1) * ldv];
                        }
                    }
                }
            } else {
                Rlaset("U", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
            }
            //
            Rgeqrf(n, nr, v, ldv, &work[(n + 1) - 1], &work[(2 * n + 1) - 1], lwork - 2 * n, ierr);
            Rlacpy("L", n, nr, v, ldv, &work[(2 * n + 1) - 1], n);
            //
            for (p = 1; p <= nr; p = p + 1) {
                Rcopy(nr - p + 1, &v[(p - 1) + (p - 1) * ldv], ldv, &u[(p - 1) + (p - 1) * ldu], 1);
            }
            //
            if (l2pert) {
                xsc = sqrt(small / epsln);
                for (q = 2; q <= nr; q = q + 1) {
                    for (p = 1; p <= q - 1; p = p + 1) {
                        temp1 = xsc * min(abs(u[(p - 1) + (p - 1) * ldu]), abs(u[(q - 1) + (q - 1) * ldu]));
                        u[(p - 1) + (q - 1) * ldu] = -sign(temp1, u[(q - 1) + (p - 1) * ldu]);
                    }
                }
            } else {
                Rlaset("U", nr - 1, nr - 1, zero, zero, &u[(2 - 1) * ldu], ldu);
            }
            //
            Rgesvj("G", "U", "V", nr, nr, u, ldu, sva, n, v, ldv, &work[(2 * n + n * nr + 1) - 1], lwork - 2 * n - n * nr, info);
            scalem = work[(2 * n + n * nr + 1) - 1];
            numrank = nint(work[(2 * n + n * nr + 2) - 1]);
            //
            if (nr < n) {
                Rlaset("A", n - nr, nr, zero, zero, &v[((nr + 1) - 1)], ldv);
                Rlaset("A", nr, n - nr, zero, zero, &v[((nr + 1) - 1) * ldv], ldv);
                Rlaset("A", n - nr, n - nr, zero, one, &v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
            }
            //
            Rormqr("L", "N", n, n, nr, &work[(2 * n + 1) - 1], n, &work[(n + 1) - 1], v, ldv, &work[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
            //
            //           Permute the rows of V using the (column) permutation from the
            //           first QRF. Also, scale the columns to make them unit in
            //           Euclidean norm. This applies to all cases.
            //
            temp1 = sqrt(castREAL(n)) * epsln;
            for (q = 1; q <= n; q = q + 1) {
                for (p = 1; p <= n; p = p + 1) {
                    work[(2 * n + n * nr + nr + iwork[p - 1]) - 1] = v[(p - 1) + (q - 1) * ldv];
                }
                for (p = 1; p <= n; p = p + 1) {
                    v[(p - 1) + (q - 1) * ldv] = work[(2 * n + n * nr + nr + p) - 1];
                }
                xsc = one / Rnrm2(n, &v[(q - 1) * ldv], 1);
                if ((xsc < (one - temp1)) || (xsc > (one + temp1))) {
                    Rscal(n, xsc, &v[(q - 1) * ldv], 1);
                }
            }
            //
            //           At this moment, V contains the right singular vectors of A.
            //           Next, assemble the left singular vector matrix U (M x N).
            //
            if (nr < m) {
                Rlaset("A", m - nr, nr, zero, zero, &u[((nr + 1) - 1)], ldu);
                if (nr < n1) {
                    Rlaset("A", nr, n1 - nr, zero, zero, &u[((nr + 1) - 1) * ldu], ldu);
                    Rlaset("A", m - nr, n1 - nr, zero, one, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                }
            }
            //
            Rormqr("Left", "No Tr", m, n1, n, a, lda, work, u, ldu, &work[(n + 1) - 1], lwork - n, ierr);
            //
            if (rowpiv) {
                Rlaswp(n1, u, ldu, 1, m - 1, &iwork[(2 * n + 1) - 1], -1);
            }
            //
        }
        if (transp) {
            //           .. swap U and V because the procedure worked on A^t
            for (p = 1; p <= n; p = p + 1) {
                Rswap(n, &u[(p - 1) * ldu], 1, &v[(p - 1) * ldv], 1);
            }
        }
        //
    }
    //     end of the full SVD
    //
    //     Undo scaling, if necessary (and possible)
    //
    if (uscal2 <= (big / sva[1 - 1]) * uscal1) {
        Rlascl("G", 0, 0, uscal1, uscal2, nr, 1, sva, n, ierr);
        uscal1 = one;
        uscal2 = one;
    }
    //
    if (nr < n) {
        for (p = nr + 1; p <= n; p = p + 1) {
            sva[p - 1] = zero;
        }
    }
    //
    work[1 - 1] = uscal2 * scalem;
    work[2 - 1] = uscal1;
    if (errest) {
        work[3 - 1] = sconda;
    }
    if (lsvec && rsvec) {
        work[4 - 1] = condr1;
        work[5 - 1] = condr2;
    }
    if (l2tran) {
        work[6 - 1] = entra;
        work[7 - 1] = entrat;
    }
    //
    iwork[1 - 1] = nr;
    iwork[2 - 1] = numrank;
    iwork[3 - 1] = warning;
    //
    //     ..
    //     .. END OF Rgejsv
    //     ..
}
