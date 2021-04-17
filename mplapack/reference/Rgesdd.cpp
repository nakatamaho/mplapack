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

void Rgesdd(const char *jobz, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *s, REAL *u, INTEGER const ldu, REAL *vt, INTEGER const ldvt, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER &info) {
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    INTEGER minmn = min(m, n);
    bool wntqa = Mlsame(jobz, "A");
    bool wntqs = Mlsame(jobz, "S");
    bool wntqas = wntqa || wntqs;
    bool wntqo = Mlsame(jobz, "O");
    bool wntqn = Mlsame(jobz, "N");
    bool lquery = (lwork == -1);
    //
    if (!(wntqa || wntqs || wntqo || wntqn)) {
        info = -1;
    } else if (m < 0) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    } else if (ldu < 1 || (wntqas && ldu < m) || (wntqo && m < n && ldu < m)) {
        info = -8;
    } else if (ldvt < 1 || (wntqa && ldvt < n) || (wntqs && ldvt < minmn) || (wntqo && m >= n && ldvt < n)) {
        info = -10;
    }
    //
    //     Compute workspace
    //       Note: Comments in the code beginning "Workspace:" describe the
    //       minimal amount of workspace allocated at that point in the code,
    //       as well as the preferred amount for good performance.
    //       NB refers to the optimal block size for the immediately
    //       following subroutine, as returned by iMlaenv.
    //
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER bdspac = 0;
    INTEGER mnthr = 0;
    REAL dum[1];
    INTEGER ierr = 0;
    INTEGER lwork_Rgebrd_mn = 0;
    INTEGER lwork_Rgebrd_nn = 0;
    INTEGER lwork_Rgeqrf_mn = 0;
    INTEGER lwork_Rorgbr_q_nn = 0;
    INTEGER lwork_Rorgqr_mm = 0;
    INTEGER lwork_Rorgqr_mn = 0;
    INTEGER lwork_Rormbr_prt_nn = 0;
    INTEGER lwork_Rormbr_qln_nn = 0;
    INTEGER lwork_Rormbr_qln_mn = 0;
    INTEGER lwork_Rormbr_qln_mm = 0;
    INTEGER wrkbl = 0;
    INTEGER lwork_Rgebrd_mm = 0;
    INTEGER lwork_Rgelqf_mn = 0;
    INTEGER lwork_Rorglq_nn = 0;
    INTEGER lwork_Rorglq_mn = 0;
    INTEGER lwork_Rorgbr_p_mm = 0;
    INTEGER lwork_Rormbr_prt_mm = 0;
    INTEGER lwork_Rormbr_prt_mn = 0;
    if (info == 0) {
        minwrk = 1;
        maxwrk = 1;
        bdspac = 0;
        mnthr = int(minmn * 11.0 / 6.0);
        if (m >= n && minmn > 0) {
            //
            //           Compute space needed for Rbdsdc
            //
            if (wntqn) {
                //              Rbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6)
                //              keep 7*N for backwards compatibility.
                bdspac = 7 * n;
            } else {
                bdspac = 3 * n * n + 4 * n;
            }
            //
            //           Compute space preferred for each routine
            Rgebrd(m, n, &dum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rgebrd_mn = castINTEGER(dum[1 - 1]);
            //
            Rgebrd(n, n, &dum[1 - 1], n, &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rgebrd_nn = castINTEGER(dum[1 - 1]);
            //
            Rgeqrf(m, n, &dum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rgeqrf_mn = castINTEGER(dum[1 - 1]);
            //
            Rorgbr("Q", n, n, n, &dum[1 - 1], n, &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rorgbr_q_nn = castINTEGER(dum[1 - 1]);
            //
            Rorgqr(m, m, n, &dum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rorgqr_mm = castINTEGER(dum[1 - 1]);
            //
            Rorgqr(m, n, n, &dum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rorgqr_mn = castINTEGER(dum[1 - 1]);
            //
            Rormbr("P", "R", "T", n, n, n, &dum[1 - 1], n, &dum[1 - 1], &dum[1 - 1], n, &dum[1 - 1], -1, ierr);
            lwork_Rormbr_prt_nn = castINTEGER(dum[1 - 1]);
            //
            Rormbr("Q", "L", "N", n, n, n, &dum[1 - 1], n, &dum[1 - 1], &dum[1 - 1], n, &dum[1 - 1], -1, ierr);
            lwork_Rormbr_qln_nn = castINTEGER(dum[1 - 1]);
            //
            Rormbr("Q", "L", "N", m, n, n, &dum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], m, &dum[1 - 1], -1, ierr);
            lwork_Rormbr_qln_mn = castINTEGER(dum[1 - 1]);
            //
            Rormbr("Q", "L", "N", m, m, n, &dum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], m, &dum[1 - 1], -1, ierr);
            lwork_Rormbr_qln_mm = castINTEGER(dum[1 - 1]);
            //
            if (m >= mnthr) {
                if (wntqn) {
                    //
                    //                 Path 1 (M >> N, JOBZ='N')
                    //
                    wrkbl = n + lwork_Rgeqrf_mn;
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd_nn);
                    maxwrk = max(wrkbl, bdspac + n);
                    minwrk = bdspac + n;
                } else if (wntqo) {
                    //
                    //                 Path 2 (M >> N, JOBZ='O')
                    //
                    wrkbl = n + lwork_Rgeqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_mn);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd_nn);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_qln_nn);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_prt_nn);
                    wrkbl = max(wrkbl, 3 * n + bdspac);
                    maxwrk = wrkbl + 2 * n * n;
                    minwrk = bdspac + 2 * n * n + 3 * n;
                } else if (wntqs) {
                    //
                    //                 Path 3 (M >> N, JOBZ='S')
                    //
                    wrkbl = n + lwork_Rgeqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_mn);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd_nn);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_qln_nn);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_prt_nn);
                    wrkbl = max(wrkbl, 3 * n + bdspac);
                    maxwrk = wrkbl + n * n;
                    minwrk = bdspac + n * n + 3 * n;
                } else if (wntqa) {
                    //
                    //                 Path 4 (M >> N, JOBZ='A')
                    //
                    wrkbl = n + lwork_Rgeqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_mm);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd_nn);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_qln_nn);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_prt_nn);
                    wrkbl = max(wrkbl, 3 * n + bdspac);
                    maxwrk = wrkbl + n * n;
                    minwrk = n * n + max(3 * n + bdspac, n + m);
                }
            } else {
                //
                //              Path 5 (M >= N, but not much larger)
                //
                wrkbl = 3 * n + lwork_Rgebrd_mn;
                if (wntqn) {
                    //                 Path 5n (M >= N, jobz='N')
                    maxwrk = max(wrkbl, 3 * n + bdspac);
                    minwrk = 3 * n + max(m, bdspac);
                } else if (wntqo) {
                    //                 Path 5o (M >= N, jobz='O')
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_prt_nn);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_qln_mn);
                    wrkbl = max(wrkbl, 3 * n + bdspac);
                    maxwrk = wrkbl + m * n;
                    minwrk = 3 * n + max(m, n * n + bdspac);
                } else if (wntqs) {
                    //                 Path 5s (M >= N, jobz='S')
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_qln_mn);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_prt_nn);
                    maxwrk = max(wrkbl, 3 * n + bdspac);
                    minwrk = 3 * n + max(m, bdspac);
                } else if (wntqa) {
                    //                 Path 5a (M >= N, jobz='A')
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_qln_mm);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rormbr_prt_nn);
                    maxwrk = max(wrkbl, 3 * n + bdspac);
                    minwrk = 3 * n + max(m, bdspac);
                }
            }
        } else if (minmn > 0) {
            //
            //           Compute space needed for Rbdsdc
            //
            if (wntqn) {
                //              Rbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6)
                //              keep 7*N for backwards compatibility.
                bdspac = 7 * m;
            } else {
                bdspac = 3 * m * m + 4 * m;
            }
            //
            //           Compute space preferred for each routine
            Rgebrd(m, n, &dum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rgebrd_mn = castINTEGER(dum[1 - 1]);
            //
            Rgebrd(m, m, a, m, s, &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rgebrd_mm = castINTEGER(dum[1 - 1]);
            //
            Rgelqf(m, n, a, m, &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rgelqf_mn = castINTEGER(dum[1 - 1]);
            //
            Rorglq(n, n, m, &dum[1 - 1], n, &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rorglq_nn = castINTEGER(dum[1 - 1]);
            //
            Rorglq(m, n, m, a, m, &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rorglq_mn = castINTEGER(dum[1 - 1]);
            //
            Rorgbr("P", m, m, m, a, n, &dum[1 - 1], &dum[1 - 1], -1, ierr);
            lwork_Rorgbr_p_mm = castINTEGER(dum[1 - 1]);
            //
            Rormbr("P", "R", "T", m, m, m, &dum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], m, &dum[1 - 1], -1, ierr);
            lwork_Rormbr_prt_mm = castINTEGER(dum[1 - 1]);
            //
            Rormbr("P", "R", "T", m, n, m, &dum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], m, &dum[1 - 1], -1, ierr);
            lwork_Rormbr_prt_mn = castINTEGER(dum[1 - 1]);
            //
            Rormbr("P", "R", "T", n, n, m, &dum[1 - 1], n, &dum[1 - 1], &dum[1 - 1], n, &dum[1 - 1], -1, ierr);
            lwork_Rormbr_prt_nn = castINTEGER(dum[1 - 1]);
            //
            Rormbr("Q", "L", "N", m, m, m, &dum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], m, &dum[1 - 1], -1, ierr);
            lwork_Rormbr_qln_mm = castINTEGER(dum[1 - 1]);
            //
            if (n >= mnthr) {
                if (wntqn) {
                    //
                    //                 Path 1t (N >> M, JOBZ='N')
                    //
                    wrkbl = m + lwork_Rgelqf_mn;
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd_mm);
                    maxwrk = max(wrkbl, bdspac + m);
                    minwrk = bdspac + m;
                } else if (wntqo) {
                    //
                    //                 Path 2t (N >> M, JOBZ='O')
                    //
                    wrkbl = m + lwork_Rgelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_mn);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd_mm);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_qln_mm);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_prt_mm);
                    wrkbl = max(wrkbl, 3 * m + bdspac);
                    maxwrk = wrkbl + 2 * m * m;
                    minwrk = bdspac + 2 * m * m + 3 * m;
                } else if (wntqs) {
                    //
                    //                 Path 3t (N >> M, JOBZ='S')
                    //
                    wrkbl = m + lwork_Rgelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_mn);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd_mm);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_qln_mm);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_prt_mm);
                    wrkbl = max(wrkbl, 3 * m + bdspac);
                    maxwrk = wrkbl + m * m;
                    minwrk = bdspac + m * m + 3 * m;
                } else if (wntqa) {
                    //
                    //                 Path 4t (N >> M, JOBZ='A')
                    //
                    wrkbl = m + lwork_Rgelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_nn);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd_mm);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_qln_mm);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_prt_mm);
                    wrkbl = max(wrkbl, 3 * m + bdspac);
                    maxwrk = wrkbl + m * m;
                    minwrk = m * m + max(3 * m + bdspac, m + n);
                }
            } else {
                //
                //              Path 5t (N > M, but not much larger)
                //
                wrkbl = 3 * m + lwork_Rgebrd_mn;
                if (wntqn) {
                    //                 Path 5tn (N > M, jobz='N')
                    maxwrk = max(wrkbl, 3 * m + bdspac);
                    minwrk = 3 * m + max(n, bdspac);
                } else if (wntqo) {
                    //                 Path 5to (N > M, jobz='O')
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_qln_mm);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_prt_mn);
                    wrkbl = max(wrkbl, 3 * m + bdspac);
                    maxwrk = wrkbl + m * n;
                    minwrk = 3 * m + max(n, m * m + bdspac);
                } else if (wntqs) {
                    //                 Path 5ts (N > M, jobz='S')
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_qln_mm);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_prt_mn);
                    maxwrk = max(wrkbl, 3 * m + bdspac);
                    minwrk = 3 * m + max(n, bdspac);
                } else if (wntqa) {
                    //                 Path 5ta (N > M, jobz='A')
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_qln_mm);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rormbr_prt_nn);
                    maxwrk = max(wrkbl, 3 * m + bdspac);
                    minwrk = 3 * m + max(n, bdspac);
                }
            }
        }
        //
        maxwrk = max(maxwrk, minwrk);
        work[1 - 1] = maxwrk;
        //
        if (lwork < minwrk && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rgesdd", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     Get machine constants
    //
    REAL eps = Rlamch("P");
    REAL smlnum = sqrt(Rlamch("S")) / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    //
    //     Scale A if max element outside range [SMLNUM,BIGNUM]
    //
    REAL anrm = Rlange("M", m, n, a, lda, dum);
    if (Risnan(anrm)) {
        info = -4;
        return;
    }
    INTEGER iscl = 0;
    const REAL zero = 0.0;
    if (anrm > zero && anrm < smlnum) {
        iscl = 1;
        Rlascl("G", 0, 0, anrm, smlnum, m, n, a, lda, ierr);
    } else if (anrm > bignum) {
        iscl = 1;
        Rlascl("G", 0, 0, anrm, bignum, m, n, a, lda, ierr);
    }
    //
    INTEGER itau = 0;
    INTEGER nwork = 0;
    INTEGER ie = 0;
    INTEGER itauq = 0;
    INTEGER itaup = 0;
    INTEGER idum[1];
    INTEGER ir = 0;
    INTEGER ldwrkr = 0;
    INTEGER iu = 0;
    INTEGER i = 0;
    INTEGER chunk = 0;
    INTEGER ldwrku = 0;
    INTEGER ivt = 0;
    INTEGER il = 0;
    INTEGER ldwrkl = 0;
    INTEGER blk = 0;
    INTEGER ldwkvt = 0;
    if (m >= n) {
        //
        //        A has at least as many rows as columns. If A has sufficiently
        //        more rows than columns, first reduce using the QR
        //        decomposition (if sufficient workspace available)
        //
        if (m >= mnthr) {
            //
            if (wntqn) {
                //
                //              Path 1 (M >> N, JOBZ='N')
                //              No singular vectors to be computed
                //
                itau = 1;
                nwork = itau + n;
                //
                //              Compute A=Q*R
                //              Workspace: need   N [tau] + N    [work]
                //              Workspace: prefer N [tau] + N*NB [work]
                //
                Rgeqrf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Zero out below R
                //
                Rlaset("L", n - 1, n - 1, zero, zero, &a[(2 - 1)], lda);
                ie = 1;
                itauq = ie + n;
                itaup = itauq + n;
                nwork = itaup + n;
                //
                //              Bidiagonalize R in A
                //              Workspace: need   3*N [e, tauq, taup] + N      [work]
                //              Workspace: prefer 3*N [e, tauq, taup] + 2*N*NB [work]
                //
                Rgebrd(n, n, a, lda, s, &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                nwork = ie + n;
                //
                //              Perform bidiagonal SVD, computing singular values only
                //              Workspace: need   N [e] + BDSPAC
                //
                Rbdsdc("U", "N", n, s, &work[ie - 1], dum, 1, dum, 1, dum, idum, &work[nwork - 1], iwork, info);
                //
            } else if (wntqo) {
                //
                //              Path 2 (M >> N, JOBZ = 'O')
                //              N left singular vectors to be overwritten on A and
                //              N right singular vectors to be computed in VT
                //
                ir = 1;
                //
                //              WORK(IR) is LDWRKR by N
                //
                if (lwork >= lda * n + n * n + 3 * n + bdspac) {
                    ldwrkr = lda;
                } else {
                    ldwrkr = (lwork - n * n - 3 * n - bdspac) / n;
                }
                itau = ir + ldwrkr * n;
                nwork = itau + n;
                //
                //              Compute A=Q*R
                //              Workspace: need   N*N [R] + N [tau] + N    [work]
                //              Workspace: prefer N*N [R] + N [tau] + N*NB [work]
                //
                Rgeqrf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy R to WORK(IR), zeroing out below it
                //
                Rlacpy("U", n, n, a, lda, &work[ir - 1], ldwrkr);
                Rlaset("L", n - 1, n - 1, zero, zero, &work[(ir + 1) - 1], ldwrkr);
                //
                //              Generate Q in A
                //              Workspace: need   N*N [R] + N [tau] + N    [work]
                //              Workspace: prefer N*N [R] + N [tau] + N*NB [work]
                //
                Rorgqr(m, n, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                ie = itau;
                itauq = ie + n;
                itaup = itauq + n;
                nwork = itaup + n;
                //
                //              Bidiagonalize R in WORK(IR)
                //              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work]
                //              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work]
                //
                Rgebrd(n, n, &work[ir - 1], ldwrkr, s, &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              WORK(IU) is N by N
                //
                iu = nwork;
                nwork = iu + n * n;
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in WORK(IU) and computing right
                //              singular vectors of bidiagonal matrix in VT
                //              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + BDSPAC
                //
                Rbdsdc("U", "I", n, s, &work[ie - 1], &work[iu - 1], n, vt, ldvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Overwrite WORK(IU) by left singular vectors of R
                //              and VT by right singular vectors of R
                //              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N    [work]
                //              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N*NB [work]
                //
                Rormbr("Q", "L", "N", n, n, n, &work[ir - 1], ldwrkr, &work[itauq - 1], &work[iu - 1], n, &work[nwork - 1], lwork - nwork + 1, ierr);
                Rormbr("P", "R", "T", n, n, n, &work[ir - 1], ldwrkr, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply Q in A by left singular vectors of R in
                //              WORK(IU), storing result in WORK(IR) and copying to A
                //              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U]
                //              Workspace: prefer M*N [R] + 3*N [e, tauq, taup] + N*N [U]
                //
                for (i = 1; i <= m; i = i + ldwrkr) {
                    chunk = min(m - i + 1, ldwrkr);
                    Rgemm("N", "N", chunk, n, n, one, &a[(i - 1)], lda, &work[iu - 1], n, zero, &work[ir - 1], ldwrkr);
                    Rlacpy("F", chunk, n, &work[ir - 1], ldwrkr, &a[(i - 1)], lda);
                }
                //
            } else if (wntqs) {
                //
                //              Path 3 (M >> N, JOBZ='S')
                //              N left singular vectors to be computed in U and
                //              N right singular vectors to be computed in VT
                //
                ir = 1;
                //
                //              WORK(IR) is N by N
                //
                ldwrkr = n;
                itau = ir + ldwrkr * n;
                nwork = itau + n;
                //
                //              Compute A=Q*R
                //              Workspace: need   N*N [R] + N [tau] + N    [work]
                //              Workspace: prefer N*N [R] + N [tau] + N*NB [work]
                //
                Rgeqrf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy R to WORK(IR), zeroing out below it
                //
                Rlacpy("U", n, n, a, lda, &work[ir - 1], ldwrkr);
                Rlaset("L", n - 1, n - 1, zero, zero, &work[(ir + 1) - 1], ldwrkr);
                //
                //              Generate Q in A
                //              Workspace: need   N*N [R] + N [tau] + N    [work]
                //              Workspace: prefer N*N [R] + N [tau] + N*NB [work]
                //
                Rorgqr(m, n, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                ie = itau;
                itauq = ie + n;
                itaup = itauq + n;
                nwork = itaup + n;
                //
                //              Bidiagonalize R in WORK(IR)
                //              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work]
                //              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work]
                //
                Rgebrd(n, n, &work[ir - 1], ldwrkr, s, &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagoal matrix in U and computing right singular
                //              vectors of bidiagonal matrix in VT
                //              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + BDSPAC
                //
                Rbdsdc("U", "I", n, s, &work[ie - 1], u, ldu, vt, ldvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Overwrite U by left singular vectors of R and VT
                //              by right singular vectors of R
                //              Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N    [work]
                //              Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*NB [work]
                //
                Rormbr("Q", "L", "N", n, n, n, &work[ir - 1], ldwrkr, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                Rormbr("P", "R", "T", n, n, n, &work[ir - 1], ldwrkr, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply Q in A by left singular vectors of R in
                //              WORK(IR), storing result in U
                //              Workspace: need   N*N [R]
                //
                Rlacpy("F", n, n, u, ldu, &work[ir - 1], ldwrkr);
                Rgemm("N", "N", m, n, n, one, a, lda, &work[ir - 1], ldwrkr, zero, u, ldu);
                //
            } else if (wntqa) {
                //
                //              Path 4 (M >> N, JOBZ='A')
                //              M left singular vectors to be computed in U and
                //              N right singular vectors to be computed in VT
                //
                iu = 1;
                //
                //              WORK(IU) is N by N
                //
                ldwrku = n;
                itau = iu + ldwrku * n;
                nwork = itau + n;
                //
                //              Compute A=Q*R, copying result to U
                //              Workspace: need   N*N [U] + N [tau] + N    [work]
                //              Workspace: prefer N*N [U] + N [tau] + N*NB [work]
                //
                Rgeqrf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                Rlacpy("L", m, n, a, lda, u, ldu);
                //
                //              Generate Q in U
                //              Workspace: need   N*N [U] + N [tau] + M    [work]
                //              Workspace: prefer N*N [U] + N [tau] + M*NB [work]
                Rorgqr(m, m, n, u, ldu, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Produce R in A, zeroing out other entries
                //
                Rlaset("L", n - 1, n - 1, zero, zero, &a[(2 - 1)], lda);
                ie = itau;
                itauq = ie + n;
                itaup = itauq + n;
                nwork = itaup + n;
                //
                //              Bidiagonalize R in A
                //              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N      [work]
                //              Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + 2*N*NB [work]
                //
                Rgebrd(n, n, a, lda, s, &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in WORK(IU) and computing right
                //              singular vectors of bidiagonal matrix in VT
                //              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + BDSPAC
                //
                Rbdsdc("U", "I", n, s, &work[ie - 1], &work[iu - 1], n, vt, ldvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Overwrite WORK(IU) by left singular vectors of R and VT
                //              by right singular vectors of R
                //              Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N    [work]
                //              Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + N*NB [work]
                //
                Rormbr("Q", "L", "N", n, n, n, a, lda, &work[itauq - 1], &work[iu - 1], ldwrku, &work[nwork - 1], lwork - nwork + 1, ierr);
                Rormbr("P", "R", "T", n, n, n, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply Q in U by left singular vectors of R in
                //              WORK(IU), storing result in A
                //              Workspace: need   N*N [U]
                //
                Rgemm("N", "N", m, n, n, one, u, ldu, &work[iu - 1], ldwrku, zero, a, lda);
                //
                //              Copy left singular vectors of A from A to U
                //
                Rlacpy("F", m, n, a, lda, u, ldu);
                //
            }
            //
        } else {
            //
            //           M .LT. MNTHR
            //
            //           Path 5 (M >= N, but not much larger)
            //           Reduce to bidiagonal form without QR decomposition
            //
            ie = 1;
            itauq = ie + n;
            itaup = itauq + n;
            nwork = itaup + n;
            //
            //           Bidiagonalize A
            //           Workspace: need   3*N [e, tauq, taup] + M        [work]
            //           Workspace: prefer 3*N [e, tauq, taup] + (M+N)*NB [work]
            //
            Rgebrd(m, n, a, lda, s, &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
            if (wntqn) {
                //
                //              Path 5n (M >= N, JOBZ='N')
                //              Perform bidiagonal SVD, only computing singular values
                //              Workspace: need   3*N [e, tauq, taup] + BDSPAC
                //
                Rbdsdc("U", "N", n, s, &work[ie - 1], dum, 1, dum, 1, dum, idum, &work[nwork - 1], iwork, info);
            } else if (wntqo) {
                //              Path 5o (M >= N, JOBZ='O')
                iu = nwork;
                if (lwork >= m * n + 3 * n + bdspac) {
                    //
                    //                 WORK( IU ) is M by N
                    //
                    ldwrku = m;
                    nwork = iu + ldwrku * n;
                    Rlaset("F", m, n, zero, zero, &work[iu - 1], ldwrku);
                    //                 IR is unused; silence compile warnings
                    ir = -1;
                } else {
                    //
                    //                 WORK( IU ) is N by N
                    //
                    ldwrku = n;
                    nwork = iu + ldwrku * n;
                    //
                    //                 WORK(IR) is LDWRKR by N
                    //
                    ir = nwork;
                    ldwrkr = (lwork - n * n - 3 * n) / n;
                }
                nwork = iu + ldwrku * n;
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in WORK(IU) and computing right
                //              singular vectors of bidiagonal matrix in VT
                //              Workspace: need   3*N [e, tauq, taup] + N*N [U] + BDSPAC
                //
                Rbdsdc("U", "I", n, s, &work[ie - 1], &work[iu - 1], ldwrku, vt, ldvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Overwrite VT by right singular vectors of A
                //              Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work]
                //              Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work]
                //
                Rormbr("P", "R", "T", n, n, n, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                if (lwork >= m * n + 3 * n + bdspac) {
                    //
                    //                 Path 5o-fast
                    //                 Overwrite WORK(IU) by left singular vectors of A
                    //                 Workspace: need   3*N [e, tauq, taup] + M*N [U] + N    [work]
                    //                 Workspace: prefer 3*N [e, tauq, taup] + M*N [U] + N*NB [work]
                    //
                    Rormbr("Q", "L", "N", m, n, n, a, lda, &work[itauq - 1], &work[iu - 1], ldwrku, &work[nwork - 1], lwork - nwork + 1, ierr);
                    //
                    //                 Copy left singular vectors of A from WORK(IU) to A
                    //
                    Rlacpy("F", m, n, &work[iu - 1], ldwrku, a, lda);
                } else {
                    //
                    //                 Path 5o-slow
                    //                 Generate Q in A
                    //                 Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work]
                    //                 Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work]
                    //
                    Rorgbr("Q", m, n, n, a, lda, &work[itauq - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                    //
                    //                 Multiply Q in A by left singular vectors of
                    //                 bidiagonal matrix in WORK(IU), storing result in
                    //                 WORK(IR) and copying to A
                    //                 Workspace: need   3*N [e, tauq, taup] + N*N [U] + NB*N [R]
                    //                 Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + M*N  [R]
                    //
                    for (i = 1; i <= m; i = i + ldwrkr) {
                        chunk = min(m - i + 1, ldwrkr);
                        Rgemm("N", "N", chunk, n, n, one, &a[(i - 1)], lda, &work[iu - 1], ldwrku, zero, &work[ir - 1], ldwrkr);
                        Rlacpy("F", chunk, n, &work[ir - 1], ldwrkr, &a[(i - 1)], lda);
                    }
                }
                //
            } else if (wntqs) {
                //
                //              Path 5s (M >= N, JOBZ='S')
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in U and computing right singular
                //              vectors of bidiagonal matrix in VT
                //              Workspace: need   3*N [e, tauq, taup] + BDSPAC
                //
                Rlaset("F", m, n, zero, zero, u, ldu);
                Rbdsdc("U", "I", n, s, &work[ie - 1], u, ldu, vt, ldvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Overwrite U by left singular vectors of A and VT
                //              by right singular vectors of A
                //              Workspace: need   3*N [e, tauq, taup] + N    [work]
                //              Workspace: prefer 3*N [e, tauq, taup] + N*NB [work]
                //
                Rormbr("Q", "L", "N", m, n, n, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                Rormbr("P", "R", "T", n, n, n, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
            } else if (wntqa) {
                //
                //              Path 5a (M >= N, JOBZ='A')
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in U and computing right singular
                //              vectors of bidiagonal matrix in VT
                //              Workspace: need   3*N [e, tauq, taup] + BDSPAC
                //
                Rlaset("F", m, m, zero, zero, u, ldu);
                Rbdsdc("U", "I", n, s, &work[ie - 1], u, ldu, vt, ldvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Set the right corner of U to identity matrix
                //
                if (m > n) {
                    Rlaset("F", m - n, m - n, zero, one, &u[((n + 1) - 1) + ((n + 1) - 1) * ldu], ldu);
                }
                //
                //              Overwrite U by left singular vectors of A and VT
                //              by right singular vectors of A
                //              Workspace: need   3*N [e, tauq, taup] + M    [work]
                //              Workspace: prefer 3*N [e, tauq, taup] + M*NB [work]
                //
                Rormbr("Q", "L", "N", m, m, n, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                Rormbr("P", "R", "T", n, n, m, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
            }
            //
        }
        //
    } else {
        //
        //        A has more columns than rows. If A has sufficiently more
        //        columns than rows, first reduce using the LQ decomposition (if
        //        sufficient workspace available)
        //
        if (n >= mnthr) {
            //
            if (wntqn) {
                //
                //              Path 1t (N >> M, JOBZ='N')
                //              No singular vectors to be computed
                //
                itau = 1;
                nwork = itau + m;
                //
                //              Compute A=L*Q
                //              Workspace: need   M [tau] + M [work]
                //              Workspace: prefer M [tau] + M*NB [work]
                //
                Rgelqf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Zero out above L
                //
                Rlaset("U", m - 1, m - 1, zero, zero, &a[(2 - 1) * lda], lda);
                ie = 1;
                itauq = ie + m;
                itaup = itauq + m;
                nwork = itaup + m;
                //
                //              Bidiagonalize L in A
                //              Workspace: need   3*M [e, tauq, taup] + M      [work]
                //              Workspace: prefer 3*M [e, tauq, taup] + 2*M*NB [work]
                //
                Rgebrd(m, m, a, lda, s, &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                nwork = ie + m;
                //
                //              Perform bidiagonal SVD, computing singular values only
                //              Workspace: need   M [e] + BDSPAC
                //
                Rbdsdc("U", "N", m, s, &work[ie - 1], dum, 1, dum, 1, dum, idum, &work[nwork - 1], iwork, info);
                //
            } else if (wntqo) {
                //
                //              Path 2t (N >> M, JOBZ='O')
                //              M right singular vectors to be overwritten on A and
                //              M left singular vectors to be computed in U
                //
                ivt = 1;
                //
                //              WORK(IVT) is M by M
                //              WORK(IL)  is M by M; it is later resized to M by chunk for gemm
                //
                il = ivt + m * m;
                if (lwork >= m * n + m * m + 3 * m + bdspac) {
                    ldwrkl = m;
                    chunk = n;
                } else {
                    ldwrkl = m;
                    chunk = (lwork - m * m) / m;
                }
                itau = il + ldwrkl * m;
                nwork = itau + m;
                //
                //              Compute A=L*Q
                //              Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
                //              Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
                //
                Rgelqf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy L to WORK(IL), zeroing about above it
                //
                Rlacpy("L", m, m, a, lda, &work[il - 1], ldwrkl);
                Rlaset("U", m - 1, m - 1, zero, zero, &work[(il + ldwrkl) - 1], ldwrkl);
                //
                //              Generate Q in A
                //              Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
                //              Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
                //
                Rorglq(m, n, m, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                ie = itau;
                itauq = ie + m;
                itaup = itauq + m;
                nwork = itaup + m;
                //
                //              Bidiagonalize L in WORK(IL)
                //              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M      [work]
                //              Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work]
                //
                Rgebrd(m, m, &work[il - 1], ldwrkl, s, &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in U, and computing right singular
                //              vectors of bidiagonal matrix in WORK(IVT)
                //              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + BDSPAC
                //
                Rbdsdc("U", "I", m, s, &work[ie - 1], u, ldu, &work[ivt - 1], m, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Overwrite U by left singular vectors of L and WORK(IVT)
                //              by right singular vectors of L
                //              Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M    [work]
                //              Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M*NB [work]
                //
                Rormbr("Q", "L", "N", m, m, m, &work[il - 1], ldwrkl, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                Rormbr("P", "R", "T", m, m, m, &work[il - 1], ldwrkl, &work[itaup - 1], &work[ivt - 1], m, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply right singular vectors of L in WORK(IVT) by Q
                //              in A, storing result in WORK(IL) and copying to A
                //              Workspace: need   M*M [VT] + M*M [L]
                //              Workspace: prefer M*M [VT] + M*N [L]
                //              At this point, L is resized as M by chunk.
                //
                for (i = 1; i <= n; i = i + chunk) {
                    blk = min(n - i + 1, chunk);
                    Rgemm("N", "N", m, blk, m, one, &work[ivt - 1], m, &a[(i - 1) * lda], lda, zero, &work[il - 1], ldwrkl);
                    Rlacpy("F", m, blk, &work[il - 1], ldwrkl, &a[(i - 1) * lda], lda);
                }
                //
            } else if (wntqs) {
                //
                //              Path 3t (N >> M, JOBZ='S')
                //              M right singular vectors to be computed in VT and
                //              M left singular vectors to be computed in U
                //
                il = 1;
                //
                //              WORK(IL) is M by M
                //
                ldwrkl = m;
                itau = il + ldwrkl * m;
                nwork = itau + m;
                //
                //              Compute A=L*Q
                //              Workspace: need   M*M [L] + M [tau] + M    [work]
                //              Workspace: prefer M*M [L] + M [tau] + M*NB [work]
                //
                Rgelqf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy L to WORK(IL), zeroing out above it
                //
                Rlacpy("L", m, m, a, lda, &work[il - 1], ldwrkl);
                Rlaset("U", m - 1, m - 1, zero, zero, &work[(il + ldwrkl) - 1], ldwrkl);
                //
                //              Generate Q in A
                //              Workspace: need   M*M [L] + M [tau] + M    [work]
                //              Workspace: prefer M*M [L] + M [tau] + M*NB [work]
                //
                Rorglq(m, n, m, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                ie = itau;
                itauq = ie + m;
                itaup = itauq + m;
                nwork = itaup + m;
                //
                //              Bidiagonalize L in WORK(IU).
                //              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M      [work]
                //              Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work]
                //
                Rgebrd(m, m, &work[il - 1], ldwrkl, s, &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in U and computing right singular
                //              vectors of bidiagonal matrix in VT
                //              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + BDSPAC
                //
                Rbdsdc("U", "I", m, s, &work[ie - 1], u, ldu, vt, ldvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Overwrite U by left singular vectors of L and VT
                //              by right singular vectors of L
                //              Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M    [work]
                //              Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + M*NB [work]
                //
                Rormbr("Q", "L", "N", m, m, m, &work[il - 1], ldwrkl, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                Rormbr("P", "R", "T", m, m, m, &work[il - 1], ldwrkl, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply right singular vectors of L in WORK(IL) by
                //              Q in A, storing result in VT
                //              Workspace: need   M*M [L]
                //
                Rlacpy("F", m, m, vt, ldvt, &work[il - 1], ldwrkl);
                Rgemm("N", "N", m, n, m, one, &work[il - 1], ldwrkl, a, lda, zero, vt, ldvt);
                //
            } else if (wntqa) {
                //
                //              Path 4t (N >> M, JOBZ='A')
                //              N right singular vectors to be computed in VT and
                //              M left singular vectors to be computed in U
                //
                ivt = 1;
                //
                //              WORK(IVT) is M by M
                //
                ldwkvt = m;
                itau = ivt + ldwkvt * m;
                nwork = itau + m;
                //
                //              Compute A=L*Q, copying result to VT
                //              Workspace: need   M*M [VT] + M [tau] + M    [work]
                //              Workspace: prefer M*M [VT] + M [tau] + M*NB [work]
                //
                Rgelqf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                Rlacpy("U", m, n, a, lda, vt, ldvt);
                //
                //              Generate Q in VT
                //              Workspace: need   M*M [VT] + M [tau] + N    [work]
                //              Workspace: prefer M*M [VT] + M [tau] + N*NB [work]
                //
                Rorglq(n, n, m, vt, ldvt, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Produce L in A, zeroing out other entries
                //
                Rlaset("U", m - 1, m - 1, zero, zero, &a[(2 - 1) * lda], lda);
                ie = itau;
                itauq = ie + m;
                itaup = itauq + m;
                nwork = itaup + m;
                //
                //              Bidiagonalize L in A
                //              Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + M      [work]
                //              Workspace: prefer M*M [VT] + 3*M [e, tauq, taup] + 2*M*NB [work]
                //
                Rgebrd(m, m, a, lda, s, &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in U and computing right singular
                //              vectors of bidiagonal matrix in WORK(IVT)
                //              Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + BDSPAC
                //
                Rbdsdc("U", "I", m, s, &work[ie - 1], u, ldu, &work[ivt - 1], ldwkvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Overwrite U by left singular vectors of L and WORK(IVT)
                //              by right singular vectors of L
                //              Workspace: need   M*M [VT] + 3*M [e, tauq, taup]+ M    [work]
                //              Workspace: prefer M*M [VT] + 3*M [e, tauq, taup]+ M*NB [work]
                //
                Rormbr("Q", "L", "N", m, m, m, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                Rormbr("P", "R", "T", m, m, m, a, lda, &work[itaup - 1], &work[ivt - 1], ldwkvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply right singular vectors of L in WORK(IVT) by
                //              Q in VT, storing result in A
                //              Workspace: need   M*M [VT]
                //
                Rgemm("N", "N", m, n, m, one, &work[ivt - 1], ldwkvt, vt, ldvt, zero, a, lda);
                //
                //              Copy right singular vectors of A from A to VT
                //
                Rlacpy("F", m, n, a, lda, vt, ldvt);
                //
            }
            //
        } else {
            //
            //           N .LT. MNTHR
            //
            //           Path 5t (N > M, but not much larger)
            //           Reduce to bidiagonal form without LQ decomposition
            //
            ie = 1;
            itauq = ie + m;
            itaup = itauq + m;
            nwork = itaup + m;
            //
            //           Bidiagonalize A
            //           Workspace: need   3*M [e, tauq, taup] + N        [work]
            //           Workspace: prefer 3*M [e, tauq, taup] + (M+N)*NB [work]
            //
            Rgebrd(m, n, a, lda, s, &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
            if (wntqn) {
                //
                //              Path 5tn (N > M, JOBZ='N')
                //              Perform bidiagonal SVD, only computing singular values
                //              Workspace: need   3*M [e, tauq, taup] + BDSPAC
                //
                Rbdsdc("L", "N", m, s, &work[ie - 1], dum, 1, dum, 1, dum, idum, &work[nwork - 1], iwork, info);
            } else if (wntqo) {
                //              Path 5to (N > M, JOBZ='O')
                ldwkvt = m;
                ivt = nwork;
                if (lwork >= m * n + 3 * m + bdspac) {
                    //
                    //                 WORK( IVT ) is M by N
                    //
                    Rlaset("F", m, n, zero, zero, &work[ivt - 1], ldwkvt);
                    nwork = ivt + ldwkvt * n;
                    //                 IL is unused; silence compile warnings
                    il = -1;
                } else {
                    //
                    //                 WORK( IVT ) is M by M
                    //
                    nwork = ivt + ldwkvt * m;
                    il = nwork;
                    //
                    //                 WORK(IL) is M by CHUNK
                    //
                    chunk = (lwork - m * m - 3 * m) / m;
                }
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in U and computing right singular
                //              vectors of bidiagonal matrix in WORK(IVT)
                //              Workspace: need   3*M [e, tauq, taup] + M*M [VT] + BDSPAC
                //
                Rbdsdc("L", "I", m, s, &work[ie - 1], u, ldu, &work[ivt - 1], ldwkvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Overwrite U by left singular vectors of A
                //              Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work]
                //              Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work]
                //
                Rormbr("Q", "L", "N", m, m, n, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                if (lwork >= m * n + 3 * m + bdspac) {
                    //
                    //                 Path 5to-fast
                    //                 Overwrite WORK(IVT) by left singular vectors of A
                    //                 Workspace: need   3*M [e, tauq, taup] + M*N [VT] + M    [work]
                    //                 Workspace: prefer 3*M [e, tauq, taup] + M*N [VT] + M*NB [work]
                    //
                    Rormbr("P", "R", "T", m, n, m, a, lda, &work[itaup - 1], &work[ivt - 1], ldwkvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                    //
                    //                 Copy right singular vectors of A from WORK(IVT) to A
                    //
                    Rlacpy("F", m, n, &work[ivt - 1], ldwkvt, a, lda);
                } else {
                    //
                    //                 Path 5to-slow
                    //                 Generate P**T in A
                    //                 Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work]
                    //                 Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work]
                    //
                    Rorgbr("P", m, n, m, a, lda, &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                    //
                    //                 Multiply Q in A by right singular vectors of
                    //                 bidiagonal matrix in WORK(IVT), storing result in
                    //                 WORK(IL) and copying to A
                    //                 Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M*NB [L]
                    //                 Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*N  [L]
                    //
                    for (i = 1; i <= n; i = i + chunk) {
                        blk = min(n - i + 1, chunk);
                        Rgemm("N", "N", m, blk, m, one, &work[ivt - 1], ldwkvt, &a[(i - 1) * lda], lda, zero, &work[il - 1], m);
                        Rlacpy("F", m, blk, &work[il - 1], m, &a[(i - 1) * lda], lda);
                    }
                }
            } else if (wntqs) {
                //
                //              Path 5ts (N > M, JOBZ='S')
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in U and computing right singular
                //              vectors of bidiagonal matrix in VT
                //              Workspace: need   3*M [e, tauq, taup] + BDSPAC
                //
                Rlaset("F", m, n, zero, zero, vt, ldvt);
                Rbdsdc("L", "I", m, s, &work[ie - 1], u, ldu, vt, ldvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Overwrite U by left singular vectors of A and VT
                //              by right singular vectors of A
                //              Workspace: need   3*M [e, tauq, taup] + M    [work]
                //              Workspace: prefer 3*M [e, tauq, taup] + M*NB [work]
                //
                Rormbr("Q", "L", "N", m, m, n, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                Rormbr("P", "R", "T", m, n, m, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
            } else if (wntqa) {
                //
                //              Path 5ta (N > M, JOBZ='A')
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in U and computing right singular
                //              vectors of bidiagonal matrix in VT
                //              Workspace: need   3*M [e, tauq, taup] + BDSPAC
                //
                Rlaset("F", n, n, zero, zero, vt, ldvt);
                Rbdsdc("L", "I", m, s, &work[ie - 1], u, ldu, vt, ldvt, dum, idum, &work[nwork - 1], iwork, info);
                //
                //              Set the right corner of VT to identity matrix
                //
                if (n > m) {
                    Rlaset("F", n - m, n - m, zero, one, &vt[((m + 1) - 1) + ((m + 1) - 1) * ldvt], ldvt);
                }
                //
                //              Overwrite U by left singular vectors of A and VT
                //              by right singular vectors of A
                //              Workspace: need   3*M [e, tauq, taup] + N    [work]
                //              Workspace: prefer 3*M [e, tauq, taup] + N*NB [work]
                //
                Rormbr("Q", "L", "N", m, m, n, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                Rormbr("P", "R", "T", n, n, m, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
            }
            //
        }
        //
    }
    //
    //     Undo scaling if necessary
    //
    if (iscl == 1) {
        if (anrm > bignum) {
            Rlascl("G", 0, 0, bignum, anrm, minmn, 1, s, minmn, ierr);
        }
        if (anrm < smlnum) {
            Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, s, minmn, ierr);
        }
    }
    //
    //     Return optimal workspace in WORK(1)
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Rgesdd
    //
}
