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

void Cgesdd(const char *jobz, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL *s, COMPLEX *u, INTEGER const ldu, COMPLEX *vt, INTEGER const ldvt, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER *iwork, INTEGER &info) {
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
    INTEGER mnthr1 = int(minmn * 17.0 / 9.0);
    INTEGER mnthr2 = int(minmn * 5.0 / 3.0);
    bool wntqa = Mlsame(jobz, "A");
    bool wntqs = Mlsame(jobz, "S");
    bool wntqas = wntqa || wntqs;
    bool wntqo = Mlsame(jobz, "O");
    bool wntqn = Mlsame(jobz, "N");
    bool lquery = (lwork == -1);
    INTEGER minwrk = 1;
    INTEGER maxwrk = 1;
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
    //       CWorkspace refers to complex workspace, and RWorkspace to
    //       real workspace. NB refers to the optimal block size for the
    //       immediately following subroutine, as returned by iMlaenv.)
    //
    COMPLEX cdum[1];
    REAL dum[1];
    INTEGER ierr = 0;
    INTEGER lwork_Cgebrd_mn = 0;
    INTEGER lwork_Cgebrd_nn = 0;
    INTEGER lwork_Cgeqrf_mn = 0;
    INTEGER lwork_Cungbr_p_nn = 0;
    INTEGER lwork_Cungbr_q_mm = 0;
    INTEGER lwork_Cungbr_q_mn = 0;
    INTEGER lwork_Cungqr_mm = 0;
    INTEGER lwork_Cungqr_mn = 0;
    INTEGER lwork_Cunmbr_prc_nn = 0;
    INTEGER lwork_Cunmbr_qln_mm = 0;
    INTEGER lwork_Cunmbr_qln_mn = 0;
    INTEGER lwork_Cunmbr_qln_nn = 0;
    INTEGER wrkbl = 0;
    INTEGER lwork_Cgebrd_mm = 0;
    INTEGER lwork_Cgelqf_mn = 0;
    INTEGER lwork_Cungbr_p_mn = 0;
    INTEGER lwork_Cunglq_mn = 0;
    INTEGER lwork_Cunglq_nn = 0;
    INTEGER lwork_Cunmbr_prc_mm = 0;
    INTEGER lwork_Cunmbr_prc_mn = 0;
    if (info == 0) {
        minwrk = 1;
        maxwrk = 1;
        if (m >= n && minmn > 0) {
            //
            //           There is no complex work space needed for bidiagonal SVD
            //           The real work space needed for bidiagonal SVD (Rbdsdc) is
            //           BDSPAC = 3*N*N + 4*N for singular values and vectors;
            //           BDSPAC = 4*N         for singular values only;
            //           not including e, RU, and RVT matrices.
            //
            //           Compute space preferred for each routine
            Cgebrd(m, n, &cdum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cgebrd_mn = castINTEGER(cdum[1 - 1].real());
            //
            Cgebrd(n, n, &cdum[1 - 1], n, &dum[1 - 1], &dum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cgebrd_nn = castINTEGER(cdum[1 - 1].real());
            //
            Cgeqrf(m, n, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cgeqrf_mn = castINTEGER(cdum[1 - 1].real());
            //
            Cungbr("P", n, n, n, &cdum[1 - 1], n, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungbr_p_nn = castINTEGER(cdum[1 - 1].real());
            //
            Cungbr("Q", m, m, n, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungbr_q_mm = castINTEGER(cdum[1 - 1].real());
            //
            Cungbr("Q", m, n, n, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungbr_q_mn = castINTEGER(cdum[1 - 1].real());
            //
            Cungqr(m, m, n, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungqr_mm = castINTEGER(cdum[1 - 1].real());
            //
            Cungqr(m, n, n, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungqr_mn = castINTEGER(cdum[1 - 1].real());
            //
            Cunmbr("P", "R", "C", n, n, n, &cdum[1 - 1], n, &cdum[1 - 1], &cdum[1 - 1], n, &cdum[1 - 1], -1, ierr);
            lwork_Cunmbr_prc_nn = castINTEGER(cdum[1 - 1].real());
            //
            Cunmbr("Q", "L", "N", m, m, n, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], m, &cdum[1 - 1], -1, ierr);
            lwork_Cunmbr_qln_mm = castINTEGER(cdum[1 - 1].real());
            //
            Cunmbr("Q", "L", "N", m, n, n, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], m, &cdum[1 - 1], -1, ierr);
            lwork_Cunmbr_qln_mn = castINTEGER(cdum[1 - 1].real());
            //
            Cunmbr("Q", "L", "N", n, n, n, &cdum[1 - 1], n, &cdum[1 - 1], &cdum[1 - 1], n, &cdum[1 - 1], -1, ierr);
            lwork_Cunmbr_qln_nn = castINTEGER(cdum[1 - 1].real());
            //
            if (m >= mnthr1) {
                if (wntqn) {
                    //
                    //                 Path 1 (M >> N, JOBZ='N')
                    //
                    maxwrk = n + lwork_Cgeqrf_mn;
                    maxwrk = max(maxwrk, 2 * n + lwork_Cgebrd_nn);
                    minwrk = 3 * n;
                } else if (wntqo) {
                    //
                    //                 Path 2 (M >> N, JOBZ='O')
                    //
                    wrkbl = n + lwork_Cgeqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_mn);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd_nn);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cunmbr_qln_nn);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cunmbr_prc_nn);
                    maxwrk = m * n + n * n + wrkbl;
                    minwrk = 2 * n * n + 3 * n;
                } else if (wntqs) {
                    //
                    //                 Path 3 (M >> N, JOBZ='S')
                    //
                    wrkbl = n + lwork_Cgeqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_mn);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd_nn);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cunmbr_qln_nn);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cunmbr_prc_nn);
                    maxwrk = n * n + wrkbl;
                    minwrk = n * n + 3 * n;
                } else if (wntqa) {
                    //
                    //                 Path 4 (M >> N, JOBZ='A')
                    //
                    wrkbl = n + lwork_Cgeqrf_mn;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_mm);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd_nn);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cunmbr_qln_nn);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cunmbr_prc_nn);
                    maxwrk = n * n + wrkbl;
                    minwrk = n * n + max(3 * n, n + m);
                }
            } else if (m >= mnthr2) {
                //
                //              Path 5 (M >> N, but not as much as MNTHR1)
                //
                maxwrk = 2 * n + lwork_Cgebrd_mn;
                minwrk = 2 * n + m;
                if (wntqo) {
                    //                 Path 5o (M >> N, JOBZ='O')
                    maxwrk = max(maxwrk, 2 * n + lwork_Cungbr_p_nn);
                    maxwrk = max(maxwrk, 2 * n + lwork_Cungbr_q_mn);
                    maxwrk += m * n;
                    minwrk += n * n;
                } else if (wntqs) {
                    //                 Path 5s (M >> N, JOBZ='S')
                    maxwrk = max(maxwrk, 2 * n + lwork_Cungbr_p_nn);
                    maxwrk = max(maxwrk, 2 * n + lwork_Cungbr_q_mn);
                } else if (wntqa) {
                    //                 Path 5a (M >> N, JOBZ='A')
                    maxwrk = max(maxwrk, 2 * n + lwork_Cungbr_p_nn);
                    maxwrk = max(maxwrk, 2 * n + lwork_Cungbr_q_mm);
                }
            } else {
                //
                //              Path 6 (M >= N, but not much larger)
                //
                maxwrk = 2 * n + lwork_Cgebrd_mn;
                minwrk = 2 * n + m;
                if (wntqo) {
                    //                 Path 6o (M >= N, JOBZ='O')
                    maxwrk = max(maxwrk, 2 * n + lwork_Cunmbr_prc_nn);
                    maxwrk = max(maxwrk, 2 * n + lwork_Cunmbr_qln_mn);
                    maxwrk += m * n;
                    minwrk += n * n;
                } else if (wntqs) {
                    //                 Path 6s (M >= N, JOBZ='S')
                    maxwrk = max(maxwrk, 2 * n + lwork_Cunmbr_qln_mn);
                    maxwrk = max(maxwrk, 2 * n + lwork_Cunmbr_prc_nn);
                } else if (wntqa) {
                    //                 Path 6a (M >= N, JOBZ='A')
                    maxwrk = max(maxwrk, 2 * n + lwork_Cunmbr_qln_mm);
                    maxwrk = max(maxwrk, 2 * n + lwork_Cunmbr_prc_nn);
                }
            }
        } else if (minmn > 0) {
            //
            //           There is no complex work space needed for bidiagonal SVD
            //           The real work space needed for bidiagonal SVD (Rbdsdc) is
            //           BDSPAC = 3*M*M + 4*M for singular values and vectors;
            //           BDSPAC = 4*M         for singular values only;
            //           not including e, RU, and RVT matrices.
            //
            //           Compute space preferred for each routine
            Cgebrd(m, n, &cdum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cgebrd_mn = castINTEGER(cdum[1 - 1].real());
            //
            Cgebrd(m, m, &cdum[1 - 1], m, &dum[1 - 1], &dum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cgebrd_mm = castINTEGER(cdum[1 - 1].real());
            //
            Cgelqf(m, n, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cgelqf_mn = castINTEGER(cdum[1 - 1].real());
            //
            Cungbr("P", m, n, m, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungbr_p_mn = castINTEGER(cdum[1 - 1].real());
            //
            Cungbr("P", n, n, m, &cdum[1 - 1], n, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungbr_p_nn = castINTEGER(cdum[1 - 1].real());
            //
            Cungbr("Q", m, m, n, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungbr_q_mm = castINTEGER(cdum[1 - 1].real());
            //
            Cunglq(m, n, m, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cunglq_mn = castINTEGER(cdum[1 - 1].real());
            //
            Cunglq(n, n, m, &cdum[1 - 1], n, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cunglq_nn = castINTEGER(cdum[1 - 1].real());
            //
            Cunmbr("P", "R", "C", m, m, m, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], m, &cdum[1 - 1], -1, ierr);
            lwork_Cunmbr_prc_mm = castINTEGER(cdum[1 - 1].real());
            //
            Cunmbr("P", "R", "C", m, n, m, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], m, &cdum[1 - 1], -1, ierr);
            lwork_Cunmbr_prc_mn = castINTEGER(cdum[1 - 1].real());
            //
            Cunmbr("P", "R", "C", n, n, m, &cdum[1 - 1], n, &cdum[1 - 1], &cdum[1 - 1], n, &cdum[1 - 1], -1, ierr);
            lwork_Cunmbr_prc_nn = castINTEGER(cdum[1 - 1].real());
            //
            Cunmbr("Q", "L", "N", m, m, m, &cdum[1 - 1], m, &cdum[1 - 1], &cdum[1 - 1], m, &cdum[1 - 1], -1, ierr);
            lwork_Cunmbr_qln_mm = castINTEGER(cdum[1 - 1].real());
            //
            if (n >= mnthr1) {
                if (wntqn) {
                    //
                    //                 Path 1t (N >> M, JOBZ='N')
                    //
                    maxwrk = m + lwork_Cgelqf_mn;
                    maxwrk = max(maxwrk, 2 * m + lwork_Cgebrd_mm);
                    minwrk = 3 * m;
                } else if (wntqo) {
                    //
                    //                 Path 2t (N >> M, JOBZ='O')
                    //
                    wrkbl = m + lwork_Cgelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_mn);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd_mm);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cunmbr_qln_mm);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cunmbr_prc_mm);
                    maxwrk = m * n + m * m + wrkbl;
                    minwrk = 2 * m * m + 3 * m;
                } else if (wntqs) {
                    //
                    //                 Path 3t (N >> M, JOBZ='S')
                    //
                    wrkbl = m + lwork_Cgelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_mn);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd_mm);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cunmbr_qln_mm);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cunmbr_prc_mm);
                    maxwrk = m * m + wrkbl;
                    minwrk = m * m + 3 * m;
                } else if (wntqa) {
                    //
                    //                 Path 4t (N >> M, JOBZ='A')
                    //
                    wrkbl = m + lwork_Cgelqf_mn;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_nn);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd_mm);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cunmbr_qln_mm);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cunmbr_prc_mm);
                    maxwrk = m * m + wrkbl;
                    minwrk = m * m + max(3 * m, m + n);
                }
            } else if (n >= mnthr2) {
                //
                //              Path 5t (N >> M, but not as much as MNTHR1)
                //
                maxwrk = 2 * m + lwork_Cgebrd_mn;
                minwrk = 2 * m + n;
                if (wntqo) {
                    //                 Path 5to (N >> M, JOBZ='O')
                    maxwrk = max(maxwrk, 2 * m + lwork_Cungbr_q_mm);
                    maxwrk = max(maxwrk, 2 * m + lwork_Cungbr_p_mn);
                    maxwrk += m * n;
                    minwrk += m * m;
                } else if (wntqs) {
                    //                 Path 5ts (N >> M, JOBZ='S')
                    maxwrk = max(maxwrk, 2 * m + lwork_Cungbr_q_mm);
                    maxwrk = max(maxwrk, 2 * m + lwork_Cungbr_p_mn);
                } else if (wntqa) {
                    //                 Path 5ta (N >> M, JOBZ='A')
                    maxwrk = max(maxwrk, 2 * m + lwork_Cungbr_q_mm);
                    maxwrk = max(maxwrk, 2 * m + lwork_Cungbr_p_nn);
                }
            } else {
                //
                //              Path 6t (N > M, but not much larger)
                //
                maxwrk = 2 * m + lwork_Cgebrd_mn;
                minwrk = 2 * m + n;
                if (wntqo) {
                    //                 Path 6to (N > M, JOBZ='O')
                    maxwrk = max(maxwrk, 2 * m + lwork_Cunmbr_qln_mm);
                    maxwrk = max(maxwrk, 2 * m + lwork_Cunmbr_prc_mn);
                    maxwrk += m * n;
                    minwrk += m * m;
                } else if (wntqs) {
                    //                 Path 6ts (N > M, JOBZ='S')
                    maxwrk = max(maxwrk, 2 * m + lwork_Cunmbr_qln_mm);
                    maxwrk = max(maxwrk, 2 * m + lwork_Cunmbr_prc_mn);
                } else if (wntqa) {
                    //                 Path 6ta (N > M, JOBZ='A')
                    maxwrk = max(maxwrk, 2 * m + lwork_Cunmbr_qln_mm);
                    maxwrk = max(maxwrk, 2 * m + lwork_Cunmbr_prc_nn);
                }
            }
        }
        maxwrk = max(maxwrk, minwrk);
    }
    if (info == 0) {
        work[1 - 1] = maxwrk;
        if (lwork < minwrk && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cgesdd", -info);
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
    REAL anrm = Clange("M", m, n, a, lda, dum);
    if (Risnan(anrm)) {
        info = -4;
        return;
    }
    INTEGER iscl = 0;
    const REAL zero = 0.0;
    if (anrm > zero && anrm < smlnum) {
        iscl = 1;
        Clascl("G", 0, 0, anrm, smlnum, m, n, a, lda, ierr);
    } else if (anrm > bignum) {
        iscl = 1;
        Clascl("G", 0, 0, anrm, bignum, m, n, a, lda, ierr);
    }
    //
    INTEGER itau = 0;
    INTEGER nwork = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER ie = 0;
    INTEGER itauq = 0;
    INTEGER itaup = 0;
    INTEGER nrwork = 0;
    INTEGER idum[1];
    INTEGER iu = 0;
    INTEGER ldwrku = 0;
    INTEGER ir = 0;
    INTEGER ldwrkr = 0;
    INTEGER iru = 0;
    INTEGER irvt = 0;
    INTEGER i = 0;
    INTEGER chunk = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER ivt = 0;
    INTEGER ldwkvt = 0;
    INTEGER il = 0;
    INTEGER ldwrkl = 0;
    INTEGER blk = 0;
    if (m >= n) {
        //
        //        A has at least as many rows as columns. If A has sufficiently
        //        more rows than columns, first reduce using the QR
        //        decomposition (if sufficient workspace available)
        //
        if (m >= mnthr1) {
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
                //              CWorkspace: need   N [tau] + N    [work]
                //              CWorkspace: prefer N [tau] + N*NB [work]
                //              RWorkspace: need   0
                //
                Cgeqrf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Zero out below R
                //
                Claset("L", n - 1, n - 1, czero, czero, &a[(2 - 1)], lda);
                ie = 1;
                itauq = 1;
                itaup = itauq + n;
                nwork = itaup + n;
                //
                //              Bidiagonalize R in A
                //              CWorkspace: need   2*N [tauq, taup] + N      [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + 2*N*NB [work]
                //              RWorkspace: need   N [e]
                //
                Cgebrd(n, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                nrwork = ie + n;
                //
                //              Perform bidiagonal SVD, compute singular values only
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + BDSPAC
                //
                Rbdsdc("U", "N", n, s, &rwork[ie - 1], dum, 1, dum, 1, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
            } else if (wntqo) {
                //
                //              Path 2 (M >> N, JOBZ='O')
                //              N left singular vectors to be overwritten on A and
                //              N right singular vectors to be computed in VT
                //
                iu = 1;
                //
                //              WORK(IU) is N by N
                //
                ldwrku = n;
                ir = iu + ldwrku * n;
                if (lwork >= m * n + n * n + 3 * n) {
                    //
                    //                 WORK(IR) is M by N
                    //
                    ldwrkr = m;
                } else {
                    ldwrkr = (lwork - n * n - 3 * n) / n;
                }
                itau = ir + ldwrkr * n;
                nwork = itau + n;
                //
                //              Compute A=Q*R
                //              CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work]
                //              CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work]
                //              RWorkspace: need   0
                //
                Cgeqrf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy R to WORK( IR ), zeroing out below it
                //
                Clacpy("U", n, n, a, lda, &work[ir - 1], ldwrkr);
                Claset("L", n - 1, n - 1, czero, czero, &work[(ir + 1) - 1], ldwrkr);
                //
                //              Generate Q in A
                //              CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work]
                //              CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work]
                //              RWorkspace: need   0
                //
                Cungqr(m, n, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                ie = 1;
                itauq = itau;
                itaup = itauq + n;
                nwork = itaup + n;
                //
                //              Bidiagonalize R in WORK(IR)
                //              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N      [work]
                //              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + 2*N*NB [work]
                //              RWorkspace: need   N [e]
                //
                Cgebrd(n, n, &work[ir - 1], ldwrkr, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of R in WORK(IRU) and computing right singular vectors
                //              of R in WORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
                //
                iru = ie + n;
                irvt = iru + n * n;
                nrwork = irvt + n * n;
                Rbdsdc("U", "I", n, s, &rwork[ie - 1], &rwork[iru - 1], n, &rwork[irvt - 1], n, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
                //              Overwrite WORK(IU) by the left singular vectors of R
                //              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", n, n, &rwork[iru - 1], n, &work[iu - 1], ldwrku);
                Cunmbr("Q", "L", "N", n, n, n, &work[ir - 1], ldwrkr, &work[itauq - 1], &work[iu - 1], ldwrku, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix VT
                //              Overwrite VT by the right singular vectors of R
                //              CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", n, n, &rwork[irvt - 1], n, vt, ldvt);
                Cunmbr("P", "R", "C", n, n, n, &work[ir - 1], ldwrkr, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply Q in A by left singular vectors of R in
                //              WORK(IU), storing result in WORK(IR) and copying to A
                //              CWorkspace: need   N*N [U] + N*N [R]
                //              CWorkspace: prefer N*N [U] + M*N [R]
                //              RWorkspace: need   0
                //
                for (i = 1; i <= m; i = i + ldwrkr) {
                    chunk = min(m - i + 1, ldwrkr);
                    Cgemm("N", "N", chunk, n, n, cone, &a[(i - 1)], lda, &work[iu - 1], ldwrku, czero, &work[ir - 1], ldwrkr);
                    Clacpy("F", chunk, n, &work[ir - 1], ldwrkr, &a[(i - 1)], lda);
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
                //              CWorkspace: need   N*N [R] + N [tau] + N    [work]
                //              CWorkspace: prefer N*N [R] + N [tau] + N*NB [work]
                //              RWorkspace: need   0
                //
                Cgeqrf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy R to WORK(IR), zeroing out below it
                //
                Clacpy("U", n, n, a, lda, &work[ir - 1], ldwrkr);
                Claset("L", n - 1, n - 1, czero, czero, &work[(ir + 1) - 1], ldwrkr);
                //
                //              Generate Q in A
                //              CWorkspace: need   N*N [R] + N [tau] + N    [work]
                //              CWorkspace: prefer N*N [R] + N [tau] + N*NB [work]
                //              RWorkspace: need   0
                //
                Cungqr(m, n, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                ie = 1;
                itauq = itau;
                itaup = itauq + n;
                nwork = itaup + n;
                //
                //              Bidiagonalize R in WORK(IR)
                //              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N      [work]
                //              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + 2*N*NB [work]
                //              RWorkspace: need   N [e]
                //
                Cgebrd(n, n, &work[ir - 1], ldwrkr, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
                //
                iru = ie + n;
                irvt = iru + n * n;
                nrwork = irvt + n * n;
                Rbdsdc("U", "I", n, s, &rwork[ie - 1], &rwork[iru - 1], n, &rwork[irvt - 1], n, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRU) to complex matrix U
                //              Overwrite U by left singular vectors of R
                //              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", n, n, &rwork[iru - 1], n, u, ldu);
                Cunmbr("Q", "L", "N", n, n, n, &work[ir - 1], ldwrkr, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix VT
                //              Overwrite VT by right singular vectors of R
                //              CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", n, n, &rwork[irvt - 1], n, vt, ldvt);
                Cunmbr("P", "R", "C", n, n, n, &work[ir - 1], ldwrkr, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply Q in A by left singular vectors of R in
                //              WORK(IR), storing result in U
                //              CWorkspace: need   N*N [R]
                //              RWorkspace: need   0
                //
                Clacpy("F", n, n, u, ldu, &work[ir - 1], ldwrkr);
                Cgemm("N", "N", m, n, n, cone, a, lda, &work[ir - 1], ldwrkr, czero, u, ldu);
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
                //              CWorkspace: need   N*N [U] + N [tau] + N    [work]
                //              CWorkspace: prefer N*N [U] + N [tau] + N*NB [work]
                //              RWorkspace: need   0
                //
                Cgeqrf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                Clacpy("L", m, n, a, lda, u, ldu);
                //
                //              Generate Q in U
                //              CWorkspace: need   N*N [U] + N [tau] + M    [work]
                //              CWorkspace: prefer N*N [U] + N [tau] + M*NB [work]
                //              RWorkspace: need   0
                //
                Cungqr(m, m, n, u, ldu, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Produce R in A, zeroing out below it
                //
                Claset("L", n - 1, n - 1, czero, czero, &a[(2 - 1)], lda);
                ie = 1;
                itauq = itau;
                itaup = itauq + n;
                nwork = itaup + n;
                //
                //              Bidiagonalize R in A
                //              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N      [work]
                //              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + 2*N*NB [work]
                //              RWorkspace: need   N [e]
                //
                Cgebrd(n, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                iru = ie + n;
                irvt = iru + n * n;
                nrwork = irvt + n * n;
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
                //
                Rbdsdc("U", "I", n, s, &rwork[ie - 1], &rwork[iru - 1], n, &rwork[irvt - 1], n, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
                //              Overwrite WORK(IU) by left singular vectors of R
                //              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", n, n, &rwork[iru - 1], n, &work[iu - 1], ldwrku);
                Cunmbr("Q", "L", "N", n, n, n, a, lda, &work[itauq - 1], &work[iu - 1], ldwrku, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix VT
                //              Overwrite VT by right singular vectors of R
                //              CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", n, n, &rwork[irvt - 1], n, vt, ldvt);
                Cunmbr("P", "R", "C", n, n, n, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply Q in U by left singular vectors of R in
                //              WORK(IU), storing result in A
                //              CWorkspace: need   N*N [U]
                //              RWorkspace: need   0
                //
                Cgemm("N", "N", m, n, n, cone, u, ldu, &work[iu - 1], ldwrku, czero, a, lda);
                //
                //              Copy left singular vectors of A from A to U
                //
                Clacpy("F", m, n, a, lda, u, ldu);
                //
            }
            //
        } else if (m >= mnthr2) {
            //
            //           MNTHR2 <= M < MNTHR1
            //
            //           Path 5 (M >> N, but not as much as MNTHR1)
            //           Reduce to bidiagonal form without QR decomposition, use
            //           Cungbr and matrix multiplication to compute singular vectors
            //
            ie = 1;
            nrwork = ie + n;
            itauq = 1;
            itaup = itauq + n;
            nwork = itaup + n;
            //
            //           Bidiagonalize A
            //           CWorkspace: need   2*N [tauq, taup] + M        [work]
            //           CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work]
            //           RWorkspace: need   N [e]
            //
            Cgebrd(m, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
            if (wntqn) {
                //
                //              Path 5n (M >> N, JOBZ='N')
                //              Compute singular values only
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + BDSPAC
                //
                Rbdsdc("U", "N", n, s, &rwork[ie - 1], dum, 1, dum, 1, dum, idum, &rwork[nrwork - 1], iwork, info);
            } else if (wntqo) {
                iu = nwork;
                iru = nrwork;
                irvt = iru + n * n;
                nrwork = irvt + n * n;
                //
                //              Path 5o (M >> N, JOBZ='O')
                //              Copy A to VT, generate P**H
                //              CWorkspace: need   2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacpy("U", n, n, a, lda, vt, ldvt);
                Cungbr("P", n, n, n, vt, ldvt, &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Generate Q in A
                //              CWorkspace: need   2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Cungbr("Q", m, n, n, a, lda, &work[itauq - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                if (lwork >= m * n + 3 * n) {
                    //
                    //                 WORK( IU ) is M by N
                    //
                    ldwrku = m;
                } else {
                    //
                    //                 WORK(IU) is LDWRKU by N
                    //
                    ldwrku = (lwork - 3 * n) / n;
                }
                nwork = iu + ldwrku * n;
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
                //
                Rbdsdc("U", "I", n, s, &rwork[ie - 1], &rwork[iru - 1], n, &rwork[irvt - 1], n, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Multiply real matrix RWORK(IRVT) by P**H in VT,
                //              storing the result in WORK(IU), copying to VT
                //              CWorkspace: need   2*N [tauq, taup] + N*N [U]
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]
                //
                Clarcm(n, n, &rwork[irvt - 1], n, vt, ldvt, &work[iu - 1], ldwrku, &rwork[nrwork - 1]);
                Clacpy("F", n, n, &work[iu - 1], ldwrku, vt, ldvt);
                //
                //              Multiply Q in A by real matrix RWORK(IRU), storing the
                //              result in WORK(IU), copying to A
                //              CWorkspace: need   2*N [tauq, taup] + N*N [U]
                //              CWorkspace: prefer 2*N [tauq, taup] + M*N [U]
                //              RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork]
                //              RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here
                //
                nrwork = irvt;
                for (i = 1; i <= m; i = i + ldwrku) {
                    chunk = min(m - i + 1, ldwrku);
                    Clacrm(chunk, n, &a[(i - 1)], lda, &rwork[iru - 1], n, &work[iu - 1], ldwrku, &rwork[nrwork - 1]);
                    Clacpy("F", chunk, n, &work[iu - 1], ldwrku, &a[(i - 1)], lda);
                }
                //
            } else if (wntqs) {
                //
                //              Path 5s (M >> N, JOBZ='S')
                //              Copy A to VT, generate P**H
                //              CWorkspace: need   2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacpy("U", n, n, a, lda, vt, ldvt);
                Cungbr("P", n, n, n, vt, ldvt, &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy A to U, generate Q
                //              CWorkspace: need   2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacpy("L", m, n, a, lda, u, ldu);
                Cungbr("Q", m, n, n, u, ldu, &work[itauq - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
                //
                iru = nrwork;
                irvt = iru + n * n;
                nrwork = irvt + n * n;
                Rbdsdc("U", "I", n, s, &rwork[ie - 1], &rwork[iru - 1], n, &rwork[irvt - 1], n, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Multiply real matrix RWORK(IRVT) by P**H in VT,
                //              storing the result in A, copying to VT
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]
                //
                Clarcm(n, n, &rwork[irvt - 1], n, vt, ldvt, a, lda, &rwork[nrwork - 1]);
                Clacpy("F", n, n, a, lda, vt, ldvt);
                //
                //              Multiply Q in U by real matrix RWORK(IRU), storing the
                //              result in A, copying to U
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here
                //
                nrwork = irvt;
                Clacrm(m, n, u, ldu, &rwork[iru - 1], n, a, lda, &rwork[nrwork - 1]);
                Clacpy("F", m, n, a, lda, u, ldu);
            } else {
                //
                //              Path 5a (M >> N, JOBZ='A')
                //              Copy A to VT, generate P**H
                //              CWorkspace: need   2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacpy("U", n, n, a, lda, vt, ldvt);
                Cungbr("P", n, n, n, vt, ldvt, &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy A to U, generate Q
                //              CWorkspace: need   2*N [tauq, taup] + M    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacpy("L", m, n, a, lda, u, ldu);
                Cungbr("Q", m, m, n, u, ldu, &work[itauq - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
                //
                iru = nrwork;
                irvt = iru + n * n;
                nrwork = irvt + n * n;
                Rbdsdc("U", "I", n, s, &rwork[ie - 1], &rwork[iru - 1], n, &rwork[irvt - 1], n, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Multiply real matrix RWORK(IRVT) by P**H in VT,
                //              storing the result in A, copying to VT
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]
                //
                Clarcm(n, n, &rwork[irvt - 1], n, vt, ldvt, a, lda, &rwork[nrwork - 1]);
                Clacpy("F", n, n, a, lda, vt, ldvt);
                //
                //              Multiply Q in U by real matrix RWORK(IRU), storing the
                //              result in A, copying to U
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here
                //
                nrwork = irvt;
                Clacrm(m, n, u, ldu, &rwork[iru - 1], n, a, lda, &rwork[nrwork - 1]);
                Clacpy("F", m, n, a, lda, u, ldu);
            }
            //
        } else {
            //
            //           M .LT. MNTHR2
            //
            //           Path 6 (M >= N, but not much larger)
            //           Reduce to bidiagonal form without QR decomposition
            //           Use Cunmbr to compute singular vectors
            //
            ie = 1;
            nrwork = ie + n;
            itauq = 1;
            itaup = itauq + n;
            nwork = itaup + n;
            //
            //           Bidiagonalize A
            //           CWorkspace: need   2*N [tauq, taup] + M        [work]
            //           CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work]
            //           RWorkspace: need   N [e]
            //
            Cgebrd(m, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
            if (wntqn) {
                //
                //              Path 6n (M >= N, JOBZ='N')
                //              Compute singular values only
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + BDSPAC
                //
                Rbdsdc("U", "N", n, s, &rwork[ie - 1], dum, 1, dum, 1, dum, idum, &rwork[nrwork - 1], iwork, info);
            } else if (wntqo) {
                iu = nwork;
                iru = nrwork;
                irvt = iru + n * n;
                nrwork = irvt + n * n;
                if (lwork >= m * n + 3 * n) {
                    //
                    //                 WORK( IU ) is M by N
                    //
                    ldwrku = m;
                } else {
                    //
                    //                 WORK( IU ) is LDWRKU by N
                    //
                    ldwrku = (lwork - 3 * n) / n;
                }
                nwork = iu + ldwrku * n;
                //
                //              Path 6o (M >= N, JOBZ='O')
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
                //
                Rbdsdc("U", "I", n, s, &rwork[ie - 1], &rwork[iru - 1], n, &rwork[irvt - 1], n, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix VT
                //              Overwrite VT by right singular vectors of A
                //              CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work]
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]
                //
                Clacp2("F", n, n, &rwork[irvt - 1], n, vt, ldvt);
                Cunmbr("P", "R", "C", n, n, n, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                if (lwork >= m * n + 3 * n) {
                    //
                    //                 Path 6o-fast
                    //                 Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
                    //                 Overwrite WORK(IU) by left singular vectors of A, copying
                    //                 to A
                    //                 CWorkspace: need   2*N [tauq, taup] + M*N [U] + N    [work]
                    //                 CWorkspace: prefer 2*N [tauq, taup] + M*N [U] + N*NB [work]
                    //                 RWorkspace: need   N [e] + N*N [RU]
                    //
                    Claset("F", m, n, czero, czero, &work[iu - 1], ldwrku);
                    Clacp2("F", n, n, &rwork[iru - 1], n, &work[iu - 1], ldwrku);
                    Cunmbr("Q", "L", "N", m, n, n, a, lda, &work[itauq - 1], &work[iu - 1], ldwrku, &work[nwork - 1], lwork - nwork + 1, ierr);
                    Clacpy("F", m, n, &work[iu - 1], ldwrku, a, lda);
                } else {
                    //
                    //                 Path 6o-slow
                    //                 Generate Q in A
                    //                 CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work]
                    //                 CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work]
                    //                 RWorkspace: need   0
                    //
                    Cungbr("Q", m, n, n, a, lda, &work[itauq - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                    //
                    //                 Multiply Q in A by real matrix RWORK(IRU), storing the
                    //                 result in WORK(IU), copying to A
                    //                 CWorkspace: need   2*N [tauq, taup] + N*N [U]
                    //                 CWorkspace: prefer 2*N [tauq, taup] + M*N [U]
                    //                 RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork]
                    //                 RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here
                    //
                    nrwork = irvt;
                    for (i = 1; i <= m; i = i + ldwrku) {
                        chunk = min(m - i + 1, ldwrku);
                        Clacrm(chunk, n, &a[(i - 1)], lda, &rwork[iru - 1], n, &work[iu - 1], ldwrku, &rwork[nrwork - 1]);
                        Clacpy("F", chunk, n, &work[iu - 1], ldwrku, &a[(i - 1)], lda);
                    }
                }
                //
            } else if (wntqs) {
                //
                //              Path 6s (M >= N, JOBZ='S')
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
                //
                iru = nrwork;
                irvt = iru + n * n;
                nrwork = irvt + n * n;
                Rbdsdc("U", "I", n, s, &rwork[ie - 1], &rwork[iru - 1], n, &rwork[irvt - 1], n, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRU) to complex matrix U
                //              Overwrite U by left singular vectors of A
                //              CWorkspace: need   2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]
                //
                Claset("F", m, n, czero, czero, u, ldu);
                Clacp2("F", n, n, &rwork[iru - 1], n, u, ldu);
                Cunmbr("Q", "L", "N", m, n, n, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix VT
                //              Overwrite VT by right singular vectors of A
                //              CWorkspace: need   2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]
                //
                Clacp2("F", n, n, &rwork[irvt - 1], n, vt, ldvt);
                Cunmbr("P", "R", "C", n, n, n, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
            } else {
                //
                //              Path 6a (M >= N, JOBZ='A')
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC
                //
                iru = nrwork;
                irvt = iru + n * n;
                nrwork = irvt + n * n;
                Rbdsdc("U", "I", n, s, &rwork[ie - 1], &rwork[iru - 1], n, &rwork[irvt - 1], n, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Set the right corner of U to identity matrix
                //
                Claset("F", m, m, czero, czero, u, ldu);
                if (m > n) {
                    Claset("F", m - n, m - n, czero, cone, &u[((n + 1) - 1) + ((n + 1) - 1) * ldu], ldu);
                }
                //
                //              Copy real matrix RWORK(IRU) to complex matrix U
                //              Overwrite U by left singular vectors of A
                //              CWorkspace: need   2*N [tauq, taup] + M    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + M*NB [work]
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]
                //
                Clacp2("F", n, n, &rwork[iru - 1], n, u, ldu);
                Cunmbr("Q", "L", "N", m, m, n, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix VT
                //              Overwrite VT by right singular vectors of A
                //              CWorkspace: need   2*N [tauq, taup] + N    [work]
                //              CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
                //              RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]
                //
                Clacp2("F", n, n, &rwork[irvt - 1], n, vt, ldvt);
                Cunmbr("P", "R", "C", n, n, n, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
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
        if (n >= mnthr1) {
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
                //              CWorkspace: need   M [tau] + M    [work]
                //              CWorkspace: prefer M [tau] + M*NB [work]
                //              RWorkspace: need   0
                //
                Cgelqf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Zero out above L
                //
                Claset("U", m - 1, m - 1, czero, czero, &a[(2 - 1) * lda], lda);
                ie = 1;
                itauq = 1;
                itaup = itauq + m;
                nwork = itaup + m;
                //
                //              Bidiagonalize L in A
                //              CWorkspace: need   2*M [tauq, taup] + M      [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + 2*M*NB [work]
                //              RWorkspace: need   M [e]
                //
                Cgebrd(m, m, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                nrwork = ie + m;
                //
                //              Perform bidiagonal SVD, compute singular values only
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + BDSPAC
                //
                Rbdsdc("U", "N", m, s, &rwork[ie - 1], dum, 1, dum, 1, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
            } else if (wntqo) {
                //
                //              Path 2t (N >> M, JOBZ='O')
                //              M right singular vectors to be overwritten on A and
                //              M left singular vectors to be computed in U
                //
                ivt = 1;
                ldwkvt = m;
                //
                //              WORK(IVT) is M by M
                //
                il = ivt + ldwkvt * m;
                if (lwork >= m * n + m * m + 3 * m) {
                    //
                    //                 WORK(IL) M by N
                    //
                    ldwrkl = m;
                    chunk = n;
                } else {
                    //
                    //                 WORK(IL) is M by CHUNK
                    //
                    ldwrkl = m;
                    chunk = (lwork - m * m - 3 * m) / m;
                }
                itau = il + ldwrkl * chunk;
                nwork = itau + m;
                //
                //              Compute A=L*Q
                //              CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
                //              CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
                //              RWorkspace: need   0
                //
                Cgelqf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy L to WORK(IL), zeroing about above it
                //
                Clacpy("L", m, m, a, lda, &work[il - 1], ldwrkl);
                Claset("U", m - 1, m - 1, czero, czero, &work[(il + ldwrkl) - 1], ldwrkl);
                //
                //              Generate Q in A
                //              CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
                //              CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
                //              RWorkspace: need   0
                //
                Cunglq(m, n, m, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                ie = 1;
                itauq = itau;
                itaup = itauq + m;
                nwork = itaup + m;
                //
                //              Bidiagonalize L in WORK(IL)
                //              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M      [work]
                //              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + 2*M*NB [work]
                //              RWorkspace: need   M [e]
                //
                Cgebrd(m, m, &work[il - 1], ldwrkl, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC
                //
                iru = ie + m;
                irvt = iru + m * m;
                nrwork = irvt + m * m;
                Rbdsdc("U", "I", m, s, &rwork[ie - 1], &rwork[iru - 1], m, &rwork[irvt - 1], m, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
                //              Overwrite WORK(IU) by the left singular vectors of L
                //              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", m, m, &rwork[iru - 1], m, u, ldu);
                Cunmbr("Q", "L", "N", m, m, m, &work[il - 1], ldwrkl, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
                //              Overwrite WORK(IVT) by the right singular vectors of L
                //              CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", m, m, &rwork[irvt - 1], m, &work[ivt - 1], ldwkvt);
                Cunmbr("P", "R", "C", m, m, m, &work[il - 1], ldwrkl, &work[itaup - 1], &work[ivt - 1], ldwkvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply right singular vectors of L in WORK(IL) by Q
                //              in A, storing result in WORK(IL) and copying to A
                //              CWorkspace: need   M*M [VT] + M*M [L]
                //              CWorkspace: prefer M*M [VT] + M*N [L]
                //              RWorkspace: need   0
                //
                for (i = 1; i <= n; i = i + chunk) {
                    blk = min(n - i + 1, chunk);
                    Cgemm("N", "N", m, blk, m, cone, &work[ivt - 1], m, &a[(i - 1) * lda], lda, czero, &work[il - 1], ldwrkl);
                    Clacpy("F", m, blk, &work[il - 1], ldwrkl, &a[(i - 1) * lda], lda);
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
                //              CWorkspace: need   M*M [L] + M [tau] + M    [work]
                //              CWorkspace: prefer M*M [L] + M [tau] + M*NB [work]
                //              RWorkspace: need   0
                //
                Cgelqf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy L to WORK(IL), zeroing out above it
                //
                Clacpy("L", m, m, a, lda, &work[il - 1], ldwrkl);
                Claset("U", m - 1, m - 1, czero, czero, &work[(il + ldwrkl) - 1], ldwrkl);
                //
                //              Generate Q in A
                //              CWorkspace: need   M*M [L] + M [tau] + M    [work]
                //              CWorkspace: prefer M*M [L] + M [tau] + M*NB [work]
                //              RWorkspace: need   0
                //
                Cunglq(m, n, m, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                ie = 1;
                itauq = itau;
                itaup = itauq + m;
                nwork = itaup + m;
                //
                //              Bidiagonalize L in WORK(IL)
                //              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M      [work]
                //              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + 2*M*NB [work]
                //              RWorkspace: need   M [e]
                //
                Cgebrd(m, m, &work[il - 1], ldwrkl, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC
                //
                iru = ie + m;
                irvt = iru + m * m;
                nrwork = irvt + m * m;
                Rbdsdc("U", "I", m, s, &rwork[ie - 1], &rwork[iru - 1], m, &rwork[irvt - 1], m, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRU) to complex matrix U
                //              Overwrite U by left singular vectors of L
                //              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", m, m, &rwork[iru - 1], m, u, ldu);
                Cunmbr("Q", "L", "N", m, m, m, &work[il - 1], ldwrkl, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix VT
                //              Overwrite VT by left singular vectors of L
                //              CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", m, m, &rwork[irvt - 1], m, vt, ldvt);
                Cunmbr("P", "R", "C", m, m, m, &work[il - 1], ldwrkl, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy VT to WORK(IL), multiply right singular vectors of L
                //              in WORK(IL) by Q in A, storing result in VT
                //              CWorkspace: need   M*M [L]
                //              RWorkspace: need   0
                //
                Clacpy("F", m, m, vt, ldvt, &work[il - 1], ldwrkl);
                Cgemm("N", "N", m, n, m, cone, &work[il - 1], ldwrkl, a, lda, czero, vt, ldvt);
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
                //              CWorkspace: need   M*M [VT] + M [tau] + M    [work]
                //              CWorkspace: prefer M*M [VT] + M [tau] + M*NB [work]
                //              RWorkspace: need   0
                //
                Cgelqf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                Clacpy("U", m, n, a, lda, vt, ldvt);
                //
                //              Generate Q in VT
                //              CWorkspace: need   M*M [VT] + M [tau] + N    [work]
                //              CWorkspace: prefer M*M [VT] + M [tau] + N*NB [work]
                //              RWorkspace: need   0
                //
                Cunglq(n, n, m, vt, ldvt, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Produce L in A, zeroing out above it
                //
                Claset("U", m - 1, m - 1, czero, czero, &a[(2 - 1) * lda], lda);
                ie = 1;
                itauq = itau;
                itaup = itauq + m;
                nwork = itaup + m;
                //
                //              Bidiagonalize L in A
                //              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M      [work]
                //              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + 2*M*NB [work]
                //              RWorkspace: need   M [e]
                //
                Cgebrd(m, m, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC
                //
                iru = ie + m;
                irvt = iru + m * m;
                nrwork = irvt + m * m;
                Rbdsdc("U", "I", m, s, &rwork[ie - 1], &rwork[iru - 1], m, &rwork[irvt - 1], m, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRU) to complex matrix U
                //              Overwrite U by left singular vectors of L
                //              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", m, m, &rwork[iru - 1], m, u, ldu);
                Cunmbr("Q", "L", "N", m, m, m, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
                //              Overwrite WORK(IVT) by right singular vectors of L
                //              CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacp2("F", m, m, &rwork[irvt - 1], m, &work[ivt - 1], ldwkvt);
                Cunmbr("P", "R", "C", m, m, m, a, lda, &work[itaup - 1], &work[ivt - 1], ldwkvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Multiply right singular vectors of L in WORK(IVT) by
                //              Q in VT, storing result in A
                //              CWorkspace: need   M*M [VT]
                //              RWorkspace: need   0
                //
                Cgemm("N", "N", m, n, m, cone, &work[ivt - 1], ldwkvt, vt, ldvt, czero, a, lda);
                //
                //              Copy right singular vectors of A from A to VT
                //
                Clacpy("F", m, n, a, lda, vt, ldvt);
                //
            }
            //
        } else if (n >= mnthr2) {
            //
            //           MNTHR2 <= N < MNTHR1
            //
            //           Path 5t (N >> M, but not as much as MNTHR1)
            //           Reduce to bidiagonal form without QR decomposition, use
            //           Cungbr and matrix multiplication to compute singular vectors
            //
            ie = 1;
            nrwork = ie + m;
            itauq = 1;
            itaup = itauq + m;
            nwork = itaup + m;
            //
            //           Bidiagonalize A
            //           CWorkspace: need   2*M [tauq, taup] + N        [work]
            //           CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work]
            //           RWorkspace: need   M [e]
            //
            Cgebrd(m, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
            //
            if (wntqn) {
                //
                //              Path 5tn (N >> M, JOBZ='N')
                //              Compute singular values only
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + BDSPAC
                //
                Rbdsdc("L", "N", m, s, &rwork[ie - 1], dum, 1, dum, 1, dum, idum, &rwork[nrwork - 1], iwork, info);
            } else if (wntqo) {
                irvt = nrwork;
                iru = irvt + m * m;
                nrwork = iru + m * m;
                ivt = nwork;
                //
                //              Path 5to (N >> M, JOBZ='O')
                //              Copy A to U, generate Q
                //              CWorkspace: need   2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacpy("L", m, m, a, lda, u, ldu);
                Cungbr("Q", m, m, n, u, ldu, &work[itauq - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Generate P**H in A
                //              CWorkspace: need   2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Cungbr("P", m, n, m, a, lda, &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                ldwkvt = m;
                if (lwork >= m * n + 3 * m) {
                    //
                    //                 WORK( IVT ) is M by N
                    //
                    nwork = ivt + ldwkvt * n;
                    chunk = n;
                } else {
                    //
                    //                 WORK( IVT ) is M by CHUNK
                    //
                    chunk = (lwork - 3 * m) / m;
                    nwork = ivt + ldwkvt * chunk;
                }
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
                //
                Rbdsdc("L", "I", m, s, &rwork[ie - 1], &rwork[iru - 1], m, &rwork[irvt - 1], m, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Multiply Q in U by real matrix RWORK(IRVT)
                //              storing the result in WORK(IVT), copying to U
                //              CWorkspace: need   2*M [tauq, taup] + M*M [VT]
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]
                //
                Clacrm(m, m, u, ldu, &rwork[iru - 1], m, &work[ivt - 1], ldwkvt, &rwork[nrwork - 1]);
                Clacpy("F", m, m, &work[ivt - 1], ldwkvt, u, ldu);
                //
                //              Multiply RWORK(IRVT) by P**H in A, storing the
                //              result in WORK(IVT), copying to A
                //              CWorkspace: need   2*M [tauq, taup] + M*M [VT]
                //              CWorkspace: prefer 2*M [tauq, taup] + M*N [VT]
                //              RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork]
                //              RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here
                //
                nrwork = iru;
                for (i = 1; i <= n; i = i + chunk) {
                    blk = min(n - i + 1, chunk);
                    Clarcm(m, blk, &rwork[irvt - 1], m, &a[(i - 1) * lda], lda, &work[ivt - 1], ldwkvt, &rwork[nrwork - 1]);
                    Clacpy("F", m, blk, &work[ivt - 1], ldwkvt, &a[(i - 1) * lda], lda);
                }
            } else if (wntqs) {
                //
                //              Path 5ts (N >> M, JOBZ='S')
                //              Copy A to U, generate Q
                //              CWorkspace: need   2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacpy("L", m, m, a, lda, u, ldu);
                Cungbr("Q", m, m, n, u, ldu, &work[itauq - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy A to VT, generate P**H
                //              CWorkspace: need   2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacpy("U", m, n, a, lda, vt, ldvt);
                Cungbr("P", m, n, m, vt, ldvt, &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
                //
                irvt = nrwork;
                iru = irvt + m * m;
                nrwork = iru + m * m;
                Rbdsdc("L", "I", m, s, &rwork[ie - 1], &rwork[iru - 1], m, &rwork[irvt - 1], m, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Multiply Q in U by real matrix RWORK(IRU), storing the
                //              result in A, copying to U
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]
                //
                Clacrm(m, m, u, ldu, &rwork[iru - 1], m, a, lda, &rwork[nrwork - 1]);
                Clacpy("F", m, m, a, lda, u, ldu);
                //
                //              Multiply real matrix RWORK(IRVT) by P**H in VT,
                //              storing the result in A, copying to VT
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here
                //
                nrwork = iru;
                Clarcm(m, n, &rwork[irvt - 1], m, vt, ldvt, a, lda, &rwork[nrwork - 1]);
                Clacpy("F", m, n, a, lda, vt, ldvt);
            } else {
                //
                //              Path 5ta (N >> M, JOBZ='A')
                //              Copy A to U, generate Q
                //              CWorkspace: need   2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   0
                //
                Clacpy("L", m, m, a, lda, u, ldu);
                Cungbr("Q", m, m, n, u, ldu, &work[itauq - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy A to VT, generate P**H
                //              CWorkspace: need   2*M [tauq, taup] + N    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + N*NB [work]
                //              RWorkspace: need   0
                //
                Clacpy("U", m, n, a, lda, vt, ldvt);
                Cungbr("P", n, n, m, vt, ldvt, &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
                //
                irvt = nrwork;
                iru = irvt + m * m;
                nrwork = iru + m * m;
                Rbdsdc("L", "I", m, s, &rwork[ie - 1], &rwork[iru - 1], m, &rwork[irvt - 1], m, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Multiply Q in U by real matrix RWORK(IRU), storing the
                //              result in A, copying to U
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]
                //
                Clacrm(m, m, u, ldu, &rwork[iru - 1], m, a, lda, &rwork[nrwork - 1]);
                Clacpy("F", m, m, a, lda, u, ldu);
                //
                //              Multiply real matrix RWORK(IRVT) by P**H in VT,
                //              storing the result in A, copying to VT
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here
                //
                nrwork = iru;
                Clarcm(m, n, &rwork[irvt - 1], m, vt, ldvt, a, lda, &rwork[nrwork - 1]);
                Clacpy("F", m, n, a, lda, vt, ldvt);
            }
            //
        } else {
            //
            //           N .LT. MNTHR2
            //
            //           Path 6t (N > M, but not much larger)
            //           Reduce to bidiagonal form without LQ decomposition
            //           Use Cunmbr to compute singular vectors
            //
            ie = 1;
            nrwork = ie + m;
            itauq = 1;
            itaup = itauq + m;
            nwork = itaup + m;
            //
            //           Bidiagonalize A
            //           CWorkspace: need   2*M [tauq, taup] + N        [work]
            //           CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work]
            //           RWorkspace: need   M [e]
            //
            Cgebrd(m, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
            if (wntqn) {
                //
                //              Path 6tn (N > M, JOBZ='N')
                //              Compute singular values only
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + BDSPAC
                //
                Rbdsdc("L", "N", m, s, &rwork[ie - 1], dum, 1, dum, 1, dum, idum, &rwork[nrwork - 1], iwork, info);
            } else if (wntqo) {
                //              Path 6to (N > M, JOBZ='O')
                ldwkvt = m;
                ivt = nwork;
                if (lwork >= m * n + 3 * m) {
                    //
                    //                 WORK( IVT ) is M by N
                    //
                    Claset("F", m, n, czero, czero, &work[ivt - 1], ldwkvt);
                    nwork = ivt + ldwkvt * n;
                } else {
                    //
                    //                 WORK( IVT ) is M by CHUNK
                    //
                    chunk = (lwork - 3 * m) / m;
                    nwork = ivt + ldwkvt * chunk;
                }
                //
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
                //
                irvt = nrwork;
                iru = irvt + m * m;
                nrwork = iru + m * m;
                Rbdsdc("L", "I", m, s, &rwork[ie - 1], &rwork[iru - 1], m, &rwork[irvt - 1], m, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRU) to complex matrix U
                //              Overwrite U by left singular vectors of A
                //              CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work]
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]
                //
                Clacp2("F", m, m, &rwork[iru - 1], m, u, ldu);
                Cunmbr("Q", "L", "N", m, m, n, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                if (lwork >= m * n + 3 * m) {
                    //
                    //                 Path 6to-fast
                    //                 Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
                    //                 Overwrite WORK(IVT) by right singular vectors of A,
                    //                 copying to A
                    //                 CWorkspace: need   2*M [tauq, taup] + M*N [VT] + M    [work]
                    //                 CWorkspace: prefer 2*M [tauq, taup] + M*N [VT] + M*NB [work]
                    //                 RWorkspace: need   M [e] + M*M [RVT]
                    //
                    Clacp2("F", m, m, &rwork[irvt - 1], m, &work[ivt - 1], ldwkvt);
                    Cunmbr("P", "R", "C", m, n, m, a, lda, &work[itaup - 1], &work[ivt - 1], ldwkvt, &work[nwork - 1], lwork - nwork + 1, ierr);
                    Clacpy("F", m, n, &work[ivt - 1], ldwkvt, a, lda);
                } else {
                    //
                    //                 Path 6to-slow
                    //                 Generate P**H in A
                    //                 CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work]
                    //                 CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work]
                    //                 RWorkspace: need   0
                    //
                    Cungbr("P", m, n, m, a, lda, &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, ierr);
                    //
                    //                 Multiply Q in A by real matrix RWORK(IRU), storing the
                    //                 result in WORK(IU), copying to A
                    //                 CWorkspace: need   2*M [tauq, taup] + M*M [VT]
                    //                 CWorkspace: prefer 2*M [tauq, taup] + M*N [VT]
                    //                 RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork]
                    //                 RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here
                    //
                    nrwork = iru;
                    for (i = 1; i <= n; i = i + chunk) {
                        blk = min(n - i + 1, chunk);
                        Clarcm(m, blk, &rwork[irvt - 1], m, &a[(i - 1) * lda], lda, &work[ivt - 1], ldwkvt, &rwork[nrwork - 1]);
                        Clacpy("F", m, blk, &work[ivt - 1], ldwkvt, &a[(i - 1) * lda], lda);
                    }
                }
            } else if (wntqs) {
                //
                //              Path 6ts (N > M, JOBZ='S')
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
                //
                irvt = nrwork;
                iru = irvt + m * m;
                nrwork = iru + m * m;
                Rbdsdc("L", "I", m, s, &rwork[ie - 1], &rwork[iru - 1], m, &rwork[irvt - 1], m, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRU) to complex matrix U
                //              Overwrite U by left singular vectors of A
                //              CWorkspace: need   2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]
                //
                Clacp2("F", m, m, &rwork[iru - 1], m, u, ldu);
                Cunmbr("Q", "L", "N", m, m, n, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix VT
                //              Overwrite VT by right singular vectors of A
                //              CWorkspace: need   2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   M [e] + M*M [RVT]
                //
                Claset("F", m, n, czero, czero, vt, ldvt);
                Clacp2("F", m, m, &rwork[irvt - 1], m, vt, ldvt);
                Cunmbr("P", "R", "C", m, n, m, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
            } else {
                //
                //              Path 6ta (N > M, JOBZ='A')
                //              Perform bidiagonal SVD, computing left singular vectors
                //              of bidiagonal matrix in RWORK(IRU) and computing right
                //              singular vectors of bidiagonal matrix in RWORK(IRVT)
                //              CWorkspace: need   0
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC
                //
                irvt = nrwork;
                iru = irvt + m * m;
                nrwork = iru + m * m;
                //
                Rbdsdc("L", "I", m, s, &rwork[ie - 1], &rwork[iru - 1], m, &rwork[irvt - 1], m, dum, idum, &rwork[nrwork - 1], iwork, info);
                //
                //              Copy real matrix RWORK(IRU) to complex matrix U
                //              Overwrite U by left singular vectors of A
                //              CWorkspace: need   2*M [tauq, taup] + M    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
                //              RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]
                //
                Clacp2("F", m, m, &rwork[iru - 1], m, u, ldu);
                Cunmbr("Q", "L", "N", m, m, n, a, lda, &work[itauq - 1], u, ldu, &work[nwork - 1], lwork - nwork + 1, ierr);
                //
                //              Set all of VT to identity matrix
                //
                Claset("F", n, n, czero, cone, vt, ldvt);
                //
                //              Copy real matrix RWORK(IRVT) to complex matrix VT
                //              Overwrite VT by right singular vectors of A
                //              CWorkspace: need   2*M [tauq, taup] + N    [work]
                //              CWorkspace: prefer 2*M [tauq, taup] + N*NB [work]
                //              RWorkspace: need   M [e] + M*M [RVT]
                //
                Clacp2("F", m, m, &rwork[irvt - 1], m, vt, ldvt);
                Cunmbr("P", "R", "C", n, n, m, a, lda, &work[itaup - 1], vt, ldvt, &work[nwork - 1], lwork - nwork + 1, ierr);
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
        if (info != 0 && anrm > bignum) {
            Rlascl("G", 0, 0, bignum, anrm, minmn - 1, 1, &rwork[ie - 1], minmn, ierr);
        }
        if (anrm < smlnum) {
            Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, s, minmn, ierr);
        }
        if (info != 0 && anrm < smlnum) {
            Rlascl("G", 0, 0, smlnum, anrm, minmn - 1, 1, &rwork[ie - 1], minmn, ierr);
        }
    }
    //
    //     Return optimal workspace in WORK(1)
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Cgesdd
    //
}
