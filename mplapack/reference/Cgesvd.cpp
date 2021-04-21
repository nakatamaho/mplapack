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

void Cgesvd(const char *jobu, const char *jobvt, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL *s, COMPLEX *u, INTEGER const ldu, COMPLEX *vt, INTEGER const ldvt, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER &info) {
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
    bool wntua = Mlsame(jobu, "A");
    bool wntus = Mlsame(jobu, "S");
    bool wntuas = wntua || wntus;
    bool wntuo = Mlsame(jobu, "O");
    bool wntun = Mlsame(jobu, "N");
    bool wntva = Mlsame(jobvt, "A");
    bool wntvs = Mlsame(jobvt, "S");
    bool wntvas = wntva || wntvs;
    bool wntvo = Mlsame(jobvt, "O");
    bool wntvn = Mlsame(jobvt, "N");
    bool lquery = (lwork == -1);
    //
    if (!(wntua || wntus || wntuo || wntun)) {
        info = -1;
    } else if (!(wntva || wntvs || wntvo || wntvn) || (wntvo && wntuo)) {
        info = -2;
    } else if (m < 0) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, m)) {
        info = -6;
    } else if (ldu < 1 || (wntuas && ldu < m)) {
        info = -9;
    } else if (ldvt < 1 || (wntva && ldvt < n) || (wntvs && ldvt < minmn)) {
        info = -11;
    }
    //
    //     Compute workspace
    //      (Note: Comments in the code beginning "Workspace:" describe the
    //       minimal amount of workspace needed at that point in the code,
    //       as well as the preferred amount for good performance.
    //       CWorkspace refers to complex workspace, and RWorkspace to
    //       real workspace. NB refers to the optimal block size for the
    //       immediately following subroutine, as returned by iMlaenv.)
    //
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER mnthr = 0;
    COMPLEX cdum[1];
    INTEGER ierr = 0;
    INTEGER lwork_Cgeqrf = 0;
    INTEGER lwork_Cungqr_n = 0;
    INTEGER lwork_Cungqr_m = 0;
    REAL dum[1];
    INTEGER lwork_Cgebrd = 0;
    INTEGER lwork_Cungbr_p = 0;
    INTEGER lwork_Cungbr_q = 0;
    INTEGER wrkbl = 0;
    INTEGER lwork_Cgelqf = 0;
    INTEGER lwork_Cunglq_n = 0;
    INTEGER lwork_Cunglq_m = 0;
    char jobu_jobvt[3];
    jobu_jobvt[0] = jobu[0];
    jobu_jobvt[1] = jobvt[0];
    jobu_jobvt[2] = '\0';
    if (info == 0) {
        minwrk = 1;
        maxwrk = 1;
        if (m >= n && minmn > 0) {
            //
            //           Space needed for Cbdsqr is BDSPAC = 5*N
            //
            mnthr = iMlaenv(6, "Cgesvd", jobu_jobvt, m, n, 0, 0);
            //           Compute space needed for Cgeqrf
            Cgeqrf(m, n, a, lda, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cgeqrf = castINTEGER(cdum[1 - 1].real());
            //           Compute space needed for Cungqr
            Cungqr(m, n, n, a, lda, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungqr_n = castINTEGER(cdum[1 - 1].real());
            Cungqr(m, m, n, a, lda, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungqr_m = castINTEGER(cdum[1 - 1].real());
            //           Compute space needed for Cgebrd
            Cgebrd(n, n, a, lda, s, &dum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cgebrd = castINTEGER(cdum[1 - 1].real());
            //           Compute space needed for Cungbr
            Cungbr("P", n, n, n, a, lda, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungbr_p = castINTEGER(cdum[1 - 1].real());
            Cungbr("Q", n, n, n, a, lda, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungbr_q = castINTEGER(cdum[1 - 1].real());
            //
            if (m >= mnthr) {
                if (wntun) {
                    //
                    //                 Path 1 (M much larger than N, JOBU='N')
                    //
                    maxwrk = n + lwork_Cgeqrf;
                    maxwrk = max(maxwrk, 2 * n + lwork_Cgebrd);
                    if (wntvo || wntvas) {
                        maxwrk = max(maxwrk, 2 * n + lwork_Cungbr_p);
                    }
                    minwrk = 3 * n;
                } else if (wntuo && wntvn) {
                    //
                    //                 Path 2 (M much larger than N, JOBU='O', JOBVT='N')
                    //
                    wrkbl = n + lwork_Cgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_n);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_q);
                    maxwrk = max(n * n + wrkbl, n * n + m * n);
                    minwrk = 2 * n + m;
                } else if (wntuo && wntvas) {
                    //
                    //                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
                    //                 'A')
                    //
                    wrkbl = n + lwork_Cgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_n);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_q);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_p);
                    maxwrk = max(n * n + wrkbl, n * n + m * n);
                    minwrk = 2 * n + m;
                } else if (wntus && wntvn) {
                    //
                    //                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
                    //
                    wrkbl = n + lwork_Cgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_n);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_q);
                    maxwrk = n * n + wrkbl;
                    minwrk = 2 * n + m;
                } else if (wntus && wntvo) {
                    //
                    //                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
                    //
                    wrkbl = n + lwork_Cgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_n);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_q);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_p);
                    maxwrk = 2 * n * n + wrkbl;
                    minwrk = 2 * n + m;
                } else if (wntus && wntvas) {
                    //
                    //                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
                    //                 'A')
                    //
                    wrkbl = n + lwork_Cgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_n);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_q);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_p);
                    maxwrk = n * n + wrkbl;
                    minwrk = 2 * n + m;
                } else if (wntua && wntvn) {
                    //
                    //                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
                    //
                    wrkbl = n + lwork_Cgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_m);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_q);
                    maxwrk = n * n + wrkbl;
                    minwrk = 2 * n + m;
                } else if (wntua && wntvo) {
                    //
                    //                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
                    //
                    wrkbl = n + lwork_Cgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_m);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_q);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_p);
                    maxwrk = 2 * n * n + wrkbl;
                    minwrk = 2 * n + m;
                } else if (wntua && wntvas) {
                    //
                    //                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
                    //                 'A')
                    //
                    wrkbl = n + lwork_Cgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Cungqr_m);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_q);
                    wrkbl = max(wrkbl, 2 * n + lwork_Cungbr_p);
                    maxwrk = n * n + wrkbl;
                    minwrk = 2 * n + m;
                }
            } else {
                //
                //              Path 10 (M at least N, but not much larger)
                //
                Cgebrd(m, n, a, lda, s, &dum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
                lwork_Cgebrd = castINTEGER(cdum[1 - 1].real());
                maxwrk = 2 * n + lwork_Cgebrd;
                if (wntus || wntuo) {
                    Cungbr("Q", m, n, n, a, lda, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
                    lwork_Cungbr_q = castINTEGER(cdum[1 - 1].real());
                    maxwrk = max(maxwrk, 2 * n + lwork_Cungbr_q);
                }
                if (wntua) {
                    Cungbr("Q", m, m, n, a, lda, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
                    lwork_Cungbr_q = castINTEGER(cdum[1 - 1].real());
                    maxwrk = max(maxwrk, 2 * n + lwork_Cungbr_q);
                }
                if (!wntvn) {
                    maxwrk = max(maxwrk, 2 * n + lwork_Cungbr_p);
                }
                minwrk = 2 * n + m;
            }
        } else if (minmn > 0) {
            //
            //           Space needed for Cbdsqr is BDSPAC = 5*M
            //
            mnthr = iMlaenv(6, "Cgesvd", jobu_jobvt, m, n, 0, 0);
            //           Compute space needed for Cgelqf
            Cgelqf(m, n, a, lda, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cgelqf = castINTEGER(cdum[1 - 1].real());
            //           Compute space needed for Cunglq
            Cunglq(n, n, m, &cdum[1 - 1], n, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cunglq_n = castINTEGER(cdum[1 - 1].real());
            Cunglq(m, n, m, a, lda, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cunglq_m = castINTEGER(cdum[1 - 1].real());
            //           Compute space needed for Cgebrd
            Cgebrd(m, m, a, lda, s, &dum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cgebrd = castINTEGER(cdum[1 - 1].real());
            //            Compute space needed for Cungbr P
            Cungbr("P", m, m, m, a, n, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungbr_p = castINTEGER(cdum[1 - 1].real());
            //           Compute space needed for Cungbr Q
            Cungbr("Q", m, m, m, a, n, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
            lwork_Cungbr_q = castINTEGER(cdum[1 - 1].real());
            if (n >= mnthr) {
                if (wntvn) {
                    //
                    //                 Path 1 (N much larger than M, JOBVT='N')
                    //
                    maxwrk = m + lwork_Cgelqf;
                    maxwrk = max(maxwrk, 2 * m + lwork_Cgebrd);
                    if (wntuo || wntuas) {
                        maxwrk = max(maxwrk, 2 * m + lwork_Cungbr_q);
                    }
                    minwrk = 3 * m;
                } else if (wntvo && wntun) {
                    //
                    //                 Path 2 (N much larger than M, JOBU='N', JOBVT='O')
                    //
                    wrkbl = m + lwork_Cgelqf;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_m);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_p);
                    maxwrk = max(m * m + wrkbl, m * m + m * n);
                    minwrk = 2 * m + n;
                } else if (wntvo && wntuas) {
                    //
                    //                 Path 3 (N much larger than M, JOBU='S' or 'A',
                    //                 JOBVT='O')
                    //
                    wrkbl = m + lwork_Cgelqf;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_m);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_p);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_q);
                    maxwrk = max(m * m + wrkbl, m * m + m * n);
                    minwrk = 2 * m + n;
                } else if (wntvs && wntun) {
                    //
                    //                 Path 4 (N much larger than M, JOBU='N', JOBVT='S')
                    //
                    wrkbl = m + lwork_Cgelqf;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_m);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_p);
                    maxwrk = m * m + wrkbl;
                    minwrk = 2 * m + n;
                } else if (wntvs && wntuo) {
                    //
                    //                 Path 5 (N much larger than M, JOBU='O', JOBVT='S')
                    //
                    wrkbl = m + lwork_Cgelqf;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_m);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_p);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_q);
                    maxwrk = 2 * m * m + wrkbl;
                    minwrk = 2 * m + n;
                } else if (wntvs && wntuas) {
                    //
                    //                 Path 6 (N much larger than M, JOBU='S' or 'A',
                    //                 JOBVT='S')
                    //
                    wrkbl = m + lwork_Cgelqf;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_m);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_p);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_q);
                    maxwrk = m * m + wrkbl;
                    minwrk = 2 * m + n;
                } else if (wntva && wntun) {
                    //
                    //                 Path 7 (N much larger than M, JOBU='N', JOBVT='A')
                    //
                    wrkbl = m + lwork_Cgelqf;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_n);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_p);
                    maxwrk = m * m + wrkbl;
                    minwrk = 2 * m + n;
                } else if (wntva && wntuo) {
                    //
                    //                 Path 8 (N much larger than M, JOBU='O', JOBVT='A')
                    //
                    wrkbl = m + lwork_Cgelqf;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_n);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_p);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_q);
                    maxwrk = 2 * m * m + wrkbl;
                    minwrk = 2 * m + n;
                } else if (wntva && wntuas) {
                    //
                    //                 Path 9 (N much larger than M, JOBU='S' or 'A',
                    //                 JOBVT='A')
                    //
                    wrkbl = m + lwork_Cgelqf;
                    wrkbl = max(wrkbl, m + lwork_Cunglq_n);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cgebrd);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_p);
                    wrkbl = max(wrkbl, 2 * m + lwork_Cungbr_q);
                    maxwrk = m * m + wrkbl;
                    minwrk = 2 * m + n;
                }
            } else {
                //
                //              Path 10 (N greater than M, but not much larger)
                //
                Cgebrd(m, n, a, lda, s, &dum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
                lwork_Cgebrd = castINTEGER(cdum[1 - 1].real());
                maxwrk = 2 * m + lwork_Cgebrd;
                if (wntvs || wntvo) {
                    //                Compute space needed for Cungbr P
                    Cungbr("P", m, n, m, a, n, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
                    lwork_Cungbr_p = castINTEGER(cdum[1 - 1].real());
                    maxwrk = max(maxwrk, 2 * m + lwork_Cungbr_p);
                }
                if (wntva) {
                    Cungbr("P", n, n, m, a, n, &cdum[1 - 1], &cdum[1 - 1], -1, ierr);
                    lwork_Cungbr_p = castINTEGER(cdum[1 - 1].real());
                    maxwrk = max(maxwrk, 2 * m + lwork_Cungbr_p);
                }
                if (!wntun) {
                    maxwrk = max(maxwrk, 2 * m + lwork_Cungbr_q);
                }
                minwrk = 2 * m + n;
            }
        }
        maxwrk = max(maxwrk, minwrk);
        work[1 - 1] = maxwrk;
        //
        if (lwork < minwrk && !lquery) {
            info = -13;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cgesvd", -info);
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
    INTEGER iwork = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER ie = 0;
    INTEGER itauq = 0;
    INTEGER itaup = 0;
    INTEGER ncvt = 0;
    INTEGER irwork = 0;
    INTEGER ir = 0;
    INTEGER ldwrku = 0;
    INTEGER ldwrkr = 0;
    INTEGER iu = 0;
    INTEGER i = 0;
    INTEGER chunk = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER ncu = 0;
    INTEGER nru = 0;
    INTEGER blk = 0;
    INTEGER nrvt = 0;
    if (m >= n) {
        //
        //        A has at least as many rows as columns. If A has sufficiently
        //        more rows than columns, first reduce using the QR
        //        decomposition (if sufficient workspace available)
        //
        if (m >= mnthr) {
            //
            if (wntun) {
                //
                //              Path 1 (M much larger than N, JOBU='N')
                //              No left singular vectors to be computed
                //
                itau = 1;
                iwork = itau + n;
                //
                //              Compute A=Q*R
                //              (CWorkspace: need 2*N, prefer N+N*NB)
                //              (RWorkspace: need 0)
                //
                Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                //
                //              Zero out below R
                //
                if (n > 1) {
                    Claset("L", n - 1, n - 1, czero, czero, &a[(2 - 1)], lda);
                }
                ie = 1;
                itauq = 1;
                itaup = itauq + n;
                iwork = itaup + n;
                //
                //              Bidiagonalize R in A
                //              (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                //              (RWorkspace: need N)
                //
                Cgebrd(n, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                ncvt = 0;
                if (wntvo || wntvas) {
                    //
                    //                 If right singular vectors desired, generate P'.
                    //                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("P", n, n, n, a, lda, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    ncvt = n;
                }
                irwork = ie + n;
                //
                //              Perform bidiagonal QR iteration, computing right
                //              singular vectors of A in A if desired
                //              (CWorkspace: 0)
                //              (RWorkspace: need BDSPAC)
                //
                Cbdsqr("U", n, ncvt, 0, 0, s, &rwork[ie - 1], a, lda, cdum, 1, cdum, 1, &rwork[irwork - 1], info);
                //
                //              If right singular vectors desired in VT, copy them there
                //
                if (wntvas) {
                    Clacpy("F", n, n, a, lda, vt, ldvt);
                }
                //
            } else if (wntuo && wntvn) {
                //
                //              Path 2 (M much larger than N, JOBU='O', JOBVT='N')
                //              N left singular vectors to be overwritten on A and
                //              no right singular vectors to be computed
                //
                if (lwork >= n * n + 3 * n) {
                    //
                    //                 Sufficient workspace for a fast algorithm
                    //
                    ir = 1;
                    if (lwork >= max(wrkbl, lda * n) + lda * n) {
                        //
                        //                    WORK(IU) is LDA by N, WORK(IR) is LDA by N
                        //
                        ldwrku = lda;
                        ldwrkr = lda;
                    } else if (lwork >= max(wrkbl, lda * n) + n * n) {
                        //
                        //                    WORK(IU) is LDA by N, WORK(IR) is N by N
                        //
                        ldwrku = lda;
                        ldwrkr = n;
                    } else {
                        //
                        //                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N
                        //
                        ldwrku = (lwork - n * n) / n;
                        ldwrkr = n;
                    }
                    itau = ir + ldwrkr * n;
                    iwork = itau + n;
                    //
                    //                 Compute A=Q*R
                    //                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy R to WORK(IR) and zero out below it
                    //
                    Clacpy("U", n, n, a, lda, &work[ir - 1], ldwrkr);
                    Claset("L", n - 1, n - 1, czero, czero, &work[(ir + 1) - 1], ldwrkr);
                    //
                    //                 Generate Q in A
                    //                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungqr(m, n, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = 1;
                    itauq = itau;
                    itaup = itauq + n;
                    iwork = itaup + n;
                    //
                    //                 Bidiagonalize R in WORK(IR)
                    //                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                    //                 (RWorkspace: need N)
                    //
                    Cgebrd(n, n, &work[ir - 1], ldwrkr, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate left vectors bidiagonalizing R
                    //                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                    //                 (RWorkspace: need 0)
                    //
                    Cungbr("Q", n, n, n, &work[ir - 1], ldwrkr, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    irwork = ie + n;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of R in WORK(IR)
                    //                 (CWorkspace: need N*N)
                    //                 (RWorkspace: need BDSPAC)
                    //
                    Cbdsqr("U", n, 0, n, 0, s, &rwork[ie - 1], cdum, 1, &work[ir - 1], ldwrkr, cdum, 1, &rwork[irwork - 1], info);
                    iu = itauq;
                    //
                    //                 Multiply Q in A by left singular vectors of R in
                    //                 WORK(IR), storing result in WORK(IU) and copying to A
                    //                 (CWorkspace: need N*N+N, prefer N*N+M*N)
                    //                 (RWorkspace: 0)
                    //
                    for (i = 1; i <= m; i = i + ldwrku) {
                        chunk = min(m - i + 1, ldwrku);
                        Cgemm("N", "N", chunk, n, n, cone, &a[(i - 1)], lda, &work[ir - 1], ldwrkr, czero, &work[iu - 1], ldwrku);
                        Clacpy("F", chunk, n, &work[iu - 1], ldwrku, &a[(i - 1)], lda);
                    }
                    //
                } else {
                    //
                    //                 Insufficient workspace for a fast algorithm
                    //
                    ie = 1;
                    itauq = 1;
                    itaup = itauq + n;
                    iwork = itaup + n;
                    //
                    //                 Bidiagonalize A
                    //                 (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
                    //                 (RWorkspace: N)
                    //
                    Cgebrd(m, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate left vectors bidiagonalizing A
                    //                 (CWorkspace: need 3*N, prefer 2*N+N*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("Q", m, n, n, a, lda, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    irwork = ie + n;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of A in A
                    //                 (CWorkspace: need 0)
                    //                 (RWorkspace: need BDSPAC)
                    //
                    Cbdsqr("U", n, 0, m, 0, s, &rwork[ie - 1], cdum, 1, a, lda, cdum, 1, &rwork[irwork - 1], info);
                    //
                }
                //
            } else if (wntuo && wntvas) {
                //
                //              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
                //              N left singular vectors to be overwritten on A and
                //              N right singular vectors to be computed in VT
                //
                if (lwork >= n * n + 3 * n) {
                    //
                    //                 Sufficient workspace for a fast algorithm
                    //
                    ir = 1;
                    if (lwork >= max(wrkbl, lda * n) + lda * n) {
                        //
                        //                    WORK(IU) is LDA by N and WORK(IR) is LDA by N
                        //
                        ldwrku = lda;
                        ldwrkr = lda;
                    } else if (lwork >= max(wrkbl, lda * n) + n * n) {
                        //
                        //                    WORK(IU) is LDA by N and WORK(IR) is N by N
                        //
                        ldwrku = lda;
                        ldwrkr = n;
                    } else {
                        //
                        //                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N
                        //
                        ldwrku = (lwork - n * n) / n;
                        ldwrkr = n;
                    }
                    itau = ir + ldwrkr * n;
                    iwork = itau + n;
                    //
                    //                 Compute A=Q*R
                    //                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy R to VT, zeroing out below it
                    //
                    Clacpy("U", n, n, a, lda, vt, ldvt);
                    if (n > 1) {
                        Claset("L", n - 1, n - 1, czero, czero, &vt[(2 - 1)], ldvt);
                    }
                    //
                    //                 Generate Q in A
                    //                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungqr(m, n, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = 1;
                    itauq = itau;
                    itaup = itauq + n;
                    iwork = itaup + n;
                    //
                    //                 Bidiagonalize R in VT, copying result to WORK(IR)
                    //                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                    //                 (RWorkspace: need N)
                    //
                    Cgebrd(n, n, vt, ldvt, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    Clacpy("L", n, n, vt, ldvt, &work[ir - 1], ldwrkr);
                    //
                    //                 Generate left vectors bidiagonalizing R in WORK(IR)
                    //                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("Q", n, n, n, &work[ir - 1], ldwrkr, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate right vectors bidiagonalizing R in VT
                    //                 (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("P", n, n, n, vt, ldvt, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    irwork = ie + n;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of R in WORK(IR) and computing right
                    //                 singular vectors of R in VT
                    //                 (CWorkspace: need N*N)
                    //                 (RWorkspace: need BDSPAC)
                    //
                    Cbdsqr("U", n, n, n, 0, s, &rwork[ie - 1], vt, ldvt, &work[ir - 1], ldwrkr, cdum, 1, &rwork[irwork - 1], info);
                    iu = itauq;
                    //
                    //                 Multiply Q in A by left singular vectors of R in
                    //                 WORK(IR), storing result in WORK(IU) and copying to A
                    //                 (CWorkspace: need N*N+N, prefer N*N+M*N)
                    //                 (RWorkspace: 0)
                    //
                    for (i = 1; i <= m; i = i + ldwrku) {
                        chunk = min(m - i + 1, ldwrku);
                        Cgemm("N", "N", chunk, n, n, cone, &a[(i - 1)], lda, &work[ir - 1], ldwrkr, czero, &work[iu - 1], ldwrku);
                        Clacpy("F", chunk, n, &work[iu - 1], ldwrku, &a[(i - 1)], lda);
                    }
                    //
                } else {
                    //
                    //                 Insufficient workspace for a fast algorithm
                    //
                    itau = 1;
                    iwork = itau + n;
                    //
                    //                 Compute A=Q*R
                    //                 (CWorkspace: need 2*N, prefer N+N*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy R to VT, zeroing out below it
                    //
                    Clacpy("U", n, n, a, lda, vt, ldvt);
                    if (n > 1) {
                        Claset("L", n - 1, n - 1, czero, czero, &vt[(2 - 1)], ldvt);
                    }
                    //
                    //                 Generate Q in A
                    //                 (CWorkspace: need 2*N, prefer N+N*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungqr(m, n, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = 1;
                    itauq = itau;
                    itaup = itauq + n;
                    iwork = itaup + n;
                    //
                    //                 Bidiagonalize R in VT
                    //                 (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                    //                 (RWorkspace: N)
                    //
                    Cgebrd(n, n, vt, ldvt, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Multiply Q in A by left vectors bidiagonalizing R
                    //                 (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cunmbr("Q", "R", "N", m, n, n, vt, ldvt, &work[itauq - 1], a, lda, &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate right vectors bidiagonalizing R in VT
                    //                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("P", n, n, n, vt, ldvt, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    irwork = ie + n;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of A in A and computing right
                    //                 singular vectors of A in VT
                    //                 (CWorkspace: 0)
                    //                 (RWorkspace: need BDSPAC)
                    //
                    Cbdsqr("U", n, n, m, 0, s, &rwork[ie - 1], vt, ldvt, a, lda, cdum, 1, &rwork[irwork - 1], info);
                    //
                }
                //
            } else if (wntus) {
                //
                if (wntvn) {
                    //
                    //                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
                    //                 N left singular vectors to be computed in U and
                    //                 no right singular vectors to be computed
                    //
                    if (lwork >= n * n + 3 * n) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        ir = 1;
                        if (lwork >= wrkbl + lda * n) {
                            //
                            //                       WORK(IR) is LDA by N
                            //
                            ldwrkr = lda;
                        } else {
                            //
                            //                       WORK(IR) is N by N
                            //
                            ldwrkr = n;
                        }
                        itau = ir + ldwrkr * n;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R
                        //                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to WORK(IR), zeroing out below it
                        //
                        Clacpy("U", n, n, a, lda, &work[ir - 1], ldwrkr);
                        Claset("L", n - 1, n - 1, czero, czero, &work[(ir + 1) - 1], ldwrkr);
                        //
                        //                    Generate Q in A
                        //                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, n, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IR)
                        //                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                        //                    (RWorkspace: need N)
                        //
                        Cgebrd(n, n, &work[ir - 1], ldwrkr, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left vectors bidiagonalizing R in WORK(IR)
                        //                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", n, n, n, &work[ir - 1], ldwrkr, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IR)
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, 0, n, 0, s, &rwork[ie - 1], cdum, 1, &work[ir - 1], ldwrkr, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply Q in A by left singular vectors of R in
                        //                    WORK(IR), storing result in U
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, n, cone, a, lda, &work[ir - 1], ldwrkr, czero, u, ldu);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (CWorkspace: need 2*N, prefer N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (CWorkspace: need 2*N, prefer N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, n, n, u, ldu, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Zero out below R in A
                        //
                        if (n > 1) {
                            Claset("L", n - 1, n - 1, czero, czero, &a[(2 - 1)], lda);
                        }
                        //
                        //                    Bidiagonalize R in A
                        //                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                        //                    (RWorkspace: need N)
                        //
                        Cgebrd(n, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left vectors bidiagonalizing R
                        //                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("Q", "R", "N", m, n, n, a, lda, &work[itauq - 1], u, ldu, &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, 0, m, 0, s, &rwork[ie - 1], cdum, 1, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                } else if (wntvo) {
                    //
                    //                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
                    //                 N left singular vectors to be computed in U and
                    //                 N right singular vectors to be overwritten on A
                    //
                    if (lwork >= 2 * n * n + 3 * n) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        iu = 1;
                        if (lwork >= wrkbl + 2 * lda * n) {
                            //
                            //                       WORK(IU) is LDA by N and WORK(IR) is LDA by N
                            //
                            ldwrku = lda;
                            ir = iu + ldwrku * n;
                            ldwrkr = lda;
                        } else if (lwork >= wrkbl + (lda + n) * n) {
                            //
                            //                       WORK(IU) is LDA by N and WORK(IR) is N by N
                            //
                            ldwrku = lda;
                            ir = iu + ldwrku * n;
                            ldwrkr = n;
                        } else {
                            //
                            //                       WORK(IU) is N by N and WORK(IR) is N by N
                            //
                            ldwrku = n;
                            ir = iu + ldwrku * n;
                            ldwrkr = n;
                        }
                        itau = ir + ldwrkr * n;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R
                        //                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to WORK(IU), zeroing out below it
                        //
                        Clacpy("U", n, n, a, lda, &work[iu - 1], ldwrku);
                        Claset("L", n - 1, n - 1, czero, czero, &work[(iu + 1) - 1], ldwrku);
                        //
                        //                    Generate Q in A
                        //                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, n, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IU), copying result to
                        //                    WORK(IR)
                        //                    (CWorkspace: need   2*N*N+3*N,
                        //                                 prefer 2*N*N+2*N+2*N*NB)
                        //                    (RWorkspace: need   N)
                        //
                        Cgebrd(n, n, &work[iu - 1], ldwrku, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", n, n, &work[iu - 1], ldwrku, &work[ir - 1], ldwrkr);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IU)
                        //                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", n, n, n, &work[iu - 1], ldwrku, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IR)
                        //                    (CWorkspace: need   2*N*N+3*N-1,
                        //                                 prefer 2*N*N+2*N+(N-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", n, n, n, &work[ir - 1], ldwrkr, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IU) and computing
                        //                    right singular vectors of R in WORK(IR)
                        //                    (CWorkspace: need 2*N*N)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, n, n, 0, s, &rwork[ie - 1], &work[ir - 1], ldwrkr, &work[iu - 1], ldwrku, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply Q in A by left singular vectors of R in
                        //                    WORK(IU), storing result in U
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, n, cone, a, lda, &work[iu - 1], ldwrku, czero, u, ldu);
                        //
                        //                    Copy right singular vectors of R to A
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: 0)
                        //
                        Clacpy("F", n, n, &work[ir - 1], ldwrkr, a, lda);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (CWorkspace: need 2*N, prefer N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (CWorkspace: need 2*N, prefer N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, n, n, u, ldu, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Zero out below R in A
                        //
                        if (n > 1) {
                            Claset("L", n - 1, n - 1, czero, czero, &a[(2 - 1)], lda);
                        }
                        //
                        //                    Bidiagonalize R in A
                        //                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                        //                    (RWorkspace: need N)
                        //
                        Cgebrd(n, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left vectors bidiagonalizing R
                        //                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("Q", "R", "N", m, n, n, a, lda, &work[itauq - 1], u, ldu, &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right vectors bidiagonalizing R in A
                        //                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", n, n, n, a, lda, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in A
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, n, m, 0, s, &rwork[ie - 1], a, lda, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                } else if (wntvas) {
                    //
                    //                 Path 6 (M much larger than N, JOBU='S', JOBVT='S'
                    //                         or 'A')
                    //                 N left singular vectors to be computed in U and
                    //                 N right singular vectors to be computed in VT
                    //
                    if (lwork >= n * n + 3 * n) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        iu = 1;
                        if (lwork >= wrkbl + lda * n) {
                            //
                            //                       WORK(IU) is LDA by N
                            //
                            ldwrku = lda;
                        } else {
                            //
                            //                       WORK(IU) is N by N
                            //
                            ldwrku = n;
                        }
                        itau = iu + ldwrku * n;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R
                        //                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to WORK(IU), zeroing out below it
                        //
                        Clacpy("U", n, n, a, lda, &work[iu - 1], ldwrku);
                        Claset("L", n - 1, n - 1, czero, czero, &work[(iu + 1) - 1], ldwrku);
                        //
                        //                    Generate Q in A
                        //                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, n, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IU), copying result to VT
                        //                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                        //                    (RWorkspace: need N)
                        //
                        Cgebrd(n, n, &work[iu - 1], ldwrku, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", n, n, &work[iu - 1], ldwrku, vt, ldvt);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IU)
                        //                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", n, n, n, &work[iu - 1], ldwrku, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in VT
                        //                    (CWorkspace: need   N*N+3*N-1,
                        //                                 prefer N*N+2*N+(N-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", n, n, n, vt, ldvt, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IU) and computing
                        //                    right singular vectors of R in VT
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, n, n, 0, s, &rwork[ie - 1], vt, ldvt, &work[iu - 1], ldwrku, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply Q in A by left singular vectors of R in
                        //                    WORK(IU), storing result in U
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, n, cone, a, lda, &work[iu - 1], ldwrku, czero, u, ldu);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (CWorkspace: need 2*N, prefer N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (CWorkspace: need 2*N, prefer N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, n, n, u, ldu, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to VT, zeroing out below it
                        //
                        Clacpy("U", n, n, a, lda, vt, ldvt);
                        if (n > 1) {
                            Claset("L", n - 1, n - 1, czero, czero, &vt[(2 - 1)], ldvt);
                        }
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in VT
                        //                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                        //                    (RWorkspace: need N)
                        //
                        Cgebrd(n, n, vt, ldvt, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left bidiagonalizing vectors
                        //                    in VT
                        //                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("Q", "R", "N", m, n, n, vt, ldvt, &work[itauq - 1], u, ldu, &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in VT
                        //                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", n, n, n, vt, ldvt, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in VT
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, n, m, 0, s, &rwork[ie - 1], vt, ldvt, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                }
                //
            } else if (wntua) {
                //
                if (wntvn) {
                    //
                    //                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
                    //                 M left singular vectors to be computed in U and
                    //                 no right singular vectors to be computed
                    //
                    if (lwork >= n * n + max(n + m, 3 * n)) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        ir = 1;
                        if (lwork >= wrkbl + lda * n) {
                            //
                            //                       WORK(IR) is LDA by N
                            //
                            ldwrkr = lda;
                        } else {
                            //
                            //                       WORK(IR) is N by N
                            //
                            ldwrkr = n;
                        }
                        itau = ir + ldwrkr * n;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Copy R to WORK(IR), zeroing out below it
                        //
                        Clacpy("U", n, n, a, lda, &work[ir - 1], ldwrkr);
                        Claset("L", n - 1, n - 1, czero, czero, &work[(ir + 1) - 1], ldwrkr);
                        //
                        //                    Generate Q in U
                        //                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, m, n, u, ldu, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IR)
                        //                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                        //                    (RWorkspace: need N)
                        //
                        Cgebrd(n, n, &work[ir - 1], ldwrkr, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IR)
                        //                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", n, n, n, &work[ir - 1], ldwrkr, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IR)
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, 0, n, 0, s, &rwork[ie - 1], cdum, 1, &work[ir - 1], ldwrkr, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply Q in U by left singular vectors of R in
                        //                    WORK(IR), storing result in A
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, n, cone, u, ldu, &work[ir - 1], ldwrkr, czero, a, lda);
                        //
                        //                    Copy left singular vectors of A from A to U
                        //
                        Clacpy("F", m, n, a, lda, u, ldu);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (CWorkspace: need 2*N, prefer N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (CWorkspace: need N+M, prefer N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, m, n, u, ldu, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Zero out below R in A
                        //
                        if (n > 1) {
                            Claset("L", n - 1, n - 1, czero, czero, &a[(2 - 1)], lda);
                        }
                        //
                        //                    Bidiagonalize R in A
                        //                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                        //                    (RWorkspace: need N)
                        //
                        Cgebrd(n, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left bidiagonalizing vectors
                        //                    in A
                        //                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("Q", "R", "N", m, n, n, a, lda, &work[itauq - 1], u, ldu, &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, 0, m, 0, s, &rwork[ie - 1], cdum, 1, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                } else if (wntvo) {
                    //
                    //                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
                    //                 M left singular vectors to be computed in U and
                    //                 N right singular vectors to be overwritten on A
                    //
                    if (lwork >= 2 * n * n + max(n + m, 3 * n)) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        iu = 1;
                        if (lwork >= wrkbl + 2 * lda * n) {
                            //
                            //                       WORK(IU) is LDA by N and WORK(IR) is LDA by N
                            //
                            ldwrku = lda;
                            ir = iu + ldwrku * n;
                            ldwrkr = lda;
                        } else if (lwork >= wrkbl + (lda + n) * n) {
                            //
                            //                       WORK(IU) is LDA by N and WORK(IR) is N by N
                            //
                            ldwrku = lda;
                            ir = iu + ldwrku * n;
                            ldwrkr = n;
                        } else {
                            //
                            //                       WORK(IU) is N by N and WORK(IR) is N by N
                            //
                            ldwrku = n;
                            ir = iu + ldwrku * n;
                            ldwrkr = n;
                        }
                        itau = ir + ldwrkr * n;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, m, n, u, ldu, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to WORK(IU), zeroing out below it
                        //
                        Clacpy("U", n, n, a, lda, &work[iu - 1], ldwrku);
                        Claset("L", n - 1, n - 1, czero, czero, &work[(iu + 1) - 1], ldwrku);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IU), copying result to
                        //                    WORK(IR)
                        //                    (CWorkspace: need   2*N*N+3*N,
                        //                                 prefer 2*N*N+2*N+2*N*NB)
                        //                    (RWorkspace: need   N)
                        //
                        Cgebrd(n, n, &work[iu - 1], ldwrku, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", n, n, &work[iu - 1], ldwrku, &work[ir - 1], ldwrkr);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IU)
                        //                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", n, n, n, &work[iu - 1], ldwrku, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IR)
                        //                    (CWorkspace: need   2*N*N+3*N-1,
                        //                                 prefer 2*N*N+2*N+(N-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", n, n, n, &work[ir - 1], ldwrkr, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IU) and computing
                        //                    right singular vectors of R in WORK(IR)
                        //                    (CWorkspace: need 2*N*N)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, n, n, 0, s, &rwork[ie - 1], &work[ir - 1], ldwrkr, &work[iu - 1], ldwrku, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply Q in U by left singular vectors of R in
                        //                    WORK(IU), storing result in A
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, n, cone, u, ldu, &work[iu - 1], ldwrku, czero, a, lda);
                        //
                        //                    Copy left singular vectors of A from A to U
                        //
                        Clacpy("F", m, n, a, lda, u, ldu);
                        //
                        //                    Copy right singular vectors of R from WORK(IR) to A
                        //
                        Clacpy("F", n, n, &work[ir - 1], ldwrkr, a, lda);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (CWorkspace: need 2*N, prefer N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (CWorkspace: need N+M, prefer N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, m, n, u, ldu, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Zero out below R in A
                        //
                        if (n > 1) {
                            Claset("L", n - 1, n - 1, czero, czero, &a[(2 - 1)], lda);
                        }
                        //
                        //                    Bidiagonalize R in A
                        //                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                        //                    (RWorkspace: need N)
                        //
                        Cgebrd(n, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left bidiagonalizing vectors
                        //                    in A
                        //                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("Q", "R", "N", m, n, n, a, lda, &work[itauq - 1], u, ldu, &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in A
                        //                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", n, n, n, a, lda, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in A
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, n, m, 0, s, &rwork[ie - 1], a, lda, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                } else if (wntvas) {
                    //
                    //                 Path 9 (M much larger than N, JOBU='A', JOBVT='S'
                    //                         or 'A')
                    //                 M left singular vectors to be computed in U and
                    //                 N right singular vectors to be computed in VT
                    //
                    if (lwork >= n * n + max(n + m, 3 * n)) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        iu = 1;
                        if (lwork >= wrkbl + lda * n) {
                            //
                            //                       WORK(IU) is LDA by N
                            //
                            ldwrku = lda;
                        } else {
                            //
                            //                       WORK(IU) is N by N
                            //
                            ldwrku = n;
                        }
                        itau = iu + ldwrku * n;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, m, n, u, ldu, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to WORK(IU), zeroing out below it
                        //
                        Clacpy("U", n, n, a, lda, &work[iu - 1], ldwrku);
                        Claset("L", n - 1, n - 1, czero, czero, &work[(iu + 1) - 1], ldwrku);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IU), copying result to VT
                        //                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                        //                    (RWorkspace: need N)
                        //
                        Cgebrd(n, n, &work[iu - 1], ldwrku, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", n, n, &work[iu - 1], ldwrku, vt, ldvt);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IU)
                        //                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", n, n, n, &work[iu - 1], ldwrku, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in VT
                        //                    (CWorkspace: need   N*N+3*N-1,
                        //                                 prefer N*N+2*N+(N-1)*NB)
                        //                    (RWorkspace: need   0)
                        //
                        Cungbr("P", n, n, n, vt, ldvt, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IU) and computing
                        //                    right singular vectors of R in VT
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, n, n, 0, s, &rwork[ie - 1], vt, ldvt, &work[iu - 1], ldwrku, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply Q in U by left singular vectors of R in
                        //                    WORK(IU), storing result in A
                        //                    (CWorkspace: need N*N)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, n, cone, u, ldu, &work[iu - 1], ldwrku, czero, a, lda);
                        //
                        //                    Copy left singular vectors of A from A to U
                        //
                        Clacpy("F", m, n, a, lda, u, ldu);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (CWorkspace: need 2*N, prefer N+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgeqrf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (CWorkspace: need N+M, prefer N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungqr(m, m, n, u, ldu, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R from A to VT, zeroing out below it
                        //
                        Clacpy("U", n, n, a, lda, vt, ldvt);
                        if (n > 1) {
                            Claset("L", n - 1, n - 1, czero, czero, &vt[(2 - 1)], ldvt);
                        }
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in VT
                        //                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                        //                    (RWorkspace: need N)
                        //
                        Cgebrd(n, n, vt, ldvt, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left bidiagonalizing vectors
                        //                    in VT
                        //                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("Q", "R", "N", m, n, n, vt, ldvt, &work[itauq - 1], u, ldu, &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in VT
                        //                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", n, n, n, vt, ldvt, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in VT
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", n, n, m, 0, s, &rwork[ie - 1], vt, ldvt, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                }
                //
            }
            //
        } else {
            //
            //           M .LT. MNTHR
            //
            //           Path 10 (M at least N, but not much larger)
            //           Reduce to bidiagonal form without QR decomposition
            //
            ie = 1;
            itauq = 1;
            itaup = itauq + n;
            iwork = itaup + n;
            //
            //           Bidiagonalize A
            //           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
            //           (RWorkspace: need N)
            //
            Cgebrd(m, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
            if (wntuas) {
                //
                //              If left singular vectors desired in U, copy result to U
                //              and generate left bidiagonalizing vectors in U
                //              (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB)
                //              (RWorkspace: 0)
                //
                Clacpy("L", m, n, a, lda, u, ldu);
                if (wntus) {
                    ncu = n;
                }
                if (wntua) {
                    ncu = m;
                }
                Cungbr("Q", m, ncu, n, u, ldu, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntvas) {
                //
                //              If right singular vectors desired in VT, copy result to
                //              VT and generate right bidiagonalizing vectors in VT
                //              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                //              (RWorkspace: 0)
                //
                Clacpy("U", n, n, a, lda, vt, ldvt);
                Cungbr("P", n, n, n, vt, ldvt, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntuo) {
                //
                //              If left singular vectors desired in A, generate left
                //              bidiagonalizing vectors in A
                //              (CWorkspace: need 3*N, prefer 2*N+N*NB)
                //              (RWorkspace: 0)
                //
                Cungbr("Q", m, n, n, a, lda, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntvo) {
                //
                //              If right singular vectors desired in A, generate right
                //              bidiagonalizing vectors in A
                //              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                //              (RWorkspace: 0)
                //
                Cungbr("P", n, n, n, a, lda, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
            }
            irwork = ie + n;
            if (wntuas || wntuo) {
                nru = m;
            }
            if (wntun) {
                nru = 0;
            }
            if (wntvas || wntvo) {
                ncvt = n;
            }
            if (wntvn) {
                ncvt = 0;
            }
            if ((!wntuo) && (!wntvo)) {
                //
                //              Perform bidiagonal QR iteration, if desired, computing
                //              left singular vectors in U and computing right singular
                //              vectors in VT
                //              (CWorkspace: 0)
                //              (RWorkspace: need BDSPAC)
                //
                Cbdsqr("U", n, ncvt, nru, 0, s, &rwork[ie - 1], vt, ldvt, u, ldu, cdum, 1, &rwork[irwork - 1], info);
            } else if ((!wntuo) && wntvo) {
                //
                //              Perform bidiagonal QR iteration, if desired, computing
                //              left singular vectors in U and computing right singular
                //              vectors in A
                //              (CWorkspace: 0)
                //              (RWorkspace: need BDSPAC)
                //
                Cbdsqr("U", n, ncvt, nru, 0, s, &rwork[ie - 1], a, lda, u, ldu, cdum, 1, &rwork[irwork - 1], info);
            } else {
                //
                //              Perform bidiagonal QR iteration, if desired, computing
                //              left singular vectors in A and computing right singular
                //              vectors in VT
                //              (CWorkspace: 0)
                //              (RWorkspace: need BDSPAC)
                //
                Cbdsqr("U", n, ncvt, nru, 0, s, &rwork[ie - 1], vt, ldvt, a, lda, cdum, 1, &rwork[irwork - 1], info);
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
            if (wntvn) {
                //
                //              Path 1t(N much larger than M, JOBVT='N')
                //              No right singular vectors to be computed
                //
                itau = 1;
                iwork = itau + m;
                //
                //              Compute A=L*Q
                //              (CWorkspace: need 2*M, prefer M+M*NB)
                //              (RWorkspace: 0)
                //
                Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                //
                //              Zero out above L
                //
                Claset("U", m - 1, m - 1, czero, czero, &a[(2 - 1) * lda], lda);
                ie = 1;
                itauq = 1;
                itaup = itauq + m;
                iwork = itaup + m;
                //
                //              Bidiagonalize L in A
                //              (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                //              (RWorkspace: need M)
                //
                Cgebrd(m, m, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                if (wntuo || wntuas) {
                    //
                    //                 If left singular vectors desired, generate Q
                    //                 (CWorkspace: need 3*M, prefer 2*M+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("Q", m, m, m, a, lda, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                }
                irwork = ie + m;
                nru = 0;
                if (wntuo || wntuas) {
                    nru = m;
                }
                //
                //              Perform bidiagonal QR iteration, computing left singular
                //              vectors of A in A if desired
                //              (CWorkspace: 0)
                //              (RWorkspace: need BDSPAC)
                //
                Cbdsqr("U", m, 0, nru, 0, s, &rwork[ie - 1], cdum, 1, a, lda, cdum, 1, &rwork[irwork - 1], info);
                //
                //              If left singular vectors desired in U, copy them there
                //
                if (wntuas) {
                    Clacpy("F", m, m, a, lda, u, ldu);
                }
                //
            } else if (wntvo && wntun) {
                //
                //              Path 2t(N much larger than M, JOBU='N', JOBVT='O')
                //              M right singular vectors to be overwritten on A and
                //              no left singular vectors to be computed
                //
                if (lwork >= m * m + 3 * m) {
                    //
                    //                 Sufficient workspace for a fast algorithm
                    //
                    ir = 1;
                    if (lwork >= max(wrkbl, lda * n) + lda * m) {
                        //
                        //                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
                        //
                        ldwrku = lda;
                        chunk = n;
                        ldwrkr = lda;
                    } else if (lwork >= max(wrkbl, lda * n) + m * m) {
                        //
                        //                    WORK(IU) is LDA by N and WORK(IR) is M by M
                        //
                        ldwrku = lda;
                        chunk = n;
                        ldwrkr = m;
                    } else {
                        //
                        //                    WORK(IU) is M by CHUNK and WORK(IR) is M by M
                        //
                        ldwrku = m;
                        chunk = (lwork - m * m) / m;
                        ldwrkr = m;
                    }
                    itau = ir + ldwrkr * m;
                    iwork = itau + m;
                    //
                    //                 Compute A=L*Q
                    //                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy L to WORK(IR) and zero out above it
                    //
                    Clacpy("L", m, m, a, lda, &work[ir - 1], ldwrkr);
                    Claset("U", m - 1, m - 1, czero, czero, &work[(ir + ldwrkr) - 1], ldwrkr);
                    //
                    //                 Generate Q in A
                    //                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cunglq(m, n, m, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = 1;
                    itauq = itau;
                    itaup = itauq + m;
                    iwork = itaup + m;
                    //
                    //                 Bidiagonalize L in WORK(IR)
                    //                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                    //                 (RWorkspace: need M)
                    //
                    Cgebrd(m, m, &work[ir - 1], ldwrkr, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate right vectors bidiagonalizing L
                    //                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("P", m, m, m, &work[ir - 1], ldwrkr, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    irwork = ie + m;
                    //
                    //                 Perform bidiagonal QR iteration, computing right
                    //                 singular vectors of L in WORK(IR)
                    //                 (CWorkspace: need M*M)
                    //                 (RWorkspace: need BDSPAC)
                    //
                    Cbdsqr("U", m, m, 0, 0, s, &rwork[ie - 1], &work[ir - 1], ldwrkr, cdum, 1, cdum, 1, &rwork[irwork - 1], info);
                    iu = itauq;
                    //
                    //                 Multiply right singular vectors of L in WORK(IR) by Q
                    //                 in A, storing result in WORK(IU) and copying to A
                    //                 (CWorkspace: need M*M+M, prefer M*M+M*N)
                    //                 (RWorkspace: 0)
                    //
                    for (i = 1; i <= n; i = i + chunk) {
                        blk = min(n - i + 1, chunk);
                        Cgemm("N", "N", m, blk, m, cone, &work[ir - 1], ldwrkr, &a[(i - 1) * lda], lda, czero, &work[iu - 1], ldwrku);
                        Clacpy("F", m, blk, &work[iu - 1], ldwrku, &a[(i - 1) * lda], lda);
                    }
                    //
                } else {
                    //
                    //                 Insufficient workspace for a fast algorithm
                    //
                    ie = 1;
                    itauq = 1;
                    itaup = itauq + m;
                    iwork = itaup + m;
                    //
                    //                 Bidiagonalize A
                    //                 (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
                    //                 (RWorkspace: need M)
                    //
                    Cgebrd(m, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate right vectors bidiagonalizing A
                    //                 (CWorkspace: need 3*M, prefer 2*M+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("P", m, n, m, a, lda, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    irwork = ie + m;
                    //
                    //                 Perform bidiagonal QR iteration, computing right
                    //                 singular vectors of A in A
                    //                 (CWorkspace: 0)
                    //                 (RWorkspace: need BDSPAC)
                    //
                    Cbdsqr("L", m, n, 0, 0, s, &rwork[ie - 1], a, lda, cdum, 1, cdum, 1, &rwork[irwork - 1], info);
                    //
                }
                //
            } else if (wntvo && wntuas) {
                //
                //              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
                //              M right singular vectors to be overwritten on A and
                //              M left singular vectors to be computed in U
                //
                if (lwork >= m * m + 3 * m) {
                    //
                    //                 Sufficient workspace for a fast algorithm
                    //
                    ir = 1;
                    if (lwork >= max(wrkbl, lda * n) + lda * m) {
                        //
                        //                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
                        //
                        ldwrku = lda;
                        chunk = n;
                        ldwrkr = lda;
                    } else if (lwork >= max(wrkbl, lda * n) + m * m) {
                        //
                        //                    WORK(IU) is LDA by N and WORK(IR) is M by M
                        //
                        ldwrku = lda;
                        chunk = n;
                        ldwrkr = m;
                    } else {
                        //
                        //                    WORK(IU) is M by CHUNK and WORK(IR) is M by M
                        //
                        ldwrku = m;
                        chunk = (lwork - m * m) / m;
                        ldwrkr = m;
                    }
                    itau = ir + ldwrkr * m;
                    iwork = itau + m;
                    //
                    //                 Compute A=L*Q
                    //                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy L to U, zeroing about above it
                    //
                    Clacpy("L", m, m, a, lda, u, ldu);
                    Claset("U", m - 1, m - 1, czero, czero, &u[(2 - 1) * ldu], ldu);
                    //
                    //                 Generate Q in A
                    //                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cunglq(m, n, m, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = 1;
                    itauq = itau;
                    itaup = itauq + m;
                    iwork = itaup + m;
                    //
                    //                 Bidiagonalize L in U, copying result to WORK(IR)
                    //                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                    //                 (RWorkspace: need M)
                    //
                    Cgebrd(m, m, u, ldu, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    Clacpy("U", m, m, u, ldu, &work[ir - 1], ldwrkr);
                    //
                    //                 Generate right vectors bidiagonalizing L in WORK(IR)
                    //                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("P", m, m, m, &work[ir - 1], ldwrkr, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate left vectors bidiagonalizing L in U
                    //                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("Q", m, m, m, u, ldu, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    irwork = ie + m;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of L in U, and computing right
                    //                 singular vectors of L in WORK(IR)
                    //                 (CWorkspace: need M*M)
                    //                 (RWorkspace: need BDSPAC)
                    //
                    Cbdsqr("U", m, m, m, 0, s, &rwork[ie - 1], &work[ir - 1], ldwrkr, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                    iu = itauq;
                    //
                    //                 Multiply right singular vectors of L in WORK(IR) by Q
                    //                 in A, storing result in WORK(IU) and copying to A
                    //                 (CWorkspace: need M*M+M, prefer M*M+M*N))
                    //                 (RWorkspace: 0)
                    //
                    for (i = 1; i <= n; i = i + chunk) {
                        blk = min(n - i + 1, chunk);
                        Cgemm("N", "N", m, blk, m, cone, &work[ir - 1], ldwrkr, &a[(i - 1) * lda], lda, czero, &work[iu - 1], ldwrku);
                        Clacpy("F", m, blk, &work[iu - 1], ldwrku, &a[(i - 1) * lda], lda);
                    }
                    //
                } else {
                    //
                    //                 Insufficient workspace for a fast algorithm
                    //
                    itau = 1;
                    iwork = itau + m;
                    //
                    //                 Compute A=L*Q
                    //                 (CWorkspace: need 2*M, prefer M+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy L to U, zeroing out above it
                    //
                    Clacpy("L", m, m, a, lda, u, ldu);
                    Claset("U", m - 1, m - 1, czero, czero, &u[(2 - 1) * ldu], ldu);
                    //
                    //                 Generate Q in A
                    //                 (CWorkspace: need 2*M, prefer M+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cunglq(m, n, m, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = 1;
                    itauq = itau;
                    itaup = itauq + m;
                    iwork = itaup + m;
                    //
                    //                 Bidiagonalize L in U
                    //                 (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                    //                 (RWorkspace: need M)
                    //
                    Cgebrd(m, m, u, ldu, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Multiply right vectors bidiagonalizing L by Q in A
                    //                 (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cunmbr("P", "L", "C", m, n, m, u, ldu, &work[itaup - 1], a, lda, &work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate left vectors bidiagonalizing L in U
                    //                 (CWorkspace: need 3*M, prefer 2*M+M*NB)
                    //                 (RWorkspace: 0)
                    //
                    Cungbr("Q", m, m, m, u, ldu, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                    irwork = ie + m;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of A in U and computing right
                    //                 singular vectors of A in A
                    //                 (CWorkspace: 0)
                    //                 (RWorkspace: need BDSPAC)
                    //
                    Cbdsqr("U", m, n, m, 0, s, &rwork[ie - 1], a, lda, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                    //
                }
                //
            } else if (wntvs) {
                //
                if (wntun) {
                    //
                    //                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
                    //                 M right singular vectors to be computed in VT and
                    //                 no left singular vectors to be computed
                    //
                    if (lwork >= m * m + 3 * m) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        ir = 1;
                        if (lwork >= wrkbl + lda * m) {
                            //
                            //                       WORK(IR) is LDA by M
                            //
                            ldwrkr = lda;
                        } else {
                            //
                            //                       WORK(IR) is M by M
                            //
                            ldwrkr = m;
                        }
                        itau = ir + ldwrkr * m;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q
                        //                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to WORK(IR), zeroing out above it
                        //
                        Clacpy("L", m, m, a, lda, &work[ir - 1], ldwrkr);
                        Claset("U", m - 1, m - 1, czero, czero, &work[(ir + ldwrkr) - 1], ldwrkr);
                        //
                        //                    Generate Q in A
                        //                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(m, n, m, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IR)
                        //                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                        //                    (RWorkspace: need M)
                        //
                        Cgebrd(m, m, &work[ir - 1], ldwrkr, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right vectors bidiagonalizing L in
                        //                    WORK(IR)
                        //                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", m, m, m, &work[ir - 1], ldwrkr, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing right
                        //                    singular vectors of L in WORK(IR)
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, m, 0, 0, s, &rwork[ie - 1], &work[ir - 1], ldwrkr, cdum, 1, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IR) by
                        //                    Q in A, storing result in VT
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, m, cone, &work[ir - 1], ldwrkr, a, lda, czero, vt, ldvt);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q
                        //                    (CWorkspace: need 2*M, prefer M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy result to VT
                        //
                        Clacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (CWorkspace: need 2*M, prefer M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(m, n, m, vt, ldvt, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Zero out above L in A
                        //
                        Claset("U", m - 1, m - 1, czero, czero, &a[(2 - 1) * lda], lda);
                        //
                        //                    Bidiagonalize L in A
                        //                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                        //                    (RWorkspace: need M)
                        //
                        Cgebrd(m, m, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right vectors bidiagonalizing L by Q in VT
                        //                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("P", "L", "C", m, n, m, a, lda, &work[itaup - 1], vt, ldvt, &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing right
                        //                    singular vectors of A in VT
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, n, 0, 0, s, &rwork[ie - 1], vt, ldvt, cdum, 1, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                } else if (wntuo) {
                    //
                    //                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
                    //                 M right singular vectors to be computed in VT and
                    //                 M left singular vectors to be overwritten on A
                    //
                    if (lwork >= 2 * m * m + 3 * m) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        iu = 1;
                        if (lwork >= wrkbl + 2 * lda * m) {
                            //
                            //                       WORK(IU) is LDA by M and WORK(IR) is LDA by M
                            //
                            ldwrku = lda;
                            ir = iu + ldwrku * m;
                            ldwrkr = lda;
                        } else if (lwork >= wrkbl + (lda + m) * m) {
                            //
                            //                       WORK(IU) is LDA by M and WORK(IR) is M by M
                            //
                            ldwrku = lda;
                            ir = iu + ldwrku * m;
                            ldwrkr = m;
                        } else {
                            //
                            //                       WORK(IU) is M by M and WORK(IR) is M by M
                            //
                            ldwrku = m;
                            ir = iu + ldwrku * m;
                            ldwrkr = m;
                        }
                        itau = ir + ldwrkr * m;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q
                        //                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to WORK(IU), zeroing out below it
                        //
                        Clacpy("L", m, m, a, lda, &work[iu - 1], ldwrku);
                        Claset("U", m - 1, m - 1, czero, czero, &work[(iu + ldwrku) - 1], ldwrku);
                        //
                        //                    Generate Q in A
                        //                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(m, n, m, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IU), copying result to
                        //                    WORK(IR)
                        //                    (CWorkspace: need   2*M*M+3*M,
                        //                                 prefer 2*M*M+2*M+2*M*NB)
                        //                    (RWorkspace: need   M)
                        //
                        Cgebrd(m, m, &work[iu - 1], ldwrku, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, m, &work[iu - 1], ldwrku, &work[ir - 1], ldwrkr);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IU)
                        //                    (CWorkspace: need   2*M*M+3*M-1,
                        //                                 prefer 2*M*M+2*M+(M-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", m, m, m, &work[iu - 1], ldwrku, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IR)
                        //                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", m, m, m, &work[ir - 1], ldwrkr, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of L in WORK(IR) and computing
                        //                    right singular vectors of L in WORK(IU)
                        //                    (CWorkspace: need 2*M*M)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, m, m, 0, s, &rwork[ie - 1], &work[iu - 1], ldwrku, &work[ir - 1], ldwrkr, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IU) by
                        //                    Q in A, storing result in VT
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, m, cone, &work[iu - 1], ldwrku, a, lda, czero, vt, ldvt);
                        //
                        //                    Copy left singular vectors of L to A
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: 0)
                        //
                        Clacpy("F", m, m, &work[ir - 1], ldwrkr, a, lda);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (CWorkspace: need 2*M, prefer M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (CWorkspace: need 2*M, prefer M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(m, n, m, vt, ldvt, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Zero out above L in A
                        //
                        Claset("U", m - 1, m - 1, czero, czero, &a[(2 - 1) * lda], lda);
                        //
                        //                    Bidiagonalize L in A
                        //                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                        //                    (RWorkspace: need M)
                        //
                        Cgebrd(m, m, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right vectors bidiagonalizing L by Q in VT
                        //                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("P", "L", "C", m, n, m, a, lda, &work[itaup - 1], vt, ldvt, &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors of L in A
                        //                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", m, m, m, a, lda, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in A and computing right
                        //                    singular vectors of A in VT
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, n, m, 0, s, &rwork[ie - 1], vt, ldvt, a, lda, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                } else if (wntuas) {
                    //
                    //                 Path 6t(N much larger than M, JOBU='S' or 'A',
                    //                         JOBVT='S')
                    //                 M right singular vectors to be computed in VT and
                    //                 M left singular vectors to be computed in U
                    //
                    if (lwork >= m * m + 3 * m) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        iu = 1;
                        if (lwork >= wrkbl + lda * m) {
                            //
                            //                       WORK(IU) is LDA by N
                            //
                            ldwrku = lda;
                        } else {
                            //
                            //                       WORK(IU) is LDA by M
                            //
                            ldwrku = m;
                        }
                        itau = iu + ldwrku * m;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q
                        //                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to WORK(IU), zeroing out above it
                        //
                        Clacpy("L", m, m, a, lda, &work[iu - 1], ldwrku);
                        Claset("U", m - 1, m - 1, czero, czero, &work[(iu + ldwrku) - 1], ldwrku);
                        //
                        //                    Generate Q in A
                        //                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(m, n, m, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IU), copying result to U
                        //                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                        //                    (RWorkspace: need M)
                        //
                        Cgebrd(m, m, &work[iu - 1], ldwrku, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, m, &work[iu - 1], ldwrku, u, ldu);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IU)
                        //                    (CWorkspace: need   M*M+3*M-1,
                        //                                 prefer M*M+2*M+(M-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", m, m, m, &work[iu - 1], ldwrku, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in U
                        //                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", m, m, m, u, ldu, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of L in U and computing right
                        //                    singular vectors of L in WORK(IU)
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, m, m, 0, s, &rwork[ie - 1], &work[iu - 1], ldwrku, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IU) by
                        //                    Q in A, storing result in VT
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, m, cone, &work[iu - 1], ldwrku, a, lda, czero, vt, ldvt);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (CWorkspace: need 2*M, prefer M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (CWorkspace: need 2*M, prefer M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(m, n, m, vt, ldvt, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to U, zeroing out above it
                        //
                        Clacpy("L", m, m, a, lda, u, ldu);
                        Claset("U", m - 1, m - 1, czero, czero, &u[(2 - 1) * ldu], ldu);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in U
                        //                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                        //                    (RWorkspace: need M)
                        //
                        Cgebrd(m, m, u, ldu, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right bidiagonalizing vectors in U by Q
                        //                    in VT
                        //                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("P", "L", "C", m, n, m, u, ldu, &work[itaup - 1], vt, ldvt, &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in U
                        //                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", m, m, m, u, ldu, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in VT
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, n, m, 0, s, &rwork[ie - 1], vt, ldvt, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                }
                //
            } else if (wntva) {
                //
                if (wntun) {
                    //
                    //                 Path 7t(N much larger than M, JOBU='N', JOBVT='A')
                    //                 N right singular vectors to be computed in VT and
                    //                 no left singular vectors to be computed
                    //
                    if (lwork >= m * m + max(n + m, 3 * m)) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        ir = 1;
                        if (lwork >= wrkbl + lda * m) {
                            //
                            //                       WORK(IR) is LDA by M
                            //
                            ldwrkr = lda;
                        } else {
                            //
                            //                       WORK(IR) is M by M
                            //
                            ldwrkr = m;
                        }
                        itau = ir + ldwrkr * m;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Copy L to WORK(IR), zeroing out above it
                        //
                        Clacpy("L", m, m, a, lda, &work[ir - 1], ldwrkr);
                        Claset("U", m - 1, m - 1, czero, czero, &work[(ir + ldwrkr) - 1], ldwrkr);
                        //
                        //                    Generate Q in VT
                        //                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(n, n, m, vt, ldvt, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IR)
                        //                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                        //                    (RWorkspace: need M)
                        //
                        Cgebrd(m, m, &work[ir - 1], ldwrkr, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IR)
                        //                    (CWorkspace: need   M*M+3*M-1,
                        //                                 prefer M*M+2*M+(M-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", m, m, m, &work[ir - 1], ldwrkr, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing right
                        //                    singular vectors of L in WORK(IR)
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, m, 0, 0, s, &rwork[ie - 1], &work[ir - 1], ldwrkr, cdum, 1, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IR) by
                        //                    Q in VT, storing result in A
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, m, cone, &work[ir - 1], ldwrkr, vt, ldvt, czero, a, lda);
                        //
                        //                    Copy right singular vectors of A from A to VT
                        //
                        Clacpy("F", m, n, a, lda, vt, ldvt);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (CWorkspace: need 2*M, prefer M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (CWorkspace: need M+N, prefer M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(n, n, m, vt, ldvt, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Zero out above L in A
                        //
                        Claset("U", m - 1, m - 1, czero, czero, &a[(2 - 1) * lda], lda);
                        //
                        //                    Bidiagonalize L in A
                        //                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                        //                    (RWorkspace: need M)
                        //
                        Cgebrd(m, m, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right bidiagonalizing vectors in A by Q
                        //                    in VT
                        //                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("P", "L", "C", m, n, m, a, lda, &work[itaup - 1], vt, ldvt, &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing right
                        //                    singular vectors of A in VT
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, n, 0, 0, s, &rwork[ie - 1], vt, ldvt, cdum, 1, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                } else if (wntuo) {
                    //
                    //                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
                    //                 N right singular vectors to be computed in VT and
                    //                 M left singular vectors to be overwritten on A
                    //
                    if (lwork >= 2 * m * m + max(n + m, 3 * m)) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        iu = 1;
                        if (lwork >= wrkbl + 2 * lda * m) {
                            //
                            //                       WORK(IU) is LDA by M and WORK(IR) is LDA by M
                            //
                            ldwrku = lda;
                            ir = iu + ldwrku * m;
                            ldwrkr = lda;
                        } else if (lwork >= wrkbl + (lda + m) * m) {
                            //
                            //                       WORK(IU) is LDA by M and WORK(IR) is M by M
                            //
                            ldwrku = lda;
                            ir = iu + ldwrku * m;
                            ldwrkr = m;
                        } else {
                            //
                            //                       WORK(IU) is M by M and WORK(IR) is M by M
                            //
                            ldwrku = m;
                            ir = iu + ldwrku * m;
                            ldwrkr = m;
                        }
                        itau = ir + ldwrkr * m;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(n, n, m, vt, ldvt, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to WORK(IU), zeroing out above it
                        //
                        Clacpy("L", m, m, a, lda, &work[iu - 1], ldwrku);
                        Claset("U", m - 1, m - 1, czero, czero, &work[(iu + ldwrku) - 1], ldwrku);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IU), copying result to
                        //                    WORK(IR)
                        //                    (CWorkspace: need   2*M*M+3*M,
                        //                                 prefer 2*M*M+2*M+2*M*NB)
                        //                    (RWorkspace: need   M)
                        //
                        Cgebrd(m, m, &work[iu - 1], ldwrku, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, m, &work[iu - 1], ldwrku, &work[ir - 1], ldwrkr);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IU)
                        //                    (CWorkspace: need   2*M*M+3*M-1,
                        //                                 prefer 2*M*M+2*M+(M-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", m, m, m, &work[iu - 1], ldwrku, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IR)
                        //                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", m, m, m, &work[ir - 1], ldwrkr, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of L in WORK(IR) and computing
                        //                    right singular vectors of L in WORK(IU)
                        //                    (CWorkspace: need 2*M*M)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, m, m, 0, s, &rwork[ie - 1], &work[iu - 1], ldwrku, &work[ir - 1], ldwrkr, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IU) by
                        //                    Q in VT, storing result in A
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, m, cone, &work[iu - 1], ldwrku, vt, ldvt, czero, a, lda);
                        //
                        //                    Copy right singular vectors of A from A to VT
                        //
                        Clacpy("F", m, n, a, lda, vt, ldvt);
                        //
                        //                    Copy left singular vectors of A from WORK(IR) to A
                        //
                        Clacpy("F", m, m, &work[ir - 1], ldwrkr, a, lda);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (CWorkspace: need 2*M, prefer M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (CWorkspace: need M+N, prefer M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(n, n, m, vt, ldvt, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Zero out above L in A
                        //
                        Claset("U", m - 1, m - 1, czero, czero, &a[(2 - 1) * lda], lda);
                        //
                        //                    Bidiagonalize L in A
                        //                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                        //                    (RWorkspace: need M)
                        //
                        Cgebrd(m, m, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right bidiagonalizing vectors in A by Q
                        //                    in VT
                        //                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("P", "L", "C", m, n, m, a, lda, &work[itaup - 1], vt, ldvt, &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in A
                        //                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", m, m, m, a, lda, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in A and computing right
                        //                    singular vectors of A in VT
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, n, m, 0, s, &rwork[ie - 1], vt, ldvt, a, lda, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                } else if (wntuas) {
                    //
                    //                 Path 9t(N much larger than M, JOBU='S' or 'A',
                    //                         JOBVT='A')
                    //                 N right singular vectors to be computed in VT and
                    //                 M left singular vectors to be computed in U
                    //
                    if (lwork >= m * m + max(n + m, 3 * m)) {
                        //
                        //                    Sufficient workspace for a fast algorithm
                        //
                        iu = 1;
                        if (lwork >= wrkbl + lda * m) {
                            //
                            //                       WORK(IU) is LDA by M
                            //
                            ldwrku = lda;
                        } else {
                            //
                            //                       WORK(IU) is M by M
                            //
                            ldwrku = m;
                        }
                        itau = iu + ldwrku * m;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(n, n, m, vt, ldvt, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to WORK(IU), zeroing out above it
                        //
                        Clacpy("L", m, m, a, lda, &work[iu - 1], ldwrku);
                        Claset("U", m - 1, m - 1, czero, czero, &work[(iu + ldwrku) - 1], ldwrku);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IU), copying result to U
                        //                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                        //                    (RWorkspace: need M)
                        //
                        Cgebrd(m, m, &work[iu - 1], ldwrku, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("L", m, m, &work[iu - 1], ldwrku, u, ldu);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IU)
                        //                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("P", m, m, m, &work[iu - 1], ldwrku, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in U
                        //                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", m, m, m, u, ldu, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of L in U and computing right
                        //                    singular vectors of L in WORK(IU)
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, m, m, 0, s, &rwork[ie - 1], &work[iu - 1], ldwrku, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IU) by
                        //                    Q in VT, storing result in A
                        //                    (CWorkspace: need M*M)
                        //                    (RWorkspace: 0)
                        //
                        Cgemm("N", "N", m, n, m, cone, &work[iu - 1], ldwrku, vt, ldvt, czero, a, lda);
                        //
                        //                    Copy right singular vectors of A from A to VT
                        //
                        Clacpy("F", m, n, a, lda, vt, ldvt);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (CWorkspace: need 2*M, prefer M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cgelqf(m, n, a, lda, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        Clacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (CWorkspace: need M+N, prefer M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunglq(n, n, m, vt, ldvt, &work[itau - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to U, zeroing out above it
                        //
                        Clacpy("L", m, m, a, lda, u, ldu);
                        Claset("U", m - 1, m - 1, czero, czero, &u[(2 - 1) * ldu], ldu);
                        ie = 1;
                        itauq = itau;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in U
                        //                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                        //                    (RWorkspace: need M)
                        //
                        Cgebrd(m, m, u, ldu, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right bidiagonalizing vectors in U by Q
                        //                    in VT
                        //                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cunmbr("P", "L", "C", m, n, m, u, ldu, &work[itaup - 1], vt, ldvt, &work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in U
                        //                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
                        //                    (RWorkspace: 0)
                        //
                        Cungbr("Q", m, m, m, u, ldu, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
                        irwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in VT
                        //                    (CWorkspace: 0)
                        //                    (RWorkspace: need BDSPAC)
                        //
                        Cbdsqr("U", m, n, m, 0, s, &rwork[ie - 1], vt, ldvt, u, ldu, cdum, 1, &rwork[irwork - 1], info);
                        //
                    }
                    //
                }
                //
            }
            //
        } else {
            //
            //           N .LT. MNTHR
            //
            //           Path 10t(N greater than M, but not much larger)
            //           Reduce to bidiagonal form without LQ decomposition
            //
            ie = 1;
            itauq = 1;
            itaup = itauq + m;
            iwork = itaup + m;
            //
            //           Bidiagonalize A
            //           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
            //           (RWorkspace: M)
            //
            Cgebrd(m, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
            if (wntuas) {
                //
                //              If left singular vectors desired in U, copy result to U
                //              and generate left bidiagonalizing vectors in U
                //              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB)
                //              (RWorkspace: 0)
                //
                Clacpy("L", m, m, a, lda, u, ldu);
                Cungbr("Q", m, m, n, u, ldu, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntvas) {
                //
                //              If right singular vectors desired in VT, copy result to
                //              VT and generate right bidiagonalizing vectors in VT
                //              (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB)
                //              (RWorkspace: 0)
                //
                Clacpy("U", m, n, a, lda, vt, ldvt);
                if (wntva) {
                    nrvt = n;
                }
                if (wntvs) {
                    nrvt = m;
                }
                Cungbr("P", nrvt, n, m, vt, ldvt, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntuo) {
                //
                //              If left singular vectors desired in A, generate left
                //              bidiagonalizing vectors in A
                //              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB)
                //              (RWorkspace: 0)
                //
                Cungbr("Q", m, m, n, a, lda, &work[itauq - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntvo) {
                //
                //              If right singular vectors desired in A, generate right
                //              bidiagonalizing vectors in A
                //              (CWorkspace: need 3*M, prefer 2*M+M*NB)
                //              (RWorkspace: 0)
                //
                Cungbr("P", m, n, m, a, lda, &work[itaup - 1], &work[iwork - 1], lwork - iwork + 1, ierr);
            }
            irwork = ie + m;
            if (wntuas || wntuo) {
                nru = m;
            }
            if (wntun) {
                nru = 0;
            }
            if (wntvas || wntvo) {
                ncvt = n;
            }
            if (wntvn) {
                ncvt = 0;
            }
            if ((!wntuo) && (!wntvo)) {
                //
                //              Perform bidiagonal QR iteration, if desired, computing
                //              left singular vectors in U and computing right singular
                //              vectors in VT
                //              (CWorkspace: 0)
                //              (RWorkspace: need BDSPAC)
                //
                Cbdsqr("L", m, ncvt, nru, 0, s, &rwork[ie - 1], vt, ldvt, u, ldu, cdum, 1, &rwork[irwork - 1], info);
            } else if ((!wntuo) && wntvo) {
                //
                //              Perform bidiagonal QR iteration, if desired, computing
                //              left singular vectors in U and computing right singular
                //              vectors in A
                //              (CWorkspace: 0)
                //              (RWorkspace: need BDSPAC)
                //
                Cbdsqr("L", m, ncvt, nru, 0, s, &rwork[ie - 1], a, lda, u, ldu, cdum, 1, &rwork[irwork - 1], info);
            } else {
                //
                //              Perform bidiagonal QR iteration, if desired, computing
                //              left singular vectors in A and computing right singular
                //              vectors in VT
                //              (CWorkspace: 0)
                //              (RWorkspace: need BDSPAC)
                //
                Cbdsqr("L", m, ncvt, nru, 0, s, &rwork[ie - 1], vt, ldvt, a, lda, cdum, 1, &rwork[irwork - 1], info);
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
    //     End of Cgesvd
    //
}
