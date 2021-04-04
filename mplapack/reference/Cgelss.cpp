/*
 * Copyright (c) 2021
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

void Cgelss(INTEGER const &m, INTEGER const &n, INTEGER const &nrhs, COMPLEX *a, INTEGER const &lda, COMPLEX *b, INTEGER const &ldb, REAL *s, REAL const &rcond, INTEGER &rank, COMPLEX *work, INTEGER const &lwork, REAL *rwork, INTEGER &info) {
    INTEGER minmn = 0;
    INTEGER maxmn = 0;
    bool lquery = false;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER mm = 0;
    INTEGER mnthr = 0;
    arr_1d<1, COMPLEX> dum(fill0);
    INTEGER lwork_Cgeqrf = 0;
    INTEGER lwork_Cunmqr = 0;
    INTEGER lwork_Cgebrd = 0;
    INTEGER lwork_Cunmbr = 0;
    INTEGER lwork_Cungbr = 0;
    INTEGER lwork_Cgelqf = 0;
    INTEGER lwork_Cunmlq = 0;
    REAL eps = 0.0;
    REAL sfmin = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL anrm = 0.0;
    INTEGER iascl = 0;
    const REAL zero = 0.0;
    const COMPLEX czero = (0.0, 0.0);
    REAL bnrm = 0.0;
    INTEGER ibscl = 0;
    INTEGER itau = 0;
    INTEGER iwork = 0;
    INTEGER ie = 0;
    INTEGER itauq = 0;
    INTEGER itaup = 0;
    INTEGER irwork = 0;
    REAL thr = 0.0;
    INTEGER i = 0;
    const COMPLEX cone = (1.0, 0.0);
    INTEGER chunk = 0;
    INTEGER bl = 0;
    INTEGER ldwork = 0;
    INTEGER il = 0;
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
    minmn = min(m, n);
    maxmn = max(m, n);
    lquery = (lwork == -1);
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (nrhs < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, m)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, maxmn)) {
        info = -7;
    }
    //
    //     Compute workspace
    //      (Note: Comments in the code beginning "Workspace:" describe the
    //       minimal amount of workspace needed at that poINTEGER in the code,
    //       as well as the preferred amount for good performance.
    //       CWorkspace refers to complex workspace, and RWorkspace refers
    //       to real workspace. NB refers to the optimal block size for the
    //       immediately following subroutine, as returned by iMlaenv.)
    //
    if (info == 0) {
        minwrk = 1;
        maxwrk = 1;
        if (minmn > 0) {
            mm = m;
            mnthr = iMlaenv[(6 - 1) + ("Cgelss" - 1) * ldiMlaenv];
            if (m >= n && m >= mnthr) {
                //
                //              Path 1a - overdetermined, with many more rows than
                //                        columns
                //
                //              Compute space needed for Cgeqrf
                Cgeqrf(m, n, a, lda, dum[1 - 1], dum[1 - 1], -1, info);
                lwork_Cgeqrf = dum[1 - 1];
                //              Compute space needed for Cunmqr
                Cunmqr("L", "C", m, nrhs, n, a, lda, dum[1 - 1], b, ldb, dum[1 - 1], -1, info);
                lwork_Cunmqr = dum[1 - 1];
                mm = n;
                maxwrk = max(maxwrk, n + n * iMlaenv[("Cgeqrf" - 1) * ldiMlaenv]);
                maxwrk = max(maxwrk, n + nrhs * iMlaenv[("Cunmqr" - 1) * ldiMlaenv]);
            }
            if (m >= n) {
                //
                //              Path 1 - overdetermined or exactly determined
                //
                //              Compute space needed for Cgebrd
                Cgebrd(mm, n, a, lda, s, s, dum[1 - 1], dum[1 - 1], dum[1 - 1], -1, info);
                lwork_Cgebrd = dum[1 - 1];
                //              Compute space needed for Cunmbr
                Cunmbr("Q", "L", "C", mm, nrhs, n, a, lda, dum[1 - 1], b, ldb, dum[1 - 1], -1, info);
                lwork_Cunmbr = dum[1 - 1];
                //              Compute space needed for Cungbr
                Cungbr("P", n, n, n, a, lda, dum[1 - 1], dum[1 - 1], -1, info);
                lwork_Cungbr = dum[1 - 1];
                //              Compute total workspace needed
                maxwrk = max(maxwrk, 2 * n + lwork_Cgebrd);
                maxwrk = max(maxwrk, 2 * n + lwork_Cunmbr);
                maxwrk = max(maxwrk, 2 * n + lwork_Cungbr);
                maxwrk = max(maxwrk, n * nrhs);
                minwrk = 2 * n + max(nrhs, m);
            }
            if (n > m) {
                minwrk = 2 * m + max(nrhs, n);
                if (n >= mnthr) {
                    //
                    //                 Path 2a - underdetermined, with many more columns
                    //                 than rows
                    //
                    //                 Compute space needed for Cgelqf
                    Cgelqf(m, n, a, lda, dum[1 - 1], dum[1 - 1], -1, info);
                    lwork_Cgelqf = dum[1 - 1];
                    //                 Compute space needed for Cgebrd
                    Cgebrd(m, m, a, lda, s, s, dum[1 - 1], dum[1 - 1], dum[1 - 1], -1, info);
                    lwork_Cgebrd = dum[1 - 1];
                    //                 Compute space needed for Cunmbr
                    Cunmbr("Q", "L", "C", m, nrhs, n, a, lda, dum[1 - 1], b, ldb, dum[1 - 1], -1, info);
                    lwork_Cunmbr = dum[1 - 1];
                    //                 Compute space needed for Cungbr
                    Cungbr("P", m, m, m, a, lda, dum[1 - 1], dum[1 - 1], -1, info);
                    lwork_Cungbr = dum[1 - 1];
                    //                 Compute space needed for Cunmlq
                    Cunmlq("L", "C", n, nrhs, m, a, lda, dum[1 - 1], b, ldb, dum[1 - 1], -1, info);
                    lwork_Cunmlq = dum[1 - 1];
                    //                 Compute total workspace needed
                    maxwrk = m + lwork_Cgelqf;
                    maxwrk = max(maxwrk, 3 * m + m * m + lwork_Cgebrd);
                    maxwrk = max(maxwrk, 3 * m + m * m + lwork_Cunmbr);
                    maxwrk = max(maxwrk, 3 * m + m * m + lwork_Cungbr);
                    if (nrhs > 1) {
                        maxwrk = max(maxwrk, m * m + m + m * nrhs);
                    } else {
                        maxwrk = max(maxwrk, m * m + 2 * m);
                    }
                    maxwrk = max(maxwrk, m + lwork_Cunmlq);
                } else {
                    //
                    //                 Path 2 - underdetermined
                    //
                    //                 Compute space needed for Cgebrd
                    Cgebrd(m, n, a, lda, s, s, dum[1 - 1], dum[1 - 1], dum[1 - 1], -1, info);
                    lwork_Cgebrd = dum[1 - 1];
                    //                 Compute space needed for Cunmbr
                    Cunmbr("Q", "L", "C", m, nrhs, m, a, lda, dum[1 - 1], b, ldb, dum[1 - 1], -1, info);
                    lwork_Cunmbr = dum[1 - 1];
                    //                 Compute space needed for Cungbr
                    Cungbr("P", m, n, m, a, lda, dum[1 - 1], dum[1 - 1], -1, info);
                    lwork_Cungbr = dum[1 - 1];
                    maxwrk = 2 * m + lwork_Cgebrd;
                    maxwrk = max(maxwrk, 2 * m + lwork_Cunmbr);
                    maxwrk = max(maxwrk, 2 * m + lwork_Cungbr);
                    maxwrk = max(maxwrk, n * nrhs);
                }
            }
            maxwrk = max(minwrk, maxwrk);
        }
        work[1 - 1] = maxwrk;
        //
        if (lwork < minwrk && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cgelss", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        rank = 0;
        return;
    }
    //
    //     Get machine parameters
    //
    eps = dlamch("P");
    sfmin = dlamch("S");
    smlnum = sfmin / eps;
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Scale A if max element outside range [SMLNUM,BIGNUM]
    //
    anrm = Clange[("M" - 1) + (m - 1) * ldClange];
    iascl = 0;
    if (anrm > zero && anrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM
        //
        Clascl("G", 0, 0, anrm, smlnum, m, n, a, lda, info);
        iascl = 1;
    } else if (anrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM
        //
        Clascl("G", 0, 0, anrm, bignum, m, n, a, lda, info);
        iascl = 2;
    } else if (anrm == zero) {
        //
        //        Matrix all zero. Return zero solution.
        //
        Claset("F", max(m, n), nrhs, czero, czero, b, ldb);
        Rlaset("F", minmn, 1, zero, zero, s, minmn);
        rank = 0;
        goto statement_70;
    }
    //
    //     Scale B if max element outside range [SMLNUM,BIGNUM]
    //
    bnrm = Clange[("M" - 1) + (m - 1) * ldClange];
    ibscl = 0;
    if (bnrm > zero && bnrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM
        //
        Clascl("G", 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info);
        ibscl = 1;
    } else if (bnrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM
        //
        Clascl("G", 0, 0, bnrm, bignum, m, nrhs, b, ldb, info);
        ibscl = 2;
    }
    //
    //     Overdetermined case
    //
    if (m >= n) {
        //
        //        Path 1 - overdetermined or exactly determined
        //
        mm = m;
        if (m >= mnthr) {
            //
            //           Path 1a - overdetermined, with many more rows than columns
            //
            mm = n;
            itau = 1;
            iwork = itau + n;
            //
            //           Compute A=Q*R
            //           (CWorkspace: need 2*N, prefer N+N*NB)
            //           (RWorkspace: none)
            //
            Cgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, info);
            //
            //           Multiply B by transpose(Q)
            //           (CWorkspace: need N+NRHS, prefer N+NRHS*NB)
            //           (RWorkspace: none)
            //
            Cunmqr("L", "C", m, nrhs, n, a, lda, work[itau - 1], b, ldb, work[iwork - 1], lwork - iwork + 1, info);
            //
            //           Zero out below R
            //
            if (n > 1) {
                Claset("L", n - 1, n - 1, czero, czero, a[(2 - 1)], lda);
            }
        }
        //
        ie = 1;
        itauq = 1;
        itaup = itauq + n;
        iwork = itaup + n;
        //
        //        Bidiagonalize R in A
        //        (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
        //        (RWorkspace: need N)
        //
        Cgebrd(mm, n, a, lda, s, rwork[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors of R
        //        (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
        //        (RWorkspace: none)
        //
        Cunmbr("Q", "L", "C", mm, nrhs, n, a, lda, work[itauq - 1], b, ldb, work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Generate right bidiagonalizing vectors of R in A
        //        (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
        //        (RWorkspace: none)
        //
        Cungbr("P", n, n, n, a, lda, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        irwork = ie + n;
        //
        //        Perform bidiagonal QR iteration
        //          multiply B by transpose of left singular vectors
        //          compute right singular vectors in A
        //        (CWorkspace: none)
        //        (RWorkspace: need BDSPAC)
        //
        Cbdsqr("U", n, n, 0, nrhs, s, rwork[ie - 1], a, lda, dum, 1, b, ldb, rwork[irwork - 1], info);
        if (info != 0) {
            goto statement_70;
        }
        //
        //        Multiply B by reciprocals of singular values
        //
        thr = max(rcond * s[1 - 1], sfmin);
        if (rcond < zero) {
            thr = max(eps * s[1 - 1], sfmin);
        }
        rank = 0;
        for (i = 1; i <= n; i = i + 1) {
            if (s[i - 1] > thr) {
                CRrscl(nrhs, s[i - 1], b[(i - 1)], ldb);
                rank++;
            } else {
                Claset("F", 1, nrhs, czero, czero, b[(i - 1)], ldb);
            }
        }
        //
        //        Multiply B by right singular vectors
        //        (CWorkspace: need N, prefer N*NRHS)
        //        (RWorkspace: none)
        //
        if (lwork >= ldb * nrhs && nrhs > 1) {
            Cgemm("C", "N", n, nrhs, n, cone, a, lda, b, ldb, czero, work, ldb);
            Clacpy("G", n, nrhs, work, ldb, b, ldb);
        } else if (nrhs > 1) {
            chunk = lwork / n;
            for (i = 1; i <= nrhs; i = i + chunk) {
                bl = min(nrhs - i + 1, chunk);
                Cgemm("C", "N", n, bl, n, cone, a, lda, b[(i - 1) * ldb], ldb, czero, work, n);
                Clacpy("G", n, bl, work, n, b[(i - 1) * ldb], ldb);
            }
        } else {
            Cgemv("C", n, n, cone, a, lda, b, 1, czero, work, 1);
            Ccopy(n, work, 1, b, 1);
        }
        //
    } else if (n >= mnthr && lwork >= 3 * m + m * m + max(m, nrhs, n - 2 * m)) {
        //
        //        Underdetermined case, M much less than N
        //
        //        Path 2a - underdetermined, with many more columns than rows
        //        and sufficient workspace for an efficient algorithm
        //
        ldwork = m;
        if (lwork >= 3 * m + m * lda + max(m, nrhs, n - 2 * m)) {
            ldwork = lda;
        }
        itau = 1;
        iwork = m + 1;
        //
        //        Compute A=L*Q
        //        (CWorkspace: need 2*M, prefer M+M*NB)
        //        (RWorkspace: none)
        //
        Cgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, info);
        il = iwork;
        //
        //        Copy L to WORK(IL), zeroing out above it
        //
        Clacpy("L", m, m, a, lda, work[il - 1], ldwork);
        Claset("U", m - 1, m - 1, czero, czero, work[(il + ldwork) - 1], ldwork);
        ie = 1;
        itauq = il + ldwork * m;
        itaup = itauq + m;
        iwork = itaup + m;
        //
        //        Bidiagonalize L in WORK(IL)
        //        (CWorkspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
        //        (RWorkspace: need M)
        //
        Cgebrd(m, m, work[il - 1], ldwork, s, rwork[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors of L
        //        (CWorkspace: need M*M+3*M+NRHS, prefer M*M+3*M+NRHS*NB)
        //        (RWorkspace: none)
        //
        Cunmbr("Q", "L", "C", m, nrhs, m, work[il - 1], ldwork, work[itauq - 1], b, ldb, work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Generate right bidiagonalizing vectors of R in WORK(IL)
        //        (CWorkspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)
        //        (RWorkspace: none)
        //
        Cungbr("P", m, m, m, work[il - 1], ldwork, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        irwork = ie + m;
        //
        //        Perform bidiagonal QR iteration, computing right singular
        //        vectors of L in WORK(IL) and multiplying B by transpose of
        //        left singular vectors
        //        (CWorkspace: need M*M)
        //        (RWorkspace: need BDSPAC)
        //
        Cbdsqr("U", m, m, 0, nrhs, s, rwork[ie - 1], work[il - 1], ldwork, a, lda, b, ldb, rwork[irwork - 1], info);
        if (info != 0) {
            goto statement_70;
        }
        //
        //        Multiply B by reciprocals of singular values
        //
        thr = max(rcond * s[1 - 1], sfmin);
        if (rcond < zero) {
            thr = max(eps * s[1 - 1], sfmin);
        }
        rank = 0;
        for (i = 1; i <= m; i = i + 1) {
            if (s[i - 1] > thr) {
                CRrscl(nrhs, s[i - 1], b[(i - 1)], ldb);
                rank++;
            } else {
                Claset("F", 1, nrhs, czero, czero, b[(i - 1)], ldb);
            }
        }
        iwork = il + m * ldwork;
        //
        //        Multiply B by right singular vectors of L in WORK(IL)
        //        (CWorkspace: need M*M+2*M, prefer M*M+M+M*NRHS)
        //        (RWorkspace: none)
        //
        if (lwork >= ldb * nrhs + iwork - 1 && nrhs > 1) {
            Cgemm("C", "N", m, nrhs, m, cone, work[il - 1], ldwork, b, ldb, czero, work[iwork - 1], ldb);
            Clacpy("G", m, nrhs, work[iwork - 1], ldb, b, ldb);
        } else if (nrhs > 1) {
            chunk = (lwork - iwork + 1) / m;
            for (i = 1; i <= nrhs; i = i + chunk) {
                bl = min(nrhs - i + 1, chunk);
                Cgemm("C", "N", m, bl, m, cone, work[il - 1], ldwork, b[(i - 1) * ldb], ldb, czero, work[iwork - 1], m);
                Clacpy("G", m, bl, work[iwork - 1], m, b[(i - 1) * ldb], ldb);
            }
        } else {
            Cgemv("C", m, m, cone, work[il - 1], ldwork, b[(1 - 1)], 1, czero, work[iwork - 1], 1);
            Ccopy(m, work[iwork - 1], 1, b[(1 - 1)], 1);
        }
        //
        //        Zero out below first M rows of B
        //
        Claset("F", n - m, nrhs, czero, czero, b[((m + 1) - 1)], ldb);
        iwork = itau + m;
        //
        //        Multiply transpose(Q) by B
        //        (CWorkspace: need M+NRHS, prefer M+NHRS*NB)
        //        (RWorkspace: none)
        //
        Cunmlq("L", "C", n, nrhs, m, a, lda, work[itau - 1], b, ldb, work[iwork - 1], lwork - iwork + 1, info);
        //
    } else {
        //
        //        Path 2 - remaining underdetermined cases
        //
        ie = 1;
        itauq = 1;
        itaup = itauq + m;
        iwork = itaup + m;
        //
        //        Bidiagonalize A
        //        (CWorkspace: need 3*M, prefer 2*M+(M+N)*NB)
        //        (RWorkspace: need N)
        //
        Cgebrd(m, n, a, lda, s, rwork[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors
        //        (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
        //        (RWorkspace: none)
        //
        Cunmbr("Q", "L", "C", m, nrhs, n, a, lda, work[itauq - 1], b, ldb, work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Generate right bidiagonalizing vectors in A
        //        (CWorkspace: need 3*M, prefer 2*M+M*NB)
        //        (RWorkspace: none)
        //
        Cungbr("P", m, n, m, a, lda, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        irwork = ie + m;
        //
        //        Perform bidiagonal QR iteration,
        //           computing right singular vectors of A in A and
        //           multiplying B by transpose of left singular vectors
        //        (CWorkspace: none)
        //        (RWorkspace: need BDSPAC)
        //
        Cbdsqr("L", m, n, 0, nrhs, s, rwork[ie - 1], a, lda, dum, 1, b, ldb, rwork[irwork - 1], info);
        if (info != 0) {
            goto statement_70;
        }
        //
        //        Multiply B by reciprocals of singular values
        //
        thr = max(rcond * s[1 - 1], sfmin);
        if (rcond < zero) {
            thr = max(eps * s[1 - 1], sfmin);
        }
        rank = 0;
        for (i = 1; i <= m; i = i + 1) {
            if (s[i - 1] > thr) {
                CRrscl(nrhs, s[i - 1], b[(i - 1)], ldb);
                rank++;
            } else {
                Claset("F", 1, nrhs, czero, czero, b[(i - 1)], ldb);
            }
        }
        //
        //        Multiply B by right singular vectors of A
        //        (CWorkspace: need N, prefer N*NRHS)
        //        (RWorkspace: none)
        //
        if (lwork >= ldb * nrhs && nrhs > 1) {
            Cgemm("C", "N", n, nrhs, m, cone, a, lda, b, ldb, czero, work, ldb);
            Clacpy("G", n, nrhs, work, ldb, b, ldb);
        } else if (nrhs > 1) {
            chunk = lwork / n;
            for (i = 1; i <= nrhs; i = i + chunk) {
                bl = min(nrhs - i + 1, chunk);
                Cgemm("C", "N", n, bl, m, cone, a, lda, b[(i - 1) * ldb], ldb, czero, work, n);
                Clacpy("F", n, bl, work, n, b[(i - 1) * ldb], ldb);
            }
        } else {
            Cgemv("C", m, n, cone, a, lda, b, 1, czero, work, 1);
            Ccopy(n, work, 1, b, 1);
        }
    }
    //
    //     Undo scaling
    //
    if (iascl == 1) {
        Clascl("G", 0, 0, anrm, smlnum, n, nrhs, b, ldb, info);
        Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, s, minmn, info);
    } else if (iascl == 2) {
        Clascl("G", 0, 0, anrm, bignum, n, nrhs, b, ldb, info);
        Rlascl("G", 0, 0, bignum, anrm, minmn, 1, s, minmn, info);
    }
    if (ibscl == 1) {
        Clascl("G", 0, 0, smlnum, bnrm, n, nrhs, b, ldb, info);
    } else if (ibscl == 2) {
        Clascl("G", 0, 0, bignum, bnrm, n, nrhs, b, ldb, info);
    }
statement_70:
    work[1 - 1] = maxwrk;
    //
    //     End of Cgelss
    //
}
