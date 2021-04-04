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

void Rgelss(INTEGER const &m, INTEGER const &n, INTEGER const &nrhs, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL *s, REAL const &rcond, INTEGER &rank, REAL *work, INTEGER const &lwork, INTEGER &info) {
    INTEGER minmn = 0;
    INTEGER maxmn = 0;
    bool lquery = false;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER mm = 0;
    INTEGER mnthr = 0;
    arr_1d<1, REAL> dum(fill0);
    INTEGER lwork_Rgeqrf = 0;
    INTEGER lwork_Rormqr = 0;
    INTEGER bdspac = 0;
    INTEGER lwork_Rgebrd = 0;
    INTEGER lwork_Rormbr = 0;
    INTEGER lwork_Rorgbr = 0;
    INTEGER lwork_Rgelqf = 0;
    INTEGER lwork_Rormlq = 0;
    REAL eps = 0.0;
    REAL sfmin = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL anrm = 0.0;
    INTEGER iascl = 0;
    const REAL zero = 0.0;
    REAL bnrm = 0.0;
    INTEGER ibscl = 0;
    INTEGER itau = 0;
    INTEGER iwork = 0;
    INTEGER ie = 0;
    INTEGER itauq = 0;
    INTEGER itaup = 0;
    REAL thr = 0.0;
    INTEGER i = 0;
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
    //       NB refers to the optimal block size for the immediately
    //       following subroutine, as returned by iMlaenv.)
    //
    if (info == 0) {
        minwrk = 1;
        maxwrk = 1;
        if (minmn > 0) {
            mm = m;
            mnthr = iMlaenv[(6 - 1) + ("Rgelss" - 1) * ldiMlaenv];
            if (m >= n && m >= mnthr) {
                //
                //              Path 1a - overdetermined, with many more rows than
                //                        columns
                //
                //              Compute space needed for Rgeqrf
                Rgeqrf(m, n, a, lda, dum[1 - 1], dum[1 - 1], -1, info);
                lwork_Rgeqrf = dum[1 - 1];
                //              Compute space needed for Rormqr
                Rormqr("L", "T", m, nrhs, n, a, lda, dum[1 - 1], b, ldb, dum[1 - 1], -1, info);
                lwork_Rormqr = dum[1 - 1];
                mm = n;
                maxwrk = max(maxwrk, n + lwork_Rgeqrf);
                maxwrk = max(maxwrk, n + lwork_Rormqr);
            }
            if (m >= n) {
                //
                //              Path 1 - overdetermined or exactly determined
                //
                //              Compute workspace needed for Rbdsqr
                //
                bdspac = max((INTEGER)1, 5 * n);
                //              Compute space needed for Rgebrd
                Rgebrd(mm, n, a, lda, s, dum[1 - 1], dum[1 - 1], dum[1 - 1], dum[1 - 1], -1, info);
                lwork_Rgebrd = dum[1 - 1];
                //              Compute space needed for Rormbr
                Rormbr("Q", "L", "T", mm, nrhs, n, a, lda, dum[1 - 1], b, ldb, dum[1 - 1], -1, info);
                lwork_Rormbr = dum[1 - 1];
                //              Compute space needed for Rorgbr
                Rorgbr("P", n, n, n, a, lda, dum[1 - 1], dum[1 - 1], -1, info);
                lwork_Rorgbr = dum[1 - 1];
                //              Compute total workspace needed
                maxwrk = max(maxwrk, 3 * n + lwork_Rgebrd);
                maxwrk = max(maxwrk, 3 * n + lwork_Rormbr);
                maxwrk = max(maxwrk, 3 * n + lwork_Rorgbr);
                maxwrk = max(maxwrk, bdspac);
                maxwrk = max(maxwrk, n * nrhs);
                minwrk = max(3 * n + mm, 3 * n + nrhs, bdspac);
                maxwrk = max(minwrk, maxwrk);
            }
            if (n > m) {
                //
                //              Compute workspace needed for Rbdsqr
                //
                bdspac = max((INTEGER)1, 5 * m);
                minwrk = max(3 * m + nrhs, 3 * m + n, bdspac);
                if (n >= mnthr) {
                    //
                    //                 Path 2a - underdetermined, with many more columns
                    //                 than rows
                    //
                    //                 Compute space needed for Rgelqf
                    Rgelqf(m, n, a, lda, dum[1 - 1], dum[1 - 1], -1, info);
                    lwork_Rgelqf = dum[1 - 1];
                    //                 Compute space needed for Rgebrd
                    Rgebrd(m, m, a, lda, s, dum[1 - 1], dum[1 - 1], dum[1 - 1], dum[1 - 1], -1, info);
                    lwork_Rgebrd = dum[1 - 1];
                    //                 Compute space needed for Rormbr
                    Rormbr("Q", "L", "T", m, nrhs, n, a, lda, dum[1 - 1], b, ldb, dum[1 - 1], -1, info);
                    lwork_Rormbr = dum[1 - 1];
                    //                 Compute space needed for Rorgbr
                    Rorgbr("P", m, m, m, a, lda, dum[1 - 1], dum[1 - 1], -1, info);
                    lwork_Rorgbr = dum[1 - 1];
                    //                 Compute space needed for Rormlq
                    Rormlq("L", "T", n, nrhs, m, a, lda, dum[1 - 1], b, ldb, dum[1 - 1], -1, info);
                    lwork_Rormlq = dum[1 - 1];
                    //                 Compute total workspace needed
                    maxwrk = m + lwork_Rgelqf;
                    maxwrk = max(maxwrk, m * m + 4 * m + lwork_Rgebrd);
                    maxwrk = max(maxwrk, m * m + 4 * m + lwork_Rormbr);
                    maxwrk = max(maxwrk, m * m + 4 * m + lwork_Rorgbr);
                    maxwrk = max(maxwrk, m * m + m + bdspac);
                    if (nrhs > 1) {
                        maxwrk = max(maxwrk, m * m + m + m * nrhs);
                    } else {
                        maxwrk = max(maxwrk, m * m + 2 * m);
                    }
                    maxwrk = max(maxwrk, m + lwork_Rormlq);
                } else {
                    //
                    //                 Path 2 - underdetermined
                    //
                    //                 Compute space needed for Rgebrd
                    Rgebrd(m, n, a, lda, s, dum[1 - 1], dum[1 - 1], dum[1 - 1], dum[1 - 1], -1, info);
                    lwork_Rgebrd = dum[1 - 1];
                    //                 Compute space needed for Rormbr
                    Rormbr("Q", "L", "T", m, nrhs, m, a, lda, dum[1 - 1], b, ldb, dum[1 - 1], -1, info);
                    lwork_Rormbr = dum[1 - 1];
                    //                 Compute space needed for Rorgbr
                    Rorgbr("P", m, n, m, a, lda, dum[1 - 1], dum[1 - 1], -1, info);
                    lwork_Rorgbr = dum[1 - 1];
                    maxwrk = 3 * m + lwork_Rgebrd;
                    maxwrk = max(maxwrk, 3 * m + lwork_Rormbr);
                    maxwrk = max(maxwrk, 3 * m + lwork_Rorgbr);
                    maxwrk = max(maxwrk, bdspac);
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
        Mxerbla("Rgelss", -info);
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
    anrm = Rlange[("M" - 1) + (m - 1) * ldRlange];
    iascl = 0;
    if (anrm > zero && anrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM
        //
        Rlascl("G", 0, 0, anrm, smlnum, m, n, a, lda, info);
        iascl = 1;
    } else if (anrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM
        //
        Rlascl("G", 0, 0, anrm, bignum, m, n, a, lda, info);
        iascl = 2;
    } else if (anrm == zero) {
        //
        //        Matrix all zero. Return zero solution.
        //
        Rlaset("F", max(m, n), nrhs, zero, zero, b, ldb);
        Rlaset("F", minmn, 1, zero, zero, s, minmn);
        rank = 0;
        goto statement_70;
    }
    //
    //     Scale B if max element outside range [SMLNUM,BIGNUM]
    //
    bnrm = Rlange[("M" - 1) + (m - 1) * ldRlange];
    ibscl = 0;
    if (bnrm > zero && bnrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM
        //
        Rlascl("G", 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info);
        ibscl = 1;
    } else if (bnrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM
        //
        Rlascl("G", 0, 0, bnrm, bignum, m, nrhs, b, ldb, info);
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
            //           (Workspace: need 2*N, prefer N+N*NB)
            //
            Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, info);
            //
            //           Multiply B by transpose(Q)
            //           (Workspace: need N+NRHS, prefer N+NRHS*NB)
            //
            Rormqr("L", "T", m, nrhs, n, a, lda, work[itau - 1], b, ldb, work[iwork - 1], lwork - iwork + 1, info);
            //
            //           Zero out below R
            //
            if (n > 1) {
                Rlaset("L", n - 1, n - 1, zero, zero, a[(2 - 1)], lda);
            }
        }
        //
        ie = 1;
        itauq = ie + n;
        itaup = itauq + n;
        iwork = itaup + n;
        //
        //        Bidiagonalize R in A
        //        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
        //
        Rgebrd(mm, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors of R
        //        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
        //
        Rormbr("Q", "L", "T", mm, nrhs, n, a, lda, work[itauq - 1], b, ldb, work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Generate right bidiagonalizing vectors of R in A
        //        (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
        //
        Rorgbr("P", n, n, n, a, lda, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        iwork = ie + n;
        //
        //        Perform bidiagonal QR iteration
        //          multiply B by transpose of left singular vectors
        //          compute right singular vectors in A
        //        (Workspace: need BDSPAC)
        //
        Rbdsqr("U", n, n, 0, nrhs, s, work[ie - 1], a, lda, dum, 1, b, ldb, work[iwork - 1], info);
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
                Rrscl(nrhs, s[i - 1], b[(i - 1)], ldb);
                rank++;
            } else {
                Rlaset("F", 1, nrhs, zero, zero, b[(i - 1)], ldb);
            }
        }
        //
        //        Multiply B by right singular vectors
        //        (Workspace: need N, prefer N*NRHS)
        //
        if (lwork >= ldb * nrhs && nrhs > 1) {
            Rgemm("T", "N", n, nrhs, n, one, a, lda, b, ldb, zero, work, ldb);
            Rlacpy("G", n, nrhs, work, ldb, b, ldb);
        } else if (nrhs > 1) {
            chunk = lwork / n;
            for (i = 1; i <= nrhs; i = i + chunk) {
                bl = min(nrhs - i + 1, chunk);
                Rgemm("T", "N", n, bl, n, one, a, lda, b[(i - 1) * ldb], ldb, zero, work, n);
                Rlacpy("G", n, bl, work, n, b[(i - 1) * ldb], ldb);
            }
        } else {
            Rgemv("T", n, n, one, a, lda, b, 1, zero, work, 1);
            Rcopy(n, work, 1, b, 1);
        }
        //
    } else if (n >= mnthr && lwork >= 4 * m + m * m + max(m, 2 * m - 4, nrhs, n - 3 * m)) {
        //
        //        Path 2a - underdetermined, with many more columns than rows
        //        and sufficient workspace for an efficient algorithm
        //
        ldwork = m;
        if (lwork >= max(4 * m + m * lda + max(m, 2 * m - 4, nrhs, n - 3 * m), m * lda + m + m * nrhs)) {
            ldwork = lda;
        }
        itau = 1;
        iwork = m + 1;
        //
        //        Compute A=L*Q
        //        (Workspace: need 2*M, prefer M+M*NB)
        //
        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, info);
        il = iwork;
        //
        //        Copy L to WORK(IL), zeroing out above it
        //
        Rlacpy("L", m, m, a, lda, work[il - 1], ldwork);
        Rlaset("U", m - 1, m - 1, zero, zero, work[(il + ldwork) - 1], ldwork);
        ie = il + ldwork * m;
        itauq = ie + m;
        itaup = itauq + m;
        iwork = itaup + m;
        //
        //        Bidiagonalize L in WORK(IL)
        //        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
        //
        Rgebrd(m, m, work[il - 1], ldwork, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors of L
        //        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
        //
        Rormbr("Q", "L", "T", m, nrhs, m, work[il - 1], ldwork, work[itauq - 1], b, ldb, work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Generate right bidiagonalizing vectors of R in WORK(IL)
        //        (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB)
        //
        Rorgbr("P", m, m, m, work[il - 1], ldwork, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        iwork = ie + m;
        //
        //        Perform bidiagonal QR iteration,
        //           computing right singular vectors of L in WORK(IL) and
        //           multiplying B by transpose of left singular vectors
        //        (Workspace: need M*M+M+BDSPAC)
        //
        Rbdsqr("U", m, m, 0, nrhs, s, work[ie - 1], work[il - 1], ldwork, a, lda, b, ldb, work[iwork - 1], info);
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
                Rrscl(nrhs, s[i - 1], b[(i - 1)], ldb);
                rank++;
            } else {
                Rlaset("F", 1, nrhs, zero, zero, b[(i - 1)], ldb);
            }
        }
        iwork = ie;
        //
        //        Multiply B by right singular vectors of L in WORK(IL)
        //        (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS)
        //
        if (lwork >= ldb * nrhs + iwork - 1 && nrhs > 1) {
            Rgemm("T", "N", m, nrhs, m, one, work[il - 1], ldwork, b, ldb, zero, work[iwork - 1], ldb);
            Rlacpy("G", m, nrhs, work[iwork - 1], ldb, b, ldb);
        } else if (nrhs > 1) {
            chunk = (lwork - iwork + 1) / m;
            for (i = 1; i <= nrhs; i = i + chunk) {
                bl = min(nrhs - i + 1, chunk);
                Rgemm("T", "N", m, bl, m, one, work[il - 1], ldwork, b[(i - 1) * ldb], ldb, zero, work[iwork - 1], m);
                Rlacpy("G", m, bl, work[iwork - 1], m, b[(i - 1) * ldb], ldb);
            }
        } else {
            Rgemv("T", m, m, one, work[il - 1], ldwork, b[(1 - 1)], 1, zero, work[iwork - 1], 1);
            Rcopy(m, work[iwork - 1], 1, b[(1 - 1)], 1);
        }
        //
        //        Zero out below first M rows of B
        //
        Rlaset("F", n - m, nrhs, zero, zero, b[((m + 1) - 1)], ldb);
        iwork = itau + m;
        //
        //        Multiply transpose(Q) by B
        //        (Workspace: need M+NRHS, prefer M+NRHS*NB)
        //
        Rormlq("L", "T", n, nrhs, m, a, lda, work[itau - 1], b, ldb, work[iwork - 1], lwork - iwork + 1, info);
        //
    } else {
        //
        //        Path 2 - remaining underdetermined cases
        //
        ie = 1;
        itauq = ie + m;
        itaup = itauq + m;
        iwork = itaup + m;
        //
        //        Bidiagonalize A
        //        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
        //
        Rgebrd(m, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors
        //        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
        //
        Rormbr("Q", "L", "T", m, nrhs, n, a, lda, work[itauq - 1], b, ldb, work[iwork - 1], lwork - iwork + 1, info);
        //
        //        Generate right bidiagonalizing vectors in A
        //        (Workspace: need 4*M, prefer 3*M+M*NB)
        //
        Rorgbr("P", m, n, m, a, lda, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, info);
        iwork = ie + m;
        //
        //        Perform bidiagonal QR iteration,
        //           computing right singular vectors of A in A and
        //           multiplying B by transpose of left singular vectors
        //        (Workspace: need BDSPAC)
        //
        Rbdsqr("L", m, n, 0, nrhs, s, work[ie - 1], a, lda, dum, 1, b, ldb, work[iwork - 1], info);
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
                Rrscl(nrhs, s[i - 1], b[(i - 1)], ldb);
                rank++;
            } else {
                Rlaset("F", 1, nrhs, zero, zero, b[(i - 1)], ldb);
            }
        }
        //
        //        Multiply B by right singular vectors of A
        //        (Workspace: need N, prefer N*NRHS)
        //
        if (lwork >= ldb * nrhs && nrhs > 1) {
            Rgemm("T", "N", n, nrhs, m, one, a, lda, b, ldb, zero, work, ldb);
            Rlacpy("F", n, nrhs, work, ldb, b, ldb);
        } else if (nrhs > 1) {
            chunk = lwork / n;
            for (i = 1; i <= nrhs; i = i + chunk) {
                bl = min(nrhs - i + 1, chunk);
                Rgemm("T", "N", n, bl, m, one, a, lda, b[(i - 1) * ldb], ldb, zero, work, n);
                Rlacpy("F", n, bl, work, n, b[(i - 1) * ldb], ldb);
            }
        } else {
            Rgemv("T", m, n, one, a, lda, b, 1, zero, work, 1);
            Rcopy(n, work, 1, b, 1);
        }
    }
    //
    //     Undo scaling
    //
    if (iascl == 1) {
        Rlascl("G", 0, 0, anrm, smlnum, n, nrhs, b, ldb, info);
        Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, s, minmn, info);
    } else if (iascl == 2) {
        Rlascl("G", 0, 0, anrm, bignum, n, nrhs, b, ldb, info);
        Rlascl("G", 0, 0, bignum, anrm, minmn, 1, s, minmn, info);
    }
    if (ibscl == 1) {
        Rlascl("G", 0, 0, smlnum, bnrm, n, nrhs, b, ldb, info);
    } else if (ibscl == 2) {
        Rlascl("G", 0, 0, bignum, bnrm, n, nrhs, b, ldb, info);
    }
//
statement_70:
    work[1 - 1] = maxwrk;
    //
    //     End of Rgelss
    //
}
