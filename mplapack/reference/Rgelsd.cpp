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

void Rgelsd(INTEGER const &m, INTEGER const &n, INTEGER const &nrhs, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL *s, REAL const &rcond, INTEGER &rank, REAL *work, INTEGER const &lwork, arr_ref<INTEGER> iwork, INTEGER &info) {
    INTEGER minmn = 0;
    INTEGER maxmn = 0;
    INTEGER mnthr = 0;
    bool lquery = false;
    INTEGER smlsiz = 0;
    INTEGER minwrk = 0;
    INTEGER liwork = 0;
    const REAL two = 2.0;
    INTEGER nlvl = 0;
    INTEGER maxwrk = 0;
    INTEGER mm = 0;
    INTEGER wlalsd = 0;
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
    INTEGER nwork = 0;
    INTEGER ie = 0;
    INTEGER itauq = 0;
    INTEGER itaup = 0;
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
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments.
    //
    info = 0;
    minmn = min(m, n);
    maxmn = max(m, n);
    mnthr = iMlaenv[(6 - 1) + ("Rgelsd" - 1) * ldiMlaenv];
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
    smlsiz = iMlaenv[(9 - 1) + ("Rgelsd" - 1) * ldiMlaenv];
    //
    //     Compute workspace.
    //     (Note: Comments in the code beginning "Workspace:" describe the
    //     minimal amount of workspace needed at that poINTEGER in the code,
    //     as well as the preferred amount for good performance.
    //     NB refers to the optimal block size for the immediately
    //     following subroutine, as returned by iMlaenv.)
    //
    minwrk = 1;
    liwork = 1;
    minmn = max((INTEGER)1, minmn);
    nlvl = max(INTEGER(log[(minmn.real() / smlsiz + 1.real()) - 1] / log[two - 1]) + 1, 0);
    //
    if (info == 0) {
        maxwrk = 0;
        liwork = 3 * minmn * nlvl + 11 * minmn;
        mm = m;
        if (m >= n && m >= mnthr) {
            //
            //           Path 1a - overdetermined, with many more rows than columns.
            //
            mm = n;
            maxwrk = max(maxwrk, n + n * iMlaenv[("Rgeqrf" - 1) * ldiMlaenv]);
            maxwrk = max(maxwrk, n + nrhs * iMlaenv[("Rormqr" - 1) * ldiMlaenv]);
        }
        if (m >= n) {
            //
            //           Path 1 - overdetermined or exactly determined.
            //
            maxwrk = max(maxwrk, 3 * n + (mm + n) * iMlaenv[("Rgebrd" - 1) * ldiMlaenv]);
            maxwrk = max(maxwrk, 3 * n + nrhs * iMlaenv[("Rormbr" - 1) * ldiMlaenv]);
            maxwrk = max(maxwrk, 3 * n + (n - 1) * iMlaenv[("Rormbr" - 1) * ldiMlaenv]);
            wlalsd = 9 * n + 2 * n * smlsiz + 8 * n * nlvl + n * nrhs + pow2((smlsiz + 1));
            maxwrk = max(maxwrk, 3 * n + wlalsd);
            minwrk = max(3 * n + mm, 3 * n + nrhs, 3 * n + wlalsd);
        }
        if (n > m) {
            wlalsd = 9 * m + 2 * m * smlsiz + 8 * m * nlvl + m * nrhs + pow2((smlsiz + 1));
            if (n >= mnthr) {
                //
                //              Path 2a - underdetermined, with many more columns
                //              than rows.
                //
                maxwrk = m + m * iMlaenv[("Rgelqf" - 1) * ldiMlaenv];
                maxwrk = max(maxwrk, m * m + 4 * m + 2 * m * iMlaenv[("Rgebrd" - 1) * ldiMlaenv]);
                maxwrk = max(maxwrk, m * m + 4 * m + nrhs * iMlaenv[("Rormbr" - 1) * ldiMlaenv]);
                maxwrk = max(maxwrk, m * m + 4 * m + (m - 1) * iMlaenv[("Rormbr" - 1) * ldiMlaenv]);
                if (nrhs > 1) {
                    maxwrk = max(maxwrk, m * m + m + m * nrhs);
                } else {
                    maxwrk = max(maxwrk, m * m + 2 * m);
                }
                maxwrk = max(maxwrk, m + nrhs * iMlaenv[("Rormlq" - 1) * ldiMlaenv]);
                maxwrk = max(maxwrk, m * m + 4 * m + wlalsd);
                //     XXX: Ensure the Path 2a case below is triggered.  The workspace
                //     calculation should use queries for all routines eventually.
                maxwrk = max(maxwrk, 4 * m + m * m + max(m, 2 * m - 4, nrhs, n - 3 * m));
            } else {
                //
                //              Path 2 - remaining underdetermined cases.
                //
                maxwrk = 3 * m + (n + m) * iMlaenv[("Rgebrd" - 1) * ldiMlaenv];
                maxwrk = max(maxwrk, 3 * m + nrhs * iMlaenv[("Rormbr" - 1) * ldiMlaenv]);
                maxwrk = max(maxwrk, 3 * m + m * iMlaenv[("Rormbr" - 1) * ldiMlaenv]);
                maxwrk = max(maxwrk, 3 * m + wlalsd);
            }
            minwrk = max(3 * m + nrhs, 3 * m + m, 3 * m + wlalsd);
        }
        minwrk = min(minwrk, maxwrk);
        work[1 - 1] = maxwrk;
        iwork[1 - 1] = liwork;
        //
        if (lwork < minwrk && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rgelsd", -info);
        return;
    } else if (lquery) {
        goto statement_10;
    }
    //
    //     Quick return if possible.
    //
    if (m == 0 || n == 0) {
        rank = 0;
        return;
    }
    //
    //     Get machine parameters.
    //
    eps = dlamch("P");
    sfmin = dlamch("S");
    smlnum = sfmin / eps;
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Scale A if max entry outside range [SMLNUM,BIGNUM].
    //
    anrm = Rlange[("M" - 1) + (m - 1) * ldRlange];
    iascl = 0;
    if (anrm > zero && anrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM.
        //
        Rlascl("G", 0, 0, anrm, smlnum, m, n, a, lda, info);
        iascl = 1;
    } else if (anrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM.
        //
        Rlascl("G", 0, 0, anrm, bignum, m, n, a, lda, info);
        iascl = 2;
    } else if (anrm == zero) {
        //
        //        Matrix all zero. Return zero solution.
        //
        Rlaset("F", max(m, n), nrhs, zero, zero, b, ldb);
        Rlaset("F", minmn, 1, zero, zero, s, 1);
        rank = 0;
        goto statement_10;
    }
    //
    //     Scale B if max entry outside range [SMLNUM,BIGNUM].
    //
    bnrm = Rlange[("M" - 1) + (m - 1) * ldRlange];
    ibscl = 0;
    if (bnrm > zero && bnrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM.
        //
        Rlascl("G", 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info);
        ibscl = 1;
    } else if (bnrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM.
        //
        Rlascl("G", 0, 0, bnrm, bignum, m, nrhs, b, ldb, info);
        ibscl = 2;
    }
    //
    //     If M < N make sure certain entries of B are zero.
    //
    if (m < n) {
        Rlaset("F", n - m, nrhs, zero, zero, b[((m + 1) - 1)], ldb);
    }
    //
    //     Overdetermined case.
    //
    if (m >= n) {
        //
        //        Path 1 - overdetermined or exactly determined.
        //
        mm = m;
        if (m >= mnthr) {
            //
            //           Path 1a - overdetermined, with many more rows than columns.
            //
            mm = n;
            itau = 1;
            nwork = itau + n;
            //
            //           Compute A=Q*R.
            //           (Workspace: need 2*N, prefer N+N*NB)
            //
            Rgeqrf(m, n, a, lda, work[itau - 1], work[nwork - 1], lwork - nwork + 1, info);
            //
            //           Multiply B by transpose(Q).
            //           (Workspace: need N+NRHS, prefer N+NRHS*NB)
            //
            Rormqr("L", "T", m, nrhs, n, a, lda, work[itau - 1], b, ldb, work[nwork - 1], lwork - nwork + 1, info);
            //
            //           Zero out below R.
            //
            if (n > 1) {
                Rlaset("L", n - 1, n - 1, zero, zero, a[(2 - 1)], lda);
            }
        }
        //
        ie = 1;
        itauq = ie + n;
        itaup = itauq + n;
        nwork = itaup + n;
        //
        //        Bidiagonalize R in A.
        //        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
        //
        Rgebrd(mm, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors of R.
        //        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
        //
        Rormbr("Q", "L", "T", mm, nrhs, n, a, lda, work[itauq - 1], b, ldb, work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Solve the bidiagonal least squares problem.
        //
        Rlalsd("U", smlsiz, n, nrhs, s, work[ie - 1], b, ldb, rcond, rank, work[nwork - 1], iwork, info);
        if (info != 0) {
            goto statement_10;
        }
        //
        //        Multiply B by right bidiagonalizing vectors of R.
        //
        Rormbr("P", "L", "N", n, nrhs, n, a, lda, work[itaup - 1], b, ldb, work[nwork - 1], lwork - nwork + 1, info);
        //
    } else if (n >= mnthr && lwork >= 4 * m + m * m + max(m, 2 * m - 4, nrhs, n - 3 * m, wlalsd)) {
        //
        //        Path 2a - underdetermined, with many more columns than rows
        //        and sufficient workspace for an efficient algorithm.
        //
        ldwork = m;
        if (lwork >= max(4 * m + m * lda + max(m, 2 * m - 4, nrhs, n - 3 * m), m * lda + m + m * nrhs, 4 * m + m * lda + wlalsd)) {
            ldwork = lda;
        }
        itau = 1;
        nwork = m + 1;
        //
        //        Compute A=L*Q.
        //        (Workspace: need 2*M, prefer M+M*NB)
        //
        Rgelqf(m, n, a, lda, work[itau - 1], work[nwork - 1], lwork - nwork + 1, info);
        il = nwork;
        //
        //        Copy L to WORK(IL), zeroing out above its diagonal.
        //
        Rlacpy("L", m, m, a, lda, work[il - 1], ldwork);
        Rlaset("U", m - 1, m - 1, zero, zero, work[(il + ldwork) - 1], ldwork);
        ie = il + ldwork * m;
        itauq = ie + m;
        itaup = itauq + m;
        nwork = itaup + m;
        //
        //        Bidiagonalize L in WORK(IL).
        //        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
        //
        Rgebrd(m, m, work[il - 1], ldwork, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors of L.
        //        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
        //
        Rormbr("Q", "L", "T", m, nrhs, m, work[il - 1], ldwork, work[itauq - 1], b, ldb, work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Solve the bidiagonal least squares problem.
        //
        Rlalsd("U", smlsiz, m, nrhs, s, work[ie - 1], b, ldb, rcond, rank, work[nwork - 1], iwork, info);
        if (info != 0) {
            goto statement_10;
        }
        //
        //        Multiply B by right bidiagonalizing vectors of L.
        //
        Rormbr("P", "L", "N", m, nrhs, m, work[il - 1], ldwork, work[itaup - 1], b, ldb, work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Zero out below first M rows of B.
        //
        Rlaset("F", n - m, nrhs, zero, zero, b[((m + 1) - 1)], ldb);
        nwork = itau + m;
        //
        //        Multiply transpose(Q) by B.
        //        (Workspace: need M+NRHS, prefer M+NRHS*NB)
        //
        Rormlq("L", "T", n, nrhs, m, a, lda, work[itau - 1], b, ldb, work[nwork - 1], lwork - nwork + 1, info);
        //
    } else {
        //
        //        Path 2 - remaining underdetermined cases.
        //
        ie = 1;
        itauq = ie + m;
        itaup = itauq + m;
        nwork = itaup + m;
        //
        //        Bidiagonalize A.
        //        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
        //
        Rgebrd(m, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors.
        //        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
        //
        Rormbr("Q", "L", "T", m, nrhs, n, a, lda, work[itauq - 1], b, ldb, work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Solve the bidiagonal least squares problem.
        //
        Rlalsd("L", smlsiz, m, nrhs, s, work[ie - 1], b, ldb, rcond, rank, work[nwork - 1], iwork, info);
        if (info != 0) {
            goto statement_10;
        }
        //
        //        Multiply B by right bidiagonalizing vectors of A.
        //
        Rormbr("P", "L", "N", n, nrhs, m, a, lda, work[itaup - 1], b, ldb, work[nwork - 1], lwork - nwork + 1, info);
        //
    }
    //
    //     Undo scaling.
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
statement_10:
    work[1 - 1] = maxwrk;
    iwork[1 - 1] = liwork;
    //
    //     End of Rgelsd
    //
}
