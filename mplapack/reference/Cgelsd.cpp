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

void Cgelsd(INTEGER const m, INTEGER const n, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL *s, REAL const rcond, INTEGER &rank, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER *iwork, INTEGER &info) {
    INTEGER minmn = 0;
    INTEGER maxmn = 0;
    bool lquery = false;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER liwork = 0;
    INTEGER lrwork = 0;
    INTEGER smlsiz = 0;
    INTEGER mnthr = 0;
    const REAL two = 2.0e+0;
    INTEGER nlvl = 0;
    INTEGER mm = 0;
    REAL eps = 0.0;
    REAL sfmin = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL anrm = 0.0;
    INTEGER iascl = 0;
    const REAL zero = 0.0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    REAL bnrm = 0.0;
    INTEGER ibscl = 0;
    INTEGER itau = 0;
    INTEGER nwork = 0;
    INTEGER itauq = 0;
    INTEGER itaup = 0;
    INTEGER ie = 0;
    INTEGER nrwork = 0;
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
    //     Compute workspace.
    //     (Note: Comments in the code beginning "Workspace:" describe the
    //     minimal amount of workspace needed at that point in the code,
    //     as well as the preferred amount for good performance.
    //     NB refers to the optimal block size for the immediately
    //     following subroutine, as returned by iMlaenv.)
    //
    if (info == 0) {
        minwrk = 1;
        maxwrk = 1;
        liwork = 1;
        lrwork = 1;
        if (minmn > 0) {
            smlsiz = iMlaenv(9, "Cgelsd", " ", 0, 0, 0, 0);
            mnthr = iMlaenv(6, "Cgelsd", " ", m, n, nrhs, -1);
            nlvl = max(castINTEGER(log(castREAL(minmn) / castREAL(smlsiz + 1)) / log(two)) + 1, (INTEGER)0);
            liwork = 3 * minmn * nlvl + 11 * minmn;
            mm = m;
            if (m >= n && m >= mnthr) {
                //
                //              Path 1a - overdetermined, with many more rows than
                //                        columns.
                //
                mm = n;
                maxwrk = max({maxwrk, n * iMlaenv(1, "Cgeqrf", " ", m, n, -1, -1)});
                maxwrk = max({maxwrk, nrhs * iMlaenv(1, "Cunmqr", "LC", m, nrhs, n, -1)});
            }
            if (m >= n) {
                //
                //              Path 1 - overdetermined or exactly determined.
                //
                lrwork = 10 * n + 2 * n * smlsiz + 8 * n * nlvl + 3 * smlsiz * nrhs + max(pow2((smlsiz + 1)), n * (1 + nrhs) + 2 * nrhs);
                maxwrk = max({maxwrk, 2 * n + (mm + n) * iMlaenv(1, "Cgebrd", " ", mm, n, -1, -1)});
                maxwrk = max({maxwrk, 2 * n + nrhs * iMlaenv(1, "Cunmbr", "QLC", mm, nrhs, n, -1)});
                maxwrk = max({maxwrk, 2 * n + (n - 1) * iMlaenv(1, "Cunmbr", "PLN", n, nrhs, n, -1)});
                maxwrk = max(maxwrk, 2 * n + n * nrhs);
                minwrk = max((INTEGER)2 * n + mm, 2 * n + n * nrhs);
            }
            if (n > m) {
                lrwork = 10 * m + 2 * m * smlsiz + 8 * m * nlvl + 3 * smlsiz * nrhs + max((smlsiz + 1) * (smlsiz + 1), INTEGER(n * (1 + nrhs) + 2 * nrhs));
                if (n >= mnthr) {
                    //
                    //                 Path 2a - underdetermined, with many more columns
                    //                           than rows.
                    //
                    maxwrk = m + m * iMlaenv(1, "Cgelqf", " ", m, n, -1, -1);
                    maxwrk = max({maxwrk, m * m + 4 * m + 2 * m * iMlaenv(1, "Cgebrd", " ", m, m, -1, -1)});
                    maxwrk = max({maxwrk, m * m + 4 * m + nrhs * iMlaenv(1, "Cunmbr", "QLC", m, nrhs, m, -1)});
                    maxwrk = max({maxwrk, m * m + 4 * m + (m - 1) * iMlaenv(1, "Cunmlq", "LC", n, nrhs, m, -1)});
                    if (nrhs > 1) {
                        maxwrk = max(maxwrk, m * m + m + m * nrhs);
                    } else {
                        maxwrk = max(maxwrk, m * m + 2 * m);
                    }
                    maxwrk = max(maxwrk, m * m + 4 * m + m * nrhs);
                    //     XXX: Ensure the Path 2a case below is triggered.  The workspace
                    //     calculation should use queries for all routines eventually.
                    maxwrk = max({maxwrk, 4 * m + m * m + max({m, 2 * m - 4, nrhs, n - 3 * m})});
                } else {
                    //
                    //                 Path 2 - underdetermined.
                    //
                    maxwrk = 2 * m + (n + m) * iMlaenv(1, "Cgebrd", " ", m, n, -1, -1);
                    maxwrk = max({maxwrk, 2 * m + nrhs * iMlaenv(1, "Cunmbr", "QLC", m, nrhs, m, -1)});
                    maxwrk = max({maxwrk, 2 * m + m * iMlaenv(1, "Cunmbr", "PLN", n, nrhs, m, -1)});
                    maxwrk = max(maxwrk, 2 * m + m * nrhs);
                }
                minwrk = max((INTEGER)2 * m + n, 2 * m + m * nrhs);
            }
        }
        minwrk = min(minwrk, maxwrk);
        work[1 - 1] = maxwrk;
        iwork[1 - 1] = liwork;
        rwork[1 - 1] = lrwork;
        //
        if (lwork < minwrk && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cgelsd", -info);
        return;
    } else if (lquery) {
        return;
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
    eps = Rlamch("P");
    sfmin = Rlamch("S");
    smlnum = sfmin / eps;
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Scale A if max entry outside range [SMLNUM,BIGNUM].
    //
    anrm = Clange("M", m, n, a, lda, rwork);
    iascl = 0;
    if (anrm > zero && anrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM
        //
        Clascl("G", 0, 0, anrm, smlnum, m, n, a, lda, info);
        iascl = 1;
    } else if (anrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM.
        //
        Clascl("G", 0, 0, anrm, bignum, m, n, a, lda, info);
        iascl = 2;
    } else if (anrm == zero) {
        //
        //        Matrix all zero. Return zero solution.
        //
        Claset("F", max(m, n), nrhs, czero, czero, b, ldb);
        Rlaset("F", minmn, 1, zero, zero, s, 1);
        rank = 0;
        goto statement_10;
    }
    //
    //     Scale B if max entry outside range [SMLNUM,BIGNUM].
    //
    bnrm = Clange("M", m, nrhs, b, ldb, rwork);
    ibscl = 0;
    if (bnrm > zero && bnrm < smlnum) {
        //
        //        Scale matrix norm up to SMLNUM.
        //
        Clascl("G", 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info);
        ibscl = 1;
    } else if (bnrm > bignum) {
        //
        //        Scale matrix norm down to BIGNUM.
        //
        Clascl("G", 0, 0, bnrm, bignum, m, nrhs, b, ldb, info);
        ibscl = 2;
    }
    //
    //     If M < N make sure B(M+1:N,:) = 0
    //
    if (m < n) {
        Claset("F", n - m, nrhs, czero, czero, &b[((m + 1) - 1)], ldb);
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
            //           Path 1a - overdetermined, with many more rows than columns
            //
            mm = n;
            itau = 1;
            nwork = itau + n;
            //
            //           Compute A=Q*R.
            //           (RWorkspace: need N)
            //           (CWorkspace: need N, prefer N*NB)
            //
            Cgeqrf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, info);
            //
            //           Multiply B by transpose(Q).
            //           (RWorkspace: need N)
            //           (CWorkspace: need NRHS, prefer NRHS*NB)
            //
            Cunmqr("L", "C", m, nrhs, n, a, lda, &work[itau - 1], b, ldb, &work[nwork - 1], lwork - nwork + 1, info);
            //
            //           Zero out below R.
            //
            if (n > 1) {
                Claset("L", n - 1, n - 1, czero, czero, &a[(2 - 1)], lda);
            }
        }
        //
        itauq = 1;
        itaup = itauq + n;
        nwork = itaup + n;
        ie = 1;
        nrwork = ie + n;
        //
        //        Bidiagonalize R in A.
        //        (RWorkspace: need N)
        //        (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
        //
        Cgebrd(mm, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors of R.
        //        (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
        //
        Cunmbr("Q", "L", "C", mm, nrhs, n, a, lda, &work[itauq - 1], b, ldb, &work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Solve the bidiagonal least squares problem.
        //
        Clalsd("U", smlsiz, n, nrhs, s, &rwork[ie - 1], b, ldb, rcond, rank, &work[nwork - 1], &rwork[nrwork - 1], iwork, info);
        if (info != 0) {
            goto statement_10;
        }
        //
        //        Multiply B by right bidiagonalizing vectors of R.
        //
        Cunmbr("P", "L", "N", n, nrhs, n, a, lda, &work[itaup - 1], b, ldb, &work[nwork - 1], lwork - nwork + 1, info);
        //
    } else if (n >= mnthr && lwork >= 4 * m + m * m + max({m, 2 * m - 4, nrhs, n - 3 * m})) {
        //
        //        Path 2a - underdetermined, with many more columns than rows
        //        and sufficient workspace for an efficient algorithm.
        //
        ldwork = m;
        if (lwork >= max({4 * m + m * lda + max({m, 2 * m - 4, nrhs, n - 3 * m}), m * lda + m + m * nrhs})) {
            ldwork = lda;
        }
        itau = 1;
        nwork = m + 1;
        //
        //        Compute A=L*Q.
        //        (CWorkspace: need 2*M, prefer M+M*NB)
        //
        Cgelqf(m, n, a, lda, &work[itau - 1], &work[nwork - 1], lwork - nwork + 1, info);
        il = nwork;
        //
        //        Copy L to WORK(IL), zeroing out above its diagonal.
        //
        Clacpy("L", m, m, a, lda, &work[il - 1], ldwork);
        Claset("U", m - 1, m - 1, czero, czero, &work[(il + ldwork) - 1], ldwork);
        itauq = il + ldwork * m;
        itaup = itauq + m;
        nwork = itaup + m;
        ie = 1;
        nrwork = ie + m;
        //
        //        Bidiagonalize L in WORK(IL).
        //        (RWorkspace: need M)
        //        (CWorkspace: need M*M+4*M, prefer M*M+4*M+2*M*NB)
        //
        Cgebrd(m, m, &work[il - 1], ldwork, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors of L.
        //        (CWorkspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
        //
        Cunmbr("Q", "L", "C", m, nrhs, m, &work[il - 1], ldwork, &work[itauq - 1], b, ldb, &work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Solve the bidiagonal least squares problem.
        //
        Clalsd("U", smlsiz, m, nrhs, s, &rwork[ie - 1], b, ldb, rcond, rank, &work[nwork - 1], &rwork[nrwork - 1], iwork, info);
        if (info != 0) {
            goto statement_10;
        }
        //
        //        Multiply B by right bidiagonalizing vectors of L.
        //
        Cunmbr("P", "L", "N", m, nrhs, m, &work[il - 1], ldwork, &work[itaup - 1], b, ldb, &work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Zero out below first M rows of B.
        //
        Claset("F", n - m, nrhs, czero, czero, &b[((m + 1) - 1)], ldb);
        nwork = itau + m;
        //
        //        Multiply transpose(Q) by B.
        //        (CWorkspace: need NRHS, prefer NRHS*NB)
        //
        Cunmlq("L", "C", n, nrhs, m, a, lda, &work[itau - 1], b, ldb, &work[nwork - 1], lwork - nwork + 1, info);
        //
    } else {
        //
        //        Path 2 - remaining underdetermined cases.
        //
        itauq = 1;
        itaup = itauq + m;
        nwork = itaup + m;
        ie = 1;
        nrwork = ie + m;
        //
        //        Bidiagonalize A.
        //        (RWorkspace: need M)
        //        (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
        //
        Cgebrd(m, n, a, lda, s, &rwork[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Multiply B by transpose of left bidiagonalizing vectors.
        //        (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
        //
        Cunmbr("Q", "L", "C", m, nrhs, n, a, lda, &work[itauq - 1], b, ldb, &work[nwork - 1], lwork - nwork + 1, info);
        //
        //        Solve the bidiagonal least squares problem.
        //
        Clalsd("L", smlsiz, m, nrhs, s, &rwork[ie - 1], b, ldb, rcond, rank, &work[nwork - 1], &rwork[nrwork - 1], iwork, info);
        if (info != 0) {
            goto statement_10;
        }
        //
        //        Multiply B by right bidiagonalizing vectors of A.
        //
        Cunmbr("P", "L", "N", n, nrhs, m, a, lda, &work[itaup - 1], b, ldb, &work[nwork - 1], lwork - nwork + 1, info);
        //
    }
    //
    //     Undo scaling.
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
//
statement_10:
    work[1 - 1] = maxwrk;
    iwork[1 - 1] = liwork;
    rwork[1 - 1] = lrwork;
    //
    //     End of Cgelsd
    //
}
