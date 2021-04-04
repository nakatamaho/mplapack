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

void Rgesvd(const char *jobu, const char *jobvt, INTEGER const &m, INTEGER const &n, REAL *a, INTEGER const &lda, REAL *s, REAL *u, INTEGER const &ldu, REAL *vt, INTEGER const &ldvt, REAL *work, INTEGER const &lwork, INTEGER &info) {
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
    //       minimal amount of workspace needed at that poINTEGER in the code,
    //       as well as the preferred amount for good performance.
    //       NB refers to the optimal block size for the immediately
    //       following subroutine, as returned by iMlaenv.)
    //
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER mnthr = 0;
    INTEGER bdspac = 0;
    arr_1d<1, REAL> dum(fill0);
    INTEGER ierr = 0;
    INTEGER lwork_Rgeqrf = 0;
    INTEGER lwork_Rorgqr_n = 0;
    INTEGER lwork_Rorgqr_m = 0;
    INTEGER lwork_Rgebrd = 0;
    INTEGER lwork_Rorgbr_p = 0;
    INTEGER lwork_Rorgbr_q = 0;
    INTEGER wrkbl = 0;
    INTEGER lwork_Rgelqf = 0;
    INTEGER lwork_Rorglq_n = 0;
    INTEGER lwork_Rorglq_m = 0;
    if (info == 0) {
        minwrk = 1;
        maxwrk = 1;
        if (m >= n && minmn > 0) {
            //
            //           Compute space needed for Rbdsqr
            //
            mnthr = iMlaenv[(6 - 1) + ("Rgesvd" - 1) * ldiMlaenv];
            bdspac = 5 * n;
            //           Compute space needed for Rgeqrf
            Rgeqrf(m, n, a, lda, dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rgeqrf = INTEGER(dum[1 - 1]);
            //           Compute space needed for Rorgqr
            Rorgqr(m, n, n, a, lda, dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rorgqr_n = INTEGER(dum[1 - 1]);
            Rorgqr(m, m, n, a, lda, dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rorgqr_m = INTEGER(dum[1 - 1]);
            //           Compute space needed for Rgebrd
            Rgebrd(n, n, a, lda, s, dum[1 - 1], dum[1 - 1], dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rgebrd = INTEGER(dum[1 - 1]);
            //           Compute space needed for Rorgbr P
            Rorgbr("P", n, n, n, a, lda, dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rorgbr_p = INTEGER(dum[1 - 1]);
            //           Compute space needed for Rorgbr Q
            Rorgbr("Q", n, n, n, a, lda, dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rorgbr_q = INTEGER(dum[1 - 1]);
            //
            if (m >= mnthr) {
                if (wntun) {
                    //
                    //                 Path 1 (M much larger than N, JOBU='N')
                    //
                    maxwrk = n + lwork_Rgeqrf;
                    maxwrk = max(maxwrk, 3 * n + lwork_Rgebrd);
                    if (wntvo || wntvas) {
                        maxwrk = max(maxwrk, 3 * n + lwork_Rorgbr_p);
                    }
                    maxwrk = max(maxwrk, bdspac);
                    minwrk = max(4 * n, bdspac);
                } else if (wntuo && wntvn) {
                    //
                    //                 Path 2 (M much larger than N, JOBU='O', JOBVT='N')
                    //
                    wrkbl = n + lwork_Rgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_n);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = max(n * n + wrkbl, n * n + m * n + n);
                    minwrk = max(3 * n + m, bdspac);
                } else if (wntuo && wntvas) {
                    //
                    //                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
                    //                 'A')
                    //
                    wrkbl = n + lwork_Rgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_n);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = max(n * n + wrkbl, n * n + m * n + n);
                    minwrk = max(3 * n + m, bdspac);
                } else if (wntus && wntvn) {
                    //
                    //                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
                    //
                    wrkbl = n + lwork_Rgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_n);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = n * n + wrkbl;
                    minwrk = max(3 * n + m, bdspac);
                } else if (wntus && wntvo) {
                    //
                    //                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
                    //
                    wrkbl = n + lwork_Rgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_n);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = 2 * n * n + wrkbl;
                    minwrk = max(3 * n + m, bdspac);
                } else if (wntus && wntvas) {
                    //
                    //                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
                    //                 'A')
                    //
                    wrkbl = n + lwork_Rgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_n);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = n * n + wrkbl;
                    minwrk = max(3 * n + m, bdspac);
                } else if (wntua && wntvn) {
                    //
                    //                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
                    //
                    wrkbl = n + lwork_Rgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_m);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = n * n + wrkbl;
                    minwrk = max(3 * n + m, bdspac);
                } else if (wntua && wntvo) {
                    //
                    //                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
                    //
                    wrkbl = n + lwork_Rgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_m);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = 2 * n * n + wrkbl;
                    minwrk = max(3 * n + m, bdspac);
                } else if (wntua && wntvas) {
                    //
                    //                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
                    //                 'A')
                    //
                    wrkbl = n + lwork_Rgeqrf;
                    wrkbl = max(wrkbl, n + lwork_Rorgqr_m);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, 3 * n + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = n * n + wrkbl;
                    minwrk = max(3 * n + m, bdspac);
                }
            } else {
                //
                //              Path 10 (M at least N, but not much larger)
                //
                Rgebrd(m, n, a, lda, s, dum[1 - 1], dum[1 - 1], dum[1 - 1], dum[1 - 1], -1, ierr);
                lwork_Rgebrd = INTEGER(dum[1 - 1]);
                maxwrk = 3 * n + lwork_Rgebrd;
                if (wntus || wntuo) {
                    Rorgbr("Q", m, n, n, a, lda, dum[1 - 1], dum[1 - 1], -1, ierr);
                    lwork_Rorgbr_q = INTEGER(dum[1 - 1]);
                    maxwrk = max(maxwrk, 3 * n + lwork_Rorgbr_q);
                }
                if (wntua) {
                    Rorgbr("Q", m, m, n, a, lda, dum[1 - 1], dum[1 - 1], -1, ierr);
                    lwork_Rorgbr_q = INTEGER(dum[1 - 1]);
                    maxwrk = max(maxwrk, 3 * n + lwork_Rorgbr_q);
                }
                if (!wntvn) {
                    maxwrk = max(maxwrk, 3 * n + lwork_Rorgbr_p);
                }
                maxwrk = max(maxwrk, bdspac);
                minwrk = max(3 * n + m, bdspac);
            }
        } else if (minmn > 0) {
            //
            //           Compute space needed for Rbdsqr
            //
            mnthr = iMlaenv[(6 - 1) + ("Rgesvd" - 1) * ldiMlaenv];
            bdspac = 5 * m;
            //           Compute space needed for Rgelqf
            Rgelqf(m, n, a, lda, dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rgelqf = INTEGER(dum[1 - 1]);
            //           Compute space needed for Rorglq
            Rorglq(n, n, m, dum[1 - 1], n, dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rorglq_n = INTEGER(dum[1 - 1]);
            Rorglq(m, n, m, a, lda, dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rorglq_m = INTEGER(dum[1 - 1]);
            //           Compute space needed for Rgebrd
            Rgebrd(m, m, a, lda, s, dum[1 - 1], dum[1 - 1], dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rgebrd = INTEGER(dum[1 - 1]);
            //            Compute space needed for Rorgbr P
            Rorgbr("P", m, m, m, a, n, dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rorgbr_p = INTEGER(dum[1 - 1]);
            //           Compute space needed for Rorgbr Q
            Rorgbr("Q", m, m, m, a, n, dum[1 - 1], dum[1 - 1], -1, ierr);
            lwork_Rorgbr_q = INTEGER(dum[1 - 1]);
            if (n >= mnthr) {
                if (wntvn) {
                    //
                    //                 Path 1t(N much larger than M, JOBVT='N')
                    //
                    maxwrk = m + lwork_Rgelqf;
                    maxwrk = max(maxwrk, 3 * m + lwork_Rgebrd);
                    if (wntuo || wntuas) {
                        maxwrk = max(maxwrk, 3 * m + lwork_Rorgbr_q);
                    }
                    maxwrk = max(maxwrk, bdspac);
                    minwrk = max(4 * m, bdspac);
                } else if (wntvo && wntun) {
                    //
                    //                 Path 2t(N much larger than M, JOBU='N', JOBVT='O')
                    //
                    wrkbl = m + lwork_Rgelqf;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_m);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = max(m * m + wrkbl, m * m + m * n + m);
                    minwrk = max(3 * m + n, bdspac);
                } else if (wntvo && wntuas) {
                    //
                    //                 Path 3t(N much larger than M, JOBU='S' or 'A',
                    //                 JOBVT='O')
                    //
                    wrkbl = m + lwork_Rgelqf;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_m);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = max(m * m + wrkbl, m * m + m * n + m);
                    minwrk = max(3 * m + n, bdspac);
                } else if (wntvs && wntun) {
                    //
                    //                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
                    //
                    wrkbl = m + lwork_Rgelqf;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_m);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = m * m + wrkbl;
                    minwrk = max(3 * m + n, bdspac);
                } else if (wntvs && wntuo) {
                    //
                    //                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
                    //
                    wrkbl = m + lwork_Rgelqf;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_m);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = 2 * m * m + wrkbl;
                    minwrk = max(3 * m + n, bdspac);
                } else if (wntvs && wntuas) {
                    //
                    //                 Path 6t(N much larger than M, JOBU='S' or 'A',
                    //                 JOBVT='S')
                    //
                    wrkbl = m + lwork_Rgelqf;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_m);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = m * m + wrkbl;
                    minwrk = max(3 * m + n, bdspac);
                } else if (wntva && wntun) {
                    //
                    //                 Path 7t(N much larger than M, JOBU='N', JOBVT='A')
                    //
                    wrkbl = m + lwork_Rgelqf;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_n);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = m * m + wrkbl;
                    minwrk = max(3 * m + n, bdspac);
                } else if (wntva && wntuo) {
                    //
                    //                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
                    //
                    wrkbl = m + lwork_Rgelqf;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_n);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = 2 * m * m + wrkbl;
                    minwrk = max(3 * m + n, bdspac);
                } else if (wntva && wntuas) {
                    //
                    //                 Path 9t(N much larger than M, JOBU='S' or 'A',
                    //                 JOBVT='A')
                    //
                    wrkbl = m + lwork_Rgelqf;
                    wrkbl = max(wrkbl, m + lwork_Rorglq_n);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rgebrd);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_p);
                    wrkbl = max(wrkbl, 3 * m + lwork_Rorgbr_q);
                    wrkbl = max(wrkbl, bdspac);
                    maxwrk = m * m + wrkbl;
                    minwrk = max(3 * m + n, bdspac);
                }
            } else {
                //
                //              Path 10t(N greater than M, but not much larger)
                //
                Rgebrd(m, n, a, lda, s, dum[1 - 1], dum[1 - 1], dum[1 - 1], dum[1 - 1], -1, ierr);
                lwork_Rgebrd = INTEGER(dum[1 - 1]);
                maxwrk = 3 * m + lwork_Rgebrd;
                if (wntvs || wntvo) {
                    //                Compute space needed for Rorgbr P
                    Rorgbr("P", m, n, m, a, n, dum[1 - 1], dum[1 - 1], -1, ierr);
                    lwork_Rorgbr_p = INTEGER(dum[1 - 1]);
                    maxwrk = max(maxwrk, 3 * m + lwork_Rorgbr_p);
                }
                if (wntva) {
                    Rorgbr("P", n, n, m, a, n, dum[1 - 1], dum[1 - 1], -1, ierr);
                    lwork_Rorgbr_p = INTEGER(dum[1 - 1]);
                    maxwrk = max(maxwrk, 3 * m + lwork_Rorgbr_p);
                }
                if (!wntun) {
                    maxwrk = max(maxwrk, 3 * m + lwork_Rorgbr_q);
                }
                maxwrk = max(maxwrk, bdspac);
                minwrk = max(3 * m + n, bdspac);
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
        Mxerbla("Rgesvd", -info);
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
    REAL eps = dlamch("P");
    REAL smlnum = sqrt(dlamch("S")) / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    //
    //     Scale A if max element outside range [SMLNUM,BIGNUM]
    //
    REAL anrm = Rlange[("M" - 1) + (m - 1) * ldRlange];
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
    INTEGER iwork = 0;
    INTEGER ie = 0;
    INTEGER itauq = 0;
    INTEGER itaup = 0;
    INTEGER ncvt = 0;
    INTEGER ir = 0;
    INTEGER ldwrku = 0;
    INTEGER ldwrkr = 0;
    INTEGER iu = 0;
    INTEGER i = 0;
    INTEGER chunk = 0;
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
                //              (Workspace: need 2*N, prefer N + N*NB)
                //
                Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                //
                //              Zero out below R
                //
                if (n > 1) {
                    Rlaset("L", n - 1, n - 1, zero, zero, a[(2 - 1)], lda);
                }
                ie = 1;
                itauq = ie + n;
                itaup = itauq + n;
                iwork = itaup + n;
                //
                //              Bidiagonalize R in A
                //              (Workspace: need 4*N, prefer 3*N + 2*N*NB)
                //
                Rgebrd(n, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                ncvt = 0;
                if (wntvo || wntvas) {
                    //
                    //                 If right singular vectors desired, generate P'.
                    //                 (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
                    //
                    Rorgbr("P", n, n, n, a, lda, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    ncvt = n;
                }
                iwork = ie + n;
                //
                //              Perform bidiagonal QR iteration, computing right
                //              singular vectors of A in A if desired
                //              (Workspace: need BDSPAC)
                //
                Rbdsqr("U", n, ncvt, 0, 0, s, work[ie - 1], a, lda, dum, 1, dum, 1, work[iwork - 1], info);
                //
                //              If right singular vectors desired in VT, copy them there
                //
                if (wntvas) {
                    Rlacpy("F", n, n, a, lda, vt, ldvt);
                }
                //
            } else if (wntuo && wntvn) {
                //
                //              Path 2 (M much larger than N, JOBU='O', JOBVT='N')
                //              N left singular vectors to be overwritten on A and
                //              no right singular vectors to be computed
                //
                if (lwork >= n * n + max(4 * n, bdspac)) {
                    //
                    //                 Sufficient workspace for a fast algorithm
                    //
                    ir = 1;
                    if (lwork >= max(wrkbl, lda * n + n) + lda * n) {
                        //
                        //                    WORK(IU) is LDA by N, WORK(IR) is LDA by N
                        //
                        ldwrku = lda;
                        ldwrkr = lda;
                    } else if (lwork >= max(wrkbl, lda * n + n) + n * n) {
                        //
                        //                    WORK(IU) is LDA by N, WORK(IR) is N by N
                        //
                        ldwrku = lda;
                        ldwrkr = n;
                    } else {
                        //
                        //                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N
                        //
                        ldwrku = (lwork - n * n - n) / n;
                        ldwrkr = n;
                    }
                    itau = ir + ldwrkr * n;
                    iwork = itau + n;
                    //
                    //                 Compute A=Q*R
                    //                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
                    //
                    Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy R to WORK(IR) and zero out below it
                    //
                    Rlacpy("U", n, n, a, lda, work[ir - 1], ldwrkr);
                    Rlaset("L", n - 1, n - 1, zero, zero, work[(ir + 1) - 1], ldwrkr);
                    //
                    //                 Generate Q in A
                    //                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
                    //
                    Rorgqr(m, n, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = itau;
                    itauq = ie + n;
                    itaup = itauq + n;
                    iwork = itaup + n;
                    //
                    //                 Bidiagonalize R in WORK(IR)
                    //                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
                    //
                    Rgebrd(n, n, work[ir - 1], ldwrkr, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate left vectors bidiagonalizing R
                    //                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
                    //
                    Rorgbr("Q", n, n, n, work[ir - 1], ldwrkr, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    iwork = ie + n;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of R in WORK(IR)
                    //                 (Workspace: need N*N + BDSPAC)
                    //
                    Rbdsqr("U", n, 0, n, 0, s, work[ie - 1], dum, 1, work[ir - 1], ldwrkr, dum, 1, work[iwork - 1], info);
                    iu = ie + n;
                    //
                    //                 Multiply Q in A by left singular vectors of R in
                    //                 WORK(IR), storing result in WORK(IU) and copying to A
                    //                 (Workspace: need N*N + 2*N, prefer N*N + M*N + N)
                    //
                    for (i = 1; i <= m; i = i + ldwrku) {
                        chunk = min(m - i + 1, ldwrku);
                        Rgemm("N", "N", chunk, n, n, one, a[(i - 1)], lda, work[ir - 1], ldwrkr, zero, work[iu - 1], ldwrku);
                        Rlacpy("F", chunk, n, work[iu - 1], ldwrku, a[(i - 1)], lda);
                    }
                    //
                } else {
                    //
                    //                 Insufficient workspace for a fast algorithm
                    //
                    ie = 1;
                    itauq = ie + n;
                    itaup = itauq + n;
                    iwork = itaup + n;
                    //
                    //                 Bidiagonalize A
                    //                 (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB)
                    //
                    Rgebrd(m, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate left vectors bidiagonalizing A
                    //                 (Workspace: need 4*N, prefer 3*N + N*NB)
                    //
                    Rorgbr("Q", m, n, n, a, lda, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    iwork = ie + n;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of A in A
                    //                 (Workspace: need BDSPAC)
                    //
                    Rbdsqr("U", n, 0, m, 0, s, work[ie - 1], dum, 1, a, lda, dum, 1, work[iwork - 1], info);
                    //
                }
                //
            } else if (wntuo && wntvas) {
                //
                //              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
                //              N left singular vectors to be overwritten on A and
                //              N right singular vectors to be computed in VT
                //
                if (lwork >= n * n + max(4 * n, bdspac)) {
                    //
                    //                 Sufficient workspace for a fast algorithm
                    //
                    ir = 1;
                    if (lwork >= max(wrkbl, lda * n + n) + lda * n) {
                        //
                        //                    WORK(IU) is LDA by N and WORK(IR) is LDA by N
                        //
                        ldwrku = lda;
                        ldwrkr = lda;
                    } else if (lwork >= max(wrkbl, lda * n + n) + n * n) {
                        //
                        //                    WORK(IU) is LDA by N and WORK(IR) is N by N
                        //
                        ldwrku = lda;
                        ldwrkr = n;
                    } else {
                        //
                        //                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N
                        //
                        ldwrku = (lwork - n * n - n) / n;
                        ldwrkr = n;
                    }
                    itau = ir + ldwrkr * n;
                    iwork = itau + n;
                    //
                    //                 Compute A=Q*R
                    //                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
                    //
                    Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy R to VT, zeroing out below it
                    //
                    Rlacpy("U", n, n, a, lda, vt, ldvt);
                    if (n > 1) {
                        Rlaset("L", n - 1, n - 1, zero, zero, vt[(2 - 1)], ldvt);
                    }
                    //
                    //                 Generate Q in A
                    //                 (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
                    //
                    Rorgqr(m, n, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = itau;
                    itauq = ie + n;
                    itaup = itauq + n;
                    iwork = itaup + n;
                    //
                    //                 Bidiagonalize R in VT, copying result to WORK(IR)
                    //                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
                    //
                    Rgebrd(n, n, vt, ldvt, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    Rlacpy("L", n, n, vt, ldvt, work[ir - 1], ldwrkr);
                    //
                    //                 Generate left vectors bidiagonalizing R in WORK(IR)
                    //                 (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
                    //
                    Rorgbr("Q", n, n, n, work[ir - 1], ldwrkr, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate right vectors bidiagonalizing R in VT
                    //                 (Workspace: need N*N + 4*N-1, prefer N*N + 3*N + (N-1)*NB)
                    //
                    Rorgbr("P", n, n, n, vt, ldvt, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    iwork = ie + n;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of R in WORK(IR) and computing right
                    //                 singular vectors of R in VT
                    //                 (Workspace: need N*N + BDSPAC)
                    //
                    Rbdsqr("U", n, n, n, 0, s, work[ie - 1], vt, ldvt, work[ir - 1], ldwrkr, dum, 1, work[iwork - 1], info);
                    iu = ie + n;
                    //
                    //                 Multiply Q in A by left singular vectors of R in
                    //                 WORK(IR), storing result in WORK(IU) and copying to A
                    //                 (Workspace: need N*N + 2*N, prefer N*N + M*N + N)
                    //
                    for (i = 1; i <= m; i = i + ldwrku) {
                        chunk = min(m - i + 1, ldwrku);
                        Rgemm("N", "N", chunk, n, n, one, a[(i - 1)], lda, work[ir - 1], ldwrkr, zero, work[iu - 1], ldwrku);
                        Rlacpy("F", chunk, n, work[iu - 1], ldwrku, a[(i - 1)], lda);
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
                    //                 (Workspace: need 2*N, prefer N + N*NB)
                    //
                    Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy R to VT, zeroing out below it
                    //
                    Rlacpy("U", n, n, a, lda, vt, ldvt);
                    if (n > 1) {
                        Rlaset("L", n - 1, n - 1, zero, zero, vt[(2 - 1)], ldvt);
                    }
                    //
                    //                 Generate Q in A
                    //                 (Workspace: need 2*N, prefer N + N*NB)
                    //
                    Rorgqr(m, n, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = itau;
                    itauq = ie + n;
                    itaup = itauq + n;
                    iwork = itaup + n;
                    //
                    //                 Bidiagonalize R in VT
                    //                 (Workspace: need 4*N, prefer 3*N + 2*N*NB)
                    //
                    Rgebrd(n, n, vt, ldvt, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Multiply Q in A by left vectors bidiagonalizing R
                    //                 (Workspace: need 3*N + M, prefer 3*N + M*NB)
                    //
                    Rormbr("Q", "R", "N", m, n, n, vt, ldvt, work[itauq - 1], a, lda, work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate right vectors bidiagonalizing R in VT
                    //                 (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
                    //
                    Rorgbr("P", n, n, n, vt, ldvt, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    iwork = ie + n;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of A in A and computing right
                    //                 singular vectors of A in VT
                    //                 (Workspace: need BDSPAC)
                    //
                    Rbdsqr("U", n, n, m, 0, s, work[ie - 1], vt, ldvt, a, lda, dum, 1, work[iwork - 1], info);
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
                    if (lwork >= n * n + max(4 * n, bdspac)) {
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
                        //                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to WORK(IR), zeroing out below it
                        //
                        Rlacpy("U", n, n, a, lda, work[ir - 1], ldwrkr);
                        Rlaset("L", n - 1, n - 1, zero, zero, work[(ir + 1) - 1], ldwrkr);
                        //
                        //                    Generate Q in A
                        //                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
                        //
                        Rorgqr(m, n, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IR)
                        //                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
                        //
                        Rgebrd(n, n, work[ir - 1], ldwrkr, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left vectors bidiagonalizing R in WORK(IR)
                        //                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
                        //
                        Rorgbr("Q", n, n, n, work[ir - 1], ldwrkr, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IR)
                        //                    (Workspace: need N*N + BDSPAC)
                        //
                        Rbdsqr("U", n, 0, n, 0, s, work[ie - 1], dum, 1, work[ir - 1], ldwrkr, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply Q in A by left singular vectors of R in
                        //                    WORK(IR), storing result in U
                        //                    (Workspace: need N*N)
                        //
                        Rgemm("N", "N", m, n, n, one, a, lda, work[ir - 1], ldwrkr, zero, u, ldu);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (Workspace: need 2*N, prefer N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (Workspace: need 2*N, prefer N + N*NB)
                        //
                        Rorgqr(m, n, n, u, ldu, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Zero out below R in A
                        //
                        if (n > 1) {
                            Rlaset("L", n - 1, n - 1, zero, zero, a[(2 - 1)], lda);
                        }
                        //
                        //                    Bidiagonalize R in A
                        //                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
                        //
                        Rgebrd(n, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left vectors bidiagonalizing R
                        //                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
                        //
                        Rormbr("Q", "R", "N", m, n, n, a, lda, work[itauq - 1], u, ldu, work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", n, 0, m, 0, s, work[ie - 1], dum, 1, u, ldu, dum, 1, work[iwork - 1], info);
                        //
                    }
                    //
                } else if (wntvo) {
                    //
                    //                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
                    //                 N left singular vectors to be computed in U and
                    //                 N right singular vectors to be overwritten on A
                    //
                    if (lwork >= 2 * n * n + max(4 * n, bdspac)) {
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
                        //                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to WORK(IU), zeroing out below it
                        //
                        Rlacpy("U", n, n, a, lda, work[iu - 1], ldwrku);
                        Rlaset("L", n - 1, n - 1, zero, zero, work[(iu + 1) - 1], ldwrku);
                        //
                        //                    Generate Q in A
                        //                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)
                        //
                        Rorgqr(m, n, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IU), copying result to
                        //                    WORK(IR)
                        //                    (Workspace: need 2*N*N + 4*N,
                        //                                prefer 2*N*N+3*N+2*N*NB)
                        //
                        Rgebrd(n, n, work[iu - 1], ldwrku, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", n, n, work[iu - 1], ldwrku, work[ir - 1], ldwrkr);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IU)
                        //                    (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB)
                        //
                        Rorgbr("Q", n, n, n, work[iu - 1], ldwrku, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IR)
                        //                    (Workspace: need 2*N*N + 4*N-1,
                        //                                prefer 2*N*N+3*N+(N-1)*NB)
                        //
                        Rorgbr("P", n, n, n, work[ir - 1], ldwrkr, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IU) and computing
                        //                    right singular vectors of R in WORK(IR)
                        //                    (Workspace: need 2*N*N + BDSPAC)
                        //
                        Rbdsqr("U", n, n, n, 0, s, work[ie - 1], work[ir - 1], ldwrkr, work[iu - 1], ldwrku, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply Q in A by left singular vectors of R in
                        //                    WORK(IU), storing result in U
                        //                    (Workspace: need N*N)
                        //
                        Rgemm("N", "N", m, n, n, one, a, lda, work[iu - 1], ldwrku, zero, u, ldu);
                        //
                        //                    Copy right singular vectors of R to A
                        //                    (Workspace: need N*N)
                        //
                        Rlacpy("F", n, n, work[ir - 1], ldwrkr, a, lda);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (Workspace: need 2*N, prefer N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (Workspace: need 2*N, prefer N + N*NB)
                        //
                        Rorgqr(m, n, n, u, ldu, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Zero out below R in A
                        //
                        if (n > 1) {
                            Rlaset("L", n - 1, n - 1, zero, zero, a[(2 - 1)], lda);
                        }
                        //
                        //                    Bidiagonalize R in A
                        //                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
                        //
                        Rgebrd(n, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left vectors bidiagonalizing R
                        //                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
                        //
                        Rormbr("Q", "R", "N", m, n, n, a, lda, work[itauq - 1], u, ldu, work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right vectors bidiagonalizing R in A
                        //                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
                        //
                        Rorgbr("P", n, n, n, a, lda, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in A
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", n, n, m, 0, s, work[ie - 1], a, lda, u, ldu, dum, 1, work[iwork - 1], info);
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
                    if (lwork >= n * n + max(4 * n, bdspac)) {
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
                        //                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to WORK(IU), zeroing out below it
                        //
                        Rlacpy("U", n, n, a, lda, work[iu - 1], ldwrku);
                        Rlaset("L", n - 1, n - 1, zero, zero, work[(iu + 1) - 1], ldwrku);
                        //
                        //                    Generate Q in A
                        //                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
                        //
                        Rorgqr(m, n, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IU), copying result to VT
                        //                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
                        //
                        Rgebrd(n, n, work[iu - 1], ldwrku, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", n, n, work[iu - 1], ldwrku, vt, ldvt);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IU)
                        //                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
                        //
                        Rorgbr("Q", n, n, n, work[iu - 1], ldwrku, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in VT
                        //                    (Workspace: need N*N + 4*N-1,
                        //                                prefer N*N+3*N+(N-1)*NB)
                        //
                        Rorgbr("P", n, n, n, vt, ldvt, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IU) and computing
                        //                    right singular vectors of R in VT
                        //                    (Workspace: need N*N + BDSPAC)
                        //
                        Rbdsqr("U", n, n, n, 0, s, work[ie - 1], vt, ldvt, work[iu - 1], ldwrku, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply Q in A by left singular vectors of R in
                        //                    WORK(IU), storing result in U
                        //                    (Workspace: need N*N)
                        //
                        Rgemm("N", "N", m, n, n, one, a, lda, work[iu - 1], ldwrku, zero, u, ldu);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (Workspace: need 2*N, prefer N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (Workspace: need 2*N, prefer N + N*NB)
                        //
                        Rorgqr(m, n, n, u, ldu, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to VT, zeroing out below it
                        //
                        Rlacpy("U", n, n, a, lda, vt, ldvt);
                        if (n > 1) {
                            Rlaset("L", n - 1, n - 1, zero, zero, vt[(2 - 1)], ldvt);
                        }
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in VT
                        //                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
                        //
                        Rgebrd(n, n, vt, ldvt, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left bidiagonalizing vectors
                        //                    in VT
                        //                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
                        //
                        Rormbr("Q", "R", "N", m, n, n, vt, ldvt, work[itauq - 1], u, ldu, work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in VT
                        //                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
                        //
                        Rorgbr("P", n, n, n, vt, ldvt, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in VT
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", n, n, m, 0, s, work[ie - 1], vt, ldvt, u, ldu, dum, 1, work[iwork - 1], info);
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
                    if (lwork >= n * n + max(n + m, 4 * n, bdspac)) {
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
                        //                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Copy R to WORK(IR), zeroing out below it
                        //
                        Rlacpy("U", n, n, a, lda, work[ir - 1], ldwrkr);
                        Rlaset("L", n - 1, n - 1, zero, zero, work[(ir + 1) - 1], ldwrkr);
                        //
                        //                    Generate Q in U
                        //                    (Workspace: need N*N + N + M, prefer N*N + N + M*NB)
                        //
                        Rorgqr(m, m, n, u, ldu, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IR)
                        //                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
                        //
                        Rgebrd(n, n, work[ir - 1], ldwrkr, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IR)
                        //                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
                        //
                        Rorgbr("Q", n, n, n, work[ir - 1], ldwrkr, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IR)
                        //                    (Workspace: need N*N + BDSPAC)
                        //
                        Rbdsqr("U", n, 0, n, 0, s, work[ie - 1], dum, 1, work[ir - 1], ldwrkr, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply Q in U by left singular vectors of R in
                        //                    WORK(IR), storing result in A
                        //                    (Workspace: need N*N)
                        //
                        Rgemm("N", "N", m, n, n, one, u, ldu, work[ir - 1], ldwrkr, zero, a, lda);
                        //
                        //                    Copy left singular vectors of A from A to U
                        //
                        Rlacpy("F", m, n, a, lda, u, ldu);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (Workspace: need 2*N, prefer N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (Workspace: need N + M, prefer N + M*NB)
                        //
                        Rorgqr(m, m, n, u, ldu, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Zero out below R in A
                        //
                        if (n > 1) {
                            Rlaset("L", n - 1, n - 1, zero, zero, a[(2 - 1)], lda);
                        }
                        //
                        //                    Bidiagonalize R in A
                        //                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
                        //
                        Rgebrd(n, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left bidiagonalizing vectors
                        //                    in A
                        //                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
                        //
                        Rormbr("Q", "R", "N", m, n, n, a, lda, work[itauq - 1], u, ldu, work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", n, 0, m, 0, s, work[ie - 1], dum, 1, u, ldu, dum, 1, work[iwork - 1], info);
                        //
                    }
                    //
                } else if (wntvo) {
                    //
                    //                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
                    //                 M left singular vectors to be computed in U and
                    //                 N right singular vectors to be overwritten on A
                    //
                    if (lwork >= 2 * n * n + max(n + m, 4 * n, bdspac)) {
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
                        //                    (Workspace: need 2*N*N + 2*N, prefer 2*N*N + N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (Workspace: need 2*N*N + N + M, prefer 2*N*N + N + M*NB)
                        //
                        Rorgqr(m, m, n, u, ldu, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to WORK(IU), zeroing out below it
                        //
                        Rlacpy("U", n, n, a, lda, work[iu - 1], ldwrku);
                        Rlaset("L", n - 1, n - 1, zero, zero, work[(iu + 1) - 1], ldwrku);
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IU), copying result to
                        //                    WORK(IR)
                        //                    (Workspace: need 2*N*N + 4*N,
                        //                                prefer 2*N*N+3*N+2*N*NB)
                        //
                        Rgebrd(n, n, work[iu - 1], ldwrku, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", n, n, work[iu - 1], ldwrku, work[ir - 1], ldwrkr);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IU)
                        //                    (Workspace: need 2*N*N + 4*N, prefer 2*N*N + 3*N + N*NB)
                        //
                        Rorgbr("Q", n, n, n, work[iu - 1], ldwrku, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IR)
                        //                    (Workspace: need 2*N*N + 4*N-1,
                        //                                prefer 2*N*N+3*N+(N-1)*NB)
                        //
                        Rorgbr("P", n, n, n, work[ir - 1], ldwrkr, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IU) and computing
                        //                    right singular vectors of R in WORK(IR)
                        //                    (Workspace: need 2*N*N + BDSPAC)
                        //
                        Rbdsqr("U", n, n, n, 0, s, work[ie - 1], work[ir - 1], ldwrkr, work[iu - 1], ldwrku, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply Q in U by left singular vectors of R in
                        //                    WORK(IU), storing result in A
                        //                    (Workspace: need N*N)
                        //
                        Rgemm("N", "N", m, n, n, one, u, ldu, work[iu - 1], ldwrku, zero, a, lda);
                        //
                        //                    Copy left singular vectors of A from A to U
                        //
                        Rlacpy("F", m, n, a, lda, u, ldu);
                        //
                        //                    Copy right singular vectors of R from WORK(IR) to A
                        //
                        Rlacpy("F", n, n, work[ir - 1], ldwrkr, a, lda);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (Workspace: need 2*N, prefer N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (Workspace: need N + M, prefer N + M*NB)
                        //
                        Rorgqr(m, m, n, u, ldu, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Zero out below R in A
                        //
                        if (n > 1) {
                            Rlaset("L", n - 1, n - 1, zero, zero, a[(2 - 1)], lda);
                        }
                        //
                        //                    Bidiagonalize R in A
                        //                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
                        //
                        Rgebrd(n, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left bidiagonalizing vectors
                        //                    in A
                        //                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
                        //
                        Rormbr("Q", "R", "N", m, n, n, a, lda, work[itauq - 1], u, ldu, work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in A
                        //                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
                        //
                        Rorgbr("P", n, n, n, a, lda, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in A
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", n, n, m, 0, s, work[ie - 1], a, lda, u, ldu, dum, 1, work[iwork - 1], info);
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
                    if (lwork >= n * n + max(n + m, 4 * n, bdspac)) {
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
                        //                    (Workspace: need N*N + 2*N, prefer N*N + N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (Workspace: need N*N + N + M, prefer N*N + N + M*NB)
                        //
                        Rorgqr(m, m, n, u, ldu, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R to WORK(IU), zeroing out below it
                        //
                        Rlacpy("U", n, n, a, lda, work[iu - 1], ldwrku);
                        Rlaset("L", n - 1, n - 1, zero, zero, work[(iu + 1) - 1], ldwrku);
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in WORK(IU), copying result to VT
                        //                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + 2*N*NB)
                        //
                        Rgebrd(n, n, work[iu - 1], ldwrku, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", n, n, work[iu - 1], ldwrku, vt, ldvt);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IU)
                        //                    (Workspace: need N*N + 4*N, prefer N*N + 3*N + N*NB)
                        //
                        Rorgbr("Q", n, n, n, work[iu - 1], ldwrku, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in VT
                        //                    (Workspace: need N*N + 4*N-1,
                        //                                prefer N*N+3*N+(N-1)*NB)
                        //
                        Rorgbr("P", n, n, n, vt, ldvt, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of R in WORK(IU) and computing
                        //                    right singular vectors of R in VT
                        //                    (Workspace: need N*N + BDSPAC)
                        //
                        Rbdsqr("U", n, n, n, 0, s, work[ie - 1], vt, ldvt, work[iu - 1], ldwrku, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply Q in U by left singular vectors of R in
                        //                    WORK(IU), storing result in A
                        //                    (Workspace: need N*N)
                        //
                        Rgemm("N", "N", m, n, n, one, u, ldu, work[iu - 1], ldwrku, zero, a, lda);
                        //
                        //                    Copy left singular vectors of A from A to U
                        //
                        Rlacpy("F", m, n, a, lda, u, ldu);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + n;
                        //
                        //                    Compute A=Q*R, copying result to U
                        //                    (Workspace: need 2*N, prefer N + N*NB)
                        //
                        Rgeqrf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, n, a, lda, u, ldu);
                        //
                        //                    Generate Q in U
                        //                    (Workspace: need N + M, prefer N + M*NB)
                        //
                        Rorgqr(m, m, n, u, ldu, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy R from A to VT, zeroing out below it
                        //
                        Rlacpy("U", n, n, a, lda, vt, ldvt);
                        if (n > 1) {
                            Rlaset("L", n - 1, n - 1, zero, zero, vt[(2 - 1)], ldvt);
                        }
                        ie = itau;
                        itauq = ie + n;
                        itaup = itauq + n;
                        iwork = itaup + n;
                        //
                        //                    Bidiagonalize R in VT
                        //                    (Workspace: need 4*N, prefer 3*N + 2*N*NB)
                        //
                        Rgebrd(n, n, vt, ldvt, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply Q in U by left bidiagonalizing vectors
                        //                    in VT
                        //                    (Workspace: need 3*N + M, prefer 3*N + M*NB)
                        //
                        Rormbr("Q", "R", "N", m, n, n, vt, ldvt, work[itauq - 1], u, ldu, work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in VT
                        //                    (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
                        //
                        Rorgbr("P", n, n, n, vt, ldvt, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + n;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in VT
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", n, n, m, 0, s, work[ie - 1], vt, ldvt, u, ldu, dum, 1, work[iwork - 1], info);
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
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //
            //           Bidiagonalize A
            //           (Workspace: need 3*N + M, prefer 3*N + (M + N)*NB)
            //
            Rgebrd(m, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
            if (wntuas) {
                //
                //              If left singular vectors desired in U, copy result to U
                //              and generate left bidiagonalizing vectors in U
                //              (Workspace: need 3*N + NCU, prefer 3*N + NCU*NB)
                //
                Rlacpy("L", m, n, a, lda, u, ldu);
                if (wntus) {
                    ncu = n;
                }
                if (wntua) {
                    ncu = m;
                }
                Rorgbr("Q", m, ncu, n, u, ldu, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntvas) {
                //
                //              If right singular vectors desired in VT, copy result to
                //              VT and generate right bidiagonalizing vectors in VT
                //              (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
                //
                Rlacpy("U", n, n, a, lda, vt, ldvt);
                Rorgbr("P", n, n, n, vt, ldvt, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntuo) {
                //
                //              If left singular vectors desired in A, generate left
                //              bidiagonalizing vectors in A
                //              (Workspace: need 4*N, prefer 3*N + N*NB)
                //
                Rorgbr("Q", m, n, n, a, lda, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntvo) {
                //
                //              If right singular vectors desired in A, generate right
                //              bidiagonalizing vectors in A
                //              (Workspace: need 4*N-1, prefer 3*N + (N-1)*NB)
                //
                Rorgbr("P", n, n, n, a, lda, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
            }
            iwork = ie + n;
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
                //              (Workspace: need BDSPAC)
                //
                Rbdsqr("U", n, ncvt, nru, 0, s, work[ie - 1], vt, ldvt, u, ldu, dum, 1, work[iwork - 1], info);
            } else if ((!wntuo) && wntvo) {
                //
                //              Perform bidiagonal QR iteration, if desired, computing
                //              left singular vectors in U and computing right singular
                //              vectors in A
                //              (Workspace: need BDSPAC)
                //
                Rbdsqr("U", n, ncvt, nru, 0, s, work[ie - 1], a, lda, u, ldu, dum, 1, work[iwork - 1], info);
            } else {
                //
                //              Perform bidiagonal QR iteration, if desired, computing
                //              left singular vectors in A and computing right singular
                //              vectors in VT
                //              (Workspace: need BDSPAC)
                //
                Rbdsqr("U", n, ncvt, nru, 0, s, work[ie - 1], vt, ldvt, a, lda, dum, 1, work[iwork - 1], info);
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
                //              (Workspace: need 2*M, prefer M + M*NB)
                //
                Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                //
                //              Zero out above L
                //
                Rlaset("U", m - 1, m - 1, zero, zero, a[(2 - 1) * lda], lda);
                ie = 1;
                itauq = ie + m;
                itaup = itauq + m;
                iwork = itaup + m;
                //
                //              Bidiagonalize L in A
                //              (Workspace: need 4*M, prefer 3*M + 2*M*NB)
                //
                Rgebrd(m, m, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                if (wntuo || wntuas) {
                    //
                    //                 If left singular vectors desired, generate Q
                    //                 (Workspace: need 4*M, prefer 3*M + M*NB)
                    //
                    Rorgbr("Q", m, m, m, a, lda, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                }
                iwork = ie + m;
                nru = 0;
                if (wntuo || wntuas) {
                    nru = m;
                }
                //
                //              Perform bidiagonal QR iteration, computing left singular
                //              vectors of A in A if desired
                //              (Workspace: need BDSPAC)
                //
                Rbdsqr("U", m, 0, nru, 0, s, work[ie - 1], dum, 1, a, lda, dum, 1, work[iwork - 1], info);
                //
                //              If left singular vectors desired in U, copy them there
                //
                if (wntuas) {
                    Rlacpy("F", m, m, a, lda, u, ldu);
                }
                //
            } else if (wntvo && wntun) {
                //
                //              Path 2t(N much larger than M, JOBU='N', JOBVT='O')
                //              M right singular vectors to be overwritten on A and
                //              no left singular vectors to be computed
                //
                if (lwork >= m * m + max(4 * m, bdspac)) {
                    //
                    //                 Sufficient workspace for a fast algorithm
                    //
                    ir = 1;
                    if (lwork >= max(wrkbl, lda * n + m) + lda * m) {
                        //
                        //                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
                        //
                        ldwrku = lda;
                        chunk = n;
                        ldwrkr = lda;
                    } else if (lwork >= max(wrkbl, lda * n + m) + m * m) {
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
                        chunk = (lwork - m * m - m) / m;
                        ldwrkr = m;
                    }
                    itau = ir + ldwrkr * m;
                    iwork = itau + m;
                    //
                    //                 Compute A=L*Q
                    //                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
                    //
                    Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy L to WORK(IR) and zero out above it
                    //
                    Rlacpy("L", m, m, a, lda, work[ir - 1], ldwrkr);
                    Rlaset("U", m - 1, m - 1, zero, zero, work[(ir + ldwrkr) - 1], ldwrkr);
                    //
                    //                 Generate Q in A
                    //                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
                    //
                    Rorglq(m, n, m, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = itau;
                    itauq = ie + m;
                    itaup = itauq + m;
                    iwork = itaup + m;
                    //
                    //                 Bidiagonalize L in WORK(IR)
                    //                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
                    //
                    Rgebrd(m, m, work[ir - 1], ldwrkr, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate right vectors bidiagonalizing L
                    //                 (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB)
                    //
                    Rorgbr("P", m, m, m, work[ir - 1], ldwrkr, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    iwork = ie + m;
                    //
                    //                 Perform bidiagonal QR iteration, computing right
                    //                 singular vectors of L in WORK(IR)
                    //                 (Workspace: need M*M + BDSPAC)
                    //
                    Rbdsqr("U", m, m, 0, 0, s, work[ie - 1], work[ir - 1], ldwrkr, dum, 1, dum, 1, work[iwork - 1], info);
                    iu = ie + m;
                    //
                    //                 Multiply right singular vectors of L in WORK(IR) by Q
                    //                 in A, storing result in WORK(IU) and copying to A
                    //                 (Workspace: need M*M + 2*M, prefer M*M + M*N + M)
                    //
                    for (i = 1; i <= n; i = i + chunk) {
                        blk = min(n - i + 1, chunk);
                        Rgemm("N", "N", m, blk, m, one, work[ir - 1], ldwrkr, a[(i - 1) * lda], lda, zero, work[iu - 1], ldwrku);
                        Rlacpy("F", m, blk, work[iu - 1], ldwrku, a[(i - 1) * lda], lda);
                    }
                    //
                } else {
                    //
                    //                 Insufficient workspace for a fast algorithm
                    //
                    ie = 1;
                    itauq = ie + m;
                    itaup = itauq + m;
                    iwork = itaup + m;
                    //
                    //                 Bidiagonalize A
                    //                 (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB)
                    //
                    Rgebrd(m, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate right vectors bidiagonalizing A
                    //                 (Workspace: need 4*M, prefer 3*M + M*NB)
                    //
                    Rorgbr("P", m, n, m, a, lda, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    iwork = ie + m;
                    //
                    //                 Perform bidiagonal QR iteration, computing right
                    //                 singular vectors of A in A
                    //                 (Workspace: need BDSPAC)
                    //
                    Rbdsqr("L", m, n, 0, 0, s, work[ie - 1], a, lda, dum, 1, dum, 1, work[iwork - 1], info);
                    //
                }
                //
            } else if (wntvo && wntuas) {
                //
                //              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
                //              M right singular vectors to be overwritten on A and
                //              M left singular vectors to be computed in U
                //
                if (lwork >= m * m + max(4 * m, bdspac)) {
                    //
                    //                 Sufficient workspace for a fast algorithm
                    //
                    ir = 1;
                    if (lwork >= max(wrkbl, lda * n + m) + lda * m) {
                        //
                        //                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
                        //
                        ldwrku = lda;
                        chunk = n;
                        ldwrkr = lda;
                    } else if (lwork >= max(wrkbl, lda * n + m) + m * m) {
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
                        chunk = (lwork - m * m - m) / m;
                        ldwrkr = m;
                    }
                    itau = ir + ldwrkr * m;
                    iwork = itau + m;
                    //
                    //                 Compute A=L*Q
                    //                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
                    //
                    Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy L to U, zeroing about above it
                    //
                    Rlacpy("L", m, m, a, lda, u, ldu);
                    Rlaset("U", m - 1, m - 1, zero, zero, u[(2 - 1) * ldu], ldu);
                    //
                    //                 Generate Q in A
                    //                 (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
                    //
                    Rorglq(m, n, m, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = itau;
                    itauq = ie + m;
                    itaup = itauq + m;
                    iwork = itaup + m;
                    //
                    //                 Bidiagonalize L in U, copying result to WORK(IR)
                    //                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
                    //
                    Rgebrd(m, m, u, ldu, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    Rlacpy("U", m, m, u, ldu, work[ir - 1], ldwrkr);
                    //
                    //                 Generate right vectors bidiagonalizing L in WORK(IR)
                    //                 (Workspace: need M*M + 4*M-1, prefer M*M + 3*M + (M-1)*NB)
                    //
                    Rorgbr("P", m, m, m, work[ir - 1], ldwrkr, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate left vectors bidiagonalizing L in U
                    //                 (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)
                    //
                    Rorgbr("Q", m, m, m, u, ldu, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    iwork = ie + m;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of L in U, and computing right
                    //                 singular vectors of L in WORK(IR)
                    //                 (Workspace: need M*M + BDSPAC)
                    //
                    Rbdsqr("U", m, m, m, 0, s, work[ie - 1], work[ir - 1], ldwrkr, u, ldu, dum, 1, work[iwork - 1], info);
                    iu = ie + m;
                    //
                    //                 Multiply right singular vectors of L in WORK(IR) by Q
                    //                 in A, storing result in WORK(IU) and copying to A
                    //                 (Workspace: need M*M + 2*M, prefer M*M + M*N + M))
                    //
                    for (i = 1; i <= n; i = i + chunk) {
                        blk = min(n - i + 1, chunk);
                        Rgemm("N", "N", m, blk, m, one, work[ir - 1], ldwrkr, a[(i - 1) * lda], lda, zero, work[iu - 1], ldwrku);
                        Rlacpy("F", m, blk, work[iu - 1], ldwrku, a[(i - 1) * lda], lda);
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
                    //                 (Workspace: need 2*M, prefer M + M*NB)
                    //
                    Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Copy L to U, zeroing out above it
                    //
                    Rlacpy("L", m, m, a, lda, u, ldu);
                    Rlaset("U", m - 1, m - 1, zero, zero, u[(2 - 1) * ldu], ldu);
                    //
                    //                 Generate Q in A
                    //                 (Workspace: need 2*M, prefer M + M*NB)
                    //
                    Rorglq(m, n, m, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    ie = itau;
                    itauq = ie + m;
                    itaup = itauq + m;
                    iwork = itaup + m;
                    //
                    //                 Bidiagonalize L in U
                    //                 (Workspace: need 4*M, prefer 3*M + 2*M*NB)
                    //
                    Rgebrd(m, m, u, ldu, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Multiply right vectors bidiagonalizing L by Q in A
                    //                 (Workspace: need 3*M + N, prefer 3*M + N*NB)
                    //
                    Rormbr("P", "L", "T", m, n, m, u, ldu, work[itaup - 1], a, lda, work[iwork - 1], lwork - iwork + 1, ierr);
                    //
                    //                 Generate left vectors bidiagonalizing L in U
                    //                 (Workspace: need 4*M, prefer 3*M + M*NB)
                    //
                    Rorgbr("Q", m, m, m, u, ldu, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                    iwork = ie + m;
                    //
                    //                 Perform bidiagonal QR iteration, computing left
                    //                 singular vectors of A in U and computing right
                    //                 singular vectors of A in A
                    //                 (Workspace: need BDSPAC)
                    //
                    Rbdsqr("U", m, n, m, 0, s, work[ie - 1], a, lda, u, ldu, dum, 1, work[iwork - 1], info);
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
                    if (lwork >= m * m + max(4 * m, bdspac)) {
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
                        //                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to WORK(IR), zeroing out above it
                        //
                        Rlacpy("L", m, m, a, lda, work[ir - 1], ldwrkr);
                        Rlaset("U", m - 1, m - 1, zero, zero, work[(ir + ldwrkr) - 1], ldwrkr);
                        //
                        //                    Generate Q in A
                        //                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
                        //
                        Rorglq(m, n, m, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IR)
                        //                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
                        //
                        Rgebrd(m, m, work[ir - 1], ldwrkr, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right vectors bidiagonalizing L in
                        //                    WORK(IR)
                        //                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB)
                        //
                        Rorgbr("P", m, m, m, work[ir - 1], ldwrkr, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing right
                        //                    singular vectors of L in WORK(IR)
                        //                    (Workspace: need M*M + BDSPAC)
                        //
                        Rbdsqr("U", m, m, 0, 0, s, work[ie - 1], work[ir - 1], ldwrkr, dum, 1, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IR) by
                        //                    Q in A, storing result in VT
                        //                    (Workspace: need M*M)
                        //
                        Rgemm("N", "N", m, n, m, one, work[ir - 1], ldwrkr, a, lda, zero, vt, ldvt);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q
                        //                    (Workspace: need 2*M, prefer M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy result to VT
                        //
                        Rlacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (Workspace: need 2*M, prefer M + M*NB)
                        //
                        Rorglq(m, n, m, vt, ldvt, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Zero out above L in A
                        //
                        Rlaset("U", m - 1, m - 1, zero, zero, a[(2 - 1) * lda], lda);
                        //
                        //                    Bidiagonalize L in A
                        //                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
                        //
                        Rgebrd(m, m, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right vectors bidiagonalizing L by Q in VT
                        //                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
                        //
                        Rormbr("P", "L", "T", m, n, m, a, lda, work[itaup - 1], vt, ldvt, work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing right
                        //                    singular vectors of A in VT
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", m, n, 0, 0, s, work[ie - 1], vt, ldvt, dum, 1, dum, 1, work[iwork - 1], info);
                        //
                    }
                    //
                } else if (wntuo) {
                    //
                    //                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
                    //                 M right singular vectors to be computed in VT and
                    //                 M left singular vectors to be overwritten on A
                    //
                    if (lwork >= 2 * m * m + max(4 * m, bdspac)) {
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
                        //                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to WORK(IU), zeroing out below it
                        //
                        Rlacpy("L", m, m, a, lda, work[iu - 1], ldwrku);
                        Rlaset("U", m - 1, m - 1, zero, zero, work[(iu + ldwrku) - 1], ldwrku);
                        //
                        //                    Generate Q in A
                        //                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)
                        //
                        Rorglq(m, n, m, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IU), copying result to
                        //                    WORK(IR)
                        //                    (Workspace: need 2*M*M + 4*M,
                        //                                prefer 2*M*M+3*M+2*M*NB)
                        //
                        Rgebrd(m, m, work[iu - 1], ldwrku, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, m, work[iu - 1], ldwrku, work[ir - 1], ldwrkr);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IU)
                        //                    (Workspace: need 2*M*M + 4*M-1,
                        //                                prefer 2*M*M+3*M+(M-1)*NB)
                        //
                        Rorgbr("P", m, m, m, work[iu - 1], ldwrku, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IR)
                        //                    (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB)
                        //
                        Rorgbr("Q", m, m, m, work[ir - 1], ldwrkr, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of L in WORK(IR) and computing
                        //                    right singular vectors of L in WORK(IU)
                        //                    (Workspace: need 2*M*M + BDSPAC)
                        //
                        Rbdsqr("U", m, m, m, 0, s, work[ie - 1], work[iu - 1], ldwrku, work[ir - 1], ldwrkr, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IU) by
                        //                    Q in A, storing result in VT
                        //                    (Workspace: need M*M)
                        //
                        Rgemm("N", "N", m, n, m, one, work[iu - 1], ldwrku, a, lda, zero, vt, ldvt);
                        //
                        //                    Copy left singular vectors of L to A
                        //                    (Workspace: need M*M)
                        //
                        Rlacpy("F", m, m, work[ir - 1], ldwrkr, a, lda);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (Workspace: need 2*M, prefer M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (Workspace: need 2*M, prefer M + M*NB)
                        //
                        Rorglq(m, n, m, vt, ldvt, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Zero out above L in A
                        //
                        Rlaset("U", m - 1, m - 1, zero, zero, a[(2 - 1) * lda], lda);
                        //
                        //                    Bidiagonalize L in A
                        //                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
                        //
                        Rgebrd(m, m, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right vectors bidiagonalizing L by Q in VT
                        //                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
                        //
                        Rormbr("P", "L", "T", m, n, m, a, lda, work[itaup - 1], vt, ldvt, work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors of L in A
                        //                    (Workspace: need 4*M, prefer 3*M + M*NB)
                        //
                        Rorgbr("Q", m, m, m, a, lda, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, compute left
                        //                    singular vectors of A in A and compute right
                        //                    singular vectors of A in VT
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", m, n, m, 0, s, work[ie - 1], vt, ldvt, a, lda, dum, 1, work[iwork - 1], info);
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
                    if (lwork >= m * m + max(4 * m, bdspac)) {
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
                        //                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to WORK(IU), zeroing out above it
                        //
                        Rlacpy("L", m, m, a, lda, work[iu - 1], ldwrku);
                        Rlaset("U", m - 1, m - 1, zero, zero, work[(iu + ldwrku) - 1], ldwrku);
                        //
                        //                    Generate Q in A
                        //                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
                        //
                        Rorglq(m, n, m, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IU), copying result to U
                        //                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
                        //
                        Rgebrd(m, m, work[iu - 1], ldwrku, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, m, work[iu - 1], ldwrku, u, ldu);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IU)
                        //                    (Workspace: need M*M + 4*M-1,
                        //                                prefer M*M+3*M+(M-1)*NB)
                        //
                        Rorgbr("P", m, m, m, work[iu - 1], ldwrku, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in U
                        //                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)
                        //
                        Rorgbr("Q", m, m, m, u, ldu, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of L in U and computing right
                        //                    singular vectors of L in WORK(IU)
                        //                    (Workspace: need M*M + BDSPAC)
                        //
                        Rbdsqr("U", m, m, m, 0, s, work[ie - 1], work[iu - 1], ldwrku, u, ldu, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IU) by
                        //                    Q in A, storing result in VT
                        //                    (Workspace: need M*M)
                        //
                        Rgemm("N", "N", m, n, m, one, work[iu - 1], ldwrku, a, lda, zero, vt, ldvt);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (Workspace: need 2*M, prefer M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (Workspace: need 2*M, prefer M + M*NB)
                        //
                        Rorglq(m, n, m, vt, ldvt, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to U, zeroing out above it
                        //
                        Rlacpy("L", m, m, a, lda, u, ldu);
                        Rlaset("U", m - 1, m - 1, zero, zero, u[(2 - 1) * ldu], ldu);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in U
                        //                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
                        //
                        Rgebrd(m, m, u, ldu, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right bidiagonalizing vectors in U by Q
                        //                    in VT
                        //                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
                        //
                        Rormbr("P", "L", "T", m, n, m, u, ldu, work[itaup - 1], vt, ldvt, work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in U
                        //                    (Workspace: need 4*M, prefer 3*M + M*NB)
                        //
                        Rorgbr("Q", m, m, m, u, ldu, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in VT
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", m, n, m, 0, s, work[ie - 1], vt, ldvt, u, ldu, dum, 1, work[iwork - 1], info);
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
                    if (lwork >= m * m + max(n + m, 4 * m, bdspac)) {
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
                        //                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Copy L to WORK(IR), zeroing out above it
                        //
                        Rlacpy("L", m, m, a, lda, work[ir - 1], ldwrkr);
                        Rlaset("U", m - 1, m - 1, zero, zero, work[(ir + ldwrkr) - 1], ldwrkr);
                        //
                        //                    Generate Q in VT
                        //                    (Workspace: need M*M + M + N, prefer M*M + M + N*NB)
                        //
                        Rorglq(n, n, m, vt, ldvt, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IR)
                        //                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
                        //
                        Rgebrd(m, m, work[ir - 1], ldwrkr, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IR)
                        //                    (Workspace: need M*M + 4*M-1,
                        //                                prefer M*M+3*M+(M-1)*NB)
                        //
                        Rorgbr("P", m, m, m, work[ir - 1], ldwrkr, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing right
                        //                    singular vectors of L in WORK(IR)
                        //                    (Workspace: need M*M + BDSPAC)
                        //
                        Rbdsqr("U", m, m, 0, 0, s, work[ie - 1], work[ir - 1], ldwrkr, dum, 1, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IR) by
                        //                    Q in VT, storing result in A
                        //                    (Workspace: need M*M)
                        //
                        Rgemm("N", "N", m, n, m, one, work[ir - 1], ldwrkr, vt, ldvt, zero, a, lda);
                        //
                        //                    Copy right singular vectors of A from A to VT
                        //
                        Rlacpy("F", m, n, a, lda, vt, ldvt);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (Workspace: need 2*M, prefer M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (Workspace: need M + N, prefer M + N*NB)
                        //
                        Rorglq(n, n, m, vt, ldvt, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Zero out above L in A
                        //
                        Rlaset("U", m - 1, m - 1, zero, zero, a[(2 - 1) * lda], lda);
                        //
                        //                    Bidiagonalize L in A
                        //                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
                        //
                        Rgebrd(m, m, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right bidiagonalizing vectors in A by Q
                        //                    in VT
                        //                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
                        //
                        Rormbr("P", "L", "T", m, n, m, a, lda, work[itaup - 1], vt, ldvt, work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing right
                        //                    singular vectors of A in VT
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", m, n, 0, 0, s, work[ie - 1], vt, ldvt, dum, 1, dum, 1, work[iwork - 1], info);
                        //
                    }
                    //
                } else if (wntuo) {
                    //
                    //                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
                    //                 N right singular vectors to be computed in VT and
                    //                 M left singular vectors to be overwritten on A
                    //
                    if (lwork >= 2 * m * m + max(n + m, 4 * m, bdspac)) {
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
                        //                    (Workspace: need 2*M*M + 2*M, prefer 2*M*M + M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (Workspace: need 2*M*M + M + N, prefer 2*M*M + M + N*NB)
                        //
                        Rorglq(n, n, m, vt, ldvt, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to WORK(IU), zeroing out above it
                        //
                        Rlacpy("L", m, m, a, lda, work[iu - 1], ldwrku);
                        Rlaset("U", m - 1, m - 1, zero, zero, work[(iu + ldwrku) - 1], ldwrku);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IU), copying result to
                        //                    WORK(IR)
                        //                    (Workspace: need 2*M*M + 4*M,
                        //                                prefer 2*M*M+3*M+2*M*NB)
                        //
                        Rgebrd(m, m, work[iu - 1], ldwrku, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, m, work[iu - 1], ldwrku, work[ir - 1], ldwrkr);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IU)
                        //                    (Workspace: need 2*M*M + 4*M-1,
                        //                                prefer 2*M*M+3*M+(M-1)*NB)
                        //
                        Rorgbr("P", m, m, m, work[iu - 1], ldwrku, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in WORK(IR)
                        //                    (Workspace: need 2*M*M + 4*M, prefer 2*M*M + 3*M + M*NB)
                        //
                        Rorgbr("Q", m, m, m, work[ir - 1], ldwrkr, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of L in WORK(IR) and computing
                        //                    right singular vectors of L in WORK(IU)
                        //                    (Workspace: need 2*M*M + BDSPAC)
                        //
                        Rbdsqr("U", m, m, m, 0, s, work[ie - 1], work[iu - 1], ldwrku, work[ir - 1], ldwrkr, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IU) by
                        //                    Q in VT, storing result in A
                        //                    (Workspace: need M*M)
                        //
                        Rgemm("N", "N", m, n, m, one, work[iu - 1], ldwrku, vt, ldvt, zero, a, lda);
                        //
                        //                    Copy right singular vectors of A from A to VT
                        //
                        Rlacpy("F", m, n, a, lda, vt, ldvt);
                        //
                        //                    Copy left singular vectors of A from WORK(IR) to A
                        //
                        Rlacpy("F", m, m, work[ir - 1], ldwrkr, a, lda);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (Workspace: need 2*M, prefer M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (Workspace: need M + N, prefer M + N*NB)
                        //
                        Rorglq(n, n, m, vt, ldvt, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Zero out above L in A
                        //
                        Rlaset("U", m - 1, m - 1, zero, zero, a[(2 - 1) * lda], lda);
                        //
                        //                    Bidiagonalize L in A
                        //                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
                        //
                        Rgebrd(m, m, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right bidiagonalizing vectors in A by Q
                        //                    in VT
                        //                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
                        //
                        Rormbr("P", "L", "T", m, n, m, a, lda, work[itaup - 1], vt, ldvt, work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in A
                        //                    (Workspace: need 4*M, prefer 3*M + M*NB)
                        //
                        Rorgbr("Q", m, m, m, a, lda, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in A and computing right
                        //                    singular vectors of A in VT
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", m, n, m, 0, s, work[ie - 1], vt, ldvt, a, lda, dum, 1, work[iwork - 1], info);
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
                    if (lwork >= m * m + max(n + m, 4 * m, bdspac)) {
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
                        //                    (Workspace: need M*M + 2*M, prefer M*M + M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (Workspace: need M*M + M + N, prefer M*M + M + N*NB)
                        //
                        Rorglq(n, n, m, vt, ldvt, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to WORK(IU), zeroing out above it
                        //
                        Rlacpy("L", m, m, a, lda, work[iu - 1], ldwrku);
                        Rlaset("U", m - 1, m - 1, zero, zero, work[(iu + ldwrku) - 1], ldwrku);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in WORK(IU), copying result to U
                        //                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + 2*M*NB)
                        //
                        Rgebrd(m, m, work[iu - 1], ldwrku, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("L", m, m, work[iu - 1], ldwrku, u, ldu);
                        //
                        //                    Generate right bidiagonalizing vectors in WORK(IU)
                        //                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + (M-1)*NB)
                        //
                        Rorgbr("P", m, m, m, work[iu - 1], ldwrku, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in U
                        //                    (Workspace: need M*M + 4*M, prefer M*M + 3*M + M*NB)
                        //
                        Rorgbr("Q", m, m, m, u, ldu, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of L in U and computing right
                        //                    singular vectors of L in WORK(IU)
                        //                    (Workspace: need M*M + BDSPAC)
                        //
                        Rbdsqr("U", m, m, m, 0, s, work[ie - 1], work[iu - 1], ldwrku, u, ldu, dum, 1, work[iwork - 1], info);
                        //
                        //                    Multiply right singular vectors of L in WORK(IU) by
                        //                    Q in VT, storing result in A
                        //                    (Workspace: need M*M)
                        //
                        Rgemm("N", "N", m, n, m, one, work[iu - 1], ldwrku, vt, ldvt, zero, a, lda);
                        //
                        //                    Copy right singular vectors of A from A to VT
                        //
                        Rlacpy("F", m, n, a, lda, vt, ldvt);
                        //
                    } else {
                        //
                        //                    Insufficient workspace for a fast algorithm
                        //
                        itau = 1;
                        iwork = itau + m;
                        //
                        //                    Compute A=L*Q, copying result to VT
                        //                    (Workspace: need 2*M, prefer M + M*NB)
                        //
                        Rgelqf(m, n, a, lda, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        Rlacpy("U", m, n, a, lda, vt, ldvt);
                        //
                        //                    Generate Q in VT
                        //                    (Workspace: need M + N, prefer M + N*NB)
                        //
                        Rorglq(n, n, m, vt, ldvt, work[itau - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Copy L to U, zeroing out above it
                        //
                        Rlacpy("L", m, m, a, lda, u, ldu);
                        Rlaset("U", m - 1, m - 1, zero, zero, u[(2 - 1) * ldu], ldu);
                        ie = itau;
                        itauq = ie + m;
                        itaup = itauq + m;
                        iwork = itaup + m;
                        //
                        //                    Bidiagonalize L in U
                        //                    (Workspace: need 4*M, prefer 3*M + 2*M*NB)
                        //
                        Rgebrd(m, m, u, ldu, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Multiply right bidiagonalizing vectors in U by Q
                        //                    in VT
                        //                    (Workspace: need 3*M + N, prefer 3*M + N*NB)
                        //
                        Rormbr("P", "L", "T", m, n, m, u, ldu, work[itaup - 1], vt, ldvt, work[iwork - 1], lwork - iwork + 1, ierr);
                        //
                        //                    Generate left bidiagonalizing vectors in U
                        //                    (Workspace: need 4*M, prefer 3*M + M*NB)
                        //
                        Rorgbr("Q", m, m, m, u, ldu, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
                        iwork = ie + m;
                        //
                        //                    Perform bidiagonal QR iteration, computing left
                        //                    singular vectors of A in U and computing right
                        //                    singular vectors of A in VT
                        //                    (Workspace: need BDSPAC)
                        //
                        Rbdsqr("U", m, n, m, 0, s, work[ie - 1], vt, ldvt, u, ldu, dum, 1, work[iwork - 1], info);
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
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //
            //           Bidiagonalize A
            //           (Workspace: need 3*M + N, prefer 3*M + (M + N)*NB)
            //
            Rgebrd(m, n, a, lda, s, work[ie - 1], work[itauq - 1], work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
            if (wntuas) {
                //
                //              If left singular vectors desired in U, copy result to U
                //              and generate left bidiagonalizing vectors in U
                //              (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB)
                //
                Rlacpy("L", m, m, a, lda, u, ldu);
                Rorgbr("Q", m, m, n, u, ldu, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntvas) {
                //
                //              If right singular vectors desired in VT, copy result to
                //              VT and generate right bidiagonalizing vectors in VT
                //              (Workspace: need 3*M + NRVT, prefer 3*M + NRVT*NB)
                //
                Rlacpy("U", m, n, a, lda, vt, ldvt);
                if (wntva) {
                    nrvt = n;
                }
                if (wntvs) {
                    nrvt = m;
                }
                Rorgbr("P", nrvt, n, m, vt, ldvt, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntuo) {
                //
                //              If left singular vectors desired in A, generate left
                //              bidiagonalizing vectors in A
                //              (Workspace: need 4*M-1, prefer 3*M + (M-1)*NB)
                //
                Rorgbr("Q", m, m, n, a, lda, work[itauq - 1], work[iwork - 1], lwork - iwork + 1, ierr);
            }
            if (wntvo) {
                //
                //              If right singular vectors desired in A, generate right
                //              bidiagonalizing vectors in A
                //              (Workspace: need 4*M, prefer 3*M + M*NB)
                //
                Rorgbr("P", m, n, m, a, lda, work[itaup - 1], work[iwork - 1], lwork - iwork + 1, ierr);
            }
            iwork = ie + m;
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
                //              (Workspace: need BDSPAC)
                //
                Rbdsqr("L", m, ncvt, nru, 0, s, work[ie - 1], vt, ldvt, u, ldu, dum, 1, work[iwork - 1], info);
            } else if ((!wntuo) && wntvo) {
                //
                //              Perform bidiagonal QR iteration, if desired, computing
                //              left singular vectors in U and computing right singular
                //              vectors in A
                //              (Workspace: need BDSPAC)
                //
                Rbdsqr("L", m, ncvt, nru, 0, s, work[ie - 1], a, lda, u, ldu, dum, 1, work[iwork - 1], info);
            } else {
                //
                //              Perform bidiagonal QR iteration, if desired, computing
                //              left singular vectors in A and computing right singular
                //              vectors in VT
                //              (Workspace: need BDSPAC)
                //
                Rbdsqr("L", m, ncvt, nru, 0, s, work[ie - 1], vt, ldvt, a, lda, dum, 1, work[iwork - 1], info);
            }
            //
        }
        //
    }
    //
    //     If Rbdsqr failed to converge, copy unconverged superdiagonals
    //     to WORK( 2:MINMN )
    //
    if (info != 0) {
        if (ie > 2) {
            for (i = 1; i <= minmn - 1; i = i + 1) {
                work[(i + 1) - 1] = work[(i + ie - 1) - 1];
            }
        }
        if (ie < 2) {
            for (i = minmn - 1; i >= 1; i = i - 1) {
                work[(i + 1) - 1] = work[(i + ie - 1) - 1];
            }
        }
    }
    //
    //     Undo scaling if necessary
    //
    if (iscl == 1) {
        if (anrm > bignum) {
            Rlascl("G", 0, 0, bignum, anrm, minmn, 1, s, minmn, ierr);
        }
        if (info != 0 && anrm > bignum) {
            Rlascl("G", 0, 0, bignum, anrm, minmn - 1, 1, work[2 - 1], minmn, ierr);
        }
        if (anrm < smlnum) {
            Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, s, minmn, ierr);
        }
        if (info != 0 && anrm < smlnum) {
            Rlascl("G", 0, 0, smlnum, anrm, minmn - 1, 1, work[2 - 1], minmn, ierr);
        }
    }
    //
    //     Return optimal workspace in WORK(1)
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Rgesvd
    //
}
