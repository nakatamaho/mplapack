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

void Rggev(const char *jobvl, const char *jobvr, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *alphar, REAL *alphai, REAL *beta, REAL *vl, INTEGER const ldvl, REAL *vr, INTEGER const ldvr, REAL *work, INTEGER const lwork, INTEGER &info) {
    INTEGER ijobvl = 0;
    bool ilvl = false;
    INTEGER ijobvr = 0;
    bool ilvr = false;
    bool ilv = false;
    bool lquery = false;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL anrm = 0.0;
    bool ilascl = false;
    const REAL zero = 0.0;
    REAL anrmto = 0.0;
    INTEGER ierr = 0;
    REAL bnrm = 0.0;
    bool ilbscl = false;
    REAL bnrmto = 0.0;
    INTEGER ileft = 0;
    INTEGER iright = 0;
    INTEGER iwrk = 0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    INTEGER irows = 0;
    INTEGER icols = 0;
    INTEGER itau = 0;
    char chtemp;
    bool ldumma;
    INTEGER in = 0;
    INTEGER jc = 0;
    REAL temp = 0.0;
    INTEGER jr = 0;
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
    //     Decode the input arguments
    //
    if (Mlsame(jobvl, "N")) {
        ijobvl = 1;
        ilvl = false;
    } else if (Mlsame(jobvl, "V")) {
        ijobvl = 2;
        ilvl = true;
    } else {
        ijobvl = -1;
        ilvl = false;
    }
    //
    if (Mlsame(jobvr, "N")) {
        ijobvr = 1;
        ilvr = false;
    } else if (Mlsame(jobvr, "V")) {
        ijobvr = 2;
        ilvr = true;
    } else {
        ijobvr = -1;
        ilvr = false;
    }
    ilv = ilvl || ilvr;
    //
    //     Test the input arguments
    //
    info = 0;
    lquery = (lwork == -1);
    if (ijobvl <= 0) {
        info = -1;
    } else if (ijobvr <= 0) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldvl < 1 || (ilvl && ldvl < n)) {
        info = -12;
    } else if (ldvr < 1 || (ilvr && ldvr < n)) {
        info = -14;
    }
    //
    //     Compute workspace
    //      (Note: Comments in the code beginning "Workspace:" describe the
    //       minimal amount of workspace needed at that poINTEGER in the code,
    //       as well as the preferred amount for good performance.
    //       NB refers to the optimal block size for the immediately
    //       following subroutine, as returned by iMlaenv. The workspace is
    //       computed assuming ILO = 1 and IHI = N, the worst case.)
    //
    if (info == 0) {
        minwrk = max((INTEGER)1, 8 * n);
        maxwrk = max((INTEGER)1, n * (7 + iMlaenv(1, "Rgeqrf", " ", n, 1, n, 0)));
        maxwrk = max(maxwrk, n * (7 + iMlaenv(1, "Rormqr", " ", n, 1, n, 0)));
        if (ilvl) {
            maxwrk = max(maxwrk, n * (7 + iMlaenv(1, "Rorgqr", " ", n, 1, n, -1)));
        }
        work[1 - 1] = maxwrk;
        //
        if (lwork < minwrk && !lquery) {
            info = -16;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rggev ", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Get machine constants
    //
    eps = Rlamch("P");
    smlnum = Rlamch("S");
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = one / smlnum;
    //
    //     Scale A if max element outside range [SMLNUM,BIGNUM]
    //
    anrm = Rlange("M", n, n, a, lda, work);
    ilascl = false;
    if (anrm > zero && anrm < smlnum) {
        anrmto = smlnum;
        ilascl = true;
    } else if (anrm > bignum) {
        anrmto = bignum;
        ilascl = true;
    }
    if (ilascl) {
        Rlascl("G", 0, 0, anrm, anrmto, n, n, a, lda, ierr);
    }
    //
    //     Scale B if max element outside range [SMLNUM,BIGNUM]
    //
    bnrm = Rlange("M", n, n, b, ldb, work);
    ilbscl = false;
    if (bnrm > zero && bnrm < smlnum) {
        bnrmto = smlnum;
        ilbscl = true;
    } else if (bnrm > bignum) {
        bnrmto = bignum;
        ilbscl = true;
    }
    if (ilbscl) {
        Rlascl("G", 0, 0, bnrm, bnrmto, n, n, b, ldb, ierr);
    }
    //
    //     Permute the matrices A, B to isolate eigenvalues if possible
    //     (Workspace: need 6*N)
    //
    ileft = 1;
    iright = n + 1;
    iwrk = iright + n;
    Rggbal("P", n, a, lda, b, ldb, ilo, ihi, &work[ileft - 1], &work[iright - 1], &work[iwrk - 1], ierr);
    //
    //     Reduce B to triangular form (QR decomposition of B)
    //     (Workspace: need N, prefer N*NB)
    //
    irows = ihi + 1 - ilo;
    if (ilv) {
        icols = n + 1 - ilo;
    } else {
        icols = irows;
    }
    itau = iwrk;
    iwrk = itau + irows;
    Rgeqrf(irows, icols, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, &work[itau - 1], &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    //     Apply the orthogonal transformation to matrix A
    //     (Workspace: need N, prefer N*NB)
    //
    Rormqr("L", "T", irows, icols, irows, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, &work[itau - 1], &a[(ilo - 1) + (ilo - 1) * lda], lda, &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    //     Initialize VL
    //     (Workspace: need N, prefer N*NB)
    //
    if (ilvl) {
        Rlaset("Full", n, n, zero, one, vl, ldvl);
        if (irows > 1) {
            Rlacpy("L", irows - 1, irows - 1, &b[((ilo + 1) - 1) + (ilo - 1) * ldb], ldb, &vl[((ilo + 1) - 1) + (ilo - 1) * ldvl], ldvl);
        }
        Rorgqr(irows, irows, irows, &vl[(ilo - 1) + (ilo - 1) * ldvl], ldvl, &work[itau - 1], &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    }
    //
    //     Initialize VR
    //
    if (ilvr) {
        Rlaset("Full", n, n, zero, one, vr, ldvr);
    }
    //
    //     Reduce to generalized Hessenberg form
    //     (Workspace: none needed)
    //
    if (ilv) {
        //
        //        Eigenvectors requested -- work on whole matrix.
        //
        Rgghrd(jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl, ldvl, vr, ldvr, ierr);
    } else {
        Rgghrd("N", "N", irows, 1, irows, &a[(ilo - 1) + (ilo - 1) * lda], lda, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, vl, ldvl, vr, ldvr, ierr);
    }
    //
    //     Perform QZ algorithm (Compute eigenvalues, and optionally, the
    //     Schur forms and Schur vectors)
    //     (Workspace: need N)
    //
    iwrk = itau;
    if (ilv) {
        chtemp = 'S';
    } else {
        chtemp = 'E';
    }
    Rhgeqz(&chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    if (ierr != 0) {
        if (ierr > 0 && ierr <= n) {
            info = ierr;
        } else if (ierr > n && ierr <= 2 * n) {
            info = ierr - n;
        } else {
            info = n + 1;
        }
        goto statement_110;
    }
    //
    //     Compute Eigenvectors
    //     (Workspace: need 6*N)
    //
    if (ilv) {
        if (ilvl) {
            if (ilvr) {
                chtemp = 'B';
            } else {
                chtemp = 'L';
            }
        } else {
            chtemp = 'R';
        }
        Rtgevc(&chtemp, "B", &ldumma, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, n, in, &work[iwrk - 1], ierr);
        if (ierr != 0) {
            info = n + 2;
            goto statement_110;
        }
        //
        //        Undo balancing on VL and VR and normalization
        //        (Workspace: none needed)
        //
        if (ilvl) {
            Rggbak("P", "L", n, ilo, ihi, &work[ileft - 1], &work[iright - 1], n, vl, ldvl, ierr);
            for (jc = 1; jc <= n; jc = jc + 1) {
                if (alphai[jc - 1] < zero) {
                    goto statement_50;
                }
                temp = zero;
                if (alphai[jc - 1] == zero) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        temp = max(temp, abs(vl[(jr - 1) + (jc - 1) * ldvl]));
                    }
                } else {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        temp = max(temp, abs(vl[(jr - 1) + (jc - 1) * ldvl]) + abs(vl[(jr - 1) + ((jc + 1) - 1) * ldvl]));
                    }
                }
                if (temp < smlnum) {
                    goto statement_50;
                }
                temp = one / temp;
                if (alphai[jc - 1] == zero) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        vl[(jr - 1) + (jc - 1) * ldvl] = vl[(jr - 1) + (jc - 1) * ldvl] * temp;
                    }
                } else {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        vl[(jr - 1) + (jc - 1) * ldvl] = vl[(jr - 1) + (jc - 1) * ldvl] * temp;
                        vl[(jr - 1) + ((jc + 1) - 1) * ldvl] = vl[(jr - 1) + ((jc + 1) - 1) * ldvl] * temp;
                    }
                }
            statement_50:;
            }
        }
        if (ilvr) {
            Rggbak("P", "R", n, ilo, ihi, &work[ileft - 1], &work[iright - 1], n, vr, ldvr, ierr);
            for (jc = 1; jc <= n; jc = jc + 1) {
                if (alphai[jc - 1] < zero) {
                    goto statement_100;
                }
                temp = zero;
                if (alphai[jc - 1] == zero) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        temp = max(temp, abs(vr[(jr - 1) + (jc - 1) * ldvr]));
                    }
                } else {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        temp = max(temp, abs(vr[(jr - 1) + (jc - 1) * ldvr]) + abs(vr[(jr - 1) + ((jc + 1) - 1) * ldvr]));
                    }
                }
                if (temp < smlnum) {
                    goto statement_100;
                }
                temp = one / temp;
                if (alphai[jc - 1] == zero) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        vr[(jr - 1) + (jc - 1) * ldvr] = vr[(jr - 1) + (jc - 1) * ldvr] * temp;
                    }
                } else {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        vr[(jr - 1) + (jc - 1) * ldvr] = vr[(jr - 1) + (jc - 1) * ldvr] * temp;
                        vr[(jr - 1) + ((jc + 1) - 1) * ldvr] = vr[(jr - 1) + ((jc + 1) - 1) * ldvr] * temp;
                    }
                }
            statement_100:;
            }
        }
        //
        //        End of eigenvector calculation
        //
    }
//
//     Undo scaling if necessary
//
statement_110:
    //
    if (ilascl) {
        Rlascl("G", 0, 0, anrmto, anrm, n, 1, alphar, n, ierr);
        Rlascl("G", 0, 0, anrmto, anrm, n, 1, alphai, n, ierr);
    }
    //
    if (ilbscl) {
        Rlascl("G", 0, 0, bnrmto, bnrm, n, 1, beta, n, ierr);
    }
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Rggev
    //
}
