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

void Cgeev(const char *jobvl, const char *jobvr, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *w, COMPLEX *vl, INTEGER const ldvl, COMPLEX *vr, INTEGER const ldvr, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER &info) {
    bool lquery = false;
    bool wantvl = false;
    bool wantvr = false;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    bool select[1];
    INTEGER nout = 0;
    INTEGER ierr = 0;
    INTEGER lwork_trevc = 0;
    INTEGER hswork = 0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL dum[1];
    REAL anrm = 0.0;
    bool scalea = false;
    const REAL zero = 0.0;
    REAL cscale = 0.0;
    INTEGER ibal = 0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    INTEGER itau = 0;
    INTEGER iwrk = 0;
    char side;
    INTEGER irwork = 0;
    INTEGER i = 0;
    REAL scl = 0.0;
    INTEGER k = 0;
    COMPLEX tmp = 0.0;
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
    lquery = (lwork == -1);
    wantvl = Mlsame(jobvl, "V");
    wantvr = Mlsame(jobvr, "V");
    if ((!wantvl) && (!Mlsame(jobvl, "N"))) {
        info = -1;
    } else if ((!wantvr) && (!Mlsame(jobvr, "N"))) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldvl < 1 || (wantvl && ldvl < n)) {
        info = -8;
    } else if (ldvr < 1 || (wantvr && ldvr < n)) {
        info = -10;
    }
    //
    //     Compute workspace
    //      (Note: Comments in the code beginning "Workspace:" describe the
    //       minimal amount of workspace needed at that point in the code,
    //       as well as the preferred amount for good performance.
    //       CWorkspace refers to complex workspace, and RWorkspace to real
    //       workspace. NB refers to the optimal block size for the
    //       immediately following subroutine, as returned by iMlaenv.
    //       HSWORK refers to the workspace preferred by Chseqr, as
    //       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
    //       the worst case.)
    //
    if (info == 0) {
        if (n == 0) {
            minwrk = 1;
            maxwrk = 1;
        } else {
            maxwrk = n + n * iMlaenv(1, "Cgehrd", " ", n, 1, n, 0);
            minwrk = 2 * n;
            if (wantvl) {
                maxwrk = max(maxwrk, n + (n - 1) * iMlaenv(1, "Cunghr", " ", n, 1, n, -1));
                Ctrevc3("L", "B", select, n, a, lda, vl, ldvl, vr, ldvr, n, nout, work, -1, rwork, -1, ierr);
                lwork_trevc = castINTEGER(work[1 - 1].real());
                maxwrk = max(maxwrk, n + lwork_trevc);
                Chseqr("S", "V", n, 1, n, a, lda, w, vl, ldvl, work, -1, info);
            } else if (wantvr) {
                maxwrk = max(maxwrk, n + (n - 1) * iMlaenv(1, "Cunghr", " ", n, 1, n, -1));
                Ctrevc3("R", "B", select, n, a, lda, vl, ldvl, vr, ldvr, n, nout, work, -1, rwork, -1, ierr);
                lwork_trevc = castINTEGER(work[1 - 1].real());
                maxwrk = max(maxwrk, n + lwork_trevc);
                Chseqr("S", "V", n, 1, n, a, lda, w, vr, ldvr, work, -1, info);
            } else {
                Chseqr("E", "N", n, 1, n, a, lda, w, vr, ldvr, work, -1, info);
            }
            hswork = castINTEGER(work[1 - 1].real());
            maxwrk = max({maxwrk, hswork, minwrk});
        }
        work[1 - 1] = maxwrk;
        //
        if (lwork < minwrk && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cgeev ", -info);
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
    anrm = Clange("M", n, n, a, lda, dum);
    scalea = false;
    if (anrm > zero && anrm < smlnum) {
        scalea = true;
        cscale = smlnum;
    } else if (anrm > bignum) {
        scalea = true;
        cscale = bignum;
    }
    if (scalea) {
        Clascl("G", 0, 0, anrm, cscale, n, n, a, lda, ierr);
    }
    //
    //     Balance the matrix
    //     (CWorkspace: none)
    //     (RWorkspace: need N)
    //
    ibal = 1;
    Cgebal("B", n, a, lda, ilo, ihi, &rwork[ibal - 1], ierr);
    //
    //     Reduce to upper Hessenberg form
    //     (CWorkspace: need 2*N, prefer N+N*NB)
    //     (RWorkspace: none)
    //
    itau = 1;
    iwrk = itau + n;
    Cgehrd(n, ilo, ihi, a, lda, &work[itau - 1], &work[iwrk - 1], lwork - iwrk + 1, ierr);
    //
    if (wantvl) {
        //
        //        Want left eigenvectors
        //        Copy Householder vectors to VL
        //
        side = 'L';
        Clacpy("L", n, n, a, lda, vl, ldvl);
        //
        //        Generate unitary matrix in VL
        //        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
        //        (RWorkspace: none)
        //
        Cunghr(n, ilo, ihi, vl, ldvl, &work[itau - 1], &work[iwrk - 1], lwork - iwrk + 1, ierr);
        //
        //        Perform QR iteration, accumulating Schur vectors in VL
        //        (CWorkspace: need 1, prefer HSWORK (see comments) )
        //        (RWorkspace: none)
        //
        iwrk = itau;
        Chseqr("S", "V", n, ilo, ihi, a, lda, w, vl, ldvl, &work[iwrk - 1], lwork - iwrk + 1, info);
        //
        if (wantvr) {
            //
            //           Want left and right eigenvectors
            //           Copy Schur vectors to VR
            //
            side = 'B';
            Clacpy("F", n, n, vl, ldvl, vr, ldvr);
        }
        //
    } else if (wantvr) {
        //
        //        Want right eigenvectors
        //        Copy Householder vectors to VR
        //
        side = 'R';
        Clacpy("L", n, n, a, lda, vr, ldvr);
        //
        //        Generate unitary matrix in VR
        //        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
        //        (RWorkspace: none)
        //
        Cunghr(n, ilo, ihi, vr, ldvr, &work[itau - 1], &work[iwrk - 1], lwork - iwrk + 1, ierr);
        //
        //        Perform QR iteration, accumulating Schur vectors in VR
        //        (CWorkspace: need 1, prefer HSWORK (see comments) )
        //        (RWorkspace: none)
        //
        iwrk = itau;
        Chseqr("S", "V", n, ilo, ihi, a, lda, w, vr, ldvr, &work[iwrk - 1], lwork - iwrk + 1, info);
        //
    } else {
        //
        //        Compute eigenvalues only
        //        (CWorkspace: need 1, prefer HSWORK (see comments) )
        //        (RWorkspace: none)
        //
        iwrk = itau;
        Chseqr("E", "N", n, ilo, ihi, a, lda, w, vr, ldvr, &work[iwrk - 1], lwork - iwrk + 1, info);
    }
    //
    //     If INFO .NE. 0 from Chseqr, then quit
    //
    if (info != 0) {
        goto statement_50;
    }
    //
    if (wantvl || wantvr) {
        //
        //        Compute left and/or right eigenvectors
        //        (CWorkspace: need 2*N, prefer N + 2*N*NB)
        //        (RWorkspace: need 2*N)
        //
        irwork = ibal + n;
        Ctrevc3(&side, "B", select, n, a, lda, vl, ldvl, vr, ldvr, n, nout, &work[iwrk - 1], lwork - iwrk + 1, &rwork[irwork - 1], n, ierr);
    }
    //
    if (wantvl) {
        //
        //        Undo balancing of left eigenvectors
        //        (CWorkspace: none)
        //        (RWorkspace: need N)
        //
        Cgebak("B", "L", n, ilo, ihi, &rwork[ibal - 1], n, vl, ldvl, ierr);
        //
        //        Normalize left eigenvectors and make largest component real
        //
        for (i = 1; i <= n; i = i + 1) {
            scl = one / RCnrm2(n, &vl[(i - 1) * ldvl], 1);
            CRscal(n, scl, &vl[(i - 1) * ldvl], 1);
            for (k = 1; k <= n; k = k + 1) {
                rwork[(irwork + k - 1) - 1] = pow2(vl[(k - 1) + (i - 1) * ldvl].real()) + pow2(vl[(k - 1) + (i - 1) * ldvl].imag());
            }
            k = iRamax(n, &rwork[irwork - 1], 1);
            tmp = conj(vl[(k - 1) + (i - 1) * ldvl]) / sqrt(rwork[(irwork + k - 1) - 1]);
            Cscal(n, tmp, &vl[(i - 1) * ldvl], 1);
            vl[(k - 1) + (i - 1) * ldvl] = COMPLEX(vl[(k - 1) + (i - 1) * ldvl].real(), zero);
        }
    }
    //
    if (wantvr) {
        //
        //        Undo balancing of right eigenvectors
        //        (CWorkspace: none)
        //        (RWorkspace: need N)
        //
        Cgebak("B", "R", n, ilo, ihi, &rwork[ibal - 1], n, vr, ldvr, ierr);
        //
        //        Normalize right eigenvectors and make largest component real
        //
        for (i = 1; i <= n; i = i + 1) {
            scl = one / RCnrm2(n, &vr[(i - 1) * ldvr], 1);
            CRscal(n, scl, &vr[(i - 1) * ldvr], 1);
            for (k = 1; k <= n; k = k + 1) {
                rwork[(irwork + k - 1) - 1] = pow2(vr[(k - 1) + (i - 1) * ldvr].real()) + pow2(vr[(k - 1) + (i - 1) * ldvr].real());
            }
            k = iRamax(n, &rwork[irwork - 1], 1);
            tmp = conj(vr[(k - 1) + (i - 1) * ldvr]) / sqrt(rwork[(irwork + k - 1) - 1]);
            Cscal(n, tmp, &vr[(i - 1) * ldvr], 1);
            vr[(k - 1) + (i - 1) * ldvr] = COMPLEX(vr[(k - 1) + (i - 1) * ldvr].real(), zero);
        }
    }
//
//     Undo scaling if necessary
//
statement_50:
    if (scalea) {
        Clascl("G", 0, 0, cscale, anrm, n - info, 1, &w[(info + 1) - 1], max(n - info, (INTEGER)1), ierr);
        if (info > 0) {
            Clascl("G", 0, 0, cscale, anrm, ilo - 1, 1, w, n, ierr);
        }
    }
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Cgeev
    //
}
