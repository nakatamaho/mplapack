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

void Rgeevx(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, INTEGER const n, REAL *a, INTEGER const lda, REAL *wr, REAL *wi, REAL *vl, INTEGER const ldvl, REAL *vr, INTEGER const ldvr, INTEGER ilo, INTEGER ihi, REAL *scale, REAL &abnrm, REAL *rconde, REAL *rcondv, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER &info) {
    bool lquery = false;
    bool wantvl = false;
    bool wantvr = false;
    bool wntsnn = false;
    bool wntsne = false;
    bool wntsnv = false;
    bool wntsnb = false;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    bool select;
    INTEGER nout = 0;
    INTEGER ierr = 0;
    INTEGER lwork_trevc = 0;
    INTEGER hswork = 0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    INTEGER icond = 0;
    REAL dum;
    REAL anrm = 0.0;
    bool scalea = false;
    const REAL zero = 0.0;
    REAL cscale = 0.0;
    INTEGER itau = 0;
    INTEGER iwrk = 0;
    char side;
    char job;
    INTEGER i = 0;
    REAL scl = 0.0;
    INTEGER k = 0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    REAL r = 0.0;
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
    wntsnn = Mlsame(sense, "N");
    wntsne = Mlsame(sense, "E");
    wntsnv = Mlsame(sense, "V");
    wntsnb = Mlsame(sense, "B");
    if (!(Mlsame(balanc, "N") || Mlsame(balanc, "S") || Mlsame(balanc, "P") || Mlsame(balanc, "B"))) {
        info = -1;
    } else if ((!wantvl) && (!Mlsame(jobvl, "N"))) {
        info = -2;
    } else if ((!wantvr) && (!Mlsame(jobvr, "N"))) {
        info = -3;
    } else if (!(wntsnn || wntsne || wntsnb || wntsnv) || ((wntsne || wntsnb) && !(wantvl && wantvr))) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (lda < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldvl < 1 || (wantvl && ldvl < n)) {
        info = -11;
    } else if (ldvr < 1 || (wantvr && ldvr < n)) {
        info = -13;
    }
    //
    //     Compute workspace
    //      (Note: Comments in the code beginning "Workspace:" describe the
    //       minimal amount of workspace needed at that point in the code,
    //       as well as the preferred amount for good performance.
    //       NB refers to the optimal block size for the immediately
    //       following subroutine, as returned by iMlaenv.
    //       HSWORK refers to the workspace preferred by Rhseqr, as
    //       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
    //       the worst case.)
    //
    if (info == 0) {
        if (n == 0) {
            minwrk = 1;
            maxwrk = 1;
        } else {
            maxwrk = n + n * iMlaenv(1, "Rgehrd", " ", n, 1, n, 0);
            //
            if (wantvl) {
                Rtrevc3("L", "B", &select, n, a, lda, vl, ldvl, vr, ldvr, n, nout, work, -1, ierr);
                lwork_trevc = castINTEGER(work[1 - 1]);
                maxwrk = max(maxwrk, n + lwork_trevc);
                Rhseqr("S", "V", n, 1, n, a, lda, wr, wi, vl, ldvl, work, -1, info);
            } else if (wantvr) {
                Rtrevc3("R", "B", &select, n, a, lda, vl, ldvl, vr, ldvr, n, nout, work, -1, ierr);
                lwork_trevc = castINTEGER(work[1 - 1]);
                maxwrk = max(maxwrk, n + lwork_trevc);
                Rhseqr("S", "V", n, 1, n, a, lda, wr, wi, vr, ldvr, work, -1, info);
            } else {
                if (wntsnn) {
                    Rhseqr("E", "N", n, 1, n, a, lda, wr, wi, vr, ldvr, work, -1, info);
                } else {
                    Rhseqr("S", "N", n, 1, n, a, lda, wr, wi, vr, ldvr, work, -1, info);
                }
            }
            hswork = castINTEGER(work[1 - 1]);
            //
            if ((!wantvl) && (!wantvr)) {
                minwrk = 2 * n;
                if (!wntsnn) {
                    minwrk = max(minwrk, n * n + 6 * n);
                }
                maxwrk = max(maxwrk, hswork);
                if (!wntsnn) {
                    maxwrk = max(maxwrk, n * n + 6 * n);
                }
            } else {
                minwrk = 3 * n;
                if ((!wntsnn) && (!wntsne)) {
                    minwrk = max(minwrk, n * n + 6 * n);
                }
                maxwrk = max(maxwrk, hswork);
                maxwrk = max(maxwrk, n + (n - 1) * iMlaenv(1, "Rorghr", " ", n, 1, n, -1));
                if ((!wntsnn) && (!wntsne)) {
                    maxwrk = max(maxwrk, n * n + 6 * n);
                }
                maxwrk = max(maxwrk, 3 * n);
            }
            maxwrk = max(maxwrk, minwrk);
        }
        work[1 - 1] = maxwrk;
        //
        if (lwork < minwrk && !lquery) {
            info = -21;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rgeevx", -info);
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
    icond = 0;
    anrm = Rlange("M", n, n, a, lda, &dum);
    scalea = false;
    if (anrm > zero && anrm < smlnum) {
        scalea = true;
        cscale = smlnum;
    } else if (anrm > bignum) {
        scalea = true;
        cscale = bignum;
    }
    if (scalea) {
        Rlascl("G", 0, 0, anrm, cscale, n, n, a, lda, ierr);
    }
    //
    //     Balance the matrix and compute ABNRM
    //
    Rgebal(balanc, n, a, lda, ilo, ihi, scale, ierr);
    abnrm = Rlange("1", n, n, a, lda, &dum);
    if (scalea) {
        dum[1 - 1] = *abnrm;
        Rlascl("G", 0, 0, cscale, anrm, 1, 1, &dum, 1, ierr);
        *abnrm = dum[1 - 1];
    }
    //
    //     Reduce to upper Hessenberg form
    //     (Workspace: need 2*N, prefer N+N*NB)
    //
    itau = 1;
    iwrk = itau + n;
    Rgehrd(n, ilo, ihi, a, lda, &work[itau - 1], &work[iwrk - 1], lwork - iwrk + 1, ierr);
    //
    if (wantvl) {
        //
        //        Want left eigenvectors
        //        Copy Householder vectors to VL
        //
        side = 'L';
        Rlacpy("L", n, n, a, lda, vl, ldvl);
        //
        //        Generate orthogonal matrix in VL
        //        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
        //
        Rorghr(n, ilo, ihi, vl, ldvl, &work[itau - 1], &work[iwrk - 1], lwork - iwrk + 1, ierr);
        //
        //        Perform QR iteration, accumulating Schur vectors in VL
        //        (Workspace: need 1, prefer HSWORK (see comments) )
        //
        iwrk = itau;
        Rhseqr("S", "V", n, ilo, ihi, a, lda, wr, wi, vl, ldvl, &work[iwrk - 1], lwork - iwrk + 1, info);
        //
        if (wantvr) {
            //
            //           Want left and right eigenvectors
            //           Copy Schur vectors to VR
            //
            side = 'B';
            Rlacpy("F", n, n, vl, ldvl, vr, ldvr);
        }
        //
    } else if (wantvr) {
        //
        //        Want right eigenvectors
        //        Copy Householder vectors to VR
        //
        side = 'R';
        Rlacpy("L", n, n, a, lda, vr, ldvr);
        //
        //        Generate orthogonal matrix in VR
        //        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
        //
        Rorghr(n, ilo, ihi, vr, ldvr, &work[itau - 1], &work[iwrk - 1], lwork - iwrk + 1, ierr);
        //
        //        Perform QR iteration, accumulating Schur vectors in VR
        //        (Workspace: need 1, prefer HSWORK (see comments) )
        //
        iwrk = itau;
        Rhseqr("S", "V", n, ilo, ihi, a, lda, wr, wi, vr, ldvr, &work[iwrk - 1], lwork - iwrk + 1, info);
        //
    } else {
        //
        //        Compute eigenvalues only
        //        If condition numbers desired, compute Schur form
        //
        if (wntsnn) {
            job = 'E';
        } else {
            job = 'S';
        }
        //
        //        (Workspace: need 1, prefer HSWORK (see comments) )
        //
        iwrk = itau;
        Rhseqr(&job, "N", n, ilo, ihi, a, lda, wr, wi, vr, ldvr, &work[iwrk - 1], lwork - iwrk + 1, info);
    }
    //
    //     If INFO .NE. 0 from Rhseqr, then quit
    //
    if (info != 0) {
        goto statement_50;
    }
    //
    if (wantvl || wantvr) {
        //
        //        Compute left and/or right eigenvectors
        //        (Workspace: need 3*N, prefer N + 2*N*NB)
        //
        Rtrevc3(&side, "B", &select, n, a, lda, vl, ldvl, vr, ldvr, n, nout, &work[iwrk - 1], lwork - iwrk + 1, ierr);
    }
    //
    //     Compute condition numbers if desired
    //     (Workspace: need N*N+6*N unless SENSE = 'E')
    //
    if (!wntsnn) {
        Rtrsna(sense, "A", &select, n, a, lda, vl, ldvl, vr, ldvr, rconde, rcondv, n, nout, &work[iwrk - 1], n, iwork, icond);
    }
    //
    if (wantvl) {
        //
        //        Undo balancing of left eigenvectors
        //
        Rgebak(balanc, "L", n, ilo, ihi, scale, n, vl, ldvl, ierr);
        //
        //        Normalize left eigenvectors and make largest component real
        //
        for (i = 1; i <= n; i = i + 1) {
            if (wi[i - 1] == zero) {
                scl = one / Rnrm2(n, &vl[(i - 1) * ldvl], 1);
                Rscal(n, scl, &vl[(i - 1) * ldvl], 1);
            } else if (wi[i - 1] > zero) {
                scl = one / Rlapy2(Rnrm2(n, &vl[(i - 1) * ldvl], 1), Rnrm2(n, &vl[((i + 1) - 1) * ldvl], 1));
                Rscal(n, scl, &vl[(i - 1) * ldvl], 1);
                Rscal(n, scl, &vl[((i + 1) - 1) * ldvl], 1);
                for (k = 1; k <= n; k = k + 1) {
                    work[k - 1] = pow2(vl[(k - 1) + (i - 1) * ldvl]) + pow2(vl[(k - 1) + ((i + 1) - 1) * ldvl]);
                }
                k = iRamax(n, work, 1);
                Rlartg(vl[(k - 1) + (i - 1) * ldvl], vl[(k - 1) + ((i + 1) - 1) * ldvl], cs, sn, r);
                Rrot(n, &vl[(i - 1) * ldvl], 1, &vl[((i + 1) - 1) * ldvl], 1, cs, sn);
                vl[(k - 1) + ((i + 1) - 1) * ldvl] = zero;
            }
        }
    }
    //
    if (wantvr) {
        //
        //        Undo balancing of right eigenvectors
        //
        Rgebak(balanc, "R", n, ilo, ihi, scale, n, vr, ldvr, ierr);
        //
        //        Normalize right eigenvectors and make largest component real
        //
        for (i = 1; i <= n; i = i + 1) {
            if (wi[i - 1] == zero) {
                scl = one / Rnrm2(n, &vr[(i - 1) * ldvr], 1);
                Rscal(n, scl, &vr[(i - 1) * ldvr], 1);
            } else if (wi[i - 1] > zero) {
                scl = one / Rlapy2(Rnrm2(n, &vr[(i - 1) * ldvr], 1), Rnrm2(n, &vr[((i + 1) - 1) * ldvr], 1));
                Rscal(n, scl, &vr[(i - 1) * ldvr], 1);
                Rscal(n, scl, &vr[((i + 1) - 1) * ldvr], 1);
                for (k = 1; k <= n; k = k + 1) {
                    work[k - 1] = pow2(vr[(k - 1) + (i - 1) * ldvr]) + pow2(vr[(k - 1) + ((i + 1) - 1) * ldvr]);
                }
                k = iRamax(n, work, 1);
                Rlartg(vr[(k - 1) + (i - 1) * ldvr], vr[(k - 1) + ((i + 1) - 1) * ldvr], cs, sn, r);
                Rrot(n, &vr[(i - 1) * ldvr], 1, &vr[((i + 1) - 1) * ldvr], 1, cs, sn);
                vr[(k - 1) + ((i + 1) - 1) * ldvr] = zero;
            }
        }
    }
//
//     Undo scaling if necessary
//
statement_50:
    if (scalea) {
        Rlascl("G", 0, 0, cscale, anrm, n - info, 1, &wr[(info + 1) - 1], max(n - info, 1), ierr);
        Rlascl("G", 0, 0, cscale, anrm, n - info, 1, &wi[(info + 1) - 1], max(n - info, 1), ierr);
        if (info == 0) {
            if ((wntsnv || wntsnb) && icond == 0) {
                Rlascl("G", 0, 0, cscale, anrm, n, 1, rcondv, n, ierr);
            }
        } else {
            Rlascl("G", 0, 0, cscale, anrm, ilo - 1, 1, wr, n, ierr);
            Rlascl("G", 0, 0, cscale, anrm, ilo - 1, 1, wi, n, ierr);
        }
    }
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Rgeevx
    //
}
