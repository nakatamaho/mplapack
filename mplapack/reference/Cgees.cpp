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

void Cgees(const char *jobvs, const char *sort, UNHANDLED_function_poINTEGERer select, INTEGER const &n, COMPLEX *a, INTEGER const &lda, INTEGER &sdim, COMPLEX *w, COMPLEX *vs, INTEGER const &ldvs, COMPLEX *work, INTEGER const &lwork, REAL *rwork, arr_ref<bool> bwork, INTEGER &info) {
    //
    //  -- LAPACK driver routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //     .. Function Arguments ..
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
    bool lquery = (lwork == -1);
    bool wantvs = Mlsame(jobvs, "V");
    bool wantst = Mlsame(sort, "S");
    if ((!wantvs) && (!Mlsame(jobvs, "N"))) {
        info = -1;
    } else if ((!wantst) && (!Mlsame(sort, "N"))) {
        info = -2;
    } else if (n < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldvs < 1 || (wantvs && ldvs < n)) {
        info = -10;
    }
    //
    //     Compute workspace
    //      (Note: Comments in the code beginning "Workspace:" describe the
    //       minimal amount of workspace needed at that poINTEGER in the code,
    //       as well as the preferred amount for good performance.
    //       CWorkspace refers to complex workspace, and RWorkspace to real
    //       workspace. NB refers to the optimal block size for the
    //       immediately following subroutine, as returned by iMlaenv.
    //       HSWORK refers to the workspace preferred by Chseqr, as
    //       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
    //       the worst case.)
    //
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER ieval = 0;
    INTEGER hswork = 0;
    if (info == 0) {
        if (n == 0) {
            minwrk = 1;
            maxwrk = 1;
        } else {
            maxwrk = n + n * iMlaenv[("Cgehrd" - 1) * ldiMlaenv];
            minwrk = 2 * n;
            //
            Chseqr("S", jobvs, n, 1, n, a, lda, w, vs, ldvs, work, -1, ieval);
            hswork = work[1 - 1];
            //
            if (!wantvs) {
                maxwrk = max(maxwrk, hswork);
            } else {
                maxwrk = max(maxwrk, n + (n - 1) * iMlaenv[("Cunghr" - 1) * ldiMlaenv]);
                maxwrk = max(maxwrk, hswork);
            }
        }
        work[1 - 1] = maxwrk;
        //
        if (lwork < minwrk && !lquery) {
            info = -12;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cgees ", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        sdim = 0;
        return;
    }
    //
    //     Get machine constants
    //
    REAL eps = dlamch("P");
    REAL smlnum = dlamch("S");
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = one / smlnum;
    //
    //     Scale A if max element outside range [SMLNUM,BIGNUM]
    //
    arr_1d<1, REAL> dum(fill0);
    REAL anrm = Clange[("M" - 1) + (n - 1) * ldClange];
    bool scalea = false;
    const REAL zero = 0.0;
    REAL cscale = 0.0;
    if (anrm > zero && anrm < smlnum) {
        scalea = true;
        cscale = smlnum;
    } else if (anrm > bignum) {
        scalea = true;
        cscale = bignum;
    }
    INTEGER ierr = 0;
    if (scalea) {
        Clascl("G", 0, 0, anrm, cscale, n, n, a, lda, ierr);
    }
    //
    //     Permute the matrix to make it more nearly triangular
    //     (CWorkspace: none)
    //     (RWorkspace: need N)
    //
    INTEGER ibal = 1;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    Cgebal("P", n, a, lda, ilo, ihi, rwork[ibal - 1], ierr);
    //
    //     Reduce to upper Hessenberg form
    //     (CWorkspace: need 2*N, prefer N+N*NB)
    //     (RWorkspace: none)
    //
    INTEGER itau = 1;
    INTEGER iwrk = n + itau;
    Cgehrd(n, ilo, ihi, a, lda, work[itau - 1], work[iwrk - 1], lwork - iwrk + 1, ierr);
    //
    if (wantvs) {
        //
        //        Copy Householder vectors to VS
        //
        Clacpy("L", n, n, a, lda, vs, ldvs);
        //
        //        Generate unitary matrix in VS
        //        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
        //        (RWorkspace: none)
        //
        Cunghr(n, ilo, ihi, vs, ldvs, work[itau - 1], work[iwrk - 1], lwork - iwrk + 1, ierr);
    }
    //
    sdim = 0;
    //
    //     Perform QR iteration, accumulating Schur vectors in VS if desired
    //     (CWorkspace: need 1, prefer HSWORK (see comments) )
    //     (RWorkspace: none)
    //
    iwrk = itau;
    Chseqr("S", jobvs, n, ilo, ihi, a, lda, w, vs, ldvs, work[iwrk - 1], lwork - iwrk + 1, ieval);
    if (ieval > 0) {
        info = ieval;
    }
    //
    //     Sort eigenvalues if desired
    //
    INTEGER i = 0;
    REAL s = 0.0;
    REAL sep = 0.0;
    INTEGER icond = 0;
    if (wantst && info == 0) {
        if (scalea) {
            Clascl("G", 0, 0, cscale, anrm, n, 1, w, n, ierr);
        }
        for (i = 1; i <= n; i = i + 1) {
            bwork[i - 1] = select[w[i - 1] - 1];
        }
        //
        //        Reorder eigenvalues and transform Schur vectors
        //        (CWorkspace: none)
        //        (RWorkspace: none)
        //
        Ctrsen("N", jobvs, bwork, n, a, lda, vs, ldvs, w, sdim, s, sep, work[iwrk - 1], lwork - iwrk + 1, icond);
    }
    //
    if (wantvs) {
        //
        //        Undo balancing
        //        (CWorkspace: none)
        //        (RWorkspace: need N)
        //
        Cgebak("P", "R", n, ilo, ihi, rwork[ibal - 1], n, vs, ldvs, ierr);
    }
    //
    if (scalea) {
        //
        //        Undo scaling for the Schur form of A
        //
        Clascl("U", 0, 0, cscale, anrm, n, n, a, lda, ierr);
        Ccopy(n, a, lda + 1, w, 1);
    }
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Cgees
    //
}
