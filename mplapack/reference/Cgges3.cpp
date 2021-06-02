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

void Cgges3(const char *jobvsl, const char *jobvsr, const char *sort, bool (*selctg)(COMPLEX, COMPLEX), INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, INTEGER &sdim, COMPLEX *alpha, COMPLEX *beta, COMPLEX *vsl, INTEGER const ldvsl, COMPLEX *vsr, INTEGER const ldvsr, COMPLEX *work, INTEGER const lwork, REAL *rwork, bool *bwork, INTEGER &info) {
    INTEGER ijobvl = 0;
    bool ilvsl = false;
    INTEGER ijobvr = 0;
    bool ilvsr = false;
    bool wantst = false;
    bool lquery = false;
    INTEGER ierr = 0;
    INTEGER lwkopt = 0;
    REAL pvsl = 0.0;
    REAL pvsr = 0.0;
    REAL dif[2];
    INTEGER idum[1];
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL anrm = 0.0;
    bool ilascl = false;
    const REAL zero = 0.0;
    REAL anrmto = 0.0;
    REAL bnrm = 0.0;
    bool ilbscl = false;
    REAL bnrmto = 0.0;
    INTEGER ileft = 0;
    INTEGER iright = 0;
    INTEGER irwrk = 0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    INTEGER irows = 0;
    INTEGER icols = 0;
    INTEGER itau = 0;
    INTEGER iwrk = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER i = 0;
    bool lastsl = false;
    bool cursl = false;
    //
    //     Decode the input arguments
    //
    if (Mlsame(jobvsl, "N")) {
        ijobvl = 1;
        ilvsl = false;
    } else if (Mlsame(jobvsl, "V")) {
        ijobvl = 2;
        ilvsl = true;
    } else {
        ijobvl = -1;
        ilvsl = false;
    }
    //
    if (Mlsame(jobvsr, "N")) {
        ijobvr = 1;
        ilvsr = false;
    } else if (Mlsame(jobvsr, "V")) {
        ijobvr = 2;
        ilvsr = true;
    } else {
        ijobvr = -1;
        ilvsr = false;
    }
    //
    wantst = Mlsame(sort, "S");
    //
    //     Test the input arguments
    //
    info = 0;
    lquery = (lwork == -1);
    if (ijobvl <= 0) {
        info = -1;
    } else if (ijobvr <= 0) {
        info = -2;
    } else if ((!wantst) && (!Mlsame(sort, "N"))) {
        info = -3;
    } else if (n < 0) {
        info = -5;
    } else if (lda < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -9;
    } else if (ldvsl < 1 || (ilvsl && ldvsl < n)) {
        info = -14;
    } else if (ldvsr < 1 || (ilvsr && ldvsr < n)) {
        info = -16;
    } else if (lwork < max((INTEGER)1, 2 * n) && !lquery) {
        info = -18;
    }
    //
    //     Compute workspace
    //
    if (info == 0) {
        Cgeqrf(n, n, b, ldb, work, work, -1, ierr);
        lwkopt = max((INTEGER)1, n + castINTEGER(work[1 - 1].real()));
        Cunmqr("L", "C", n, n, n, b, ldb, work, a, lda, work, -1, ierr);
        lwkopt = max(lwkopt, n + castINTEGER(work[1 - 1].real()));
        if (ilvsl) {
            Cungqr(n, n, n, vsl, ldvsl, work, work, -1, ierr);
            lwkopt = max(lwkopt, n + castINTEGER(work[1 - 1].real()));
        }
        Cgghd3(jobvsl, jobvsr, n, 1, n, a, lda, b, ldb, vsl, ldvsl, vsr, ldvsr, work, -1, ierr);
        lwkopt = max(lwkopt, n + castINTEGER(work[1 - 1].real()));
        Chgeqz("S", jobvsl, jobvsr, n, 1, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, -1, rwork, ierr);
        lwkopt = max(lwkopt, castINTEGER(work[1 - 1].real()));
        if (wantst) {
            Ctgsen(0, ilvsl, ilvsr, bwork, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, sdim, pvsl, pvsr, dif, work, -1, idum, 1, ierr);
            lwkopt = max(lwkopt, castINTEGER(work[1 - 1].real()));
        }
        work[1 - 1] = COMPLEX(lwkopt);
    }
    //
    if (info != 0) {
        Mxerbla("Cgges3", -info);
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
    eps = Rlamch("P");
    smlnum = Rlamch("S");
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = one / smlnum;
    //
    //     Scale A if max element outside range [SMLNUM,BIGNUM]
    //
    anrm = Clange("M", n, n, a, lda, rwork);
    ilascl = false;
    if (anrm > zero && anrm < smlnum) {
        anrmto = smlnum;
        ilascl = true;
    } else if (anrm > bignum) {
        anrmto = bignum;
        ilascl = true;
    }
    //
    if (ilascl) {
        Clascl("G", 0, 0, anrm, anrmto, n, n, a, lda, ierr);
    }
    //
    //     Scale B if max element outside range [SMLNUM,BIGNUM]
    //
    bnrm = Clange("M", n, n, b, ldb, rwork);
    ilbscl = false;
    if (bnrm > zero && bnrm < smlnum) {
        bnrmto = smlnum;
        ilbscl = true;
    } else if (bnrm > bignum) {
        bnrmto = bignum;
        ilbscl = true;
    }
    //
    if (ilbscl) {
        Clascl("G", 0, 0, bnrm, bnrmto, n, n, b, ldb, ierr);
    }
    //
    //     Permute the matrix to make it more nearly triangular
    //
    ileft = 1;
    iright = n + 1;
    irwrk = iright + n;
    Cggbal("P", n, a, lda, b, ldb, ilo, ihi, &rwork[ileft - 1], &rwork[iright - 1], &rwork[irwrk - 1], ierr);
    //
    //     Reduce B to triangular form (QR decomposition of B)
    //
    irows = ihi + 1 - ilo;
    icols = n + 1 - ilo;
    itau = 1;
    iwrk = itau + irows;
    Cgeqrf(irows, icols, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, &work[itau - 1], &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    //     Apply the orthogonal transformation to matrix A
    //
    Cunmqr("L", "C", irows, icols, irows, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, &work[itau - 1], &a[(ilo - 1) + (ilo - 1) * lda], lda, &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    //     Initialize VSL
    //
    if (ilvsl) {
        Claset("Full", n, n, czero, cone, vsl, ldvsl);
        if (irows > 1) {
            Clacpy("L", irows - 1, irows - 1, &b[((ilo + 1) - 1) + (ilo - 1) * ldb], ldb, &vsl[((ilo + 1) - 1) + (ilo - 1) * ldvsl], ldvsl);
        }
        Cungqr(irows, irows, irows, &vsl[(ilo - 1) + (ilo - 1) * ldvsl], ldvsl, &work[itau - 1], &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    }
    //
    //     Initialize VSR
    //
    if (ilvsr) {
        Claset("Full", n, n, czero, cone, vsr, ldvsr);
    }
    //
    //     Reduce to generalized Hessenberg form
    //
    Cgghd3(jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, vsl, ldvsl, vsr, ldvsr, &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    sdim = 0;
    //
    //     Perform QZ algorithm, computing Schur vectors if desired
    //
    iwrk = itau;
    Chgeqz("S", jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, &work[iwrk - 1], lwork + 1 - iwrk, &rwork[irwrk - 1], ierr);
    if (ierr != 0) {
        if (ierr > 0 && ierr <= n) {
            info = ierr;
        } else if (ierr > n && ierr <= 2 * n) {
            info = ierr - n;
        } else {
            info = n + 1;
        }
        goto statement_30;
    }
    //
    //     Sort eigenvalues ALPHA/BETA if desired
    //
    if (wantst) {
        //
        //        Undo scaling on eigenvalues before selecting
        //
        if (ilascl) {
            Clascl("G", 0, 0, anrm, anrmto, n, 1, alpha, n, ierr);
        }
        if (ilbscl) {
            Clascl("G", 0, 0, bnrm, bnrmto, n, 1, beta, n, ierr);
        }
        //
        //        Select eigenvalues
        //
        for (i = 1; i <= n; i = i + 1) {
            bwork[i - 1] = selctg(alpha[i - 1], beta[i - 1]);
        }
        //
        Ctgsen(0, ilvsl, ilvsr, bwork, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, sdim, pvsl, pvsr, dif, &work[iwrk - 1], lwork - iwrk + 1, idum, 1, ierr);
        if (ierr == 1) {
            info = n + 3;
        }
        //
    }
    //
    //     Apply back-permutation to VSL and VSR
    //
    if (ilvsl) {
        Cggbak("P", "L", n, ilo, ihi, &rwork[ileft - 1], &rwork[iright - 1], n, vsl, ldvsl, ierr);
    }
    if (ilvsr) {
        Cggbak("P", "R", n, ilo, ihi, &rwork[ileft - 1], &rwork[iright - 1], n, vsr, ldvsr, ierr);
    }
    //
    //     Undo scaling
    //
    if (ilascl) {
        Clascl("U", 0, 0, anrmto, anrm, n, n, a, lda, ierr);
        Clascl("G", 0, 0, anrmto, anrm, n, 1, alpha, n, ierr);
    }
    //
    if (ilbscl) {
        Clascl("U", 0, 0, bnrmto, bnrm, n, n, b, ldb, ierr);
        Clascl("G", 0, 0, bnrmto, bnrm, n, 1, beta, n, ierr);
    }
    //
    if (wantst) {
        //
        //        Check if reordering is correct
        //
        lastsl = true;
        sdim = 0;
        for (i = 1; i <= n; i = i + 1) {
            cursl = selctg(alpha[i - 1], beta[i - 1]);
            if (cursl) {
                sdim++;
            }
            if (cursl && !lastsl) {
                info = n + 2;
            }
            lastsl = cursl;
        }
        //
    }
//
statement_30:
    //
    work[1 - 1] = COMPLEX(lwkopt);
    //
    //     End of Cgges3
    //
}
