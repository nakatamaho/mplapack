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

void Cggesx(const char *jobvsl, const char *jobvsr, const char *sort, bool (*selctg)(COMPLEX, COMPLEX), const char *sense, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, INTEGER &sdim, COMPLEX *alpha, COMPLEX *beta, COMPLEX *vsl, INTEGER const ldvsl, COMPLEX *vsr, INTEGER const ldvsr, REAL *rconde, REAL *rcondv, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER *iwork, INTEGER const liwork, bool *bwork, INTEGER &info) {
    INTEGER ijobvl = 0;
    bool ilvsl = false;
    INTEGER ijobvr = 0;
    bool ilvsr = false;
    bool wantst = false;
    bool wantsn = false;
    bool wantse = false;
    bool wantsv = false;
    bool wantsb = false;
    bool lquery = false;
    INTEGER ijob = 0;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER lwrk = 0;
    INTEGER liwmin = 0;
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
    REAL pl = 0.0;
    REAL pr = 0.0;
    REAL dif[2];
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
    wantsn = Mlsame(sense, "N");
    wantse = Mlsame(sense, "E");
    wantsv = Mlsame(sense, "V");
    wantsb = Mlsame(sense, "B");
    lquery = (lwork == -1 || liwork == -1);
    if (wantsn) {
        ijob = 0;
    } else if (wantse) {
        ijob = 1;
    } else if (wantsv) {
        ijob = 2;
    } else if (wantsb) {
        ijob = 4;
    }
    //
    //     Test the input arguments
    //
    info = 0;
    if (ijobvl <= 0) {
        info = -1;
    } else if (ijobvr <= 0) {
        info = -2;
    } else if ((!wantst) && (!Mlsame(sort, "N"))) {
        info = -3;
    } else if (!(wantsn || wantse || wantsv || wantsb) || (!wantst && !wantsn)) {
        info = -5;
    } else if (n < 0) {
        info = -6;
    } else if (lda < max((INTEGER)1, n)) {
        info = -8;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -10;
    } else if (ldvsl < 1 || (ilvsl && ldvsl < n)) {
        info = -15;
    } else if (ldvsr < 1 || (ilvsr && ldvsr < n)) {
        info = -17;
    }
    //
    //     Compute workspace
    //      (Note: Comments in the code beginning "Workspace:" describe the
    //       minimal amount of workspace needed at that point in the code,
    //       as well as the preferred amount for good performance.
    //       NB refers to the optimal block size for the immediately
    //       following subroutine, as returned by iMlaenv.)
    //
    if (info == 0) {
        if (n > 0) {
            minwrk = 2 * n;
            maxwrk = n * (1 + iMlaenv(1, "Cgeqrf", " ", n, 1, n, 0));
            maxwrk = max({maxwrk, n * (1 + iMlaenv(1, "Cunmqr", " ", n, 1, n, -1))});
            if (ilvsl) {
                maxwrk = max({maxwrk, n * (1 + iMlaenv(1, "Cungqr", " ", n, 1, n, -1))});
            }
            lwrk = maxwrk;
            if (ijob >= 1) {
                lwrk = max(lwrk, n * n / 2);
            }
        } else {
            minwrk = 1;
            maxwrk = 1;
            lwrk = 1;
        }
        work[1 - 1] = lwrk;
        if (wantsn || n == 0) {
            liwmin = 1;
        } else {
            liwmin = n + 2;
        }
        iwork[1 - 1] = liwmin;
        //
        if (lwork < minwrk && !lquery) {
            info = -21;
        } else if (liwork < liwmin && !lquery) {
            info = -24;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cggesx", -info);
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
    if (ilbscl) {
        Clascl("G", 0, 0, bnrm, bnrmto, n, n, b, ldb, ierr);
    }
    //
    //     Permute the matrix to make it more nearly triangular
    //     (Real Workspace: need 6*N)
    //
    ileft = 1;
    iright = n + 1;
    irwrk = iright + n;
    Cggbal("P", n, a, lda, b, ldb, ilo, ihi, &rwork[ileft - 1], &rwork[iright - 1], &rwork[irwrk - 1], ierr);
    //
    //     Reduce B to triangular form (QR decomposition of B)
    //     (Complex Workspace: need N, prefer N*NB)
    //
    irows = ihi + 1 - ilo;
    icols = n + 1 - ilo;
    itau = 1;
    iwrk = itau + irows;
    Cgeqrf(irows, icols, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, &work[itau - 1], &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    //     Apply the unitary transformation to matrix A
    //     (Complex Workspace: need N, prefer N*NB)
    //
    Cunmqr("L", "C", irows, icols, irows, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, &work[itau - 1], &a[(ilo - 1) + (ilo - 1) * lda], lda, &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    //     Initialize VSL
    //     (Complex Workspace: need N, prefer N*NB)
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
    //     (Workspace: none needed)
    //
    Cgghrd(jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, vsl, ldvsl, vsr, ldvsr, ierr);
    //
    sdim = 0;
    //
    //     Perform QZ algorithm, computing Schur vectors if desired
    //     (Complex Workspace: need N)
    //     (Real Workspace:    need N)
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
        goto statement_40;
    }
    //
    //     Sort eigenvalues ALPHA/BETA and compute the reciprocal of
    //     condition number(s)
    //
    if (wantst) {
        //
        //        Undo scaling on eigenvalues before SELCTGing
        //
        if (ilascl) {
            Clascl("G", 0, 0, anrmto, anrm, n, 1, alpha, n, ierr);
        }
        if (ilbscl) {
            Clascl("G", 0, 0, bnrmto, bnrm, n, 1, beta, n, ierr);
        }
        //
        //        Select eigenvalues
        //
        for (i = 1; i <= n; i = i + 1) {
            bwork[i - 1] = selctg(alpha[i - 1], beta[i - 1]);
        }
        //
        //        Reorder eigenvalues, transform Generalized Schur vectors, and
        //        compute reciprocal condition numbers
        //        (Complex Workspace: If IJOB >= 1, need MAX(1, 2*SDIM*(N-SDIM))
        //                            otherwise, need 1 )
        //
        Ctgsen(ijob, ilvsl, ilvsr, bwork, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, sdim, pl, pr, dif, &work[iwrk - 1], lwork - iwrk + 1, iwork, liwork, ierr);
        //
        if (ijob >= 1) {
            maxwrk = max(maxwrk, 2 * sdim * (n - sdim));
        }
        if (ierr == -21) {
            //
            //            not enough complex workspace
            //
            info = -21;
        } else {
            if (ijob == 1 || ijob == 4) {
                rconde[1 - 1] = pl;
                rconde[2 - 1] = pr;
            }
            if (ijob == 2 || ijob == 4) {
                rcondv[1 - 1] = dif[1 - 1];
                rcondv[2 - 1] = dif[2 - 1];
            }
            if (ierr == 1) {
                info = n + 3;
            }
        }
        //
    }
    //
    //     Apply permutation to VSL and VSR
    //     (Workspace: none needed)
    //
    if (ilvsl) {
        Cggbak("P", "L", n, ilo, ihi, &rwork[ileft - 1], &rwork[iright - 1], n, vsl, ldvsl, ierr);
    }
    //
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
statement_40:
    //
    work[1 - 1] = maxwrk;
    iwork[1 - 1] = liwmin;
    //
    //     End of Cggesx
    //
}
