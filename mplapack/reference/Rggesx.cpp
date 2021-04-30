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

void Rggesx(const char *jobvsl, const char *jobvsr, const char *sort, bool (*selctg)(REAL, REAL, REAL), const char *sense, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, INTEGER &sdim, REAL *alphar, REAL *alphai, REAL *beta, REAL *vsl, INTEGER const ldvsl, REAL *vsr, INTEGER const ldvsr, REAL *rconde, REAL *rcondv, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, bool *bwork, INTEGER &info) {
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
    REAL safmin = 0.0;
    const REAL one = 1.0;
    REAL safmax = 0.0;
    REAL smlnum = 0.0;
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
    INTEGER i = 0;
    REAL pl = 0.0;
    REAL pr = 0.0;
    REAL dif[2];
    bool lastsl = false;
    bool lst2sl = false;
    INTEGER ip = 0;
    bool cursl = false;
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
        info = -16;
    } else if (ldvsr < 1 || (ilvsr && ldvsr < n)) {
        info = -18;
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
            minwrk = max(8 * n, 6 * n + 16);
            maxwrk = minwrk - n + n * iMlaenv(1, "Rgeqrf", " ", n, 1, n, 0);
            maxwrk = max({maxwrk, minwrk - n + n * iMlaenv(1, "Rormqr", " ", n, 1, n, -1)});
            if (ilvsl) {
                maxwrk = max({maxwrk, minwrk - n + n * iMlaenv(1, "Rorgqr", " ", n, 1, n, -1)});
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
            liwmin = n + 6;
        }
        iwork[1 - 1] = liwmin;
        //
        if (lwork < minwrk && !lquery) {
            info = -22;
        } else if (liwork < liwmin && !lquery) {
            info = -24;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rggesx", -info);
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
    safmin = Rlamch("S");
    safmax = one / safmin;
    Rlabad(safmin, safmax);
    smlnum = sqrt(safmin) / eps;
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
    //     Permute the matrix to make it more nearly triangular
    //     (Workspace: need 6*N + 2*N for permutation parameters)
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
    icols = n + 1 - ilo;
    itau = iwrk;
    iwrk = itau + irows;
    Rgeqrf(irows, icols, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, &work[itau - 1], &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    //     Apply the orthogonal transformation to matrix A
    //     (Workspace: need N, prefer N*NB)
    //
    Rormqr("L", "T", irows, icols, irows, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, &work[itau - 1], &a[(ilo - 1) + (ilo - 1) * lda], lda, &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    //     Initialize VSL
    //     (Workspace: need N, prefer N*NB)
    //
    if (ilvsl) {
        Rlaset("Full", n, n, zero, one, vsl, ldvsl);
        if (irows > 1) {
            Rlacpy("L", irows - 1, irows - 1, &b[((ilo + 1) - 1) + (ilo - 1) * ldb], ldb, &vsl[((ilo + 1) - 1) + (ilo - 1) * ldvsl], ldvsl);
        }
        Rorgqr(irows, irows, irows, &vsl[(ilo - 1) + (ilo - 1) * ldvsl], ldvsl, &work[itau - 1], &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    }
    //
    //     Initialize VSR
    //
    if (ilvsr) {
        Rlaset("Full", n, n, zero, one, vsr, ldvsr);
    }
    //
    //     Reduce to generalized Hessenberg form
    //     (Workspace: none needed)
    //
    Rgghrd(jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, vsl, ldvsl, vsr, ldvsr, ierr);
    //
    sdim = 0;
    //
    //     Perform QZ algorithm, computing Schur vectors if desired
    //     (Workspace: need N)
    //
    iwrk = itau;
    Rhgeqz("S", jobvsl, jobvsr, n, ilo, ihi, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, &work[iwrk - 1], lwork + 1 - iwrk, ierr);
    if (ierr != 0) {
        if (ierr > 0 && ierr <= n) {
            info = ierr;
        } else if (ierr > n && ierr <= 2 * n) {
            info = ierr - n;
        } else {
            info = n + 1;
        }
        goto statement_60;
    }
    //
    //     Sort eigenvalues ALPHA/BETA and compute the reciprocal of
    //     condition number(s)
    //     (Workspace: If IJOB >= 1, need MAX( 8*(N+1), 2*SDIM*(N-SDIM) )
    //                 otherwise, need 8*(N+1) )
    //
    if (wantst) {
        //
        //        Undo scaling on eigenvalues before SELCTGing
        //
        if (ilascl) {
            Rlascl("G", 0, 0, anrmto, anrm, n, 1, alphar, n, ierr);
            Rlascl("G", 0, 0, anrmto, anrm, n, 1, alphai, n, ierr);
        }
        if (ilbscl) {
            Rlascl("G", 0, 0, bnrmto, bnrm, n, 1, beta, n, ierr);
        }
        //
        //        Select eigenvalues
        //
        for (i = 1; i <= n; i = i + 1) {
            bwork[i - 1] = selctg(alphar[i - 1], alphai[i - 1], beta[i - 1]);
        }
        //
        //        Reorder eigenvalues, transform Generalized Schur vectors, and
        //        compute reciprocal condition numbers
        //
        Rtgsen(ijob, ilvsl, ilvsr, bwork, n, a, lda, b, ldb, alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, sdim, pl, pr, dif, &work[iwrk - 1], lwork - iwrk + 1, iwork, liwork, ierr);
        //
        if (ijob >= 1) {
            maxwrk = max(maxwrk, 2 * sdim * (n - sdim));
        }
        if (ierr == -22) {
            //
            //            not enough real workspace
            //
            info = -22;
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
        Rggbak("P", "L", n, ilo, ihi, &work[ileft - 1], &work[iright - 1], n, vsl, ldvsl, ierr);
    }
    //
    if (ilvsr) {
        Rggbak("P", "R", n, ilo, ihi, &work[ileft - 1], &work[iright - 1], n, vsr, ldvsr, ierr);
    }
    //
    //     Check if unscaling would cause over/underflow, if so, rescale
    //     (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
    //     B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)
    //
    if (ilascl) {
        for (i = 1; i <= n; i = i + 1) {
            if (alphai[i - 1] != zero) {
                if ((alphar[i - 1] / safmax) > (anrmto / anrm) || (safmin / alphar[i - 1]) > (anrm / anrmto)) {
                    work[1 - 1] = abs(a[(i - 1) + (i - 1) * lda] / alphar[i - 1]);
                    beta[i - 1] = beta[i - 1] * work[1 - 1];
                    alphar[i - 1] = alphar[i - 1] * work[1 - 1];
                    alphai[i - 1] = alphai[i - 1] * work[1 - 1];
                } else if ((alphai[i - 1] / safmax) > (anrmto / anrm) || (safmin / alphai[i - 1]) > (anrm / anrmto)) {
                    work[1 - 1] = abs(a[(i - 1) + ((i + 1) - 1) * lda] / alphai[i - 1]);
                    beta[i - 1] = beta[i - 1] * work[1 - 1];
                    alphar[i - 1] = alphar[i - 1] * work[1 - 1];
                    alphai[i - 1] = alphai[i - 1] * work[1 - 1];
                }
            }
        }
    }
    //
    if (ilbscl) {
        for (i = 1; i <= n; i = i + 1) {
            if (alphai[i - 1] != zero) {
                if ((beta[i - 1] / safmax) > (bnrmto / bnrm) || (safmin / beta[i - 1]) > (bnrm / bnrmto)) {
                    work[1 - 1] = abs(b[(i - 1) + (i - 1) * ldb] / beta[i - 1]);
                    beta[i - 1] = beta[i - 1] * work[1 - 1];
                    alphar[i - 1] = alphar[i - 1] * work[1 - 1];
                    alphai[i - 1] = alphai[i - 1] * work[1 - 1];
                }
            }
        }
    }
    //
    //     Undo scaling
    //
    if (ilascl) {
        Rlascl("H", 0, 0, anrmto, anrm, n, n, a, lda, ierr);
        Rlascl("G", 0, 0, anrmto, anrm, n, 1, alphar, n, ierr);
        Rlascl("G", 0, 0, anrmto, anrm, n, 1, alphai, n, ierr);
    }
    //
    if (ilbscl) {
        Rlascl("U", 0, 0, bnrmto, bnrm, n, n, b, ldb, ierr);
        Rlascl("G", 0, 0, bnrmto, bnrm, n, 1, beta, n, ierr);
    }
    //
    if (wantst) {
        //
        //        Check if reordering is correct
        //
        lastsl = true;
        lst2sl = true;
        sdim = 0;
        ip = 0;
        for (i = 1; i <= n; i = i + 1) {
            cursl = selctg(alphar[i - 1], alphai[i - 1], beta[i - 1]);
            if (alphai[i - 1] == zero) {
                if (cursl) {
                    sdim++;
                }
                ip = 0;
                if (cursl && !lastsl) {
                    info = n + 2;
                }
            } else {
                if (ip == 1) {
                    //
                    //                 Last eigenvalue of conjugate pair
                    //
                    cursl = cursl || lastsl;
                    lastsl = cursl;
                    if (cursl) {
                        sdim += 2;
                    }
                    ip = -1;
                    if (cursl && !lst2sl) {
                        info = n + 2;
                    }
                } else {
                    //
                    //                 First eigenvalue of conjugate pair
                    //
                    ip = 1;
                }
            }
            lst2sl = lastsl;
            lastsl = cursl;
        }
        //
    }
//
statement_60:
    //
    work[1 - 1] = maxwrk;
    iwork[1 - 1] = liwmin;
    //
    //     End of Rggesx
    //
}
