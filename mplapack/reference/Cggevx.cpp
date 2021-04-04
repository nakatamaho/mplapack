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

void Cggevx(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, INTEGER const &n, COMPLEX *a, INTEGER const &lda, COMPLEX *b, INTEGER const &ldb, COMPLEX *alpha, COMPLEX *beta, COMPLEX *vl, INTEGER const &ldvl, COMPLEX *vr, INTEGER const &ldvr, INTEGER const &ilo, INTEGER const &ihi, REAL *lscale, REAL *rscale, REAL &abnrm, REAL &bbnrm, REAL *rconde, REAL *rcondv, COMPLEX *work, INTEGER const &lwork, REAL *rwork, INTEGER *iwork, arr_ref<bool> bwork, INTEGER &info) {
    COMPLEX x = 0.0;
    INTEGER ijobvl = 0;
    bool ilvl = false;
    INTEGER ijobvr = 0;
    bool ilvr = false;
    bool ilv = false;
    bool noscl = false;
    bool wantsn = false;
    bool wantse = false;
    bool wantsv = false;
    bool wantsb = false;
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
    INTEGER irows = 0;
    INTEGER icols = 0;
    INTEGER itau = 0;
    INTEGER iwrk = 0;
    const COMPLEX czero = (0.0, 0.0);
    const COMPLEX cone = (1.0, 0.0);
    str<1> chtemp = char0;
    arr_1d<1, bool> ldumma(fill0);
    INTEGER in = 0;
    INTEGER i = 0;
    INTEGER j = 0;
    INTEGER iwrk1 = 0;
    INTEGER m = 0;
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
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    abs1[x - 1] = abs(x.real()) + abs(x.imag());
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
    noscl = Mlsame(balanc, "N") || Mlsame(balanc, "P");
    wantsn = Mlsame(sense, "N");
    wantse = Mlsame(sense, "E");
    wantsv = Mlsame(sense, "V");
    wantsb = Mlsame(sense, "B");
    //
    //     Test the input arguments
    //
    info = 0;
    lquery = (lwork == -1);
    if (!(noscl || Mlsame(balanc, "S") || Mlsame(balanc, "B"))) {
        info = -1;
    } else if (ijobvl <= 0) {
        info = -2;
    } else if (ijobvr <= 0) {
        info = -3;
    } else if (!(wantsn || wantse || wantsb || wantsv)) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (lda < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -9;
    } else if (ldvl < 1 || (ilvl && ldvl < n)) {
        info = -13;
    } else if (ldvr < 1 || (ilvr && ldvr < n)) {
        info = -15;
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
        if (n == 0) {
            minwrk = 1;
            maxwrk = 1;
        } else {
            minwrk = 2 * n;
            if (wantse) {
                minwrk = 4 * n;
            } else if (wantsv || wantsb) {
                minwrk = 2 * n * (n + 1);
            }
            maxwrk = minwrk;
            maxwrk = max(maxwrk, n + n * iMlaenv[("Cgeqrf" - 1) * ldiMlaenv]);
            maxwrk = max(maxwrk, n + n * iMlaenv[("Cunmqr" - 1) * ldiMlaenv]);
            if (ilvl) {
                maxwrk = max(maxwrk, n + n * iMlaenv[("Cungqr" - 1) * ldiMlaenv]);
            }
        }
        work[1 - 1] = maxwrk;
        //
        if (lwork < minwrk && !lquery) {
            info = -25;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cggevx", -info);
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
    eps = dlamch("P");
    smlnum = dlamch("S");
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = one / smlnum;
    //
    //     Scale A if max element outside range [SMLNUM,BIGNUM]
    //
    anrm = Clange[("M" - 1) + (n - 1) * ldClange];
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
    bnrm = Clange[("M" - 1) + (n - 1) * ldClange];
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
    //     Permute and/or balance the matrix pair (A,B)
    //     (Real Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)
    //
    Cggbal(balanc, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, rwork, ierr);
    //
    //     Compute ABNRM and BBNRM
    //
    abnrm = Clange[("1" - 1) + (n - 1) * ldClange];
    if (ilascl) {
        rwork[1 - 1] = abnrm;
        Rlascl("G", 0, 0, anrmto, anrm, 1, 1, rwork[1 - 1], 1, ierr);
        abnrm = rwork[1 - 1];
    }
    //
    bbnrm = Clange[("1" - 1) + (n - 1) * ldClange];
    if (ilbscl) {
        rwork[1 - 1] = bbnrm;
        Rlascl("G", 0, 0, bnrmto, bnrm, 1, 1, rwork[1 - 1], 1, ierr);
        bbnrm = rwork[1 - 1];
    }
    //
    //     Reduce B to triangular form (QR decomposition of B)
    //     (Complex Workspace: need N, prefer N*NB )
    //
    irows = ihi + 1 - ilo;
    if (ilv || !wantsn) {
        icols = n + 1 - ilo;
    } else {
        icols = irows;
    }
    itau = 1;
    iwrk = itau + irows;
    Cgeqrf(irows, icols, b[(ilo - 1) + (ilo - 1) * ldb], ldb, work[itau - 1], work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    //     Apply the unitary transformation to A
    //     (Complex Workspace: need N, prefer N*NB)
    //
    Cunmqr("L", "C", irows, icols, irows, b[(ilo - 1) + (ilo - 1) * ldb], ldb, work[itau - 1], a[(ilo - 1) + (ilo - 1) * lda], lda, work[iwrk - 1], lwork + 1 - iwrk, ierr);
    //
    //     Initialize VL and/or VR
    //     (Workspace: need N, prefer N*NB)
    //
    if (ilvl) {
        Claset("Full", n, n, czero, cone, vl, ldvl);
        if (irows > 1) {
            Clacpy("L", irows - 1, irows - 1, b[((ilo + 1) - 1) + (ilo - 1) * ldb], ldb, vl[((ilo + 1) - 1) + (ilo - 1) * ldvl], ldvl);
        }
        Cungqr(irows, irows, irows, vl[(ilo - 1) + (ilo - 1) * ldvl], ldvl, work[itau - 1], work[iwrk - 1], lwork + 1 - iwrk, ierr);
    }
    //
    if (ilvr) {
        Claset("Full", n, n, czero, cone, vr, ldvr);
    }
    //
    //     Reduce to generalized Hessenberg form
    //     (Workspace: none needed)
    //
    if (ilv || !wantsn) {
        //
        //        Eigenvectors requested -- work on whole matrix.
        //
        Cgghrd(jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, vl, ldvl, vr, ldvr, ierr);
    } else {
        Cgghrd("N", "N", irows, 1, irows, a[(ilo - 1) + (ilo - 1) * lda], lda, b[(ilo - 1) + (ilo - 1) * ldb], ldb, vl, ldvl, vr, ldvr, ierr);
    }
    //
    //     Perform QZ algorithm (Compute eigenvalues, and optionally, the
    //     Schur forms and Schur vectors)
    //     (Complex Workspace: need N)
    //     (Real Workspace: need N)
    //
    iwrk = itau;
    if (ilv || !wantsn) {
        chtemp = "S";
    } else {
        chtemp = "E";
    }
    //
    Chgeqz(chtemp, jobvl, jobvr, n, ilo, ihi, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work[iwrk - 1], lwork + 1 - iwrk, rwork, ierr);
    if (ierr != 0) {
        if (ierr > 0 && ierr <= n) {
            info = ierr;
        } else if (ierr > n && ierr <= 2 * n) {
            info = ierr - n;
        } else {
            info = n + 1;
        }
        goto statement_90;
    }
    //
    //     Compute Eigenvectors and estimate condition numbers if desired
    //     Ctgevc: (Complex Workspace: need 2*N )
    //             (Real Workspace:    need 2*N )
    //     Ctgsna: (Complex Workspace: need 2*N*N if SENSE='V' or 'B')
    //             (Integer Workspace: need N+2 )
    //
    if (ilv || !wantsn) {
        if (ilv) {
            if (ilvl) {
                if (ilvr) {
                    chtemp = "B";
                } else {
                    chtemp = "L";
                }
            } else {
                chtemp = "R";
            }
            //
            Ctgevc(chtemp, "B", ldumma, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, n, in, work[iwrk - 1], rwork, ierr);
            if (ierr != 0) {
                info = n + 2;
                goto statement_90;
            }
        }
        //
        if (!wantsn) {
            //
            //           compute eigenvectors (Rtgevc) and estimate condition
            //           numbers (Rtgsna). Note that the definition of the condition
            //           number is not invariant under transformation (u,v) to
            //           (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
            //           Schur form (S,T), Q and Z are orthogonal matrices. In order
            //           to avoid using extra 2*N*N workspace, we have to
            //           re-calculate eigenvectors and estimate the condition numbers
            //           one at a time.
            //
            for (i = 1; i <= n; i = i + 1) {
                //
                for (j = 1; j <= n; j = j + 1) {
                    bwork[j - 1] = false;
                }
                bwork[i - 1] = true;
                //
                iwrk = n + 1;
                iwrk1 = iwrk + n;
                //
                if (wantse || wantsb) {
                    Ctgevc("B", "S", bwork, n, a, lda, b, ldb, work[1 - 1], n, work[iwrk - 1], n, 1, m, work[iwrk1 - 1], rwork, ierr);
                    if (ierr != 0) {
                        info = n + 2;
                        goto statement_90;
                    }
                }
                //
                Ctgsna(sense, "S", bwork, n, a, lda, b, ldb, work[1 - 1], n, work[iwrk - 1], n, rconde[i - 1], rcondv[i - 1], 1, m, work[iwrk1 - 1], lwork - iwrk1 + 1, iwork, ierr);
                //
            }
        }
    }
    //
    //     Undo balancing on VL and VR and normalization
    //     (Workspace: none needed)
    //
    if (ilvl) {
        Cggbak(balanc, "L", n, ilo, ihi, lscale, rscale, n, vl, ldvl, ierr);
        //
        for (jc = 1; jc <= n; jc = jc + 1) {
            temp = zero;
            for (jr = 1; jr <= n; jr = jr + 1) {
                temp = max(temp, abs1[vl[(jr - 1) + (jc - 1) * ldvl] - 1]);
            }
            if (temp < smlnum) {
                goto statement_50;
            }
            temp = one / temp;
            for (jr = 1; jr <= n; jr = jr + 1) {
                vl[(jr - 1) + (jc - 1) * ldvl] = vl[(jr - 1) + (jc - 1) * ldvl] * temp;
            }
        statement_50:;
        }
    }
    //
    if (ilvr) {
        Cggbak(balanc, "R", n, ilo, ihi, lscale, rscale, n, vr, ldvr, ierr);
        for (jc = 1; jc <= n; jc = jc + 1) {
            temp = zero;
            for (jr = 1; jr <= n; jr = jr + 1) {
                temp = max(temp, abs1[vr[(jr - 1) + (jc - 1) * ldvr] - 1]);
            }
            if (temp < smlnum) {
                goto statement_80;
            }
            temp = one / temp;
            for (jr = 1; jr <= n; jr = jr + 1) {
                vr[(jr - 1) + (jc - 1) * ldvr] = vr[(jr - 1) + (jc - 1) * ldvr] * temp;
            }
        statement_80:;
        }
    }
//
//     Undo scaling if necessary
//
statement_90:
    //
    if (ilascl) {
        Clascl("G", 0, 0, anrmto, anrm, n, 1, alpha, n, ierr);
    }
    //
    if (ilbscl) {
        Clascl("G", 0, 0, bnrmto, bnrm, n, 1, beta, n, ierr);
    }
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Cggevx
    //
}
