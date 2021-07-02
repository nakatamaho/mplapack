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

void Rgees(const char *jobvs, const char *sort, bool (*select)(REAL, REAL), INTEGER const n, REAL *a, INTEGER const lda, INTEGER &sdim, REAL *wr, REAL *wi, REAL *vs, INTEGER const ldvs, REAL *work, INTEGER const lwork, bool *bwork, INTEGER &info) {
    bool lquery = false;
    bool wantvs = false;
    bool wantst = false;
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER ieval = 0;
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
    INTEGER ierr = 0;
    INTEGER ibal = 0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    INTEGER itau = 0;
    INTEGER iwrk = 0;
    INTEGER i = 0;
    REAL s = 0.0;
    REAL sep = 0.0;
    INTEGER idum[1];
    INTEGER icond = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER inxt = 0;
    bool lastsl = false;
    bool lst2sl = false;
    INTEGER ip = 0;
    bool cursl = false;
    //
    //     Test the input arguments
    //
    info = 0;
    lquery = (lwork == -1);
    wantvs = Mlsame(jobvs, "V");
    wantst = Mlsame(sort, "S");
    if ((!wantvs) && (!Mlsame(jobvs, "N"))) {
        info = -1;
    } else if ((!wantst) && (!Mlsame(sort, "N"))) {
        info = -2;
    } else if (n < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldvs < 1 || (wantvs && ldvs < n)) {
        info = -11;
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
            maxwrk = 2 * n + n * iMlaenv(1, "Rgehrd", " ", n, 1, n, 0);
            minwrk = 3 * n;
            //
            Rhseqr("S", jobvs, n, 1, n, a, lda, wr, wi, vs, ldvs, work, -1, ieval);
            hswork = castINTEGER(work[1 - 1]);
            //
            if (!wantvs) {
                maxwrk = max(maxwrk, n + hswork);
            } else {
                maxwrk = max(maxwrk, 2 * n + (n - 1) * iMlaenv(1, "Rorghr", " ", n, 1, n, -1));
                maxwrk = max(maxwrk, n + hswork);
            }
        }
        work[1 - 1] = maxwrk;
        //
        if (lwork < minwrk && !lquery) {
            info = -13;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rgees", -info);
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
    anrm = Rlange("M", n, n, a, lda, dum);
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
    //     Permute the matrix to make it more nearly triangular
    //     (Workspace: need N)
    //
    ibal = 1;
    Rgebal("P", n, a, lda, ilo, ihi, &work[ibal - 1], ierr);
    //
    //     Reduce to upper Hessenberg form
    //     (Workspace: need 3*N, prefer 2*N+N*NB)
    //
    itau = n + ibal;
    iwrk = n + itau;
    Rgehrd(n, ilo, ihi, a, lda, &work[itau - 1], &work[iwrk - 1], lwork - iwrk + 1, ierr);
    //
    if (wantvs) {
        //
        //        Copy Householder vectors to VS
        //
        Rlacpy("L", n, n, a, lda, vs, ldvs);
        //
        //        Generate orthogonal matrix in VS
        //        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
        //
        Rorghr(n, ilo, ihi, vs, ldvs, &work[itau - 1], &work[iwrk - 1], lwork - iwrk + 1, ierr);
    }
    //
    sdim = 0;
    //
    //     Perform QR iteration, accumulating Schur vectors in VS if desired
    //     (Workspace: need N+1, prefer N+HSWORK (see comments) )
    //
    iwrk = itau;
    Rhseqr("S", jobvs, n, ilo, ihi, a, lda, wr, wi, vs, ldvs, &work[iwrk - 1], lwork - iwrk + 1, ieval);
    if (ieval > 0) {
        info = ieval;
    }
    //
    //     Sort eigenvalues if desired
    //
    if (wantst && info == 0) {
        if (scalea) {
            Rlascl("G", 0, 0, cscale, anrm, n, 1, wr, n, ierr);
            Rlascl("G", 0, 0, cscale, anrm, n, 1, wi, n, ierr);
        }
        for (i = 1; i <= n; i = i + 1) {
            bwork[i - 1] = select(wr[i - 1], wi[i - 1]);
        }
        //
        //        Reorder eigenvalues and transform Schur vectors
        //        (Workspace: none needed)
        //
        Rtrsen("N", jobvs, bwork, n, a, lda, vs, ldvs, wr, wi, sdim, s, sep, &work[iwrk - 1], lwork - iwrk + 1, idum, 1, icond);
        if (icond > 0) {
            info = n + icond;
        }
    }
    //
    if (wantvs) {
        //
        //        Undo balancing
        //        (Workspace: need N)
        //
        Rgebak("P", "R", n, ilo, ihi, &work[ibal - 1], n, vs, ldvs, ierr);
    }
    //
    if (scalea) {
        //
        //        Undo scaling for the Schur form of A
        //
        Rlascl("H", 0, 0, cscale, anrm, n, n, a, lda, ierr);
        Rcopy(n, a, lda + 1, wr, 1);
        if (cscale == smlnum) {
            //
            //           If scaling back towards underflow, adjust WI if an
            //           offdiagonal element of a 2-by-2 block in the Schur form
            //           underflows.
            //
            if (ieval > 0) {
                i1 = ieval + 1;
                i2 = ihi - 1;
                Rlascl("G", 0, 0, cscale, anrm, ilo - 1, 1, wi, max(ilo - 1, (INTEGER)1), ierr);
            } else if (wantst) {
                i1 = 1;
                i2 = n - 1;
            } else {
                i1 = ilo;
                i2 = ihi - 1;
            }
            inxt = i1 - 1;
            for (i = i1; i <= i2; i = i + 1) {
                if (i < inxt) {
                    goto statement_20;
                }
                if (wi[i - 1] == zero) {
                    inxt = i + 1;
                } else {
                    if (a[((i + 1) - 1) + (i - 1) * lda] == zero) {
                        wi[i - 1] = zero;
                        wi[(i + 1) - 1] = zero;
                    } else if (a[((i + 1) - 1) + (i - 1) * lda] != zero && a[(i - 1) + ((i + 1) - 1) * lda] == zero) {
                        wi[i - 1] = zero;
                        wi[(i + 1) - 1] = zero;
                        if (i > 1) {
                            Rswap(i - 1, &a[(i - 1) * lda], 1, &a[((i + 1) - 1) * lda], 1);
                        }
                        if (n > i + 1) {
                            Rswap(n - i - 1, &a[(i - 1) + ((i + 2) - 1) * lda], lda, &a[((i + 1) - 1) + ((i + 2) - 1) * lda], lda);
                        }
                        if (wantvs) {
                            Rswap(n, &vs[(i - 1) * ldvs], 1, &vs[((i + 1) - 1) * ldvs], 1);
                        }
                        a[(i - 1) + ((i + 1) - 1) * lda] = a[((i + 1) - 1) + (i - 1) * lda];
                        a[((i + 1) - 1) + (i - 1) * lda] = zero;
                    }
                    inxt = i + 2;
                }
            statement_20:;
            }
        }
        //
        //        Undo scaling for the imaginary part of the eigenvalues
        //
        Rlascl("G", 0, 0, cscale, anrm, n - ieval, 1, &wi[(ieval + 1) - 1], max(n - ieval, (INTEGER)1), ierr);
    }
    //
    if (wantst && info == 0) {
        //
        //        Check if reordering successful
        //
        lastsl = true;
        lst2sl = true;
        sdim = 0;
        ip = 0;
        for (i = 1; i <= n; i = i + 1) {
            cursl = select(wr[i - 1], wi[i - 1]);
            if (wi[i - 1] == zero) {
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
    }
    //
    work[1 - 1] = maxwrk;
    //
    //     End of Rgees
    //
}
