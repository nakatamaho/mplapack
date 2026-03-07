/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ZLAQZ0.
// Original LAPACK authors:
//   Thijs Steel, KU Leuven

#include <mpblas.h>
#include <mplapack.h>

void Claqz0(const char *wants, const char *wantq, const char *wantz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *alpha, COMPLEX *beta, COMPLEX *q, INTEGER const ldq, COMPLEX *z, INTEGER const ldz, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER const rec, INTEGER &info) {
    bool ilschur = false;
    INTEGER iwants = 0;
    bool ilq = false;
    INTEGER iwantq = 0;
    bool ilz = false;
    INTEGER iwantz = 0;
    char jbcmpz[3];
    INTEGER nmin = 0;
    INTEGER nwr = 0;
    INTEGER nibble = 0;
    INTEGER nsr = 0;
    INTEGER rcost = 0;
    INTEGER itemp1 = 0;
    INTEGER nbr = 0;
    INTEGER nw = 0;
    INTEGER n_undeflated = 0;
    INTEGER n_deflated = 0;
    INTEGER aed_info = 0;
    INTEGER sweep_info = 0;
    INTEGER itemp2 = 0;
    INTEGER lworkreq = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    REAL safmin = 0.0;
    const REAL one = 1.0;
    REAL safmax = 0.0;
    REAL ulp = 0.0;
    REAL smlnum = 0.0;
    REAL bnorm = 0.0;
    REAL btol = 0.0;
    INTEGER istart = 0;
    INTEGER istop = 0;
    INTEGER maxit = 0;
    INTEGER ld = 0;
    INTEGER iiter = 0;
    COMPLEX eshift = 0.0;
    INTEGER istart2 = 0;
    INTEGER k = 0;
    INTEGER istartm = 0;
    INTEGER istopm = 0;
    INTEGER k2 = 0;
    REAL c1 = 0.0;
    COMPLEX s1 = 0.0;
    COMPLEX temp = 0.0;
    INTEGER nshifts = 0;
    INTEGER nblock = 0;
    INTEGER ns = 0;
    INTEGER shiftpos = 0;
    INTEGER norm_info = 0;
    //
    // Arguments
    //
    // Parameters
    //
    // Local scalars
    //
    // External Functions
    //
    // Decode wantS,wantQ,wantZ
    //
    if (Mlsame(wants, "E")) {
        ilschur = false;
        iwants = 1;
    } else if (Mlsame(wants, "S")) {
        ilschur = true;
        iwants = 2;
    } else {
        iwants = 0;
    }
    //
    if (Mlsame(wantq, "N")) {
        ilq = false;
        iwantq = 1;
    } else if (Mlsame(wantq, "V")) {
        ilq = true;
        iwantq = 2;
    } else if (Mlsame(wantq, "I")) {
        ilq = true;
        iwantq = 3;
    } else {
        iwantq = 0;
    }
    //
    if (Mlsame(wantz, "N")) {
        ilz = false;
        iwantz = 1;
    } else if (Mlsame(wantz, "V")) {
        ilz = true;
        iwantz = 2;
    } else if (Mlsame(wantz, "I")) {
        ilz = true;
        iwantz = 3;
    } else {
        iwantz = 0;
    }
    //
    // Check Argument Values
    //
    info = 0;
    if (iwants == 0) {
        info = -1;
    } else if (iwantq == 0) {
        info = -2;
    } else if (iwantz == 0) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (ilo < 1) {
        info = -5;
    } else if (ihi > n || ihi < ilo - 1) {
        info = -6;
    } else if (lda < n) {
        info = -8;
    } else if (ldb < n) {
        info = -10;
    } else if (ldq < 1 || (ilq && ldq < n)) {
        info = -15;
    } else if (ldz < 1 || (ilz && ldz < n)) {
        info = -17;
    }
    if (info != 0) {
        Mxerbla("Claqz0", -info);
        return;
    }
    //
    // Quick return if possible
    //
    if (n <= 0) {
        work[1 - 1] = castREAL(1);
        return;
    }
    //
    // Get the parameters
    //
    jbcmpz[0] = *wants;
    jbcmpz[1] = *wantq;
    jbcmpz[2] = *wantz;
    //
    nmin = iMlaenv(12, "Claqz0", jbcmpz, n, ilo, ihi, lwork);
    //
    nwr = iMlaenv(13, "Claqz0", jbcmpz, n, ilo, ihi, lwork);
    nwr = max((INTEGER)2, nwr);
    nwr = min(ihi - ilo + 1, (n - 1) / 3, nwr);
    //
    nibble = iMlaenv(14, "Claqz0", jbcmpz, n, ilo, ihi, lwork);
    //
    nsr = iMlaenv(15, "Claqz0", jbcmpz, n, ilo, ihi, lwork);
    nsr = min(nsr, (n + 6) / 9, ihi - ilo);
    nsr = max((INTEGER)2, nsr - mod(nsr, 2));
    //
    rcost = iMlaenv(17, "Claqz0", jbcmpz, n, ilo, ihi, lwork);
    itemp1 = castINTEGER(nsr / sqrt(1 + 2 * nsr / (castREAL(rcost) / 100 * n)));
    itemp1 = ((itemp1 - 1) / 4) * 4 + 4;
    nbr = nsr + itemp1;
    //
    if (n < nmin || rec >= 2) {
        Chgeqz(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, info);
        return;
    }
    //
    // Find out required workspace
    //
    // Workspace query to Claqz2
    nw = max(nwr, nmin);
    Claqz2(ilschur, ilq, ilz, n, ilo, ihi, nw, a, lda, b, ldb, q, ldq, z, ldz, n_undeflated, n_deflated, alpha, beta, work, nw, work, nw, work, -1, rwork, rec, aed_info);
    itemp1 = castINTEGER(work[1 - 1].real());
    // Workspace query to Claqz3
    Claqz3(ilschur, ilq, ilz, n, ilo, ihi, nsr, nbr, alpha, beta, a, lda, b, ldb, q, ldq, z, ldz, work, nbr, work, nbr, work, -1, sweep_info);
    itemp2 = castINTEGER(work[1 - 1].real());
    //
    lworkreq = max(itemp1 + 2 * pow2(nw), itemp2 + 2 * pow2(nbr));
    if (lwork == -1) {
        work[1 - 1] = castREAL(lworkreq);
        return;
    } else if (lwork < lworkreq) {
        info = -19;
    }
    if (info != 0) {
        Mxerbla("Claqz0", info);
        return;
    }
    //
    // Initialize Q and Z
    //
    if (iwantq == 3) {
        Claset("FULL", n, n, czero, cone, q, ldq);
    }
    if (iwantz == 3) {
        Claset("FULL", n, n, czero, cone, z, ldz);
    }
    //
    // Get machine constants
    safmin = Rlamch("SAFE MINIMUM");
    safmax = one / safmin;
    ulp = Rlamch("PRECISION");
    smlnum = safmin * (castREAL(n) / ulp);
    //
    bnorm = Clanhs("F", ihi - ilo + 1, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, rwork);
    btol = max(safmin, ulp * bnorm);
    //
    istart = ilo;
    istop = ihi;
    maxit = 30 * (ihi - ilo + 1);
    ld = 0;
    //
    for (iiter = 1; iiter <= maxit; iiter = iiter + 1) {
        if (iiter >= maxit) {
            info = istop + 1;
            goto statement_80;
        }
        if (istart + 1 >= istop) {
            istop = istart;
            break;
        }
        //
        // Check deflations at the end
        if (abs(a[(istop - 1) + ((istop - 1) - 1) * lda]) <= max(smlnum, ulp * (abs(a[(istop - 1) + (istop - 1) * lda]) + abs(a[((istop - 1) - 1) + ((istop - 1) - 1) * lda])))) {
            a[(istop - 1) + ((istop - 1) - 1) * lda] = czero;
            istop = istop - 1;
            ld = 0;
            eshift = czero;
        }
        // Check deflations at the start
        if (abs(a[((istart + 1) - 1) + (istart - 1) * lda]) <= max(smlnum, ulp * (abs(a[(istart - 1) + (istart - 1) * lda]) + abs(a[((istart + 1) - 1) + ((istart + 1) - 1) * lda])))) {
            a[((istart + 1) - 1) + (istart - 1) * lda] = czero;
            istart++;
            ld = 0;
            eshift = czero;
        }
        //
        if (istart + 1 >= istop) {
            break;
        }
        //
        // Check interior deflations
        istart2 = istart;
        for (k = istop; k >= istart + 1; k = k - 1) {
            if (abs(a[(k - 1) + ((k - 1) - 1) * lda]) <= max(smlnum, ulp * (abs(a[(k - 1) + (k - 1) * lda]) + abs(a[((k - 1) - 1) + ((k - 1) - 1) * lda])))) {
                a[(k - 1) + ((k - 1) - 1) * lda] = czero;
                istart2 = k;
                break;
            }
        }
        //
        // Get range to apply rotations to
        if (ilschur) {
            istartm = 1;
            istopm = n;
        } else {
            istartm = istart2;
            istopm = istop;
        }
        //
        // Check infinite eigenvalues, this is done without blocking so might
        // slow down the method when many infinite eigenvalues are present
        k = istop;
        while (k >= istart2) {
            //
            if (abs(b[(k - 1) + (k - 1) * ldb]) < btol) {
                // A diagonal element of B is negligible, move it
                // to the top and deflate it
                //
                for (k2 = k; k2 >= istart2 + 1; k2 = k2 - 1) {
                    Clartg(b[((k2 - 1) - 1) + (k2 - 1) * ldb], b[((k2 - 1) - 1) + ((k2 - 1) - 1) * ldb], c1, s1, temp);
                    b[((k2 - 1) - 1) + (k2 - 1) * ldb] = temp;
                    b[((k2 - 1) - 1) + ((k2 - 1) - 1) * ldb] = czero;
                    //
                    Crot(k2 - 2 - istartm + 1, &b[(istartm - 1) + (k2 - 1) * ldb], 1, &b[(istartm - 1) + ((k2 - 1) - 1) * ldb], 1, c1, s1);
                    Crot(min(k2 + 1, istop) - istartm + 1, &a[(istartm - 1) + (k2 - 1) * lda], 1, &a[(istartm - 1) + ((k2 - 1) - 1) * lda], 1, c1, s1);
                    if (ilz) {
                        Crot(n, &z[(k2 - 1) * ldz], 1, &z[((k2 - 1) - 1) * ldz], 1, c1, s1);
                    }
                    //
                    if (k2 < istop) {
                        Clartg(a[(k2 - 1) + ((k2 - 1) - 1) * lda], a[((k2 + 1) - 1) + ((k2 - 1) - 1) * lda], c1, s1, temp);
                        a[(k2 - 1) + ((k2 - 1) - 1) * lda] = temp;
                        a[((k2 + 1) - 1) + ((k2 - 1) - 1) * lda] = czero;
                        //
                        Crot(istopm - k2 + 1, &a[(k2 - 1) + (k2 - 1) * lda], lda, &a[((k2 + 1) - 1) + (k2 - 1) * lda], lda, c1, s1);
                        Crot(istopm - k2 + 1, &b[(k2 - 1) + (k2 - 1) * ldb], ldb, &b[((k2 + 1) - 1) + (k2 - 1) * ldb], ldb, c1, s1);
                        if (ilq) {
                            Crot(n, &q[(k2 - 1) * ldq], 1, &q[((k2 + 1) - 1) * ldq], 1, c1, conj(s1));
                        }
                    }
                    //
                }
                //
                if (istart2 < istop) {
                    Clartg(a[(istart2 - 1) + (istart2 - 1) * lda], a[((istart2 + 1) - 1) + (istart2 - 1) * lda], c1, s1, temp);
                    a[(istart2 - 1) + (istart2 - 1) * lda] = temp;
                    a[((istart2 + 1) - 1) + (istart2 - 1) * lda] = czero;
                    //
                    Crot(istopm - (istart2 + 1) + 1, &a[(istart2 - 1) + ((istart2 + 1) - 1) * lda], lda, &a[((istart2 + 1) - 1) + ((istart2 + 1) - 1) * lda], lda, c1, s1);
                    Crot(istopm - (istart2 + 1) + 1, &b[(istart2 - 1) + ((istart2 + 1) - 1) * ldb], ldb, &b[((istart2 + 1) - 1) + ((istart2 + 1) - 1) * ldb], ldb, c1, s1);
                    if (ilq) {
                        Crot(n, &q[(istart2 - 1) * ldq], 1, &q[((istart2 + 1) - 1) * ldq], 1, c1, conj(s1));
                    }
                }
                //
                istart2++;
                //
            }
            k = k - 1;
        }
        //
        // istart2 now points to the top of the bottom right
        // unreduced Hessenberg block
        if (istart2 >= istop) {
            istop = istart2 - 1;
            ld = 0;
            eshift = czero;
            continue;
        }
        //
        nw = nwr;
        nshifts = nsr;
        nblock = nbr;
        //
        if (istop - istart2 + 1 < nmin) {
            // Setting nw to the size of the subblock will make AED deflate
            // all the eigenvalues. This is slightly more efficient than just
            // using qz_small because the off diagonal part gets updated via BLAS.
            if (istop - istart + 1 < nmin) {
                nw = istop - istart + 1;
                istart2 = istart;
            } else {
                nw = istop - istart2 + 1;
            }
        }
        //
        // Time for AED
        //
        Claqz2(ilschur, ilq, ilz, n, istart2, istop, nw, a, lda, b, ldb, q, ldq, z, ldz, n_undeflated, n_deflated, alpha, beta, work, nw, &work[(pow2(nw) + 1) - 1], nw, &work[(2 * pow2(nw) + 1) - 1], lwork - 2 * pow2(nw), rwork, rec, aed_info);
        //
        if (n_deflated > 0) {
            istop = istop - n_deflated;
            ld = 0;
            eshift = czero;
        }
        //
        if (100 * n_deflated > nibble * (n_deflated + n_undeflated) || istop - istart2 + 1 < nmin) {
            // AED has uncovered many eigenvalues. Skip a QZ sweep and run
            // AED again.
            continue;
        }
        //
        ld++;
        //
        ns = min(nshifts, istop - istart2);
        ns = min(ns, n_undeflated);
        shiftpos = istop - n_undeflated + 1;
        //
        if (mod(ld, 6) == 0) {
            //
            // Exceptional shift.  Chosen for no particularly good reason.
            //
            if ((castREAL(maxit) * safmin) * abs(a[(istop - 1) + ((istop - 1) - 1) * lda]) < abs(a[((istop - 1) - 1) + ((istop - 1) - 1) * lda])) {
                eshift = a[(istop - 1) + ((istop - 1) - 1) * lda] / b[((istop - 1) - 1) + ((istop - 1) - 1) * ldb];
            } else {
                eshift += cone / (safmin * castREAL(maxit));
            }
            alpha[shiftpos - 1] = cone;
            beta[shiftpos - 1] = eshift;
            ns = 1;
        }
        //
        // Time for a QZ sweep
        //
        Claqz3(ilschur, ilq, ilz, n, istart2, istop, ns, nblock, &alpha[shiftpos - 1], &beta[shiftpos - 1], a, lda, b, ldb, q, ldq, z, ldz, work, nblock, &work[(pow2(nblock) + 1) - 1], nblock, &work[(2 * pow2(nblock) + 1) - 1], lwork - 2 * pow2(nblock), sweep_info);
        //
    }
//
// Call Chgeqz to normalize the eigenvalue blocks and set the eigenvalues
// If all the eigenvalues have been found, Chgeqz will not do any iterations
// and only normalize the blocks. In case of a rare convergence failure,
// the single shift might perform better.
//
statement_80:
    Chgeqz(wants, wantq, wantz, n, ilo, ihi, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz, work, lwork, rwork, norm_info);
    //
    info = norm_info;
    //
}
