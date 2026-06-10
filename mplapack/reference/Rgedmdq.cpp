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

// Derived from LAPACK routine DGEDMDQ.
// Original LAPACK authors:
//   Univ. of Tennessee,
//   Univ. of California Berkeley,
//   Univ. of Colorado Denver,
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Rgedmdq(const char *jobs, const char *jobz, const char *jobr, const char *jobq, const char *jobt, const char *jobf, INTEGER const whtsvd, INTEGER const m, INTEGER const n, REAL *f, INTEGER const ldf, REAL *x, INTEGER const ldx, REAL *y, INTEGER const ldy, INTEGER const nrnk, REAL const tol, INTEGER &k, REAL *reig, REAL *imeig, REAL *z, INTEGER const ldz, REAL *res, REAL *b, INTEGER const ldb, REAL *v, INTEGER const ldv, REAL *s, INTEGER const lds, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
    //
    //
    // Test the input arguments
    bool wntres = Mlsame(jobr, "R");
    bool sccolx = Mlsame(jobs, "S") || Mlsame(jobs, "C");
    bool sccoly = Mlsame(jobs, "Y");
    bool wntvec = Mlsame(jobz, "V");
    bool wntvcf = Mlsame(jobz, "F");
    bool wntvcq = Mlsame(jobz, "Q");
    bool wntref = Mlsame(jobf, "R");
    bool wntex = Mlsame(jobf, "E");
    bool wantq = Mlsame(jobq, "Q");
    bool wnttrf = Mlsame(jobt, "R");
    INTEGER minmn = min(m, n);
    info = 0;
    bool lquery = ((lwork == -1) || (liwork == -1));
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if (!(sccolx || sccoly || Mlsame(jobs, "N"))) {
        info = -1;
    } else if (!(wntvec || wntvcf || wntvcq || Mlsame(jobz, "N"))) {
        info = -2;
    } else if (!(wntres || Mlsame(jobr, "N")) || (wntres && Mlsame(jobz, "N"))) {
        info = -3;
    } else if (!(wantq || Mlsame(jobq, "N"))) {
        info = -4;
    } else if (!(wnttrf || Mlsame(jobt, "N"))) {
        info = -5;
    } else if (!(wntref || wntex || Mlsame(jobf, "N"))) {
        info = -6;
    } else if (!((whtsvd == 1) || (whtsvd == 2) || (whtsvd == 3) || (whtsvd == 4))) {
        info = -7;
    } else if (m < 0) {
        info = -8;
    } else if ((n < 0) || (n > m + 1)) {
        info = -9;
    } else if (ldf < m) {
        info = -11;
    } else if (ldx < minmn) {
        info = -13;
    } else if (ldy < minmn) {
        info = -15;
    } else if (!((nrnk == -2) || (nrnk == -1) || ((nrnk >= 1) && (nrnk <= n)))) {
        info = -16;
    } else if ((tol < zero) || (tol >= one)) {
        info = -17;
    } else if (ldz < m) {
        info = -22;
    } else if ((wntref || wntex) && (ldb < minmn)) {
        info = -25;
    } else if (ldv < n - 1) {
        info = -27;
    } else if (lds < n - 1) {
        info = -29;
    }
    //
    char jobvl;
    if (wntvec || wntvcf || wntvcq) {
        jobvl = 'V';
    } else {
        jobvl = 'N';
    }
    INTEGER mlwqr = 0;
    INTEGER mlwork = 0;
    REAL rdummy[2];
    INTEGER info1 = 0;
    INTEGER olwqr = 0;
    INTEGER olwork = 0;
    INTEGER mlwdmd = 0;
    INTEGER iminwr = 0;
    INTEGER olwdmd = 0;
    INTEGER mlwmqr = 0;
    INTEGER olwmqr = 0;
    INTEGER mlwgqr = 0;
    INTEGER olwgqr = 0;
    if (info == 0) {
        // Compute the minimal and the optimal workspace
        // requirements. Simulate running the code and
        // determine minimal and optimal sizes of the
        // workspace at any moment of the run.
        if ((n == 0) || (n == 1)) {
            // All output except K is void. INFO=1 signals
            // the void input. In case of a workspace query,
            // the minimal workspace lengths are returned.
            if (lquery) {
                iwork[1 - 1] = 1;
                work[1 - 1] = 2.0;
                work[2 - 1] = 2.0;
            } else {
                k = 0;
            }
            info = 1;
            return;
        }
        // Minimal workspace length for Rgeqrf.
        mlwqr = max((INTEGER)1, n);
        mlwork = minmn + mlwqr;
        if (lquery) {
            Rgeqrf(m, n, f, ldf, work, rdummy, -1, info1);
            olwqr = castINTEGER(rdummy[1 - 1]);
            olwork = min(m, n) + olwqr;
        }
        Rgedmd(jobs, &jobvl, jobr, jobf, whtsvd, minmn, n - 1, x, ldx, y, ldy, nrnk, tol, k, reig, imeig, z, ldz, res, b, ldb, v, ldv, s, lds, work, -1, iwork, liwork, info1);
        mlwdmd = castINTEGER(work[1 - 1]);
        mlwork = max(mlwork, minmn + mlwdmd);
        iminwr = iwork[1 - 1];
        if (lquery) {
            olwdmd = castINTEGER(work[2 - 1]);
            olwork = max(olwork, minmn + olwdmd);
        }
        if (wntvec || wntvcf) {
            mlwmqr = max((INTEGER)1, n);
            mlwork = max(mlwork, minmn + n - 1 + mlwmqr);
            if (lquery) {
                Rormqr("L", "N", m, n, minmn, f, ldf, work, z, ldz, work, -1, info1);
                olwmqr = castINTEGER(work[1 - 1]);
                olwork = max(olwork, minmn + n - 1 + olwmqr);
            }
        }
        if (wantq) {
            mlwgqr = n;
            mlwork = max(mlwork, minmn + n - 1 + mlwgqr);
            if (lquery) {
                Rorgqr(m, minmn, minmn, f, ldf, work, work, -1, info1);
                olwgqr = castINTEGER(work[1 - 1]);
                olwork = max(olwork, minmn + n - 1 + olwgqr);
            }
        }
        iminwr = max((INTEGER)1, iminwr);
        mlwork = max((INTEGER)2, mlwork);
        if (lwork < mlwork && (!lquery)) {
            info = -31;
        }
        if (liwork < iminwr && (!lquery)) {
            info = -33;
        }
    }
    if (info != 0) {
        Mxerbla("Rgedmdq", -info);
        return;
    } else if (lquery) {
        // Return minimal and optimal workspace sizes
        iwork[1 - 1] = iminwr;
        work[1 - 1] = mlwork;
        work[2 - 1] = olwork;
        return;
    }
    // .....
    // Initial QR factorization that is used to represent the
    // snapshots as elements of lower dimensional subspace.
    // For large scale computation with M >>N , at this place
    // one can use an out of core QRF.
    //
    Rgeqrf(m, n, f, ldf, work, &work[(minmn + 1) - 1], lwork - minmn, info1);
    //
    // Define X and Y as the snapshots representations in the
    // orthogonal basis computed in the QR factorization.
    // X corresponds to the leading N-1 and Y to the trailing
    // N-1 snapshots.
    Rlaset("L", minmn, n - 1, zero, zero, x, ldx);
    Rlacpy("U", minmn, n - 1, f, ldf, x, ldx);
    Rlacpy("A", minmn, n - 1, &f[(2 - 1) * ldf], ldf, y, ldy);
    if (m >= 3) {
        Rlaset("L", minmn - 2, n - 2, zero, zero, &y[(3 - 1)], ldy);
    }
    //
    // Compute the DMD of the projected snapshot pairs (X,Y)
    Rgedmd(jobs, &jobvl, jobr, jobf, whtsvd, minmn, n - 1, x, ldx, y, ldy, nrnk, tol, k, reig, imeig, z, ldz, res, b, ldb, v, ldv, s, lds, &work[(minmn + 1) - 1], lwork - minmn, iwork, liwork, info1);
    if (info1 == 2 || info1 == 3) {
        // Return with error code. See Rgedmd for details.
        info = info1;
        return;
    } else {
        info = info1;
    }
    //
    // The Ritz vectors (Koopman modes) can be explicitly
    // formed or returned in factored form.
    if (wntvec) {
        // Compute the eigenvectors explicitly.
        if (m > minmn) {
            Rlaset("A", m - minmn, k, zero, zero, &z[((minmn + 1) - 1)], ldz);
        }
        Rormqr("L", "N", m, k, minmn, f, ldf, work, z, ldz, &work[(minmn + n) - 1], lwork - (minmn + n - 1), info1);
    } else if (wntvcf) {
        // Return the Ritz vectors (eigenvectors) in factored
        // form Z*V, where Z contains orthonormal matrix (the
        // product of Q from the initial QR factorization and
        // the SVD/POD_basis returned by Rgedmd in X) and the
        // second factor (the eigenvectors of the Rayleigh
        // quotient) is in the array V, as returned by Rgedmd.
        Rlacpy("A", n, k, x, ldx, z, ldz);
        if (m > n) {
            Rlaset("A", m - n, k, zero, zero, &z[((n + 1) - 1)], ldz);
        }
        Rormqr("L", "N", m, k, minmn, f, ldf, work, z, ldz, &work[(minmn + n) - 1], lwork - (minmn + n - 1), info1);
    }
    //
    // Some optional output variables:
    //
    // The upper triangular factor R in the initial QR
    // factorization is optionally returned in the array Y.
    // This is useful if this call to Rgedmdq is to be
    // followed by a streaming DMD that is implemented in a
    // QR compressed form.
    // Return the upper triangular R in Y
    if (wnttrf) {
        Rlaset("A", minmn, n, zero, zero, y, ldy);
        Rlacpy("U", minmn, n, f, ldf, y, ldy);
    }
    //
    // The orthonormal/orthogonal factor Q in the initial QR
    // factorization is optionally returned in the array F.
    // Same as with the triangular factor above, this is
    // useful in a streaming DMD.
    // Q overwrites F
    if (wantq) {
        Rorgqr(m, minmn, minmn, f, ldf, work, &work[(minmn + n) - 1], lwork - (minmn + n - 1), info1);
    }
    //
}
