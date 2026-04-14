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

// Derived from LAPACK routine ZLAQZ2.
// Original LAPACK authors:
//   Thijs Steel, KU Leuven

#include <mpblas.h>
#include <mplapack.h>

void Claqz2(bool const ilschur, bool const ilq, bool const ilz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, INTEGER const nw, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *q, INTEGER const ldq, COMPLEX *z, INTEGER const ldz, INTEGER &ns, INTEGER &nd, COMPLEX *alpha, COMPLEX *beta, COMPLEX *qc, INTEGER const ldqc, COMPLEX *zc, INTEGER const ldzc, COMPLEX *work, INTEGER const lwork, REAL *rwork, INTEGER const rec, INTEGER &info) {
    //
    // Arguments
    //
    // Parameters
    //
    // Local Scalars
    //
    // External Functions
    //
    info = 0;
    //
    // Set up deflation window
    INTEGER jw = min(nw, ihi - ilo + 1);
    INTEGER kwtop = ihi - jw + 1;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    COMPLEX s = 0.0;
    if (kwtop == ilo) {
        s = czero;
    } else {
        s = a[(kwtop - 1) + ((kwtop - 1) - 1) * lda];
    }
    //
    // Determine required workspace
    INTEGER ifst = 1;
    INTEGER ilst = jw;
    INTEGER qz_small_info = 0;
    Claqz0("S", "V", "V", jw, 1, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb, alpha, beta, qc, ldqc, zc, ldzc, work, -1, rwork, rec + 1, qz_small_info);
    INTEGER lworkreq = castINTEGER(work[1 - 1].real()) + 2 * pow2(jw);
    lworkreq = max(lworkreq, n * nw, 2 * pow2(nw) + n);
    if (lwork == -1) {
        // workspace query, quick return
        work[1 - 1] = lworkreq;
        return;
    } else if (lwork < lworkreq) {
        info = -26;
    }
    //
    if (info != 0) {
        Mxerbla("Claqz2", -info);
        return;
    }
    //
    // Get machine constants
    REAL safmin = Rlamch("SAFE MINIMUM");
    const REAL one = 1.0;
    REAL safmax = one / safmin;
    REAL ulp = Rlamch("PRECISION");
    REAL smlnum = safmin * (castREAL(n) / ulp);
    //
    if (ihi == kwtop) {
        // 1 by 1 deflation window, just try a regular deflation
        alpha[kwtop - 1] = a[(kwtop - 1) + (kwtop - 1) * lda];
        beta[kwtop - 1] = b[(kwtop - 1) + (kwtop - 1) * ldb];
        ns = 1;
        nd = 0;
        if (abs(s) <= max(smlnum, ulp * abs(a[(kwtop - 1) + (kwtop - 1) * lda]))) {
            ns = 0;
            nd = 1;
            if (kwtop > ilo) {
                a[(kwtop - 1) + ((kwtop - 1) - 1) * lda] = czero;
            }
        }
    }
    //
    // Store window in case of convergence failure
    Clacpy("ALL", jw, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda, work, jw);
    Clacpy("ALL", jw, jw, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb, &work[(pow2(jw) + 1) - 1], jw);
    //
    // Transform window to real schur form
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    Claset("FULL", jw, jw, czero, cone, qc, ldqc);
    Claset("FULL", jw, jw, czero, cone, zc, ldzc);
    Claqz0("S", "V", "V", jw, 1, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb, alpha, beta, qc, ldqc, zc, ldzc, &work[(2 * pow2(jw) + 1) - 1], lwork - 2 * pow2(jw), rwork, rec + 1, qz_small_info);
    //
    if (qz_small_info != 0) {
        // Convergence failure, restore the window and exit
        nd = 0;
        ns = jw - qz_small_info;
        Clacpy("ALL", jw, jw, work, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda);
        Clacpy("ALL", jw, jw, &work[(pow2(jw) + 1) - 1], jw, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb);
        return;
    }
    //
    // Deflation detection loop
    INTEGER kwbot = 0;
    INTEGER k = 0;
    INTEGER k2 = 0;
    REAL tempr = 0.0;
    const REAL zero = 0.0;
    INTEGER ztgexc_info = 0;
    if (kwtop == ilo || s == czero) {
        kwbot = kwtop - 1;
    } else {
        kwbot = ihi;
        k = 1;
        k2 = 1;
        while (k <= jw) {
            // Try to deflate eigenvalue
            tempr = abs(a[(kwbot - 1) + (kwbot - 1) * lda]);
            if (tempr == zero) {
                tempr = abs(s);
            }
            if ((abs(s * qc[((kwbot - kwtop + 1) - 1) * ldqc])) <= max(ulp * tempr, smlnum)) {
                // Deflatable
                kwbot = kwbot - 1;
            } else {
                // Not deflatable, move out of the way
                ifst = kwbot - kwtop + 1;
                ilst = k2;
                Ctgexc(true, true, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb, qc, ldqc, zc, ldzc, ifst, ilst, ztgexc_info);
                k2++;
            }
            //
            k++;
        }
    }
    //
    // Store eigenvalues
    nd = ihi - kwbot;
    ns = jw - nd;
    k = kwtop;
    while (k <= ihi) {
        alpha[k - 1] = a[(k - 1) + (k - 1) * lda];
        beta[k - 1] = b[(k - 1) + (k - 1) * ldb];
        k++;
    }
    //
    REAL c1 = 0.0;
    COMPLEX s1 = 0.0;
    COMPLEX temp = 0.0;
    INTEGER istartm = 0;
    INTEGER istopm = 0;
    if (kwtop != ilo && s != czero) {
        // Reflect spike back, this will create optimally packed bulges
        const COMPLEX spike = a[(kwtop - 1) + ((kwtop - 1) - 1) * lda];
        for (INTEGER i_ = kwtop; i_ <= kwbot; i_++) {
            a[(i_ - 1) + ((kwtop - 1) - 1) * lda] = spike * conj(qc[(i_ - kwtop) * ldqc]);
        }
        for (k = kwbot - 1; k >= kwtop; k = k - 1) {
            Clartg(a[(k - 1) + ((kwtop - 1) - 1) * lda], a[((k + 1) - 1) + ((kwtop - 1) - 1) * lda], c1, s1, temp);
            a[(k - 1) + ((kwtop - 1) - 1) * lda] = temp;
            a[((k + 1) - 1) + ((kwtop - 1) - 1) * lda] = czero;
            k2 = max(kwtop, k - 1);
            Crot(ihi - k2 + 1, &a[(k - 1) + (k2 - 1) * lda], lda, &a[((k + 1) - 1) + (k2 - 1) * lda], lda, c1, s1);
            Crot(ihi - (k - 1) + 1, &b[(k - 1) + ((k - 1) - 1) * ldb], ldb, &b[((k + 1) - 1) + ((k - 1) - 1) * ldb], ldb, c1, s1);
            Crot(jw, &qc[((k - kwtop + 1) - 1) * ldqc], 1, &qc[((k + 1 - kwtop + 1) - 1) * ldqc], 1, c1, conj(s1));
        }
        //
        // Chase bulges down
        istartm = kwtop;
        istopm = ihi;
        k = kwbot - 1;
        while (k >= kwtop) {
            //
            // Move bulge down and remove it
            for (k2 = k; k2 <= kwbot - 1; k2 = k2 + 1) {
                Claqz1(true, true, k2, kwtop, kwtop + jw - 1, kwbot, a, lda, b, ldb, jw, kwtop, qc, ldqc, jw, kwtop, zc, ldzc);
            }
            //
            k = k - 1;
        }
        //
    }
    //
    // Apply Qc and Zc to rest of the matrix
    if (ilschur) {
        istartm = 1;
        istopm = n;
    } else {
        istartm = ilo;
        istopm = ihi;
    }
    //
    if (istopm - ihi > 0) {
        Cgemm("C", "N", jw, istopm - ihi, jw, cone, qc, ldqc, &a[(kwtop - 1) + ((ihi + 1) - 1) * lda], lda, czero, work, jw);
        Clacpy("ALL", jw, istopm - ihi, work, jw, &a[(kwtop - 1) + ((ihi + 1) - 1) * lda], lda);
        Cgemm("C", "N", jw, istopm - ihi, jw, cone, qc, ldqc, &b[(kwtop - 1) + ((ihi + 1) - 1) * ldb], ldb, czero, work, jw);
        Clacpy("ALL", jw, istopm - ihi, work, jw, &b[(kwtop - 1) + ((ihi + 1) - 1) * ldb], ldb);
    }
    if (ilq) {
        Cgemm("N", "N", n, jw, jw, cone, &q[(kwtop - 1) * ldq], ldq, qc, ldqc, czero, work, n);
        Clacpy("ALL", n, jw, work, n, &q[(kwtop - 1) * ldq], ldq);
    }
    //
    if (kwtop - 1 - istartm + 1 > 0) {
        Cgemm("N", "N", kwtop - istartm, jw, jw, cone, &a[(istartm - 1) + (kwtop - 1) * lda], lda, zc, ldzc, czero, work, kwtop - istartm);
        Clacpy("ALL", kwtop - istartm, jw, work, kwtop - istartm, &a[(istartm - 1) + (kwtop - 1) * lda], lda);
        Cgemm("N", "N", kwtop - istartm, jw, jw, cone, &b[(istartm - 1) + (kwtop - 1) * ldb], ldb, zc, ldzc, czero, work, kwtop - istartm);
        Clacpy("ALL", kwtop - istartm, jw, work, kwtop - istartm, &b[(istartm - 1) + (kwtop - 1) * ldb], ldb);
    }
    if (ilz) {
        Cgemm("N", "N", n, jw, jw, cone, &z[(kwtop - 1) * ldz], ldz, zc, ldzc, czero, work, n);
        Clacpy("ALL", n, jw, work, n, &z[(kwtop - 1) * ldz], ldz);
    }
    //
}
