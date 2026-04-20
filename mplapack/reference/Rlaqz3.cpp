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

// Derived from LAPACK routine DLAQZ3.
// Original LAPACK authors:
//   Thijs Steel, KU Leuven

#include <mpblas.h>
#include <mplapack.h>

void Rlaqz3(bool const ilschur, bool const ilq, bool const ilz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, INTEGER const nw, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *q, INTEGER const ldq, REAL *z, INTEGER const ldz, INTEGER &ns, INTEGER &nd, REAL *alphar, REAL *alphai, REAL *beta, REAL *qc, INTEGER const ldqc, REAL *zc, INTEGER const ldzc, REAL *work, INTEGER const lwork, INTEGER const rec, INTEGER &info) {
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
    const REAL zero = 0.0;
    REAL s = 0.0;
    if (kwtop == ilo) {
        s = zero;
    } else {
        s = a[(kwtop - 1) + ((kwtop - 1) - 1) * lda];
    }
    //
    // Determine required workspace
    INTEGER ifst = 1;
    INTEGER ilst = jw;
    INTEGER dtgexc_info = 0;
    Rtgexc(true, true, jw, a, lda, b, ldb, qc, ldqc, zc, ldzc, ifst, ilst, work, -1, dtgexc_info);
    INTEGER lworkreq = castINTEGER(work[1 - 1]);
    INTEGER qz_small_info = 0;
    Rlaqz0("S", "V", "V", jw, 1, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb, alphar, alphai, beta, qc, ldqc, zc, ldzc, work, -1, rec + 1, qz_small_info);
    lworkreq = max(lworkreq, castINTEGER(work[1 - 1]) + 2 * pow2(jw));
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
        Mxerbla("Rlaqz3", -info);
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
        alphar[kwtop - 1] = a[(kwtop - 1) + (kwtop - 1) * lda];
        alphai[kwtop - 1] = zero;
        beta[kwtop - 1] = b[(kwtop - 1) + (kwtop - 1) * ldb];
        ns = 1;
        nd = 0;
        if (abs(s) <= max(smlnum, ulp * abs(a[(kwtop - 1) + (kwtop - 1) * lda]))) {
            ns = 0;
            nd = 1;
            if (kwtop > ilo) {
                a[(kwtop - 1) + ((kwtop - 1) - 1) * lda] = zero;
            }
        }
    }
    //
    // Store window in case of convergence failure
    Rlacpy("ALL", jw, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda, work, jw);
    Rlacpy("ALL", jw, jw, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb, &work[(pow2(jw) + 1) - 1], jw);
    //
    // Transform window to real schur form
    Rlaset("FULL", jw, jw, zero, one, qc, ldqc);
    Rlaset("FULL", jw, jw, zero, one, zc, ldzc);
    Rlaqz0("S", "V", "V", jw, 1, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb, alphar, alphai, beta, qc, ldqc, zc, ldzc, &work[(2 * pow2(jw) + 1) - 1], lwork - 2 * pow2(jw), rec + 1, qz_small_info);
    //
    if (qz_small_info != 0) {
        // Convergence failure, restore the window and exit
        nd = 0;
        ns = jw - qz_small_info;
        Rlacpy("ALL", jw, jw, work, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda);
        Rlacpy("ALL", jw, jw, &work[(pow2(jw) + 1) - 1], jw, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb);
        return;
    }
    //
    // Deflation detection loop
    INTEGER kwbot = 0;
    INTEGER k = 0;
    INTEGER k2 = 0;
    bool bulge = false;
    REAL temp = 0.0;
    if (kwtop == ilo || s == zero) {
        kwbot = kwtop - 1;
    } else {
        kwbot = ihi;
        k = 1;
        k2 = 1;
        while (k <= jw) {
            bulge = false;
            if (kwbot - kwtop + 1 >= 2) {
                bulge = a[(kwbot - 1) + ((kwbot - 1) - 1) * lda] != zero;
            }
            if (bulge) {
                //
                // Try to deflate complex conjugate eigenvalue pair
                temp = abs(a[(kwbot - 1) + (kwbot - 1) * lda]) + sqrt(abs(a[(kwbot - 1) + ((kwbot - 1) - 1) * lda])) * sqrt(abs(a[((kwbot - 1) - 1) + (kwbot - 1) * lda]));
                if (temp == zero) {
                    temp = abs(s);
                }
                if (max(abs(s * qc[((kwbot - kwtop) - 1) * ldqc]), abs(s * qc[((kwbot - kwtop + 1) - 1) * ldqc])) <= max(smlnum, ulp * temp)) {
                    // Deflatable
                    kwbot = kwbot - 2;
                } else {
                    // Not deflatable, move out of the way
                    ifst = kwbot - kwtop + 1;
                    ilst = k2;
                    Rtgexc(true, true, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb, qc, ldqc, zc, ldzc, ifst, ilst, work, lwork, dtgexc_info);
                    k2 += 2;
                }
                k += 2;
            } else {
                //
                // Try to deflate real eigenvalue
                temp = abs(a[(kwbot - 1) + (kwbot - 1) * lda]);
                if (temp == zero) {
                    temp = abs(s);
                }
                if ((abs(s * qc[((kwbot - kwtop + 1) - 1) * ldqc])) <= max(ulp * temp, smlnum)) {
                    // Deflatable
                    kwbot = kwbot - 1;
                } else {
                    // Not deflatable, move out of the way
                    ifst = kwbot - kwtop + 1;
                    ilst = k2;
                    Rtgexc(true, true, jw, &a[(kwtop - 1) + (kwtop - 1) * lda], lda, &b[(kwtop - 1) + (kwtop - 1) * ldb], ldb, qc, ldqc, zc, ldzc, ifst, ilst, work, lwork, dtgexc_info);
                    k2++;
                }
                //
                k++;
                //
            }
        }
    }
    //
    // Store eigenvalues
    nd = ihi - kwbot;
    ns = jw - nd;
    k = kwtop;
    while (k <= ihi) {
        bulge = false;
        if (k < ihi) {
            if (a[((k + 1) - 1) + (k - 1) * lda] != zero) {
                bulge = true;
            }
        }
        if (bulge) {
            // 2x2 eigenvalue block
            Rlag2(&a[(k - 1) + (k - 1) * lda], lda, &b[(k - 1) + (k - 1) * ldb], ldb, safmin, beta[k - 1], beta[(k + 1) - 1], alphar[k - 1], alphar[(k + 1) - 1], alphai[k - 1]);
            alphai[(k + 1) - 1] = -alphai[k - 1];
            k += 2;
        } else {
            // 1x1 eigenvalue block
            alphar[k - 1] = a[(k - 1) + (k - 1) * lda];
            alphai[k - 1] = zero;
            beta[k - 1] = b[(k - 1) + (k - 1) * ldb];
            k++;
        }
    }
    //
    REAL c1 = 0.0;
    REAL s1 = 0.0;
    INTEGER istartm = 0;
    INTEGER istopm = 0;
    if (kwtop != ilo && s != zero) {
        // Reflect spike back, this will create optimally packed bulges
        const REAL spike = a[(kwtop - 1) + ((kwtop - 1) - 1) * lda];
        for (INTEGER i_ = kwtop; i_ <= kwbot; i_++) {
            a[(i_ - 1) + ((kwtop - 1) - 1) * lda] = spike * qc[(i_ - kwtop) * ldqc];
        }
        for (k = kwbot - 1; k >= kwtop; k = k - 1) {
            Rlartg(a[(k - 1) + ((kwtop - 1) - 1) * lda], a[((k + 1) - 1) + ((kwtop - 1) - 1) * lda], c1, s1, temp);
            a[(k - 1) + ((kwtop - 1) - 1) * lda] = temp;
            a[((k + 1) - 1) + ((kwtop - 1) - 1) * lda] = zero;
            k2 = max(kwtop, k - 1);
            Rrot(ihi - k2 + 1, &a[(k - 1) + (k2 - 1) * lda], lda, &a[((k + 1) - 1) + (k2 - 1) * lda], lda, c1, s1);
            Rrot(ihi - (k - 1) + 1, &b[(k - 1) + ((k - 1) - 1) * ldb], ldb, &b[((k + 1) - 1) + ((k - 1) - 1) * ldb], ldb, c1, s1);
            Rrot(jw, &qc[((k - kwtop + 1) - 1) * ldqc], 1, &qc[((k + 1 - kwtop + 1) - 1) * ldqc], 1, c1, s1);
        }
        //
        // Chase bulges down
        istartm = kwtop;
        istopm = ihi;
        k = kwbot - 1;
        while (k >= kwtop) {
            if ((k >= kwtop + 1) && a[((k + 1) - 1) + ((k - 1) - 1) * lda] != zero) {
                //
                // Move double pole block down and remove it
                for (k2 = k - 1; k2 <= kwbot - 2; k2 = k2 + 1) {
                    Rlaqz2(true, true, k2, kwtop, kwtop + jw - 1, kwbot, a, lda, b, ldb, jw, kwtop, qc, ldqc, jw, kwtop, zc, ldzc);
                }
                //
                k = k - 2;
            } else {
                //
                // k points to single shift
                for (k2 = k; k2 <= kwbot - 2; k2 = k2 + 1) {
                    //
                    // Move shift down
                    Rlartg(b[((k2 + 1) - 1) + ((k2 + 1) - 1) * ldb], b[((k2 + 1) - 1) + (k2 - 1) * ldb], c1, s1, temp);
                    b[((k2 + 1) - 1) + ((k2 + 1) - 1) * ldb] = temp;
                    b[((k2 + 1) - 1) + (k2 - 1) * ldb] = zero;
                    Rrot(k2 + 2 - istartm + 1, &a[(istartm - 1) + ((k2 + 1) - 1) * lda], 1, &a[(istartm - 1) + (k2 - 1) * lda], 1, c1, s1);
                    Rrot(k2 - istartm + 1, &b[(istartm - 1) + ((k2 + 1) - 1) * ldb], 1, &b[(istartm - 1) + (k2 - 1) * ldb], 1, c1, s1);
                    Rrot(jw, &zc[((k2 + 1 - kwtop + 1) - 1) * ldzc], 1, &zc[((k2 - kwtop + 1) - 1) * ldzc], 1, c1, s1);
                    //
                    Rlartg(a[((k2 + 1) - 1) + (k2 - 1) * lda], a[((k2 + 2) - 1) + (k2 - 1) * lda], c1, s1, temp);
                    a[((k2 + 1) - 1) + (k2 - 1) * lda] = temp;
                    a[((k2 + 2) - 1) + (k2 - 1) * lda] = zero;
                    Rrot(istopm - k2, &a[((k2 + 1) - 1) + ((k2 + 1) - 1) * lda], lda, &a[((k2 + 2) - 1) + ((k2 + 1) - 1) * lda], lda, c1, s1);
                    Rrot(istopm - k2, &b[((k2 + 1) - 1) + ((k2 + 1) - 1) * ldb], ldb, &b[((k2 + 2) - 1) + ((k2 + 1) - 1) * ldb], ldb, c1, s1);
                    Rrot(jw, &qc[((k2 + 1 - kwtop + 1) - 1) * ldqc], 1, &qc[((k2 + 2 - kwtop + 1) - 1) * ldqc], 1, c1, s1);
                    //
                }
                //
                // Remove the shift
                Rlartg(b[(kwbot - 1) + (kwbot - 1) * ldb], b[(kwbot - 1) + ((kwbot - 1) - 1) * ldb], c1, s1, temp);
                b[(kwbot - 1) + (kwbot - 1) * ldb] = temp;
                b[(kwbot - 1) + ((kwbot - 1) - 1) * ldb] = zero;
                Rrot(kwbot - istartm, &b[(istartm - 1) + (kwbot - 1) * ldb], 1, &b[(istartm - 1) + ((kwbot - 1) - 1) * ldb], 1, c1, s1);
                Rrot(kwbot - istartm + 1, &a[(istartm - 1) + (kwbot - 1) * lda], 1, &a[(istartm - 1) + ((kwbot - 1) - 1) * lda], 1, c1, s1);
                Rrot(jw, &zc[((kwbot - kwtop + 1) - 1) * ldzc], 1, &zc[((kwbot - 1 - kwtop + 1) - 1) * ldzc], 1, c1, s1);
                //
                k = k - 1;
            }
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
        Rgemm("T", "N", jw, istopm - ihi, jw, one, qc, ldqc, &a[(kwtop - 1) + ((ihi + 1) - 1) * lda], lda, zero, work, jw);
        Rlacpy("ALL", jw, istopm - ihi, work, jw, &a[(kwtop - 1) + ((ihi + 1) - 1) * lda], lda);
        Rgemm("T", "N", jw, istopm - ihi, jw, one, qc, ldqc, &b[(kwtop - 1) + ((ihi + 1) - 1) * ldb], ldb, zero, work, jw);
        Rlacpy("ALL", jw, istopm - ihi, work, jw, &b[(kwtop - 1) + ((ihi + 1) - 1) * ldb], ldb);
    }
    if (ilq) {
        Rgemm("N", "N", n, jw, jw, one, &q[(kwtop - 1) * ldq], ldq, qc, ldqc, zero, work, n);
        Rlacpy("ALL", n, jw, work, n, &q[(kwtop - 1) * ldq], ldq);
    }
    //
    if (kwtop - 1 - istartm + 1 > 0) {
        Rgemm("N", "N", kwtop - istartm, jw, jw, one, &a[(istartm - 1) + (kwtop - 1) * lda], lda, zc, ldzc, zero, work, kwtop - istartm);
        Rlacpy("ALL", kwtop - istartm, jw, work, kwtop - istartm, &a[(istartm - 1) + (kwtop - 1) * lda], lda);
        Rgemm("N", "N", kwtop - istartm, jw, jw, one, &b[(istartm - 1) + (kwtop - 1) * ldb], ldb, zc, ldzc, zero, work, kwtop - istartm);
        Rlacpy("ALL", kwtop - istartm, jw, work, kwtop - istartm, &b[(istartm - 1) + (kwtop - 1) * ldb], ldb);
    }
    if (ilz) {
        Rgemm("N", "N", n, jw, jw, one, &z[(kwtop - 1) * ldz], ldz, zc, ldzc, zero, work, n);
        Rlacpy("ALL", n, jw, work, n, &z[(kwtop - 1) * ldz], ldz);
    }
    //
}
