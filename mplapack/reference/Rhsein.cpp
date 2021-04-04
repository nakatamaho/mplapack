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

void Rhsein(const char *side, const char *eigsrc, const char *initv, arr_ref<bool> select, INTEGER const &n, REAL *h, INTEGER const &ldh, REAL *wr, REAL *wi, REAL *vl, INTEGER const &ldvl, REAL *vr, INTEGER const &ldvr, INTEGER const &mm, INTEGER &m, REAL *work, arr_ref<INTEGER> ifaill, arr_ref<INTEGER> ifailr, INTEGER &info) {
    bool bothv = false;
    bool rightv = false;
    bool leftv = false;
    bool fromqr = false;
    bool noinit = false;
    bool pair = false;
    INTEGER k = 0;
    const REAL zero = 0.0;
    REAL unfl = 0.0;
    REAL ulp = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    INTEGER ldwork = 0;
    INTEGER kl = 0;
    INTEGER kln = 0;
    INTEGER kr = 0;
    INTEGER ksr = 0;
    INTEGER i = 0;
    REAL hnorm = 0.0;
    REAL eps3 = 0.0;
    REAL wkr = 0.0;
    REAL wki = 0.0;
    INTEGER ksi = 0;
    INTEGER iinfo = 0;
    //
    //  -- LAPACK computational routine --
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Decode and test the input parameters.
    //
    bothv = Mlsame(side, "B");
    rightv = Mlsame(side, "R") || bothv;
    leftv = Mlsame(side, "L") || bothv;
    //
    fromqr = Mlsame(eigsrc, "Q");
    //
    noinit = Mlsame(initv, "N");
    //
    //     Set M to the number of columns required to store the selected
    //     eigenvectors, and standardize the array SELECT.
    //
    m = 0;
    pair = false;
    for (k = 1; k <= n; k = k + 1) {
        if (pair) {
            pair = false;
            select[k - 1] = false;
        } else {
            if (wi[k - 1] == zero) {
                if (select[k - 1]) {
                    m++;
                }
            } else {
                pair = true;
                if (select[k - 1] || select[(k + 1) - 1]) {
                    select[k - 1] = true;
                    m += 2;
                }
            }
        }
    }
    //
    info = 0;
    if (!rightv && !leftv) {
        info = -1;
    } else if (!fromqr && !Mlsame(eigsrc, "N")) {
        info = -2;
    } else if (!noinit && !Mlsame(initv, "U")) {
        info = -3;
    } else if (n < 0) {
        info = -5;
    } else if (ldh < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldvl < 1 || (leftv && ldvl < n)) {
        info = -11;
    } else if (ldvr < 1 || (rightv && ldvr < n)) {
        info = -13;
    } else if (mm < m) {
        info = -14;
    }
    if (info != 0) {
        Mxerbla("Rhsein", -info);
        return;
    }
    //
    //     Quick return if possible.
    //
    if (n == 0) {
        return;
    }
    //
    //     Set machine-dependent constants.
    //
    unfl = dlamch("Safe minimum");
    ulp = dlamch("Precision");
    smlnum = unfl * (n / ulp);
    bignum = (one - ulp) / smlnum;
    //
    ldwork = n + 1;
    //
    kl = 1;
    kln = 0;
    if (fromqr) {
        kr = 0;
    } else {
        kr = n;
    }
    ksr = 1;
    //
    for (k = 1; k <= n; k = k + 1) {
        if (select[k - 1]) {
            //
            //           Compute eigenvector(s) corresponding to W(K).
            //
            if (fromqr) {
                //
                //              If affiliation of eigenvalues is known, check whether
                //              the matrix splits.
                //
                //              Determine KL and KR such that 1 <= KL <= K <= KR <= N
                //              and H(KL,KL-1) and H(KR+1,KR) are zero (or KL = 1 or
                //              KR = N).
                //
                //              Then inverse iteration can be performed with the
                //              submatrix H(KL:N,KL:N) for a left eigenvector, and with
                //              the submatrix H(1:KR,1:KR) for a right eigenvector.
                //
                for (i = k; i >= kl + 1; i = i - 1) {
                    if (h[(i - 1) + ((i - 1) - 1) * ldh] == zero) {
                        goto statement_30;
                    }
                }
            statement_30:
                kl = i;
                if (k > kr) {
                    for (i = k; i <= n - 1; i = i + 1) {
                        if (h[((i + 1) - 1) + (i - 1) * ldh] == zero) {
                            goto statement_50;
                        }
                    }
                statement_50:
                    kr = i;
                }
            }
            //
            if (kl != kln) {
                kln = kl;
                //
                //              Compute infinity-norm of submatrix H(KL:KR,KL:KR) if it
                //              has not ben computed before.
                //
                hnorm = Rlanhs[("I" - 1) + ((kr - kl + 1) - 1) * ldRlanhs];
                if (Risnan(hnorm)) {
                    info = -6;
                    return;
                } else if (hnorm > zero) {
                    eps3 = hnorm * ulp;
                } else {
                    eps3 = smlnum;
                }
            }
            //
            //           Perturb eigenvalue if it is close to any previous
            //           selected eigenvalues affiliated to the submatrix
            //           H(KL:KR,KL:KR). Close roots are modified by EPS3.
            //
            wkr = wr[k - 1];
            wki = wi[k - 1];
        statement_60:
            for (i = k - 1; i >= kl; i = i - 1) {
                if (select[i - 1] && abs(wr[i - 1] - wkr) + abs(wi[i - 1] - wki) < eps3) {
                    wkr += eps3;
                    goto statement_60;
                }
            }
            wr[k - 1] = wkr;
            //
            pair = wki != zero;
            if (pair) {
                ksi = ksr + 1;
            } else {
                ksi = ksr;
            }
            if (leftv) {
                //
                //              Compute left eigenvector.
                //
                Rlaein(false, noinit, n - kl + 1, h[(kl - 1) + (kl - 1) * ldh], ldh, wkr, wki, vl[(kl - 1) + (ksr - 1) * ldvl], vl[(kl - 1) + (ksi - 1) * ldvl], work, ldwork, work[(n * n + n + 1) - 1], eps3, smlnum, bignum, iinfo);
                if (iinfo > 0) {
                    if (pair) {
                        info += 2;
                    } else {
                        info++;
                    }
                    ifaill[ksr - 1] = k;
                    ifaill[ksi - 1] = k;
                } else {
                    ifaill[ksr - 1] = 0;
                    ifaill[ksi - 1] = 0;
                }
                for (i = 1; i <= kl - 1; i = i + 1) {
                    vl[(i - 1) + (ksr - 1) * ldvl] = zero;
                }
                if (pair) {
                    for (i = 1; i <= kl - 1; i = i + 1) {
                        vl[(i - 1) + (ksi - 1) * ldvl] = zero;
                    }
                }
            }
            if (rightv) {
                //
                //              Compute right eigenvector.
                //
                Rlaein(true, noinit, kr, h, ldh, wkr, wki, vr[(ksr - 1) * ldvr], vr[(ksi - 1) * ldvr], work, ldwork, work[(n * n + n + 1) - 1], eps3, smlnum, bignum, iinfo);
                if (iinfo > 0) {
                    if (pair) {
                        info += 2;
                    } else {
                        info++;
                    }
                    ifailr[ksr - 1] = k;
                    ifailr[ksi - 1] = k;
                } else {
                    ifailr[ksr - 1] = 0;
                    ifailr[ksi - 1] = 0;
                }
                for (i = kr + 1; i <= n; i = i + 1) {
                    vr[(i - 1) + (ksr - 1) * ldvr] = zero;
                }
                if (pair) {
                    for (i = kr + 1; i <= n; i = i + 1) {
                        vr[(i - 1) + (ksi - 1) * ldvr] = zero;
                    }
                }
            }
            //
            if (pair) {
                ksr += 2;
            } else {
                ksr++;
            }
        }
    }
    //
    //     End of Rhsein
    //
}
