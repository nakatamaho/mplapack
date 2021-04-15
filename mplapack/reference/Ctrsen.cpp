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

void Ctrsen(const char *job, const char *compq, bool *select, INTEGER const n, COMPLEX *t, INTEGER const ldt, COMPLEX *q, INTEGER const ldq, COMPLEX *w, INTEGER &m, REAL &s, REAL &sep, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
    bool wantbh = false;
    bool wants = false;
    bool wantsp = false;
    bool wantq = false;
    INTEGER k = 0;
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER nn = 0;
    bool lquery = false;
    INTEGER lwmin = 0;
    const REAL one = 1.0;
    REAL rwork[1];
    INTEGER ks = 0;
    INTEGER ierr = 0;
    REAL scale = 0.0;
    REAL rnorm = 0.0;
    const REAL zero = 0.0;
    REAL est = 0.0;
    INTEGER kase = 0;
    INTEGER isave[3];
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
    //     .. Local Arrays ..
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
    wantbh = Mlsame(job, "B");
    wants = Mlsame(job, "E") || wantbh;
    wantsp = Mlsame(job, "V") || wantbh;
    wantq = Mlsame(compq, "V");
    //
    //     Set M to the number of selected eigenvalues.
    //
    m = 0;
    for (k = 1; k <= n; k = k + 1) {
        if (select[k - 1]) {
            m++;
        }
    }
    //
    n1 = m;
    n2 = n - m;
    nn = n1 * n2;
    //
    info = 0;
    lquery = (lwork == -1);
    //
    if (wantsp) {
        lwmin = max((INTEGER)1, 2 * nn);
    } else if (Mlsame(job, "N")) {
        lwmin = 1;
    } else if (Mlsame(job, "E")) {
        lwmin = max((INTEGER)1, nn);
    }
    //
    if (!Mlsame(job, "N") && !wants && !wantsp) {
        info = -1;
    } else if (!Mlsame(compq, "N") && !wantq) {
        info = -2;
    } else if (n < 0) {
        info = -4;
    } else if (ldt < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldq < 1 || (wantq && ldq < n)) {
        info = -8;
    } else if (lwork < lwmin && !lquery) {
        info = -14;
    }
    //
    if (info == 0) {
        work[1 - 1] = lwmin;
    }
    //
    if (info != 0) {
        Mxerbla("Ctrsen", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == n || m == 0) {
        if (wants) {
            s = one;
        }
        if (wantsp) {
            sep = Clange("1", n, n, t, ldt, rwork);
        }
        goto statement_40;
    }
    //
    //     Collect the selected eigenvalues at the top left corner of T.
    //
    ks = 0;
    for (k = 1; k <= n; k = k + 1) {
        if (select[k - 1]) {
            ks++;
            //
            //           Swap the K-th eigenvalue to position KS.
            //
            if (k != ks) {
                Ctrexc(compq, n, t, ldt, q, ldq, k, ks, ierr);
            }
        }
    }
    //
    if (wants) {
        //
        //        Solve the Sylvester equation for R:
        //
        //           T11*R - R*T22 = scale*T12
        //
        Clacpy("F", n1, n2, &t[((n1 + 1) - 1) * ldt], ldt, work, n1);
        Ctrsyl("N", "N", -1, n1, n2, t, ldt, &t[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], ldt, work, n1, scale, ierr);
        //
        //        Estimate the reciprocal of the condition number of the cluster
        //        of eigenvalues.
        //
        rnorm = Clange("F", n1, n2, work, n1, rwork);
        if (rnorm == zero) {
            s = one;
        } else {
            s = scale / (sqrt(scale * scale / rnorm + rnorm) * sqrt(rnorm));
        }
    }
    //
    if (wantsp) {
        //
        //        Estimate sep(T11,T22).
        //
        est = zero;
        kase = 0;
    statement_30:
        Clacn2(nn, &work[(nn + 1) - 1], work, est, kase, isave);
        if (kase != 0) {
            if (kase == 1) {
                //
                //              Solve T11*R - R*T22 = scale*X.
                //
                Ctrsyl("N", "N", -1, n1, n2, t, ldt, &t[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], ldt, work, n1, scale, ierr);
            } else {
                //
                //              Solve T11**H*R - R*T22**H = scale*X.
                //
                Ctrsyl("C", "C", -1, n1, n2, t, ldt, &t[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], ldt, work, n1, scale, ierr);
            }
            goto statement_30;
        }
        //
        sep = scale / est;
    }
//
statement_40:
    //
    //     Copy reordered eigenvalues to W.
    //
    for (k = 1; k <= n; k = k + 1) {
        w[k - 1] = t[(k - 1) + (k - 1) * ldt];
    }
    //
    work[1 - 1] = lwmin;
    //
    //     End of Ctrsen
    //
}
