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

void Rggsvd3(const char *jobu, const char *jobv, const char *jobq, INTEGER const &m, INTEGER const &n, INTEGER const &p, INTEGER const &k, INTEGER const &l, REAL *a, INTEGER const &lda, REAL *b, INTEGER const &ldb, REAL *alpha, REAL *beta, REAL *u, INTEGER const &ldu, REAL *v, INTEGER const &ldv, REAL *q, INTEGER const &ldq, REAL *work, INTEGER const &lwork, arr_ref<INTEGER> iwork, INTEGER &info) {
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
    //     Decode and test the input parameters
    //
    bool wantu = Mlsame(jobu, "U");
    bool wantv = Mlsame(jobv, "V");
    bool wantq = Mlsame(jobq, "Q");
    bool lquery = (lwork == -1);
    INTEGER lwkopt = 1;
    //
    //     Test the input arguments
    //
    info = 0;
    if (!(wantu || Mlsame(jobu, "N"))) {
        info = -1;
    } else if (!(wantv || Mlsame(jobv, "N"))) {
        info = -2;
    } else if (!(wantq || Mlsame(jobq, "N"))) {
        info = -3;
    } else if (m < 0) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (p < 0) {
        info = -6;
    } else if (lda < max((INTEGER)1, m)) {
        info = -10;
    } else if (ldb < max((INTEGER)1, p)) {
        info = -12;
    } else if (ldu < 1 || (wantu && ldu < m)) {
        info = -16;
    } else if (ldv < 1 || (wantv && ldv < p)) {
        info = -18;
    } else if (ldq < 1 || (wantq && ldq < n)) {
        info = -20;
    } else if (lwork < 1 && !lquery) {
        info = -24;
    }
    //
    //     Compute workspace
    //
    REAL tola = 0.0;
    REAL tolb = 0.0;
    if (info == 0) {
        Rggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, work, work, -1, info);
        lwkopt = n + INTEGER(work[1 - 1]);
        lwkopt = max(2 * n, lwkopt);
        lwkopt = max((INTEGER)1, lwkopt);
        work[1 - 1] = lwkopt.real();
    }
    //
    if (info != 0) {
        Mxerbla("Rggsvd3", -info);
        return;
    }
    if (lquery) {
        return;
    }
    //
    //     Compute the Frobenius norm of matrices A and B
    //
    REAL anorm = Rlange[("1" - 1) + (m - 1) * ldRlange];
    REAL bnorm = Rlange[("1" - 1) + (p - 1) * ldRlange];
    //
    //     Get machine precision and set up threshold for determining
    //     the effective numerical rank of the matrices A and B.
    //
    REAL ulp = dlamch("Precision");
    REAL unfl = dlamch("Safe Minimum");
    tola = max(m, n) * max(anorm, unfl) * ulp;
    tolb = max(p, n) * max(bnorm, unfl) * ulp;
    //
    //     Preprocessing
    //
    Rggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, work, work[(n + 1) - 1], lwork - n, info);
    //
    //     Compute the GSVD of two upper "triangular" matrices
    //
    INTEGER ncycle = 0;
    Rtgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work, ncycle, info);
    //
    //     Sort the singular values and store the pivot indices in IWORK
    //     Copy ALPHA to WORK, then sort ALPHA in WORK
    //
    Rcopy(n, alpha, 1, work, 1);
    INTEGER ibnd = min(l, m - k);
    INTEGER i = 0;
    INTEGER isub = 0;
    REAL smax = 0.0;
    INTEGER j = 0;
    REAL temp = 0.0;
    for (i = 1; i <= ibnd; i = i + 1) {
        //
        //        Scan for largest ALPHA(K+I)
        //
        isub = i;
        smax = work[(k + i) - 1];
        for (j = i + 1; j <= ibnd; j = j + 1) {
            temp = work[(k + j) - 1];
            if (temp > smax) {
                isub = j;
                smax = temp;
            }
        }
        if (isub != i) {
            work[(k + isub) - 1] = work[(k + i) - 1];
            work[(k + i) - 1] = smax;
            iwork[(k + i) - 1] = k + isub;
        } else {
            iwork[(k + i) - 1] = k + i;
        }
    }
    //
    work[1 - 1] = lwkopt.real();
    //
    //     End of Rggsvd3
    //
}
