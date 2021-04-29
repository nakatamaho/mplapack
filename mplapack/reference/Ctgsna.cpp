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

void Ctgsna(const char *job, const char *howmny, bool *select, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *vl, INTEGER const ldvl, COMPLEX *vr, INTEGER const ldvr, REAL *s, REAL *dif, INTEGER const mm, INTEGER &m, COMPLEX *work, INTEGER const lwork, INTEGER *iwork, INTEGER &info) {
    bool wantbh = false;
    bool wants = false;
    bool wantdf = false;
    bool somcon = false;
    bool lquery = false;
    INTEGER k = 0;
    INTEGER lwmin = 0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    INTEGER ks = 0;
    REAL rnrm = 0.0;
    REAL lnrm = 0.0;
    const REAL zero = 0.0;
    COMPLEX yhax = 0.0;
    COMPLEX yhbx = 0.0;
    REAL cond = 0.0;
    INTEGER ifst = 0;
    INTEGER ilst = 0;
    COMPLEX dummy[1];
    COMPLEX dummy1[1];
    INTEGER ierr = 0;
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER i = 0;
    const INTEGER idifjb = 3;
    REAL scale = 0.0;
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
    //     Decode and test the input parameters
    //
    wantbh = Mlsame(job, "B");
    wants = Mlsame(job, "E") || wantbh;
    wantdf = Mlsame(job, "V") || wantbh;
    //
    somcon = Mlsame(howmny, "S");
    //
    info = 0;
    lquery = (lwork == -1);
    //
    if (!wants && !wantdf) {
        info = -1;
    } else if (!Mlsame(howmny, "A") && !somcon) {
        info = -2;
    } else if (n < 0) {
        info = -4;
    } else if (lda < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -8;
    } else if (wants && ldvl < n) {
        info = -10;
    } else if (wants && ldvr < n) {
        info = -12;
    } else {
        //
        //        Set M to the number of eigenpairs for which condition numbers
        //        are required, and test MM.
        //
        if (somcon) {
            m = 0;
            for (k = 1; k <= n; k = k + 1) {
                if (select[k - 1]) {
                    m++;
                }
            }
        } else {
            m = n;
        }
        //
        if (n == 0) {
            lwmin = 1;
        } else if (Mlsame(job, "V") || Mlsame(job, "B")) {
            lwmin = 2 * n * n;
        } else {
            lwmin = n;
        }
        work[1 - 1] = lwmin;
        //
        if (mm < m) {
            info = -15;
        } else if (lwork < lwmin && !lquery) {
            info = -18;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Ctgsna", -info);
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
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    ks = 0;
    for (k = 1; k <= n; k = k + 1) {
        //
        //        Determine whether condition numbers are required for the k-th
        //        eigenpair.
        //
        if (somcon) {
            if (!select[k - 1]) {
                goto statement_20;
            }
        }
        //
        ks++;
        //
        if (wants) {
            //
            //           Compute the reciprocal condition number of the k-th
            //           eigenvalue.
            //
            rnrm = RCnrm2(n, &vr[(ks - 1) * ldvr], 1);
            lnrm = RCnrm2(n, &vl[(ks - 1) * ldvl], 1);
            Cgemv("N", n, n, COMPLEX(one, zero), a, lda, &vr[(ks - 1) * ldvr], 1, COMPLEX(zero, zero), work, 1);
            yhax = Cdotc(n, work, 1, &vl[(ks - 1) * ldvl], 1);
            Cgemv("N", n, n, COMPLEX(one, zero), b, ldb, &vr[(ks - 1) * ldvr], 1, COMPLEX(zero, zero), work, 1);
            yhbx = Cdotc(n, work, 1, &vl[(ks - 1) * ldvl], 1);
            cond = Rlapy2(abs(yhax), abs(yhbx));
            if (cond == zero) {
                s[ks - 1] = -one;
            } else {
                s[ks - 1] = cond / (rnrm * lnrm);
            }
        }
        //
        if (wantdf) {
            if (n == 1) {
                dif[ks - 1] = Rlapy2(abs(a[(1 - 1)]), abs(b[(1 - 1)]));
            } else {
                //
                //              Estimate the reciprocal condition number of the k-th
                //              eigenvectors.
                //
                //              Copy the matrix (A, B) to the array WORK and move the
                //              (k,k)th pair to the (1,1) position.
                //
                Clacpy("Full", n, n, a, lda, work, n);
                Clacpy("Full", n, n, b, ldb, &work[(n * n + 1) - 1], n);
                ifst = k;
                ilst = 1;
                //
                Ctgexc(false, false, n, work, n, &work[(n * n + 1) - 1], n, dummy, 1, dummy1, 1, ifst, ilst, ierr);
                //
                if (ierr > 0) {
                    //
                    //                 Ill-conditioned problem - swap rejected.
                    //
                    dif[ks - 1] = zero;
                } else {
                    //
                    //                 Reordering successful, solve generalized Sylvester
                    //                 equation for R and L,
                    //                            A22 * R - L * A11 = A12
                    //                            B22 * R - L * B11 = B12,
                    //                 and compute estimate of Difl[(A11,B11), (A22, B22)].
                    //
                    n1 = 1;
                    n2 = n - n1;
                    i = n * n + 1;
                    Ctgsyl("N", idifjb, n2, n1, &work[(n * n1 + n1 + 1) - 1], n, work, n, &work[(n1 + 1) - 1], n, &work[(n * n1 + n1 + i) - 1], n, &work[i - 1], n, &work[(n1 + i) - 1], n, scale, dif[ks - 1], dummy, 1, iwork, ierr);
                }
            }
        }
    //
    statement_20:;
    }
    work[1 - 1] = lwmin;
    //
    //     End of Ctgsna
    //
}
