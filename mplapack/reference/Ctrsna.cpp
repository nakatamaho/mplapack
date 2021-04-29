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

void Ctrsna(const char *job, const char *howmny, bool *select, INTEGER const n, COMPLEX *t, INTEGER const ldt, COMPLEX *vl, INTEGER const ldvl, COMPLEX *vr, INTEGER const ldvr, REAL *s, REAL *sep, INTEGER const mm, INTEGER &m, COMPLEX *work, INTEGER const ldwork, REAL *rwork, INTEGER &info) {
    COMPLEX cdum = 0.0;
    bool wantbh = false;
    bool wants = false;
    bool wantsp = false;
    bool somcon = false;
    INTEGER j = 0;
    const REAL one = 1.0 + 0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    REAL bignum = 0.0;
    INTEGER ks = 0;
    INTEGER k = 0;
    COMPLEX prod = 0.0;
    REAL rnrm = 0.0;
    REAL lnrm = 0.0;
    arr_1d<1, COMPLEX> dummy(fill0);
    INTEGER ierr = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    REAL est = 0.0;
    INTEGER kase = 0;
    char normin = char0;
    INTEGER isave[3];
    REAL scale = 0.0;
    INTEGER ix = 0;
    REAL xnorm = 0.0;
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
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    abs1(cdum) = abs(cdum.real()) + abs(cdum.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     Decode and test the input parameters
    //
    wantbh = Mlsame(job, "B");
    wants = Mlsame(job, "E") || wantbh;
    wantsp = Mlsame(job, "V") || wantbh;
    //
    somcon = Mlsame(howmny, "S");
    //
    //     Set M to the number of eigenpairs for which condition numbers are
    //     to be computed.
    //
    if (somcon) {
        m = 0;
        for (j = 1; j <= n; j = j + 1) {
            if (select(j)) {
                m++;
            }
        }
    } else {
        m = n;
    }
    //
    info = 0;
    if (!wants && !wantsp) {
        info = -1;
    } else if (!Mlsame(howmny, "A") && !somcon) {
        info = -2;
    } else if (n < 0) {
        info = -4;
    } else if (ldt < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldvl < 1 || (wants && ldvl < n)) {
        info = -8;
    } else if (ldvr < 1 || (wants && ldvr < n)) {
        info = -10;
    } else if (mm < m) {
        info = -13;
    } else if (ldwork < 1 || (wantsp && ldwork < n)) {
        info = -16;
    }
    if (info != 0) {
        Mxerbla("Ctrsna", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    if (n == 1) {
        if (somcon) {
            if (!select(1)) {
                return;
            }
        }
        if (wants) {
            s[1 - 1] = one;
        }
        if (wantsp) {
            sep[1 - 1] = abs(t[(1 - 1)]);
        }
        return;
    }
    //
    //     Get machine constants
    //
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    ks = 1;
    for (k = 1; k <= n; k = k + 1) {
        //
        if (somcon) {
            if (!select(k)) {
                goto statement_50;
            }
        }
        //
        if (wants) {
            //
            //           Compute the reciprocal condition number of the k-th
            //           eigenvalue.
            //
            prod = Cdotc(n, vr[(ks - 1) * ldvr], 1, vl[(ks - 1) * ldvl], 1);
            rnrm = RCnrm2(n, vr[(ks - 1) * ldvr], 1);
            lnrm = RCnrm2(n, vl[(ks - 1) * ldvl], 1);
            s[ks - 1] = abs(prod) / (rnrm * lnrm);
            //
        }
        //
        if (wantsp) {
            //
            //           Estimate the reciprocal condition number of the k-th
            //           eigenvector.
            //
            //           Copy the matrix T to the array WORK and swap the k-th
            //           diagonal element to the (1,1) position.
            //
            Clacpy("Full", n, n, t, ldt, work, ldwork);
            Ctrexc("No Q", n, work, ldwork, dummy, 1, k, 1, ierr);
            //
            //           Form  C = T22 - lambda*I in WORK(2:N,2:N).
            //
            for (i = 2; i <= n; i = i + 1) {
                work[(i - 1) + (i - 1) * ldwork] = work[(i - 1) + (i - 1) * ldwork] - work[(1 - 1)];
            }
            //
            //           Estimate a lower bound for the 1-norm of inv(C**H). The 1st
            //           and (N+1)th columns of WORK are used to store work vectors.
            //
            sep[ks - 1] = zero;
            est = zero;
            kase = 0;
            normin = "N";
        statement_30:
            Clacn2(n - 1, &work[((n + 1) - 1) * ldwork], work, est, kase, isave);
            //
            if (kase != 0) {
                if (kase == 1) {
                    //
                    //                 Solve C**H*x = scale*b
                    //
                    Clatrs("Upper", "Conjugate transpose", "Nonunit", normin, n - 1, &work[(2 - 1) + (2 - 1) * ldwork], ldwork, work, scale, rwork, ierr);
                } else {
                    //
                    //                 Solve C*x = scale*b
                    //
                    Clatrs("Upper", "No transpose", "Nonunit", normin, n - 1, &work[(2 - 1) + (2 - 1) * ldwork], ldwork, work, scale, rwork, ierr);
                }
                normin = "Y";
                if (scale != one) {
                    //
                    //                 Multiply by 1/SCALE if doing so will not cause
                    //                 overflow.
                    //
                    ix = iCamax(n - 1, work, 1);
                    xnorm = abs1(work[(ix - 1)]);
                    if (scale < xnorm * smlnum || scale == zero) {
                        goto statement_40;
                    }
                    CRrscl(n, scale, work, 1);
                }
                goto statement_30;
            }
            //
            sep[ks - 1] = one / max(est, smlnum);
        }
    //
    statement_40:
        ks++;
    statement_50:;
    }
    //
    //     End of Ctrsna
    //
}
