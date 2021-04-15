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

void Rtrsna(const char *job, const char *howmny, bool *select, INTEGER const n, REAL *t, INTEGER const ldt, REAL *vl, INTEGER const ldvl, REAL *vr, INTEGER const ldvr, REAL *s, REAL *sep, INTEGER const mm, INTEGER &m, REAL *work, INTEGER const ldwork, INTEGER *iwork, INTEGER &info) {
    bool wantbh = false;
    bool wants = false;
    bool wantsp = false;
    bool somcon = false;
    bool pair = false;
    INTEGER k = 0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    REAL bignum = 0.0;
    INTEGER ks = 0;
    REAL prod = 0.0;
    REAL rnrm = 0.0;
    REAL lnrm = 0.0;
    REAL prod1 = 0.0;
    REAL prod2 = 0.0;
    REAL cond = 0.0;
    INTEGER ifst = 0;
    INTEGER ilst = 0;
    arr_1d<1, REAL> dummy(fill0);
    INTEGER ierr = 0;
    REAL scale = 0.0;
    REAL est = 0.0;
    INTEGER i = 0;
    INTEGER n2 = 0;
    INTEGER nn = 0;
    REAL mu = 0.0;
    REAL delta = 0.0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    INTEGER j = 0;
    const REAL two = 2.0e+0;
    INTEGER kase = 0;
    arr_1d<3, int> isave(fill0);
    REAL dumm = 0.0;
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
    wantsp = Mlsame(job, "V") || wantbh;
    //
    somcon = Mlsame(howmny, "S");
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
    } else {
        //
        //        Set M to the number of eigenpairs for which condition numbers
        //        are required, and test MM.
        //
        if (somcon) {
            m = 0;
            pair = false;
            for (k = 1; k <= n; k = k + 1) {
                if (pair) {
                    pair = false;
                } else {
                    if (k < n) {
                        if (t[((k + 1) - 1) + (k - 1) * ldt] == zero) {
                            if (select(k)) {
                                m++;
                            }
                        } else {
                            pair = true;
                            if (select(k) || select(k + 1)) {
                                m += 2;
                            }
                        }
                    } else {
                        if (select(n)) {
                            m++;
                        }
                    }
                }
            }
        } else {
            m = n;
        }
        //
        if (mm < m) {
            info = -13;
        } else if (ldwork < 1 || (wantsp && ldwork < n)) {
            info = -16;
        }
    }
    if (info != 0) {
        Mxerbla("Rtrsna", -info);
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
    ks = 0;
    pair = false;
    for (k = 1; k <= n; k = k + 1) {
        //
        //        Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block.
        //
        if (pair) {
            pair = false;
            goto statement_60;
        } else {
            if (k < n) {
                pair = t[((k + 1) - 1) + (k - 1) * ldt] != zero;
            }
        }
        //
        //        Determine whether condition numbers are required for the k-th
        //        eigenpair.
        //
        if (somcon) {
            if (pair) {
                if (!select(k) && !select(k + 1)) {
                    goto statement_60;
                }
            } else {
                if (!select(k)) {
                    goto statement_60;
                }
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
            if (!pair) {
                //
                //              Real eigenvalue.
                //
                prod = Rdot(n, vr[(ks - 1) * ldvr], 1, vl[(ks - 1) * ldvl], 1);
                rnrm = Rnrm2(n, vr[(ks - 1) * ldvr], 1);
                lnrm = Rnrm2(n, vl[(ks - 1) * ldvl], 1);
                s[ks - 1] = abs(prod) / (rnrm * lnrm);
            } else {
                //
                //              Complex eigenvalue.
                //
                prod1 = Rdot(n, vr[(ks - 1) * ldvr], 1, vl[(ks - 1) * ldvl], 1);
                prod1 += Rdot(n, vr[((ks + 1) - 1) * ldvr], 1, vl[((ks + 1) - 1) * ldvl], 1);
                prod2 = Rdot(n, vl[(ks - 1) * ldvl], 1, vr[((ks + 1) - 1) * ldvr], 1);
                prod2 = prod2 - Rdot(n, vl[((ks + 1) - 1) * ldvl], 1, vr[(ks - 1) * ldvr], 1);
                rnrm = Rlapy2(Rnrm2(n, vr[(ks - 1) * ldvr], 1), Rnrm2(n, vr[((ks + 1) - 1) * ldvr], 1));
                lnrm = Rlapy2(Rnrm2(n, vl[(ks - 1) * ldvl], 1), Rnrm2(n, vl[((ks + 1) - 1) * ldvl], 1));
                cond = Rlapy2(prod1, prod2) / (rnrm * lnrm);
                s[ks - 1] = cond;
                s[(ks + 1) - 1] = cond;
            }
        }
        //
        if (wantsp) {
            //
            //           Estimate the reciprocal condition number of the k-th
            //           eigenvector.
            //
            //           Copy the matrix T to the array WORK and swap the diagonal
            //           block beginning at T(k,k) to the (1,1) position.
            //
            Rlacpy("Full", n, n, t, ldt, work, ldwork);
            ifst = k;
            ilst = 1;
            Rtrexc("No Q", n, work, ldwork, dummy, 1, ifst, ilst, &work[((n + 1) - 1) * ldwork], ierr);
            //
            if (ierr == 1 || ierr == 2) {
                //
                //              Could not swap because blocks not well separated
                //
                scale = one;
                est = bignum;
            } else {
                //
                //              Reordering successful
                //
                if (work[(2 - 1)] == zero) {
                    //
                    //                 Form C = T22 - lambda*I in WORK(2:N,2:N).
                    //
                    for (i = 2; i <= n; i = i + 1) {
                        work[(i - 1) + (i - 1) * ldwork] = work[(i - 1) + (i - 1) * ldwork] - work[(1 - 1)];
                    }
                    n2 = 1;
                    nn = n - 1;
                } else {
                    //
                    //                 Triangularize the 2 by 2 block by unitary
                    //                 transformation U = [  cs   i*ss ]
                    //                                    [ i*ss   cs  ].
                    //                 such that the (1,1) position of WORK is complex
                    //                 eigenvalue lambda with positive imaginary part. (2,2)
                    //                 position of WORK is the complex eigenvalue lambda
                    //                 with negative imaginary  part.
                    //
                    mu = sqrt(abs(work[(2 - 1) * ldwork])) * sqrt(abs(work[(2 - 1)]));
                    delta = Rlapy2(mu, &work[(2 - 1)]);
                    cs = mu / delta;
                    sn = -work[(2 - 1)] / delta;
                    //
                    //                 Form
                    //
                    //                 C**T = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ]
                    //                                          [   mu                     ]
                    //                                          [         ..               ]
                    //                                          [             ..           ]
                    //                                          [                  mu      ]
                    //                 where C**T is transpose of matrix C,
                    //                 and RWORK is stored starting in the N+1-st column of
                    //                 WORK.
                    //
                    for (j = 3; j <= n; j = j + 1) {
                        work[(2 - 1) + (j - 1) * ldwork] = cs * work[(2 - 1) + (j - 1) * ldwork];
                        work[(j - 1) + (j - 1) * ldwork] = work[(j - 1) + (j - 1) * ldwork] - work[(1 - 1)];
                    }
                    work[(2 - 1) + (2 - 1) * ldwork] = zero;
                    //
                    work[((n + 1) - 1) * ldwork] = two * mu;
                    for (i = 2; i <= n - 1; i = i + 1) {
                        work[(i - 1) + ((n + 1) - 1) * ldwork] = sn * work[((i + 1) - 1) * ldwork];
                    }
                    n2 = 2;
                    nn = 2 * (n - 1);
                }
                //
                //              Estimate norm(inv(C**T))
                //
                est = zero;
                kase = 0;
            statement_50:
                Rlacn2(nn, &work[((n + 2) - 1) * ldwork], &work[((n + 4) - 1) * ldwork], iwork, est, kase, isave);
                if (kase != 0) {
                    if (kase == 1) {
                        if (n2 == 1) {
                            //
                            //                       Real eigenvalue: solve C**T*x = scale*c.
                            //
                            Rlaqtr(true, true, n - 1, &work[(2 - 1) + (2 - 1) * ldwork], ldwork, dummy, dumm, scale, &work[((n + 4) - 1) * ldwork], &work[((n + 6) - 1) * ldwork], ierr);
                        } else {
                            //
                            //                       Complex eigenvalue: solve
                            //                       C**T*(p+iq) = scale*(c+id) in real arithmetic.
                            //
                            Rlaqtr(true, false, n - 1, &work[(2 - 1) + (2 - 1) * ldwork], ldwork, &work[((n + 1) - 1) * ldwork], mu, scale, &work[((n + 4) - 1) * ldwork], &work[((n + 6) - 1) * ldwork], ierr);
                        }
                    } else {
                        if (n2 == 1) {
                            //
                            //                       Real eigenvalue: solve C*x = scale*c.
                            //
                            Rlaqtr(false, true, n - 1, &work[(2 - 1) + (2 - 1) * ldwork], ldwork, dummy, dumm, scale, &work[((n + 4) - 1) * ldwork], &work[((n + 6) - 1) * ldwork], ierr);
                        } else {
                            //
                            //                       Complex eigenvalue: solve
                            //                       C*(p+iq) = scale*(c+id) in real arithmetic.
                            //
                            Rlaqtr(false, false, n - 1, &work[(2 - 1) + (2 - 1) * ldwork], ldwork, &work[((n + 1) - 1) * ldwork], mu, scale, &work[((n + 4) - 1) * ldwork], &work[((n + 6) - 1) * ldwork], ierr);
                            //
                        }
                    }
                    //
                    goto statement_50;
                }
            }
            //
            sep[ks - 1] = scale / max(est, smlnum);
            if (pair) {
                sep[(ks + 1) - 1] = sep[ks - 1];
            }
        }
        //
        if (pair) {
            ks++;
        }
    //
    statement_60:;
    }
    //
    //     End of Rtrsna
    //
}
