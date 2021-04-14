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

void Rtrsen(const char *job, const char *compq, bool *select, INTEGER const n, REAL *t, INTEGER const ldt, REAL *q, INTEGER const ldq, REAL *wr, REAL *wi, INTEGER &m, REAL &s, REAL &sep, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
    bool wantbh = false;
    bool wants = false;
    bool wantsp = false;
    bool wantq = false;
    bool lquery = false;
    bool pair = false;
    INTEGER k = 0;
    const REAL zero = 0.0;
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER nn = 0;
    INTEGER lwmin = 0;
    INTEGER liwmin = 0;
    const REAL one = 1.0;
    INTEGER ks = 0;
    bool swap = false;
    INTEGER ierr = 0;
    INTEGER kk = 0;
    REAL scale = 0.0;
    REAL rnorm = 0.0;
    REAL est = 0.0;
    INTEGER kase = 0;
    arr_1d<3, int> isave(fill0);
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
    wantq = Mlsame(compq, "V");
    //
    info = 0;
    lquery = (lwork == -1);
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
    } else {
        //
        //        and test LWORK and LIWORK.
        //
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
        //
        n1 = m;
        n2 = n - m;
        nn = n1 * n2;
        //
        if (wantsp) {
            lwmin = max((INTEGER)1, 2 * nn);
            liwmin = max((INTEGER)1, nn);
        } else if (Mlsame(job, "N")) {
            lwmin = max((INTEGER)1, n);
            liwmin = 1;
        } else if (Mlsame(job, "E")) {
            lwmin = max((INTEGER)1, nn);
            liwmin = 1;
        }
        //
        if (lwork < lwmin && !lquery) {
            info = -15;
        } else if (liwork < liwmin && !lquery) {
            info = -17;
        }
    }
    //
    if (info == 0) {
        work[1 - 1] = lwmin;
        iwork[1 - 1] = liwmin;
    }
    //
    if (info != 0) {
        Mxerbla("Rtrsen", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible.
    //
    if (m == n || m == 0) {
        if (wants) {
            s = one;
        }
        if (wantsp) {
            sep = Rlange("1", n, n, t, ldt, work);
        }
        goto statement_40;
    }
    //
    //     Collect the selected blocks at the top-left corner of T.
    //
    ks = 0;
    pair = false;
    for (k = 1; k <= n; k = k + 1) {
        if (pair) {
            pair = false;
        } else {
            swap = select(k);
            if (k < n) {
                if (t[((k + 1) - 1) + (k - 1) * ldt] != zero) {
                    pair = true;
                    swap = swap || select(k + 1);
                }
            }
            if (swap) {
                ks++;
                //
                //              Swap the K-th block to position KS.
                //
                ierr = 0;
                kk = k;
                if (k != ks) {
                    Rtrexc(compq, n, t, ldt, q, ldq, kk, ks, work, ierr);
                }
                if (ierr == 1 || ierr == 2) {
                    //
                    //                 Blocks too close to swap: exit.
                    //
                    info = 1;
                    if (wants) {
                        s = zero;
                    }
                    if (wantsp) {
                        sep = zero;
                    }
                    goto statement_40;
                }
                if (pair) {
                    ks++;
                }
            }
        }
    }
    //
    if (wants) {
        //
        //        Solve Sylvester equation for R:
        //
        //           T11*R - R*T22 = scale*T12
        //
        Rlacpy("F", n1, n2, &t[((n1 + 1) - 1) * ldt], ldt, work, n1);
        Rtrsyl("N", "N", -1, n1, n2, t, ldt, &t[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], ldt, work, n1, scale, ierr);
        //
        //        Estimate the reciprocal of the condition number of the cluster
        //        of eigenvalues.
        //
        rnorm = Rlange("F", n1, n2, work, n1, work);
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
        Rlacn2(nn, &work[(nn + 1) - 1], work, iwork, est, kase, isave);
        if (kase != 0) {
            if (kase == 1) {
                //
                //              Solve  T11*R - R*T22 = scale*X.
                //
                Rtrsyl("N", "N", -1, n1, n2, t, ldt, &t[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], ldt, work, n1, scale, ierr);
            } else {
                //
                //              Solve T11**T*R - R*T22**T = scale*X.
                //
                Rtrsyl("T", "T", -1, n1, n2, t, ldt, &t[((n1 + 1) - 1) + ((n1 + 1) - 1) * ldt], ldt, work, n1, scale, ierr);
            }
            goto statement_30;
        }
        //
        sep = scale / est;
    }
//
statement_40:
    //
    //     Store the output eigenvalues in WR and WI.
    //
    for (k = 1; k <= n; k = k + 1) {
        wr[k - 1] = t[(k - 1) + (k - 1) * ldt];
        wi[k - 1] = zero;
    }
    for (k = 1; k <= n - 1; k = k + 1) {
        if (t[((k + 1) - 1) + (k - 1) * ldt] != zero) {
            wi[k - 1] = sqrt(abs(t[(k - 1) + ((k + 1) - 1) * ldt])) * sqrt(abs(t[((k + 1) - 1) + (k - 1) * ldt]));
            wi[(k + 1) - 1] = -wi[k - 1];
        }
    }
    //
    work[1 - 1] = lwmin;
    iwork[1 - 1] = liwmin;
    //
    //     End of Rtrsen
    //
}
