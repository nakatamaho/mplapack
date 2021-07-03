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

void Rtgsen(INTEGER const ijob, bool const wantq, bool const wantz, bool *select, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *alphar, REAL *alphai, REAL *beta, REAL *q, INTEGER const ldq, REAL *z, INTEGER const ldz, INTEGER &m, REAL &pl, REAL &pr, REAL *dif, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
    bool lquery = false;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    INTEGER ierr = 0;
    bool wantp = false;
    bool wantd1 = false;
    bool wantd2 = false;
    bool wantd = false;
    bool pair = false;
    INTEGER k = 0;
    const REAL zero = 0.0;
    INTEGER lwmin = 0;
    INTEGER liwmin = 0;
    const REAL one = 1.0;
    REAL Rscale = 0.0;
    REAL dsum = 0.0;
    INTEGER i = 0;
    INTEGER ks = 0;
    bool swap = false;
    INTEGER kk = 0;
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER ijb = 0;
    REAL rRscal = 0.0;
    const INTEGER idifjb = 3;
    INTEGER kase = 0;
    INTEGER mn2 = 0;
    INTEGER isave[3];
    //
    //     Decode and test the input parameters
    //
    info = 0;
    lquery = (lwork == -1 || liwork == -1);
    //
    if (ijob < 0 || ijob > 5) {
        info = -1;
    } else if (n < 0) {
        info = -5;
    } else if (lda < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -9;
    } else if (ldq < 1 || (wantq && ldq < n)) {
        info = -14;
    } else if (ldz < 1 || (wantz && ldz < n)) {
        info = -16;
    }
    //
    if (info != 0) {
        Mxerbla("Rtgsen", -info);
        return;
    }
    //
    //     Get machine constants
    //
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    ierr = 0;
    //
    wantp = ijob == 1 || ijob >= 4;
    wantd1 = ijob == 2 || ijob == 4;
    wantd2 = ijob == 3 || ijob == 5;
    wantd = wantd1 || wantd2;
    //
    //     subspaces.
    //
    m = 0;
    pair = false;
    if (!lquery || ijob != 0) {
        for (k = 1; k <= n; k = k + 1) {
            if (pair) {
                pair = false;
            } else {
                if (k < n) {
                    if (a[((k + 1) - 1) + (k - 1) * lda] == zero) {
                        if (select[k - 1]) {
                            m++;
                        }
                    } else {
                        pair = true;
                        if (select[k - 1] || select[(k + 1) - 1]) {
                            m += 2;
                        }
                    }
                } else {
                    if (select[n - 1]) {
                        m++;
                    }
                }
            }
        }
    }
    //
    if (ijob == 1 || ijob == 2 || ijob == 4) {
        lwmin = max({(INTEGER)1, 4 * n + 16, 2 * m * (n - m)});
        liwmin = max((INTEGER)1, n + 6);
    } else if (ijob == 3 || ijob == 5) {
        lwmin = max({(INTEGER)1, 4 * n + 16, 4 * m * (n - m)});
        liwmin = max({(INTEGER)1, 2 * m * (n - m), n + 6});
    } else {
        lwmin = max((INTEGER)1, 4 * n + 16);
        liwmin = 1;
    }
    //
    work[1 - 1] = lwmin;
    iwork[1 - 1] = liwmin;
    //
    if (lwork < lwmin && !lquery) {
        info = -22;
    } else if (liwork < liwmin && !lquery) {
        info = -24;
    }
    //
    if (info != 0) {
        Mxerbla("Rtgsen", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible.
    //
    if (m == n || m == 0) {
        if (wantp) {
            pl = one;
            pr = one;
        }
        if (wantd) {
            Rscale = zero;
            dsum = one;
            for (i = 1; i <= n; i = i + 1) {
                Rlassq(n, &a[(i - 1) * lda], 1, Rscale, dsum);
                Rlassq(n, &b[(i - 1) * ldb], 1, Rscale, dsum);
            }
            dif[1 - 1] = Rscale * sqrt(dsum);
            dif[2 - 1] = dif[1 - 1];
        }
        goto statement_60;
    }
    //
    //     Collect the selected blocks at the top-left corner of (A, B).
    //
    ks = 0;
    pair = false;
    for (k = 1; k <= n; k = k + 1) {
        if (pair) {
            pair = false;
        } else {
            //
            swap = select[k - 1];
            if (k < n) {
                if (a[((k + 1) - 1) + (k - 1) * lda] != zero) {
                    pair = true;
                    swap = swap || select[(k + 1) - 1];
                }
            }
            //
            if (swap) {
                ks++;
                //
                //              Swap the K-th block to position KS.
                //              Perform the reordering of diagonal blocks in (A, B)
                //              by orthogonal transformation matrices and update
                //              Q and Z accordingly (if requested):
                //
                kk = k;
                if (k != ks) {
                    Rtgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, kk, ks, work, lwork, ierr);
                }
                //
                if (ierr > 0) {
                    //
                    //                 Swap is rejected: exit.
                    //
                    info = 1;
                    if (wantp) {
                        pl = zero;
                        pr = zero;
                    }
                    if (wantd) {
                        dif[1 - 1] = zero;
                        dif[2 - 1] = zero;
                    }
                    goto statement_60;
                }
                //
                if (pair) {
                    ks++;
                }
            }
        }
    }
    if (wantp) {
        //
        //        Solve generalized Sylvester equation for R and L
        //        and compute PL and PR.
        //
        n1 = m;
        n2 = n - m;
        i = n1 + 1;
        ijb = 0;
        Rlacpy("Full", n1, n2, &a[(i - 1) * lda], lda, work, n1);
        Rlacpy("Full", n1, n2, &b[(i - 1) * ldb], ldb, &work[(n1 * n2 + 1) - 1], n1);
        Rtgsyl("N", ijb, n1, n2, a, lda, &a[(i - 1) + (i - 1) * lda], lda, work, n1, b, ldb, &b[(i - 1) + (i - 1) * ldb], ldb, &work[(n1 * n2 + 1) - 1], n1, Rscale, dif[1 - 1], &work[(n1 * n2 * 2 + 1) - 1], lwork - 2 * n1 * n2, iwork, ierr);
        //
        //        Estimate the reciprocal of norms of "projections" onto left
        //        and right eigenspaces.
        //
        rRscal = zero;
        dsum = one;
        Rlassq(n1 * n2, work, 1, rRscal, dsum);
        pl = rRscal * sqrt(dsum);
        if (pl == zero) {
            pl = one;
        } else {
            pl = Rscale / (sqrt(Rscale * Rscale / pl + pl) * sqrt(pl));
        }
        rRscal = zero;
        dsum = one;
        Rlassq(n1 * n2, &work[(n1 * n2 + 1) - 1], 1, rRscal, dsum);
        pr = rRscal * sqrt(dsum);
        if (pr == zero) {
            pr = one;
        } else {
            pr = Rscale / (sqrt(Rscale * Rscale / pr + pr) * sqrt(pr));
        }
    }
    //
    if (wantd) {
        //
        //        Compute estimates of Difu and Difl.
        //
        if (wantd1) {
            n1 = m;
            n2 = n - m;
            i = n1 + 1;
            ijb = idifjb;
            //
            //           Frobenius norm-based Difu-estimate.
            //
            Rtgsyl("N", ijb, n1, n2, a, lda, &a[(i - 1) + (i - 1) * lda], lda, work, n1, b, ldb, &b[(i - 1) + (i - 1) * ldb], ldb, &work[(n1 * n2 + 1) - 1], n1, Rscale, dif[1 - 1], &work[(2 * n1 * n2 + 1) - 1], lwork - 2 * n1 * n2, iwork, ierr);
            //
            //           Frobenius norm-based Difl-estimate.
            //
            Rtgsyl("N", ijb, n2, n1, &a[(i - 1) + (i - 1) * lda], lda, a, lda, work, n2, &b[(i - 1) + (i - 1) * ldb], ldb, b, ldb, &work[(n1 * n2 + 1) - 1], n2, Rscale, dif[2 - 1], &work[(2 * n1 * n2 + 1) - 1], lwork - 2 * n1 * n2, iwork, ierr);
        } else {
            //
            //           Compute 1-norm-based estimates of Difu and Difl using
            //           reversed communication with Rlacn2. In each step a
            //           generalized Sylvester equation or a transposed variant
            //           is solved.
            //
            kase = 0;
            n1 = m;
            n2 = n - m;
            i = n1 + 1;
            ijb = 0;
            mn2 = 2 * n1 * n2;
        //
        //           1-norm-based estimate of Difu.
        //
        statement_40:
            Rlacn2(mn2, &work[(mn2 + 1) - 1], work, iwork, dif[1 - 1], kase, isave);
            if (kase != 0) {
                if (kase == 1) {
                    //
                    //                 Solve generalized Sylvester equation.
                    //
                    Rtgsyl("N", ijb, n1, n2, a, lda, &a[(i - 1) + (i - 1) * lda], lda, work, n1, b, ldb, &b[(i - 1) + (i - 1) * ldb], ldb, &work[(n1 * n2 + 1) - 1], n1, Rscale, dif[1 - 1], &work[(2 * n1 * n2 + 1) - 1], lwork - 2 * n1 * n2, iwork, ierr);
                } else {
                    //
                    //                 Solve the transposed variant.
                    //
                    Rtgsyl("T", ijb, n1, n2, a, lda, &a[(i - 1) + (i - 1) * lda], lda, work, n1, b, ldb, &b[(i - 1) + (i - 1) * ldb], ldb, &work[(n1 * n2 + 1) - 1], n1, Rscale, dif[1 - 1], &work[(2 * n1 * n2 + 1) - 1], lwork - 2 * n1 * n2, iwork, ierr);
                }
                goto statement_40;
            }
            dif[1 - 1] = Rscale / dif[1 - 1];
        //
        //           1-norm-based estimate of Difl.
        //
        statement_50:
            Rlacn2(mn2, &work[(mn2 + 1) - 1], work, iwork, dif[2 - 1], kase, isave);
            if (kase != 0) {
                if (kase == 1) {
                    //
                    //                 Solve generalized Sylvester equation.
                    //
                    Rtgsyl("N", ijb, n2, n1, &a[(i - 1) + (i - 1) * lda], lda, a, lda, work, n2, &b[(i - 1) + (i - 1) * ldb], ldb, b, ldb, &work[(n1 * n2 + 1) - 1], n2, Rscale, dif[2 - 1], &work[(2 * n1 * n2 + 1) - 1], lwork - 2 * n1 * n2, iwork, ierr);
                } else {
                    //
                    //                 Solve the transposed variant.
                    //
                    Rtgsyl("T", ijb, n2, n1, &a[(i - 1) + (i - 1) * lda], lda, a, lda, work, n2, &b[(i - 1) + (i - 1) * ldb], ldb, b, ldb, &work[(n1 * n2 + 1) - 1], n2, Rscale, dif[2 - 1], &work[(2 * n1 * n2 + 1) - 1], lwork - 2 * n1 * n2, iwork, ierr);
                }
                goto statement_50;
            }
            dif[2 - 1] = Rscale / dif[2 - 1];
            //
        }
    }
//
statement_60:
    //
    //     Compute generalized eigenvalues of reordered pair (A, B) and
    //     normalize the generalized Schur form.
    //
    pair = false;
    for (k = 1; k <= n; k = k + 1) {
        if (pair) {
            pair = false;
        } else {
            //
            if (k < n) {
                if (a[((k + 1) - 1) + (k - 1) * lda] != zero) {
                    pair = true;
                }
            }
            //
            if (pair) {
                //
                //             Compute the eigenvalue(s) at position K.
                //
                work[1 - 1] = a[(k - 1) + (k - 1) * lda];
                work[2 - 1] = a[((k + 1) - 1) + (k - 1) * lda];
                work[3 - 1] = a[(k - 1) + ((k + 1) - 1) * lda];
                work[4 - 1] = a[((k + 1) - 1) + ((k + 1) - 1) * lda];
                work[5 - 1] = b[(k - 1) + (k - 1) * ldb];
                work[6 - 1] = b[((k + 1) - 1) + (k - 1) * ldb];
                work[7 - 1] = b[(k - 1) + ((k + 1) - 1) * ldb];
                work[8 - 1] = b[((k + 1) - 1) + ((k + 1) - 1) * ldb];
                Rlag2(work, 2, &work[5 - 1], 2, smlnum * eps, beta[k - 1], beta[(k + 1) - 1], alphar[k - 1], alphar[(k + 1) - 1], alphai[k - 1]);
                alphai[(k + 1) - 1] = -alphai[k - 1];
                //
            } else {
                //
                if (sign(one, b[(k - 1) + (k - 1) * ldb]) < zero) {
                    //
                    //                 If B(K,K) is negative, make it positive
                    //
                    for (i = 1; i <= n; i = i + 1) {
                        a[(k - 1) + (i - 1) * lda] = -a[(k - 1) + (i - 1) * lda];
                        b[(k - 1) + (i - 1) * ldb] = -b[(k - 1) + (i - 1) * ldb];
                        if (wantq) {
                            q[(i - 1) + (k - 1) * ldq] = -q[(i - 1) + (k - 1) * ldq];
                        }
                    }
                }
                //
                alphar[k - 1] = a[(k - 1) + (k - 1) * lda];
                alphai[k - 1] = zero;
                beta[k - 1] = b[(k - 1) + (k - 1) * ldb];
                //
            }
        }
    }
    //
    work[1 - 1] = lwmin;
    iwork[1 - 1] = liwmin;
    //
    //     End of Rtgsen
    //
}
