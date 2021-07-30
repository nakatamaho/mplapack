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

void Rtgsna(const char *job, const char *howmny, bool *select, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *vl, INTEGER const ldvl, REAL *vr, INTEGER const ldvr, REAL *s, REAL *dif, INTEGER const mm, INTEGER &m, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER &info) {
    bool wantbh = false;
    bool wants = false;
    bool wantdf = false;
    bool somcon = false;
    bool lquery = false;
    bool pair = false;
    INTEGER k = 0;
    const REAL zero = 0.0;
    INTEGER lwmin = 0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    INTEGER ks = 0;
    REAL rnrm = 0.0;
    REAL lnrm = 0.0;
    const REAL one = 1.0;
    REAL tmprr = 0.0;
    REAL tmpri = 0.0;
    REAL tmpii = 0.0;
    REAL tmpir = 0.0;
    REAL uhav = 0.0;
    REAL uhavi = 0.0;
    REAL uhbv = 0.0;
    REAL uhbvi = 0.0;
    REAL cond = 0.0;
    REAL beta = 0.0;
    REAL dummy1[1];
    REAL alphar = 0.0;
    REAL dummy[1];
    REAL alphai = 0.0;
    REAL alprqt = 0.0;
    const REAL two = 2.0e+0;
    REAL c1 = 0.0;
    const REAL four = 4.0e+0;
    REAL c2 = 0.0;
    REAL root1 = 0.0;
    REAL root2 = 0.0;
    INTEGER ifst = 0;
    INTEGER ilst = 0;
    INTEGER ierr = 0;
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER i = 0;
    INTEGER iz = 0;
    const INTEGER difdri = 3;
    REAL scale = 0.0;
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
            pair = false;
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
        } else {
            m = n;
        }
        //
        if (n == 0) {
            lwmin = 1;
        } else if (Mlsame(job, "V") || Mlsame(job, "B")) {
            lwmin = 2 * n * (n + 2) + 16;
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
        Mxerbla("Rtgsna", -info);
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
    ks = 0;
    pair = false;
    //
    for (k = 1; k <= n; k = k + 1) {
        //
        //        Determine whether A(k,k) begins a 1-by-1 or 2-by-2 block.
        //
        if (pair) {
            pair = false;
            goto statement_20;
        } else {
            if (k < n) {
                pair = a[((k + 1) - 1) + (k - 1) * lda] != zero;
            }
        }
        //
        //        Determine whether condition numbers are required for the k-th
        //        eigenpair.
        //
        if (somcon) {
            if (pair) {
                if (!select[k - 1] && !select[(k + 1) - 1]) {
                    goto statement_20;
                }
            } else {
                if (!select[k - 1]) {
                    goto statement_20;
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
            if (pair) {
                //
                //              Complex eigenvalue pair.
                //
                rnrm = Rlapy2(Rnrm2(n, &vr[(ks - 1) * ldvr], 1), Rnrm2(n, &vr[((ks + 1) - 1) * ldvr], 1));
                lnrm = Rlapy2(Rnrm2(n, &vl[(ks - 1) * ldvl], 1), Rnrm2(n, &vl[((ks + 1) - 1) * ldvl], 1));
                Rgemv("N", n, n, one, a, lda, &vr[(ks - 1) * ldvr], 1, zero, work, 1);
                tmprr = Rdot(n, work, 1, &vl[(ks - 1) * ldvl], 1);
                tmpri = Rdot(n, work, 1, &vl[((ks + 1) - 1) * ldvl], 1);
                Rgemv("N", n, n, one, a, lda, &vr[((ks + 1) - 1) * ldvr], 1, zero, work, 1);
                tmpii = Rdot(n, work, 1, &vl[((ks + 1) - 1) * ldvl], 1);
                tmpir = Rdot(n, work, 1, &vl[(ks - 1) * ldvl], 1);
                uhav = tmprr + tmpii;
                uhavi = tmpir - tmpri;
                Rgemv("N", n, n, one, b, ldb, &vr[(ks - 1) * ldvr], 1, zero, work, 1);
                tmprr = Rdot(n, work, 1, &vl[(ks - 1) * ldvl], 1);
                tmpri = Rdot(n, work, 1, &vl[((ks + 1) - 1) * ldvl], 1);
                Rgemv("N", n, n, one, b, ldb, &vr[((ks + 1) - 1) * ldvr], 1, zero, work, 1);
                tmpii = Rdot(n, work, 1, &vl[((ks + 1) - 1) * ldvl], 1);
                tmpir = Rdot(n, work, 1, &vl[(ks - 1) * ldvl], 1);
                uhbv = tmprr + tmpii;
                uhbvi = tmpir - tmpri;
                uhav = Rlapy2(uhav, uhavi);
                uhbv = Rlapy2(uhbv, uhbvi);
                cond = Rlapy2(uhav, uhbv);
                s[ks - 1] = cond / (rnrm * lnrm);
                s[(ks + 1) - 1] = s[ks - 1];
                //
            } else {
                //
                //              Real eigenvalue.
                //
                rnrm = Rnrm2(n, &vr[(ks - 1) * ldvr], 1);
                lnrm = Rnrm2(n, &vl[(ks - 1) * ldvl], 1);
                Rgemv("N", n, n, one, a, lda, &vr[(ks - 1) * ldvr], 1, zero, work, 1);
                uhav = Rdot(n, work, 1, &vl[(ks - 1) * ldvl], 1);
                Rgemv("N", n, n, one, b, ldb, &vr[(ks - 1) * ldvr], 1, zero, work, 1);
                uhbv = Rdot(n, work, 1, &vl[(ks - 1) * ldvl], 1);
                cond = Rlapy2(uhav, uhbv);
                if (cond == zero) {
                    s[ks - 1] = -one;
                } else {
                    s[ks - 1] = cond / (rnrm * lnrm);
                }
            }
        }
        //
        if (wantdf) {
            if (n == 1) {
                dif[ks - 1] = Rlapy2(a[(1 - 1)], b[(1 - 1)]);
                goto statement_20;
            }
            //
            //           Estimate the reciprocal condition number of the k-th
            //           eigenvectors.
            if (pair) {
                //
                //              Copy the  2-by 2 pencil beginning at (A(k,k), B(k, k)).
                //              Compute the eigenvalue(s) at position K.
                //
                work[1 - 1] = a[(k - 1) + (k - 1) * lda];
                work[2 - 1] = a[((k + 1) - 1) + (k - 1) * lda];
                work[3 - 1] = a[(k - 1) + ((k + 1) - 1) * lda];
                work[4 - 1] = a[((k + 1) - 1) + ((k + 1) - 1) * lda];
                work[5 - 1] = b[(k - 1) + (k - 1) * ldb];
                work[6 - 1] = b[((k + 1) - 1) + (k - 1) * ldb];
                work[7 - 1] = b[(k - 1) + ((k + 1) - 1) * ldb];
                work[8 - 1] = b[((k + 1) - 1) + ((k + 1) - 1) * ldb];
                Rlag2(work, 2, &work[5 - 1], 2, smlnum * eps, beta, dummy1[1 - 1], alphar, dummy[1 - 1], alphai);
                alprqt = one;
                c1 = two * (alphar * alphar + alphai * alphai + beta * beta);
                c2 = four * beta * beta * alphai * alphai;
                root1 = c1 + sqrt(c1 * c1 - 4.0 * c2);
                root2 = c2 / root1;
                root1 = root1 / two;
                cond = min(sqrt(root1), sqrt(root2));
            }
            //
            //           Copy the matrix (A, B) to the array WORK and swap the
            //           diagonal block beginning at A(k,k) to the (1,1) position.
            //
            Rlacpy("Full", n, n, a, lda, work, n);
            Rlacpy("Full", n, n, b, ldb, &work[(n * n + 1) - 1], n);
            ifst = k;
            ilst = 1;
            //
            Rtgexc(false, false, n, work, n, &work[(n * n + 1) - 1], n, dummy, 1, dummy1, 1, ifst, ilst, &work[(n * n * 2 + 1) - 1], lwork - 2 * n * n, ierr);
            //
            if (ierr > 0) {
                //
                //              Ill-conditioned problem - swap rejected.
                //
                dif[ks - 1] = zero;
            } else {
                //
                //              Reordering successful, solve generalized Sylvester
                //              equation for R and L,
                //                         A22 * R - L * A11 = A12
                //                         B22 * R - L * B11 = B12,
                //              and compute estimate of Difl((A11,B11), (A22, B22)).
                //
                n1 = 1;
                if (work[2 - 1] != zero) {
                    n1 = 2;
                }
                n2 = n - n1;
                if (n2 == 0) {
                    dif[ks - 1] = cond;
                } else {
                    i = n * n + 1;
                    iz = 2 * n * n + 1;
                    Rtgsyl("N", difdri, n2, n1, &work[(n * n1 + n1 + 1) - 1], n, work, n, &work[(n1 + 1) - 1], n, &work[(n * n1 + n1 + i) - 1], n, &work[i - 1], n, &work[(n1 + i) - 1], n, scale, dif[ks - 1], &work[(iz + 1) - 1], lwork - 2 * n * n, iwork, ierr);
                    //
                    if (pair) {
                        dif[ks - 1] = min(REAL(max(one, alprqt) * dif[ks - 1]), cond);
                    }
                }
            }
            if (pair) {
                dif[(ks + 1) - 1] = dif[ks - 1];
            }
        }
        if (pair) {
            ks++;
        }
    //
    statement_20:;
    }
    work[1 - 1] = lwmin;
    //
    //     End of Rtgsna
    //
}
