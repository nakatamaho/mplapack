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

void Rgesvdq(const char *joba, const char *jobp, const char *jobr, const char *jobu, const char *jobv, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *s, REAL *u, INTEGER const ldu, REAL *v, INTEGER const ldv, INTEGER &numrank, INTEGER *iwork, INTEGER const liwork, REAL *work, INTEGER const lwork, REAL *rwork, INTEGER const lrwork, INTEGER &info) {
    bool wntus = false;
    bool wntur = false;
    bool wntua = false;
    bool wntuf = false;
    bool lsvc0 = false;
    bool lsvec = false;
    bool dntwu = false;
    bool wntvr = false;
    bool wntva = false;
    bool rsvec = false;
    bool dntwv = false;
    bool accla = false;
    bool acclm = false;
    bool conda = false;
    bool acclh = false;
    bool rowprm = false;
    bool rtrans = false;
    INTEGER iminwrk = 0;
    INTEGER rminwrk = 0;
    bool lquery = false;
    INTEGER lwqp3 = 0;
    INTEGER lworq = 0;
    INTEGER lwcon = 0;
    INTEGER lwsvd = 0;
    REAL rdummy[1];
    INTEGER ierr = 0;
    INTEGER lwrk_Rgeqp3 = 0;
    INTEGER lwrk_Rormqr = 0;
    INTEGER minwrk = 0;
    INTEGER optwrk = 0;
    INTEGER lwrk_Rgesvd = 0;
    INTEGER lwqrf = 0;
    INTEGER lwsvd2 = 0;
    INTEGER lworq2 = 0;
    INTEGER minwrk2 = 0;
    INTEGER lwlqf = 0;
    INTEGER lworlq = 0;
    INTEGER lwrk_Rgeqrf = 0;
    INTEGER lwrk_Rgesvd2 = 0;
    INTEGER lwrk_Rormqr2 = 0;
    INTEGER optwrk2 = 0;
    INTEGER lwrk_Rgelqf = 0;
    INTEGER lwrk_Rormlq = 0;
    REAL big = 0.0;
    bool ascaled = false;
    INTEGER iwoff = 0;
    INTEGER p = 0;
    const REAL zero = 0.0;
    INTEGER q = 0;
    REAL rtmp = 0.0;
    const REAL one = 1.0;
    REAL epsln = 0.0;
    REAL sfmin = 0.0;
    INTEGER nr = 0;
    REAL sconda = 0.0;
    INTEGER n1 = 0;
    INTEGER optratio = 0;
    //
    //     Test the input arguments
    //
    wntus = Mlsame(jobu, "S") || Mlsame(jobu, "U");
    wntur = Mlsame(jobu, "R");
    wntua = Mlsame(jobu, "A");
    wntuf = Mlsame(jobu, "F");
    lsvc0 = wntus || wntur || wntua;
    lsvec = lsvc0 || wntuf;
    dntwu = Mlsame(jobu, "N");
    //
    wntvr = Mlsame(jobv, "R");
    wntva = Mlsame(jobv, "A") || Mlsame(jobv, "V");
    rsvec = wntvr || wntva;
    dntwv = Mlsame(jobv, "N");
    //
    accla = Mlsame(joba, "A");
    acclm = Mlsame(joba, "M");
    conda = Mlsame(joba, "E");
    acclh = Mlsame(joba, "H") || conda;
    //
    rowprm = Mlsame(jobp, "P");
    rtrans = Mlsame(jobr, "T");
    //
    if (rowprm) {
        if (conda) {
            iminwrk = max((INTEGER)1, n + m - 1 + n);
        } else {
            iminwrk = max((INTEGER)1, n + m - 1);
        }
        rminwrk = max((INTEGER)2, m);
    } else {
        if (conda) {
            iminwrk = max((INTEGER)1, n + n);
        } else {
            iminwrk = max((INTEGER)1, n);
        }
        rminwrk = 2;
    }
    lquery = (liwork == -1 || lwork == -1 || lrwork == -1);
    info = 0;
    if (!(accla || acclm || acclh)) {
        info = -1;
    } else if (!(rowprm || Mlsame(jobp, "N"))) {
        info = -2;
    } else if (!(rtrans || Mlsame(jobr, "N"))) {
        info = -3;
    } else if (!(lsvec || dntwu)) {
        info = -4;
    } else if (wntur && wntva) {
        info = -5;
    } else if (!(rsvec || dntwv)) {
        info = -5;
    } else if (m < 0) {
        info = -6;
    } else if ((n < 0) || (n > m)) {
        info = -7;
    } else if (lda < max((INTEGER)1, m)) {
        info = -9;
    } else if (ldu < 1 || (lsvc0 && ldu < m) || (wntuf && ldu < n)) {
        info = -12;
    } else if (ldv < 1 || (rsvec && ldv < n) || (conda && ldv < n)) {
        info = -14;
    } else if (liwork < iminwrk && !lquery) {
        info = -17;
    }
    //
    if (info == 0) {
        //        .. compute the minimal and the optimal workspace lengths
        //        [[The expressions for computing the minimal and the optimal
        //        values of LWORK are written with a lot of redundancy and
        //        can be simplified. However, this detailed form is easier for
        //        maintenance and modifications of the code.]]
        //
        //        .. minimal workspace length for Rgeqp3 of an M x N matrix
        lwqp3 = 3 * n + 1;
        //        .. minimal workspace length for Rormqr to build left singular vectors
        if (wntus || wntur) {
            lworq = max(n, (INTEGER)1);
        } else if (wntua) {
            lworq = max(m, (INTEGER)1);
        }
        //        .. minimal workspace length for Rpocon of an N x N matrix
        lwcon = 3 * n;
        //        .. Rgesvd of an N x N matrix
        lwsvd = max(5 * n, (INTEGER)1);
        if (lquery) {
            Rgeqp3(m, n, a, lda, iwork, rdummy, rdummy, -1, ierr);
            lwrk_Rgeqp3 = castINTEGER(rdummy[1 - 1]);
            if (wntus || wntur) {
                Rormqr("L", "N", m, n, n, a, lda, rdummy, u, ldu, rdummy, -1, ierr);
                lwrk_Rormqr = castINTEGER(rdummy[1 - 1]);
            } else if (wntua) {
                Rormqr("L", "N", m, m, n, a, lda, rdummy, u, ldu, rdummy, -1, ierr);
                lwrk_Rormqr = castINTEGER(rdummy[1 - 1]);
            } else {
                lwrk_Rormqr = 0;
            }
        }
        minwrk = 2;
        optwrk = 2;
        if (!(lsvec || rsvec)) {
            //            .. minimal and optimal sizes of the workspace if
            //            only the singular values are requested
            if (conda) {
                minwrk = max({n + lwqp3, lwcon, lwsvd});
            } else {
                minwrk = max(n + lwqp3, lwsvd);
            }
            if (lquery) {
                Rgesvd("N", "N", n, n, a, lda, s, u, ldu, v, ldv, rdummy, -1, ierr);
                lwrk_Rgesvd = castINTEGER(rdummy[1 - 1]);
                if (conda) {
                    optwrk = max({n + lwrk_Rgeqp3, n + lwcon, lwrk_Rgesvd});
                } else {
                    optwrk = max(n + lwrk_Rgeqp3, lwrk_Rgesvd);
                }
            }
        } else if (lsvec && (!rsvec)) {
            //            .. minimal and optimal sizes of the workspace if the
            //            singular values and the left singular vectors are requested
            if (conda) {
                minwrk = n + max({lwqp3, lwcon, lwsvd, lworq});
            } else {
                minwrk = n + max({lwqp3, lwsvd, lworq});
            }
            if (lquery) {
                if (rtrans) {
                    Rgesvd("N", "O", n, n, a, lda, s, u, ldu, v, ldv, rdummy, -1, ierr);
                } else {
                    Rgesvd("O", "N", n, n, a, lda, s, u, ldu, v, ldv, rdummy, -1, ierr);
                }
                lwrk_Rgesvd = castINTEGER(rdummy[1 - 1]);
                if (conda) {
                    optwrk = n + max({lwrk_Rgeqp3, lwcon, lwrk_Rgesvd, lwrk_Rormqr});
                } else {
                    optwrk = n + max({lwrk_Rgeqp3, lwrk_Rgesvd, lwrk_Rormqr});
                }
            }
        } else if (rsvec && (!lsvec)) {
            //            .. minimal and optimal sizes of the workspace if the
            //            singular values and the right singular vectors are requested
            if (conda) {
                minwrk = n + max({lwqp3, lwcon, lwsvd});
            } else {
                minwrk = n + max(lwqp3, lwsvd);
            }
            if (lquery) {
                if (rtrans) {
                    Rgesvd("O", "N", n, n, a, lda, s, u, ldu, v, ldv, rdummy, -1, ierr);
                } else {
                    Rgesvd("N", "O", n, n, a, lda, s, u, ldu, v, ldv, rdummy, -1, ierr);
                }
                lwrk_Rgesvd = castINTEGER(rdummy[1 - 1]);
                if (conda) {
                    optwrk = n + max({lwrk_Rgeqp3, lwcon, lwrk_Rgesvd});
                } else {
                    optwrk = n + max(lwrk_Rgeqp3, lwrk_Rgesvd);
                }
            }
        } else {
            //            .. minimal and optimal sizes of the workspace if the
            //            full SVD is requested
            if (rtrans) {
                minwrk = max({lwqp3, lwsvd, lworq});
                if (conda) {
                    minwrk = max(minwrk, lwcon);
                }
                minwrk += n;
                if (wntva) {
                    //                   .. minimal workspace length for N x N/2 Rgeqrf
                    lwqrf = max(n / 2, (INTEGER)1);
                    //                   .. minimal workspace length for N/2 x N/2 Rgesvd
                    lwsvd2 = max(5 * (n / 2), (INTEGER)1);
                    lworq2 = max(n, (INTEGER)1);
                    minwrk2 = max({lwqp3, n / 2 + lwqrf, n / 2 + lwsvd2, n / 2 + lworq2, lworq});
                    if (conda) {
                        minwrk2 = max(minwrk2, lwcon);
                    }
                    minwrk2 += n;
                    minwrk = max(minwrk, minwrk2);
                }
            } else {
                minwrk = max({lwqp3, lwsvd, lworq});
                if (conda) {
                    minwrk = max(minwrk, lwcon);
                }
                minwrk += n;
                if (wntva) {
                    //                   .. minimal workspace length for N/2 x N Rgelqf
                    lwlqf = max(n / 2, (INTEGER)1);
                    lwsvd2 = max(5 * (n / 2), (INTEGER)1);
                    lworlq = max(n, (INTEGER)1);
                    minwrk2 = max({lwqp3, n / 2 + lwlqf, n / 2 + lwsvd2, n / 2 + lworlq, lworq});
                    if (conda) {
                        minwrk2 = max(minwrk2, lwcon);
                    }
                    minwrk2 += n;
                    minwrk = max(minwrk, minwrk2);
                }
            }
            if (lquery) {
                if (rtrans) {
                    Rgesvd("O", "A", n, n, a, lda, s, u, ldu, v, ldv, rdummy, -1, ierr);
                    lwrk_Rgesvd = castINTEGER(rdummy[1 - 1]);
                    optwrk = max({lwrk_Rgeqp3, lwrk_Rgesvd, lwrk_Rormqr});
                    if (conda) {
                        optwrk = max(optwrk, lwcon);
                    }
                    optwrk += n;
                    if (wntva) {
                        Rgeqrf(n, n / 2, u, ldu, rdummy, rdummy, -1, ierr);
                        lwrk_Rgeqrf = castINTEGER(rdummy[1 - 1]);
                        Rgesvd("S", "O", n / 2, n / 2, v, ldv, s, u, ldu, v, ldv, rdummy, -1, ierr);
                        lwrk_Rgesvd2 = castINTEGER(rdummy[1 - 1]);
                        Rormqr("R", "C", n, n, n / 2, u, ldu, rdummy, v, ldv, rdummy, -1, ierr);
                        lwrk_Rormqr2 = castINTEGER(rdummy[1 - 1]);
                        optwrk2 = max({lwrk_Rgeqp3, n / 2 + lwrk_Rgeqrf, n / 2 + lwrk_Rgesvd2, n / 2 + lwrk_Rormqr2});
                        if (conda) {
                            optwrk2 = max(optwrk2, lwcon);
                        }
                        optwrk2 += n;
                        optwrk = max(optwrk, optwrk2);
                    }
                } else {
                    Rgesvd("S", "O", n, n, a, lda, s, u, ldu, v, ldv, rdummy, -1, ierr);
                    lwrk_Rgesvd = castINTEGER(rdummy[1 - 1]);
                    optwrk = max({lwrk_Rgeqp3, lwrk_Rgesvd, lwrk_Rormqr});
                    if (conda) {
                        optwrk = max(optwrk, lwcon);
                    }
                    optwrk += n;
                    if (wntva) {
                        Rgelqf(n / 2, n, u, ldu, rdummy, rdummy, -1, ierr);
                        lwrk_Rgelqf = castINTEGER(rdummy[1 - 1]);
                        Rgesvd("S", "O", n / 2, n / 2, v, ldv, s, u, ldu, v, ldv, rdummy, -1, ierr);
                        lwrk_Rgesvd2 = castINTEGER(rdummy[1 - 1]);
                        Rormlq("R", "N", n, n, n / 2, u, ldu, rdummy, v, ldv, rdummy, -1, ierr);
                        lwrk_Rormlq = castINTEGER(rdummy[1 - 1]);
                        optwrk2 = max({lwrk_Rgeqp3, n / 2 + lwrk_Rgelqf, n / 2 + lwrk_Rgesvd2, n / 2 + lwrk_Rormlq});
                        if (conda) {
                            optwrk2 = max(optwrk2, lwcon);
                        }
                        optwrk2 += n;
                        optwrk = max(optwrk, optwrk2);
                    }
                }
            }
        }
        //
        minwrk = max((INTEGER)2, minwrk);
        optwrk = max((INTEGER)2, optwrk);
        if (lwork < minwrk && (!lquery)) {
            info = -19;
        }
        //
    }
    //
    if (info == 0 && lrwork < rminwrk && !lquery) {
        info = -21;
    }
    if (info != 0) {
        Mxerbla("Rgesvdq", -info);
        return;
    } else if (lquery) {
        //
        //     Return optimal workspace
        //
        iwork[1 - 1] = iminwrk;
        work[1 - 1] = optwrk;
        work[2 - 1] = minwrk;
        rwork[1 - 1] = rminwrk;
        return;
    }
    //
    //     Quick return if the matrix is void.
    //
    if ((m == 0) || (n == 0)) {
        //     .. all output is void.
        return;
    }
    //
    big = Rlamch("O");
    ascaled = false;
    iwoff = 1;
    if (rowprm) {
        iwoff = m;
        //           .. reordering the rows in decreasing sequence in the
        //           ell-infinity norm - this enhances numerical robustness in
        //           the case of differently scaled rows.
        for (p = 1; p <= m; p = p + 1) {
            //               RWORK(p) = ABS( A(p,ICAMAX(N,A(p,1),LDA)) )
            //               [[Rlange will return NaN if an entry of the p-th row is Nan]]
            rwork[p - 1] = Rlange("M", 1, n, &a[(p - 1)], lda, rdummy);
            //               .. check for NaN's and Inf's
            if ((rwork[p - 1] != rwork[p - 1]) || ((rwork[p - 1] * zero) != zero)) {
                info = -8;
                Mxerbla("Rgesvdq", -info);
                return;
            }
        }
        for (p = 1; p <= m - 1; p = p + 1) {
            q = iRamax(m - p + 1, &rwork[p - 1], 1) + p - 1;
            iwork[(n + p) - 1] = q;
            if (p != q) {
                rtmp = rwork[p - 1];
                rwork[p - 1] = rwork[q - 1];
                rwork[q - 1] = rtmp;
            }
        }
        //
        if (rwork[1 - 1] == zero) {
            //              Quick return: A is the M x N zero matrix.
            numrank = 0;
            Rlaset("G", n, 1, zero, zero, s, n);
            if (wntus) {
                Rlaset("G", m, n, zero, one, u, ldu);
            }
            if (wntua) {
                Rlaset("G", m, m, zero, one, u, ldu);
            }
            if (wntva) {
                Rlaset("G", n, n, zero, one, v, ldv);
            }
            if (wntuf) {
                Rlaset("G", n, 1, zero, zero, work, n);
                Rlaset("G", m, n, zero, one, u, ldu);
            }
            for (p = 1; p <= n; p = p + 1) {
                iwork[p - 1] = p;
            }
            if (rowprm) {
                for (p = n + 1; p <= n + m - 1; p = p + 1) {
                    iwork[p - 1] = p - n;
                }
            }
            if (conda) {
                rwork[1 - 1] = -1;
            }
            rwork[2 - 1] = -1;
            return;
        }
        //
        if (rwork[1 - 1] > big / sqrt(castREAL(m))) {
            //               .. to prevent overflow in the QR factorization, scale the
            //               matrix by 1/sqrt(M) if too large entry detected
            Rlascl("G", 0, 0, sqrt(castREAL(m)), one, m, n, a, lda, ierr);
            ascaled = true;
        }
        Rlaswp(n, a, lda, 1, m - 1, &iwork[(n + 1) - 1], 1);
    }
    //
    //    .. At this stage, preemptive scaling is done only to avoid column
    //    norms overflows during the QR factorization. The SVD procedure should
    //    have its own scaling to save the singular values from overflows and
    //    underflows. That depends on the SVD procedure.
    //
    if (!rowprm) {
        rtmp = Rlange("M", m, n, a, lda, rdummy);
        if ((rtmp != rtmp) || ((rtmp * zero) != zero)) {
            info = -8;
            Mxerbla("Rgesvdq", -info);
            return;
        }
        if (rtmp > big / sqrt(castREAL(m))) {
            //             .. to prevent overflow in the QR factorization, scale the
            //             matrix by 1/sqrt(M) if too large entry detected
            Rlascl("G", 0, 0, sqrt(castREAL(m)), one, m, n, a, lda, ierr);
            ascaled = true;
        }
    }
    //
    //     .. QR factorization with column pivoting
    //
    //     A * P = Q * [ R ]
    //                 [ 0 ]
    //
    for (p = 1; p <= n; p = p + 1) {
        //        .. all columns are free columns
        iwork[p - 1] = 0;
    }
    Rgeqp3(m, n, a, lda, iwork, work, &work[(n + 1) - 1], lwork - n, ierr);
    //
    //    If the user requested accuracy level allows truncation in the
    //    computed upper triangular factor, the matrix R is examined and,
    //    if possible, replaced with its leading upper trapezoidal part.
    //
    epsln = Rlamch("E");
    sfmin = Rlamch("S");
    //     SMALL = SFMIN / EPSLN
    nr = n;
    //
    if (accla) {
        //
        //        Standard absolute error bound suffices. All sigma_i with
        //        sigma_i < N*EPS*||A||_F are flushed to zero. This is an
        //        aggressive enforcement of lower numerical rank by introducing a
        //        backward error of the order of N*EPS*||A||_F.
        nr = 1;
        rtmp = sqrt(castREAL(n)) * epsln;
        for (p = 2; p <= n; p = p + 1) {
            if (abs(a[(p - 1) + (p - 1) * lda]) < (rtmp * abs(a[(1 - 1)]))) {
                goto statement_3002;
            }
            nr++;
        }
    statement_3002:;
        //
    } else if (acclm) {
        //        .. similarly as above, only slightly more gentle (less aggressive).
        //        Sudden drop on the diagonal of R is used as the criterion for being
        //        close-to-rank-deficient. The threshold is set to EPSLN=DLAMCH('E').
        //        [[This can be made more flexible by replacing this hard-coded value
        //        with a user specified threshold.]] Also, the values that underflow
        //        will be truncated.
        nr = 1;
        for (p = 2; p <= n; p = p + 1) {
            if ((abs(a[(p - 1) + (p - 1) * lda]) < (epsln * abs(a[((p - 1) - 1) + ((p - 1) - 1) * lda]))) || (abs(a[(p - 1) + (p - 1) * lda]) < sfmin)) {
                goto statement_3402;
            }
            nr++;
        }
    statement_3402:;
        //
    } else {
        //        .. RRQR not authorized to determine numerical rank except in the
        //        obvious case of zero pivots.
        //        .. inspect R for exact zeros on the diagonal;
        //        R(i,i)=0 => R(i:N,i:N)=0.
        nr = 1;
        for (p = 2; p <= n; p = p + 1) {
            if (abs(a[(p - 1) + (p - 1) * lda]) == zero) {
                goto statement_3502;
            }
            nr++;
        }
    statement_3502:
        //
        if (conda) {
            //           Estimate the scaled condition number of A. Use the fact that it is
            //           the same as the scaled condition number of R.
            //              .. V is used as workspace
            Rlacpy("U", n, n, a, lda, v, ldv);
            //              Only the leading NR x NR submatrix of the triangular factor
            //              is considered. Only if NR=N will this give a reliable error
            //              bound. However, even for NR < N, this can be used on an
            //              expert level and obtain useful information in the sense of
            //              perturbation theory.
            for (p = 1; p <= nr; p = p + 1) {
                rtmp = Rnrm2(p, &v[(p - 1) * ldv], 1);
                Rscal(p, one / rtmp, &v[(p - 1) * ldv], 1);
            }
            if (!(lsvec || rsvec)) {
                Rpocon("U", nr, v, ldv, one, rtmp, work, &iwork[(n + iwoff) - 1], ierr);
            } else {
                Rpocon("U", nr, v, ldv, one, rtmp, &work[(n + 1) - 1], &iwork[(n + iwoff) - 1], ierr);
            }
            sconda = one / sqrt(rtmp);
            //           For NR=N, SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1),
            //           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
            //           See the reference [1] for more details.
        }
        //
    }
    //
    if (wntur) {
        n1 = nr;
    } else if (wntus || wntuf) {
        n1 = n;
    } else if (wntua) {
        n1 = m;
    }
    //
    if (!(rsvec || lsvec)) {
        //.......................................................................
        //        .. only the singular values are requested
        //.......................................................................
        if (rtrans) {
            //
            //         .. compute the singular values of R**T = [A](1:NR,1:N)**T
            //           .. set the lower triangle of [A] to [A](1:NR,1:N)**T and
            //           the upper triangle of [A] to zero.
            for (p = 1; p <= min(n, nr); p = p + 1) {
                for (q = p + 1; q <= n; q = q + 1) {
                    a[(q - 1) + (p - 1) * lda] = a[(p - 1) + (q - 1) * lda];
                    if (q <= nr) {
                        a[(p - 1) + (q - 1) * lda] = zero;
                    }
                }
            }
            //
            Rgesvd("N", "N", n, nr, a, lda, s, u, ldu, v, ldv, work, lwork, info);
            //
        } else {
            //
            //           .. compute the singular values of R = [A](1:NR,1:N)
            //
            if (nr > 1) {
                Rlaset("L", nr - 1, nr - 1, zero, zero, &a[(2 - 1)], lda);
            }
            Rgesvd("N", "N", nr, n, a, lda, s, u, ldu, v, ldv, work, lwork, info);
            //
        }
        //
    } else if (lsvec && (!rsvec)) {
        //.......................................................................
        //       .. the singular values and the left singular vectors requested
        //.......................................................................""""""""
        if (rtrans) {
            //            .. apply Rgesvd to R**T
            //            .. copy R**T into [U] and overwrite [U] with the right singular
            //            vectors of R
            for (p = 1; p <= nr; p = p + 1) {
                for (q = p; q <= n; q = q + 1) {
                    u[(q - 1) + (p - 1) * ldu] = a[(p - 1) + (q - 1) * lda];
                }
            }
            if (nr > 1) {
                Rlaset("U", nr - 1, nr - 1, zero, zero, &u[(2 - 1) * ldu], ldu);
            }
            //           .. the left singular vectors not computed, the NR right singular
            //           vectors overwrite [U](1:NR,1:NR) as transposed. These
            //           will be pre-multiplied by Q to build the left singular vectors of A.
            Rgesvd("N", "O", n, nr, u, ldu, s, u, ldu, u, ldu, &work[(n + 1) - 1], lwork - n, info);
            //
            for (p = 1; p <= nr; p = p + 1) {
                for (q = p + 1; q <= nr; q = q + 1) {
                    rtmp = u[(q - 1) + (p - 1) * ldu];
                    u[(q - 1) + (p - 1) * ldu] = u[(p - 1) + (q - 1) * ldu];
                    u[(p - 1) + (q - 1) * ldu] = rtmp;
                }
            }
            //
        } else {
            //            .. apply Rgesvd to R
            //            .. copy R into [U] and overwrite [U] with the left singular vectors
            Rlacpy("U", nr, n, a, lda, u, ldu);
            if (nr > 1) {
                Rlaset("L", nr - 1, nr - 1, zero, zero, &u[(2 - 1)], ldu);
            }
            //            .. the right singular vectors not computed, the NR left singular
            //            vectors overwrite [U](1:NR,1:NR)
            Rgesvd("O", "N", nr, n, u, ldu, s, u, ldu, v, ldv, &work[(n + 1) - 1], lwork - n, info);
            //               .. now [U](1:NR,1:NR) contains the NR left singular vectors of
            //               R. These will be pre-multiplied by Q to build the left singular
            //               vectors of A.
        }
        //
        //              (M x NR) or (M x N) or (M x M).
        if ((nr < m) && (!wntuf)) {
            Rlaset("A", m - nr, nr, zero, zero, &u[((nr + 1) - 1)], ldu);
            if (nr < n1) {
                Rlaset("A", nr, n1 - nr, zero, zero, &u[((nr + 1) - 1) * ldu], ldu);
                Rlaset("A", m - nr, n1 - nr, zero, one, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
            }
        }
        //
        //           The Q matrix from the first QRF is built into the left singular
        //           vectors matrix U.
        //
        if (!wntuf) {
            Rormqr("L", "N", m, n1, n, a, lda, work, u, ldu, &work[(n + 1) - 1], lwork - n, ierr);
        }
        if (rowprm && !wntuf) {
            Rlaswp(n1, u, ldu, 1, m - 1, &iwork[(n + 1) - 1], -1);
        }
        //
    } else if (rsvec && (!lsvec)) {
        //.......................................................................
        //       .. the singular values and the right singular vectors requested
        //.......................................................................
        if (rtrans) {
            //            .. apply Rgesvd to R**T
            //            .. copy R**T into V and overwrite V with the left singular vectors
            for (p = 1; p <= nr; p = p + 1) {
                for (q = p; q <= n; q = q + 1) {
                    v[(q - 1) + (p - 1) * ldv] = (a[(p - 1) + (q - 1) * lda]);
                }
            }
            if (nr > 1) {
                Rlaset("U", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
            }
            //           .. the left singular vectors of R**T overwrite V, the right singular
            //           vectors not computed
            if (wntvr || (nr == n)) {
                Rgesvd("O", "N", n, nr, v, ldv, s, u, ldu, u, ldu, &work[(n + 1) - 1], lwork - n, info);
                //
                for (p = 1; p <= nr; p = p + 1) {
                    for (q = p + 1; q <= nr; q = q + 1) {
                        rtmp = v[(q - 1) + (p - 1) * ldv];
                        v[(q - 1) + (p - 1) * ldv] = v[(p - 1) + (q - 1) * ldv];
                        v[(p - 1) + (q - 1) * ldv] = rtmp;
                    }
                }
                //
                if (nr < n) {
                    for (p = 1; p <= nr; p = p + 1) {
                        for (q = nr + 1; q <= n; q = q + 1) {
                            v[(p - 1) + (q - 1) * ldv] = v[(q - 1) + (p - 1) * ldv];
                        }
                    }
                }
                Rlapmt(false, nr, n, v, ldv, iwork);
            } else {
                //               .. need all N right singular vectors and NR < N
                //               [!] This is simple implementation that augments [V](1:N,1:NR)
                //               by padding a zero block. In the case NR << N, a more efficient
                //               way is to first use the QR factorization. For more details
                //               how to implement this, see the " FULL SVD " branch.
                Rlaset("G", n, n - nr, zero, zero, &v[((nr + 1) - 1) * ldv], ldv);
                Rgesvd("O", "N", n, n, v, ldv, s, u, ldu, u, ldu, &work[(n + 1) - 1], lwork - n, info);
                //
                for (p = 1; p <= n; p = p + 1) {
                    for (q = p + 1; q <= n; q = q + 1) {
                        rtmp = v[(q - 1) + (p - 1) * ldv];
                        v[(q - 1) + (p - 1) * ldv] = v[(p - 1) + (q - 1) * ldv];
                        v[(p - 1) + (q - 1) * ldv] = rtmp;
                    }
                }
                Rlapmt(false, n, n, v, ldv, iwork);
            }
            //
        } else {
            //            .. aply Rgesvd to R
            //            .. copy R into V and overwrite V with the right singular vectors
            Rlacpy("U", nr, n, a, lda, v, ldv);
            if (nr > 1) {
                Rlaset("L", nr - 1, nr - 1, zero, zero, &v[(2 - 1)], ldv);
            }
            //            .. the right singular vectors overwrite V, the NR left singular
            //            vectors stored in U(1:NR,1:NR)
            if (wntvr || (nr == n)) {
                Rgesvd("N", "O", nr, n, v, ldv, s, u, ldu, v, ldv, &work[(n + 1) - 1], lwork - n, info);
                Rlapmt(false, nr, n, v, ldv, iwork);
                //               .. now [V](1:NR,1:N) contains V(1:N,1:NR)**T
            } else {
                //               .. need all N right singular vectors and NR < N
                //               [!] This is simple implementation that augments [V](1:NR,1:N)
                //               by padding a zero block. In the case NR << N, a more efficient
                //               way is to first use the LQ factorization. For more details
                //               how to implement this, see the " FULL SVD " branch.
                Rlaset("G", n - nr, n, zero, zero, &v[((nr + 1) - 1)], ldv);
                Rgesvd("N", "O", n, n, v, ldv, s, u, ldu, v, ldv, &work[(n + 1) - 1], lwork - n, info);
                Rlapmt(false, n, n, v, ldv, iwork);
            }
            //            .. now [V] contains the transposed matrix of the right singular
            //            vectors of A.
        }
        //
    } else {
        //.......................................................................
        //       .. FULL SVD requested
        //.......................................................................
        if (rtrans) {
            //
            //            .. apply Rgesvd to R**T [[this option is left for R&D&T]]
            //
            if (wntvr || (nr == n)) {
                //            .. copy R**T into [V] and overwrite [V] with the left singular
                //            vectors of R**T
                for (p = 1; p <= nr; p = p + 1) {
                    for (q = p; q <= n; q = q + 1) {
                        v[(q - 1) + (p - 1) * ldv] = a[(p - 1) + (q - 1) * lda];
                    }
                }
                if (nr > 1) {
                    Rlaset("U", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
                }
                //
                //           .. the left singular vectors of R**T overwrite [V], the NR right
                //           singular vectors of R**T stored in [U](1:NR,1:NR) as transposed
                Rgesvd("O", "A", n, nr, v, ldv, s, v, ldv, u, ldu, &work[(n + 1) - 1], lwork - n, info);
                //              .. assemble V
                for (p = 1; p <= nr; p = p + 1) {
                    for (q = p + 1; q <= nr; q = q + 1) {
                        rtmp = v[(q - 1) + (p - 1) * ldv];
                        v[(q - 1) + (p - 1) * ldv] = v[(p - 1) + (q - 1) * ldv];
                        v[(p - 1) + (q - 1) * ldv] = rtmp;
                    }
                }
                if (nr < n) {
                    for (p = 1; p <= nr; p = p + 1) {
                        for (q = nr + 1; q <= n; q = q + 1) {
                            v[(p - 1) + (q - 1) * ldv] = v[(q - 1) + (p - 1) * ldv];
                        }
                    }
                }
                Rlapmt(false, nr, n, v, ldv, iwork);
                //
                for (p = 1; p <= nr; p = p + 1) {
                    for (q = p + 1; q <= nr; q = q + 1) {
                        rtmp = u[(q - 1) + (p - 1) * ldu];
                        u[(q - 1) + (p - 1) * ldu] = u[(p - 1) + (q - 1) * ldu];
                        u[(p - 1) + (q - 1) * ldu] = rtmp;
                    }
                }
                //
                if ((nr < m) && !(wntuf)) {
                    Rlaset("A", m - nr, nr, zero, zero, &u[((nr + 1) - 1)], ldu);
                    if (nr < n1) {
                        Rlaset("A", nr, n1 - nr, zero, zero, &u[((nr + 1) - 1) * ldu], ldu);
                        Rlaset("A", m - nr, n1 - nr, zero, one, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                    }
                }
                //
            } else {
                //               .. need all N right singular vectors and NR < N
                //            .. copy R**T into [V] and overwrite [V] with the left singular
                //            vectors of R**T
                //               [[The optimal ratio N/NR for using QRF instead of padding
                //                 with zeros. Here hard coded to 2; it must be at least
                //                 two due to work space constraints.]]
                //               OPTRATIO = iMlaenv(6, 'Rgesvd', 'S' // 'O', NR,N,0,0)
                //               OPTRATIO = MAX( OPTRATIO, 2 )
                optratio = 2;
                if (optratio * nr > n) {
                    for (p = 1; p <= nr; p = p + 1) {
                        for (q = p; q <= n; q = q + 1) {
                            v[(q - 1) + (p - 1) * ldv] = a[(p - 1) + (q - 1) * lda];
                        }
                    }
                    if (nr > 1) {
                        Rlaset("U", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
                    }
                    //
                    Rlaset("A", n, n - nr, zero, zero, &v[((nr + 1) - 1) * ldv], ldv);
                    Rgesvd("O", "A", n, n, v, ldv, s, v, ldv, u, ldu, &work[(n + 1) - 1], lwork - n, info);
                    //
                    for (p = 1; p <= n; p = p + 1) {
                        for (q = p + 1; q <= n; q = q + 1) {
                            rtmp = v[(q - 1) + (p - 1) * ldv];
                            v[(q - 1) + (p - 1) * ldv] = v[(p - 1) + (q - 1) * ldv];
                            v[(p - 1) + (q - 1) * ldv] = rtmp;
                        }
                    }
                    Rlapmt(false, n, n, v, ldv, iwork);
                    //              (M x N1), i.e. (M x N) or (M x M).
                    //
                    for (p = 1; p <= n; p = p + 1) {
                        for (q = p + 1; q <= n; q = q + 1) {
                            rtmp = u[(q - 1) + (p - 1) * ldu];
                            u[(q - 1) + (p - 1) * ldu] = u[(p - 1) + (q - 1) * ldu];
                            u[(p - 1) + (q - 1) * ldu] = rtmp;
                        }
                    }
                    //
                    if ((n < m) && !(wntuf)) {
                        Rlaset("A", m - n, n, zero, zero, &u[((n + 1) - 1)], ldu);
                        if (n < n1) {
                            Rlaset("A", n, n1 - n, zero, zero, &u[((n + 1) - 1) * ldu], ldu);
                            Rlaset("A", m - n, n1 - n, zero, one, &u[((n + 1) - 1) + ((n + 1) - 1) * ldu], ldu);
                        }
                    }
                } else {
                    //                  .. copy R**T into [U] and overwrite [U] with the right
                    //                  singular vectors of R
                    for (p = 1; p <= nr; p = p + 1) {
                        for (q = p; q <= n; q = q + 1) {
                            u[(q - 1) + ((nr + p) - 1) * ldu] = a[(p - 1) + (q - 1) * lda];
                        }
                    }
                    if (nr > 1) {
                        Rlaset("U", nr - 1, nr - 1, zero, zero, &u[((nr + 2) - 1) * ldu], ldu);
                    }
                    Rgeqrf(n, nr, &u[((nr + 1) - 1) * ldu], ldu, &work[(n + 1) - 1], &work[(n + nr + 1) - 1], lwork - n - nr, ierr);
                    for (p = 1; p <= nr; p = p + 1) {
                        for (q = 1; q <= n; q = q + 1) {
                            v[(q - 1) + (p - 1) * ldv] = u[(p - 1) + ((nr + q) - 1) * ldu];
                        }
                    }
                    Rlaset("U", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
                    Rgesvd("S", "O", nr, nr, v, ldv, s, u, ldu, v, ldv, &work[(n + nr + 1) - 1], lwork - n - nr, info);
                    Rlaset("A", n - nr, nr, zero, zero, &v[((nr + 1) - 1)], ldv);
                    Rlaset("A", nr, n - nr, zero, zero, &v[((nr + 1) - 1) * ldv], ldv);
                    Rlaset("A", n - nr, n - nr, zero, one, &v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
                    Rormqr("R", "C", n, n, nr, &u[((nr + 1) - 1) * ldu], ldu, &work[(n + 1) - 1], v, ldv, &work[(n + nr + 1) - 1], lwork - n - nr, ierr);
                    Rlapmt(false, n, n, v, ldv, iwork);
                    //                 (M x NR) or (M x N) or (M x M).
                    if ((nr < m) && !(wntuf)) {
                        Rlaset("A", m - nr, nr, zero, zero, &u[((nr + 1) - 1)], ldu);
                        if (nr < n1) {
                            Rlaset("A", nr, n1 - nr, zero, zero, &u[((nr + 1) - 1) * ldu], ldu);
                            Rlaset("A", m - nr, n1 - nr, zero, one, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                        }
                    }
                }
            }
            //
        } else {
            //
            //            .. apply Rgesvd to R [[this is the recommended option]]
            //
            if (wntvr || (nr == n)) {
                //                .. copy R into [V] and overwrite V with the right singular vectors
                Rlacpy("U", nr, n, a, lda, v, ldv);
                if (nr > 1) {
                    Rlaset("L", nr - 1, nr - 1, zero, zero, &v[(2 - 1)], ldv);
                }
                //               .. the right singular vectors of R overwrite [V], the NR left
                //               singular vectors of R stored in [U](1:NR,1:NR)
                Rgesvd("S", "O", nr, n, v, ldv, s, u, ldu, v, ldv, &work[(n + 1) - 1], lwork - n, info);
                Rlapmt(false, nr, n, v, ldv, iwork);
                //               .. now [V](1:NR,1:N) contains V(1:N,1:NR)**T
                //              (M x NR) or (M x N) or (M x M).
                if ((nr < m) && !(wntuf)) {
                    Rlaset("A", m - nr, nr, zero, zero, &u[((nr + 1) - 1)], ldu);
                    if (nr < n1) {
                        Rlaset("A", nr, n1 - nr, zero, zero, &u[((nr + 1) - 1) * ldu], ldu);
                        Rlaset("A", m - nr, n1 - nr, zero, one, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                    }
                }
                //
            } else {
                //              .. need all N right singular vectors and NR < N
                //              .. the requested number of the left singular vectors
                //               is then N1 (N or M)
                //               [[The optimal ratio N/NR for using LQ instead of padding
                //                 with zeros. Here hard coded to 2; it must be at least
                //                 two due to work space constraints.]]
                //               OPTRATIO = iMlaenv(6, 'Rgesvd', 'S' // 'O', NR,N,0,0)
                //               OPTRATIO = MAX( OPTRATIO, 2 )
                optratio = 2;
                if (optratio * nr > n) {
                    Rlacpy("U", nr, n, a, lda, v, ldv);
                    if (nr > 1) {
                        Rlaset("L", nr - 1, nr - 1, zero, zero, &v[(2 - 1)], ldv);
                    }
                    //              .. the right singular vectors of R overwrite [V], the NR left
                    //                 singular vectors of R stored in [U](1:NR,1:NR)
                    Rlaset("A", n - nr, n, zero, zero, &v[((nr + 1) - 1)], ldv);
                    Rgesvd("S", "O", n, n, v, ldv, s, u, ldu, v, ldv, &work[(n + 1) - 1], lwork - n, info);
                    Rlapmt(false, n, n, v, ldv, iwork);
                    //                 .. now [V] contains the transposed matrix of the right
                    //                 singular vectors of A. The leading N left singular vectors
                    //                 are in [U](1:N,1:N)
                    //                 (M x N1), i.e. (M x N) or (M x M).
                    if ((n < m) && !(wntuf)) {
                        Rlaset("A", m - n, n, zero, zero, &u[((n + 1) - 1)], ldu);
                        if (n < n1) {
                            Rlaset("A", n, n1 - n, zero, zero, &u[((n + 1) - 1) * ldu], ldu);
                            Rlaset("A", m - n, n1 - n, zero, one, &u[((n + 1) - 1) + ((n + 1) - 1) * ldu], ldu);
                        }
                    }
                } else {
                    Rlacpy("U", nr, n, a, lda, &u[((nr + 1) - 1)], ldu);
                    if (nr > 1) {
                        Rlaset("L", nr - 1, nr - 1, zero, zero, &u[((nr + 2) - 1)], ldu);
                    }
                    Rgelqf(nr, n, &u[((nr + 1) - 1)], ldu, &work[(n + 1) - 1], &work[(n + nr + 1) - 1], lwork - n - nr, ierr);
                    Rlacpy("L", nr, nr, &u[((nr + 1) - 1)], ldu, v, ldv);
                    if (nr > 1) {
                        Rlaset("U", nr - 1, nr - 1, zero, zero, &v[(2 - 1) * ldv], ldv);
                    }
                    Rgesvd("S", "O", nr, nr, v, ldv, s, u, ldu, v, ldv, &work[(n + nr + 1) - 1], lwork - n - nr, info);
                    Rlaset("A", n - nr, nr, zero, zero, &v[((nr + 1) - 1)], ldv);
                    Rlaset("A", nr, n - nr, zero, zero, &v[((nr + 1) - 1) * ldv], ldv);
                    Rlaset("A", n - nr, n - nr, zero, one, &v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
                    Rormlq("R", "N", n, n, nr, &u[((nr + 1) - 1)], ldu, &work[(n + 1) - 1], v, ldv, &work[(n + nr + 1) - 1], lwork - n - nr, ierr);
                    Rlapmt(false, n, n, v, ldv, iwork);
                    //              (M x NR) or (M x N) or (M x M).
                    if ((nr < m) && !(wntuf)) {
                        Rlaset("A", m - nr, nr, zero, zero, &u[((nr + 1) - 1)], ldu);
                        if (nr < n1) {
                            Rlaset("A", nr, n1 - nr, zero, zero, &u[((nr + 1) - 1) * ldu], ldu);
                            Rlaset("A", m - nr, n1 - nr, zero, one, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                        }
                    }
                }
            }
            //        .. end of the "R**T or R" branch
        }
        //
        //           The Q matrix from the first QRF is built into the left singular
        //           vectors matrix U.
        //
        if (!wntuf) {
            Rormqr("L", "N", m, n1, n, a, lda, work, u, ldu, &work[(n + 1) - 1], lwork - n, ierr);
        }
        if (rowprm && !wntuf) {
            Rlaswp(n1, u, ldu, 1, m - 1, &iwork[(n + 1) - 1], -1);
        }
        //
        //     ... end of the "full SVD" branch
    }
    //
    //     Check whether some singular values are returned as zeros, e.g.
    //     due to underflow, and update the numerical rank.
    p = nr;
    for (q = p; q >= 1; q = q - 1) {
        if (s[q - 1] > zero) {
            goto statement_4002;
        }
        nr = nr - 1;
    }
statement_4002:
    //
    //     .. if numerical rank deficiency is detected, the truncated
    //     singular values are set to zero.
    if (nr < n) {
        Rlaset("G", n - nr, 1, zero, zero, &s[(nr + 1) - 1], n);
    }
    //     .. undo scaling; this may cause overflow in the largest singular
    //     values.
    if (ascaled) {
        Rlascl("G", 0, 0, one, sqrt(castREAL(m)), nr, 1, s, n, ierr);
    }
    if (conda) {
        rwork[1 - 1] = sconda;
    }
    rwork[2 - 1] = p - nr;
    //     .. p-NR is the number of singular values that are computed as
    //     exact zeros in Rgesvd() applied to the (possibly truncated)
    //     full row rank triangular (trapezoidal) factor of A.
    numrank = nr;
    //
    //     End of Rgesvdq
    //
}
