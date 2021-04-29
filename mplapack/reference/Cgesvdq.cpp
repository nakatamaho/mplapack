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

void Cgesvdq(const char *joba, const char *jobp, const char *jobr, const char *jobu, const char *jobv, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL *s, COMPLEX *u, INTEGER const ldu, COMPLEX *v, INTEGER const ldv, INTEGER &numrank, INTEGER *iwork, INTEGER const liwork, COMPLEX *cwork, INTEGER const lcwork, REAL *rwork, INTEGER const lrwork, INTEGER &info) {
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
    INTEGER lwunq = 0;
    INTEGER lwcon = 0;
    INTEGER lwsvd = 0;
    COMPLEX cdummy[1];
    REAL rdummy[1];
    INTEGER ierr = 0;
    INTEGER lwrk_Cgeqp3 = 0;
    INTEGER lwrk_Cunmqr = 0;
    INTEGER minwrk = 0;
    INTEGER optwrk = 0;
    INTEGER lwrk_Cgesvd = 0;
    INTEGER lwqrf = 0;
    INTEGER lwsvd2 = 0;
    INTEGER lwunq2 = 0;
    INTEGER minwrk2 = 0;
    INTEGER lwlqf = 0;
    INTEGER lwunlq = 0;
    INTEGER lwrk_Cgeqrf = 0;
    INTEGER lwrk_Cgesvd2 = 0;
    INTEGER lwrk_Cunmqr2 = 0;
    INTEGER optwrk2 = 0;
    INTEGER lwrk_Cgelqf = 0;
    INTEGER lwrk_Cunmlq = 0;
    REAL big = 0.0;
    bool ascaled = false;
    INTEGER p = 0;
    const REAL zero = 0.0;
    INTEGER q = 0;
    REAL rtmp = 0.0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const REAL one = 1.0;
    REAL epsln = 0.0;
    REAL sfmin = 0.0;
    INTEGER nr = 0;
    REAL sconda = 0.0;
    INTEGER n1 = 0;
    COMPLEX ctmp = 0.0;
    INTEGER optratio = 0;
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays
    //     ..
    //     .. External Subroutines (BLAS, LAPACK)
    //     ..
    //     .. External Functions (BLAS, LAPACK)
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
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
        iminwrk = max((INTEGER)1, n + m - 1);
        rminwrk = max({(INTEGER)2, m, 5 * n});
    } else {
        iminwrk = max((INTEGER)1, n);
        rminwrk = max((INTEGER)2, 5 * n);
    }
    lquery = (liwork == -1 || lcwork == -1 || lrwork == -1);
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
        //        values of LCWORK are written with a lot of redundancy and
        //        can be simplified. However, this detailed form is easier for
        //        maintenance and modifications of the code.]]
        //
        //        .. minimal workspace length for Cgeqp3 of an M x N matrix
        lwqp3 = n + 1;
        //        .. minimal workspace length for Cunmqr to build left singular vectors
        if (wntus || wntur) {
            lwunq = max(n, 1);
        } else if (wntua) {
            lwunq = max(m, 1);
        }
        //        .. minimal workspace length for Cpocon of an N x N matrix
        lwcon = 2 * n;
        //        .. Cgesvd of an N x N matrix
        lwsvd = max(3 * n, 1);
        if (lquery) {
            Cgeqp3(m, n, a, lda, iwork, cdummy, cdummy, -1, rdummy, ierr);
            lwrk_Cgeqp3 = castINTEGER(cdummy[1 - 1].real());
            if (wntus || wntur) {
                Cunmqr("L", "N", m, n, n, a, lda, cdummy, u, ldu, cdummy, -1, ierr);
                lwrk_Cunmqr = castINTEGER(cdummy[1 - 1].real());
            } else if (wntua) {
                Cunmqr("L", "N", m, m, n, a, lda, cdummy, u, ldu, cdummy, -1, ierr);
                lwrk_Cunmqr = castINTEGER(cdummy[1 - 1].real());
            } else {
                lwrk_Cunmqr = 0;
            }
        }
        minwrk = 2;
        optwrk = 2;
        if (!(lsvec || rsvec)) {
            //            .. minimal and optimal sizes of the complex workspace if
            //            only the singular values are requested
            if (conda) {
                minwrk = max({n + lwqp3, lwcon, lwsvd});
            } else {
                minwrk = max(n + lwqp3, lwsvd);
            }
            if (lquery) {
                Cgesvd("N", "N", n, n, a, lda, s, u, ldu, v, ldv, cdummy, -1, rdummy, ierr);
                lwrk_Cgesvd = castINTEGER(cdummy[1 - 1].real());
                if (conda) {
                    optwrk = max({n + lwrk_Cgeqp3, n + lwcon, lwrk_Cgesvd});
                } else {
                    optwrk = max(n + lwrk_Cgeqp3, lwrk_Cgesvd);
                }
            }
        } else if (lsvec && (!rsvec)) {
            //            .. minimal and optimal sizes of the complex workspace if the
            //            singular values and the left singular vectors are requested
            if (conda) {
                minwrk = n + max({lwqp3, lwcon, lwsvd, lwunq});
            } else {
                minwrk = n + max({lwqp3, lwsvd, lwunq});
            }
            if (lquery) {
                if (rtrans) {
                    Cgesvd("N", "O", n, n, a, lda, s, u, ldu, v, ldv, cdummy, -1, rdummy, ierr);
                } else {
                    Cgesvd("O", "N", n, n, a, lda, s, u, ldu, v, ldv, cdummy, -1, rdummy, ierr);
                }
                lwrk_Cgesvd = castINTEGER(cdummy[1 - 1].real());
                if (conda) {
                    optwrk = n + max({lwrk_Cgeqp3, lwcon, lwrk_Cgesvd, lwrk_Cunmqr});
                } else {
                    optwrk = n + max({lwrk_Cgeqp3, lwrk_Cgesvd, lwrk_Cunmqr});
                }
            }
        } else if (rsvec && (!lsvec)) {
            //            .. minimal and optimal sizes of the complex workspace if the
            //            singular values and the right singular vectors are requested
            if (conda) {
                minwrk = n + max({lwqp3, lwcon, lwsvd});
            } else {
                minwrk = n + max(lwqp3, lwsvd);
            }
            if (lquery) {
                if (rtrans) {
                    Cgesvd("O", "N", n, n, a, lda, s, u, ldu, v, ldv, cdummy, -1, rdummy, ierr);
                } else {
                    Cgesvd("N", "O", n, n, a, lda, s, u, ldu, v, ldv, cdummy, -1, rdummy, ierr);
                }
                lwrk_Cgesvd = castINTEGER(cdummy[1 - 1].real());
                if (conda) {
                    optwrk = n + max({lwrk_Cgeqp3, lwcon, lwrk_Cgesvd});
                } else {
                    optwrk = n + max(lwrk_Cgeqp3, lwrk_Cgesvd);
                }
            }
        } else {
            //            .. minimal and optimal sizes of the complex workspace if the
            //            full SVD is requested
            if (rtrans) {
                minwrk = max({lwqp3, lwsvd, lwunq});
                if (conda) {
                    minwrk = max(minwrk, lwcon);
                }
                minwrk += n;
                if (wntva) {
                    //                   .. minimal workspace length for N x N/2 Cgeqrf
                    lwqrf = max(n / 2, 1);
                    //                   .. minimal workspace length for N/2 x N/2 Cgesvd
                    lwsvd2 = max(3 * (n / 2), 1);
                    lwunq2 = max(n, 1);
                    minwrk2 = max({lwqp3, n / 2 + lwqrf, n / 2 + lwsvd2, n / 2 + lwunq2, lwunq});
                    if (conda) {
                        minwrk2 = max(minwrk2, lwcon);
                    }
                    minwrk2 += n;
                    minwrk = max(minwrk, minwrk2);
                }
            } else {
                minwrk = max({lwqp3, lwsvd, lwunq});
                if (conda) {
                    minwrk = max(minwrk, lwcon);
                }
                minwrk += n;
                if (wntva) {
                    //                   .. minimal workspace length for N/2 x N Cgelqf
                    lwlqf = max(n / 2, 1);
                    lwsvd2 = max(3 * (n / 2), 1);
                    lwunlq = max(n, 1);
                    minwrk2 = max({lwqp3, n / 2 + lwlqf, n / 2 + lwsvd2, n / 2 + lwunlq, lwunq});
                    if (conda) {
                        minwrk2 = max(minwrk2, lwcon);
                    }
                    minwrk2 += n;
                    minwrk = max(minwrk, minwrk2);
                }
            }
            if (lquery) {
                if (rtrans) {
                    Cgesvd("O", "A", n, n, a, lda, s, u, ldu, v, ldv, cdummy, -1, rdummy, ierr);
                    lwrk_Cgesvd = castINTEGER(cdummy[1 - 1].real());
                    optwrk = max({lwrk_Cgeqp3, lwrk_Cgesvd, lwrk_Cunmqr});
                    if (conda) {
                        optwrk = max(optwrk, lwcon);
                    }
                    optwrk += n;
                    if (wntva) {
                        Cgeqrf(n, n / 2, u, ldu, cdummy, cdummy, -1, ierr);
                        lwrk_Cgeqrf = castINTEGER(cdummy[1 - 1].real());
                        Cgesvd("S", "O", n / 2, n / 2, v, ldv, s, u, ldu, v, ldv, cdummy, -1, rdummy, ierr);
                        lwrk_Cgesvd2 = castINTEGER(cdummy[1 - 1].real());
                        Cunmqr("R", "C", n, n, n / 2, u, ldu, cdummy, v, ldv, cdummy, -1, ierr);
                        lwrk_Cunmqr2 = castINTEGER(cdummy[1 - 1].real());
                        optwrk2 = max({lwrk_Cgeqp3, n / 2 + lwrk_Cgeqrf, n / 2 + lwrk_Cgesvd2, n / 2 + lwrk_Cunmqr2});
                        if (conda) {
                            optwrk2 = max(optwrk2, lwcon);
                        }
                        optwrk2 += n;
                        optwrk = max(optwrk, optwrk2);
                    }
                } else {
                    Cgesvd("S", "O", n, n, a, lda, s, u, ldu, v, ldv, cdummy, -1, rdummy, ierr);
                    lwrk_Cgesvd = castINTEGER(cdummy[1 - 1].real());
                    optwrk = max({lwrk_Cgeqp3, lwrk_Cgesvd, lwrk_Cunmqr});
                    if (conda) {
                        optwrk = max(optwrk, lwcon);
                    }
                    optwrk += n;
                    if (wntva) {
                        Cgelqf(n / 2, n, u, ldu, cdummy, cdummy, -1, ierr);
                        lwrk_Cgelqf = castINTEGER(cdummy[1 - 1].real());
                        Cgesvd("S", "O", n / 2, n / 2, v, ldv, s, u, ldu, v, ldv, cdummy, -1, rdummy, ierr);
                        lwrk_Cgesvd2 = castINTEGER(cdummy[1 - 1].real());
                        Cunmlq("R", "N", n, n, n / 2, u, ldu, cdummy, v, ldv, cdummy, -1, ierr);
                        lwrk_Cunmlq = castINTEGER(cdummy[1 - 1].real());
                        optwrk2 = max({lwrk_Cgeqp3, n / 2 + lwrk_Cgelqf, n / 2 + lwrk_Cgesvd2, n / 2 + lwrk_Cunmlq});
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
        if (lcwork < minwrk && (!lquery)) {
            info = -19;
        }
        //
    }
    //
    if (info == 0 && lrwork < rminwrk && !lquery) {
        info = -21;
    }
    if (info != 0) {
        Mxerbla("Cgesvdq", -info);
        return;
    } else if (lquery) {
        //
        //     Return optimal workspace
        //
        iwork[1 - 1] = iminwrk;
        cwork[1 - 1] = optwrk;
        cwork[2 - 1] = minwrk;
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
    if (rowprm) {
        //           .. reordering the rows in decreasing sequence in the
        //           ell-infinity norm - this enhances numerical robustness in
        //           the case of differently scaled rows.
        for (p = 1; p <= m; p = p + 1) {
            //               RWORK(p) = ABS( A(p,iCamax(N,A(p,1),LDA)) )
            //               [[Clange will return NaN if an entry of the p-th row is Nan]]
            rwork[p - 1] = Clange("M", 1, n, &a[(p - 1)], lda, rdummy);
            //               .. check for NaN's and Inf's
            if ((rwork[p - 1] != rwork[p - 1]) || ((rwork[p - 1] * zero) != zero)) {
                info = -8;
                Mxerbla("Cgesvdq", -info);
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
                Claset("G", m, n, czero, cone, u, ldu);
            }
            if (wntua) {
                Claset("G", m, m, czero, cone, u, ldu);
            }
            if (wntva) {
                Claset("G", n, n, czero, cone, v, ldv);
            }
            if (wntuf) {
                Claset("G", n, 1, czero, czero, cwork, n);
                Claset("G", m, n, czero, cone, u, ldu);
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
            Clascl("G", 0, 0, sqrt(castREAL(m)), one, m, n, a, lda, ierr);
            ascaled = true;
        }
        Claswp(n, a, lda, 1, m - 1, &iwork[(n + 1) - 1], 1);
    }
    //
    //    .. At this stage, preemptive scaling is done only to avoid column
    //    norms overflows during the QR factorization. The SVD procedure should
    //    have its own scaling to save the singular values from overflows and
    //    underflows. That depends on the SVD procedure.
    //
    if (!rowprm) {
        rtmp = Clange("M", m, n, a, lda, rwork);
        if ((rtmp != rtmp) || ((rtmp * zero) != zero)) {
            info = -8;
            Mxerbla("Cgesvdq", -info);
            return;
        }
        if (rtmp > big / sqrt(castREAL(m))) {
            //             .. to prevent overflow in the QR factorization, scale the
            //             matrix by 1/sqrt(M) if too large entry detected
            Clascl("G", 0, 0, sqrt(castREAL(m)), one, m, n, a, lda, ierr);
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
    Cgeqp3(m, n, a, lda, iwork, cwork, &cwork[(n + 1) - 1], lcwork - n, rwork, ierr);
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
            Clacpy("U", n, n, a, lda, v, ldv);
            //              Only the leading NR x NR submatrix of the triangular factor
            //              is considered. Only if NR=N will this give a reliable error
            //              bound. However, even for NR < N, this can be used on an
            //              expert level and obtain useful information in the sense of
            //              perturbation theory.
            for (p = 1; p <= nr; p = p + 1) {
                rtmp = RCnrm2(p, &v[(p - 1) * ldv], 1);
                CRscal(p, one / rtmp, &v[(p - 1) * ldv], 1);
            }
            if (!(lsvec || rsvec)) {
                Cpocon("U", nr, v, ldv, one, rtmp, cwork, rwork, ierr);
            } else {
                Cpocon("U", nr, v, ldv, one, rtmp, &cwork[(n + 1) - 1], rwork, ierr);
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
            //         .. compute the singular values of R**H = [A](1:NR,1:N)**H
            //           .. set the lower triangle of [A] to [A](1:NR,1:N)**H and
            //           the upper triangle of [A] to zero.
            for (p = 1; p <= min(n, nr); p = p + 1) {
                a[(p - 1) + (p - 1) * lda] = conj(a[(p - 1) + (p - 1) * lda]);
                for (q = p + 1; q <= n; q = q + 1) {
                    a[(q - 1) + (p - 1) * lda] = conj(a[(p - 1) + (q - 1) * lda]);
                    if (q <= nr) {
                        a[(p - 1) + (q - 1) * lda] = czero;
                    }
                }
            }
            //
            Cgesvd("N", "N", n, nr, a, lda, s, u, ldu, v, ldv, cwork, lcwork, rwork, info);
            //
        } else {
            //
            //           .. compute the singular values of R = [A](1:NR,1:N)
            //
            if (nr > 1) {
                Claset("L", nr - 1, nr - 1, czero, czero, &a[(2 - 1)], lda);
            }
            Cgesvd("N", "N", nr, n, a, lda, s, u, ldu, v, ldv, cwork, lcwork, rwork, info);
            //
        }
        //
    } else if (lsvec && (!rsvec)) {
        //.......................................................................
        //       .. the singular values and the left singular vectors requested
        //.......................................................................""""""""
        if (rtrans) {
            //            .. apply Cgesvd to R**H
            //            .. copy R**H into [U] and overwrite [U] with the right singular
            //            vectors of R
            for (p = 1; p <= nr; p = p + 1) {
                for (q = p; q <= n; q = q + 1) {
                    u[(q - 1) + (p - 1) * ldu] = conj(a[(p - 1) + (q - 1) * lda]);
                }
            }
            if (nr > 1) {
                Claset("U", nr - 1, nr - 1, czero, czero, &u[(2 - 1) * ldu], ldu);
            }
            //           .. the left singular vectors not computed, the NR right singular
            //           vectors overwrite [U](1:NR,1:NR) as conjugate transposed. These
            //           will be pre-multiplied by Q to build the left singular vectors of A.
            Cgesvd("N", "O", n, nr, u, ldu, s, u, ldu, u, ldu, &cwork[(n + 1) - 1], lcwork - n, rwork, info);
            //
            for (p = 1; p <= nr; p = p + 1) {
                u[(p - 1) + (p - 1) * ldu] = conj(u[(p - 1) + (p - 1) * ldu]);
                for (q = p + 1; q <= nr; q = q + 1) {
                    ctmp = conj(u[(q - 1) + (p - 1) * ldu]);
                    u[(q - 1) + (p - 1) * ldu] = conj(u[(p - 1) + (q - 1) * ldu]);
                    u[(p - 1) + (q - 1) * ldu] = ctmp;
                }
            }
            //
        } else {
            //            .. apply Cgesvd to R
            //            .. copy R into [U] and overwrite [U] with the left singular vectors
            Clacpy("U", nr, n, a, lda, u, ldu);
            if (nr > 1) {
                Claset("L", nr - 1, nr - 1, czero, czero, &u[(2 - 1)], ldu);
            }
            //            .. the right singular vectors not computed, the NR left singular
            //            vectors overwrite [U](1:NR,1:NR)
            Cgesvd("O", "N", nr, n, u, ldu, s, u, ldu, v, ldv, &cwork[(n + 1) - 1], lcwork - n, rwork, info);
            //               .. now [U](1:NR,1:NR) contains the NR left singular vectors of
            //               R. These will be pre-multiplied by Q to build the left singular
            //               vectors of A.
        }
        //
        //              (M x NR) or (M x N) or (M x M).
        if ((nr < m) && (!wntuf)) {
            Claset("A", m - nr, nr, czero, czero, &u[((nr + 1) - 1)], ldu);
            if (nr < n1) {
                Claset("A", nr, n1 - nr, czero, czero, &u[((nr + 1) - 1) * ldu], ldu);
                Claset("A", m - nr, n1 - nr, czero, cone, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
            }
        }
        //
        //           The Q matrix from the first QRF is built into the left singular
        //           vectors matrix U.
        //
        if (!wntuf) {
            Cunmqr("L", "N", m, n1, n, a, lda, cwork, u, ldu, &cwork[(n + 1) - 1], lcwork - n, ierr);
        }
        if (rowprm && !wntuf) {
            Claswp(n1, u, ldu, 1, m - 1, &iwork[(n + 1) - 1], -1);
        }
        //
    } else if (rsvec && (!lsvec)) {
        //.......................................................................
        //       .. the singular values and the right singular vectors requested
        //.......................................................................
        if (rtrans) {
            //            .. apply Cgesvd to R**H
            //            .. copy R**H into V and overwrite V with the left singular vectors
            for (p = 1; p <= nr; p = p + 1) {
                for (q = p; q <= n; q = q + 1) {
                    v[(q - 1) + (p - 1) * ldv] = conj(a[(p - 1) + (q - 1) * lda]);
                }
            }
            if (nr > 1) {
                Claset("U", nr - 1, nr - 1, czero, czero, &v[(2 - 1) * ldv], ldv);
            }
            //           .. the left singular vectors of R**H overwrite V, the right singular
            //           vectors not computed
            if (wntvr || (nr == n)) {
                Cgesvd("O", "N", n, nr, v, ldv, s, u, ldu, u, ldu, &cwork[(n + 1) - 1], lcwork - n, rwork, info);
                //
                for (p = 1; p <= nr; p = p + 1) {
                    v[(p - 1) + (p - 1) * ldv] = conj(v[(p - 1) + (p - 1) * ldv]);
                    for (q = p + 1; q <= nr; q = q + 1) {
                        ctmp = conj(v[(q - 1) + (p - 1) * ldv]);
                        v[(q - 1) + (p - 1) * ldv] = conj(v[(p - 1) + (q - 1) * ldv]);
                        v[(p - 1) + (q - 1) * ldv] = ctmp;
                    }
                }
                //
                if (nr < n) {
                    for (p = 1; p <= nr; p = p + 1) {
                        for (q = nr + 1; q <= n; q = q + 1) {
                            v[(p - 1) + (q - 1) * ldv] = conj(v[(q - 1) + (p - 1) * ldv]);
                        }
                    }
                }
                Clapmt(false, nr, n, v, ldv, iwork);
            } else {
                //               .. need all N right singular vectors and NR < N
                //               [!] This is simple implementation that augments [V](1:N,1:NR)
                //               by padding a zero block. In the case NR << N, a more efficient
                //               way is to first use the QR factorization. For more details
                //               how to implement this, see the " FULL SVD " branch.
                Claset("G", n, n - nr, czero, czero, &v[((nr + 1) - 1) * ldv], ldv);
                Cgesvd("O", "N", n, n, v, ldv, s, u, ldu, u, ldu, &cwork[(n + 1) - 1], lcwork - n, rwork, info);
                //
                for (p = 1; p <= n; p = p + 1) {
                    v[(p - 1) + (p - 1) * ldv] = conj(v[(p - 1) + (p - 1) * ldv]);
                    for (q = p + 1; q <= n; q = q + 1) {
                        ctmp = conj(v[(q - 1) + (p - 1) * ldv]);
                        v[(q - 1) + (p - 1) * ldv] = conj(v[(p - 1) + (q - 1) * ldv]);
                        v[(p - 1) + (q - 1) * ldv] = ctmp;
                    }
                }
                Clapmt(false, n, n, v, ldv, iwork);
            }
            //
        } else {
            //            .. aply Cgesvd to R
            //            .. copy R into V and overwrite V with the right singular vectors
            Clacpy("U", nr, n, a, lda, v, ldv);
            if (nr > 1) {
                Claset("L", nr - 1, nr - 1, czero, czero, &v[(2 - 1)], ldv);
            }
            //            .. the right singular vectors overwrite V, the NR left singular
            //            vectors stored in U(1:NR,1:NR)
            if (wntvr || (nr == n)) {
                Cgesvd("N", "O", nr, n, v, ldv, s, u, ldu, v, ldv, &cwork[(n + 1) - 1], lcwork - n, rwork, info);
                Clapmt(false, nr, n, v, ldv, iwork);
                //               .. now [V](1:NR,1:N) contains V(1:N,1:NR)**H
            } else {
                //               .. need all N right singular vectors and NR < N
                //               [!] This is simple implementation that augments [V](1:NR,1:N)
                //               by padding a zero block. In the case NR << N, a more efficient
                //               way is to first use the LQ factorization. For more details
                //               how to implement this, see the " FULL SVD " branch.
                Claset("G", n - nr, n, czero, czero, &v[((nr + 1) - 1)], ldv);
                Cgesvd("N", "O", n, n, v, ldv, s, u, ldu, v, ldv, &cwork[(n + 1) - 1], lcwork - n, rwork, info);
                Clapmt(false, n, n, v, ldv, iwork);
            }
            //            .. now [V] contains the adjoINTEGER of the matrix of the right singular
            //            vectors of A.
        }
        //
    } else {
        //.......................................................................
        //       .. FULL SVD requested
        //.......................................................................
        if (rtrans) {
            //
            //            .. apply Cgesvd to R**H [[this option is left for R&D&T]]
            //
            if (wntvr || (nr == n)) {
                //            .. copy R**H into [V] and overwrite [V] with the left singular
                //            vectors of R**H
                for (p = 1; p <= nr; p = p + 1) {
                    for (q = p; q <= n; q = q + 1) {
                        v[(q - 1) + (p - 1) * ldv] = conj(a[(p - 1) + (q - 1) * lda]);
                    }
                }
                if (nr > 1) {
                    Claset("U", nr - 1, nr - 1, czero, czero, &v[(2 - 1) * ldv], ldv);
                }
                //
                //           .. the left singular vectors of R**H overwrite [V], the NR right
                //           singular vectors of R**H stored in [U](1:NR,1:NR) as conjugate
                //           transposed
                Cgesvd("O", "A", n, nr, v, ldv, s, v, ldv, u, ldu, &cwork[(n + 1) - 1], lcwork - n, rwork, info);
                //              .. assemble V
                for (p = 1; p <= nr; p = p + 1) {
                    v[(p - 1) + (p - 1) * ldv] = conj(v[(p - 1) + (p - 1) * ldv]);
                    for (q = p + 1; q <= nr; q = q + 1) {
                        ctmp = conj(v[(q - 1) + (p - 1) * ldv]);
                        v[(q - 1) + (p - 1) * ldv] = conj(v[(p - 1) + (q - 1) * ldv]);
                        v[(p - 1) + (q - 1) * ldv] = ctmp;
                    }
                }
                if (nr < n) {
                    for (p = 1; p <= nr; p = p + 1) {
                        for (q = nr + 1; q <= n; q = q + 1) {
                            v[(p - 1) + (q - 1) * ldv] = conj(v[(q - 1) + (p - 1) * ldv]);
                        }
                    }
                }
                Clapmt(false, nr, n, v, ldv, iwork);
                //
                for (p = 1; p <= nr; p = p + 1) {
                    u[(p - 1) + (p - 1) * ldu] = conj(u[(p - 1) + (p - 1) * ldu]);
                    for (q = p + 1; q <= nr; q = q + 1) {
                        ctmp = conj(u[(q - 1) + (p - 1) * ldu]);
                        u[(q - 1) + (p - 1) * ldu] = conj(u[(p - 1) + (q - 1) * ldu]);
                        u[(p - 1) + (q - 1) * ldu] = ctmp;
                    }
                }
                //
                if ((nr < m) && !(wntuf)) {
                    Claset("A", m - nr, nr, czero, czero, &u[((nr + 1) - 1)], ldu);
                    if (nr < n1) {
                        Claset("A", nr, n1 - nr, czero, czero, &u[((nr + 1) - 1) * ldu], ldu);
                        Claset("A", m - nr, n1 - nr, czero, cone, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                    }
                }
                //
            } else {
                //               .. need all N right singular vectors and NR < N
                //            .. copy R**H into [V] and overwrite [V] with the left singular
                //            vectors of R**H
                //               [[The optimal ratio N/NR for using QRF instead of padding
                //                 with zeros. Here hard coded to 2; it must be at least
                //                 two due to work space constraints.]]
                //               OPTRATIO = iMlaenv(6, 'Cgesvd', 'S' // 'O', NR,N,0,0)
                //               OPTRATIO = MAX( OPTRATIO, 2 )
                optratio = 2;
                if (optratio * nr > n) {
                    for (p = 1; p <= nr; p = p + 1) {
                        for (q = p; q <= n; q = q + 1) {
                            v[(q - 1) + (p - 1) * ldv] = conj(a[(p - 1) + (q - 1) * lda]);
                        }
                    }
                    if (nr > 1) {
                        Claset("U", nr - 1, nr - 1, czero, czero, &v[(2 - 1) * ldv], ldv);
                    }
                    //
                    Claset("A", n, n - nr, czero, czero, &v[((nr + 1) - 1) * ldv], ldv);
                    Cgesvd("O", "A", n, n, v, ldv, s, v, ldv, u, ldu, &cwork[(n + 1) - 1], lcwork - n, rwork, info);
                    //
                    for (p = 1; p <= n; p = p + 1) {
                        v[(p - 1) + (p - 1) * ldv] = conj(v[(p - 1) + (p - 1) * ldv]);
                        for (q = p + 1; q <= n; q = q + 1) {
                            ctmp = conj(v[(q - 1) + (p - 1) * ldv]);
                            v[(q - 1) + (p - 1) * ldv] = conj(v[(p - 1) + (q - 1) * ldv]);
                            v[(p - 1) + (q - 1) * ldv] = ctmp;
                        }
                    }
                    Clapmt(false, n, n, v, ldv, iwork);
                    //              (M x N1), i.e. (M x N) or (M x M).
                    //
                    for (p = 1; p <= n; p = p + 1) {
                        u[(p - 1) + (p - 1) * ldu] = conj(u[(p - 1) + (p - 1) * ldu]);
                        for (q = p + 1; q <= n; q = q + 1) {
                            ctmp = conj(u[(q - 1) + (p - 1) * ldu]);
                            u[(q - 1) + (p - 1) * ldu] = conj(u[(p - 1) + (q - 1) * ldu]);
                            u[(p - 1) + (q - 1) * ldu] = ctmp;
                        }
                    }
                    //
                    if ((n < m) && !(wntuf)) {
                        Claset("A", m - n, n, czero, czero, &u[((n + 1) - 1)], ldu);
                        if (n < n1) {
                            Claset("A", n, n1 - n, czero, czero, &u[((n + 1) - 1) * ldu], ldu);
                            Claset("A", m - n, n1 - n, czero, cone, &u[((n + 1) - 1) + ((n + 1) - 1) * ldu], ldu);
                        }
                    }
                } else {
                    //                  .. copy R**H into [U] and overwrite [U] with the right
                    //                  singular vectors of R
                    for (p = 1; p <= nr; p = p + 1) {
                        for (q = p; q <= n; q = q + 1) {
                            u[(q - 1) + ((nr + p) - 1) * ldu] = conj(a[(p - 1) + (q - 1) * lda]);
                        }
                    }
                    if (nr > 1) {
                        Claset("U", nr - 1, nr - 1, czero, czero, &u[((nr + 2) - 1) * ldu], ldu);
                    }
                    Cgeqrf(n, nr, &u[((nr + 1) - 1) * ldu], ldu, &cwork[(n + 1) - 1], &cwork[(n + nr + 1) - 1], lcwork - n - nr, ierr);
                    for (p = 1; p <= nr; p = p + 1) {
                        for (q = 1; q <= n; q = q + 1) {
                            v[(q - 1) + (p - 1) * ldv] = conj(u[(p - 1) + ((nr + q) - 1) * ldu]);
                        }
                    }
                    Claset("U", nr - 1, nr - 1, czero, czero, &v[(2 - 1) * ldv], ldv);
                    Cgesvd("S", "O", nr, nr, v, ldv, s, u, ldu, v, ldv, &cwork[(n + nr + 1) - 1], lcwork - n - nr, rwork, info);
                    Claset("A", n - nr, nr, czero, czero, &v[((nr + 1) - 1)], ldv);
                    Claset("A", nr, n - nr, czero, czero, &v[((nr + 1) - 1) * ldv], ldv);
                    Claset("A", n - nr, n - nr, czero, cone, &v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
                    Cunmqr("R", "C", n, n, nr, &u[((nr + 1) - 1) * ldu], ldu, &cwork[(n + 1) - 1], v, ldv, &cwork[(n + nr + 1) - 1], lcwork - n - nr, ierr);
                    Clapmt(false, n, n, v, ldv, iwork);
                    //                 (M x NR) or (M x N) or (M x M).
                    if ((nr < m) && !(wntuf)) {
                        Claset("A", m - nr, nr, czero, czero, &u[((nr + 1) - 1)], ldu);
                        if (nr < n1) {
                            Claset("A", nr, n1 - nr, czero, czero, &u[((nr + 1) - 1) * ldu], ldu);
                            Claset("A", m - nr, n1 - nr, czero, cone, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                        }
                    }
                }
            }
            //
        } else {
            //
            //            .. apply Cgesvd to R [[this is the recommended option]]
            //
            if (wntvr || (nr == n)) {
                //                .. copy R into [V] and overwrite V with the right singular vectors
                Clacpy("U", nr, n, a, lda, v, ldv);
                if (nr > 1) {
                    Claset("L", nr - 1, nr - 1, czero, czero, &v[(2 - 1)], ldv);
                }
                //               .. the right singular vectors of R overwrite [V], the NR left
                //               singular vectors of R stored in [U](1:NR,1:NR)
                Cgesvd("S", "O", nr, n, v, ldv, s, u, ldu, v, ldv, &cwork[(n + 1) - 1], lcwork - n, rwork, info);
                Clapmt(false, nr, n, v, ldv, iwork);
                //               .. now [V](1:NR,1:N) contains V(1:N,1:NR)**H
                //              (M x NR) or (M x N) or (M x M).
                if ((nr < m) && !(wntuf)) {
                    Claset("A", m - nr, nr, czero, czero, &u[((nr + 1) - 1)], ldu);
                    if (nr < n1) {
                        Claset("A", nr, n1 - nr, czero, czero, &u[((nr + 1) - 1) * ldu], ldu);
                        Claset("A", m - nr, n1 - nr, czero, cone, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
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
                //               OPTRATIO = iMlaenv(6, 'Cgesvd', 'S' // 'O', NR,N,0,0)
                //               OPTRATIO = MAX( OPTRATIO, 2 )
                optratio = 2;
                if (optratio * nr > n) {
                    Clacpy("U", nr, n, a, lda, v, ldv);
                    if (nr > 1) {
                        Claset("L", nr - 1, nr - 1, czero, czero, &v[(2 - 1)], ldv);
                    }
                    //              .. the right singular vectors of R overwrite [V], the NR left
                    //                 singular vectors of R stored in [U](1:NR,1:NR)
                    Claset("A", n - nr, n, czero, czero, &v[((nr + 1) - 1)], ldv);
                    Cgesvd("S", "O", n, n, v, ldv, s, u, ldu, v, ldv, &cwork[(n + 1) - 1], lcwork - n, rwork, info);
                    Clapmt(false, n, n, v, ldv, iwork);
                    //                 .. now [V] contains the adjoINTEGER of the matrix of the right
                    //                 singular vectors of A. The leading N left singular vectors
                    //                 are in [U](1:N,1:N)
                    //                 (M x N1), i.e. (M x N) or (M x M).
                    if ((n < m) && !(wntuf)) {
                        Claset("A", m - n, n, czero, czero, &u[((n + 1) - 1)], ldu);
                        if (n < n1) {
                            Claset("A", n, n1 - n, czero, czero, &u[((n + 1) - 1) * ldu], ldu);
                            Claset("A", m - n, n1 - n, czero, cone, &u[((n + 1) - 1) + ((n + 1) - 1) * ldu], ldu);
                        }
                    }
                } else {
                    Clacpy("U", nr, n, a, lda, &u[((nr + 1) - 1)], ldu);
                    if (nr > 1) {
                        Claset("L", nr - 1, nr - 1, czero, czero, &u[((nr + 2) - 1)], ldu);
                    }
                    Cgelqf(nr, n, &u[((nr + 1) - 1)], ldu, &cwork[(n + 1) - 1], &cwork[(n + nr + 1) - 1], lcwork - n - nr, ierr);
                    Clacpy("L", nr, nr, &u[((nr + 1) - 1)], ldu, v, ldv);
                    if (nr > 1) {
                        Claset("U", nr - 1, nr - 1, czero, czero, &v[(2 - 1) * ldv], ldv);
                    }
                    Cgesvd("S", "O", nr, nr, v, ldv, s, u, ldu, v, ldv, &cwork[(n + nr + 1) - 1], lcwork - n - nr, rwork, info);
                    Claset("A", n - nr, nr, czero, czero, &v[((nr + 1) - 1)], ldv);
                    Claset("A", nr, n - nr, czero, czero, &v[((nr + 1) - 1) * ldv], ldv);
                    Claset("A", n - nr, n - nr, czero, cone, &v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
                    Cunmlq("R", "N", n, n, nr, &u[((nr + 1) - 1)], ldu, &cwork[(n + 1) - 1], v, ldv, &cwork[(n + nr + 1) - 1], lcwork - n - nr, ierr);
                    Clapmt(false, n, n, v, ldv, iwork);
                    //              (M x NR) or (M x N) or (M x M).
                    if ((nr < m) && !(wntuf)) {
                        Claset("A", m - nr, nr, czero, czero, &u[((nr + 1) - 1)], ldu);
                        if (nr < n1) {
                            Claset("A", nr, n1 - nr, czero, czero, &u[((nr + 1) - 1) * ldu], ldu);
                            Claset("A", m - nr, n1 - nr, czero, cone, &u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                        }
                    }
                }
            }
            //        .. end of the "R**H or R" branch
        }
        //
        //           The Q matrix from the first QRF is built into the left singular
        //           vectors matrix U.
        //
        if (!wntuf) {
            Cunmqr("L", "N", m, n1, n, a, lda, cwork, u, ldu, &cwork[(n + 1) - 1], lcwork - n, ierr);
        }
        if (rowprm && !wntuf) {
            Claswp(n1, u, ldu, 1, m - 1, &iwork[(n + 1) - 1], -1);
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
    //     exact zeros in Cgesvd() applied to the (possibly truncated)
    //     full row rank triangular (trapezoidal) factor of A.
    numrank = nr;
    //
    //     End of Cgesvdq
    //
}
