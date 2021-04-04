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

void Rtrevc(const char *side, const char *howmny, arr_ref<bool> select, INTEGER const &n, REAL *t, INTEGER const &ldt, REAL *vl, INTEGER const &ldvl, REAL *vr, INTEGER const &ldvr, INTEGER const &mm, INTEGER &m, REAL *work, INTEGER &info) {
    bool bothv = false;
    bool rightv = false;
    bool leftv = false;
    bool allv = false;
    bool over = false;
    bool somev = false;
    bool pair = false;
    INTEGER j = 0;
    const REAL zero = 0.0;
    REAL unfl = 0.0;
    const REAL one = 1.0;
    REAL ovfl = 0.0;
    REAL ulp = 0.0;
    REAL smlnum = 0.0;
    REAL bignum = 0.0;
    INTEGER i = 0;
    INTEGER n2 = 0;
    INTEGER ip = 0;
    INTEGER is = 0;
    INTEGER ki = 0;
    REAL wr = 0.0;
    REAL wi = 0.0;
    REAL smin = 0.0;
    INTEGER k = 0;
    INTEGER jnxt = 0;
    INTEGER j1 = 0;
    INTEGER j2 = 0;
    arr_2d<2, 2, REAL> x(fill0);
    REAL scale = 0.0;
    REAL xnorm = 0.0;
    INTEGER ierr = 0;
    REAL beta = 0.0;
    INTEGER ii = 0;
    REAL remax = 0.0;
    REAL rec = 0.0;
    REAL emax = 0.0;
    REAL vmax = 0.0;
    REAL vcrit = 0.0;
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Decode and test the input parameters
    //
    bothv = Mlsame(side, "B");
    rightv = Mlsame(side, "R") || bothv;
    leftv = Mlsame(side, "L") || bothv;
    //
    allv = Mlsame(howmny, "A");
    over = Mlsame(howmny, "B");
    somev = Mlsame(howmny, "S");
    //
    info = 0;
    if (!rightv && !leftv) {
        info = -1;
    } else if (!allv && !over && !somev) {
        info = -2;
    } else if (n < 0) {
        info = -4;
    } else if (ldt < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldvl < 1 || (leftv && ldvl < n)) {
        info = -8;
    } else if (ldvr < 1 || (rightv && ldvr < n)) {
        info = -10;
    } else {
        //
        //        Set M to the number of columns required to store the selected
        //        eigenvectors, standardize the array SELECT if necessary, and
        //        test MM.
        //
        if (somev) {
            m = 0;
            pair = false;
            for (j = 1; j <= n; j = j + 1) {
                if (pair) {
                    pair = false;
                    select[j - 1] = false;
                } else {
                    if (j < n) {
                        if (t[((j + 1) - 1) + (j - 1) * ldt] == zero) {
                            if (select[j - 1]) {
                                m++;
                            }
                        } else {
                            pair = true;
                            if (select[j - 1] || select[(j + 1) - 1]) {
                                select[j - 1] = true;
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
        if (mm < m) {
            info = -11;
        }
    }
    if (info != 0) {
        Mxerbla("Rtrevc", -info);
        return;
    }
    //
    //     Quick return if possible.
    //
    if (n == 0) {
        return;
    }
    //
    //     Set the constants to control overflow.
    //
    unfl = dlamch("Safe minimum");
    ovfl = one / unfl;
    Rlabad(unfl, ovfl);
    ulp = dlamch("Precision");
    smlnum = unfl * (n / ulp);
    bignum = (one - ulp) / smlnum;
    //
    //     Compute 1-norm of each column of strictly upper triangular
    //     part of T to control overflow in triangular solver.
    //
    work[1 - 1] = zero;
    for (j = 2; j <= n; j = j + 1) {
        work[j - 1] = zero;
        for (i = 1; i <= j - 1; i = i + 1) {
            work[j - 1] += abs(t[(i - 1) + (j - 1) * ldt]);
        }
    }
    //
    //     Index IP is used to specify the real or complex eigenvalue:
    //       IP = 0, real eigenvalue,
    //            1, first of conjugate complex pair: (wr,wi)
    //           -1, second of conjugate complex pair: (wr,wi)
    //
    n2 = 2 * n;
    //
    if (rightv) {
        //
        //        Compute right eigenvectors.
        //
        ip = 0;
        is = m;
        for (ki = n; ki >= 1; ki = ki - 1) {
            //
            if (ip == 1) {
                goto statement_130;
            }
            if (ki == 1) {
                goto statement_40;
            }
            if (t[(ki - 1) + ((ki - 1) - 1) * ldt] == zero) {
                goto statement_40;
            }
            ip = -1;
        //
        statement_40:
            if (somev) {
                if (ip == 0) {
                    if (!select[ki - 1]) {
                        goto statement_130;
                    }
                } else {
                    if (!select[(ki - 1) - 1]) {
                        goto statement_130;
                    }
                }
            }
            //
            //           Compute the KI-th eigenvalue (WR,WI).
            //
            wr = t[(ki - 1) + (ki - 1) * ldt];
            wi = zero;
            if (ip != 0) {
                wi = sqrt(abs(t[(ki - 1) + ((ki - 1) - 1) * ldt])) * sqrt(abs(t[((ki - 1) - 1) + (ki - 1) * ldt]));
            }
            smin = max(ulp * (abs(wr) + abs(wi)), smlnum);
            //
            if (ip == 0) {
                //
                //              Real right eigenvector
                //
                work[(ki + n) - 1] = one;
                //
                //              Form right-hand side
                //
                for (k = 1; k <= ki - 1; k = k + 1) {
                    work[(k + n) - 1] = -t[(k - 1) + (ki - 1) * ldt];
                }
                //
                //              Solve the upper quasi-triangular system:
                //                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
                //
                jnxt = ki - 1;
                for (j = ki - 1; j >= 1; j = j - 1) {
                    if (j > jnxt) {
                        goto statement_60;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j - 1;
                    if (j > 1) {
                        if (t[(j - 1) + ((j - 1) - 1) * ldt] != zero) {
                            j1 = j - 1;
                            jnxt = j - 2;
                        }
                    }
                    //
                    if (j1 == j2) {
                        //
                        //                    1-by-1 diagonal block
                        //
                        Rlaln2(false, 1, 1, smin, one, t[(j - 1) + (j - 1) * ldt], ldt, one, one, work[(j + n) - 1], n, wr, zero, x, 2, scale, xnorm, ierr);
                        //
                        //                    Scale X(1,1) to avoid overflow when updating
                        //                    the right-hand side.
                        //
                        if (xnorm > one) {
                            if (work[j - 1] > bignum / xnorm) {
                                x[(1 - 1)] = x[(1 - 1)] / xnorm;
                                scale = scale / xnorm;
                            }
                        }
                        //
                        //                    Scale if necessary
                        //
                        if (scale != one) {
                            Rscal(ki, scale, work[(1 + n) - 1], 1);
                        }
                        work[(j + n) - 1] = x[(1 - 1)];
                        //
                        //                    Update right-hand side
                        //
                        Raxpy(j - 1, -x[(1 - 1)], t[(j - 1) * ldt], 1, work[(1 + n) - 1], 1);
                        //
                    } else {
                        //
                        //                    2-by-2 diagonal block
                        //
                        Rlaln2(false, 2, 1, smin, one, t[((j - 1) - 1) + ((j - 1) - 1) * ldt], ldt, one, one, work[(j - 1 + n) - 1], n, wr, zero, x, 2, scale, xnorm, ierr);
                        //
                        //                    Scale X(1,1) and X(2,1) to avoid overflow when
                        //                    updating the right-hand side.
                        //
                        if (xnorm > one) {
                            beta = max(work[(j - 1) - 1], work[j - 1]);
                            if (beta > bignum / xnorm) {
                                x[(1 - 1)] = x[(1 - 1)] / xnorm;
                                x[(2 - 1)] = x[(2 - 1)] / xnorm;
                                scale = scale / xnorm;
                            }
                        }
                        //
                        //                    Scale if necessary
                        //
                        if (scale != one) {
                            Rscal(ki, scale, work[(1 + n) - 1], 1);
                        }
                        work[(j - 1 + n) - 1] = x[(1 - 1)];
                        work[(j + n) - 1] = x[(2 - 1)];
                        //
                        //                    Update right-hand side
                        //
                        Raxpy(j - 2, -x[(1 - 1)], t[((j - 1) - 1) * ldt], 1, work[(1 + n) - 1], 1);
                        Raxpy(j - 2, -x[(2 - 1)], t[(j - 1) * ldt], 1, work[(1 + n) - 1], 1);
                    }
                statement_60:;
                }
                //
                //              Copy the vector x or Q*x to VR and normalize.
                //
                if (!over) {
                    Rcopy(ki, work[(1 + n) - 1], 1, vr[(is - 1) * ldvr], 1);
                    //
                    ii = iRamax[(ki - 1) + (vr[(is - 1) * ldvr] - 1) * ldiRamax];
                    remax = one / abs(vr[(ii - 1) + (is - 1) * ldvr]);
                    Rscal(ki, remax, vr[(is - 1) * ldvr], 1);
                    //
                    for (k = ki + 1; k <= n; k = k + 1) {
                        vr[(k - 1) + (is - 1) * ldvr] = zero;
                    }
                } else {
                    if (ki > 1) {
                        Rgemv("N", n, ki - 1, one, vr, ldvr, work[(1 + n) - 1], 1, work[(ki + n) - 1], vr[(ki - 1) * ldvr], 1);
                    }
                    //
                    ii = iRamax[(n - 1) + (vr[(ki - 1) * ldvr] - 1) * ldiRamax];
                    remax = one / abs(vr[(ii - 1) + (ki - 1) * ldvr]);
                    Rscal(n, remax, vr[(ki - 1) * ldvr], 1);
                }
                //
            } else {
                //
                //              Complex right eigenvector.
                //
                //              Initial solve
                //                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
                //                [ (T(KI,KI-1)   T(KI,KI)   )               ]
                //
                if (abs(t[((ki - 1) - 1) + (ki - 1) * ldt]) >= abs(t[(ki - 1) + ((ki - 1) - 1) * ldt])) {
                    work[(ki - 1 + n) - 1] = one;
                    work[(ki + n2) - 1] = wi / t[((ki - 1) - 1) + (ki - 1) * ldt];
                } else {
                    work[(ki - 1 + n) - 1] = -wi / t[(ki - 1) + ((ki - 1) - 1) * ldt];
                    work[(ki + n2) - 1] = one;
                }
                work[(ki + n) - 1] = zero;
                work[(ki - 1 + n2) - 1] = zero;
                //
                //              Form right-hand side
                //
                for (k = 1; k <= ki - 2; k = k + 1) {
                    work[(k + n) - 1] = -work[(ki - 1 + n) - 1] * t[(k - 1) + ((ki - 1) - 1) * ldt];
                    work[(k + n2) - 1] = -work[(ki + n2) - 1] * t[(k - 1) + (ki - 1) * ldt];
                }
                //
                //              Solve upper quasi-triangular system:
                //              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
                //
                jnxt = ki - 2;
                for (j = ki - 2; j >= 1; j = j - 1) {
                    if (j > jnxt) {
                        goto statement_90;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j - 1;
                    if (j > 1) {
                        if (t[(j - 1) + ((j - 1) - 1) * ldt] != zero) {
                            j1 = j - 1;
                            jnxt = j - 2;
                        }
                    }
                    //
                    if (j1 == j2) {
                        //
                        //                    1-by-1 diagonal block
                        //
                        Rlaln2(false, 1, 2, smin, one, t[(j - 1) + (j - 1) * ldt], ldt, one, one, work[(j + n) - 1], n, wr, wi, x, 2, scale, xnorm, ierr);
                        //
                        //                    Scale X(1,1) and X(1,2) to avoid overflow when
                        //                    updating the right-hand side.
                        //
                        if (xnorm > one) {
                            if (work[j - 1] > bignum / xnorm) {
                                x[(1 - 1)] = x[(1 - 1)] / xnorm;
                                x[(2 - 1) * ldx] = x[(2 - 1) * ldx] / xnorm;
                                scale = scale / xnorm;
                            }
                        }
                        //
                        //                    Scale if necessary
                        //
                        if (scale != one) {
                            Rscal(ki, scale, work[(1 + n) - 1], 1);
                            Rscal(ki, scale, work[(1 + n2) - 1], 1);
                        }
                        work[(j + n) - 1] = x[(1 - 1)];
                        work[(j + n2) - 1] = x[(2 - 1) * ldx];
                        //
                        //                    Update the right-hand side
                        //
                        Raxpy(j - 1, -x[(1 - 1)], t[(j - 1) * ldt], 1, work[(1 + n) - 1], 1);
                        Raxpy(j - 1, -x[(2 - 1) * ldx], t[(j - 1) * ldt], 1, work[(1 + n2) - 1], 1);
                        //
                    } else {
                        //
                        //                    2-by-2 diagonal block
                        //
                        Rlaln2(false, 2, 2, smin, one, t[((j - 1) - 1) + ((j - 1) - 1) * ldt], ldt, one, one, work[(j - 1 + n) - 1], n, wr, wi, x, 2, scale, xnorm, ierr);
                        //
                        //                    Scale X to avoid overflow when updating
                        //                    the right-hand side.
                        //
                        if (xnorm > one) {
                            beta = max(work[(j - 1) - 1], work[j - 1]);
                            if (beta > bignum / xnorm) {
                                rec = one / xnorm;
                                x[(1 - 1)] = x[(1 - 1)] * rec;
                                x[(2 - 1) * ldx] = x[(2 - 1) * ldx] * rec;
                                x[(2 - 1)] = x[(2 - 1)] * rec;
                                x[(2 - 1) + (2 - 1) * ldx] = x[(2 - 1) + (2 - 1) * ldx] * rec;
                                scale = scale * rec;
                            }
                        }
                        //
                        //                    Scale if necessary
                        //
                        if (scale != one) {
                            Rscal(ki, scale, work[(1 + n) - 1], 1);
                            Rscal(ki, scale, work[(1 + n2) - 1], 1);
                        }
                        work[(j - 1 + n) - 1] = x[(1 - 1)];
                        work[(j + n) - 1] = x[(2 - 1)];
                        work[(j - 1 + n2) - 1] = x[(2 - 1) * ldx];
                        work[(j + n2) - 1] = x[(2 - 1) + (2 - 1) * ldx];
                        //
                        //                    Update the right-hand side
                        //
                        Raxpy(j - 2, -x[(1 - 1)], t[((j - 1) - 1) * ldt], 1, work[(1 + n) - 1], 1);
                        Raxpy(j - 2, -x[(2 - 1)], t[(j - 1) * ldt], 1, work[(1 + n) - 1], 1);
                        Raxpy(j - 2, -x[(2 - 1) * ldx], t[((j - 1) - 1) * ldt], 1, work[(1 + n2) - 1], 1);
                        Raxpy(j - 2, -x[(2 - 1) + (2 - 1) * ldx], t[(j - 1) * ldt], 1, work[(1 + n2) - 1], 1);
                    }
                statement_90:;
                }
                //
                //              Copy the vector x or Q*x to VR and normalize.
                //
                if (!over) {
                    Rcopy(ki, work[(1 + n) - 1], 1, vr[((is - 1) - 1) * ldvr], 1);
                    Rcopy(ki, work[(1 + n2) - 1], 1, vr[(is - 1) * ldvr], 1);
                    //
                    emax = zero;
                    for (k = 1; k <= ki; k = k + 1) {
                        emax = max(emax, abs(vr[(k - 1) + ((is - 1) - 1) * ldvr]) + abs(vr[(k - 1) + (is - 1) * ldvr]));
                    }
                    //
                    remax = one / emax;
                    Rscal(ki, remax, vr[((is - 1) - 1) * ldvr], 1);
                    Rscal(ki, remax, vr[(is - 1) * ldvr], 1);
                    //
                    for (k = ki + 1; k <= n; k = k + 1) {
                        vr[(k - 1) + ((is - 1) - 1) * ldvr] = zero;
                        vr[(k - 1) + (is - 1) * ldvr] = zero;
                    }
                    //
                } else {
                    //
                    if (ki > 2) {
                        Rgemv("N", n, ki - 2, one, vr, ldvr, work[(1 + n) - 1], 1, work[(ki - 1 + n) - 1], vr[((ki - 1) - 1) * ldvr], 1);
                        Rgemv("N", n, ki - 2, one, vr, ldvr, work[(1 + n2) - 1], 1, work[(ki + n2) - 1], vr[(ki - 1) * ldvr], 1);
                    } else {
                        Rscal(n, work[(ki - 1 + n) - 1], vr[((ki - 1) - 1) * ldvr], 1);
                        Rscal(n, work[(ki + n2) - 1], vr[(ki - 1) * ldvr], 1);
                    }
                    //
                    emax = zero;
                    for (k = 1; k <= n; k = k + 1) {
                        emax = max(emax, abs(vr[(k - 1) + ((ki - 1) - 1) * ldvr]) + abs(vr[(k - 1) + (ki - 1) * ldvr]));
                    }
                    remax = one / emax;
                    Rscal(n, remax, vr[((ki - 1) - 1) * ldvr], 1);
                    Rscal(n, remax, vr[(ki - 1) * ldvr], 1);
                }
            }
            //
            is = is - 1;
            if (ip != 0) {
                is = is - 1;
            }
        statement_130:
            if (ip == 1) {
                ip = 0;
            }
            if (ip == -1) {
                ip = 1;
            }
        }
    }
    //
    if (leftv) {
        //
        //        Compute left eigenvectors.
        //
        ip = 0;
        is = 1;
        for (ki = 1; ki <= n; ki = ki + 1) {
            //
            if (ip == -1) {
                goto statement_250;
            }
            if (ki == n) {
                goto statement_150;
            }
            if (t[((ki + 1) - 1) + (ki - 1) * ldt] == zero) {
                goto statement_150;
            }
            ip = 1;
        //
        statement_150:
            if (somev) {
                if (!select[ki - 1]) {
                    goto statement_250;
                }
            }
            //
            //           Compute the KI-th eigenvalue (WR,WI).
            //
            wr = t[(ki - 1) + (ki - 1) * ldt];
            wi = zero;
            if (ip != 0) {
                wi = sqrt(abs(t[(ki - 1) + ((ki + 1) - 1) * ldt])) * sqrt(abs(t[((ki + 1) - 1) + (ki - 1) * ldt]));
            }
            smin = max(ulp * (abs(wr) + abs(wi)), smlnum);
            //
            if (ip == 0) {
                //
                //              Real left eigenvector.
                //
                work[(ki + n) - 1] = one;
                //
                //              Form right-hand side
                //
                for (k = ki + 1; k <= n; k = k + 1) {
                    work[(k + n) - 1] = -t[(ki - 1) + (k - 1) * ldt];
                }
                //
                //              Solve the quasi-triangular system:
                //                 (T(KI+1:N,KI+1:N) - WR)**T*X = SCALE*WORK
                //
                vmax = one;
                vcrit = bignum;
                //
                jnxt = ki + 1;
                for (j = ki + 1; j <= n; j = j + 1) {
                    if (j < jnxt) {
                        goto statement_170;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j + 1;
                    if (j < n) {
                        if (t[((j + 1) - 1) + (j - 1) * ldt] != zero) {
                            j2 = j + 1;
                            jnxt = j + 2;
                        }
                    }
                    //
                    if (j1 == j2) {
                        //
                        //                    1-by-1 diagonal block
                        //
                        //                    Scale if necessary to avoid overflow when forming
                        //                    the right-hand side.
                        //
                        if (work[j - 1] > vcrit) {
                            rec = one / vmax;
                            Rscal(n - ki + 1, rec, work[(ki + n) - 1], 1);
                            vmax = one;
                            vcrit = bignum;
                        }
                        //
                        work[(j + n) - 1] = work[(j + n) - 1] - Rdot(j - ki - 1, t[((ki + 1) - 1) + (j - 1) * ldt], 1, work[(ki + 1 + n) - 1], 1);
                        //
                        //                    Solve (T(J,J)-WR)**T*X = WORK
                        //
                        Rlaln2(false, 1, 1, smin, one, t[(j - 1) + (j - 1) * ldt], ldt, one, one, work[(j + n) - 1], n, wr, zero, x, 2, scale, xnorm, ierr);
                        //
                        //                    Scale if necessary
                        //
                        if (scale != one) {
                            Rscal(n - ki + 1, scale, work[(ki + n) - 1], 1);
                        }
                        work[(j + n) - 1] = x[(1 - 1)];
                        vmax = max(abs(work[(j + n) - 1]), vmax);
                        vcrit = bignum / vmax;
                        //
                    } else {
                        //
                        //                    2-by-2 diagonal block
                        //
                        //                    Scale if necessary to avoid overflow when forming
                        //                    the right-hand side.
                        //
                        beta = max(work[j - 1], work[(j + 1) - 1]);
                        if (beta > vcrit) {
                            rec = one / vmax;
                            Rscal(n - ki + 1, rec, work[(ki + n) - 1], 1);
                            vmax = one;
                            vcrit = bignum;
                        }
                        //
                        work[(j + n) - 1] = work[(j + n) - 1] - Rdot(j - ki - 1, t[((ki + 1) - 1) + (j - 1) * ldt], 1, work[(ki + 1 + n) - 1], 1);
                        //
                        work[(j + 1 + n) - 1] = work[(j + 1 + n) - 1] - Rdot(j - ki - 1, t[((ki + 1) - 1) + ((j + 1) - 1) * ldt], 1, work[(ki + 1 + n) - 1], 1);
                        //
                        //                    Solve
                        //                      [T(J,J)-WR   T(J,J+1)     ]**T * X = SCALE*( WORK1 )
                        //                      [T(J+1,J)    T(J+1,J+1)-WR]                ( WORK2 )
                        //
                        Rlaln2(true, 2, 1, smin, one, t[(j - 1) + (j - 1) * ldt], ldt, one, one, work[(j + n) - 1], n, wr, zero, x, 2, scale, xnorm, ierr);
                        //
                        //                    Scale if necessary
                        //
                        if (scale != one) {
                            Rscal(n - ki + 1, scale, work[(ki + n) - 1], 1);
                        }
                        work[(j + n) - 1] = x[(1 - 1)];
                        work[(j + 1 + n) - 1] = x[(2 - 1)];
                        //
                        vmax = max(abs(work[(j + n) - 1]), abs(work[(j + 1 + n) - 1]), vmax);
                        vcrit = bignum / vmax;
                        //
                    }
                statement_170:;
                }
                //
                //              Copy the vector x or Q*x to VL and normalize.
                //
                if (!over) {
                    Rcopy(n - ki + 1, work[(ki + n) - 1], 1, vl[(ki - 1) + (is - 1) * ldvl], 1);
                    //
                    ii = iRamax[((n - ki + 1) - 1) + (vl[(ki - 1) + (is - 1) * ldvl] - 1) * ldiRamax] + ki - 1;
                    remax = one / abs(vl[(ii - 1) + (is - 1) * ldvl]);
                    Rscal(n - ki + 1, remax, vl[(ki - 1) + (is - 1) * ldvl], 1);
                    //
                    for (k = 1; k <= ki - 1; k = k + 1) {
                        vl[(k - 1) + (is - 1) * ldvl] = zero;
                    }
                    //
                } else {
                    //
                    if (ki < n) {
                        Rgemv("N", n, n - ki, one, vl[((ki + 1) - 1) * ldvl], ldvl, work[(ki + 1 + n) - 1], 1, work[(ki + n) - 1], vl[(ki - 1) * ldvl], 1);
                    }
                    //
                    ii = iRamax[(n - 1) + (vl[(ki - 1) * ldvl] - 1) * ldiRamax];
                    remax = one / abs(vl[(ii - 1) + (ki - 1) * ldvl]);
                    Rscal(n, remax, vl[(ki - 1) * ldvl], 1);
                    //
                }
                //
            } else {
                //
                //              Complex left eigenvector.
                //
                //               Initial solve:
                //                 ((T(KI,KI)    T(KI,KI+1) )**T - (WR - I* WI))*X = 0.
                //                 ((T(KI+1,KI) T(KI+1,KI+1))                )
                //
                if (abs(t[(ki - 1) + ((ki + 1) - 1) * ldt]) >= abs(t[((ki + 1) - 1) + (ki - 1) * ldt])) {
                    work[(ki + n) - 1] = wi / t[(ki - 1) + ((ki + 1) - 1) * ldt];
                    work[(ki + 1 + n2) - 1] = one;
                } else {
                    work[(ki + n) - 1] = one;
                    work[(ki + 1 + n2) - 1] = -wi / t[((ki + 1) - 1) + (ki - 1) * ldt];
                }
                work[(ki + 1 + n) - 1] = zero;
                work[(ki + n2) - 1] = zero;
                //
                //              Form right-hand side
                //
                for (k = ki + 2; k <= n; k = k + 1) {
                    work[(k + n) - 1] = -work[(ki + n) - 1] * t[(ki - 1) + (k - 1) * ldt];
                    work[(k + n2) - 1] = -work[(ki + 1 + n2) - 1] * t[((ki + 1) - 1) + (k - 1) * ldt];
                }
                //
                //              Solve complex quasi-triangular system:
                //              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
                //
                vmax = one;
                vcrit = bignum;
                //
                jnxt = ki + 2;
                for (j = ki + 2; j <= n; j = j + 1) {
                    if (j < jnxt) {
                        goto statement_200;
                    }
                    j1 = j;
                    j2 = j;
                    jnxt = j + 1;
                    if (j < n) {
                        if (t[((j + 1) - 1) + (j - 1) * ldt] != zero) {
                            j2 = j + 1;
                            jnxt = j + 2;
                        }
                    }
                    //
                    if (j1 == j2) {
                        //
                        //                    1-by-1 diagonal block
                        //
                        //                    Scale if necessary to avoid overflow when
                        //                    forming the right-hand side elements.
                        //
                        if (work[j - 1] > vcrit) {
                            rec = one / vmax;
                            Rscal(n - ki + 1, rec, work[(ki + n) - 1], 1);
                            Rscal(n - ki + 1, rec, work[(ki + n2) - 1], 1);
                            vmax = one;
                            vcrit = bignum;
                        }
                        //
                        work[(j + n) - 1] = work[(j + n) - 1] - Rdot(j - ki - 2, t[((ki + 2) - 1) + (j - 1) * ldt], 1, work[(ki + 2 + n) - 1], 1);
                        work[(j + n2) - 1] = work[(j + n2) - 1] - Rdot(j - ki - 2, t[((ki + 2) - 1) + (j - 1) * ldt], 1, work[(ki + 2 + n2) - 1], 1);
                        //
                        //                    Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
                        //
                        Rlaln2(false, 1, 2, smin, one, t[(j - 1) + (j - 1) * ldt], ldt, one, one, work[(j + n) - 1], n, wr, -wi, x, 2, scale, xnorm, ierr);
                        //
                        //                    Scale if necessary
                        //
                        if (scale != one) {
                            Rscal(n - ki + 1, scale, work[(ki + n) - 1], 1);
                            Rscal(n - ki + 1, scale, work[(ki + n2) - 1], 1);
                        }
                        work[(j + n) - 1] = x[(1 - 1)];
                        work[(j + n2) - 1] = x[(2 - 1) * ldx];
                        vmax = max(abs(work[(j + n) - 1]), abs(work[(j + n2) - 1]), vmax);
                        vcrit = bignum / vmax;
                        //
                    } else {
                        //
                        //                    2-by-2 diagonal block
                        //
                        //                    Scale if necessary to avoid overflow when forming
                        //                    the right-hand side elements.
                        //
                        beta = max(work[j - 1], work[(j + 1) - 1]);
                        if (beta > vcrit) {
                            rec = one / vmax;
                            Rscal(n - ki + 1, rec, work[(ki + n) - 1], 1);
                            Rscal(n - ki + 1, rec, work[(ki + n2) - 1], 1);
                            vmax = one;
                            vcrit = bignum;
                        }
                        //
                        work[(j + n) - 1] = work[(j + n) - 1] - Rdot(j - ki - 2, t[((ki + 2) - 1) + (j - 1) * ldt], 1, work[(ki + 2 + n) - 1], 1);
                        //
                        work[(j + n2) - 1] = work[(j + n2) - 1] - Rdot(j - ki - 2, t[((ki + 2) - 1) + (j - 1) * ldt], 1, work[(ki + 2 + n2) - 1], 1);
                        //
                        work[(j + 1 + n) - 1] = work[(j + 1 + n) - 1] - Rdot(j - ki - 2, t[((ki + 2) - 1) + ((j + 1) - 1) * ldt], 1, work[(ki + 2 + n) - 1], 1);
                        //
                        work[(j + 1 + n2) - 1] = work[(j + 1 + n2) - 1] - Rdot(j - ki - 2, t[((ki + 2) - 1) + ((j + 1) - 1) * ldt], 1, work[(ki + 2 + n2) - 1], 1);
                        //
                        //                    Solve 2-by-2 complex linear equation
                        //                      ([T(j,j)   T(j,j+1)  ]**T-(wr-i*wi)*I)*X = SCALE*B
                        //                      ([T(j+1,j) T(j+1,j+1)]               )
                        //
                        Rlaln2(true, 2, 2, smin, one, t[(j - 1) + (j - 1) * ldt], ldt, one, one, work[(j + n) - 1], n, wr, -wi, x, 2, scale, xnorm, ierr);
                        //
                        //                    Scale if necessary
                        //
                        if (scale != one) {
                            Rscal(n - ki + 1, scale, work[(ki + n) - 1], 1);
                            Rscal(n - ki + 1, scale, work[(ki + n2) - 1], 1);
                        }
                        work[(j + n) - 1] = x[(1 - 1)];
                        work[(j + n2) - 1] = x[(2 - 1) * ldx];
                        work[(j + 1 + n) - 1] = x[(2 - 1)];
                        work[(j + 1 + n2) - 1] = x[(2 - 1) + (2 - 1) * ldx];
                        vmax = max(abs(x[(1 - 1)]), abs(x[(2 - 1) * ldx]), abs(x[(2 - 1)]), abs(x[(2 - 1) + (2 - 1) * ldx]), vmax);
                        vcrit = bignum / vmax;
                        //
                    }
                statement_200:;
                }
                //
                //              Copy the vector x or Q*x to VL and normalize.
                //
                if (!over) {
                    Rcopy(n - ki + 1, work[(ki + n) - 1], 1, vl[(ki - 1) + (is - 1) * ldvl], 1);
                    Rcopy(n - ki + 1, work[(ki + n2) - 1], 1, vl[(ki - 1) + ((is + 1) - 1) * ldvl], 1);
                    //
                    emax = zero;
                    for (k = ki; k <= n; k = k + 1) {
                        emax = max(emax, abs(vl[(k - 1) + (is - 1) * ldvl]) + abs(vl[(k - 1) + ((is + 1) - 1) * ldvl]));
                    }
                    remax = one / emax;
                    Rscal(n - ki + 1, remax, vl[(ki - 1) + (is - 1) * ldvl], 1);
                    Rscal(n - ki + 1, remax, vl[(ki - 1) + ((is + 1) - 1) * ldvl], 1);
                    //
                    for (k = 1; k <= ki - 1; k = k + 1) {
                        vl[(k - 1) + (is - 1) * ldvl] = zero;
                        vl[(k - 1) + ((is + 1) - 1) * ldvl] = zero;
                    }
                } else {
                    if (ki < n - 1) {
                        Rgemv("N", n, n - ki - 1, one, vl[((ki + 2) - 1) * ldvl], ldvl, work[(ki + 2 + n) - 1], 1, work[(ki + n) - 1], vl[(ki - 1) * ldvl], 1);
                        Rgemv("N", n, n - ki - 1, one, vl[((ki + 2) - 1) * ldvl], ldvl, work[(ki + 2 + n2) - 1], 1, work[(ki + 1 + n2) - 1], vl[((ki + 1) - 1) * ldvl], 1);
                    } else {
                        Rscal(n, work[(ki + n) - 1], vl[(ki - 1) * ldvl], 1);
                        Rscal(n, work[(ki + 1 + n2) - 1], vl[((ki + 1) - 1) * ldvl], 1);
                    }
                    //
                    emax = zero;
                    for (k = 1; k <= n; k = k + 1) {
                        emax = max(emax, abs(vl[(k - 1) + (ki - 1) * ldvl]) + abs(vl[(k - 1) + ((ki + 1) - 1) * ldvl]));
                    }
                    remax = one / emax;
                    Rscal(n, remax, vl[(ki - 1) * ldvl], 1);
                    Rscal(n, remax, vl[((ki + 1) - 1) * ldvl], 1);
                    //
                }
                //
            }
            //
            is++;
            if (ip != 0) {
                is++;
            }
        statement_250:
            if (ip == -1) {
                ip = 0;
            }
            if (ip == 1) {
                ip = -1;
            }
            //
        }
        //
    }
    //
    //     End of Rtrevc
    //
}
