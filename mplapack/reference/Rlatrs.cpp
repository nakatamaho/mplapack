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

void Rlatrs(const char *uplo, const char *trans, const char *diag, const char *normin, INTEGER const &n, REAL *a, INTEGER const &lda, REAL *x, REAL &scale, REAL *cnorm, INTEGER &info) {
    bool upper = false;
    bool notran = false;
    bool nounit = false;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    INTEGER imax = 0;
    REAL tmax = 0.0;
    REAL tscal = 0.0;
    REAL xmax = 0.0;
    REAL xbnd = 0.0;
    INTEGER jfirst = 0;
    INTEGER jlast = 0;
    INTEGER jinc = 0;
    REAL grow = 0.0;
    REAL tjj = 0.0;
    REAL xj = 0.0;
    REAL tjjs = 0.0;
    REAL rec = 0.0;
    INTEGER i = 0;
    const REAL half = 0.5e+0;
    REAL uscal = 0.0;
    REAL sumj = 0.0;
    //
    //  -- LAPACK auxiliary routine --
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
    //     .. Executable Statements ..
    //
    info = 0;
    upper = Mlsame(uplo, "U");
    notran = Mlsame(trans, "N");
    nounit = Mlsame(diag, "N");
    //
    //     Test the input parameters.
    //
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (!notran && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
        info = -2;
    } else if (!nounit && !Mlsame(diag, "U")) {
        info = -3;
    } else if (!Mlsame(normin, "Y") && !Mlsame(normin, "N")) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (lda < max((INTEGER)1, n)) {
        info = -7;
    }
    if (info != 0) {
        Mxerbla("Rlatrs", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Determine machine dependent parameters to control overflow.
    //
    smlnum = dlamch("Safe minimum") / dlamch("Precision");
    bignum = one / smlnum;
    scale = one;
    //
    if (Mlsame(normin, "N")) {
        //
        //        Compute the 1-norm of each column, not including the diagonal.
        //
        if (upper) {
            //
            //           A is upper triangular.
            //
            for (j = 1; j <= n; j = j + 1) {
                cnorm[j - 1] = Rasum[((j - 1) - 1) + (a[(j - 1) * lda] - 1) * ldRasum];
            }
        } else {
            //
            //           A is lower triangular.
            //
            for (j = 1; j <= n - 1; j = j + 1) {
                cnorm[j - 1] = Rasum[((n - j) - 1) + ((a[((j + 1) - 1) + (j - 1) * lda]) - 1) * ldRasum];
            }
            cnorm[n - 1] = zero;
        }
    }
    //
    //     Scale the column norms by TSCAL if the maximum element in CNORM is
    //     greater than BIGNUM.
    //
    imax = iRamax[(n - 1) + (cnorm - 1) * ldiRamax];
    tmax = cnorm[imax - 1];
    if (tmax <= bignum) {
        tscal = one;
    } else {
        tscal = one / (smlnum * tmax);
        Rscal(n, tscal, cnorm, 1);
    }
    //
    //     Compute a bound on the computed solution vector to see if the
    //     Level 2 BLAS routine Rtrsv can be used.
    //
    j = iRamax[(n - 1) + (x - 1) * ldiRamax];
    xmax = abs(x[j - 1]);
    xbnd = xmax;
    if (notran) {
        //
        //        Compute the growth in A * x = b.
        //
        if (upper) {
            jfirst = n;
            jlast = 1;
            jinc = -1;
        } else {
            jfirst = 1;
            jlast = n;
            jinc = 1;
        }
        //
        if (tscal != one) {
            grow = zero;
            goto statement_50;
        }
        //
        if (nounit) {
            //
            //           A is non-unit triangular.
            //
            //           Compute GROW = 1/G(j) and XBND = 1/M(j).
            //           Initially, G(0) = max{x(i), i=1,...,n}.
            //
            grow = one / max(xbnd, smlnum);
            xbnd = grow;
            for (j = jfirst; j <= jlast; j = j + jinc) {
                //
                //              Exit the loop if the growth factor is too small.
                //
                if (grow <= smlnum) {
                    goto statement_50;
                }
                //
                //              M(j) = G(j-1) / abs(A(j,j))
                //
                tjj = abs(a[(j - 1) + (j - 1) * lda]);
                xbnd = min(xbnd, min(one, tjj) * grow);
                if (tjj + cnorm[j - 1] >= smlnum) {
                    //
                    //                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
                    //
                    grow = grow * (tjj / (tjj + cnorm[j - 1]));
                } else {
                    //
                    //                 G(j) could overflow, set GROW to 0.
                    //
                    grow = zero;
                }
            }
            grow = xbnd;
        } else {
            //
            //           A is unit triangular.
            //
            //           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
            //
            grow = min(one, one / max(xbnd, smlnum));
            for (j = jfirst; j <= jlast; j = j + jinc) {
                //
                //              Exit the loop if the growth factor is too small.
                //
                if (grow <= smlnum) {
                    goto statement_50;
                }
                //
                //              G(j) = G(j-1)*( 1 + CNORM(j) )
                //
                grow = grow * (one / (one + cnorm[j - 1]));
            }
        }
    statement_50:;
        //
    } else {
        //
        //        Compute the growth in A**T * x = b.
        //
        if (upper) {
            jfirst = 1;
            jlast = n;
            jinc = 1;
        } else {
            jfirst = n;
            jlast = 1;
            jinc = -1;
        }
        //
        if (tscal != one) {
            grow = zero;
            goto statement_80;
        }
        //
        if (nounit) {
            //
            //           A is non-unit triangular.
            //
            //           Compute GROW = 1/G(j) and XBND = 1/M(j).
            //           Initially, M(0) = max{x(i), i=1,...,n}.
            //
            grow = one / max(xbnd, smlnum);
            xbnd = grow;
            for (j = jfirst; j <= jlast; j = j + jinc) {
                //
                //              Exit the loop if the growth factor is too small.
                //
                if (grow <= smlnum) {
                    goto statement_80;
                }
                //
                //              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
                //
                xj = one + cnorm[j - 1];
                grow = min(grow, xbnd / xj);
                //
                //              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
                //
                tjj = abs(a[(j - 1) + (j - 1) * lda]);
                if (xj > tjj) {
                    xbnd = xbnd * (tjj / xj);
                }
            }
            grow = min(grow, xbnd);
        } else {
            //
            //           A is unit triangular.
            //
            //           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
            //
            grow = min(one, one / max(xbnd, smlnum));
            for (j = jfirst; j <= jlast; j = j + jinc) {
                //
                //              Exit the loop if the growth factor is too small.
                //
                if (grow <= smlnum) {
                    goto statement_80;
                }
                //
                //              G(j) = ( 1 + CNORM(j) )*G(j-1)
                //
                xj = one + cnorm[j - 1];
                grow = grow / xj;
            }
        }
    statement_80:;
    }
    //
    if ((grow * tscal) > smlnum) {
        //
        //        Use the Level 2 BLAS solve if the reciprocal of the bound on
        //        elements of X is not too small.
        //
        Rtrsv(uplo, trans, diag, n, a, lda, x, 1);
    } else {
        //
        //        Use a Level 1 BLAS solve, scaling INTEGERermediate results.
        //
        if (xmax > bignum) {
            //
            //           Scale X so that its components are less than or equal to
            //           BIGNUM in absolute value.
            //
            scale = bignum / xmax;
            Rscal(n, scale, x, 1);
            xmax = bignum;
        }
        //
        if (notran) {
            //
            //           Solve A * x = b
            //
            for (j = jfirst; j <= jlast; j = j + jinc) {
                //
                //              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
                //
                xj = abs(x[j - 1]);
                if (nounit) {
                    tjjs = a[(j - 1) + (j - 1) * lda] * tscal;
                } else {
                    tjjs = tscal;
                    if (tscal == one) {
                        goto statement_100;
                    }
                }
                tjj = abs(tjjs);
                if (tjj > smlnum) {
                    //
                    //                    abs(A(j,j)) > SMLNUM:
                    //
                    if (tjj < one) {
                        if (xj > tjj * bignum) {
                            //
                            //                          Scale x by 1/b(j).
                            //
                            rec = one / xj;
                            Rscal(n, rec, x, 1);
                            scale = scale * rec;
                            xmax = xmax * rec;
                        }
                    }
                    x[j - 1] = x[j - 1] / tjjs;
                    xj = abs(x[j - 1]);
                } else if (tjj > zero) {
                    //
                    //                    0 < abs(A(j,j)) <= SMLNUM:
                    //
                    if (xj > tjj * bignum) {
                        //
                        //                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                        //                       to avoid overflow when dividing by A(j,j).
                        //
                        rec = (tjj * bignum) / xj;
                        if (cnorm[j - 1] > one) {
                            //
                            //                          Scale by 1/CNORM(j) to avoid overflow when
                            //                          multiplying x(j) times column j.
                            //
                            rec = rec / cnorm[j - 1];
                        }
                        Rscal(n, rec, x, 1);
                        scale = scale * rec;
                        xmax = xmax * rec;
                    }
                    x[j - 1] = x[j - 1] / tjjs;
                    xj = abs(x[j - 1]);
                } else {
                    //
                    //                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                    //                    scale = 0, and compute a solution to A*x = 0.
                    //
                    for (i = 1; i <= n; i = i + 1) {
                        x[i - 1] = zero;
                    }
                    x[j - 1] = one;
                    xj = one;
                    scale = zero;
                    xmax = zero;
                }
            statement_100:
                //
                //              Scale x if necessary to avoid overflow when adding a
                //              multiple of column j of A.
                //
                if (xj > one) {
                    rec = one / xj;
                    if (cnorm[j - 1] > (bignum - xmax) * rec) {
                        //
                        //                    Scale x by 1/(2*abs(x(j))).
                        //
                        rec = rec * half;
                        Rscal(n, rec, x, 1);
                        scale = scale * rec;
                    }
                } else if (xj * cnorm[j - 1] > (bignum - xmax)) {
                    //
                    //                 Scale x by 1/2.
                    //
                    Rscal(n, half, x, 1);
                    scale = scale * half;
                }
                //
                if (upper) {
                    if (j > 1) {
                        //
                        //                    Compute the update
                        //                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
                        //
                        Raxpy(j - 1, -x[j - 1] * tscal, a[(j - 1) * lda], 1, x, 1);
                        i = iRamax[((j - 1) - 1) + (x - 1) * ldiRamax];
                        xmax = abs(x[i - 1]);
                    }
                } else {
                    if (j < n) {
                        //
                        //                    Compute the update
                        //                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
                        //
                        Raxpy(n - j, -x[j - 1] * tscal, a[((j + 1) - 1) + (j - 1) * lda], 1, x[(j + 1) - 1], 1);
                        i = j + iRamax[((n - j) - 1) + ((x[(j + 1) - 1]) - 1) * ldiRamax];
                        xmax = abs(x[i - 1]);
                    }
                }
            }
            //
        } else {
            //
            //           Solve A**T * x = b
            //
            for (j = jfirst; j <= jlast; j = j + jinc) {
                //
                //              Compute x(j) = b(j) - sum A(k,j)*x(k).
                //                                    k<>j
                //
                xj = abs(x[j - 1]);
                uscal = tscal;
                rec = one / max(xmax, one);
                if (cnorm[j - 1] > (bignum - xj) * rec) {
                    //
                    //                 If x(j) could overflow, scale x by 1/(2*XMAX).
                    //
                    rec = rec * half;
                    if (nounit) {
                        tjjs = a[(j - 1) + (j - 1) * lda] * tscal;
                    } else {
                        tjjs = tscal;
                    }
                    tjj = abs(tjjs);
                    if (tjj > one) {
                        //
                        //                       Divide by A(j,j) when scaling x if A(j,j) > 1.
                        //
                        rec = min(one, rec * tjj);
                        uscal = uscal / tjjs;
                    }
                    if (rec < one) {
                        Rscal(n, rec, x, 1);
                        scale = scale * rec;
                        xmax = xmax * rec;
                    }
                }
                //
                sumj = zero;
                if (uscal == one) {
                    //
                    //                 If the scaling needed for A in the dot product is 1,
                    //                 call Rdot to perform the dot product.
                    //
                    if (upper) {
                        sumj = Rdot(j - 1, a[(j - 1) * lda], 1, x, 1);
                    } else if (j < n) {
                        sumj = Rdot(n - j, a[((j + 1) - 1) + (j - 1) * lda], 1, x[(j + 1) - 1], 1);
                    }
                } else {
                    //
                    //                 Otherwise, use in-line code for the dot product.
                    //
                    if (upper) {
                        for (i = 1; i <= j - 1; i = i + 1) {
                            sumj += (a[(i - 1) + (j - 1) * lda] * uscal) * x[i - 1];
                        }
                    } else if (j < n) {
                        for (i = j + 1; i <= n; i = i + 1) {
                            sumj += (a[(i - 1) + (j - 1) * lda] * uscal) * x[i - 1];
                        }
                    }
                }
                //
                if (uscal == tscal) {
                    //
                    //                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
                    //                 was not used to scale the dotproduct.
                    //
                    x[j - 1] = x[j - 1] - sumj;
                    xj = abs(x[j - 1]);
                    if (nounit) {
                        tjjs = a[(j - 1) + (j - 1) * lda] * tscal;
                    } else {
                        tjjs = tscal;
                        if (tscal == one) {
                            goto statement_150;
                        }
                    }
                    //
                    //                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
                    //
                    tjj = abs(tjjs);
                    if (tjj > smlnum) {
                        //
                        //                       abs(A(j,j)) > SMLNUM:
                        //
                        if (tjj < one) {
                            if (xj > tjj * bignum) {
                                //
                                //                             Scale X by 1/abs(x(j)).
                                //
                                rec = one / xj;
                                Rscal(n, rec, x, 1);
                                scale = scale * rec;
                                xmax = xmax * rec;
                            }
                        }
                        x[j - 1] = x[j - 1] / tjjs;
                    } else if (tjj > zero) {
                        //
                        //                       0 < abs(A(j,j)) <= SMLNUM:
                        //
                        if (xj > tjj * bignum) {
                            //
                            //                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
                            //
                            rec = (tjj * bignum) / xj;
                            Rscal(n, rec, x, 1);
                            scale = scale * rec;
                            xmax = xmax * rec;
                        }
                        x[j - 1] = x[j - 1] / tjjs;
                    } else {
                        //
                        //                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        //                       scale = 0, and compute a solution to A**T*x = 0.
                        //
                        for (i = 1; i <= n; i = i + 1) {
                            x[i - 1] = zero;
                        }
                        x[j - 1] = one;
                        scale = zero;
                        xmax = zero;
                    }
                statement_150:;
                } else {
                    //
                    //                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot
                    //                 product has already been divided by 1/A(j,j).
                    //
                    x[j - 1] = x[j - 1] / tjjs - sumj;
                }
                xmax = max(xmax, abs(x[j - 1]));
            }
        }
        scale = scale / tscal;
    }
    //
    //     Scale the column norms by 1/TSCAL for return.
    //
    if (tscal != one) {
        Rscal(n, one / tscal, cnorm, 1);
    }
    //
    //     End of Rlatrs
    //
}
