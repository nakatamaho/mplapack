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

void Chetf2_rk(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *e, INTEGER *ipiv, INTEGER &info) {
    COMPLEX z = 0.0;
    bool upper = false;
    const REAL one = 1.0;
    const REAL sevten = 17.0e+0;
    const REAL eight = 8.0e+0;
    REAL alpha = 0.0;
    REAL sfmin = 0.0;
    const COMPLEX czero = (0.0, 0.0);
    INTEGER k = 0;
    INTEGER kstep = 0;
    INTEGER p = 0;
    REAL absakk = 0.0;
    INTEGER imax = 0;
    REAL colmax = 0.0;
    const REAL zero = 0.0;
    INTEGER kp = 0;
    bool done = false;
    INTEGER jmax = 0;
    REAL rowmax = 0.0;
    INTEGER itemp = 0;
    REAL dtemp = 0.0;
    INTEGER kk = 0;
    INTEGER j = 0;
    COMPLEX t = 0.0;
    REAL r1 = 0.0;
    REAL d11 = 0.0;
    INTEGER ii = 0;
    REAL d = 0.0;
    REAL d22 = 0.0;
    COMPLEX d12 = 0.0;
    REAL tt = 0.0;
    COMPLEX wkm1 = 0.0;
    COMPLEX wk = 0.0;
    INTEGER i = 0;
    COMPLEX d21 = 0.0;
    COMPLEX wkp1 = 0.0;
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
    //  ======================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    abs1[z - 1] = abs(z.real()) + abs(z.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Chetf2_rk", -info);
        return;
    }
    //
    //     Initialize ALPHA for use in choosing pivot block size.
    //
    alpha = (one + sqrt(sevten)) / eight;
    //
    //     Compute machine safe minimum
    //
    sfmin = Rlamch("S");
    //
    if (upper) {
        //
        //        Factorize A as U*D*U**H using the upper triangle of A
        //
        //        Initialize the first entry of array E, where superdiagonal
        //        elements of D are stored
        //
        e[1 - 1] = czero;
        //
        //        K is the main loop index, decreasing from N to 1 in steps of
        //        1 or 2
        //
        k = n;
    statement_10:
        //
        //        If K < 1, exit from loop
        //
        if (k < 1) {
            goto statement_34;
        }
        kstep = 1;
        p = k;
        //
        //        Determine rows and columns to be interchanged and whether
        //        a 1-by-1 or 2-by-2 pivot block will be used
        //
        absakk = abs(a[(k - 1) + (k - 1) * lda].real());
        //
        //        IMAX is the row-index of the largest off-diagonal element in
        //        column K, and COLMAX is its absolute value.
        //        Determine both COLMAX and IMAX.
        //
        if (k > 1) {
            imax = iCamax(k - 1, &a[(k - 1) * lda], 1);
            colmax = abs1[a[(imax - 1) + (k - 1) * lda] - 1];
        } else {
            colmax = zero;
        }
        //
        if ((max(absakk, colmax) == zero)) {
            //
            //           Column K is zero or underflow: set INFO and continue
            //
            if (info == 0) {
                info = k;
            }
            kp = k;
            a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda].real();
            //
            //           Set E( K ) to zero
            //
            if (k > 1) {
                e[k - 1] = czero;
            }
            //
        } else {
            //
            //           ============================================================
            //
            //           BEGIN pivot search
            //
            //           Case(1)
            //           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX
            //           (used to handle NaN and Inf)
            //
            if (!(absakk < alpha * colmax)) {
                //
                //              no interchange, use 1-by-1 pivot block
                //
                kp = k;
                //
            } else {
                //
                done = false;
            //
            //              Loop until pivot found
            //
            statement_12:
                //
                //                 BEGIN pivot search loop body
                //
                //                 JMAX is the column-index of the largest off-diagonal
                //                 element in row IMAX, and ROWMAX is its absolute value.
                //                 Determine both ROWMAX and JMAX.
                //
                if (imax != k) {
                    jmax = imax + iCamax(k - imax, &a[(imax - 1) + ((imax + 1) - 1) * lda], lda);
                    rowmax = abs1[a[(imax - 1) + (jmax - 1) * lda] - 1];
                } else {
                    rowmax = zero;
                }
                //
                if (imax > 1) {
                    itemp = iCamax(imax - 1, &a[(imax - 1) * lda], 1);
                    dtemp = abs1[a[(itemp - 1) + (imax - 1) * lda] - 1];
                    if (dtemp > rowmax) {
                        rowmax = dtemp;
                        jmax = itemp;
                    }
                }
                //
                //                 Case(2)
                //                 Equivalent to testing for
                //                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX
                //                 (used to handle NaN and Inf)
                //
                if (!(abs(a[(imax - 1) + (imax - 1) * lda].real()) < alpha * rowmax)) {
                    //
                    //                    interchange rows and columns K and IMAX,
                    //                    use 1-by-1 pivot block
                    //
                    kp = imax;
                    done = true;
                    //
                    //                 Case(3)
                    //                 Equivalent to testing for ROWMAX.EQ.COLMAX,
                    //                 (used to handle NaN and Inf)
                    //
                } else if ((p == jmax) || (rowmax <= colmax)) {
                    //
                    //                    interchange rows and columns K-1 and IMAX,
                    //                    use 2-by-2 pivot block
                    //
                    kp = imax;
                    kstep = 2;
                    done = true;
                    //
                    //                 Case(4)
                } else {
                    //
                    //                    Pivot not found: set params and repeat
                    //
                    p = imax;
                    colmax = rowmax;
                    imax = jmax;
                }
                //
                //                 END pivot search loop body
                //
                if (!done) {
                    goto statement_12;
                }
                //
            }
            //
            //           END pivot search
            //
            //           ============================================================
            //
            //           KK is the column of A where pivoting step stopped
            //
            kk = k - kstep + 1;
            //
            //           For only a 2x2 pivot, interchange rows and columns K and P
            //           in the leading submatrix A(1:k,1:k)
            //
            if ((kstep == 2) && (p != k)) {
                //              (1) Swap columnar parts
                if (p > 1) {
                    Cswap(p - 1, &a[(k - 1) * lda], 1, &a[(p - 1) * lda], 1);
                }
                //              (2) Swap and conjugate middle parts
                for (j = p + 1; j <= k - 1; j = j + 1) {
                    t = conj(a[(j - 1) + (k - 1) * lda]);
                    a[(j - 1) + (k - 1) * lda] = conj(a[(p - 1) + (j - 1) * lda]);
                    a[(p - 1) + (j - 1) * lda] = t;
                }
                //              (3) Swap and conjugate corner elements at row-col interserction
                a[(p - 1) + (k - 1) * lda] = conj(a[(p - 1) + (k - 1) * lda]);
                //              (4) Swap diagonal elements at row-col intersection
                r1 = a[(k - 1) + (k - 1) * lda].real();
                a[(k - 1) + (k - 1) * lda] = a[(p - 1) + (p - 1) * lda].real();
                a[(p - 1) + (p - 1) * lda] = r1;
                //
                //              Convert upper triangle of A into U form by applying
                //              the interchanges in columns k+1:N.
                //
                if (k < n) {
                    Cswap(n - k, &a[(k - 1) + ((k + 1) - 1) * lda], lda, &a[(p - 1) + ((k + 1) - 1) * lda], lda);
                }
                //
            }
            //
            //           For both 1x1 and 2x2 pivots, interchange rows and
            //           columns KK and KP in the leading submatrix A(1:k,1:k)
            //
            if (kp != kk) {
                //              (1) Swap columnar parts
                if (kp > 1) {
                    Cswap(kp - 1, &a[(kk - 1) * lda], 1, &a[(kp - 1) * lda], 1);
                }
                //              (2) Swap and conjugate middle parts
                for (j = kp + 1; j <= kk - 1; j = j + 1) {
                    t = conj(a[(j - 1) + (kk - 1) * lda]);
                    a[(j - 1) + (kk - 1) * lda] = conj(a[(kp - 1) + (j - 1) * lda]);
                    a[(kp - 1) + (j - 1) * lda] = t;
                }
                //              (3) Swap and conjugate corner elements at row-col interserction
                a[(kp - 1) + (kk - 1) * lda] = conj(a[(kp - 1) + (kk - 1) * lda]);
                //              (4) Swap diagonal elements at row-col intersection
                r1 = a[(kk - 1) + (kk - 1) * lda].real();
                a[(kk - 1) + (kk - 1) * lda] = a[(kp - 1) + (kp - 1) * lda].real();
                a[(kp - 1) + (kp - 1) * lda] = r1;
                //
                if (kstep == 2) {
                    //                 (*) Make sure that diagonal element of pivot is real
                    a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda].real();
                    //                 (5) Swap row elements
                    t = a[((k - 1) - 1) + (k - 1) * lda];
                    a[((k - 1) - 1) + (k - 1) * lda] = a[(kp - 1) + (k - 1) * lda];
                    a[(kp - 1) + (k - 1) * lda] = t;
                }
                //
                //              Convert upper triangle of A into U form by applying
                //              the interchanges in columns k+1:N.
                //
                if (k < n) {
                    Cswap(n - k, &a[(kk - 1) + ((k + 1) - 1) * lda], lda, &a[(kp - 1) + ((k + 1) - 1) * lda], lda);
                }
                //
            } else {
                //              (*) Make sure that diagonal element of pivot is real
                a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda].real();
                if (kstep == 2) {
                    a[((k - 1) - 1) + ((k - 1) - 1) * lda] = a[((k - 1) - 1) + ((k - 1) - 1) * lda].real();
                }
            }
            //
            //           Update the leading submatrix
            //
            if (kstep == 1) {
                //
                //              1-by-1 pivot block D(k): column k now holds
                //
                //              W(k) = U(k)*D(k)
                //
                //              where U(k) is the k-th column of U
                //
                if (k > 1) {
                    //
                    //                 Perform a rank-1 update of A(1:k-1,1:k-1) and
                    //                 store U(k) in column k
                    //
                    if (abs(a[(k - 1) + (k - 1) * lda].real()) >= sfmin) {
                        //
                        //                    Perform a rank-1 update of A(1:k-1,1:k-1) as
                        //                    A := A - U(k)*D(k)*U(k)**T
                        //                       = A - W(k)*1/D(k)*W(k)**T
                        //
                        d11 = one / a[(k - 1) + (k - 1) * lda].real();
                        Cher(uplo, k - 1, -d11, &a[(k - 1) * lda], 1, a, lda);
                        //
                        //                    Store U(k) in column k
                        //
                        CRscal(k - 1, d11, &a[(k - 1) * lda], 1);
                    } else {
                        //
                        //                    Store L(k) in column K
                        //
                        d11 = a[(k - 1) + (k - 1) * lda].real();
                        for (ii = 1; ii <= k - 1; ii = ii + 1) {
                            a[(ii - 1) + (k - 1) * lda] = a[(ii - 1) + (k - 1) * lda] / d11;
                        }
                        //
                        //                    Perform a rank-1 update of A(k+1:n,k+1:n) as
                        //                    A := A - U(k)*D(k)*U(k)**T
                        //                       = A - W(k)*(1/D(k))*W(k)**T
                        //                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
                        //
                        Cher(uplo, k - 1, -d11, &a[(k - 1) * lda], 1, a, lda);
                    }
                    //
                    //                 Store the superdiagonal element of D in array E
                    //
                    e[k - 1] = czero;
                    //
                }
                //
            } else {
                //
                //              2-by-2 pivot block D(k): columns k and k-1 now hold
                //
                //              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
                //
                //              where U(k) and U(k-1) are the k-th and (k-1)-th columns
                //              of U
                //
                //              Perform a rank-2 update of A(1:k-2,1:k-2) as
                //
                //              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
                //                 = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T
                //
                //              and store L(k) and L(k+1) in columns k and k+1
                //
                if (k > 2) {
                    //                 D = |A12|
                    d = Rlapy2(a[((k - 1) - 1) + (k - 1) * lda].real(), &a[((k - 1) - 1) + (k - 1) * lda].imag());
                    d11 = a[(k - 1) + (k - 1) * lda] / d;
                    d22 = a[((k - 1) - 1) + ((k - 1) - 1) * lda] / d;
                    d12 = a[((k - 1) - 1) + (k - 1) * lda] / d;
                    tt = one / (d11 * d22 - one);
                    //
                    for (j = k - 2; j >= 1; j = j - 1) {
                        //
                        //                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J
                        //
                        wkm1 = tt * (d11 * a[(j - 1) + ((k - 1) - 1) * lda] - conj(d12) * a[(j - 1) + (k - 1) * lda]);
                        wk = tt * (d22 * a[(j - 1) + (k - 1) * lda] - d12 * a[(j - 1) + ((k - 1) - 1) * lda]);
                        //
                        //                    Perform a rank-2 update of A(1:k-2,1:k-2)
                        //
                        for (i = j; i >= 1; i = i - 1) {
                            a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - (a[(i - 1) + (k - 1) * lda] / d) * conj(wk) - (a[(i - 1) + ((k - 1) - 1) * lda] / d) * conj(wkm1);
                        }
                        //
                        //                    Store U(k) and U(k-1) in cols k and k-1 for row J
                        //
                        a[(j - 1) + (k - 1) * lda] = wk / d;
                        a[(j - 1) + ((k - 1) - 1) * lda] = wkm1 / d;
                        //                    (*) Make sure that diagonal element of pivot is real
                        a[(j - 1) + (j - 1) * lda] = COMPLEX(a[(j - 1) + (j - 1) * lda].real(), zero);
                        //
                    }
                    //
                }
                //
                //              Copy superdiagonal elements of D(K) to E(K) and
                //              ZERO out superdiagonal entry of A
                //
                e[k - 1] = a[((k - 1) - 1) + (k - 1) * lda];
                e[(k - 1) - 1] = czero;
                a[((k - 1) - 1) + (k - 1) * lda] = czero;
                //
            }
            //
            //           End column K is nonsingular
            //
        }
        //
        //        Store details of the interchanges in IPIV
        //
        if (kstep == 1) {
            ipiv[k - 1] = kp;
        } else {
            ipiv[k - 1] = -p;
            ipiv[(k - 1) - 1] = -kp;
        }
        //
        //        Decrease K and return to the start of the main loop
        //
        k = k - kstep;
        goto statement_10;
    //
    statement_34:;
        //
    } else {
        //
        //        Factorize A as L*D*L**H using the lower triangle of A
        //
        //        Initialize the unused last entry of the subdiagonal array E.
        //
        e[n - 1] = czero;
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        1 or 2
        //
        k = 1;
    statement_40:
        //
        //        If K > N, exit from loop
        //
        if (k > n) {
            goto statement_64;
        }
        kstep = 1;
        p = k;
        //
        //        Determine rows and columns to be interchanged and whether
        //        a 1-by-1 or 2-by-2 pivot block will be used
        //
        absakk = abs(a[(k - 1) + (k - 1) * lda].real());
        //
        //        IMAX is the row-index of the largest off-diagonal element in
        //        column K, and COLMAX is its absolute value.
        //        Determine both COLMAX and IMAX.
        //
        if (k < n) {
            imax = k + iCamax(n - k, &a[((k + 1) - 1) + (k - 1) * lda], 1);
            colmax = abs1[a[(imax - 1) + (k - 1) * lda] - 1];
        } else {
            colmax = zero;
        }
        //
        if (max(absakk, colmax) == zero) {
            //
            //           Column K is zero or underflow: set INFO and continue
            //
            if (info == 0) {
                info = k;
            }
            kp = k;
            a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda].real();
            //
            //           Set E( K ) to zero
            //
            if (k < n) {
                e[k - 1] = czero;
            }
            //
        } else {
            //
            //           ============================================================
            //
            //           BEGIN pivot search
            //
            //           Case(1)
            //           Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX
            //           (used to handle NaN and Inf)
            //
            if (!(absakk < alpha * colmax)) {
                //
                //              no interchange, use 1-by-1 pivot block
                //
                kp = k;
                //
            } else {
                //
                done = false;
            //
            //              Loop until pivot found
            //
            statement_42:
                //
                //                 BEGIN pivot search loop body
                //
                //                 JMAX is the column-index of the largest off-diagonal
                //                 element in row IMAX, and ROWMAX is its absolute value.
                //                 Determine both ROWMAX and JMAX.
                //
                if (imax != k) {
                    jmax = k - 1 + iCamax(imax - k, &a[(imax - 1) + (k - 1) * lda], lda);
                    rowmax = abs1[a[(imax - 1) + (jmax - 1) * lda] - 1];
                } else {
                    rowmax = zero;
                }
                //
                if (imax < n) {
                    itemp = imax + iCamax(n - imax, &a[((imax + 1) - 1) + (imax - 1) * lda], 1);
                    dtemp = abs1[a[(itemp - 1) + (imax - 1) * lda] - 1];
                    if (dtemp > rowmax) {
                        rowmax = dtemp;
                        jmax = itemp;
                    }
                }
                //
                //                 Case(2)
                //                 Equivalent to testing for
                //                 ABS( REAL( W( IMAX,KW-1 ) ) ).GE.ALPHA*ROWMAX
                //                 (used to handle NaN and Inf)
                //
                if (!(abs(a[(imax - 1) + (imax - 1) * lda].real()) < alpha * rowmax)) {
                    //
                    //                    interchange rows and columns K and IMAX,
                    //                    use 1-by-1 pivot block
                    //
                    kp = imax;
                    done = true;
                    //
                    //                 Case(3)
                    //                 Equivalent to testing for ROWMAX.EQ.COLMAX,
                    //                 (used to handle NaN and Inf)
                    //
                } else if ((p == jmax) || (rowmax <= colmax)) {
                    //
                    //                    interchange rows and columns K+1 and IMAX,
                    //                    use 2-by-2 pivot block
                    //
                    kp = imax;
                    kstep = 2;
                    done = true;
                    //
                    //                 Case(4)
                } else {
                    //
                    //                    Pivot not found: set params and repeat
                    //
                    p = imax;
                    colmax = rowmax;
                    imax = jmax;
                }
                //
                //                 END pivot search loop body
                //
                if (!done) {
                    goto statement_42;
                }
                //
            }
            //
            //           END pivot search
            //
            //           ============================================================
            //
            //           KK is the column of A where pivoting step stopped
            //
            kk = k + kstep - 1;
            //
            //           For only a 2x2 pivot, interchange rows and columns K and P
            //           in the trailing submatrix A(k:n,k:n)
            //
            if ((kstep == 2) && (p != k)) {
                //              (1) Swap columnar parts
                if (p < n) {
                    Cswap(n - p, &a[((p + 1) - 1) + (k - 1) * lda], 1, &a[((p + 1) - 1) + (p - 1) * lda], 1);
                }
                //              (2) Swap and conjugate middle parts
                for (j = k + 1; j <= p - 1; j = j + 1) {
                    t = conj(a[(j - 1) + (k - 1) * lda]);
                    a[(j - 1) + (k - 1) * lda] = conj(a[(p - 1) + (j - 1) * lda]);
                    a[(p - 1) + (j - 1) * lda] = t;
                }
                //              (3) Swap and conjugate corner elements at row-col interserction
                a[(p - 1) + (k - 1) * lda] = conj(a[(p - 1) + (k - 1) * lda]);
                //              (4) Swap diagonal elements at row-col intersection
                r1 = a[(k - 1) + (k - 1) * lda].real();
                a[(k - 1) + (k - 1) * lda] = a[(p - 1) + (p - 1) * lda].real();
                a[(p - 1) + (p - 1) * lda] = r1;
                //
                //              Convert lower triangle of A into L form by applying
                //              the interchanges in columns 1:k-1.
                //
                if (k > 1) {
                    Cswap(k - 1, &a[(k - 1)], lda, &a[(p - 1)], lda);
                }
                //
            }
            //
            //           For both 1x1 and 2x2 pivots, interchange rows and
            //           columns KK and KP in the trailing submatrix A(k:n,k:n)
            //
            if (kp != kk) {
                //              (1) Swap columnar parts
                if (kp < n) {
                    Cswap(n - kp, &a[((kp + 1) - 1) + (kk - 1) * lda], 1, &a[((kp + 1) - 1) + (kp - 1) * lda], 1);
                }
                //              (2) Swap and conjugate middle parts
                for (j = kk + 1; j <= kp - 1; j = j + 1) {
                    t = conj(a[(j - 1) + (kk - 1) * lda]);
                    a[(j - 1) + (kk - 1) * lda] = conj(a[(kp - 1) + (j - 1) * lda]);
                    a[(kp - 1) + (j - 1) * lda] = t;
                }
                //              (3) Swap and conjugate corner elements at row-col interserction
                a[(kp - 1) + (kk - 1) * lda] = conj(a[(kp - 1) + (kk - 1) * lda]);
                //              (4) Swap diagonal elements at row-col intersection
                r1 = a[(kk - 1) + (kk - 1) * lda].real();
                a[(kk - 1) + (kk - 1) * lda] = a[(kp - 1) + (kp - 1) * lda].real();
                a[(kp - 1) + (kp - 1) * lda] = r1;
                //
                if (kstep == 2) {
                    //                 (*) Make sure that diagonal element of pivot is real
                    a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda].real();
                    //                 (5) Swap row elements
                    t = a[((k + 1) - 1) + (k - 1) * lda];
                    a[((k + 1) - 1) + (k - 1) * lda] = a[(kp - 1) + (k - 1) * lda];
                    a[(kp - 1) + (k - 1) * lda] = t;
                }
                //
                //              Convert lower triangle of A into L form by applying
                //              the interchanges in columns 1:k-1.
                //
                if (k > 1) {
                    Cswap(k - 1, &a[(kk - 1)], lda, &a[(kp - 1)], lda);
                }
                //
            } else {
                //              (*) Make sure that diagonal element of pivot is real
                a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda].real();
                if (kstep == 2) {
                    a[((k + 1) - 1) + ((k + 1) - 1) * lda] = a[((k + 1) - 1) + ((k + 1) - 1) * lda].real();
                }
            }
            //
            //           Update the trailing submatrix
            //
            if (kstep == 1) {
                //
                //              1-by-1 pivot block D(k): column k of A now holds
                //
                //              W(k) = L(k)*D(k),
                //
                //              where L(k) is the k-th column of L
                //
                if (k < n) {
                    //
                    //                 Perform a rank-1 update of A(k+1:n,k+1:n) and
                    //                 store L(k) in column k
                    //
                    //                 Handle division by a small number
                    //
                    if (abs(a[(k - 1) + (k - 1) * lda].real()) >= sfmin) {
                        //
                        //                    Perform a rank-1 update of A(k+1:n,k+1:n) as
                        //                    A := A - L(k)*D(k)*L(k)**T
                        //                       = A - W(k)*(1/D(k))*W(k)**T
                        //
                        d11 = one / a[(k - 1) + (k - 1) * lda].real();
                        Cher(uplo, n - k, -d11, &a[((k + 1) - 1) + (k - 1) * lda], 1, &a[((k + 1) - 1) + ((k + 1) - 1) * lda], lda);
                        //
                        //                    Store L(k) in column k
                        //
                        CRscal(n - k, d11, &a[((k + 1) - 1) + (k - 1) * lda], 1);
                    } else {
                        //
                        //                    Store L(k) in column k
                        //
                        d11 = a[(k - 1) + (k - 1) * lda].real();
                        for (ii = k + 1; ii <= n; ii = ii + 1) {
                            a[(ii - 1) + (k - 1) * lda] = a[(ii - 1) + (k - 1) * lda] / d11;
                        }
                        //
                        //                    Perform a rank-1 update of A(k+1:n,k+1:n) as
                        //                    A := A - L(k)*D(k)*L(k)**T
                        //                       = A - W(k)*(1/D(k))*W(k)**T
                        //                       = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
                        //
                        Cher(uplo, n - k, -d11, &a[((k + 1) - 1) + (k - 1) * lda], 1, &a[((k + 1) - 1) + ((k + 1) - 1) * lda], lda);
                    }
                    //
                    //                 Store the subdiagonal element of D in array E
                    //
                    e[k - 1] = czero;
                    //
                }
                //
            } else {
                //
                //              2-by-2 pivot block D(k): columns k and k+1 now hold
                //
                //              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
                //
                //              where L(k) and L(k+1) are the k-th and (k+1)-th columns
                //              of L
                //
                //              Perform a rank-2 update of A(k+2:n,k+2:n) as
                //
                //              A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T
                //                 = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T
                //
                //              and store L(k) and L(k+1) in columns k and k+1
                //
                if (k < n - 1) {
                    //                 D = |A21|
                    d = Rlapy2(a[((k + 1) - 1) + (k - 1) * lda].real(), &a[((k + 1) - 1) + (k - 1) * lda].imag());
                    d11 = a[((k + 1) - 1) + ((k + 1) - 1) * lda].real() / d;
                    d22 = a[(k - 1) + (k - 1) * lda].real() / d;
                    d21 = a[((k + 1) - 1) + (k - 1) * lda] / d;
                    tt = one / (d11 * d22 - one);
                    //
                    for (j = k + 2; j <= n; j = j + 1) {
                        //
                        //                    Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J
                        //
                        wk = tt * (d11 * a[(j - 1) + (k - 1) * lda] - d21 * a[(j - 1) + ((k + 1) - 1) * lda]);
                        wkp1 = tt * (d22 * a[(j - 1) + ((k + 1) - 1) * lda] - conj(d21) * a[(j - 1) + (k - 1) * lda]);
                        //
                        //                    Perform a rank-2 update of A(k+2:n,k+2:n)
                        //
                        for (i = j; i <= n; i = i + 1) {
                            a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - (a[(i - 1) + (k - 1) * lda] / d) * conj(wk) - (a[(i - 1) + ((k + 1) - 1) * lda] / d) * conj(wkp1);
                        }
                        //
                        //                    Store L(k) and L(k+1) in cols k and k+1 for row J
                        //
                        a[(j - 1) + (k - 1) * lda] = wk / d;
                        a[(j - 1) + ((k + 1) - 1) * lda] = wkp1 / d;
                        //                    (*) Make sure that diagonal element of pivot is real
                        a[(j - 1) + (j - 1) * lda] = COMPLEX(a[(j - 1) + (j - 1) * lda].real(), zero);
                        //
                    }
                    //
                }
                //
                //              Copy subdiagonal elements of D(K) to E(K) and
                //              ZERO out subdiagonal entry of A
                //
                e[k - 1] = a[((k + 1) - 1) + (k - 1) * lda];
                e[(k + 1) - 1] = czero;
                a[((k + 1) - 1) + (k - 1) * lda] = czero;
                //
            }
            //
            //           End column K is nonsingular
            //
        }
        //
        //        Store details of the interchanges in IPIV
        //
        if (kstep == 1) {
            ipiv[k - 1] = kp;
        } else {
            ipiv[k - 1] = -p;
            ipiv[(k + 1) - 1] = -kp;
        }
        //
        //        Increase K and return to the start of the main loop
        //
        k += kstep;
        goto statement_40;
    //
    statement_64:;
        //
    }
    //
    //     End of Chetf2_rk
    //
}
