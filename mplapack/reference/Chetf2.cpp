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

inline REAL abs1(COMPLEX ff) { return max(abs(ff.real()), abs(ff.imag())); }

void Chetf2(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info) {
    COMPLEX zdum = 0.0;
    bool upper = false;
    const REAL one = 1.0;
    const REAL sevten = 17.0e+0;
    const REAL eight = 8.0e+0;
    REAL alpha = 0.0;
    INTEGER k = 0;
    INTEGER kstep = 0;
    REAL absakk = 0.0;
    INTEGER imax = 0;
    REAL colmax = 0.0;
    const REAL zero = 0.0;
    INTEGER kp = 0;
    INTEGER jmax = 0;
    REAL rowmax = 0.0;
    INTEGER kk = 0;
    INTEGER j = 0;
    COMPLEX t = 0.0;
    REAL r1 = 0.0;
    REAL d = 0.0;
    REAL d22 = 0.0;
    REAL d11 = 0.0;
    REAL tt = 0.0;
    COMPLEX d12 = 0.0;
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
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
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
        Mxerbla("Chetf2", -info);
        return;
    }
    //
    //     Initialize ALPHA for use in choosing pivot block size.
    //
    alpha = (one + sqrt(sevten)) / eight;
    //
    if (upper) {
        //
        //        Factorize A as U*D*U**H using the upper triangle of A
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
            goto statement_90;
        }
        kstep = 1;
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
            colmax = abs1(a[(imax - 1) + (k - 1) * lda]);
        } else {
            colmax = zero;
        }
        //
        if ((max(absakk, colmax) == zero) || Risnan(absakk)) {
            //
            //           Column K is zero or underflow, or contains a NaN:
            //           set INFO and continue
            //
            if (info == 0) {
                info = k;
            }
            kp = k;
            a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda].real();
        } else {
            //
            //           ============================================================
            //
            //           Test for interchange
            //
            if (absakk >= alpha * colmax) {
                //
                //              no interchange, use 1-by-1 pivot block
                //
                kp = k;
            } else {
                //
                //              JMAX is the column-index of the largest off-diagonal
                //              element in row IMAX, and ROWMAX is its absolute value.
                //              Determine only ROWMAX.
                //
                jmax = imax + iCamax(k - imax, &a[(imax - 1) + ((imax + 1) - 1) * lda], lda);
                rowmax = abs1(a[(imax - 1) + (jmax - 1) * lda] - 1);
                if (imax > 1) {
                    jmax = iCamax(imax - 1, &a[(imax - 1) * lda], 1);
                    rowmax = max(rowmax, abs1(a[(jmax - 1) + (imax - 1) * lda]));
                }
                //
                if (absakk >= alpha * colmax * (colmax / rowmax)) {
                    //
                    //                 no interchange, use 1-by-1 pivot block
                    //
                    kp = k;
                    //
                } else if (abs(a[(imax - 1) + (imax - 1) * lda].real()) >= alpha * rowmax) {
                    //
                    //                 interchange rows and columns K and IMAX, use 1-by-1
                    //                 pivot block
                    //
                    kp = imax;
                } else {
                    //
                    //                 interchange rows and columns K-1 and IMAX, use 2-by-2
                    //                 pivot block
                    //
                    kp = imax;
                    kstep = 2;
                }
                //
            }
            //
            //           ============================================================
            //
            kk = k - kstep + 1;
            if (kp != kk) {
                //
                //              Interchange rows and columns KK and KP in the leading
                //              submatrix A(1:k,1:k)
                //
                Cswap(kp - 1, &a[(kk - 1) * lda], 1, &a[(kp - 1) * lda], 1);
                for (j = kp + 1; j <= kk - 1; j = j + 1) {
                    t = conj(a[(j - 1) + (kk - 1) * lda]);
                    a[(j - 1) + (kk - 1) * lda] = conj(a[(kp - 1) + (j - 1) * lda]);
                    a[(kp - 1) + (j - 1) * lda] = t;
                }
                a[(kp - 1) + (kk - 1) * lda] = conj(a[(kp - 1) + (kk - 1) * lda]);
                r1 = a[(kk - 1) + (kk - 1) * lda].real();
                a[(kk - 1) + (kk - 1) * lda] = a[(kp - 1) + (kp - 1) * lda].real();
                a[(kp - 1) + (kp - 1) * lda] = r1;
                if (kstep == 2) {
                    a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda].real();
                    t = a[((k - 1) - 1) + (k - 1) * lda];
                    a[((k - 1) - 1) + (k - 1) * lda] = a[(kp - 1) + (k - 1) * lda];
                    a[(kp - 1) + (k - 1) * lda] = t;
                }
            } else {
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
                //              Perform a rank-1 update of A(1:k-1,1:k-1) as
                //
                //              A := A - U(k)*D(k)*U(k)**H = A - W(k)*1/D(k)*W(k)**H
                //
                r1 = one / a[(k - 1) + (k - 1) * lda].real();
                Cher(uplo, k - 1, -r1, &a[(k - 1) * lda], 1, a, lda);
                //
                //              Store U(k) in column k
                //
                CRscal(k - 1, r1, &a[(k - 1) * lda], 1);
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
                //              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**H
                //                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**H
                //
                if (k > 2) {
                    //
                    d = Rlapy2(a[((k - 1) - 1) + (k - 1) * lda].real(), a[((k - 1) - 1) + (k - 1) * lda].imag());
                    d22 = a[((k - 1) - 1) + ((k - 1) - 1) * lda].real() / d;
                    d11 = a[(k - 1) + (k - 1) * lda].real() / d;
                    tt = one / (d11 * d22 - one);
                    d12 = a[((k - 1) - 1) + (k - 1) * lda] / d;
                    d = tt / d;
                    //
                    for (j = k - 2; j >= 1; j = j - 1) {
                        wkm1 = d * (d11 * a[(j - 1) + ((k - 1) - 1) * lda] - conj(d12) * a[(j - 1) + (k - 1) * lda]);
                        wk = d * (d22 * a[(j - 1) + (k - 1) * lda] - d12 * a[(j - 1) + ((k - 1) - 1) * lda]);
                        for (i = j; i >= 1; i = i - 1) {
                            a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - a[(i - 1) + (k - 1) * lda] * conj(wk) - a[(i - 1) + ((k - 1) - 1) * lda] * conj(wkm1);
                        }
                        a[(j - 1) + (k - 1) * lda] = wk;
                        a[(j - 1) + ((k - 1) - 1) * lda] = wkm1;
                        a[(j - 1) + (j - 1) * lda] = COMPLEX(a[(j - 1) + (j - 1) * lda].real(), 0.0);
                    }
                    //
                }
                //
            }
        }
        //
        //        Store details of the interchanges in IPIV
        //
        if (kstep == 1) {
            ipiv[k - 1] = kp;
        } else {
            ipiv[k - 1] = -kp;
            ipiv[(k - 1) - 1] = -kp;
        }
        //
        //        Decrease K and return to the start of the main loop
        //
        k = k - kstep;
        goto statement_10;
        //
    } else {
        //
        //        Factorize A as L*D*L**H using the lower triangle of A
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        1 or 2
        //
        k = 1;
    statement_50:
        //
        //        If K > N, exit from loop
        //
        if (k > n) {
            goto statement_90;
        }
        kstep = 1;
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
            colmax = abs1(a[(imax - 1) + (k - 1) * lda]);
        } else {
            colmax = zero;
        }
        //
        if ((max(absakk, colmax) == zero) || Risnan(absakk)) {
            //
            //           Column K is zero or underflow, or contains a NaN:
            //           set INFO and continue
            //
            if (info == 0) {
                info = k;
            }
            kp = k;
            a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda].real();
        } else {
            //
            //           ============================================================
            //
            //           Test for interchange
            //
            if (absakk >= alpha * colmax) {
                //
                //              no interchange, use 1-by-1 pivot block
                //
                kp = k;
            } else {
                //
                //              JMAX is the column-index of the largest off-diagonal
                //              element in row IMAX, and ROWMAX is its absolute value.
                //              Determine only ROWMAX.
                //
                jmax = k - 1 + iCamax(imax - k, &a[(imax - 1) + (k - 1) * lda], lda);
                rowmax = abs1(a[(imax - 1) + (jmax - 1) * lda]);
                if (imax < n) {
                    jmax = imax + iCamax(n - imax, &a[((imax + 1) - 1) + (imax - 1) * lda], 1);
                    rowmax = max(rowmax, abs1(a[(jmax - 1) + (imax - 1) * lda]));
                }
                //
                if (absakk >= alpha * colmax * (colmax / rowmax)) {
                    //
                    //                 no interchange, use 1-by-1 pivot block
                    //
                    kp = k;
                    //
                } else if (abs(a[(imax - 1) + (imax - 1) * lda].real()) >= alpha * rowmax) {
                    //
                    //                 interchange rows and columns K and IMAX, use 1-by-1
                    //                 pivot block
                    //
                    kp = imax;
                } else {
                    //
                    //                 interchange rows and columns K+1 and IMAX, use 2-by-2
                    //                 pivot block
                    //
                    kp = imax;
                    kstep = 2;
                }
                //
            }
            //
            //           ============================================================
            //
            kk = k + kstep - 1;
            if (kp != kk) {
                //
                //              Interchange rows and columns KK and KP in the trailing
                //              submatrix A(k:n,k:n)
                //
                if (kp < n) {
                    Cswap(n - kp, &a[((kp + 1) - 1) + (kk - 1) * lda], 1, &a[((kp + 1) - 1) + (kp - 1) * lda], 1);
                }
                for (j = kk + 1; j <= kp - 1; j = j + 1) {
                    t = conj(a[(j - 1) + (kk - 1) * lda]);
                    a[(j - 1) + (kk - 1) * lda] = conj(a[(kp - 1) + (j - 1) * lda]);
                    a[(kp - 1) + (j - 1) * lda] = t;
                }
                a[(kp - 1) + (kk - 1) * lda] = conj(a[(kp - 1) + (kk - 1) * lda]);
                r1 = a[(kk - 1) + (kk - 1) * lda].real();
                a[(kk - 1) + (kk - 1) * lda] = a[(kp - 1) + (kp - 1) * lda].real();
                a[(kp - 1) + (kp - 1) * lda] = r1;
                if (kstep == 2) {
                    a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda].real();
                    t = a[((k + 1) - 1) + (k - 1) * lda];
                    a[((k + 1) - 1) + (k - 1) * lda] = a[(kp - 1) + (k - 1) * lda];
                    a[(kp - 1) + (k - 1) * lda] = t;
                }
            } else {
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
                //              1-by-1 pivot block D(k): column k now holds
                //
                //              W(k) = L(k)*D(k)
                //
                //              where L(k) is the k-th column of L
                //
                if (k < n) {
                    //
                    //                 Perform a rank-1 update of A(k+1:n,k+1:n) as
                    //
                    //                 A := A - L(k)*D(k)*L(k)**H = A - W(k)*(1/D(k))*W(k)**H
                    //
                    r1 = one / a[(k - 1) + (k - 1) * lda].real();
                    Cher(uplo, n - k, -r1, &a[((k + 1) - 1) + (k - 1) * lda], 1, &a[((k + 1) - 1) + ((k + 1) - 1) * lda], lda);
                    //
                    //                 Store L(k) in column K
                    //
                    CRscal(n - k, r1, &a[((k + 1) - 1) + (k - 1) * lda], 1);
                }
            } else {
                //
                //              2-by-2 pivot block D(k)
                //
                if (k < n - 1) {
                    //
                    //                 Perform a rank-2 update of A(k+2:n,k+2:n) as
                    //
                    //                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**H
                    //                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**H
                    //
                    //                 where L(k) and L(k+1) are the k-th and (k+1)-th
                    //                 columns of L
                    //
                    d = Rlapy2(a[((k + 1) - 1) + (k - 1) * lda].real(), a[((k + 1) - 1) + (k - 1) * lda].imag());
                    d11 = a[((k + 1) - 1) + ((k + 1) - 1) * lda].real() / d;
                    d22 = a[(k - 1) + (k - 1) * lda].real() / d;
                    tt = one / (d11 * d22 - one);
                    d21 = a[((k + 1) - 1) + (k - 1) * lda] / d;
                    d = tt / d;
                    //
                    for (j = k + 2; j <= n; j = j + 1) {
                        wk = d * (d11 * a[(j - 1) + (k - 1) * lda] - d21 * a[(j - 1) + ((k + 1) - 1) * lda]);
                        wkp1 = d * (d22 * a[(j - 1) + ((k + 1) - 1) * lda] - conj(d21) * a[(j - 1) + (k - 1) * lda]);
                        for (i = j; i <= n; i = i + 1) {
                            a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] - a[(i - 1) + (k - 1) * lda] * conj(wk) - a[(i - 1) + ((k + 1) - 1) * lda] * conj(wkp1);
                        }
                        a[(j - 1) + (k - 1) * lda] = wk;
                        a[(j - 1) + ((k + 1) - 1) * lda] = wkp1;
                        a[(j - 1) + (j - 1) * lda] = COMPLEX(a[(j - 1) + (j - 1) * lda].real(), 0.0);
                    }
                }
            }
        }
        //
        //        Store details of the interchanges in IPIV
        //
        if (kstep == 1) {
            ipiv[k - 1] = kp;
        } else {
            ipiv[k - 1] = -kp;
            ipiv[(k + 1) - 1] = -kp;
        }
        //
        //        Increase K and return to the start of the main loop
        //
        k += kstep;
        goto statement_50;
        //
    }
//
statement_90:;
    //
    //     End of Chetf2
    //
}
