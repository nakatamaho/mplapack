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

void Rsptrf(const char *uplo, INTEGER const n, REAL *ap, INTEGER *ipiv, INTEGER &info) {
    bool upper = false;
    const REAL one = 1.0;
    const REAL sevten = 17.0e+0;
    const REAL eight = 8.0e+0;
    REAL alpha = 0.0;
    INTEGER k = 0;
    INTEGER kc = 0;
    INTEGER knc = 0;
    INTEGER kstep = 0;
    REAL absakk = 0.0;
    INTEGER imax = 0;
    REAL colmax = 0.0;
    const REAL zero = 0.0;
    INTEGER kp = 0;
    REAL rowmax = 0.0;
    INTEGER jmax = 0;
    INTEGER kx = 0;
    INTEGER j = 0;
    INTEGER kpc = 0;
    INTEGER kk = 0;
    REAL t = 0.0;
    REAL r1 = 0.0;
    REAL d12 = 0.0;
    REAL d22 = 0.0;
    REAL d11 = 0.0;
    REAL wkm1 = 0.0;
    REAL wk = 0.0;
    INTEGER i = 0;
    INTEGER npp = 0;
    REAL d21 = 0.0;
    REAL wkp1 = 0.0;
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
    }
    if (info != 0) {
        Mxerbla("Rsptrf", -info);
        return;
    }
    //
    //     Initialize ALPHA for use in choosing pivot block size.
    //
    alpha = (one + sqrt(sevten)) / eight;
    //
    if (upper) {
        //
        //        Factorize A as U*D*U**T using the upper triangle of A
        //
        //        K is the main loop index, decreasing from N to 1 in steps of
        //        1 or 2
        //
        k = n;
        kc = (n - 1) * n / 2 + 1;
    statement_10:
        knc = kc;
        //
        //        If K < 1, exit from loop
        //
        if (k < 1) {
            goto statement_110;
        }
        kstep = 1;
        //
        //        Determine rows and columns to be interchanged and whether
        //        a 1-by-1 or 2-by-2 pivot block will be used
        //
        absakk = abs(ap[(kc + k - 1) - 1]);
        //
        //        IMAX is the row-index of the largest off-diagonal element in
        //        column K, and COLMAX is its absolute value
        //
        if (k > 1) {
            imax = iRamax(k - 1, &ap[kc - 1], 1);
            colmax = abs(ap[(kc + imax - 1) - 1]);
        } else {
            colmax = zero;
        }
        //
        if (max(absakk, colmax) == zero) {
            //
            //           Column K is zero: set INFO and continue
            //
            if (info == 0) {
                info = k;
            }
            kp = k;
        } else {
            if (absakk >= alpha * colmax) {
                //
                //              no interchange, use 1-by-1 pivot block
                //
                kp = k;
            } else {
                //
                rowmax = zero;
                jmax = imax;
                kx = imax * (imax + 1) / 2 + imax;
                for (j = imax + 1; j <= k; j = j + 1) {
                    if (abs(ap[kx - 1]) > rowmax) {
                        rowmax = abs(ap[kx - 1]);
                        jmax = j;
                    }
                    kx += j;
                }
                kpc = (imax - 1) * imax / 2 + 1;
                if (imax > 1) {
                    jmax = iRamax(imax - 1, &ap[kpc - 1], 1);
                    rowmax = max(rowmax, REAL(abs(ap[(kpc + jmax - 1) - 1])));
                }
                //
                if (absakk >= alpha * colmax * (colmax / rowmax)) {
                    //
                    //                 no interchange, use 1-by-1 pivot block
                    //
                    kp = k;
                } else if (abs(ap[(kpc + imax - 1) - 1]) >= alpha * rowmax) {
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
            }
            //
            kk = k - kstep + 1;
            if (kstep == 2) {
                knc = knc - k + 1;
            }
            if (kp != kk) {
                //
                //              Interchange rows and columns KK and KP in the leading
                //              submatrix A(1:k,1:k)
                //
                Rswap(kp - 1, &ap[knc - 1], 1, &ap[kpc - 1], 1);
                kx = kpc + kp - 1;
                for (j = kp + 1; j <= kk - 1; j = j + 1) {
                    kx += j - 1;
                    t = ap[(knc + j - 1) - 1];
                    ap[(knc + j - 1) - 1] = ap[kx - 1];
                    ap[kx - 1] = t;
                }
                t = ap[(knc + kk - 1) - 1];
                ap[(knc + kk - 1) - 1] = ap[(kpc + kp - 1) - 1];
                ap[(kpc + kp - 1) - 1] = t;
                if (kstep == 2) {
                    t = ap[(kc + k - 2) - 1];
                    ap[(kc + k - 2) - 1] = ap[(kc + kp - 1) - 1];
                    ap[(kc + kp - 1) - 1] = t;
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
                //              A := A - U(k)*D(k)*U(k)**T = A - W(k)*1/D(k)*W(k)**T
                //
                r1 = one / ap[(kc + k - 1) - 1];
                Rspr(uplo, k - 1, -r1, &ap[kc - 1], 1, ap);
                //
                //              Store U(k) in column k
                //
                Rscal(k - 1, r1, &ap[kc - 1], 1);
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
                //                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )**T
                //
                if (k > 2) {
                    //
                    d12 = ap[(k - 1 + (k - 1) * k / 2) - 1];
                    d22 = ap[(k - 1 + (k - 2) * (k - 1) / 2) - 1] / d12;
                    d11 = ap[(k + (k - 1) * k / 2) - 1] / d12;
                    t = one / (d11 * d22 - one);
                    d12 = t / d12;
                    //
                    for (j = k - 2; j >= 1; j = j - 1) {
                        wkm1 = d12 * (d11 * ap[(j + (k - 2) * (k - 1) / 2) - 1] - ap[(j + (k - 1) * k / 2) - 1]);
                        wk = d12 * (d22 * ap[(j + (k - 1) * k / 2) - 1] - ap[(j + (k - 2) * (k - 1) / 2) - 1]);
                        for (i = j; i >= 1; i = i - 1) {
                            ap[(i + (j - 1) * j / 2) - 1] = ap[(i + (j - 1) * j / 2) - 1] - ap[(i + (k - 1) * k / 2) - 1] * wk - ap[(i + (k - 2) * (k - 1) / 2) - 1] * wkm1;
                        }
                        ap[(j + (k - 1) * k / 2) - 1] = wk;
                        ap[(j + (k - 2) * (k - 1) / 2) - 1] = wkm1;
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
        kc = knc - k;
        goto statement_10;
        //
    } else {
        //
        //        Factorize A as L*D*L**T using the lower triangle of A
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        1 or 2
        //
        k = 1;
        kc = 1;
        npp = n * (n + 1) / 2;
    statement_60:
        knc = kc;
        //
        //        If K > N, exit from loop
        //
        if (k > n) {
            goto statement_110;
        }
        kstep = 1;
        //
        //        Determine rows and columns to be interchanged and whether
        //        a 1-by-1 or 2-by-2 pivot block will be used
        //
        absakk = abs(ap[kc - 1]);
        //
        //        IMAX is the row-index of the largest off-diagonal element in
        //        column K, and COLMAX is its absolute value
        //
        if (k < n) {
            imax = k + iRamax(n - k, &ap[(kc + 1) - 1], 1);
            colmax = abs(ap[(kc + imax - k) - 1]);
        } else {
            colmax = zero;
        }
        //
        if (max(absakk, colmax) == zero) {
            //
            //           Column K is zero: set INFO and continue
            //
            if (info == 0) {
                info = k;
            }
            kp = k;
        } else {
            if (absakk >= alpha * colmax) {
                //
                //              no interchange, use 1-by-1 pivot block
                //
                kp = k;
            } else {
                //
                //              JMAX is the column-index of the largest off-diagonal
                //              element in row IMAX, and ROWMAX is its absolute value
                //
                rowmax = zero;
                kx = kc + imax - k;
                for (j = k; j <= imax - 1; j = j + 1) {
                    if (abs(ap[kx - 1]) > rowmax) {
                        rowmax = abs(ap[kx - 1]);
                        jmax = j;
                    }
                    kx += n - j;
                }
                kpc = npp - (n - imax + 1) * (n - imax + 2) / 2 + 1;
                if (imax < n) {
                    jmax = imax + iRamax(n - imax, &ap[(kpc + 1) - 1], 1);
                    rowmax = max(rowmax, REAL(abs(ap[(kpc + jmax - imax) - 1])));
                }
                //
                if (absakk >= alpha * colmax * (colmax / rowmax)) {
                    //
                    //                 no interchange, use 1-by-1 pivot block
                    //
                    kp = k;
                } else if (abs(ap[kpc - 1]) >= alpha * rowmax) {
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
            }
            //
            kk = k + kstep - 1;
            if (kstep == 2) {
                knc += n - k + 1;
            }
            if (kp != kk) {
                //
                //              Interchange rows and columns KK and KP in the trailing
                //              submatrix A(k:n,k:n)
                //
                if (kp < n) {
                    Rswap(n - kp, &ap[(knc + kp - kk + 1) - 1], 1, &ap[(kpc + 1) - 1], 1);
                }
                kx = knc + kp - kk;
                for (j = kk + 1; j <= kp - 1; j = j + 1) {
                    kx += n - j + 1;
                    t = ap[(knc + j - kk) - 1];
                    ap[(knc + j - kk) - 1] = ap[kx - 1];
                    ap[kx - 1] = t;
                }
                t = ap[knc - 1];
                ap[knc - 1] = ap[kpc - 1];
                ap[kpc - 1] = t;
                if (kstep == 2) {
                    t = ap[(kc + 1) - 1];
                    ap[(kc + 1) - 1] = ap[(kc + kp - k) - 1];
                    ap[(kc + kp - k) - 1] = t;
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
                    //                 A := A - L(k)*D(k)*L(k)**T = A - W(k)*(1/D(k))*W(k)**T
                    //
                    r1 = one / ap[kc - 1];
                    Rspr(uplo, n - k, -r1, &ap[(kc + 1) - 1], 1, &ap[(kc + n - k + 1) - 1]);
                    //
                    //                 Store L(k) in column K
                    //
                    Rscal(n - k, r1, &ap[(kc + 1) - 1], 1);
                }
            } else {
                //
                //              2-by-2 pivot block D(k): columns K and K+1 now hold
                //
                //              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
                //
                //              where L(k) and L(k+1) are the k-th and (k+1)-th columns
                //              of L
                //
                if (k < n - 1) {
                    //
                    //                 Perform a rank-2 update of A(k+2:n,k+2:n) as
                    //
                    //                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )**T
                    //                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )**T
                    //
                    //                 where L(k) and L(k+1) are the k-th and (k+1)-th
                    //                 columns of L
                    //
                    d21 = ap[(k + 1 + (k - 1) * (2 * n - k) / 2) - 1];
                    d11 = ap[(k + 1 + k * (2 * n - k - 1) / 2) - 1] / d21;
                    d22 = ap[(k + (k - 1) * (2 * n - k) / 2) - 1] / d21;
                    t = one / (d11 * d22 - one);
                    d21 = t / d21;
                    //
                    for (j = k + 2; j <= n; j = j + 1) {
                        wk = d21 * (d11 * ap[(j + (k - 1) * (2 * n - k) / 2) - 1] - ap[(j + k * (2 * n - k - 1) / 2) - 1]);
                        wkp1 = d21 * (d22 * ap[(j + k * (2 * n - k - 1) / 2) - 1] - ap[(j + (k - 1) * (2 * n - k) / 2) - 1]);
                        //
                        for (i = j; i <= n; i = i + 1) {
                            ap[(i + (j - 1) * (2 * n - j) / 2) - 1] = ap[(i + (j - 1) * (2 * n - j) / 2) - 1] - ap[(i + (k - 1) * (2 * n - k) / 2) - 1] * wk - ap[(i + k * (2 * n - k - 1) / 2) - 1] * wkp1;
                        }
                        //
                        ap[(j + (k - 1) * (2 * n - k) / 2) - 1] = wk;
                        ap[(j + k * (2 * n - k - 1) / 2) - 1] = wkp1;
                        //
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
        kc = knc + n - k + 2;
        goto statement_60;
        //
    }
//
statement_110:;
    //
    //     End of Rsptrf
    //
}
