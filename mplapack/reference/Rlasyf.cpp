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

void Rlasyf(const char *uplo, INTEGER const n, INTEGER const nb, INTEGER &kb, REAL *a, INTEGER const lda, INTEGER *ipiv, REAL *w, INTEGER const ldw, INTEGER &info) {
    const REAL one = 1.0;
    const REAL sevten = 17.0e+0;
    const REAL eight = 8.0e+0;
    REAL alpha = 0.0;
    INTEGER k = 0;
    INTEGER kw = 0;
    INTEGER kstep = 0;
    REAL absakk = 0.0;
    INTEGER imax = 0;
    REAL colmax = 0.0;
    const REAL zero = 0.0;
    INTEGER kp = 0;
    INTEGER jmax = 0;
    REAL rowmax = 0.0;
    INTEGER kk = 0;
    INTEGER kkw = 0;
    REAL r1 = 0.0;
    REAL d21 = 0.0;
    REAL d11 = 0.0;
    REAL d22 = 0.0;
    REAL t = 0.0;
    INTEGER j = 0;
    INTEGER jb = 0;
    INTEGER jj = 0;
    INTEGER jp = 0;
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
    info = 0;
    //
    //     Initialize ALPHA for use in choosing pivot block size.
    //
    alpha = (one + sqrt(sevten)) / eight;
    //
    if (Mlsame(uplo, "U")) {
        //
        //        Factorize the trailing columns of A using the upper triangle
        //        of A and working backwards, and compute the matrix W = U12*D
        //        for use in updating A11
        //
        //        K is the main loop index, decreasing from N in steps of 1 or 2
        //
        //        KW is the column of W which corresponds to column K of A
        //
        k = n;
    statement_10:
        kw = nb + k - n;
        //
        //        Exit from loop
        //
        if ((k <= n - nb + 1 && nb < n) || k < 1) {
            goto statement_30;
        }
        //
        //        Copy column K of A to column KW of W and update it
        //
        Rcopy(k, &a[(k - 1) * lda], 1, &w[(kw - 1) * ldw], 1);
        if (k < n) {
            Rgemv("No transpose", k, n - k, -one, &a[((k + 1) - 1) * lda], lda, &w[(k - 1) + ((kw + 1) - 1) * ldw], ldw, one, &w[(kw - 1) * ldw], 1);
        }
        //
        kstep = 1;
        //
        //        Determine rows and columns to be INTEGERerchanged and whether
        //        a 1-by-1 or 2-by-2 pivot block will be used
        //
        absakk = abs(w[(k - 1) + (kw - 1) * ldw]);
        //
        //        IMAX is the row-index of the largest off-diagonal element in
        //        column K, and COLMAX is its absolute value.
        //        Determine both COLMAX and IMAX.
        //
        if (k > 1) {
            imax = iRamax(k - 1, &w[(kw - 1) * ldw], 1);
            colmax = abs(w[(imax - 1) + (kw - 1) * ldw]);
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
        } else {
            if (absakk >= alpha * colmax) {
                //
                //              no INTEGERerchange, use 1-by-1 pivot block
                //
                kp = k;
            } else {
                //
                //              Copy column IMAX to column KW-1 of W and update it
                //
                Rcopy(imax, &a[(imax - 1) * lda], 1, &w[((kw - 1) - 1) * ldw], 1);
                Rcopy(k - imax, &a[(imax - 1) + ((imax + 1) - 1) * lda], lda, &w[((imax + 1) - 1) + ((kw - 1) - 1) * ldw], 1);
                if (k < n) {
                    Rgemv("No transpose", k, n - k, -one, &a[((k + 1) - 1) * lda], lda, &w[(imax - 1) + ((kw + 1) - 1) * ldw], ldw, one, &w[((kw - 1) - 1) * ldw], 1);
                }
                //
                //              JMAX is the column-index of the largest off-diagonal
                //              element in row IMAX, and ROWMAX is its absolute value
                //
                jmax = imax + iRamax(k - imax, &w[((imax + 1) - 1) + ((kw - 1) - 1) * ldw], 1);
                rowmax = abs(w[(jmax - 1) + ((kw - 1) - 1) * ldw]);
                if (imax > 1) {
                    jmax = iRamax(imax - 1, &w[((kw - 1) - 1) * ldw], 1);
                    rowmax = max(rowmax, abs(w[(jmax - 1) + ((kw - 1) - 1) * ldw]));
                }
                //
                if (absakk >= alpha * colmax * (colmax / rowmax)) {
                    //
                    //                 no INTEGERerchange, use 1-by-1 pivot block
                    //
                    kp = k;
                } else if (abs(w[(imax - 1) + ((kw - 1) - 1) * ldw]) >= alpha * rowmax) {
                    //
                    //                 INTEGERerchange rows and columns K and IMAX, use 1-by-1
                    //                 pivot block
                    //
                    kp = imax;
                    //
                    //                 copy column KW-1 of W to column KW of W
                    //
                    Rcopy(k, &w[((kw - 1) - 1) * ldw], 1, &w[(kw - 1) * ldw], 1);
                } else {
                    //
                    //                 INTEGERerchange rows and columns K-1 and IMAX, use 2-by-2
                    //                 pivot block
                    //
                    kp = imax;
                    kstep = 2;
                }
            }
            //
            //           ============================================================
            //
            //           KK is the column of A where pivoting step stopped
            //
            kk = k - kstep + 1;
            //
            //           KKW is the column of W which corresponds to column KK of A
            //
            kkw = nb + kk - n;
            //
            //           Interchange rows and columns KP and KK.
            //           Updated column KP is already stored in column KKW of W.
            //
            if (kp != kk) {
                //
                //              Copy non-updated column KK to column KP of submatrix A
                //              at step K. No need to copy element INTEGERo column K
                //              (or K and K-1 for 2-by-2 pivot) of A, since these columns
                //              will be later overwritten.
                //
                a[(kp - 1) + (kp - 1) * lda] = a[(kk - 1) + (kk - 1) * lda];
                Rcopy(kk - 1 - kp, &a[((kp + 1) - 1) + (kk - 1) * lda], 1, &a[(kp - 1) + ((kp + 1) - 1) * lda], lda);
                if (kp > 1) {
                    Rcopy(kp - 1, &a[(kk - 1) * lda], 1, &a[(kp - 1) * lda], 1);
                }
                //
                //              Interchange rows KK and KP in last K+1 to N columns of A
                //              (columns K (or K and K-1 for 2-by-2 pivot) of A will be
                //              later overwritten). Interchange rows KK and KP
                //              in last KKW to NB columns of W.
                //
                if (k < n) {
                    Rswap(n - k, &a[(kk - 1) + ((k + 1) - 1) * lda], lda, &a[(kp - 1) + ((k + 1) - 1) * lda], lda);
                }
                Rswap(n - kk + 1, &w[(kk - 1) + (kkw - 1) * ldw], ldw, &w[(kp - 1) + (kkw - 1) * ldw], ldw);
            }
            //
            if (kstep == 1) {
                //
                //              1-by-1 pivot block D(k): column kw of W now holds
                //
                //              W(kw) = U(k)*D(k),
                //
                //              where U(k) is the k-th column of U
                //
                //              Store subdiag. elements of column U(k)
                //              and 1-by-1 block D(k) in column k of A.
                //              NOTE: Diagonal element U(k,k) is a UNIT element
                //              and not stored.
                //                 A(k,k) := D(k,k) = W(k,kw)
                //                 A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k)
                //
                Rcopy(k, &w[(kw - 1) * ldw], 1, &a[(k - 1) * lda], 1);
                r1 = one / a[(k - 1) + (k - 1) * lda];
                Rscal(k - 1, r1, &a[(k - 1) * lda], 1);
                //
            } else {
                //
                //              2-by-2 pivot block D(k): columns kw and kw-1 of W now hold
                //
                //              ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k)
                //
                //              where U(k) and U(k-1) are the k-th and (k-1)-th columns
                //              of U
                //
                //              Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2
                //              block D(k-1:k,k-1:k) in columns k-1 and k of A.
                //              NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT
                //              block and not stored.
                //                 A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw)
                //                 A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) =
                //                 = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) )
                //
                if (k > 2) {
                    //
                    //                 Compose the columns of the inverse of 2-by-2 pivot
                    //                 block D in the following way to reduce the number
                    //                 of FLOPS when we myltiply panel ( W(kw-1) W(kw) ) by
                    //                 this inverse
                    //
                    //                 D**(-1) = ( d11 d21 )**(-1) =
                    //                           ( d21 d22 )
                    //
                    //                 = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                    //                                        ( (-d21 ) ( d11 ) )
                    //
                    //                 = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                    //
                    //                   * ( ( d22/d21 ) (      -1 ) ) =
                    //                     ( (      -1 ) ( d11/d21 ) )
                    //
                    //                 = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) =
                    //                                           ( ( -1  ) ( D22 ) )
                    //
                    //                 = 1/d21 * T * ( ( D11 ) (  -1 ) )
                    //                               ( (  -1 ) ( D22 ) )
                    //
                    //                 = D21 * ( ( D11 ) (  -1 ) )
                    //                         ( (  -1 ) ( D22 ) )
                    //
                    d21 = w[((k - 1) - 1) + (kw - 1) * ldw];
                    d11 = w[(k - 1) + (kw - 1) * ldw] / d21;
                    d22 = w[((k - 1) - 1) + ((kw - 1) - 1) * ldw] / d21;
                    t = one / (d11 * d22 - one);
                    d21 = t / d21;
                    //
                    //                 Update elements in columns A(k-1) and A(k) as
                    //                 dot products of rows of ( W(kw-1) W(kw) ) and columns
                    //                 of D**(-1)
                    //
                    for (j = 1; j <= k - 2; j = j + 1) {
                        a[(j - 1) + ((k - 1) - 1) * lda] = d21 * (d11 * w[(j - 1) + ((kw - 1) - 1) * ldw] - w[(j - 1) + (kw - 1) * ldw]);
                        a[(j - 1) + (k - 1) * lda] = d21 * (d22 * w[(j - 1) + (kw - 1) * ldw] - w[(j - 1) + ((kw - 1) - 1) * ldw]);
                    }
                }
                //
                //              Copy D(k) to A
                //
                a[((k - 1) - 1) + ((k - 1) - 1) * lda] = w[((k - 1) - 1) + ((kw - 1) - 1) * ldw];
                a[((k - 1) - 1) + (k - 1) * lda] = w[((k - 1) - 1) + (kw - 1) * ldw];
                a[(k - 1) + (k - 1) * lda] = w[(k - 1) + (kw - 1) * ldw];
                //
            }
            //
        }
        //
        //        Store details of the INTEGERerchanges in IPIV
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
    statement_30:
        //
        //        Update the upper triangle of A11 (= A(1:k,1:k)) as
        //
        //        A11 := A11 - U12*D*U12**T = A11 - U12*W**T
        //
        //        computing blocks of NB columns at a time
        //
        for (j = ((k - 1) / nb) * nb + 1; j <= 1; j = j + -nb) {
            jb = min(nb, k - j + 1);
            //
            //           Update the upper triangle of the diagonal block
            //
            for (jj = j; jj <= j + jb - 1; jj = jj + 1) {
                Rgemv("No transpose", jj - j + 1, n - k, -one, &a[(j - 1) + ((k + 1) - 1) * lda], lda, &w[(jj - 1) + ((kw + 1) - 1) * ldw], ldw, one, &a[(j - 1) + (jj - 1) * lda], 1);
            }
            //
            //           Update the rectangular superdiagonal block
            //
            Rgemm("No transpose", "Transpose", j - 1, jb, n - k, -one, &a[((k + 1) - 1) * lda], lda, &w[(j - 1) + ((kw + 1) - 1) * ldw], ldw, one, &a[(j - 1) * lda], lda);
        }
        //
        //        Put U12 in standard form by partially undoing the INTEGERerchanges
        //        in columns k+1:n looping backwards from k+1 to n
        //
        j = k + 1;
    statement_60:
        //
        //           Undo the INTEGERerchanges (if any) of rows JJ and JP at each
        //           step J
        //
        //           (Here, J is a diagonal index)
        jj = j;
        jp = ipiv[j - 1];
        if (jp < 0) {
            jp = -jp;
            //              (Here, J is a diagonal index)
            j++;
        }
        //           (NOTE: Here, J is used to determine row length. Length N-J+1
        //           of the rows to swap back doesn't include diagonal element)
        j++;
        if (jp != jj && j <= n) {
            Rswap(n - j + 1, &a[(jp - 1) + (j - 1) * lda], lda, &a[(jj - 1) + (j - 1) * lda], lda);
        }
        if (j < n) {
            goto statement_60;
        }
        //
        //        Set KB to the number of columns factorized
        //
        kb = n - k;
        //
    } else {
        //
        //        Factorize the leading columns of A using the lower triangle
        //        of A and working forwards, and compute the matrix W = L21*D
        //        for use in updating A22
        //
        //        K is the main loop index, increasing from 1 in steps of 1 or 2
        //
        k = 1;
    statement_70:
        //
        //        Exit from loop
        //
        if ((k >= nb && nb < n) || k > n) {
            goto statement_90;
        }
        //
        //        Copy column K of A to column K of W and update it
        //
        Rcopy(n - k + 1, &a[(k - 1) + (k - 1) * lda], 1, &w[(k - 1) + (k - 1) * ldw], 1);
        Rgemv("No transpose", n - k + 1, k - 1, -one, &a[(k - 1)], lda, &w[(k - 1)], ldw, one, &w[(k - 1) + (k - 1) * ldw], 1);
        //
        kstep = 1;
        //
        //        Determine rows and columns to be INTEGERerchanged and whether
        //        a 1-by-1 or 2-by-2 pivot block will be used
        //
        absakk = abs(w[(k - 1) + (k - 1) * ldw]);
        //
        //        IMAX is the row-index of the largest off-diagonal element in
        //        column K, and COLMAX is its absolute value.
        //        Determine both COLMAX and IMAX.
        //
        if (k < n) {
            imax = k + iRamax(n - k, &w[((k + 1) - 1) + (k - 1) * ldw], 1);
            colmax = abs(w[(imax - 1) + (k - 1) * ldw]);
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
        } else {
            if (absakk >= alpha * colmax) {
                //
                //              no INTEGERerchange, use 1-by-1 pivot block
                //
                kp = k;
            } else {
                //
                //              Copy column IMAX to column K+1 of W and update it
                //
                Rcopy(imax - k, &a[(imax - 1) + (k - 1) * lda], lda, &w[(k - 1) + ((k + 1) - 1) * ldw], 1);
                Rcopy(n - imax + 1, &a[(imax - 1) + (imax - 1) * lda], 1, &w[(imax - 1) + ((k + 1) - 1) * ldw], 1);
                Rgemv("No transpose", n - k + 1, k - 1, -one, &a[(k - 1)], lda, &w[(imax - 1)], ldw, one, &w[(k - 1) + ((k + 1) - 1) * ldw], 1);
                //
                //              JMAX is the column-index of the largest off-diagonal
                //              element in row IMAX, and ROWMAX is its absolute value
                //
                jmax = k - 1 + iRamax(imax - k, &w[(k - 1) + ((k + 1) - 1) * ldw], 1);
                rowmax = abs(w[(jmax - 1) + ((k + 1) - 1) * ldw]);
                if (imax < n) {
                    jmax = imax + iRamax(n - imax, &w[((imax + 1) - 1) + ((k + 1) - 1) * ldw], 1);
                    rowmax = max(rowmax, abs(w[(jmax - 1) + ((k + 1) - 1) * ldw]));
                }
                //
                if (absakk >= alpha * colmax * (colmax / rowmax)) {
                    //
                    //                 no INTEGERerchange, use 1-by-1 pivot block
                    //
                    kp = k;
                } else if (abs(w[(imax - 1) + ((k + 1) - 1) * ldw]) >= alpha * rowmax) {
                    //
                    //                 INTEGERerchange rows and columns K and IMAX, use 1-by-1
                    //                 pivot block
                    //
                    kp = imax;
                    //
                    //                 copy column K+1 of W to column K of W
                    //
                    Rcopy(n - k + 1, &w[(k - 1) + ((k + 1) - 1) * ldw], 1, &w[(k - 1) + (k - 1) * ldw], 1);
                } else {
                    //
                    //                 INTEGERerchange rows and columns K+1 and IMAX, use 2-by-2
                    //                 pivot block
                    //
                    kp = imax;
                    kstep = 2;
                }
            }
            //
            //           ============================================================
            //
            //           KK is the column of A where pivoting step stopped
            //
            kk = k + kstep - 1;
            //
            //           Interchange rows and columns KP and KK.
            //           Updated column KP is already stored in column KK of W.
            //
            if (kp != kk) {
                //
                //              Copy non-updated column KK to column KP of submatrix A
                //              at step K. No need to copy element INTEGERo column K
                //              (or K and K+1 for 2-by-2 pivot) of A, since these columns
                //              will be later overwritten.
                //
                a[(kp - 1) + (kp - 1) * lda] = a[(kk - 1) + (kk - 1) * lda];
                Rcopy(kp - kk - 1, &a[((kk + 1) - 1) + (kk - 1) * lda], 1, &a[(kp - 1) + ((kk + 1) - 1) * lda], lda);
                if (kp < n) {
                    Rcopy(n - kp, &a[((kp + 1) - 1) + (kk - 1) * lda], 1, &a[((kp + 1) - 1) + (kp - 1) * lda], 1);
                }
                //
                //              Interchange rows KK and KP in first K-1 columns of A
                //              (columns K (or K and K+1 for 2-by-2 pivot) of A will be
                //              later overwritten). Interchange rows KK and KP
                //              in first KK columns of W.
                //
                if (k > 1) {
                    Rswap(k - 1, &a[(kk - 1)], lda, &a[(kp - 1)], lda);
                }
                Rswap(kk, &w[(kk - 1)], ldw, &w[(kp - 1)], ldw);
            }
            //
            if (kstep == 1) {
                //
                //              1-by-1 pivot block D(k): column k of W now holds
                //
                //              W(k) = L(k)*D(k),
                //
                //              where L(k) is the k-th column of L
                //
                //              Store subdiag. elements of column L(k)
                //              and 1-by-1 block D(k) in column k of A.
                //              (NOTE: Diagonal element L(k,k) is a UNIT element
                //              and not stored)
                //                 A(k,k) := D(k,k) = W(k,k)
                //                 A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k)
                //
                Rcopy(n - k + 1, &w[(k - 1) + (k - 1) * ldw], 1, &a[(k - 1) + (k - 1) * lda], 1);
                if (k < n) {
                    r1 = one / a[(k - 1) + (k - 1) * lda];
                    Rscal(n - k, r1, &a[((k + 1) - 1) + (k - 1) * lda], 1);
                }
                //
            } else {
                //
                //              2-by-2 pivot block D(k): columns k and k+1 of W now hold
                //
                //              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
                //
                //              where L(k) and L(k+1) are the k-th and (k+1)-th columns
                //              of L
                //
                //              Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2
                //              block D(k:k+1,k:k+1) in columns k and k+1 of A.
                //              (NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT
                //              block and not stored)
                //                 A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1)
                //                 A(k+2:N,k:k+1) := L(k+2:N,k:k+1) =
                //                 = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) )
                //
                if (k < n - 1) {
                    //
                    //                 Compose the columns of the inverse of 2-by-2 pivot
                    //                 block D in the following way to reduce the number
                    //                 of FLOPS when we myltiply panel ( W(k) W(k+1) ) by
                    //                 this inverse
                    //
                    //                 D**(-1) = ( d11 d21 )**(-1) =
                    //                           ( d21 d22 )
                    //
                    //                 = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                    //                                        ( (-d21 ) ( d11 ) )
                    //
                    //                 = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                    //
                    //                   * ( ( d22/d21 ) (      -1 ) ) =
                    //                     ( (      -1 ) ( d11/d21 ) )
                    //
                    //                 = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) =
                    //                                           ( ( -1  ) ( D22 ) )
                    //
                    //                 = 1/d21 * T * ( ( D11 ) (  -1 ) )
                    //                               ( (  -1 ) ( D22 ) )
                    //
                    //                 = D21 * ( ( D11 ) (  -1 ) )
                    //                         ( (  -1 ) ( D22 ) )
                    //
                    d21 = w[((k + 1) - 1) + (k - 1) * ldw];
                    d11 = w[((k + 1) - 1) + ((k + 1) - 1) * ldw] / d21;
                    d22 = w[(k - 1) + (k - 1) * ldw] / d21;
                    t = one / (d11 * d22 - one);
                    d21 = t / d21;
                    //
                    //                 Update elements in columns A(k) and A(k+1) as
                    //                 dot products of rows of ( W(k) W(k+1) ) and columns
                    //                 of D**(-1)
                    //
                    for (j = k + 2; j <= n; j = j + 1) {
                        a[(j - 1) + (k - 1) * lda] = d21 * (d11 * w[(j - 1) + (k - 1) * ldw] - w[(j - 1) + ((k + 1) - 1) * ldw]);
                        a[(j - 1) + ((k + 1) - 1) * lda] = d21 * (d22 * w[(j - 1) + ((k + 1) - 1) * ldw] - w[(j - 1) + (k - 1) * ldw]);
                    }
                }
                //
                //              Copy D(k) to A
                //
                a[(k - 1) + (k - 1) * lda] = w[(k - 1) + (k - 1) * ldw];
                a[((k + 1) - 1) + (k - 1) * lda] = w[((k + 1) - 1) + (k - 1) * ldw];
                a[((k + 1) - 1) + ((k + 1) - 1) * lda] = w[((k + 1) - 1) + ((k + 1) - 1) * ldw];
                //
            }
            //
        }
        //
        //        Store details of the INTEGERerchanges in IPIV
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
        goto statement_70;
    //
    statement_90:
        //
        //        Update the lower triangle of A22 (= A(k:n,k:n)) as
        //
        //        A22 := A22 - L21*D*L21**T = A22 - L21*W**T
        //
        //        computing blocks of NB columns at a time
        //
        for (j = k; j <= n; j = j + nb) {
            jb = min(nb, n - j + 1);
            //
            //           Update the lower triangle of the diagonal block
            //
            for (jj = j; jj <= j + jb - 1; jj = jj + 1) {
                Rgemv("No transpose", j + jb - jj, k - 1, -one, &a[(jj - 1)], lda, &w[(jj - 1)], ldw, one, &a[(jj - 1) + (jj - 1) * lda], 1);
            }
            //
            //           Update the rectangular subdiagonal block
            //
            if (j + jb <= n) {
                Rgemm("No transpose", "Transpose", n - j - jb + 1, jb, k - 1, -one, &a[((j + jb) - 1)], lda, &w[(j - 1)], ldw, one, &a[((j + jb) - 1) + (j - 1) * lda], lda);
            }
        }
        //
        //        Put L21 in standard form by partially undoing the INTEGERerchanges
        //        of rows in columns 1:k-1 looping backwards from k-1 to 1
        //
        j = k - 1;
    statement_120:
        //
        //           Undo the INTEGERerchanges (if any) of rows JJ and JP at each
        //           step J
        //
        //           (Here, J is a diagonal index)
        jj = j;
        jp = ipiv[j - 1];
        if (jp < 0) {
            jp = -jp;
            //              (Here, J is a diagonal index)
            j = j - 1;
        }
        //           (NOTE: Here, J is used to determine row length. Length J
        //           of the rows to swap back doesn't include diagonal element)
        j = j - 1;
        if (jp != jj && j >= 1) {
            Rswap(j, &a[(jp - 1)], lda, &a[(jj - 1)], lda);
        }
        if (j > 1) {
            goto statement_120;
        }
        //
        //        Set KB to the number of columns factorized
        //
        kb = k - 1;
        //
    }
    //
    //     End of Rlasyf
    //
}
