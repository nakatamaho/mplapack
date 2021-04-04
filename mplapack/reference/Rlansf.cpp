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

REAL Rlansf(const char *norm, const char *transr, const char *uplo, INTEGER const &n, REAL *a, REAL *work) {
    REAL return_value = 0.0;
    a(dim1(0, star));
    work(dim1(0, star));
    const REAL zero = 0.0;
    INTEGER noe = 0;
    INTEGER ifm = 0;
    INTEGER ilu = 0;
    INTEGER lda = 0;
    INTEGER k = 0;
    REAL value = 0.0;
    INTEGER j = 0;
    INTEGER i = 0;
    REAL temp = 0.0;
    REAL s = 0.0;
    REAL aa = 0.0;
    INTEGER l = 0;
    INTEGER n1 = 0;
    REAL scale = 0.0;
    const REAL one = 1.0;
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
    if (n == 0) {
        return_value = zero;
        return return_value;
    } else if (n == 1) {
        return_value = abs(a[0 - 1]);
        return return_value;
    }
    //
    //     set noe = 1 if n is odd. if n is even set noe=0
    //
    noe = 1;
    if (mod(n, 2) == 0) {
        noe = 0;
    }
    //
    //     set ifm = 0 when form='T or 't' and 1 otherwise
    //
    ifm = 1;
    if (Mlsame(transr, "T")) {
        ifm = 0;
    }
    //
    //     set ilu = 0 when uplo='U or 'u' and 1 otherwise
    //
    ilu = 1;
    if (Mlsame(uplo, "U")) {
        ilu = 0;
    }
    //
    //     set lda = (n+1)/2 when ifm = 0
    //     set lda = n when ifm = 1 and noe = 1
    //     set lda = n+1 when ifm = 1 and noe = 0
    //
    if (ifm == 1) {
        if (noe == 1) {
            lda = n;
        } else {
            //           noe=0
            lda = n + 1;
        }
    } else {
        //        ifm=0
        lda = (n + 1) / 2;
    }
    //
    if (Mlsame(norm, "M")) {
        //
        //       Find max(abs(A(i,j))).
        //
        k = (n + 1) / 2;
        value = zero;
        if (noe == 1) {
            //           n is odd
            if (ifm == 1) {
                //           A is n by k
                for (j = 0; j <= k - 1; j = j + 1) {
                    for (i = 0; i <= n - 1; i = i + 1) {
                        temp = abs(a[(i + j * lda) - 1]);
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                }
            } else {
                //              xpose case; A is k by n
                for (j = 0; j <= n - 1; j = j + 1) {
                    for (i = 0; i <= k - 1; i = i + 1) {
                        temp = abs(a[(i + j * lda) - 1]);
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                }
            }
        } else {
            //           n is even
            if (ifm == 1) {
                //              A is n+1 by k
                for (j = 0; j <= k - 1; j = j + 1) {
                    for (i = 0; i <= n; i = i + 1) {
                        temp = abs(a[(i + j * lda) - 1]);
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                }
            } else {
                //              xpose case; A is k by n+1
                for (j = 0; j <= n; j = j + 1) {
                    for (i = 0; i <= k - 1; i = i + 1) {
                        temp = abs(a[(i + j * lda) - 1]);
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                }
            }
        }
    } else if ((Mlsame(norm, "I")) || (Mlsame(norm, "O")) || (norm == "1")) {
        //
        //        Find normI(A) ( = norm1(A), since A is symmetric).
        //
        if (ifm == 1) {
            k = n / 2;
            if (noe == 1) {
                //              n is odd
                if (ilu == 0) {
                    for (i = 0; i <= k - 1; i = i + 1) {
                        work[i - 1] = zero;
                    }
                    for (j = 0; j <= k; j = j + 1) {
                        s = zero;
                        for (i = 0; i <= k + j - 1; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       -> A(i,j+k)
                            s += aa;
                            work[i - 1] += aa;
                        }
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    -> A(j+k,j+k)
                        work[(j + k) - 1] = s + aa;
                        if (i == k + k) {
                            goto statement_10;
                        }
                        i++;
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    -> A(j,j)
                        work[j - 1] += aa;
                        s = zero;
                        for (l = j + 1; l <= k - 1; l = l + 1) {
                            i++;
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       -> A(l,j)
                            s += aa;
                            work[l - 1] += aa;
                        }
                        work[j - 1] += s;
                    }
                statement_10:
                    value = work[0 - 1];
                    for (i = 1; i <= n - 1; i = i + 1) {
                        temp = work[i - 1];
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                } else {
                    //                 ilu = 1
                    k++;
                    //                 k=(n+1)/2 for n odd and ilu=1
                    for (i = k; i <= n - 1; i = i + 1) {
                        work[i - 1] = zero;
                    }
                    for (j = k - 1; j >= 0; j = j - 1) {
                        s = zero;
                        for (i = 0; i <= j - 2; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       -> A(j+k,i+k)
                            s += aa;
                            work[(i + k) - 1] += aa;
                        }
                        if (j > 0) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       -> A(j+k,j+k)
                            s += aa;
                            work[(i + k) - 1] += s;
                            //                       i=j
                            i++;
                        }
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    -> A(j,j)
                        work[j - 1] = aa;
                        s = zero;
                        for (l = j + 1; l <= n - 1; l = l + 1) {
                            i++;
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       -> A(l,j)
                            s += aa;
                            work[l - 1] += aa;
                        }
                        work[j - 1] += s;
                    }
                    value = work[0 - 1];
                    for (i = 1; i <= n - 1; i = i + 1) {
                        temp = work[i - 1];
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                }
            } else {
                //              n is even
                if (ilu == 0) {
                    for (i = 0; i <= k - 1; i = i + 1) {
                        work[i - 1] = zero;
                    }
                    for (j = 0; j <= k - 1; j = j + 1) {
                        s = zero;
                        for (i = 0; i <= k + j - 1; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       -> A(i,j+k)
                            s += aa;
                            work[i - 1] += aa;
                        }
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    -> A(j+k,j+k)
                        work[(j + k) - 1] = s + aa;
                        i++;
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    -> A(j,j)
                        work[j - 1] += aa;
                        s = zero;
                        for (l = j + 1; l <= k - 1; l = l + 1) {
                            i++;
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       -> A(l,j)
                            s += aa;
                            work[l - 1] += aa;
                        }
                        work[j - 1] += s;
                    }
                    value = work[0 - 1];
                    for (i = 1; i <= n - 1; i = i + 1) {
                        temp = work[i - 1];
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                } else {
                    //                 ilu = 1
                    for (i = k; i <= n - 1; i = i + 1) {
                        work[i - 1] = zero;
                    }
                    for (j = k - 1; j >= 0; j = j - 1) {
                        s = zero;
                        for (i = 0; i <= j - 1; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       -> A(j+k,i+k)
                            s += aa;
                            work[(i + k) - 1] += aa;
                        }
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    -> A(j+k,j+k)
                        s += aa;
                        work[(i + k) - 1] += s;
                        //                    i=j
                        i++;
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    -> A(j,j)
                        work[j - 1] = aa;
                        s = zero;
                        for (l = j + 1; l <= n - 1; l = l + 1) {
                            i++;
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       -> A(l,j)
                            s += aa;
                            work[l - 1] += aa;
                        }
                        work[j - 1] += s;
                    }
                    value = work[0 - 1];
                    for (i = 1; i <= n - 1; i = i + 1) {
                        temp = work[i - 1];
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                }
            }
        } else {
            //           ifm=0
            k = n / 2;
            if (noe == 1) {
                //              n is odd
                if (ilu == 0) {
                    n1 = k;
                    //                 n/2
                    k++;
                    //                 k is the row size and lda
                    for (i = n1; i <= n - 1; i = i + 1) {
                        work[i - 1] = zero;
                    }
                    for (j = 0; j <= n1 - 1; j = j + 1) {
                        s = zero;
                        for (i = 0; i <= k - 1; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(j,n1+i)
                            work[(i + n1) - 1] += aa;
                            s += aa;
                        }
                        work[j - 1] = s;
                    }
                    //                 j=n1=k-1 is special
                    s = abs(a[(0 + j * lda) - 1]);
                    //                 A(k-1,k-1)
                    for (i = 1; i <= k - 1; i = i + 1) {
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    A(k-1,i+n1)
                        work[(i + n1) - 1] += aa;
                        s += aa;
                    }
                    work[j - 1] += s;
                    for (j = k; j <= n - 1; j = j + 1) {
                        s = zero;
                        for (i = 0; i <= j - k - 1; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(i,j-k)
                            work[i - 1] += aa;
                            s += aa;
                        }
                        //                    i=j-k
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    A(j-k,j-k)
                        s += aa;
                        work[(j - k) - 1] += s;
                        i++;
                        s = abs(a[(i + j * lda) - 1]);
                        //                    A(j,j)
                        for (l = j + 1; l <= n - 1; l = l + 1) {
                            i++;
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(j,l)
                            work[l - 1] += aa;
                            s += aa;
                        }
                        work[j - 1] += s;
                    }
                    value = work[0 - 1];
                    for (i = 1; i <= n - 1; i = i + 1) {
                        temp = work[i - 1];
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                } else {
                    //                 ilu=1
                    k++;
                    //                 k=(n+1)/2 for n odd and ilu=1
                    for (i = k; i <= n - 1; i = i + 1) {
                        work[i - 1] = zero;
                    }
                    for (j = 0; j <= k - 2; j = j + 1) {
                        //                    process
                        s = zero;
                        for (i = 0; i <= j - 1; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(j,i)
                            work[i - 1] += aa;
                            s += aa;
                        }
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    i=j so process of A(j,j)
                        s += aa;
                        work[j - 1] = s;
                        //                    is initialised here
                        i++;
                        //                    i=j process A(j+k,j+k)
                        aa = abs(a[(i + j * lda) - 1]);
                        s = aa;
                        for (l = k + j + 1; l <= n - 1; l = l + 1) {
                            i++;
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(l,k+j)
                            s += aa;
                            work[l - 1] += aa;
                        }
                        work[(k + j) - 1] += s;
                    }
                    //                 j=k-1 is special :process col A(k-1,0:k-1)
                    s = zero;
                    for (i = 0; i <= k - 2; i = i + 1) {
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    A(k,i)
                        work[i - 1] += aa;
                        s += aa;
                    }
                    //                 i=k-1
                    aa = abs(a[(i + j * lda) - 1]);
                    //                 A(k-1,k-1)
                    s += aa;
                    work[i - 1] = s;
                    //                 done with col j=k+1
                    for (j = k; j <= n - 1; j = j + 1) {
                        //                    process col j of A = A(j,0:k-1)
                        s = zero;
                        for (i = 0; i <= k - 1; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(j,i)
                            work[i - 1] += aa;
                            s += aa;
                        }
                        work[j - 1] += s;
                    }
                    value = work[0 - 1];
                    for (i = 1; i <= n - 1; i = i + 1) {
                        temp = work[i - 1];
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                }
            } else {
                //              n is even
                if (ilu == 0) {
                    for (i = k; i <= n - 1; i = i + 1) {
                        work[i - 1] = zero;
                    }
                    for (j = 0; j <= k - 1; j = j + 1) {
                        s = zero;
                        for (i = 0; i <= k - 1; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(j,i+k)
                            work[(i + k) - 1] += aa;
                            s += aa;
                        }
                        work[j - 1] = s;
                    }
                    //                 j=k
                    aa = abs(a[(0 + j * lda) - 1]);
                    //                 A(k,k)
                    s = aa;
                    for (i = 1; i <= k - 1; i = i + 1) {
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    A(k,k+i)
                        work[(i + k) - 1] += aa;
                        s += aa;
                    }
                    work[j - 1] += s;
                    for (j = k + 1; j <= n - 1; j = j + 1) {
                        s = zero;
                        for (i = 0; i <= j - 2 - k; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(i,j-k-1)
                            work[i - 1] += aa;
                            s += aa;
                        }
                        //                     i=j-1-k
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    A(j-k-1,j-k-1)
                        s += aa;
                        work[(j - k - 1) - 1] += s;
                        i++;
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    A(j,j)
                        s = aa;
                        for (l = j + 1; l <= n - 1; l = l + 1) {
                            i++;
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(j,l)
                            work[l - 1] += aa;
                            s += aa;
                        }
                        work[j - 1] += s;
                    }
                    //                 j=n
                    s = zero;
                    for (i = 0; i <= k - 2; i = i + 1) {
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    A(i,k-1)
                        work[i - 1] += aa;
                        s += aa;
                    }
                    //                 i=k-1
                    aa = abs(a[(i + j * lda) - 1]);
                    //                 A(k-1,k-1)
                    s += aa;
                    work[i - 1] += s;
                    value = work[0 - 1];
                    for (i = 1; i <= n - 1; i = i + 1) {
                        temp = work[i - 1];
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                } else {
                    //                 ilu=1
                    for (i = k; i <= n - 1; i = i + 1) {
                        work[i - 1] = zero;
                    }
                    //                 j=0 is special :process col A(k:n-1,k)
                    s = abs(a[0 - 1]);
                    //                 A(k,k)
                    for (i = 1; i <= k - 1; i = i + 1) {
                        aa = abs(a[i - 1]);
                        //                    A(k+i,k)
                        work[(i + k) - 1] += aa;
                        s += aa;
                    }
                    work[k - 1] += s;
                    for (j = 1; j <= k - 1; j = j + 1) {
                        //                    process
                        s = zero;
                        for (i = 0; i <= j - 2; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(j-1,i)
                            work[i - 1] += aa;
                            s += aa;
                        }
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    i=j-1 so process of A(j-1,j-1)
                        s += aa;
                        work[(j - 1) - 1] = s;
                        //                    is initialised here
                        i++;
                        //                    i=j process A(j+k,j+k)
                        aa = abs(a[(i + j * lda) - 1]);
                        s = aa;
                        for (l = k + j + 1; l <= n - 1; l = l + 1) {
                            i++;
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(l,k+j)
                            s += aa;
                            work[l - 1] += aa;
                        }
                        work[(k + j) - 1] += s;
                    }
                    //                 j=k is special :process col A(k,0:k-1)
                    s = zero;
                    for (i = 0; i <= k - 2; i = i + 1) {
                        aa = abs(a[(i + j * lda) - 1]);
                        //                    A(k,i)
                        work[i - 1] += aa;
                        s += aa;
                    }
                    //                 i=k-1
                    aa = abs(a[(i + j * lda) - 1]);
                    //                 A(k-1,k-1)
                    s += aa;
                    work[i - 1] = s;
                    //                 done with col j=k+1
                    for (j = k + 1; j <= n; j = j + 1) {
                        //                    process col j-1 of A = A(j-1,0:k-1)
                        s = zero;
                        for (i = 0; i <= k - 1; i = i + 1) {
                            aa = abs(a[(i + j * lda) - 1]);
                            //                       A(j-1,i)
                            work[i - 1] += aa;
                            s += aa;
                        }
                        work[(j - 1) - 1] += s;
                    }
                    value = work[0 - 1];
                    for (i = 1; i <= n - 1; i = i + 1) {
                        temp = work[i - 1];
                        if (value < temp || Risnan(temp)) {
                            value = temp;
                        }
                    }
                }
            }
        }
    } else if ((Mlsame(norm, "F")) || (Mlsame(norm, "E"))) {
        //
        //       Find normF(A).
        //
        k = (n + 1) / 2;
        scale = zero;
        s = one;
        if (noe == 1) {
            //           n is odd
            if (ifm == 1) {
                //              A is normal
                if (ilu == 0) {
                    //                 A is upper
                    for (j = 0; j <= k - 3; j = j + 1) {
                        Rlassq(k - j - 2, a[(k + j + 1 + j * lda) - 1], 1, scale, s);
                        //                    L at A(k,0)
                    }
                    for (j = 0; j <= k - 1; j = j + 1) {
                        Rlassq(k + j - 1, a[(0 + j * lda) - 1], 1, scale, s);
                        //                    trap U at A(0,0)
                    }
                    s += s;
                    //                 REAL s for the off diagonal elements
                    Rlassq(k - 1, a[k - 1], lda + 1, scale, s);
                    //                 tri L at A(k,0)
                    Rlassq(k, a[(k - 1) - 1], lda + 1, scale, s);
                    //                 tri U at A(k-1,0)
                } else {
                    //                 ilu=1 & A is lower
                    for (j = 0; j <= k - 1; j = j + 1) {
                        Rlassq(n - j - 1, a[(j + 1 + j * lda) - 1], 1, scale, s);
                        //                    trap L at A(0,0)
                    }
                    for (j = 0; j <= k - 2; j = j + 1) {
                        Rlassq(j, a[(0 + (1 + j) * lda) - 1], 1, scale, s);
                        //                    U at A(0,1)
                    }
                    s += s;
                    //                 REAL s for the off diagonal elements
                    Rlassq(k, a[0 - 1], lda + 1, scale, s);
                    //                 tri L at A(0,0)
                    Rlassq(k - 1, a[(0 + lda) - 1], lda + 1, scale, s);
                    //                 tri U at A(0,1)
                }
            } else {
                //              A is xpose
                if (ilu == 0) {
                    //                 A**T is upper
                    for (j = 1; j <= k - 2; j = j + 1) {
                        Rlassq(j, a[(0 + (k + j) * lda) - 1], 1, scale, s);
                        //                    U at A(0,k)
                    }
                    for (j = 0; j <= k - 2; j = j + 1) {
                        Rlassq(k, a[(0 + j * lda) - 1], 1, scale, s);
                        //                    k by k-1 rect. at A(0,0)
                    }
                    for (j = 0; j <= k - 2; j = j + 1) {
                        Rlassq(k - j - 1, a[(j + 1 + (j + k - 1) * lda) - 1], 1, scale, s);
                        //                    L at A(0,k-1)
                    }
                    s += s;
                    //                 REAL s for the off diagonal elements
                    Rlassq(k - 1, a[(0 + k * lda) - 1], lda + 1, scale, s);
                    //                 tri U at A(0,k)
                    Rlassq(k, a[(0 + (k - 1) * lda) - 1], lda + 1, scale, s);
                    //                 tri L at A(0,k-1)
                } else {
                    //                 A**T is lower
                    for (j = 1; j <= k - 1; j = j + 1) {
                        Rlassq(j, a[(0 + j * lda) - 1], 1, scale, s);
                        //                    U at A(0,0)
                    }
                    for (j = k; j <= n - 1; j = j + 1) {
                        Rlassq(k, a[(0 + j * lda) - 1], 1, scale, s);
                        //                    k by k-1 rect. at A(0,k)
                    }
                    for (j = 0; j <= k - 3; j = j + 1) {
                        Rlassq(k - j - 2, a[(j + 2 + j * lda) - 1], 1, scale, s);
                        //                    L at A(1,0)
                    }
                    s += s;
                    //                 REAL s for the off diagonal elements
                    Rlassq(k, a[0 - 1], lda + 1, scale, s);
                    //                 tri U at A(0,0)
                    Rlassq(k - 1, a[1 - 1], lda + 1, scale, s);
                    //                 tri L at A(1,0)
                }
            }
        } else {
            //           n is even
            if (ifm == 1) {
                //              A is normal
                if (ilu == 0) {
                    //                 A is upper
                    for (j = 0; j <= k - 2; j = j + 1) {
                        Rlassq(k - j - 1, a[(k + j + 2 + j * lda) - 1], 1, scale, s);
                        //                    L at A(k+1,0)
                    }
                    for (j = 0; j <= k - 1; j = j + 1) {
                        Rlassq(k + j, a[(0 + j * lda) - 1], 1, scale, s);
                        //                    trap U at A(0,0)
                    }
                    s += s;
                    //                 REAL s for the off diagonal elements
                    Rlassq(k, a[(k + 1) - 1], lda + 1, scale, s);
                    //                 tri L at A(k+1,0)
                    Rlassq(k, a[k - 1], lda + 1, scale, s);
                    //                 tri U at A(k,0)
                } else {
                    //                 ilu=1 & A is lower
                    for (j = 0; j <= k - 1; j = j + 1) {
                        Rlassq(n - j - 1, a[(j + 2 + j * lda) - 1], 1, scale, s);
                        //                    trap L at A(1,0)
                    }
                    for (j = 1; j <= k - 1; j = j + 1) {
                        Rlassq(j, a[(0 + j * lda) - 1], 1, scale, s);
                        //                    U at A(0,0)
                    }
                    s += s;
                    //                 REAL s for the off diagonal elements
                    Rlassq(k, a[1 - 1], lda + 1, scale, s);
                    //                 tri L at A(1,0)
                    Rlassq(k, a[0 - 1], lda + 1, scale, s);
                    //                 tri U at A(0,0)
                }
            } else {
                //              A is xpose
                if (ilu == 0) {
                    //                 A**T is upper
                    for (j = 1; j <= k - 1; j = j + 1) {
                        Rlassq(j, a[(0 + (k + 1 + j) * lda) - 1], 1, scale, s);
                        //                    U at A(0,k+1)
                    }
                    for (j = 0; j <= k - 1; j = j + 1) {
                        Rlassq(k, a[(0 + j * lda) - 1], 1, scale, s);
                        //                    k by k rect. at A(0,0)
                    }
                    for (j = 0; j <= k - 2; j = j + 1) {
                        Rlassq(k - j - 1, a[(j + 1 + (j + k) * lda) - 1], 1, scale, s);
                        //                    L at A(0,k)
                    }
                    s += s;
                    //                 REAL s for the off diagonal elements
                    Rlassq(k, a[(0 + (k + 1) * lda) - 1], lda + 1, scale, s);
                    //                 tri U at A(0,k+1)
                    Rlassq(k, a[(0 + k * lda) - 1], lda + 1, scale, s);
                    //                 tri L at A(0,k)
                } else {
                    //                 A**T is lower
                    for (j = 1; j <= k - 1; j = j + 1) {
                        Rlassq(j, a[(0 + (j + 1) * lda) - 1], 1, scale, s);
                        //                    U at A(0,1)
                    }
                    for (j = k + 1; j <= n; j = j + 1) {
                        Rlassq(k, a[(0 + j * lda) - 1], 1, scale, s);
                        //                    k by k rect. at A(0,k+1)
                    }
                    for (j = 0; j <= k - 2; j = j + 1) {
                        Rlassq(k - j - 1, a[(j + 1 + j * lda) - 1], 1, scale, s);
                        //                    L at A(0,0)
                    }
                    s += s;
                    //                 REAL s for the off diagonal elements
                    Rlassq(k, a[lda - 1], lda + 1, scale, s);
                    //                 tri L at A(0,1)
                    Rlassq(k, a[0 - 1], lda + 1, scale, s);
                    //                 tri U at A(0,0)
                }
            }
        }
        value = scale * sqrt(s);
    }
    //
    return_value = value;
    return return_value;
    //
    //     End of Rlansf
    //
}
