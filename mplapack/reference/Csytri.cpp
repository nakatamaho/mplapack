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

void Csytri(const char *uplo, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, COMPLEX *work, INTEGER &info) {
    bool upper = false;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    INTEGER k = 0;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    INTEGER kstep = 0;
    COMPLEX t = 0.0;
    COMPLEX ak = 0.0;
    COMPLEX akp1 = 0.0;
    COMPLEX akkp1 = 0.0;
    COMPLEX d = 0.0;
    INTEGER kp = 0;
    COMPLEX temp = 0.0;
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
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Csytri", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Check that the diagonal matrix D is nonsingular.
    //
    if (upper) {
        //
        //        Upper triangular storage: examine D from bottom to top
        //
        for (info = n; info >= 1; info = info - 1) {
            if (ipiv[info - 1] > 0 && a[(info - 1) + (info - 1) * lda] == zero) {
                return;
            }
        }
    } else {
        //
        //        Lower triangular storage: examine D from top to bottom.
        //
        for (info = 1; info <= n; info = info + 1) {
            if (ipiv[info - 1] > 0 && a[(info - 1) + (info - 1) * lda] == zero) {
                return;
            }
        }
    }
    info = 0;
    //
    if (upper) {
        //
        //        Compute inv(A) from the factorization A = U*D*U**T.
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        1 or 2, depending on the size of the diagonal blocks.
        //
        k = 1;
    statement_30:
        //
        //        If K > N, exit from loop.
        //
        if (k > n) {
            goto statement_40;
        }
        //
        if (ipiv[k - 1] > 0) {
            //
            //           1 x 1 diagonal block
            //
            //           Invert the diagonal block.
            //
            a[(k - 1) + (k - 1) * lda] = one / a[(k - 1) + (k - 1) * lda];
            //
            //           Compute column K of the inverse.
            //
            if (k > 1) {
                Ccopy(k - 1, &a[(k - 1) * lda], 1, work, 1);
                Csymv(uplo, k - 1, -one, a, lda, work, 1, zero, &a[(k - 1) * lda], 1);
                a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda] - Cdotu(k - 1, work, 1, &a[(k - 1) * lda], 1);
            }
            kstep = 1;
        } else {
            //
            //           2 x 2 diagonal block
            //
            //           Invert the diagonal block.
            //
            t = a[(k - 1) + ((k + 1) - 1) * lda];
            ak = a[(k - 1) + (k - 1) * lda] / t;
            akp1 = a[((k + 1) - 1) + ((k + 1) - 1) * lda] / t;
            akkp1 = a[(k - 1) + ((k + 1) - 1) * lda] / t;
            d = t * (ak * akp1 - one);
            a[(k - 1) + (k - 1) * lda] = akp1 / d;
            a[((k + 1) - 1) + ((k + 1) - 1) * lda] = ak / d;
            a[(k - 1) + ((k + 1) - 1) * lda] = -akkp1 / d;
            //
            //           Compute columns K and K+1 of the inverse.
            //
            if (k > 1) {
                Ccopy(k - 1, &a[(k - 1) * lda], 1, work, 1);
                Csymv(uplo, k - 1, -one, a, lda, work, 1, zero, &a[(k - 1) * lda], 1);
                a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda] - Cdotu(k - 1, work, 1, &a[(k - 1) * lda], 1);
                a[(k - 1) + ((k + 1) - 1) * lda] = a[(k - 1) + ((k + 1) - 1) * lda] - Cdotu(k - 1, &a[(k - 1) * lda], 1, &a[((k + 1) - 1) * lda], 1);
                Ccopy(k - 1, &a[((k + 1) - 1) * lda], 1, work, 1);
                Csymv(uplo, k - 1, -one, a, lda, work, 1, zero, &a[((k + 1) - 1) * lda], 1);
                a[((k + 1) - 1) + ((k + 1) - 1) * lda] = a[((k + 1) - 1) + ((k + 1) - 1) * lda] - Cdotu(k - 1, work, 1, &a[((k + 1) - 1) * lda], 1);
            }
            kstep = 2;
        }
        //
        kp = abs(ipiv[k - 1]);
        if (kp != k) {
            //
            //           Interchange rows and columns K and KP in the leading
            //           submatrix A(1:k+1,1:k+1)
            //
            Cswap(kp - 1, &a[(k - 1) * lda], 1, &a[(kp - 1) * lda], 1);
            Cswap(k - kp - 1, &a[((kp + 1) - 1) + (k - 1) * lda], 1, &a[(kp - 1) + ((kp + 1) - 1) * lda], lda);
            temp = a[(k - 1) + (k - 1) * lda];
            a[(k - 1) + (k - 1) * lda] = a[(kp - 1) + (kp - 1) * lda];
            a[(kp - 1) + (kp - 1) * lda] = temp;
            if (kstep == 2) {
                temp = a[(k - 1) + ((k + 1) - 1) * lda];
                a[(k - 1) + ((k + 1) - 1) * lda] = a[(kp - 1) + ((k + 1) - 1) * lda];
                a[(kp - 1) + ((k + 1) - 1) * lda] = temp;
            }
        }
        //
        k += kstep;
        goto statement_30;
    statement_40:;
        //
    } else {
        //
        //        Compute inv(A) from the factorization A = L*D*L**T.
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        1 or 2, depending on the size of the diagonal blocks.
        //
        k = n;
    statement_50:
        //
        //        If K < 1, exit from loop.
        //
        if (k < 1) {
            goto statement_60;
        }
        //
        if (ipiv[k - 1] > 0) {
            //
            //           1 x 1 diagonal block
            //
            //           Invert the diagonal block.
            //
            a[(k - 1) + (k - 1) * lda] = one / a[(k - 1) + (k - 1) * lda];
            //
            //           Compute column K of the inverse.
            //
            if (k < n) {
                Ccopy(n - k, &a[((k + 1) - 1) + (k - 1) * lda], 1, work, 1);
                Csymv(uplo, n - k, -one, &a[((k + 1) - 1) + ((k + 1) - 1) * lda], lda, work, 1, zero, &a[((k + 1) - 1) + (k - 1) * lda], 1);
                a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda] - Cdotu(n - k, work, 1, &a[((k + 1) - 1) + (k - 1) * lda], 1);
            }
            kstep = 1;
        } else {
            //
            //           2 x 2 diagonal block
            //
            //           Invert the diagonal block.
            //
            t = a[(k - 1) + ((k - 1) - 1) * lda];
            ak = a[((k - 1) - 1) + ((k - 1) - 1) * lda] / t;
            akp1 = a[(k - 1) + (k - 1) * lda] / t;
            akkp1 = a[(k - 1) + ((k - 1) - 1) * lda] / t;
            d = t * (ak * akp1 - one);
            a[((k - 1) - 1) + ((k - 1) - 1) * lda] = akp1 / d;
            a[(k - 1) + (k - 1) * lda] = ak / d;
            a[(k - 1) + ((k - 1) - 1) * lda] = -akkp1 / d;
            //
            //           Compute columns K-1 and K of the inverse.
            //
            if (k < n) {
                Ccopy(n - k, &a[((k + 1) - 1) + (k - 1) * lda], 1, work, 1);
                Csymv(uplo, n - k, -one, &a[((k + 1) - 1) + ((k + 1) - 1) * lda], lda, work, 1, zero, &a[((k + 1) - 1) + (k - 1) * lda], 1);
                a[(k - 1) + (k - 1) * lda] = a[(k - 1) + (k - 1) * lda] - Cdotu(n - k, work, 1, &a[((k + 1) - 1) + (k - 1) * lda], 1);
                a[(k - 1) + ((k - 1) - 1) * lda] = a[(k - 1) + ((k - 1) - 1) * lda] - Cdotu(n - k, &a[((k + 1) - 1) + (k - 1) * lda], 1, &a[((k + 1) - 1) + ((k - 1) - 1) * lda], 1);
                Ccopy(n - k, &a[((k + 1) - 1) + ((k - 1) - 1) * lda], 1, work, 1);
                Csymv(uplo, n - k, -one, &a[((k + 1) - 1) + ((k + 1) - 1) * lda], lda, work, 1, zero, &a[((k + 1) - 1) + ((k - 1) - 1) * lda], 1);
                a[((k - 1) - 1) + ((k - 1) - 1) * lda] = a[((k - 1) - 1) + ((k - 1) - 1) * lda] - Cdotu(n - k, work, 1, &a[((k + 1) - 1) + ((k - 1) - 1) * lda], 1);
            }
            kstep = 2;
        }
        //
        kp = abs(ipiv[k - 1]);
        if (kp != k) {
            //
            //           Interchange rows and columns K and KP in the trailing
            //           submatrix A(k-1:n,k-1:n)
            //
            if (kp < n) {
                Cswap(n - kp, &a[((kp + 1) - 1) + (k - 1) * lda], 1, &a[((kp + 1) - 1) + (kp - 1) * lda], 1);
            }
            Cswap(kp - k - 1, &a[((k + 1) - 1) + (k - 1) * lda], 1, &a[(kp - 1) + ((k + 1) - 1) * lda], lda);
            temp = a[(k - 1) + (k - 1) * lda];
            a[(k - 1) + (k - 1) * lda] = a[(kp - 1) + (kp - 1) * lda];
            a[(kp - 1) + (kp - 1) * lda] = temp;
            if (kstep == 2) {
                temp = a[(k - 1) + ((k - 1) - 1) * lda];
                a[(k - 1) + ((k - 1) - 1) * lda] = a[(kp - 1) + ((k - 1) - 1) * lda];
                a[(kp - 1) + ((k - 1) - 1) * lda] = temp;
            }
        }
        //
        k = k - kstep;
        goto statement_50;
    statement_60:;
    }
    //
    //     End of Csytri
    //
}
