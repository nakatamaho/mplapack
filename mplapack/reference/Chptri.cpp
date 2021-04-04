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

void Chptri(const char *uplo, INTEGER const &n, COMPLEX *ap, INTEGER *ipiv, COMPLEX *work, INTEGER &info) {
    bool upper = false;
    INTEGER kp = 0;
    const COMPLEX zero = (0.0, 0.0);
    INTEGER k = 0;
    INTEGER kc = 0;
    INTEGER kcnext = 0;
    const REAL one = 1.0;
    const COMPLEX cone = (1.0, 0.0);
    INTEGER kstep = 0;
    REAL t = 0.0;
    REAL ak = 0.0;
    REAL akp1 = 0.0;
    COMPLEX akkp1 = 0.0;
    REAL d = 0.0;
    INTEGER kpc = 0;
    INTEGER kx = 0;
    INTEGER j = 0;
    COMPLEX temp = 0.0;
    INTEGER npp = 0;
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
        Mxerbla("Chptri", -info);
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
        kp = n * (n + 1) / 2;
        for (info = n; info >= 1; info = info - 1) {
            if (ipiv[info - 1] > 0 && ap[kp - 1] == zero) {
                return;
            }
            kp = kp - info;
        }
    } else {
        //
        //        Lower triangular storage: examine D from top to bottom.
        //
        kp = 1;
        for (info = 1; info <= n; info = info + 1) {
            if (ipiv[info - 1] > 0 && ap[kp - 1] == zero) {
                return;
            }
            kp += n - info + 1;
        }
    }
    info = 0;
    //
    if (upper) {
        //
        //        Compute inv(A) from the factorization A = U*D*U**H.
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        1 or 2, depending on the size of the diagonal blocks.
        //
        k = 1;
        kc = 1;
    statement_30:
        //
        //        If K > N, exit from loop.
        //
        if (k > n) {
            goto statement_50;
        }
        //
        kcnext = kc + k;
        if (ipiv[k - 1] > 0) {
            //
            //           1 x 1 diagonal block
            //
            //           Invert the diagonal block.
            //
            ap[(kc + k - 1) - 1] = one / ap[(kc + k - 1) - 1].real();
            //
            //           Compute column K of the inverse.
            //
            if (k > 1) {
                Ccopy(k - 1, ap[kc - 1], 1, work, 1);
                Chpmv(uplo, k - 1, -cone, ap, work, 1, zero, ap[kc - 1], 1);
                ap[(kc + k - 1) - 1] = ap[(kc + k - 1) - 1] - Cdotc[((k - 1) - 1) + (work - 1) * ldCdotc].real();
            }
            kstep = 1;
        } else {
            //
            //           2 x 2 diagonal block
            //
            //           Invert the diagonal block.
            //
            t = abs(ap[(kcnext + k - 1) - 1]);
            ak = ap[(kc + k - 1) - 1].real() / t;
            akp1 = ap[(kcnext + k) - 1].real() / t;
            akkp1 = ap[(kcnext + k - 1) - 1] / t;
            d = t * (ak * akp1 - one);
            ap[(kc + k - 1) - 1] = akp1 / d;
            ap[(kcnext + k) - 1] = ak / d;
            ap[(kcnext + k - 1) - 1] = -akkp1 / d;
            //
            //           Compute columns K and K+1 of the inverse.
            //
            if (k > 1) {
                Ccopy(k - 1, ap[kc - 1], 1, work, 1);
                Chpmv(uplo, k - 1, -cone, ap, work, 1, zero, ap[kc - 1], 1);
                ap[(kc + k - 1) - 1] = ap[(kc + k - 1) - 1] - Cdotc[((k - 1) - 1) + (work - 1) * ldCdotc].real();
                ap[(kcnext + k - 1) - 1] = ap[(kcnext + k - 1) - 1] - Cdotc[((k - 1) - 1) + (ap[kc - 1] - 1) * ldCdotc];
                Ccopy(k - 1, ap[kcnext - 1], 1, work, 1);
                Chpmv(uplo, k - 1, -cone, ap, work, 1, zero, ap[kcnext - 1], 1);
                ap[(kcnext + k) - 1] = ap[(kcnext + k) - 1] - Cdotc[((k - 1) - 1) + (work - 1) * ldCdotc].real();
            }
            kstep = 2;
            kcnext += k + 1;
        }
        //
        kp = abs(ipiv[k - 1]);
        if (kp != k) {
            //
            //           Interchange rows and columns K and KP in the leading
            //           submatrix A(1:k+1,1:k+1)
            //
            kpc = (kp - 1) * kp / 2 + 1;
            Cswap(kp - 1, ap[kc - 1], 1, ap[kpc - 1], 1);
            kx = kpc + kp - 1;
            for (j = kp + 1; j <= k - 1; j = j + 1) {
                kx += j - 1;
                temp = conj(ap[(kc + j - 1) - 1]);
                ap[(kc + j - 1) - 1] = conj(ap[kx - 1]);
                ap[kx - 1] = temp;
            }
            ap[(kc + kp - 1) - 1] = conj(ap[(kc + kp - 1) - 1]);
            temp = ap[(kc + k - 1) - 1];
            ap[(kc + k - 1) - 1] = ap[(kpc + kp - 1) - 1];
            ap[(kpc + kp - 1) - 1] = temp;
            if (kstep == 2) {
                temp = ap[(kc + k + k - 1) - 1];
                ap[(kc + k + k - 1) - 1] = ap[(kc + k + kp - 1) - 1];
                ap[(kc + k + kp - 1) - 1] = temp;
            }
        }
        //
        k += kstep;
        kc = kcnext;
        goto statement_30;
    statement_50:;
        //
    } else {
        //
        //        Compute inv(A) from the factorization A = L*D*L**H.
        //
        //        K is the main loop index, increasing from 1 to N in steps of
        //        1 or 2, depending on the size of the diagonal blocks.
        //
        npp = n * (n + 1) / 2;
        k = n;
        kc = npp;
    statement_60:
        //
        //        If K < 1, exit from loop.
        //
        if (k < 1) {
            goto statement_80;
        }
        //
        kcnext = kc - (n - k + 2);
        if (ipiv[k - 1] > 0) {
            //
            //           1 x 1 diagonal block
            //
            //           Invert the diagonal block.
            //
            ap[kc - 1] = one / ap[kc - 1].real();
            //
            //           Compute column K of the inverse.
            //
            if (k < n) {
                Ccopy(n - k, ap[(kc + 1) - 1], 1, work, 1);
                Chpmv(uplo, n - k, -cone, ap[(kc + n - k + 1) - 1], work, 1, zero, ap[(kc + 1) - 1], 1);
                ap[kc - 1] = ap[kc - 1] - Cdotc[((n - k) - 1) + (work - 1) * ldCdotc].real();
            }
            kstep = 1;
        } else {
            //
            //           2 x 2 diagonal block
            //
            //           Invert the diagonal block.
            //
            t = abs(ap[(kcnext + 1) - 1]);
            ak = ap[kcnext - 1].real() / t;
            akp1 = ap[kc - 1].real() / t;
            akkp1 = ap[(kcnext + 1) - 1] / t;
            d = t * (ak * akp1 - one);
            ap[kcnext - 1] = akp1 / d;
            ap[kc - 1] = ak / d;
            ap[(kcnext + 1) - 1] = -akkp1 / d;
            //
            //           Compute columns K-1 and K of the inverse.
            //
            if (k < n) {
                Ccopy(n - k, ap[(kc + 1) - 1], 1, work, 1);
                Chpmv(uplo, n - k, -cone, ap[(kc + (n - k + 1)) - 1], work, 1, zero, ap[(kc + 1) - 1], 1);
                ap[kc - 1] = ap[kc - 1] - Cdotc[((n - k) - 1) + (work - 1) * ldCdotc].real();
                ap[(kcnext + 1) - 1] = ap[(kcnext + 1) - 1] - Cdotc[((n - k) - 1) + ((ap[(kc + 1) - 1]) - 1) * ldCdotc];
                Ccopy(n - k, ap[(kcnext + 2) - 1], 1, work, 1);
                Chpmv(uplo, n - k, -cone, ap[(kc + (n - k + 1)) - 1], work, 1, zero, ap[(kcnext + 2) - 1], 1);
                ap[kcnext - 1] = ap[kcnext - 1] - Cdotc[((n - k) - 1) + (work - 1) * ldCdotc].real();
            }
            kstep = 2;
            kcnext = kcnext - (n - k + 3);
        }
        //
        kp = abs(ipiv[k - 1]);
        if (kp != k) {
            //
            //           Interchange rows and columns K and KP in the trailing
            //           submatrix A(k-1:n,k-1:n)
            //
            kpc = npp - (n - kp + 1) * (n - kp + 2) / 2 + 1;
            if (kp < n) {
                Cswap(n - kp, ap[(kc + kp - k + 1) - 1], 1, ap[(kpc + 1) - 1], 1);
            }
            kx = kc + kp - k;
            for (j = k + 1; j <= kp - 1; j = j + 1) {
                kx += n - j + 1;
                temp = conj(ap[(kc + j - k) - 1]);
                ap[(kc + j - k) - 1] = conj(ap[kx - 1]);
                ap[kx - 1] = temp;
            }
            ap[(kc + kp - k) - 1] = conj(ap[(kc + kp - k) - 1]);
            temp = ap[kc - 1];
            ap[kc - 1] = ap[kpc - 1];
            ap[kpc - 1] = temp;
            if (kstep == 2) {
                temp = ap[(kc - n + k - 1) - 1];
                ap[(kc - n + k - 1) - 1] = ap[(kc - n + kp - 1) - 1];
                ap[(kc - n + kp - 1) - 1] = temp;
            }
        }
        //
        k = k - kstep;
        kc = kcnext;
        goto statement_60;
    statement_80:;
    }
    //
    //     End of Chptri
    //
}
