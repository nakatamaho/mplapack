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

void Cgbt01(INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, COMPLEX *a, INTEGER const lda, COMPLEX *afac, INTEGER const ldafac, INTEGER *ipiv, COMPLEX *work, REAL &resid) {
    //
    //  -- LAPACK test routine --
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
    //     Quick exit if M = 0 or N = 0.
    //
    const REAL zero = 0.0;
    resid = zero;
    if (m <= 0 || n <= 0) {
        return;
    }
    //
    //     Determine EPS and the norm of A.
    //
    REAL eps = Rlamch("Epsilon");
    INTEGER kd = ku + 1;
    REAL anorm = zero;
    INTEGER j = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    for (j = 1; j <= n; j = j + 1) {
        i1 = max(kd + 1 - j, 1);
        i2 = min(kd + m - j, kl + kd);
        if (i2 >= i1) {
            anorm = max(anorm, RCasum(i2 - i1 + 1, &a[(i1 - 1) + (j - 1) * lda], 1));
        }
    }
    //
    //     Compute one column at a time of L*U - A.
    //
    kd = kl + ku + 1;
    INTEGER ju = 0;
    INTEGER jl = 0;
    INTEGER lenj = 0;
    INTEGER i = 0;
    INTEGER il = 0;
    INTEGER iw = 0;
    COMPLEX t = 0.0;
    INTEGER ip = 0;
    INTEGER jua = 0;
    const REAL one = 1.0;
    for (j = 1; j <= n; j = j + 1) {
        //
        //        Copy the J-th column of U to WORK.
        //
        ju = min(kl + ku, j - 1);
        jl = min(kl, m - j);
        lenj = min(m, j) - j + ju + 1;
        if (lenj > 0) {
            Ccopy(lenj, &afac[((kd - ju) - 1) + (j - 1) * ldafac], 1, work, 1);
            for (i = lenj + 1; i <= ju + jl + 1; i = i + 1) {
                work[i - 1] = zero;
            }
            //
            //           Multiply by the unit lower triangular matrix L.  Note that L
            //           is stored as a product of transformations and permutations.
            //
            for (i = min(m - 1, j); i >= j - ju; i = i - 1) {
                il = min(kl, m - i);
                if (il > 0) {
                    iw = i - j + ju + 1;
                    t = work[iw - 1];
                    Caxpy(il, t, &afac[((kd + 1) - 1) + (i - 1) * ldafac], 1, &work[(iw + 1) - 1], 1);
                    ip = ipiv[i - 1];
                    if (i != ip) {
                        ip = ip - j + ju + 1;
                        work[iw - 1] = work[ip - 1];
                        work[ip - 1] = t;
                    }
                }
            }
            //
            //           Subtract the corresponding column of A.
            //
            jua = min(ju, ku);
            if (jua + jl + 1 > 0) {
                Caxpy(jua + jl + 1, -COMPLEX(one), &a[((ku + 1 - jua) - 1) + (j - 1) * lda], 1, &work[(ju + 1 - jua) - 1], 1);
            }
            //
            //           Compute the 1-norm of the column.
            //
            resid = max(resid, RCasum(ju + jl + 1, work, 1));
        }
    }
    //
    //     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
    //
    if (anorm <= zero) {
        if (resid != zero) {
            resid = one / eps;
        }
    } else {
        resid = ((resid / n.real()) / anorm) / eps;
    }
    //
    //     End of Cgbt01
    //
}
