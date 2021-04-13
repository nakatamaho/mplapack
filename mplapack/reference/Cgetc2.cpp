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

void Cgetc2(INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, INTEGER *jpiv, INTEGER &info) {
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
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Set constants to control overflow
    //
    REAL eps = Rlamch("P");
    REAL smlnum = Rlamch("S") / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Handle the case N=1 by itself
    //
    const REAL zero = 0.0;
    if (n == 1) {
        ipiv[1 - 1] = 1;
        jpiv[1 - 1] = 1;
        if (abs(a[(1 - 1)]) < smlnum) {
            info = 1;
            a[(1 - 1)] = COMPLEX(smlnum, zero);
        }
        return;
    }
    //
    //     Factorize A using complete pivoting.
    //     Set pivots less than SMIN to SMIN
    //
    INTEGER i = 0;
    REAL xmax = 0.0;
    INTEGER ip = 0;
    INTEGER jp = 0;
    INTEGER ipv = 0;
    INTEGER jpv = 0;
    REAL smin = 0.0;
    INTEGER j = 0;
    for (i = 1; i <= n - 1; i = i + 1) {
        //
        //        Find max element in matrix A
        //
        xmax = zero;
        for (ip = i; ip <= n; ip = ip + 1) {
            for (jp = i; jp <= n; jp = jp + 1) {
                if (abs(a[(ip - 1) + (jp - 1) * lda]) >= xmax) {
                    xmax = abs(a[(ip - 1) + (jp - 1) * lda]);
                    ipv = ip;
                    jpv = jp;
                }
            }
        }
        if (i == 1) {
            smin = max(REAL(eps * xmax), smlnum);
        }
        //
        //        Swap rows
        //
        if (ipv != i) {
            Cswap(n, &a[(ipv - 1)], lda, &a[(i - 1)], lda);
        }
        ipiv[i - 1] = ipv;
        //
        //        Swap columns
        //
        if (jpv != i) {
            Cswap(n, &a[(jpv - 1) * lda], 1, &a[(i - 1) * lda], 1);
        }
        jpiv[i - 1] = jpv;
        //
        //        Check for singularity
        //
        if (abs(a[(i - 1) + (i - 1) * lda]) < smin) {
            info = i;
            a[(i - 1) + (i - 1) * lda] = COMPLEX(smin, zero);
        }
        for (j = i + 1; j <= n; j = j + 1) {
            a[(j - 1) + (i - 1) * lda] = a[(j - 1) + (i - 1) * lda] / a[(i - 1) + (i - 1) * lda];
        }
        Cgeru(n - i, n - i, -COMPLEX(one), &a[((i + 1) - 1) + (i - 1) * lda], 1, &a[(i - 1) + ((i + 1) - 1) * lda], lda, &a[((i + 1) - 1) + ((i + 1) - 1) * lda], lda);
    }
    //
    if (abs(a[(n - 1) + (n - 1) * lda]) < smin) {
        info = n;
        a[(n - 1) + (n - 1) * lda] = COMPLEX(smin, zero);
    }
    //
    //     Set last pivots to N
    //
    ipiv[n - 1] = n;
    jpiv[n - 1] = n;
    //
    //     End of Cgetc2
    //
}
