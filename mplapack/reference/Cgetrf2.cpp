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

void Cgetrf2(INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER *ipiv, INTEGER &info) {
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
    //     Test the input parameters
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, m)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Cgetrf2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    const COMPLEX zero = (0.0, 0.0);
    REAL sfmin = 0.0;
    INTEGER i = 0;
    COMPLEX temp = 0.0;
    const COMPLEX one = (1.0, 0.0);
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER iinfo = 0;
    if (m == 1) {
        //
        //        Use unblocked code for one row case
        //        Just need to handle IPIV and INFO
        //
        ipiv[1 - 1] = 1;
        if (a[(1 - 1)] == zero) {
            info = 1;
        }
        //
    } else if (n == 1) {
        //
        //        Use unblocked code for one column case
        //
        //        Compute machine safe minimum
        //
        sfmin = Rlamch("S");
        //
        //        Find pivot and test for singularity
        //
        i = iCamax(m, &a[(1 - 1)], 1);
        ipiv[1 - 1] = i;
        if (a[(i - 1)] != zero) {
            //
            //           Apply the INTEGERerchange
            //
            if (i != 1) {
                temp = a[(1 - 1)];
                a[(1 - 1)] = a[(i - 1)];
                a[(i - 1)] = temp;
            }
            //
            //           Compute elements 2:M of the column
            //
            if (abs(a[(1 - 1)]) >= sfmin) {
                Cscal(m - 1, one / a[(1 - 1)], &a[(2 - 1)], 1);
            } else {
                for (i = 1; i <= m - 1; i = i + 1) {
                    a[((1 + i) - 1)] = a[((1 + i) - 1)] / a[(1 - 1)];
                }
            }
            //
        } else {
            info = 1;
        }
        //
    } else {
        //
        //        Use recursive code
        //
        n1 = min(m, n) / 2;
        n2 = n - n1;
        //
        //               [ A11 ]
        //        Factor [ --- ]
        //               [ A21 ]
        //
        Cgetrf2(m, n1, a, lda, ipiv, iinfo);
        //
        if (info == 0 && iinfo > 0) {
            info = iinfo;
        }
        //
        //                              [ A12 ]
        //        Apply INTEGERerchanges to [ --- ]
        //                              [ A22 ]
        //
        Claswp(n2, &a[((n1 + 1) - 1) * lda], lda, 1, n1, ipiv, 1);
        //
        //        Solve A12
        //
        Ctrsm("L", "L", "N", "U", n1, n2, one, a, lda, &a[((n1 + 1) - 1) * lda], lda);
        //
        //        Update A22
        //
        Cgemm("N", "N", m - n1, n2, n1, -one, &a[((n1 + 1) - 1)], lda, &a[((n1 + 1) - 1) * lda], lda, one, &a[((n1 + 1) - 1) + ((n1 + 1) - 1) * lda], lda);
        //
        //        Factor A22
        //
        Cgetrf2(m - n1, n2, &a[((n1 + 1) - 1) + ((n1 + 1) - 1) * lda], lda, &ipiv[(n1 + 1) - 1], iinfo);
        //
        //        Adjust INFO and the pivot indices
        //
        if (info == 0 && iinfo > 0) {
            info = iinfo + n1;
        }
        for (i = n1 + 1; i <= min(m, n); i = i + 1) {
            ipiv[i - 1] += n1;
        }
        //
        //        Apply INTEGERerchanges to A21
        //
        Claswp(n1, &a[(1 - 1)], lda, n1 + 1, min(m, n), ipiv, 1);
        //
    }
    //
    //     End of Cgetrf2
    //
}
