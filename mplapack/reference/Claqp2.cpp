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

void Claqp2(INTEGER const m, INTEGER const n, INTEGER const offset, COMPLEX *a, INTEGER const lda, INTEGER *jpvt, COMPLEX *tau, REAL *vn1, REAL *vn2, COMPLEX *work) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER mn = min(m - offset, n);
    REAL tol3z = sqrt(Rlamch("Epsilon"));
    //
    //     Compute factorization.
    //
    INTEGER i = 0;
    INTEGER offpi = 0;
    INTEGER pvt = 0;
    INTEGER itemp = 0;
    COMPLEX aii = 0.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    INTEGER j = 0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL temp = 0.0;
    REAL temp2 = 0.0;
    for (i = 1; i <= mn; i = i + 1) {
        //
        offpi = offset + i;
        //
        //        Determine ith pivot column and swap if necessary.
        //
        pvt = (i - 1) + iRamax(n - i + 1, vn1[i - 1], 1);
        //
        if (pvt != i) {
            Cswap(m, &a[(pvt - 1) * lda], 1, &a[(i - 1) * lda], 1);
            itemp = jpvt[pvt - 1];
            jpvt[pvt - 1] = jpvt[i - 1];
            jpvt[i - 1] = itemp;
            vn1[pvt - 1] = vn1[i - 1];
            vn2[pvt - 1] = vn2[i - 1];
        }
        //
        //        Generate elementary reflector H(i).
        //
        if (offpi < m) {
            Clarfg(m - offpi + 1, &a[(offpi - 1) + (i - 1) * lda], &a[((offpi + 1) - 1) + (i - 1) * lda], 1, &tau[i - 1]);
        } else {
            Clarfg(1, &a[(m - 1) + (i - 1) * lda], &a[(m - 1) + (i - 1) * lda], 1, &tau[i - 1]);
        }
        //
        if (i < n) {
            //
            //           Apply H(i)**H to A(offset+i:m,i+1:n) from the left.
            //
            aii = a[(offpi - 1) + (i - 1) * lda];
            a[(offpi - 1) + (i - 1) * lda] = cone;
            Clarf("Left", m - offpi + 1, n - i, &a[(offpi - 1) + (i - 1) * lda], 1, conj(tau[i - 1]), &a[(offpi - 1) + ((i + 1) - 1) * lda], lda, &work[1 - 1]);
            a[(offpi - 1) + (i - 1) * lda] = aii;
        }
        //
        //        Update partial column norms.
        //
        for (j = i + 1; j <= n; j = j + 1) {
            if (vn1[j - 1] != zero) {
                //
                //              NOTE: The following 4 lines follow from the analysis in
                //              Lapack Working Note 176.
                //
                temp = one - pow2((abs(a[(offpi - 1) + (j - 1) * lda]) / vn1[j - 1]));
                temp = max(temp, zero);
                temp2 = temp * pow2((vn1[j - 1] / vn2[j - 1]));
                if (temp2 <= tol3z) {
                    if (offpi < m) {
                        vn1[j - 1] = RCnrm2(m - offpi, &a[((offpi + 1) - 1) + (j - 1) * lda], 1);
                        vn2[j - 1] = vn1[j - 1];
                    } else {
                        vn1[j - 1] = zero;
                        vn2[j - 1] = zero;
                    }
                } else {
                    vn1[j - 1] = vn1[j - 1] * sqrt(temp);
                }
            }
        }
        //
    }
    //
    //     End of Claqp2
    //
}
