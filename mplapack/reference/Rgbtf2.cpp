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

void Rgbtf2(INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, REAL *ab, INTEGER const ldab, INTEGER *ipiv, INTEGER &info) {
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
    //     KV is the number of superdiagonals in the factor U, allowing for
    //     fill-in.
    //
    INTEGER kv = ku + kl;
    //
    //     Test the input parameters.
    //
    info = 0;
    if (m < 0) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (kl < 0) {
        info = -3;
    } else if (ku < 0) {
        info = -4;
    } else if (ldab < kl + kv + 1) {
        info = -6;
    }
    if (info != 0) {
        Mxerbla("Rgbtf2", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     Gaussian elimination with partial pivoting
    //
    //     Set fill-in elements in columns KU+2 to KV to zero.
    //
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    for (j = ku + 2; j <= min(kv, n); j = j + 1) {
        for (i = kv - j + 2; i <= kl; i = i + 1) {
            ab[(i - 1) + (j - 1) * ldab] = zero;
        }
    }
    //
    //     JU is the index of the last column affected by the current stage
    //     of the factorization.
    //
    INTEGER ju = 1;
    //
    INTEGER km = 0;
    INTEGER jp = 0;
    const REAL one = 1.0;
    for (j = 1; j <= min(m, n); j = j + 1) {
        //
        //        Set fill-in elements in column J+KV to zero.
        //
        if (j + kv <= n) {
            for (i = 1; i <= kl; i = i + 1) {
                ab[(i - 1) + ((j + kv) - 1) * ldab] = zero;
            }
        }
        //
        //        Find pivot and test for singularity. KM is the number of
        //        subdiagonal elements in the current column.
        //
        km = min(kl, m - j);
        jp = iRamax(km + 1, ab[((kv + 1) - 1) + (j - 1) * ldab], 1);
        ipiv[j - 1] = jp + j - 1;
        if (ab[((kv + jp) - 1) + (j - 1) * ldab] != zero) {
            ju = max(ju, min(j + ku + jp - 1, n));
            //
            //           Apply interchange to columns J to JU.
            //
            if (jp != 1) {
                Rswap(ju - j + 1, ab[((kv + jp) - 1) + (j - 1) * ldab], ldab - 1, ab[((kv + 1) - 1) + (j - 1) * ldab], ldab - 1);
            }
            //
            if (km > 0) {
                //
                //              Compute multipliers.
                //
                Rscal(km, one / ab[((kv + 1) - 1) + (j - 1) * ldab], ab[((kv + 2) - 1) + (j - 1) * ldab], 1);
                //
                //              Update trailing submatrix within the band.
                //
                if (ju > j) {
                    Rger(km, ju - j, -one, ab[((kv + 2) - 1) + (j - 1) * ldab], 1, ab[(kv - 1) + ((j + 1) - 1) * ldab], ldab - 1, ab[((kv + 1) - 1) + ((j + 1) - 1) * ldab], ldab - 1);
                }
            }
        } else {
            //
            //           If pivot is zero, set INFO to the index of the pivot
            //           unless a zero pivot has already been found.
            //
            if (info == 0) {
                info = j;
            }
        }
    }
    //
    //     End of Rgbtf2
    //
}
