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

REAL Rla_syrpvgrw(const char *uplo, INTEGER const &n, INTEGER const &info, REAL *a, INTEGER const &lda, REAL *af, INTEGER const &ldaf, INTEGER *ipiv, REAL *work) {
    REAL return_value = 0.0;
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
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    bool upper = Mlsame("Upper", uplo);
    INTEGER ncols = 0;
    if (info == 0) {
        if (upper) {
            ncols = 1;
        } else {
            ncols = n;
        }
    } else {
        ncols = info;
    }
    //
    REAL rpvgrw = 1.0;
    INTEGER i = 0;
    for (i = 1; i <= 2 * n; i = i + 1) {
        work[i - 1] = 0.0;
    }
    //
    //     Find the max magnitude entry of each column of A.  Compute the max
    //     for all N columns so we can apply the pivot permutation while
    //     looping below.  Assume a full factorization is the common case.
    //
    INTEGER j = 0;
    if (upper) {
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= j; i = i + 1) {
                work[(n + i) - 1] = max(abs(a[(i - 1) + (j - 1) * lda]), work[(n + i) - 1]);
                work[(n + j) - 1] = max(abs(a[(i - 1) + (j - 1) * lda]), work[(n + j) - 1]);
            }
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            for (i = j; i <= n; i = i + 1) {
                work[(n + i) - 1] = max(abs(a[(i - 1) + (j - 1) * lda]), work[(n + i) - 1]);
                work[(n + j) - 1] = max(abs(a[(i - 1) + (j - 1) * lda]), work[(n + j) - 1]);
            }
        }
    }
    //
    //     Now find the max magnitude entry of each column of U or L.  Also
    //     permute the magnitudes of A above so they're in the same order as
    //     the factor.
    //
    //     The iteration orders and permutations were copied from Rsytrs.
    //     Calls to SSWAP would be severe overkill.
    //
    INTEGER k = 0;
    INTEGER kp = 0;
    REAL tmp = 0.0;
    if (upper) {
        k = n;
        while (k < ncols && k > 0) {
            if (ipiv[k - 1] > 0) {
                //              1x1 pivot
                kp = ipiv[k - 1];
                if (kp != k) {
                    tmp = work[(n + k) - 1];
                    work[(n + k) - 1] = work[(n + kp) - 1];
                    work[(n + kp) - 1] = tmp;
                }
                for (i = 1; i <= k; i = i + 1) {
                    work[k - 1] = max(abs(af[(i - 1) + (k - 1) * ldaf]), work[k - 1]);
                }
                k = k - 1;
            } else {
                //              2x2 pivot
                kp = -ipiv[k - 1];
                tmp = work[(n + k - 1) - 1];
                work[(n + k - 1) - 1] = work[(n + kp) - 1];
                work[(n + kp) - 1] = tmp;
                for (i = 1; i <= k - 1; i = i + 1) {
                    work[k - 1] = max(abs(af[(i - 1) + (k - 1) * ldaf]), work[k - 1]);
                    work[(k - 1) - 1] = max(abs(af[(i - 1) + ((k - 1) - 1) * ldaf]), work[(k - 1) - 1]);
                }
                work[k - 1] = max(abs(af[(k - 1) + (k - 1) * ldaf]), work[k - 1]);
                k = k - 2;
            }
        }
        k = ncols;
        while (k <= n) {
            if (ipiv[k - 1] > 0) {
                kp = ipiv[k - 1];
                if (kp != k) {
                    tmp = work[(n + k) - 1];
                    work[(n + k) - 1] = work[(n + kp) - 1];
                    work[(n + kp) - 1] = tmp;
                }
                k++;
            } else {
                kp = -ipiv[k - 1];
                tmp = work[(n + k) - 1];
                work[(n + k) - 1] = work[(n + kp) - 1];
                work[(n + kp) - 1] = tmp;
                k += 2;
            }
        }
    } else {
        k = 1;
        while (k <= ncols) {
            if (ipiv[k - 1] > 0) {
                //              1x1 pivot
                kp = ipiv[k - 1];
                if (kp != k) {
                    tmp = work[(n + k) - 1];
                    work[(n + k) - 1] = work[(n + kp) - 1];
                    work[(n + kp) - 1] = tmp;
                }
                for (i = k; i <= n; i = i + 1) {
                    work[k - 1] = max(abs(af[(i - 1) + (k - 1) * ldaf]), work[k - 1]);
                }
                k++;
            } else {
                //              2x2 pivot
                kp = -ipiv[k - 1];
                tmp = work[(n + k + 1) - 1];
                work[(n + k + 1) - 1] = work[(n + kp) - 1];
                work[(n + kp) - 1] = tmp;
                for (i = k + 1; i <= n; i = i + 1) {
                    work[k - 1] = max(abs(af[(i - 1) + (k - 1) * ldaf]), work[k - 1]);
                    work[(k + 1) - 1] = max(abs(af[(i - 1) + ((k + 1) - 1) * ldaf]), work[(k + 1) - 1]);
                }
                work[k - 1] = max(abs(af[(k - 1) + (k - 1) * ldaf]), work[k - 1]);
                k += 2;
            }
        }
        k = ncols;
        while (k >= 1) {
            if (ipiv[k - 1] > 0) {
                kp = ipiv[k - 1];
                if (kp != k) {
                    tmp = work[(n + k) - 1];
                    work[(n + k) - 1] = work[(n + kp) - 1];
                    work[(n + kp) - 1] = tmp;
                }
                k = k - 1;
            } else {
                kp = -ipiv[k - 1];
                tmp = work[(n + k) - 1];
                work[(n + k) - 1] = work[(n + kp) - 1];
                work[(n + kp) - 1] = tmp;
                k = k - 2;
            }
        }
    }
    //
    //     Compute the *inverse* of the max element growth factor.  Dividing
    //     by zero would imply the largest entry of the factor's column is
    //     zero.  Than can happen when either the column of A is zero or
    //     massive pivots made the factor underflow to zero.  Neither counts
    //     as growth in itself, so simply ignore terms with zero
    //     denominators.
    //
    REAL umax = 0.0;
    REAL amax = 0.0;
    if (upper) {
        for (i = ncols; i <= n; i = i + 1) {
            umax = work[i - 1];
            amax = work[(n + i) - 1];
            if (umax != 0.0) {
                rpvgrw = min(amax / umax, rpvgrw);
            }
        }
    } else {
        for (i = 1; i <= ncols; i = i + 1) {
            umax = work[i - 1];
            amax = work[(n + i) - 1];
            if (umax != 0.0) {
                rpvgrw = min(amax / umax, rpvgrw);
            }
        }
    }
    //
    return_value = rpvgrw;
    return return_value;
}
