/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine DLANTP.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

REAL Rlantp(const char *norm, const char *uplo, const char *diag, INTEGER const n, REAL *ap, REAL *work) {
    REAL return_value = 0.0;
    //
    const REAL zero = 0.0;
    REAL value = 0.0;
    INTEGER k = 0;
    const REAL one = 1.0;
    INTEGER j = 0;
    INTEGER i = 0;
    REAL sum = 0.0;
    bool udiag = false;
    REAL scale = 0.0;
    if (n == 0) {
        value = zero;
    } else if (Mlsame(norm, "M")) {
        //
        // Find max(abs(A(i,j))).
        //
        k = 1;
        if (Mlsame(diag, "U")) {
            value = one;
            if (Mlsame(uplo, "U")) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = k; i <= k + j - 2; i = i + 1) {
                        sum = abs(ap[i - 1]);
                        if (value < sum || Risnan(sum)) {
                            value = sum;
                        }
                    }
                    k += j;
                }
            } else {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = k + 1; i <= k + n - j; i = i + 1) {
                        sum = abs(ap[i - 1]);
                        if (value < sum || Risnan(sum)) {
                            value = sum;
                        }
                    }
                    k += n - j + 1;
                }
            }
        } else {
            value = zero;
            if (Mlsame(uplo, "U")) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = k; i <= k + j - 1; i = i + 1) {
                        sum = abs(ap[i - 1]);
                        if (value < sum || Risnan(sum)) {
                            value = sum;
                        }
                    }
                    k += j;
                }
            } else {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = k; i <= k + n - j; i = i + 1) {
                        sum = abs(ap[i - 1]);
                        if (value < sum || Risnan(sum)) {
                            value = sum;
                        }
                    }
                    k += n - j + 1;
                }
            }
        }
    } else if ((Mlsame(norm, "O")) || (Mlsame(norm, "1"))) {
        //
        // Find norm1(A).
        //
        value = zero;
        k = 1;
        udiag = Mlsame(diag, "U");
        if (Mlsame(uplo, "U")) {
            for (j = 1; j <= n; j = j + 1) {
                if (udiag) {
                    sum = one;
                    for (i = k; i <= k + j - 2; i = i + 1) {
                        sum += abs(ap[i - 1]);
                    }
                } else {
                    sum = zero;
                    for (i = k; i <= k + j - 1; i = i + 1) {
                        sum += abs(ap[i - 1]);
                    }
                }
                k += j;
                if (value < sum || Risnan(sum)) {
                    value = sum;
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                if (udiag) {
                    sum = one;
                    for (i = k + 1; i <= k + n - j; i = i + 1) {
                        sum += abs(ap[i - 1]);
                    }
                } else {
                    sum = zero;
                    for (i = k; i <= k + n - j; i = i + 1) {
                        sum += abs(ap[i - 1]);
                    }
                }
                k += n - j + 1;
                if (value < sum || Risnan(sum)) {
                    value = sum;
                }
            }
        }
    } else if (Mlsame(norm, "I")) {
        //
        // Find normI(A).
        //
        k = 1;
        if (Mlsame(uplo, "U")) {
            if (Mlsame(diag, "U")) {
                for (i = 1; i <= n; i = i + 1) {
                    work[i - 1] = one;
                }
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= j - 1; i = i + 1) {
                        work[i - 1] += abs(ap[k - 1]);
                        k++;
                    }
                    k++;
                }
            } else {
                for (i = 1; i <= n; i = i + 1) {
                    work[i - 1] = zero;
                }
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= j; i = i + 1) {
                        work[i - 1] += abs(ap[k - 1]);
                        k++;
                    }
                }
            }
        } else {
            if (Mlsame(diag, "U")) {
                for (i = 1; i <= n; i = i + 1) {
                    work[i - 1] = one;
                }
                for (j = 1; j <= n; j = j + 1) {
                    k++;
                    for (i = j + 1; i <= n; i = i + 1) {
                        work[i - 1] += abs(ap[k - 1]);
                        k++;
                    }
                }
            } else {
                for (i = 1; i <= n; i = i + 1) {
                    work[i - 1] = zero;
                }
                for (j = 1; j <= n; j = j + 1) {
                    for (i = j; i <= n; i = i + 1) {
                        work[i - 1] += abs(ap[k - 1]);
                        k++;
                    }
                }
            }
        }
        value = zero;
        for (i = 1; i <= n; i = i + 1) {
            sum = work[i - 1];
            if (value < sum || Risnan(sum)) {
                value = sum;
            }
        }
    } else if ((Mlsame(norm, "F")) || (Mlsame(norm, "E"))) {
        //
        // Find normF(A).
        //
        if (Mlsame(uplo, "U")) {
            if (Mlsame(diag, "U")) {
                scale = one;
                sum = n;
                k = 2;
                for (j = 2; j <= n; j = j + 1) {
                    Rlassq(j - 1, &ap[k - 1], 1, scale, sum);
                    k += j;
                }
            } else {
                scale = zero;
                sum = one;
                k = 1;
                for (j = 1; j <= n; j = j + 1) {
                    Rlassq(j, &ap[k - 1], 1, scale, sum);
                    k += j;
                }
            }
        } else {
            if (Mlsame(diag, "U")) {
                scale = one;
                sum = n;
                k = 2;
                for (j = 1; j <= n - 1; j = j + 1) {
                    Rlassq(n - j, &ap[k - 1], 1, scale, sum);
                    k += n - j + 1;
                }
            } else {
                scale = zero;
                sum = one;
                k = 1;
                for (j = 1; j <= n; j = j + 1) {
                    Rlassq(n - j + 1, &ap[k - 1], 1, scale, sum);
                    k += n - j + 1;
                }
            }
        }
        value = scale * sqrt(sum);
    }
    //
    return_value = value;
    return return_value;
    //
    // End of Rlantp
    //
}
