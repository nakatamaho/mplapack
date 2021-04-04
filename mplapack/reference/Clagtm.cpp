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

void Clagtm(const char *trans, INTEGER const &n, INTEGER const &nrhs, REAL const &alpha, COMPLEX *dl, COMPLEX *d, COMPLEX *du, COMPLEX *x, INTEGER const &ldx, REAL const &beta, COMPLEX *b, INTEGER const &ldb) {
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    if (n == 0) {
        return;
    }
    //
    //     Multiply B by BETA if BETA.NE.1.
    //
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL one = 1.0;
    if (beta == zero) {
        for (j = 1; j <= nrhs; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = zero;
            }
        }
    } else if (beta == -one) {
        for (j = 1; j <= nrhs; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = -b[(i - 1) + (j - 1) * ldb];
            }
        }
    }
    //
    if (alpha == one) {
        if (Mlsame(trans, "N")) {
            //
            //           Compute B := B + A*X
            //
            for (j = 1; j <= nrhs; j = j + 1) {
                if (n == 1) {
                    b[(j - 1) * ldb] += d[1 - 1] * x[(j - 1) * ldx];
                } else {
                    b[(j - 1) * ldb] += d[1 - 1] * x[(j - 1) * ldx] + du[1 - 1] * x[(2 - 1) + (j - 1) * ldx];
                    b[(n - 1) + (j - 1) * ldb] += dl[(n - 1) - 1] * x[((n - 1) - 1) + (j - 1) * ldx] + d[n - 1] * x[(n - 1) + (j - 1) * ldx];
                    for (i = 2; i <= n - 1; i = i + 1) {
                        b[(i - 1) + (j - 1) * ldb] += dl[(i - 1) - 1] * x[((i - 1) - 1) + (j - 1) * ldx] + d[i - 1] * x[(i - 1) + (j - 1) * ldx] + du[i - 1] * x[((i + 1) - 1) + (j - 1) * ldx];
                    }
                }
            }
        } else if (Mlsame(trans, "T")) {
            //
            //           Compute B := B + A**T * X
            //
            for (j = 1; j <= nrhs; j = j + 1) {
                if (n == 1) {
                    b[(j - 1) * ldb] += d[1 - 1] * x[(j - 1) * ldx];
                } else {
                    b[(j - 1) * ldb] += d[1 - 1] * x[(j - 1) * ldx] + dl[1 - 1] * x[(2 - 1) + (j - 1) * ldx];
                    b[(n - 1) + (j - 1) * ldb] += du[(n - 1) - 1] * x[((n - 1) - 1) + (j - 1) * ldx] + d[n - 1] * x[(n - 1) + (j - 1) * ldx];
                    for (i = 2; i <= n - 1; i = i + 1) {
                        b[(i - 1) + (j - 1) * ldb] += du[(i - 1) - 1] * x[((i - 1) - 1) + (j - 1) * ldx] + d[i - 1] * x[(i - 1) + (j - 1) * ldx] + dl[i - 1] * x[((i + 1) - 1) + (j - 1) * ldx];
                    }
                }
            }
        } else if (Mlsame(trans, "C")) {
            //
            //           Compute B := B + A**H * X
            //
            for (j = 1; j <= nrhs; j = j + 1) {
                if (n == 1) {
                    b[(j - 1) * ldb] += conj(d[1 - 1]) * x[(j - 1) * ldx];
                } else {
                    b[(j - 1) * ldb] += conj(d[1 - 1]) * x[(j - 1) * ldx] + conj(dl[1 - 1]) * x[(2 - 1) + (j - 1) * ldx];
                    b[(n - 1) + (j - 1) * ldb] += conj(du[(n - 1) - 1]) * x[((n - 1) - 1) + (j - 1) * ldx] + conj(d[n - 1]) * x[(n - 1) + (j - 1) * ldx];
                    for (i = 2; i <= n - 1; i = i + 1) {
                        b[(i - 1) + (j - 1) * ldb] += conj(du[(i - 1) - 1]) * x[((i - 1) - 1) + (j - 1) * ldx] + conj(d[i - 1]) * x[(i - 1) + (j - 1) * ldx] + conj(dl[i - 1]) * x[((i + 1) - 1) + (j - 1) * ldx];
                    }
                }
            }
        }
    } else if (alpha == -one) {
        if (Mlsame(trans, "N")) {
            //
            //           Compute B := B - A*X
            //
            for (j = 1; j <= nrhs; j = j + 1) {
                if (n == 1) {
                    b[(j - 1) * ldb] = b[(j - 1) * ldb] - d[1 - 1] * x[(j - 1) * ldx];
                } else {
                    b[(j - 1) * ldb] = b[(j - 1) * ldb] - d[1 - 1] * x[(j - 1) * ldx] - du[1 - 1] * x[(2 - 1) + (j - 1) * ldx];
                    b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] - dl[(n - 1) - 1] * x[((n - 1) - 1) + (j - 1) * ldx] - d[n - 1] * x[(n - 1) + (j - 1) * ldx];
                    for (i = 2; i <= n - 1; i = i + 1) {
                        b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - dl[(i - 1) - 1] * x[((i - 1) - 1) + (j - 1) * ldx] - d[i - 1] * x[(i - 1) + (j - 1) * ldx] - du[i - 1] * x[((i + 1) - 1) + (j - 1) * ldx];
                    }
                }
            }
        } else if (Mlsame(trans, "T")) {
            //
            //           Compute B := B - A**T *X
            //
            for (j = 1; j <= nrhs; j = j + 1) {
                if (n == 1) {
                    b[(j - 1) * ldb] = b[(j - 1) * ldb] - d[1 - 1] * x[(j - 1) * ldx];
                } else {
                    b[(j - 1) * ldb] = b[(j - 1) * ldb] - d[1 - 1] * x[(j - 1) * ldx] - dl[1 - 1] * x[(2 - 1) + (j - 1) * ldx];
                    b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] - du[(n - 1) - 1] * x[((n - 1) - 1) + (j - 1) * ldx] - d[n - 1] * x[(n - 1) + (j - 1) * ldx];
                    for (i = 2; i <= n - 1; i = i + 1) {
                        b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - du[(i - 1) - 1] * x[((i - 1) - 1) + (j - 1) * ldx] - d[i - 1] * x[(i - 1) + (j - 1) * ldx] - dl[i - 1] * x[((i + 1) - 1) + (j - 1) * ldx];
                    }
                }
            }
        } else if (Mlsame(trans, "C")) {
            //
            //           Compute B := B - A**H *X
            //
            for (j = 1; j <= nrhs; j = j + 1) {
                if (n == 1) {
                    b[(j - 1) * ldb] = b[(j - 1) * ldb] - conj(d[1 - 1]) * x[(j - 1) * ldx];
                } else {
                    b[(j - 1) * ldb] = b[(j - 1) * ldb] - conj(d[1 - 1]) * x[(j - 1) * ldx] - conj(dl[1 - 1]) * x[(2 - 1) + (j - 1) * ldx];
                    b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] - conj(du[(n - 1) - 1]) * x[((n - 1) - 1) + (j - 1) * ldx] - conj(d[n - 1]) * x[(n - 1) + (j - 1) * ldx];
                    for (i = 2; i <= n - 1; i = i + 1) {
                        b[(i - 1) + (j - 1) * ldb] = b[(i - 1) + (j - 1) * ldb] - conj(du[(i - 1) - 1]) * x[((i - 1) - 1) + (j - 1) * ldx] - conj(d[i - 1]) * x[(i - 1) + (j - 1) * ldx] - conj(dl[i - 1]) * x[((i + 1) - 1) + (j - 1) * ldx];
                    }
                }
            }
        }
    }
    //
    //     End of Clagtm
    //
}
