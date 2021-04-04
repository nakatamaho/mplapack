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

void Rgtsv(INTEGER const &n, INTEGER const &nrhs, REAL *dl, REAL *d, REAL *du, REAL *b, INTEGER const &ldb, INTEGER &info) {
    INTEGER i = 0;
    const REAL zero = 0.0;
    REAL fact = 0.0;
    REAL temp = 0.0;
    INTEGER j = 0;
    //
    //  -- LAPACK driver routine --
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    if (n < 0) {
        info = -1;
    } else if (nrhs < 0) {
        info = -2;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -7;
    }
    if (info != 0) {
        Mxerbla("Rgtsv ", -info);
        return;
    }
    //
    if (n == 0) {
        return;
    }
    //
    if (nrhs == 1) {
        for (i = 1; i <= n - 2; i = i + 1) {
            if (abs(d[i - 1]) >= abs(dl[i - 1])) {
                //
                //              No row INTEGERerchange required
                //
                if (d[i - 1] != zero) {
                    fact = dl[i - 1] / d[i - 1];
                    d[(i + 1) - 1] = d[(i + 1) - 1] - fact * du[i - 1];
                    b[((i + 1) - 1)] = b[((i + 1) - 1)] - fact * b[(i - 1)];
                } else {
                    info = i;
                    return;
                }
                dl[i - 1] = zero;
            } else {
                //
                //              Interchange rows I and I+1
                //
                fact = d[i - 1] / dl[i - 1];
                d[i - 1] = dl[i - 1];
                temp = d[(i + 1) - 1];
                d[(i + 1) - 1] = du[i - 1] - fact * temp;
                dl[i - 1] = du[(i + 1) - 1];
                du[(i + 1) - 1] = -fact * dl[i - 1];
                du[i - 1] = temp;
                temp = b[(i - 1)];
                b[(i - 1)] = b[((i + 1) - 1)];
                b[((i + 1) - 1)] = temp - fact * b[((i + 1) - 1)];
            }
        }
        if (n > 1) {
            i = n - 1;
            if (abs(d[i - 1]) >= abs(dl[i - 1])) {
                if (d[i - 1] != zero) {
                    fact = dl[i - 1] / d[i - 1];
                    d[(i + 1) - 1] = d[(i + 1) - 1] - fact * du[i - 1];
                    b[((i + 1) - 1)] = b[((i + 1) - 1)] - fact * b[(i - 1)];
                } else {
                    info = i;
                    return;
                }
            } else {
                fact = d[i - 1] / dl[i - 1];
                d[i - 1] = dl[i - 1];
                temp = d[(i + 1) - 1];
                d[(i + 1) - 1] = du[i - 1] - fact * temp;
                du[i - 1] = temp;
                temp = b[(i - 1)];
                b[(i - 1)] = b[((i + 1) - 1)];
                b[((i + 1) - 1)] = temp - fact * b[((i + 1) - 1)];
            }
        }
        if (d[n - 1] == zero) {
            info = n;
            return;
        }
    } else {
        for (i = 1; i <= n - 2; i = i + 1) {
            if (abs(d[i - 1]) >= abs(dl[i - 1])) {
                //
                //              No row INTEGERerchange required
                //
                if (d[i - 1] != zero) {
                    fact = dl[i - 1] / d[i - 1];
                    d[(i + 1) - 1] = d[(i + 1) - 1] - fact * du[i - 1];
                    for (j = 1; j <= nrhs; j = j + 1) {
                        b[((i + 1) - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb] - fact * b[(i - 1) + (j - 1) * ldb];
                    }
                } else {
                    info = i;
                    return;
                }
                dl[i - 1] = zero;
            } else {
                //
                //              Interchange rows I and I+1
                //
                fact = d[i - 1] / dl[i - 1];
                d[i - 1] = dl[i - 1];
                temp = d[(i + 1) - 1];
                d[(i + 1) - 1] = du[i - 1] - fact * temp;
                dl[i - 1] = du[(i + 1) - 1];
                du[(i + 1) - 1] = -fact * dl[i - 1];
                du[i - 1] = temp;
                for (j = 1; j <= nrhs; j = j + 1) {
                    temp = b[(i - 1) + (j - 1) * ldb];
                    b[(i - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb];
                    b[((i + 1) - 1) + (j - 1) * ldb] = temp - fact * b[((i + 1) - 1) + (j - 1) * ldb];
                }
            }
        }
        if (n > 1) {
            i = n - 1;
            if (abs(d[i - 1]) >= abs(dl[i - 1])) {
                if (d[i - 1] != zero) {
                    fact = dl[i - 1] / d[i - 1];
                    d[(i + 1) - 1] = d[(i + 1) - 1] - fact * du[i - 1];
                    for (j = 1; j <= nrhs; j = j + 1) {
                        b[((i + 1) - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb] - fact * b[(i - 1) + (j - 1) * ldb];
                    }
                } else {
                    info = i;
                    return;
                }
            } else {
                fact = d[i - 1] / dl[i - 1];
                d[i - 1] = dl[i - 1];
                temp = d[(i + 1) - 1];
                d[(i + 1) - 1] = du[i - 1] - fact * temp;
                du[i - 1] = temp;
                for (j = 1; j <= nrhs; j = j + 1) {
                    temp = b[(i - 1) + (j - 1) * ldb];
                    b[(i - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb];
                    b[((i + 1) - 1) + (j - 1) * ldb] = temp - fact * b[((i + 1) - 1) + (j - 1) * ldb];
                }
            }
        }
        if (d[n - 1] == zero) {
            info = n;
            return;
        }
    }
    //
    //     Back solve with the matrix U from the factorization.
    //
    if (nrhs <= 2) {
        j = 1;
    statement_70:
        b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] / d[n - 1];
        if (n > 1) {
            b[((n - 1) - 1) + (j - 1) * ldb] = (b[((n - 1) - 1) + (j - 1) * ldb] - du[(n - 1) - 1] * b[(n - 1) + (j - 1) * ldb]) / d[(n - 1) - 1];
        }
        for (i = n - 2; i >= 1; i = i - 1) {
            b[(i - 1) + (j - 1) * ldb] = (b[(i - 1) + (j - 1) * ldb] - du[i - 1] * b[((i + 1) - 1) + (j - 1) * ldb] - dl[i - 1] * b[((i + 2) - 1) + (j - 1) * ldb]) / d[i - 1];
        }
        if (j < nrhs) {
            j++;
            goto statement_70;
        }
    } else {
        for (j = 1; j <= nrhs; j = j + 1) {
            b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] / d[n - 1];
            if (n > 1) {
                b[((n - 1) - 1) + (j - 1) * ldb] = (b[((n - 1) - 1) + (j - 1) * ldb] - du[(n - 1) - 1] * b[(n - 1) + (j - 1) * ldb]) / d[(n - 1) - 1];
            }
            for (i = n - 2; i >= 1; i = i - 1) {
                b[(i - 1) + (j - 1) * ldb] = (b[(i - 1) + (j - 1) * ldb] - du[i - 1] * b[((i + 1) - 1) + (j - 1) * ldb] - dl[i - 1] * b[((i + 2) - 1) + (j - 1) * ldb]) / d[i - 1];
            }
        }
    }
    //
    //     End of Rgtsv
    //
}
