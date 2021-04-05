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

void Rtpmv(const char *uplo, const char *trans, const char *diag, INTEGER const n, REAL *ap, REAL *x, INTEGER const incx) {
    //
    //  -- Reference BLAS level2 routine --
    //  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
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
    //
    //     Test the input parameters.
    //
    INTEGER info = 0;
    if (!Mlsame(uplo, "U") && !Mlsame(uplo, "L")) {
        info = 1;
    } else if (!Mlsame(trans, "N") && !Mlsame(trans, "T") && !Mlsame(trans, "C")) {
        info = 2;
    } else if (!Mlsame(diag, "U") && !Mlsame(diag, "N")) {
        info = 3;
    } else if (n < 0) {
        info = 4;
    } else if (incx == 0) {
        info = 7;
    }
    if (info != 0) {
        Mxerbla("Rtpmv ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    if (n == 0) {
        return;
    }
    //
    bool nounit = Mlsame(diag, "N");
    //
    //     Set up the start poINTEGER in X if the increment is not unity. This
    //     will be  ( N - 1 )*INCX  too small for descending loops.
    //
    INTEGER kx = 0;
    if (incx <= 0) {
        kx = 1 - (n - 1) * incx;
    } else if (incx != 1) {
        kx = 1;
    }
    //
    //     Start the operations. In this version the elements of AP are
    //     accessed sequentially with one pass through AP.
    //
    INTEGER kk = 0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    REAL temp = 0.0;
    INTEGER k = 0;
    INTEGER i = 0;
    INTEGER jx = 0;
    INTEGER ix = 0;
    if (Mlsame(trans, "N")) {
        //
        //        Form  x:= A*x.
        //
        if (Mlsame(uplo, "U")) {
            kk = 1;
            if (incx == 1) {
                for (j = 1; j <= n; j = j + 1) {
                    if (x[j - 1] != zero) {
                        temp = x[j - 1];
                        k = kk;
                        for (i = 1; i <= j - 1; i = i + 1) {
                            x[i - 1] += temp * ap[k - 1];
                            k++;
                        }
                        if (nounit) {
                            x[j - 1] = x[j - 1] * ap[(kk + j - 1) - 1];
                        }
                    }
                    kk += j;
                }
            } else {
                jx = kx;
                for (j = 1; j <= n; j = j + 1) {
                    if (x[jx - 1] != zero) {
                        temp = x[jx - 1];
                        ix = kx;
                        for (k = kk; k <= kk + j - 2; k = k + 1) {
                            x[ix - 1] += temp * ap[k - 1];
                            ix += incx;
                        }
                        if (nounit) {
                            x[jx - 1] = x[jx - 1] * ap[(kk + j - 1) - 1];
                        }
                    }
                    jx += incx;
                    kk += j;
                }
            }
        } else {
            kk = (n * (n + 1)) / 2;
            if (incx == 1) {
                for (j = n; j >= 1; j = j - 1) {
                    if (x[j - 1] != zero) {
                        temp = x[j - 1];
                        k = kk;
                        for (i = n; i >= j + 1; i = i - 1) {
                            x[i - 1] += temp * ap[k - 1];
                            k = k - 1;
                        }
                        if (nounit) {
                            x[j - 1] = x[j - 1] * ap[(kk - n + j) - 1];
                        }
                    }
                    kk = kk - (n - j + 1);
                }
            } else {
                kx += (n - 1) * incx;
                jx = kx;
                for (j = n; j >= 1; j = j - 1) {
                    if (x[jx - 1] != zero) {
                        temp = x[jx - 1];
                        ix = kx;
                        for (k = kk; k >= kk - (n - (j + 1)); k = k - 1) {
                            x[ix - 1] += temp * ap[k - 1];
                            ix = ix - incx;
                        }
                        if (nounit) {
                            x[jx - 1] = x[jx - 1] * ap[(kk - n + j) - 1];
                        }
                    }
                    jx = jx - incx;
                    kk = kk - (n - j + 1);
                }
            }
        }
    } else {
        //
        //        Form  x := A**T*x.
        //
        if (Mlsame(uplo, "U")) {
            kk = (n * (n + 1)) / 2;
            if (incx == 1) {
                for (j = n; j >= 1; j = j - 1) {
                    temp = x[j - 1];
                    if (nounit) {
                        temp = temp * ap[kk - 1];
                    }
                    k = kk - 1;
                    for (i = j - 1; i >= 1; i = i - 1) {
                        temp += ap[k - 1] * x[i - 1];
                        k = k - 1;
                    }
                    x[j - 1] = temp;
                    kk = kk - j;
                }
            } else {
                jx = kx + (n - 1) * incx;
                for (j = n; j >= 1; j = j - 1) {
                    temp = x[jx - 1];
                    ix = jx;
                    if (nounit) {
                        temp = temp * ap[kk - 1];
                    }
                    for (k = kk - 1; k >= kk - j + 1; k = k - 1) {
                        ix = ix - incx;
                        temp += ap[k - 1] * x[ix - 1];
                    }
                    x[jx - 1] = temp;
                    jx = jx - incx;
                    kk = kk - j;
                }
            }
        } else {
            kk = 1;
            if (incx == 1) {
                for (j = 1; j <= n; j = j + 1) {
                    temp = x[j - 1];
                    if (nounit) {
                        temp = temp * ap[kk - 1];
                    }
                    k = kk + 1;
                    for (i = j + 1; i <= n; i = i + 1) {
                        temp += ap[k - 1] * x[i - 1];
                        k++;
                    }
                    x[j - 1] = temp;
                    kk += (n - j + 1);
                }
            } else {
                jx = kx;
                for (j = 1; j <= n; j = j + 1) {
                    temp = x[jx - 1];
                    ix = jx;
                    if (nounit) {
                        temp = temp * ap[kk - 1];
                    }
                    for (k = kk + 1; k <= kk + n - j; k = k + 1) {
                        ix += incx;
                        temp += ap[k - 1] * x[ix - 1];
                    }
                    x[jx - 1] = temp;
                    jx += incx;
                    kk += (n - j + 1);
                }
            }
        }
    }
    //
    //     End of Rtpmv .
    //
}
