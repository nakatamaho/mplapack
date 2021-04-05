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

void Ctbmv(const char *uplo, const char *trans, const char *diag, INTEGER const n, INTEGER const k, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const incx) {
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
    //     .. Intrinsic Functions ..
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
    } else if (k < 0) {
        info = 5;
    } else if (lda < (k + 1)) {
        info = 7;
    } else if (incx == 0) {
        info = 9;
    }
    if (info != 0) {
        Mxerbla("Ctbmv ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    if (n == 0) {
        return;
    }
    //
    bool noconj = Mlsame(trans, "T");
    bool nounit = Mlsame(diag, "N");
    //
    //     Set up the start poINTEGER in X if the increment is not unity. This
    //     will be  ( N - 1 )*INCX   too small for descending loops.
    //
    INTEGER kx = 0;
    if (incx <= 0) {
        kx = 1 - (n - 1) * incx;
    } else if (incx != 1) {
        kx = 1;
    }
    //
    //     Start the operations. In this version the elements of A are
    //     accessed sequentially with one pass through A.
    //
    INTEGER kplus1 = 0;
    INTEGER j = 0;
    const COMPLEX zero = (0.0, 0.0);
    COMPLEX temp = 0.0;
    INTEGER l = 0;
    INTEGER i = 0;
    INTEGER jx = 0;
    INTEGER ix = 0;
    if (Mlsame(trans, "N")) {
        //
        //         Form  x := A*x.
        //
        if (Mlsame(uplo, "U")) {
            kplus1 = k + 1;
            if (incx == 1) {
                for (j = 1; j <= n; j = j + 1) {
                    if (x[j - 1] != zero) {
                        temp = x[j - 1];
                        l = kplus1 - j;
                        for (i = max((INTEGER)1, j - k); i <= j - 1; i = i + 1) {
                            x[i - 1] += temp * a[((l + i) - 1) + (j - 1) * lda];
                        }
                        if (nounit) {
                            x[j - 1] = x[j - 1] * a[(kplus1 - 1) + (j - 1) * lda];
                        }
                    }
                }
            } else {
                jx = kx;
                for (j = 1; j <= n; j = j + 1) {
                    if (x[jx - 1] != zero) {
                        temp = x[jx - 1];
                        ix = kx;
                        l = kplus1 - j;
                        for (i = max((INTEGER)1, j - k); i <= j - 1; i = i + 1) {
                            x[ix - 1] += temp * a[((l + i) - 1) + (j - 1) * lda];
                            ix += incx;
                        }
                        if (nounit) {
                            x[jx - 1] = x[jx - 1] * a[(kplus1 - 1) + (j - 1) * lda];
                        }
                    }
                    jx += incx;
                    if (j > k) {
                        kx += incx;
                    }
                }
            }
        } else {
            if (incx == 1) {
                for (j = n; j >= 1; j = j - 1) {
                    if (x[j - 1] != zero) {
                        temp = x[j - 1];
                        l = 1 - j;
                        for (i = min(n, j + k); i >= j + 1; i = i - 1) {
                            x[i - 1] += temp * a[((l + i) - 1) + (j - 1) * lda];
                        }
                        if (nounit) {
                            x[j - 1] = x[j - 1] * a[(j - 1) * lda];
                        }
                    }
                }
            } else {
                kx += (n - 1) * incx;
                jx = kx;
                for (j = n; j >= 1; j = j - 1) {
                    if (x[jx - 1] != zero) {
                        temp = x[jx - 1];
                        ix = kx;
                        l = 1 - j;
                        for (i = min(n, j + k); i >= j + 1; i = i - 1) {
                            x[ix - 1] += temp * a[((l + i) - 1) + (j - 1) * lda];
                            ix = ix - incx;
                        }
                        if (nounit) {
                            x[jx - 1] = x[jx - 1] * a[(j - 1) * lda];
                        }
                    }
                    jx = jx - incx;
                    if ((n - j) >= k) {
                        kx = kx - incx;
                    }
                }
            }
        }
    } else {
        //
        //        Form  x := A**T*x  or  x := A**H*x.
        //
        if (Mlsame(uplo, "U")) {
            kplus1 = k + 1;
            if (incx == 1) {
                for (j = n; j >= 1; j = j - 1) {
                    temp = x[j - 1];
                    l = kplus1 - j;
                    if (noconj) {
                        if (nounit) {
                            temp = temp * a[(kplus1 - 1) + (j - 1) * lda];
                        }
                        for (i = j - 1; i >= max((INTEGER)1, j - k); i = i - 1) {
                            temp += a[((l + i) - 1) + (j - 1) * lda] * x[i - 1];
                        }
                    } else {
                        if (nounit) {
                            temp = temp * conj(a[(kplus1 - 1) + (j - 1) * lda]);
                        }
                        for (i = j - 1; i >= max((INTEGER)1, j - k); i = i - 1) {
                            temp += conj(a[((l + i) - 1) + (j - 1) * lda]) * x[i - 1];
                        }
                    }
                    x[j - 1] = temp;
                }
            } else {
                kx += (n - 1) * incx;
                jx = kx;
                for (j = n; j >= 1; j = j - 1) {
                    temp = x[jx - 1];
                    kx = kx - incx;
                    ix = kx;
                    l = kplus1 - j;
                    if (noconj) {
                        if (nounit) {
                            temp = temp * a[(kplus1 - 1) + (j - 1) * lda];
                        }
                        for (i = j - 1; i >= max((INTEGER)1, j - k); i = i - 1) {
                            temp += a[((l + i) - 1) + (j - 1) * lda] * x[ix - 1];
                            ix = ix - incx;
                        }
                    } else {
                        if (nounit) {
                            temp = temp * conj(a[(kplus1 - 1) + (j - 1) * lda]);
                        }
                        for (i = j - 1; i >= max((INTEGER)1, j - k); i = i - 1) {
                            temp += conj(a[((l + i) - 1) + (j - 1) * lda]) * x[ix - 1];
                            ix = ix - incx;
                        }
                    }
                    x[jx - 1] = temp;
                    jx = jx - incx;
                }
            }
        } else {
            if (incx == 1) {
                for (j = 1; j <= n; j = j + 1) {
                    temp = x[j - 1];
                    l = 1 - j;
                    if (noconj) {
                        if (nounit) {
                            temp = temp * a[(j - 1) * lda];
                        }
                        for (i = j + 1; i <= min(n, j + k); i = i + 1) {
                            temp += a[((l + i) - 1) + (j - 1) * lda] * x[i - 1];
                        }
                    } else {
                        if (nounit) {
                            temp = temp * conj(a[(j - 1) * lda]);
                        }
                        for (i = j + 1; i <= min(n, j + k); i = i + 1) {
                            temp += conj(a[((l + i) - 1) + (j - 1) * lda]) * x[i - 1];
                        }
                    }
                    x[j - 1] = temp;
                }
            } else {
                jx = kx;
                for (j = 1; j <= n; j = j + 1) {
                    temp = x[jx - 1];
                    kx += incx;
                    ix = kx;
                    l = 1 - j;
                    if (noconj) {
                        if (nounit) {
                            temp = temp * a[(j - 1) * lda];
                        }
                        for (i = j + 1; i <= min(n, j + k); i = i + 1) {
                            temp += a[((l + i) - 1) + (j - 1) * lda] * x[ix - 1];
                            ix += incx;
                        }
                    } else {
                        if (nounit) {
                            temp = temp * conj(a[(j - 1) * lda]);
                        }
                        for (i = j + 1; i <= min(n, j + k); i = i + 1) {
                            temp += conj(a[((l + i) - 1) + (j - 1) * lda]) * x[ix - 1];
                            ix += incx;
                        }
                    }
                    x[jx - 1] = temp;
                    jx += incx;
                }
            }
        }
    }
    //
    //     End of Ctbmv .
    //
}
