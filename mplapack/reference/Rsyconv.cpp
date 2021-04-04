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

void Rsyconv(const char *uplo, const char *way, INTEGER const &n, REAL *a, INTEGER const &lda, INTEGER *ipiv, REAL *e, INTEGER &info) {
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
    //     .. External Functions ..
    //
    //     .. External Subroutines ..
    //     .. Local Scalars ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    bool upper = Mlsame(uplo, "U");
    bool convert = Mlsame(way, "C");
    if (!upper && !Mlsame(uplo, "L")) {
        info = -1;
    } else if (!convert && !Mlsame(way, "R")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
        //
    }
    if (info != 0) {
        Mxerbla("Rsyconv", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    INTEGER i = 0;
    const REAL zero = 0.0;
    INTEGER ip = 0;
    INTEGER j = 0;
    REAL temp = 0.0;
    if (upper) {
        //
        //      A is UPPER
        //
        //      Convert A (A is upper)
        //
        //        Convert VALUE
        //
        if (convert) {
            i = n;
            e[1 - 1] = zero;
            while (i > 1) {
                if (ipiv[i - 1] < 0) {
                    e[i - 1] = a[((i - 1) - 1) + (i - 1) * lda];
                    e[(i - 1) - 1] = zero;
                    a[((i - 1) - 1) + (i - 1) * lda] = zero;
                    i = i - 1;
                } else {
                    e[i - 1] = zero;
                }
                i = i - 1;
            }
            //
            //        Convert PERMUTATIONS
            //
            i = n;
            while (i >= 1) {
                if (ipiv[i - 1] > 0) {
                    ip = ipiv[i - 1];
                    if (i < n) {
                        for (j = i + 1; j <= n; j = j + 1) {
                            temp = a[(ip - 1) + (j - 1) * lda];
                            a[(ip - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = temp;
                        }
                    }
                } else {
                    ip = -ipiv[i - 1];
                    if (i < n) {
                        for (j = i + 1; j <= n; j = j + 1) {
                            temp = a[(ip - 1) + (j - 1) * lda];
                            a[(ip - 1) + (j - 1) * lda] = a[((i - 1) - 1) + (j - 1) * lda];
                            a[((i - 1) - 1) + (j - 1) * lda] = temp;
                        }
                    }
                    i = i - 1;
                }
                i = i - 1;
            }
            //
        } else {
            //
            //      Revert A (A is upper)
            //
            //        Revert PERMUTATIONS
            //
            i = 1;
            while (i <= n) {
                if (ipiv[i - 1] > 0) {
                    ip = ipiv[i - 1];
                    if (i < n) {
                        for (j = i + 1; j <= n; j = j + 1) {
                            temp = a[(ip - 1) + (j - 1) * lda];
                            a[(ip - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = temp;
                        }
                    }
                } else {
                    ip = -ipiv[i - 1];
                    i++;
                    if (i < n) {
                        for (j = i + 1; j <= n; j = j + 1) {
                            temp = a[(ip - 1) + (j - 1) * lda];
                            a[(ip - 1) + (j - 1) * lda] = a[((i - 1) - 1) + (j - 1) * lda];
                            a[((i - 1) - 1) + (j - 1) * lda] = temp;
                        }
                    }
                }
                i++;
            }
            //
            //        Revert VALUE
            //
            i = n;
            while (i > 1) {
                if (ipiv[i - 1] < 0) {
                    a[((i - 1) - 1) + (i - 1) * lda] = e[i - 1];
                    i = i - 1;
                }
                i = i - 1;
            }
        }
    } else {
        //
        //      A is LOWER
        //
        if (convert) {
            //
            //      Convert A (A is lower)
            //
            //        Convert VALUE
            //
            i = 1;
            e[n - 1] = zero;
            while (i <= n) {
                if (i < n && ipiv[i - 1] < 0) {
                    e[i - 1] = a[((i + 1) - 1) + (i - 1) * lda];
                    e[(i + 1) - 1] = zero;
                    a[((i + 1) - 1) + (i - 1) * lda] = zero;
                    i++;
                } else {
                    e[i - 1] = zero;
                }
                i++;
            }
            //
            //        Convert PERMUTATIONS
            //
            i = 1;
            while (i <= n) {
                if (ipiv[i - 1] > 0) {
                    ip = ipiv[i - 1];
                    if (i > 1) {
                        for (j = 1; j <= i - 1; j = j + 1) {
                            temp = a[(ip - 1) + (j - 1) * lda];
                            a[(ip - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = temp;
                        }
                    }
                } else {
                    ip = -ipiv[i - 1];
                    if (i > 1) {
                        for (j = 1; j <= i - 1; j = j + 1) {
                            temp = a[(ip - 1) + (j - 1) * lda];
                            a[(ip - 1) + (j - 1) * lda] = a[((i + 1) - 1) + (j - 1) * lda];
                            a[((i + 1) - 1) + (j - 1) * lda] = temp;
                        }
                    }
                    i++;
                }
                i++;
            }
        } else {
            //
            //      Revert A (A is lower)
            //
            //        Revert PERMUTATIONS
            //
            i = n;
            while (i >= 1) {
                if (ipiv[i - 1] > 0) {
                    ip = ipiv[i - 1];
                    if (i > 1) {
                        for (j = 1; j <= i - 1; j = j + 1) {
                            temp = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = a[(ip - 1) + (j - 1) * lda];
                            a[(ip - 1) + (j - 1) * lda] = temp;
                        }
                    }
                } else {
                    ip = -ipiv[i - 1];
                    i = i - 1;
                    if (i > 1) {
                        for (j = 1; j <= i - 1; j = j + 1) {
                            temp = a[((i + 1) - 1) + (j - 1) * lda];
                            a[((i + 1) - 1) + (j - 1) * lda] = a[(ip - 1) + (j - 1) * lda];
                            a[(ip - 1) + (j - 1) * lda] = temp;
                        }
                    }
                }
                i = i - 1;
            }
            //
            //        Revert VALUE
            //
            i = 1;
            while (i <= n - 1) {
                if (ipiv[i - 1] < 0) {
                    a[((i + 1) - 1) + (i - 1) * lda] = e[i - 1];
                    i++;
                }
                i++;
            }
        }
    }
    //
    //     End of Rsyconv
    //
}
