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

void Chemm(const char *side, const char *uplo, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc) {
    //
    //  -- Reference BLAS level3 routine --
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Parameters ..
    //     ..
    //
    //     Set NROWA as the number of rows of A.
    //
    INTEGER nrowa = 0;
    if (Mlsame(side, "L")) {
        nrowa = m;
    } else {
        nrowa = n;
    }
    bool upper = Mlsame(uplo, "U");
    //
    //     Test the input parameters.
    //
    INTEGER info = 0;
    if ((!Mlsame(side, "L")) && (!Mlsame(side, "R"))) {
        info = 1;
    } else if ((!upper) && (!Mlsame(uplo, "L"))) {
        info = 2;
    } else if (m < 0) {
        info = 3;
    } else if (n < 0) {
        info = 4;
    } else if (lda < max((INTEGER)1, nrowa)) {
        info = 7;
    } else if (ldb < max((INTEGER)1, m)) {
        info = 9;
    } else if (ldc < max((INTEGER)1, m)) {
        info = 12;
    }
    if (info != 0) {
        Mxerbla("Chemm ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    const COMPLEX one = COMPLEX(1.0, 0.0);
    if ((m == 0) || (n == 0) || ((alpha == zero) && (beta == one))) {
        return;
    }
    //
    //     And when  alpha.eq.zero.
    //
    INTEGER j = 0;
    INTEGER i = 0;
    if (alpha == zero) {
        if (beta == zero) {
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= m; i = i + 1) {
                    c[(i - 1) + (j - 1) * ldc] = zero;
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= m; i = i + 1) {
                    c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                }
            }
        }
        return;
    }
    //
    //     Start the operations.
    //
    COMPLEX temp1 = 0.0;
    COMPLEX temp2 = 0.0;
    INTEGER k = 0;
    if (Mlsame(side, "L")) {
        //
        //        Form  C := alpha*A*B + beta*C.
        //
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= m; i = i + 1) {
                    temp1 = alpha * b[(i - 1) + (j - 1) * ldb];
                    temp2 = zero;
                    for (k = 1; k <= i - 1; k = k + 1) {
                        c[(k - 1) + (j - 1) * ldc] += temp1 * a[(k - 1) + (i - 1) * lda];
                        temp2 += b[(k - 1) + (j - 1) * ldb] * conj(a[(k - 1) + (i - 1) * lda]);
                    }
                    if (beta == zero) {
                        c[(i - 1) + (j - 1) * ldc] = temp1 * (a[(i - 1) + (i - 1) * lda]).real() + alpha * temp2;
                    } else {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc] + temp1 * (a[(i - 1) + (i - 1) * lda]).real() + alpha * temp2;
                    }
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                for (i = m; i >= 1; i = i - 1) {
                    temp1 = alpha * b[(i - 1) + (j - 1) * ldb];
                    temp2 = zero;
                    for (k = i + 1; k <= m; k = k + 1) {
                        c[(k - 1) + (j - 1) * ldc] += temp1 * a[(k - 1) + (i - 1) * lda];
                        temp2 += b[(k - 1) + (j - 1) * ldb] * conj(a[(k - 1) + (i - 1) * lda]);
                    }
                    if (beta == zero) {
                        c[(i - 1) + (j - 1) * ldc] = temp1 * (a[(i - 1) + (i - 1) * lda]).real() + alpha * temp2;
                    } else {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc] + temp1 * (a[(i - 1) + (i - 1) * lda]).real() + alpha * temp2;
                    }
                }
            }
        }
    } else {
        //
        //        Form  C := alpha*B*A + beta*C.
        //
        for (j = 1; j <= n; j = j + 1) {
            temp1 = alpha * (a[(j - 1) + (j - 1) * lda]).real();
            if (beta == zero) {
                for (i = 1; i <= m; i = i + 1) {
                    c[(i - 1) + (j - 1) * ldc] = temp1 * b[(i - 1) + (j - 1) * ldb];
                }
            } else {
                for (i = 1; i <= m; i = i + 1) {
                    c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc] + temp1 * b[(i - 1) + (j - 1) * ldb];
                }
            }
            for (k = 1; k <= j - 1; k = k + 1) {
                if (upper) {
                    temp1 = alpha * a[(k - 1) + (j - 1) * lda];
                } else {
                    temp1 = alpha * conj(a[(j - 1) + (k - 1) * lda]);
                }
                for (i = 1; i <= m; i = i + 1) {
                    c[(i - 1) + (j - 1) * ldc] += temp1 * b[(i - 1) + (k - 1) * ldb];
                }
            }
            for (k = j + 1; k <= n; k = k + 1) {
                if (upper) {
                    temp1 = alpha * conj(a[(j - 1) + (k - 1) * lda]);
                } else {
                    temp1 = alpha * a[(k - 1) + (j - 1) * lda];
                }
                for (i = 1; i <= m; i = i + 1) {
                    c[(i - 1) + (j - 1) * ldc] += temp1 * b[(i - 1) + (k - 1) * ldb];
                }
            }
        }
    }
    //
    //     End of Chemm .
    //
}
