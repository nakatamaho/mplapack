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

void Cgemm(const char *transa, const char *transb, INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX const alpha, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX const beta, COMPLEX *c, INTEGER const ldc) {
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
    //     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    //     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
    //     B  respectively are to be  transposed but  not conjugated  and set
    //     NROWA and NROWB  as the number of rows  of  A  and  B  respectively.
    //
    bool nota = Mlsame(transa, "N");
    bool notb = Mlsame(transb, "N");
    bool conja = Mlsame(transa, "C");
    bool conjb = Mlsame(transb, "C");
    INTEGER nrowa = 0;
    if (nota) {
        nrowa = m;
    } else {
        nrowa = k;
    }
    INTEGER nrowb = 0;
    if (notb) {
        nrowb = k;
    } else {
        nrowb = n;
    }
    //
    //     Test the input parameters.
    //
    INTEGER info = 0;
    if ((!nota) && (!conja) && (!Mlsame(transa, "T"))) {
        info = 1;
    } else if ((!notb) && (!conjb) && (!Mlsame(transb, "T"))) {
        info = 2;
    } else if (m < 0) {
        info = 3;
    } else if (n < 0) {
        info = 4;
    } else if (k < 0) {
        info = 5;
    } else if (lda < max((INTEGER)1, nrowa)) {
        info = 8;
    } else if (ldb < max((INTEGER)1, nrowb)) {
        info = 10;
    } else if (ldc < max((INTEGER)1, m)) {
        info = 13;
    }
    if (info != 0) {
        Mxerbla("Cgemm ", info);
        return;
    }
    //
    //     Quick return if possible.
    //
    const COMPLEX zero = (0.0, 0.0);
    const COMPLEX one = (1.0, 0.0);
    if ((m == 0) || (n == 0) || (((alpha == zero) || (k == 0)) && (beta == one))) {
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
    INTEGER l = 0;
    COMPLEX temp = 0.0;
    if (notb) {
        if (nota) {
            //
            //           Form  C := alpha*A*B + beta*C.
            //
            for (j = 1; j <= n; j = j + 1) {
                if (beta == zero) {
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = zero;
                    }
                } else if (beta != one) {
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
                for (l = 1; l <= k; l = l + 1) {
                    temp = alpha * b[(l - 1) + (j - 1) * ldb];
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] += temp * a[(i - 1) + (l - 1) * lda];
                    }
                }
            }
        } else if (conja) {
            //
            //           Form  C := alpha*A**H*B + beta*C.
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= m; i = i + 1) {
                    temp = zero;
                    for (l = 1; l <= k; l = l + 1) {
                        temp += conj(a[(l - 1) + (i - 1) * lda]) * b[(l - 1) + (j - 1) * ldb];
                    }
                    if (beta == zero) {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp;
                    } else {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp + beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
            }
        } else {
            //
            //           Form  C := alpha*A**T*B + beta*C
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= m; i = i + 1) {
                    temp = zero;
                    for (l = 1; l <= k; l = l + 1) {
                        temp += a[(l - 1) + (i - 1) * lda] * b[(l - 1) + (j - 1) * ldb];
                    }
                    if (beta == zero) {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp;
                    } else {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp + beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
            }
        }
    } else if (nota) {
        if (conjb) {
            //
            //           Form  C := alpha*A*B**H + beta*C.
            //
            for (j = 1; j <= n; j = j + 1) {
                if (beta == zero) {
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = zero;
                    }
                } else if (beta != one) {
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
                for (l = 1; l <= k; l = l + 1) {
                    temp = alpha * conj(b[(j - 1) + (l - 1) * ldb]);
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] += temp * a[(i - 1) + (l - 1) * lda];
                    }
                }
            }
        } else {
            //
            //           Form  C := alpha*A*B**T + beta*C
            //
            for (j = 1; j <= n; j = j + 1) {
                if (beta == zero) {
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = zero;
                    }
                } else if (beta != one) {
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
                for (l = 1; l <= k; l = l + 1) {
                    temp = alpha * b[(j - 1) + (l - 1) * ldb];
                    for (i = 1; i <= m; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] += temp * a[(i - 1) + (l - 1) * lda];
                    }
                }
            }
        }
    } else if (conja) {
        if (conjb) {
            //
            //           Form  C := alpha*A**H*B**H + beta*C.
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= m; i = i + 1) {
                    temp = zero;
                    for (l = 1; l <= k; l = l + 1) {
                        temp += conj(a[(l - 1) + (i - 1) * lda]) * conj(b[(j - 1) + (l - 1) * ldb]);
                    }
                    if (beta == zero) {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp;
                    } else {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp + beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
            }
        } else {
            //
            //           Form  C := alpha*A**H*B**T + beta*C
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= m; i = i + 1) {
                    temp = zero;
                    for (l = 1; l <= k; l = l + 1) {
                        temp += conj(a[(l - 1) + (i - 1) * lda]) * b[(j - 1) + (l - 1) * ldb];
                    }
                    if (beta == zero) {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp;
                    } else {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp + beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
            }
        }
    } else {
        if (conjb) {
            //
            //           Form  C := alpha*A**T*B**H + beta*C
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= m; i = i + 1) {
                    temp = zero;
                    for (l = 1; l <= k; l = l + 1) {
                        temp += a[(l - 1) + (i - 1) * lda] * conj(b[(j - 1) + (l - 1) * ldb]);
                    }
                    if (beta == zero) {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp;
                    } else {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp + beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
            }
        } else {
            //
            //           Form  C := alpha*A**T*B**T + beta*C
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= m; i = i + 1) {
                    temp = zero;
                    for (l = 1; l <= k; l = l + 1) {
                        temp += a[(l - 1) + (i - 1) * lda] * b[(j - 1) + (l - 1) * ldb];
                    }
                    if (beta == zero) {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp;
                    } else {
                        c[(i - 1) + (j - 1) * ldc] = alpha * temp + beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
            }
        }
    }
    //
    //     End of Cgemm .
    //
}
