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

// Derived from BLAS routine DGEMMTR.
// Original BLAS authors:
//   Martin Koehler

#include <mpblas.h>

void Rgemmtr(const char *uplo, const char *transa, const char *transb, INTEGER const n, INTEGER const k, REAL const alpha, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL const beta, REAL *c, INTEGER const ldc) {
    //
    // Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    // transposed and set  NROWA and NROWB  as the number of rows of  A
    // and  B  respectively.
    //
    bool nota = Mlsame(transa, "N");
    bool notb = Mlsame(transb, "N");
    INTEGER nrowa = 0;
    if (nota) {
        nrowa = n;
    } else {
        nrowa = k;
    }
    INTEGER nrowb = 0;
    if (notb) {
        nrowb = k;
    } else {
        nrowb = n;
    }
    bool upper = Mlsame(uplo, "U");
    //
    // Test the input parameters.
    //
    INTEGER info = 0;
    if ((!upper) && (!Mlsame(uplo, "L"))) {
        info = 1;
    } else if ((!nota) && (!Mlsame(transa, "C")) && (!Mlsame(transa, "T"))) {
        info = 2;
    } else if ((!notb) && (!Mlsame(transb, "C")) && (!Mlsame(transb, "T"))) {
        info = 3;
    } else if (n < 0) {
        info = 4;
    } else if (k < 0) {
        info = 5;
    } else if (lda < max((INTEGER)1, nrowa)) {
        info = 8;
    } else if (ldb < max((INTEGER)1, nrowb)) {
        info = 10;
    } else if (ldc < max((INTEGER)1, n)) {
        info = 13;
    }
    if (info != 0) {
        Mxerbla("Rgemmtr", info);
        return;
    }
    //
    // Quick return if possible.
    //
    if (n == 0) {
        return;
    }
    //
    // And if  alpha.eq.zero.
    //
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER istart = 0;
    INTEGER istop = 0;
    INTEGER i = 0;
    if (alpha == zero) {
        if (beta == zero) {
            for (j = 1; j <= n; j = j + 1) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = n;
                }
                //
                for (i = istart; i <= istop; i = i + 1) {
                    c[(i - 1) + (j - 1) * ldc] = zero;
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = n;
                }
                //
                for (i = istart; i <= istop; i = i + 1) {
                    c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                }
            }
        }
        return;
    }
    //
    // Start the operations.
    //
    const REAL one = 1.0;
    INTEGER l = 0;
    REAL temp = 0.0;
    if (notb) {
        if (nota) {
            //
            // Form  C := alpha*A*B + beta*C.
            //
            for (j = 1; j <= n; j = j + 1) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = n;
                }
                if (beta == zero) {
                    for (i = istart; i <= istop; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = zero;
                    }
                } else if (beta != one) {
                    for (i = istart; i <= istop; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
                for (l = 1; l <= k; l = l + 1) {
                    temp = alpha * b[(l - 1) + (j - 1) * ldb];
                    for (i = istart; i <= istop; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] += temp * a[(i - 1) + (l - 1) * lda];
                    }
                }
            }
        } else {
            //
            // Form  C := alpha*A**T*B + beta*C
            //
            for (j = 1; j <= n; j = j + 1) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = n;
                }
                //
                for (i = istart; i <= istop; i = i + 1) {
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
    } else {
        if (nota) {
            //
            // Form  C := alpha*A*B**T + beta*C
            //
            for (j = 1; j <= n; j = j + 1) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = n;
                }
                //
                if (beta == zero) {
                    for (i = istart; i <= istop; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = zero;
                    }
                } else if (beta != one) {
                    for (i = istart; i <= istop; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] = beta * c[(i - 1) + (j - 1) * ldc];
                    }
                }
                for (l = 1; l <= k; l = l + 1) {
                    temp = alpha * b[(j - 1) + (l - 1) * ldb];
                    for (i = istart; i <= istop; i = i + 1) {
                        c[(i - 1) + (j - 1) * ldc] += temp * a[(i - 1) + (l - 1) * lda];
                    }
                }
            }
        } else {
            //
            // Form  C := alpha*A**T*B**T + beta*C
            //
            for (j = 1; j <= n; j = j + 1) {
                if (upper) {
                    istart = 1;
                    istop = j;
                } else {
                    istart = j;
                    istop = n;
                }
                //
                for (i = istart; i <= istop; i = i + 1) {
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
    // End of Rgemm
    //
}
