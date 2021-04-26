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
#include <mplapack.h>

void Rlatm5(INTEGER const prtype, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *c, INTEGER const ldc, REAL *d, INTEGER const ldd, REAL *e, INTEGER const lde, REAL *f, INTEGER const ldf, REAL *r, INTEGER const ldr, REAL *l, INTEGER const ldl, REAL const alpha, INTEGER &qblcka, INTEGER &qblckb) {
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
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER i = 0;
    INTEGER j = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    const REAL half = 0.5e+0;
    const REAL twenty = 2.0e+1;
    const REAL two = 2.0e+0;
    INTEGER k = 0;
    REAL reeps = 0.0;
    REAL imeps = 0.0;
    if (prtype == 1) {
        for (i = 1; i <= m; i = i + 1) {
            for (j = 1; j <= m; j = j + 1) {
                if (i == j) {
                    a[(i - 1) + (j - 1) * lda] = one;
                    d[(i - 1) + (j - 1) * ldd] = one;
                } else if (i == j - 1) {
                    a[(i - 1) + (j - 1) * lda] = -one;
                    d[(i - 1) + (j - 1) * ldd] = zero;
                } else {
                    a[(i - 1) + (j - 1) * lda] = zero;
                    d[(i - 1) + (j - 1) * ldd] = zero;
                }
            }
        }
        //
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (i == j) {
                    b[(i - 1) + (j - 1) * ldb] = one - alpha;
                    e[(i - 1) + (j - 1) * lde] = one;
                } else if (i == j - 1) {
                    b[(i - 1) + (j - 1) * ldb] = one;
                    e[(i - 1) + (j - 1) * lde] = zero;
                } else {
                    b[(i - 1) + (j - 1) * ldb] = zero;
                    e[(i - 1) + (j - 1) * lde] = zero;
                }
            }
        }
        //
        for (i = 1; i <= m; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                r[(i - 1) + (j - 1) * ldr] = (half - sin(castREAL(i / j))) * twenty;
                l[(i - 1) + (j - 1) * ldl] = r[(i - 1) + (j - 1) * ldr];
            }
        }
        //
    } else if (prtype == 2 || prtype == 3) {
        for (i = 1; i <= m; i = i + 1) {
            for (j = 1; j <= m; j = j + 1) {
                if (i <= j) {
                    a[(i - 1) + (j - 1) * lda] = (half - sin(castREAL(i))) * two;
                    d[(i - 1) + (j - 1) * ldd] = (half - sin(castREAL(i * j))) * two;
                } else {
                    a[(i - 1) + (j - 1) * lda] = zero;
                    d[(i - 1) + (j - 1) * ldd] = zero;
                }
            }
        }
        //
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (i <= j) {
                    b[(i - 1) + (j - 1) * ldb] = (half - sin(castREAL(i + j))) * two;
                    e[(i - 1) + (j - 1) * lde] = (half - sin(castREAL(j))) * two;
                } else {
                    b[(i - 1) + (j - 1) * ldb] = zero;
                    e[(i - 1) + (j - 1) * lde] = zero;
                }
            }
        }
        //
        for (i = 1; i <= m; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                r[(i - 1) + (j - 1) * ldr] = (half - sin(castREAL(i * j))) * twenty;
                l[(i - 1) + (j - 1) * ldl] = (half - sin(castREAL(i + j))) * twenty;
            }
        }
        //
        if (prtype == 3) {
            if (qblcka <= 1) {
                qblcka = 2;
            }
            for (k = 1; k <= m - 1; k = k + qblcka) {
                a[((k + 1) - 1) + ((k + 1) - 1) * lda] = a[(k - 1) + (k - 1) * lda];
                a[((k + 1) - 1) + (k - 1) * lda] = -sin(a[(k - 1) + ((k + 1) - 1) * lda]);
            }
            //
            if (qblckb <= 1) {
                qblckb = 2;
            }
            for (k = 1; k <= n - 1; k = k + qblckb) {
                b[((k + 1) - 1) + ((k + 1) - 1) * ldb] = b[(k - 1) + (k - 1) * ldb];
                b[((k + 1) - 1) + (k - 1) * ldb] = -sin(b[(k - 1) + ((k + 1) - 1) * ldb]);
            }
        }
        //
    } else if (prtype == 4) {
        for (i = 1; i <= m; i = i + 1) {
            for (j = 1; j <= m; j = j + 1) {
                a[(i - 1) + (j - 1) * lda] = (half - sin(castREAL(i * j))) * twenty;
                d[(i - 1) + (j - 1) * ldd] = (half - sin(castREAL(i + j))) * two;
            }
        }
        //
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                b[(i - 1) + (j - 1) * ldb] = (half - castREAL(sin((i + j)))) * twenty;
                e[(i - 1) + (j - 1) * lde] = (half - castREAL(sin((i * j)))) * two;
            }
        }
        //
        for (i = 1; i <= m; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                r[(i - 1) + (j - 1) * ldr] = (half - sin(castREAL(j / i))) * twenty;
                l[(i - 1) + (j - 1) * ldl] = (half - sin(castREAL(i * j))) * two;
            }
        }
        //
    } else if (prtype >= 5) {
        reeps = half * two * twenty / alpha;
        imeps = (half - two) / alpha;
        for (i = 1; i <= m; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                r[(i - 1) + (j - 1) * ldr] = (half - sin(castREAL(i * j))) * alpha / twenty;
                l[(i - 1) + (j - 1) * ldl] = (half - sin(castREAL(i + j))) * alpha / twenty;
            }
        }
        //
        for (i = 1; i <= m; i = i + 1) {
            d[(i - 1) + (i - 1) * ldd] = one;
        }
        //
        for (i = 1; i <= m; i = i + 1) {
            if (i <= 4) {
                a[(i - 1) + (i - 1) * lda] = one;
                if (i > 2) {
                    a[(i - 1) + (i - 1) * lda] = one + reeps;
                }
                if (mod(i, 2) != 0 && i < m) {
                    a[(i - 1) + ((i + 1) - 1) * lda] = imeps;
                } else if (i > 1) {
                    a[(i - 1) + ((i - 1) - 1) * lda] = -imeps;
                }
            } else if (i <= 8) {
                if (i <= 6) {
                    a[(i - 1) + (i - 1) * lda] = reeps;
                } else {
                    a[(i - 1) + (i - 1) * lda] = -reeps;
                }
                if (mod(i, 2) != 0 && i < m) {
                    a[(i - 1) + ((i + 1) - 1) * lda] = one;
                } else if (i > 1) {
                    a[(i - 1) + ((i - 1) - 1) * lda] = -one;
                }
            } else {
                a[(i - 1) + (i - 1) * lda] = one;
                if (mod(i, 2) != 0 && i < m) {
                    a[(i - 1) + ((i + 1) - 1) * lda] = imeps * 2;
                } else if (i > 1) {
                    a[(i - 1) + ((i - 1) - 1) * lda] = -imeps * 2;
                }
            }
        }
        //
        for (i = 1; i <= n; i = i + 1) {
            e[(i - 1) + (i - 1) * lde] = one;
            if (i <= 4) {
                b[(i - 1) + (i - 1) * ldb] = -one;
                if (i > 2) {
                    b[(i - 1) + (i - 1) * ldb] = one - reeps;
                }
                if (mod(i, 2) != 0 && i < n) {
                    b[(i - 1) + ((i + 1) - 1) * ldb] = imeps;
                } else if (i > 1) {
                    b[(i - 1) + ((i - 1) - 1) * ldb] = -imeps;
                }
            } else if (i <= 8) {
                if (i <= 6) {
                    b[(i - 1) + (i - 1) * ldb] = reeps;
                } else {
                    b[(i - 1) + (i - 1) * ldb] = -reeps;
                }
                if (mod(i, 2) != 0 && i < n) {
                    b[(i - 1) + ((i + 1) - 1) * ldb] = one + imeps;
                } else if (i > 1) {
                    b[(i - 1) + ((i - 1) - 1) * ldb] = -one - imeps;
                }
            } else {
                b[(i - 1) + (i - 1) * ldb] = one - reeps;
                if (mod(i, 2) != 0 && i < n) {
                    b[(i - 1) + ((i + 1) - 1) * ldb] = imeps * 2;
                } else if (i > 1) {
                    b[(i - 1) + ((i - 1) - 1) * ldb] = -imeps * 2;
                }
            }
        }
    }
    //
    //     Compute rhs (C, F)
    //
    Rgemm("N", "N", m, n, m, one, a, lda, r, ldr, zero, c, ldc);
    Rgemm("N", "N", m, n, n, -one, l, ldl, b, ldb, one, c, ldc);
    Rgemm("N", "N", m, n, m, one, d, ldd, r, ldr, zero, f, ldf);
    Rgemm("N", "N", m, n, n, -one, l, ldl, e, lde, one, f, ldf);
    //
    //     End of Rlatm5
    //
}
