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

void Clasr(const char *side, const char *pivot, const char *direct, INTEGER const m, INTEGER const n, REAL *c, REAL *s, COMPLEX *a, INTEGER const lda) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    INTEGER info = 0;
    if (!(Mlsame(side, "L") || Mlsame(side, "R"))) {
        info = 1;
    } else if (!(Mlsame(pivot, "V") || Mlsame(pivot, "T") || Mlsame(pivot, "B"))) {
        info = 2;
    } else if (!(Mlsame(direct, "F") || Mlsame(direct, "B"))) {
        info = 3;
    } else if (m < 0) {
        info = 4;
    } else if (n < 0) {
        info = 5;
    } else if (lda < max((INTEGER)1, m)) {
        info = 9;
    }
    if (info != 0) {
        Mxerbla("Clasr ", info);
        return;
    }
    //
    //     Quick return if possible
    //
    if ((m == 0) || (n == 0)) {
        return;
    }
    INTEGER j = 0;
    REAL ctemp = 0.0;
    REAL stemp = 0.0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    INTEGER i = 0;
    COMPLEX temp = 0.0;
    if (Mlsame(side, "L")) {
        //
        //        Form  P * A
        //
        if (Mlsame(pivot, "V")) {
            if (Mlsame(direct, "F")) {
                for (j = 1; j <= m - 1; j = j + 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[((j + 1) - 1) + (i - 1) * lda];
                            a[((j + 1) - 1) + (i - 1) * lda] = ctemp * temp - stemp * a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = stemp * temp + ctemp * a[(j - 1) + (i - 1) * lda];
                        }
                    }
                }
            } else if (Mlsame(direct, "B")) {
                for (j = m - 1; j >= 1; j = j - 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[((j + 1) - 1) + (i - 1) * lda];
                            a[((j + 1) - 1) + (i - 1) * lda] = ctemp * temp - stemp * a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = stemp * temp + ctemp * a[(j - 1) + (i - 1) * lda];
                        }
                    }
                }
            }
        } else if (Mlsame(pivot, "T")) {
            if (Mlsame(direct, "F")) {
                for (j = 2; j <= m; j = j + 1) {
                    ctemp = c[(j - 1) - 1];
                    stemp = s[(j - 1) - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = ctemp * temp - stemp * a[(i - 1) * lda];
                            a[(i - 1) * lda] = stemp * temp + ctemp * a[(i - 1) * lda];
                        }
                    }
                }
            } else if (Mlsame(direct, "B")) {
                for (j = m; j >= 2; j = j - 1) {
                    ctemp = c[(j - 1) - 1];
                    stemp = s[(j - 1) - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = ctemp * temp - stemp * a[(i - 1) * lda];
                            a[(i - 1) * lda] = stemp * temp + ctemp * a[(i - 1) * lda];
                        }
                    }
                }
            }
        } else if (Mlsame(pivot, "B")) {
            if (Mlsame(direct, "F")) {
                for (j = 1; j <= m - 1; j = j + 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = stemp * a[(m - 1) + (i - 1) * lda] + ctemp * temp;
                            a[(m - 1) + (i - 1) * lda] = ctemp * a[(m - 1) + (i - 1) * lda] - stemp * temp;
                        }
                    }
                }
            } else if (Mlsame(direct, "B")) {
                for (j = m - 1; j >= 1; j = j - 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= n; i = i + 1) {
                            temp = a[(j - 1) + (i - 1) * lda];
                            a[(j - 1) + (i - 1) * lda] = stemp * a[(m - 1) + (i - 1) * lda] + ctemp * temp;
                            a[(m - 1) + (i - 1) * lda] = ctemp * a[(m - 1) + (i - 1) * lda] - stemp * temp;
                        }
                    }
                }
            }
        }
    } else if (Mlsame(side, "R")) {
        //
        //        Form A * P**T
        //
        if (Mlsame(pivot, "V")) {
            if (Mlsame(direct, "F")) {
                for (j = 1; j <= n - 1; j = j + 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + ((j + 1) - 1) * lda];
                            a[(i - 1) + ((j + 1) - 1) * lda] = ctemp * temp - stemp * a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = stemp * temp + ctemp * a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
            } else if (Mlsame(direct, "B")) {
                for (j = n - 1; j >= 1; j = j - 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + ((j + 1) - 1) * lda];
                            a[(i - 1) + ((j + 1) - 1) * lda] = ctemp * temp - stemp * a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = stemp * temp + ctemp * a[(i - 1) + (j - 1) * lda];
                        }
                    }
                }
            }
        } else if (Mlsame(pivot, "T")) {
            if (Mlsame(direct, "F")) {
                for (j = 2; j <= n; j = j + 1) {
                    ctemp = c[(j - 1) - 1];
                    stemp = s[(j - 1) - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = ctemp * temp - stemp * a[(i - 1)];
                            a[(i - 1)] = stemp * temp + ctemp * a[(i - 1)];
                        }
                    }
                }
            } else if (Mlsame(direct, "B")) {
                for (j = n; j >= 2; j = j - 1) {
                    ctemp = c[(j - 1) - 1];
                    stemp = s[(j - 1) - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = ctemp * temp - stemp * a[(i - 1)];
                            a[(i - 1)] = stemp * temp + ctemp * a[(i - 1)];
                        }
                    }
                }
            }
        } else if (Mlsame(pivot, "B")) {
            if (Mlsame(direct, "F")) {
                for (j = 1; j <= n - 1; j = j + 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = stemp * a[(i - 1) + (n - 1) * lda] + ctemp * temp;
                            a[(i - 1) + (n - 1) * lda] = ctemp * a[(i - 1) + (n - 1) * lda] - stemp * temp;
                        }
                    }
                }
            } else if (Mlsame(direct, "B")) {
                for (j = n - 1; j >= 1; j = j - 1) {
                    ctemp = c[j - 1];
                    stemp = s[j - 1];
                    if ((ctemp != one) || (stemp != zero)) {
                        for (i = 1; i <= m; i = i + 1) {
                            temp = a[(i - 1) + (j - 1) * lda];
                            a[(i - 1) + (j - 1) * lda] = stemp * a[(i - 1) + (n - 1) * lda] + ctemp * temp;
                            a[(i - 1) + (n - 1) * lda] = ctemp * a[(i - 1) + (n - 1) * lda] - stemp * temp;
                        }
                    }
                }
            }
        }
    }
    //
    //     End of Clasr
    //
}
