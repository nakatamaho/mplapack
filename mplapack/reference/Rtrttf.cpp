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

void Rtrttf(const char *transr, const char *uplo, INTEGER const n, REAL *a, INTEGER const lda, REAL *arf, INTEGER &info) {
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
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    bool normaltransr = Mlsame(transr, "N");
    bool lower = Mlsame(uplo, "L");
    if (!normaltransr && !Mlsame(transr, "T")) {
        info = -1;
    } else if (!lower && !Mlsame(uplo, "U")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    }
    if (info != 0) {
        Mxerbla("Rtrttf", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        if (n == 1) {
            arf[0 - 1] = a[(0 - 1) + (0 - 1) * lda];
        }
        return;
    }
    //
    //     Size of array ARF(0:nt-1)
    //
    INTEGER nt = n * (n + 1) / 2;
    //
    //     Set N1 and N2 depending on LOWER: for N even N1=N2=K
    //
    INTEGER n2 = 0;
    INTEGER n1 = 0;
    if (lower) {
        n2 = n / 2;
        n1 = n - n2;
    } else {
        n1 = n / 2;
        n2 = n - n1;
    }
    //
    //     If N is odd, set NISODD = .TRUE., LDA=N+1 and A is (N+1)--by--K2.
    //     If N is even, set K = N/2 and NISODD = .FALSE., LDA=N and A is
    //     N--by--(N+1)/2.
    //
    INTEGER k = 0;
    bool nisodd = false;
    INTEGER np1x2 = 0;
    INTEGER nx2 = 0;
    if (mod(n, 2) == 0) {
        k = n / 2;
        nisodd = false;
        if (!lower) {
            np1x2 = n + n + 2;
        }
    } else {
        nisodd = true;
        if (!lower) {
            nx2 = n + n;
        }
    }
    //
    INTEGER ij = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    INTEGER l = 0;
    if (nisodd) {
        //
        //        N is odd
        //
        if (normaltransr) {
            //
            //           N is odd and TRANSR = 'N'
            //
            if (lower) {
                //
                //              N is odd, TRANSR = 'N', and UPLO = 'L'
                //
                ij = 0;
                for (j = 0; j <= n2; j = j + 1) {
                    for (i = n1; i <= n2 + j; i = i + 1) {
                        arf[ij - 1] = a[((n2 + j) - 1) + (i - 1) * lda];
                        ij++;
                    }
                    for (i = j; i <= n - 1; i = i + 1) {
                        arf[ij - 1] = a[(i - 1) + (j - 1) * lda];
                        ij++;
                    }
                }
                //
            } else {
                //
                //              N is odd, TRANSR = 'N', and UPLO = 'U'
                //
                ij = nt - n;
                for (j = n - 1; j >= n1; j = j - 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        arf[ij - 1] = a[(i - 1) + (j - 1) * lda];
                        ij++;
                    }
                    for (l = j - n1; l <= n1 - 1; l = l + 1) {
                        arf[ij - 1] = a[((j - n1) - 1) + (l - 1) * lda];
                        ij++;
                    }
                    ij = ij - nx2;
                }
                //
            }
            //
        } else {
            //
            //           N is odd and TRANSR = 'T'
            //
            if (lower) {
                //
                //              N is odd, TRANSR = 'T', and UPLO = 'L'
                //
                ij = 0;
                for (j = 0; j <= n2 - 1; j = j + 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        arf[ij - 1] = a[(j - 1) + (i - 1) * lda];
                        ij++;
                    }
                    for (i = n1 + j; i <= n - 1; i = i + 1) {
                        arf[ij - 1] = a[(i - 1) + ((n1 + j) - 1) * lda];
                        ij++;
                    }
                }
                for (j = n2; j <= n - 1; j = j + 1) {
                    for (i = 0; i <= n1 - 1; i = i + 1) {
                        arf[ij - 1] = a[(j - 1) + (i - 1) * lda];
                        ij++;
                    }
                }
                //
            } else {
                //
                //              N is odd, TRANSR = 'T', and UPLO = 'U'
                //
                ij = 0;
                for (j = 0; j <= n1; j = j + 1) {
                    for (i = n1; i <= n - 1; i = i + 1) {
                        arf[ij - 1] = a[(j - 1) + (i - 1) * lda];
                        ij++;
                    }
                }
                for (j = 0; j <= n1 - 1; j = j + 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        arf[ij - 1] = a[(i - 1) + (j - 1) * lda];
                        ij++;
                    }
                    for (l = n2 + j; l <= n - 1; l = l + 1) {
                        arf[ij - 1] = a[((n2 + j) - 1) + (l - 1) * lda];
                        ij++;
                    }
                }
                //
            }
            //
        }
        //
    } else {
        //
        //        N is even
        //
        if (normaltransr) {
            //
            //           N is even and TRANSR = 'N'
            //
            if (lower) {
                //
                //              N is even, TRANSR = 'N', and UPLO = 'L'
                //
                ij = 0;
                for (j = 0; j <= k - 1; j = j + 1) {
                    for (i = k; i <= k + j; i = i + 1) {
                        arf[ij - 1] = a[((k + j) - 1) + (i - 1) * lda];
                        ij++;
                    }
                    for (i = j; i <= n - 1; i = i + 1) {
                        arf[ij - 1] = a[(i - 1) + (j - 1) * lda];
                        ij++;
                    }
                }
                //
            } else {
                //
                //              N is even, TRANSR = 'N', and UPLO = 'U'
                //
                ij = nt - n - 1;
                for (j = n - 1; j >= k; j = j - 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        arf[ij - 1] = a[(i - 1) + (j - 1) * lda];
                        ij++;
                    }
                    for (l = j - k; l <= k - 1; l = l + 1) {
                        arf[ij - 1] = a[((j - k) - 1) + (l - 1) * lda];
                        ij++;
                    }
                    ij = ij - np1x2;
                }
                //
            }
            //
        } else {
            //
            //           N is even and TRANSR = 'T'
            //
            if (lower) {
                //
                //              N is even, TRANSR = 'T', and UPLO = 'L'
                //
                ij = 0;
                j = k;
                for (i = k; i <= n - 1; i = i + 1) {
                    arf[ij - 1] = a[(i - 1) + (j - 1) * lda];
                    ij++;
                }
                for (j = 0; j <= k - 2; j = j + 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        arf[ij - 1] = a[(j - 1) + (i - 1) * lda];
                        ij++;
                    }
                    for (i = k + 1 + j; i <= n - 1; i = i + 1) {
                        arf[ij - 1] = a[(i - 1) + ((k + 1 + j) - 1) * lda];
                        ij++;
                    }
                }
                for (j = k - 1; j <= n - 1; j = j + 1) {
                    for (i = 0; i <= k - 1; i = i + 1) {
                        arf[ij - 1] = a[(j - 1) + (i - 1) * lda];
                        ij++;
                    }
                }
                //
            } else {
                //
                //              N is even, TRANSR = 'T', and UPLO = 'U'
                //
                ij = 0;
                for (j = 0; j <= k; j = j + 1) {
                    for (i = k; i <= n - 1; i = i + 1) {
                        arf[ij - 1] = a[(j - 1) + (i - 1) * lda];
                        ij++;
                    }
                }
                for (j = 0; j <= k - 2; j = j + 1) {
                    for (i = 0; i <= j; i = i + 1) {
                        arf[ij - 1] = a[(i - 1) + (j - 1) * lda];
                        ij++;
                    }
                    for (l = k + 1 + j; l <= n - 1; l = l + 1) {
                        arf[ij - 1] = a[((k + 1 + j) - 1) + (l - 1) * lda];
                        ij++;
                    }
                }
                //              Note that here, on exit of the loop, J = K-1
                for (i = 0; i <= j; i = i + 1) {
                    arf[ij - 1] = a[(i - 1) + (j - 1) * lda];
                    ij++;
                }
                //
            }
            //
        }
        //
    }
    //
    //     End of Rtrttf
    //
}
