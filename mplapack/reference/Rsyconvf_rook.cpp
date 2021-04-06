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

void Rsyconvf_rook(const char *uplo, const char *way, INTEGER const n, REAL *a, INTEGER const lda, REAL *e, INTEGER *ipiv, INTEGER &info) {
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
        Mxerbla("Rsyconvf_rook", -info);
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
    INTEGER ip2 = 0;
    if (upper) {
        //
        //        Begin A is UPPER
        //
        if (convert) {
            //
            //           Convert A (A is upper)
            //
            //           Convert VALUE
            //
            //           Assign superdiagonal entries of D to array E and zero out
            //           corresponding entries in input storage A
            //
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
            //           Convert PERMUTATIONS
            //
            //           Apply permutations to submatrices of upper part of A
            //           in factorization order where i decreases from N to 1
            //
            i = n;
            while (i >= 1) {
                if (ipiv[i - 1] > 0) {
                    //
                    //                 1-by-1 pivot interchange
                    //
                    //                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
                    //
                    ip = ipiv[i - 1];
                    if (i < n) {
                        if (ip != i) {
                            Rswap(n - i, &a[(i - 1) + ((i + 1) - 1) * lda], lda, &a[(ip - 1) + ((i + 1) - 1) * lda], lda);
                        }
                    }
                    //
                } else {
                    //
                    //                 2-by-2 pivot interchange
                    //
                    //                 Swap rows i and IPIV(i) and i-1 and IPIV(i-1)
                    //                 in A(1:i,N-i:N)
                    //
                    ip = -ipiv[i - 1];
                    ip2 = -ipiv[(i - 1) - 1];
                    if (i < n) {
                        if (ip != i) {
                            Rswap(n - i, &a[(i - 1) + ((i + 1) - 1) * lda], lda, &a[(ip - 1) + ((i + 1) - 1) * lda], lda);
                        }
                        if (ip2 != (i - 1)) {
                            Rswap(n - i, &a[((i - 1) - 1) + ((i + 1) - 1) * lda], lda, &a[(ip2 - 1) + ((i + 1) - 1) * lda], lda);
                        }
                    }
                    i = i - 1;
                    //
                }
                i = i - 1;
            }
            //
        } else {
            //
            //           Revert A (A is upper)
            //
            //           Revert PERMUTATIONS
            //
            //           Apply permutations to submatrices of upper part of A
            //           in reverse factorization order where i increases from 1 to N
            //
            i = 1;
            while (i <= n) {
                if (ipiv[i - 1] > 0) {
                    //
                    //                 1-by-1 pivot interchange
                    //
                    //                 Swap rows i and IPIV(i) in A(1:i,N-i:N)
                    //
                    ip = ipiv[i - 1];
                    if (i < n) {
                        if (ip != i) {
                            Rswap(n - i, &a[(ip - 1) + ((i + 1) - 1) * lda], lda, &a[(i - 1) + ((i + 1) - 1) * lda], lda);
                        }
                    }
                    //
                } else {
                    //
                    //                 2-by-2 pivot interchange
                    //
                    //                 Swap rows i-1 and IPIV(i-1) and i and IPIV(i)
                    //                 in A(1:i,N-i:N)
                    //
                    i++;
                    ip = -ipiv[i - 1];
                    ip2 = -ipiv[(i - 1) - 1];
                    if (i < n) {
                        if (ip2 != (i - 1)) {
                            Rswap(n - i, &a[(ip2 - 1) + ((i + 1) - 1) * lda], lda, &a[((i - 1) - 1) + ((i + 1) - 1) * lda], lda);
                        }
                        if (ip != i) {
                            Rswap(n - i, &a[(ip - 1) + ((i + 1) - 1) * lda], lda, &a[(i - 1) + ((i + 1) - 1) * lda], lda);
                        }
                    }
                    //
                }
                i++;
            }
            //
            //           Revert VALUE
            //           Assign superdiagonal entries of D from array E to
            //           superdiagonal entries of A.
            //
            i = n;
            while (i > 1) {
                if (ipiv[i - 1] < 0) {
                    a[((i - 1) - 1) + (i - 1) * lda] = e[i - 1];
                    i = i - 1;
                }
                i = i - 1;
            }
            //
            //        End A is UPPER
            //
        }
        //
    } else {
        //
        //        Begin A is LOWER
        //
        if (convert) {
            //
            //           Convert A (A is lower)
            //
            //           Convert VALUE
            //           Assign subdiagonal entries of D to array E and zero out
            //           corresponding entries in input storage A
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
            //           Convert PERMUTATIONS
            //
            //           Apply permutations to submatrices of lower part of A
            //           in factorization order where i increases from 1 to N
            //
            i = 1;
            while (i <= n) {
                if (ipiv[i - 1] > 0) {
                    //
                    //                 1-by-1 pivot interchange
                    //
                    //                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
                    //
                    ip = ipiv[i - 1];
                    if (i > 1) {
                        if (ip != i) {
                            Rswap(i - 1, &a[(i - 1)], lda, &a[(ip - 1)], lda);
                        }
                    }
                    //
                } else {
                    //
                    //                 2-by-2 pivot interchange
                    //
                    //                 Swap rows i and IPIV(i) and i+1 and IPIV(i+1)
                    //                 in A(i:N,1:i-1)
                    //
                    ip = -ipiv[i - 1];
                    ip2 = -ipiv[(i + 1) - 1];
                    if (i > 1) {
                        if (ip != i) {
                            Rswap(i - 1, &a[(i - 1)], lda, &a[(ip - 1)], lda);
                        }
                        if (ip2 != (i + 1)) {
                            Rswap(i - 1, &a[((i + 1) - 1)], lda, &a[(ip2 - 1)], lda);
                        }
                    }
                    i++;
                    //
                }
                i++;
            }
            //
        } else {
            //
            //           Revert A (A is lower)
            //
            //           Revert PERMUTATIONS
            //
            //           Apply permutations to submatrices of lower part of A
            //           in reverse factorization order where i decreases from N to 1
            //
            i = n;
            while (i >= 1) {
                if (ipiv[i - 1] > 0) {
                    //
                    //                 1-by-1 pivot interchange
                    //
                    //                 Swap rows i and IPIV(i) in A(i:N,1:i-1)
                    //
                    ip = ipiv[i - 1];
                    if (i > 1) {
                        if (ip != i) {
                            Rswap(i - 1, &a[(ip - 1)], lda, &a[(i - 1)], lda);
                        }
                    }
                    //
                } else {
                    //
                    //                 2-by-2 pivot interchange
                    //
                    //                 Swap rows i+1 and IPIV(i+1) and i and IPIV(i)
                    //                 in A(i:N,1:i-1)
                    //
                    i = i - 1;
                    ip = -ipiv[i - 1];
                    ip2 = -ipiv[(i + 1) - 1];
                    if (i > 1) {
                        if (ip2 != (i + 1)) {
                            Rswap(i - 1, &a[(ip2 - 1)], lda, &a[((i + 1) - 1)], lda);
                        }
                        if (ip != i) {
                            Rswap(i - 1, &a[(ip - 1)], lda, &a[(i - 1)], lda);
                        }
                    }
                    //
                }
                i = i - 1;
            }
            //
            //           Revert VALUE
            //           Assign subdiagonal entries of D from array E to
            //           subgiagonal entries of A.
            //
            i = 1;
            while (i <= n - 1) {
                if (ipiv[i - 1] < 0) {
                    a[((i + 1) - 1) + (i - 1) * lda] = e[i - 1];
                    i++;
                }
                i++;
            }
            //
        }
        //
        //        End A is LOWER
        //
    }
    //
    //     End of Rsyconvf_rook
    //
}
