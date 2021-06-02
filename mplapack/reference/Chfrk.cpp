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

void Chfrk(const char *transr, const char *uplo, const char *trans, INTEGER const n, INTEGER const k, REAL const alpha, COMPLEX *a, INTEGER const lda, REAL const beta, COMPLEX *c) {
    //
    //     Test the input parameters.
    //
    INTEGER info = 0;
    bool normaltransr = Mlsame(transr, "N");
    bool lower = Mlsame(uplo, "L");
    bool notrans = Mlsame(trans, "N");
    //
    INTEGER nrowa = 0;
    if (notrans) {
        nrowa = n;
    } else {
        nrowa = k;
    }
    //
    if (!normaltransr && !Mlsame(transr, "C")) {
        info = -1;
    } else if (!lower && !Mlsame(uplo, "U")) {
        info = -2;
    } else if (!notrans && !Mlsame(trans, "C")) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (k < 0) {
        info = -5;
    } else if (lda < max((INTEGER)1, nrowa)) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Chfrk", -info);
        return;
    }
    //
    //     Quick return if possible.
    //
    //     The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not
    //     done (it is in Cherk for example) and left in the general case.
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if ((n == 0) || (((alpha == zero) || (k == 0)) && (beta == one))) {
        return;
    }
    //
    INTEGER j = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    if ((alpha == zero) && (beta == zero)) {
        for (j = 1; j <= ((n * (n + 1)) / 2); j = j + 1) {
            c[j - 1] = czero;
        }
        return;
    }
    //
    COMPLEX calpha = COMPLEX(alpha, zero);
    COMPLEX cbeta = COMPLEX(beta, zero);
    //
    //     C is N-by-N.
    //     If N is odd, set NISODD = .TRUE., and N1 and N2.
    //     If N is even, NISODD = .FALSE., and NK.
    //
    bool nisodd = false;
    INTEGER nk = 0;
    INTEGER n2 = 0;
    INTEGER n1 = 0;
    if (mod(n, 2) == 0) {
        nisodd = false;
        nk = n / 2;
    } else {
        nisodd = true;
        if (lower) {
            n2 = n / 2;
            n1 = n - n2;
        } else {
            n1 = n / 2;
            n2 = n - n1;
        }
    }
    //
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
                if (notrans) {
                    //
                    //                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N'
                    //
                    Cherk("L", "N", n1, k, alpha, &a[(1 - 1)], lda, beta, &c[1 - 1], n);
                    Cherk("U", "N", n2, k, alpha, &a[((n1 + 1) - 1)], lda, beta, &c[(n + 1) - 1], n);
                    Cgemm("N", "C", n2, n1, k, calpha, &a[((n1 + 1) - 1)], lda, &a[(1 - 1)], lda, cbeta, &c[(n1 + 1) - 1], n);
                    //
                } else {
                    //
                    //                 N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'C'
                    //
                    Cherk("L", "C", n1, k, alpha, &a[(1 - 1)], lda, beta, &c[1 - 1], n);
                    Cherk("U", "C", n2, k, alpha, &a[((n1 + 1) - 1) * lda], lda, beta, &c[(n + 1) - 1], n);
                    Cgemm("C", "N", n2, n1, k, calpha, &a[((n1 + 1) - 1) * lda], lda, &a[(1 - 1)], lda, cbeta, &c[(n1 + 1) - 1], n);
                    //
                }
                //
            } else {
                //
                //              N is odd, TRANSR = 'N', and UPLO = 'U'
                //
                if (notrans) {
                    //
                    //                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'
                    //
                    Cherk("L", "N", n1, k, alpha, &a[(1 - 1)], lda, beta, &c[(n2 + 1) - 1], n);
                    Cherk("U", "N", n2, k, alpha, &a[(n2 - 1)], lda, beta, &c[(n1 + 1) - 1], n);
                    Cgemm("N", "C", n1, n2, k, calpha, &a[(1 - 1)], lda, &a[(n2 - 1)], lda, cbeta, &c[1 - 1], n);
                    //
                } else {
                    //
                    //                 N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'C'
                    //
                    Cherk("L", "C", n1, k, alpha, &a[(1 - 1)], lda, beta, &c[(n2 + 1) - 1], n);
                    Cherk("U", "C", n2, k, alpha, &a[(n2 - 1) * lda], lda, beta, &c[(n1 + 1) - 1], n);
                    Cgemm("C", "N", n1, n2, k, calpha, &a[(1 - 1)], lda, &a[(n2 - 1) * lda], lda, cbeta, &c[1 - 1], n);
                    //
                }
                //
            }
            //
        } else {
            //
            //           N is odd, and TRANSR = 'C'
            //
            if (lower) {
                //
                //              N is odd, TRANSR = 'C', and UPLO = 'L'
                //
                if (notrans) {
                    //
                    //                 N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'N'
                    //
                    Cherk("U", "N", n1, k, alpha, &a[(1 - 1)], lda, beta, &c[1 - 1], n1);
                    Cherk("L", "N", n2, k, alpha, &a[((n1 + 1) - 1)], lda, beta, &c[2 - 1], n1);
                    Cgemm("N", "C", n1, n2, k, calpha, &a[(1 - 1)], lda, &a[((n1 + 1) - 1)], lda, cbeta, &c[(n1 * n1 + 1) - 1], n1);
                    //
                } else {
                    //
                    //                 N is odd, TRANSR = 'C', UPLO = 'L', and TRANS = 'C'
                    //
                    Cherk("U", "C", n1, k, alpha, &a[(1 - 1)], lda, beta, &c[1 - 1], n1);
                    Cherk("L", "C", n2, k, alpha, &a[((n1 + 1) - 1) * lda], lda, beta, &c[2 - 1], n1);
                    Cgemm("C", "N", n1, n2, k, calpha, &a[(1 - 1)], lda, &a[((n1 + 1) - 1) * lda], lda, cbeta, &c[(n1 * n1 + 1) - 1], n1);
                    //
                }
                //
            } else {
                //
                //              N is odd, TRANSR = 'C', and UPLO = 'U'
                //
                if (notrans) {
                    //
                    //                 N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'N'
                    //
                    Cherk("U", "N", n1, k, alpha, &a[(1 - 1)], lda, beta, &c[(n2 * n2 + 1) - 1], n2);
                    Cherk("L", "N", n2, k, alpha, &a[((n1 + 1) - 1)], lda, beta, &c[(n1 * n2 + 1) - 1], n2);
                    Cgemm("N", "C", n2, n1, k, calpha, &a[((n1 + 1) - 1)], lda, &a[(1 - 1)], lda, cbeta, &c[1 - 1], n2);
                    //
                } else {
                    //
                    //                 N is odd, TRANSR = 'C', UPLO = 'U', and TRANS = 'C'
                    //
                    Cherk("U", "C", n1, k, alpha, &a[(1 - 1)], lda, beta, &c[(n2 * n2 + 1) - 1], n2);
                    Cherk("L", "C", n2, k, alpha, &a[((n1 + 1) - 1) * lda], lda, beta, &c[(n1 * n2 + 1) - 1], n2);
                    Cgemm("C", "N", n2, n1, k, calpha, &a[((n1 + 1) - 1) * lda], lda, &a[(1 - 1)], lda, cbeta, &c[1 - 1], n2);
                    //
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
                if (notrans) {
                    //
                    //                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N'
                    //
                    Cherk("L", "N", nk, k, alpha, &a[(1 - 1)], lda, beta, &c[2 - 1], n + 1);
                    Cherk("U", "N", nk, k, alpha, &a[((nk + 1) - 1)], lda, beta, &c[1 - 1], n + 1);
                    Cgemm("N", "C", nk, nk, k, calpha, &a[((nk + 1) - 1)], lda, &a[(1 - 1)], lda, cbeta, &c[(nk + 2) - 1], n + 1);
                    //
                } else {
                    //
                    //                 N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'C'
                    //
                    Cherk("L", "C", nk, k, alpha, &a[(1 - 1)], lda, beta, &c[2 - 1], n + 1);
                    Cherk("U", "C", nk, k, alpha, &a[((nk + 1) - 1) * lda], lda, beta, &c[1 - 1], n + 1);
                    Cgemm("C", "N", nk, nk, k, calpha, &a[((nk + 1) - 1) * lda], lda, &a[(1 - 1)], lda, cbeta, &c[(nk + 2) - 1], n + 1);
                    //
                }
                //
            } else {
                //
                //              N is even, TRANSR = 'N', and UPLO = 'U'
                //
                if (notrans) {
                    //
                    //                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'
                    //
                    Cherk("L", "N", nk, k, alpha, &a[(1 - 1)], lda, beta, &c[(nk + 2) - 1], n + 1);
                    Cherk("U", "N", nk, k, alpha, &a[((nk + 1) - 1)], lda, beta, &c[(nk + 1) - 1], n + 1);
                    Cgemm("N", "C", nk, nk, k, calpha, &a[(1 - 1)], lda, &a[((nk + 1) - 1)], lda, cbeta, &c[1 - 1], n + 1);
                    //
                } else {
                    //
                    //                 N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'C'
                    //
                    Cherk("L", "C", nk, k, alpha, &a[(1 - 1)], lda, beta, &c[(nk + 2) - 1], n + 1);
                    Cherk("U", "C", nk, k, alpha, &a[((nk + 1) - 1) * lda], lda, beta, &c[(nk + 1) - 1], n + 1);
                    Cgemm("C", "N", nk, nk, k, calpha, &a[(1 - 1)], lda, &a[((nk + 1) - 1) * lda], lda, cbeta, &c[1 - 1], n + 1);
                    //
                }
                //
            }
            //
        } else {
            //
            //           N is even, and TRANSR = 'C'
            //
            if (lower) {
                //
                //              N is even, TRANSR = 'C', and UPLO = 'L'
                //
                if (notrans) {
                    //
                    //                 N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'N'
                    //
                    Cherk("U", "N", nk, k, alpha, &a[(1 - 1)], lda, beta, &c[(nk + 1) - 1], nk);
                    Cherk("L", "N", nk, k, alpha, &a[((nk + 1) - 1)], lda, beta, &c[1 - 1], nk);
                    Cgemm("N", "C", nk, nk, k, calpha, &a[(1 - 1)], lda, &a[((nk + 1) - 1)], lda, cbeta, &c[(((nk + 1) * nk) + 1) - 1], nk);
                    //
                } else {
                    //
                    //                 N is even, TRANSR = 'C', UPLO = 'L', and TRANS = 'C'
                    //
                    Cherk("U", "C", nk, k, alpha, &a[(1 - 1)], lda, beta, &c[(nk + 1) - 1], nk);
                    Cherk("L", "C", nk, k, alpha, &a[((nk + 1) - 1) * lda], lda, beta, &c[1 - 1], nk);
                    Cgemm("C", "N", nk, nk, k, calpha, &a[(1 - 1)], lda, &a[((nk + 1) - 1) * lda], lda, cbeta, &c[(((nk + 1) * nk) + 1) - 1], nk);
                    //
                }
                //
            } else {
                //
                //              N is even, TRANSR = 'C', and UPLO = 'U'
                //
                if (notrans) {
                    //
                    //                 N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'N'
                    //
                    Cherk("U", "N", nk, k, alpha, &a[(1 - 1)], lda, beta, &c[(nk * (nk + 1) + 1) - 1], nk);
                    Cherk("L", "N", nk, k, alpha, &a[((nk + 1) - 1)], lda, beta, &c[(nk * nk + 1) - 1], nk);
                    Cgemm("N", "C", nk, nk, k, calpha, &a[((nk + 1) - 1)], lda, &a[(1 - 1)], lda, cbeta, &c[1 - 1], nk);
                    //
                } else {
                    //
                    //                 N is even, TRANSR = 'C', UPLO = 'U', and TRANS = 'C'
                    //
                    Cherk("U", "C", nk, k, alpha, &a[(1 - 1)], lda, beta, &c[(nk * (nk + 1) + 1) - 1], nk);
                    Cherk("L", "C", nk, k, alpha, &a[((nk + 1) - 1) * lda], lda, beta, &c[(nk * nk + 1) - 1], nk);
                    Cgemm("C", "N", nk, nk, k, calpha, &a[((nk + 1) - 1) * lda], lda, &a[(1 - 1)], lda, cbeta, &c[1 - 1], nk);
                    //
                }
                //
            }
            //
        }
        //
    }
    //
    //     End of Chfrk
    //
}
