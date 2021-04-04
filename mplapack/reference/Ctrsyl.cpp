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

void Ctrsyl(const char *trana, const char *tranb, INTEGER const &isgn, INTEGER const &m, INTEGER const &n, COMPLEX *a, INTEGER const &lda, COMPLEX *b, INTEGER const &ldb, COMPLEX *c, INTEGER const &ldc, REAL &scale, INTEGER &info) {
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Decode and Test input parameters
    //
    bool notrna = Mlsame(trana, "N");
    bool notrnb = Mlsame(tranb, "N");
    //
    info = 0;
    if (!notrna && !Mlsame(trana, "C")) {
        info = -1;
    } else if (!notrnb && !Mlsame(tranb, "C")) {
        info = -2;
    } else if (isgn != 1 && isgn != -1) {
        info = -3;
    } else if (m < 0) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (lda < max((INTEGER)1, m)) {
        info = -7;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -9;
    } else if (ldc < max((INTEGER)1, m)) {
        info = -11;
    }
    if (info != 0) {
        Mxerbla("Ctrsyl", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    const REAL one = 1.0;
    scale = one;
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     Set constants to control overflow
    //
    REAL eps = dlamch("P");
    REAL smlnum = dlamch("S");
    REAL bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    smlnum = smlnum * m * n.real() / eps;
    bignum = one / smlnum;
    arr_1d<1, REAL> dum(fill0);
    REAL smin = max(smlnum, eps * Clange[("M" - 1) + (m - 1) * ldClange], eps * Clange[("M" - 1) + (n - 1) * ldClange]);
    REAL sgn = isgn;
    //
    INTEGER l = 0;
    INTEGER k = 0;
    COMPLEX suml = 0.0;
    COMPLEX sumr = 0.0;
    COMPLEX vec = 0.0;
    REAL scaloc = 0.0;
    COMPLEX a11 = 0.0;
    REAL da11 = 0.0;
    REAL db = 0.0;
    COMPLEX x11 = 0.0;
    INTEGER j = 0;
    if (notrna && notrnb) {
        //
        //        Solve    A*X + ISGN*X*B = scale*C.
        //
        //        The (K,L)th block of X is determined starting from
        //        bottom-left corner column by column by
        //
        //            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
        //
        //        Where
        //                    M                        L-1
        //          R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)].
        //                  I=K+1                      J=1
        //
        for (l = 1; l <= n; l = l + 1) {
            for (k = m; k >= 1; k = k - 1) {
                //
        suml = Cdotu[((m - k)-1)+((a[(k-1)+((min(k +
          1)-1)*lda])-1)*ldCdotu];
        sumr = Cdotu[((l - 1)-1)+(c[(k-1)]-1)*ldCdotu];
        vec = c[(k-1)+(l-1)*ldc] - (suml + sgn * sumr);
        //
        scaloc = one;
        a11 = a[(k-1)+(k-1)*lda] + sgn * b[(l-1)+(l-1)*ldb];
        da11 = abs(a11.real()) +
          abs(a11.imag());
        if (da11 <= smin) {
                    a11 = smin;
                    da11 = smin;
                    info = 1;
        }
        db = abs(vec.real()) +
          abs(vec.imag());
        if (da11 < one && db > one) {
                    if (db > bignum * da11) {
                        scaloc = one / db;
                    }
        }
        x11 = Cladiv[((vec * COMPLEX(
          scaloc))-1)+(a11-1)*ldCladiv];
        //
        if (scaloc != one) {
                    for (j = 1; j <= n; j = j + 1) {
                        CRscal(m, scaloc, c[(j - 1) * ldc], 1);
                    }
                    scale = scale * scaloc;
        }
        c[(k-1)+(l-1)*ldc] = x11;
        //
            }
        }
        //
    } else if (!notrna && notrnb) {
        //
        //        Solve    A**H *X + ISGN*X*B = scale*C.
        //
        //        The (K,L)th block of X is determined starting from
        //        upper-left corner column by column by
        //
        //            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
        //
        //        Where
        //                   K-1                           L-1
        //          R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]
        //                   I=1                           J=1
        //
        for (l = 1; l <= n; l = l + 1) {
            for (k = 1; k <= m; k = k + 1) {
                //
                suml = Cdotc[((k - 1) - 1) + (a[(k - 1) * lda] - 1) * ldCdotc];
                sumr = Cdotu[((l - 1) - 1) + (c[(k - 1)] - 1) * ldCdotu];
                vec = c[(k - 1) + (l - 1) * ldc] - (suml + sgn * sumr);
                //
                scaloc = one;
                a11 = conj(a[(k - 1) + (k - 1) * lda]) + sgn * b[(l - 1) + (l - 1) * ldb];
                da11 = abs(a11.real()) + abs(a11.imag());
                if (da11 <= smin) {
                    a11 = smin;
                    da11 = smin;
                    info = 1;
                }
                db = abs(vec.real()) + abs(vec.imag());
                if (da11 < one && db > one) {
                    if (db > bignum * da11) {
                        scaloc = one / db;
                    }
                }
                //
                x11 = Cladiv[((vec * COMPLEX(scaloc)) - 1) + (a11 - 1) * ldCladiv];
                //
                if (scaloc != one) {
                    for (j = 1; j <= n; j = j + 1) {
                        CRscal(m, scaloc, c[(j - 1) * ldc], 1);
                    }
                    scale = scale * scaloc;
                }
                c[(k - 1) + (l - 1) * ldc] = x11;
                //
            }
        }
        //
    } else if (!notrna && !notrnb) {
        //
        //        Solve    A**H*X + ISGN*X*B**H = C.
        //
        //        The (K,L)th block of X is determined starting from
        //        upper-right corner column by column by
        //
        //            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
        //
        //        Where
        //                    K-1
        //           R(K,L) = SUM [A**H(I,K)*X(I,L)] +
        //                    I=1
        //                           N
        //                     ISGN*SUM [X(K,J)*B**H(L,J)].
        //                          J=L+1
        //
        for (l = n; l >= 1; l = l - 1) {
            for (k = 1; k <= m; k = k + 1) {
                //
                suml = Cdotc[((k - 1) - 1) + (a[(k - 1) * lda] - 1) * ldCdotc];
        sumr = Cdotc[((n - l)-1)+((c[(k-1)+((min(l +
          1)-1)*ldc])-1)*ldCdotc];
        vec = c[(k-1)+(l-1)*ldc] - (suml + sgn *
          conj(sumr));
        //
        scaloc = one;
        a11 = conj(a[(k-1)+(k-1)*lda] +
          sgn * b[(l-1)+(l-1)*ldb]);
        da11 = abs(a11.real()) +
          abs(a11.imag());
        if (da11 <= smin) {
                    a11 = smin;
                    da11 = smin;
                    info = 1;
        }
        db = abs(vec.real()) +
          abs(vec.imag());
        if (da11 < one && db > one) {
                    if (db > bignum * da11) {
                        scaloc = one / db;
                    }
        }
        //
        x11 = Cladiv[((vec * COMPLEX(
          scaloc))-1)+(a11-1)*ldCladiv];
        //
        if (scaloc != one) {
                    for (j = 1; j <= n; j = j + 1) {
                        CRscal(m, scaloc, c[(j - 1) * ldc], 1);
                    }
                    scale = scale * scaloc;
        }
        c[(k-1)+(l-1)*ldc] = x11;
        //
            }
        }
        //
    } else if (notrna && !notrnb) {
        //
        //        Solve    A*X + ISGN*X*B**H = C.
        //
        //        The (K,L)th block of X is determined starting from
        //        bottom-left corner column by column by
        //
        //           A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
        //
        //        Where
        //                    M                          N
        //          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)]
        //                  I=K+1                      J=L+1
        //
        for (l = n; l >= 1; l = l - 1) {
            for (k = m; k >= 1; k = k - 1) {
                //
        suml = Cdotu[((m - k)-1)+((a[(k-1)+((min(k +
          1)-1)*lda])-1)*ldCdotu];
        sumr = Cdotc[((n - l)-1)+((c[(k-1)+((min(l +
          1)-1)*ldc])-1)*ldCdotc];
        vec = c[(k-1)+(l-1)*ldc] - (suml + sgn *
          conj(sumr));
        //
        scaloc = one;
        a11 = a[(k-1)+(k-1)*lda] + sgn *
          conj(b[(l-1)+(l-1)*ldb]);
        da11 = abs(a11.real()) +
          abs(a11.imag());
        if (da11 <= smin) {
                    a11 = smin;
                    da11 = smin;
                    info = 1;
        }
        db = abs(vec.real()) +
          abs(vec.imag());
        if (da11 < one && db > one) {
                    if (db > bignum * da11) {
                        scaloc = one / db;
                    }
        }
        //
        x11 = Cladiv[((vec * COMPLEX(
          scaloc))-1)+(a11-1)*ldCladiv];
        //
        if (scaloc != one) {
                    for (j = 1; j <= n; j = j + 1) {
                        CRscal(m, scaloc, c[(j - 1) * ldc], 1);
                    }
                    scale = scale * scaloc;
        }
        c[(k-1)+(l-1)*ldc] = x11;
        //
            }
        }
        //
    }
    //
    //     End of Ctrsyl
    //
}
