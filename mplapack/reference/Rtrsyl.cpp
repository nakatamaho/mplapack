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

void Rtrsyl(const char *trana, const char *tranb, INTEGER const isgn, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *c, INTEGER const ldc, REAL &scale, INTEGER &info) {
    bool notrna = false;
    bool notrnb = false;
    const REAL one = 1.0;
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    REAL bignum = 0.0;
    REAL dum[1];
    REAL smin = 0.0;
    REAL sgn = 0.0;
    INTEGER lnext = 0;
    INTEGER l = 0;
    INTEGER l1 = 0;
    INTEGER l2 = 0;
    const REAL zero = 0.0;
    INTEGER knext = 0;
    INTEGER k = 0;
    INTEGER k1 = 0;
    INTEGER k2 = 0;
    REAL suml = 0.0;
    REAL sumr = 0.0;
    REAL vec[2 * 2];
    INTEGER ldvec = 2;
    REAL scaloc = 0.0;
    REAL a11 = 0.0;
    REAL da11 = 0.0;
    REAL db = 0.0;
    REAL x[4];
    INTEGER ldx = 2;
    INTEGER j = 0;
    REAL xnorm = 0.0;
    INTEGER ierr = 0;
    //
    //     Decode and Test input parameters
    //
    notrna = Mlsame(trana, "N");
    notrnb = Mlsame(tranb, "N");
    //
    info = 0;
    if (!notrna && !Mlsame(trana, "T") && !Mlsame(trana, "C")) {
        info = -1;
    } else if (!notrnb && !Mlsame(tranb, "T") && !Mlsame(tranb, "C")) {
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
        Mxerbla("Rtrsyl", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    scale = one;
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     Set constants to control overflow
    //
    eps = Rlamch("P");
    smlnum = Rlamch("S");
    bignum = one / smlnum;
    smlnum = smlnum * castREAL(m * n) / eps;
    bignum = one / smlnum;
    //
    smin = max({smlnum, REAL(eps * Rlange("M", m, m, a, lda, dum)), REAL(eps * Rlange("M", n, n, b, ldb, dum))});
    //
    sgn = isgn;
    //
    if (notrna && notrnb) {
        //
        //        Solve    A*X + ISGN*X*B = scale*C.
        //
        //        The (K,L)th block of X is determined starting from
        //        bottom-left corner column by column by
        //
        //         A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
        //
        //        Where
        //                  M                         L-1
        //        R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
        //                I=K+1                       J=1
        //
        //        Start column loop (index = L)
        //        L1 (L2) : column index of the first (first) row of X(K,L).
        //
        lnext = 1;
        for (l = 1; l <= n; l = l + 1) {
            if (l < lnext) {
                goto statement_60;
            }
            if (l == n) {
                l1 = l;
                l2 = l;
            } else {
                if (b[((l + 1) - 1) + (l - 1) * ldb] != zero) {
                    l1 = l;
                    l2 = l + 1;
                    lnext = l + 2;
                } else {
                    l1 = l;
                    l2 = l;
                    lnext = l + 1;
                }
            }
            //
            //           Start row loop (index = K)
            //           K1 (K2): row index of the first (last) row of X(K,L).
            //
            knext = m;
            for (k = m; k >= 1; k = k - 1) {
                if (k > knext) {
                    goto statement_50;
                }
                if (k == 1) {
                    k1 = k;
                    k2 = k;
                } else {
                    if (a[(k - 1) + ((k - 1) - 1) * lda] != zero) {
                        k1 = k - 1;
                        k2 = k;
                        knext = k - 2;
                    } else {
                        k1 = k;
                        k2 = k;
                        knext = k - 1;
                    }
                }
                //
                if (l1 == l2 && k1 == k2) {
                    suml = Rdot(m - k1, &a[(k1 - 1) + (min(k1 + 1, m) - 1) * lda], lda, &c[(min(k1 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[0] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    scaloc = one;
                    //
                    a11 = a[(k1 - 1) + (k1 - 1) * lda] + sgn * b[(l1 - 1) + (l1 - 1) * ldb];
                    da11 = abs(a11);
                    if (da11 <= smin) {
                        a11 = smin;
                        da11 = smin;
                        info = 1;
                    }
                    db = abs(vec[(1 - 1)]);
                    if (da11 < one && db > one) {
                        if (db > bignum * da11) {
                            scaloc = one / db;
                        }
                    }
                    x[(1 - 1)] = (vec[(1 - 1)] * scaloc) / a11;
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    //
                } else if (l1 == l2 && k1 != k2) {
                    //
                    suml = Rdot(m - k2, &a[(k1 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(m - k2, &a[(k2 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k2 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(2 - 1)] = c[(k2 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    Rlaln2(false, 2, 1, smin, one, &a[(k1 - 1) + (k1 - 1) * lda], lda, one, one, vec, 2, -sgn * b[(l1 - 1) + (l1 - 1) * ldb], zero, x, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k2 - 1) + (l1 - 1) * ldc] = x[(2 - 1)];
                    //
                } else if (l1 != l2 && k1 == k2) {
                    //
                    suml = Rdot(m - k1, &a[(k1 - 1) + (min(k1 + 1, m) - 1) * lda], lda, &c[(min(k1 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(1 - 1)] = sgn * (c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr));
                    //
                    suml = Rdot(m - k1, &a[(k1 - 1) + (min(k1 + 1, m) - 1) * lda], lda, &c[(min(k1 + 1, m) - 1) + (l2 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l2 - 1) * ldb], 1);
                    vec[(2 - 1)] = sgn * (c[(k1 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr));
                    //
                    Rlaln2(true, 2, 1, smin, one, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, one, one, vec, 2, -sgn * a[(k1 - 1) + (k1 - 1) * lda], zero, x, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k1 - 1) + (l2 - 1) * ldc] = x[(2 - 1)];
                    //
                } else if (l1 != l2 && k1 != k2) {
                    //
                    suml = Rdot(m - k2, &a[(k1 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(m - k2, &a[(k1 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l2 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l2 - 1) * ldb], 1);
                    vec[(2 - 1) * ldvec] = c[(k1 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(m - k2, &a[(k2 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k2 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(2 - 1)] = c[(k2 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(m - k2, &a[(k2 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l2 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k2 - 1)], ldc, &b[(l2 - 1) * ldb], 1);
                    vec[(2 - 1) + (2 - 1) * ldvec] = c[(k2 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    Rlasy2(false, false, isgn, 2, 2, &a[(k1 - 1) + (k1 - 1) * lda], lda, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, vec, 2, scaloc, x, 2, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k1 - 1) + (l2 - 1) * ldc] = x[(2 - 1) * ldx];
                    c[(k2 - 1) + (l1 - 1) * ldc] = x[(2 - 1)];
                    c[(k2 - 1) + (l2 - 1) * ldc] = x[(2 - 1) + (2 - 1) * ldx];
                }
            //
            statement_50:;
            }
        //
        statement_60:;
        }
        //
    } else if (!notrna && notrnb) {
        //
        //        Solve    A**T *X + ISGN*X*B = scale*C.
        //
        //        The (K,L)th block of X is determined starting from
        //        upper-left corner column by column by
        //
        //          A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
        //
        //        Where
        //                   K-1        T                    L-1
        //          R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
        //                   I=1                          J=1
        //
        //        Start column loop (index = L)
        //        L1 (L2): column index of the first (last) row of X(K,L)
        //
        lnext = 1;
        for (l = 1; l <= n; l = l + 1) {
            if (l < lnext) {
                goto statement_120;
            }
            if (l == n) {
                l1 = l;
                l2 = l;
            } else {
                if (b[((l + 1) - 1) + (l - 1) * ldb] != zero) {
                    l1 = l;
                    l2 = l + 1;
                    lnext = l + 2;
                } else {
                    l1 = l;
                    l2 = l;
                    lnext = l + 1;
                }
            }
            //
            //           Start row loop (index = K)
            //           K1 (K2): row index of the first (last) row of X(K,L)
            //
            knext = 1;
            for (k = 1; k <= m; k = k + 1) {
                if (k < knext) {
                    goto statement_110;
                }
                if (k == m) {
                    k1 = k;
                    k2 = k;
                } else {
                    if (a[((k + 1) - 1) + (k - 1) * lda] != zero) {
                        k1 = k;
                        k2 = k + 1;
                        knext = k + 2;
                    } else {
                        k1 = k;
                        k2 = k;
                        knext = k + 1;
                    }
                }
                //
                if (l1 == l2 && k1 == k2) {
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    scaloc = one;
                    //
                    a11 = a[(k1 - 1) + (k1 - 1) * lda] + sgn * b[(l1 - 1) + (l1 - 1) * ldb];
                    da11 = abs(a11);
                    if (da11 <= smin) {
                        a11 = smin;
                        da11 = smin;
                        info = 1;
                    }
                    db = abs(vec[(1 - 1)]);
                    if (da11 < one && db > one) {
                        if (db > bignum * da11) {
                            scaloc = one / db;
                        }
                    }
                    x[(1 - 1)] = (vec[(1 - 1)] * scaloc) / a11;
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    //
                } else if (l1 == l2 && k1 != k2) {
                    //
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(k1 - 1, &a[(k2 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k2 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(2 - 1)] = c[(k2 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    Rlaln2(true, 2, 1, smin, one, &a[(k1 - 1) + (k1 - 1) * lda], lda, one, one, vec, 2, -sgn * b[(l1 - 1) + (l1 - 1) * ldb], zero, x, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k2 - 1) + (l1 - 1) * ldc] = x[(2 - 1)];
                    //
                } else if (l1 != l2 && k1 == k2) {
                    //
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(1 - 1)] = sgn * (c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr));
                    //
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l2 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l2 - 1) * ldb], 1);
                    vec[(2 - 1)] = sgn * (c[(k1 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr));
                    //
                    Rlaln2(true, 2, 1, smin, one, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, one, one, vec, 2, -sgn * a[(k1 - 1) + (k1 - 1) * lda], zero, x, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k1 - 1) + (l2 - 1) * ldc] = x[(2 - 1)];
                    //
                } else if (l1 != l2 && k1 != k2) {
                    //
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l2 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k1 - 1)], ldc, &b[(l2 - 1) * ldb], 1);
                    vec[(2 - 1) * ldvec] = c[(k1 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(k1 - 1, &a[(k2 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k2 - 1)], ldc, &b[(l1 - 1) * ldb], 1);
                    vec[(2 - 1)] = c[(k2 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(k1 - 1, &a[(k2 - 1) * lda], 1, &c[(l2 - 1) * ldc], 1);
                    sumr = Rdot(l1 - 1, &c[(k2 - 1)], ldc, &b[(l2 - 1) * ldb], 1);
                    vec[(2 - 1) + (2 - 1) * ldvec] = c[(k2 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    Rlasy2(true, false, isgn, 2, 2, &a[(k1 - 1) + (k1 - 1) * lda], lda, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, vec, 2, scaloc, x, 2, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k1 - 1) + (l2 - 1) * ldc] = x[(2 - 1) * ldx];
                    c[(k2 - 1) + (l1 - 1) * ldc] = x[(2 - 1)];
                    c[(k2 - 1) + (l2 - 1) * ldc] = x[(2 - 1) + (2 - 1) * ldx];
                }
            //
            statement_110:;
            }
        statement_120:;
        }
        //
    } else if (!notrna && !notrnb) {
        //
        //        Solve    A**T*X + ISGN*X*B**T = scale*C.
        //
        //        The (K,L)th block of X is determined starting from
        //        top-right corner column by column by
        //
        //           A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)
        //
        //        Where
        //                     K-1                            N
        //            R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
        //                     I=1                          J=L+1
        //
        //        Start column loop (index = L)
        //        L1 (L2): column index of the first (last) row of X(K,L)
        //
        lnext = n;
        for (l = n; l >= 1; l = l - 1) {
            if (l > lnext) {
                goto statement_180;
            }
            if (l == 1) {
                l1 = l;
                l2 = l;
            } else {
                if (b[(l - 1) + ((l - 1) - 1) * ldb] != zero) {
                    l1 = l - 1;
                    l2 = l;
                    lnext = l - 2;
                } else {
                    l1 = l;
                    l2 = l;
                    lnext = l - 1;
                }
            }
            //
            //           Start row loop (index = K)
            //           K1 (K2): row index of the first (last) row of X(K,L)
            //
            knext = 1;
            for (k = 1; k <= m; k = k + 1) {
                if (k < knext) {
                    goto statement_170;
                }
                if (k == m) {
                    k1 = k;
                    k2 = k;
                } else {
                    if (a[((k + 1) - 1) + (k - 1) * lda] != zero) {
                        k1 = k;
                        k2 = k + 1;
                        knext = k + 2;
                    } else {
                        k1 = k;
                        k2 = k;
                        knext = k + 1;
                    }
                }
                //
                if (l1 == l2 && k1 == k2) {
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l1, &c[(k1 - 1) + (min(l1 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l1 + 1, n) - 1) * ldb], ldb);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    scaloc = one;
                    //
                    a11 = a[(k1 - 1) + (k1 - 1) * lda] + sgn * b[(l1 - 1) + (l1 - 1) * ldb];
                    da11 = abs(a11);
                    if (da11 <= smin) {
                        a11 = smin;
                        da11 = smin;
                        info = 1;
                    }
                    db = abs(vec[(1 - 1)]);
                    if (da11 < one && db > one) {
                        if (db > bignum * da11) {
                            scaloc = one / db;
                        }
                    }
                    x[(1 - 1)] = (vec[(1 - 1)] * scaloc) / a11;
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    //
                } else if (l1 == l2 && k1 != k2) {
                    //
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k1 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(k1 - 1, &a[(k2 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k2 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(2 - 1)] = c[(k2 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    Rlaln2(true, 2, 1, smin, one, &a[(k1 - 1) + (k1 - 1) * lda], lda, one, one, vec, 2, -sgn * b[(l1 - 1) + (l1 - 1) * ldb], zero, x, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k2 - 1) + (l1 - 1) * ldc] = x[(2 - 1)];
                    //
                } else if (l1 != l2 && k1 == k2) {
                    //
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k1 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(1 - 1)] = sgn * (c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr));
                    //
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l2 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k1 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l2 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(2 - 1)] = sgn * (c[(k1 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr));
                    //
                    Rlaln2(false, 2, 1, smin, one, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, one, one, vec, 2, -sgn * a[(k1 - 1) + (k1 - 1) * lda], zero, x, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k1 - 1) + (l2 - 1) * ldc] = x[(2 - 1)];
                    //
                } else if (l1 != l2 && k1 != k2) {
                    //
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k1 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(k1 - 1, &a[(k1 - 1) * lda], 1, &c[(l2 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k1 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l2 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(2 - 1) * ldvec] = c[(k1 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(k1 - 1, &a[(k2 - 1) * lda], 1, &c[(l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k2 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(2 - 1)] = c[(k2 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(k1 - 1, &a[(k2 - 1) * lda], 1, &c[(l2 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k2 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l2 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(2 - 1) + (2 - 1) * ldvec] = c[(k2 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    Rlasy2(true, true, isgn, 2, 2, &a[(k1 - 1) + (k1 - 1) * lda], lda, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, vec, 2, scaloc, x, 2, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k1 - 1) + (l2 - 1) * ldc] = x[(2 - 1) * ldx];
                    c[(k2 - 1) + (l1 - 1) * ldc] = x[(2 - 1)];
                    c[(k2 - 1) + (l2 - 1) * ldc] = x[(2 - 1) + (2 - 1) * ldx];
                }
            //
            statement_170:;
            }
        statement_180:;
        }
        //
    } else if (notrna && !notrnb) {
        //
        //        Solve    A*X + ISGN*X*B**T = scale*C.
        //
        //        The (K,L)th block of X is determined starting from
        //        bottom-right corner column by column by
        //
        //            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)
        //
        //        Where
        //                      M                          N
        //            R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
        //                    I=K+1                      J=L+1
        //
        //        Start column loop (index = L)
        //        L1 (L2): column index of the first (last) row of X(K,L)
        //
        lnext = n;
        for (l = n; l >= 1; l = l - 1) {
            if (l > lnext) {
                goto statement_240;
            }
            if (l == 1) {
                l1 = l;
                l2 = l;
            } else {
                if (b[(l - 1) + ((l - 1) - 1) * ldb] != zero) {
                    l1 = l - 1;
                    l2 = l;
                    lnext = l - 2;
                } else {
                    l1 = l;
                    l2 = l;
                    lnext = l - 1;
                }
            }
            //
            //           Start row loop (index = K)
            //           K1 (K2): row index of the first (last) row of X(K,L)
            //
            knext = m;
            for (k = m; k >= 1; k = k - 1) {
                if (k > knext) {
                    goto statement_230;
                }
                if (k == 1) {
                    k1 = k;
                    k2 = k;
                } else {
                    if (a[(k - 1) + ((k - 1) - 1) * lda] != zero) {
                        k1 = k - 1;
                        k2 = k;
                        knext = k - 2;
                    } else {
                        k1 = k;
                        k2 = k;
                        knext = k - 1;
                    }
                }
                //
                if (l1 == l2 && k1 == k2) {
                    suml = Rdot(m - k1, &a[(k1 - 1) + (min(k1 + 1, m) - 1) * lda], lda, &c[(min(k1 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l1, &c[(k1 - 1) + (min(l1 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l1 + 1, n) - 1) * ldb], ldb);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    scaloc = one;
                    //
                    a11 = a[(k1 - 1) + (k1 - 1) * lda] + sgn * b[(l1 - 1) + (l1 - 1) * ldb];
                    da11 = abs(a11);
                    if (da11 <= smin) {
                        a11 = smin;
                        da11 = smin;
                        info = 1;
                    }
                    db = abs(vec[(1 - 1)]);
                    if (da11 < one && db > one) {
                        if (db > bignum * da11) {
                            scaloc = one / db;
                        }
                    }
                    x[(1 - 1)] = (vec[(1 - 1)] * scaloc) / a11;
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    //
                } else if (l1 == l2 && k1 != k2) {
                    //
                    suml = Rdot(m - k2, &a[(k1 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k1 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(m - k2, &a[(k2 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k2 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(2 - 1)] = c[(k2 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    Rlaln2(false, 2, 1, smin, one, &a[(k1 - 1) + (k1 - 1) * lda], lda, one, one, vec, 2, -sgn * b[(l1 - 1) + (l1 - 1) * ldb], zero, x, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k2 - 1) + (l1 - 1) * ldc] = x[(2 - 1)];
                    //
                } else if (l1 != l2 && k1 == k2) {
                    //
                    suml = Rdot(m - k1, &a[(k1 - 1) + (min(k1 + 1, m) - 1) * lda], lda, &c[(min(k1 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k1 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(1 - 1)] = sgn * (c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr));
                    //
                    suml = Rdot(m - k1, &a[(k1 - 1) + (min(k1 + 1, m) - 1) * lda], lda, &c[(min(k1 + 1, m) - 1) + (l2 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k1 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l2 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(2 - 1)] = sgn * (c[(k1 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr));
                    //
                    Rlaln2(false, 2, 1, smin, one, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, one, one, vec, 2, -sgn * a[(k1 - 1) + (k1 - 1) * lda], zero, x, 2, scaloc, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k1 - 1) + (l2 - 1) * ldc] = x[(2 - 1)];
                    //
                } else if (l1 != l2 && k1 != k2) {
                    //
                    suml = Rdot(m - k2, &a[(k1 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k1 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(1 - 1)] = c[(k1 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(m - k2, &a[(k1 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l2 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k1 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l2 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(2 - 1) * ldvec] = c[(k1 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(m - k2, &a[(k2 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l1 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k2 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l1 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(2 - 1)] = c[(k2 - 1) + (l1 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    suml = Rdot(m - k2, &a[(k2 - 1) + (min(k2 + 1, m) - 1) * lda], lda, &c[(min(k2 + 1, m) - 1) + (l2 - 1) * ldc], 1);
                    sumr = Rdot(n - l2, &c[(k2 - 1) + (min(l2 + 1, n) - 1) * ldc], ldc, &b[(l2 - 1) + (min(l2 + 1, n) - 1) * ldb], ldb);
                    vec[(2 - 1) + (2 - 1) * ldvec] = c[(k2 - 1) + (l2 - 1) * ldc] - (suml + sgn * sumr);
                    //
                    Rlasy2(false, true, isgn, 2, 2, &a[(k1 - 1) + (k1 - 1) * lda], lda, &b[(l1 - 1) + (l1 - 1) * ldb], ldb, vec, 2, scaloc, x, 2, xnorm, ierr);
                    if (ierr != 0) {
                        info = 1;
                    }
                    //
                    if (scaloc != one) {
                        for (j = 1; j <= n; j = j + 1) {
                            Rscal(m, scaloc, &c[(j - 1) * ldc], 1);
                        }
                        scale = scale * scaloc;
                    }
                    c[(k1 - 1) + (l1 - 1) * ldc] = x[(1 - 1)];
                    c[(k1 - 1) + (l2 - 1) * ldc] = x[(2 - 1) * ldx];
                    c[(k2 - 1) + (l1 - 1) * ldc] = x[(2 - 1)];
                    c[(k2 - 1) + (l2 - 1) * ldc] = x[(2 - 1) + (2 - 1) * ldx];
                }
            //
            statement_230:;
            }
        statement_240:;
        }
        //
    }
    //
    //     End of Rtrsyl
    //
}
