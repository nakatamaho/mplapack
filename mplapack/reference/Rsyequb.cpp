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

void Rsyequb(const char *uplo, INTEGER const &n, REAL *a, INTEGER const &lda, REAL *s, REAL &scond, REAL &amax, REAL *work, INTEGER &info) {
    bool up = false;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    INTEGER i = 0;
    INTEGER j = 0;
    REAL tol = 0.0;
    INTEGER iter = 0;
    const INTEGER max_iter = 100;
    REAL scale = 0.0;
    REAL sumsq = 0.0;
    REAL avg = 0.0;
    REAL std = 0.0;
    REAL t = 0.0;
    REAL si = 0.0;
    REAL c2 = 0.0;
    REAL c1 = 0.0;
    REAL c0 = 0.0;
    REAL d = 0.0;
    REAL u = 0.0;
    REAL smlnum = 0.0;
    REAL bignum = 0.0;
    REAL smin = 0.0;
    REAL smax = 0.0;
    REAL base = 0.0;
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
    if (!(Mlsame(uplo, "U") || Mlsame(uplo, "L"))) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Rsyequb", -info);
        return;
    }
    //
    up = Mlsame(uplo, "U");
    amax = zero;
    //
    //     Quick return if possible.
    //
    if (n == 0) {
        scond = one;
        return;
    }
    //
    for (i = 1; i <= n; i = i + 1) {
        s[i - 1] = zero;
    }
    //
    amax = zero;
    if (up) {
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= j - 1; i = i + 1) {
                s[i - 1] = max(s[i - 1], abs(a[(i - 1) + (j - 1) * lda]));
                s[j - 1] = max(s[j - 1], abs(a[(i - 1) + (j - 1) * lda]));
                amax = max(amax, abs(a[(i - 1) + (j - 1) * lda]));
            }
            s[j - 1] = max(s[j - 1], abs(a[(j - 1) + (j - 1) * lda]));
            amax = max(amax, abs(a[(j - 1) + (j - 1) * lda]));
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            s[j - 1] = max(s[j - 1], abs(a[(j - 1) + (j - 1) * lda]));
            amax = max(amax, abs(a[(j - 1) + (j - 1) * lda]));
            for (i = j + 1; i <= n; i = i + 1) {
                s[i - 1] = max(s[i - 1], abs(a[(i - 1) + (j - 1) * lda]));
                s[j - 1] = max(s[j - 1], abs(a[(i - 1) + (j - 1) * lda]));
                amax = max(amax, abs(a[(i - 1) + (j - 1) * lda]));
            }
        }
    }
    for (j = 1; j <= n; j = j + 1) {
        s[j - 1] = 1.0 / s[j - 1];
    }
    //
    tol = one / sqrt(2.0 * n);
    //
    for (iter = 1; iter <= max_iter; iter = iter + 1) {
        scale = 0.0;
        sumsq = 0.0;
        //        beta = |A|s
        for (i = 1; i <= n; i = i + 1) {
            work[i - 1] = zero;
        }
        if (up) {
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j - 1; i = i + 1) {
                    work[i - 1] += abs(a[(i - 1) + (j - 1) * lda]) * s[j - 1];
                    work[j - 1] += abs(a[(i - 1) + (j - 1) * lda]) * s[i - 1];
                }
                work[j - 1] += abs(a[(j - 1) + (j - 1) * lda]) * s[j - 1];
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                work[j - 1] += abs(a[(j - 1) + (j - 1) * lda]) * s[j - 1];
                for (i = j + 1; i <= n; i = i + 1) {
                    work[i - 1] += abs(a[(i - 1) + (j - 1) * lda]) * s[j - 1];
                    work[j - 1] += abs(a[(i - 1) + (j - 1) * lda]) * s[i - 1];
                }
            }
        }
        //
        //        avg = s^T beta / n
        avg = 0.0;
        for (i = 1; i <= n; i = i + 1) {
            avg += s[i - 1] * work[i - 1];
        }
        avg = avg / n;
        //
        std = 0.0;
        for (i = n + 1; i <= 2 * n; i = i + 1) {
            work[i - 1] = s[(i - n) - 1] * work[(i - n) - 1] - avg;
        }
        Rlassq(n, work[(n + 1) - 1], 1, scale, sumsq);
        std = scale * sqrt(sumsq / n);
        //
        if (std < tol * avg) {
            goto statement_999;
        }
        //
        for (i = 1; i <= n; i = i + 1) {
            t = abs(a[(i - 1) + (i - 1) * lda]);
            si = s[i - 1];
            c2 = (n - 1) * t;
            c1 = (n - 2) * (work[i - 1] - t * si);
            c0 = -(t * si) * si + 2 * work[i - 1] * si - n * avg;
            d = c1 * c1 - 4 * c0 * c2;
            //
            if (d <= 0) {
                info = -1;
                return;
            }
            si = -2 * c0 / (c1 + sqrt(d));
            //
            d = si - s[i - 1];
            u = zero;
            if (up) {
                for (j = 1; j <= i; j = j + 1) {
                    t = abs(a[(j - 1) + (i - 1) * lda]);
                    u += s[j - 1] * t;
                    work[j - 1] += d * t;
                }
                for (j = i + 1; j <= n; j = j + 1) {
                    t = abs(a[(i - 1) + (j - 1) * lda]);
                    u += s[j - 1] * t;
                    work[j - 1] += d * t;
                }
            } else {
                for (j = 1; j <= i; j = j + 1) {
                    t = abs(a[(i - 1) + (j - 1) * lda]);
                    u += s[j - 1] * t;
                    work[j - 1] += d * t;
                }
                for (j = i + 1; j <= n; j = j + 1) {
                    t = abs(a[(j - 1) + (i - 1) * lda]);
                    u += s[j - 1] * t;
                    work[j - 1] += d * t;
                }
            }
            //
            avg += (u + work[i - 1]) * d / n;
            s[i - 1] = si;
        }
    }
//
statement_999:
    //
    smlnum = dlamch("SAFEMIN");
    bignum = one / smlnum;
    smin = bignum;
    smax = zero;
    t = one / sqrt(avg);
    base = dlamch("B");
    u = one / log[base - 1];
    for (i = 1; i <= n; i = i + 1) {
        s[i - 1] = pow(base, INTEGER(u * log[(s[i - 1] * t) - 1]));
        smin = min(smin, s[i - 1]);
        smax = max(smax, s[i - 1]);
    }
    scond = max(smin, smlnum) / min(smax, bignum);
    //
}
