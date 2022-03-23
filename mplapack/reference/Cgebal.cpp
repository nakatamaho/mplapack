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

void Cgebal(const char *job, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER &ilo, INTEGER &ihi, REAL *scale, INTEGER &info) {
    INTEGER k = 0;
    INTEGER l = 0;
    INTEGER i = 0;
    const REAL one = 1.0;
    INTEGER m = 0;
    INTEGER j = 0;
    INTEGER iexc = 0;
    const REAL zero = 0.0;
    REAL sfmin1 = 0.0;
    REAL sfmax1 = 0.0;
    const REAL sclfac = 2.0e+0;
    REAL sfmin2 = 0.0;
    REAL sfmax2 = 0.0;
    bool noconv = false;
    REAL c = 0.0;
    REAL r = 0.0;
    INTEGER ica = 0;
    REAL ca = 0.0;
    INTEGER ira = 0;
    REAL ra = 0.0;
    REAL g = 0.0;
    REAL f = 0.0;
    REAL s = 0.0;
    const REAL factor = 0.95e+0;
    //
    //     Test the input parameters
    //
    info = 0;
    if (!Mlsame(job, "N") && !Mlsame(job, "P") && !Mlsame(job, "S") && !Mlsame(job, "B")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (lda < max((INTEGER)1, n)) {
        info = -4;
    }
    if (info != 0) {
        Mxerbla("Cgebal", -info);
        return;
    }
    //
    k = 1;
    l = n;
    //
    if (n == 0) {
        goto statement_210;
    }
    //
    if (Mlsame(job, "N")) {
        for (i = 1; i <= n; i = i + 1) {
            scale[i - 1] = one;
        }
        goto statement_210;
    }
    //
    if (Mlsame(job, "S")) {
        goto statement_120;
    }
    //
    //     Permutation to isolate eigenvalues if possible
    //
    goto statement_50;
//
//     Row and column exchange.
//
statement_20:
    scale[m - 1] = j;
    if (j == m) {
        goto statement_30;
    }
    //
    Cswap(l, &a[(j - 1) * lda], 1, &a[(m - 1) * lda], 1);
    Cswap(n - k + 1, &a[(j - 1) + (k - 1) * lda], lda, &a[(m - 1) + (k - 1) * lda], lda);
//
statement_30:
    switch (iexc) {
    case 1:
        goto statement_40;
    case 2:
        goto statement_80;
    default:
        break;
    }
//
//     Search for rows isolating an eigenvalue and push them down.
//
statement_40:
    if (l == 1) {
        goto statement_210;
    }
    l = l - 1;
//
statement_50:
    for (j = l; j >= 1; j = j - 1) {
        //
        for (i = 1; i <= l; i = i + 1) {
            if (i == j) {
                goto statement_60;
            }
            if (a[(j - 1) + (i - 1) * lda].real() != zero || a[(j - 1) + (i - 1) * lda].imag() != zero) {
                goto statement_70;
            }
        statement_60:;
        }
        //
        m = l;
        iexc = 1;
        goto statement_20;
    statement_70:;
    }
    //
    goto statement_90;
//
//     Search for columns isolating an eigenvalue and push them left.
//
statement_80:
    k++;
//
statement_90:
    for (j = k; j <= l; j = j + 1) {
        //
        for (i = k; i <= l; i = i + 1) {
            if (i == j) {
                goto statement_100;
            }
            if (a[(i - 1) + (j - 1) * lda].real() != zero || a[(i - 1) + (j - 1) * lda].imag() != zero) {
                goto statement_110;
            }
        statement_100:;
        }
        //
        m = k;
        iexc = 2;
        goto statement_20;
    statement_110:;
    }
//
statement_120:
    for (i = k; i <= l; i = i + 1) {
        scale[i - 1] = one;
    }
    //
    if (Mlsame(job, "P")) {
        goto statement_210;
    }
    //
    //     Balance the submatrix in rows K to L.
    //
    //     Iterative loop for norm reduction
    //
    sfmin1 = Rlamch("S") / Rlamch("P");
    sfmax1 = one / sfmin1;
    sfmin2 = sfmin1 * sclfac;
    sfmax2 = one / sfmin2;
statement_140:
    noconv = false;
    //
    for (i = k; i <= l; i = i + 1) {
        //
        c = RCnrm2(l - k + 1, &a[(k - 1) + (i - 1) * lda], 1);
        r = RCnrm2(l - k + 1, &a[(i - 1) + (k - 1) * lda], lda);
        ica = iCamax(l, &a[(i - 1) * lda], 1);
        ca = abs(a[(ica - 1) + (i - 1) * lda]);
        ira = iCamax(n - k + 1, &a[(i - 1) + (k - 1) * lda], lda);
        ra = abs(a[(i - 1) + ((ira + k - 1) - 1) * lda]);
        //
        //        Guard against zero C or R due to underflow.
        //
        if (c == zero || r == zero) {
            goto statement_200;
        }
        g = r / sclfac;
        f = one;
        s = c + r;
    statement_160:
        if (c >= g || max({f, c, ca}) >= sfmax2 || min({r, g, ra}) <= sfmin2) {
            goto statement_170;
        }
        if (Risnan(c + f + ca + r + g + ra)) {
            //
            //           Exit if NaN to avoid infinite loop
            //
            info = -3;
            Mxerbla("Cgebal", -info);
            return;
        }
        f = f * sclfac;
        c = c * sclfac;
        ca = ca * sclfac;
        r = r / sclfac;
        g = g / sclfac;
        ra = ra / sclfac;
        goto statement_160;
    //
    statement_170:
        g = c / sclfac;
    statement_180:
        if (g < r || max(r, ra) >= sfmax2 || min({f, c, g, ca}) <= sfmin2) {
            goto statement_190;
        }
        f = f / sclfac;
        c = c / sclfac;
        g = g / sclfac;
        ca = ca / sclfac;
        r = r * sclfac;
        ra = ra * sclfac;
        goto statement_180;
    //
    //        Now balance.
    //
    statement_190:
        if ((c + r) >= factor * s) {
            goto statement_200;
        }
        if (f < one && scale[i - 1] < one) {
            if (f * scale[i - 1] <= sfmin1) {
                goto statement_200;
            }
        }
        if (f > one && scale[i - 1] > one) {
            if (scale[i - 1] >= sfmax1 / f) {
                goto statement_200;
            }
        }
        g = one / f;
        scale[i - 1] = scale[i - 1] * f;
        noconv = true;
        //
        CRscal(n - k + 1, g, &a[(i - 1) + (k - 1) * lda], lda);
        CRscal(l, f, &a[(i - 1) * lda], 1);
    //
    statement_200:;
    }
    //
    if (noconv) {
        goto statement_140;
    }
//
statement_210:
    ilo = k;
    ihi = l;
    //
    //     End of Cgebal
    //
}
