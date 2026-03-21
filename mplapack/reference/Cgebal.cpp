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

// Derived from LAPACK routine ZGEBAL.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Cgebal(const char *job, INTEGER const n, COMPLEX *a, INTEGER const lda, INTEGER &ilo, INTEGER &ihi, REAL *scale, INTEGER &info) {
    //
    // Test the input parameters
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
    // Quick returns.
    //
    if (n == 0) {
        ilo = 1;
        ihi = 0;
        return;
    }
    //
    INTEGER i = 0;
    const REAL one = 1.0;
    if (Mlsame(job, "N")) {
        for (i = 1; i <= n; i = i + 1) {
            scale[i - 1] = one;
        }
        ilo = 1;
        ihi = n;
        return;
    }
    //
    // Permutation to isolate eigenvalues if possible.
    //
    INTEGER k = 1;
    INTEGER l = n;
    //
    bool noconv = false;
    bool canswap = false;
    INTEGER j = 0;
    const REAL zero = 0.0;
    if (!Mlsame(job, "S")) {
        //
        // Row and column exchange.
        //
        noconv = true;
        while (noconv) {
            //
            // Search for rows isolating an eigenvalue and push them down.
            //
            noconv = false;
            for (i = l; i >= 1; i = i - 1) {
                canswap = true;
                for (j = 1; j <= l; j = j + 1) {
                    if (i != j && (a[(i - 1) + (j - 1) * lda].real() != zero || a[(i - 1) + (j - 1) * lda].imag() != zero)) {
                        canswap = false;
                        break;
                    }
                }
                //
                if (canswap) {
                    scale[l - 1] = i;
                    if (i != l) {
                        Cswap(l, &a[(i - 1) * lda], 1, &a[(l - 1) * lda], 1);
                        Cswap(n - k + 1, &a[(i - 1) + (k - 1) * lda], lda, &a[(l - 1) + (k - 1) * lda], lda);
                    }
                    noconv = true;
                    //
                    if (l == 1) {
                        ilo = 1;
                        ihi = 1;
                        return;
                    }
                    //
                    l = l - 1;
                }
            }
            //
        }
        //
        noconv = true;
        while (noconv) {
            //
            // Search for columns isolating an eigenvalue and push them left.
            //
            noconv = false;
            for (j = k; j <= l; j = j + 1) {
                canswap = true;
                for (i = k; i <= l; i = i + 1) {
                    if (i != j && (a[(i - 1) + (j - 1) * lda].real() != zero || a[(i - 1) + (j - 1) * lda].imag() != zero)) {
                        canswap = false;
                        break;
                    }
                }
                //
                if (canswap) {
                    scale[k - 1] = j;
                    if (j != k) {
                        Cswap(l, &a[(j - 1) * lda], 1, &a[(k - 1) * lda], 1);
                        Cswap(n - k + 1, &a[(j - 1) + (k - 1) * lda], lda, &a[(k - 1) + (k - 1) * lda], lda);
                    }
                    noconv = true;
                    //
                    k++;
                }
            }
            //
        }
        //
    }
    //
    // Initialize SCALE for non-permuted submatrix.
    //
    for (i = k; i <= l; i = i + 1) {
        scale[i - 1] = one;
    }
    //
    // If we only had to permute, we are done.
    //
    if (Mlsame(job, "P")) {
        ilo = k;
        ihi = l;
        return;
    }
    //
    // Balance the submatrix in rows K to L.
    //
    // Iterative loop for norm reduction.
    //
    REAL sfmin1 = Rlamch("S") / Rlamch("P");
    REAL sfmax1 = one / sfmin1;
    const REAL sclfac = 2.0;
    REAL sfmin2 = sfmin1 * sclfac;
    REAL sfmax2 = one / sfmin2;
    //
    noconv = true;
    REAL c = 0.0;
    REAL r = 0.0;
    INTEGER ica = 0;
    REAL ca = 0.0;
    INTEGER ira = 0;
    REAL ra = 0.0;
    REAL g = 0.0;
    REAL f = 0.0;
    REAL s = 0.0;
    const REAL factor = 0.95;
    while (noconv) {
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
            // Guard against zero C or R due to underflow.
            //
            if (c == zero || r == zero) {
                continue;
            }
            //
            // Exit if NaN to avoid infinite loop
            //
            if (Risnan(c + ca + r + ra)) {
                info = -3;
                Mxerbla("Cgebal", -info);
                return;
            }
            //
            g = r / sclfac;
            f = one;
            s = c + r;
            //
            while (c < g && max(f, c, ca) < sfmax2 && min(r, g, ra) > sfmin2) {
                f = f * sclfac;
                c = c * sclfac;
                ca = ca * sclfac;
                r = r / sclfac;
                g = g / sclfac;
                ra = ra / sclfac;
            }
            //
            g = c / sclfac;
            //
            while (g >= r && max(r, ra) < sfmax2 && min(f, c, g, ca) > sfmin2) {
                f = f / sclfac;
                c = c / sclfac;
                g = g / sclfac;
                ca = ca / sclfac;
                r = r * sclfac;
                ra = ra * sclfac;
            }
            //
            // Now balance.
            //
            if ((c + r) >= factor * s) {
                continue;
            }
            if (f < one && scale[i - 1] < one) {
                if (f * scale[i - 1] <= sfmin1) {
                    continue;
                }
            }
            if (f > one && scale[i - 1] > one) {
                if (scale[i - 1] >= sfmax1 / f) {
                    continue;
                }
            }
            g = one / f;
            scale[i - 1] = scale[i - 1] * f;
            noconv = true;
            //
            CRscal(n - k + 1, g, &a[(i - 1) + (k - 1) * lda], lda);
            CRscal(l, f, &a[(i - 1) * lda], 1);
            //
        }
        //
    }
    //
    ilo = k;
    ihi = l;
    //
    // End of Cgebal
    //
}
