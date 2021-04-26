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

#include <mplapack_matgen.h>

void Rlatme(INTEGER const n, const char *dist, INTEGER *iseed, REAL *d, INTEGER const mode, REAL const cond, REAL const dmax, const char *ei, const char *rsign, const char *upper, const char *sim, REAL *ds, INTEGER const modes, REAL const conds, INTEGER const kl, INTEGER const ku, REAL const anorm, REAL *a, INTEGER const lda, REAL *work, INTEGER &info) {
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
    //     1)      Decode and Test the input parameters.
    //             Initialize flags & seed.
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    //     Decode DIST
    //
    INTEGER idist = 0;
    if (Mlsame(dist, "U")) {
        idist = 1;
    } else if (Mlsame(dist, "S")) {
        idist = 2;
    } else if (Mlsame(dist, "N")) {
        idist = 3;
    } else {
        idist = -1;
    }
    //
    //     Check EI
    //
    bool useei = true;
    bool badei = false;
    INTEGER j = 0;
    if (Mlsame(&ei[1 - 1], " ") || mode != 0) {
        useei = false;
    } else {
        if (Mlsame(&ei[1 - 1], "R")) {
            for (j = 2; j <= n; j = j + 1) {
                if (Mlsame(&ei[j - 1], "I")) {
                    if (Mlsame(&ei[(j - 1) - 1], "I")) {
                        badei = true;
                    }
                } else {
                    if (!Mlsame(&ei[j - 1], "R")) {
                        badei = true;
                    }
                }
            }
        } else {
            badei = true;
        }
    }
    //
    //     Decode RSIGN
    //
    INTEGER irsign = 0;
    if (Mlsame(rsign, "T")) {
        irsign = 1;
    } else if (Mlsame(rsign, "F")) {
        irsign = 0;
    } else {
        irsign = -1;
    }
    //
    //     Decode UPPER
    //
    INTEGER iupper = 0;
    if (Mlsame(upper, "T")) {
        iupper = 1;
    } else if (Mlsame(upper, "F")) {
        iupper = 0;
    } else {
        iupper = -1;
    }
    //
    //     Decode SIM
    //
    INTEGER isim = 0;
    if (Mlsame(sim, "T")) {
        isim = 1;
    } else if (Mlsame(sim, "F")) {
        isim = 0;
    } else {
        isim = -1;
    }
    //
    //     Check DS, if MODES=0 and ISIM=1
    //
    bool bads = false;
    const REAL zero = 0.0;
    if (modes == 0 && isim == 1) {
        for (j = 1; j <= n; j = j + 1) {
            if (ds[j - 1] == zero) {
                bads = true;
            }
        }
    }
    //
    //     Set INFO if an error
    //
    const REAL one = 1.0;
    if (n < 0) {
        info = -1;
    } else if (idist == -1) {
        info = -2;
    } else if (abs(mode) > 6) {
        info = -5;
    } else if ((mode != 0 && abs(mode) != 6) && cond < one) {
        info = -6;
    } else if (badei) {
        info = -8;
    } else if (irsign == -1) {
        info = -9;
    } else if (iupper == -1) {
        info = -10;
    } else if (isim == -1) {
        info = -11;
    } else if (bads) {
        info = -12;
    } else if (isim == 1 && abs(modes) > 5) {
        info = -13;
    } else if (isim == 1 && modes != 0 && conds < one) {
        info = -14;
    } else if (kl < 1) {
        info = -15;
    } else if (ku < 1 || (ku < n - 1 && kl < n - 1)) {
        info = -16;
    } else if (lda < max((INTEGER)1, n)) {
        info = -19;
    }
    //
    if (info != 0) {
        Mxerbla("Rlatme", -info);
        return;
    }
    //
    //     Initialize random number generator
    //
    INTEGER i = 0;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = mod(abs(iseed[i - 1]), 4096);
    }
    //
    if (mod(iseed[4 - 1], 2) != 1) {
        iseed[4 - 1]++;
    }
    //
    //     2)      Set up diagonal of A
    //
    //             Compute D according to COND and MODE
    //
    INTEGER iinfo = 0;
    Rlatm1(mode, cond, irsign, idist, iseed, d, n, iinfo);
    if (iinfo != 0) {
        info = 1;
        return;
    }
    REAL temp = 0.0;
    REAL alpha = 0.0;
    if (mode != 0 && abs(mode) != 6) {
        //
        //        Scale by DMAX
        //
        temp = abs(d[1 - 1]);
        for (i = 2; i <= n; i = i + 1) {
            temp = max(temp, abs(d[i - 1]));
        }
        //
        if (temp > zero) {
            alpha = dmax / temp;
        } else if (dmax != zero) {
            info = 2;
            return;
        } else {
            alpha = zero;
        }
        //
        Rscal(n, alpha, d, 1);
        //
    }
    //
    Rlaset("Full", n, n, zero, zero, a, lda);
    Rcopy(n, d, 1, a, lda + 1);
    //
    //     Set up complex conjugate pairs
    //
    const REAL half = 1.0 / 2.0;
    if (mode == 0) {
        if (useei) {
            for (j = 2; j <= n; j = j + 1) {
                if (Mlsame(&ei[j - 1], "I")) {
                    a[((j - 1) - 1) + (j - 1) * lda] = a[(j - 1) + (j - 1) * lda];
                    a[(j - 1) + ((j - 1) - 1) * lda] = -a[(j - 1) + (j - 1) * lda];
                    a[(j - 1) + (j - 1) * lda] = a[((j - 1) - 1) + ((j - 1) - 1) * lda];
                }
            }
        }
        //
    } else if (abs(mode) == 5) {
        //
        for (j = 2; j <= n; j = j + 2) {
            if (Rlaran(iseed) > half) {
                a[((j - 1) - 1) + (j - 1) * lda] = a[(j - 1) + (j - 1) * lda];
                a[(j - 1) + ((j - 1) - 1) * lda] = -a[(j - 1) + (j - 1) * lda];
                a[(j - 1) + (j - 1) * lda] = a[((j - 1) - 1) + ((j - 1) - 1) * lda];
            }
        }
    }
    //
    //     3)      If UPPER='T', set upper triangle of A to random numbers.
    //             (but don't modify the corners of 2x2 blocks.)
    //
    INTEGER jc = 0;
    INTEGER jr = 0;
    if (iupper != 0) {
        for (jc = 2; jc <= n; jc = jc + 1) {
            if (a[((jc - 1) - 1) + (jc - 1) * lda] != zero) {
                jr = jc - 2;
            } else {
                jr = jc - 1;
            }
            Rlarnv(idist, iseed, jr, &a[(jc - 1) * lda]);
        }
    }
    //
    //     4)      If SIM='T', apply similarity transformation.
    //
    //                                -1
    //             Transform is  X A X  , where X = U S V, thus
    //
    //             it is  U S V A V' (1/S) U'
    //
    if (isim != 0) {
        //
        //        Compute S (singular values of the eigenvector matrix)
        //        according to CONDS and MODES
        //
        Rlatm1(modes, conds, 0, 0, iseed, ds, n, iinfo);
        if (iinfo != 0) {
            info = 3;
            return;
        }
        //
        //        Multiply by V and V'
        //
        Rlarge(n, a, lda, iseed, work, iinfo);
        if (iinfo != 0) {
            info = 4;
            return;
        }
        //
        //        Multiply by S and (1/S)
        //
        for (j = 1; j <= n; j = j + 1) {
            Rscal(n, ds[j - 1], &a[(j - 1)], lda);
            if (ds[j - 1] != zero) {
                Rscal(n, one / ds[j - 1], &a[(j - 1) * lda], 1);
            } else {
                info = 5;
                return;
            }
        }
        //
        //        Multiply by U and U'
        //
        Rlarge(n, a, lda, iseed, work, iinfo);
        if (iinfo != 0) {
            info = 4;
            return;
        }
    }
    //
    //     5)      Reduce the bandwidth.
    //
    INTEGER jcr = 0;
    INTEGER ic = 0;
    INTEGER irows = 0;
    INTEGER icols = 0;
    REAL xnorms = 0.0;
    REAL tau = 0.0;
    INTEGER ir = 0;
    if (kl < n - 1) {
        //
        //        Reduce bandwidth -- kill column
        //
        for (jcr = kl + 1; jcr <= n - 1; jcr = jcr + 1) {
            ic = jcr - kl;
            irows = n + 1 - jcr;
            icols = n + kl - jcr;
            //
            Rcopy(irows, &a[(jcr - 1) + (ic - 1) * lda], 1, work, 1);
            xnorms = work[1 - 1];
            Rlarfg(irows, xnorms, &work[2 - 1], 1, tau);
            work[1 - 1] = one;
            //
            Rgemv("T", irows, icols, one, &a[(jcr - 1) + ((ic + 1) - 1) * lda], lda, work, 1, zero, &work[(irows + 1) - 1], 1);
            Rger(irows, icols, -tau, work, 1, &work[(irows + 1) - 1], 1, &a[(jcr - 1) + ((ic + 1) - 1) * lda], lda);
            //
            Rgemv("N", n, irows, one, &a[(jcr - 1) * lda], lda, work, 1, zero, &work[(irows + 1) - 1], 1);
            Rger(n, irows, -tau, &work[(irows + 1) - 1], 1, work, 1, &a[(jcr - 1) * lda], lda);
            //
            a[(jcr - 1) + (ic - 1) * lda] = xnorms;
            Rlaset("Full", irows - 1, 1, zero, zero, &a[((jcr + 1) - 1) + (ic - 1) * lda], lda);
        }
    } else if (ku < n - 1) {
        //
        //        Reduce upper bandwidth -- kill a row at a time.
        //
        for (jcr = ku + 1; jcr <= n - 1; jcr = jcr + 1) {
            ir = jcr - ku;
            irows = n + ku - jcr;
            icols = n + 1 - jcr;
            //
            Rcopy(icols, &a[(ir - 1) + (jcr - 1) * lda], lda, work, 1);
            xnorms = work[1 - 1];
            Rlarfg(icols, xnorms, &work[2 - 1], 1, tau);
            work[1 - 1] = one;
            //
            Rgemv("N", irows, icols, one, &a[((ir + 1) - 1) + (jcr - 1) * lda], lda, work, 1, zero, &work[(icols + 1) - 1], 1);
            Rger(irows, icols, -tau, &work[(icols + 1) - 1], 1, work, 1, &a[((ir + 1) - 1) + (jcr - 1) * lda], lda);
            //
            Rgemv("C", icols, n, one, &a[(jcr - 1)], lda, work, 1, zero, &work[(icols + 1) - 1], 1);
            Rger(icols, n, -tau, work, 1, &work[(icols + 1) - 1], 1, &a[(jcr - 1)], lda);
            //
            a[(ir - 1) + (jcr - 1) * lda] = xnorms;
            Rlaset("Full", 1, icols - 1, zero, zero, &a[(ir - 1) + ((jcr + 1) - 1) * lda], lda);
        }
    }
    //
    //     Scale the matrix to have norm ANORM
    //
    REAL tempa[1];
    if (anorm >= zero) {
        temp = Rlange("M", n, n, a, lda, tempa);
        if (temp > zero) {
            alpha = anorm / temp;
            for (j = 1; j <= n; j = j + 1) {
                Rscal(n, alpha, &a[(j - 1) * lda], 1);
            }
        }
    }
    //
    //     End of Rlatme
    //
}
