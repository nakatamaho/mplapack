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

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_eig.h>

#include <mplapack_debug.h>

void Rlatm4(INTEGER const itype, INTEGER const n, INTEGER const nz1, INTEGER const nz2, INTEGER const isign, REAL const amagn, REAL const rcond, REAL const triang, INTEGER const idist, INTEGER *iseed, REAL *a, INTEGER const lda) {
    iseed([4]);
    a([lda * star]);
    const REAL zero = 0.0;
    INTEGER kbeg = 0;
    INTEGER kend = 0;
    INTEGER klen = 0;
    INTEGER isdb = 0;
    INTEGER isde = 0;
    INTEGER jd = 0;
    const REAL one = 1.0;
    INTEGER k = 0;
    REAL alpha = 0.0;
    INTEGER i = 0;
    const REAL two = 2.0;
    const REAL half = one / two;
    REAL temp = 0.0;
    REAL safmin = 0.0;
    REAL cl = 0.0;
    REAL sl = 0.0;
    REAL cr = 0.0;
    REAL sr = 0.0;
    REAL sv1 = 0.0;
    REAL sv2 = 0.0;
    INTEGER ioff = 0;
    INTEGER jr = 0;
    INTEGER jc = 0;
    //
    //  -- LAPACK test routine --
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
    if (n <= 0) {
        return;
    }
    Rlaset("Full", n, n, zero, zero, a, lda);
    //
    //     Insure a correct ISEED
    //
    if (mod(iseed[4 - 1], 2) != 1) {
        iseed[4 - 1]++;
    }
    //
    //     Compute diagonal and subdiagonal according to ITYPE, NZ1, NZ2,
    //     and RCOND
    //
    if (itype != 0) {
        if (abs(itype) >= 4) {
            kbeg = max({(INTEGER)1, min(n, nz1 + 1)});
            kend = max({kbeg, min(n, n - nz2)});
            klen = kend + 1 - kbeg;
        } else {
            kbeg = 1;
            kend = n;
            klen = n;
        }
        isdb = 1;
        isde = 0;
        switch (abs(itype)) {
        case 1:
            goto statement_10;
        case 2:
            goto statement_30;
        case 3:
            goto statement_50;
        case 4:
            goto statement_80;
        case 5:
            goto statement_100;
        case 6:
            goto statement_120;
        case 7:
            goto statement_140;
        case 8:
            goto statement_160;
        case 9:
            goto statement_180;
        case 10:
            goto statement_200;
        default:
            break;
        }
    //
    //        abs(ITYPE) = 1: Identity
    //
    statement_10:
        for (jd = 1; jd <= n; jd = jd + 1) {
            a[(jd - 1) + (jd - 1) * lda] = one;
        }
        goto statement_220;
    //
    //        abs(ITYPE) = 2: Transposed Jordan block
    //
    statement_30:
        for (jd = 1; jd <= n - 1; jd = jd + 1) {
            a[((jd + 1) - 1) + (jd - 1) * lda] = one;
        }
        isdb = 1;
        isde = n - 1;
        goto statement_220;
    //
    //        abs(ITYPE) = 3: Transposed Jordan block, followed by the
    //                        identity.
    //
    statement_50:
        k = (n - 1) / 2;
        for (jd = 1; jd <= k; jd = jd + 1) {
            a[((jd + 1) - 1) + (jd - 1) * lda] = one;
        }
        isdb = 1;
        isde = k;
        for (jd = k + 2; jd <= 2 * k + 1; jd = jd + 1) {
            a[(jd - 1) + (jd - 1) * lda] = one;
        }
        goto statement_220;
    //
    //        abs(ITYPE) = 4: 1,...,k
    //
    statement_80:
        for (jd = kbeg; jd <= kend; jd = jd + 1) {
            a[(jd - 1) + (jd - 1) * lda] = (jd - nz1).real();
        }
        goto statement_220;
    //
    //        abs(ITYPE) = 5: One large D value:
    //
    statement_100:
        for (jd = kbeg + 1; jd <= kend; jd = jd + 1) {
            a[(jd - 1) + (jd - 1) * lda] = rcond;
        }
        a[(kbeg - 1) + (kbeg - 1) * lda] = one;
        goto statement_220;
    //
    //        abs(ITYPE) = 6: One small D value:
    //
    statement_120:
        for (jd = kbeg; jd <= kend - 1; jd = jd + 1) {
            a[(jd - 1) + (jd - 1) * lda] = one;
        }
        a[(kend - 1) + (kend - 1) * lda] = rcond;
        goto statement_220;
    //
    //        abs(ITYPE) = 7: Exponentially distributed D values:
    //
    statement_140:
        a[(kbeg - 1) + (kbeg - 1) * lda] = one;
        if (klen > 1) {
            alpha = pow(rcond, [(one / (klen - 1).real()) - 1]);
            for (i = 2; i <= klen; i = i + 1) {
                a[((nz1 + i) - 1) + ((nz1 + i) - 1) * lda] = pow(alpha, (i - 1).real());
            }
        }
        goto statement_220;
    //
    //        abs(ITYPE) = 8: Arithmetically distributed D values:
    //
    statement_160:
        a[(kbeg - 1) + (kbeg - 1) * lda] = one;
        if (klen > 1) {
            alpha = (one - rcond) / (klen - 1).real();
            for (i = 2; i <= klen; i = i + 1) {
                a[((nz1 + i) - 1) + ((nz1 + i) - 1) * lda] = (klen - i).real() * alpha + rcond;
            }
        }
        goto statement_220;
    //
    //        abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):
    //
    statement_180:
        alpha = log(rcond);
        for (jd = kbeg; jd <= kend; jd = jd + 1) {
            a[(jd - 1) + (jd - 1) * lda] = exp(alpha * dlaran(iseed));
        }
        goto statement_220;
    //
    //        abs(ITYPE) = 10: Randomly distributed D values from DIST
    //
    statement_200:
        for (jd = kbeg; jd <= kend; jd = jd + 1) {
            a[(jd - 1) + (jd - 1) * lda] = dlarnd(idist, iseed);
        }
    //
    statement_220:
        //
        //        Scale by AMAGN
        //
        for (jd = kbeg; jd <= kend; jd = jd + 1) {
            a[(jd - 1) + (jd - 1) * lda] = amagn * a[(jd - 1) + (jd - 1) * lda].real();
        }
        for (jd = isdb; jd <= isde; jd = jd + 1) {
            a[((jd + 1) - 1) + (jd - 1) * lda] = amagn * a[((jd + 1) - 1) + (jd - 1) * lda].real();
        }
        //
        //        If ISIGN = 1 or 2, assign random signs to diagonal and
        //        subdiagonal
        //
        if (isign > 0) {
            for (jd = kbeg; jd <= kend; jd = jd + 1) {
                if (a[(jd - 1) + (jd - 1) * lda].real() != zero) {
                    if (dlaran(iseed) > half) {
                        a[(jd - 1) + (jd - 1) * lda] = -a[(jd - 1) + (jd - 1) * lda];
                    }
                }
            }
            for (jd = isdb; jd <= isde; jd = jd + 1) {
                if (a[((jd + 1) - 1) + (jd - 1) * lda].real() != zero) {
                    if (dlaran(iseed) > half) {
                        a[((jd + 1) - 1) + (jd - 1) * lda] = -a[((jd + 1) - 1) + (jd - 1) * lda];
                    }
                }
            }
        }
        //
        //        Reverse if ITYPE < 0
        //
        if (itype < 0) {
            for (jd = kbeg; jd <= (kbeg + kend - 1) / 2; jd = jd + 1) {
                temp = a[(jd - 1) + (jd - 1) * lda];
                a[(jd - 1) + (jd - 1) * lda] = a[((kbeg + kend - jd) - 1) + ((kbeg + kend - jd) - 1) * lda];
                a[((kbeg + kend - jd) - 1) + ((kbeg + kend - jd) - 1) * lda] = temp;
            }
            for (jd = 1; jd <= (n - 1) / 2; jd = jd + 1) {
                temp = a[((jd + 1) - 1) + (jd - 1) * lda];
                a[((jd + 1) - 1) + (jd - 1) * lda] = a[((n + 1 - jd) - 1) + ((n - jd) - 1) * lda];
                a[((n + 1 - jd) - 1) + ((n - jd) - 1) * lda] = temp;
            }
        }
        //
        //        If ISIGN = 2, and no subdiagonals already, then apply
        //        random rotations to make 2x2 blocks.
        //
        if (isign == 2 && itype != 2 && itype != 3) {
            safmin = Rlamch("S");
            for (jd = kbeg; jd <= kend - 1; jd = jd + 2) {
                if (dlaran(iseed) > half) {
                    //
                    //                 Rotation on left.
                    //
                    cl = two * dlaran(iseed) - one;
                    sl = two * dlaran(iseed) - one;
                    temp = one / max(safmin, sqrt(pow2(cl) + pow2(sl)));
                    cl = cl * temp;
                    sl = sl * temp;
                    //
                    //                 Rotation on right.
                    //
                    cr = two * dlaran(iseed) - one;
                    sr = two * dlaran(iseed) - one;
                    temp = one / max(safmin, sqrt(pow2(cr) + pow2(sr)));
                    cr = cr * temp;
                    sr = sr * temp;
                    //
                    //                 Apply
                    //
                    sv1 = a[(jd - 1) + (jd - 1) * lda];
                    sv2 = a[((jd + 1) - 1) + ((jd + 1) - 1) * lda];
                    a[(jd - 1) + (jd - 1) * lda] = cl * cr * sv1 + sl * sr * sv2;
                    a[((jd + 1) - 1) + (jd - 1) * lda] = -sl * cr * sv1 + cl * sr * sv2;
                    a[(jd - 1) + ((jd + 1) - 1) * lda] = -cl * sr * sv1 + sl * cr * sv2;
                    a[((jd + 1) - 1) + ((jd + 1) - 1) * lda] = sl * sr * sv1 + cl * cr * sv2;
                }
            }
        }
        //
    }
    //
    //     Fill in upper triangle (except for 2x2 blocks)
    //
    if (triang != zero) {
        if (isign != 2 || itype == 2 || itype == 3) {
            ioff = 1;
        } else {
            ioff = 2;
            for (jr = 1; jr <= n - 1; jr = jr + 1) {
                if (a[((jr + 1) - 1) + (jr - 1) * lda] == zero) {
                    a[(jr - 1) + ((jr + 1) - 1) * lda] = triang * dlarnd(idist, iseed);
                }
            }
        }
        //
        for (jc = 2; jc <= n; jc = jc + 1) {
            for (jr = 1; jr <= jc - ioff; jr = jr + 1) {
                a[(jr - 1) + (jc - 1) * lda] = triang * dlarnd(idist, iseed);
            }
        }
    }
    //
    //     End of Rlatm4
    //
}
