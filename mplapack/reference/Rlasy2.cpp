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

void Rlasy2(bool const ltranl, bool const ltranr, INTEGER const isgn, INTEGER const n1, INTEGER const n2, REAL *tl, INTEGER const ldtl, REAL *tr, INTEGER const ldtr, REAL *b, INTEGER const ldb, REAL &scale, REAL *x, INTEGER const ldx, REAL &xnorm, INTEGER &info) {
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    REAL sgn = 0.0;
    INTEGER k = 0;
    REAL tau1 = 0.0;
    REAL bet = 0.0;
    const REAL one = 1.0;
    REAL gam = 0.0;
    REAL smin = 0.0;
    REAL tmp[4];
    REAL btmp[4];
    INTEGER ipiv = 0;
    REAL u11 = 0.0;
    REAL u12 = 0.0;
    REAL l21 = 0.0;
    REAL u22 = 0.0;
    bool xswap = false;
    bool bswap = false;
    REAL temp = 0.0;
    const REAL two = 2.0e+0;
    const REAL half = 0.5e+0;
    REAL x2[4];
    const REAL zero = 0.0;
    REAL t16[4 * 4];
    INTEGER ldt16 = 4;
    INTEGER i = 0;
    REAL xmax = 0.0;
    INTEGER ip = 0;
    INTEGER jp = 0;
    INTEGER ipsv = 0;
    INTEGER jpsv = 0;
    INTEGER j = 0;
    INTEGER jpiv[4];
    static INTEGER locl12[4] = {3, 4, 1, 2};
    static INTEGER locu21[4] = {2, 1, 4, 3};
    static INTEGER locu22[4] = {4, 3, 2, 1};
    static bool xswpiv[4] = {false, false, true, true};
    static bool bswpiv[4] = {false, true, false, true};
    const REAL eight = 8.0e+0;
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    // =====================================================================
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
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Do not check the input parameters for errors
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (n1 == 0 || n2 == 0) {
        return;
    }
    //
    //     Set constants to control overflow
    //
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    sgn = isgn;
    //
    k = n1 + n1 + n2 - 2;
    switch (k) {
    case 1:
        goto statement_10;
    case 2:
        goto statement_20;
    case 3:
        goto statement_30;
    case 4:
        goto statement_50;
    default:
        break;
    }
//
//     1 by 1: TL11*X + SGN*X*TR11 = B11
//
statement_10:
    tau1 = tl[(1 - 1)] + sgn * tr[(1 - 1)];
    bet = abs(tau1);
    if (bet <= smlnum) {
        tau1 = smlnum;
        bet = smlnum;
        info = 1;
    }
    //
    scale = one;
    gam = abs(b[(1 - 1)]);
    if (smlnum * gam > bet) {
        scale = one / gam;
    }
    //
    x[(1 - 1)] = (b[(1 - 1)] * scale) / tau1;
    xnorm = abs(x[(1 - 1)]);
    return;
//
//     1 by 2:
//     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
//                                       [TR21 TR22]
//
statement_20:
    //
    smin = max(eps * max({abs(tl[(1 - 1)]), abs(tr[(1 - 1)]), abs(tr[(2 - 1) * ldtr]), abs(tr[(2 - 1)]), abs(tr[(2 - 1) + (2 - 1) * ldtr])}), smlnum);
    tmp[1 - 1] = tl[(1 - 1)] + sgn * tr[(1 - 1)];
    tmp[4 - 1] = tl[(1 - 1)] + sgn * tr[(2 - 1) + (2 - 1) * ldtr];
    if (ltranr) {
        tmp[2 - 1] = sgn * tr[(2 - 1)];
        tmp[3 - 1] = sgn * tr[(2 - 1) * ldtr];
    } else {
        tmp[2 - 1] = sgn * tr[(2 - 1) * ldtr];
        tmp[3 - 1] = sgn * tr[(2 - 1)];
    }
    btmp[1 - 1] = b[(1 - 1)];
    btmp[2 - 1] = b[(2 - 1) * ldb];
    goto statement_40;
//
//     2 by 1:
//          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
//            [TL21 TL22] [X21]         [X21]         [B21]
//
statement_30:
    smin = max(eps * max({abs(tr[(1 - 1)]), abs(tl[(1 - 1)]), abs(tl[(2 - 1) * ldtl]), abs(tl[(2 - 1)]), abs(tl[(2 - 1) + (2 - 1) * ldtl])}), smlnum);
    tmp[1 - 1] = tl[(1 - 1)] + sgn * tr[(1 - 1)];
    tmp[4 - 1] = tl[(2 - 1) + (2 - 1) * ldtl] + sgn * tr[(1 - 1)];
    if (ltranl) {
        tmp[2 - 1] = tl[(2 - 1) * ldtl];
        tmp[3 - 1] = tl[(2 - 1)];
    } else {
        tmp[2 - 1] = tl[(2 - 1)];
        tmp[3 - 1] = tl[(2 - 1) * ldtl];
    }
    btmp[1 - 1] = b[(1 - 1)];
    btmp[2 - 1] = b[(2 - 1)];
statement_40:
    //
    //     Solve 2 by 2 system using complete pivoting.
    //     Set pivots less than SMIN to SMIN.
    //
    ipiv = iRamax(4, tmp, 1);
    u11 = tmp[ipiv - 1];
    if (abs(u11) <= smin) {
        info = 1;
        u11 = smin;
    }
    u12 = tmp[locu12[ipiv - 1] - 1];
    l21 = tmp[locl21[ipiv - 1] - 1] / u11;
    u22 = tmp[locu22[ipiv - 1] - 1] - u12 * l21;
    xswap = xswpiv[ipiv - 1];
    bswap = bswpiv[ipiv - 1];
    if (abs(u22) <= smin) {
        info = 1;
        u22 = smin;
    }
    if (bswap) {
        temp = btmp[2 - 1];
        btmp[2 - 1] = btmp[1 - 1] - l21 * temp;
        btmp[1 - 1] = temp;
    } else {
        btmp[2 - 1] = btmp[2 - 1] - l21 * btmp[1 - 1];
    }
    scale = one;
    if ((two * smlnum) * abs(btmp[2 - 1]) > abs(u22) || (two * smlnum) * abs(btmp[1 - 1]) > abs(u11)) {
        scale = half / max(abs(btmp[1 - 1]), abs(btmp[2 - 1]));
        btmp[1 - 1] = btmp[1 - 1] * scale;
        btmp[2 - 1] = btmp[2 - 1] * scale;
    }
    x2[2 - 1] = btmp[2 - 1] / u22;
    x2[1 - 1] = btmp[1 - 1] / u11 - (u12 / u11) * x2[2 - 1];
    if (xswap) {
        temp = x2[2 - 1];
        x2[2 - 1] = x2[1 - 1];
        x2[1 - 1] = temp;
    }
    x[(1 - 1)] = x2[1 - 1];
    if (n1 == 1) {
        x[(2 - 1) * ldx] = x2[2 - 1];
        xnorm = abs(x[(1 - 1)]) + abs(x[(2 - 1) * ldx]);
    } else {
        x[(2 - 1)] = x2[2 - 1];
        xnorm = max(abs(x[(1 - 1)]), abs(x[(2 - 1)]));
    }
    return;
//
//     2 by 2:
//     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
//       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
//
//     Solve equivalent 4 by 4 system using complete pivoting.
//     Set pivots less than SMIN to SMIN.
//
statement_50:
    smin = max({abs(tr[(1 - 1)]), abs(tr[(2 - 1) * ldtr]), abs(tr[(2 - 1)]), abs(tr[(2 - 1) + (2 - 1) * ldtr])});
    smin = max({smin, abs(tl[(1 - 1)]), abs(tl[(2 - 1) * ldtl]), abs(tl[(2 - 1)]), abs(tl[(2 - 1) + (2 - 1) * ldtl])});
    smin = max(eps * smin, smlnum);
    btmp[1 - 1] = zero;
    Rcopy(16, btmp, 0, t16, 1);
    t16[(1 - 1)] = tl[(1 - 1)] + sgn * tr[(1 - 1)];
    t16[(2 - 1) + (2 - 1) * ldt16] = tl[(2 - 1) + (2 - 1) * ldtl] + sgn * tr[(1 - 1)];
    t16[(3 - 1) + (3 - 1) * ldt16] = tl[(1 - 1)] + sgn * tr[(2 - 1) + (2 - 1) * ldtr];
    t16[(4 - 1) + (4 - 1) * ldt16] = tl[(2 - 1) + (2 - 1) * ldtl] + sgn * tr[(2 - 1) + (2 - 1) * ldtr];
    if (ltranl) {
        t16[(2 - 1) * ldt16] = tl[(2 - 1)];
        t16[(2 - 1)] = tl[(2 - 1) * ldtl];
        t16[(3 - 1) + (4 - 1) * ldt16] = tl[(2 - 1)];
        t16[(4 - 1) + (3 - 1) * ldt16] = tl[(2 - 1) * ldtl];
    } else {
        t16[(2 - 1) * ldt16] = tl[(2 - 1) * ldtl];
        t16[(2 - 1)] = tl[(2 - 1)];
        t16[(3 - 1) + (4 - 1) * ldt16] = tl[(2 - 1) * ldtl];
        t16[(4 - 1) + (3 - 1) * ldt16] = tl[(2 - 1)];
    }
    if (ltranr) {
        t16[(3 - 1) * ldt16] = sgn * tr[(2 - 1) * ldtr];
        t16[(2 - 1) + (4 - 1) * ldt16] = sgn * tr[(2 - 1) * ldtr];
        t16[(3 - 1)] = sgn * tr[(2 - 1)];
        t16[(4 - 1) + (2 - 1) * ldt16] = sgn * tr[(2 - 1)];
    } else {
        t16[(3 - 1) * ldt16] = sgn * tr[(2 - 1)];
        t16[(2 - 1) + (4 - 1) * ldt16] = sgn * tr[(2 - 1)];
        t16[(3 - 1)] = sgn * tr[(2 - 1) * ldtr];
        t16[(4 - 1) + (2 - 1) * ldt16] = sgn * tr[(2 - 1) * ldtr];
    }
    btmp[1 - 1] = b[(1 - 1)];
    btmp[2 - 1] = b[(2 - 1)];
    btmp[3 - 1] = b[(2 - 1) * ldb];
    btmp[4 - 1] = b[(2 - 1) + (2 - 1) * ldb];
    //
    //     Perform elimination
    //
    for (i = 1; i <= 3; i = i + 1) {
        xmax = zero;
        for (ip = i; ip <= 4; ip = ip + 1) {
            for (jp = i; jp <= 4; jp = jp + 1) {
                if (abs(t16[(ip - 1) + (jp - 1) * ldt16]) >= xmax) {
                    xmax = abs(t16[(ip - 1) + (jp - 1) * ldt16]);
                    ipsv = ip;
                    jpsv = jp;
                }
            }
        }
        if (ipsv != i) {
            Rswap(4, &t16[(ipsv - 1)], 4, &t16[(i - 1)], 4);
            temp = btmp[i - 1];
            btmp[i - 1] = btmp[ipsv - 1];
            btmp[ipsv - 1] = temp;
        }
        if (jpsv != i) {
            Rswap(4, &t16[(jpsv - 1) * ldt16], 1, &t16[(i - 1) * ldt16], 1);
        }
        jpiv[i - 1] = jpsv;
        if (abs(t16[(i - 1) + (i - 1) * ldt16]) < smin) {
            info = 1;
            t16[(i - 1) + (i - 1) * ldt16] = smin;
        }
        for (j = i + 1; j <= 4; j = j + 1) {
            t16[(j - 1) + (i - 1) * ldt16] = t16[(j - 1) + (i - 1) * ldt16] / t16[(i - 1) + (i - 1) * ldt16];
            btmp[j - 1] = btmp[j - 1] - t16[(j - 1) + (i - 1) * ldt16] * btmp[i - 1];
            for (k = i + 1; k <= 4; k = k + 1) {
                t16[(j - 1) + (k - 1) * ldt16] = t16[(j - 1) + (k - 1) * ldt16] - t16[(j - 1) + (i - 1) * ldt16] * t16[(i - 1) + (k - 1) * ldt16];
            }
        }
    }
    if (abs(t16[(4 - 1) + (4 - 1) * ldt16]) < smin) {
        info = 1;
        t16[(4 - 1) + (4 - 1) * ldt16] = smin;
    }
    scale = one;
    if ((eight * smlnum) * abs(btmp[1 - 1]) > abs(t16[(1 - 1)]) || (eight * smlnum) * abs(btmp[2 - 1]) > abs(t16[(2 - 1) + (2 - 1) * ldt16]) || (eight * smlnum) * abs(btmp[3 - 1]) > abs(t16[(3 - 1) + (3 - 1) * ldt16]) || (eight * smlnum) * abs(btmp[4 - 1]) > abs(t16[(4 - 1) + (4 - 1) * ldt16])) {
        scale = (one / eight) / max({abs(btmp[1 - 1]), abs(btmp[2 - 1]), abs(btmp[3 - 1]), abs(btmp[4 - 1])});
        btmp[1 - 1] = btmp[1 - 1] * scale;
        btmp[2 - 1] = btmp[2 - 1] * scale;
        btmp[3 - 1] = btmp[3 - 1] * scale;
        btmp[4 - 1] = btmp[4 - 1] * scale;
    }
    for (i = 1; i <= 4; i = i + 1) {
        k = 5 - i;
        temp = one / t16[(k - 1) + (k - 1) * ldt16];
        tmp[k - 1] = btmp[k - 1] * temp;
        for (j = k + 1; j <= 4; j = j + 1) {
            tmp[k - 1] = tmp[k - 1] - (temp * t16[(k - 1) + (j - 1) * ldt16]) * tmp[j - 1];
        }
    }
    for (i = 1; i <= 3; i = i + 1) {
        if (jpiv[(4 - i) - 1] != 4 - i) {
            temp = tmp[(4 - i) - 1];
            tmp[(4 - i) - 1] = tmp[(jpiv[(4 - i) - 1]) - 1];
            tmp[(jpiv[(4 - i) - 1]) - 1] = temp;
        }
    }
    x[(1 - 1)] = tmp[1 - 1];
    x[(2 - 1)] = tmp[2 - 1];
    x[(2 - 1) * ldx] = tmp[3 - 1];
    x[(2 - 1) + (2 - 1) * ldx] = tmp[4 - 1];
    xnorm = max(abs(tmp[1 - 1]) + abs(tmp[3 - 1]), abs(tmp[2 - 1]) + abs(tmp[4 - 1]));
    //
    //     End of Rlasy2
    //
}
