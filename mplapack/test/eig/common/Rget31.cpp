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

void Rget31(REAL &rmax, INTEGER &lmax, INTEGER *ninfo, INTEGER &knt) {
    FEM_CMN_SVE(Rget31);
    ninfo([2]);
    // SAVE
    bool *ltrans(sve.ltrans, dim1(0, 1));
    //
    if (is_called_first_time) {
        static const bool values[] = {false, true};
        data_of_type<bool>(FEM_VALUES_AND_SIZE), ltrans;
    }
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
    //     Get machine parameters
    //
    REAL eps = Rlamch("P");
    REAL unfl = Rlamch("U");
    REAL smlnum = Rlamch("S") / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Set up test case parameters
    //
    REAL vsmin[4];
    vsmin[1 - 1] = smlnum;
    vsmin[2 - 1] = eps;
    const REAL ten = 10.0;
    vsmin[3 - 1] = one / (ten * ten);
    vsmin[4 - 1] = one / eps;
    REAL vab[3];
    vab[1 - 1] = sqrt(smlnum);
    vab[2 - 1] = one;
    vab[3 - 1] = sqrt(bignum);
    const REAL zero = 0.0;
    REAL vwr[4];
    vwr[1 - 1] = zero;
    const REAL half = 0.5e0;
    vwr[2 - 1] = half;
    const REAL two = 2.0;
    vwr[3 - 1] = two;
    vwr[4 - 1] = one;
    REAL vwi[4];
    vwi[1 - 1] = smlnum;
    vwi[2 - 1] = eps;
    vwi[3 - 1] = one;
    vwi[4 - 1] = two;
    REAL vdd[4];
    vdd[1 - 1] = sqrt(smlnum);
    vdd[2 - 1] = one;
    vdd[3 - 1] = two;
    vdd[4 - 1] = sqrt(bignum);
    REAL vca[5];
    vca[1 - 1] = zero;
    vca[2 - 1] = sqrt(smlnum);
    vca[3 - 1] = eps;
    vca[4 - 1] = half;
    vca[5 - 1] = one;
    //
    knt = 0;
    ninfo[1 - 1] = 0;
    ninfo[2 - 1] = 0;
    lmax = 0;
    rmax = zero;
    //
    //     Begin test loop
    //
    INTEGER id1 = 0;
    REAL d1 = 0.0;
    INTEGER id2 = 0;
    REAL d2 = 0.0;
    INTEGER ica = 0;
    REAL ca = 0.0;
    INTEGER itrans = 0;
    INTEGER ismin = 0;
    REAL smin = 0.0;
    INTEGER na = 0;
    INTEGER nw = 0;
    INTEGER ia = 0;
    INTEGER &a[(1 - 1) + (1 - 1) * lda] = 0;
    INTEGER ib = 0;
    INTEGER &b[(1 - 1) + (1 - 1) * ldb] = 0;
    INTEGER iwr = 0;
    REAL wr = 0.0;
    REAL wi = 0.0;
    REAL a[2 * 2];
    REAL b[2 * 2];
    REAL x[2 * 2];
    REAL scale = 0.0;
    REAL xnorm = 0.0;
    INTEGER info = 0;
    REAL res = 0.0;
    REAL den = 0.0;
    INTEGER &b[(1 - 1) + (2 - 1) * ldb] = 0;
    INTEGER iwi = 0;
    const REAL three = 3.0;
    INTEGER &a[(1 - 1) + (2 - 1) * lda] = 0;
    const REAL seven = 7.0;
    const REAL twnone = 21.0;
    REAL tmp = 0.0;
    const REAL four = 4.0;
    for (id1 = 1; id1 <= 4; id1 = id1 + 1) {
        d1 = vdd[id1 - 1];
        for (id2 = 1; id2 <= 4; id2 = id2 + 1) {
            d2 = vdd[id2 - 1];
            for (ica = 1; ica <= 5; ica = ica + 1) {
                ca = vca[ica - 1];
                for (itrans = 0; itrans <= 1; itrans = itrans + 1) {
                    for (ismin = 1; ismin <= 4; ismin = ismin + 1) {
                        smin = vsmin[ismin - 1];
                        //
                        na = 1;
                        nw = 1;
                        for (ia = 1; ia <= 3; ia = ia + 1) {
                            &a[(1 - 1) + (1 - 1) * lda] = vab[ia - 1];
                            for (ib = 1; ib <= 3; ib = ib + 1) {
                                &b[(1 - 1) + (1 - 1) * ldb] = vab[ib - 1];
                                for (iwr = 1; iwr <= 4; iwr = iwr + 1) {
                                    if (d1 == one && d2 == one && ca == one) {
                                        wr = vwr[iwr - 1] * &a[(1 - 1) + (1 - 1) * lda];
                                    } else {
                                        wr = vwr[iwr - 1];
                                    }
                                    wi = zero;
                                    Rlaln2(ltrans[itrans - 1], na, nw, smin, ca, a, 2, d1, d2, b, 2, wr, wi, x, 2, scale, xnorm, info);
                                    if (info < 0) {
                                        ninfo[1 - 1]++;
                                    }
                                    if (info > 0) {
                                        ninfo[2 - 1]++;
                                    }
                                    res = abs((ca * &a[(1 - 1) + (1 - 1) * lda] - wr * d1) * x[(1 - 1)] - scale * &b[(1 - 1) + (1 - 1) * ldb]);
                                    if (info == 0) {
                                        den = max(eps * (abs((ca * a[(1 - 1)] - wr * d1) * x[(1 - 1)])), smlnum);
                                    } else {
                                        den = max(smin * abs(x[(1 - 1)]), smlnum);
                                    }
                                    res = res / den;
                                    if (abs(x[(1 - 1)]) < unfl && abs(&b[(1 - 1) + (1 - 1) * ldb]) <= smlnum * abs(ca * &a[(1 - 1) + (1 - 1) * lda] - wr * d1)) {
                                        res = zero;
                                    }
                                    if (scale > one) {
                                        res += one / eps;
                                    }
                                    res += abs(xnorm - abs(x[(1 - 1)])) / max(smlnum, xnorm) / eps;
                                    if (info != 0 && info != 1) {
                                        res += one / eps;
                                    }
                                    knt++;
                                    if (res > rmax) {
                                        lmax = knt;
                                        rmax = res;
                                    }
                                }
                            }
                        }
                        //
                        na = 1;
                        nw = 2;
                        for (ia = 1; ia <= 3; ia = ia + 1) {
                            &a[(1 - 1) + (1 - 1) * lda] = vab[ia - 1];
                            for (ib = 1; ib <= 3; ib = ib + 1) {
                                &b[(1 - 1) + (1 - 1) * ldb] = vab[ib - 1];
                                &b[(1 - 1) + (2 - 1) * ldb] = -half * vab[ib - 1];
                                for (iwr = 1; iwr <= 4; iwr = iwr + 1) {
                                    if (d1 == one && d2 == one && ca == one) {
                                        wr = vwr[iwr - 1] * &a[(1 - 1) + (1 - 1) * lda];
                                    } else {
                                        wr = vwr[iwr - 1];
                                    }
                                    for (iwi = 1; iwi <= 4; iwi = iwi + 1) {
                                        if (d1 == one && d2 == one && ca == one) {
                                            wi = vwi[iwi - 1] * &a[(1 - 1) + (1 - 1) * lda];
                                        } else {
                                            wi = vwi[iwi - 1];
                                        }
                                        Rlaln2(ltrans[itrans - 1], na, nw, smin, ca, a, 2, d1, d2, b, 2, wr, wi, x, 2, scale, xnorm, info);
                                        if (info < 0) {
                                            ninfo[1 - 1]++;
                                        }
                                        if (info > 0) {
                                            ninfo[2 - 1]++;
                                        }
                                        res = abs((ca * &a[(1 - 1) + (1 - 1) * lda] - wr * d1) * x[(1 - 1)] + (wi * d1) * x[(2 - 1) * ldx] - scale * &b[(1 - 1) + (1 - 1) * ldb]);
                                        res += abs((-wi * d1) * x[(1 - 1)] + (ca * &a[(1 - 1) + (1 - 1) * lda] - wr * d1) * x[(2 - 1) * ldx] - scale * &b[(1 - 1) + (2 - 1) * ldb]);
                                        if (info == 0) {
                                            den = max({eps * (max(abs(ca * a[(1 - 1)] - wr * d1), abs(d1 * wi)) * (abs(x[(1 - 1)]) + abs(x[(2 - 1) * ldx]))), smlnum});
                                        } else {
                                            den = max(smin * (abs(x[(1 - 1)]) + abs(x[(2 - 1) * ldx])), smlnum);
                                        }
                                        res = res / den;
                                        if (abs(x[(1 - 1)]) < unfl && abs(x[(2 - 1) * ldx]) < unfl && abs(&b[(1 - 1) + (1 - 1) * ldb]) <= smlnum * abs(ca * &a[(1 - 1) + (1 - 1) * lda] - wr * d1)) {
                                            res = zero;
                                        }
                                        if (scale > one) {
                                            res += one / eps;
                                        }
                                        res += abs(xnorm - abs(x[(1 - 1)]) - abs(x[(2 - 1) * ldx])) / max(smlnum, xnorm) / eps;
                                        if (info != 0 && info != 1) {
                                            res += one / eps;
                                        }
                                        knt++;
                                        if (res > rmax) {
                                            lmax = knt;
                                            rmax = res;
                                        }
                                    }
                                }
                            }
                        }
                        //
                        na = 2;
                        nw = 1;
                        for (ia = 1; ia <= 3; ia = ia + 1) {
                            &a[(1 - 1) + (1 - 1) * lda] = vab[ia - 1];
                            &a[(1 - 1) + (2 - 1) * lda] = -three * vab[ia - 1];
                            a[(2 - 1)] = -seven * vab[ia - 1];
                            a[(2 - 1) + (2 - 1) * lda] = twnone * vab[ia - 1];
                            for (ib = 1; ib <= 3; ib = ib + 1) {
                                &b[(1 - 1) + (1 - 1) * ldb] = vab[ib - 1];
                                b[(2 - 1)] = -two * vab[ib - 1];
                                for (iwr = 1; iwr <= 4; iwr = iwr + 1) {
                                    if (d1 == one && d2 == one && ca == one) {
                                        wr = vwr[iwr - 1] * &a[(1 - 1) + (1 - 1) * lda];
                                    } else {
                                        wr = vwr[iwr - 1];
                                    }
                                    wi = zero;
                                    Rlaln2(ltrans[itrans - 1], na, nw, smin, ca, a, 2, d1, d2, b, 2, wr, wi, x, 2, scale, xnorm, info);
                                    if (info < 0) {
                                        ninfo[1 - 1]++;
                                    }
                                    if (info > 0) {
                                        ninfo[2 - 1]++;
                                    }
                                    if (itrans == 1) {
                                        tmp = &a[(1 - 1) + (2 - 1) * lda];
                                        &a[(1 - 1) + (2 - 1) * lda] = a[(2 - 1)];
                                        a[(2 - 1)] = tmp;
                                    }
                                    res = abs((ca * &a[(1 - 1) + (1 - 1) * lda] - wr * d1) * x[(1 - 1)] + (ca * &a[(1 - 1) + (2 - 1) * lda]) * x[(2 - 1)] - scale * &b[(1 - 1) + (1 - 1) * ldb]);
                                    res += abs((ca * a[(2 - 1)]) * x[(1 - 1)] + (ca * a[(2 - 1) + (2 - 1) * lda] - wr * d2) * x[(2 - 1)] - scale * b[(2 - 1)]);
                                    if (info == 0) {
                                        den = max({eps * (max(abs(ca * a[(1 - 1)] - wr * d1) + abs(ca * &a[(1 - 1) + (2 - 1) * lda]), abs(ca * a[(2 - 1)]) + abs(ca * a[(2 - 1) + (2 - 1) * lda] - wr * d2)) * max(abs(x[(1 - 1)]), abs(x[(2 - 1)]))), smlnum});
                                    } else {
                                        den = max({eps * (max({smin / eps, max(abs(ca * a[(1 - 1)] - wr * d1) + abs(ca * &a[(1 - 1) + (2 - 1) * lda]), abs(ca * a[(2 - 1)]) + abs(ca * a[(2 - 1) + (2 - 1) * lda] - wr * d2))}) * max(abs(x[(1 - 1)]), abs(x[(2 - 1)]))), smlnum});
                                    }
                                    res = res / den;
                                    if (abs(x[(1 - 1)]) < unfl && abs(x[(2 - 1)]) < unfl && abs(&b[(1 - 1) + (1 - 1) * ldb]) + abs(b[(2 - 1)]) <= smlnum * (abs(ca * a[(1 - 1)] - wr * d1) + abs(ca * a[(2 - 1) * lda]) + abs(ca * a[(2 - 1)]) + abs(ca * a[(2 - 1) + (2 - 1) * lda] - wr * d2))) {
                                        res = zero;
                                    }
                                    if (scale > one) {
                                        res += one / eps;
                                    }
                                    res += abs(xnorm - max(abs(x[(1 - 1)]), abs(x[(2 - 1)]))) / max(smlnum, xnorm) / eps;
                                    if (info != 0 && info != 1) {
                                        res += one / eps;
                                    }
                                    knt++;
                                    if (res > rmax) {
                                        lmax = knt;
                                        rmax = res;
                                    }
                                }
                            }
                        }
                        //
                        na = 2;
                        nw = 2;
                        for (ia = 1; ia <= 3; ia = ia + 1) {
                            &a[(1 - 1) + (1 - 1) * lda] = vab[ia - 1] * two;
                            &a[(1 - 1) + (2 - 1) * lda] = -three * vab[ia - 1];
                            a[(2 - 1)] = -seven * vab[ia - 1];
                            a[(2 - 1) + (2 - 1) * lda] = twnone * vab[ia - 1];
                            for (ib = 1; ib <= 3; ib = ib + 1) {
                                &b[(1 - 1) + (1 - 1) * ldb] = vab[ib - 1];
                                b[(2 - 1)] = -two * vab[ib - 1];
                                &b[(1 - 1) + (2 - 1) * ldb] = four * vab[ib - 1];
                                b[(2 - 1) + (2 - 1) * ldb] = -seven * vab[ib - 1];
                                for (iwr = 1; iwr <= 4; iwr = iwr + 1) {
                                    if (d1 == one && d2 == one && ca == one) {
                                        wr = vwr[iwr - 1] * &a[(1 - 1) + (1 - 1) * lda];
                                    } else {
                                        wr = vwr[iwr - 1];
                                    }
                                    for (iwi = 1; iwi <= 4; iwi = iwi + 1) {
                                        if (d1 == one && d2 == one && ca == one) {
                                            wi = vwi[iwi - 1] * &a[(1 - 1) + (1 - 1) * lda];
                                        } else {
                                            wi = vwi[iwi - 1];
                                        }
                                        Rlaln2(ltrans[itrans - 1], na, nw, smin, ca, a, 2, d1, d2, b, 2, wr, wi, x, 2, scale, xnorm, info);
                                        if (info < 0) {
                                            ninfo[1 - 1]++;
                                        }
                                        if (info > 0) {
                                            ninfo[2 - 1]++;
                                        }
                                        if (itrans == 1) {
                                            tmp = &a[(1 - 1) + (2 - 1) * lda];
                                            &a[(1 - 1) + (2 - 1) * lda] = a[(2 - 1)];
                                            a[(2 - 1)] = tmp;
                                        }
                                        res = abs((ca * &a[(1 - 1) + (1 - 1) * lda] - wr * d1) * x[(1 - 1)] + (ca * &a[(1 - 1) + (2 - 1) * lda]) * x[(2 - 1)] + (wi * d1) * x[(2 - 1) * ldx] - scale * &b[(1 - 1) + (1 - 1) * ldb]);
                                        res += abs((ca * a[(1 - 1)] - wr * d1) * x[(2 - 1) * ldx] + (ca * &a[(1 - 1) + (2 - 1) * lda]) * x[(2 - 1) + (2 - 1) * ldx] - (wi * d1) * x[(1 - 1)] - scale * &b[(1 - 1) + (2 - 1) * ldb]);
                                        res += abs((ca * a[(2 - 1)]) * x[(1 - 1)] + (ca * a[(2 - 1) + (2 - 1) * lda] - wr * d2) * x[(2 - 1)] + (wi * d2) * x[(2 - 1) + (2 - 1) * ldx] - scale * b[(2 - 1)]);
                                        res += abs((ca * a[(2 - 1)]) * x[(2 - 1) * ldx] + (ca * a[(2 - 1) + (2 - 1) * lda] - wr * d2) * x[(2 - 1) + (2 - 1) * ldx] - (wi * d2) * x[(2 - 1)] - scale * b[(2 - 1) + (2 - 1) * ldb]);
                                        if (info == 0) {
                                            den = max({eps * (max(abs(ca * a[(1 - 1)] - wr * d1) + abs(ca * a[(2 - 1) * lda]) + abs(wi * d1), abs(ca * a[(2 - 1)]) + abs(ca * a[(2 - 1) + (2 - 1) * lda] - wr * d2) + abs(wi * d2)) * max(abs(x[(1 - 1)]) + abs(x[(2 - 1)]), abs(x[(2 - 1) * ldx]) + abs(x[(2 - 1) + (2 - 1) * ldx]))), smlnum});
                                        } else {
                                            den = max({eps * (max({smin / eps, max(abs(ca * a[(1 - 1)] - wr * d1) + abs(ca * a[(2 - 1) * lda]) + abs(wi * d1), abs(ca * a[(2 - 1)]) + abs(ca * a[(2 - 1) + (2 - 1) * lda] - wr * d2) + abs(wi * d2))}) * max(abs(x[(1 - 1)]) + abs(x[(2 - 1)]), abs(x[(2 - 1) * ldx]) + abs(x[(2 - 1) + (2 - 1) * ldx]))), smlnum});
                                        }
                                        res = res / den;
                                        if (abs(x[(1 - 1)]) < unfl && abs(x[(2 - 1)]) < unfl && abs(x[(2 - 1) * ldx]) < unfl && abs(x[(2 - 1) + (2 - 1) * ldx]) < unfl && abs(&b[(1 - 1) + (1 - 1) * ldb]) + abs(b[(2 - 1)]) <= smlnum * (abs(ca * &a[(1 - 1) + (1 - 1) * lda] - wr * d1) + abs(ca * &a[(1 - 1) + (2 - 1) * lda]) + abs(ca * a[(2 - 1)]) + abs(ca * a[(2 - 1) + (2 - 1) * lda] - wr * d2) + abs(wi * d2) + abs(wi * d1))) {
                                            res = zero;
                                        }
                                        if (scale > one) {
                                            res += one / eps;
                                        }
                                        res += abs(xnorm - max(abs(x[(1 - 1)]) + abs(x[(2 - 1) * ldx]), abs(x[(2 - 1)]) + abs(x[(2 - 1) + (2 - 1) * ldx]))) / max(smlnum, xnorm) / eps;
                                        if (info != 0 && info != 1) {
                                            res += one / eps;
                                        }
                                        knt++;
                                        if (res > rmax) {
                                            lmax = knt;
                                            rmax = res;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //
    //     End of Rget31
    //
}
