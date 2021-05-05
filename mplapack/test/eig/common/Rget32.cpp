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

void Rget32(REAL &rmax, INTEGER &lmax, INTEGER &ninfo, INTEGER &knt) {
    FEM_CMN_SVE(Rget32);
    // SAVE
    arr_ref<int, 3> itval(sve.itval, [2 * 2 * 8]);
    //
    if (is_called_first_time) {
        static const INTEGER values[] = {8, 4, 2, 1, 4, 8, 1, 2, 2, 1, 8, 4, 1, 2, 4, 8, 9, 4, 2, 1, 4, 9, 1, 2, 2, 1, 9, 4, 1, 2, 4, 9};
        data_of_type<int>(FEM_VALUES_AND_SIZE), itval;
    }
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
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
    REAL smlnum = Rlamch("S") / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     Set up test case parameters
    //
    REAL val[3];
    val[1 - 1] = sqrt(smlnum);
    val[2 - 1] = one;
    val[3 - 1] = sqrt(bignum);
    //
    knt = 0;
    ninfo = 0;
    lmax = 0;
    const REAL zero = 0.0;
    rmax = zero;
    //
    //     Begin test loop
    //
    INTEGER itranl = 0;
    INTEGER itranr = 0;
    INTEGER isgn = 0;
    REAL sgn = 0.0;
    bool ltranl = false;
    bool ltranr = false;
    INTEGER n1 = 0;
    INTEGER n2 = 0;
    INTEGER itl = 0;
    INTEGER itr = 0;
    INTEGER ib = 0;
    REAL tl[2 * 2];
    REAL tr[2 * 2];
    INTEGER &b[(1 - 1) + (1 - 1) * ldb] = 0;
    REAL b[2 * 2];
    REAL scale = 0.0;
    REAL x[2 * 2];
    REAL xnorm = 0.0;
    INTEGER info = 0;
    REAL res = 0.0;
    REAL den = 0.0;
    INTEGER itlscl = 0;
    INTEGER ib1 = 0;
    INTEGER ib2 = 0;
    const REAL four = 4.0;
    REAL tmp = 0.0;
    REAL tnrm = 0.0;
    REAL xnrm = 0.0;
    INTEGER itrscl = 0;
    const REAL two = 2.0;
    INTEGER &b[(1 - 1) + (2 - 1) * ldb] = 0;
    INTEGER ib3 = 0;
    const REAL eight = 8.0;
    for (itranl = 0; itranl <= 1; itranl = itranl + 1) {
        for (itranr = 0; itranr <= 1; itranr = itranr + 1) {
            for (isgn = -1; isgn <= 1; isgn = isgn + 2) {
                sgn = isgn;
                ltranl = itranl == 1;
                ltranr = itranr == 1;
                //
                n1 = 1;
                n2 = 1;
                for (itl = 1; itl <= 3; itl = itl + 1) {
                    for (itr = 1; itr <= 3; itr = itr + 1) {
                        for (ib = 1; ib <= 3; ib = ib + 1) {
                            tl[(1 - 1)] = val[itl - 1];
                            tr[(1 - 1)] = val[itr - 1];
                            &b[(1 - 1) + (1 - 1) * ldb] = val[ib - 1];
                            knt++;
                            Rlasy2(ltranl, ltranr, isgn, n1, n2, tl, 2, tr, 2, b, 2, scale, x, 2, xnorm, info);
                            if (info != 0) {
                                ninfo++;
                            }
                            res = abs((tl[(1 - 1)] + sgn * tr[(1 - 1)]) * x[(1 - 1)] - scale * &b[(1 - 1) + (1 - 1) * ldb]);
                            if (info == 0) {
                                den = max(eps * ((abs(tr[(1 - 1)]) + abs(tl[(1 - 1)])) * abs(x[(1 - 1)])), smlnum);
                            } else {
                                den = smlnum * max(abs(x[(1 - 1)]), one);
                            }
                            res = res / den;
                            if (scale > one) {
                                res += one / eps;
                            }
                            res += abs(xnorm - abs(x[(1 - 1)])) / max(smlnum, xnorm) / eps;
                            if (info != 0 && info != 1) {
                                res += one / eps;
                            }
                            if (res > rmax) {
                                lmax = knt;
                                rmax = res;
                            }
                        }
                    }
                }
                //
                n1 = 2;
                n2 = 1;
                for (itl = 1; itl <= 8; itl = itl + 1) {
                    for (itlscl = 1; itlscl <= 3; itlscl = itlscl + 1) {
                        for (itr = 1; itr <= 3; itr = itr + 1) {
                            for (ib1 = 1; ib1 <= 3; ib1 = ib1 + 1) {
                                for (ib2 = 1; ib2 <= 3; ib2 = ib2 + 1) {
                                    &b[(1 - 1) + (1 - 1) * ldb] = val[ib1 - 1];
                                    b[(2 - 1)] = -four * val[ib2 - 1];
                                    tl[(1 - 1)] = itval1, 1, itl * val[itlscl - 1];
                                    tl[(2 - 1)] = itval2, 1, itl * val[itlscl - 1];
                                    tl[(2 - 1) * ldtl] = itval1, 2, itl * val[itlscl - 1];
                                    tl[(2 - 1) + (2 - 1) * ldtl] = itval2, 2, itl * val[itlscl - 1];
                                    tr[(1 - 1)] = val[itr - 1];
                                    knt++;
                                    Rlasy2(ltranl, ltranr, isgn, n1, n2, tl, 2, tr, 2, b, 2, scale, x, 2, xnorm, info);
                                    if (info != 0) {
                                        ninfo++;
                                    }
                                    if (ltranl) {
                                        tmp = tl[(2 - 1) * ldtl];
                                        tl[(2 - 1) * ldtl] = tl[(2 - 1)];
                                        tl[(2 - 1)] = tmp;
                                    }
                                    res = abs((tl[(1 - 1)] + sgn * tr[(1 - 1)]) * x[(1 - 1)] + tl[(2 - 1) * ldtl] * x[(2 - 1)] - scale * &b[(1 - 1) + (1 - 1) * ldb]);
                                    res += abs((tl[(2 - 1) + (2 - 1) * ldtl] + sgn * tr[(1 - 1)]) * x[(2 - 1)] + tl[(2 - 1)] * x[(1 - 1)] - scale * b[(2 - 1)]);
                                    tnrm = abs(tr[(1 - 1)]) + abs(tl[(1 - 1)]) + abs(tl[(2 - 1) * ldtl]) + abs(tl[(2 - 1)]) + abs(tl[(2 - 1) + (2 - 1) * ldtl]);
                                    xnrm = max(abs(x[(1 - 1)]), abs(x[(2 - 1)]));
                                    den = max({smlnum, smlnum * xnrm, (tnrm * eps) * xnrm});
                                    res = res / den;
                                    if (scale > one) {
                                        res += one / eps;
                                    }
                                    res += abs(xnorm - xnrm) / max(smlnum, xnorm) / eps;
                                    if (res > rmax) {
                                        lmax = knt;
                                        rmax = res;
                                    }
                                }
                            }
                        }
                    }
                }
                //
                n1 = 1;
                n2 = 2;
                for (itr = 1; itr <= 8; itr = itr + 1) {
                    for (itrscl = 1; itrscl <= 3; itrscl = itrscl + 1) {
                        for (itl = 1; itl <= 3; itl = itl + 1) {
                            for (ib1 = 1; ib1 <= 3; ib1 = ib1 + 1) {
                                for (ib2 = 1; ib2 <= 3; ib2 = ib2 + 1) {
                                    &b[(1 - 1) + (1 - 1) * ldb] = val[ib1 - 1];
                                    &b[(1 - 1) + (2 - 1) * ldb] = -two * val[ib2 - 1];
                                    tr[(1 - 1)] = itval1, 1, itr * val[itrscl - 1];
                                    tr[(2 - 1)] = itval2, 1, itr * val[itrscl - 1];
                                    tr[(2 - 1) * ldtr] = itval1, 2, itr * val[itrscl - 1];
                                    tr[(2 - 1) + (2 - 1) * ldtr] = itval2, 2, itr * val[itrscl - 1];
                                    tl[(1 - 1)] = val[itl - 1];
                                    knt++;
                                    Rlasy2(ltranl, ltranr, isgn, n1, n2, tl, 2, tr, 2, b, 2, scale, x, 2, xnorm, info);
                                    if (info != 0) {
                                        ninfo++;
                                    }
                                    if (ltranr) {
                                        tmp = tr[(2 - 1) * ldtr];
                                        tr[(2 - 1) * ldtr] = tr[(2 - 1)];
                                        tr[(2 - 1)] = tmp;
                                    }
                                    tnrm = abs(tl[(1 - 1)]) + abs(tr[(1 - 1)]) + abs(tr[(2 - 1) * ldtr]) + abs(tr[(2 - 1) + (2 - 1) * ldtr]) + abs(tr[(2 - 1)]);
                                    xnrm = abs(x[(1 - 1)]) + abs(x[(2 - 1) * ldx]);
                                    res = abs(((tl[(1 - 1)] + sgn * tr[(1 - 1)])) * (x[(1 - 1)]) + (sgn * tr[(2 - 1)]) * (x[(2 - 1) * ldx]) - (scale * &b[(1 - 1) + (1 - 1) * ldb]));
                                    res += abs(((tl[(1 - 1)] + sgn * tr[(2 - 1) + (2 - 1) * ldtr])) * (x[(2 - 1) * ldx]) + (sgn * tr[(2 - 1) * ldtr]) * (x[(1 - 1)]) - (scale * &b[(1 - 1) + (2 - 1) * ldb]));
                                    den = max({smlnum, smlnum * xnrm, (tnrm * eps) * xnrm});
                                    res = res / den;
                                    if (scale > one) {
                                        res += one / eps;
                                    }
                                    res += abs(xnorm - xnrm) / max(smlnum, xnorm) / eps;
                                    if (res > rmax) {
                                        lmax = knt;
                                        rmax = res;
                                    }
                                }
                            }
                        }
                    }
                }
                //
                n1 = 2;
                n2 = 2;
                for (itr = 1; itr <= 8; itr = itr + 1) {
                    for (itrscl = 1; itrscl <= 3; itrscl = itrscl + 1) {
                        for (itl = 1; itl <= 8; itl = itl + 1) {
                            for (itlscl = 1; itlscl <= 3; itlscl = itlscl + 1) {
                                for (ib1 = 1; ib1 <= 3; ib1 = ib1 + 1) {
                                    for (ib2 = 1; ib2 <= 3; ib2 = ib2 + 1) {
                                        for (ib3 = 1; ib3 <= 3; ib3 = ib3 + 1) {
                                            &b[(1 - 1) + (1 - 1) * ldb] = val[ib1 - 1];
                                            b[(2 - 1)] = -four * val[ib2 - 1];
                                            &b[(1 - 1) + (2 - 1) * ldb] = -two * val[ib3 - 1];
                                            b[(2 - 1) + (2 - 1) * ldb] = eight * min({val[ib1 - 1], val[ib2 - 1], val[ib3 - 1]});
                                            tr[(1 - 1)] = itval1, 1, itr * val[itrscl - 1];
                                            tr[(2 - 1)] = itval2, 1, itr * val[itrscl - 1];
                                            tr[(2 - 1) * ldtr] = itval1, 2, itr * val[itrscl - 1];
                                            tr[(2 - 1) + (2 - 1) * ldtr] = itval2, 2, itr * val[itrscl - 1];
                                            tl[(1 - 1)] = itval1, 1, itl * val[itlscl - 1];
                                            tl[(2 - 1)] = itval2, 1, itl * val[itlscl - 1];
                                            tl[(2 - 1) * ldtl] = itval1, 2, itl * val[itlscl - 1];
                                            tl[(2 - 1) + (2 - 1) * ldtl] = itval2, 2, itl * val[itlscl - 1];
                                            knt++;
                                            Rlasy2(ltranl, ltranr, isgn, n1, n2, tl, 2, tr, 2, b, 2, scale, x, 2, xnorm, info);
                                            if (info != 0) {
                                                ninfo++;
                                            }
                                            if (ltranr) {
                                                tmp = tr[(2 - 1) * ldtr];
                                                tr[(2 - 1) * ldtr] = tr[(2 - 1)];
                                                tr[(2 - 1)] = tmp;
                                            }
                                            if (ltranl) {
                                                tmp = tl[(2 - 1) * ldtl];
                                                tl[(2 - 1) * ldtl] = tl[(2 - 1)];
                                                tl[(2 - 1)] = tmp;
                                            }
                                            tnrm = abs(tr[(1 - 1)]) + abs(tr[(2 - 1)]) + abs(tr[(2 - 1) * ldtr]) + abs(tr[(2 - 1) + (2 - 1) * ldtr]) + abs(tl[(1 - 1)]) + abs(tl[(2 - 1)]) + abs(tl[(2 - 1) * ldtl]) + abs(tl[(2 - 1) + (2 - 1) * ldtl]);
                                            xnrm = max(abs(x[(1 - 1)]) + abs(x[(2 - 1) * ldx]), abs(x[(2 - 1)]) + abs(x[(2 - 1) + (2 - 1) * ldx]));
                                            res = abs(((tl[(1 - 1)] + sgn * tr[(1 - 1)])) * (x[(1 - 1)]) + (sgn * tr[(2 - 1)]) * (x[(2 - 1) * ldx]) + (tl[(2 - 1) * ldtl]) * (x[(2 - 1)]) - (scale * &b[(1 - 1) + (1 - 1) * ldb]));
                                            res += abs((tl[(1 - 1)]) * (x[(2 - 1) * ldx]) + (sgn * tr[(2 - 1) * ldtr]) * (x[(1 - 1)]) + (sgn * tr[(2 - 1) + (2 - 1) * ldtr]) * (x[(2 - 1) * ldx]) + (tl[(2 - 1) * ldtl]) * (x[(2 - 1) + (2 - 1) * ldx]) - (scale * &b[(1 - 1) + (2 - 1) * ldb]));
                                            res += abs((tl[(2 - 1)]) * (x[(1 - 1)]) + (sgn * tr[(1 - 1)]) * (x[(2 - 1)]) + (sgn * tr[(2 - 1)]) * (x[(2 - 1) + (2 - 1) * ldx]) + (tl[(2 - 1) + (2 - 1) * ldtl]) * (x[(2 - 1)]) - (scale * b[(2 - 1)]));
                                            res += abs(((tl[(2 - 1) + (2 - 1) * ldtl] + sgn * tr[(2 - 1) + (2 - 1) * ldtr])) * (x[(2 - 1) + (2 - 1) * ldx]) + (sgn * tr[(2 - 1) * ldtr]) * (x[(2 - 1)]) + (tl[(2 - 1)]) * (x[(2 - 1) * ldx]) - (scale * b[(2 - 1) + (2 - 1) * ldb]));
                                            den = max({smlnum, smlnum * xnrm, (tnrm * eps) * xnrm});
                                            res = res / den;
                                            if (scale > one) {
                                                res += one / eps;
                                            }
                                            res += abs(xnorm - xnrm) / max(smlnum, xnorm) / eps;
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
    }
    //
    //     End of Rget32
    //
}
