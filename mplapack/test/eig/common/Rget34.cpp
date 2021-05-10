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

void Rget34(REAL &rmax, INTEGER &lmax, INTEGER *ninfo, INTEGER &knt) {
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
    const REAL zero = 0.0;
    REAL val[9];
    val[1 - 1] = zero;
    val[2 - 1] = sqrt(smlnum);
    val[3 - 1] = one;
    const REAL two = 2.0;
    val[4 - 1] = two;
    val[5 - 1] = sqrt(bignum);
    val[6 - 1] = -sqrt(smlnum);
    val[7 - 1] = -one;
    val[8 - 1] = -two;
    val[9 - 1] = -sqrt(bignum);
    REAL vm[2];
    vm[1 - 1] = one;
    vm[2 - 1] = one + two * eps;
    REAL t[4 * 4];
    INTEGER ldt = 4;
    Rcopy(16, &val[4 - 1], 0, &t[(1 - 1)], 1);
    //
    ninfo[1 - 1] = 0;
    ninfo[2 - 1] = 0;
    knt = 0;
    lmax = 0;
    rmax = zero;
    //
    //     Begin test loop
    //
    INTEGER ia = 0;
    INTEGER iam = 0;
    INTEGER ib = 0;
    INTEGER ic = 0;
    REAL tnrm = 0.0;
    REAL t1[4 * 4];
    INTEGER ldt1 = 4;
    REAL q[4 * 4];
    INTEGER ldq = 4;
    const INTEGER lwork = 32;
    REAL work[lwork];
    INTEGER info = 0;
    REAL result[2];
    REAL res = 0.0;
    for (ia = 1; ia <= 9; ia = ia + 1) {
        for (iam = 1; iam <= 2; iam = iam + 1) {
            for (ib = 1; ib <= 9; ib = ib + 1) {
                for (ic = 1; ic <= 9; ic = ic + 1) {
                    t[(1 - 1)] = val[ia - 1] * vm[iam - 1];
                    t[(2 - 1) + (2 - 1) * ldt] = val[ic - 1];
                    t[(2 - 1) * ldt] = val[ib - 1];
                    t[(2 - 1)] = zero;
                    tnrm = max({abs(t[(1 - 1)]), abs(t[(2 - 1) + (2 - 1) * ldt]), abs(t[(2 - 1) * ldt])});
                    Rcopy(16, t, 1, t1, 1);
                    Rcopy(16, &val[1 - 1], 0, q, 1);
                    Rcopy(4, &val[3 - 1], 0, q, 5);
                    Rlaexc(true, 2, t, 4, q, 4, 1, 1, 1, work, info);
                    if (info != 0) {
                        ninfo[info - 1]++;
                    }
                    Rhst01(2, 1, 2, t1, 4, t, 4, q, 4, work, lwork, result);
                    res = result[1 - 1] + result[2 - 1];
                    if (info != 0) {
                        res += one / eps;
                    }
                    if (t[(1 - 1)] != t1[(2 - 1) + (2 - 1) * ldt1]) {
                        res += one / eps;
                    }
                    if (t[(2 - 1) + (2 - 1) * ldt] != t1[(1 - 1)]) {
                        res += one / eps;
                    }
                    if (t[(2 - 1)] != zero) {
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
    INTEGER ic11 = 0;
    INTEGER ic12 = 0;
    INTEGER ic21 = 0;
    INTEGER ic22 = 0;
    for (ia = 1; ia <= 5; ia = ia + 1) {
        for (iam = 1; iam <= 2; iam = iam + 1) {
            for (ib = 1; ib <= 5; ib = ib + 1) {
                for (ic11 = 1; ic11 <= 5; ic11 = ic11 + 1) {
                    for (ic12 = 2; ic12 <= 5; ic12 = ic12 + 1) {
                        for (ic21 = 2; ic21 <= 4; ic21 = ic21 + 1) {
                            for (ic22 = -1; ic22 <= 1; ic22 = ic22 + 2) {
                                t[(1 - 1)] = val[ia - 1] * vm[iam - 1];
                                t[(2 - 1) * ldt] = val[ib - 1];
                                t[(3 - 1) * ldt] = -two * val[ib - 1];
                                t[(2 - 1)] = zero;
                                t[(2 - 1) + (2 - 1) * ldt] = val[ic11 - 1];
                                t[(2 - 1) + (3 - 1) * ldt] = val[ic12 - 1];
                                t[(3 - 1)] = zero;
                                t[(3 - 1) + (2 - 1) * ldt] = -val[ic21 - 1];
                                t[(3 - 1) + (3 - 1) * ldt] = val[ic11 - 1] * castREAL(ic22);
                                tnrm = max({abs(t[(1 - 1)]), abs(t[(2 - 1) * ldt]), abs(t[(3 - 1) * ldt]), abs(t[(2 - 1) + (2 - 1) * ldt]), abs(t[(2 - 1) + (3 - 1) * ldt]), abs(t[(3 - 1) + (2 - 1) * ldt]), abs(t[(3 - 1) + (3 - 1) * ldt])});
                                Rcopy(16, t, 1, t1, 1);
                                Rcopy(16, &val[1 - 1], 0, q, 1);
                                Rcopy(4, &val[3 - 1], 0, q, 5);
                                Rlaexc(true, 3, t, 4, q, 4, 1, 1, 2, work, info);
                                if (info != 0) {
                                    ninfo[info - 1]++;
                                }
                                Rhst01(3, 1, 3, t1, 4, t, 4, q, 4, work, lwork, result);
                                res = result[1 - 1] + result[2 - 1];
                                if (info == 0) {
                                    if (t1[(1 - 1)] != t[(3 - 1) + (3 - 1) * ldt]) {
                                        res += one / eps;
                                    }
                                    if (t[(3 - 1)] != zero) {
                                        res += one / eps;
                                    }
                                    if (t[(3 - 1) + (2 - 1) * ldt] != zero) {
                                        res += one / eps;
                                    }
                                    if (t[(2 - 1)] != 0 && (t[(1 - 1)] != t[(2 - 1) + (2 - 1) * ldt] || sign(one, t[(2 - 1) * ldt]) == sign(one, t[(2 - 1)]))) {
                                        res += one / eps;
                                    }
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
    //
    INTEGER ia11 = 0;
    INTEGER ia12 = 0;
    INTEGER ia21 = 0;
    INTEGER ia22 = 0;
    INTEGER icm = 0;
    for (ia11 = 1; ia11 <= 5; ia11 = ia11 + 1) {
        for (ia12 = 2; ia12 <= 5; ia12 = ia12 + 1) {
            for (ia21 = 2; ia21 <= 4; ia21 = ia21 + 1) {
                for (ia22 = -1; ia22 <= 1; ia22 = ia22 + 2) {
                    for (icm = 1; icm <= 2; icm = icm + 1) {
                        for (ib = 1; ib <= 5; ib = ib + 1) {
                            for (ic = 1; ic <= 5; ic = ic + 1) {
                                t[(1 - 1)] = val[ia11 - 1];
                                t[(2 - 1) * ldt] = val[ia12 - 1];
                                t[(3 - 1) * ldt] = -two * val[ib - 1];
                                t[(2 - 1)] = -val[ia21 - 1];
                                t[(2 - 1) + (2 - 1) * ldt] = val[ia11 - 1] * castREAL(ia22);
                                t[(2 - 1) + (3 - 1) * ldt] = val[ib - 1];
                                t[(3 - 1)] = zero;
                                t[(3 - 1) + (2 - 1) * ldt] = zero;
                                t[(3 - 1) + (3 - 1) * ldt] = val[ic - 1] * vm[icm - 1];
                                tnrm = max({abs(t[(1 - 1)]), abs(t[(2 - 1) * ldt]), abs(t[(3 - 1) * ldt]), abs(t[(2 - 1) + (2 - 1) * ldt]), abs(t[(2 - 1) + (3 - 1) * ldt]), abs(t[(3 - 1) + (2 - 1) * ldt]), abs(t[(3 - 1) + (3 - 1) * ldt])});
                                Rcopy(16, t, 1, t1, 1);
                                Rcopy(16, &val[1 - 1], 0, q, 1);
                                Rcopy(4, &val[3 - 1], 0, q, 5);
                                Rlaexc(true, 3, t, 4, q, 4, 1, 2, 1, work, info);
                                if (info != 0) {
                                    ninfo[info - 1]++;
                                }
                                Rhst01(3, 1, 3, t1, 4, t, 4, q, 4, work, lwork, result);
                                res = result[1 - 1] + result[2 - 1];
                                if (info == 0) {
                                    if (t1[(3 - 1) + (3 - 1) * ldt1] != t[(1 - 1)]) {
                                        res += one / eps;
                                    }
                                    if (t[(2 - 1)] != zero) {
                                        res += one / eps;
                                    }
                                    if (t[(3 - 1)] != zero) {
                                        res += one / eps;
                                    }
                                    if (t[(3 - 1) + (2 - 1) * ldt] != 0 && (t[(2 - 1) + (2 - 1) * ldt] != t[(3 - 1) + (3 - 1) * ldt] || sign(one, t[(2 - 1) + (3 - 1) * ldt]) == sign(one, t[(3 - 1) + (2 - 1) * ldt]))) {
                                        res += one / eps;
                                    }
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
    //
    const REAL half = 0.5e0;
    const REAL three = 3.0;
    INTEGER i = 0;
    INTEGER j = 0;
    for (ia11 = 1; ia11 <= 5; ia11 = ia11 + 1) {
        for (ia12 = 2; ia12 <= 5; ia12 = ia12 + 1) {
            for (ia21 = 2; ia21 <= 4; ia21 = ia21 + 1) {
                for (ia22 = -1; ia22 <= 1; ia22 = ia22 + 2) {
                    for (ib = 1; ib <= 5; ib = ib + 1) {
                        for (ic11 = 3; ic11 <= 4; ic11 = ic11 + 1) {
                            for (ic12 = 3; ic12 <= 4; ic12 = ic12 + 1) {
                                for (ic21 = 3; ic21 <= 4; ic21 = ic21 + 1) {
                                    for (ic22 = -1; ic22 <= 1; ic22 = ic22 + 2) {
                                        for (icm = 5; icm <= 7; icm = icm + 1) {
                                            iam = 1;
                                            t[(1 - 1)] = val[ia11 - 1] * vm[iam - 1];
                                            t[(2 - 1) * ldt] = val[ia12 - 1] * vm[iam - 1];
                                            t[(3 - 1) * ldt] = -two * val[ib - 1];
                                            t[(4 - 1) * ldt] = half * val[ib - 1];
                                            t[(2 - 1)] = -t[(2 - 1) * ldt] * val[ia21 - 1];
                                            t[(2 - 1) + (2 - 1) * ldt] = val[ia11 - 1] * castREAL(ia22) * vm[iam - 1];
                                            t[(2 - 1) + (3 - 1) * ldt] = val[ib - 1];
                                            t[(2 - 1) + (4 - 1) * ldt] = three * val[ib - 1];
                                            t[(3 - 1)] = zero;
                                            t[(3 - 1) + (2 - 1) * ldt] = zero;
                                            t[(3 - 1) + (3 - 1) * ldt] = val[ic11 - 1] * abs(val[icm - 1]);
                                            t[(3 - 1) + (4 - 1) * ldt] = val[ic12 - 1] * abs(val[icm - 1]);
                                            t[(4 - 1)] = zero;
                                            t[(4 - 1) + (2 - 1) * ldt] = zero;
                                            t[(4 - 1) + (3 - 1) * ldt] = -t[(3 - 1) + (4 - 1) * ldt] * val[ic21 - 1] * abs(val[icm - 1]);
                                            t[(4 - 1) + (4 - 1) * ldt] = val[ic11 - 1] * castREAL(ic22) * abs(val[icm - 1]);
                                            tnrm = zero;
                                            for (i = 1; i <= 4; i = i + 1) {
                                                for (j = 1; j <= 4; j = j + 1) {
                                                    tnrm = max(tnrm, abs(t[(i - 1) + (j - 1) * ldt]));
                                                }
                                            }
                                            Rcopy(16, t, 1, t1, 1);
                                            Rcopy(16, &val[1 - 1], 0, q, 1);
                                            Rcopy(4, &val[3 - 1], 0, q, 5);
                                            Rlaexc(true, 4, t, 4, q, 4, 1, 2, 2, work, info);
                                            if (info != 0) {
                                                ninfo[info - 1]++;
                                            }
                                            Rhst01(4, 1, 4, t1, 4, t, 4, q, 4, work, lwork, result);
                                            res = result[1 - 1] + result[2 - 1];
                                            if (info == 0) {
                                                if (t[(3 - 1)] != zero) {
                                                    res += one / eps;
                                                }
                                                if (t[(4 - 1)] != zero) {
                                                    res += one / eps;
                                                }
                                                if (t[(3 - 1) + (2 - 1) * ldt] != zero) {
                                                    res += one / eps;
                                                }
                                                if (t[(4 - 1) + (2 - 1) * ldt] != zero) {
                                                    res += one / eps;
                                                }
                                                if (t[(2 - 1)] != 0 && (t[(1 - 1)] != t[(2 - 1) + (2 - 1) * ldt] || sign(one, t[(2 - 1) * ldt]) == sign(one, t[(2 - 1)]))) {
                                                    res += one / eps;
                                                }
                                                if (t[(4 - 1) + (3 - 1) * ldt] != 0 && (t[(3 - 1) + (3 - 1) * ldt] != t[(4 - 1) + (4 - 1) * ldt] || sign(one, t[(3 - 1) + (4 - 1) * ldt]) == sign(one, t[(4 - 1) + (3 - 1) * ldt]))) {
                                                    res += one / eps;
                                                }
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
    }
    //
    //     End of Rget34
    //
}
