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

void Rget33(REAL &rmax, INTEGER &lmax, INTEGER &ninfo, INTEGER &knt) {
    //
    //     Get machine parameters
    //
    REAL eps = Rlamch("P");
    REAL smlnum = Rlamch("S") / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    //
    //     Set up test case parameters
    //
    REAL val[4];
    val[1 - 1] = one;
    const REAL two = 2.0;
    val[2 - 1] = one + two * eps;
    val[3 - 1] = two;
    const REAL four = 4.0;
    val[4 - 1] = two - four * eps;
    REAL vm[3];
    vm[1 - 1] = smlnum;
    vm[2 - 1] = one;
    vm[3 - 1] = bignum;
    //
    knt = 0;
    ninfo = 0;
    lmax = 0;
    const REAL zero = 0.0;
    rmax = zero;
    //
    //     Begin test loop
    //
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER i3 = 0;
    INTEGER i4 = 0;
    INTEGER im1 = 0;
    INTEGER im2 = 0;
    INTEGER im3 = 0;
    INTEGER im4 = 0;
    REAL t[2 * 2];
    INTEGER ldt = 2;
    REAL tnrm = 0.0;
    REAL t1[2 * 2];
    INTEGER ldt1 = 2;
    REAL q[2 * 2];
    INTEGER ldq = 2;
    REAL wr1 = 0.0;
    REAL wi1 = 0.0;
    REAL wr2 = 0.0;
    REAL wi2 = 0.0;
    REAL cs = 0.0;
    REAL sn = 0.0;
    INTEGER j1 = 0;
    REAL res = 0.0;
    INTEGER j2 = 0;
    REAL t2[2 * 2];
    INTEGER ldt2 = 2;
    INTEGER j3 = 0;
    REAL sum = 0.0;
    for (i1 = 1; i1 <= 4; i1 = i1 + 1) {
        for (i2 = 1; i2 <= 4; i2 = i2 + 1) {
            for (i3 = 1; i3 <= 4; i3 = i3 + 1) {
                for (i4 = 1; i4 <= 4; i4 = i4 + 1) {
                    for (im1 = 1; im1 <= 3; im1 = im1 + 1) {
                        for (im2 = 1; im2 <= 3; im2 = im2 + 1) {
                            for (im3 = 1; im3 <= 3; im3 = im3 + 1) {
                                for (im4 = 1; im4 <= 3; im4 = im4 + 1) {
                                    t[(1 - 1)] = val[i1 - 1] * vm[im1 - 1];
                                    t[(2 - 1) * ldt] = val[i2 - 1] * vm[im2 - 1];
                                    t[(2 - 1)] = -val[i3 - 1] * vm[im3 - 1];
                                    t[(2 - 1) + (2 - 1) * ldt] = val[i4 - 1] * vm[im4 - 1];
                                    tnrm = max({abs(t[(1 - 1)]), abs(t[(2 - 1) * ldt]), abs(t[(2 - 1)]), abs(t[(2 - 1) + (2 - 1) * ldt])});
                                    t1[(1 - 1)] = t[(1 - 1)];
                                    t1[(2 - 1) * ldt1] = t[(2 - 1) * ldt];
                                    t1[(2 - 1)] = t[(2 - 1)];
                                    t1[(2 - 1) + (2 - 1) * ldt1] = t[(2 - 1) + (2 - 1) * ldt];
                                    q[(1 - 1)] = one;
                                    q[(2 - 1) * ldq] = zero;
                                    q[(2 - 1)] = zero;
                                    q[(2 - 1) + (2 - 1) * ldq] = one;
                                    //
                                    Rlanv2(t[(1 - 1)], t[(2 - 1) * ldt], t[(2 - 1)], t[(2 - 1) + (2 - 1) * ldt], wr1, wi1, wr2, wi2, cs, sn);
                                    for (j1 = 1; j1 <= 2; j1 = j1 + 1) {
                                        res = q[(j1 - 1)] * cs + q[(j1 - 1) + (2 - 1) * ldq] * sn;
                                        q[(j1 - 1) + (2 - 1) * ldq] = -q[(j1 - 1)] * sn + q[(j1 - 1) + (2 - 1) * ldq] * cs;
                                        q[(j1 - 1)] = res;
                                    }
                                    //
                                    res = zero;
                                    res += abs(pow2(q[(1 - 1)]) + pow2(q[(2 - 1) * ldq]) - one) / eps;
                                    res += abs(pow2(q[(2 - 1) + (2 - 1) * ldq]) + pow2(q[(2 - 1)]) - one) / eps;
                                    res += abs(q[(1 - 1)] * q[(2 - 1)] + q[(2 - 1) * ldq] * q[(2 - 1) + (2 - 1) * ldq]) / eps;
                                    for (j1 = 1; j1 <= 2; j1 = j1 + 1) {
                                        for (j2 = 1; j2 <= 2; j2 = j2 + 1) {
                                            t2[(j1 - 1) + (j2 - 1) * ldt2] = zero;
                                            for (j3 = 1; j3 <= 2; j3 = j3 + 1) {
                                                t2[(j1 - 1) + (j2 - 1) * ldt2] += t1[(j1 - 1) + (j3 - 1) * ldt1] * q[(j3 - 1) + (j2 - 1) * ldq];
                                            }
                                        }
                                    }
                                    for (j1 = 1; j1 <= 2; j1 = j1 + 1) {
                                        for (j2 = 1; j2 <= 2; j2 = j2 + 1) {
                                            sum = t[(j1 - 1) + (j2 - 1) * ldt];
                                            for (j3 = 1; j3 <= 2; j3 = j3 + 1) {
                                                sum = sum - q[(j3 - 1) + (j1 - 1) * ldq] * t2[(j3 - 1) + (j2 - 1) * ldt2];
                                            }
                                            res += abs(sum) / eps / tnrm;
                                        }
                                    }
                                    if (t[(2 - 1)] != zero && (t[(1 - 1)] != t[(2 - 1) + (2 - 1) * ldt] || sign(one, t[(2 - 1) * ldt]) * sign(one, t[(2 - 1)]) > zero)) {
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
    //
    //     End of Rget33
    //
}
