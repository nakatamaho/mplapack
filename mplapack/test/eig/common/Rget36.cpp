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

void Rget36(REAL &rmax, INTEGER &lmax, INTEGER *ninfo, INTEGER &knt, INTEGER const nin) {
    common_read read(cmn);
    REAL eps = 0.0;
    const REAL zero = 0.0;
    INTEGER n = 0;
    INTEGER ifst = 0;
    INTEGER ilst = 0;
    INTEGER i = 0;
    const INTEGER ldt = 10;
    arr_2d<ldt, ldt, REAL> tmp;
    INTEGER j = 0;
    arr_2d<ldt, ldt, REAL> t1;
    arr_2d<ldt, ldt, REAL> t2;
    INTEGER ifstsv = 0;
    INTEGER ilstsv = 0;
    INTEGER ifst1 = 0;
    INTEGER ilst1 = 0;
    INTEGER ifst2 = 0;
    INTEGER ilst2 = 0;
    REAL res = 0.0;
    const REAL one = 1.0;
    arr_2d<ldt, ldt, REAL> q;
    const INTEGER lwork = 2 * ldt * ldt;
    arr_1d<lwork, REAL> work;
    INTEGER info1 = 0;
    INTEGER info2 = 0;
    arr_1d<2, REAL> result;
    INTEGER loc = 0;
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
    eps = Rlamch("P");
    rmax = zero;
    lmax = 0;
    knt = 0;
    ninfo[1 - 1] = 0;
    ninfo[2 - 1] = 0;
    ninfo[3 - 1] = 0;
//
//     Read input data until N=0
//
statement_10:
    read(nin, star), n, ifst, ilst;
    if (n == 0) {
        return;
    }
    knt++;
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, tmp(i, j);
            }
        }
    }
    Rlacpy("F", n, n, tmp, ldt, t1, ldt);
    Rlacpy("F", n, n, tmp, ldt, t2, ldt);
    ifstsv = ifst;
    ilstsv = ilst;
    ifst1 = ifst;
    ilst1 = ilst;
    ifst2 = ifst;
    ilst2 = ilst;
    res = zero;
    //
    //     Test without accumulating Q
    //
    Rlaset("Full", n, n, zero, one, q, ldt);
    Rtrexc("N", n, t1, ldt, q, ldt, ifst1, ilst1, work, info1);
    for (i = 1; i <= n; i = i + 1) {
        for (j = 1; j <= n; j = j + 1) {
            if (i == j && q[(i - 1) + (j - 1) * ldq] != one) {
                res += one / eps;
            }
            if (i != j && q[(i - 1) + (j - 1) * ldq] != zero) {
                res += one / eps;
            }
        }
    }
    //
    //     Test with accumulating Q
    //
    Rlaset("Full", n, n, zero, one, q, ldt);
    Rtrexc("V", n, t2, ldt, q, ldt, ifst2, ilst2, work, info2);
    //
    //     Compare T1 with T2
    //
    for (i = 1; i <= n; i = i + 1) {
        for (j = 1; j <= n; j = j + 1) {
            if (t1[(i - 1) + (j - 1) * ldt1] != t2[(i - 1) + (j - 1) * ldt2]) {
                res += one / eps;
            }
        }
    }
    if (ifst1 != ifst2) {
        res += one / eps;
    }
    if (ilst1 != ilst2) {
        res += one / eps;
    }
    if (info1 != info2) {
        res += one / eps;
    }
    //
    //     Test for successful reordering of T2
    //
    if (info2 != 0) {
        ninfo[info2 - 1]++;
    } else {
        if (abs(ifst2 - ifstsv) > 1) {
            res += one / eps;
        }
        if (abs(ilst2 - ilstsv) > 1) {
            res += one / eps;
        }
    }
    //
    //     Test for small residual, and orthogonality of Q
    //
    Rhst01(n, 1, n, tmp, ldt, t2, ldt, q, ldt, work, lwork, result);
    res += result[1 - 1] + result[2 - 1];
    //
    //     Test for T2 being in Schur form
    //
    loc = 1;
statement_70:
    if (t2[((loc + 1) - 1) + (loc - 1) * ldt2] != zero) {
        //
        //        2 by 2 block
        //
        if (t2[(loc - 1) + ((loc + 1) - 1) * ldt2] == zero || t2[(loc - 1) + (loc - 1) * ldt2] != t2[((loc + 1) - 1) + ((loc + 1) - 1) * ldt2] || sign(one, t2[(loc - 1) + ((loc + 1) - 1) * ldt2]) == sign(one, t2[((loc + 1) - 1) + (loc - 1) * ldt2])) {
            res += one / eps;
        }
        for (i = loc + 2; i <= n; i = i + 1) {
            if (t2[(i - 1) + (loc - 1) * ldt2] != zero) {
                res += one / res;
            }
            if (t2[(i - 1) + ((loc + 1) - 1) * ldt2] != zero) {
                res += one / res;
            }
        }
        loc += 2;
    } else {
        //
        //        1 by 1 block
        //
        for (i = loc + 1; i <= n; i = i + 1) {
            if (t2[(i - 1) + (loc - 1) * ldt2] != zero) {
                res += one / res;
            }
        }
        loc++;
    }
    if (loc < n) {
        goto statement_70;
    }
    if (res > rmax) {
        rmax = res;
        lmax = knt;
    }
    goto statement_10;
    //
    //     End of Rget36
    //
}
