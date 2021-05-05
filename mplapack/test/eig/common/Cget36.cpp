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

void Cget36(REAL &rmax, INTEGER &lmax, INTEGER &ninfo, INTEGER &knt, INTEGER const nin) {
    common_read read(cmn);
    REAL eps = 0.0;
    const REAL zero = 0.0;
    INTEGER n = 0;
    INTEGER ifst = 0;
    INTEGER ilst = 0;
    INTEGER i = 0;
    const INTEGER ldt = 10;
    COMPLEX tmp[ldt * ldt];
    INTEGER j = 0;
    COMPLEX t1[ldt * ldt];
    COMPLEX t2[ldt * ldt];
    REAL res = 0.0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    COMPLEX q[ldt * ldt];
    INTEGER info1 = 0;
    const REAL one = 1.0;
    INTEGER info2 = 0;
    COMPLEX diag[ldt];
    COMPLEX ctemp = 0.0;
    const INTEGER lwork = 2 * ldt * ldt;
    COMPLEX work[lwork];
    REAL rwork[ldt];
    REAL result[2];
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
    //     .. Executable Statements ..
    //
    eps = Rlamch("P");
    rmax = zero;
    lmax = 0;
    knt = 0;
    ninfo = 0;
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
    Clacpy("F", n, n, tmp, ldt, t1, ldt);
    Clacpy("F", n, n, tmp, ldt, t2, ldt);
    res = zero;
    //
    //     Test without accumulating Q
    //
    Claset("Full", n, n, czero, cone, q, ldt);
    Ctrexc("N", n, t1, ldt, q, ldt, ifst, ilst, info1);
    for (i = 1; i <= n; i = i + 1) {
        for (j = 1; j <= n; j = j + 1) {
            if (i == j && q[(i - 1) + (j - 1) * ldq] != cone) {
                res += one / eps;
            }
            if (i != j && q[(i - 1) + (j - 1) * ldq] != czero) {
                res += one / eps;
            }
        }
    }
    //
    //     Test with accumulating Q
    //
    Claset("Full", n, n, czero, cone, q, ldt);
    Ctrexc("V", n, t2, ldt, q, ldt, ifst, ilst, info2);
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
    if (info1 != 0 || info2 != 0) {
        ninfo++;
    }
    if (info1 != info2) {
        res += one / eps;
    }
    //
    //     Test for successful reordering of T2
    //
    Ccopy(n, tmp, ldt + 1, diag, 1);
    if (ifst < ilst) {
        for (i = ifst + 1; i <= ilst; i = i + 1) {
            ctemp = diag[i - 1];
            diag[i - 1] = diag[(i - 1) - 1];
            diag[(i - 1) - 1] = ctemp;
        }
    } else if (ifst > ilst) {
        for (i = ifst - 1; i >= ilst; i = i - 1) {
            ctemp = diag[(i + 1) - 1];
            diag[(i + 1) - 1] = diag[i - 1];
            diag[i - 1] = ctemp;
        }
    }
    for (i = 1; i <= n; i = i + 1) {
        if (t2[(i - 1) + (i - 1) * ldt2] != diag[i - 1]) {
            res += one / eps;
        }
    }
    //
    //     Test for small residual, and orthogonality of Q
    //
    Chst01(n, 1, n, tmp, ldt, t2, ldt, q, ldt, work, lwork, rwork, result);
    res += result[1 - 1] + result[2 - 1];
    //
    //     Test for T2 being in Schur form
    //
    for (j = 1; j <= n - 1; j = j + 1) {
        for (i = j + 1; i <= n; i = i + 1) {
            if (t2[(i - 1) + (j - 1) * ldt2] != czero) {
                res += one / eps;
            }
        }
    }
    if (res > rmax) {
        rmax = res;
        lmax = knt;
    }
    goto statement_10;
    //
    //     End of Cget36
    //
}
