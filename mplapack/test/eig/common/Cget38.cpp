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

void Cget38(REAL *rmax, INTEGER *lmax, INTEGER *ninfo, INTEGER &knt, INTEGER const nin) {
    common_read read(cmn);
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    const REAL epsin = 5.9605e-8;
    const REAL zero = 0.0;
    arr_1d<3, REAL> val;
    INTEGER n = 0;
    INTEGER ndim = 0;
    INTEGER isrt = 0;
    const INTEGER ldt = 20;
    arr_1d<ldt, int> iselec;
    INTEGER i = 0;
    INTEGER j = 0;
    REAL sin = 0.0;
    REAL sepin = 0.0;
    arr_1d<ldt, REAL> rwork;
    REAL tnrm = 0.0;
    INTEGER iscl = 0;
    REAL vmul = 0.0;
    const INTEGER lwork = 2 * ldt * (10 + ldt);
    INTEGER info = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    arr_1d<ldt, COMPLEX> w;
    arr_1d<ldt, int> ipnt;
    arr_1d<ldt, bool> select;
    arr_1d<ldt, REAL> wsrt;
    INTEGER kmin = 0;
    REAL vmin = 0.0;
    INTEGER itmp = 0;
    arr_1d<ldt, COMPLEX> wtmp;
    INTEGER m = 0;
    REAL s = 0.0;
    REAL sep = 0.0;
    REAL septmp = 0.0;
    REAL stmp = 0.0;
    arr_1d<2, REAL> result;
    REAL vmax = 0.0;
    const REAL two = 2.0e+0;
    REAL v = 0.0;
    REAL tol = 0.0;
    REAL tolin = 0.0;
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
    smlnum = Rlamch("S") / eps;
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    //
    //     EPSIN = 2**(-24) = precision to which input data computed
    //
    eps = max(eps, epsin);
    rmax[1 - 1] = zero;
    rmax[2 - 1] = zero;
    rmax[3 - 1] = zero;
    lmax[1 - 1] = 0;
    lmax[2 - 1] = 0;
    lmax[3 - 1] = 0;
    knt = 0;
    ninfo[1 - 1] = 0;
    ninfo[2 - 1] = 0;
    ninfo[3 - 1] = 0;
    val[1 - 1] = sqrt(smlnum);
    val[2 - 1] = one;
    val[3 - 1] = sqrt(sqrt(bignum));
//
//     Read input data until N=0.  Assume input eigenvalues are sorted
//     lexicographically (increasing by real part, then decreasing by
//     imaginary part)
//
statement_10:
    read(nin, star), n, ndim, isrt;
    if (n == 0) {
        return;
    }
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= ndim; i = i + 1) {
            rloop, iselec(i);
        }
    }
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, tmp(i, j);
            }
        }
    }
    read(nin, star), sin, sepin;
    //
    tnrm = Clange("M", n, n, tmp, ldt, rwork);
    for (iscl = 1; iscl <= 3; iscl = iscl + 1) {
        //
        //        Scale input matrix
        //
        knt++;
        Clacpy("F", n, n, tmp, ldt, t, ldt);
        vmul = val[iscl - 1];
        for (i = 1; i <= n; i = i + 1) {
            CRscal(n, vmul, &t[(i - 1) * ldt], 1);
        }
        if (tnrm == zero) {
            vmul = one;
        }
        Clacpy("F", n, n, t, ldt, tsav, ldt);
        //
        //        Compute Schur form
        //
        Cgehrd(n, 1, n, t, ldt, &work[1 - 1], &work[(n + 1) - 1], lwork - n, info);
        if (info != 0) {
            lmax[1 - 1] = knt;
            ninfo[1 - 1]++;
            goto statement_200;
        }
        //
        //        Generate unitary matrix
        //
        Clacpy("L", n, n, t, ldt, q, ldt);
        Cunghr(n, 1, n, q, ldt, &work[1 - 1], &work[(n + 1) - 1], lwork - n, info);
        //
        //        Compute Schur form
        //
        for (j = 1; j <= n - 2; j = j + 1) {
            for (i = j + 2; i <= n; i = i + 1) {
                t[(i - 1) + (j - 1) * ldt] = czero;
            }
        }
        Chseqr("S", "V", n, 1, n, t, ldt, w, q, ldt, work, lwork, info);
        if (info != 0) {
            lmax[2 - 1] = knt;
            ninfo[2 - 1]++;
            goto statement_200;
        }
        //
        //        Sort, select eigenvalues
        //
        for (i = 1; i <= n; i = i + 1) {
            ipnt[i - 1] = i;
            select[i - 1] = false;
        }
        if (isrt == 0) {
            for (i = 1; i <= n; i = i + 1) {
                wsrt[i - 1] = w[i - 1].real();
            }
        } else {
            for (i = 1; i <= n; i = i + 1) {
                wsrt[i - 1] = w[i - 1].imag();
            }
        }
        for (i = 1; i <= n - 1; i = i + 1) {
            kmin = i;
            vmin = wsrt[i - 1];
            for (j = i + 1; j <= n; j = j + 1) {
                if (wsrt[j - 1] < vmin) {
                    kmin = j;
                    vmin = wsrt[j - 1];
                }
            }
            wsrt[kmin - 1] = wsrt[i - 1];
            wsrt[i - 1] = vmin;
            itmp = ipnt[i - 1];
            ipnt[i - 1] = ipnt[kmin - 1];
            ipnt[kmin - 1] = itmp;
        }
        for (i = 1; i <= ndim; i = i + 1) {
            select[ipnt[iselec[i - 1] - 1] - 1] = true;
        }
        //
        //        Compute condition numbers
        //
        Clacpy("F", n, n, q, ldt, qsav, ldt);
        Clacpy("F", n, n, t, ldt, tsav1, ldt);
        Ctrsen("B", "V", select, n, t, ldt, q, ldt, wtmp, m, s, sep, work, lwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_200;
        }
        septmp = sep / vmul;
        stmp = s;
        //
        //        Compute residuals
        //
        Chst01(n, 1, n, tsav, ldt, t, ldt, q, ldt, work, lwork, rwork, result);
        vmax = max(result[1 - 1], result[2 - 1]);
        if (vmax > rmax[1 - 1]) {
            rmax[1 - 1] = vmax;
            if (ninfo[1 - 1] == 0) {
                lmax[1 - 1] = knt;
            }
        }
        //
        //        Compare condition number for eigenvalue cluster
        //        taking its condition number into account
        //
        v = max(two * n.real() * eps * tnrm, smlnum);
        if (tnrm == zero) {
            v = one;
        }
        if (v > septmp) {
            tol = one;
        } else {
            tol = v / septmp;
        }
        if (v > sepin) {
            tolin = one;
        } else {
            tolin = v / sepin;
        }
        tol = max(tol, smlnum / eps);
        tolin = max(tolin, smlnum / eps);
        if (eps * (sin - tolin) > stmp + tol) {
            vmax = one / eps;
        } else if (sin - tolin > stmp + tol) {
            vmax = (sin - tolin) / (stmp + tol);
        } else if (sin + tolin < eps * (stmp - tol)) {
            vmax = one / eps;
        } else if (sin + tolin < stmp - tol) {
            vmax = (stmp - tol) / (sin + tolin);
        } else {
            vmax = one;
        }
        if (vmax > rmax[2 - 1]) {
            rmax[2 - 1] = vmax;
            if (ninfo[2 - 1] == 0) {
                lmax[2 - 1] = knt;
            }
        }
        //
        //        Compare condition numbers for invariant subspace
        //        taking its condition number into account
        //
        if (v > septmp * stmp) {
            tol = septmp;
        } else {
            tol = v / stmp;
        }
        if (v > sepin * sin) {
            tolin = sepin;
        } else {
            tolin = v / sin;
        }
        tol = max(tol, smlnum / eps);
        tolin = max(tolin, smlnum / eps);
        if (eps * (sepin - tolin) > septmp + tol) {
            vmax = one / eps;
        } else if (sepin - tolin > septmp + tol) {
            vmax = (sepin - tolin) / (septmp + tol);
        } else if (sepin + tolin < eps * (septmp - tol)) {
            vmax = one / eps;
        } else if (sepin + tolin < septmp - tol) {
            vmax = (septmp - tol) / (sepin + tolin);
        } else {
            vmax = one;
        }
        if (vmax > rmax[2 - 1]) {
            rmax[2 - 1] = vmax;
            if (ninfo[2 - 1] == 0) {
                lmax[2 - 1] = knt;
            }
        }
        //
        //        Compare condition number for eigenvalue cluster
        //        without taking its condition number into account
        //
        if (sin <= (2 * n).real() * eps && stmp <= (2 * n).real() * eps) {
            vmax = one;
        } else if (eps * sin > stmp) {
            vmax = one / eps;
        } else if (sin > stmp) {
            vmax = sin / stmp;
        } else if (sin < eps * stmp) {
            vmax = one / eps;
        } else if (sin < stmp) {
            vmax = stmp / sin;
        } else {
            vmax = one;
        }
        if (vmax > rmax[3 - 1]) {
            rmax[3 - 1] = vmax;
            if (ninfo[3 - 1] == 0) {
                lmax[3 - 1] = knt;
            }
        }
        //
        //        Compare condition numbers for invariant subspace
        //        without taking its condition number into account
        //
        if (sepin <= v && septmp <= v) {
            vmax = one;
        } else if (eps * sepin > septmp) {
            vmax = one / eps;
        } else if (sepin > septmp) {
            vmax = sepin / septmp;
        } else if (sepin < eps * septmp) {
            vmax = one / eps;
        } else if (sepin < septmp) {
            vmax = septmp / sepin;
        } else {
            vmax = one;
        }
        if (vmax > rmax[3 - 1]) {
            rmax[3 - 1] = vmax;
            if (ninfo[3 - 1] == 0) {
                lmax[3 - 1] = knt;
            }
        }
        //
        //        Compute eigenvalue condition number only and compare
        //        Update Q
        //
        vmax = zero;
        Clacpy("F", n, n, tsav1, ldt, ttmp, ldt);
        Clacpy("F", n, n, qsav, ldt, qtmp, ldt);
        septmp = -one;
        stmp = -one;
        Ctrsen("E", "V", select, n, ttmp, ldt, qtmp, ldt, wtmp, m, stmp, septmp, work, lwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_200;
        }
        if (s != stmp) {
            vmax = one / eps;
        }
        if (-one != septmp) {
            vmax = one / eps;
        }
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (ttmp[(i - 1) + (j - 1) * ldttmp] != t[(i - 1) + (j - 1) * ldt]) {
                    vmax = one / eps;
                }
                if (qtmp[(i - 1) + (j - 1) * ldqtmp] != q[(i - 1) + (j - 1) * ldq]) {
                    vmax = one / eps;
                }
            }
        }
        //
        //        Compute invariant subspace condition number only and compare
        //        Update Q
        //
        Clacpy("F", n, n, tsav1, ldt, ttmp, ldt);
        Clacpy("F", n, n, qsav, ldt, qtmp, ldt);
        septmp = -one;
        stmp = -one;
        Ctrsen("V", "V", select, n, ttmp, ldt, qtmp, ldt, wtmp, m, stmp, septmp, work, lwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_200;
        }
        if (-one != stmp) {
            vmax = one / eps;
        }
        if (sep != septmp) {
            vmax = one / eps;
        }
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (ttmp[(i - 1) + (j - 1) * ldttmp] != t[(i - 1) + (j - 1) * ldt]) {
                    vmax = one / eps;
                }
                if (qtmp[(i - 1) + (j - 1) * ldqtmp] != q[(i - 1) + (j - 1) * ldq]) {
                    vmax = one / eps;
                }
            }
        }
        //
        //        Compute eigenvalue condition number only and compare
        //        Do not update Q
        //
        Clacpy("F", n, n, tsav1, ldt, ttmp, ldt);
        Clacpy("F", n, n, qsav, ldt, qtmp, ldt);
        septmp = -one;
        stmp = -one;
        Ctrsen("E", "N", select, n, ttmp, ldt, qtmp, ldt, wtmp, m, stmp, septmp, work, lwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_200;
        }
        if (s != stmp) {
            vmax = one / eps;
        }
        if (-one != septmp) {
            vmax = one / eps;
        }
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (ttmp[(i - 1) + (j - 1) * ldttmp] != t[(i - 1) + (j - 1) * ldt]) {
                    vmax = one / eps;
                }
                if (qtmp[(i - 1) + (j - 1) * ldqtmp] != qsav[(i - 1) + (j - 1) * ldqsav]) {
                    vmax = one / eps;
                }
            }
        }
        //
        //        Compute invariant subspace condition number only and compare
        //        Do not update Q
        //
        Clacpy("F", n, n, tsav1, ldt, ttmp, ldt);
        Clacpy("F", n, n, qsav, ldt, qtmp, ldt);
        septmp = -one;
        stmp = -one;
        Ctrsen("V", "N", select, n, ttmp, ldt, qtmp, ldt, wtmp, m, stmp, septmp, work, lwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_200;
        }
        if (-one != stmp) {
            vmax = one / eps;
        }
        if (sep != septmp) {
            vmax = one / eps;
        }
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                if (ttmp[(i - 1) + (j - 1) * ldttmp] != t[(i - 1) + (j - 1) * ldt]) {
                    vmax = one / eps;
                }
                if (qtmp[(i - 1) + (j - 1) * ldqtmp] != qsav[(i - 1) + (j - 1) * ldqsav]) {
                    vmax = one / eps;
                }
            }
        }
        if (vmax > rmax[1 - 1]) {
            rmax[1 - 1] = vmax;
            if (ninfo[1 - 1] == 0) {
                lmax[1 - 1] = knt;
            }
        }
    statement_200:;
    }
    goto statement_10;
    //
    //     End of Cget38
    //
}
