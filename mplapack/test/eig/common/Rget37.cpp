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

void Rget37(REAL *rmax, INTEGER *lmax, INTEGER *ninfo, INTEGER &knt, INTEGER const nin) {
    rmax([3]);
    lmax([3]);
    ninfo([3]);
    common_read read(cmn);
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    const REAL epsin = 5.9605e-8;
    const REAL zero = 0.0;
    REAL val[3];
    INTEGER n = 0;
    INTEGER i = 0;
    const INTEGER ldt = 20;
    REAL tmp[ldt * ldt];
    INTEGER j = 0;
    REAL wrin[ldt];
    REAL wiin[ldt];
    REAL sin[ldt];
    REAL sepin[ldt];
    const INTEGER lwork = 2 * ldt * (10 + ldt);
    REAL work[lwork];
    REAL tnrm = 0.0;
    INTEGER iscl = 0;
    REAL t[ldt * ldt];
    REAL vmul = 0.0;
    INTEGER info = 0;
    REAL wr[ldt];
    REAL wi[ldt];
    REAL dum[1];
    bool select[ldt];
    REAL le[ldt * ldt];
    REAL re[ldt * ldt];
    INTEGER m = 0;
    REAL s[ldt];
    REAL sep[ldt];
    INTEGER iwork[2 * ldt];
    REAL wrtmp[ldt];
    REAL witmp[ldt];
    REAL stmp[ldt];
    REAL septmp[ldt];
    INTEGER kmin = 0;
    REAL vrmin = 0.0;
    REAL vimin = 0.0;
    const REAL two = 2.0;
    REAL v = 0.0;
    REAL tol = 0.0;
    REAL tolin = 0.0;
    REAL vmax = 0.0;
    INTEGER lcmp[3];
    INTEGER ifnd = 0;
    INTEGER icmp = 0;
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
    //
    val[1 - 1] = sqrt(smlnum);
    val[2 - 1] = one;
    val[3 - 1] = sqrt(bignum);
//
//     Read input data until N=0.  Assume input eigenvalues are sorted
//     lexicographically (increasing by real part, then decreasing by
//     imaginary part)
//
statement_10:
    read(nin, star), n;
    if (n == 0) {
        return;
    }
    for (i = 1; i <= n; i = i + 1) {
        {
            read_loop rloop(cmn, nin, star);
            for (j = 1; j <= n; j = j + 1) {
                rloop, tmp(i, j);
            }
        }
    }
    for (i = 1; i <= n; i = i + 1) {
        read(nin, star), wrin(i), wiin(i), sin(i), sepin(i);
    }
    tnrm = Rlange("M", n, n, tmp, ldt, work);
    //
    //     Begin test
    //
    for (iscl = 1; iscl <= 3; iscl = iscl + 1) {
        //
        //        Scale input matrix
        //
        knt++;
        Rlacpy("F", n, n, tmp, ldt, t, ldt);
        vmul = val[iscl - 1];
        for (i = 1; i <= n; i = i + 1) {
            Rscal(n, vmul, &t[(i - 1) * ldt], 1);
        }
        if (tnrm == zero) {
            vmul = one;
        }
        //
        //        Compute eigenvalues and eigenvectors
        //
        Rgehrd(n, 1, n, t, ldt, &work[1 - 1], &work[(n + 1) - 1], lwork - n, info);
        if (info != 0) {
            lmax[1 - 1] = knt;
            ninfo[1 - 1]++;
            goto statement_240;
        }
        for (j = 1; j <= n - 2; j = j + 1) {
            for (i = j + 2; i <= n; i = i + 1) {
                t[(i - 1) + (j - 1) * ldt] = zero;
            }
        }
        //
        //        Compute Schur form
        //
        Rhseqr("S", "N", n, 1, n, t, ldt, wr, wi, dum, 1, work, lwork, info);
        if (info != 0) {
            lmax[2 - 1] = knt;
            ninfo[2 - 1]++;
            goto statement_240;
        }
        //
        //        Compute eigenvectors
        //
        Rtrevc("Both", "All", select, n, t, ldt, le, ldt, re, ldt, n, m, work, info);
        //
        //        Compute condition numbers
        //
        Rtrsna("Both", "All", select, n, t, ldt, le, ldt, re, ldt, s, sep, n, m, work, n, iwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_240;
        }
        //
        //        Sort eigenvalues and condition numbers lexicographically
        //        to compare with inputs
        //
        Rcopy(n, wr, 1, wrtmp, 1);
        Rcopy(n, wi, 1, witmp, 1);
        Rcopy(n, s, 1, stmp, 1);
        Rcopy(n, sep, 1, septmp, 1);
        Rscal(n, one / vmul, septmp, 1);
        for (i = 1; i <= n - 1; i = i + 1) {
            kmin = i;
            vrmin = wrtmp[i - 1];
            vimin = witmp[i - 1];
            for (j = i + 1; j <= n; j = j + 1) {
                if (wrtmp[j - 1] < vrmin) {
                    kmin = j;
                    vrmin = wrtmp[j - 1];
                    vimin = witmp[j - 1];
                }
            }
            wrtmp[kmin - 1] = wrtmp[i - 1];
            witmp[kmin - 1] = witmp[i - 1];
            wrtmp[i - 1] = vrmin;
            witmp[i - 1] = vimin;
            vrmin = stmp[kmin - 1];
            stmp[kmin - 1] = stmp[i - 1];
            stmp[i - 1] = vrmin;
            vrmin = septmp[kmin - 1];
            septmp[kmin - 1] = septmp[i - 1];
            septmp[i - 1] = vrmin;
        }
        //
        //        Compare condition numbers for eigenvalues
        //        taking their condition numbers into account
        //
        v = max(two * n.real() * eps * tnrm, smlnum);
        if (tnrm == zero) {
            v = one;
        }
        for (i = 1; i <= n; i = i + 1) {
            if (v > septmp[i - 1]) {
                tol = one;
            } else {
                tol = v / septmp[i - 1];
            }
            if (v > sepin[i - 1]) {
                tolin = one;
            } else {
                tolin = v / sepin[i - 1];
            }
            tol = max(tol, smlnum / eps);
            tolin = max(tolin, smlnum / eps);
            if (eps * (sin(i) - tolin) > stmp[i - 1] + tol) {
                vmax = one / eps;
            } else if (sin(i) - tolin > stmp[i - 1] + tol) {
                vmax = (sin(i) - tolin) / (stmp[i - 1] + tol);
            } else if (sin(i) + tolin < eps * (stmp[i - 1] - tol)) {
                vmax = one / eps;
            } else if (sin(i) + tolin < stmp[i - 1] - tol) {
                vmax = (stmp[i - 1] - tol) / (sin(i) + tolin);
            } else {
                vmax = one;
            }
            if (vmax > rmax[2 - 1]) {
                rmax[2 - 1] = vmax;
                if (ninfo[2 - 1] == 0) {
                    lmax[2 - 1] = knt;
                }
            }
        }
        //
        //        Compare condition numbers for eigenvectors
        //        taking their condition numbers into account
        //
        for (i = 1; i <= n; i = i + 1) {
            if (v > septmp[i - 1] * stmp[i - 1]) {
                tol = septmp[i - 1];
            } else {
                tol = v / stmp[i - 1];
            }
            if (v > sepin[i - 1] * sin(i)) {
                tolin = sepin[i - 1];
            } else {
                tolin = v / sin(i);
            }
            tol = max(tol, smlnum / eps);
            tolin = max(tolin, smlnum / eps);
            if (eps * (sepin[i - 1] - tolin) > septmp[i - 1] + tol) {
                vmax = one / eps;
            } else if (sepin[i - 1] - tolin > septmp[i - 1] + tol) {
                vmax = (sepin[i - 1] - tolin) / (septmp[i - 1] + tol);
            } else if (sepin[i - 1] + tolin < eps * (septmp[i - 1] - tol)) {
                vmax = one / eps;
            } else if (sepin[i - 1] + tolin < septmp[i - 1] - tol) {
                vmax = (septmp[i - 1] - tol) / (sepin[i - 1] + tolin);
            } else {
                vmax = one;
            }
            if (vmax > rmax[2 - 1]) {
                rmax[2 - 1] = vmax;
                if (ninfo[2 - 1] == 0) {
                    lmax[2 - 1] = knt;
                }
            }
        }
        //
        //        Compare condition numbers for eigenvalues
        //        without taking their condition numbers into account
        //
        for (i = 1; i <= n; i = i + 1) {
            if (sin(i) <= (2 * n).real() * eps && stmp[i - 1] <= (2 * n).real() * eps) {
                vmax = one;
            } else if (eps * sin(i) > stmp[i - 1]) {
                vmax = one / eps;
            } else if (sin(i) > stmp[i - 1]) {
                vmax = sin(i) / stmp[i - 1];
            } else if (sin(i) < eps * stmp[i - 1]) {
                vmax = one / eps;
            } else if (sin(i) < stmp[i - 1]) {
                vmax = stmp[i - 1] / sin(i);
            } else {
                vmax = one;
            }
            if (vmax > rmax[3 - 1]) {
                rmax[3 - 1] = vmax;
                if (ninfo[3 - 1] == 0) {
                    lmax[3 - 1] = knt;
                }
            }
        }
        //
        //        Compare condition numbers for eigenvectors
        //        without taking their condition numbers into account
        //
        for (i = 1; i <= n; i = i + 1) {
            if (sepin[i - 1] <= v && septmp[i - 1] <= v) {
                vmax = one;
            } else if (eps * sepin[i - 1] > septmp[i - 1]) {
                vmax = one / eps;
            } else if (sepin[i - 1] > septmp[i - 1]) {
                vmax = sepin[i - 1] / septmp[i - 1];
            } else if (sepin[i - 1] < eps * septmp[i - 1]) {
                vmax = one / eps;
            } else if (sepin[i - 1] < septmp[i - 1]) {
                vmax = septmp[i - 1] / sepin[i - 1];
            } else {
                vmax = one;
            }
            if (vmax > rmax[3 - 1]) {
                rmax[3 - 1] = vmax;
                if (ninfo[3 - 1] == 0) {
                    lmax[3 - 1] = knt;
                }
            }
        }
        //
        //        Compute eigenvalue condition numbers only and compare
        //
        vmax = zero;
        dum[1 - 1] = -one;
        Rcopy(n, dum, 0, stmp, 1);
        Rcopy(n, dum, 0, septmp, 1);
        Rtrsna("Eigcond", "All", select, n, t, ldt, le, ldt, re, ldt, stmp, septmp, n, m, work, n, iwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_240;
        }
        for (i = 1; i <= n; i = i + 1) {
            if (stmp[i - 1] != s[i - 1]) {
                vmax = one / eps;
            }
            if (septmp[i - 1] != dum[1 - 1]) {
                vmax = one / eps;
            }
        }
        //
        //        Compute eigenvector condition numbers only and compare
        //
        Rcopy(n, dum, 0, stmp, 1);
        Rcopy(n, dum, 0, septmp, 1);
        Rtrsna("Veccond", "All", select, n, t, ldt, le, ldt, re, ldt, stmp, septmp, n, m, work, n, iwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_240;
        }
        for (i = 1; i <= n; i = i + 1) {
            if (stmp[i - 1] != dum[1 - 1]) {
                vmax = one / eps;
            }
            if (septmp[i - 1] != sep[i - 1]) {
                vmax = one / eps;
            }
        }
        //
        //        Compute all condition numbers using SELECT and compare
        //
        for (i = 1; i <= n; i = i + 1) {
            select[i - 1] = true;
        }
        Rcopy(n, dum, 0, stmp, 1);
        Rcopy(n, dum, 0, septmp, 1);
        Rtrsna("Bothcond", "Some", select, n, t, ldt, le, ldt, re, ldt, stmp, septmp, n, m, work, n, iwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_240;
        }
        for (i = 1; i <= n; i = i + 1) {
            if (septmp[i - 1] != sep[i - 1]) {
                vmax = one / eps;
            }
            if (stmp[i - 1] != s[i - 1]) {
                vmax = one / eps;
            }
        }
        //
        //        Compute eigenvalue condition numbers using SELECT and compare
        //
        Rcopy(n, dum, 0, stmp, 1);
        Rcopy(n, dum, 0, septmp, 1);
        Rtrsna("Eigcond", "Some", select, n, t, ldt, le, ldt, re, ldt, stmp, septmp, n, m, work, n, iwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_240;
        }
        for (i = 1; i <= n; i = i + 1) {
            if (stmp[i - 1] != s[i - 1]) {
                vmax = one / eps;
            }
            if (septmp[i - 1] != dum[1 - 1]) {
                vmax = one / eps;
            }
        }
        //
        //        Compute eigenvector condition numbers using SELECT and compare
        //
        Rcopy(n, dum, 0, stmp, 1);
        Rcopy(n, dum, 0, septmp, 1);
        Rtrsna("Veccond", "Some", select, n, t, ldt, le, ldt, re, ldt, stmp, septmp, n, m, work, n, iwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_240;
        }
        for (i = 1; i <= n; i = i + 1) {
            if (stmp[i - 1] != dum[1 - 1]) {
                vmax = one / eps;
            }
            if (septmp[i - 1] != sep[i - 1]) {
                vmax = one / eps;
            }
        }
        if (vmax > rmax[1 - 1]) {
            rmax[1 - 1] = vmax;
            if (ninfo[1 - 1] == 0) {
                lmax[1 - 1] = knt;
            }
        }
        //
        //        Select first real and first complex eigenvalue
        //
        if (wi[1 - 1] == zero) {
            lcmp[1 - 1] = 1;
            ifnd = 0;
            for (i = 2; i <= n; i = i + 1) {
                if (ifnd == 1 || wi[i - 1] == zero) {
                    select[i - 1] = false;
                } else {
                    ifnd = 1;
                    lcmp[2 - 1] = i;
                    lcmp[3 - 1] = i + 1;
                    Rcopy(n, re[(i - 1) * ldre], 1, re[(2 - 1) * ldre], 1);
                    Rcopy(n, re[((i + 1) - 1) * ldre], 1, re[(3 - 1) * ldre], 1);
                    Rcopy(n, le[(i - 1) * ldle], 1, le[(2 - 1) * ldle], 1);
                    Rcopy(n, le[((i + 1) - 1) * ldle], 1, le[(3 - 1) * ldle], 1);
                }
            }
            if (ifnd == 0) {
                icmp = 1;
            } else {
                icmp = 3;
            }
        } else {
            lcmp[1 - 1] = 1;
            lcmp[2 - 1] = 2;
            ifnd = 0;
            for (i = 3; i <= n; i = i + 1) {
                if (ifnd == 1 || wi[i - 1] != zero) {
                    select[i - 1] = false;
                } else {
                    lcmp[3 - 1] = i;
                    ifnd = 1;
                    Rcopy(n, re[(i - 1) * ldre], 1, re[(3 - 1) * ldre], 1);
                    Rcopy(n, le[(i - 1) * ldle], 1, le[(3 - 1) * ldle], 1);
                }
            }
            if (ifnd == 0) {
                icmp = 2;
            } else {
                icmp = 3;
            }
        }
        //
        //        Compute all selected condition numbers
        //
        Rcopy(icmp, dum, 0, stmp, 1);
        Rcopy(icmp, dum, 0, septmp, 1);
        Rtrsna("Bothcond", "Some", select, n, t, ldt, le, ldt, re, ldt, stmp, septmp, n, m, work, n, iwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_240;
        }
        for (i = 1; i <= icmp; i = i + 1) {
            j = lcmp[i - 1];
            if (septmp[i - 1] != sep[j - 1]) {
                vmax = one / eps;
            }
            if (stmp[i - 1] != s[j - 1]) {
                vmax = one / eps;
            }
        }
        //
        //        Compute selected eigenvalue condition numbers
        //
        Rcopy(icmp, dum, 0, stmp, 1);
        Rcopy(icmp, dum, 0, septmp, 1);
        Rtrsna("Eigcond", "Some", select, n, t, ldt, le, ldt, re, ldt, stmp, septmp, n, m, work, n, iwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_240;
        }
        for (i = 1; i <= icmp; i = i + 1) {
            j = lcmp[i - 1];
            if (stmp[i - 1] != s[j - 1]) {
                vmax = one / eps;
            }
            if (septmp[i - 1] != dum[1 - 1]) {
                vmax = one / eps;
            }
        }
        //
        //        Compute selected eigenvector condition numbers
        //
        Rcopy(icmp, dum, 0, stmp, 1);
        Rcopy(icmp, dum, 0, septmp, 1);
        Rtrsna("Veccond", "Some", select, n, t, ldt, le, ldt, re, ldt, stmp, septmp, n, m, work, n, iwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_240;
        }
        for (i = 1; i <= icmp; i = i + 1) {
            j = lcmp[i - 1];
            if (stmp[i - 1] != dum[1 - 1]) {
                vmax = one / eps;
            }
            if (septmp[i - 1] != sep[j - 1]) {
                vmax = one / eps;
            }
        }
        if (vmax > rmax[1 - 1]) {
            rmax[1 - 1] = vmax;
            if (ninfo[1 - 1] == 0) {
                lmax[1 - 1] = knt;
            }
        }
    statement_240:;
    }
    goto statement_10;
    //
    //     End of Rget37
    //
}
