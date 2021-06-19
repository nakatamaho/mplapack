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

#include <mplapack_matgen.h>
#include <mplapack_eig.h>

#include <mplapack_debug.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>

using namespace std;
using std::regex;
using std::regex_replace;

void Rget38(REAL *rmax, INTEGER *lmax, INTEGER *ninfo, INTEGER &knt, INTEGER const nin) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    double dtmp;
    char buf[1024];
    REAL eps = 0.0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    const REAL epsin = 5.9605e-8;
    const REAL zero = 0.0;
    REAL val[3];
    INTEGER n = 0;
    INTEGER ndim = 0;
    const INTEGER ldt = 20;
    INTEGER iselec[ldt];
    INTEGER i = 0;
    REAL tmp[ldt * ldt];
    INTEGER ldtmp = ldt;
    INTEGER j = 0;
    REAL sin = 0.0;
    REAL sepin = 0.0;
    const INTEGER lwork = 2 * ldt * (10 + ldt);
    REAL work[lwork];
    REAL tnrm = 0.0;
    INTEGER iscl = 0;
    REAL t[ldt * ldt];
    REAL vmul = 0.0;
    REAL tsav[ldt * ldt];
    INTEGER info = 0;
    REAL q[ldt * ldt];
    REAL wr[ldt];
    REAL wi[ldt];
    INTEGER ipnt[ldt];
    bool select[ldt];
    REAL wrtmp[ldt];
    REAL witmp[ldt];
    INTEGER kmin = 0;
    REAL vrmin = 0.0;
    REAL vimin = 0.0;
    INTEGER itmp = 0;
    REAL qsav[ldt * ldt];
    REAL tsav1[ldt * ldt];
    INTEGER ldqsav = ldt;
    INTEGER ldtsav1 = ldt;
    INTEGER m = 0;
    REAL s = 0.0;
    REAL sep = 0.0;
    const INTEGER liwork = ldt * ldt;
    INTEGER iwork[liwork];
    REAL septmp = 0.0;
    REAL stmp = 0.0;
    REAL result[2];
    REAL vmax = 0.0;
    const REAL two = 2.0;
    REAL v = 0.0;
    REAL tol = 0.0;
    REAL tolin = 0.0;
    REAL ttmp[ldt * ldt];
    REAL qtmp[ldt * ldt];
    INTEGER ldttmp = ldt;
    INTEGER ldq = ldt;
    INTEGER ldqtmp = ldt;
    //
    //     .. Executable Statements ..
    //
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    bignum = one / smlnum;
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
    val[3 - 1] = sqrt(sqrt(bignum));
//
    string str;
    istringstream iss;
//
//     Read input data until N=0.  Assume input eigenvalues are sorted
//     lexicographically (increasing by real part, then decreasing by
//     imaginary part)
//
statement_10:
    getline(cin, str);
    iss.clear();
    iss.str(str);
    iss >> n;
    iss >> ndim;    
    if (n == 0) {
        return;
    }
    getline(cin, str);
    string _r = regex_replace(str, regex("D\\+"), "e+");
    str = regex_replace(_r, regex("D\\-"), "e-");
    iss.clear();
    iss.str(str);
    for (i = 1; i <= ndim; i = i + 1) {
      iss >> itmp;
      iselec[i - 1] = itmp;
    }
    for (i = 1; i <= n; i = i + 1) {
        getline(cin, str);
        _r = regex_replace(str, regex("D\\+"), "e+");
        str = regex_replace(_r, regex("D\\-"), "e-");
        iss.clear();
        iss.str(str);
        for (j = 1; j <= n; j = j + 1) {
            iss >> dtmp;
            tmp[(i - 1) + (j - 1) * ldtmp] = dtmp;
        }
    }
    getline(cin, str);
    _r = regex_replace(str, regex("D\\+"), "e+");
    str = regex_replace(_r, regex("D\\-"), "e-");
    iss.clear();
    iss.str(str);
    iss >> dtmp;  dtmp = sin;
    iss >> dtmp;  dtmp = sepin;
    //
    tnrm = Rlange("M", n, n, tmp, ldt, work);
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
        Rlacpy("F", n, n, t, ldt, tsav, ldt);
        //
        //        Compute Schur form
        //
        Rgehrd(n, 1, n, t, ldt, &work[1 - 1], &work[(n + 1) - 1], lwork - n, info);
        if (info != 0) {
            lmax[1 - 1] = knt;
            ninfo[1 - 1]++;
            goto statement_160;
        }
        //
        //        Generate orthogonal matrix
        //
        Rlacpy("L", n, n, t, ldt, q, ldt);
        Rorghr(n, 1, n, q, ldt, &work[1 - 1], &work[(n + 1) - 1], lwork - n, info);
        //
        //        Compute Schur form
        //
        Rhseqr("S", "V", n, 1, n, t, ldt, wr, wi, q, ldt, work, lwork, info);
        if (info != 0) {
            lmax[2 - 1] = knt;
            ninfo[2 - 1]++;
            goto statement_160;
        }
        //
        //        Sort, select eigenvalues
        //
        for (i = 1; i <= n; i = i + 1) {
            ipnt[i - 1] = i;
            select[i - 1] = false;
        }
        Rcopy(n, wr, 1, wrtmp, 1);
        Rcopy(n, wi, 1, witmp, 1);
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
        Rlacpy("F", n, n, q, ldt, qsav, ldt);
        Rlacpy("F", n, n, t, ldt, tsav1, ldt);
        Rtrsen("B", "V", select, n, t, ldt, q, ldt, wrtmp, witmp, m, s, sep, work, lwork, iwork, liwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_160;
        }
        septmp = sep / vmul;
        stmp = s;
        //
        //        Compute residuals
        //
        Rhst01(n, 1, n, tsav, ldt, t, ldt, q, ldt, work, lwork, result);
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
        v = max(two * castREAL(n) * eps * tnrm, smlnum);
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
        if (sin <= castREAL(2 * n) * eps && stmp <= castREAL(2 * n) * eps) {
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
        Rlacpy("F", n, n, tsav1, ldt, ttmp, ldt);
        Rlacpy("F", n, n, qsav, ldt, qtmp, ldt);
        septmp = -one;
        stmp = -one;
        Rtrsen("E", "V", select, n, ttmp, ldt, qtmp, ldt, wrtmp, witmp, m, stmp, septmp, work, lwork, iwork, liwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_160;
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
        Rlacpy("F", n, n, tsav1, ldt, ttmp, ldt);
        Rlacpy("F", n, n, qsav, ldt, qtmp, ldt);
        septmp = -one;
        stmp = -one;
        Rtrsen("V", "V", select, n, ttmp, ldt, qtmp, ldt, wrtmp, witmp, m, stmp, septmp, work, lwork, iwork, liwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_160;
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
        Rlacpy("F", n, n, tsav1, ldt, ttmp, ldt);
        Rlacpy("F", n, n, qsav, ldt, qtmp, ldt);
        septmp = -one;
        stmp = -one;
        Rtrsen("E", "N", select, n, ttmp, ldt, qtmp, ldt, wrtmp, witmp, m, stmp, septmp, work, lwork, iwork, liwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_160;
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
        Rlacpy("F", n, n, tsav1, ldt, ttmp, ldt);
        Rlacpy("F", n, n, qsav, ldt, qtmp, ldt);
        septmp = -one;
        stmp = -one;
        Rtrsen("V", "N", select, n, ttmp, ldt, qtmp, ldt, wrtmp, witmp, m, stmp, septmp, work, lwork, iwork, liwork, info);
        if (info != 0) {
            lmax[3 - 1] = knt;
            ninfo[3 - 1]++;
            goto statement_160;
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
    statement_160:;
    }
    goto statement_10;
    //
    //     End of Rget38
    //
}
