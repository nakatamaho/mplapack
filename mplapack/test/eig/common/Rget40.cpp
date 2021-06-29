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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>

using namespace std;
using std::regex;
using std::regex_replace;

void Rget40(REAL &rmax, INTEGER &lmax, INTEGER *ninfo, INTEGER &knt, INTEGER const nin) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    double dtmp;
    char buf[1024];
    REAL eps = 0.0;
    const REAL zero = 0.0;
    INTEGER n = 0;
    INTEGER ifst = 0;
    INTEGER ilst = 0;
    INTEGER i = 0;
    const INTEGER ldt = 10;
    REAL tmp[ldt * ldt];
    INTEGER j = 0;
    REAL t[ldt * ldt];
    REAL t1[ldt * ldt];
    REAL t2[ldt * ldt];
    REAL s[ldt * ldt];
    REAL s1[ldt * ldt];
    REAL s2[ldt * ldt];
    INTEGER ldtmp = ldt;
    INTEGER ldt1 = ldt;
    INTEGER ldt2 = ldt;
    INTEGER lds = ldt;
    INTEGER lds1 = ldt;
    INTEGER lds2 = ldt;
    INTEGER ifstsv = 0;
    INTEGER ilstsv = 0;
    INTEGER ifst1 = 0;
    INTEGER ilst1 = 0;
    INTEGER ifst2 = 0;
    INTEGER ilst2 = 0;
    REAL res = 0.0;
    const REAL one = 1.0;
    REAL q[ldt * ldt];
    REAL z[ldt * ldt];
    INTEGER ldq = ldt;
    INTEGER ldz = ldt;
    const INTEGER lwork = 100 + 4 * ldt + 16;
    REAL work[lwork];
    INTEGER info1 = 0;
    INTEGER info2 = 0;
    //
    REAL result[4];
    string str;
    char line[1024];
    //
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
    getline(cin, str);
    stringstream ss(str);
    ss >> n;
    ss >> ifst;
    ss >> ilst;
    if (n == 0) {
        return;
    }
    knt++;
    for (i = 1; i <= n; i = i + 1) {
        getline(cin, str);
        string _r = regex_replace(str, regex("D\\+"), "e+");
        str = regex_replace(_r, regex("D\\-"), "e-");
        istringstream iss(str);
        for (j = 1; j <= n; j = j + 1) {
            iss >> dtmp;
            tmp[(i - 1) + (j - 1) * ldtmp] = dtmp;
        }
    }
    Rlacpy("F", n, n, tmp, ldt, t, ldt);
    Rlacpy("F", n, n, tmp, ldt, t1, ldt);
    Rlacpy("F", n, n, tmp, ldt, t2, ldt);
    for (i = 1; i <= n; i = i + 1) {
        getline(cin, str);
        string _r = regex_replace(str, regex("D\\+"), "e+");
        str = regex_replace(_r, regex("D\\-"), "e-");
        istringstream iss(str);
        for (j = 1; j <= n; j = j + 1) {
            iss >> dtmp;
            tmp[(i - 1) + (j - 1) * ldtmp] = dtmp;
        }
    }
    Rlacpy("F", n, n, tmp, ldt, s, ldt);
    Rlacpy("F", n, n, tmp, ldt, s1, ldt);
    Rlacpy("F", n, n, tmp, ldt, s2, ldt);
    ifstsv = ifst;
    ilstsv = ilst;
    ifst1 = ifst;
    ilst1 = ilst;
    ifst2 = ifst;
    ilst2 = ilst;
    res = zero;
    //
    //     Test without accumulating Q and Z
    //
    Rlaset("Full", n, n, zero, one, q, ldt);
    Rlaset("Full", n, n, zero, one, z, ldt);
    Rtgexc(false, false, n, t1, ldt, s1, ldt, q, ldt, z, ldt, ifst1, ilst1, work, lwork, info1);
    for (i = 1; i <= n; i = i + 1) {
        for (j = 1; j <= n; j = j + 1) {
            if (i == j && q[(i - 1) + (j - 1) * ldq] != one) {
                res += one / eps;
            }
            if (i != j && q[(i - 1) + (j - 1) * ldq] != zero) {
                res += one / eps;
            }
            if (i == j && z[(i - 1) + (j - 1) * ldz] != one) {
                res += one / eps;
            }
            if (i != j && z[(i - 1) + (j - 1) * ldz] != zero) {
                res += one / eps;
            }
        }
    }
    //
    //     Test with accumulating Q
    //
    Rlaset("Full", n, n, zero, one, q, ldt);
    Rlaset("Full", n, n, zero, one, z, ldt);
    Rtgexc(true, true, n, t2, ldt, s2, ldt, q, ldt, z, ldt, ifst2, ilst2, work, lwork, info2);
    //
    //     Compare T1 with T2 and S1 with S2
    //
    for (i = 1; i <= n; i = i + 1) {
        for (j = 1; j <= n; j = j + 1) {
            if (t1[(i - 1) + (j - 1) * ldt1] != t2[(i - 1) + (j - 1) * ldt2]) {
                res += one / eps;
            }
            if (s1[(i - 1) + (j - 1) * lds1] != s2[(i - 1) + (j - 1) * lds2]) {
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
    //     Test orthogonality of Q and Z and backward error on T2 and S2
    //
    Rget51(1, n, t, ldt, t2, ldt, q, ldt, z, ldt, work, result[1 - 1]);
    Rget51(1, n, s, ldt, s2, ldt, q, ldt, z, ldt, work, result[2 - 1]);
    Rget51(3, n, t, ldt, t2, ldt, q, ldt, q, ldt, work, result[3 - 1]);
    Rget51(3, n, t, ldt, t2, ldt, z, ldt, z, ldt, work, result[4 - 1]);
    res += result[1 - 1] + result[2 - 1] + result[3 - 1] + result[4 - 1];
    //
    //     Read next matrix pair
    //
    goto statement_10;
    //
    //     End of Rget40
    //
}
