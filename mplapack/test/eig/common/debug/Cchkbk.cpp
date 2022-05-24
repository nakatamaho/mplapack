/*
 * Copyright (c) 2021-2022
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

inline REAL abs1(COMPLEX cdum) { return abs(cdum.real()) + abs(cdum.imag()); }

void Cchkbk(INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    char buf[1024];
    COMPLEX cdum = 0.0;
    INTEGER lmax[2];
    INTEGER ninfo = 0;
    INTEGER knt = 0;
    const REAL zero = 0.0;
    REAL rmax = 0.0;
    REAL eps = 0.0;
    REAL safmin = 0.0;
    INTEGER n = 0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    const INTEGER lde = 20;
    REAL scale[lde];
    INTEGER i = 0;
    COMPLEX e[lde * lde];
    INTEGER j = 0;
    COMPLEX ein[lde * lde];
    INTEGER ldein = lde;
    INTEGER info = 0;
    REAL vmax = 0.0;
    REAL x = 0.0;
    double dtmp;
    //
    lmax[1 - 1] = 0;
    lmax[2 - 1] = 0;
    ninfo = 0;
    knt = 0;
    rmax = zero;
    eps = Rlamch("E");
    safmin = Rlamch("S");
    string str;
    istringstream iss;
    double dtmp_r;
    double dtmp_i;
    //
    while (getline(cin, str)) {
        stringstream ss(str);
        ss >> n;
        ss >> ilo;
        ss >> ihi;
        // printf("n is %d, iss is %d, ihi is %d\n", (int)n, (int)ilo, (int)ihi);
        if (n == 0)
            break;
        //
        getline(cin, str);
        string _r = regex_replace(str, regex("D\\+"), "e+");
        str = regex_replace(_r, regex("D\\-"), "e-");
        iss.clear();
        iss.str(str);
        for (i = 1; i <= n; i = i + 1) {
            iss >> dtmp;
            scale[i - 1] = dtmp;
        }
        printf("scale=");printvec(scale,n);printf("\n");
        getline(cin, str); // ignore blank line
        for (i = 1; i <= n; i = i + 1) {
            getline(cin, str);
            string ____r = regex_replace(str, regex(","), " ");
            string ___r = regex_replace(____r, regex("\\)"), " ");
            string __r = regex_replace(___r, regex("\\("), " ");
            string _r = regex_replace(__r, regex("D\\+"), "e+");
            str = regex_replace(_r, regex("D\\-"), "e-");
            iss.clear();
            iss.str(str);
            for (j = 1; j <= n; j = j + 1) {
                iss >> dtmp_r;
                iss >> dtmp_i;
                e[(i - 1) + (j - 1) * lde] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
        printf("e=");printmat(n,n,e,lde);printf("\n");
        getline(cin, str); // ignore blank line
        for (i = 1; i <= n; i = i + 1) {
            getline(cin, str);
            string ____r = regex_replace(str, regex(","), " ");
            string ___r = regex_replace(____r, regex("\\)"), " ");
            string __r = regex_replace(___r, regex("\\("), " ");
            string _r = regex_replace(__r, regex("D\\+"), "e+");
            str = regex_replace(_r, regex("D\\-"), "e-");
            iss.clear();
            iss.str(str);
            for (j = 1; j <= n; j = j + 1) {
                iss >> dtmp_r;
                iss >> dtmp_i;
                ein[(i - 1) + (j - 1) * lde] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
        printf("ein=");printmat(n,n,ein,lde);printf("\n");
        //
        knt++;
        Cgebak("B", "R", n, ilo, ihi, scale, n, e, lde, info);

        printf("eout=");printmat(n,n,e,lde);printf("\n");
        printf("info=%d\n",(int)info);
        // printf("ein-eout\n");
        //
        if (info != 0) {
            ninfo++;
            lmax[1 - 1] = knt;
        }
        //
        getline(cin, str); // ignore blank line
        vmax = zero;
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                x = abs1(e[(i - 1) + (j - 1) * lde] - ein[(i - 1) + (j - 1) * ldein]) / eps;
                if (abs1(e[(i - 1) + (j - 1) * lde]) > safmin) {
                    x = x / abs1(e[(i - 1) + (j - 1) * lde]);
                }
                vmax = max(vmax, x);
            }
        }
	printf("vmax=");printnum(vmax);printf("\n");
        //
        if (vmax > rmax) {
            lmax[2 - 1] = knt;
            rmax = vmax;
        }
        //
	printf("\n");
    }
//
statement_60:
    //
    write(nout, "(1x,'.. test output of Cgebak .. ')");
    //
    sprintnum_short(buf, rmax);
    write(nout, "(1x,'value of largest test error             = ',a)"), buf;
    write(nout, "(1x,'example number where info is not zero   = ',i4)"), lmax[1 - 1];
    write(nout, "(1x,'example number having largest error     = ',i4)"), lmax[2 - 1];
    write(nout, "(1x,'number of examples where info is not 0  = ',i4)"), ninfo;
    write(nout, "(1x,'total number of examples tested         = ',i4)"), knt;
    //
    //     End of Cchkbk
    //
}
