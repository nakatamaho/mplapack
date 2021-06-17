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

void Rchkbk(INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    char buf[1024];
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
    REAL e[lde * lde];
    INTEGER j = 0;
    REAL ein[lde * lde];
    INTEGER ldein = lde;
    INTEGER info = 0;
    REAL vmax = 0.0;
    REAL x = 0.0;
    //
    lmax[1 - 1] = 0;
    lmax[2 - 1] = 0;
    ninfo = 0;
    knt = 0;
    rmax = zero;
    safmin = Rlamch("S");
    // following should be double of Rlamch("E") since input data is at most in double prec.
    eps = 1.1102230246251565E-016; // Rlamch("E");
    string str;
    char line[1024];
    double dtmp;
    //
    while (getline(cin, str)) {
        stringstream ss(str);
        ss >> n;
        ss >> ilo;
        ss >> ihi;
        if (n == 0)
            break;
        //
        // printf("%d %d %d\n", (int)n, (int)ilo, (int)ihi);
        getline(cin, str);
        string _r = regex_replace(str, regex("D\\+"), "e+");
        str = regex_replace(_r, regex("D\\-"), "e-");
        istringstream iss(str);
        for (i = 1; i <= n; i = i + 1) {
            iss >> dtmp;
            scale[i - 1] = dtmp;
        }
        // printf("scale=");printvec(scale,n);printf("\n");
        getline(cin, str); // ignore blank line
        for (i = 1; i <= n; i = i + 1) {
            getline(cin, str);
            string _r = regex_replace(str, regex("D\\+"), "e+");
            str = regex_replace(_r, regex("D\\-"), "e-");
            istringstream iss(str);
            for (j = 1; j <= n; j = j + 1) {
                iss >> dtmp;
                e[(i - 1) + (j - 1) * lde] = dtmp;
            }
        }
        //
        // printf("e=");printmat(n,n,e,lde);printf("\n");
        getline(cin, str); // ignore blank line
        for (i = 1; i <= n; i = i + 1) {
            getline(cin, str);
            string _r = regex_replace(str, regex("D\\+"), "e+");
            str = regex_replace(_r, regex("D\\-"), "e-");
            istringstream iss(str);
            for (j = 1; j <= n; j = j + 1) {
                iss >> dtmp;
                ein[(i - 1) + (j - 1) * ldein] = dtmp;
            }
        }
        //
        knt++;
        // printf("ein=");printmat(n,n,ein,lde);printf("\n");
        Rgebak("B", "R", n, ilo, ihi, scale, n, e, lde, info);
        // printf("eout=");printmat(n,n,e,lde);printf("\n");
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
                x = abs(e[(i - 1) + (j - 1) * lde] - ein[(i - 1) + (j - 1) * ldein]) / eps;
                if (abs(e[(i - 1) + (j - 1) * lde]) > safmin) {
                    x = x / abs(e[(i - 1) + (j - 1) * lde]);
                }
                vmax = max(vmax, x);
            }
        }
        //
        if (vmax > rmax) {
            lmax[2 - 1] = knt;
            rmax = vmax;
        }
        //
    }
    //
    //
    write(nout, "(1x,'.. test output of Rgebak .. ')");
    //
    sprintnum_short(buf, rmax);
    write(nout, "(1x,'value of largest test error             = ',a)"), buf;
    write(nout, "(1x,'example number where info is not zero   = ',i4)"), lmax[1 - 1];
    write(nout, "(1x,'example number having largest error     = ',i4)"), lmax[2 - 1];
    write(nout, "(1x,'number of examples where info is not 0  = ',i4)"), ninfo;
    write(nout, "(1x,'total number of examples tested         = ',i4)"), knt;
    //
    //     End of Rchkbk
    //
}
