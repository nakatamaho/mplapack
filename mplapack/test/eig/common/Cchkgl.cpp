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

void Cchkgl(INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    char buf[1024];
    INTEGER lmax[3];
    INTEGER ninfo = 0;
    INTEGER knt = 0;
    const REAL zero = 0.0;
    REAL rmax = 0.0;
    REAL eps = 0.0;
    INTEGER n = 0;
    INTEGER i = 0;
    const INTEGER lda = 20;
    COMPLEX a[lda * lda];
    INTEGER j = 0;
    const INTEGER ldb = 20;
    COMPLEX b[ldb * ldb];
    INTEGER iloin = 0;
    INTEGER ihiin = 0;
    COMPLEX ain[lda * lda];
    COMPLEX bin[ldb * ldb];
    INTEGER ldain = lda;
    INTEGER ldbin = ldb;
    REAL lsclin[lda];
    REAL rsclin[lda];
    const INTEGER lwork = 6 * lda;
    REAL work[lwork];
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL lscale[lda];
    REAL rscale[lda];
    INTEGER info = 0;
    REAL vmax = 0.0;
    //
    lmax[1 - 1] = 0;
    lmax[2 - 1] = 0;
    lmax[3 - 1] = 0;
    ninfo = 0;
    knt = 0;
    rmax = zero;
    double dtmp;
    std::complex<double> ctmp;
    //
    eps = Rlamch("Precision");
    string str;
    istringstream iss;
    double dtmp_r;
    double dtmp_i;
    //
    while (getline(cin, str)) {
        getline(cin, str);
        stringstream ss(str);
        ss >> n;
        // printf("n is %d\n", (int)n);
        if (n == 0)
            break;
        for (i = 1; i <= n; i = i + 1) {
            getline(cin, str);
            string ____r = regex_replace(str, regex(","), " ");
            string ___r = regex_replace(____r, regex("\\)"), " ");
            string __r = regex_replace(___r, regex("\\("), " ");
            string _r = regex_replace(__r, regex("D\\+"), "e+");
            str = regex_replace(_r, regex("D\\-"), "e-");
            iss.clear();
            iss.str(str);
            // cout << str << "\n";
            for (j = 1; j <= n; j = j + 1) {
                iss >> dtmp_r;
                iss >> dtmp_i;
                a[(i - 1) + (j - 1) * lda] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
        // printf("a="); printmat(n, n, a, lda); printf("\n");
        getline(cin, str);
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
                b[(i - 1) + (j - 1) * lda] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
        // printf("b="); printmat(n, n, b, lda); printf("\n");
        //
        getline(cin, str);
        getline(cin, str);
        istringstream iss(str);
        iss >> iloin;
        iss >> ihiin;
        // printf("iloin is %d ihiin is %d\n", (int)iloin, (int)ihiin);  
        getline(cin, str);
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
                ain[(i - 1) + (j - 1) * lda] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
        getline(cin, str);
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
                bin[(i - 1) + (j - 1) * lda] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
        //
        getline(cin, str);
        getline(cin, str);
        string _r = regex_replace(str, regex("D\\+"), "e+");
        str = regex_replace(_r, regex("D\\-"), "e-");
        iss.clear();
        iss.str(str);
        for (i = 1; i <= n; i = i + 1) {
            iss >> dtmp;
            lsclin[i - 1] = dtmp;
        }
        // printf("lsclin=");printvec(lsclin,n);printf("\n");
        getline(cin, str);
        getline(cin, str);
        _r = regex_replace(str, regex("D\\+"), "e+");
        str = regex_replace(_r, regex("D\\-"), "e-");
        iss.clear();
        iss.str(str);
        for (i = 1; i <= n; i = i + 1) {
            iss >> dtmp;
            rsclin[i - 1] = dtmp;
        }
        // printf("rsclin=");printvec(rsclin,n);printf("\n");
        //
        anorm = Clange("M", n, n, a, lda, work);
        bnorm = Clange("M", n, n, b, ldb, work);
        //
        knt++;
        //
        Cggbal("B", n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info);
        //
        if (info != 0) {
            ninfo++;
            lmax[1 - 1] = knt;
        }
        //
        if (ilo != iloin || ihi != ihiin) {
            ninfo++;
            lmax[2 - 1] = knt;
        }
        //
        vmax = zero;
        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                vmax = max(vmax, REAL(abs(a[(i - 1) + (j - 1) * lda] - ain[(i - 1) + (j - 1) * ldain])));
                vmax = max(vmax, REAL(abs(b[(i - 1) + (j - 1) * ldb] - bin[(i - 1) + (j - 1) * ldbin])));
            }
        }
        //
        for (i = 1; i <= n; i = i + 1) {
            vmax = max(vmax, REAL(abs(lscale[i - 1] - lsclin[i - 1])));
            vmax = max(vmax, REAL(abs(rscale[i - 1] - rsclin[i - 1])));
        }
        //
        vmax = vmax / (eps * max(anorm, bnorm));
        //
        if (vmax > rmax) {
            lmax[3 - 1] = knt;
            rmax = vmax;
        }
        //
    }
    //
    write(nout, "(' .. test output of Cggbal .. ')");
    //
    sprintnum_short(buf, rmax);
    write(nout, "(1x,'value of largest test error            = ',a)"), buf;
    write(nout, "(' example number where info is not zero    = ',i4)"), lmax[1 - 1];
    write(nout, "(' example number where ILO or IHI is wrong = ',i4)"), lmax[2 - 1];
    write(nout, "(' example number having largest error      = ',i4)"), lmax[3 - 1];
    write(nout, "(' number of examples where info is not 0   = ',i4)"), ninfo;
    write(nout, "(' total number of examples tested          = ',i4)"), knt;
    //
    //     End of Cchkgl
    //
}
