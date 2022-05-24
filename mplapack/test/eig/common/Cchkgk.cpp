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

inline REAL cabs1(COMPLEX cdum) { return abs(cdum.real()) + abs(cdum.imag()); }

void Cchkgk(INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    double dtmp;
    std::complex<double> ctmp;
    char buf[1024];
    COMPLEX cdum = 0.0;
    INTEGER lmax[4];
    INTEGER ninfo = 0;
    INTEGER knt = 0;
    const REAL zero = 0.0;
    REAL rmax = 0.0;
    REAL eps = 0.0;
    INTEGER n = 0;
    INTEGER m = 0;
    INTEGER i = 0;
    INTEGER j = 0;
    const INTEGER lda = 50;
    const INTEGER ldb = 50;
    const INTEGER ldvl = 50;
    const INTEGER ldvr = 50;
    const INTEGER lde = 50;
    const INTEGER ldf = 50;
    const INTEGER lrwork = 6 * 50;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const INTEGER ldwork = 50;
    INTEGER ldvlf = ldvl;
    INTEGER ldvrf = ldvr;
    COMPLEX a[lda * lda];
    COMPLEX b[ldb * ldb];
    COMPLEX vl[ldvl * ldvl];
    COMPLEX vr[ldvr * ldvr];
    REAL rwork[lrwork];
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    COMPLEX af[lda * lda];
    COMPLEX bf[ldb * ldb];
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL lscale[lda];
    REAL rscale[lda];
    INTEGER info = 0;
    COMPLEX vlf[ldvl * ldvl];
    COMPLEX vrf[ldvr * ldvr];
    COMPLEX work[ldwork * ldwork];
    COMPLEX e[lde * lde];
    COMPLEX f[ldf * ldf];
    REAL vmax = 0.0;
    //
    lmax[1 - 1] = 0;
    lmax[2 - 1] = 0;
    lmax[3 - 1] = 0;
    lmax[4 - 1] = 0;
    ninfo = 0;
    knt = 0;
    rmax = zero;
    //
    eps = Rlamch("Precision");
    string str;
    istringstream iss;
    double dtmp_r;
    double dtmp_i;
    char line[1024];
    //
    while (getline(cin, str)) {
        stringstream ss(str);
        ss >> n;
        ss >> m;
        if (n == 0)
            break;
        //
        string _r = regex_replace(str, regex("D\\+"), "e+");
        str = regex_replace(_r, regex("D\\-"), "e-");
        iss.clear();
        iss.str(str);
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
                a[(i - 1) + (j - 1) * lda] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
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
                b[(i - 1) + (j - 1) * ldb] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
        getline(cin, str); // ignore blank line
        //
        for (i = 1; i <= n; i = i + 1) {
            getline(cin, str);
            string ____r = regex_replace(str, regex(","), " ");
            string ___r = regex_replace(____r, regex("\\)"), " ");
            string __r = regex_replace(___r, regex("\\("), " ");
            string _r = regex_replace(__r, regex("D\\+"), "e+");
            str = regex_replace(_r, regex("D\\-"), "e-");
            iss.clear();
            iss.str(str);
            for (j = 1; j <= m; j = j + 1) {
                iss >> dtmp_r;
                iss >> dtmp_i;
                vl[(i - 1) + (j - 1) * ldvl] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
        //
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
            for (j = 1; j <= m; j = j + 1) {
                iss >> dtmp_r;
                iss >> dtmp_i;
                vr[(i - 1) + (j - 1) * ldvr] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
        getline(cin, str); // ignore blank line
        //
        knt++;
        //
        anorm = Clange("M", n, n, a, lda, rwork);
        bnorm = Clange("M", n, n, b, ldb, rwork);
        //
        Clacpy("FULL", n, n, a, lda, af, lda);
        Clacpy("FULL", n, n, b, ldb, bf, ldb);
        //
        Cggbal("B", n, a, lda, b, ldb, ilo, ihi, lscale, rscale, rwork, info);
        if (info != 0) {
            ninfo++;
            lmax[1 - 1] = knt;
        }
        //
        Clacpy("FULL", n, m, vl, ldvl, vlf, ldvl);
        Clacpy("FULL", n, m, vr, ldvr, vrf, ldvr);
        //
        Cggbak("B", "L", n, ilo, ihi, lscale, rscale, m, vl, ldvl, info);
        if (info != 0) {
            ninfo++;
            lmax[2 - 1] = knt;
        }
        //
        Cggbak("B", "R", n, ilo, ihi, lscale, rscale, m, vr, ldvr, info);
        if (info != 0) {
            ninfo++;
            lmax[3 - 1] = knt;
        }
        //
        //     Test of Cggbak
        //
        //     Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
        //     where tilde(A) denotes the transformed matrix.
        //
        Cgemm("N", "N", n, m, n, cone, af, lda, vr, ldvr, czero, work, ldwork);
        Cgemm("C", "N", m, m, n, cone, vl, ldvl, work, ldwork, czero, e, lde);
        //
        Cgemm("N", "N", n, m, n, cone, a, lda, vrf, ldvr, czero, work, ldwork);
        Cgemm("C", "N", m, m, n, cone, vlf, ldvl, work, ldwork, czero, f, ldf);
        //
        vmax = zero;
        for (j = 1; j <= m; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                vmax = max(vmax, cabs1(e[(i - 1) + (j - 1) * lde] - f[(i - 1) + (j - 1) * ldf]));
            }
        }
        vmax = vmax / (eps * max(anorm, bnorm));
        if (vmax > rmax) {
            lmax[4 - 1] = knt;
            rmax = vmax;
        }
        //
        //     Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR
        //
        Cgemm("N", "N", n, m, n, cone, bf, ldb, vr, ldvr, czero, work, ldwork);
        Cgemm("C", "N", m, m, n, cone, vl, ldvl, work, ldwork, czero, e, lde);
        //
        Cgemm("n", "n", n, m, n, cone, b, ldb, vrf, ldvr, czero, work, ldwork);
        Cgemm("C", "N", m, m, n, cone, vlf, ldvl, work, ldwork, czero, f, ldf);
        //
        vmax = zero;
        for (j = 1; j <= m; j = j + 1) {
            for (i = 1; i <= m; i = i + 1) {
                vmax = max(vmax, cabs1(e[(i - 1) + (j - 1) * lde] - f[(i - 1) + (j - 1) * ldf]));
            }
        }
        vmax = vmax / (eps * max(anorm, bnorm));
        if (vmax > rmax) {
            lmax[4 - 1] = knt;
            rmax = vmax;
        }
        //
    }
    //
    write(nout, "(1x,'.. test output of Cggbak .. ')");
    //
    sprintnum_short(buf, rmax);
    write(nout, "(' value of largest test error                  =',a)"), buf;
    write(nout, "(' example number where Cggbal info is not 0    =',i4)"), lmax[1 - 1];
    write(nout, "(' example number where Cggbak(L) info is not 0 =',i4)"), lmax[2 - 1];
    write(nout, "(' example number where Cggbak(R) info is not 0 =',i4)"), lmax[3 - 1];
    write(nout, "(' example number having largest error          =',i4)"), lmax[4 - 1];
    write(nout, "(' number of examples where info is not 0       =',i4)"), ninfo;
    write(nout, "(' total number of examples tested              =',i4)"), knt;
    //
    //     End of Cchkgk
    //
}
