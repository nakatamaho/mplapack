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
#include <mplapack_lin.h>

#include <time.h>

#include <iostream>
#include <sstream>
#include <vector>

using namespace std;
void Cchkrfp(void) {
    common cmn;
    common_write write(cmn);
    static const char *format_9991 = "(' Relative machine ',a,' is taken to be',d16.6)";
    static const char *format_9993 = "(4x,a4,':  ',10i6,/,11x,10i6)";
    static const char *format_9995 = "(' !! Invalid input value: ',a4,'=',i6,'; must be <=',i6)";
    static const char *format_9996 = "(' !! Invalid input value: ',a4,'=',i6,'; must be >=',i6)";
    time_t s1;
    time_t s2;
    bool fatal = false;

    const INTEGER nin = 5;
    const INTEGER nout = 6;
    const INTEGER maxin = 12;
    const INTEGER maxrhs = 16;
    const INTEGER nmax = 50;
    const INTEGER ntypes = 9;
    INTEGER nn = 0;
    INTEGER i = 0;
    INTEGER nval[maxin];
    INTEGER nns = 0;
    INTEGER nsval[maxin];
    INTEGER nnt = 0;
    char buf[1024];
    char tsterr_str[1];
    bool tsterr = false;
    REAL eps = 0.0;
    std::string str;
    s1 = time(NULL);
    fatal = false;
    //
    //     Read a dummy line.
    //
    getline(cin, str);
    //
    //     Report LAPACK version tag (e.g. LAPACK-3.2.0)
    //
    INTEGER vers_major = 0;
    INTEGER vers_minor = 0;
    INTEGER vers_patch = 0;
    INTEGER mplapack_vers_major = 0;
    INTEGER mplapack_vers_minor = 0;
    INTEGER mplapack_vers_patch = 0;
    INTEGER lapack_vers_major = 0;
    INTEGER lapack_vers_minor = 0;
    INTEGER lapack_vers_patch = 0;
    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
    write(nout, "(' Tests of the Multiple precision version of LAPACK MPLAPACK VERSION ',i1,'.',i1,'.',i1,/, "
                "' Based on original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, 'The following parameter values will be used:')"),
        mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
    //
    //     Read the values of N
    //
    stringstream ss(str);
    ss.str("");
    getline(cin, str);
    ss.str(str);
    ss >> nn;
    if (nn < 1) {
        write(nout, format_9996), " NN ", nn, 1;
        nn = 0;
        fatal = true;
    } else if (nn > maxin) {
        write(nout, format_9995), " NN ", nn, maxin;
        nn = 0;
        fatal = true;
    }
    ss.str("");
    getline(cin, str);
    ss.str(str);
    for (i = 1; i <= nn; i = i + 1) {
        ss >> nval[i - 1];
    }
    for (i = 1; i <= nn; i = i + 1) {
        if (nval[i - 1] < 0) {
            write(nout, format_9996), " N  ", nval[i - 1], 0;
            fatal = true;
        } else if (nval[i - 1] > nmax) {
            write(nout, format_9995), " N  ", nval[i - 1], nmax;
            fatal = true;
        }
    }
    if (nn > 0) {
        printf("     %4s  :", "N");
        for (i = 1; i <= nn; i = i + 1) {
            printf("%6ld", nval[i - 1]);
        }
        printf("\n");
    }
    //
    //     Read the values of NRHS
    //
    ss.str("");
    getline(cin, str);
    ss.str(str);
    ss >> nns;
    if (nns < 1) {
        write(nout, format_9996), " NNS", nns, 1;
        nns = 0;
        fatal = true;
    } else if (nns > maxin) {
        write(nout, format_9995), " NNS", nns, maxin;
        nns = 0;
        fatal = true;
    }
    ss.str("");
    getline(cin, str);
    ss.str(str);
    for (i = 1; i <= nns; i = i + 1)
        ss >> nsval[i - 1];
    for (i = 1; i <= nns; i = i + 1) {
        if (nsval[i - 1] < 0) {
            write(nout, format_9996), "NRHS", nsval[i - 1], 0;
            fatal = true;
        } else if (nsval[i - 1] > maxrhs) {
            write(nout, format_9995), "NRHS", nsval[i - 1], maxrhs;
            fatal = true;
        }
    }
    if (nns > 0) {
        printf("     %4s  :", "NRHS");
        for (i = 1; i <= nns; i = i + 1) {
            printf("%6ld", nsval[i - 1]);
        }
    }
    printf("\n");
    //
    //     Read the matrix types
    //
    getline(cin, str);
    ss >> nnt;
    if (nnt < 1) {
        write(nout, format_9996), " NMA", nnt, 1;
        nnt = 0;
        fatal = true;
    } else if (nnt > ntypes) {
        write(nout, format_9995), " NMA", nnt, ntypes;
        nnt = 0;
        fatal = true;
    }
    ss.str("");
    getline(cin, str);
    ss.str(str);
    INTEGER ntval[ntypes];
    for (i = 1; i <= nnt; i = i + 1) {
        ss >> ntval[i - 1];
    }
    for (i = 1; i <= nnt; i = i + 1) {
        if (ntval[i - 1] < 0) {
            write(nout, format_9996), " M  ", ntval[i - 1], 0;
            fatal = true;
        } else if (ntval[i - 1] > ntypes) {
            write(nout, format_9995), "TYPE", ntval[i - 1], ntypes;
            fatal = true;
        }
    }
    if (nnt > 0) {
        printf("     %4s  :", "M");
        for (i = 1; i <= nnt; i = i + 1) {
            printf("%6ld", ntval[i - 1]);
        }
    }
    //
    //     Read the threshold value for the test ratios.
    //
    REAL thresh = 0.0;
    ss.str("");
    getline(cin, str);
    ss.str(str);
    double dtmp;
    ss >> dtmp;
    thresh = dtmp;
    write(nout, "(/,' Routines pass computational tests if test ratio is ','less than',"
                "f8.2,/)"),
        cast2double(thresh);
    //
    //     Read the flag that indicates whether to test the error exits.
    //
    ss.str("");
    getline(cin, str);
    ss.str(str);
    ss >> tsterr;
    if (Mlsame(tsterr_str, "T"))
        tsterr = true;
    else
        tsterr = false;
    //
    if (fatal) {
        write(nout, "(/,' Execution not attempted due to input errors')");
        exit(0);
    }
    //
    //     Calculate and print the machine dependent constants.
    //
    eps = Rlamch("Underflow threshold");
    sprintnum_short(buf, eps);
    cout << " Relative machine underflow is taken to be : " << buf << endl;
    eps = Rlamch("Overflow threshold");
    cout << " Relative machine overflow  is taken to be : " << buf << endl;
    eps = Rlamch("Epsilon");
    cout << " Relative machine precision is taken to be : " << buf << endl;
    //
    //     Test the error exit of:
    //
    if (tsterr) {
        Cerrrfp(nout);
    }
    //
    //    Test the routines: Cpftrf, Cpftri, Cpftrs (as in Cdrvpo).
    //    This also tests the routines: Ctfsm, Ctftri, Ctfttr, Ctrttf.
    //
    COMPLEX worka[nmax * nmax];
    COMPLEX workasav[nmax * nmax];
    COMPLEX workafac[nmax * nmax];
    COMPLEX workainv[nmax * nmax];
    COMPLEX workb[nmax * maxrhs];
    COMPLEX workbsav[nmax * maxrhs];
    COMPLEX workxact[nmax * maxrhs];
    COMPLEX workx[nmax * maxrhs];
    COMPLEX workarf[(nmax * (nmax + 1)) / 2];
    COMPLEX workarfinv[(nmax * (nmax + 1)) / 2];
    COMPLEX z_work_Clatms[3 * nmax];
    COMPLEX z_work_Cpot02[nmax * maxrhs];
    COMPLEX z_work_Cpot03[nmax * nmax];
    REAL d_work_Clatms[nmax];
    REAL d_work_Clanhe[nmax];
    REAL d_work_Cpot01[nmax];
    REAL d_work_Cpot02[nmax];
    REAL d_work_Cpot03[nmax];
    Cdrvrfp(nout, nn, nval, nns, nsval, nnt, ntval, thresh, worka, workasav, workafac, workainv, workb, workbsav, workxact, workx, workarf, workarfinv, z_work_Clatms, z_work_Cpot02, z_work_Cpot03, d_work_Clatms, d_work_Clanhe, d_work_Cpot01, d_work_Cpot02, d_work_Cpot03);
    //
    //    Test the routine: Clanhf
    //
    Cdrvrf1(nout, nn, nval, thresh, worka, nmax, workarf, d_work_Clanhe);
    //
    //    Test the conversion routines:
    //       zhfttp, ztpthf, Ctfttr, Ctrttf, Ctrttp and Ctpttr.
    //
    COMPLEX workap[(nmax * (nmax + 1)) / 2];
    Cdrvrf2(nout, nn, nval, worka, nmax, workarf, workap, workasav);
    //
    //    Test the routine: Ctfsm
    //
    Cdrvrf3(nout, nn, nval, thresh, worka, nmax, workarf, workainv, workafac, d_work_Clanhe, z_work_Cpot03, z_work_Cpot02);
    //
    //    Test the routine: Chfrk
    //
    Cdrvrf4(nout, nn, nval, thresh, worka, workafac, nmax, workarf, workainv, nmax, d_work_Clanhe);
    //
    s2 = time(NULL);
    write(nout, "(/,' End of tests')");
    write(nout, "(' Total time used = ',f12.2,' seconds',/)"), int(s2 - s1);
    //
    //     End of Cchkrfp
    //
}

int main(int argc, char const *argv[]) { Cchkrfp(); }
