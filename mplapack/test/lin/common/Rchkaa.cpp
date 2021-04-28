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

void Rchkaa(void) {
    static REAL threq = 2.0;
    time_t s1;
    const INTEGER nmax = 132;
    bool fatal = false;
    const INTEGER nin = 5;
    int vers_major = 0;
    int vers_minor = 0;
    int vers_patch = 0;
    const INTEGER nout = 6;
    INTEGER nm = 0;
    const INTEGER maxin = 12;
    INTEGER mval[maxin];
    INTEGER i = 0;
    INTEGER nn = 0;
    INTEGER nval[maxin];
    INTEGER nns = 0;
    INTEGER nsval[maxin];
    const INTEGER maxrhs = 16;
    INTEGER nnb = 0;
    INTEGER nbval[maxin];
    INTEGER nnb2 = 0;
    INTEGER nb = 0;
    INTEGER j = 0;
    INTEGER nbval2[maxin];
    INTEGER nxval[maxin];
    INTEGER nrank = 0;
    INTEGER rankval[maxin];
    REAL thresh = 0.0;
    bool tstchk = false;
    bool tstdrv = false;
    bool tsterr = false;
    char tstchk_str[1];
    char tstdrv_str[1];
    char tsterr_str[1];
    REAL eps = 0.0;
    char aline[72];
    char path[3];
    const INTEGER matmax = 30;
    INTEGER nmats = 0;
    char c1[1];
    INTEGER k = 0;
    INTEGER ic = 0;
    char c2[2];
    INTEGER nrhs = 0;
    INTEGER ntypes = 0;
    bool dotype[matmax];
    INTEGER iwork[25 * nmax];
    const INTEGER kdmax = nmax + (nmax + 1) / 4;
    INTEGER la = 0;
    INTEGER lafac = 0;
    INTEGER piv[nmax];
    REAL a[((kdmax + 1) * nmax) * 7];
    INTEGER lda = ((kdmax + 1) * nmax);
    REAL b[(nmax * maxrhs) * 4];
    INTEGER ldb = (nmax * maxrhs);
    REAL e[nmax];
    REAL rwork[5 * nmax + 2 * maxrhs];
    REAL work[(nmax) * (3 * nmax + maxrhs + 30)];
    INTEGER lwork = 3 * nmax + maxrhs + 30;
    REAL s[2 * nmax];

    time_t s2;
    std::string str;

    s1 = time(NULL);
    lda = nmax;
    fatal = false;
    //
    //     Read a dummy line.
    //
    getline(cin, str);
    //
    //     Report values of parameters.
    //
    iMlaver(vers_major, vers_minor, vers_patch);
    printf("Tests of the DOUBLE PRECISION MPLAPACK routines \n");
    printf("MPLAPACK VERSION  %d.%d.%d\n", vers_major, vers_minor, vers_patch);
    printf("The following parameter values will be used:\n");
    //
    //     Read the values of M
    //
    getline(cin, str);
    stringstream ss(str);
    ss >> nm;
    if (nm < 1) {
        printf(" Invalid input value: %4s = %6ld ; must be >= %6d", "NM", nm, 1);
        nm = 0;
        fatal = true;
    } else if (nm > maxin) {
        printf(" Invalid input value: %4s = %6ld ; must be <= %6ld", "NM", nm, maxin);
        nm = 0;
        fatal = true;
    }
    ss.str("");
    getline(cin, str);
    ss.str(str);
    for (i = 1; i <= nm; i = i + 1) {
        ss >> mval[i - 1];
    }
    for (i = 1; i <= nm; i = i + 1) {
        if (mval[i - 1] < 0) {
            printf(" Invalid input value: %4s = %6ld; must be >= %6d\n", "M", mval[i - 1], 0);
            fatal = true;
        } else if (mval[i - 1] > nmax) {
            printf(" Invalid input value: %4s = %6ld; must be <= %6ld\n", "M", mval[i - 1], nmax);
            fatal = true;
        }
    }
    if (nm > 0) {
        printf("     %4s  :", "M");
        for (i = 1; i <= nm; i = i + 1) {
            printf("%6ld", mval[i - 1]);
        }
        printf("\n");
    }

    //
    //     Read the values of N
    //
    ss.str("");
    getline(cin, str);
    ss.str(str);
    ss >> nn;
    if (nn < 1) {
        printf(" Invalid input value: %4s = %6ld; must be >= %6d\n", " NN ", nn, 1);
        nn = 0;
        fatal = true;
    } else if (nn > maxin) {
        printf(" Invalid input value: %4s = %6ld; must be <= %6ld\n", " NN ", nn, maxin);
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
            printf(" Invalid input value: %4s = %6ld; must be >= %6d\n", " N ", nval[i - 1], 0);
            fatal = true;
        } else if (nval[i - 1] > nmax) {
            printf(" Invalid input value: %4s = %6ld; must be <= %6ld\n", " N ", nval[i - 1], nmax);
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
        printf(" Invalid input value: %4s = %6ld; must be >= %6d\n", "NNS", nns, 1);
        nns = 0;
        fatal = true;
    } else if (nns > maxin) {
        printf(" Invalid input value: %4s = %6ld; must be <= %6ld\n", "NNS", nns, maxin);
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
            printf(" Invalid input value: %4s = %6ld; must be >= %6d\n", "NRHS", nsval[i - 1], 0);
            fatal = true;
        } else if (nsval[i - 1] > maxrhs) {
            printf(" Invalid input value: %4s = %6ld; must be <= %6ld\n", "NRHS", nsval[i - 1], maxrhs);
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
    //     Read the values of NB
    //
    ss.str("");
    getline(cin, str);
    ss.str(str);
    ss >> nnb;
    if (nnb < 1) {
        printf(" Invalid input value: %4s = %6ld; must be >= %6d\n", "NNB ", nnb, 1);
        nnb = 0;
        fatal = true;
    } else if (nnb > maxin) {
        printf(" Invalid input value: %4s = %6ld; must be <= %6ld\n", "NNB ", nnb, maxin);
        nnb = 0;
        fatal = true;
    }
    ss.str("");
    getline(cin, str);
    ss.str(str);
    for (i = 1; i <= nnb; i = i + 1) {
        ss >> nbval[i - 1];
    }

    for (i = 1; i <= nnb; i = i + 1) {
        if (nbval[i - 1] < 0) {
            printf(" Invalid input value: %4s = %6ld; must be >= %6d\n", " NB ", nbval[i - 1], 0);
            fatal = true;
        }
    }
    if (nnb > 0) {
        printf("     %4s  :", "NB");
        for (i = 1; i <= nnb; i = i + 1) {
            printf("%6ld", nbval[i - 1]);
        }
    }
    printf("\n");
    //
    //     Set NBVAL2 to be the set of unique values of NB
    //
    nnb2 = 0;
    for (i = 1; i <= nnb; i = i + 1) {
        nb = nbval[i - 1];
        for (j = 1; j <= nnb2; j = j + 1) {
            if (nb == nbval2[j - 1]) {
                break;
            }
        }
        nnb2++;
        nbval2[nnb2 - 1] = nb;
    }
    //
    //     Read the values of NX
    //
    ss.str("");
    getline(cin, str);
    ss.str(str);
    for (i = 1; i <= nnb; i = i + 1) {
        ss >> nxval[i - 1];
    }
    for (i = 1; i <= nnb; i = i + 1) {
        if (nxval[i - 1] < 0) {
            printf(" Invalid input value: %4s = %6ld; must be >= %6d\n", "NX", nxval[i - 1], 0);
            fatal = true;
        }
    }
    if (nnb > 0) {
        printf("     %4s  :", "NX");
        for (i = 1; i <= nnb; i = i + 1)
            printf("%6ld", nxval[i - 1]);
    }
    printf("\n");
    //
    //     Read the values of RANKVAL
    //
    ss.str("");
    getline(cin, str);
    ss.str(str);
    ss >> nrank;
    if (nn < 1) {
        printf(" Invalid input value: %4s = %6ld ; must be >= %6d", "NRANK", nrank, 1);
        nrank = 0;
        fatal = true;
    } else if (nn > maxin) {
        printf(" Invalid input value: %4s = %6ld ; must be <= %6ld", "NRANK", nrank, maxin);
        nrank = 0;
        fatal = true;
    }
    ss.str("");
    getline(cin, str);
    ss.str(str);
    for (i = 1; i <= nrank; i = i + 1) {
        ss >> rankval[i - 1];
    }

    for (i = 1; i <= nrank; i = i + 1) {
        if (rankval[i - 1] < 0) {
            printf(" Invalid input value: %4s = %6ld ; must be >= %6d", "RANK", rankval[i - 1], 0);
            fatal = true;
        } else if (rankval[i - 1] > 100) {
            printf(" Invalid input value: %4s = %6ld ; must be <= %6d", "RANK", rankval[i - 1], 100);
            fatal = true;
        }
    }
    if (nrank > 0) {
        printf("     %4s  :", "RANK");
        for (i = 1; i <= nrank; i = i + 1)
            printf("%6ld", rankval[i - 1]);
    }
    printf("\n");
    //
    //     Read the threshold value for the test ratios.
    //
    ss.str("");
    getline(cin, str);
    ss.str(str);
    ss >> thresh;
    printf(" Routines pass computational tests if test ratio is less than %.2f \n", double(thresh));

    //
    //     Read the flag that indicates whether to test the LAPACK routines.
    //
    ss.str("");
    getline(cin, str);
    ss.str(str);
    ss >> tstchk_str;
    if (Mlsame(tstchk_str, "T"))
        tstchk = true;
    else
        tstchk = false;
    //
    //     Read the flag that indicates whether to test the driver routines.
    //
    ss.str("");
    getline(cin, str);
    ss.str(str);
    ss >> tstdrv_str;
    if (Mlsame(tstdrv_str, "T"))
        tstdrv = true;
    else
        tstdrv = false;
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
        printf(" Execution not attempted due to input errors \n");
        exit(0);
    }
    //
    //     Calculate and print the machine dependent constants.
    //
    eps = Rlamch("Underflow threshold");
    cout << " Relative machine underflow is taken to be : " << eps << endl;
    eps = Rlamch("Overflow threshold");
    cout << " Relative machine overflow  is taken to be : " << eps << endl;
    eps = Rlamch("Epsilon");
    cout << " Relative machine precision is taken to be : " << eps << endl;
    //
    //     Read a test path and the number of matrix types to use.
    //
    while (getline(cin, str)) {
        istringstream iss(str);
        vector<string> result;
        for (string s; iss >> s;)
            result.push_back(s);
        INTEGER n = result.size();
        if (n >= 1) {
            if (result[0].length() == 3) {
                path[0] = result[0][0];
                path[1] = result[0][1];
                path[2] = result[0][2];
                c1[0] = result[0][0];
                c2[0] = result[0][1];
                c2[1] = result[0][2];
            } else {
                printf("wrong three letters\n");
                exit(0);
            }
            if (n >= 2) {
                nmats = stoi(result[1]);
            }
        }
        nrhs = nsval[1 - 1];
        //
        //     Check first character for correct precision.
        //
        if (!Mlsame(c1, "Double precision")) {
            printf(" %3s :  Unrecognized path name\n", path);
        } else if (nmats <= 0) {
            //
            //        Check for a positive number of tests requested.
            //
            printf(" %3s   routines were not tested\n", path);
            //
        } else if (Mlsamen(2, c2, "GE")) {
            //
            //        GE:  general matrices
            //
            ntypes = 11;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkge(dotype, nm, mval, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                printf(" %3s   routines were not tested\n", path);
            }
#ifdef DODRV
            if (tstdrv) {
                Rdrvge(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], s, work, rwork, iwork, nout);
            } else {
                printf(" %3s   driver routines were not tested\n", path);
            }
#endif
        }
#ifdef NOTYET
        else if (Mlsamen2, c2, "GB") {
            //
            //        GB:  general banded matrices
            //
            la = (2 * kdmax + 1) * nmax;
            lafac = (3 * kdmax + 1) * nmax;
            ntypes = 8;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkgb(dotype, nm, mval, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, &a[(1 - 1)], la, &a[(3 - 1) * lda], lafac, &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvgb(dotype, nn, nval, nrhs, thresh, tsterr, &a[(1 - 1)], la, &a[(3 - 1) * lda], lafac, &a[(6 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], s, work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "GT") {
            //
            //        GT:  general tridiagonal matrices
            //
            ntypes = 12;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkgt(dotype, nn, nval, nns, nsval, thresh, tsterr, &a[(1 - 1)], &a[(2 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvgt(dotype, nn, nval, nrhs, thresh, tsterr, &a[(1 - 1)], &a[(2 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "PO") {
            //
            //        PO:  positive definite matrices
            //
            ntypes = 9;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkpo(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvpo(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], s, work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "PS") {
            //
            //        PS:  positive semi-definite matrices
            //
            ntypes = 9;
            //
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkps(dotype, nn, nval, nnb2, nbval2, nrank, rankval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], piv, work, rwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "PP") {
            //
            //        PP:  positive definite packed matrices
            //
            ntypes = 9;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkpp(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvpp(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], s, work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "PB") {
            //
            //        PB:  positive definite banded matrices
            //
            ntypes = 8;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkpb(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvpb(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], s, work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "PT") {
            //
            //        PT:  positive definite tridiagonal matrices
            //
            ntypes = 12;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkpt(dotype, nn, nval, nns, nsval, thresh, tsterr, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvpt(dotype, nn, nval, nrhs, thresh, tsterr, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "SY") {
            //
            //        SY:  symmetric indefinite matrices,
            //             with partial (Bunch-Kaufman) pivoting algorithm
            //
            ntypes = 10;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchksy(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvsy(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "SR") {
            //
            //        SR:  symmetric indefinite matrices,
            //             with bounded Bunch-Kaufman (rook) pivoting algorithm
            //
            ntypes = 10;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchksy_rook(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvsy_rook(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "SK") {
            //
            //        SK:  symmetric indefinite matrices,
            //             with bounded Bunch-Kaufman (rook) pivoting algorithm,
            //             different matrix storage format than SR path version.
            //
            ntypes = 10;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchksy_rk(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], e, &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvsy_rk(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], e, &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "SA") {
            //
            //        SA:  symmetric indefinite matrices,
            //             with partial (Aasen's) pivoting algorithm
            //
            ntypes = 10;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchksy_aa(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvsy_aa(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "S2") {
            //
            //        SA:  symmetric indefinite matrices,
            //             with partial (Aasen's) pivoting algorithm
            //
            ntypes = 10;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchksy_aa_2stage(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvsy_aa_2stage(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "SP") {
            //
            //        SP:  symmetric indefinite packed matrices,
            //             with partial (Bunch-Kaufman) pivoting algorithm
            //
            ntypes = 10;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchksp(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
            if (tstdrv) {
                Rdrvsp(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "TR") {
            //
            //        TR:  triangular matrices
            //
            ntypes = 18;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchktr(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "TP") {
            //
            //        TP:  triangular packed matrices
            //
            ntypes = 18;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchktp(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "TB") {
            //
            //        TB:  triangular banded matrices
            //
            ntypes = 17;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchktb(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "QR") {
            //
            //        QR:  QR factorization
            //
            ntypes = 8;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkqr(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "LQ") {
            //
            //        LQ:  LQ factorization
            //
            ntypes = 8;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchklq(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, rwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "QL") {
            //
            //        QL:  QL factorization
            //
            ntypes = 8;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkql(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, rwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "RQ") {
            //
            //        RQ:  RQ factorization
            //
            ntypes = 8;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkrq(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, rwork, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "QP") {
            //
            //        QP:  QR factorization with pivoting
            //
            ntypes = 6;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchkq3(dotype, nm, mval, nn, nval, nnb, nbval, nxval, thresh, &a[(1 - 1)], &a[(2 - 1) * lda], &b[(1 - 1)], &b[(3 - 1) * ldb], work, iwork, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "TZ") {
            //
            //        TZ:  Trapezoidal matrix
            //
            ntypes = 3;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstchk) {
                Rchktz(dotype, nm, mval, nn, nval, thresh, tsterr, &a[(1 - 1)], &a[(2 - 1) * lda], &b[(1 - 1)], &b[(3 - 1) * ldb], work, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "LS") {
            //
            //        LS:  Least squares drivers
            //
            ntypes = 6;
            Alareq(path, nmats, dotype, ntypes, nin, nout);
            //
            if (tstdrv) {
                Rdrvls(dotype, nm, mval, nn, nval, nns, nsval, nnb, nbval, nxval, thresh, tsterr, &a[(1 - 1)], &a[(2 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], rwork, &rwork[(nmax + 1) - 1], nout);
            } else {
                write(nout, format_9988), path;
            }
            //
        } else if (Mlsamen2, c2, "EQ") {
            //
            //        EQ:  Equilibration routines for general and positive definite
            //             matrices (THREQ should be between 2 and 10)
            //
            if (tstchk) {
                Rchkeq(threq, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "QT") {
            //
            //        QT:  QRT routines for general matrices
            //
            if (tstchk) {
                Rchkqrt(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "QX") {
            //
            //        QX:  QRT routines for triangular-pentagonal matrices
            //
            if (tstchk) {
                Rchkqrtp(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "TQ") {
            //
            //        TQ:  LQT routines for general matrices
            //
            if (tstchk) {
                Rchklqt(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "XQ") {
            //
            //        XQ:  LQT routines for triangular-pentagonal matrices
            //
            if (tstchk) {
                Rchklqtp(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "TS") {
            //
            //        TS:  QR routines for tall-skinny matrices
            //
            if (tstchk) {
                Rchktsqr(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else if (Mlsamen2, c2, "HH") {
            //
            //        HH:  Householder reconstruction for tall-skinny matrices
            //
            if (tstchk) {
                Rchkorhr_col(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
            } else {
                write(nout, format_9989), path;
            }
            //
        } else {
            //
            write(nout, format_9990), path;
        }
#endif
        //
        //     Go back to get another input line.
        //
    }
    //
    //     Branch to this line when the last record is read.
    //
    s2 = time(NULL);
    printf(" End of tests\n");
    printf(" Total time used =  %d seconds\n", int(s2 - s1));
    //
    //     End of Rchkaa
    //
}

int main(int argc, char const *argv[]) { Rchkaa(); }
