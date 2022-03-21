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

#include <time.h>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

void Cchkee(void) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);

    const INTEGER nmax = 132;
    const INTEGER ncmax = 20;
    const INTEGER need = 14;
    const INTEGER lwork = nmax * (5 * nmax + 20);
    const INTEGER liwork = nmax * (nmax + 20);
    const INTEGER maxin = 20;
    const INTEGER maxt = 30;
    const INTEGER nout = 6;
    const INTEGER nin = 5;

    bool nep = false;
    bool sep = false;
    bool svd = false;
    bool zev = false;
    bool zes = false;
    bool zvx = false;
    bool zsx = false;
    bool zgg = false;
    bool zgs = false;
    bool zgx = false;
    bool zgv = false;
    bool zxv = false;
    bool zhb = false;
    bool zbb = false;
    bool glm = false;
    bool gqr = false;
    bool gsv = false;
    bool csd = false;
    bool lse = false;
    bool zbl = false;
    bool zbk = false;
    bool zgl = false;
    bool zgk = false;
    bool fatal = false;
    bool tstchk = false;
    bool tstdif = false;
    bool tstdrv = false;
    bool tsterr = false;

    char c1[2];
    char c3[4];
    char path[4];
    char vname[32];
    char line[1024];

    INTEGER i = 0;
    INTEGER i1 = 0;
    INTEGER ic = 0;
    INTEGER info = 0;
    INTEGER itmp = 0;
    INTEGER k = 0;
    INTEGER lenp = 0;
    INTEGER maxtyp = 0;
    INTEGER newsd = 0;
    INTEGER nk = 0;
    INTEGER nn = 0;
    INTEGER nparms = 0;
    INTEGER nrhs = 0;
    INTEGER ntypes = 0;

    INTEGER mplapack_vers_major = 0;
    INTEGER mplapack_vers_minor = 0;
    INTEGER mplapack_vers_patch = 0;
    INTEGER lapack_vers_major = 0;
    INTEGER lapack_vers_minor = 0;
    INTEGER lapack_vers_patch = 0;
    INTEGER n_threads = 0;

    REAL eps = 0.0;
    time_t s1 = 0.0;
    time_t s2 = 0.0;
    REAL thresh = 0.0;
    REAL thrshn = 0.0;

    bool dotype[maxt];
    bool logwrk[nmax];

    INTEGER ioldsd[4];
    INTEGER iseed[4];
    INTEGER iwork[liwork];
    INTEGER kval[maxin];
    INTEGER mval[maxin];
    INTEGER mxbval[maxin];
    INTEGER nbcol[maxin];
    INTEGER nbmin[maxin];
    INTEGER nbval[maxin];
    INTEGER nsval[maxin];
    INTEGER nval[maxin];
    INTEGER nxval[maxin];
    INTEGER pval[maxin];
    INTEGER inmin[maxin];
    INTEGER inwin[maxin];
    INTEGER inibl[maxin];
    INTEGER ishfts[maxin];
    INTEGER iacc22[maxin];

    REAL alpha[nmax];
    REAL beta[nmax];
    REAL dr[nmax * 12];
    REAL result[500];

    COMPLEX dc[nmax * 6];
    COMPLEX taua[nmax];
    COMPLEX taub[nmax];
    COMPLEX x[5 * nmax];

    COMPLEX *s = new COMPLEX[nmax * nmax];
    COMPLEX *a = new COMPLEX[nmax * nmax * need];
    COMPLEX *b = new COMPLEX[nmax * nmax * 5];
    COMPLEX *c = new COMPLEX[ncmax * ncmax * ncmax * ncmax];
    REAL *rwork = new REAL[lwork];
    COMPLEX *work = new COMPLEX[lwork];

    for (int ppp = 0; ppp < nmax * nmax; ppp++)
        s[ppp] = 0.0;
    for (int ppp = 0; ppp < nmax * nmax * need; ppp++)
        a[ppp] = 0.0;
    for (int ppp = 0; ppp < nmax * nmax * 5; ppp++)
        b[ppp] = 0.0;
    for (int ppp = 0; ppp < ncmax * ncmax * ncmax * ncmax; ppp++)
        c[ppp] = 0.0;
    for (int ppp = 0; ppp < lwork; ppp++)
        rwork[ppp] = 0.0;
    for (int ppp = 0; ppp < lwork; ppp++)
        work[ppp] = 0.0;
    static const char *format_9973 = "(/,1x,71('-'))";
    static const char *format_9980 = "(' *** Error code from ',a,' = ',i4)";
    static const char *format_9981 = "(' Relative machine ',a,' is taken to be',d16.6)";
    static const char *format_9983 = "(4x,a,10i6,/,10x,10i6)";
    static const char *format_9988 = "(' Invalid input value: ',a,'=',i6,'; must be <=',i6)";
    static const char *format_9989 = "(' Invalid input value: ',a,'=',i6,'; must be >=',i6)";
    static const char *format_9990 = "(/,/,1x,a3,' routines were not tested')";
    static const char *format_9992 = "(1x,a3,':  Unrecognized path name')";
    static const char *format_9997 = "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4)";

    std::string str;
    stringstream ss(str);
    s1 = time(NULL);
    fatal = false;
//
//     Return to here to read multiple sets of data
//
statement_10:
    //
    //     Read the first line and set the 3-character test path
    //
    getline(cin, str);
    char line[1024];
    ss >> line;
    path[0] = line[0];
    path[1] = line[1];
    path[2] = line[2];
    nep = Mlsamen(3, path, "NEP") || Mlsamen(3, path, "ZHS");
    sep = Mlsamen(3, path, "SEP") || Mlsamen(3, path, "ZST") || Mlsamen(3, path, "ZSG") || Mlsamen(3, path, "SE2");
    svd = Mlsamen(3, path, "SVD") || Mlsamen(3, path, "ZBD");
    zev = Mlsamen(3, path, "ZEV");
    zes = Mlsamen(3, path, "ZES");
    zvx = Mlsamen(3, path, "ZVX");
    zsx = Mlsamen(3, path, "ZSX");
    zgg = Mlsamen(3, path, "ZGG");
    zgs = Mlsamen(3, path, "ZGS");
    zgx = Mlsamen(3, path, "ZGX");
    zgv = Mlsamen(3, path, "ZGV");
    zxv = Mlsamen(3, path, "ZXV");
    zhb = Mlsamen(3, path, "ZHB");
    zbb = Mlsamen(3, path, "ZBB");
    glm = Mlsamen(3, path, "GLM");
    gqr = Mlsamen(3, path, "GQR") || Mlsamen(3, path, "GRQ");
    gsv = Mlsamen(3, path, "GSV");
    csd = Mlsamen(3, path, "CSD");
    lse = Mlsamen(3, path, "LSE");
    zbl = Mlsamen(3, path, "ZBL");
    zbk = Mlsamen(3, path, "ZBK");
    zgl = Mlsamen(3, path, "ZGL");
    zgk = Mlsamen(3, path, "ZGK");
    //
    //     Report values of parameters.
    //
    if (path == "   ") {
        goto statement_10;
    } else if (nep) {
        write(nout, "(' Tests of the Nonsymmetric Eigenvalue Problem routines')");
    } else if (sep) {
        write(nout, "(' Tests of the Hermitian Eigenvalue Problem routines')");
    } else if (svd) {
        write(nout, "(' Tests of the Singular Value Decomposition routines')");
    } else if (zev) {
        write(nout, "(/,' Tests of the Nonsymmetric Eigenvalue Problem Driver',/,"
                    "'    Cgeev (eigenvalues and eigevectors)')");
    } else if (zes) {
        write(nout, "(/,' Tests of the Nonsymmetric Eigenvalue Problem Driver',/,"
                    "'    Cgees (Schur form)')");
    } else if (zvx) {
        write(nout, "(/,' Tests of the Nonsymmetric Eigenvalue Problem Expert',' Driver',/,"
                    "'    Cgeevx (eigenvalues, eigenvectors and',' condition numbers)')");
    } else if (zsx) {
        write(nout, "(/,' Tests of the Nonsymmetric Eigenvalue Problem Expert',' Driver',/,"
                    "'    Cgeesx (Schur form and condition',' numbers)')");
    } else if (zgg) {
        write(nout, "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                    "'Problem routines')");
    } else if (zgs) {
        write(nout, "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                    "'Problem Driver Cgges')");
    } else if (zgx) {
        write(nout, "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                    "'Problem Expert Driver Cggesx')");
    } else if (zgv) {
        write(nout, "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                    "'Problem Driver Cggev')");
    } else if (zxv) {
        write(nout, "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                    "'Problem Expert Driver Cggevx')");
    } else if (zhb) {
        write(nout, "(' Tests of Chbtrd',/,' (reduction of a Hermitian band ',"
                    "'matrix to real tridiagonal form)')");
    } else if (zbb) {
        write(nout, "(' Tests of Cgbbrd',/,' (reduction of a general band ',"
                    "'matrix to real bidiagonal form)')");
    } else if (glm) {
        write(nout, "(/,' Tests of the Generalized Linear Regression Model ','routines')");
    } else if (gqr) {
        write(nout, "(/,' Tests of the Generalized QR and RQ routines')");
    } else if (gsv) {
        write(nout, "(/,' Tests of the Generalized Singular Value',"
                    "' Decomposition routines')");
    } else if (csd) {
        write(nout, "(/,' Tests of the CS Decomposition routines')");
    } else if (lse) {
        write(nout, "(/,' Tests of the Linear Least Squares routines')");
    } else if (zbl) {
        //
        //        Cgebal:  Balancing
        //
        Cchkbl(nin, nout);
        goto statement_380;
    } else if (zbk) {
        //
        //        Cgebak:  Back transformation
        //
        Cchkbk(nin, nout);
        goto statement_380;
    } else if (zgl) {
        //
        //        Cggbal:  Balancing
        //
        Cchkgl(nin, nout);
        goto statement_380;
    } else if (zgk) {
        //
        //        Cggbak:  Back transformation
        //
        Cchkgk(nin, nout);
        goto statement_380;
    } else if (Mlsamen(3, path, "ZEC")) {
        //
        //        ZEC:  Eigencondition estimation
        //
        read(nin, star), thresh;
        xlaenv(1, 1);
        xlaenv(12, 1);
        tsterr = true;
        Cchkec(thresh, tsterr, nin, nout);
        goto statement_380;
    } else {
        write(nout, format_9992), path;
        goto statement_380;
    }
    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
    write(nout, "(' Tests of the Multiple precision version of LAPACK MPLAPACK VERSION ',i1,'.',i1,'.',i1,/, "
                "' Based on original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, 'The following parameter values will be used:')"),
        mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
    //
    //     Read the number of values of M, P, and N.
    //
    read(nin, star), nn;
    if (nn < 0) {
        write(nout, format_9989), "   NN ", nn, 1;
        nn = 0;
        fatal = true;
    } else if (nn > maxin) {
        write(nout, format_9988), "   NN ", nn, maxin;
        nn = 0;
        fatal = true;
    }
    //
    //     Read the values of M
    //
    if (!(zgx || zxv)) {
        {
            read_loop rloop(cmn, nin, star);
            for (i = 1; i <= nn; i = i + 1) {
                rloop, mval[i - 1];
            }
        }
        if (svd) {
            strncmp(vname, "    M ", 8);
        } else {
            strncmp(vname, "    N ", 8);
        }
        for (i = 1; i <= nn; i = i + 1) {
            if (mval[i - 1] < 0) {
                write(nout, format_9989), vname, mval[i - 1], 0;
                fatal = true;
            } else if (mval[i - 1] > nmax) {
                write(nout, format_9988), vname, mval[i - 1], nmax;
                fatal = true;
            }
        }
        {
            write_loop wloop(cmn, nout, format_9983);
            wloop, "M:    ";
            for (i = 1; i <= nn; i = i + 1) {
                wloop, mval[i - 1];
            }
        }
    }
    //
    //     Read the values of P
    //
    if (glm || gqr || gsv || csd || lse) {
        {
            read_loop rloop(cmn, nin, star);
            for (i = 1; i <= nn; i = i + 1) {
                rloop, pval[i - 1];
            }
        }
        for (i = 1; i <= nn; i = i + 1) {
            if (pval[i - 1] < 0) {
                write(nout, format_9989), " P  ", pval[i - 1], 0;
                fatal = true;
            } else if (pval[i - 1] > nmax) {
                write(nout, format_9988), " P  ", pval[i - 1], nmax;
                fatal = true;
            }
        }
        {
            write_loop wloop(cmn, nout, format_9983);
            wloop, "P:    ";
            for (i = 1; i <= nn; i = i + 1) {
                wloop, pval[i - 1];
            }
        }
    }
    //
    //     Read the values of N
    //
    if (svd || zbb || glm || gqr || gsv || csd || lse) {
        {
            read_loop rloop(cmn, nin, star);
            for (i = 1; i <= nn; i = i + 1) {
                rloop, nval[i - 1];
            }
        }
        for (i = 1; i <= nn; i = i + 1) {
            if (nval[i - 1] < 0) {
                write(nout, format_9989), "    N ", nval[i - 1], 0;
                fatal = true;
            } else if (nval[i - 1] > nmax) {
                write(nout, format_9988), "    N ", nval[i - 1], nmax;
                fatal = true;
            }
        }
    } else {
        for (i = 1; i <= nn; i = i + 1) {
            nval[i - 1] = mval[i - 1];
        }
    }
    if (!(zgx || zxv)) {
        {
            write_loop wloop(cmn, nout, format_9983);
            wloop, "N:    ";
            for (i = 1; i <= nn; i = i + 1) {
                wloop, nval[i - 1];
            }
        }
    } else {
        write(nout, format_9983), "N:    ", nn;
    }
    //
    //     Read the number of values of K, followed by the values of K
    //
    if (zhb || zbb) {
        read(nin, star), nk;
        {
            read_loop rloop(cmn, nin, star);
            for (i = 1; i <= nk; i = i + 1) {
                rloop, kval[i - 1];
            }
        }
        for (i = 1; i <= nk; i = i + 1) {
            if (kval[i - 1] < 0) {
                write(nout, format_9989), "    K ", kval[i - 1], 0;
                fatal = true;
            } else if (kval[i - 1] > nmax) {
                write(nout, format_9988), "    K ", kval[i - 1], nmax;
                fatal = true;
            }
        }
        {
            write_loop wloop(cmn, nout, format_9983);
            wloop, "K:    ";
            for (i = 1; i <= nk; i = i + 1) {
                wloop, kval[i - 1];
            }
        }
    }
    //
    if (zev || zes || zvx || zsx) {
        //
        //        For the nonsymmetric QR driver routines, only one set of
        //        parameters is allowed.
        //
        read(nin, star), nbval[1 - 1], nbmin[1 - 1], nxval[1 - 1], inmin[1 - 1], inwin[1 - 1], inibl[1 - 1], ishfts[1 - 1], iacc22[1 - 1];
        if (nbval[1 - 1] < 1) {
            write(nout, format_9989), "   NB ", nbval[1 - 1], 1;
            fatal = true;
        } else if (nbmin[1 - 1] < 1) {
            write(nout, format_9989), "NBMIN ", nbmin[1 - 1], 1;
            fatal = true;
        } else if (nxval[1 - 1] < 1) {
            write(nout, format_9989), "   NX ", nxval[1 - 1], 1;
            fatal = true;
        } else if (inmin[1 - 1] < 1) {
            write(nout, format_9989), "   INMIN ", inmin[1 - 1], 1;
            fatal = true;
        } else if (inwin[1 - 1] < 1) {
            write(nout, format_9989), "   INWIN ", inwin[1 - 1], 1;
            fatal = true;
        } else if (inibl[1 - 1] < 1) {
            write(nout, format_9989), "   INIBL ", inibl[1 - 1], 1;
            fatal = true;
        } else if (ishfts[1 - 1] < 1) {
            write(nout, format_9989), "   ISHFTS ", ishfts[1 - 1], 1;
            fatal = true;
        } else if (iacc22[1 - 1] < 0) {
            write(nout, format_9989), "   IACC22 ", iacc22[1 - 1], 0;
            fatal = true;
        }
        xlaenv(1, nbval[1 - 1]);
        xlaenv(2, nbmin[1 - 1]);
        xlaenv(3, nxval[1 - 1]);
        xlaenv(12, max((INTEGER)11, inmin[1 - 1]));
        xlaenv(13, inwin[1 - 1]);
        xlaenv(14, inibl[1 - 1]);
        xlaenv(15, ishfts[1 - 1]);
        xlaenv(16, iacc22[1 - 1]);
        write(nout, format_9983), "NB:   ", nbval[1 - 1];
        write(nout, format_9983), "NBMIN:", nbmin[1 - 1];
        write(nout, format_9983), "NX:   ", nxval[1 - 1];
        write(nout, format_9983), "INMIN:   ", inmin[1 - 1];
        write(nout, format_9983), "INWIN: ", inwin[1 - 1];
        write(nout, format_9983), "INIBL: ", inibl[1 - 1];
        write(nout, format_9983), "ISHFTS: ", ishfts[1 - 1];
        write(nout, format_9983), "IACC22: ", iacc22[1 - 1];
        //
    } else if (zgs || zgx || zgv || zxv) {
        //
        //        For the nonsymmetric generalized driver routines, only one set of
        //        parameters is allowed.
        //
        read(nin, star), nbval[1 - 1], nbmin[1 - 1], nxval[1 - 1], nsval[1 - 1], mxbval[1 - 1];
        if (nbval[1 - 1] < 1) {
            write(nout, format_9989), "   NB ", nbval[1 - 1], 1;
            fatal = true;
        } else if (nbmin[1 - 1] < 1) {
            write(nout, format_9989), "NBMIN ", nbmin[1 - 1], 1;
            fatal = true;
        } else if (nxval[1 - 1] < 1) {
            write(nout, format_9989), "   NX ", nxval[1 - 1], 1;
            fatal = true;
        } else if (nsval[1 - 1] < 2) {
            write(nout, format_9989), "   NS ", nsval[1 - 1], 2;
            fatal = true;
        } else if (mxbval[1 - 1] < 1) {
            write(nout, format_9989), " MAXB ", mxbval[1 - 1], 1;
            fatal = true;
        }
        xlaenv(1, nbval[1 - 1]);
        xlaenv(2, nbmin[1 - 1]);
        xlaenv(3, nxval[1 - 1]);
        xlaenv(4, nsval[1 - 1]);
        xlaenv(8, mxbval[1 - 1]);
        write(nout, format_9983), "NB:   ", nbval[1 - 1];
        write(nout, format_9983), "NBMIN:", nbmin[1 - 1];
        write(nout, format_9983), "NX:   ", nxval[1 - 1];
        write(nout, format_9983), "NS:   ", nsval[1 - 1];
        write(nout, format_9983), "MAXB: ", mxbval[1 - 1];
    }
#ifdef NOTYET
    else if (!zhb && !glm && !gqr && !gsv && !csd && !lse) {
        //
        //        For the other paths, the number of parameters can be varied
        //        from the input file.  Read the number of parameter values.
        //
        read(nin, star), nparms;
        if (nparms < 1) {
            write(nout, format_9989), "NPARMS", nparms, 1;
            nparms = 0;
            fatal = true;
        } else if (nparms > maxin) {
            write(nout, format_9988), "NPARMS", nparms, maxin;
            nparms = 0;
            fatal = true;
        }
        //
        //        Read the values of NB
        //
        if (!zbb) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, nbval[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (nbval[i - 1] < 0) {
                    write(nout, format_9989), "   NB ", nbval[i - 1], 0;
                    fatal = true;
                } else if (nbval[i - 1] > nmax) {
                    write(nout, format_9988), "   NB ", nbval[i - 1], nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NB:   ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nbval[i - 1];
                }
            }
        }
        //
        //        Read the values of NBMIN
        //
        if (nep || sep || svd || zgg) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, nbmin[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (nbmin[i - 1] < 0) {
                    write(nout, format_9989), "NBMIN ", nbmin[i - 1], 0;
                    fatal = true;
                } else if (nbmin[i - 1] > nmax) {
                    write(nout, format_9988), "NBMIN ", nbmin[i - 1], nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NBMIN:";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nbmin[i - 1];
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                nbmin[i - 1] = 1;
            }
        }
        //
        //        Read the values of NX
        //
        if (nep || sep || svd) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, nxval[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (nxval[i - 1] < 0) {
                    write(nout, format_9989), "   NX ", nxval[i - 1], 0;
                    fatal = true;
                } else if (nxval[i - 1] > nmax) {
                    write(nout, format_9988), "   NX ", nxval[i - 1], nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NX:   ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nxval[i - 1];
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                nxval[i - 1] = 1;
            }
        }
        //
        //        Read the values of NSHIFT (if ZGG) or NRHS (if SVD
        //        or ZBB).
        //
        if (svd || zbb || zgg) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, nsval[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (nsval[i - 1] < 0) {
                    write(nout, format_9989), "   NS ", nsval[i - 1], 0;
                    fatal = true;
                } else if (nsval[i - 1] > nmax) {
                    write(nout, format_9988), "   NS ", nsval[i - 1], nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NS:   ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nsval[i - 1];
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                nsval[i - 1] = 1;
            }
        }
        //
        //        Read the values for MAXB.
        //
        if (zgg) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, mxbval[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (mxbval[i - 1] < 0) {
                    write(nout, format_9989), " MAXB ", mxbval[i - 1], 0;
                    fatal = true;
                } else if (mxbval[i - 1] > nmax) {
                    write(nout, format_9988), " MAXB ", mxbval[i - 1], nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "MAXB: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, mxbval[i - 1];
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                mxbval[i - 1] = 1;
            }
        }
        //
        //        Read the values for INMIN.
        //
        if (nep) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, inmin[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (inmin[i - 1] < 0) {
                    write(nout, format_9989), " INMIN ", inmin[i - 1], 0;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "INMIN: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, inmin[i - 1];
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                inmin[i - 1] = 1;
            }
        }
        //
        //        Read the values for INWIN.
        //
        if (nep) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, inwin[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (inwin[i - 1] < 0) {
                    write(nout, format_9989), " INWIN ", inwin[i - 1], 0;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "INWIN: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, inwin[i - 1];
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                inwin[i - 1] = 1;
            }
        }
        //
        //        Read the values for INIBL.
        //
        if (nep) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, inibl[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (inibl[i - 1] < 0) {
                    write(nout, format_9989), " INIBL ", inibl[i - 1], 0;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "INIBL: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, inibl[i - 1];
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                inibl[i - 1] = 1;
            }
        }
        //
        //        Read the values for ISHFTS.
        //
        if (nep) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, ishfts[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (ishfts[i - 1] < 0) {
                    write(nout, format_9989), " ISHFTS ", ishfts[i - 1], 0;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "ISHFTS: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, ishfts[i - 1];
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                ishfts[i - 1] = 1;
            }
        }
        //
        //        Read the values for IACC22.
        //
        if (nep || zgg) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, iacc22[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (iacc22[i - 1] < 0) {
                    write(nout, format_9989), " IACC22 ", iacc22[i - 1], 0;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "IACC22: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, iacc22[i - 1];
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                iacc22[i - 1] = 1;
            }
        }
        //
        //        Read the values for NBCOL.
        //
        if (zgg) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, nbcol[i - 1];
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (nbcol[i - 1] < 0) {
                    write(nout, format_9989), "NBCOL ", nbcol[i - 1], 0;
                    fatal = true;
                } else if (nbcol[i - 1] > nmax) {
                    write(nout, format_9988), "NBCOL ", nbcol[i - 1], nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NBCOL:";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nbcol[i - 1];
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                nbcol[i - 1] = 1;
            }
        }
    }
    //
    //     Calculate and print the machine dependent constants.
    //
    write(nout, star);
    eps = Rlamch("Underflow threshold");
    write(nout, format_9981), "underflow", eps;
    eps = Rlamch("Overflow threshold");
    write(nout, format_9981), "overflow ", eps;
    eps = Rlamch("Epsilon");
    write(nout, format_9981), "precision", eps;
    //
    //     Read the threshold value for the test ratios.
    //
    read(nin, star), thresh;
    write(nout, "(/,' Routines pass computational tests if test ratio is ','less than',"
                "a,/)"),
        thresh;
    if (sep || svd || zgg) {
        //
        //        Read the flag that indicates whether to test LAPACK routines.
        //
        read(nin, star), tstchk;
        //
        //        Read the flag that indicates whether to test driver routines.
        //
        read(nin, star), tstdrv;
    }
    //
    //     Read the flag that indicates whether to test the error exits.
    //
    read(nin, star), tsterr;
    //
    //     Read the code describing how to set the random number seed.
    //
    read(nin, star), newsd;
    //
    //     If NEWSD = 2, read another line with 4 integers for the seed.
    //
    if (newsd == 2) {
        {
            read_loop rloop(cmn, nin, star);
            for (i = 1; i <= 4; i = i + 1) {
                rloop, ioldsd[i - 1];
            }
        }
    }
    //
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = ioldsd[i - 1];
    }
    //
    if (fatal) {
        write(nout, "(/,' Execution not attempted due to input errors')");
        FEM_STOP(0);
    }
//
//     Read the input lines indicating the test path and its parameters.
//     The first three characters indicate the test path, and the number
//     of test matrix types must be the first nonblank item in columns
//     4-80.
//
statement_190:
    //
    if (!(zgx || zxv)) {
    //
    statement_200:
        try {
            read(nin, "(a80)"), line;
        } catch (read_end const) {
            goto statement_380;
        }
        c3 = line[(3 - 1) * ldline];
        lenp = len[line - 1];
        i = 3;
        itmp = 0;
        i1 = 0;
    statement_210:
        i++;
        if (i > lenp) {
            if (i1 > 0) {
                goto statement_240;
            } else {
                ntypes = maxt;
                goto statement_240;
            }
        }
        if (line[(i - 1) + (i - 1) * ldline] != " " && line[(i - 1) + (i - 1) * ldline] != ",") {
            i1 = i;
            c1 = line[(i1 - 1) + (i1 - 1) * ldline];
            //
            //        Check that a valid integer was read
            //
            for (k = 1; k <= 10; k = k + 1) {
                if (c1 == intstr[(k - 1) + (k - 1) * ldintstr]) {
                    ic = k - 1;
                    goto statement_230;
                }
            }
            write(nout, "(/,/,' *** Invalid integer value in column ',i2,' of input',' line:',"
                        "/,a79)"),
                i, line;
            goto statement_200;
        statement_230:
            itmp = 10 * itmp + ic;
            goto statement_210;
        } else if (i1 > 0) {
            goto statement_240;
        } else {
            goto statement_210;
        }
    statement_240:
        ntypes = itmp;
        //
        //     Skip the tests if NTYPES is <= 0.
        //
        if (!(zev || zes || zvx || zsx || zgv || zgs) && ntypes <= 0) {
            write(nout, format_9990), c3;
            goto statement_200;
        }
        //
    } else {
        if (zgx) {
            c3 = "ZGX";
        }
        if (zxv) {
            c3 = "ZXV";
        }
    }
    //
    //     Reset the random number seed.
    //
    if (newsd == 0) {
        for (k = 1; k <= 4; k = k + 1) {
            iseed[k - 1] = ioldsd[k - 1];
        }
    }
    //
    if (Mlsamen(3, c3, "ZHS") || Mlsamen(3, c3, "NEP")) {
        //
        //        -------------------------------------
        //        NEP:  Nonsymmetric Eigenvalue Problem
        //        -------------------------------------
        //        Vary the parameters
        //           NB    = block size
        //           NBMIN = minimum block size
        //           NX    = crossover point
        //           NS    = number of shifts
        //           MAXB  = minimum submatrix size
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        xlaenv(1, 1);
        if (tsterr) {
            Cerrhs("Chseqr", nout);
        }
        for (i = 1; i <= nparms; i = i + 1) {
            xlaenv(1, nbval[i - 1]);
            xlaenv(2, nbmin[i - 1]);
            xlaenv(3, nxval[i - 1]);
            xlaenv(12, max((INTEGER)11, inmin[i - 1]));
            xlaenv(13, inwin[i - 1]);
            xlaenv(14, inibl[i - 1]);
            xlaenv(15, ishfts[i - 1]);
            xlaenv(16, iacc22[i - 1]);
            //
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4,', INMIN=',i4,"
                        "', INWIN =',i4,', INIBL =',i4,', ISHFTS =',i4,', IACC22 =',i4)"),
                c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], max((INTEGER)11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
            Cchkhs(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], nmax, &a[(1 - 1) + (6 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], &a[(1 - 1) + (8 - 1) * lda], &a[(1 - 1) + (9 - 1) * lda], &a[(1 - 1) + (10 - 1) * lda], &a[(1 - 1) + (11 - 1) * lda], &a[(1 - 1) + (12 - 1) * lda], d & c[(1 - 1) + (3 - 1) * ldc], work, lwork, rwork, iwork, logwrk, result, info);
            if (info != 0) {
                write(nout, format_9980), "Cchkhs", info;
            }
        }
        //
    } else if (Mlsamen(3, c3, "ZST") || Mlsamen(3, c3, "SEP") || Mlsamen(3, c3, "SE2")) {
        //
        //        ----------------------------------
        //        SEP:  Symmetric Eigenvalue Problem
        //        ----------------------------------
        //        Vary the parameters
        //           NB    = block size
        //           NBMIN = minimum block size
        //           NX    = crossover point
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        xlaenv(1, 1);
        xlaenv(9, 25);
        if (tsterr) {
            Cerrst("ZST", nout);
        }
        for (i = 1; i <= nparms; i = i + 1) {
            xlaenv(1, nbval[i - 1]);
            xlaenv(2, nbmin[i - 1]);
            xlaenv(3, nxval[i - 1]);
            //
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, format_9997), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1];
            if (tstchk) {
                if (Mlsamen(3, c3, "SE2")) {
                    Cchkst2stg(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], dr[(1 - 1)], dr[(2 - 1) * lddr], dr[(3 - 1) * lddr], dr[(4 - 1) * lddr], dr[(5 - 1) * lddr], dr[(6 - 1) * lddr], dr[(7 - 1) * lddr], dr[(8 - 1) * lddr], dr[(9 - 1) * lddr], dr[(10 - 1) * lddr], dr[(11 - 1) * lddr], &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], &a[(1 - 1) + (6 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
                } else {
                    Cchkst(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], dr[(1 - 1)], dr[(2 - 1) * lddr], dr[(3 - 1) * lddr], dr[(4 - 1) * lddr], dr[(5 - 1) * lddr], dr[(6 - 1) * lddr], dr[(7 - 1) * lddr], dr[(8 - 1) * lddr], dr[(9 - 1) * lddr], dr[(10 - 1) * lddr], dr[(11 - 1) * lddr], &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], &a[(1 - 1) + (6 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
                }
                if (info != 0) {
                    write(nout, format_9980), "Cchkst", info;
                }
            }
            if (tstdrv) {
                if (Mlsamen(3, c3, "SE2")) {
                    Cdrvst2stg(nn, nval, 18, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, dr[(3 - 1) * lddr], dr[(4 - 1) * lddr], dr[(5 - 1) * lddr], dr[(8 - 1) * lddr], dr[(9 - 1) * lddr], dr[(10 - 1) * lddr], &a[(1 - 1) + (2 - 1) * lda], nmax, &a[(1 - 1) + (3 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], &a[(1 - 1) + (4 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
                } else {
                    Cdrvst(nn, nval, 18, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, dr[(3 - 1) * lddr], dr[(4 - 1) * lddr], dr[(5 - 1) * lddr], dr[(8 - 1) * lddr], dr[(9 - 1) * lddr], dr[(10 - 1) * lddr], &a[(1 - 1) + (2 - 1) * lda], nmax, &a[(1 - 1) + (3 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], &a[(1 - 1) + (4 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
                }
                if (info != 0) {
                    write(nout, format_9980), "Cdrvst", info;
                }
            }
        }
        //
    } else if (Mlsamen(3, c3, "ZSG")) {
        //
        //        ----------------------------------------------
        //        ZSG:  Hermitian Generalized Eigenvalue Problem
        //        ----------------------------------------------
        //        Vary the parameters
        //           NB    = block size
        //           NBMIN = minimum block size
        //           NX    = crossover point
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        xlaenv(9, 25);
        for (i = 1; i <= nparms; i = i + 1) {
            xlaenv(1, nbval[i - 1]);
            xlaenv(2, nbmin[i - 1]);
            xlaenv(3, nxval[i - 1]);
            //
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, format_9997), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1];
            if (tstchk) {
                //               CALL Cdrvsg( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
                //     $                      NOUT, &a[(1-1)+(1-1)*lda], NMAX, &a[(1-1)+(2-1)*lda], NMAX,
                //     $                      DR( 1, 3 ), &a[(1-1)+(3-1)*lda], NMAX, &a[(1-1)+(4-1)*lda],
                //     $                      &a[(1-1)+(5-1)*lda], &a[(1-1)+(6-1)*lda], &a[(1-1)+(7-1)*lda], WORK,
                //     $                      LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT,
                //     $                      INFO )
                Cdrvsg2stg(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], nmax, dr[(3 - 1) * lddr], dr[(4 - 1) * lddr], &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
                if (info != 0) {
                    write(nout, format_9980), "Cdrvsg", info;
                }
            }
        }
        //
    } else if (Mlsamen(3, c3, "ZBD") || Mlsamen(3, c3, "SVD")) {
        //
        //        ----------------------------------
        //        SVD:  Singular Value Decomposition
        //        ----------------------------------
        //        Vary the parameters
        //           NB    = block size
        //           NBMIN = minimum block size
        //           NX    = crossover point
        //           NRHS  = number of right hand sides
        //
        maxtyp = 16;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        xlaenv(9, 25);
        //
        //        Test the error exits
        //
        xlaenv(1, 1);
        if (tsterr && tstchk) {
            Cerrbd("ZBD", nout);
        }
        if (tsterr && tstdrv) {
            Cerred("ZBD", nout);
        }
        //
        for (i = 1; i <= nparms; i = i + 1) {
            nrhs = nsval[i - 1];
            xlaenv(1, nbval[i - 1]);
            xlaenv(2, nbmin[i - 1]);
            xlaenv(3, nxval[i - 1]);
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4,', NRHS =',i4)"), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], nrhs;
            if (tstchk) {
                Cchkbd(nn, mval, nval, maxtyp, dotype, nrhs, iseed, thresh, &a[(1 - 1) + (1 - 1) * lda], nmax, dr[(1 - 1)], dr[(2 - 1) * lddr], dr[(3 - 1) * lddr], dr[(4 - 1) * lddr], &a[(1 - 1) + (2 - 1) * lda], nmax, &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], nmax, &a[(1 - 1) + (6 - 1) * lda], nmax, &a[(1 - 1) + (7 - 1) * lda], &a[(1 - 1) + (8 - 1) * lda], work, lwork, rwork, nout, info);
                if (info != 0) {
                    write(nout, format_9980), "Cchkbd", info;
                }
            }
            if (tstdrv) {
                Cdrvbd(nn, mval, nval, maxtyp, dotype, iseed, thresh, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], nmax, &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], dr[(1 - 1)], dr[(2 - 1) * lddr], dr[(3 - 1) * lddr], work, lwork, rwork, iwork, nout, info);
            }
        }
        //
    } else if (Mlsamen(3, c3, "ZEV")) {
        //
        //        --------------------------------------------
        //        ZEV:  Nonsymmetric Eigenvalue Problem Driver
        //              Cgeev (eigenvalues and eigenvectors)
        //        --------------------------------------------
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        if (ntypes <= 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerred(c3, nout);
            }
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrvev(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], nmax, &a[(1 - 1) + (5 - 1) * lda], nmax, result, work, lwork, rwork, iwork, info);
            if (info != 0) {
                write(nout, format_9980), "Cgeev", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3, "ZES")) {
        //
        //        --------------------------------------------
        //        ZES:  Nonsymmetric Eigenvalue Problem Driver
        //              Cgees (Schur form)
        //        --------------------------------------------
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        if (ntypes <= 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerred(c3, nout);
            }
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrves(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], &a[(1 - 1) + (4 - 1) * lda], nmax, result, work, lwork, rwork, iwork, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Cgees", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3, "ZVX")) {
        //
        //        --------------------------------------------------------------
        //        ZVX:  Nonsymmetric Eigenvalue Problem Expert Driver
        //              Cgeevx (eigenvalues, eigenvectors and condition numbers)
        //        --------------------------------------------------------------
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        if (ntypes < 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerred(c3, nout);
            }
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrvvx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], nmax, &a[(1 - 1) + (5 - 1) * lda], nmax, dr[(1 - 1)], dr[(2 - 1) * lddr], dr[(3 - 1) * lddr], dr[(4 - 1) * lddr], dr[(5 - 1) * lddr], dr[(6 - 1) * lddr], dr[(7 - 1) * lddr], dr[(8 - 1) * lddr], result, work, lwork, rwork, info);
            if (info != 0) {
                write(nout, format_9980), "Cgeevx", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3, "ZSX")) {
        //
        //        ---------------------------------------------------
        //        ZSX:  Nonsymmetric Eigenvalue Problem Expert Driver
        //              Cgeesx (Schur form and condition numbers)
        //        ---------------------------------------------------
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        if (ntypes < 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerred(c3, nout);
            }
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrvsx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], d & c[(1 - 1) + (3 - 1) * ldc], &a[(1 - 1) + (4 - 1) * lda], nmax, &a[(1 - 1) + (5 - 1) * lda], result, work, lwork, rwork, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Cgeesx", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3, "ZGG")) {
        //
        //        -------------------------------------------------
        //        ZGG:  Generalized Nonsymmetric Eigenvalue Problem
        //        -------------------------------------------------
        //        Vary the parameters
        //           NB    = block size
        //           NBMIN = minimum block size
        //           NS    = number of shifts
        //           MAXB  = minimum submatrix size
        //           IACC22: structured matrix multiply
        //
        maxtyp = 26;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        xlaenv(1, 1);
        if (tstchk && tsterr) {
            Cerrgg(c3, nout);
        }
        for (i = 1; i <= nparms; i = i + 1) {
            xlaenv(1, nbval[i - 1]);
            xlaenv(2, nbmin[i - 1]);
            xlaenv(4, nsval[i - 1]);
            xlaenv(8, mxbval[i - 1]);
            xlaenv(16, iacc22[i - 1]);
            xlaenv(5, nbcol[i - 1]);
            //
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NS =',i4,', MAXB =',i4,"
                        "', IACC22 =',i4,', NBCOL =',i4)"),
                c3, nbval[i - 1], nbmin[i - 1], nsval[i - 1], mxbval[i - 1], iacc22[i - 1], nbcol[i - 1];
            tstdif = false;
            thrshn = 10.0;
            if (tstchk) {
                Cchkgg(nn, nval, maxtyp, dotype, iseed, thresh, tstdif, thrshn, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], &a[(1 - 1) + (8 - 1) * lda], &a[(1 - 1) + (9 - 1) * lda], nmax, &a[(1 - 1) + (10 - 1) * lda], &a[(1 - 1) + (11 - 1) * lda], &a[(1 - 1) + (12 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], d & c[(1 - 1) + (3 - 1) * ldc], d & c[(1 - 1) + (4 - 1) * ldc], &a[(1 - 1) + (13 - 1) * lda], &a[(1 - 1) + (14 - 1) * lda], work, lwork, rwork, logwrk, result, info);
                if (info != 0) {
                    write(nout, format_9980), "Cchkgg", info;
                }
            }
        }
        //
    } else if (Mlsamen(3, c3, "ZGS")) {
        //
        //        -------------------------------------------------
        //        ZGS:  Generalized Nonsymmetric Eigenvalue Problem
        //              Cgges (Schur form)
        //        -------------------------------------------------
        //
        maxtyp = 26;
        ntypes = min(maxtyp, ntypes);
        if (ntypes <= 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerrgg(c3, nout);
            }
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrges(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], nmax, &a[(1 - 1) + (8 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], work, lwork, rwork, result, logwrk, info);
            //
            if (info != 0) {
                write(nout, format_9980), "Cdrges", info;
            }
            //
            // Blocked version
            //
            Cdrges3(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], nmax, &a[(1 - 1) + (8 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], work, lwork, rwork, result, logwrk, info);
            //
            if (info != 0) {
                write(nout, format_9980), "Cdrges3", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (zgx) {
        //
        //        -------------------------------------------------
        //        ZGX  Generalized Nonsymmetric Eigenvalue Problem
        //              Cggesx (Schur form and condition numbers)
        //        -------------------------------------------------
        //
        maxtyp = 5;
        ntypes = maxtyp;
        if (nn < 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerrgg(c3, nout);
            }
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            xlaenv(5, 2);
            Cdrgsx(nn, ncmax, thresh, nin, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], c, ncmax * ncmax, s, work, lwork, rwork, iwork, liwork, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Cdrgsx", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3, "ZGV")) {
        //
        //        -------------------------------------------------
        //        ZGV:  Generalized Nonsymmetric Eigenvalue Problem
        //              Cggev (Eigenvalue/vector form)
        //        -------------------------------------------------
        //
        maxtyp = 26;
        ntypes = min(maxtyp, ntypes);
        if (ntypes <= 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerrgg(c3, nout);
            }
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrgev(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], nmax, &a[(1 - 1) + (8 - 1) * lda], &a[(1 - 1) + (9 - 1) * lda], nmax, d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], d & c[(1 - 1) + (3 - 1) * ldc], d & c[(1 - 1) + (4 - 1) * ldc], work, lwork, rwork, result, info);
            if (info != 0) {
                write(nout, format_9980), "Cdrgev", info;
            }
            //
            // Blocked version
            //
            xlaenv(16, 2);
            Cdrgev3(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], nmax, &a[(1 - 1) + (8 - 1) * lda], &a[(1 - 1) + (9 - 1) * lda], nmax, d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], d & c[(1 - 1) + (3 - 1) * ldc], d & c[(1 - 1) + (4 - 1) * ldc], work, lwork, rwork, result, info);
            if (info != 0) {
                write(nout, format_9980), "Cdrgev3", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (zxv) {
        //
        //        -------------------------------------------------
        //        ZXV:  Generalized Nonsymmetric Eigenvalue Problem
        //              Cggevx (eigenvalue/vector with condition numbers)
        //        -------------------------------------------------
        //
        maxtyp = 2;
        ntypes = maxtyp;
        if (nn < 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerrgg(c3, nout);
            }
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrgvx(nn, thresh, nin, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], d & c[(1 - 1) + (1 - 1) * ldc], d & c[(1 - 1) + (2 - 1) * ldc], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], &iwork[1 - 1], &iwork[2 - 1], dr[(1 - 1)], dr[(2 - 1) * lddr], dr[(3 - 1) * lddr], dr[(4 - 1) * lddr], dr[(5 - 1) * lddr], dr[(6 - 1) * lddr], work, lwork, rwork, &iwork[3 - 1], liwork - 2, result, logwrk, info);
            //
            if (info != 0) {
                write(nout, format_9980), "Cdrgvx", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3, "ZHB")) {
        //
        //        ------------------------------
        //        ZHB:  Hermitian Band Reduction
        //        ------------------------------
        //
        maxtyp = 15;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        if (tsterr) {
            Cerrst("ZHB", nout);
        }
        //         CALL Cchkhb( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
        //     $                NOUT, &a[(1-1)+(1-1)*lda], NMAX, DR( 1, 1 ), DR( 1, 2 ),
        //     $                &a[(1-1)+(2-1)*lda], NMAX, WORK, LWORK, RWORK, RESULT,
        //     $                INFO )
        Cchkhb2stg(nn, nval, nk, kval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, dr[(1 - 1)], dr[(2 - 1) * lddr], dr[(3 - 1) * lddr], dr[(4 - 1) * lddr], dr[(5 - 1) * lddr], &a[(1 - 1) + (2 - 1) * lda], nmax, work, lwork, rwork, result, info);
        if (info != 0) {
            write(nout, format_9980), "Cchkhb", info;
        }
        //
    } else if (Mlsamen(3, c3, "ZBB")) {
        //
        //        ------------------------------
        //        ZBB:  General Band Reduction
        //        ------------------------------
        //
        maxtyp = 15;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        for (i = 1; i <= nparms; i = i + 1) {
            nrhs = nsval[i - 1];
            //
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, "(/,/,1x,a3,':  NRHS =',i4)"), c3, nrhs;
            Cchkbb(nn, mval, nval, nk, kval, maxtyp, dotype, nrhs, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], 2 * nmax, dr[(1 - 1)], dr[(2 - 1) * lddr], &a[(1 - 1) + (4 - 1) * lda], nmax, &a[(1 - 1) + (5 - 1) * lda], nmax, &a[(1 - 1) + (6 - 1) * lda], nmax, &a[(1 - 1) + (7 - 1) * lda], work, lwork, rwork, result, info);
            if (info != 0) {
                write(nout, format_9980), "Cchkbb", info;
            }
        }
        //
    } else if (Mlsamen(3, c3, "GLM")) {
        //
        //        -----------------------------------------
        //        GLM:  Generalized Linear Regression Model
        //        -----------------------------------------
        //
        xlaenv(1, 1);
        if (tsterr) {
            Cerrgg("GLM", nout);
        }
        Cckglm(nn, nval, mval, pval, ntypes, iseed, thresh, nmax, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &b[(1 - 1) + (1 - 1) * ldb], &b[(1 - 1) + (2 - 1) * ldb], x, work, dr[(1 - 1)], nin, nout, info);
        if (info != 0) {
            write(nout, format_9980), "Cckglm", info;
        }
        //
    } else if (Mlsamen(3, c3, "GQR")) {
        //
        //        ------------------------------------------
        //        GQR:  Generalized QR and RQ factorizations
        //        ------------------------------------------
        //
        xlaenv(1, 1);
        if (tsterr) {
            Cerrgg("GQR", nout);
        }
        Cckgqr(nn, mval, nn, pval, nn, nval, ntypes, iseed, thresh, nmax, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], taua, &b[(1 - 1) + (1 - 1) * ldb], &b[(1 - 1) + (2 - 1) * ldb], &b[(1 - 1) + (3 - 1) * ldb], &b[(1 - 1) + (4 - 1) * ldb], &b[(1 - 1) + (5 - 1) * ldb], taub, work, dr[(1 - 1)], nin, nout, info);
        if (info != 0) {
            write(nout, format_9980), "Cckgqr", info;
        }
        //
    } else if (Mlsamen(3, c3, "GSV")) {
        //
        //        ----------------------------------------------
        //        GSV:  Generalized Singular Value Decomposition
        //        ----------------------------------------------
        //
        xlaenv(1, 1);
        if (tsterr) {
            Cerrgg("GSV", nout);
        }
        Cckgsv(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &b[(1 - 1) + (1 - 1) * ldb], &b[(1 - 1) + (2 - 1) * ldb], &a[(1 - 1) + (3 - 1) * lda], &b[(1 - 1) + (3 - 1) * ldb], &a[(1 - 1) + (4 - 1) * lda], alpha, beta, &b[(1 - 1) + (4 - 1) * ldb], iwork, work, dr[(1 - 1)], nin, nout, info);
        if (info != 0) {
            write(nout, format_9980), "Cckgsv", info;
        }
        //
    } else if (Mlsamen(3, c3, "CSD")) {
        //
        //        ----------------------------------------------
        //        CSD:  CS Decomposition
        //        ----------------------------------------------
        //
        xlaenv(1, 1);
        if (tsterr) {
            Cerrgg("CSD", nout);
        }
        Cckcsd(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], rwork, iwork, work, dr[(1 - 1)], nin, nout, info);
        if (info != 0) {
            write(nout, format_9980), "Cckcsd", info;
        }
        //
    } else if (Mlsamen(3, c3, "LSE")) {
        //
        //        --------------------------------------
        //        LSE:  Constrained Linear Least Squares
        //        --------------------------------------
        //
        xlaenv(1, 1);
        if (tsterr) {
            Cerrgg("LSE", nout);
        }
        Ccklse(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &b[(1 - 1) + (1 - 1) * ldb], &b[(1 - 1) + (2 - 1) * ldb], x, work, dr[(1 - 1)], nin, nout, info);
        if (info != 0) {
            write(nout, format_9980), "Ccklse", info;
        }
    } else {
        write(nout, star);
        write(nout, star);
        write(nout, format_9992), c3;
    }
    if (!(zgx || zxv)) {
        goto statement_190;
    }
#endif
statement_380:
    write(nout, "(/,/,' End of tests')");
    s2 = time(NULL);
    write(nout, "(' Total time used = ',f12.2,' seconds',/)"), double(s2 - s1);
    //
    //
    //     End of Cchkee
    //
}

int main(int argc, char const *argv[]) { Cchkee(); }
