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
#include <string>
#include <vector>
using namespace std;

void Rchkee(void) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    INTEGER allocatestatus = 0;
    INTEGER need = 14;
    const INTEGER nmax = 132;
    REAL d[nmax * 12];
    INTEGER ldd = nmax;
    time_t s1;
    bool fatal = false;
    const INTEGER nout = 6;
    const INTEGER nin = 5;
    char path[3];
    bool nep = false;
    bool sep = false;
    bool svd = false;
    bool dev = false;
    bool des = false;
    bool dvx = false;
    bool dsx = false;
    bool dgg = false;
    bool dgs = false;
    bool dgx = false;
    bool dgv = false;
    bool dxv = false;
    bool dsb = false;
    bool dbb = false;
    bool glm = false;
    bool gqr = false;
    bool gsv = false;
    bool csd = false;
    bool lse = false;
    bool dbl = false;
    bool dbk = false;
    bool dgl = false;
    bool dgk = false;
    REAL thresh = 0.0;
    bool tsterr = false;
    INTEGER vers_major = 0;
    INTEGER vers_minor = 0;
    INTEGER vers_patch = 0;
    INTEGER nn = 0;
    const INTEGER maxin = 20;
    INTEGER mval[maxin];
    INTEGER i = 0;
    char vname[32];
    INTEGER pval[maxin];
    INTEGER nval[maxin];
    INTEGER nk = 0;
    INTEGER kval[maxin];
    INTEGER nbval[maxin];
    INTEGER nbmin[maxin];
    INTEGER nxval[maxin];
    INTEGER inmin[maxin];
    INTEGER inwin[maxin];
    INTEGER inibl[maxin];
    INTEGER ishfts[maxin];
    INTEGER iacc22[maxin];
    INTEGER nsval[maxin];
    INTEGER mxbval[maxin];
    INTEGER nparms = 0;
    INTEGER nbcol[maxin];
    REAL eps = 0.0;
    bool tstchk = false;
    bool tstdrv = false;
    INTEGER newsd = 0;
    INTEGER iseed[4];
    INTEGER ioldsd[4];
    char c3[3];
    INTEGER lenp = 0;
    INTEGER i1 = 0;
    const INTEGER maxt = 30;
    INTEGER ntypes = 0;
    char c1;
    INTEGER k = 0;
    INTEGER ic = 0;
    INTEGER maxtyp = 0;
    bool dotype[maxt];
    const INTEGER lwork = nmax * (5 * nmax + 5) + 1;
    const INTEGER liwork = nmax * (5 * nmax + 20);
    INTEGER iwork[liwork];
    bool logwrk[nmax];
    INTEGER ldb = nmax;
    REAL result[500];
    INTEGER info = 0;
    INTEGER nrhs = 0;
    bool tstdif = false;
    REAL thrshn = 0.0;
    const INTEGER ncmax = 20;
    REAL x[5 * nmax];
    REAL taua[nmax];
    REAL taub[nmax];
    time_t s2;
    INTEGER mplapack_vers_major = 0;
    INTEGER mplapack_vers_minor = 0;
    INTEGER mplapack_vers_patch = 0;
    INTEGER lapack_vers_major = 0;
    INTEGER lapack_vers_minor = 0;
    INTEGER lapack_vers_patch = 0;
    static const char *format_9973 = "(/,1x,71('-'))";
    static const char *format_9980 = "(' *** Error code from ',a,' = ',i4)";
    static const char *format_9981 = "(' Relative machine ',a,' is taken to be',a)";
    static const char *format_9983 = "(4x,a,10i6,/,10x,10i6)";
    static const char *format_9988 = "(' Invalid input value: ',a,'=',i6,'; must be <=',i6)";
    static const char *format_9989 = "(' Invalid input value: ',a,'=',i6,'; must be >=',i6)";
    static const char *format_9990 = "(/,/,1x,a3,' routines were not tested')";
    static const char *format_9992 = "(1x,a3,':  Unrecognized path name')";
    static const char *format_9997 = "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4)";
    //
    //     .. Allocate memory dynamically ..
    //
    REAL *a = new REAL[nmax * nmax * need];
    INTEGER lda = nmax;
    REAL *b = new REAL[nmax * nmax * 5];
    REAL *c = new REAL[ncmax * ncmax * ncmax * ncmax];
    INTEGER ldc = ncmax;
    REAL *work = new REAL[lwork];
    //
    for (int ppp = 0; ppp < nmax * nmax * need; ppp++)
        a[ppp] = 0.0;
    for (int ppp = 0; ppp < nmax * nmax * 5; ppp++)
        b[ppp] = 0.0;
    for (int ppp = 0; ppp < ncmax * ncmax * ncmax * ncmax; ppp++)
        c[ppp] = 0.0;
    for (int ppp = 0; ppp < lwork; ppp++)
        work[ppp] = 0.0;
    //
    string str;
    istringstream iss;
    char line[1024];
    char buf0[1024];
    double dtmp;
    int itmp;
    //
    s1 = time(NULL);
    fatal = false;
    //
    //     Return to here to read multiple sets of data
    //
    while (1) {
        getline(cin, str);
        if (cin.bad() || cin.eof()) {
            break;
        }
        //
        //     Read the first line and set the 3-character test path
        //
        stringstream ss(str);
        ss >> line;
        path[0] = line[0];
        path[1] = line[1];
        path[2] = line[2];
        nep = Mlsamen(3, path, "NEP") || Mlsamen(3, path, "DHS");
        sep = Mlsamen(3, path, "SEP") || Mlsamen(3, path, "DST") || Mlsamen(3, path, "DSG") || Mlsamen(3, path, "SE2");
        svd = Mlsamen(3, path, "SVD") || Mlsamen(3, path, "DBD");
        dev = Mlsamen(3, path, "DEV");
        des = Mlsamen(3, path, "DES");
        dvx = Mlsamen(3, path, "DVX");
        dsx = Mlsamen(3, path, "DSX");
        dgg = Mlsamen(3, path, "DGG");
        dgs = Mlsamen(3, path, "DGS");
        dgx = Mlsamen(3, path, "DGX");
        dgv = Mlsamen(3, path, "DGV");
        dxv = Mlsamen(3, path, "DXV");
        dsb = Mlsamen(3, path, "DSB");
        dbb = Mlsamen(3, path, "DBB");
        glm = Mlsamen(3, path, "GLM");
        gqr = Mlsamen(3, path, "GQR") || Mlsamen(3, path, "GRQ");
        gsv = Mlsamen(3, path, "GSV");
        csd = Mlsamen(3, path, "CSD");
        lse = Mlsamen(3, path, "LSE");
        dbl = Mlsamen(3, path, "DBL");
        dbk = Mlsamen(3, path, "DBK");
        dgl = Mlsamen(3, path, "DGL");
        dgk = Mlsamen(3, path, "DGK");
        //
        //     Report values of parameters.
        //
        if (path == "   ") {
            continue;
        } else if (nep) {
            write(nout, "(' Tests of the Nonsymmetric Eigenvalue Problem routines')");
        } else if (sep) {
            write(nout, "(' Tests of the Symmetric Eigenvalue Problem routines')");
        } else if (svd) {
            write(nout, "(' Tests of the Singular Value Decomposition routines')");
        } else if (dev) {
            write(nout, "(/,' Tests of the Nonsymmetric Eigenvalue Problem Driver',/,"
                        "'    Rgeev (eigenvalues and eigevectors)')");
        } else if (des) {
            write(nout, "(/,' Tests of the Nonsymmetric Eigenvalue Problem Driver',/,"
                        "'    Rgees (Schur form)')");
        } else if (dvx) {
            write(nout, "(/,' Tests of the Nonsymmetric Eigenvalue Problem Expert',' Driver',/,"
                        "'    Rgeevx (eigenvalues, eigenvectors and',' condition numbers)')");
        } else if (dsx) {
            write(nout, "(/,' Tests of the Nonsymmetric Eigenvalue Problem Expert',' Driver',/,"
                        "'    Rgeesx (Schur form and condition',' numbers)')");
        } else if (dgg) {
            write(nout, "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                        "'Problem routines')");
        } else if (dgs) {
            write(nout, "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                        "'Problem Driver Rgges')");
        } else if (dgx) {
            write(nout, "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                        "'Problem Expert Driver Rggesx')");
        } else if (dgv) {
            write(nout, "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                        "'Problem Driver Rggev')");
        } else if (dxv) {
            write(nout, "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                        "'Problem Expert Driver Rggevx')");
        } else if (dsb) {
            write(nout, "(' Tests of Rsbtrd',/,' (reduction of a symmetric band ',"
                        "'matrix to tridiagonal form)')");
        } else if (dbb) {
            write(nout, "(' Tests of Rgbbrd',/,' (reduction of a general band ',"
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
        } else if (dbl) {
            //
            //        Rgebal:  Balancing
            //
            Rchkbl(nin, nout);
            continue;
        } else if (dbk) {
            //
            //        Rgebak:  Back transformation
            //
            Rchkbk(nin, nout);
            continue;
        } else if (dgl) {
            //
            //        Rggbal:  Balancing
            //
            Rchkgl(nin, nout);
            continue;
        } else if (dgk) {
            //
            //        Rggbak:  Back transformation
            //
            Rchkgk(nin, nout);
            continue;
        } else if (Mlsamen(3, path, "DEC")) {
            //
            //        DEC:  Eigencondition estimation
            //
            getline(cin, str);
            ss.str(str);
            ss >> thresh;
            xlaenv(1, 1);
            xlaenv(12, 11);
            xlaenv(13, 2);
            xlaenv(14, 0);
            xlaenv(15, 2);
            xlaenv(16, 2);
            tsterr = true;
            Rchkec(thresh, tsterr, nin, nout);
            continue;
        } else {
            write(nout, format_9992), path;
            continue;
        }
        iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
        write(nout, "(' Tests of the Multiple precision version of LAPACK MPLAPACK VERSION ',i1,'.',i1,'.',i1,/, "
                    "' Based on original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, 'The following parameter values will be used:')"),
            mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
        //
        //     Read the number of values of M, P, and N.
        //
        getline(cin, str);
        iss.clear();
        iss.str(str);
        iss >> nn;
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
        if (!(dgx || dxv)) {
            getline(cin, str);
            iss.clear();
            iss.str(str);
            for (i = 1; i <= nn; i = i + 1) {
                iss >> mval[i - 1];
            }
            if (svd) {
                strncpy(vname, "    M ", 8);
            } else {
                strncpy(vname, "    N ", 8);
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
            write_loop wloop(cmn, nout, format_9983);
            wloop, "M:    ";
            for (i = 1; i <= nn; i = i + 1) {
                wloop, mval[i - 1];
            }
        }
        //
        //     Read the values of P
        //
        if (glm || gqr || gsv || csd || lse) {
            getline(cin, str);
            iss.clear();
            iss.str(str);
            for (i = 1; i <= nn; i = i + 1) {
                iss >> itmp;
                pval[i - 1] = itmp;
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
            write_loop wloop(cmn, nout, format_9983);
            wloop, "P:    ";
            for (i = 1; i <= nn; i = i + 1) {
                wloop, pval[i - 1];
            }
        }
        //
        //     Read the values of N
        //
        if (svd || dbb || glm || gqr || gsv || csd || lse) {
            getline(cin, str);
            iss.clear();
            iss.str(str);
            for (i = 1; i <= nn; i = i + 1) {
                iss >> itmp;
                pval[i - 1] = itmp;
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
        if (!(dgx || dxv)) {
            write_loop wloop(cmn, nout, format_9983);
            wloop, "N:    ";
            for (i = 1; i <= nn; i = i + 1) {
                wloop, nval[i - 1];
            }
        } else {
            write(nout, format_9983), "N:    ", nn;
        }
        //
        //     Read the number of values of K, followed by the values of K
        //
        if (dsb || dbb) {
            getline(cin, str);
            iss.clear();
            iss.str(str);
            iss >> nk;
            getline(cin, str);
            iss.clear();
            iss.str(str);
            for (i = 1; i <= nk; i = i + 1) {
                iss >> itmp;
                kval[i - 1] = itmp;
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
            write_loop wloop(cmn, nout, format_9983);
            wloop, "K:    ";
            for (i = 1; i <= nk; i = i + 1) {
                wloop, kval[i - 1];
            }
        }
        //
        if (dev || des || dvx || dsx) {
            //
            //        For the nonsymmetric QR driver routines, only one set of
            //        parameters is allowed.
            //
            getline(cin, str);
            iss.clear();
            iss.str(str);
            iss >> nbval[1 - 1];
            iss >> nbmin[1 - 1];
            iss >> nxval[1 - 1];
            iss >> inmin[1 - 1];
            iss >> inwin[1 - 1];
            iss >> inibl[1 - 1];
            iss >> ishfts[1 - 1];
            iss >> iacc22[1 - 1];
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
        } else if (dgs || dgx || dgv || dxv) {
            //
            //        For the nonsymmetric generalized driver routines, only one set
            //        of parameters is allowed.
            //
            getline(cin, str);
            iss.clear();
            iss.str(str);
            iss >> nbval[1 - 1];
            iss >> nbmin[1 - 1];
            iss >> nxval[1 - 1];
            iss >> nsval[1 - 1];
            iss >> mxbval[1 - 1];
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
        } else if (!dsb && !glm && !gqr && !gsv && !csd && !lse) {
            //
            //        For the other paths, the number of parameters can be varied
            //        from the input file.  Read the number of parameter values.
            //
            getline(cin, str);
            iss.clear();
            iss.str(str);
            iss >> nparms;
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
            if (!dbb) {
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> itmp;
                    nbval[i - 1] = itmp;
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
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NB:   ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nbval[i - 1];
                }
            }
            //
            //        Read the values of NBMIN
            //
            if (nep || sep || svd || dgg) {
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> nbmin[i - 1];
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
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NBMIN:";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nbmin[i - 1];
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
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> nxval[i - 1];
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
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NX:   ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nxval[i - 1];
                }
            } else {
                for (i = 1; i <= nparms; i = i + 1) {
                    nxval[i - 1] = 1;
                }
            }
            //
            //        Read the values of NSHIFT (if DGG) or NRHS (if SVD
            //        or DBB).
            //
            if (svd || dbb || dgg) {
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> nsval[i - 1];
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
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NS:   ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nsval[i - 1];
                }
            } else {
                for (i = 1; i <= nparms; i = i + 1) {
                    nsval[i - 1] = 1;
                }
            }
            //
            //        Read the values for MAXB.
            //
            if (dgg) {
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> mxbval[i - 1];
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
                write_loop wloop(cmn, nout, format_9983);
                wloop, "MAXB: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, mxbval[i - 1];
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
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> inmin[i - 1];
                }
                for (i = 1; i <= nparms; i = i + 1) {
                    if (inmin[i - 1] < 0) {
                        write(nout, format_9989), " INMIN ", inmin[i - 1], 0;
                        fatal = true;
                    }
                }
                write_loop wloop(cmn, nout, format_9983);
                wloop, "INMIN: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, inmin[i - 1];
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
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> inwin[i - 1];
                }
                for (i = 1; i <= nparms; i = i + 1) {
                    if (inwin[i - 1] < 0) {
                        write(nout, format_9989), " INWIN ", inwin[i - 1], 0;
                        fatal = true;
                    }
                }
                write_loop wloop(cmn, nout, format_9983);
                wloop, "INWIN: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, inwin[i - 1];
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
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> inibl[i - 1];
                }
                for (i = 1; i <= nparms; i = i + 1) {
                    if (inibl[i - 1] < 0) {
                        write(nout, format_9989), " INIBL ", inibl[i - 1], 0;
                        fatal = true;
                    }
                }
                write_loop wloop(cmn, nout, format_9983);
                wloop, "INIBL: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, inibl[i - 1];
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
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> ishfts[i - 1];
                }
                for (i = 1; i <= nparms; i = i + 1) {
                    if (ishfts[i - 1] < 0) {
                        write(nout, format_9989), " ISHFTS ", ishfts[i - 1], 0;
                        fatal = true;
                    }
                }
                write_loop wloop(cmn, nout, format_9983);
                wloop, "ISHFTS: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, ishfts[i - 1];
                }
            } else {
                for (i = 1; i <= nparms; i = i + 1) {
                    ishfts[i - 1] = 1;
                }
            }
            //
            //        Read the values for IACC22.
            //
            if (nep || dgg) {
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> iacc22[i - 1];
                }
                for (i = 1; i <= nparms; i = i + 1) {
                    if (iacc22[i - 1] < 0) {
                        write(nout, format_9989), " IACC22 ", iacc22[i - 1], 0;
                        fatal = true;
                    }
                }
                write_loop wloop(cmn, nout, format_9983);
                wloop, "IACC22: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, iacc22[i - 1];
                }
            } else {
                for (i = 1; i <= nparms; i = i + 1) {
                    iacc22[i - 1] = 1;
                }
            }
            //
            //        Read the values for NBCOL.
            //
            if (dgg) {
                getline(cin, str);
                iss.clear();
                iss.str(str);
                for (i = 1; i <= nparms; i = i + 1) {
                    iss >> nbcol[i - 1];
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
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NBCOL:";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nbcol[i - 1];
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
        sprintnum_short(buf0, eps);
        write(nout, format_9981), "underflow", buf0;
        eps = Rlamch("Overflow threshold");
        sprintnum_short(buf0, eps);
        write(nout, format_9981), "overflow ", buf0;
        eps = Rlamch("Epsilon");
        sprintnum_short(buf0, eps);
        write(nout, format_9981), "precision", buf0;
        //
        //     Read the threshold value for the test ratios.
        //
        getline(cin, str);
        iss.clear();
        iss.str(str);
        iss >> dtmp;
        thresh = dtmp;
        sprintnum_short(buf0, thresh);
        write(nout, "(/,' Routines pass computational tests if test ratio is ','less than',a,/)"), buf0;
        if (sep || svd || dgg) {
            //
            //        Read the flag that indicates whether to test LAPACK routines.
            //
            getline(cin, str);
            ss.clear();
            ss.str(str);
            ss >> line;
            if (Mlsame(line, "T"))
                tstchk = true;
            else
                tstchk = false;
            //
            //        Read the flag that indicates whether to test driver routines.
            //
            getline(cin, str);
            ss.clear();
            ss.str(str);
            ss >> line;
            if (Mlsame(line, "T"))
                tstdrv = true;
            else
                tstdrv = false;
        }
        //
        //     Read the flag that indicates whether to test the error exits.
        //
        getline(cin, str);
        ss.clear();
        ss.str(str);
        ss >> line;
        if (Mlsame(line, "T"))
            tsterr = true;
        else
            tsterr = false;

        //
        //     Read the code describing how to set the random number seed.
        //
        getline(cin, str);
        ss.clear();
        ss.str(str);
        ss >> newsd;
        //
        //     If NEWSD = 2, read another line with 4 integers for the seed.
        //
        if (newsd == 2) {
            getline(cin, str);
            iss.clear();
            iss.str(str);
            for (i = 1; i <= 4; i = i + 1) {
                iss >> ioldsd[i - 1];
            }
        }
        //
        for (i = 1; i <= 4; i = i + 1) {
            iseed[i - 1] = ioldsd[i - 1];
        }
        //
        if (fatal) {
            write(nout, "(/,' Execution not attempted due to input errors')");
            exit(1);
        }
        //
        //     Read the input lines indicating the test path and its parameters.
        //     The first three characters indicate the test path, and the number
        //     of test matrix types must be the first nonblank item in columns
        //     4-80.
        //
        //
        if (!(dgx || dxv)) {
            //
            string _str;
            getline(cin, str);
            c3[0] = str[0];
            c3[1] = str[1];
            c3[2] = str[2];
            iss.clear();
            iss.str(str);
            iss >> _str; // dummy read
            iss >> ntypes;
            //
            //     Skip the tests if NTYPES is <= 0.
            //
            if (!(dev || des || dvx || dsx || dgv || dgs) && ntypes <= 0) {
                write(nout, format_9990), c3;
                continue;
            } //
        } else {
            if (dxv) {
                strncpy(c3, "DXV", 3);
            }
            if (dgx) {
                strncpy(c3, "DGX", 3);
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
        if (Mlsamen(3, c3, "DHS") || Mlsamen(3, c3, "NEP")) {
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
                Rerrhs("Rhseqr", nout);
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
                Rchkhs(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], &a[nmax * nmax * 2], &a[nmax * nmax * 3], &a[nmax * nmax * 4], nmax, &a[nmax * nmax * 5], &a[nmax * nmax * 6], &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &d[nmax * 4], &d[nmax * 5], &a[nmax * nmax * 7], &a[nmax * nmax * 8], &a[nmax * nmax * 9], &a[nmax * nmax * 10], &a[nmax * nmax * 11], &d[nmax * 6], work, lwork, iwork, logwrk, result, info);
                if (info != 0) {
                    write(nout, format_9980), "Rchkhs", info;
                }
            }
            //
        } else if (Mlsamen(3, c3, "DST") || Mlsamen(3, c3, "SEP") || Mlsamen(3, c3, "SE2")) {
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
                Rerrst("DST", nout);
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
                        Rchkst2stg(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &d[nmax * 4], &d[nmax * 5], &d[nmax * 6], &d[nmax * 7], &d[nmax * 8], &d[nmax * 9], &d[nmax * 10], &a[nmax * nmax * 2], nmax, &a[nmax * nmax * 3], &a[nmax * nmax * 4], &d[nmax * 11], &a[nmax * nmax * 5], work, lwork, iwork, liwork, result, info);
                    } else {
                        Rchkst(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &d[nmax * 4], &d[nmax * 5], &d[nmax * 6], &d[nmax * 7], &d[nmax * 8], &d[nmax * 9], &d[nmax * 10], &a[nmax * nmax * 2], nmax, &a[nmax * nmax * 3], &a[nmax * nmax * 4], &d[nmax * 11], &a[nmax * nmax * 5], work, lwork, iwork, liwork, result, info);
                    }
                    if (info != 0) {
                        write(nout, format_9980), "Rchkst", info;
                    }
                }
                if (tstdrv) {
                    if (Mlsamen(3, c3, "SE2")) {
                        Rdrvst2stg(nn, nval, 18, dotype, iseed, thresh, nout, &a[0], nmax, &d[nmax * 2], &d[nmax * 3], &d[nmax * 4], &d[nmax * 5], &d[nmax * 7], &d[nmax * 8], &d[nmax * 9], &d[nmax * 10], &a[nmax * nmax], nmax, &a[nmax * nmax * 2], &d[nmax * 11], &a[nmax * nmax * 3], work, lwork, iwork, liwork, result, info);
                    } else {
                        Rdrvst(nn, nval, 18, dotype, iseed, thresh, nout, &a[0], nmax, &d[nmax * 2], &d[nmax * 3], &d[nmax * 4], &d[nmax * 5], &d[nmax * 7], &d[nmax * 8], &d[nmax * 9], &d[nmax * 10], &a[nmax * nmax], nmax, &a[nmax * nmax * 2], &d[nmax * 11], &a[nmax * nmax * 3], work, lwork, iwork, liwork, result, info);
                    }
                    if (info != 0) {
                        write(nout, format_9980), "Rdrvst", info;
                    }
                }
            }
            //
            break;
        } else if (Mlsamen(3, c3, "DSG")) {
            //
            //        ----------------------------------------------
            //        DSG:  Symmetric Generalized Eigenvalue Problem
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
                    //               CALL Rdrvsg( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
                    //     $                      NOUT, &a[(1-1)+(1-1)*lda], NMAX, &a[(1-1)+(2-1)*lda], NMAX,
                    //     $                      D( 1, 3 ), &a[(1-1)+(3-1)*lda], NMAX, &a[(1-1)+(4-1)*lda],
                    //     $                      &a[(1-1)+(5-1)*lda], &a[(1-1)+(6-1)*lda], &a[(1-1)+(7-1)*lda], WORK,
                    //     $                      LWORK, IWORK, LIWORK, RESULT, INFO )
                    Rdrvsg2stg(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], nmax, &d[nmax * 2], &d[nmax * 2], &a[nmax * nmax * 2], nmax, &a[nmax * nmax * 3], &a[nmax * nmax * 4], &a[nmax * nmax * 5], &a[nmax * nmax * 6], work, lwork, iwork, liwork, result, info);
                    if (info != 0) {
                        write(nout, format_9980), "Rdrvsg", info;
                    }
                }
            }
            //
        } else if (Mlsamen(3, c3, "DBD") || Mlsamen(3, c3, "SVD")) {
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
            xlaenv(1, 1);
            xlaenv(9, 25);
            //
            //        Test the error exits
            //
            if (tsterr && tstchk) {
                Rerrbd("DBD", nout);
            }
            if (tsterr && tstdrv) {
                Rerred("DBD", nout);
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
                    Rchkbd(nn, mval, nval, maxtyp, dotype, nrhs, iseed, thresh, &a[0], nmax, &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &a[nmax * nmax], nmax, &a[nmax * nmax * 2], &a[nmax * nmax * 3], &a[nmax * nmax * 4], nmax, &a[nmax * nmax * 5], nmax, &a[nmax * nmax * 6], &a[nmax * nmax * 7], work, lwork, iwork, nout, info);
                    if (info != 0) {
                        write(nout, format_9980), "Rchkbd", info;
                    }
                }
                if (tstdrv) {
                    Rdrvbd(nn, mval, nval, maxtyp, dotype, iseed, thresh, &a[0], nmax, &a[nmax * nmax], nmax, &a[nmax * nmax * 2], nmax, &a[nmax * nmax * 3], &a[nmax * nmax * 4], &a[nmax * nmax * 5], &d[0], &d[nmax], &d[nmax * 2], work, lwork, iwork, nout, info);
                }
            }
            //
        } else if (Mlsamen(3, c3, "DEV")) {
            //
            //        --------------------------------------------
            //        DEV:  Nonsymmetric Eigenvalue Problem Driver
            //              Rgeev (eigenvalues and eigenvectors)
            //        --------------------------------------------
            //
            maxtyp = 21;
            ntypes = min(maxtyp, ntypes);
            if (ntypes <= 0) {
                write(nout, format_9990), c3;
            } else {
                if (tsterr) {
                    Rerred(c3, nout);
                }
                Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
                Rdrvev(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &a[nmax * nmax * 2], nmax, &a[nmax * nmax * 3], nmax, &a[nmax * nmax * 4], nmax, result, work, lwork, iwork, info);
                if (info != 0) {
                    write(nout, format_9980), "Rgeev", info;
                }
            }
            write(nout, format_9973);
            continue;
            //
        } else if (Mlsamen(3, c3, "DES")) {
            //
            //        --------------------------------------------
            //        DES:  Nonsymmetric Eigenvalue Problem Driver
            //              Rgees (Schur form)
            //        --------------------------------------------
            //
            maxtyp = 21;
            ntypes = min(maxtyp, ntypes);
            if (ntypes <= 0) {
                write(nout, format_9990), c3;
            } else {
                if (tsterr) {
                    Rerred(c3, nout);
                }
                Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
                Rdrves(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], &a[nmax * nmax * 2], &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &a[nmax * nmax * 3], nmax, result, work, lwork, iwork, logwrk, info);
                if (info != 0) {
                    write(nout, format_9980), "Rgees", info;
                }
            }
            write(nout, format_9973);
            continue;
            //
        } else if (Mlsamen(3, c3, "DVX")) {
            //
            //        --------------------------------------------------------------
            //        DVX:  Nonsymmetric Eigenvalue Problem Expert Driver
            //              Rgeevx (eigenvalues, eigenvectors and condition numbers)
            //        --------------------------------------------------------------
            //
            maxtyp = 21;
            ntypes = min(maxtyp, ntypes);
            if (ntypes < 0) {
                write(nout, format_9990), c3;
            } else {
                if (tsterr) {
                    Rerred(c3, nout);
                }
                Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
                Rdrvvx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[0], nmax, &a[nmax * nmax], &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &a[nmax * nmax * 2], nmax, &a[nmax * nmax * 3], nmax, &a[nmax * nmax * 4], nmax, &d[nmax * 4], &d[nmax * 5], &d[nmax * 6], &d[nmax * 7], &d[nmax * 8], &d[nmax * 9], &d[nmax * 10], &d[nmax * 11], result, work, lwork, iwork, info);
                if (info != 0) {
                    write(nout, format_9980), "Rgeevx", info;
                }
            }
            write(nout, format_9973);
            continue;
            //
        } else if (Mlsamen(3, c3, "DSX")) {
            //
            //        ---------------------------------------------------
            //        DSX:  Nonsymmetric Eigenvalue Problem Expert Driver
            //              Rgeesx (Schur form and condition numbers)
            //        ---------------------------------------------------
            //
            maxtyp = 21;
            ntypes = min(maxtyp, ntypes);
            if (ntypes < 0) {
                write(nout, format_9990), c3;
            } else {
                if (tsterr) {
                    Rerred(c3, nout);
                }
                Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
                Rdrvsx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[0], nmax, &a[nmax * nmax], &a[nmax * nmax * 2], &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &d[nmax * 4], &d[nmax * 5], &a[nmax * nmax * 3], nmax, &a[nmax * nmax * 4], result, work, lwork, iwork, logwrk, info);
                if (info != 0) {
                    write(nout, format_9980), "Rgeesx", info;
                }
            }
            write(nout, format_9973);
            continue;
            //
        } else if (Mlsamen(3, c3, "DGG")) {
            //
            //        -------------------------------------------------
            //        DGG:  Generalized Nonsymmetric Eigenvalue Problem
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
                Rerrgg(c3, nout);
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
                    Rchkgg(nn, nval, maxtyp, dotype, iseed, thresh, tstdif, thrshn, nout, &a[0], nmax, &a[nmax * nmax], &a[nmax * nmax * 2], &a[nmax * nmax * 3], &a[nmax * nmax * 4], &a[nmax * nmax * 5], &a[nmax * nmax * 6], &a[nmax * nmax * 7], &a[nmax * nmax * 8], nmax, &a[nmax * nmax * 9], &a[nmax * nmax * 10], &a[nmax * nmax * 11], &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &d[nmax * 4], &d[nmax * 5], &a[(1 - 1) + (13 - 1) * lda], &a[(1 - 1) + (14 - 1) * lda], work, lwork, logwrk, result, info);
                    if (info != 0) {
                        write(nout, format_9980), "Rchkgg", info;
                    }
                }
            }
            //
        } else if (Mlsamen(3, c3, "DGS")) {
            //
            //        -------------------------------------------------
            //        DGS:  Generalized Nonsymmetric Eigenvalue Problem
            //              Rgges (Schur form)
            //        -------------------------------------------------
            //
            maxtyp = 26;
            ntypes = min(maxtyp, ntypes);
            if (ntypes <= 0) {
                write(nout, format_9990), c3;
            } else {
                if (tsterr) {
                    Rerrgg(c3, nout);
                }
                Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
                Rdrges(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], &a[nmax * nmax * 2], &a[nmax * nmax * 3], &a[nmax * nmax * 6], nmax, &a[nmax * nmax * 7], &d[0], &d[nmax], &d[nmax * 2], work, lwork, result, logwrk, info);
                if (info != 0) {
                    write(nout, format_9980), "Rdrges", info;
                }
                //
                //     Blocked version
                //
                xlaenv(16, 2);
                Rdrges3(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], &a[nmax * nmax * 2], &a[nmax * nmax * 3], &a[nmax * nmax * 6], nmax, &a[nmax * nmax * 7], &d[0], &d[nmax], &d[nmax * 2], work, lwork, result, logwrk, info);
                if (info != 0) {
                    write(nout, format_9980), "Rdrges3", info;
                }
            }
            write(nout, format_9973);
            continue;
            //
        } else if (dgx) {
            //
            //        -------------------------------------------------
            //        DGX:  Generalized Nonsymmetric Eigenvalue Problem
            //              Rggesx (Schur form and condition numbers)
            //        -------------------------------------------------
            //
            maxtyp = 5;
            ntypes = maxtyp;
            if (nn < 0) {
                write(nout, format_9990), c3;
            } else {
                if (tsterr) {
                    Rerrgg(c3, nout);
                }
                Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
                xlaenv(5, 2);
                Rdrgsx(nn, ncmax, thresh, nin, nout, &a[0], nmax, &a[nmax * nmax], &a[nmax * nmax * 2], &a[nmax * nmax * 3], &a[nmax * nmax * 4], &a[nmax * nmax * 5], &d[0], &d[nmax], &d[nmax * 2], &c[0], ncmax * ncmax, &a[nmax * nmax * 11], work, lwork, iwork, liwork, logwrk, info);
                if (info != 0) {
                    write(nout, format_9980), "Rdrgsx", info;
                }
            }
            write(nout, format_9973);
            continue;
            //
        } else if (Mlsamen(3, c3, "DGV")) {
            //
            //        -------------------------------------------------
            //        DGV:  Generalized Nonsymmetric Eigenvalue Problem
            //              Rggev (Eigenvalue/vector form)
            //        -------------------------------------------------
            //
            maxtyp = 26;
            ntypes = min(maxtyp, ntypes);
            if (ntypes <= 0) {
                write(nout, format_9990), c3;
            } else {
                if (tsterr) {
                    Rerrgg(c3, nout);
                }
                Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
                Rdrgev(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], &a[nmax * nmax * 2], &a[nmax * nmax * 3], &a[nmax * nmax * 6], nmax, &a[nmax * nmax * 7], &a[nmax * nmax * 8], nmax, &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &d[nmax * 4], &d[nmax * 5], work, lwork, result, info);
                if (info != 0) {
                    write(nout, format_9980), "Rdrgev", info;
                }
                //
                //     Blocked version
                //
                Rdrgev3(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], &a[nmax * nmax * 2], &a[nmax * nmax * 3], &a[nmax * nmax * 6], nmax, &a[nmax * nmax * 7], &a[nmax * nmax * 8], nmax, &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &d[nmax * 4], &d[nmax * 5], work, lwork, result, info);
                if (info != 0) {
                    write(nout, format_9980), "Rdrgev3", info;
                }
            }
            write(nout, format_9973);
            continue;
            //
        } else if (dxv) {
            //
            //        -------------------------------------------------
            //        DXV:  Generalized Nonsymmetric Eigenvalue Problem
            //              Rggevx (eigenvalue/vector with condition numbers)
            //        -------------------------------------------------
            //
            maxtyp = 2;
            ntypes = maxtyp;
            if (nn < 0) {
                write(nout, format_9990), c3;
            } else {
                if (tsterr) {
                    Rerrgg(c3, nout);
                }
                Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
                Rdrgvx(nn, thresh, nin, nout, &a[0], nmax, &a[nmax * nmax], &a[nmax * nmax * 2], &a[nmax * nmax * 3], &d[0], &d[nmax], &d[nmax * 2], &a[nmax * nmax * 4], &a[nmax * nmax * 5], iwork[1 - 1], iwork[2 - 1], &d[nmax * 3], &d[nmax * 4], &d[nmax * 5], &d[nmax * 6], &d[nmax * 7], &d[nmax * 8], work, lwork, &iwork[3 - 1], liwork - 2, result, logwrk, info);
                //
                if (info != 0) {
                    write(nout, format_9980), "Rdrgvx", info;
                }
            }
            write(nout, format_9973);
            continue;
            //
        } else if (Mlsamen(3, c3, "DSB")) {
            //
            //        ------------------------------
            //        DSB:  Symmetric Band Reduction
            //        ------------------------------
            //
            maxtyp = 15;
            ntypes = min(maxtyp, ntypes);
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            if (tsterr) {
                Rerrst("DSB", nout);
            }
            //         CALL Rchksb( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
            //     $                NOUT, &a[(1-1)+(1-1)*lda], NMAX, D( 1, 1 ), D( 1, 2 ),
            //     $                &a[(1-1)+(2-1)*lda], NMAX, WORK, LWORK, RESULT, INFO )
            Rchksb2stg(nn, nval, nk, kval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &d[0], &d[nmax], &d[nmax * 2], &d[nmax * 3], &d[nmax * 4], &a[nmax * nmax], nmax, work, lwork, result, info);
            if (info != 0) {
                write(nout, format_9980), "Rchksb", info;
            }
            //
        } else if (Mlsamen(3, c3, "DBB")) {
            //
            //        ------------------------------
            //        DBB:  General Band Reduction
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
                Rchkbb(nn, mval, nval, nk, kval, maxtyp, dotype, nrhs, iseed, thresh, nout, &a[0], nmax, &a[nmax * nmax], 2 * nmax, &d[0], &d[nmax], &a[nmax * nmax * 3], nmax, &a[nmax * nmax * 4], nmax, &a[nmax * nmax * 5], nmax, &a[nmax * nmax * 6], work, lwork, result, info);
                if (info != 0) {
                    write(nout, format_9980), "Rchkbb", info;
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
                Rerrgg("GLM", nout);
            }
            Rckglm(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[0], &a[nmax * nmax], &b[0], &b[nmax * nmax], x, work, &d[0], nin, nout, info);
            if (info != 0) {
                write(nout, format_9980), "Rckglm", info;
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
                Rerrgg("GQR", nout);
            }
            Rckgqr(nn, mval, nn, pval, nn, nval, ntypes, iseed, thresh, nmax, &a[0], &a[nmax * nmax], &a[nmax * nmax * 2], &a[nmax * nmax * 3], taua, &b[0], &b[nmax * nmax], &b[nmax * nmax * 2], &b[nmax * nmax * 3], &b[nmax * nmax * 4], taub, work, &d[0], nin, nout, info);
            if (info != 0) {
                write(nout, format_9980), "Rckgqr", info;
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
                Rerrgg("GSV", nout);
            }
            Rckgsv(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[0], &a[nmax * nmax], &b[0], &b[nmax * nmax], &a[nmax * nmax * 2], &b[nmax * nmax * 2], &a[nmax * nmax * 3], taua, taub, &b[nmax * nmax * 3], iwork, work, &d[0], nin, nout, info);
            if (info != 0) {
                write(nout, format_9980), "Rckgsv", info;
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
                Rerrgg("CSD", nout);
            }
            Rckcsd(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[0], &a[nmax * nmax], &a[nmax * nmax * 2], &a[nmax * nmax * 3], &a[nmax * nmax * 4], &a[nmax * nmax * 5], &a[nmax * nmax * 6], iwork, work, &d[0], nin, nout, info);
            if (info != 0) {
                write(nout, format_9980), "Rckcsd", info;
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
                Rerrgg("LSE", nout);
            }
            Rcklse(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[0], &a[nmax * nmax], &b[0], &b[nmax * nmax], x, work, &d[0], nin, nout, info);
            if (info != 0) {
                write(nout, format_9980), "Rcklse", info;
            }
            //
        } else {
            write(nout, star);
            write(nout, star);
            write(nout, format_9992), c3;
        }
        if (!(dgx || dxv)) {
            continue;
        }
    }
    write(nout, "(/,/,' End of tests')");
    s2 = time(NULL);
    write(nout, "(' Total time used = ',f12.2,' seconds',/)"), double(s2 - s1);
    //
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] work;
    //
    //     End of Rchkee
    //
}

int main(int argc, char const *argv[]) { Rchkee(); }
