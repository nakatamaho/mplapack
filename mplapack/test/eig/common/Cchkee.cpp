/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ZCHKEE.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_eig.h>
#include <memory>

void Cchkee(void) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    static INTEGER ioldsd[4] = {0, 0, 0, 1};
    static fem::str<10> intstr = "0123456789";
    const INTEGER nmax = 132;
    std::unique_ptr<REAL[]> s_storage;
    REAL *s = nullptr;
    const INTEGER need = 14;
    std::unique_ptr<COMPLEX[]> a_storage;
    COMPLEX *a = nullptr;
    std::unique_ptr<COMPLEX[]> b_storage;
    COMPLEX *b = nullptr;
    const INTEGER ncmax = 20;
    std::unique_ptr<COMPLEX[]> c_storage;
    COMPLEX *c = nullptr;
    const INTEGER lwork = nmax * (5 * nmax + 20);
    std::unique_ptr<REAL[]> rwork_storage;
    REAL *rwork = nullptr;
    std::unique_ptr<COMPLEX[]> work_storage;
    COMPLEX *work = nullptr;
    auto dc_storage = std::make_unique<COMPLEX[]>(std::max<INTEGER>(1, nmax * 6));
    COMPLEX *dc = dc_storage.get();
    REAL s1 = 0.0;
    bool fatal = false;
    const INTEGER nout = 6;
    const INTEGER nin = 5;
    fem::str<80> line;
    fem::str<3> path;
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
    REAL thresh = 0.0;
    REAL thresh_org = 0.0;
    bool tsterr = false;
    INTEGER mplapack_vers_major = 0;
    INTEGER mplapack_vers_minor = 0;
    INTEGER mplapack_vers_patch = 0;
    INTEGER lapack_vers_major = 0;
    INTEGER lapack_vers_minor = 0;
    INTEGER lapack_vers_patch = 0;
    INTEGER nn = 0;
    const INTEGER maxin = 20;
    INTEGER mval[maxin];
    INTEGER i = 0;
    fem::str<32> vname;
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
    fem::str<3> c3;
    INTEGER lenp = 0;
    INTEGER itmp = 0;
    INTEGER i1 = 0;
    const INTEGER maxt = 30;
    INTEGER ntypes = 0;
    fem::str<1> c1;
    INTEGER k = 0;
    INTEGER ic = 0;
    INTEGER maxtyp = 0;
    bool dotype[maxt];
    const INTEGER liwork = nmax * (nmax + 20);
    INTEGER iwork[liwork];
    bool logwrk[nmax];
    REAL result[500];
    INTEGER info = 0;
    auto dr_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * 12));
    REAL *dr = dr_storage.get();
    INTEGER nrhs = 0;
    bool tstdif = false;
    REAL thrshn = 0.0;
    COMPLEX x[5 * nmax];
    COMPLEX taua[nmax];
    COMPLEX taub[nmax];
    REAL alpha[nmax];
    REAL beta[nmax];
    REAL s2 = 0.0;
    INTEGER lda = nmax * nmax;
    INTEGER ldb = nmax * nmax;
    //
    static const char *format_9999 = "(/,' Execution not attempted due to input errors')";
    static const char *format_9997 = "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4)";
    static const char *format_9996 = "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NS =',i4,', MAXB =',i4,"
                                     "', IACC22 =',i4,', NBCOL =',i4)";
    static const char *format_9995 = "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4,', NRHS =',i4)";
    static const char *format_9994 = "(/,/,' End of tests')";
    static const char *format_9993 = "(' Total time used = ',f12.2,' seconds',/)";
    static const char *format_9992 = "(1x,a3,':  Unrecognized path name')";
    static const char *format_9991 = "(/,/,' *** Invalid integer value in column ',i2,' of input',' line:',/,"
                                     "a79)";
    static const char *format_9990 = "(/,/,1x,a3,' routines were not tested')";
    static const char *format_9989 = "(' Invalid input value: ',a,'=',i6,'; must be >=',i6)";
    static const char *format_9988 = "(' Invalid input value: ',a,'=',i6,'; must be <=',i6)";
    static const char *format_9987 = "(' Tests of the Nonsymmetric Eigenvalue Problem routines')";
    static const char *format_9986 = "(' Tests of the Hermitian Eigenvalue Problem routines')";
    static const char *format_9985 = "(' Tests of the Singular Value Decomposition routines')";
    static const char *format_9984 = "(/,' The following parameter values will be used:')";
    static const char *format_9983 = "(4x,a,10i6,/,10x,10i6)";
    static const char *format_9982 = "(/,' Routines pass computational tests if test ratio is ','less than',"
                                     "f8.2,/)";
    static const char *format_9981 = "(' Relative machine ',a,' is taken to be',d16.6)";
    static const char *format_9980 = "(' *** Error code from ',a,' = ',i4)";
    static const char *format_9979 = "(/,' Tests of the Nonsymmetric Eigenvalue Problem Driver',/,"
                                     "'    Cgeev (eigenvalues and eigevectors)')";
    static const char *format_9978 = "(/,' Tests of the Nonsymmetric Eigenvalue Problem Driver',/,"
                                     "'    Cgees (Schur form)')";
    static const char *format_9977 = "(/,' Tests of the Nonsymmetric Eigenvalue Problem Expert',' Driver',/,"
                                     "'    Cgeevx (eigenvalues, eigenvectors and',' condition numbers)')";
    static const char *format_9976 = "(/,' Tests of the Nonsymmetric Eigenvalue Problem Expert',' Driver',/,"
                                     "'    Cgeesx (Schur form and condition',' numbers)')";
    static const char *format_9975 = "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                                     "'Problem routines')";
    static const char *format_9974 = "(' Tests of Chbtrd',/,' (reduction of a Hermitian band ',"
                                     "'matrix to real tridiagonal form)')";
    static const char *format_9973 = "(/,1x,71('-'))";
    static const char *format_9972 = "(' Tests of the Multiple precision version of LAPACK ',i1,'.',i1,'.',i1,/, "
                                     "' Based on the original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, "
                                     "'The following parameter values will be used:')";
    static const char *format_9971 = "(/,' Tests of the Generalized Linear Regression Model ','routines')";
    static const char *format_9970 = "(/,' Tests of the Generalized QR and RQ routines')";
    static const char *format_9969 = "(/,' Tests of the Generalized Singular Value',' Decomposition routines')";
    static const char *format_9968 = "(/,' Tests of the Linear Least Squares routines')";
    static const char *format_9967 = "(' Tests of Cgbbrd',/,' (reduction of a general band ',"
                                     "'matrix to real bidiagonal form)')";
    static const char *format_9966 = "(/,/,1x,a3,':  NRHS =',i4)";
    static const char *format_9965 = "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                                     "'Problem Expert Driver Cggesx')";
    static const char *format_9964 = "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                                     "'Problem Driver Cgges')";
    static const char *format_9963 = "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                                     "'Problem Driver Cggev')";
    static const char *format_9962 = "(/,' Tests of the Generalized Nonsymmetric Eigenvalue ',"
                                     "'Problem Expert Driver Cggevx')";
    static const char *format_9961 = "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4,', INMIN=',i4,"
                                     "', INWIN =',i4,', INIBL =',i4,', ISHFTS =',i4,', IACC22 =',i4)";
    static const char *format_9960 = "(/,' Tests of the CS Decomposition routines')";
    //
    s_storage = std::make_unique<REAL[]>(max((INTEGER)1, (nmax * nmax)));
    s = s_storage.get();
    a_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, (nmax * nmax) * need));
    a = a_storage.get();
    b_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, (nmax * nmax) * 5));
    b = b_storage.get();
    c_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, (ncmax * ncmax) * (ncmax * ncmax)));
    c = c_storage.get();
    rwork_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
    rwork = rwork_storage.get();
    work_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lwork));
    work = work_storage.get();
    //
    const COMPLEX zero = COMPLEX(0.0);
    std::fill_n(a, nmax * nmax * need, zero);
    std::fill_n(b, nmax * nmax * 5, zero);
    std::fill_n(c, ncmax * ncmax * ncmax * ncmax, zero);
    std::fill_n(dc, nmax * 6, zero);
    s1 = dsecnd();
    fatal = false;
    nunit = nout;
//
// Return to here to read multiple sets of data
//
statement_10:
    //
    // Read the first line and set the 3-character test path
    //
    try {
        read(nin, "(a80)"), line;
    } catch (fem::read_end const &) {
        goto statement_380;
    }
    path = line(1, 3);
    nep = Mlsamen(3, path.elems, "NEP") || Mlsamen(3, path.elems, "ZHS");
    sep = Mlsamen(3, path.elems, "SEP") || Mlsamen(3, path.elems, "ZST") || Mlsamen(3, path.elems, "ZSG") || Mlsamen(3, path.elems, "SE2");
    svd = Mlsamen(3, path.elems, "SVD") || Mlsamen(3, path.elems, "ZBD");
    zev = Mlsamen(3, path.elems, "ZEV");
    zes = Mlsamen(3, path.elems, "ZES");
    zvx = Mlsamen(3, path.elems, "ZVX");
    zsx = Mlsamen(3, path.elems, "ZSX");
    zgg = Mlsamen(3, path.elems, "ZGG");
    zgs = Mlsamen(3, path.elems, "ZGS");
    zgx = Mlsamen(3, path.elems, "ZGX");
    zgv = Mlsamen(3, path.elems, "ZGV");
    zxv = Mlsamen(3, path.elems, "ZXV");
    zhb = Mlsamen(3, path.elems, "ZHB");
    zbb = Mlsamen(3, path.elems, "ZBB");
    glm = Mlsamen(3, path.elems, "GLM");
    gqr = Mlsamen(3, path.elems, "GQR") || Mlsamen(3, path.elems, "GRQ");
    gsv = Mlsamen(3, path.elems, "GSV");
    csd = Mlsamen(3, path.elems, "CSD");
    lse = Mlsamen(3, path.elems, "LSE");
    zbl = Mlsamen(3, path.elems, "ZBL");
    zbk = Mlsamen(3, path.elems, "ZBK");
    zgl = Mlsamen(3, path.elems, "ZGL");
    zgk = Mlsamen(3, path.elems, "ZGK");
    //
    // Report values of parameters.
    //
    if (path == "   ") {
        goto statement_10;
    } else if (nep) {
        write(nout, format_9987);
    } else if (sep) {
        write(nout, format_9986);
    } else if (svd) {
        write(nout, format_9985);
    } else if (zev) {
        write(nout, format_9979);
    } else if (zes) {
        write(nout, format_9978);
    } else if (zvx) {
        write(nout, format_9977);
    } else if (zsx) {
        write(nout, format_9976);
    } else if (zgg) {
        write(nout, format_9975);
    } else if (zgs) {
        write(nout, format_9964);
    } else if (zgx) {
        write(nout, format_9965);
    } else if (zgv) {
        write(nout, format_9963);
    } else if (zxv) {
        write(nout, format_9962);
    } else if (zhb) {
        write(nout, format_9974);
    } else if (zbb) {
        write(nout, format_9967);
    } else if (glm) {
        write(nout, format_9971);
    } else if (gqr) {
        write(nout, format_9970);
    } else if (gsv) {
        write(nout, format_9969);
    } else if (csd) {
        write(nout, format_9960);
    } else if (lse) {
        write(nout, format_9968);
    } else if (zbl) {
        //
        // Cgebal:  Balancing
        //
        Cchkbl(nin, nout);
        goto statement_380;
    } else if (zbk) {
        //
        // Cgebak:  Back transformation
        //
        Cchkbk(nin, nout);
        goto statement_380;
    } else if (zgl) {
        //
        // Cggbal:  Balancing
        //
        Cchkgl(nin, nout);
        goto statement_380;
    } else if (zgk) {
        //
        // Cggbak:  Back transformation
        //
        Cchkgk(nin, nout);
        goto statement_380;
    } else if (Mlsamen(3, path.elems, "ZEC")) {
        //
        // ZEC:  Eigencondition estimation
        //
        read(nin, star), thresh;
        Mxlaenv(1, 1);
        Mxlaenv(12, 1);
        tsterr = true;
        Cchkec(thresh, tsterr, nin, nout);
        goto statement_380;
    } else {
        write(nout, format_9992), path;
        goto statement_380;
    }
    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
    write(nout, format_9972), mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
    write(nout, format_9984);
    //
    // Read the number of values of M, P, and N.
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
    // Read the values of M
    //
    if (!(zgx || zxv)) {
        {
            read_loop rloop(cmn, nin, star);
            for (i = 1; i <= nn; i = i + 1) {
                rloop, mval[i - 1];
            }
        }
        if (svd) {
            vname = "    M ";
        } else {
            vname = "    N ";
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
    // Read the values of P
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
    // Read the values of N
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
    // Read the number of values of K, followed by the values of K
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
        // For the nonsymmetric QR driver routines, only one set of
        // parameters is allowed.
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
        Mxlaenv(1, nbval[1 - 1]);
        Mxlaenv(2, nbmin[1 - 1]);
        Mxlaenv(3, nxval[1 - 1]);
        Mxlaenv(12, max((INTEGER)11, inmin[1 - 1]));
        Mxlaenv(13, inwin[1 - 1]);
        Mxlaenv(14, inibl[1 - 1]);
        Mxlaenv(15, ishfts[1 - 1]);
        Mxlaenv(16, iacc22[1 - 1]);
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
        // For the nonsymmetric generalized driver routines, only one set of
        // parameters is allowed.
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
        Mxlaenv(1, nbval[1 - 1]);
        Mxlaenv(2, nbmin[1 - 1]);
        Mxlaenv(3, nxval[1 - 1]);
        Mxlaenv(4, nsval[1 - 1]);
        Mxlaenv(8, mxbval[1 - 1]);
        write(nout, format_9983), "NB:   ", nbval[1 - 1];
        write(nout, format_9983), "NBMIN:", nbmin[1 - 1];
        write(nout, format_9983), "NX:   ", nxval[1 - 1];
        write(nout, format_9983), "NS:   ", nsval[1 - 1];
        write(nout, format_9983), "MAXB: ", mxbval[1 - 1];
    } else if (!zhb && !glm && !gqr && !gsv && !csd && !lse) {
        //
        // For the other paths, the number of parameters can be varied
        // from the input file.  Read the number of parameter values.
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
        // Read the values of NB
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
        // Read the values of NBMIN
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
        // Read the values of NX
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
        // Read the values of NSHIFT (if ZGG) or NRHS (if SVD
        // or ZBB).
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
        // Read the values for MAXB.
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
        // Read the values for INMIN.
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
        // Read the values for INWIN.
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
        // Read the values for INIBL.
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
        // Read the values for ISHFTS.
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
        // Read the values for IACC22.
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
        // Read the values for NBCOL.
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
    // Calculate and print the machine dependent constants.
    //
    write(nout, star);
    eps = Rlamch("Underflow threshold");
    write(nout, format_9981), "underflow", eps;
    eps = Rlamch("Overflow threshold");
    write(nout, format_9981), "overflow ", eps;
    eps = Rlamch("Epsilon");
    write(nout, format_9981), "precision", eps;
    //
    // Read the threshold value for the test ratios.
    //
    read(nin, star), thresh;
    write(nout, format_9982), thresh;
    if (sep || svd || zgg) {
        //
        // Read the flag that indicates whether to test LAPACK routines.
        //
        read(nin, star), tstchk;
        //
        // Read the flag that indicates whether to test driver routines.
        //
        read(nin, star), tstdrv;
    }
    //
    // Read the flag that indicates whether to test the error exits.
    //
    read(nin, star), tsterr;
    //
    // Read the code describing how to set the random number seed.
    //
    read(nin, star), newsd;
    //
    // If NEWSD = 2, read another line with 4 integers for the seed.
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
        write(nout, format_9999);
        FEM_STOP(0);
    }
//
// Read the input lines indicating the test path and its parameters.
// The first three characters indicate the test path, and the number
// of test matrix types must be the first nonblank item in columns
// 4-80.
//
statement_190:
    //
    if (!(zgx || zxv)) {
    //
    statement_200:
        try {
            read(nin, "(a80)"), line;
        } catch (fem::read_end const &) {
            goto statement_380;
        }
        c3 = line(1, 3);
        lenp = fem::len(line);
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
        if (line(i, i) != " " && line(i, i) != ",") {
            i1 = i;
            c1 = line(i1, i1);
            //
            // Check that a valid integer was read
            //
            for (k = 1; k <= 10; k = k + 1) {
                if (c1 == intstr(k, k)) {
                    ic = k - 1;
                    goto statement_230;
                }
            }
            write(nout, format_9991), i, line;
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
        // Skip the tests if NTYPES is <= 0.
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
    // Reset the random number seed.
    //
    if (newsd == 0) {
        for (k = 1; k <= 4; k = k + 1) {
            iseed[k - 1] = ioldsd[k - 1];
        }
    }
    //
    if (Mlsamen(3, c3.elems, "ZHS") || Mlsamen(3, c3.elems, "NEP")) {
        //
        // -------------------------------------
        // NEP:  Nonsymmetric Eigenvalue Problem
        // -------------------------------------
        // Vary the parameters
        // NB    = block size
        // NBMIN = minimum block size
        // NX    = crossover point
        // NS    = number of shifts
        // MAXB  = minimum submatrix size
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        Mxlaenv(1, 1);
        if (tsterr) {
            Cerrhs("ZHSEQR", nout);
        }
        for (i = 1; i <= nparms; i = i + 1) {
            Mxlaenv(1, nbval[i - 1]);
            Mxlaenv(2, nbmin[i - 1]);
            Mxlaenv(3, nxval[i - 1]);
            Mxlaenv(12, max((INTEGER)11, inmin[i - 1]));
            Mxlaenv(13, inwin[i - 1]);
            Mxlaenv(14, inibl[i - 1]);
            Mxlaenv(15, ishfts[i - 1]);
            Mxlaenv(16, iacc22[i - 1]);
            //
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, format_9961), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], max((INTEGER)11, inmin[i - 1]), inwin[i - 1], inibl[i - 1], ishfts[i - 1], iacc22[i - 1];
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh_org = thresh;
            thresh = thresh * 5.0;
            printf("Warning! Threshold has been lifted to: ");
            printnum_short(thresh);
            printf(" for GMP\n");
#endif
            Cchkhs(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], &a[(7 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(8 - 1) * lda], &a[(9 - 1) * lda], &a[(10 - 1) * lda], &a[(11 - 1) * lda], &a[(12 - 1) * lda], &dc[(3 - 1) * nmax], work, lwork, rwork, iwork, logwrk, result, info);
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh = thresh_org;
#endif
            if (info != 0) {
                write(nout, format_9980), "Cchkhs", info;
            }
        }
        //
    } else if (Mlsamen(3, c3.elems, "ZST") || Mlsamen(3, c3.elems, "SEP") || Mlsamen(3, c3.elems, "SE2")) {
        //
        // ----------------------------------
        // SEP:  Symmetric Eigenvalue Problem
        // ----------------------------------
        // Vary the parameters
        // NB    = block size
        // NBMIN = minimum block size
        // NX    = crossover point
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        Mxlaenv(1, 1);
        Mxlaenv(9, 25);
        if (tsterr) {
            Cerrst("ZST", nout);
        }
        for (i = 1; i <= nparms; i = i + 1) {
            Mxlaenv(1, nbval[i - 1]);
            Mxlaenv(2, nbmin[i - 1]);
            Mxlaenv(3, nxval[i - 1]);
            //
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, format_9997), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1];
            if (tstchk) {
                if (Mlsamen(3, c3.elems, "SE2")) {
                    Cchkst2stg(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &dr[0], &dr[(2 - 1) * nmax], &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &dr[(5 - 1) * nmax], &dr[(6 - 1) * nmax], &dr[(7 - 1) * nmax], &dr[(8 - 1) * nmax], &dr[(9 - 1) * nmax], &dr[(10 - 1) * nmax], &dr[(11 - 1) * nmax], &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], &a[(5 - 1) * lda], &dc[0], &a[(6 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
                } else {
                    Cchkst(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &dr[0], &dr[(2 - 1) * nmax], &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &dr[(5 - 1) * nmax], &dr[(6 - 1) * nmax], &dr[(7 - 1) * nmax], &dr[(8 - 1) * nmax], &dr[(9 - 1) * nmax], &dr[(10 - 1) * nmax], &dr[(11 - 1) * nmax], &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], &a[(5 - 1) * lda], &dc[0], &a[(6 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
                }
                if (info != 0) {
                    write(nout, format_9980), "Cchkst", info;
                }
            }
            if (tstdrv) {
                if (Mlsamen(3, c3.elems, "SE2")) {
                    Cdrvst2stg(nn, nval, 18, dotype, iseed, thresh, nout, &a[0], nmax, &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &dr[(5 - 1) * nmax], &dr[(8 - 1) * nmax], &dr[(9 - 1) * nmax], &dr[(10 - 1) * nmax], &a[(2 - 1) * lda], nmax, &a[(3 - 1) * lda], &dc[0], &a[(4 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
                } else {
                    Cdrvst(nn, nval, 18, dotype, iseed, thresh, nout, &a[0], nmax, &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &dr[(5 - 1) * nmax], &dr[(8 - 1) * nmax], &dr[(9 - 1) * nmax], &dr[(10 - 1) * nmax], &a[(2 - 1) * lda], nmax, &a[(3 - 1) * lda], &dc[0], &a[(4 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
                }
                if (info != 0) {
                    write(nout, format_9980), "Cdrvst", info;
                }
            }
        }
        //
    } else if (Mlsamen(3, c3.elems, "ZSG")) {
        //
        // ----------------------------------------------
        // ZSG:  Hermitian Generalized Eigenvalue Problem
        // ----------------------------------------------
        // Vary the parameters
        // NB    = block size
        // NBMIN = minimum block size
        // NX    = crossover point
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        Mxlaenv(9, 25);
        for (i = 1; i <= nparms; i = i + 1) {
            Mxlaenv(1, nbval[i - 1]);
            Mxlaenv(2, nbmin[i - 1]);
            Mxlaenv(3, nxval[i - 1]);
            //
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, format_9997), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1];
            if (tstchk) {
                // CALL Cdrvsg( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
                // $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
                // $                      DR( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
                // $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
                // $                      LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT,
                // $                      INFO )
                Cdrvsg2stg(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], nmax, &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], &a[(5 - 1) * lda], &a[(6 - 1) * lda], &a[(7 - 1) * lda], work, lwork, rwork, lwork, iwork, liwork, result, info);
                if (info != 0) {
                    write(nout, format_9980), "Cdrvsg", info;
                }
            }
        }
        //
    } else if (Mlsamen(3, c3.elems, "ZBD") || Mlsamen(3, c3.elems, "SVD")) {
        //
        // ----------------------------------
        // SVD:  Singular Value Decomposition
        // ----------------------------------
        // Vary the parameters
        // NB    = block size
        // NBMIN = minimum block size
        // NX    = crossover point
        // NRHS  = number of right hand sides
        //
        maxtyp = 16;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        Mxlaenv(9, 25);
        //
        // Test the error exits
        //
        Mxlaenv(1, 1);
        if (tsterr && tstchk) {
            Cerrbd("ZBD", nout);
        }
        if (tsterr && tstdrv) {
            Cerred("ZBD", nout);
        }
        //
        for (i = 1; i <= nparms; i = i + 1) {
            nrhs = nsval[i - 1];
            Mxlaenv(1, nbval[i - 1]);
            Mxlaenv(2, nbmin[i - 1]);
            Mxlaenv(3, nxval[i - 1]);
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, format_9995), c3, nbval[i - 1], nbmin[i - 1], nxval[i - 1], nrhs;
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh_org = thresh;
            thresh = thresh * 4.0;
            printf("Warning! Threshold has been lifted to: ");
            printnum_short(thresh);
            printf(" for GMP\n");
#endif
            if (tstchk) {
                Cchkbd(nn, mval, nval, maxtyp, dotype, nrhs, iseed, thresh, &a[0], nmax, &dr[0], &dr[(2 - 1) * nmax], &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &a[(2 - 1) * lda], nmax, &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], nmax, &a[(7 - 1) * lda], &a[(8 - 1) * lda], work, lwork, rwork, nout, info);
                if (info != 0) {
                    write(nout, format_9980), "Cchkbd", info;
                }
            }
            if (tstdrv) {
                Cdrvbd(nn, mval, nval, maxtyp, dotype, iseed, thresh, &a[0], nmax, &a[(2 - 1) * lda], nmax, &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], &a[(5 - 1) * lda], &a[(6 - 1) * lda], &dr[0], &dr[(2 - 1) * nmax], &dr[(3 - 1) * nmax], work, lwork, rwork, iwork, nout, info);
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh = thresh_org;
#endif
        }
        //
    } else if (Mlsamen(3, c3.elems, "ZEV")) {
        //
        // --------------------------------------------
        // ZEV:  Nonsymmetric Eigenvalue Problem Driver
        // Cgeev (eigenvalues and eigenvectors)
        // --------------------------------------------
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        if (ntypes <= 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerred(c3, nout);
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh_org = thresh;
            thresh = thresh * 5.0;
            printf("Warning! Threshold has been lifted to: ");
            printnum_short(thresh);
            printf(" for GMP\n");
#endif
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrvev(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], nmax, &a[(5 - 1) * lda], nmax, result, work, lwork, rwork, iwork, info);
            if (info != 0) {
                write(nout, format_9980), "Cgeev", info;
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh = thresh_org;
#endif
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3.elems, "ZES")) {
        //
        // --------------------------------------------
        // ZES:  Nonsymmetric Eigenvalue Problem Driver
        // Cgees (Schur form)
        // --------------------------------------------
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        if (ntypes <= 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerred(c3, nout);
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh_org = thresh;
            thresh = thresh * 5.0;
            printf("Warning! Threshold has been lifted to: ");
            printnum_short(thresh);
            printf(" for GMP\n");
#endif
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrves(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(4 - 1) * lda], nmax, result, work, lwork, rwork, iwork, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Cgees", info;
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh = thresh_org;
#endif
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3.elems, "ZVX")) {
        //
        // --------------------------------------------------------------
        // ZVX:  Nonsymmetric Eigenvalue Problem Expert Driver
        // Cgeevx (eigenvalues, eigenvectors and condition numbers)
        // --------------------------------------------------------------
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        if (ntypes < 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerred(c3, nout);
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh_org = thresh;
            thresh = thresh * 4.0;
            printf("Warning! Threshold has been lifted to: ");
            printnum_short(thresh);
            printf(" for GMP\n");
#endif
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrvvx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[0], nmax, &a[(2 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(3 - 1) * lda], nmax, &a[(4 - 1) * lda], nmax, &a[(5 - 1) * lda], nmax, &dr[0], &dr[(2 - 1) * nmax], &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &dr[(5 - 1) * nmax], &dr[(6 - 1) * nmax], &dr[(7 - 1) * nmax], &dr[(8 - 1) * nmax], result, work, lwork, rwork, info);
            if (info != 0) {
                write(nout, format_9980), "Cgeevx", info;
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh = thresh_org;
#endif
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3.elems, "ZSX")) {
        //
        // ---------------------------------------------------
        // ZSX:  Nonsymmetric Eigenvalue Problem Expert Driver
        // Cgeesx (Schur form and condition numbers)
        // ---------------------------------------------------
        //
        maxtyp = 21;
        ntypes = min(maxtyp, ntypes);
        if (ntypes < 0) {
            write(nout, format_9990), c3;
        } else {
            if (tsterr) {
                Cerred(c3, nout);
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh_org = thresh;
            thresh = thresh * 5.0;
            printf("Warning! Threshold has been lifted to: ");
            printnum_short(thresh);
            printf(" for GMP\n");
#endif
            Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
            Cdrvsx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &dc[(3 - 1) * nmax], &a[(4 - 1) * lda], nmax, &a[(5 - 1) * lda], result, work, lwork, rwork, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Cgeesx", info;
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh = thresh_org;
#endif
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3.elems, "ZGG")) {
        //
        // -------------------------------------------------
        // ZGG:  Generalized Nonsymmetric Eigenvalue Problem
        // -------------------------------------------------
        // Vary the parameters
        // NB    = block size
        // NBMIN = minimum block size
        // NS    = number of shifts
        // MAXB  = minimum submatrix size
        // IACC22: structured matrix multiply
        // NBCOL = minimum column dimension for blocks
        //
        maxtyp = 26;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        Mxlaenv(1, 1);
        if (tstchk && tsterr) {
            Cerrgg(c3, nout);
        }
        for (i = 1; i <= nparms; i = i + 1) {
            Mxlaenv(1, nbval[i - 1]);
            Mxlaenv(2, nbmin[i - 1]);
            Mxlaenv(4, nsval[i - 1]);
            Mxlaenv(8, mxbval[i - 1]);
            Mxlaenv(16, iacc22[i - 1]);
            Mxlaenv(5, nbcol[i - 1]);
            //
            if (newsd == 0) {
                for (k = 1; k <= 4; k = k + 1) {
                    iseed[k - 1] = ioldsd[k - 1];
                }
            }
            write(nout, format_9996), c3, nbval[i - 1], nbmin[i - 1], nsval[i - 1], mxbval[i - 1], iacc22[i - 1], nbcol[i - 1];
            tstdif = false;
            thrshn = 10.0;
            if (tstchk) {
                Cchkgg(nn, nval, maxtyp, dotype, iseed, thresh, tstdif, thrshn, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], &a[(6 - 1) * lda], &a[(7 - 1) * lda], &a[(8 - 1) * lda], &a[(9 - 1) * lda], nmax, &a[(10 - 1) * lda], &a[(11 - 1) * lda], &a[(12 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &dc[(3 - 1) * nmax], &dc[(4 - 1) * nmax], &a[(13 - 1) * lda], &a[(14 - 1) * lda], work, lwork, rwork, logwrk, result, info);
                if (info != 0) {
                    write(nout, format_9980), "Cchkgg", info;
                }
            }
        }
        //
    } else if (Mlsamen(3, c3.elems, "ZGS")) {
        //
        // -------------------------------------------------
        // ZGS:  Generalized Nonsymmetric Eigenvalue Problem
        // Cgges (Schur form)
        // -------------------------------------------------
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
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh = thresh * 5.0;
            printf("Warning! Threshold has been lifted to: ");
            printnum_short(thresh);
            printf(" for GMP\n");
#endif
            Cdrges(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(7 - 1) * lda], nmax, &a[(8 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], work, lwork, rwork, result, logwrk, info);
            //
            if (info != 0) {
                write(nout, format_9980), "Cdrges", info;
            }
            //
            // Blocked version
            //
            Cdrges3(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(7 - 1) * lda], nmax, &a[(8 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], work, lwork, rwork, result, logwrk, info);
            //
            if (info != 0) {
                write(nout, format_9980), "Cdrges3", info;
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh = thresh_org;
#endif
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (zgx) {
        //
        // -------------------------------------------------
        // ZGX  Generalized Nonsymmetric Eigenvalue Problem
        // Cggesx (Schur form and condition numbers)
        // -------------------------------------------------
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
            Mxlaenv(5, 2);
            Cdrgsx(nn, ncmax, thresh, nin, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], &a[(6 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], c, ncmax * ncmax, s, work, lwork, rwork, iwork, liwork, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Cdrgsx", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3.elems, "ZGV")) {
        //
        // -------------------------------------------------
        // ZGV:  Generalized Nonsymmetric Eigenvalue Problem
        // Cggev (Eigenvalue/vector form)
        // -------------------------------------------------
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
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh_org = thresh;
            thresh = thresh * 4.0;
            printf("Warning! Threshold has been lifted to: ");
            printnum_short(thresh);
            printf(" for GMP\n");
#endif
            Cdrgev(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(7 - 1) * lda], nmax, &a[(8 - 1) * lda], &a[(9 - 1) * lda], nmax, &dc[0], &dc[(2 - 1) * nmax], &dc[(3 - 1) * nmax], &dc[(4 - 1) * nmax], work, lwork, rwork, result, info);
            if (info != 0) {
                write(nout, format_9980), "Cdrgev", info;
            }
            //
            // Blocked version
            //
            Mxlaenv(16, 2);
            Cdrgev3(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(7 - 1) * lda], nmax, &a[(8 - 1) * lda], &a[(9 - 1) * lda], nmax, &dc[0], &dc[(2 - 1) * nmax], &dc[(3 - 1) * nmax], &dc[(4 - 1) * nmax], work, lwork, rwork, result, info);
            if (info != 0) {
                write(nout, format_9980), "Cdrgev3", info;
            }
#if defined ___MPLAPACK_BUILD_WITH_GMP___
            thresh = thresh_org;
#endif
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (zxv) {
        //
        // -------------------------------------------------
        // ZXV:  Generalized Nonsymmetric Eigenvalue Problem
        // Cggevx (eigenvalue/vector with condition numbers)
        // -------------------------------------------------
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
            Cdrgvx(nn, thresh, nin, nout, &a[0], nmax, &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &dc[0], &dc[(2 - 1) * nmax], &a[(5 - 1) * lda], &a[(6 - 1) * lda], iwork[1 - 1], iwork[2 - 1], &dr[0], &dr[(2 - 1) * nmax], &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &dr[(5 - 1) * nmax], &dr[(6 - 1) * nmax], work, lwork, rwork, &iwork[3 - 1], liwork - 2, result, logwrk, info);
            //
            if (info != 0) {
                write(nout, format_9980), "Cdrgvx", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
        //
    } else if (Mlsamen(3, c3.elems, "ZHB")) {
        //
        // ------------------------------
        // ZHB:  Hermitian Band Reduction
        // ------------------------------
        //
        maxtyp = 15;
        ntypes = min(maxtyp, ntypes);
        Alareq(c3, ntypes, dotype, maxtyp, nin, nout);
        if (tsterr) {
            Cerrst("ZHB", nout);
        }
        // CALL Cchkhb( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
        // $                NOUT, A( 1, 1 ), NMAX, DR( 1, 1 ), DR( 1, 2 ),
        // $                A( 1, 2 ), NMAX, WORK, LWORK, RWORK, RESULT,
        // $                INFO )
        Cchkhb2stg(nn, nval, nk, kval, maxtyp, dotype, iseed, thresh, nout, &a[0], nmax, &dr[0], &dr[(2 - 1) * nmax], &dr[(3 - 1) * nmax], &dr[(4 - 1) * nmax], &dr[(5 - 1) * nmax], &a[(2 - 1) * lda], nmax, work, lwork, rwork, result, info);
        if (info != 0) {
            write(nout, format_9980), "Cchkhb", info;
        }
        //
    } else if (Mlsamen(3, c3.elems, "ZBB")) {
        //
        // ------------------------------
        // ZBB:  General Band Reduction
        // ------------------------------
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
            write(nout, format_9966), c3, nrhs;
            Cchkbb(nn, mval, nval, nk, kval, maxtyp, dotype, nrhs, iseed, thresh, nout, &a[0], nmax, &a[(2 - 1) * lda], 2 * nmax, &dr[0], &dr[(2 - 1) * nmax], &a[(4 - 1) * lda], nmax, &a[(5 - 1) * lda], nmax, &a[(6 - 1) * lda], nmax, &a[(7 - 1) * lda], work, lwork, rwork, result, info);
            if (info != 0) {
                write(nout, format_9980), "Cchkbb", info;
            }
        }
        //
    } else if (Mlsamen(3, c3.elems, "GLM")) {
        //
        // -----------------------------------------
        // GLM:  Generalized Linear Regression Model
        // -----------------------------------------
        //
        Mxlaenv(1, 1);
        if (tsterr) {
            Cerrgg("GLM", nout);
        }
        Cckglm(nn, nval, mval, pval, ntypes, iseed, thresh, nmax, &a[0], &a[(2 - 1) * lda], &b[0], &b[(2 - 1) * ldb], x, work, &dr[0], nin, nout, info);
        if (info != 0) {
            write(nout, format_9980), "Cckglm", info;
        }
        //
    } else if (Mlsamen(3, c3.elems, "GQR")) {
        //
        // ------------------------------------------
        // GQR:  Generalized QR and RQ factorizations
        // ------------------------------------------
        //
        Mxlaenv(1, 1);
        if (tsterr) {
            Cerrgg("GQR", nout);
        }
        Cckgqr(nn, mval, nn, pval, nn, nval, ntypes, iseed, thresh, nmax, &a[0], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], taua, &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], &b[(5 - 1) * ldb], taub, work, &dr[0], nin, nout, info);
        if (info != 0) {
            write(nout, format_9980), "Cckgqr", info;
        }
        //
    } else if (Mlsamen(3, c3.elems, "GSV")) {
        //
        // ----------------------------------------------
        // GSV:  Generalized Singular Value Decomposition
        // ----------------------------------------------
        //
        Mxlaenv(1, 1);
        if (tsterr) {
            Cerrgg("GSV", nout);
        }
        Cckgsv(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[0], &a[(2 - 1) * lda], &b[0], &b[(2 - 1) * ldb], &a[(3 - 1) * lda], &b[(3 - 1) * ldb], &a[(4 - 1) * lda], alpha, beta, &b[(4 - 1) * ldb], iwork, work, &dr[0], nin, nout, info);
        if (info != 0) {
            write(nout, format_9980), "Cckgsv", info;
        }
        //
    } else if (Mlsamen(3, c3.elems, "CSD")) {
        //
        // ----------------------------------------------
        // CSD:  CS Decomposition
        // ----------------------------------------------
        //
        Mxlaenv(1, 1);
        if (tsterr) {
            Cerrgg("CSD", nout);
        }
        Cckcsd(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[0], &a[(2 - 1) * lda], &a[(3 - 1) * lda], &a[(4 - 1) * lda], &a[(5 - 1) * lda], &a[(6 - 1) * lda], rwork, iwork, work, &dr[0], nin, nout, info);
        if (info != 0) {
            write(nout, format_9980), "Cckcsd", info;
        }
        //
    } else if (Mlsamen(3, c3.elems, "LSE")) {
        //
        // --------------------------------------
        // LSE:  Constrained Linear Least Squares
        // --------------------------------------
        //
        Mxlaenv(1, 1);
        if (tsterr) {
            Cerrgg("LSE", nout);
        }
        Ccklse(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[0], &a[(2 - 1) * lda], &b[0], &b[(2 - 1) * ldb], x, work, &dr[0], nin, nout, info);
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
statement_380:
    write(nout, format_9994);
    s2 = dsecnd();
    write(nout, format_9993), cast2double(s2 - s1);
    //
    // End of Cchkee
    //
}

int main(int argc, char const *argv[]) { Cchkee(); }
