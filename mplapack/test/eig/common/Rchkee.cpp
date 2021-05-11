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

void program_Rchkee(INTEGER argc, char const *argv[]) {
    common cmn(argc, argv);
    FEM_CMN_SVE(program_Rchkee);
    common_read read(cmn);
    common_write write(cmn);
    char &intstr = sve.intstr;
    INTEGER *ioldsd(sve.ioldsd, [4]);
    if (is_called_first_time) {
        intstr = "0123456789";
        {
            static const INTEGER values[] = {0, 0, 0, 1};
            data_of_type<int>(FEM_VALUES_AND_SIZE), ioldsd;
        }
    }
    INTEGER allocatestatus = 0;
    REAL a = 0.0;
    REAL b = 0.0;
    REAL c = 0.0;
    const INTEGER nmax = 132;
    REAL d[nmax * 12];
    REAL s1 = 0.0;
    bool fatal = false;
    const INTEGER nout = 6;
    const INTEGER nin = 5;
    char line[80];
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
    char c3[3];
    INTEGER lenp = 0;
    INTEGER itmp = 0;
    INTEGER i1 = 0;
    const INTEGER maxt = 30;
    INTEGER ntypes = 0;
    char c1;
    INTEGER k = 0;
    INTEGER ic = 0;
    INTEGER maxtyp = 0;
    bool dotype[maxt];
    INTEGER &a[(1 - 1) + (1 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (2 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (3 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (4 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (5 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (6 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (7 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (8 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (9 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (10 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (11 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (12 - 1) * lda] = 0;
    REAL work = 0.0;
    const INTEGER lwork = nmax * (5 * nmax + 5) + 1;
    const INTEGER liwork = nmax * (5 * nmax + 20);
    INTEGER iwork[liwork];
    bool logwrk[nmax];
    REAL result[500];
    INTEGER info = 0;
    INTEGER nrhs = 0;
    bool tstdif = false;
    REAL thrshn = 0.0;
    INTEGER &a[(1 - 1) + (13 - 1) * lda] = 0;
    INTEGER &a[(1 - 1) + (14 - 1) * lda] = 0;
    const INTEGER ncmax = 20;
    INTEGER &c[(1 - 1) + (1 - 1) * ldc] = 0;
    INTEGER &b[(1 - 1) + (1 - 1) * ldb] = 0;
    INTEGER &b[(1 - 1) + (2 - 1) * ldb] = 0;
    REAL x[5 * nmax];
    REAL taua[nmax];
    INTEGER &b[(1 - 1) + (3 - 1) * ldb] = 0;
    INTEGER &b[(1 - 1) + (4 - 1) * ldb] = 0;
    INTEGER &b[(1 - 1) + (5 - 1) * ldb] = 0;
    REAL taub[nmax];
    REAL s2 = 0.0;
    static const char *format_9973 = "(/,1x,71('-'))";
    static const char *format_9980 = "(' *** Error code from ',a,' = ',i4)";
    static const char *format_9981 = "(' Relative machine ',a,' is taken to be',d16.6)";
    static const char *format_9983 = "(4x,a,10i6,/,10x,10i6)";
    static const char *format_9988 = "(' Invalid input value: ',a,'=',i6,'; must be <=',i6)";
    static const char *format_9989 = "(' Invalid input value: ',a,'=',i6,'; must be >=',i6)";
    static const char *format_9990 = "(/,/,1x,a3,' routines were not tested')";
    static const char *format_9992 = "(1x,a3,':  Unrecognized path name')";
    static const char *format_9997 = "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4)";
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. Allocatable Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Arrays in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Allocate memory dynamically ..
    //
    FEM_THROW_UNHANDLED("executable allocate: allocate(a(nmax*nmax,need),stat=allocatestatus)");
    if (allocatestatus != 0) {
        FEM_STOP("*** Not enough memory ***");
    }
    FEM_THROW_UNHANDLED("executable allocate: allocate(b(nmax*nmax,5),stat=allocatestatus)");
    if (allocatestatus != 0) {
        FEM_STOP("*** Not enough memory ***");
    }
    FEM_THROW_UNHANDLED("executable allocate: allocate(c(ncmax*ncmax,ncmax*ncmax),stat=allocatesta"
                        "tus)");
    if (allocatestatus != 0) {
        FEM_STOP("*** Not enough memory ***");
    }
    FEM_THROW_UNHANDLED("executable allocate: allocate(work(lwork),stat=allocatestatus)");
    if (allocatestatus != 0) {
        FEM_STOP("*** Not enough memory ***");
    }
    //     ..
    //     .. Executable Statements ..
    //
    a = 0.0f;
    b = 0.0f;
    c = 0.0f;
    d = 0.0f;
    s1 = dsecnd[-1];
    fatal = false;
    cmn.nunit = nout;
//
//     Return to here to read multiple sets of data
//
statement_10:
    //
    //     Read the first line and set the 3-character test path
    //
    try {
        read(nin, "(a80)"), line;
    } catch (read_end const) {
        goto statement_380;
    }
    path = line[(3 - 1) * ldline];
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
        goto statement_10;
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
        goto statement_10;
    } else if (dbk) {
        //
        //        Rgebak:  Back transformation
        //
        Rchkbk(nin, nout);
        goto statement_10;
    } else if (dgl) {
        //
        //        Rggbal:  Balancing
        //
        Rchkgl(nin, nout);
        goto statement_10;
    } else if (dgk) {
        //
        //        Rggbak:  Back transformation
        //
        Rchkgk(nin, nout);
        goto statement_10;
    } else if (Mlsamen(3, path, "DEC")) {
        //
        //        DEC:  Eigencondition estimation
        //
        read(nin, star), thresh;
        xlaenv(1, 1);
        xlaenv(12, 11);
        xlaenv(13, 2);
        xlaenv(14, 0);
        xlaenv(15, 2);
        xlaenv(16, 2);
        tsterr = true;
        Rchkec(thresh, tsterr, nin, nout);
        goto statement_10;
    } else {
        write(nout, format_9992), path;
        goto statement_10;
    }
    ilaver(vers_major, vers_minor, vers_patch);
    write(nout, "(/,' LAPACK VERSION ',i1,'.',i1,'.',i1)"), vers_major, vers_minor, vers_patch;
    write(nout, "(/,' The following parameter values will be used:')");
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
    if (!(dgx || dxv)) {
        {
            read_loop rloop(cmn, nin, star);
            for (i = 1; i <= nn; i = i + 1) {
                rloop, mval(i);
            }
        }
        if (svd) {
            vname = "    M ";
        } else {
            vname = "    N ";
        }
        for (i = 1; i <= nn; i = i + 1) {
            if (mval[i - 1] < 0) {
                write(nout, format_9989), vname, mval(i), 0;
                fatal = true;
            } else if (mval[i - 1] > nmax) {
                write(nout, format_9988), vname, mval(i), nmax;
                fatal = true;
            }
        }
        {
            write_loop wloop(cmn, nout, format_9983);
            wloop, "M:    ";
            for (i = 1; i <= nn; i = i + 1) {
                wloop, mval(i);
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
                rloop, pval(i);
            }
        }
        for (i = 1; i <= nn; i = i + 1) {
            if (pval[i - 1] < 0) {
                write(nout, format_9989), " P  ", pval(i), 0;
                fatal = true;
            } else if (pval[i - 1] > nmax) {
                write(nout, format_9988), " P  ", pval(i), nmax;
                fatal = true;
            }
        }
        {
            write_loop wloop(cmn, nout, format_9983);
            wloop, "P:    ";
            for (i = 1; i <= nn; i = i + 1) {
                wloop, pval(i);
            }
        }
    }
    //
    //     Read the values of N
    //
    if (svd || dbb || glm || gqr || gsv || csd || lse) {
        {
            read_loop rloop(cmn, nin, star);
            for (i = 1; i <= nn; i = i + 1) {
                rloop, nval(i);
            }
        }
        for (i = 1; i <= nn; i = i + 1) {
            if (nval[i - 1] < 0) {
                write(nout, format_9989), "    N ", nval(i), 0;
                fatal = true;
            } else if (nval[i - 1] > nmax) {
                write(nout, format_9988), "    N ", nval(i), nmax;
                fatal = true;
            }
        }
    } else {
        for (i = 1; i <= nn; i = i + 1) {
            nval[i - 1] = mval[i - 1];
        }
    }
    if (!(dgx || dxv)) {
        {
            write_loop wloop(cmn, nout, format_9983);
            wloop, "N:    ";
            for (i = 1; i <= nn; i = i + 1) {
                wloop, nval(i);
            }
        }
    } else {
        write(nout, format_9983), "N:    ", nn;
    }
    //
    //     Read the number of values of K, followed by the values of K
    //
    if (dsb || dbb) {
        read(nin, star), nk;
        {
            read_loop rloop(cmn, nin, star);
            for (i = 1; i <= nk; i = i + 1) {
                rloop, kval(i);
            }
        }
        for (i = 1; i <= nk; i = i + 1) {
            if (kval[i - 1] < 0) {
                write(nout, format_9989), "    K ", kval(i), 0;
                fatal = true;
            } else if (kval[i - 1] > nmax) {
                write(nout, format_9988), "    K ", kval(i), nmax;
                fatal = true;
            }
        }
        {
            write_loop wloop(cmn, nout, format_9983);
            wloop, "K:    ";
            for (i = 1; i <= nk; i = i + 1) {
                wloop, kval(i);
            }
        }
    }
    //
    if (dev || des || dvx || dsx) {
        //
        //        For the nonsymmetric QR driver routines, only one set of
        //        parameters is allowed.
        //
        read(nin, star), nbval(1), nbmin(1), nxval(1), inmin(1), inwin(1), inibl(1), ishfts(1), iacc22(1);
        if (nbval[1 - 1] < 1) {
            write(nout, format_9989), "   NB ", nbval(1), 1;
            fatal = true;
        } else if (nbmin[1 - 1] < 1) {
            write(nout, format_9989), "NBMIN ", nbmin(1), 1;
            fatal = true;
        } else if (nxval[1 - 1] < 1) {
            write(nout, format_9989), "   NX ", nxval(1), 1;
            fatal = true;
        } else if (inmin[1 - 1] < 1) {
            write(nout, format_9989), "   INMIN ", inmin(1), 1;
            fatal = true;
        } else if (inwin[1 - 1] < 1) {
            write(nout, format_9989), "   INWIN ", inwin(1), 1;
            fatal = true;
        } else if (inibl[1 - 1] < 1) {
            write(nout, format_9989), "   INIBL ", inibl(1), 1;
            fatal = true;
        } else if (ishfts[1 - 1] < 1) {
            write(nout, format_9989), "   ISHFTS ", ishfts(1), 1;
            fatal = true;
        } else if (iacc22[1 - 1] < 0) {
            write(nout, format_9989), "   IACC22 ", iacc22(1), 0;
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
        write(nout, format_9983), "NB:   ", nbval(1);
        write(nout, format_9983), "NBMIN:", nbmin(1);
        write(nout, format_9983), "NX:   ", nxval(1);
        write(nout, format_9983), "INMIN:   ", inmin(1);
        write(nout, format_9983), "INWIN: ", inwin(1);
        write(nout, format_9983), "INIBL: ", inibl(1);
        write(nout, format_9983), "ISHFTS: ", ishfts(1);
        write(nout, format_9983), "IACC22: ", iacc22(1);
        //
    } else if (dgs || dgx || dgv || dxv) {
        //
        //        For the nonsymmetric generalized driver routines, only one set
        //        of parameters is allowed.
        //
        read(nin, star), nbval(1), nbmin(1), nxval(1), nsval(1), mxbval(1);
        if (nbval[1 - 1] < 1) {
            write(nout, format_9989), "   NB ", nbval(1), 1;
            fatal = true;
        } else if (nbmin[1 - 1] < 1) {
            write(nout, format_9989), "NBMIN ", nbmin(1), 1;
            fatal = true;
        } else if (nxval[1 - 1] < 1) {
            write(nout, format_9989), "   NX ", nxval(1), 1;
            fatal = true;
        } else if (nsval[1 - 1] < 2) {
            write(nout, format_9989), "   NS ", nsval(1), 2;
            fatal = true;
        } else if (mxbval[1 - 1] < 1) {
            write(nout, format_9989), " MAXB ", mxbval(1), 1;
            fatal = true;
        }
        xlaenv(1, nbval[1 - 1]);
        xlaenv(2, nbmin[1 - 1]);
        xlaenv(3, nxval[1 - 1]);
        xlaenv(4, nsval[1 - 1]);
        xlaenv(8, mxbval[1 - 1]);
        write(nout, format_9983), "NB:   ", nbval(1);
        write(nout, format_9983), "NBMIN:", nbmin(1);
        write(nout, format_9983), "NX:   ", nxval(1);
        write(nout, format_9983), "NS:   ", nsval(1);
        write(nout, format_9983), "MAXB: ", mxbval(1);
        //
    } else if (!dsb && !glm && !gqr && !gsv && !csd && !lse) {
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
        if (!dbb) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, nbval(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (nbval[i - 1] < 0) {
                    write(nout, format_9989), "   NB ", nbval(i), 0;
                    fatal = true;
                } else if (nbval[i - 1] > nmax) {
                    write(nout, format_9988), "   NB ", nbval(i), nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NB:   ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nbval(i);
                }
            }
        }
        //
        //        Read the values of NBMIN
        //
        if (nep || sep || svd || dgg) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, nbmin(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (nbmin[i - 1] < 0) {
                    write(nout, format_9989), "NBMIN ", nbmin(i), 0;
                    fatal = true;
                } else if (nbmin[i - 1] > nmax) {
                    write(nout, format_9988), "NBMIN ", nbmin(i), nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NBMIN:";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nbmin(i);
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
                    rloop, nxval(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (nxval[i - 1] < 0) {
                    write(nout, format_9989), "   NX ", nxval(i), 0;
                    fatal = true;
                } else if (nxval[i - 1] > nmax) {
                    write(nout, format_9988), "   NX ", nxval(i), nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NX:   ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nxval(i);
                }
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
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, nsval(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (nsval[i - 1] < 0) {
                    write(nout, format_9989), "   NS ", nsval(i), 0;
                    fatal = true;
                } else if (nsval[i - 1] > nmax) {
                    write(nout, format_9988), "   NS ", nsval(i), nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NS:   ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nsval(i);
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
        if (dgg) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, mxbval(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (mxbval[i - 1] < 0) {
                    write(nout, format_9989), " MAXB ", mxbval(i), 0;
                    fatal = true;
                } else if (mxbval[i - 1] > nmax) {
                    write(nout, format_9988), " MAXB ", mxbval(i), nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "MAXB: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, mxbval(i);
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
                    rloop, inmin(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (inmin[i - 1] < 0) {
                    write(nout, format_9989), " INMIN ", inmin(i), 0;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "INMIN: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, inmin(i);
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
                    rloop, inwin(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (inwin[i - 1] < 0) {
                    write(nout, format_9989), " INWIN ", inwin(i), 0;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "INWIN: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, inwin(i);
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
                    rloop, inibl(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (inibl[i - 1] < 0) {
                    write(nout, format_9989), " INIBL ", inibl(i), 0;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "INIBL: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, inibl(i);
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
                    rloop, ishfts(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (ishfts[i - 1] < 0) {
                    write(nout, format_9989), " ISHFTS ", ishfts(i), 0;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "ISHFTS: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, ishfts(i);
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
        if (nep || dgg) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, iacc22(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (iacc22[i - 1] < 0) {
                    write(nout, format_9989), " IACC22 ", iacc22(i), 0;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "IACC22: ";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, iacc22(i);
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
        if (dgg) {
            {
                read_loop rloop(cmn, nin, star);
                for (i = 1; i <= nparms; i = i + 1) {
                    rloop, nbcol(i);
                }
            }
            for (i = 1; i <= nparms; i = i + 1) {
                if (nbcol[i - 1] < 0) {
                    write(nout, format_9989), "NBCOL ", nbcol(i), 0;
                    fatal = true;
                } else if (nbcol[i - 1] > nmax) {
                    write(nout, format_9988), "NBCOL ", nbcol(i), nmax;
                    fatal = true;
                }
            }
            {
                write_loop wloop(cmn, nout, format_9983);
                wloop, "NBCOL:";
                for (i = 1; i <= nparms; i = i + 1) {
                    wloop, nbcol(i);
                }
            }
        } else {
            for (i = 1; i <= nparms; i = i + 1) {
                nbcol[i - 1] = 1;
            }
        }
    }
    //
    //     Calculate and prINTEGER the machine dependent constants.
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
    if (sep || svd || dgg) {
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
                rloop, ioldsd(i);
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
    if (!(dgx || dxv)) {
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
        if (!(dev || des || dvx || dsx || dgv || dgs) && ntypes <= 0) {
            write(nout, format_9990), c3;
            goto statement_200;
        }
        //
    } else {
        if (dxv) {
            c3 = "DXV";
        }
        if (dgx) {
            c3 = "DGX";
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
                c3, nbval(i), nbmin(i), nxval(i), max((INTEGER)11, inmin(i)), inwin(i), inibl(i), ishfts(i), iacc22(i);
            Rchkhs(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], nmax, &a[(1 - 1) + (6 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], &a[(1 - 1) + (8 - 1) * lda], &a[(1 - 1) + (9 - 1) * lda], &a[(1 - 1) + (10 - 1) * lda], &a[(1 - 1) + (11 - 1) * lda], &a[(1 - 1) + (12 - 1) * lda], &d[(7 - 1) * ldd], work, lwork, iwork, logwrk, result, info);
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
            write(nout, format_9997), c3, nbval(i), nbmin(i), nxval(i);
            if (tstchk) {
                if (Mlsamen(3, c3, "SE2")) {
                    Rchkst2stg(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], &d[(7 - 1) * ldd], &d[(8 - 1) * ldd], &d[(9 - 1) * ldd], &d[(10 - 1) * ldd], &d[(11 - 1) * ldd], &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &d[(12 - 1) * ldd], &a[(1 - 1) + (6 - 1) * lda], work, lwork, iwork, liwork, result, info);
                } else {
                    Rchkst(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], &d[(7 - 1) * ldd], &d[(8 - 1) * ldd], &d[(9 - 1) * ldd], &d[(10 - 1) * ldd], &d[(11 - 1) * ldd], &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &d[(12 - 1) * ldd], &a[(1 - 1) + (6 - 1) * lda], work, lwork, iwork, liwork, result, info);
                }
                if (info != 0) {
                    write(nout, format_9980), "Rchkst", info;
                }
            }
            if (tstdrv) {
                if (Mlsamen(3, c3, "SE2")) {
                    Rdrvst2stg(nn, nval, 18, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], &d[(8 - 1) * ldd], &d[(9 - 1) * ldd], &d[(10 - 1) * ldd], &d[(11 - 1) * ldd], &a[(1 - 1) + (2 - 1) * lda], nmax, &a[(1 - 1) + (3 - 1) * lda], &d[(12 - 1) * ldd], &a[(1 - 1) + (4 - 1) * lda], work, lwork, iwork, liwork, result, info);
                } else {
                    Rdrvst(nn, nval, 18, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], &d[(8 - 1) * ldd], &d[(9 - 1) * ldd], &d[(10 - 1) * ldd], &d[(11 - 1) * ldd], &a[(1 - 1) + (2 - 1) * lda], nmax, &a[(1 - 1) + (3 - 1) * lda], &d[(12 - 1) * ldd], &a[(1 - 1) + (4 - 1) * lda], work, lwork, iwork, liwork, result, info);
                }
                if (info != 0) {
                    write(nout, format_9980), "Rdrvst", info;
                }
            }
        }
        //
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
            write(nout, format_9997), c3, nbval(i), nbmin(i), nxval(i);
            if (tstchk) {
                //               CALL Rdrvsg( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
                //     $                      NOUT, &a[(1-1)+(1-1)*lda], NMAX, &a[(1-1)+(2-1)*lda], NMAX,
                //     $                      D( 1, 3 ), &a[(1-1)+(3-1)*lda], NMAX, &a[(1-1)+(4-1)*lda],
                //     $                      &a[(1-1)+(5-1)*lda], &a[(1-1)+(6-1)*lda], &a[(1-1)+(7-1)*lda], WORK,
                //     $                      LWORK, IWORK, LIWORK, RESULT, INFO )
                Rdrvsg2stg(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], nmax, &d[(3 - 1) * ldd], &d[(3 - 1) * ldd], &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], work, lwork, iwork, liwork, result, info);
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
            write(nout, "(/,/,1x,a3,':  NB =',i4,', NBMIN =',i4,', NX =',i4,', NRHS =',i4)"), c3, nbval(i), nbmin(i), nxval(i), nrhs;
            if (tstchk) {
                Rchkbd(nn, mval, nval, maxtyp, dotype, nrhs, iseed, thresh, &a[(1 - 1) + (1 - 1) * lda], nmax, &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &a[(1 - 1) + (2 - 1) * lda], nmax, &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], nmax, &a[(1 - 1) + (6 - 1) * lda], nmax, &a[(1 - 1) + (7 - 1) * lda], &a[(1 - 1) + (8 - 1) * lda], work, lwork, iwork, nout, info);
                if (info != 0) {
                    write(nout, format_9980), "Rchkbd", info;
                }
            }
            if (tstdrv) {
                Rdrvbd(nn, mval, nval, maxtyp, dotype, iseed, thresh, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], nmax, &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], work, lwork, iwork, nout, info);
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
            Rdrvev(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], nmax, &a[(1 - 1) + (5 - 1) * lda], nmax, result, work, lwork, iwork, info);
            if (info != 0) {
                write(nout, format_9980), "Rgeev", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
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
            Rdrves(nn, nval, ntypes, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &a[(1 - 1) + (4 - 1) * lda], nmax, result, work, lwork, iwork, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Rgees", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
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
            Rdrvvx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &a[(1 - 1) + (3 - 1) * lda], nmax, &a[(1 - 1) + (4 - 1) * lda], nmax, &a[(1 - 1) + (5 - 1) * lda], nmax, &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], &d[(7 - 1) * ldd], &d[(8 - 1) * ldd], &d[(9 - 1) * ldd], &d[(10 - 1) * ldd], &d[(11 - 1) * ldd], &d[(12 - 1) * ldd], result, work, lwork, iwork, info);
            if (info != 0) {
                write(nout, format_9980), "Rgeevx", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
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
            Rdrvsx(nn, nval, ntypes, dotype, iseed, thresh, nin, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], &a[(1 - 1) + (4 - 1) * lda], nmax, &a[(1 - 1) + (5 - 1) * lda], result, work, lwork, iwork, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Rgeesx", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
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
                c3, nbval(i), nbmin(i), nsval(i), mxbval(i), iacc22(i), nbcol(i);
            tstdif = false;
            thrshn = 10.0;
            if (tstchk) {
                Rchkgg(nn, nval, maxtyp, dotype, iseed, thresh, tstdif, thrshn, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], &a[(1 - 1) + (8 - 1) * lda], &a[(1 - 1) + (9 - 1) * lda], nmax, &a[(1 - 1) + (10 - 1) * lda], &a[(1 - 1) + (11 - 1) * lda], &a[(1 - 1) + (12 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], &a[(1 - 1) + (13 - 1) * lda], &a[(1 - 1) + (14 - 1) * lda], work, lwork, logwrk, result, info);
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
            Rdrges(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], nmax, &a[(1 - 1) + (8 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], work, lwork, result, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Rdrges", info;
            }
            //
            //     Blocked version
            //
            xlaenv(16, 2);
            Rdrges3(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], nmax, &a[(1 - 1) + (8 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], work, lwork, result, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Rdrges3", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
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
            Rdrgsx(nn, ncmax, thresh, nin, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &c[(1 - 1) + (1 - 1) * ldc], ncmax * ncmax, &a[(1 - 1) + (12 - 1) * lda], work, lwork, iwork, liwork, logwrk, info);
            if (info != 0) {
                write(nout, format_9980), "Rdrgsx", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
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
            Rdrgev(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], nmax, &a[(1 - 1) + (8 - 1) * lda], &a[(1 - 1) + (9 - 1) * lda], nmax, &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], work, lwork, result, info);
            if (info != 0) {
                write(nout, format_9980), "Rdrgev", info;
            }
            //
            //     Blocked version
            //
            Rdrgev3(nn, nval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], nmax, &a[(1 - 1) + (8 - 1) * lda], &a[(1 - 1) + (9 - 1) * lda], nmax, &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], work, lwork, result, info);
            if (info != 0) {
                write(nout, format_9980), "Rdrgev3", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
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
            Rdrgvx(nn, thresh, nin, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], &iwork[1 - 1], &iwork[2 - 1], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &d[(6 - 1) * ldd], &d[(7 - 1) * ldd], &d[(8 - 1) * ldd], &d[(9 - 1) * ldd], work, lwork, &iwork[3 - 1], liwork - 2, result, logwrk, info);
            //
            if (info != 0) {
                write(nout, format_9980), "Rdrgvx", info;
            }
        }
        write(nout, format_9973);
        goto statement_10;
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
        Rchksb2stg(nn, nval, nk, kval, maxtyp, dotype, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &d[(1 - 1)], &d[(2 - 1) * ldd], &d[(3 - 1) * ldd], &d[(4 - 1) * ldd], &d[(5 - 1) * ldd], &a[(1 - 1) + (2 - 1) * lda], nmax, work, lwork, result, info);
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
            Rchkbb(nn, mval, nval, nk, kval, maxtyp, dotype, nrhs, iseed, thresh, nout, &a[(1 - 1) + (1 - 1) * lda], nmax, &a[(1 - 1) + (2 - 1) * lda], 2 * nmax, &d[(1 - 1)], &d[(2 - 1) * ldd], &a[(1 - 1) + (4 - 1) * lda], nmax, &a[(1 - 1) + (5 - 1) * lda], nmax, &a[(1 - 1) + (6 - 1) * lda], nmax, &a[(1 - 1) + (7 - 1) * lda], work, lwork, result, info);
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
        Rckglm(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &b[(1 - 1) + (1 - 1) * ldb], &b[(1 - 1) + (2 - 1) * ldb], x, work, &d[(1 - 1)], nin, nout, info);
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
        Rckgqr(nn, mval, nn, pval, nn, nval, ntypes, iseed, thresh, nmax, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], taua, &b[(1 - 1) + (1 - 1) * ldb], &b[(1 - 1) + (2 - 1) * ldb], &b[(1 - 1) + (3 - 1) * ldb], &b[(1 - 1) + (4 - 1) * ldb], &b[(1 - 1) + (5 - 1) * ldb], taub, work, &d[(1 - 1)], nin, nout, info);
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
        Rckgsv(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &b[(1 - 1) + (1 - 1) * ldb], &b[(1 - 1) + (2 - 1) * ldb], &a[(1 - 1) + (3 - 1) * lda], &b[(1 - 1) + (3 - 1) * ldb], &a[(1 - 1) + (4 - 1) * lda], taua, taub, &b[(1 - 1) + (4 - 1) * ldb], iwork, work, &d[(1 - 1)], nin, nout, info);
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
        Rckcsd(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &a[(1 - 1) + (4 - 1) * lda], &a[(1 - 1) + (5 - 1) * lda], &a[(1 - 1) + (6 - 1) * lda], &a[(1 - 1) + (7 - 1) * lda], iwork, work, &d[(1 - 1)], nin, nout, info);
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
        Rcklse(nn, mval, pval, nval, ntypes, iseed, thresh, nmax, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &b[(1 - 1) + (1 - 1) * ldb], &b[(1 - 1) + (2 - 1) * ldb], x, work, &d[(1 - 1)], nin, nout, info);
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
        goto statement_190;
    }
statement_380:
    write(nout, "(/,/,' End of tests')");
    s2 = dsecnd[-1];
    write(nout, "(' Total time used = ',f12.2,' seconds',/)"), s2 - s1;
    //
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,stat=allocatestatus)");
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(b,stat=allocatestatus)");
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(c,stat=allocatestatus)");
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(work,stat=allocatestatus)");
    //
    //     End of Rchkee
    //
}

INTEGER main(INTEGER argc, char const *argv[]) { return main_with_catch(argc, argv, placeholder_please_replace::program_Rchkee); }
