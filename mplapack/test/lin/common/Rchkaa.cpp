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

// Derived from LAPACK routine DCHKAA.
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
#include <mplapack_lin.h>

#include <memory>

void Rchkaa(void) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    fem::str<10> intstr = "0123456789";
    REAL threq = 2.0;
    const INTEGER nmax = 132;
    const INTEGER kdmax = nmax + (nmax + 1) / 4;
    REAL s1 = 0.0;
    INTEGER lda = 0;
    bool fatal = false;
    const INTEGER nin = 5;
    const INTEGER nout = 6;
    INTEGER mplapack_vers_major = 0;
    INTEGER mplapack_vers_minor = 0;
    INTEGER mplapack_vers_patch = 0;
    INTEGER lapack_vers_major = 0;
    INTEGER lapack_vers_minor = 0;
    INTEGER lapack_vers_patch = 0;
    INTEGER nm = 0;
    const INTEGER maxin = 12;
    const INTEGER maxrhs = 16;
    INTEGER mval[maxin];
    INTEGER i = 0;
    INTEGER nn = 0;
    INTEGER nval[maxin];
    INTEGER nns = 0;
    INTEGER nsval[maxin];
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
    REAL eps = 0.0;
    fem::str<72> aline;
    fem::str<3> path;
    const INTEGER matmax = 30;
    INTEGER nmats = 0;
    fem::str<1> c1;
    INTEGER k = 0;
    INTEGER ic = 0;
    fem::str<2> c2;
    INTEGER nrhs = 0;
    INTEGER ntypes = 0;
    bool dotype[matmax];
    auto a_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, ((kdmax + 1) * nmax) * 7));
    REAL *a = a_storage.get();
    auto b_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, (nmax * maxrhs) * 4));
    REAL *b = b_storage.get();
    auto work_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * (3 * nmax + maxrhs + 30)));
    REAL *work = work_storage.get();
    auto rwork_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, 5 * nmax + 2 * maxrhs));
    REAL *rwork = rwork_storage.get();
    INTEGER iwork[25 * nmax];
    REAL s[2 * nmax];
    INTEGER la = 0;
    INTEGER lafac = 0;
    INTEGER piv[nmax];
    REAL e[nmax];
    REAL s2 = 0.0;
    INTEGER ldaw = (kdmax + 1) * nmax;
    INTEGER ldb = nmax * maxrhs;
    //
    static const char *format_9999 = "(/,' Execution not attempted due to input errors')";
    static const char *format_9998 = "(/,' End of tests')";
    static const char *format_9997 = "(' Total time used = ',f12.2,' seconds',/)";
    static const char *format_9996 = "(' Invalid input value: ',a4,'=',i6,'; must be >=',i6)";
    static const char *format_9995 = "(' Invalid input value: ',a4,'=',i6,'; must be <=',i6)";
    static const char *format_9994 = "(' Tests of the Multiple precision version of LAPACK ',i1,'.',i1,'.',i1,/, "
                                     "' Based on the original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, "
                                     "' The following parameter values will be used:')";
    static const char *format_9993 = "(4x,a4,':  ',10i6,/,11x,10i6)";
    static const char *format_9992 = "(/,' Routines pass computational tests if test ratio is ','less than',"
                                     "f8.2,/)";
    static const char *format_9991 = "(' Relative machine ',a,' is taken to be',d16.6)";
    static const char *format_9990 = "(/,1x,a3,':  Unrecognized path name')";
    static const char *format_9989 = "(/,1x,a3,' routines were not tested')";
    static const char *format_9988 = "(/,1x,a3,' driver routines were not tested')";
    //
    s1 = dsecnd();
    lda = nmax;
    fatal = false;
    //
    // Read a dummy line.
    //
    read(nin, star);
    //
    // Report values of parameters.
    //
    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
    write(nout, format_9994), mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
    //
    // Read the values of M
    //
    read(nin, star), nm;
    if (nm < 1) {
        write(nout, format_9996), " NM ", nm, 1;
        nm = 0;
        fatal = true;
    } else if (nm > maxin) {
        write(nout, format_9995), " NM ", nm, maxin;
        nm = 0;
        fatal = true;
    }
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= nm; i = i + 1) {
            rloop, mval[i - 1];
        }
    }
    for (i = 1; i <= nm; i = i + 1) {
        if (mval[i - 1] < 0) {
            write(nout, format_9996), " M  ", mval[i - 1], 0;
            fatal = true;
        } else if (mval[i - 1] > nmax) {
            write(nout, format_9995), " M  ", mval[i - 1], nmax;
            fatal = true;
        }
    }
    if (nm > 0) {
        {
            write_loop wloop(cmn, nout, format_9993);
            wloop, "M   ";
            for (i = 1; i <= nm; i = i + 1) {
                wloop, mval[i - 1];
            }
        }
    }
    //
    // Read the values of N
    //
    read(nin, star), nn;
    if (nn < 1) {
        write(nout, format_9996), " NN ", nn, 1;
        nn = 0;
        fatal = true;
    } else if (nn > maxin) {
        write(nout, format_9995), " NN ", nn, maxin;
        nn = 0;
        fatal = true;
    }
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= nn; i = i + 1) {
            rloop, nval[i - 1];
        }
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
        {
            write_loop wloop(cmn, nout, format_9993);
            wloop, "N   ";
            for (i = 1; i <= nn; i = i + 1) {
                wloop, nval[i - 1];
            }
        }
    }
    //
    // Read the values of NRHS
    //
    read(nin, star), nns;
    if (nns < 1) {
        write(nout, format_9996), " NNS", nns, 1;
        nns = 0;
        fatal = true;
    } else if (nns > maxin) {
        write(nout, format_9995), " NNS", nns, maxin;
        nns = 0;
        fatal = true;
    }
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= nns; i = i + 1) {
            rloop, nsval[i - 1];
        }
    }
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
        {
            write_loop wloop(cmn, nout, format_9993);
            wloop, "NRHS";
            for (i = 1; i <= nns; i = i + 1) {
                wloop, nsval[i - 1];
            }
        }
    }
    //
    // Read the values of NB
    //
    read(nin, star), nnb;
    if (nnb < 1) {
        write(nout, format_9996), "NNB ", nnb, 1;
        nnb = 0;
        fatal = true;
    } else if (nnb > maxin) {
        write(nout, format_9995), "NNB ", nnb, maxin;
        nnb = 0;
        fatal = true;
    }
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= nnb; i = i + 1) {
            rloop, nbval[i - 1];
        }
    }
    for (i = 1; i <= nnb; i = i + 1) {
        if (nbval[i - 1] < 0) {
            write(nout, format_9996), " NB ", nbval[i - 1], 0;
            fatal = true;
        }
    }
    if (nnb > 0) {
        {
            write_loop wloop(cmn, nout, format_9993);
            wloop, "NB  ";
            for (i = 1; i <= nnb; i = i + 1) {
                wloop, nbval[i - 1];
            }
        }
    }
    //
    // Set NBVAL2 to be the set of unique values of NB
    //
    nnb2 = 0;
    for (i = 1; i <= nnb; i = i + 1) {
        nb = nbval[i - 1];
        for (j = 1; j <= nnb2; j = j + 1) {
            if (nb == nbval2[j - 1]) {
                goto statement_60;
            }
        }
        nnb2++;
        nbval2[nnb2 - 1] = nb;
    statement_60:;
    }
    //
    // Read the values of NX
    //
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= nnb; i = i + 1) {
            rloop, nxval[i - 1];
        }
    }
    for (i = 1; i <= nnb; i = i + 1) {
        if (nxval[i - 1] < 0) {
            write(nout, format_9996), " NX ", nxval[i - 1], 0;
            fatal = true;
        }
    }
    if (nnb > 0) {
        {
            write_loop wloop(cmn, nout, format_9993);
            wloop, "NX  ";
            for (i = 1; i <= nnb; i = i + 1) {
                wloop, nxval[i - 1];
            }
        }
    }
    //
    // Read the values of RANKVAL
    //
    read(nin, star), nrank;
    if (nn < 1) {
        write(nout, format_9996), " NRANK ", nrank, 1;
        nrank = 0;
        fatal = true;
    } else if (nn > maxin) {
        write(nout, format_9995), " NRANK ", nrank, maxin;
        nrank = 0;
        fatal = true;
    }
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= nrank; i = i + 1) {
            rloop, rankval[i - 1];
        }
    }
    for (i = 1; i <= nrank; i = i + 1) {
        if (rankval[i - 1] < 0) {
            write(nout, format_9996), " RANK  ", rankval[i - 1], 0;
            fatal = true;
        } else if (rankval[i - 1] > 100) {
            write(nout, format_9995), " RANK  ", rankval[i - 1], 100;
            fatal = true;
        }
    }
    if (nrank > 0) {
        {
            write_loop wloop(cmn, nout, format_9993);
            wloop, "RANK % OF N";
            for (i = 1; i <= nrank; i = i + 1) {
                wloop, rankval[i - 1];
            }
        }
    }
    //
    // Read the threshold value for the test ratios.
    //
    read(nin, star), thresh;
    write(nout, format_9992), thresh;
    //
    // Read the flag that indicates whether to test the LAPACK routines.
    //
    read(nin, star), tstchk;
    //
    // Read the flag that indicates whether to test the driver routines.
    //
    read(nin, star), tstdrv;
    //
    // Read the flag that indicates whether to test the error exits.
    //
    read(nin, star), tsterr;
    //
    if (fatal) {
        write(nout, format_9999);
        FEM_STOP(0);
    }
    //
    // Calculate and print the machine dependent constants.
    //
    eps = Rlamch("Underflow threshold");
    write(nout, format_9991), "underflow", eps;
    eps = Rlamch("Overflow threshold");
    write(nout, format_9991), "overflow ", eps;
    eps = Rlamch("Epsilon");
    write(nout, format_9991), "precision", eps;
    write(nout, star);
//
statement_80:
    //
    // Read a test path and the number of matrix types to use.
    //
    try {
        read(nin, "(a72)"), aline;
    } catch (fem::read_end const &) {
        goto statement_140;
    }
    path = aline(1, 3);
    nmats = matmax;
    i = 3;
statement_90:
    i++;
    if (i > 72) {
        nmats = matmax;
        goto statement_130;
    }
    if (aline(i, i) == " ") {
        goto statement_90;
    }
    nmats = 0;
statement_100:
    c1 = aline(i, i);
    for (k = 1; k <= 10; k = k + 1) {
        if (c1 == intstr(k, k)) {
            ic = k - 1;
            goto statement_120;
        }
    }
    goto statement_130;
statement_120:
    nmats = nmats * 10 + ic;
    i++;
    if (i > 72) {
        goto statement_130;
    }
    goto statement_100;
statement_130:
    c1 = path(1, 1);
    c2 = path(2, 3);
    nrhs = nsval[1 - 1];
    //
    // Check first character for correct precision.
    //
    if (!Mlsame(c1.elems, "Double precision") && !Mlsame(c1.elems, "R")) {
        write(nout, format_9990), path;
        //
    } else if (nmats <= 0) {
        //
        // Check for a positive number of tests requested.
        //
        write(nout, format_9989), path;
        //
    } else if (Mlsamen(2, c2.elems, "GE")) {
        //
        // GE:  general matrices
        //
        ntypes = 11;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkge(dotype, nm, mval, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvge(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], s, work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "GB")) {
        //
        // GB:  general banded matrices
        //
        la = (2 * kdmax + 1) * nmax;
        lafac = (3 * kdmax + 1) * nmax;
        ntypes = 8;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkgb(dotype, nm, mval, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, &a[0], la, &a[(3 - 1) * ldaw], lafac, &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvgb(dotype, nn, nval, nrhs, thresh, tsterr, &a[0], la, &a[(3 - 1) * ldaw], lafac, &a[(6 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], s, work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "GT")) {
        //
        // GT:  general tridiagonal matrices
        //
        ntypes = 12;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkgt(dotype, nn, nval, nns, nsval, thresh, tsterr, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvgt(dotype, nn, nval, nrhs, thresh, tsterr, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "PO")) {
        //
        // PO:  positive definite matrices
        //
        ntypes = 9;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkpo(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvpo(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], s, work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "PS")) {
        //
        // PS:  positive semi-definite matrices
        //
        ntypes = 9;
        //
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkps(dotype, nn, nval, nnb2, nbval2, nrank, rankval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], piv, work, rwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "PP")) {
        //
        // PP:  positive definite packed matrices
        //
        ntypes = 9;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkpp(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvpp(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], s, work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "PB")) {
        //
        // PB:  positive definite banded matrices
        //
        ntypes = 8;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkpb(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvpb(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], s, work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "PT")) {
        //
        // PT:  positive definite tridiagonal matrices
        //
        ntypes = 12;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkpt(dotype, nn, nval, nns, nsval, thresh, tsterr, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvpt(dotype, nn, nval, nrhs, thresh, tsterr, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "SY")) {
        //
        // SY:  symmetric indefinite matrices,
        // with partial (Bunch-Kaufman) pivoting algorithm
        //
        ntypes = 10;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchksy(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvsy(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "SR")) {
        //
        // SR:  symmetric indefinite matrices,
        // with bounded Bunch-Kaufman (rook) pivoting algorithm
        //
        ntypes = 10;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchksy_rook(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvsy_rook(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "SK")) {
        //
        // SK:  symmetric indefinite matrices,
        // with bounded Bunch-Kaufman (rook) pivoting algorithm,
        // different matrix storage format than SR path version.
        //
        ntypes = 10;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchksy_rk(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], e, &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvsy_rk(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], e, &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "SA")) {
        //
        // SA:  symmetric indefinite matrices,
        // with partial (Aasen's) pivoting algorithm
        //
        ntypes = 10;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchksy_aa(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvsy_aa(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "S2")) {
        //
        // SA:  symmetric indefinite matrices,
        // with partial (Aasen's) pivoting algorithm
        //
        ntypes = 10;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchksy_aa_2stage(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvsy_aa_2stage(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "SP")) {
        //
        // SP:  symmetric indefinite packed matrices,
        // with partial (Bunch-Kaufman) pivoting algorithm
        //
        ntypes = 10;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchksp(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
        if (tstdrv) {
            Rdrvsp(dotype, nn, nval, nrhs, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "TR")) {
        //
        // TR:  triangular matrices
        //
        ntypes = 18;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchktr(dotype, nn, nval, nnb2, nbval2, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "TP")) {
        //
        // TP:  triangular packed matrices
        //
        ntypes = 18;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchktp(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "TB")) {
        //
        // TB:  triangular banded matrices
        //
        ntypes = 17;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchktb(dotype, nn, nval, nns, nsval, thresh, tsterr, lda, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "QR")) {
        //
        // QR:  QR factorization
        //
        ntypes = 8;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkqr(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &a[(4 - 1) * ldaw], &a[(5 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "LQ")) {
        //
        // LQ:  LQ factorization
        //
        ntypes = 8;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchklq(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &a[(4 - 1) * ldaw], &a[(5 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, rwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "QL")) {
        //
        // QL:  QL factorization
        //
        ntypes = 8;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkql(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &a[(4 - 1) * ldaw], &a[(5 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, rwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "RQ")) {
        //
        // RQ:  RQ factorization
        //
        ntypes = 8;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkrq(dotype, nm, mval, nn, nval, nnb, nbval, nxval, nrhs, thresh, tsterr, nmax, &a[0], &a[(2 - 1) * ldaw], &a[(3 - 1) * ldaw], &a[(4 - 1) * ldaw], &a[(5 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, rwork, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "QP")) {
        //
        // QP:  QR factorization with pivoting
        //
        ntypes = 6;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkq3(dotype, nm, mval, nn, nval, nnb, nbval, nxval, thresh, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(3 - 1) * ldb], work, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "QK")) {
        //
        // QK: truncated QR factorization with pivoting
        //
        ntypes = 19;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkqp3rk(dotype, nm, mval, nn, nval, nns, nsval, nnb, nbval, nxval, thresh, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "QK")) {
        //
        // QK: truncated QR factorization with pivoting
        //
        ntypes = 19;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchkqp3rk(dotype, nm, mval, nn, nval, nns, nsval, nnb, nbval, nxval, thresh, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], &b[(4 - 1) * ldb], work, iwork, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "TZ")) {
        //
        // TZ:  Trapezoidal matrix
        //
        ntypes = 3;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstchk) {
            Rchktz(dotype, nm, mval, nn, nval, thresh, tsterr, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(3 - 1) * ldb], work, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "LS")) {
        //
        // LS:  Least squares drivers
        //
        ntypes = 6;
        Alareq(path, nmats, dotype, ntypes, nin, nout);
        //
        if (tstdrv) {
            Rdrvls(dotype, nm, mval, nn, nval, nns, nsval, nnb, nbval, nxval, thresh, tsterr, &a[0], &a[(2 - 1) * ldaw], &b[0], &b[(2 - 1) * ldb], &b[(3 - 1) * ldb], rwork, &rwork[(nmax + 1) - 1], nout);
        } else {
            write(nout, format_9988), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "EQ")) {
        //
        // EQ:  Equilibration routines for general and positive definite
        // matrices (THREQ should be between 2 and 10)
        //
        if (tstchk) {
            Rchkeq(threq, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "QT")) {
        //
        // QT:  QRT routines for general matrices
        //
        if (tstchk) {
            Rchkqrt(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "QX")) {
        //
        // QX:  QRT routines for triangular-pentagonal matrices
        //
        if (tstchk) {
            Rchkqrtp(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "TQ")) {
        //
        // TQ:  LQT routines for general matrices
        //
        if (tstchk) {
            Rchklqt(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "XQ")) {
        //
        // XQ:  LQT routines for triangular-pentagonal matrices
        //
        if (tstchk) {
            Rchklqtp(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "TS")) {
        //
        // TS:  QR routines for tall-skinny matrices
        //
        if (tstchk) {
            Rchktsqr(thresh, tsterr, nm, mval, nn, nval, nnb, nbval, nout);
        } else {
            write(nout, format_9989), path;
        }
        //
    } else if (Mlsamen(2, c2.elems, "HH")) {
        //
        // HH:  Householder reconstruction for tall-skinny matrices
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
    //
    // Go back to get another input line.
    //
    goto statement_80;
//
// Branch to this line when the last record is read.
//
statement_140:
    cmn.io.close(nin);
    s2 = dsecnd();
    write(nout, format_9998);
    write(nout, format_9997), cast2double(s2 - s1);
    //
    // End of Rchkaa
    //
}

int main(int argc, char const *argv[]) { Rchkaa(); }
