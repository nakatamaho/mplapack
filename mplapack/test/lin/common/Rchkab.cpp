/*
 * Copyright (c) 2008-2021
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

void program_Rchkab(INTEGER argc, char const *argv[]) {
    common cmn(argc, argv);
    FEM_CMN_SVE(program_Rchkab);
    common_read read(cmn);
    common_write write(cmn);
    str<10> &intstr = sve.intstr;
    if (is_called_first_time) {
        intstr = "0123456789";
    }
    REAL s1 = 0.0;
    const INTEGER nmax = 132;
    INTEGER lda = 0;
    bool fatal = false;
    const INTEGER nin = 5;
    INTEGER vers_major = 0;
    INTEGER vers_minor = 0;
    INTEGER vers_patch = 0;
    const INTEGER nout = 6;
    INTEGER nm = 0;
    const INTEGER maxin = 12;
    arr_1d<maxin, int> mval(fill0);
    INTEGER i = 0;
    INTEGER nns = 0;
    arr_1d<maxin, int> nsval(fill0);
    const INTEGER maxrhs = 16;
    REAL thresh = 0.0;
    bool tstdrv = false;
    bool tsterr = false;
    float seps = float0;
    REAL eps = 0.0;
    str<72> aline = char0;
    str<3> path = char0;
    const INTEGER matmax = 30;
    INTEGER nmats = 0;
    char c1 = char0;
    INTEGER k = 0;
    INTEGER ic = 0;
    str<2> c2 = char0;
    INTEGER nrhs = 0;
    INTEGER ntypes = 0;
    arr_1d<matmax, bool> dotype(fill0);
    const INTEGER ldamax = nmax;
    arr_1d<nmax, REAL> rwork(fill0);
    arr_1d<nmax, int> iwork(fill0);
    REAL s2 = 0.0;
    static const char *format_9989 = "(/,1x,a6,' driver routines were not tested')";
    static const char *format_9991 = "(' Relative machine ',a,' is taken to be',d16.6)";
    static const char *format_9993 = "(4x,a4,':  ',10i6,/,11x,10i6)";
    static const char *format_9995 = "(' Invalid input value: ',a4,'=',i6,'; must be <=',i6)";
    static const char *format_9996 = "(' Invalid input value: ',a4,'=',i6,'; must be >=',i6)";
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    s1 = dsecnd[-1];
    lda = nmax;
    fatal = false;
    //
    //     Read a dummy line.
    //
    read(nin, star);
    //
    //     Report values of parameters.
    //
    ilaver(vers_major, vers_minor, vers_patch);
    write(nout, "(' Tests of the DOUBLE PRECISION LAPACK DSGESV/DSPOSV',' routines ',/,"
                "' LAPACK VERSION ',i1,'.',i1,'.',i1,/,/,"
                "' The following parameter values will be used:')"),
        vers_major, vers_minor, vers_patch;
    //
    //     Read the values of M
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
            rloop, mval(i);
        }
    }
    for (i = 1; i <= nm; i = i + 1) {
        if (mval[i - 1] < 0) {
            write(nout, format_9996), " M  ", mval(i), 0;
            fatal = true;
        } else if (mval[i - 1] > nmax) {
            write(nout, format_9995), " M  ", mval(i), nmax;
            fatal = true;
        }
    }
    if (nm > 0) {
        {
            write_loop wloop(cmn, nout, format_9993);
            wloop, "M   ";
            for (i = 1; i <= nm; i = i + 1) {
                wloop, mval(i);
            }
        }
    }
    //
    //     Read the values of NRHS
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
            rloop, nsval(i);
        }
    }
    for (i = 1; i <= nns; i = i + 1) {
        if (nsval[i - 1] < 0) {
            write(nout, format_9996), "NRHS", nsval(i), 0;
            fatal = true;
        } else if (nsval[i - 1] > maxrhs) {
            write(nout, format_9995), "NRHS", nsval(i), maxrhs;
            fatal = true;
        }
    }
    if (nns > 0) {
        {
            write_loop wloop(cmn, nout, format_9993);
            wloop, "NRHS";
            for (i = 1; i <= nns; i = i + 1) {
                wloop, nsval(i);
            }
        }
    }
    //
    //     Read the threshold value for the test ratios.
    //
    read(nin, star), thresh;
    write(nout, "(/,' Routines pass computational tests if test ratio is ','less than',"
                "f8.2,/)"),
        thresh;
    //
    //     Read the flag that indicates whether to test the driver routine.
    //
    read(nin, star), tstdrv;
    //
    //     Read the flag that indicates whether to test the error exits.
    //
    read(nin, star), tsterr;
    //
    if (fatal) {
        write(nout, "(/,' Execution not attempted due to input errors')");
        FEM_STOP(0);
    }
    //
    //     Calculate and prINTEGER the machine dependent constants.
    //
    seps = slamch[("Underflow threshold") - 1];
    write(nout, format_9991), "(single precision) underflow", seps;
    seps = slamch[("Overflow threshold") - 1];
    write(nout, format_9991), "(single precision) overflow ", seps;
    seps = slamch["Epsilon" - 1];
    write(nout, format_9991), "(single precision) precision", seps;
    write(nout, star);
    //
    eps = Rlamch("Underflow threshold");
    write(nout, format_9991), "(REAL precision) underflow", eps;
    eps = Rlamch("Overflow threshold");
    write(nout, format_9991), "(REAL precision) overflow ", eps;
    eps = Rlamch("Epsilon");
    write(nout, format_9991), "(REAL precision) precision", eps;
    write(nout, star);
//
statement_80:
    //
    //     Read a test path and the number of matrix types to use.
    //
    try {
        read(nin, "(a72)"), aline;
    } catch (read_end const ) {
        goto statement_140;
    }
    path = aline[(3 - 1) * ldaline];
    nmats = matmax;
    i = 3;
statement_90:
    i++;
    if (i > 72) {
        nmats = matmax;
        goto statement_130;
    }
    if (aline[(i - 1) + (i - 1) * ldaline] == " ") {
        goto statement_90;
    }
    nmats = 0;
statement_100:
    c1 = aline[(i - 1) + (i - 1) * ldaline];
    for (k = 1; k <= 10; k = k + 1) {
        if (c1 == intstr[(k - 1) + (k - 1) * ldintstr]) {
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
    c1 = path[(1 - 1)];
    c2 = path[(2 - 1) + (3 - 1) * ldpath];
    nrhs = nsval[1 - 1];
    //
    //     Check first character for correct precision.
    //
    if (!Mlsame(c1, "Double precision")) {
        write(nout, "(/,1x,a6,' routines were not tested')"), path;
        //
    } else if (nmats <= 0) {
        //
        //        Check for a positive number of tests requested.
        //
        write(nout, format_9989), path;
        goto statement_140;
        //
    } else if (Mlsamen2, c2, "GE") {
        //
        //        GE:  general matrices
        //
        ntypes = 11;
        Alareq("DGE", nmats, dotype, ntypes, nin, nout);
        //
        //        Test the error exits
        //
        if (tsterr) {
            Rerrab(nout);
        }
        //
        if (tstdrv) {
            Rdrvab(dotype, nm, mval, nns, nsval, thresh, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], work, rwork, swork, iwork, nout);
        } else {
            write(nout, format_9989), "DSGESV";
        }
        //
    } else if (Mlsamen2, c2, "PO") {
        //
        //        PO:  positive definite matrices
        //
        ntypes = 9;
        Alareq("DPO", nmats, dotype, ntypes, nin, nout);
        //
        if (tsterr) {
            Rerrac(nout);
        }
        //
        if (tstdrv) {
            Rdrvac(dotype, nm, mval, nns, nsval, thresh, lda, &a[(1 - 1)], &a[(2 - 1) * lda], &b[(1 - 1)], &b[(2 - 1) * ldb], work, rwork, swork, nout);
        } else {
            write(nout, format_9989), path;
        }
    } else {
        //
    }
    //
    //     Go back to get another input line.
    //
    goto statement_80;
//
//     Branch to this line when the last record is read.
//
statement_140:
    cmn.io.close(nin);
    s2 = dsecnd[-1];
    write(nout, "(/,' End of tests')");
    write(nout, "(' Total time used = ',f12.2,' seconds',/)"), s2 - s1;
    //
    //     End of Rchkab
    //
}

INTEGER main(INTEGER argc, char const *argv[]) { return main_with_catch(argc, argv, placeholder_please_replace::program_Rchkab); }
