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

// Derived from LAPACK routine DCHKRFP.
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

void Rchkrfp(void) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    //
    static const char *format_9999 = "(/,' Execution not attempted due to input errors')";
    static const char *format_9998 = "(/,' End of tests')";
    static const char *format_9997 = "(' Total time used = ',f12.2,' seconds',/)";
    static const char *format_9996 = "(' !! Invalid input value: ',a4,'=',i6,'; must be >=',i6)";
    static const char *format_9995 = "(' !! Invalid input value: ',a4,'=',i6,'; must be <=',i6)";
    static const char *format_9994 = "(' Tests of the Multiple precision version of LAPACK RFP routines',i1,'.',i1,'.',i1,/, "
                                     "' Based on original LAPACK VERSION ',i1,'.',i1,'.',i1,/,/, "
                                     "'The following parameter values will be used:')";
    static const char *format_9993 = "(4x,a4,':  ',10i6,/,11x,10i6)";
    static const char *format_9992 = "(/,' Routines pass computational tests if test ratio is ','less than',"
                                     "f8.2,/)";
    static const char *format_9991 = "(' Relative machine ',a,' is taken to be',d16.6)";
    //
    //
    REAL s1 = dsecnd();
    bool fatal = false;
    //
    // Read a dummy line.
    //
    const INTEGER nin = 5;
    const INTEGER nout = 6;
    read(nin, star);
    //
    // Report LAPACK version tag (e.g. LAPACK-3.2.0)
    //
    INTEGER mplapack_vers_major = 0;
    INTEGER mplapack_vers_minor = 0;
    INTEGER mplapack_vers_patch = 0;
    INTEGER lapack_vers_major = 0;
    INTEGER lapack_vers_minor = 0;
    INTEGER lapack_vers_patch = 0;
    iMlaver(mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch);
    write(nout, format_9994), mplapack_vers_major, mplapack_vers_minor, mplapack_vers_patch, lapack_vers_major, lapack_vers_minor, lapack_vers_patch;
    //
    // Read the values of N
    //
    INTEGER nn = 0;
    read(nin, star), nn;
    const INTEGER maxin = 12;
    if (nn < 1) {
        write(nout, format_9996), " NN ", nn, 1;
        nn = 0;
        fatal = true;
    } else if (nn > maxin) {
        write(nout, format_9995), " NN ", nn, maxin;
        nn = 0;
        fatal = true;
    }
    INTEGER nval[maxin];
    INTEGER i = 0;
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= nn; i = i + 1) {
            rloop, nval[i - 1];
        }
    }
    const INTEGER nmax = 50;
    for (i = 1; i <= nn; i = i + 1) {
        if (nval[i - 1] < 0) {
            write(nout, format_9996), " M  ", nval[i - 1], 0;
            fatal = true;
        } else if (nval[i - 1] > nmax) {
            write(nout, format_9995), " M  ", nval[i - 1], nmax;
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
    INTEGER nns = 0;
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
    INTEGER nsval[maxin];
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= nns; i = i + 1) {
            rloop, nsval[i - 1];
        }
    }
    const INTEGER maxrhs = 16;
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
    // Read the matrix types
    //
    INTEGER nnt = 0;
    read(nin, star), nnt;
    const INTEGER ntypes = 9;
    if (nnt < 1) {
        write(nout, format_9996), " NMA", nnt, 1;
        nnt = 0;
        fatal = true;
    } else if (nnt > ntypes) {
        write(nout, format_9995), " NMA", nnt, ntypes;
        nnt = 0;
        fatal = true;
    }
    INTEGER ntval[ntypes];
    {
        read_loop rloop(cmn, nin, star);
        for (i = 1; i <= nnt; i = i + 1) {
            rloop, ntval[i - 1];
        }
    }
    for (i = 1; i <= nnt; i = i + 1) {
        if (ntval[i - 1] < 0) {
            write(nout, format_9996), "TYPE", ntval[i - 1], 0;
            fatal = true;
        } else if (ntval[i - 1] > ntypes) {
            write(nout, format_9995), "TYPE", ntval[i - 1], ntypes;
            fatal = true;
        }
    }
    if (nnt > 0) {
        {
            write_loop wloop(cmn, nout, format_9993);
            wloop, "TYPE";
            for (i = 1; i <= nnt; i = i + 1) {
                wloop, ntval[i - 1];
            }
        }
    }
    //
    // Read the threshold value for the test ratios.
    //
    REAL thresh = 0.0;
    read(nin, star), thresh;
    write(nout, format_9992), thresh;
    //
    // Read the flag that indicates whether to test the error exits.
    //
    bool tsterr = false;
    read(nin, star), tsterr;
    //
    if (fatal) {
        write(nout, format_9999);
        FEM_STOP(0);
    }
    //
    // Calculate and print the machine dependent constants.
    //
    REAL eps = Rlamch("Underflow threshold");
    write(nout, format_9991), "underflow", eps;
    eps = Rlamch("Overflow threshold");
    write(nout, format_9991), "overflow ", eps;
    eps = Rlamch("Epsilon");
    write(nout, format_9991), "precision", eps;
    write(nout, star);
    //
    // Test the error exit of:
    //
    if (tsterr) {
        Rerrrfp(nout);
    }
    //
    // Test the routines: dpftrf, dpftri, dpftrs (as in Rdrvpo).
    // This also tests the routines: dtfsm, dtftri, dtfttr, dtrttf.
    //
    auto worka_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * nmax));
    REAL *worka = worka_storage.get();
    auto workasav_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * nmax));
    REAL *workasav = workasav_storage.get();
    auto workafac_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * nmax));
    REAL *workafac = workafac_storage.get();
    auto workainv_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * nmax));
    REAL *workainv = workainv_storage.get();
    auto workb_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs));
    REAL *workb = workb_storage.get();
    auto workbsav_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs));
    REAL *workbsav = workbsav_storage.get();
    auto workxact_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs));
    REAL *workxact = workxact_storage.get();
    auto workx_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs));
    REAL *workx = workx_storage.get();
    auto workarf_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, (nmax * (nmax + 1)) / 2));
    REAL *workarf = workarf_storage.get();
    auto workarfinv_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, (nmax * (nmax + 1)) / 2));
    REAL *workarfinv = workarfinv_storage.get();
    REAL d_work_dlatms[3 * nmax];
    REAL d_work_dpot01[nmax];
    auto d_temp_dpot02_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * maxrhs));
    REAL *d_temp_dpot02 = d_temp_dpot02_storage.get();
    auto d_temp_dpot03_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, nmax * nmax));
    REAL *d_temp_dpot03 = d_temp_dpot03_storage.get();
    REAL d_work_dlansy[nmax];
    REAL d_work_dpot02[nmax];
    REAL d_work_dpot03[nmax];
    Rdrvrfp(nout, nn, nval, nns, nsval, nnt, ntval, thresh, worka, workasav, workafac, workainv, workb, workbsav, workxact, workx, workarf, workarfinv, d_work_dlatms, d_work_dpot01, d_temp_dpot02, d_temp_dpot03, d_work_dlansy, d_work_dpot02, d_work_dpot03);
    //
    // Test the routine: dlansf
    //
    Rdrvrf1(nout, nn, nval, thresh, worka, nmax, workarf, d_work_dlansy);
    //
    // Test the conversion routines:
    // dtfttp, dtpttf, dtfttr, dtrttf, dtrttp and dtpttr.
    //
    auto workap_storage = std::make_unique<REAL[]>(std::max<INTEGER>(1, (nmax * (nmax + 1)) / 2));
    REAL *workap = workap_storage.get();
    Rdrvrf2(nout, nn, nval, worka, nmax, workarf, workap, workasav);
    //
    // Test the routine: dtfsm
    //
    Rdrvrf3(nout, nn, nval, thresh, worka, nmax, workarf, workainv, workafac, d_work_dlansy, d_work_dpot03, d_work_dpot01);
    //
    // Test the routine: dsfrk
    //
    Rdrvrf4(nout, nn, nval, thresh, worka, workafac, nmax, workarf, workainv, nmax, d_work_dlansy);
    //
    cmn.io.close(nin);
    REAL s2 = dsecnd();
    write(nout, format_9998);
    write(nout, format_9997), cast2double(s2 - s1);
    //
    // End of Rchkrfp
    //
}

int main(int argc, char const *argv[]) { Rchkrfp(); }
