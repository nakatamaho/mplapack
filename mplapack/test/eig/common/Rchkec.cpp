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

// Derived from LAPACK routine DCHKEC.
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

void Rchkec(REAL const thresh, bool const tsterr, INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    static const char *format_9999 = "(' Error in Rlaln2: RMAX =',d12.3,/,' LMAX = ',i8,' N','INFO=',2i8,"
                                     "' KNT=',i8)";
    static const char *format_9998 = "(' Error in Rlasy2: RMAX =',d12.3,/,' LMAX = ',i8,' N','INFO=',i8,"
                                     "' KNT=',i8)";
    static const char *format_9997 = "(' Error in Rlanv2: RMAX =',d12.3,/,' LMAX = ',i8,' N','INFO=',i8,"
                                     "' KNT=',i8)";
    static const char *format_9996 = "(' Error in Rlaexc: RMAX =',d12.3,/,' LMAX = ',i8,' N','INFO=',2i8,"
                                     "' KNT=',i8)";
    static const char *format_9995 = "(' Error in Rtrsyl: RMAX =',d12.3,/,' LMAX = ',i8,' N','INFO=',i8,"
                                     "' KNT=',i8)";
    static const char *format_9994 = "(' Error in Rtrexc: RMAX =',d12.3,/,' LMAX = ',i8,' N','INFO=',3i8,"
                                     "' KNT=',i8)";
    static const char *format_9993 = "(' Error in Rtrsna: RMAX =',3d12.3,/,' LMAX = ',3i8,' NINFO=',3i8,"
                                     "' KNT=',i8)";
    static const char *format_9992 = "(' Error in Rtrsen: RMAX =',3d12.3,/,' LMAX = ',3i8,' NINFO=',3i8,"
                                     "' KNT=',i8)";
    static const char *format_9991 = "(' Error in Rlaqtr: RMAX =',d12.3,/,' LMAX = ',i8,' N','INFO=',i8,"
                                     "' KNT=',i8)";
    static const char *format_9990 = "(/,1x,'All tests for ',a3,' routines passed the thresh','old ( ',i6,"
                                     "' tests run)')";
    static const char *format_9989 = "(' Tests of the Nonsymmetric eigenproblem condition estim',"
                                     "'ation routines',/,' Rlaln2, Rlasy2, Rlanv2, Rlaexc, DTRS',"
                                     "'YL, Rtrexc, Rtrsna, Rtrsen, Rlaqtr, Rtgexc',/)";
    static const char *format_9988 = "(' Relative machine precision (EPS) = ',d16.6,/,' Safe ',"
                                     "'minimum (SFMIN)             = ',d16.6,/)";
    static const char *format_9987 = "(' Routines pass computational tests if test ratio is les','s than',f8.2,"
                                     "/,/)";
    static const char *format_9986 = "(' Error in Rtgexc: RMAX =',d12.3,/,' LMAX = ',i8,' N','INFO=',2i8,"
                                     "' KNT=',i8)";
    static const char *format_9972 = "('Rtrsyl and Rtrsyl3 compute an inconsistent result ','factor in ',i8,"
                                     "' tests.')";
    static const char *format_9971 = "('Error in Rtrsyl3: ',i8,' tests fail the threshold.',/,"
                                     "'Maximum test ratio =',d12.3,' threshold =',d12.3)";
    static const char *format_9970 = "('Error in Rtrsyl: ',i8,' tests fail the threshold.',/,"
                                     "'Maximum test ratio =',d12.3,' threshold =',d12.3)";
    //
    fem::str<3> path = "Double precision";
    path(2, 3) = "EC";
    REAL eps = Rlamch("P");
    REAL sfmin = Rlamch("S");
    //
    // Print header information
    //
    write(nout, format_9989);
    write(nout, format_9988), eps, sfmin;
    write(nout, format_9987), thresh;
    //
    // Test error exits if TSTERR is .TRUE.
    //
    if (tsterr) {
        Rerrec(path, nout);
    }
    //
    bool ok = true;
    REAL rlaln2 = 0.0;
    INTEGER llaln2 = 0;
    INTEGER nlaln2[2];
    INTEGER klaln2 = 0;
    Rget31(rlaln2, llaln2, nlaln2, klaln2);
    if (rlaln2 > thresh || nlaln2[1 - 1] != 0) {
        ok = false;
        write(nout, format_9999), rlaln2, llaln2, nlaln2, klaln2;
    }
    //
    REAL rlasy2 = 0.0;
    INTEGER llasy2 = 0;
    INTEGER nlasy2 = 0;
    INTEGER klasy2 = 0;
    Rget32(rlasy2, llasy2, nlasy2, klasy2);
    if (rlasy2 > thresh) {
        ok = false;
        write(nout, format_9998), rlasy2, llasy2, nlasy2, klasy2;
    }
    //
    REAL rlanv2 = 0.0;
    INTEGER llanv2 = 0;
    INTEGER nlanv2 = 0;
    INTEGER klanv2 = 0;
    Rget33(rlanv2, llanv2, nlanv2, klanv2);
    if (rlanv2 > thresh || nlanv2 != 0) {
        ok = false;
        write(nout, format_9997), rlanv2, llanv2, nlanv2, klanv2;
    }
    //
    REAL rlaexc = 0.0;
    INTEGER llaexc = 0;
    INTEGER nlaexc[2];
    INTEGER klaexc = 0;
    Rget34(rlaexc, llaexc, nlaexc, klaexc);
    if (rlaexc > thresh || nlaexc[2 - 1] != 0) {
        ok = false;
        write(nout, format_9996), rlaexc, llaexc, nlaexc, klaexc;
    }
    //
    REAL rtrsyl[2];
    INTEGER ltrsyl = 0;
    INTEGER ntrsyl = 0;
    INTEGER ktrsyl = 0;
    Rget35(rtrsyl[1 - 1], ltrsyl, ntrsyl, ktrsyl);
    if (rtrsyl[1 - 1] > thresh) {
        ok = false;
        write(nout, format_9995), rtrsyl[1 - 1], ltrsyl, ntrsyl, ktrsyl;
    }
    //
    INTEGER ftrsyl[3];
    INTEGER itrsyl[2];
    INTEGER ktrsyl3 = 0;
    Rsyl01(thresh, ftrsyl, rtrsyl, itrsyl, ktrsyl3);
    if (ftrsyl[1 - 1] > 0) {
        ok = false;
        write(nout, format_9970), ftrsyl[1 - 1], rtrsyl[1 - 1], thresh;
    }
    if (ftrsyl[2 - 1] > 0) {
        ok = false;
        write(nout, format_9971), ftrsyl[2 - 1], rtrsyl[2 - 1], thresh;
    }
    if (ftrsyl[3 - 1] > 0) {
        ok = false;
        write(nout, format_9972), ftrsyl[3 - 1];
    }
    //
    REAL rtrexc = 0.0;
    INTEGER ltrexc = 0;
    INTEGER ntrexc[3];
    INTEGER ktrexc = 0;
    Rget36(rtrexc, ltrexc, ntrexc, ktrexc, nin);
    if (rtrexc > thresh || ntrexc[3 - 1] > 0) {
        ok = false;
        write(nout, format_9994), rtrexc, ltrexc, ntrexc, ktrexc;
    }
    //
    REAL rtrsna[3];
    INTEGER ltrsna[3];
    INTEGER ntrsna[3];
    INTEGER ktrsna = 0;
    Rget37(rtrsna, ltrsna, ntrsna, ktrsna, nin);
    if (rtrsna[1 - 1] > thresh || rtrsna[2 - 1] > thresh || ntrsna[1 - 1] != 0 || ntrsna[2 - 1] != 0 || ntrsna[3 - 1] != 0) {
        ok = false;
        write(nout, format_9993), rtrsna, ltrsna, ntrsna, ktrsna;
    }
    //
    REAL rtrsen[3];
    INTEGER ltrsen[3];
    INTEGER ntrsen[3];
    INTEGER ktrsen = 0;
    Rget38(rtrsen, ltrsen, ntrsen, ktrsen, nin);
    if (rtrsen[1 - 1] > thresh || rtrsen[2 - 1] > thresh || ntrsen[1 - 1] != 0 || ntrsen[2 - 1] != 0 || ntrsen[3 - 1] != 0) {
        ok = false;
        write(nout, format_9992), rtrsen, ltrsen, ntrsen, ktrsen;
    }
    //
    REAL rlaqtr = 0.0;
    INTEGER llaqtr = 0;
    INTEGER nlaqtr = 0;
    INTEGER klaqtr = 0;
    Rget39(rlaqtr, llaqtr, nlaqtr, klaqtr);
    if (rlaqtr > thresh) {
        ok = false;
        write(nout, format_9991), rlaqtr, llaqtr, nlaqtr, klaqtr;
    }
    //
    REAL rtgexc = 0.0;
    INTEGER ltgexc = 0;
    INTEGER ntgexc[2];
    INTEGER ktgexc = 0;
    Rget40(rtgexc, ltgexc, ntgexc, ktgexc, nin);
    if (rtgexc > thresh) {
        ok = false;
        write(nout, format_9986), rtgexc, ltgexc, ntgexc, ktgexc;
    }
    //
    INTEGER ntests = klaln2 + klasy2 + klanv2 + klaexc + ktrsyl + ktrexc + ktrsna + ktrsen + klaqtr + ktgexc;
    if (ok) {
        write(nout, format_9990), path, ntests;
    }
    //
    // End of Rchkec
    //
}
