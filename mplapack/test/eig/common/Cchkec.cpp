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

// Derived from LAPACK routine ZCHKEC.
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

void Cchkec(REAL const thresh, bool const tsterr, INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    //
    static const char *format_9999 = "(' Error in Ctrsyl: RMAX =',d12.3,/,' LMAX = ',i8,' NINFO=',i8,' KNT=',"
                                     "i8)";
    static const char *format_9998 = "(' Error in Ctrexc: RMAX =',d12.3,/,' LMAX = ',i8,' NINFO=',i8,' KNT=',"
                                     "i8)";
    static const char *format_9997 = "(' Error in Ctrsna: RMAX =',3d12.3,/,' LMAX = ',3i8,' NINFO=',3i8,"
                                     "' KNT=',i8)";
    static const char *format_9996 = "(' Error in Ctrsen: RMAX =',3d12.3,/,' LMAX = ',3i8,' NINFO=',3i8,"
                                     "' KNT=',i8)";
    static const char *format_9995 = "(/,1x,'All tests for ',a3,' routines passed the threshold ( ',i6,"
                                     "' tests run)')";
    static const char *format_9994 = "(' Tests of the Nonsymmetric eigenproblem condition',"
                                     "' estimation routines',/,' Ctrsyl, Ctrexc, Ctrsna, Ctrsen',/)";
    static const char *format_9993 = "(' Relative machine precision (EPS) = ',d16.6,/,"
                                     "' Safe minimum (SFMIN)             = ',d16.6,/)";
    static const char *format_9992 = "(' Routines pass computational tests if test ratio is ','less than',f8.2,"
                                     "/,/)";
    static const char *format_9970 = "('Error in Ctrsyl: ',i8,' tests fail the threshold.',/,"
                                     "'Maximum test ratio =',d12.3,' threshold =',d12.3)";
    static const char *format_9971 = "('Error in Ctrsyl3: ',i8,' tests fail the threshold.',/,"
                                     "'Maximum test ratio =',d12.3,' threshold =',d12.3)";
    static const char *format_9972 = "('Ctrsyl and Ctrsyl3 compute an inconsistent scale ','factor in ',i8,"
                                     "' tests.')";
    //
    fem::str<3> path = "Zomplex precision";
    path(2, 3) = "EC";
    REAL eps = Rlamch("P");
    REAL sfmin = Rlamch("S");
    write(nout, format_9994);
    write(nout, format_9993), eps, sfmin;
    write(nout, format_9992), thresh;
    //
    // Test error exits if TSTERR is .TRUE.
    //
    if (tsterr) {
        Cerrec(path, nout);
    }
    //
    bool ok = true;
    REAL rtrsyl[2];
    INTEGER ltrsyl = 0;
    INTEGER ntrsyl = 0;
    INTEGER ktrsyl = 0;
    Cget35(rtrsyl[1 - 1], ltrsyl, ntrsyl, ktrsyl, nin);
    if (rtrsyl[1 - 1] > thresh) {
        ok = false;
        write(nout, format_9999), rtrsyl[1 - 1], ltrsyl, ntrsyl, ktrsyl;
    }
    //
    INTEGER ftrsyl[3];
    INTEGER itrsyl[2];
    INTEGER ktrsyl3 = 0;
    Csyl01(thresh, ftrsyl, rtrsyl, itrsyl, ktrsyl3);
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
    INTEGER ntrexc = 0;
    INTEGER ktrexc = 0;
    Cget36(rtrexc, ltrexc, ntrexc, ktrexc, nin);
    if (rtrexc > thresh || ntrexc > 0) {
        ok = false;
        write(nout, format_9998), rtrexc, ltrexc, ntrexc, ktrexc;
    }
    //
    REAL rtrsna[3];
    INTEGER ltrsna[3];
    INTEGER ntrsna[3];
    INTEGER ktrsna = 0;
    Cget37(rtrsna, ltrsna, ntrsna, ktrsna, nin);
    if (rtrsna[1 - 1] > thresh || rtrsna[2 - 1] > thresh || ntrsna[1 - 1] != 0 || ntrsna[2 - 1] != 0 || ntrsna[3 - 1] != 0) {
        ok = false;
        write(nout, format_9997), rtrsna, ltrsna, ntrsna, ktrsna;
    }
    //
    REAL rtrsen[3];
    INTEGER ltrsen[3];
    INTEGER ntrsen[3];
    INTEGER ktrsen = 0;
    Cget38(rtrsen, ltrsen, ntrsen, ktrsen, nin);
    if (rtrsen[1 - 1] > thresh || rtrsen[2 - 1] > thresh || ntrsen[1 - 1] != 0 || ntrsen[2 - 1] != 0 || ntrsen[3 - 1] != 0) {
        ok = false;
        write(nout, format_9996), rtrsen, ltrsen, ntrsen, ktrsen;
    }
    //
    INTEGER ntests = ktrsyl + ktrsyl3 + ktrexc + ktrsna + ktrsen;
    if (ok) {
        write(nout, format_9995), path, ntests;
    }
    //
    // End of Cchkec
    //
}
