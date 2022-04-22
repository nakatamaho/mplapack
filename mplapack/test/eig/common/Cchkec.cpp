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

void Cchkec(REAL const thresh, bool const tsterr, INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    //
    char path[3];
    char buf0[1024];
    char buf1[1024];
    char buf2[1024];
    path[0] = 'C';
    path[1] = 'E';
    path[2] = 'C';
    REAL eps = Rlamch("P");
    REAL sfmin = Rlamch("S");
    //
    //     Print header information
    //
    sprintnum_short(buf1, eps);
    sprintnum_short(buf2, sfmin);
    write(nout, "(' Tests of the Nonsymmetric eigenproblem condition',"
                "' estimation routines',/,' Ctrsyl, Ctrexc, Ctrsna, Ctrsen',/)");
    write(nout, "(' Relative machine precision (EPS) = ',a,/,"
                "' Safe minimum (SFMIN)             = ',a,/)"),
        buf1, buf2;
    sprintnum_short(buf1, thresh);
    write(nout, "(' Routines pass computational tests if test ratio is ','less than',a,"
                "/,/)"),
        buf1;
    //
    //     Test error exits if TSTERR is .TRUE.
    //
    if (tsterr) {
        Cerrec(path, nout);
    }
    //
    ok = true;
    REAL rtrsyl = 0.0;
    INTEGER ltrsyl = 0;
    INTEGER ntrsyl = 0;
    INTEGER ktrsyl = 0;
    Cget35(rtrsyl, ltrsyl, ntrsyl, ktrsyl, nin);
    if (rtrsyl > thresh) {
        ok = false;
        sprintnum_short(buf1, rtrsyl);
        write(nout, "(' Error in Ctrsyl: RMAX =',a,/,' LMAX = ',i8,' NINFO=',i8,' KNT=',"
                    "i8)"),
            buf1, ltrsyl, ntrsyl, ktrsyl;
    }
    //
    REAL rtrexc = 0.0;
    INTEGER ltrexc = 0;
    INTEGER ntrexc = 0;
    INTEGER ktrexc = 0;
    Cget36(rtrexc, ltrexc, ntrexc, ktrexc, nin);
    if (rtrexc > thresh || ntrexc > 0) {
        ok = false;
        sprintnum_short(buf1, rtrexc);
        write(nout, "(' Error in Ctrexc: RMAX =',a,/,' LMAX = ',i8,' NINFO=',i8,' KNT=',"
                    "i8)"),
            buf1, ltrexc, ntrexc, ktrexc;
    }
    //
    REAL rtrsna[3];
    INTEGER ltrsna[3];
    INTEGER ntrsna[3];
    INTEGER ktrsna = 0;
    Cget37(rtrsna, ltrsna, ntrsna, ktrsna, nin);
    if (rtrsna[1 - 1] > thresh || rtrsna[2 - 1] > thresh || ntrsna[1 - 1] != 0 || ntrsna[2 - 1] != 0 || ntrsna[3 - 1] != 0) {
        ok = false;
        sprintnum_short(buf0, rtrsna[0]);
        sprintnum_short(buf1, rtrsna[1]);
        sprintnum_short(buf2, rtrsna[2]);
        write(nout, "(' Error in Ctrsna: RMAX =',3a,/,' LMAX = ',3i8,' NINFO=',3i8,"
                    "' KNT=',i8)"),
            buf0, buf1, buf2, ltrsna, ntrsna, ktrsna;
    }
    //
    REAL rtrsen[3];
    INTEGER ltrsen[3];
    INTEGER ntrsen[3];
    INTEGER ktrsen = 0;
    Cget38(rtrsen, ltrsen, ntrsen, ktrsen, nin);
    if (rtrsen[1 - 1] > thresh || rtrsen[2 - 1] > thresh || ntrsen[1 - 1] != 0 || ntrsen[2 - 1] != 0 || ntrsen[3 - 1] != 0) {
        ok = false;
        sprintnum_short(buf0, rtrsen[0]);
        sprintnum_short(buf1, rtrsen[1]);
        sprintnum_short(buf2, rtrsen[2]);
        write(nout, "(' Error in Ctrsen: RMAX =',3a,/,' LMAX = ',3i8,' NINFO=',3i8,"
                    "' KNT=',i8)"),
            buf0, buf1, buf2, ltrsen, ntrsen, ktrsen;
    }
    //
    INTEGER ntests = ktrsyl + ktrexc + ktrsna + ktrsen;
    if (ok) {
        write(nout, "(/,1x,'All tests for ',a3,' routines passed the threshold ( ',i6,"
                    "' tests run)')"),
            path, ntests;
    }
    //
    //     End of Cchkec
    //
}
