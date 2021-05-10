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

void Rchkec(REAL const thresh, bool const tsterr, INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    char path[3];
    path[0] = 'D';
    path[1] = 'E';
    path[2] = 'C';
    char buf0[1024];
    char buf1[1024];
    char buf2[1024];
    REAL eps = Rlamch("P");
    REAL sfmin = Rlamch("S");
    //
    //     Print header information
    //

    sprintnum_short(buf0, eps);
    sprintnum_short(buf1, sfmin);
    write(nout, "(' Tests of the Nonsymmetric eigenproblem condition estim',"
                "'ation routines',/,' Rlaln2, Rlasy2, Rlanv2, Rlaexc, DTRS',"
                "'YL, Rtrexc, Rtrsna, Rtrsen, Rlaqtr, Rtgexc',/)");
    write(nout, "(' Relative machine precision (EPS) = ',a,/,' Safe ',"
                "'minimum (SFMIN)             = ',a,/)"),
        buf0, buf1;
    sprintnum_short(buf2, thresh);
    write(nout, "(' Routines pass computational tests if test ratio is les','s than',a,"
                "/,/)"),
        buf2;
    //
    //     Test error exits if TSTERR is .TRUE.
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
        sprintnum_short(buf0, rlaln2);
        write(nout, "(' Error in Rlaln2: RMAX =',a,/,' LMAX = ',i8,' N','INFO=',2i8,"
                    "' KNT=',i8)"),
            buf0, llaln2, nlaln2, klaln2;
    }
    //
    REAL rlasy2 = 0.0;
    INTEGER llasy2 = 0;
    INTEGER nlasy2 = 0;
    INTEGER klasy2 = 0;
    Rget32(rlasy2, llasy2, nlasy2, klasy2);
    if (rlasy2 > thresh) {
        ok = false;
        sprintnum_short(buf0, rlasy2);
        write(nout, "(' Error in Rlasy2: RMAX =',a,/,' LMAX = ',i8,' N','INFO=',i8,"
                    "' KNT=',i8)"),
            buf0, llasy2, nlasy2, klasy2;
    }
    //
    REAL rlanv2 = 0.0;
    INTEGER llanv2 = 0;
    INTEGER nlanv2 = 0;
    INTEGER klanv2 = 0;
    Rget33(rlanv2, llanv2, nlanv2, klanv2);
    if (rlanv2 > thresh || nlanv2 != 0) {
        ok = false;
        sprintnum_short(buf0, rlanv2);
        write(nout, "(' Error in Rlanv2: RMAX =',a,/,' LMAX = ',i8,' N','INFO=',i8,"
                    "' KNT=',i8)"),
            buf0, llanv2, nlanv2, klanv2;
    }
    //
    REAL rlaexc = 0.0;
    INTEGER llaexc = 0;
    INTEGER nlaexc[2];
    INTEGER klaexc = 0;
    Rget34(rlaexc, llaexc, nlaexc, klaexc);
    if (rlaexc > thresh || nlaexc[2 - 1] != 0) {
        ok = false;
        sprintnum_short(buf0, rlaexc);
        write(nout, "(' Error in Rlaexc: RMAX =',a,/,' LMAX = ',i8,' N','INFO=',2i8,"
                    "' KNT=',i8)"),
            buf0, llaexc, nlaexc, klaexc;
    }
    //
    REAL rtrsyl = 0.0;
    INTEGER ltrsyl = 0;
    INTEGER ntrsyl = 0;
    INTEGER ktrsyl = 0;
    Rget35(rtrsyl, ltrsyl, ntrsyl, ktrsyl);
    if (rtrsyl > thresh) {
        ok = false;
        sprintnum_short(buf0, rtrsyl);
        write(nout, "(' Error in Rtrsyl: RMAX =',a,/,' LMAX = ',i8,' N','INFO=',i8,"
                    "' KNT=',i8)"),
            buf0, ltrsyl, ntrsyl, ktrsyl;
    }
    //
    REAL rtrexc = 0.0;
    INTEGER ltrexc = 0;
    INTEGER ntrexc[3];
    INTEGER ktrexc = 0;
    Rget36(rtrexc, ltrexc, ntrexc, ktrexc, nin);
    if (rtrexc > thresh || ntrexc[3 - 1] > 0) {
        ok = false;
        sprintnum_short(buf0, rtrexc);
        write(nout, "(' Error in Rtrexc: RMAX =',a,/,' LMAX = ',i8,' N','INFO=',3i8,"
                    "' KNT=',i8)"),
            buf0, ltrexc, ntrexc, ktrexc;
    }
    //
    REAL rtrsna[3];
    INTEGER ltrsna[3];
    INTEGER ntrsna[3];
    INTEGER ktrsna = 0;
    Rget37(rtrsna, ltrsna, ntrsna, ktrsna, nin);
    if (rtrsna[1 - 1] > thresh || rtrsna[2 - 1] > thresh || ntrsna[1 - 1] != 0 || ntrsna[2 - 1] != 0 || ntrsna[3 - 1] != 0) {
        ok = false;
        sprintnum_short(buf0, rtrsna[0]);
        sprintnum_short(buf1, rtrsna[1]);
        sprintnum_short(buf2, rtrsna[2]);
        write(nout, "(' Error in Rtrsna: RMAX =',3a,/,' LMAX = ',3i8,' NINFO=',3i8,"
                    "' KNT=',i8)"),
            buf0, buf1, buf2, ltrsna, ntrsna, ktrsna;
    }
    //
    REAL rtrsen[3];
    INTEGER ltrsen[3];
    INTEGER ntrsen[3];
    INTEGER ktrsen = 0;
    Rget38(rtrsen, ltrsen, ntrsen, ktrsen, nin);
    if (rtrsen[1 - 1] > thresh || rtrsen[2 - 1] > thresh || ntrsen[1 - 1] != 0 || ntrsen[2 - 1] != 0 || ntrsen[3 - 1] != 0) {
        ok = false;
        sprintnum_short(buf0, rtrsen[0]);
        sprintnum_short(buf1, rtrsen[1]);
        sprintnum_short(buf2, rtrsen[2]);
        write(nout, "(' Error in Rtrsen: RMAX =',3a,/,' LMAX = ',3i8,' NINFO=',3i8,"
                    "' KNT=',i8)"),
            buf0, buf1, buf2, ltrsen, ntrsen, ktrsen;
    }
    //
    REAL rlaqtr = 0.0;
    INTEGER llaqtr = 0;
    INTEGER nlaqtr = 0;
    INTEGER klaqtr = 0;
    Rget39(rlaqtr, llaqtr, nlaqtr, klaqtr);
    if (rlaqtr > thresh) {
        ok = false;
        sprintnum_short(buf0, rlaqtr);
        write(nout, "(' Error in Rlaqtr: RMAX =',a,/,' LMAX = ',i8,' N','INFO=',i8,"
                    "' KNT=',i8)"),
            buf0, llaqtr, nlaqtr, klaqtr;
    }
    //
    REAL rtgexc = 0.0;
    INTEGER ltgexc = 0;
    INTEGER ntgexc = 0;
    INTEGER ktgexc = 0;
    Rget40(rtgexc, ltgexc, &ntgexc, ktgexc, nin);
    if (rtgexc > thresh) {
        ok = false;
        sprintnum_short(buf0, rtgexc);
        write(nout, "(' Error in Rtgexc: RMAX =',a,/,' LMAX = ',i8,' N','INFO=',i8,"
                    "' KNT=',i8)"),
            buf0, ltgexc, ntgexc, ktgexc;
    }
    //
    INTEGER ntests = klaln2 + klasy2 + klanv2 + klaexc + ktrsyl + ktrexc + ktrsna + ktrsen + klaqtr + ktgexc;
    if (ok) {
        write(nout, "(/,1x,'All tests for ',a3,' routines passed the thresh','old ( ',i6,"
                    "' tests run)')"),
            path, ntests;
    }
    //
    //     End of Rchkec
    //
}
