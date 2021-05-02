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
#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;
#include <mplapack_lin.h>
#include <mplapack.h>

void Rerrls(const char *path, INTEGER const nunit) {
    common_write write(cmn);
    // COMMON infoc
    INTEGER &infot = cmn.infot;
    INTEGER &nout = cmn.nout;
    bool &ok = cmn.ok;
    bool &lerr = cmn.lerr;
    // COMMON srnamc
    char[32] &srnamt = cmn.srnamt;
    //
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
    //     .. Executable Statements ..
    //
    nout = nunit;
    write(nout, star);
    char[2] c2 = path[(2 - 1) + (3 - 1) * ldpath];
    const INTEGER nmax = 2;
    REAL a[nmax * nmax];
    a[(1 - 1)] = 1.0;
    a[(2 - 1) * lda] = 2.0e+0;
    a[(2 - 1) + (2 - 1) * lda] = 3.0e+0;
    a[(2 - 1)] = 4.0e+0;
    ok = true;
    //
    REAL b[nmax * nmax];
    REAL w[nmax];
    INTEGER info = 0;
    REAL s[nmax];
    REAL rcond = 0.0;
    INTEGER irnk = 0;
    INTEGER ip[nmax];
    if (Mlsamen(2, c2, "LS")) {
        //
        //        Test error exits for the least squares driver routines.
        //
        //        Rgels
        //
        srnamt = "Rgels ";
        infot = 1;
        Rgels("/", 0, 0, 0, a, 1, b, 1, w, 1, info);
        chkxer("Rgels ", infot, nout, lerr, ok);
        infot = 2;
        Rgels("N", -1, 0, 0, a, 1, b, 1, w, 1, info);
        chkxer("Rgels ", infot, nout, lerr, ok);
        infot = 3;
        Rgels("N", 0, -1, 0, a, 1, b, 1, w, 1, info);
        chkxer("Rgels ", infot, nout, lerr, ok);
        infot = 4;
        Rgels("N", 0, 0, -1, a, 1, b, 1, w, 1, info);
        chkxer("Rgels ", infot, nout, lerr, ok);
        infot = 6;
        Rgels("N", 2, 0, 0, a, 1, b, 2, w, 2, info);
        chkxer("Rgels ", infot, nout, lerr, ok);
        infot = 8;
        Rgels("N", 2, 0, 0, a, 2, b, 1, w, 2, info);
        chkxer("Rgels ", infot, nout, lerr, ok);
        infot = 10;
        Rgels("N", 1, 1, 0, a, 1, b, 1, w, 1, info);
        chkxer("Rgels ", infot, nout, lerr, ok);
        //
        //        Rgelss
        //
        srnamt = "Rgelss";
        infot = 1;
        Rgelss(-1, 0, 0, a, 1, b, 1, s, rcond, irnk, w, 1, info);
        chkxer("Rgelss", infot, nout, lerr, ok);
        infot = 2;
        Rgelss(0, -1, 0, a, 1, b, 1, s, rcond, irnk, w, 1, info);
        chkxer("Rgelss", infot, nout, lerr, ok);
        infot = 3;
        Rgelss(0, 0, -1, a, 1, b, 1, s, rcond, irnk, w, 1, info);
        chkxer("Rgelss", infot, nout, lerr, ok);
        infot = 5;
        Rgelss(2, 0, 0, a, 1, b, 2, s, rcond, irnk, w, 2, info);
        chkxer("Rgelss", infot, nout, lerr, ok);
        infot = 7;
        Rgelss(2, 0, 0, a, 2, b, 1, s, rcond, irnk, w, 2, info);
        chkxer("Rgelss", infot, nout, lerr, ok);
        //
        //        Rgelsy
        //
        srnamt = "Rgelsy";
        infot = 1;
        Rgelsy(-1, 0, 0, a, 1, b, 1, ip, rcond, irnk, w, 10, info);
        chkxer("Rgelsy", infot, nout, lerr, ok);
        infot = 2;
        Rgelsy(0, -1, 0, a, 1, b, 1, ip, rcond, irnk, w, 10, info);
        chkxer("Rgelsy", infot, nout, lerr, ok);
        infot = 3;
        Rgelsy(0, 0, -1, a, 1, b, 1, ip, rcond, irnk, w, 10, info);
        chkxer("Rgelsy", infot, nout, lerr, ok);
        infot = 5;
        Rgelsy(2, 0, 0, a, 1, b, 2, ip, rcond, irnk, w, 10, info);
        chkxer("Rgelsy", infot, nout, lerr, ok);
        infot = 7;
        Rgelsy(2, 0, 0, a, 2, b, 1, ip, rcond, irnk, w, 10, info);
        chkxer("Rgelsy", infot, nout, lerr, ok);
        infot = 12;
        Rgelsy(2, 2, 1, a, 2, b, 2, ip, rcond, irnk, w, 1, info);
        chkxer("Rgelsy", infot, nout, lerr, ok);
        //
        //        Rgelsd
        //
        srnamt = "Rgelsd";
        infot = 1;
        Rgelsd(-1, 0, 0, a, 1, b, 1, s, rcond, irnk, w, 10, ip, info);
        chkxer("Rgelsd", infot, nout, lerr, ok);
        infot = 2;
        Rgelsd(0, -1, 0, a, 1, b, 1, s, rcond, irnk, w, 10, ip, info);
        chkxer("Rgelsd", infot, nout, lerr, ok);
        infot = 3;
        Rgelsd(0, 0, -1, a, 1, b, 1, s, rcond, irnk, w, 10, ip, info);
        chkxer("Rgelsd", infot, nout, lerr, ok);
        infot = 5;
        Rgelsd(2, 0, 0, a, 1, b, 2, s, rcond, irnk, w, 10, ip, info);
        chkxer("Rgelsd", infot, nout, lerr, ok);
        infot = 7;
        Rgelsd(2, 0, 0, a, 2, b, 1, s, rcond, irnk, w, 10, ip, info);
        chkxer("Rgelsd", infot, nout, lerr, ok);
        infot = 12;
        Rgelsd(2, 2, 1, a, 2, b, 2, s, rcond, irnk, w, 1, ip, info);
        chkxer("Rgelsd", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrls
    //
}
