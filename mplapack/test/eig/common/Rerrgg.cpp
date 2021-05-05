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

void Rerrgg(const char *path, INTEGER const nunit) {
    common_write write(cmn);
    // COMMON infoc
    INTEGER &infot = cmn.infot;
    INTEGER &nout = cmn.nout;
    bool &ok = cmn.ok;
    bool &lerr = cmn.lerr;
    // COMMON srnamc
    char &srnamt = cmn.srnamt;
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
    char c2[2] = path[(2 - 1) + (3 - 1) * ldpath];
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 3;
    bool sel[nmax];
    INTEGER i = 0;
    const REAL zero = 0.0;
    REAL a[nmax * nmax];
    REAL b[nmax * nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        sel[j - 1] = true;
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
            b[(i - 1) + (j - 1) * ldb] = zero;
        }
    }
    const REAL one = 1.0;
    for (i = 1; i <= nmax; i = i + 1) {
        a[(i - 1) + (i - 1) * lda] = one;
        b[(i - 1) + (i - 1) * ldb] = one;
    }
    ok = true;
    REAL tola = 1.0;
    REAL tolb = 1.0;
    INTEGER ifst = 1;
    INTEGER ilst = 1;
    INTEGER nt = 0;
    INTEGER lwork = 1;
    //
    //     Test error exits for the GG path.
    //
    REAL q[nmax * nmax];
    REAL z[nmax * nmax];
    INTEGER info = 0;
    const INTEGER lw = 6 * nmax;
    REAL w[lw];
    REAL r1[nmax];
    REAL r2[nmax];
    REAL r3[nmax];
    INTEGER m = 0;
    INTEGER dummyk = 0;
    INTEGER dummyl = 0;
    REAL u[nmax * nmax];
    REAL v[nmax * nmax];
    INTEGER idum[nmax];
    INTEGER iw[nmax];
    REAL tau[nmax];
    INTEGER ncycle = 0;
    INTEGER sdim = 0;
    bool bw[nmax];
    REAL rce[2];
    REAL rcv[2];
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL ls[nmax];
    REAL rs[nmax];
    REAL anrm = 0.0;
    REAL bnrm = 0.0;
    REAL scale = 0.0;
    REAL dif = 0.0;
    if (Mlsamen(2, c2, "GG")) {
        //
        //        Rgghrd
        //
        srnamt = "Rgghrd";
        infot = 1;
        Rgghrd("/", "N", 0, 1, 0, a, 1, b, 1, q, 1, z, 1, info);
        chkxer("Rgghrd", infot, nout, lerr, ok);
        infot = 2;
        Rgghrd("N", "/", 0, 1, 0, a, 1, b, 1, q, 1, z, 1, info);
        chkxer("Rgghrd", infot, nout, lerr, ok);
        infot = 3;
        Rgghrd("N", "N", -1, 0, 0, a, 1, b, 1, q, 1, z, 1, info);
        chkxer("Rgghrd", infot, nout, lerr, ok);
        infot = 4;
        Rgghrd("N", "N", 0, 0, 0, a, 1, b, 1, q, 1, z, 1, info);
        chkxer("Rgghrd", infot, nout, lerr, ok);
        infot = 5;
        Rgghrd("N", "N", 0, 1, 1, a, 1, b, 1, q, 1, z, 1, info);
        chkxer("Rgghrd", infot, nout, lerr, ok);
        infot = 7;
        Rgghrd("N", "N", 2, 1, 1, a, 1, b, 2, q, 1, z, 1, info);
        chkxer("Rgghrd", infot, nout, lerr, ok);
        infot = 9;
        Rgghrd("N", "N", 2, 1, 1, a, 2, b, 1, q, 1, z, 1, info);
        chkxer("Rgghrd", infot, nout, lerr, ok);
        infot = 11;
        Rgghrd("V", "N", 2, 1, 1, a, 2, b, 2, q, 1, z, 1, info);
        chkxer("Rgghrd", infot, nout, lerr, ok);
        infot = 13;
        Rgghrd("N", "V", 2, 1, 1, a, 2, b, 2, q, 1, z, 1, info);
        chkxer("Rgghrd", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Rgghd3
        //
        srnamt = "Rgghd3";
        infot = 1;
        Rgghd3("/", "N", 0, 1, 0, a, 1, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Rgghd3", infot, nout, lerr, ok);
        infot = 2;
        Rgghd3("N", "/", 0, 1, 0, a, 1, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Rgghd3", infot, nout, lerr, ok);
        infot = 3;
        Rgghd3("N", "N", -1, 0, 0, a, 1, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Rgghd3", infot, nout, lerr, ok);
        infot = 4;
        Rgghd3("N", "N", 0, 0, 0, a, 1, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Rgghd3", infot, nout, lerr, ok);
        infot = 5;
        Rgghd3("N", "N", 0, 1, 1, a, 1, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Rgghd3", infot, nout, lerr, ok);
        infot = 7;
        Rgghd3("N", "N", 2, 1, 1, a, 1, b, 2, q, 1, z, 1, w, lw, info);
        chkxer("Rgghd3", infot, nout, lerr, ok);
        infot = 9;
        Rgghd3("N", "N", 2, 1, 1, a, 2, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Rgghd3", infot, nout, lerr, ok);
        infot = 11;
        Rgghd3("V", "N", 2, 1, 1, a, 2, b, 2, q, 1, z, 1, w, lw, info);
        chkxer("Rgghd3", infot, nout, lerr, ok);
        infot = 13;
        Rgghd3("N", "V", 2, 1, 1, a, 2, b, 2, q, 1, z, 1, w, lw, info);
        chkxer("Rgghd3", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Rhgeqz
        //
        srnamt = "Rhgeqz";
        infot = 1;
        Rhgeqz("/", "N", "N", 0, 1, 0, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, w, lw, info);
        chkxer("Rhgeqz", infot, nout, lerr, ok);
        infot = 2;
        Rhgeqz("E", "/", "N", 0, 1, 0, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, w, lw, info);
        chkxer("Rhgeqz", infot, nout, lerr, ok);
        infot = 3;
        Rhgeqz("E", "N", "/", 0, 1, 0, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, w, lw, info);
        chkxer("Rhgeqz", infot, nout, lerr, ok);
        infot = 4;
        Rhgeqz("E", "N", "N", -1, 0, 0, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, w, lw, info);
        chkxer("Rhgeqz", infot, nout, lerr, ok);
        infot = 5;
        Rhgeqz("E", "N", "N", 0, 0, 0, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, w, lw, info);
        chkxer("Rhgeqz", infot, nout, lerr, ok);
        infot = 6;
        Rhgeqz("E", "N", "N", 0, 1, 1, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, w, lw, info);
        chkxer("Rhgeqz", infot, nout, lerr, ok);
        infot = 8;
        Rhgeqz("E", "N", "N", 2, 1, 1, a, 1, b, 2, r1, r2, r3, q, 1, z, 1, w, lw, info);
        chkxer("Rhgeqz", infot, nout, lerr, ok);
        infot = 10;
        Rhgeqz("E", "N", "N", 2, 1, 1, a, 2, b, 1, r1, r2, r3, q, 1, z, 1, w, lw, info);
        chkxer("Rhgeqz", infot, nout, lerr, ok);
        infot = 15;
        Rhgeqz("E", "V", "N", 2, 1, 1, a, 2, b, 2, r1, r2, r3, q, 1, z, 1, w, lw, info);
        chkxer("Rhgeqz", infot, nout, lerr, ok);
        infot = 17;
        Rhgeqz("E", "N", "V", 2, 1, 1, a, 2, b, 2, r1, r2, r3, q, 1, z, 1, w, lw, info);
        chkxer("Rhgeqz", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Rtgevc
        //
        srnamt = "Rtgevc";
        infot = 1;
        Rtgevc("/", "A", sel, 0, a, 1, b, 1, q, 1, z, 1, 0, m, w, info);
        chkxer("Rtgevc", infot, nout, lerr, ok);
        infot = 2;
        Rtgevc("R", "/", sel, 0, a, 1, b, 1, q, 1, z, 1, 0, m, w, info);
        chkxer("Rtgevc", infot, nout, lerr, ok);
        infot = 4;
        Rtgevc("R", "A", sel, -1, a, 1, b, 1, q, 1, z, 1, 0, m, w, info);
        chkxer("Rtgevc", infot, nout, lerr, ok);
        infot = 6;
        Rtgevc("R", "A", sel, 2, a, 1, b, 2, q, 1, z, 2, 0, m, w, info);
        chkxer("Rtgevc", infot, nout, lerr, ok);
        infot = 8;
        Rtgevc("R", "A", sel, 2, a, 2, b, 1, q, 1, z, 2, 0, m, w, info);
        chkxer("Rtgevc", infot, nout, lerr, ok);
        infot = 10;
        Rtgevc("L", "A", sel, 2, a, 2, b, 2, q, 1, z, 1, 0, m, w, info);
        chkxer("Rtgevc", infot, nout, lerr, ok);
        infot = 12;
        Rtgevc("R", "A", sel, 2, a, 2, b, 2, q, 1, z, 1, 0, m, w, info);
        chkxer("Rtgevc", infot, nout, lerr, ok);
        infot = 13;
        Rtgevc("R", "A", sel, 2, a, 2, b, 2, q, 1, z, 2, 1, m, w, info);
        chkxer("Rtgevc", infot, nout, lerr, ok);
        nt += 8;
        //
        //     Test error exits for the GSV path.
        //
    } else if (Mlsamen(3, path, "GSV")) {
        //
        //        Rggsvd3
        //
        srnamt = "Rggsvd3";
        infot = 1;
        Rggsvd3("/", "N", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        infot = 2;
        Rggsvd3("N", "/", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        infot = 3;
        Rggsvd3("N", "N", "/", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        infot = 4;
        Rggsvd3("N", "N", "N", -1, 0, 0, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        infot = 5;
        Rggsvd3("N", "N", "N", 0, -1, 0, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        infot = 6;
        Rggsvd3("N", "N", "N", 0, 0, -1, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        infot = 10;
        Rggsvd3("N", "N", "N", 2, 1, 1, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        infot = 12;
        Rggsvd3("N", "N", "N", 1, 1, 2, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        infot = 16;
        Rggsvd3("U", "N", "N", 2, 2, 2, dummyk, dummyl, a, 2, b, 2, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        infot = 18;
        Rggsvd3("N", "V", "N", 1, 1, 2, dummyk, dummyl, a, 1, b, 2, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        infot = 20;
        Rggsvd3("N", "N", "Q", 1, 2, 1, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, idum, info);
        chkxer("Rggsvd3", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Rggsvp3
        //
        srnamt = "Rggsvp3";
        infot = 1;
        Rggsvp3("/", "N", "N", 0, 0, 0, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        infot = 2;
        Rggsvp3("N", "/", "N", 0, 0, 0, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        infot = 3;
        Rggsvp3("N", "N", "/", 0, 0, 0, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        infot = 4;
        Rggsvp3("N", "N", "N", -1, 0, 0, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        infot = 5;
        Rggsvp3("N", "N", "N", 0, -1, 0, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        infot = 6;
        Rggsvp3("N", "N", "N", 0, 0, -1, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        infot = 8;
        Rggsvp3("N", "N", "N", 2, 1, 1, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        infot = 10;
        Rggsvp3("N", "N", "N", 1, 2, 1, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        infot = 16;
        Rggsvp3("U", "N", "N", 2, 2, 2, a, 2, b, 2, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        infot = 18;
        Rggsvp3("N", "V", "N", 1, 2, 1, a, 1, b, 2, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        infot = 20;
        Rggsvp3("N", "N", "Q", 1, 1, 2, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, tau, w, lwork, info);
        chkxer("Rggsvp3", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Rtgsja
        //
        srnamt = "Rtgsja";
        infot = 1;
        Rtgsja("/", "N", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        infot = 2;
        Rtgsja("N", "/", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        infot = 3;
        Rtgsja("N", "N", "/", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        infot = 4;
        Rtgsja("N", "N", "N", -1, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        infot = 5;
        Rtgsja("N", "N", "N", 0, -1, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        infot = 6;
        Rtgsja("N", "N", "N", 0, 0, -1, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        infot = 10;
        Rtgsja("N", "N", "N", 0, 0, 0, dummyk, dummyl, a, 0, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        infot = 12;
        Rtgsja("N", "N", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 0, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        infot = 18;
        Rtgsja("U", "N", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 0, v, 1, q, 1, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        infot = 20;
        Rtgsja("N", "V", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 0, q, 1, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        infot = 22;
        Rtgsja("N", "N", "Q", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 0, w, ncycle, info);
        chkxer("Rtgsja", infot, nout, lerr, ok);
        nt += 11;
        //
        //     Test error exits for the GLM path.
        //
    } else if (Mlsamen(3, path, "GLM")) {
        //
        //        Rggglm
        //
        srnamt = "Rggglm";
        infot = 1;
        Rggglm(-1, 0, 0, a, 1, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rggglm", infot, nout, lerr, ok);
        infot = 2;
        Rggglm(0, -1, 0, a, 1, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rggglm", infot, nout, lerr, ok);
        infot = 2;
        Rggglm(0, 1, 0, a, 1, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rggglm", infot, nout, lerr, ok);
        infot = 3;
        Rggglm(0, 0, -1, a, 1, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rggglm", infot, nout, lerr, ok);
        infot = 3;
        Rggglm(1, 0, 0, a, 1, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rggglm", infot, nout, lerr, ok);
        infot = 5;
        Rggglm(0, 0, 0, a, 0, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rggglm", infot, nout, lerr, ok);
        infot = 7;
        Rggglm(0, 0, 0, a, 1, b, 0, r1, r2, r3, w, lw, info);
        chkxer("Rggglm", infot, nout, lerr, ok);
        infot = 12;
        Rggglm(1, 1, 1, a, 1, b, 1, r1, r2, r3, w, 1, info);
        chkxer("Rggglm", infot, nout, lerr, ok);
        nt += 8;
        //
        //     Test error exits for the LSE path.
        //
    } else if (Mlsamen(3, path, "LSE")) {
        //
        //        Rgglse
        //
        srnamt = "Rgglse";
        infot = 1;
        Rgglse(-1, 0, 0, a, 1, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rgglse", infot, nout, lerr, ok);
        infot = 2;
        Rgglse(0, -1, 0, a, 1, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rgglse", infot, nout, lerr, ok);
        infot = 3;
        Rgglse(0, 0, -1, a, 1, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rgglse", infot, nout, lerr, ok);
        infot = 3;
        Rgglse(0, 0, 1, a, 1, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rgglse", infot, nout, lerr, ok);
        infot = 3;
        Rgglse(0, 1, 0, a, 1, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rgglse", infot, nout, lerr, ok);
        infot = 5;
        Rgglse(0, 0, 0, a, 0, b, 1, r1, r2, r3, w, lw, info);
        chkxer("Rgglse", infot, nout, lerr, ok);
        infot = 7;
        Rgglse(0, 0, 0, a, 1, b, 0, r1, r2, r3, w, lw, info);
        chkxer("Rgglse", infot, nout, lerr, ok);
        infot = 12;
        Rgglse(1, 1, 1, a, 1, b, 1, r1, r2, r3, w, 1, info);
        chkxer("Rgglse", infot, nout, lerr, ok);
        nt += 8;
        //
        //     Test error exits for the CSD path.
        //
    } else if (Mlsamen(3, path, "CSD")) {
        //
        //        Rorcsd
        //
        srnamt = "Rorcsd";
        infot = 7;
        Rorcsd("Y", "Y", "Y", "Y", "N", "N", -1, 0, 0, a, 1, a, 1, a, 1, a, 1, a, a, 1, a, 1, a, 1, a, 1, w, lw, iw, info);
        chkxer("Rorcsd", infot, nout, lerr, ok);
        infot = 8;
        Rorcsd("Y", "Y", "Y", "Y", "N", "N", 1, -1, 0, a, 1, a, 1, a, 1, a, 1, a, a, 1, a, 1, a, 1, a, 1, w, lw, iw, info);
        chkxer("Rorcsd", infot, nout, lerr, ok);
        infot = 9;
        Rorcsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, -1, a, 1, a, 1, a, 1, a, 1, a, a, 1, a, 1, a, 1, a, 1, w, lw, iw, info);
        chkxer("Rorcsd", infot, nout, lerr, ok);
        infot = 11;
        Rorcsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, 1, a, -1, a, 1, a, 1, a, 1, a, a, 1, a, 1, a, 1, a, 1, w, lw, iw, info);
        chkxer("Rorcsd", infot, nout, lerr, ok);
        infot = 20;
        Rorcsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, 1, a, 1, a, 1, a, 1, a, 1, a, a, -1, a, 1, a, 1, a, 1, w, lw, iw, info);
        chkxer("Rorcsd", infot, nout, lerr, ok);
        infot = 22;
        Rorcsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, 1, a, 1, a, 1, a, 1, a, 1, a, a, 1, a, -1, a, 1, a, 1, w, lw, iw, info);
        chkxer("Rorcsd", infot, nout, lerr, ok);
        infot = 24;
        Rorcsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, 1, a, 1, a, 1, a, 1, a, 1, a, a, 1, a, 1, a, -1, a, 1, w, lw, iw, info);
        chkxer("Rorcsd", infot, nout, lerr, ok);
        infot = 26;
        Rorcsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, 1, a, 1, a, 1, a, 1, a, 1, a, a, 1, a, 1, a, 1, a, -1, w, lw, iw, info);
        chkxer("Rorcsd", infot, nout, lerr, ok);
        nt += 8;
        //
        //     Test error exits for the GQR path.
        //
    } else if (Mlsamen(3, path, "GQR")) {
        //
        //        Rggqrf
        //
        srnamt = "Rggqrf";
        infot = 1;
        Rggqrf(-1, 0, 0, a, 1, r1, b, 1, r2, w, lw, info);
        chkxer("Rggqrf", infot, nout, lerr, ok);
        infot = 2;
        Rggqrf(0, -1, 0, a, 1, r1, b, 1, r2, w, lw, info);
        chkxer("Rggqrf", infot, nout, lerr, ok);
        infot = 3;
        Rggqrf(0, 0, -1, a, 1, r1, b, 1, r2, w, lw, info);
        chkxer("Rggqrf", infot, nout, lerr, ok);
        infot = 5;
        Rggqrf(0, 0, 0, a, 0, r1, b, 1, r2, w, lw, info);
        chkxer("Rggqrf", infot, nout, lerr, ok);
        infot = 8;
        Rggqrf(0, 0, 0, a, 1, r1, b, 0, r2, w, lw, info);
        chkxer("Rggqrf", infot, nout, lerr, ok);
        infot = 11;
        Rggqrf(1, 1, 2, a, 1, r1, b, 1, r2, w, 1, info);
        chkxer("Rggqrf", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Rggrqf
        //
        srnamt = "Rggrqf";
        infot = 1;
        Rggrqf(-1, 0, 0, a, 1, r1, b, 1, r2, w, lw, info);
        chkxer("Rggrqf", infot, nout, lerr, ok);
        infot = 2;
        Rggrqf(0, -1, 0, a, 1, r1, b, 1, r2, w, lw, info);
        chkxer("Rggrqf", infot, nout, lerr, ok);
        infot = 3;
        Rggrqf(0, 0, -1, a, 1, r1, b, 1, r2, w, lw, info);
        chkxer("Rggrqf", infot, nout, lerr, ok);
        infot = 5;
        Rggrqf(0, 0, 0, a, 0, r1, b, 1, r2, w, lw, info);
        chkxer("Rggrqf", infot, nout, lerr, ok);
        infot = 8;
        Rggrqf(0, 0, 0, a, 1, r1, b, 0, r2, w, lw, info);
        chkxer("Rggrqf", infot, nout, lerr, ok);
        infot = 11;
        Rggrqf(1, 1, 2, a, 1, r1, b, 1, r2, w, 1, info);
        chkxer("Rggrqf", infot, nout, lerr, ok);
        nt += 6;
        //
        //     Test error exits for the DGS, DGV, DGX, and DXV paths.
        //
    } else if (Mlsamen(3, path, "DGS") || Mlsamen(3, path, "DGV") || Mlsamen(3, path, "DGX") || Mlsamen(3, path, "DXV")) {
        //
        //        Rgges
        //
        srnamt = "Rgges ";
        infot = 1;
        Rgges("/", "N", "S", Rlctes, 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        infot = 2;
        Rgges("N", "/", "S", Rlctes, 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        infot = 3;
        Rgges("N", "V", "/", Rlctes, 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        infot = 5;
        Rgges("N", "V", "S", Rlctes, -1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        infot = 7;
        Rgges("N", "V", "S", Rlctes, 1, a, 0, b, 1, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        infot = 9;
        Rgges("N", "V", "S", Rlctes, 1, a, 1, b, 0, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        infot = 15;
        Rgges("N", "V", "S", Rlctes, 1, a, 1, b, 1, sdim, r1, r2, r3, q, 0, u, 1, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        infot = 15;
        Rgges("V", "V", "S", Rlctes, 2, a, 2, b, 2, sdim, r1, r2, r3, q, 1, u, 2, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        infot = 17;
        Rgges("N", "V", "S", Rlctes, 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 0, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        infot = 17;
        Rgges("V", "V", "S", Rlctes, 2, a, 2, b, 2, sdim, r1, r2, r3, q, 2, u, 1, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        infot = 19;
        Rgges("V", "V", "S", Rlctes, 2, a, 2, b, 2, sdim, r1, r2, r3, q, 2, u, 2, w, 1, bw, info);
        chkxer("Rgges ", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Rgges3
        //
        srnamt = "Rgges3 ";
        infot = 1;
        Rgges3("/", "N", "S", Rlctes, 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        infot = 2;
        Rgges3("N", "/", "S", Rlctes, 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        infot = 3;
        Rgges3("N", "V", "/", Rlctes, 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        infot = 5;
        Rgges3("N", "V", "S", Rlctes, -1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        infot = 7;
        Rgges3("N", "V", "S", Rlctes, 1, a, 0, b, 1, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        infot = 9;
        Rgges3("N", "V", "S", Rlctes, 1, a, 1, b, 0, sdim, r1, r2, r3, q, 1, u, 1, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        infot = 15;
        Rgges3("N", "V", "S", Rlctes, 1, a, 1, b, 1, sdim, r1, r2, r3, q, 0, u, 1, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        infot = 15;
        Rgges3("V", "V", "S", Rlctes, 2, a, 2, b, 2, sdim, r1, r2, r3, q, 1, u, 2, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        infot = 17;
        Rgges3("N", "V", "S", Rlctes, 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 0, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        infot = 17;
        Rgges3("V", "V", "S", Rlctes, 2, a, 2, b, 2, sdim, r1, r2, r3, q, 2, u, 1, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        infot = 19;
        Rgges3("V", "V", "S", Rlctes, 2, a, 2, b, 2, sdim, r1, r2, r3, q, 2, u, 2, w, 1, bw, info);
        chkxer("Rgges3 ", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Rggesx
        //
        srnamt = "Rggesx";
        infot = 1;
        Rggesx("/", "N", "S", Rlctsx, "N", 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 2;
        Rggesx("N", "/", "S", Rlctsx, "N", 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 3;
        Rggesx("V", "V", "/", Rlctsx, "N", 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 5;
        Rggesx("V", "V", "S", Rlctsx, "/", 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 6;
        Rggesx("V", "V", "S", Rlctsx, "B", -1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 8;
        Rggesx("V", "V", "S", Rlctsx, "B", 1, a, 0, b, 1, sdim, r1, r2, r3, q, 1, u, 1, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 10;
        Rggesx("V", "V", "S", Rlctsx, "B", 1, a, 1, b, 0, sdim, r1, r2, r3, q, 1, u, 1, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 16;
        Rggesx("V", "V", "S", Rlctsx, "B", 1, a, 1, b, 1, sdim, r1, r2, r3, q, 0, u, 1, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 16;
        Rggesx("V", "V", "S", Rlctsx, "B", 2, a, 2, b, 2, sdim, r1, r2, r3, q, 1, u, 1, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 18;
        Rggesx("V", "V", "S", Rlctsx, "B", 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 0, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 18;
        Rggesx("V", "V", "S", Rlctsx, "B", 2, a, 2, b, 2, sdim, r1, r2, r3, q, 2, u, 1, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 22;
        Rggesx("V", "V", "S", Rlctsx, "B", 2, a, 2, b, 2, sdim, r1, r2, r3, q, 2, u, 2, rce, rcv, w, 1, iw, 1, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        infot = 24;
        Rggesx("V", "V", "S", Rlctsx, "V", 1, a, 1, b, 1, sdim, r1, r2, r3, q, 1, u, 1, rce, rcv, w, 32, iw, 0, bw, info);
        chkxer("Rggesx", infot, nout, lerr, ok);
        nt += 13;
        //
        //        Rggev
        //
        srnamt = "Rggev ";
        infot = 1;
        Rggev("/", "N", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev ", infot, nout, lerr, ok);
        infot = 2;
        Rggev("N", "/", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev ", infot, nout, lerr, ok);
        infot = 3;
        Rggev("V", "V", -1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev ", infot, nout, lerr, ok);
        infot = 5;
        Rggev("V", "V", 1, a, 0, b, 1, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev ", infot, nout, lerr, ok);
        infot = 7;
        Rggev("V", "V", 1, a, 1, b, 0, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev ", infot, nout, lerr, ok);
        infot = 12;
        Rggev("N", "V", 1, a, 1, b, 1, r1, r2, r3, q, 0, u, 1, w, 1, info);
        chkxer("Rggev ", infot, nout, lerr, ok);
        infot = 12;
        Rggev("V", "V", 2, a, 2, b, 2, r1, r2, r3, q, 1, u, 2, w, 1, info);
        chkxer("Rggev ", infot, nout, lerr, ok);
        infot = 14;
        Rggev("V", "N", 2, a, 2, b, 2, r1, r2, r3, q, 2, u, 0, w, 1, info);
        chkxer("Rggev ", infot, nout, lerr, ok);
        infot = 14;
        Rggev("V", "V", 2, a, 2, b, 2, r1, r2, r3, q, 2, u, 1, w, 1, info);
        chkxer("Rggev ", infot, nout, lerr, ok);
        infot = 16;
        Rggev("V", "V", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev ", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Rggev3
        //
        srnamt = "Rggev3 ";
        infot = 1;
        Rggev3("/", "N", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev3 ", infot, nout, lerr, ok);
        infot = 2;
        Rggev3("N", "/", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev3 ", infot, nout, lerr, ok);
        infot = 3;
        Rggev3("V", "V", -1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev3 ", infot, nout, lerr, ok);
        infot = 5;
        Rggev3("V", "V", 1, a, 0, b, 1, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev3 ", infot, nout, lerr, ok);
        infot = 7;
        Rggev3("V", "V", 1, a, 1, b, 0, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev3 ", infot, nout, lerr, ok);
        infot = 12;
        Rggev3("N", "V", 1, a, 1, b, 1, r1, r2, r3, q, 0, u, 1, w, 1, info);
        chkxer("Rggev3 ", infot, nout, lerr, ok);
        infot = 12;
        Rggev3("V", "V", 2, a, 2, b, 2, r1, r2, r3, q, 1, u, 2, w, 1, info);
        chkxer("Rggev3 ", infot, nout, lerr, ok);
        infot = 14;
        Rggev3("V", "N", 2, a, 2, b, 2, r1, r2, r3, q, 2, u, 0, w, 1, info);
        chkxer("Rggev3 ", infot, nout, lerr, ok);
        infot = 14;
        Rggev3("V", "V", 2, a, 2, b, 2, r1, r2, r3, q, 2, u, 1, w, 1, info);
        chkxer("Rggev3 ", infot, nout, lerr, ok);
        infot = 16;
        Rggev3("V", "V", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, w, 1, info);
        chkxer("Rggev3 ", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Rggevx
        //
        srnamt = "Rggevx";
        infot = 1;
        Rggevx("/", "N", "N", "N", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 2;
        Rggevx("N", "/", "N", "N", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 3;
        Rggevx("N", "N", "/", "N", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 4;
        Rggevx("N", "N", "N", "/", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 5;
        Rggevx("N", "N", "N", "N", -1, a, 1, b, 1, r1, r2, r3, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 7;
        Rggevx("N", "N", "N", "N", 1, a, 0, b, 1, r1, r2, r3, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 9;
        Rggevx("N", "N", "N", "N", 1, a, 1, b, 0, r1, r2, r3, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 14;
        Rggevx("N", "N", "N", "N", 1, a, 1, b, 1, r1, r2, r3, q, 0, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 14;
        Rggevx("N", "V", "N", "N", 2, a, 2, b, 2, r1, r2, r3, q, 1, u, 2, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 16;
        Rggevx("N", "N", "N", "N", 1, a, 1, b, 1, r1, r2, r3, q, 1, u, 0, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 16;
        Rggevx("N", "N", "V", "N", 2, a, 2, b, 2, r1, r2, r3, q, 2, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        infot = 26;
        Rggevx("N", "N", "V", "N", 2, a, 2, b, 2, r1, r2, r3, q, 2, u, 2, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, iw, bw, info);
        chkxer("Rggevx", infot, nout, lerr, ok);
        nt += 12;
        //
        //        Rtgexc
        //
        srnamt = "Rtgexc";
        infot = 3;
        Rtgexc(true, true, -1, a, 1, b, 1, q, 1, z, 1, ifst, ilst, w, 1, info);
        chkxer("Rtgexc", infot, nout, lerr, ok);
        infot = 5;
        Rtgexc(true, true, 1, a, 0, b, 1, q, 1, z, 1, ifst, ilst, w, 1, info);
        chkxer("Rtgexc", infot, nout, lerr, ok);
        infot = 7;
        Rtgexc(true, true, 1, a, 1, b, 0, q, 1, z, 1, ifst, ilst, w, 1, info);
        chkxer("Rtgexc", infot, nout, lerr, ok);
        infot = 9;
        Rtgexc(false, true, 1, a, 1, b, 1, q, 0, z, 1, ifst, ilst, w, 1, info);
        chkxer("Rtgexc", infot, nout, lerr, ok);
        infot = 9;
        Rtgexc(true, true, 1, a, 1, b, 1, q, 0, z, 1, ifst, ilst, w, 1, info);
        chkxer("Rtgexc", infot, nout, lerr, ok);
        infot = 11;
        Rtgexc(true, false, 1, a, 1, b, 1, q, 1, z, 0, ifst, ilst, w, 1, info);
        chkxer("Rtgexc", infot, nout, lerr, ok);
        infot = 11;
        Rtgexc(true, true, 1, a, 1, b, 1, q, 1, z, 0, ifst, ilst, w, 1, info);
        chkxer("Rtgexc", infot, nout, lerr, ok);
        infot = 15;
        Rtgexc(true, true, 1, a, 1, b, 1, q, 1, z, 1, ifst, ilst, w, 0, info);
        chkxer("Rtgexc", infot, nout, lerr, ok);
        nt += 8;
        //
        //        Rtgsen
        //
        srnamt = "Rtgsen";
        infot = 1;
        Rtgsen(-1, true, true, sel, 1, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 5;
        Rtgsen(1, true, true, sel, -1, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 7;
        Rtgsen(1, true, true, sel, 1, a, 0, b, 1, r1, r2, r3, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 9;
        Rtgsen(1, true, true, sel, 1, a, 1, b, 0, r1, r2, r3, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 14;
        Rtgsen(1, true, true, sel, 1, a, 1, b, 1, r1, r2, r3, q, 0, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 16;
        Rtgsen(1, true, true, sel, 1, a, 1, b, 1, r1, r2, r3, q, 1, z, 0, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 22;
        Rtgsen(0, true, true, sel, 1, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 22;
        Rtgsen(1, true, true, sel, 1, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 22;
        Rtgsen(2, true, true, sel, 1, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 24;
        Rtgsen(0, true, true, sel, 1, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, m, tola, tolb, rcv, w, 20, iw, 0, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 24;
        Rtgsen(1, true, true, sel, 1, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, m, tola, tolb, rcv, w, 20, iw, 0, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        infot = 24;
        Rtgsen(2, true, true, sel, 1, a, 1, b, 1, r1, r2, r3, q, 1, z, 1, m, tola, tolb, rcv, w, 20, iw, 1, info);
        chkxer("Rtgsen", infot, nout, lerr, ok);
        nt += 12;
        //
        //        Rtgsna
        //
        srnamt = "Rtgsna";
        infot = 1;
        Rtgsna("/", "A", sel, 1, a, 1, b, 1, q, 1, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Rtgsna", infot, nout, lerr, ok);
        infot = 2;
        Rtgsna("B", "/", sel, 1, a, 1, b, 1, q, 1, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Rtgsna", infot, nout, lerr, ok);
        infot = 4;
        Rtgsna("B", "A", sel, -1, a, 1, b, 1, q, 1, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Rtgsna", infot, nout, lerr, ok);
        infot = 6;
        Rtgsna("B", "A", sel, 1, a, 0, b, 1, q, 1, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Rtgsna", infot, nout, lerr, ok);
        infot = 8;
        Rtgsna("B", "A", sel, 1, a, 1, b, 0, q, 1, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Rtgsna", infot, nout, lerr, ok);
        infot = 10;
        Rtgsna("E", "A", sel, 1, a, 1, b, 1, q, 0, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Rtgsna", infot, nout, lerr, ok);
        infot = 12;
        Rtgsna("E", "A", sel, 1, a, 1, b, 1, q, 1, u, 0, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Rtgsna", infot, nout, lerr, ok);
        infot = 15;
        Rtgsna("E", "A", sel, 1, a, 1, b, 1, q, 1, u, 1, r1, r2, 0, m, w, 1, iw, info);
        chkxer("Rtgsna", infot, nout, lerr, ok);
        infot = 18;
        Rtgsna("E", "A", sel, 1, a, 1, b, 1, q, 1, u, 1, r1, r2, 1, m, w, 0, iw, info);
        chkxer("Rtgsna", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Rtgsyl
        //
        srnamt = "Rtgsyl";
        infot = 1;
        Rtgsyl("/", 0, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 2;
        Rtgsyl("N", -1, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 3;
        Rtgsyl("N", 0, 0, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 4;
        Rtgsyl("N", 0, 1, 0, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 6;
        Rtgsyl("N", 0, 1, 1, a, 0, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 8;
        Rtgsyl("N", 0, 1, 1, a, 1, b, 0, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 10;
        Rtgsyl("N", 0, 1, 1, a, 1, b, 1, q, 0, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 12;
        Rtgsyl("N", 0, 1, 1, a, 1, b, 1, q, 1, u, 0, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 14;
        Rtgsyl("N", 0, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 0, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 16;
        Rtgsyl("N", 0, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 0, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 20;
        Rtgsyl("N", 1, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        infot = 20;
        Rtgsyl("N", 2, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Rtgsyl", infot, nout, lerr, ok);
        nt += 12;
    }
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a3,' routines passed the tests of the error exits (',i3,"
                    "' tests done)')"),
            path, nt;
    } else {
        write(nout, "(' *** ',a3,' routines failed the tests of the error ','exits ***')"), path;
    }
    //
    //     End of Rerrgg
    //
}
