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
#include <mplapack_lin.h>

#include <mplapack_debug.h>

void Rerrtsqr(const char *path, INTEGER const nunit) {
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
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    nout = nunit;
    write(nout, star);
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 2;
    INTEGER i = 0;
    REAL a[nmax * nmax];
    REAL c[nmax * nmax];
    REAL t[nmax * nmax];
    INTEGER lda = nmax;
    INTEGER ldc = nmax;
    INTEGER ldt = nmax;
    REAL w[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / castREAL(i + j);
            c[(i - 1) + (j - 1) * ldc] = 1.0 / castREAL(i + j);
            t[(i - 1) + (j - 1) * ldt] = 1.0 / castREAL(i + j);
        }
        w[j - 1] = 0.0;
    }
    ok = true;
    //
    //     Error exits for TS factorization
    //
    //     Rgeqr
    //
    infot = 1;
    REAL tau[nmax * 2];
    INTEGER info = 0;
    strncpy(srnamt, "Rgeqr", srnamt_len);
    Rgeqr(-1, 0, a, 1, tau, 1, w, 1, info);
    chkxer("Rgeqr", infot, nout, lerr, ok);
    infot = 2;
    Rgeqr(0, -1, a, 1, tau, 1, w, 1, info);
    chkxer("Rgeqr", infot, nout, lerr, ok);
    infot = 4;
    Rgeqr(1, 1, a, 0, tau, 1, w, 1, info);
    chkxer("Rgeqr", infot, nout, lerr, ok);
    infot = 6;
    Rgeqr(3, 2, a, 3, tau, 1, w, 1, info);
    chkxer("Rgeqr", infot, nout, lerr, ok);
    infot = 8;
    Rgeqr(3, 2, a, 3, tau, 7, w, 0, info);
    chkxer("Rgeqr", infot, nout, lerr, ok);
    //
    //     Rgemqr
    //
    tau[1 - 1] = 1;
    tau[2 - 1] = 1;
    tau[3 - 1] = 1;
    tau[4 - 1] = 1;
    INTEGER nb = 1;
    infot = 1;
    strncpy(srnamt, "Rgemqr", srnamt_len);
    Rgemqr("/", "N", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    infot = 2;
    Rgemqr("L", "/", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    infot = 3;
    Rgemqr("L", "N", -1, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    infot = 4;
    Rgemqr("L", "N", 0, -1, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    infot = 5;
    Rgemqr("L", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    infot = 5;
    Rgemqr("R", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    infot = 7;
    Rgemqr("L", "N", 2, 1, 0, a, 0, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    infot = 9;
    Rgemqr("R", "N", 2, 2, 1, a, 2, tau, 0, c, 1, w, 1, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    infot = 9;
    Rgemqr("L", "N", 2, 2, 1, a, 2, tau, 0, c, 1, w, 1, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    infot = 11;
    Rgemqr("L", "N", 2, 1, 1, a, 2, tau, 6, c, 0, w, 1, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    infot = 13;
    Rgemqr("L", "N", 2, 2, 1, a, 2, tau, 6, c, 2, w, 0, info);
    chkxer("Rgemqr", infot, nout, lerr, ok);
    //
    //     Rgelq
    //
    infot = 1;
    strncpy(srnamt, "Rgelq", srnamt_len);
    Rgelq(-1, 0, a, 1, tau, 1, w, 1, info);
    chkxer("Rgelq", infot, nout, lerr, ok);
    infot = 2;
    Rgelq(0, -1, a, 1, tau, 1, w, 1, info);
    chkxer("Rgelq", infot, nout, lerr, ok);
    infot = 4;
    Rgelq(1, 1, a, 0, tau, 1, w, 1, info);
    chkxer("Rgelq", infot, nout, lerr, ok);
    infot = 6;
    Rgelq(2, 3, a, 3, tau, 1, w, 1, info);
    chkxer("Rgelq", infot, nout, lerr, ok);
    infot = 8;
    Rgelq(2, 3, a, 3, tau, 7, w, 0, info);
    chkxer("Rgelq", infot, nout, lerr, ok);
    //
    //     Rgemlq
    //
    tau[1 - 1] = 1;
    tau[2 - 1] = 1;
    nb = 1;
    infot = 1;
    strncpy(srnamt, "Rgemlq", srnamt_len);
    Rgemlq("/", "N", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    infot = 2;
    Rgemlq("L", "/", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    infot = 3;
    Rgemlq("L", "N", -1, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    infot = 4;
    Rgemlq("L", "N", 0, -1, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    infot = 5;
    Rgemlq("L", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    infot = 5;
    Rgemlq("R", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    infot = 7;
    Rgemlq("L", "N", 1, 2, 0, a, 0, tau, 1, c, 1, w, 1, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    infot = 9;
    Rgemlq("R", "N", 2, 2, 1, a, 1, tau, 0, c, 1, w, 1, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    infot = 9;
    Rgemlq("L", "N", 2, 2, 1, a, 1, tau, 0, c, 1, w, 1, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    infot = 11;
    Rgemlq("L", "N", 1, 2, 1, a, 1, tau, 6, c, 0, w, 1, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    infot = 13;
    Rgemlq("L", "N", 2, 2, 1, a, 2, tau, 6, c, 2, w, 0, info);
    chkxer("Rgemlq", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrtsqr
    //
}
