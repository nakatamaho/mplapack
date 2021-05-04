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

void Rerrec(const char *path, INTEGER const nunit) {
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
    //     .. External Subroutines ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Executable Statements ..
    //
    nout = nunit;
    ok = true;
    INTEGER nt = 0;
    //
    //     Initialize A, B and SEL
    //
    INTEGER j = 0;
    const INTEGER nmax = 4;
    INTEGER i = 0;
    const REAL zero = 0.0;
    arr_2d<nmax, nmax, REAL> a;
    arr_2d<nmax, nmax, REAL> b;
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
            b[(i - 1) + (j - 1) * ldb] = zero;
        }
    }
    const REAL one = 1.0;
    arr_1d<nmax, bool> sel;
    for (i = 1; i <= nmax; i = i + 1) {
        a[(i - 1) + (i - 1) * lda] = one;
        sel[i - 1] = true;
    }
    //
    //     Test Rtrsyl
    //
    srnamt = "Rtrsyl";
    infot = 1;
    arr_2d<nmax, nmax, REAL> c;
    REAL scale = 0.0;
    INTEGER info = 0;
    Rtrsyl("X", "N", 1, 0, 0, a, 1, b, 1, c, 1, scale, info);
    chkxer("Rtrsyl", infot, nout, lerr, ok);
    infot = 2;
    Rtrsyl("N", "X", 1, 0, 0, a, 1, b, 1, c, 1, scale, info);
    chkxer("Rtrsyl", infot, nout, lerr, ok);
    infot = 3;
    Rtrsyl("N", "N", 0, 0, 0, a, 1, b, 1, c, 1, scale, info);
    chkxer("Rtrsyl", infot, nout, lerr, ok);
    infot = 4;
    Rtrsyl("N", "N", 1, -1, 0, a, 1, b, 1, c, 1, scale, info);
    chkxer("Rtrsyl", infot, nout, lerr, ok);
    infot = 5;
    Rtrsyl("N", "N", 1, 0, -1, a, 1, b, 1, c, 1, scale, info);
    chkxer("Rtrsyl", infot, nout, lerr, ok);
    infot = 7;
    Rtrsyl("N", "N", 1, 2, 0, a, 1, b, 1, c, 2, scale, info);
    chkxer("Rtrsyl", infot, nout, lerr, ok);
    infot = 9;
    Rtrsyl("N", "N", 1, 0, 2, a, 1, b, 1, c, 1, scale, info);
    chkxer("Rtrsyl", infot, nout, lerr, ok);
    infot = 11;
    Rtrsyl("N", "N", 1, 2, 0, a, 2, b, 1, c, 1, scale, info);
    chkxer("Rtrsyl", infot, nout, lerr, ok);
    nt += 8;
    //
    //     Test Rtrexc
    //
    srnamt = "Rtrexc";
    INTEGER ifst = 1;
    INTEGER ilst = 1;
    infot = 1;
    arr_1d<nmax, REAL> work;
    Rtrexc("X", 1, a, 1, b, 1, ifst, ilst, work, info);
    chkxer("Rtrexc", infot, nout, lerr, ok);
    infot = 2;
    Rtrexc("N", -1, a, 1, b, 1, ifst, ilst, work, info);
    chkxer("Rtrexc", infot, nout, lerr, ok);
    infot = 4;
    ilst = 2;
    Rtrexc("N", 2, a, 1, b, 1, ifst, ilst, work, info);
    chkxer("Rtrexc", infot, nout, lerr, ok);
    infot = 6;
    Rtrexc("V", 2, a, 2, b, 1, ifst, ilst, work, info);
    chkxer("Rtrexc", infot, nout, lerr, ok);
    infot = 7;
    ifst = 0;
    ilst = 1;
    Rtrexc("V", 1, a, 1, b, 1, ifst, ilst, work, info);
    chkxer("Rtrexc", infot, nout, lerr, ok);
    infot = 7;
    ifst = 2;
    Rtrexc("V", 1, a, 1, b, 1, ifst, ilst, work, info);
    chkxer("Rtrexc", infot, nout, lerr, ok);
    infot = 8;
    ifst = 1;
    ilst = 0;
    Rtrexc("V", 1, a, 1, b, 1, ifst, ilst, work, info);
    chkxer("Rtrexc", infot, nout, lerr, ok);
    infot = 8;
    ilst = 2;
    Rtrexc("V", 1, a, 1, b, 1, ifst, ilst, work, info);
    chkxer("Rtrexc", infot, nout, lerr, ok);
    nt += 8;
    //
    //     Test Rtrsna
    //
    srnamt = "Rtrsna";
    infot = 1;
    arr_1d<nmax, REAL> s;
    arr_1d<nmax, REAL> sep;
    INTEGER m = 0;
    arr_1d<nmax, int> iwork;
    Rtrsna("X", "A", sel, 0, a, 1, b, 1, c, 1, s, sep, 1, m, work, 1, iwork, info);
    chkxer("Rtrsna", infot, nout, lerr, ok);
    infot = 2;
    Rtrsna("B", "X", sel, 0, a, 1, b, 1, c, 1, s, sep, 1, m, work, 1, iwork, info);
    chkxer("Rtrsna", infot, nout, lerr, ok);
    infot = 4;
    Rtrsna("B", "A", sel, -1, a, 1, b, 1, c, 1, s, sep, 1, m, work, 1, iwork, info);
    chkxer("Rtrsna", infot, nout, lerr, ok);
    infot = 6;
    Rtrsna("V", "A", sel, 2, a, 1, b, 1, c, 1, s, sep, 2, m, work, 2, iwork, info);
    chkxer("Rtrsna", infot, nout, lerr, ok);
    infot = 8;
    Rtrsna("B", "A", sel, 2, a, 2, b, 1, c, 2, s, sep, 2, m, work, 2, iwork, info);
    chkxer("Rtrsna", infot, nout, lerr, ok);
    infot = 10;
    Rtrsna("B", "A", sel, 2, a, 2, b, 2, c, 1, s, sep, 2, m, work, 2, iwork, info);
    chkxer("Rtrsna", infot, nout, lerr, ok);
    infot = 13;
    Rtrsna("B", "A", sel, 1, a, 1, b, 1, c, 1, s, sep, 0, m, work, 1, iwork, info);
    chkxer("Rtrsna", infot, nout, lerr, ok);
    infot = 13;
    Rtrsna("B", "S", sel, 2, a, 2, b, 2, c, 2, s, sep, 1, m, work, 2, iwork, info);
    chkxer("Rtrsna", infot, nout, lerr, ok);
    infot = 16;
    Rtrsna("B", "A", sel, 2, a, 2, b, 2, c, 2, s, sep, 2, m, work, 1, iwork, info);
    chkxer("Rtrsna", infot, nout, lerr, ok);
    nt += 9;
    //
    //     Test Rtrsen
    //
    sel[1 - 1] = false;
    srnamt = "Rtrsen";
    infot = 1;
    arr_1d<nmax, REAL> wr;
    arr_1d<nmax, REAL> wi;
    Rtrsen("X", "N", sel, 0, a, 1, b, 1, wr, wi, m, s[1 - 1], sep[1 - 1], work, 1, iwork, 1, info);
    chkxer("Rtrsen", infot, nout, lerr, ok);
    infot = 2;
    Rtrsen("N", "X", sel, 0, a, 1, b, 1, wr, wi, m, s[1 - 1], sep[1 - 1], work, 1, iwork, 1, info);
    chkxer("Rtrsen", infot, nout, lerr, ok);
    infot = 4;
    Rtrsen("N", "N", sel, -1, a, 1, b, 1, wr, wi, m, s[1 - 1], sep[1 - 1], work, 1, iwork, 1, info);
    chkxer("Rtrsen", infot, nout, lerr, ok);
    infot = 6;
    Rtrsen("N", "N", sel, 2, a, 1, b, 1, wr, wi, m, s[1 - 1], sep[1 - 1], work, 2, iwork, 1, info);
    chkxer("Rtrsen", infot, nout, lerr, ok);
    infot = 8;
    Rtrsen("N", "V", sel, 2, a, 2, b, 1, wr, wi, m, s[1 - 1], sep[1 - 1], work, 1, iwork, 1, info);
    chkxer("Rtrsen", infot, nout, lerr, ok);
    infot = 15;
    Rtrsen("N", "V", sel, 2, a, 2, b, 2, wr, wi, m, s[1 - 1], sep[1 - 1], work, 0, iwork, 1, info);
    chkxer("Rtrsen", infot, nout, lerr, ok);
    infot = 15;
    Rtrsen("E", "V", sel, 3, a, 3, b, 3, wr, wi, m, s[1 - 1], sep[1 - 1], work, 1, iwork, 1, info);
    chkxer("Rtrsen", infot, nout, lerr, ok);
    infot = 15;
    Rtrsen("V", "V", sel, 3, a, 3, b, 3, wr, wi, m, s[1 - 1], sep[1 - 1], work, 3, iwork, 2, info);
    chkxer("Rtrsen", infot, nout, lerr, ok);
    infot = 17;
    Rtrsen("E", "V", sel, 2, a, 2, b, 2, wr, wi, m, s[1 - 1], sep[1 - 1], work, 1, iwork, 0, info);
    chkxer("Rtrsen", infot, nout, lerr, ok);
    infot = 17;
    Rtrsen("V", "V", sel, 3, a, 3, b, 3, wr, wi, m, s[1 - 1], sep[1 - 1], work, 4, iwork, 1, info);
    chkxer("Rtrsen", infot, nout, lerr, ok);
    nt += 10;
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a3,' routines passed the tests of the error exits (',i3,"
                    "' tests done)')"),
            path, nt;
    } else {
        write(nout, "(' *** ',a3,' routines failed the tests of the error ex','its ***')"), path;
    }
    //
    //     End of Rerrec
    //
}
