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

void Cerrec(const char *path, INTEGER const nunit) {
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
    arr_2d<nmax, nmax, COMPLEX> a;
    arr_2d<nmax, nmax, COMPLEX> b;
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
    //     Test Ctrsyl
    //
    srnamt = "Ctrsyl";
    infot = 1;
    arr_2d<nmax, nmax, COMPLEX> c;
    REAL scale = 0.0;
    INTEGER info = 0;
    Ctrsyl("X", "N", 1, 0, 0, a, 1, b, 1, c, 1, scale, info);
    chkxer("Ctrsyl", infot, nout, lerr, ok);
    infot = 2;
    Ctrsyl("N", "X", 1, 0, 0, a, 1, b, 1, c, 1, scale, info);
    chkxer("Ctrsyl", infot, nout, lerr, ok);
    infot = 3;
    Ctrsyl("N", "N", 0, 0, 0, a, 1, b, 1, c, 1, scale, info);
    chkxer("Ctrsyl", infot, nout, lerr, ok);
    infot = 4;
    Ctrsyl("N", "N", 1, -1, 0, a, 1, b, 1, c, 1, scale, info);
    chkxer("Ctrsyl", infot, nout, lerr, ok);
    infot = 5;
    Ctrsyl("N", "N", 1, 0, -1, a, 1, b, 1, c, 1, scale, info);
    chkxer("Ctrsyl", infot, nout, lerr, ok);
    infot = 7;
    Ctrsyl("N", "N", 1, 2, 0, a, 1, b, 1, c, 2, scale, info);
    chkxer("Ctrsyl", infot, nout, lerr, ok);
    infot = 9;
    Ctrsyl("N", "N", 1, 0, 2, a, 1, b, 1, c, 1, scale, info);
    chkxer("Ctrsyl", infot, nout, lerr, ok);
    infot = 11;
    Ctrsyl("N", "N", 1, 2, 0, a, 2, b, 1, c, 1, scale, info);
    chkxer("Ctrsyl", infot, nout, lerr, ok);
    nt += 8;
    //
    //     Test Ctrexc
    //
    srnamt = "Ctrexc";
    INTEGER ifst = 1;
    INTEGER ilst = 1;
    infot = 1;
    Ctrexc("X", 1, a, 1, b, 1, ifst, ilst, info);
    chkxer("Ctrexc", infot, nout, lerr, ok);
    infot = 2;
    Ctrexc("N", -1, a, 1, b, 1, ifst, ilst, info);
    chkxer("Ctrexc", infot, nout, lerr, ok);
    infot = 4;
    ilst = 2;
    Ctrexc("N", 2, a, 1, b, 1, ifst, ilst, info);
    chkxer("Ctrexc", infot, nout, lerr, ok);
    infot = 6;
    Ctrexc("V", 2, a, 2, b, 1, ifst, ilst, info);
    chkxer("Ctrexc", infot, nout, lerr, ok);
    infot = 7;
    ifst = 0;
    ilst = 1;
    Ctrexc("V", 1, a, 1, b, 1, ifst, ilst, info);
    chkxer("Ctrexc", infot, nout, lerr, ok);
    infot = 7;
    ifst = 2;
    Ctrexc("V", 1, a, 1, b, 1, ifst, ilst, info);
    chkxer("Ctrexc", infot, nout, lerr, ok);
    infot = 8;
    ifst = 1;
    ilst = 0;
    Ctrexc("V", 1, a, 1, b, 1, ifst, ilst, info);
    chkxer("Ctrexc", infot, nout, lerr, ok);
    infot = 8;
    ilst = 2;
    Ctrexc("V", 1, a, 1, b, 1, ifst, ilst, info);
    chkxer("Ctrexc", infot, nout, lerr, ok);
    nt += 8;
    //
    //     Test Ctrsna
    //
    srnamt = "Ctrsna";
    infot = 1;
    arr_1d<nmax, REAL> s;
    arr_1d<nmax, REAL> sep;
    INTEGER m = 0;
    const INTEGER lw = nmax * (nmax + 2);
    arr_1d<lw, COMPLEX> work;
    arr_1d<lw, REAL> rw;
    Ctrsna("X", "A", sel, 0, a, 1, b, 1, c, 1, s, sep, 1, m, work, 1, rw, info);
    chkxer("Ctrsna", infot, nout, lerr, ok);
    infot = 2;
    Ctrsna("B", "X", sel, 0, a, 1, b, 1, c, 1, s, sep, 1, m, work, 1, rw, info);
    chkxer("Ctrsna", infot, nout, lerr, ok);
    infot = 4;
    Ctrsna("B", "A", sel, -1, a, 1, b, 1, c, 1, s, sep, 1, m, work, 1, rw, info);
    chkxer("Ctrsna", infot, nout, lerr, ok);
    infot = 6;
    Ctrsna("V", "A", sel, 2, a, 1, b, 1, c, 1, s, sep, 2, m, work, 2, rw, info);
    chkxer("Ctrsna", infot, nout, lerr, ok);
    infot = 8;
    Ctrsna("B", "A", sel, 2, a, 2, b, 1, c, 2, s, sep, 2, m, work, 2, rw, info);
    chkxer("Ctrsna", infot, nout, lerr, ok);
    infot = 10;
    Ctrsna("B", "A", sel, 2, a, 2, b, 2, c, 1, s, sep, 2, m, work, 2, rw, info);
    chkxer("Ctrsna", infot, nout, lerr, ok);
    infot = 13;
    Ctrsna("B", "A", sel, 1, a, 1, b, 1, c, 1, s, sep, 0, m, work, 1, rw, info);
    chkxer("Ctrsna", infot, nout, lerr, ok);
    infot = 13;
    Ctrsna("B", "S", sel, 2, a, 2, b, 2, c, 2, s, sep, 1, m, work, 1, rw, info);
    chkxer("Ctrsna", infot, nout, lerr, ok);
    infot = 16;
    Ctrsna("B", "A", sel, 2, a, 2, b, 2, c, 2, s, sep, 2, m, work, 1, rw, info);
    chkxer("Ctrsna", infot, nout, lerr, ok);
    nt += 9;
    //
    //     Test Ctrsen
    //
    sel[1 - 1] = false;
    srnamt = "Ctrsen";
    infot = 1;
    arr_1d<nmax, COMPLEX> x;
    Ctrsen("X", "N", sel, 0, a, 1, b, 1, x, m, s[1 - 1], sep[1 - 1], work, 1, info);
    chkxer("Ctrsen", infot, nout, lerr, ok);
    infot = 2;
    Ctrsen("N", "X", sel, 0, a, 1, b, 1, x, m, s[1 - 1], sep[1 - 1], work, 1, info);
    chkxer("Ctrsen", infot, nout, lerr, ok);
    infot = 4;
    Ctrsen("N", "N", sel, -1, a, 1, b, 1, x, m, s[1 - 1], sep[1 - 1], work, 1, info);
    chkxer("Ctrsen", infot, nout, lerr, ok);
    infot = 6;
    Ctrsen("N", "N", sel, 2, a, 1, b, 1, x, m, s[1 - 1], sep[1 - 1], work, 2, info);
    chkxer("Ctrsen", infot, nout, lerr, ok);
    infot = 8;
    Ctrsen("N", "V", sel, 2, a, 2, b, 1, x, m, s[1 - 1], sep[1 - 1], work, 1, info);
    chkxer("Ctrsen", infot, nout, lerr, ok);
    infot = 14;
    Ctrsen("N", "V", sel, 2, a, 2, b, 2, x, m, s[1 - 1], sep[1 - 1], work, 0, info);
    chkxer("Ctrsen", infot, nout, lerr, ok);
    infot = 14;
    Ctrsen("E", "V", sel, 3, a, 3, b, 3, x, m, s[1 - 1], sep[1 - 1], work, 1, info);
    chkxer("Ctrsen", infot, nout, lerr, ok);
    infot = 14;
    Ctrsen("V", "V", sel, 3, a, 3, b, 3, x, m, s[1 - 1], sep[1 - 1], work, 3, info);
    chkxer("Ctrsen", infot, nout, lerr, ok);
    nt += 8;
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
    //     End of Cerrec
    //
}
