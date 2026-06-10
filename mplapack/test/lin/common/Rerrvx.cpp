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

// Derived from LAPACK routine DERRVX.
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

void Rerrvx(fem::str_cref path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    static const char *format_9999 = "(1x,a3,' drivers passed the tests of the error exits')";
    static const char *format_9998 = "(' *** ',a3,' drivers failed the tests of the error ','exits ***')";
    //
    nout = nunit;
    write(nout, star);
    fem::str<2> c2 = path(2, 3);
    //
    // Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 4;
    INTEGER i = 0;
    REAL a[nmax * nmax];
    REAL af[nmax * nmax];
    REAL b[nmax];
    REAL e[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    REAL w[2 * nmax];
    REAL x[nmax];
    REAL c[nmax];
    REAL r[nmax];
    INTEGER ip[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * nmax] = 1.0 / castREAL(i + j);
            af[(i - 1) + (j - 1) * nmax] = 1.0 / castREAL(i + j);
        }
        b[j - 1] = 0.0;
        e[j - 1] = 0.0;
        r1[j - 1] = 0.0;
        r2[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
        c[j - 1] = 0.0;
        r[j - 1] = 0.0;
        ip[j - 1] = j;
    }
    fem::str<1> eq = " ";
    ok = true;
    //
    INTEGER info = 0;
    REAL rcond = 0.0;
    INTEGER iw[nmax];
    if (Mlsamen(2, c2.elems, "GE")) {
        //
        // Rgesv
        //
        srnamt = "Rgesv";
        infot = 1;
        Rgesv(-1, 0, a, 1, ip, b, 1, info);
        Chkxer("Rgesv", infot, nout, lerr, ok);
        infot = 2;
        Rgesv(0, -1, a, 1, ip, b, 1, info);
        Chkxer("Rgesv", infot, nout, lerr, ok);
        infot = 4;
        Rgesv(2, 1, a, 1, ip, b, 2, info);
        Chkxer("Rgesv", infot, nout, lerr, ok);
        infot = 7;
        Rgesv(2, 1, a, 2, ip, b, 1, info);
        Chkxer("Rgesv", infot, nout, lerr, ok);
        //
        // Rgesvx
        //
        srnamt = "Rgesvx";
        infot = 1;
        Rgesvx("/", "N", 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 2;
        Rgesvx("N", "/", 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 3;
        Rgesvx("N", "N", -1, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 4;
        Rgesvx("N", "N", 0, -1, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 6;
        Rgesvx("N", "N", 2, 1, a, 1, af, 2, ip, eq.elems, r, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 8;
        Rgesvx("N", "N", 2, 1, a, 2, af, 1, ip, eq.elems, r, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 10;
        eq = "/";
        Rgesvx("F", "N", 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 11;
        eq = "R";
        Rgesvx("F", "N", 1, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 12;
        eq = "C";
        Rgesvx("F", "N", 1, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 14;
        Rgesvx("N", "N", 2, 1, a, 2, af, 2, ip, eq.elems, r, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 16;
        Rgesvx("N", "N", 2, 1, a, 2, af, 2, ip, eq.elems, r, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgesvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "GB")) {
        //
        // Rgbsv
        //
        srnamt = "Rgbsv";
        infot = 1;
        Rgbsv(-1, 0, 0, 0, a, 1, ip, b, 1, info);
        Chkxer("Rgbsv", infot, nout, lerr, ok);
        infot = 2;
        Rgbsv(1, -1, 0, 0, a, 1, ip, b, 1, info);
        Chkxer("Rgbsv", infot, nout, lerr, ok);
        infot = 3;
        Rgbsv(1, 0, -1, 0, a, 1, ip, b, 1, info);
        Chkxer("Rgbsv", infot, nout, lerr, ok);
        infot = 4;
        Rgbsv(0, 0, 0, -1, a, 1, ip, b, 1, info);
        Chkxer("Rgbsv", infot, nout, lerr, ok);
        infot = 6;
        Rgbsv(1, 1, 1, 0, a, 3, ip, b, 1, info);
        Chkxer("Rgbsv", infot, nout, lerr, ok);
        infot = 9;
        Rgbsv(2, 0, 0, 0, a, 1, ip, b, 1, info);
        Chkxer("Rgbsv", infot, nout, lerr, ok);
        //
        // Rgbsvx
        //
        srnamt = "Rgbsvx";
        infot = 1;
        Rgbsvx("/", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 2;
        Rgbsvx("N", "/", 0, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 3;
        Rgbsvx("N", "N", -1, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 4;
        Rgbsvx("N", "N", 1, -1, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 5;
        Rgbsvx("N", "N", 1, 0, -1, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 6;
        Rgbsvx("N", "N", 0, 0, 0, -1, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 8;
        Rgbsvx("N", "N", 1, 1, 1, 0, a, 2, af, 4, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 10;
        Rgbsvx("N", "N", 1, 1, 1, 0, a, 3, af, 3, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 12;
        eq = "/";
        Rgbsvx("F", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 13;
        eq = "R";
        Rgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 14;
        eq = "C";
        Rgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 16;
        Rgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 18;
        Rgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgbsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "GT")) {
        //
        // Rgtsv
        //
        srnamt = "Rgtsv";
        infot = 1;
        Rgtsv(-1, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], b, 1, info);
        Chkxer("Rgtsv", infot, nout, lerr, ok);
        infot = 2;
        Rgtsv(0, -1, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], b, 1, info);
        Chkxer("Rgtsv", infot, nout, lerr, ok);
        infot = 7;
        Rgtsv(2, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], b, 1, info);
        Chkxer("Rgtsv", infot, nout, lerr, ok);
        //
        // Rgtsvx
        //
        srnamt = "Rgtsvx";
        infot = 1;
        Rgtsvx("/", "N", 0, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgtsvx", infot, nout, lerr, ok);
        infot = 2;
        Rgtsvx("N", "/", 0, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgtsvx", infot, nout, lerr, ok);
        infot = 3;
        Rgtsvx("N", "N", -1, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgtsvx", infot, nout, lerr, ok);
        infot = 4;
        Rgtsvx("N", "N", 0, -1, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgtsvx", infot, nout, lerr, ok);
        infot = 14;
        Rgtsvx("N", "N", 2, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rgtsvx", infot, nout, lerr, ok);
        infot = 16;
        Rgtsvx("N", "N", 2, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rgtsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "PO")) {
        //
        // Rposv
        //
        srnamt = "Rposv";
        infot = 1;
        Rposv("/", 0, 0, a, 1, b, 1, info);
        Chkxer("Rposv", infot, nout, lerr, ok);
        infot = 2;
        Rposv("U", -1, 0, a, 1, b, 1, info);
        Chkxer("Rposv", infot, nout, lerr, ok);
        infot = 3;
        Rposv("U", 0, -1, a, 1, b, 1, info);
        Chkxer("Rposv", infot, nout, lerr, ok);
        infot = 5;
        Rposv("U", 2, 0, a, 1, b, 2, info);
        Chkxer("Rposv", infot, nout, lerr, ok);
        infot = 7;
        Rposv("U", 2, 0, a, 2, b, 1, info);
        Chkxer("Rposv", infot, nout, lerr, ok);
        //
        // Rposvx
        //
        srnamt = "Rposvx";
        infot = 1;
        Rposvx("/", "U", 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 2;
        Rposvx("N", "/", 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 3;
        Rposvx("N", "U", -1, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 4;
        Rposvx("N", "U", 0, -1, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 6;
        Rposvx("N", "U", 2, 0, a, 1, af, 2, eq.elems, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 8;
        Rposvx("N", "U", 2, 0, a, 2, af, 1, eq.elems, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 9;
        eq = "/";
        Rposvx("F", "U", 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 10;
        eq = "Y";
        Rposvx("F", "U", 1, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 12;
        Rposvx("N", "U", 2, 0, a, 2, af, 2, eq.elems, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 14;
        Rposvx("N", "U", 2, 0, a, 2, af, 2, eq.elems, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rposvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "PP")) {
        //
        // Rppsv
        //
        srnamt = "Rppsv";
        infot = 1;
        Rppsv("/", 0, 0, a, b, 1, info);
        Chkxer("Rppsv", infot, nout, lerr, ok);
        infot = 2;
        Rppsv("U", -1, 0, a, b, 1, info);
        Chkxer("Rppsv", infot, nout, lerr, ok);
        infot = 3;
        Rppsv("U", 0, -1, a, b, 1, info);
        Chkxer("Rppsv", infot, nout, lerr, ok);
        infot = 6;
        Rppsv("U", 2, 0, a, b, 1, info);
        Chkxer("Rppsv", infot, nout, lerr, ok);
        //
        // Rppsvx
        //
        srnamt = "Rppsvx";
        infot = 1;
        Rppsvx("/", "U", 0, 0, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 2;
        Rppsvx("N", "/", 0, 0, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 3;
        Rppsvx("N", "U", -1, 0, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 4;
        Rppsvx("N", "U", 0, -1, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 7;
        eq = "/";
        Rppsvx("F", "U", 0, 0, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 8;
        eq = "Y";
        Rppsvx("F", "U", 1, 0, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 10;
        Rppsvx("N", "U", 2, 0, a, af, eq.elems, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 12;
        Rppsvx("N", "U", 2, 0, a, af, eq.elems, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rppsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "PB")) {
        //
        // Rpbsv
        //
        srnamt = "Rpbsv";
        infot = 1;
        Rpbsv("/", 0, 0, 0, a, 1, b, 1, info);
        Chkxer("Rpbsv", infot, nout, lerr, ok);
        infot = 2;
        Rpbsv("U", -1, 0, 0, a, 1, b, 1, info);
        Chkxer("Rpbsv", infot, nout, lerr, ok);
        infot = 3;
        Rpbsv("U", 1, -1, 0, a, 1, b, 1, info);
        Chkxer("Rpbsv", infot, nout, lerr, ok);
        infot = 4;
        Rpbsv("U", 0, 0, -1, a, 1, b, 1, info);
        Chkxer("Rpbsv", infot, nout, lerr, ok);
        infot = 6;
        Rpbsv("U", 1, 1, 0, a, 1, b, 2, info);
        Chkxer("Rpbsv", infot, nout, lerr, ok);
        infot = 8;
        Rpbsv("U", 2, 0, 0, a, 1, b, 1, info);
        Chkxer("Rpbsv", infot, nout, lerr, ok);
        //
        // Rpbsvx
        //
        srnamt = "Rpbsvx";
        infot = 1;
        Rpbsvx("/", "U", 0, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 2;
        Rpbsvx("N", "/", 0, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 3;
        Rpbsvx("N", "U", -1, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 4;
        Rpbsvx("N", "U", 1, -1, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 5;
        Rpbsvx("N", "U", 0, 0, -1, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 7;
        Rpbsvx("N", "U", 1, 1, 0, a, 1, af, 2, eq.elems, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 9;
        Rpbsvx("N", "U", 1, 1, 0, a, 2, af, 1, eq.elems, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 10;
        eq = "/";
        Rpbsvx("F", "U", 0, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 11;
        eq = "Y";
        Rpbsvx("F", "U", 1, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 13;
        Rpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 15;
        Rpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, eq.elems, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rpbsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "PT")) {
        //
        // Rptsv
        //
        srnamt = "Rptsv";
        infot = 1;
        Rptsv(-1, 0, &a[0], &a[(2 - 1) * nmax], b, 1, info);
        Chkxer("Rptsv", infot, nout, lerr, ok);
        infot = 2;
        Rptsv(0, -1, &a[0], &a[(2 - 1) * nmax], b, 1, info);
        Chkxer("Rptsv", infot, nout, lerr, ok);
        infot = 6;
        Rptsv(2, 0, &a[0], &a[(2 - 1) * nmax], b, 1, info);
        Chkxer("Rptsv", infot, nout, lerr, ok);
        //
        // Rptsvx
        //
        srnamt = "Rptsvx";
        infot = 1;
        Rptsvx("/", 0, 0, &a[0], &a[(2 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], b, 1, x, 1, rcond, r1, r2, w, info);
        Chkxer("Rptsvx", infot, nout, lerr, ok);
        infot = 2;
        Rptsvx("N", -1, 0, &a[0], &a[(2 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], b, 1, x, 1, rcond, r1, r2, w, info);
        Chkxer("Rptsvx", infot, nout, lerr, ok);
        infot = 3;
        Rptsvx("N", 0, -1, &a[0], &a[(2 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], b, 1, x, 1, rcond, r1, r2, w, info);
        Chkxer("Rptsvx", infot, nout, lerr, ok);
        infot = 9;
        Rptsvx("N", 2, 0, &a[0], &a[(2 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], b, 1, x, 2, rcond, r1, r2, w, info);
        Chkxer("Rptsvx", infot, nout, lerr, ok);
        infot = 11;
        Rptsvx("N", 2, 0, &a[0], &a[(2 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], b, 2, x, 1, rcond, r1, r2, w, info);
        Chkxer("Rptsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SY")) {
        //
        // Rsysv
        //
        srnamt = "Rsysv";
        infot = 1;
        Rsysv("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("Rsysv", infot, nout, lerr, ok);
        infot = 2;
        Rsysv("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("Rsysv", infot, nout, lerr, ok);
        infot = 3;
        Rsysv("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        Chkxer("Rsysv", infot, nout, lerr, ok);
        infot = 5;
        Rsysv("U", 2, 0, a, 1, ip, b, 2, w, 1, info);
        Chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 8;
        Rsysv("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        Chkxer("Rsysv", infot, nout, lerr, ok);
        infot = 10;
        Rsysv("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        Chkxer("Rsysv", infot, nout, lerr, ok);
        infot = 10;
        Rsysv("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        Chkxer("Rsysv", infot, nout, lerr, ok);
        //
        // Rsysvx
        //
        srnamt = "Rsysvx";
        infot = 1;
        Rsysvx("/", "U", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        Chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 2;
        Rsysvx("N", "/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        Chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 3;
        Rsysvx("N", "U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        Chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 4;
        Rsysvx("N", "U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        Chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 6;
        Rsysvx("N", "U", 2, 0, a, 1, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 4, iw, info);
        Chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 8;
        Rsysvx("N", "U", 2, 0, a, 2, af, 1, ip, b, 2, x, 2, rcond, r1, r2, w, 4, iw, info);
        Chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 11;
        Rsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 1, x, 2, rcond, r1, r2, w, 4, iw, info);
        Chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 13;
        Rsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 1, rcond, r1, r2, w, 4, iw, info);
        Chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 18;
        Rsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 3, iw, info);
        Chkxer("Rsysvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SR")) {
        //
        // Rsysv_rook
        //
        srnamt = "Rsysv_rook";
        infot = 1;
        Rsysv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 2;
        Rsysv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 3;
        Rsysv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 5;
        Rsysv_rook("U", 2, 0, a, 1, ip, b, 2, w, 1, info);
        Chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 8;
        Rsysv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 10;
        Rsysv_rook("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        Chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 10;
        Rsysv_rook("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        Chkxer("Rsysv_rook", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SK")) {
        //
        // Rsysv_rk
        //
        // Test error exits of the driver that uses factorization
        // of a symmetric indefinite matrix with rook
        // (bounded Bunch-Kaufman) pivoting with the new storage
        // format for factors L ( or U) and D.
        //
        // L (or U) is stored in A, diagonal of D is stored on the
        // diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        srnamt = "Rsysv_rk";
        infot = 1;
        Rsysv_rk("/", 0, 0, a, 1, e, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 2;
        Rsysv_rk("U", -1, 0, a, 1, e, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 3;
        Rsysv_rk("U", 0, -1, a, 1, e, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 5;
        Rsysv_rk("U", 2, 0, a, 1, e, ip, b, 2, w, 1, info);
        Chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 9;
        Rsysv_rk("U", 2, 0, a, 2, e, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 11;
        Rsysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, 0, info);
        Chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 11;
        Rsysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, -2, info);
        Chkxer("Rsysv_rk", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SA")) {
        //
        // DSYSV_AASEN
        //
        srnamt = "Rsysv_aa";
        infot = 1;
        Rsysv_aa("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_aa", infot, nout, lerr, ok);
        infot = 2;
        Rsysv_aa("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_aa", infot, nout, lerr, ok);
        infot = 3;
        Rsysv_aa("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_aa", infot, nout, lerr, ok);
        infot = 5;
        Rsysv_aa("U", 2, 0, a, 1, ip, b, 2, w, 1, info);
        Chkxer("Rsysv_aa", infot, nout, lerr, ok);
        infot = 8;
        Rsysv_aa("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_aa", infot, nout, lerr, ok);
        infot = 10;
        Rsysv_aa("U", 3, 1, a, 3, ip, b, 3, w, 6, info);
        Chkxer("Rsysv_aa", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "S2")) {
        //
        // DSYSV_AASEN_2STAGE
        //
        srnamt = "Rsysv_aa_2stage";
        infot = 1;
        Rsysv_aa_2stage("/", 0, 0, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsysv_aa_2stage("U", -1, 0, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 3;
        Rsysv_aa_2stage("U", 0, -1, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 5;
        Rsysv_aa_2stage("U", 2, 1, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 7;
        Rsysv_aa_2stage("U", 2, 1, a, 2, a, 1, ip, ip, b, 2, w, 1, info);
        Chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 11;
        Rsysv_aa_2stage("U", 2, 1, a, 2, a, 8, ip, ip, b, 1, w, 1, info);
        Chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 13;
        Rsysv_aa_2stage("U", 2, 1, a, 2, a, 8, ip, ip, b, 2, w, 1, info);
        Chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SP")) {
        //
        // Rspsv
        //
        srnamt = "Rspsv";
        infot = 1;
        Rspsv("/", 0, 0, a, ip, b, 1, info);
        Chkxer("Rspsv", infot, nout, lerr, ok);
        infot = 2;
        Rspsv("U", -1, 0, a, ip, b, 1, info);
        Chkxer("Rspsv", infot, nout, lerr, ok);
        infot = 3;
        Rspsv("U", 0, -1, a, ip, b, 1, info);
        Chkxer("Rspsv", infot, nout, lerr, ok);
        infot = 7;
        Rspsv("U", 2, 0, a, ip, b, 1, info);
        Chkxer("Rspsv", infot, nout, lerr, ok);
        //
        // Rspsvx
        //
        srnamt = "Rspsvx";
        infot = 1;
        Rspsvx("/", "U", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rspsvx", infot, nout, lerr, ok);
        infot = 2;
        Rspsvx("N", "/", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rspsvx", infot, nout, lerr, ok);
        infot = 3;
        Rspsvx("N", "U", -1, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rspsvx", infot, nout, lerr, ok);
        infot = 4;
        Rspsvx("N", "U", 0, -1, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rspsvx", infot, nout, lerr, ok);
        infot = 9;
        Rspsvx("N", "U", 2, 0, a, af, ip, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        Chkxer("Rspsvx", infot, nout, lerr, ok);
        infot = 11;
        Rspsvx("N", "U", 2, 0, a, af, ip, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        Chkxer("Rspsvx", infot, nout, lerr, ok);
    }
    //
    // Print a summary line.
    //
    if (ok) {
        write(nout, format_9999), path;
    } else {
        write(nout, format_9998), path;
    }
    //
    // End of Rerrvx
    //
}
