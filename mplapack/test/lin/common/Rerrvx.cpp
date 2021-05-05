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

void Rerrvx(const char *path, INTEGER const nunit) {
    common_write write(cmn);
    // COMMON infoc
    INTEGER &infot = cmn.infot;
    INTEGER &nout = cmn.nout;
    bool &ok = cmn.ok;
    bool &lerr = cmn.lerr;
    // COMMON srnamc
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
    //     .. Intrinsic Functions ..
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
            a[(i - 1) + (j - 1) * lda] = 1.0 / (i + j).real();
            af[(i - 1) + (j - 1) * ldaf] = 1.0 / (i + j).real();
        }
        b[j - 1] = 0.e+0;
        e[j - 1] = 0.e+0;
        r1[j - 1] = 0.e+0;
        r2[j - 1] = 0.e+0;
        w[j - 1] = 0.e+0;
        x[j - 1] = 0.e+0;
        c[j - 1] = 0.e+0;
        r[j - 1] = 0.e+0;
        ip[j - 1] = j;
    }
    char eq;
    ok = true;
    //
    INTEGER info = 0;
    REAL rcond = 0.0;
    INTEGER iw[nmax];
    if (Mlsamen(2, c2, "GE")) {
        //
        //        Rgesv
        //
        infot = 1;
        Rgesv(-1, 0, a, 1, ip, b, 1, info);
        chkxer("Rgesv ", infot, nout, lerr, ok);
        infot = 2;
        Rgesv(0, -1, a, 1, ip, b, 1, info);
        chkxer("Rgesv ", infot, nout, lerr, ok);
        infot = 4;
        Rgesv(2, 1, a, 1, ip, b, 2, info);
        chkxer("Rgesv ", infot, nout, lerr, ok);
        infot = 7;
        Rgesv(2, 1, a, 2, ip, b, 1, info);
        chkxer("Rgesv ", infot, nout, lerr, ok);
        //
        //        Rgesvx
        //
        infot = 1;
        Rgesvx("/", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 2;
        Rgesvx("N", "/", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 3;
        Rgesvx("N", "N", -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 4;
        Rgesvx("N", "N", 0, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 6;
        Rgesvx("N", "N", 2, 1, a, 1, af, 2, ip, eq, r, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 8;
        Rgesvx("N", "N", 2, 1, a, 2, af, 1, ip, eq, r, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 10;
        eq = "/";
        Rgesvx("F", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 11;
        eq = "R";
        Rgesvx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 12;
        eq = "C";
        Rgesvx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 14;
        Rgesvx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        infot = 16;
        Rgesvx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgesvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "GB")) {
        //
        //        Rgbsv
        //
        infot = 1;
        Rgbsv(-1, 0, 0, 0, a, 1, ip, b, 1, info);
        chkxer("Rgbsv ", infot, nout, lerr, ok);
        infot = 2;
        Rgbsv(1, -1, 0, 0, a, 1, ip, b, 1, info);
        chkxer("Rgbsv ", infot, nout, lerr, ok);
        infot = 3;
        Rgbsv(1, 0, -1, 0, a, 1, ip, b, 1, info);
        chkxer("Rgbsv ", infot, nout, lerr, ok);
        infot = 4;
        Rgbsv(0, 0, 0, -1, a, 1, ip, b, 1, info);
        chkxer("Rgbsv ", infot, nout, lerr, ok);
        infot = 6;
        Rgbsv(1, 1, 1, 0, a, 3, ip, b, 1, info);
        chkxer("Rgbsv ", infot, nout, lerr, ok);
        infot = 9;
        Rgbsv(2, 0, 0, 0, a, 1, ip, b, 1, info);
        chkxer("Rgbsv ", infot, nout, lerr, ok);
        //
        //        Rgbsvx
        //
        infot = 1;
        Rgbsvx("/", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 2;
        Rgbsvx("N", "/", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 3;
        Rgbsvx("N", "N", -1, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 4;
        Rgbsvx("N", "N", 1, -1, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 5;
        Rgbsvx("N", "N", 1, 0, -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 6;
        Rgbsvx("N", "N", 0, 0, 0, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 8;
        Rgbsvx("N", "N", 1, 1, 1, 0, a, 2, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 10;
        Rgbsvx("N", "N", 1, 1, 1, 0, a, 3, af, 3, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 12;
        eq = "/";
        Rgbsvx("F", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 13;
        eq = "R";
        Rgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 14;
        eq = "C";
        Rgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 16;
        Rgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        infot = 18;
        Rgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgbsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "GT")) {
        //
        //        Rgtsv
        //
        infot = 1;
        Rgtsv(-1, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], b, 1, info);
        chkxer("Rgtsv ", infot, nout, lerr, ok);
        infot = 2;
        Rgtsv(0, -1, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], b, 1, info);
        chkxer("Rgtsv ", infot, nout, lerr, ok);
        infot = 7;
        Rgtsv(2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], b, 1, info);
        chkxer("Rgtsv ", infot, nout, lerr, ok);
        //
        //        Rgtsvx
        //
        infot = 1;
        Rgtsvx("/", "N", 0, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgtsvx", infot, nout, lerr, ok);
        infot = 2;
        Rgtsvx("N", "/", 0, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgtsvx", infot, nout, lerr, ok);
        infot = 3;
        Rgtsvx("N", "N", -1, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgtsvx", infot, nout, lerr, ok);
        infot = 4;
        Rgtsvx("N", "N", 0, -1, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgtsvx", infot, nout, lerr, ok);
        infot = 14;
        Rgtsvx("N", "N", 2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rgtsvx", infot, nout, lerr, ok);
        infot = 16;
        Rgtsvx("N", "N", 2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rgtsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PO")) {
        //
        //        Rposv
        //
        infot = 1;
        Rposv("/", 0, 0, a, 1, b, 1, info);
        chkxer("Rposv ", infot, nout, lerr, ok);
        infot = 2;
        Rposv("U", -1, 0, a, 1, b, 1, info);
        chkxer("Rposv ", infot, nout, lerr, ok);
        infot = 3;
        Rposv("U", 0, -1, a, 1, b, 1, info);
        chkxer("Rposv ", infot, nout, lerr, ok);
        infot = 5;
        Rposv("U", 2, 0, a, 1, b, 2, info);
        chkxer("Rposv ", infot, nout, lerr, ok);
        infot = 7;
        Rposv("U", 2, 0, a, 2, b, 1, info);
        chkxer("Rposv ", infot, nout, lerr, ok);
        //
        //        Rposvx
        //
        infot = 1;
        Rposvx("/", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 2;
        Rposvx("N", "/", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 3;
        Rposvx("N", "U", -1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 4;
        Rposvx("N", "U", 0, -1, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 6;
        Rposvx("N", "U", 2, 0, a, 1, af, 2, eq, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 8;
        Rposvx("N", "U", 2, 0, a, 2, af, 1, eq, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 9;
        eq = "/";
        Rposvx("F", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 10;
        eq = "Y";
        Rposvx("F", "U", 1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 12;
        Rposvx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rposvx", infot, nout, lerr, ok);
        infot = 14;
        Rposvx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rposvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PP")) {
        //
        //        Rppsv
        //
        infot = 1;
        Rppsv("/", 0, 0, a, b, 1, info);
        chkxer("Rppsv ", infot, nout, lerr, ok);
        infot = 2;
        Rppsv("U", -1, 0, a, b, 1, info);
        chkxer("Rppsv ", infot, nout, lerr, ok);
        infot = 3;
        Rppsv("U", 0, -1, a, b, 1, info);
        chkxer("Rppsv ", infot, nout, lerr, ok);
        infot = 6;
        Rppsv("U", 2, 0, a, b, 1, info);
        chkxer("Rppsv ", infot, nout, lerr, ok);
        //
        //        Rppsvx
        //
        infot = 1;
        Rppsvx("/", "U", 0, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 2;
        Rppsvx("N", "/", 0, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 3;
        Rppsvx("N", "U", -1, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 4;
        Rppsvx("N", "U", 0, -1, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 7;
        eq = "/";
        Rppsvx("F", "U", 0, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 8;
        eq = "Y";
        Rppsvx("F", "U", 1, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 10;
        Rppsvx("N", "U", 2, 0, a, af, eq, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rppsvx", infot, nout, lerr, ok);
        infot = 12;
        Rppsvx("N", "U", 2, 0, a, af, eq, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rppsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PB")) {
        //
        //        Rpbsv
        //
        infot = 1;
        Rpbsv("/", 0, 0, 0, a, 1, b, 1, info);
        chkxer("Rpbsv ", infot, nout, lerr, ok);
        infot = 2;
        Rpbsv("U", -1, 0, 0, a, 1, b, 1, info);
        chkxer("Rpbsv ", infot, nout, lerr, ok);
        infot = 3;
        Rpbsv("U", 1, -1, 0, a, 1, b, 1, info);
        chkxer("Rpbsv ", infot, nout, lerr, ok);
        infot = 4;
        Rpbsv("U", 0, 0, -1, a, 1, b, 1, info);
        chkxer("Rpbsv ", infot, nout, lerr, ok);
        infot = 6;
        Rpbsv("U", 1, 1, 0, a, 1, b, 2, info);
        chkxer("Rpbsv ", infot, nout, lerr, ok);
        infot = 8;
        Rpbsv("U", 2, 0, 0, a, 1, b, 1, info);
        chkxer("Rpbsv ", infot, nout, lerr, ok);
        //
        //        Rpbsvx
        //
        infot = 1;
        Rpbsvx("/", "U", 0, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 2;
        Rpbsvx("N", "/", 0, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 3;
        Rpbsvx("N", "U", -1, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 4;
        Rpbsvx("N", "U", 1, -1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 5;
        Rpbsvx("N", "U", 0, 0, -1, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 7;
        Rpbsvx("N", "U", 1, 1, 0, a, 1, af, 2, eq, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 9;
        Rpbsvx("N", "U", 1, 1, 0, a, 2, af, 1, eq, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 10;
        eq = "/";
        Rpbsvx("F", "U", 0, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 11;
        eq = "Y";
        Rpbsvx("F", "U", 1, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 13;
        Rpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        infot = 15;
        Rpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, eq, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rpbsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PT")) {
        //
        //        Rptsv
        //
        infot = 1;
        Rptsv(-1, 0, &a[(1 - 1)], &a[(2 - 1) * lda], b, 1, info);
        chkxer("Rptsv ", infot, nout, lerr, ok);
        infot = 2;
        Rptsv(0, -1, &a[(1 - 1)], &a[(2 - 1) * lda], b, 1, info);
        chkxer("Rptsv ", infot, nout, lerr, ok);
        infot = 6;
        Rptsv(2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], b, 1, info);
        chkxer("Rptsv ", infot, nout, lerr, ok);
        //
        //        Rptsvx
        //
        infot = 1;
        Rptsvx("/", 0, 0, &a[(1 - 1)], &a[(2 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], b, 1, x, 1, rcond, r1, r2, w, info);
        chkxer("Rptsvx", infot, nout, lerr, ok);
        infot = 2;
        Rptsvx("N", -1, 0, &a[(1 - 1)], &a[(2 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], b, 1, x, 1, rcond, r1, r2, w, info);
        chkxer("Rptsvx", infot, nout, lerr, ok);
        infot = 3;
        Rptsvx("N", 0, -1, &a[(1 - 1)], &a[(2 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], b, 1, x, 1, rcond, r1, r2, w, info);
        chkxer("Rptsvx", infot, nout, lerr, ok);
        infot = 9;
        Rptsvx("N", 2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], b, 1, x, 2, rcond, r1, r2, w, info);
        chkxer("Rptsvx", infot, nout, lerr, ok);
        infot = 11;
        Rptsvx("N", 2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], b, 2, x, 1, rcond, r1, r2, w, info);
        chkxer("Rptsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SY")) {
        //
        //        Rsysv
        //
        infot = 1;
        Rsysv("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsysv ", infot, nout, lerr, ok);
        infot = 2;
        Rsysv("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsysv ", infot, nout, lerr, ok);
        infot = 3;
        Rsysv("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsysv ", infot, nout, lerr, ok);
        infot = 5;
        Rsysv("U", 2, 0, a, 1, ip, b, 2, w, 1, info);
        chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 8;
        Rsysv("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("Rsysv ", infot, nout, lerr, ok);
        infot = 10;
        Rsysv("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("Rsysv ", infot, nout, lerr, ok);
        infot = 10;
        Rsysv("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        chkxer("Rsysv ", infot, nout, lerr, ok);
        //
        //        Rsysvx
        //
        infot = 1;
        Rsysvx("/", "U", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 2;
        Rsysvx("N", "/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 3;
        Rsysvx("N", "U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 4;
        Rsysvx("N", "U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 6;
        Rsysvx("N", "U", 2, 0, a, 1, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 4, iw, info);
        chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 8;
        Rsysvx("N", "U", 2, 0, a, 2, af, 1, ip, b, 2, x, 2, rcond, r1, r2, w, 4, iw, info);
        chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 11;
        Rsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 1, x, 2, rcond, r1, r2, w, 4, iw, info);
        chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 13;
        Rsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 1, rcond, r1, r2, w, 4, iw, info);
        chkxer("Rsysvx", infot, nout, lerr, ok);
        infot = 18;
        Rsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 3, iw, info);
        chkxer("Rsysvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SR")) {
        //
        //        Rsysv_rook
        //
        infot = 1;
        Rsysv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 2;
        Rsysv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 3;
        Rsysv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 5;
        Rsysv_rook("U", 2, 0, a, 1, ip, b, 2, w, 1, info);
        chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 8;
        Rsysv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 10;
        Rsysv_rook("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("Rsysv_rook", infot, nout, lerr, ok);
        infot = 10;
        Rsysv_rook("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        chkxer("Rsysv_rook", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SK")) {
        //
        //        Rsysv_rk
        //
        //        Test error exits of the driver that uses factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting with the new storage
        //        format for factors L ( or U) and D.
        //
        //        L (or U) is stored in A, diagonal of D is stored on the
        //        diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        infot = 1;
        Rsysv_rk("/", 0, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 2;
        Rsysv_rk("U", -1, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 3;
        Rsysv_rk("U", 0, -1, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 5;
        Rsysv_rk("U", 2, 0, a, 1, e, ip, b, 2, w, 1, info);
        chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 9;
        Rsysv_rk("U", 2, 0, a, 2, e, ip, b, 1, w, 1, info);
        chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 11;
        Rsysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, 0, info);
        chkxer("Rsysv_rk", infot, nout, lerr, ok);
        infot = 11;
        Rsysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, -2, info);
        chkxer("Rsysv_rk", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SA")) {
        //
        //        Rsysv_aa
        //
        infot = 1;
        Rsysv_aa("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsysv_aa", infot, nout, lerr, ok);
        infot = 2;
        Rsysv_aa("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsysv_aa", infot, nout, lerr, ok);
        infot = 3;
        Rsysv_aa("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsysv_aa", infot, nout, lerr, ok);
        infot = 8;
        Rsysv_aa("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("Rsysv_aa", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "S2")) {
        //
        //        Rsysv_aaSEN_2STAGE
        //
        infot = 1;
        Rsysv_aa_2stage("/", 0, 0, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsysv_aa_2stage("U", -1, 0, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 3;
        Rsysv_aa_2stage("U", 0, -1, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 5;
        Rsysv_aa_2stage("U", 2, 1, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 11;
        Rsysv_aa_2stage("U", 2, 1, a, 2, a, 8, ip, ip, b, 1, w, 1, info);
        chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        infot = 7;
        Rsysv_aa_2stage("U", 2, 1, a, 2, a, 1, ip, ip, b, 2, w, 1, info);
        chkxer("Rsysv_aa_2stage", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SP")) {
        //
        //        Rspsv
        //
        infot = 1;
        Rspsv("/", 0, 0, a, ip, b, 1, info);
        chkxer("Rspsv ", infot, nout, lerr, ok);
        infot = 2;
        Rspsv("U", -1, 0, a, ip, b, 1, info);
        chkxer("Rspsv ", infot, nout, lerr, ok);
        infot = 3;
        Rspsv("U", 0, -1, a, ip, b, 1, info);
        chkxer("Rspsv ", infot, nout, lerr, ok);
        infot = 7;
        Rspsv("U", 2, 0, a, ip, b, 1, info);
        chkxer("Rspsv ", infot, nout, lerr, ok);
        //
        //        Rspsvx
        //
        infot = 1;
        Rspsvx("/", "U", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rspsvx", infot, nout, lerr, ok);
        infot = 2;
        Rspsvx("N", "/", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rspsvx", infot, nout, lerr, ok);
        infot = 3;
        Rspsvx("N", "U", -1, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rspsvx", infot, nout, lerr, ok);
        infot = 4;
        Rspsvx("N", "U", 0, -1, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rspsvx", infot, nout, lerr, ok);
        infot = 9;
        Rspsvx("N", "U", 2, 0, a, af, ip, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("Rspsvx", infot, nout, lerr, ok);
        infot = 11;
        Rspsvx("N", "U", 2, 0, a, af, ip, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("Rspsvx", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a3,' drivers passed the tests of the error exits')"), path;
    } else {
        write(nout, "(' *** ',a3,' drivers failed the tests of the error ','exits ***')"), path;
    }
    //
    //     End of Rerrvx
    //
}
