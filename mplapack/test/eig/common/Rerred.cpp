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
#include <mplapack_eig.h>

#include <mplapack_debug.h>

bool _Rslect(REAL dummy1, REAL dummy2) { return true; }

void Rerred(const char *path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    INTEGER infot;
    INTEGER nout;
    bool ok;
    bool lerr;
    //
    static const char *format_9998 = "(' *** ',a,' failed the tests of the error exits ***')";
    static const char *format_9999 = "(1x,a,' passed the tests of the error exits (',i3,' tests done)')";
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Arrays in Common ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Executable Statements ..
    //
    nout = nunit;
    char c2[2];
    c2[0] = path[1];
    c2[1] = path[2];
    char srnamt[32];
    memset(srnamt, 0, sizeof(srnamt));
    //
    //     Initialize A
    //
    INTEGER j = 0;
    const INTEGER nmax = 4;
    INTEGER i = 0;
    const REAL zero = 0.0;
    REAL a[nmax * nmax];
    INTEGER lda = nmax;
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
        }
    }
    const REAL one = 1.0;
    for (i = 1; i <= nmax; i = i + 1) {
        a[(i - 1) + (i - 1) * lda] = one;
    }
    ok = true;
    INTEGER nt = 0;
    //
    REAL wr[nmax];
    REAL wi[nmax];
    REAL vl[nmax * nmax];
    REAL vr[nmax * nmax];
    REAL w[10 * nmax];
    INTEGER info = 0;
    INTEGER sdim = 0;
    bool b[nmax];
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL s[nmax];
    REAL abnrm = 0.0;
    REAL r1[nmax];
    REAL r2[nmax];
    INTEGER iw[2 * nmax];
    REAL u[nmax * nmax];
    REAL vt[nmax * nmax];
    INTEGER ns = 0;
    if (Mlsamen(2, c2, "EV")) {
        //
        //        Test Rgeev
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Rgeev ", 16);
        infot = 1;
        Rgeev("X", "N", 0, a, 1, wr, wi, vl, 1, vr, 1, w, 1, info);
        chkxer("Rgeev ", infot, nout, lerr, ok);
        infot = 2;
        Rgeev("N", "X", 0, a, 1, wr, wi, vl, 1, vr, 1, w, 1, info);
        chkxer("Rgeev ", infot, nout, lerr, ok);
        infot = 3;
        Rgeev("N", "N", -1, a, 1, wr, wi, vl, 1, vr, 1, w, 1, info);
        chkxer("Rgeev ", infot, nout, lerr, ok);
        infot = 5;
        Rgeev("N", "N", 2, a, 1, wr, wi, vl, 1, vr, 1, w, 6, info);
        chkxer("Rgeev ", infot, nout, lerr, ok);
        infot = 9;
        Rgeev("V", "N", 2, a, 2, wr, wi, vl, 1, vr, 1, w, 8, info);
        chkxer("Rgeev ", infot, nout, lerr, ok);
        infot = 11;
        Rgeev("N", "V", 2, a, 2, wr, wi, vl, 1, vr, 1, w, 8, info);
        chkxer("Rgeev ", infot, nout, lerr, ok);
        infot = 13;
        Rgeev("V", "V", 1, a, 1, wr, wi, vl, 1, vr, 1, w, 3, info);
        chkxer("Rgeev ", infot, nout, lerr, ok);
        nt += 7;
        //
    } else if (Mlsamen(2, c2, "ES")) {
        //
        //        Test Rgees
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Rgees ", 16);
        infot = 1;
        Rgees("X", "N", _Rslect, 0, a, 1, sdim, wr, wi, vl, 1, w, 1, b, info);
        chkxer("Rgees ", infot, nout, lerr, ok);
        infot = 2;
        Rgees("N", "X", _Rslect, 0, a, 1, sdim, wr, wi, vl, 1, w, 1, b, info);
        chkxer("Rgees ", infot, nout, lerr, ok);
        infot = 4;
        Rgees("N", "S", _Rslect, -1, a, 1, sdim, wr, wi, vl, 1, w, 1, b, info);
        chkxer("Rgees ", infot, nout, lerr, ok);
        infot = 6;
        Rgees("N", "S", _Rslect, 2, a, 1, sdim, wr, wi, vl, 1, w, 6, b, info);
        chkxer("Rgees ", infot, nout, lerr, ok);
        infot = 11;
        Rgees("V", "S", _Rslect, 2, a, 2, sdim, wr, wi, vl, 1, w, 6, b, info);
        chkxer("Rgees ", infot, nout, lerr, ok);
        infot = 13;
        Rgees("N", "S", _Rslect, 1, a, 1, sdim, wr, wi, vl, 1, w, 2, b, info);
        chkxer("Rgees ", infot, nout, lerr, ok);
        nt += 6;
        //
    } else if (Mlsamen(2, c2, "VX")) {
        //
        //        Test Rgeevx
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Rgeevx", 16);
        infot = 1;
        Rgeevx("X", "N", "N", "N", 0, a, 1, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        infot = 2;
        Rgeevx("N", "X", "N", "N", 0, a, 1, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        infot = 3;
        Rgeevx("N", "N", "X", "N", 0, a, 1, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        infot = 4;
        Rgeevx("N", "N", "N", "X", 0, a, 1, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        infot = 5;
        Rgeevx("N", "N", "N", "N", -1, a, 1, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        infot = 7;
        Rgeevx("N", "N", "N", "N", 2, a, 1, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        infot = 11;
        Rgeevx("N", "V", "N", "N", 2, a, 2, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 6, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        infot = 13;
        Rgeevx("N", "N", "V", "N", 2, a, 2, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 6, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        infot = 21;
        Rgeevx("N", "N", "N", "N", 1, a, 1, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        infot = 21;
        Rgeevx("N", "V", "N", "N", 1, a, 1, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 2, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        infot = 21;
        Rgeevx("N", "N", "V", "V", 1, a, 1, wr, wi, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 3, iw, info);
        chkxer("Rgeevx", infot, nout, lerr, ok);
        nt += 11;
        //
    } else if (Mlsamen(2, c2, "SX")) {
        //
        //        Test Rgeesx
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Rgeesx", 16);
        infot = 1;
        Rgeesx("X", "N", _Rslect, "N", 0, a, 1, sdim, wr, wi, vl, 1, r1[1 - 1], r2[1 - 1], w, 1, iw, 1, b, info);
        chkxer("Rgeesx", infot, nout, lerr, ok);
        infot = 2;
        Rgeesx("N", "X", _Rslect, "N", 0, a, 1, sdim, wr, wi, vl, 1, r1[1 - 1], r2[1 - 1], w, 1, iw, 1, b, info);
        chkxer("Rgeesx", infot, nout, lerr, ok);
        infot = 4;
        Rgeesx("N", "N", _Rslect, "X", 0, a, 1, sdim, wr, wi, vl, 1, r1[1 - 1], r2[1 - 1], w, 1, iw, 1, b, info);
        chkxer("Rgeesx", infot, nout, lerr, ok);
        infot = 5;
        Rgeesx("N", "N", _Rslect, "N", -1, a, 1, sdim, wr, wi, vl, 1, r1[1 - 1], r2[1 - 1], w, 1, iw, 1, b, info);
        chkxer("Rgeesx", infot, nout, lerr, ok);
        infot = 7;
        Rgeesx("N", "N", _Rslect, "N", 2, a, 1, sdim, wr, wi, vl, 1, r1[1 - 1], r2[1 - 1], w, 6, iw, 1, b, info);
        chkxer("Rgeesx", infot, nout, lerr, ok);
        infot = 12;
        Rgeesx("V", "N", _Rslect, "N", 2, a, 2, sdim, wr, wi, vl, 1, r1[1 - 1], r2[1 - 1], w, 6, iw, 1, b, info);
        chkxer("Rgeesx", infot, nout, lerr, ok);
        infot = 16;
        Rgeesx("N", "N", _Rslect, "N", 1, a, 1, sdim, wr, wi, vl, 1, r1[1 - 1], r2[1 - 1], w, 2, iw, 1, b, info);
        chkxer("Rgeesx", infot, nout, lerr, ok);
        nt += 7;
        //
    } else if (Mlsamen(2, c2, "BD")) {
        //
        //        Test Rgesvd
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Rgesvd", 16);
        infot = 1;
        Rgesvd("X", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, info);
        chkxer("Rgesvd", infot, nout, lerr, ok);
        infot = 2;
        Rgesvd("N", "X", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, info);
        chkxer("Rgesvd", infot, nout, lerr, ok);
        infot = 2;
        Rgesvd("O", "O", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, info);
        chkxer("Rgesvd", infot, nout, lerr, ok);
        infot = 3;
        Rgesvd("N", "N", -1, 0, a, 1, s, u, 1, vt, 1, w, 1, info);
        chkxer("Rgesvd", infot, nout, lerr, ok);
        infot = 4;
        Rgesvd("N", "N", 0, -1, a, 1, s, u, 1, vt, 1, w, 1, info);
        chkxer("Rgesvd", infot, nout, lerr, ok);
        infot = 6;
        Rgesvd("N", "N", 2, 1, a, 1, s, u, 1, vt, 1, w, 5, info);
        chkxer("Rgesvd", infot, nout, lerr, ok);
        infot = 9;
        Rgesvd("A", "N", 2, 1, a, 2, s, u, 1, vt, 1, w, 5, info);
        chkxer("Rgesvd", infot, nout, lerr, ok);
        infot = 11;
        Rgesvd("N", "A", 1, 2, a, 1, s, u, 1, vt, 1, w, 5, info);
        chkxer("Rgesvd", infot, nout, lerr, ok);
        nt = 8;
        if (ok) {
            write(nout, format_9999), srnamt, nt;
        } else {
            write(nout, format_9998);
        }
        //
        //        Test Rgesdd
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Rgesdd", 16);
        infot = 1;
        Rgesdd("X", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesdd", infot, nout, lerr, ok);
        infot = 2;
        Rgesdd("N", -1, 0, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesdd", infot, nout, lerr, ok);
        infot = 3;
        Rgesdd("N", 0, -1, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesdd", infot, nout, lerr, ok);
        infot = 5;
        Rgesdd("N", 2, 1, a, 1, s, u, 1, vt, 1, w, 5, iw, info);
        chkxer("Rgesdd", infot, nout, lerr, ok);
        infot = 8;
        Rgesdd("A", 2, 1, a, 2, s, u, 1, vt, 1, w, 5, iw, info);
        chkxer("Rgesdd", infot, nout, lerr, ok);
        infot = 10;
        Rgesdd("A", 1, 2, a, 1, s, u, 1, vt, 1, w, 5, iw, info);
        chkxer("Rgesdd", infot, nout, lerr, ok);
        nt = 6;
        if (ok) {
            write(nout, format_9999), srnamt, nt;
        } else {
            write(nout, format_9998);
        }
        //
        //        Test Rgejsv
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Rgejsv", 16);
        infot = 1;
        Rgejsv("X", "U", "V", "R", "N", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        infot = 2;
        Rgejsv("G", "X", "V", "R", "N", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        infot = 3;
        Rgejsv("G", "U", "X", "R", "N", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        infot = 4;
        Rgejsv("G", "U", "V", "X", "N", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        infot = 5;
        Rgejsv("G", "U", "V", "R", "X", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        infot = 6;
        Rgejsv("G", "U", "V", "R", "N", "X", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        infot = 7;
        Rgejsv("G", "U", "V", "R", "N", "N", -1, 0, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        infot = 8;
        Rgejsv("G", "U", "V", "R", "N", "N", 0, -1, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        infot = 10;
        Rgejsv("G", "U", "V", "R", "N", "N", 2, 1, a, 1, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        infot = 13;
        Rgejsv("G", "U", "V", "R", "N", "N", 2, 2, a, 2, s, u, 1, vt, 2, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        infot = 15;
        Rgejsv("G", "U", "V", "R", "N", "N", 2, 2, a, 2, s, u, 2, vt, 1, w, 1, iw, info);
        chkxer("Rgejsv", infot, nout, lerr, ok);
        nt = 11;
        if (ok) {
            write(nout, format_9999), srnamt, nt;
        } else {
            write(nout, format_9998);
        }
        //
        //        Test Rgesvdx
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Rgesvdx", 16);
        infot = 1;
        Rgesvdx("X", "N", "A", 0, 0, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 2;
        Rgesvdx("N", "X", "A", 0, 0, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 3;
        Rgesvdx("N", "N", "X", 0, 0, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 4;
        Rgesvdx("N", "N", "A", -1, 0, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 5;
        Rgesvdx("N", "N", "A", 0, -1, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 7;
        Rgesvdx("N", "N", "A", 2, 1, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 8;
        Rgesvdx("N", "N", "V", 2, 1, a, 2, -one, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 9;
        Rgesvdx("N", "N", "V", 2, 1, a, 2, one, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 10;
        Rgesvdx("N", "N", "I", 2, 2, a, 2, zero, zero, 0, 1, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 11;
        Rgesvdx("V", "N", "I", 2, 2, a, 2, zero, zero, 1, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 15;
        Rgesvdx("V", "N", "A", 2, 2, a, 2, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        infot = 17;
        Rgesvdx("N", "V", "A", 2, 2, a, 2, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, iw, info);
        chkxer("Rgesvdx", infot, nout, lerr, ok);
        nt = 12;
        if (ok) {
            write(nout, format_9999), srnamt, nt;
        } else {
            write(nout, format_9998);
        }
        //
        //        Test Rgesvdq
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Rgesvdq", 16);
        infot = 1;
        Rgesvdq("X", "P", "T", "A", "A", 0, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        infot = 2;
        Rgesvdq("A", "X", "T", "A", "A", 0, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        infot = 3;
        Rgesvdq("A", "P", "X", "A", "A", 0, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        infot = 4;
        Rgesvdq("A", "P", "T", "X", "A", 0, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        infot = 5;
        Rgesvdq("A", "P", "T", "A", "X", 0, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        infot = 6;
        Rgesvdq("A", "P", "T", "A", "A", -1, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        infot = 7;
        Rgesvdq("A", "P", "T", "A", "A", 0, 1, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        infot = 9;
        Rgesvdq("A", "P", "T", "A", "A", 1, 1, a, 0, s, u, 0, vt, 0, ns, iw, 1, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        infot = 12;
        Rgesvdq("A", "P", "T", "A", "A", 1, 1, a, 1, s, u, -1, vt, 0, ns, iw, 1, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        infot = 14;
        Rgesvdq("A", "P", "T", "A", "A", 1, 1, a, 1, s, u, 1, vt, -1, ns, iw, 1, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        infot = 17;
        Rgesvdq("A", "P", "T", "A", "A", 1, 1, a, 1, s, u, 1, vt, 1, ns, iw, -5, w, 1, w, 1, info);
        chkxer("Rgesvdq", infot, nout, lerr, ok);
        nt = 11;
        if (ok) {
            write(nout, format_9999), srnamt, nt;
        } else {
            write(nout, format_9998);
        }
    }
    //
    //     Print a summary line.
    //
    if (!Mlsamen(2, c2, "BD")) {
        if (ok) {
            write(nout, format_9999), srnamt, nt;
        } else {
            write(nout, format_9998);
        }
    }
    //
    //     End of Rerred
}
