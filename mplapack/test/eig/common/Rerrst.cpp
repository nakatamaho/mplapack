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

void Rerrst(const char *path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    nout = nunit;
    char c2[2];
    c2[0] = path[1];
    c2[1] = path[2];
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 3;
    INTEGER i = 0;
    REAL a[nmax * nmax];
    INTEGER lda = nmax;
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / castREAL(i + j);
        }
    }
    REAL d[nmax];
    REAL e[nmax];
    INTEGER i1[nmax];
    INTEGER i2[nmax];
    REAL tau[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        d[j - 1] = castREAL(j);
        e[j - 1] = 0.0;
        i1[j - 1] = j;
        i2[j - 1] = j;
        tau[j - 1] = 1.0;
    }
    ok = true;
    INTEGER nt = 0;
    //
    //     Test error exits for the ST path.
    //
    const INTEGER lw = 20 * nmax;
    REAL w[lw];
    INTEGER info = 0;
    REAL c[nmax * nmax];
    REAL z[nmax * nmax];
    INTEGER m = 0;
    INTEGER nsplit = 0;
    REAL x[nmax];
    const INTEGER liw = 12 * nmax;
    INTEGER iw[liw];
    INTEGER i3[nmax];
    INTEGER n = 0;
    REAL r[nmax];
    REAL q[nmax * nmax];
    if (Mlsamen(2, c2, "ST")) {
        //
        //        Rsytrd
        //
        infot = 1;
        strncpy(srnamt, "Rsytrd", srnamt_len);
        Rsytrd("/", 0, a, 1, d, e, tau, w, 1, info);
        chkxer("Rsytrd", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd("U", -1, a, 1, d, e, tau, w, 1, info);
        chkxer("Rsytrd", infot, nout, lerr, ok);
        infot = 4;
        Rsytrd("U", 2, a, 1, d, e, tau, w, 1, info);
        chkxer("Rsytrd", infot, nout, lerr, ok);
        infot = 9;
        Rsytrd("U", 0, a, 1, d, e, tau, w, 0, info);
        chkxer("Rsytrd", infot, nout, lerr, ok);
        nt += 4;
        //
        //        Rsytrd_2stage
        //
        infot = 1;
        strncpy(srnamt, "Rsytrd_2stage", srnamt_len);
        Rsytrd_2stage("/", "U", 0, a, 1, d, e, tau, c, 1, w, 1, info);
        chkxer("Rsytrd_2stage", infot, nout, lerr, ok);
        infot = 1;
        Rsytrd_2stage("H", "U", 0, a, 1, d, e, tau, c, 1, w, 1, info);
        chkxer("Rsytrd_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_2stage("N", "/", 0, a, 1, d, e, tau, c, 1, w, 1, info);
        chkxer("Rsytrd_2stage", infot, nout, lerr, ok);
        infot = 3;
        Rsytrd_2stage("N", "U", -1, a, 1, d, e, tau, c, 1, w, 1, info);
        chkxer("Rsytrd_2stage", infot, nout, lerr, ok);
        infot = 5;
        Rsytrd_2stage("N", "U", 2, a, 1, d, e, tau, c, 1, w, 1, info);
        chkxer("Rsytrd_2stage", infot, nout, lerr, ok);
        infot = 10;
        Rsytrd_2stage("N", "U", 0, a, 1, d, e, tau, c, (INTEGER)0, w, 1, info);
        chkxer("Rsytrd_2stage", infot, nout, lerr, ok);
        infot = 12;
        Rsytrd_2stage("N", "U", 0, a, 1, d, e, tau, c, 1, w, 0, info);
        chkxer("Rsytrd_2stage", infot, nout, lerr, ok);
        nt += 7;
        //
        //        Rsytrd_sy2sb
        //
        infot = 1;
        strncpy(srnamt, "Rsytrd_sy2sb", srnamt_len);
        Rsytrd_sy2sb("/", 0, 0, a, 1, c, 1, tau, w, 1, info);
        chkxer("Rsytrd_sy2sb", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_sy2sb("U", -1, 0, a, 1, c, 1, tau, w, 1, info);
        chkxer("Rsytrd_sy2sb", infot, nout, lerr, ok);
        infot = 3;
        Rsytrd_sy2sb("U", 0, -1, a, 1, c, 1, tau, w, 1, info);
        chkxer("Rsytrd_sy2sb", infot, nout, lerr, ok);
        infot = 5;
        Rsytrd_sy2sb("U", 2, 0, a, 1, c, 1, tau, w, 1, info);
        chkxer("Rsytrd_sy2sb", infot, nout, lerr, ok);
        infot = 7;
        Rsytrd_sy2sb("U", 0, 2, a, 1, c, 1, tau, w, 1, info);
        chkxer("Rsytrd_sy2sb", infot, nout, lerr, ok);
        infot = 10;
        Rsytrd_sy2sb("U", 0, 0, a, 1, c, 1, tau, w, 0, info);
        chkxer("Rsytrd_sy2sb", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Rsytrd_sb2st
        //
        infot = 1;
        strncpy(srnamt, "Rsytrd_sb2st", srnamt_len);
        Rsytrd_sb2st("/", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_sb2st("Y", "/", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_sb2st("Y", "H", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 3;
        Rsytrd_sb2st("Y", "N", "/", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 4;
        Rsytrd_sb2st("Y", "N", "U", -1, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 5;
        Rsytrd_sb2st("Y", "N", "U", 0, -1, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 7;
        Rsytrd_sb2st("Y", "N", "U", 0, 1, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 11;
        Rsytrd_sb2st("Y", "N", "U", 0, 0, a, 1, d, e, c, 0, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 13;
        Rsytrd_sb2st("Y", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 0, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Rorgtr
        //
        infot = 1;
        strncpy(srnamt, "Rorgtr", srnamt_len);
        Rorgtr("/", 0, a, 1, tau, w, 1, info);
        chkxer("Rorgtr", infot, nout, lerr, ok);
        infot = 2;
        Rorgtr("U", -1, a, 1, tau, w, 1, info);
        chkxer("Rorgtr", infot, nout, lerr, ok);
        infot = 4;
        Rorgtr("U", 2, a, 1, tau, w, 1, info);
        chkxer("Rorgtr", infot, nout, lerr, ok);
        infot = 7;
        Rorgtr("U", 3, a, 3, tau, w, 1, info);
        chkxer("Rorgtr", infot, nout, lerr, ok);
        nt += 4;
        //
        //        Rormtr
        //
        infot = 1;
        strncpy(srnamt, "Rormtr", srnamt_len);
        Rormtr("/", "U", "N", 0, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Rormtr", infot, nout, lerr, ok);
        infot = 2;
        Rormtr("L", "/", "N", 0, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Rormtr", infot, nout, lerr, ok);
        infot = 3;
        Rormtr("L", "U", "/", 0, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Rormtr", infot, nout, lerr, ok);
        infot = 4;
        Rormtr("L", "U", "N", -1, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Rormtr", infot, nout, lerr, ok);
        infot = 5;
        Rormtr("L", "U", "N", 0, -1, a, 1, tau, c, 1, w, 1, info);
        chkxer("Rormtr", infot, nout, lerr, ok);
        infot = 7;
        Rormtr("L", "U", "N", 2, 0, a, 1, tau, c, 2, w, 1, info);
        chkxer("Rormtr", infot, nout, lerr, ok);
        infot = 7;
        Rormtr("R", "U", "N", 0, 2, a, 1, tau, c, 1, w, 1, info);
        chkxer("Rormtr", infot, nout, lerr, ok);
        infot = 10;
        Rormtr("L", "U", "N", 2, 0, a, 2, tau, c, 1, w, 1, info);
        chkxer("Rormtr", infot, nout, lerr, ok);
        infot = 12;
        Rormtr("L", "U", "N", 0, 2, a, 1, tau, c, 1, w, 1, info);
        chkxer("Rormtr", infot, nout, lerr, ok);
        infot = 12;
        Rormtr("R", "U", "N", 2, 0, a, 1, tau, c, 2, w, 1, info);
        chkxer("Rormtr", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Rsptrd
        //
        infot = 1;
        strncpy(srnamt, "Rsptrd", srnamt_len);
        Rsptrd("/", 0, a, d, e, tau, info);
        chkxer("Rsptrd", infot, nout, lerr, ok);
        infot = 2;
        Rsptrd("U", -1, a, d, e, tau, info);
        chkxer("Rsptrd", infot, nout, lerr, ok);
        nt += 2;
        //
        //        Ropgtr
        //
        infot = 1;
        strncpy(srnamt, "Ropgtr", srnamt_len);
        Ropgtr("/", 0, a, tau, z, 1, w, info);
        chkxer("Ropgtr", infot, nout, lerr, ok);
        infot = 2;
        Ropgtr("U", -1, a, tau, z, 1, w, info);
        chkxer("Ropgtr", infot, nout, lerr, ok);
        infot = 6;
        Ropgtr("U", 2, a, tau, z, 1, w, info);
        chkxer("Ropgtr", infot, nout, lerr, ok);
        nt += 3;
        //
        //        Ropmtr
        //
        infot = 1;
        strncpy(srnamt, "Ropmtr", srnamt_len);
        Ropmtr("/", "U", "N", 0, 0, a, tau, c, 1, w, info);
        chkxer("Ropmtr", infot, nout, lerr, ok);
        infot = 2;
        Ropmtr("L", "/", "N", 0, 0, a, tau, c, 1, w, info);
        chkxer("Ropmtr", infot, nout, lerr, ok);
        infot = 3;
        Ropmtr("L", "U", "/", 0, 0, a, tau, c, 1, w, info);
        chkxer("Ropmtr", infot, nout, lerr, ok);
        infot = 4;
        Ropmtr("L", "U", "N", -1, 0, a, tau, c, 1, w, info);
        chkxer("Ropmtr", infot, nout, lerr, ok);
        infot = 5;
        Ropmtr("L", "U", "N", 0, -1, a, tau, c, 1, w, info);
        chkxer("Ropmtr", infot, nout, lerr, ok);
        infot = 9;
        Ropmtr("L", "U", "N", 2, 0, a, tau, c, 1, w, info);
        chkxer("Ropmtr", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Rpteqr
        //
        infot = 1;
        strncpy(srnamt, "Rpteqr", srnamt_len);
        Rpteqr("/", 0, d, e, z, 1, w, info);
        chkxer("Rpteqr", infot, nout, lerr, ok);
        infot = 2;
        Rpteqr("N", -1, d, e, z, 1, w, info);
        chkxer("Rpteqr", infot, nout, lerr, ok);
        infot = 6;
        Rpteqr("V", 2, d, e, z, 1, w, info);
        chkxer("Rpteqr", infot, nout, lerr, ok);
        nt += 3;
        //
        //        Rstebz
        //
        infot = 1;
        strncpy(srnamt, "Rstebz", srnamt_len);
        Rstebz("/", "E", 0, 0.0, 1.0, 1, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        chkxer("Rstebz", infot, nout, lerr, ok);
        infot = 2;
        Rstebz("A", "/", 0, 0.0, 0.0, 0, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        chkxer("Rstebz", infot, nout, lerr, ok);
        infot = 3;
        Rstebz("A", "E", -1, 0.0, 0.0, 0, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        chkxer("Rstebz", infot, nout, lerr, ok);
        infot = 5;
        Rstebz("V", "E", 0, 0.0, 0.0, 0, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        chkxer("Rstebz", infot, nout, lerr, ok);
        infot = 6;
        Rstebz("I", "E", 0, 0.0, 0.0, 0, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        chkxer("Rstebz", infot, nout, lerr, ok);
        infot = 6;
        Rstebz("I", "E", 1, 0.0, 0.0, 2, 1, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        chkxer("Rstebz", infot, nout, lerr, ok);
        infot = 7;
        Rstebz("I", "E", 1, 0.0, 0.0, 1, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        chkxer("Rstebz", infot, nout, lerr, ok);
        infot = 7;
        Rstebz("I", "E", 1, 0.0, 0.0, 1, 2, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        chkxer("Rstebz", infot, nout, lerr, ok);
        nt += 8;
        //
        //        Rstein
        //
        infot = 1;
        strncpy(srnamt, "Rstein", srnamt_len);
        Rstein(-1, d, e, 0, x, i1, i2, z, 1, w, iw, i3, info);
        chkxer("Rstein", infot, nout, lerr, ok);
        infot = 4;
        Rstein(0, d, e, -1, x, i1, i2, z, 1, w, iw, i3, info);
        chkxer("Rstein", infot, nout, lerr, ok);
        infot = 4;
        Rstein(0, d, e, 1, x, i1, i2, z, 1, w, iw, i3, info);
        chkxer("Rstein", infot, nout, lerr, ok);
        infot = 9;
        Rstein(2, d, e, 0, x, i1, i2, z, 1, w, iw, i3, info);
        chkxer("Rstein", infot, nout, lerr, ok);
        nt += 4;
        //
        //        Rsteqr
        //
        infot = 1;
        strncpy(srnamt, "Rsteqr", srnamt_len);
        Rsteqr("/", 0, d, e, z, 1, w, info);
        chkxer("Rsteqr", infot, nout, lerr, ok);
        infot = 2;
        Rsteqr("N", -1, d, e, z, 1, w, info);
        chkxer("Rsteqr", infot, nout, lerr, ok);
        infot = 6;
        Rsteqr("V", 2, d, e, z, 1, w, info);
        chkxer("Rsteqr", infot, nout, lerr, ok);
        nt += 3;
        //
        //        Rsterf
        //
        infot = 1;
        strncpy(srnamt, "Rsterf", srnamt_len);
        Rsterf(-1, d, e, info);
        chkxer("Rsterf", infot, nout, lerr, ok);
        nt++;
        //
        //        Rstedc
        //
        infot = 1;
        strncpy(srnamt, "Rstedc", srnamt_len);
        Rstedc("/", 0, d, e, z, 1, w, 1, iw, 1, info);
        chkxer("Rstedc", infot, nout, lerr, ok);
        infot = 2;
        Rstedc("N", -1, d, e, z, 1, w, 1, iw, 1, info);
        chkxer("Rstedc", infot, nout, lerr, ok);
        infot = 6;
        Rstedc("V", 2, d, e, z, 1, w, 23, iw, 28, info);
        chkxer("Rstedc", infot, nout, lerr, ok);
        infot = 8;
        Rstedc("N", 1, d, e, z, 1, w, 0, iw, 1, info);
        chkxer("Rstedc", infot, nout, lerr, ok);
        infot = 8;
        Rstedc("I", 2, d, e, z, 2, w, 0, iw, 12, info);
        chkxer("Rstedc", infot, nout, lerr, ok);
        infot = 8;
        Rstedc("V", 2, d, e, z, 2, w, 0, iw, 28, info);
        chkxer("Rstedc", infot, nout, lerr, ok);
        infot = 10;
        Rstedc("N", 1, d, e, z, 1, w, 1, iw, 0, info);
        chkxer("Rstedc", infot, nout, lerr, ok);
        infot = 10;
        Rstedc("I", 2, d, e, z, 2, w, 19, iw, 0, info);
        chkxer("Rstedc", infot, nout, lerr, ok);
        infot = 10;
        Rstedc("V", 2, d, e, z, 2, w, 23, iw, 0, info);
        chkxer("Rstedc", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Rstevd
        //
        infot = 1;
        strncpy(srnamt, "Rstevd", srnamt_len);
        Rstevd("/", 0, d, e, z, 1, w, 1, iw, 1, info);
        chkxer("Rstevd", infot, nout, lerr, ok);
        infot = 2;
        Rstevd("N", -1, d, e, z, 1, w, 1, iw, 1, info);
        chkxer("Rstevd", infot, nout, lerr, ok);
        infot = 6;
        Rstevd("V", 2, d, e, z, 1, w, 19, iw, 12, info);
        chkxer("Rstevd", infot, nout, lerr, ok);
        infot = 8;
        Rstevd("N", 1, d, e, z, 1, w, 0, iw, 1, info);
        chkxer("Rstevd", infot, nout, lerr, ok);
        infot = 8;
        Rstevd("V", 2, d, e, z, 2, w, 12, iw, 12, info);
        chkxer("Rstevd", infot, nout, lerr, ok);
        infot = 10;
        Rstevd("N", 0, d, e, z, 1, w, 1, iw, 0, info);
        chkxer("Rstevd", infot, nout, lerr, ok);
        infot = 10;
        Rstevd("V", 2, d, e, z, 2, w, 19, iw, 11, info);
        chkxer("Rstevd", infot, nout, lerr, ok);
        nt += 7;
        //
        //        Rstev
        //
        infot = 1;
        strncpy(srnamt, "Rstev", srnamt_len);
        Rstev("/", 0, d, e, z, 1, w, info);
        chkxer("Rstev ", infot, nout, lerr, ok);
        infot = 2;
        Rstev("N", -1, d, e, z, 1, w, info);
        chkxer("Rstev ", infot, nout, lerr, ok);
        infot = 6;
        Rstev("V", 2, d, e, z, 1, w, info);
        chkxer("Rstev ", infot, nout, lerr, ok);
        nt += 3;
        //
        //        Rstevx
        //
        infot = 1;
        strncpy(srnamt, "Rstevx", srnamt_len);
        Rstevx("/", "A", 0, d, e, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rstevx", infot, nout, lerr, ok);
        infot = 2;
        Rstevx("N", "/", 0, d, e, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rstevx", infot, nout, lerr, ok);
        infot = 3;
        Rstevx("N", "A", -1, d, e, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rstevx", infot, nout, lerr, ok);
        infot = 7;
        Rstevx("N", "V", 1, d, e, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rstevx", infot, nout, lerr, ok);
        infot = 8;
        Rstevx("N", "I", 1, d, e, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rstevx", infot, nout, lerr, ok);
        infot = 8;
        Rstevx("N", "I", 1, d, e, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rstevx", infot, nout, lerr, ok);
        infot = 9;
        Rstevx("N", "I", 2, d, e, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rstevx", infot, nout, lerr, ok);
        infot = 9;
        Rstevx("N", "I", 1, d, e, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rstevx", infot, nout, lerr, ok);
        infot = 14;
        Rstevx("V", "A", 2, d, e, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rstevx", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Rstevr
        //
        n = 1;
        infot = 1;
        strncpy(srnamt, "Rstevr", srnamt_len);
        Rstevr("/", "A", 0, d, e, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rstevr", infot, nout, lerr, ok);
        infot = 2;
        Rstevr("V", "/", 0, d, e, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rstevr", infot, nout, lerr, ok);
        infot = 3;
        Rstevr("V", "A", -1, d, e, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rstevr", infot, nout, lerr, ok);
        infot = 7;
        Rstevr("V", "V", 1, d, e, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rstevr", infot, nout, lerr, ok);
        infot = 8;
        Rstevr("V", "I", 1, d, e, 0.0, 0.0, 0, 1, 0.0, m, w, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rstevr", infot, nout, lerr, ok);
        infot = 9;
        n = 2;
        Rstevr("V", "I", 2, d, e, 0.0, 0.0, 2, 1, 0.0, m, w, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rstevr", infot, nout, lerr, ok);
        infot = 14;
        n = 1;
        Rstevr("V", "I", 1, d, e, 0.0, 0.0, 1, 1, 0.0, m, w, z, 0, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rstevr", infot, nout, lerr, ok);
        infot = 17;
        Rstevr("V", "I", 1, d, e, 0.0, 0.0, 1, 1, 0.0, m, w, z, 1, iw, x, 20 * n - 1, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rstevr", infot, nout, lerr, ok);
        infot = 19;
        Rstevr("V", "I", 1, d, e, 0.0, 0.0, 1, 1, 0.0, m, w, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n - 1, info);
        chkxer("Rstevr", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Rsyevd
        //
        infot = 1;
        strncpy(srnamt, "Rsyevd", srnamt_len);
        Rsyevd("/", "U", 0, a, 1, x, w, 1, iw, 1, info);
        chkxer("Rsyevd", infot, nout, lerr, ok);
        infot = 2;
        Rsyevd("N", "/", 0, a, 1, x, w, 1, iw, 1, info);
        chkxer("Rsyevd", infot, nout, lerr, ok);
        infot = 3;
        Rsyevd("N", "U", -1, a, 1, x, w, 1, iw, 1, info);
        chkxer("Rsyevd", infot, nout, lerr, ok);
        infot = 5;
        Rsyevd("N", "U", 2, a, 1, x, w, 3, iw, 1, info);
        chkxer("Rsyevd", infot, nout, lerr, ok);
        infot = 8;
        Rsyevd("N", "U", 1, a, 1, x, w, 0, iw, 1, info);
        chkxer("Rsyevd", infot, nout, lerr, ok);
        infot = 8;
        Rsyevd("N", "U", 2, a, 2, x, w, 4, iw, 1, info);
        chkxer("Rsyevd", infot, nout, lerr, ok);
        infot = 8;
        Rsyevd("V", "U", 2, a, 2, x, w, 20, iw, 12, info);
        chkxer("Rsyevd", infot, nout, lerr, ok);
        infot = 10;
        Rsyevd("N", "U", 1, a, 1, x, w, 1, iw, 0, info);
        chkxer("Rsyevd", infot, nout, lerr, ok);
        infot = 10;
        Rsyevd("N", "U", 2, a, 2, x, w, 5, iw, 0, info);
        chkxer("Rsyevd", infot, nout, lerr, ok);
        infot = 10;
        Rsyevd("V", "U", 2, a, 2, x, w, 27, iw, 11, info);
        chkxer("Rsyevd", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Rsyevd_2stage
        //
        infot = 1;
        strncpy(srnamt, "Rsyevd_2stage", srnamt_len);
        Rsyevd_2stage("/", "U", 0, a, 1, x, w, 1, iw, 1, info);
        chkxer("Rsyevd_2stage", infot, nout, lerr, ok);
        infot = 1;
        Rsyevd_2stage("V", "U", 0, a, 1, x, w, 1, iw, 1, info);
        chkxer("Rsyevd_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsyevd_2stage("N", "/", 0, a, 1, x, w, 1, iw, 1, info);
        chkxer("Rsyevd_2stage", infot, nout, lerr, ok);
        infot = 3;
        Rsyevd_2stage("N", "U", -1, a, 1, x, w, 1, iw, 1, info);
        chkxer("Rsyevd_2stage", infot, nout, lerr, ok);
        infot = 5;
        Rsyevd_2stage("N", "U", 2, a, 1, x, w, 3, iw, 1, info);
        chkxer("Rsyevd_2stage", infot, nout, lerr, ok);
        infot = 8;
        Rsyevd_2stage("N", "U", 1, a, 1, x, w, 0, iw, 1, info);
        chkxer("Rsyevd_2stage", infot, nout, lerr, ok);
        infot = 8;
        Rsyevd_2stage("N", "U", 2, a, 2, x, w, 4, iw, 1, info);
        chkxer("Rsyevd_2stage", infot, nout, lerr, ok);
        //         INFOT = 8
        //         CALL Rsyevd_2stage( 'V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO )
        //         CALL CHKXER( 'Rsyevd_2stage', INFOT, NOUT, LERR, OK )
        infot = 10;
        strncpy(srnamt, "Rsyevd_2stage", srnamt_len);
        Rsyevd_2stage("N", "U", 1, a, 1, x, w, 1, iw, 0, info);
        chkxer("Rsyevd_2stage", infot, nout, lerr, ok);
        infot = 10;
        Rsyevd_2stage("N", "U", 2, a, 2, x, w, 25, iw, 0, info);
        chkxer("Rsyevd_2stage", infot, nout, lerr, ok);
        //         INFOT = 10
        //         CALL Rsyevd_2stage( 'V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO )
        //         CALL CHKXER( 'Rsyevd_2stage', INFOT, NOUT, LERR, OK )
        nt += 9;
        //
        //        Rsyevr
        //
        n = 1;
        infot = 1;
        strncpy(srnamt, "Rsyevr", srnamt_len);
        Rsyevr("/", "A", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        infot = 2;
        Rsyevr("V", "/", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        infot = 3;
        Rsyevr("V", "A", "/", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        infot = 4;
        Rsyevr("V", "A", "U", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        infot = 6;
        Rsyevr("V", "A", "U", 2, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        infot = 8;
        Rsyevr("V", "V", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        infot = 9;
        Rsyevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 0, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        infot = 10;
        //
        Rsyevr("V", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        infot = 15;
        Rsyevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 0, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        infot = 18;
        Rsyevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n - 1, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        infot = 20;
        Rsyevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n - 1, info);
        chkxer("Rsyevr", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Rsyevr_2stage
        //
        n = 1;
        infot = 1;
        strncpy(srnamt, "Rsyevr_2stage", srnamt_len);
        Rsyevr_2stage("/", "A", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 1;
        Rsyevr_2stage("V", "A", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsyevr_2stage("N", "/", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 3;
        Rsyevr_2stage("N", "A", "/", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 4;
        Rsyevr_2stage("N", "A", "U", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 6;
        Rsyevr_2stage("N", "A", "U", 2, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 8;
        Rsyevr_2stage("N", "V", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 9;
        Rsyevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 0, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 10;
        Rsyevr_2stage("N", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 15;
        Rsyevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 0, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 18;
        Rsyevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 0, &iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        infot = 20;
        Rsyevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 0, info);
        chkxer("Rsyevr_2stage", infot, nout, lerr, ok);
        nt += 12;
        //
        //        Rsyev
        //
        infot = 1;
        strncpy(srnamt, "Rsyev", srnamt_len);
        Rsyev("/", "U", 0, a, 1, x, w, 1, info);
        chkxer("Rsyev ", infot, nout, lerr, ok);
        infot = 2;
        Rsyev("N", "/", 0, a, 1, x, w, 1, info);
        chkxer("Rsyev ", infot, nout, lerr, ok);
        infot = 3;
        Rsyev("N", "U", -1, a, 1, x, w, 1, info);
        chkxer("Rsyev ", infot, nout, lerr, ok);
        infot = 5;
        Rsyev("N", "U", 2, a, 1, x, w, 3, info);
        chkxer("Rsyev ", infot, nout, lerr, ok);
        infot = 8;
        Rsyev("N", "U", 1, a, 1, x, w, 1, info);
        chkxer("Rsyev ", infot, nout, lerr, ok);
        nt += 5;
        //
        //        Rsyev_2stage
        //
        infot = 1;
        strncpy(srnamt, "Rsyev_2stage", srnamt_len);
        Rsyev_2stage("/", "U", 0, a, 1, x, w, 1, info);
        chkxer("Rsyev_2stage", infot, nout, lerr, ok);
        infot = 1;
        Rsyev_2stage("V", "U", 0, a, 1, x, w, 1, info);
        chkxer("Rsyev_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsyev_2stage("N", "/", 0, a, 1, x, w, 1, info);
        chkxer("Rsyev_2stage", infot, nout, lerr, ok);
        infot = 3;
        Rsyev_2stage("N", "U", -1, a, 1, x, w, 1, info);
        chkxer("Rsyev_2stage", infot, nout, lerr, ok);
        infot = 5;
        Rsyev_2stage("N", "U", 2, a, 1, x, w, 3, info);
        chkxer("Rsyev_2stage", infot, nout, lerr, ok);
        infot = 8;
        Rsyev_2stage("N", "U", 1, a, 1, x, w, 1, info);
        chkxer("Rsyev_2stage", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Rsyevx
        //
        infot = 1;
        strncpy(srnamt, "Rsyevx", srnamt_len);
        Rsyevx("/", "A", "U", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        infot = 2;
        Rsyevx("N", "/", "U", 0, a, 1, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        infot = 3;
        Rsyevx("N", "A", "/", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        infot = 4;
        Rsyevx("N", "A", "U", -1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        infot = 6;
        Rsyevx("N", "A", "U", 2, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 16, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        infot = 8;
        Rsyevx("N", "V", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        infot = 9;
        Rsyevx("N", "I", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        infot = 9;
        Rsyevx("N", "I", "U", 1, a, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        infot = 10;
        Rsyevx("N", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 16, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        infot = 10;
        Rsyevx("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        infot = 15;
        Rsyevx("V", "A", "U", 2, a, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 16, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        infot = 17;
        Rsyevx("V", "A", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsyevx", infot, nout, lerr, ok);
        nt += 12;
        //
        //        Rsyevx_2stage
        //
        infot = 1;
        strncpy(srnamt, "Rsyevx_2stage", srnamt_len);
        Rsyevx_2stage("/", "A", "U", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 1;
        Rsyevx_2stage("V", "A", "U", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsyevx_2stage("N", "/", "U", 0, a, 1, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 3;
        Rsyevx_2stage("N", "A", "/", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        infot = 4;
        Rsyevx_2stage("N", "A", "U", -1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 6;
        Rsyevx_2stage("N", "A", "U", 2, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 16, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 8;
        Rsyevx_2stage("N", "V", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 9;
        Rsyevx_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 9;
        Rsyevx_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 10;
        Rsyevx_2stage("N", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 16, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 10;
        Rsyevx_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 15;
        Rsyevx_2stage("N", "A", "U", 2, a, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 0, w, 16, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        infot = 17;
        Rsyevx_2stage("N", "A", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsyevx_2stage", infot, nout, lerr, ok);
        nt += 13;
        //
        //        Rspevd
        //
        infot = 1;
        strncpy(srnamt, "Rspevd", srnamt_len);
        Rspevd("/", "U", 0, a, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rspevd", infot, nout, lerr, ok);
        infot = 2;
        Rspevd("N", "/", 0, a, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rspevd", infot, nout, lerr, ok);
        infot = 3;
        Rspevd("N", "U", -1, a, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rspevd", infot, nout, lerr, ok);
        infot = 7;
        Rspevd("V", "U", 2, a, x, z, 1, w, 23, iw, 12, info);
        chkxer("Rspevd", infot, nout, lerr, ok);
        infot = 9;
        Rspevd("N", "U", 1, a, x, z, 1, w, 0, iw, 1, info);
        chkxer("Rspevd", infot, nout, lerr, ok);
        infot = 9;
        Rspevd("N", "U", 2, a, x, z, 1, w, 3, iw, 1, info);
        chkxer("Rspevd", infot, nout, lerr, ok);
        infot = 9;
        Rspevd("V", "U", 2, a, x, z, 2, w, 16, iw, 12, info);
        chkxer("Rspevd", infot, nout, lerr, ok);
        infot = 11;
        Rspevd("N", "U", 1, a, x, z, 1, w, 1, iw, 0, info);
        chkxer("Rspevd", infot, nout, lerr, ok);
        infot = 11;
        Rspevd("N", "U", 2, a, x, z, 1, w, 4, iw, 0, info);
        chkxer("Rspevd", infot, nout, lerr, ok);
        infot = 11;
        Rspevd("V", "U", 2, a, x, z, 2, w, 23, iw, 11, info);
        chkxer("Rspevd", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Rspev
        //
        infot = 1;
        strncpy(srnamt, "Rspev", srnamt_len);
        Rspev("/", "U", 0, a, w, z, 1, x, info);
        chkxer("Rspev ", infot, nout, lerr, ok);
        infot = 2;
        Rspev("N", "/", 0, a, w, z, 1, x, info);
        chkxer("Rspev ", infot, nout, lerr, ok);
        infot = 3;
        Rspev("N", "U", -1, a, w, z, 1, x, info);
        chkxer("Rspev ", infot, nout, lerr, ok);
        infot = 7;
        Rspev("V", "U", 2, a, w, z, 1, x, info);
        chkxer("Rspev ", infot, nout, lerr, ok);
        nt += 4;
        //
        //        Rspevx
        //
        infot = 1;
        strncpy(srnamt, "Rspevx", srnamt_len);
        Rspevx("/", "A", "U", 0, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rspevx", infot, nout, lerr, ok);
        infot = 2;
        Rspevx("N", "/", "U", 0, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rspevx", infot, nout, lerr, ok);
        infot = 3;
        Rspevx("N", "A", "/", 0, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        infot = 4;
        Rspevx("N", "A", "U", -1, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rspevx", infot, nout, lerr, ok);
        infot = 7;
        Rspevx("N", "V", "U", 1, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rspevx", infot, nout, lerr, ok);
        infot = 8;
        Rspevx("N", "I", "U", 1, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rspevx", infot, nout, lerr, ok);
        infot = 8;
        Rspevx("N", "I", "U", 1, a, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rspevx", infot, nout, lerr, ok);
        infot = 9;
        Rspevx("N", "I", "U", 2, a, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rspevx", infot, nout, lerr, ok);
        infot = 9;
        Rspevx("N", "I", "U", 1, a, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rspevx", infot, nout, lerr, ok);
        infot = 14;
        Rspevx("V", "A", "U", 2, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rspevx", infot, nout, lerr, ok);
        nt += 10;
        //
        //     Test error exits for the SB path.
        //
    } else if (Mlsamen(2, c2, "SB")) {
        //
        //        Rsbtrd
        //
        infot = 1;
        strncpy(srnamt, "Rsbtrd", srnamt_len);
        Rsbtrd("/", "U", 0, 0, a, 1, d, e, z, 1, w, info);
        chkxer("Rsbtrd", infot, nout, lerr, ok);
        infot = 2;
        Rsbtrd("N", "/", 0, 0, a, 1, d, e, z, 1, w, info);
        chkxer("Rsbtrd", infot, nout, lerr, ok);
        infot = 3;
        Rsbtrd("N", "U", -1, 0, a, 1, d, e, z, 1, w, info);
        chkxer("Rsbtrd", infot, nout, lerr, ok);
        infot = 4;
        Rsbtrd("N", "U", 0, -1, a, 1, d, e, z, 1, w, info);
        chkxer("Rsbtrd", infot, nout, lerr, ok);
        infot = 6;
        Rsbtrd("N", "U", 1, 1, a, 1, d, e, z, 1, w, info);
        chkxer("Rsbtrd", infot, nout, lerr, ok);
        infot = 10;
        Rsbtrd("V", "U", 2, 0, a, 1, d, e, z, 1, w, info);
        chkxer("Rsbtrd", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Rsytrd_sb2st
        //
        infot = 1;
        strncpy(srnamt, "Rsytrd_sb2st", srnamt_len);
        Rsytrd_sb2st("/", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_sb2st("N", "/", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_sb2st("N", "H", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 3;
        Rsytrd_sb2st("N", "N", "/", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 4;
        Rsytrd_sb2st("N", "N", "U", -1, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 5;
        Rsytrd_sb2st("N", "N", "U", 0, -1, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 7;
        Rsytrd_sb2st("N", "N", "U", 0, 1, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 11;
        Rsytrd_sb2st("N", "N", "U", 0, 0, a, 1, d, e, c, 0, w, 1, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        infot = 13;
        Rsytrd_sb2st("N", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 0, info);
        chkxer("Rsytrd_sb2st", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Rsbevd
        //
        infot = 1;
        strncpy(srnamt, "Rsbevd", srnamt_len);
        Rsbevd("/", "U", 0, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        infot = 2;
        Rsbevd("N", "/", 0, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        infot = 3;
        Rsbevd("N", "U", -1, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        infot = 4;
        Rsbevd("N", "U", 0, -1, a, 1, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        infot = 6;
        Rsbevd("N", "U", 2, 1, a, 1, x, z, 1, w, 4, iw, 1, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        infot = 9;
        Rsbevd("V", "U", 2, 1, a, 2, x, z, 1, w, 25, iw, 12, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        infot = 11;
        Rsbevd("N", "U", 1, 0, a, 1, x, z, 1, w, 0, iw, 1, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        infot = 11;
        Rsbevd("N", "U", 2, 0, a, 1, x, z, 1, w, 3, iw, 1, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        infot = 11;
        Rsbevd("V", "U", 2, 0, a, 1, x, z, 2, w, 18, iw, 12, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        infot = 13;
        Rsbevd("N", "U", 1, 0, a, 1, x, z, 1, w, 1, iw, 0, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        infot = 13;
        Rsbevd("V", "U", 2, 0, a, 1, x, z, 2, w, 25, iw, 11, info);
        chkxer("Rsbevd", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Rsbevd_2stage
        //
        infot = 1;
        strncpy(srnamt, "Rsbevd_2stage", srnamt_len);
        Rsbevd_2stage("/", "U", 0, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rsbevd_2stage", infot, nout, lerr, ok);
        infot = 1;
        Rsbevd_2stage("V", "U", 0, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rsbevd_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsbevd_2stage("N", "/", 0, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rsbevd_2stage", infot, nout, lerr, ok);
        infot = 3;
        Rsbevd_2stage("N", "U", -1, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rsbevd_2stage", infot, nout, lerr, ok);
        infot = 4;
        Rsbevd_2stage("N", "U", 0, -1, a, 1, x, z, 1, w, 1, iw, 1, info);
        chkxer("Rsbevd_2stage", infot, nout, lerr, ok);
        infot = 6;
        Rsbevd_2stage("N", "U", 2, 1, a, 1, x, z, 1, w, 4, iw, 1, info);
        chkxer("Rsbevd_2stage", infot, nout, lerr, ok);
        //         INFOT = 9
        //         CALL Rsbevd_2stage( 'V', 'U', 2, 1, A, 2, X, Z, 1, W,
        //     $                                      25, IW, 12, INFO )
        //         CALL CHKXER( 'Rsbevd_2stage', INFOT, NOUT, LERR, OK )
        infot = 11;
        Rsbevd_2stage("N", "U", 1, 0, a, 1, x, z, 1, w, 0, iw, 1, info);
        chkxer("Rsbevd_2stage", infot, nout, lerr, ok);
        infot = 11;
        Rsbevd_2stage("N", "U", 2, 0, a, 1, x, z, 1, w, 3, iw, 1, info);
        chkxer("Rsbevd_2stage", infot, nout, lerr, ok);
        //         INFOT = 11
        //         CALL Rsbevd_2stage( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
        //     $                                      18, IW, 12, INFO )
        //         CALL CHKXER( 'Rsbevd_2stage', INFOT, NOUT, LERR, OK )
        infot = 13;
        Rsbevd_2stage("N", "U", 1, 0, a, 1, x, z, 1, w, 1, iw, 0, info);
        chkxer("Rsbevd_2stage", infot, nout, lerr, ok);
        //         INFOT = 13
        //         CALL Rsbevd_2stage( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
        //     $                                      25, IW, 11, INFO )
        //         CALL CHKXER( 'Rsbevd_2stage', INFOT, NOUT, LERR, OK )
        //         NT = NT + 12
        nt += 9;
        //
        //        Rsbev
        //
        infot = 1;
        strncpy(srnamt, "Rsbev", srnamt_len);
        Rsbev("/", "U", 0, 0, a, 1, x, z, 1, w, info);
        chkxer("Rsbev ", infot, nout, lerr, ok);
        infot = 2;
        Rsbev("N", "/", 0, 0, a, 1, x, z, 1, w, info);
        chkxer("Rsbev ", infot, nout, lerr, ok);
        infot = 3;
        Rsbev("N", "U", -1, 0, a, 1, x, z, 1, w, info);
        chkxer("Rsbev ", infot, nout, lerr, ok);
        infot = 4;
        Rsbev("N", "U", 0, -1, a, 1, x, z, 1, w, info);
        chkxer("Rsbev ", infot, nout, lerr, ok);
        infot = 6;
        Rsbev("N", "U", 2, 1, a, 1, x, z, 1, w, info);
        chkxer("Rsbev ", infot, nout, lerr, ok);
        infot = 9;
        Rsbev("V", "U", 2, 0, a, 1, x, z, 1, w, info);
        chkxer("Rsbev ", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Rsbev_2stage
        //
        infot = 1;
        strncpy(srnamt, "Rsbev_2stage", srnamt_len);
        Rsbev_2stage("/", "U", 0, 0, a, 1, x, z, 1, w, 0, info);
        chkxer("Rsbev_2stage ", infot, nout, lerr, ok);
        infot = 1;
        Rsbev_2stage("V", "U", 0, 0, a, 1, x, z, 1, w, 0, info);
        chkxer("Rsbev_2stage ", infot, nout, lerr, ok);
        infot = 2;
        Rsbev_2stage("N", "/", 0, 0, a, 1, x, z, 1, w, 0, info);
        chkxer("Rsbev_2stage ", infot, nout, lerr, ok);
        infot = 3;
        Rsbev_2stage("N", "U", -1, 0, a, 1, x, z, 1, w, 0, info);
        chkxer("Rsbev_2stage ", infot, nout, lerr, ok);
        infot = 4;
        Rsbev_2stage("N", "U", 0, -1, a, 1, x, z, 1, w, 0, info);
        chkxer("Rsbev_2stage ", infot, nout, lerr, ok);
        infot = 6;
        Rsbev_2stage("N", "U", 2, 1, a, 1, x, z, 1, w, 0, info);
        chkxer("Rsbev_2stage ", infot, nout, lerr, ok);
        infot = 9;
        Rsbev_2stage("N", "U", 2, 0, a, 1, x, z, 0, w, 0, info);
        chkxer("Rsbev_2stage ", infot, nout, lerr, ok);
        infot = 11;
        Rsbev_2stage("N", "U", 0, 0, a, 1, x, z, 1, w, 0, info);
        chkxer("Rsbev_2stage ", infot, nout, lerr, ok);
        nt += 8;
        //
        //        Rsbevx
        //
        infot = 1;
        strncpy(srnamt, "Rsbevx", srnamt_len);
        Rsbevx("/", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 2;
        Rsbevx("N", "/", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 3;
        Rsbevx("N", "A", "/", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 4;
        Rsbevx("N", "A", "U", -1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 5;
        Rsbevx("N", "A", "U", 0, -1, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 7;
        Rsbevx("N", "A", "U", 2, 1, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 9;
        Rsbevx("V", "A", "U", 2, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 2, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 11;
        Rsbevx("N", "V", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 12;
        Rsbevx("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 12;
        Rsbevx("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 13;
        Rsbevx("N", "I", "U", 2, 0, a, 1, q, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 13;
        Rsbevx("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        infot = 18;
        Rsbevx("V", "A", "U", 2, 0, a, 1, q, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        chkxer("Rsbevx", infot, nout, lerr, ok);
        nt += 13;
        //
        //        Rsbevx_2stage
        //
        infot = 1;
        strncpy(srnamt, "Rsbevx_2stage", srnamt_len);
        Rsbevx_2stage("/", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        infot = 1;
        Rsbevx_2stage("V", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsbevx_2stage("N", "/", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        infot = 3;
        Rsbevx_2stage("N", "A", "/", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        infot = 4;
        Rsbevx_2stage("N", "A", "U", -1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        infot = 5;
        Rsbevx_2stage("N", "A", "U", 0, -1, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        infot = 7;
        Rsbevx_2stage("N", "A", "U", 2, 1, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        //         INFOT = 9
        //         CALL Rsbevx_2stage( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0D0,
        //     $          0.0D0, 0, 0, 0.0D0, M, X, Z, 2, W, 0, IW, I3, INFO )
        //         CALL CHKXER( 'Rsbevx_2stage', INFOT, NOUT, LERR, OK )
        infot = 11;
        Rsbevx_2stage("N", "V", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        infot = 12;
        Rsbevx_2stage("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        infot = 12;
        Rsbevx_2stage("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        infot = 13;
        Rsbevx_2stage("N", "I", "U", 2, 0, a, 1, q, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        infot = 13;
        Rsbevx_2stage("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        //         INFOT = 18
        //         CALL Rsbevx_2stage( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0D0,
        //     $          0.0D0, 0, 0, 0.0D0, M, X, Z, 1, W, 0, IW, I3, INFO )
        //         CALL CHKXER( 'Rsbevx_2stage', INFOT, NOUT, LERR, OK )
        infot = 20;
        Rsbevx_2stage("N", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        chkxer("Rsbevx_2stage", infot, nout, lerr, ok);
        //         NT = NT + 15
        nt += 13;
    }
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a3,' routines passed the tests of the error exits',' (',i3,"
                    "' tests done)')"),
            path, nt;
    } else {
        write(nout, "(' *** ',a3,' routines failed the tests of the error ','exits ***')"), path;
    }
    //
    //     End of Rerrst
    //
}
