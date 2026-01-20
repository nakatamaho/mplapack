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

// Derived from LAPACK routine DERRST.
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
#include <mplapack_eig.h>

void Rerrst(fem::str_cref path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    nout = nunit;
    write(nout, star);
    fem::str<2> c2 = path(2, 3);
    //
    // Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 3;
    INTEGER i = 0;
    REAL a[nmax * nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * nmax] = 1.0 / castREAL(i + j);
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
    // Test error exits for the ST path.
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
    if (Mlsamen(2, c2.elems, "ST")) {
        //
        // Rsytrd
        //
        srnamt = "DSYTRD";
        infot = 1;
        Rsytrd("/", 0, a, 1, d, e, tau, w, 1, info);
        Chkxer("DSYTRD", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd("U", -1, a, 1, d, e, tau, w, 1, info);
        Chkxer("DSYTRD", infot, nout, lerr, ok);
        infot = 4;
        Rsytrd("U", 2, a, 1, d, e, tau, w, 1, info);
        Chkxer("DSYTRD", infot, nout, lerr, ok);
        infot = 9;
        Rsytrd("U", 0, a, 1, d, e, tau, w, 0, info);
        Chkxer("DSYTRD", infot, nout, lerr, ok);
        nt += 4;
        //
        // Rsytrd_2stage
        //
        srnamt = "DSYTRD_2STAGE";
        infot = 1;
        Rsytrd_2stage("/", "U", 0, a, 1, d, e, tau, c, 1, w, 1, info);
        Chkxer("DSYTRD_2STAGE", infot, nout, lerr, ok);
        infot = 1;
        Rsytrd_2stage("H", "U", 0, a, 1, d, e, tau, c, 1, w, 1, info);
        Chkxer("DSYTRD_2STAGE", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_2stage("N", "/", 0, a, 1, d, e, tau, c, 1, w, 1, info);
        Chkxer("DSYTRD_2STAGE", infot, nout, lerr, ok);
        infot = 3;
        Rsytrd_2stage("N", "U", -1, a, 1, d, e, tau, c, 1, w, 1, info);
        Chkxer("DSYTRD_2STAGE", infot, nout, lerr, ok);
        infot = 5;
        Rsytrd_2stage("N", "U", 2, a, 1, d, e, tau, c, 1, w, 1, info);
        Chkxer("DSYTRD_2STAGE", infot, nout, lerr, ok);
        infot = 10;
        Rsytrd_2stage("N", "U", 0, a, 1, d, e, tau, c, 0, w, 1, info);
        Chkxer("DSYTRD_2STAGE", infot, nout, lerr, ok);
        infot = 12;
        Rsytrd_2stage("N", "U", 0, a, 1, d, e, tau, c, 1, w, 0, info);
        Chkxer("DSYTRD_2STAGE", infot, nout, lerr, ok);
        nt += 7;
        //
        // Rsytrd_sy2sb
        //
        srnamt = "DSYTRD_SY2SB";
        infot = 1;
        Rsytrd_sy2sb("/", 0, 0, a, 1, c, 1, tau, w, 1, info);
        Chkxer("DSYTRD_SY2SB", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_sy2sb("U", -1, 0, a, 1, c, 1, tau, w, 1, info);
        Chkxer("DSYTRD_SY2SB", infot, nout, lerr, ok);
        infot = 3;
        Rsytrd_sy2sb("U", 0, -1, a, 1, c, 1, tau, w, 1, info);
        Chkxer("DSYTRD_SY2SB", infot, nout, lerr, ok);
        infot = 5;
        Rsytrd_sy2sb("U", 2, 0, a, 1, c, 1, tau, w, 1, info);
        Chkxer("DSYTRD_SY2SB", infot, nout, lerr, ok);
        infot = 7;
        Rsytrd_sy2sb("U", 0, 2, a, 1, c, 1, tau, w, 1, info);
        Chkxer("DSYTRD_SY2SB", infot, nout, lerr, ok);
        infot = 10;
        Rsytrd_sy2sb("U", 0, 0, a, 1, c, 1, tau, w, 0, info);
        Chkxer("DSYTRD_SY2SB", infot, nout, lerr, ok);
        nt += 6;
        //
        // Rsytrd_sb2st
        //
        srnamt = "DSYTRD_SB2ST";
        infot = 1;
        Rsytrd_sb2st("/", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_sb2st("Y", "/", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_sb2st("Y", "H", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 3;
        Rsytrd_sb2st("Y", "N", "/", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 4;
        Rsytrd_sb2st("Y", "N", "U", -1, 0, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 5;
        Rsytrd_sb2st("Y", "N", "U", 0, -1, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 7;
        Rsytrd_sb2st("Y", "N", "U", 0, 1, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 11;
        Rsytrd_sb2st("Y", "N", "U", 0, 0, a, 1, d, e, c, 0, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 13;
        Rsytrd_sb2st("Y", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 0, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        nt += 9;
        //
        // Rorgtr
        //
        srnamt = "DORGTR";
        infot = 1;
        Rorgtr("/", 0, a, 1, tau, w, 1, info);
        Chkxer("DORGTR", infot, nout, lerr, ok);
        infot = 2;
        Rorgtr("U", -1, a, 1, tau, w, 1, info);
        Chkxer("DORGTR", infot, nout, lerr, ok);
        infot = 4;
        Rorgtr("U", 2, a, 1, tau, w, 1, info);
        Chkxer("DORGTR", infot, nout, lerr, ok);
        infot = 7;
        Rorgtr("U", 3, a, 3, tau, w, 1, info);
        Chkxer("DORGTR", infot, nout, lerr, ok);
        nt += 4;
        //
        // Rormtr
        //
        srnamt = "DORMTR";
        infot = 1;
        Rormtr("/", "U", "N", 0, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("DORMTR", infot, nout, lerr, ok);
        infot = 2;
        Rormtr("L", "/", "N", 0, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("DORMTR", infot, nout, lerr, ok);
        infot = 3;
        Rormtr("L", "U", "/", 0, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("DORMTR", infot, nout, lerr, ok);
        infot = 4;
        Rormtr("L", "U", "N", -1, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("DORMTR", infot, nout, lerr, ok);
        infot = 5;
        Rormtr("L", "U", "N", 0, -1, a, 1, tau, c, 1, w, 1, info);
        Chkxer("DORMTR", infot, nout, lerr, ok);
        infot = 7;
        Rormtr("L", "U", "N", 2, 0, a, 1, tau, c, 2, w, 1, info);
        Chkxer("DORMTR", infot, nout, lerr, ok);
        infot = 7;
        Rormtr("R", "U", "N", 0, 2, a, 1, tau, c, 1, w, 1, info);
        Chkxer("DORMTR", infot, nout, lerr, ok);
        infot = 10;
        Rormtr("L", "U", "N", 2, 0, a, 2, tau, c, 1, w, 1, info);
        Chkxer("DORMTR", infot, nout, lerr, ok);
        infot = 12;
        Rormtr("L", "U", "N", 0, 2, a, 1, tau, c, 1, w, 1, info);
        Chkxer("DORMTR", infot, nout, lerr, ok);
        infot = 12;
        Rormtr("R", "U", "N", 2, 0, a, 1, tau, c, 2, w, 1, info);
        Chkxer("DORMTR", infot, nout, lerr, ok);
        nt += 10;
        //
        // Rsptrd
        //
        srnamt = "DSPTRD";
        infot = 1;
        Rsptrd("/", 0, a, d, e, tau, info);
        Chkxer("DSPTRD", infot, nout, lerr, ok);
        infot = 2;
        Rsptrd("U", -1, a, d, e, tau, info);
        Chkxer("DSPTRD", infot, nout, lerr, ok);
        nt += 2;
        //
        // Ropgtr
        //
        srnamt = "DOPGTR";
        infot = 1;
        Ropgtr("/", 0, a, tau, z, 1, w, info);
        Chkxer("DOPGTR", infot, nout, lerr, ok);
        infot = 2;
        Ropgtr("U", -1, a, tau, z, 1, w, info);
        Chkxer("DOPGTR", infot, nout, lerr, ok);
        infot = 6;
        Ropgtr("U", 2, a, tau, z, 1, w, info);
        Chkxer("DOPGTR", infot, nout, lerr, ok);
        nt += 3;
        //
        // Ropmtr
        //
        srnamt = "DOPMTR";
        infot = 1;
        Ropmtr("/", "U", "N", 0, 0, a, tau, c, 1, w, info);
        Chkxer("DOPMTR", infot, nout, lerr, ok);
        infot = 2;
        Ropmtr("L", "/", "N", 0, 0, a, tau, c, 1, w, info);
        Chkxer("DOPMTR", infot, nout, lerr, ok);
        infot = 3;
        Ropmtr("L", "U", "/", 0, 0, a, tau, c, 1, w, info);
        Chkxer("DOPMTR", infot, nout, lerr, ok);
        infot = 4;
        Ropmtr("L", "U", "N", -1, 0, a, tau, c, 1, w, info);
        Chkxer("DOPMTR", infot, nout, lerr, ok);
        infot = 5;
        Ropmtr("L", "U", "N", 0, -1, a, tau, c, 1, w, info);
        Chkxer("DOPMTR", infot, nout, lerr, ok);
        infot = 9;
        Ropmtr("L", "U", "N", 2, 0, a, tau, c, 1, w, info);
        Chkxer("DOPMTR", infot, nout, lerr, ok);
        nt += 6;
        //
        // Rpteqr
        //
        srnamt = "DPTEQR";
        infot = 1;
        Rpteqr("/", 0, d, e, z, 1, w, info);
        Chkxer("DPTEQR", infot, nout, lerr, ok);
        infot = 2;
        Rpteqr("N", -1, d, e, z, 1, w, info);
        Chkxer("DPTEQR", infot, nout, lerr, ok);
        infot = 6;
        Rpteqr("V", 2, d, e, z, 1, w, info);
        Chkxer("DPTEQR", infot, nout, lerr, ok);
        nt += 3;
        //
        // Rstebz
        //
        srnamt = "DSTEBZ";
        infot = 1;
        Rstebz("/", "E", 0, 0.0, 1.0, 1, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        Chkxer("DSTEBZ", infot, nout, lerr, ok);
        infot = 2;
        Rstebz("A", "/", 0, 0.0, 0.0, 0, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        Chkxer("DSTEBZ", infot, nout, lerr, ok);
        infot = 3;
        Rstebz("A", "E", -1, 0.0, 0.0, 0, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        Chkxer("DSTEBZ", infot, nout, lerr, ok);
        infot = 5;
        Rstebz("V", "E", 0, 0.0, 0.0, 0, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        Chkxer("DSTEBZ", infot, nout, lerr, ok);
        infot = 6;
        Rstebz("I", "E", 0, 0.0, 0.0, 0, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        Chkxer("DSTEBZ", infot, nout, lerr, ok);
        infot = 6;
        Rstebz("I", "E", 1, 0.0, 0.0, 2, 1, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        Chkxer("DSTEBZ", infot, nout, lerr, ok);
        infot = 7;
        Rstebz("I", "E", 1, 0.0, 0.0, 1, 0, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        Chkxer("DSTEBZ", infot, nout, lerr, ok);
        infot = 7;
        Rstebz("I", "E", 1, 0.0, 0.0, 1, 2, 0.0, d, e, m, nsplit, x, i1, i2, w, iw, info);
        Chkxer("DSTEBZ", infot, nout, lerr, ok);
        nt += 8;
        //
        // Rstein
        //
        srnamt = "DSTEIN";
        infot = 1;
        Rstein(-1, d, e, 0, x, i1, i2, z, 1, w, iw, i3, info);
        Chkxer("DSTEIN", infot, nout, lerr, ok);
        infot = 4;
        Rstein(0, d, e, -1, x, i1, i2, z, 1, w, iw, i3, info);
        Chkxer("DSTEIN", infot, nout, lerr, ok);
        infot = 4;
        Rstein(0, d, e, 1, x, i1, i2, z, 1, w, iw, i3, info);
        Chkxer("DSTEIN", infot, nout, lerr, ok);
        infot = 9;
        Rstein(2, d, e, 0, x, i1, i2, z, 1, w, iw, i3, info);
        Chkxer("DSTEIN", infot, nout, lerr, ok);
        nt += 4;
        //
        // Rsteqr
        //
        srnamt = "DSTEQR";
        infot = 1;
        Rsteqr("/", 0, d, e, z, 1, w, info);
        Chkxer("DSTEQR", infot, nout, lerr, ok);
        infot = 2;
        Rsteqr("N", -1, d, e, z, 1, w, info);
        Chkxer("DSTEQR", infot, nout, lerr, ok);
        infot = 6;
        Rsteqr("V", 2, d, e, z, 1, w, info);
        Chkxer("DSTEQR", infot, nout, lerr, ok);
        nt += 3;
        //
        // Rsterf
        //
        srnamt = "DSTERF";
        infot = 1;
        Rsterf(-1, d, e, info);
        Chkxer("DSTERF", infot, nout, lerr, ok);
        nt++;
        //
        // Rstedc
        //
        srnamt = "DSTEDC";
        infot = 1;
        Rstedc("/", 0, d, e, z, 1, w, 1, iw, 1, info);
        Chkxer("DSTEDC", infot, nout, lerr, ok);
        infot = 2;
        Rstedc("N", -1, d, e, z, 1, w, 1, iw, 1, info);
        Chkxer("DSTEDC", infot, nout, lerr, ok);
        infot = 6;
        Rstedc("V", 2, d, e, z, 1, w, 23, iw, 28, info);
        Chkxer("DSTEDC", infot, nout, lerr, ok);
        infot = 8;
        Rstedc("N", 1, d, e, z, 1, w, 0, iw, 1, info);
        Chkxer("DSTEDC", infot, nout, lerr, ok);
        infot = 8;
        Rstedc("I", 2, d, e, z, 2, w, 0, iw, 12, info);
        Chkxer("DSTEDC", infot, nout, lerr, ok);
        infot = 8;
        Rstedc("V", 2, d, e, z, 2, w, 0, iw, 28, info);
        Chkxer("DSTEDC", infot, nout, lerr, ok);
        infot = 10;
        Rstedc("N", 1, d, e, z, 1, w, 1, iw, 0, info);
        Chkxer("DSTEDC", infot, nout, lerr, ok);
        infot = 10;
        Rstedc("I", 2, d, e, z, 2, w, 19, iw, 0, info);
        Chkxer("DSTEDC", infot, nout, lerr, ok);
        infot = 10;
        Rstedc("V", 2, d, e, z, 2, w, 23, iw, 0, info);
        Chkxer("DSTEDC", infot, nout, lerr, ok);
        nt += 9;
        //
        // Rstevd
        //
        srnamt = "DSTEVD";
        infot = 1;
        Rstevd("/", 0, d, e, z, 1, w, 1, iw, 1, info);
        Chkxer("DSTEVD", infot, nout, lerr, ok);
        infot = 2;
        Rstevd("N", -1, d, e, z, 1, w, 1, iw, 1, info);
        Chkxer("DSTEVD", infot, nout, lerr, ok);
        infot = 6;
        Rstevd("V", 2, d, e, z, 1, w, 19, iw, 12, info);
        Chkxer("DSTEVD", infot, nout, lerr, ok);
        infot = 8;
        Rstevd("N", 1, d, e, z, 1, w, 0, iw, 1, info);
        Chkxer("DSTEVD", infot, nout, lerr, ok);
        infot = 8;
        Rstevd("V", 2, d, e, z, 2, w, 12, iw, 12, info);
        Chkxer("DSTEVD", infot, nout, lerr, ok);
        infot = 10;
        Rstevd("N", 0, d, e, z, 1, w, 1, iw, 0, info);
        Chkxer("DSTEVD", infot, nout, lerr, ok);
        infot = 10;
        Rstevd("V", 2, d, e, z, 2, w, 19, iw, 11, info);
        Chkxer("DSTEVD", infot, nout, lerr, ok);
        nt += 7;
        //
        // Rstev
        //
        srnamt = "DSTEV ";
        infot = 1;
        Rstev("/", 0, d, e, z, 1, w, info);
        Chkxer("DSTEV ", infot, nout, lerr, ok);
        infot = 2;
        Rstev("N", -1, d, e, z, 1, w, info);
        Chkxer("DSTEV ", infot, nout, lerr, ok);
        infot = 6;
        Rstev("V", 2, d, e, z, 1, w, info);
        Chkxer("DSTEV ", infot, nout, lerr, ok);
        nt += 3;
        //
        // Rstevx
        //
        srnamt = "DSTEVX";
        infot = 1;
        Rstevx("/", "A", 0, d, e, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSTEVX", infot, nout, lerr, ok);
        infot = 2;
        Rstevx("N", "/", 0, d, e, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSTEVX", infot, nout, lerr, ok);
        infot = 3;
        Rstevx("N", "A", -1, d, e, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSTEVX", infot, nout, lerr, ok);
        infot = 7;
        Rstevx("N", "V", 1, d, e, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSTEVX", infot, nout, lerr, ok);
        infot = 8;
        Rstevx("N", "I", 1, d, e, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSTEVX", infot, nout, lerr, ok);
        infot = 8;
        Rstevx("N", "I", 1, d, e, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSTEVX", infot, nout, lerr, ok);
        infot = 9;
        Rstevx("N", "I", 2, d, e, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSTEVX", infot, nout, lerr, ok);
        infot = 9;
        Rstevx("N", "I", 1, d, e, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSTEVX", infot, nout, lerr, ok);
        infot = 14;
        Rstevx("V", "A", 2, d, e, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSTEVX", infot, nout, lerr, ok);
        nt += 9;
        //
        // Rstevr
        //
        n = 1;
        srnamt = "DSTEVR";
        infot = 1;
        Rstevr("/", "A", 0, d, e, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSTEVR", infot, nout, lerr, ok);
        infot = 2;
        Rstevr("V", "/", 0, d, e, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSTEVR", infot, nout, lerr, ok);
        infot = 3;
        Rstevr("V", "A", -1, d, e, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSTEVR", infot, nout, lerr, ok);
        infot = 7;
        Rstevr("V", "V", 1, d, e, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSTEVR", infot, nout, lerr, ok);
        infot = 8;
        Rstevr("V", "I", 1, d, e, 0.0, 0.0, 0, 1, 0.0, m, w, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSTEVR", infot, nout, lerr, ok);
        infot = 9;
        n = 2;
        Rstevr("V", "I", 2, d, e, 0.0, 0.0, 2, 1, 0.0, m, w, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSTEVR", infot, nout, lerr, ok);
        infot = 14;
        n = 1;
        Rstevr("V", "I", 1, d, e, 0.0, 0.0, 1, 1, 0.0, m, w, z, 0, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSTEVR", infot, nout, lerr, ok);
        infot = 17;
        Rstevr("V", "I", 1, d, e, 0.0, 0.0, 1, 1, 0.0, m, w, z, 1, iw, x, 20 * n - 1, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSTEVR", infot, nout, lerr, ok);
        infot = 19;
        Rstevr("V", "I", 1, d, e, 0.0, 0.0, 1, 1, 0.0, m, w, z, 1, iw, x, 20 * n, &iw[(2 * n + 1) - 1], 10 * n - 1, info);
        Chkxer("DSTEVR", infot, nout, lerr, ok);
        nt += 9;
        //
        // Rsyevd
        //
        srnamt = "DSYEVD";
        infot = 1;
        Rsyevd("/", "U", 0, a, 1, x, w, 1, iw, 1, info);
        Chkxer("DSYEVD", infot, nout, lerr, ok);
        infot = 2;
        Rsyevd("N", "/", 0, a, 1, x, w, 1, iw, 1, info);
        Chkxer("DSYEVD", infot, nout, lerr, ok);
        infot = 3;
        Rsyevd("N", "U", -1, a, 1, x, w, 1, iw, 1, info);
        Chkxer("DSYEVD", infot, nout, lerr, ok);
        infot = 5;
        Rsyevd("N", "U", 2, a, 1, x, w, 3, iw, 1, info);
        Chkxer("DSYEVD", infot, nout, lerr, ok);
        infot = 8;
        Rsyevd("N", "U", 1, a, 1, x, w, 0, iw, 1, info);
        Chkxer("DSYEVD", infot, nout, lerr, ok);
        infot = 8;
        Rsyevd("N", "U", 2, a, 2, x, w, 4, iw, 1, info);
        Chkxer("DSYEVD", infot, nout, lerr, ok);
        infot = 8;
        Rsyevd("V", "U", 2, a, 2, x, w, 20, iw, 12, info);
        Chkxer("DSYEVD", infot, nout, lerr, ok);
        infot = 10;
        Rsyevd("N", "U", 1, a, 1, x, w, 1, iw, 0, info);
        Chkxer("DSYEVD", infot, nout, lerr, ok);
        infot = 10;
        Rsyevd("N", "U", 2, a, 2, x, w, 5, iw, 0, info);
        Chkxer("DSYEVD", infot, nout, lerr, ok);
        infot = 10;
        Rsyevd("V", "U", 2, a, 2, x, w, 27, iw, 11, info);
        Chkxer("DSYEVD", infot, nout, lerr, ok);
        nt += 10;
        //
        // Rsyevd_2stage
        //
        srnamt = "DSYEVD_2STAGE";
        infot = 1;
        Rsyevd_2stage("/", "U", 0, a, 1, x, w, 1, iw, 1, info);
        Chkxer("DSYEVD_2STAGE", infot, nout, lerr, ok);
        infot = 1;
        Rsyevd_2stage("V", "U", 0, a, 1, x, w, 1, iw, 1, info);
        Chkxer("DSYEVD_2STAGE", infot, nout, lerr, ok);
        infot = 2;
        Rsyevd_2stage("N", "/", 0, a, 1, x, w, 1, iw, 1, info);
        Chkxer("DSYEVD_2STAGE", infot, nout, lerr, ok);
        infot = 3;
        Rsyevd_2stage("N", "U", -1, a, 1, x, w, 1, iw, 1, info);
        Chkxer("DSYEVD_2STAGE", infot, nout, lerr, ok);
        infot = 5;
        Rsyevd_2stage("N", "U", 2, a, 1, x, w, 3, iw, 1, info);
        Chkxer("DSYEVD_2STAGE", infot, nout, lerr, ok);
        infot = 8;
        Rsyevd_2stage("N", "U", 1, a, 1, x, w, 0, iw, 1, info);
        Chkxer("DSYEVD_2STAGE", infot, nout, lerr, ok);
        infot = 8;
        Rsyevd_2stage("N", "U", 2, a, 2, x, w, 4, iw, 1, info);
        Chkxer("DSYEVD_2STAGE", infot, nout, lerr, ok);
        // INFOT = 8
        // CALL Rsyevd_2stage( 'V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO )
        // CALL Chkxer( 'Rsyevd_2stage', INFOT, NOUT, LERR, OK )
        infot = 10;
        Rsyevd_2stage("N", "U", 1, a, 1, x, w, 1, iw, 0, info);
        Chkxer("DSYEVD_2STAGE", infot, nout, lerr, ok);
        infot = 10;
        Rsyevd_2stage("N", "U", 2, a, 2, x, w, 25, iw, 0, info);
        Chkxer("DSYEVD_2STAGE", infot, nout, lerr, ok);
        // INFOT = 10
        // CALL Rsyevd_2stage( 'V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO )
        // CALL Chkxer( 'Rsyevd_2stage', INFOT, NOUT, LERR, OK )
        nt += 9;
        //
        // Rsyevr
        //
        srnamt = "DSYEVR";
        n = 1;
        infot = 1;
        Rsyevr("/", "A", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        infot = 2;
        Rsyevr("V", "/", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        infot = 3;
        Rsyevr("V", "A", "/", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        infot = 4;
        Rsyevr("V", "A", "U", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        infot = 6;
        Rsyevr("V", "A", "U", 2, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        infot = 8;
        Rsyevr("V", "V", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        infot = 9;
        Rsyevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 0, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        infot = 10;
        //
        Rsyevr("V", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        infot = 15;
        Rsyevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 0, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        infot = 18;
        Rsyevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n - 1, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        infot = 20;
        Rsyevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n - 1, info);
        Chkxer("DSYEVR", infot, nout, lerr, ok);
        nt += 11;
        //
        // Rsyevr_2stage
        //
        srnamt = "DSYEVR_2STAGE";
        n = 1;
        infot = 1;
        Rsyevr_2stage("/", "A", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 1;
        Rsyevr_2stage("V", "A", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 2;
        Rsyevr_2stage("N", "/", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 3;
        Rsyevr_2stage("N", "A", "/", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 4;
        Rsyevr_2stage("N", "A", "U", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 6;
        Rsyevr_2stage("N", "A", "U", 2, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 8;
        Rsyevr_2stage("N", "V", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 9;
        Rsyevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 0, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 10;
        Rsyevr_2stage("N", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 15;
        Rsyevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 0, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 18;
        Rsyevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 0, &iw[(2 * n + 1) - 1], 10 * n, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        infot = 20;
        Rsyevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, &iw[(2 * n + 1) - 1], 0, info);
        Chkxer("DSYEVR_2STAGE", infot, nout, lerr, ok);
        nt += 12;
        //
        // Rsyev
        //
        srnamt = "DSYEV ";
        infot = 1;
        Rsyev("/", "U", 0, a, 1, x, w, 1, info);
        Chkxer("DSYEV ", infot, nout, lerr, ok);
        infot = 2;
        Rsyev("N", "/", 0, a, 1, x, w, 1, info);
        Chkxer("DSYEV ", infot, nout, lerr, ok);
        infot = 3;
        Rsyev("N", "U", -1, a, 1, x, w, 1, info);
        Chkxer("DSYEV ", infot, nout, lerr, ok);
        infot = 5;
        Rsyev("N", "U", 2, a, 1, x, w, 3, info);
        Chkxer("DSYEV ", infot, nout, lerr, ok);
        infot = 8;
        Rsyev("N", "U", 1, a, 1, x, w, 1, info);
        Chkxer("DSYEV ", infot, nout, lerr, ok);
        nt += 5;
        //
        // Rsyev_2stage
        //
        srnamt = "DSYEV_2STAGE ";
        infot = 1;
        Rsyev_2stage("/", "U", 0, a, 1, x, w, 1, info);
        Chkxer("DSYEV_2STAGE ", infot, nout, lerr, ok);
        infot = 1;
        Rsyev_2stage("V", "U", 0, a, 1, x, w, 1, info);
        Chkxer("DSYEV_2STAGE ", infot, nout, lerr, ok);
        infot = 2;
        Rsyev_2stage("N", "/", 0, a, 1, x, w, 1, info);
        Chkxer("DSYEV_2STAGE ", infot, nout, lerr, ok);
        infot = 3;
        Rsyev_2stage("N", "U", -1, a, 1, x, w, 1, info);
        Chkxer("DSYEV_2STAGE ", infot, nout, lerr, ok);
        infot = 5;
        Rsyev_2stage("N", "U", 2, a, 1, x, w, 3, info);
        Chkxer("DSYEV_2STAGE ", infot, nout, lerr, ok);
        infot = 8;
        Rsyev_2stage("N", "U", 1, a, 1, x, w, 1, info);
        Chkxer("DSYEV_2STAGE ", infot, nout, lerr, ok);
        nt += 6;
        //
        // Rsyevx
        //
        srnamt = "DSYEVX";
        infot = 1;
        Rsyevx("/", "A", "U", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        infot = 2;
        Rsyevx("N", "/", "U", 0, a, 1, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        infot = 3;
        Rsyevx("N", "A", "/", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        infot = 4;
        Rsyevx("N", "A", "U", -1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        infot = 6;
        Rsyevx("N", "A", "U", 2, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 16, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        infot = 8;
        Rsyevx("N", "V", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        infot = 9;
        Rsyevx("N", "I", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        infot = 9;
        Rsyevx("N", "I", "U", 1, a, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        infot = 10;
        Rsyevx("N", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 16, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        infot = 10;
        Rsyevx("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        infot = 15;
        Rsyevx("V", "A", "U", 2, a, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 16, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        infot = 17;
        Rsyevx("V", "A", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSYEVX", infot, nout, lerr, ok);
        nt += 12;
        //
        // Rsyevx_2stage
        //
        srnamt = "DSYEVX_2STAGE";
        infot = 1;
        Rsyevx_2stage("/", "A", "U", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 1;
        Rsyevx_2stage("V", "A", "U", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 2;
        Rsyevx_2stage("N", "/", "U", 0, a, 1, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 3;
        Rsyevx_2stage("N", "A", "/", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        infot = 4;
        Rsyevx_2stage("N", "A", "U", -1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 6;
        Rsyevx_2stage("N", "A", "U", 2, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 16, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 8;
        Rsyevx_2stage("N", "V", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 9;
        Rsyevx_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 9;
        Rsyevx_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 10;
        Rsyevx_2stage("N", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 16, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 10;
        Rsyevx_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, 8, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 15;
        Rsyevx_2stage("N", "A", "U", 2, a, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 0, w, 16, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        infot = 17;
        Rsyevx_2stage("N", "A", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSYEVX_2STAGE", infot, nout, lerr, ok);
        nt += 13;
        //
        // Rspevd
        //
        srnamt = "DSPEVD";
        infot = 1;
        Rspevd("/", "U", 0, a, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSPEVD", infot, nout, lerr, ok);
        infot = 2;
        Rspevd("N", "/", 0, a, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSPEVD", infot, nout, lerr, ok);
        infot = 3;
        Rspevd("N", "U", -1, a, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSPEVD", infot, nout, lerr, ok);
        infot = 7;
        Rspevd("V", "U", 2, a, x, z, 1, w, 23, iw, 12, info);
        Chkxer("DSPEVD", infot, nout, lerr, ok);
        infot = 9;
        Rspevd("N", "U", 1, a, x, z, 1, w, 0, iw, 1, info);
        Chkxer("DSPEVD", infot, nout, lerr, ok);
        infot = 9;
        Rspevd("N", "U", 2, a, x, z, 1, w, 3, iw, 1, info);
        Chkxer("DSPEVD", infot, nout, lerr, ok);
        infot = 9;
        Rspevd("V", "U", 2, a, x, z, 2, w, 16, iw, 12, info);
        Chkxer("DSPEVD", infot, nout, lerr, ok);
        infot = 11;
        Rspevd("N", "U", 1, a, x, z, 1, w, 1, iw, 0, info);
        Chkxer("DSPEVD", infot, nout, lerr, ok);
        infot = 11;
        Rspevd("N", "U", 2, a, x, z, 1, w, 4, iw, 0, info);
        Chkxer("DSPEVD", infot, nout, lerr, ok);
        infot = 11;
        Rspevd("V", "U", 2, a, x, z, 2, w, 23, iw, 11, info);
        Chkxer("DSPEVD", infot, nout, lerr, ok);
        nt += 10;
        //
        // Rspev
        //
        srnamt = "DSPEV ";
        infot = 1;
        Rspev("/", "U", 0, a, w, z, 1, x, info);
        Chkxer("DSPEV ", infot, nout, lerr, ok);
        infot = 2;
        Rspev("N", "/", 0, a, w, z, 1, x, info);
        Chkxer("DSPEV ", infot, nout, lerr, ok);
        infot = 3;
        Rspev("N", "U", -1, a, w, z, 1, x, info);
        Chkxer("DSPEV ", infot, nout, lerr, ok);
        infot = 7;
        Rspev("V", "U", 2, a, w, z, 1, x, info);
        Chkxer("DSPEV ", infot, nout, lerr, ok);
        nt += 4;
        //
        // Rspevx
        //
        srnamt = "DSPEVX";
        infot = 1;
        Rspevx("/", "A", "U", 0, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSPEVX", infot, nout, lerr, ok);
        infot = 2;
        Rspevx("N", "/", "U", 0, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSPEVX", infot, nout, lerr, ok);
        infot = 3;
        Rspevx("N", "A", "/", 0, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        infot = 4;
        Rspevx("N", "A", "U", -1, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSPEVX", infot, nout, lerr, ok);
        infot = 7;
        Rspevx("N", "V", "U", 1, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSPEVX", infot, nout, lerr, ok);
        infot = 8;
        Rspevx("N", "I", "U", 1, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSPEVX", infot, nout, lerr, ok);
        infot = 8;
        Rspevx("N", "I", "U", 1, a, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSPEVX", infot, nout, lerr, ok);
        infot = 9;
        Rspevx("N", "I", "U", 2, a, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSPEVX", infot, nout, lerr, ok);
        infot = 9;
        Rspevx("N", "I", "U", 1, a, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSPEVX", infot, nout, lerr, ok);
        infot = 14;
        Rspevx("V", "A", "U", 2, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSPEVX", infot, nout, lerr, ok);
        nt += 10;
        //
        // Test error exits for the SB path.
        //
    } else if (Mlsamen(2, c2.elems, "SB")) {
        //
        // Rsbtrd
        //
        srnamt = "DSBTRD";
        infot = 1;
        Rsbtrd("/", "U", 0, 0, a, 1, d, e, z, 1, w, info);
        Chkxer("DSBTRD", infot, nout, lerr, ok);
        infot = 2;
        Rsbtrd("N", "/", 0, 0, a, 1, d, e, z, 1, w, info);
        Chkxer("DSBTRD", infot, nout, lerr, ok);
        infot = 3;
        Rsbtrd("N", "U", -1, 0, a, 1, d, e, z, 1, w, info);
        Chkxer("DSBTRD", infot, nout, lerr, ok);
        infot = 4;
        Rsbtrd("N", "U", 0, -1, a, 1, d, e, z, 1, w, info);
        Chkxer("DSBTRD", infot, nout, lerr, ok);
        infot = 6;
        Rsbtrd("N", "U", 1, 1, a, 1, d, e, z, 1, w, info);
        Chkxer("DSBTRD", infot, nout, lerr, ok);
        infot = 10;
        Rsbtrd("V", "U", 2, 0, a, 1, d, e, z, 1, w, info);
        Chkxer("DSBTRD", infot, nout, lerr, ok);
        nt += 6;
        //
        // Rsytrd_sb2st
        //
        srnamt = "DSYTRD_SB2ST";
        infot = 1;
        Rsytrd_sb2st("/", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_sb2st("N", "/", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 2;
        Rsytrd_sb2st("N", "H", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 3;
        Rsytrd_sb2st("N", "N", "/", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 4;
        Rsytrd_sb2st("N", "N", "U", -1, 0, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 5;
        Rsytrd_sb2st("N", "N", "U", 0, -1, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 7;
        Rsytrd_sb2st("N", "N", "U", 0, 1, a, 1, d, e, c, 1, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 11;
        Rsytrd_sb2st("N", "N", "U", 0, 0, a, 1, d, e, c, 0, w, 1, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        infot = 13;
        Rsytrd_sb2st("N", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 0, info);
        Chkxer("DSYTRD_SB2ST", infot, nout, lerr, ok);
        nt += 9;
        //
        // Rsbevd
        //
        srnamt = "DSBEVD";
        infot = 1;
        Rsbevd("/", "U", 0, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        infot = 2;
        Rsbevd("N", "/", 0, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        infot = 3;
        Rsbevd("N", "U", -1, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        infot = 4;
        Rsbevd("N", "U", 0, -1, a, 1, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        infot = 6;
        Rsbevd("N", "U", 2, 1, a, 1, x, z, 1, w, 4, iw, 1, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        infot = 9;
        Rsbevd("V", "U", 2, 1, a, 2, x, z, 1, w, 25, iw, 12, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        infot = 11;
        Rsbevd("N", "U", 1, 0, a, 1, x, z, 1, w, 0, iw, 1, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        infot = 11;
        Rsbevd("N", "U", 2, 0, a, 1, x, z, 1, w, 3, iw, 1, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        infot = 11;
        Rsbevd("V", "U", 2, 0, a, 1, x, z, 2, w, 18, iw, 12, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        infot = 13;
        Rsbevd("N", "U", 1, 0, a, 1, x, z, 1, w, 1, iw, 0, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        infot = 13;
        Rsbevd("V", "U", 2, 0, a, 1, x, z, 2, w, 25, iw, 11, info);
        Chkxer("DSBEVD", infot, nout, lerr, ok);
        nt += 11;
        //
        // Rsbevd_2stage
        //
        srnamt = "DSBEVD_2STAGE";
        infot = 1;
        Rsbevd_2stage("/", "U", 0, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSBEVD_2STAGE", infot, nout, lerr, ok);
        infot = 1;
        Rsbevd_2stage("V", "U", 0, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSBEVD_2STAGE", infot, nout, lerr, ok);
        infot = 2;
        Rsbevd_2stage("N", "/", 0, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSBEVD_2STAGE", infot, nout, lerr, ok);
        infot = 3;
        Rsbevd_2stage("N", "U", -1, 0, a, 1, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSBEVD_2STAGE", infot, nout, lerr, ok);
        infot = 4;
        Rsbevd_2stage("N", "U", 0, -1, a, 1, x, z, 1, w, 1, iw, 1, info);
        Chkxer("DSBEVD_2STAGE", infot, nout, lerr, ok);
        infot = 6;
        Rsbevd_2stage("N", "U", 2, 1, a, 1, x, z, 1, w, 4, iw, 1, info);
        Chkxer("DSBEVD_2STAGE", infot, nout, lerr, ok);
        // INFOT = 9
        // CALL Rsbevd_2stage( 'V', 'U', 2, 1, A, 2, X, Z, 1, W,
        // $                                      25, IW, 12, INFO )
        // CALL Chkxer( 'Rsbevd_2stage', INFOT, NOUT, LERR, OK )
        infot = 11;
        Rsbevd_2stage("N", "U", 1, 0, a, 1, x, z, 1, w, 0, iw, 1, info);
        Chkxer("DSBEVD_2STAGE", infot, nout, lerr, ok);
        infot = 11;
        Rsbevd_2stage("N", "U", 2, 0, a, 1, x, z, 1, w, 3, iw, 1, info);
        Chkxer("DSBEVD_2STAGE", infot, nout, lerr, ok);
        // INFOT = 11
        // CALL Rsbevd_2stage( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
        // $                                      18, IW, 12, INFO )
        // CALL Chkxer( 'Rsbevd_2stage', INFOT, NOUT, LERR, OK )
        infot = 13;
        Rsbevd_2stage("N", "U", 1, 0, a, 1, x, z, 1, w, 1, iw, 0, info);
        Chkxer("DSBEVD_2STAGE", infot, nout, lerr, ok);
        // INFOT = 13
        // CALL Rsbevd_2stage( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
        // $                                      25, IW, 11, INFO )
        // CALL Chkxer( 'Rsbevd_2stage', INFOT, NOUT, LERR, OK )
        // NT = NT + 12
        nt += 9;
        //
        // Rsbev
        //
        srnamt = "DSBEV ";
        infot = 1;
        Rsbev("/", "U", 0, 0, a, 1, x, z, 1, w, info);
        Chkxer("DSBEV ", infot, nout, lerr, ok);
        infot = 2;
        Rsbev("N", "/", 0, 0, a, 1, x, z, 1, w, info);
        Chkxer("DSBEV ", infot, nout, lerr, ok);
        infot = 3;
        Rsbev("N", "U", -1, 0, a, 1, x, z, 1, w, info);
        Chkxer("DSBEV ", infot, nout, lerr, ok);
        infot = 4;
        Rsbev("N", "U", 0, -1, a, 1, x, z, 1, w, info);
        Chkxer("DSBEV ", infot, nout, lerr, ok);
        infot = 6;
        Rsbev("N", "U", 2, 1, a, 1, x, z, 1, w, info);
        Chkxer("DSBEV ", infot, nout, lerr, ok);
        infot = 9;
        Rsbev("V", "U", 2, 0, a, 1, x, z, 1, w, info);
        Chkxer("DSBEV ", infot, nout, lerr, ok);
        nt += 6;
        //
        // Rsbev_2stage
        //
        srnamt = "DSBEV_2STAGE ";
        infot = 1;
        Rsbev_2stage("/", "U", 0, 0, a, 1, x, z, 1, w, 0, info);
        Chkxer("DSBEV_2STAGE ", infot, nout, lerr, ok);
        infot = 1;
        Rsbev_2stage("V", "U", 0, 0, a, 1, x, z, 1, w, 0, info);
        Chkxer("DSBEV_2STAGE ", infot, nout, lerr, ok);
        infot = 2;
        Rsbev_2stage("N", "/", 0, 0, a, 1, x, z, 1, w, 0, info);
        Chkxer("DSBEV_2STAGE ", infot, nout, lerr, ok);
        infot = 3;
        Rsbev_2stage("N", "U", -1, 0, a, 1, x, z, 1, w, 0, info);
        Chkxer("DSBEV_2STAGE ", infot, nout, lerr, ok);
        infot = 4;
        Rsbev_2stage("N", "U", 0, -1, a, 1, x, z, 1, w, 0, info);
        Chkxer("DSBEV_2STAGE ", infot, nout, lerr, ok);
        infot = 6;
        Rsbev_2stage("N", "U", 2, 1, a, 1, x, z, 1, w, 0, info);
        Chkxer("DSBEV_2STAGE ", infot, nout, lerr, ok);
        infot = 9;
        Rsbev_2stage("N", "U", 2, 0, a, 1, x, z, 0, w, 0, info);
        Chkxer("DSBEV_2STAGE ", infot, nout, lerr, ok);
        infot = 11;
        Rsbev_2stage("N", "U", 0, 0, a, 1, x, z, 1, w, 0, info);
        Chkxer("DSBEV_2STAGE ", infot, nout, lerr, ok);
        nt += 8;
        //
        // Rsbevx
        //
        srnamt = "DSBEVX";
        infot = 1;
        Rsbevx("/", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 2;
        Rsbevx("N", "/", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 3;
        Rsbevx("N", "A", "/", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 4;
        Rsbevx("N", "A", "U", -1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 5;
        Rsbevx("N", "A", "U", 0, -1, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 7;
        Rsbevx("N", "A", "U", 2, 1, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 9;
        Rsbevx("V", "A", "U", 2, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 2, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 11;
        Rsbevx("N", "V", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 12;
        Rsbevx("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 12;
        Rsbevx("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 13;
        Rsbevx("N", "I", "U", 2, 0, a, 1, q, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 13;
        Rsbevx("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        infot = 18;
        Rsbevx("V", "A", "U", 2, 0, a, 1, q, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, iw, i3, info);
        Chkxer("DSBEVX", infot, nout, lerr, ok);
        nt += 13;
        //
        // Rsbevx_2stage
        //
        srnamt = "DSBEVX_2STAGE";
        infot = 1;
        Rsbevx_2stage("/", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        infot = 1;
        Rsbevx_2stage("V", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        infot = 2;
        Rsbevx_2stage("N", "/", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        infot = 3;
        Rsbevx_2stage("N", "A", "/", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        infot = 4;
        Rsbevx_2stage("N", "A", "U", -1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        infot = 5;
        Rsbevx_2stage("N", "A", "U", 0, -1, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        infot = 7;
        Rsbevx_2stage("N", "A", "U", 2, 1, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        // INFOT = 9
        // CALL Rsbevx_2stage( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0D0,
        // $          0.0D0, 0, 0, 0.0D0, M, X, Z, 2, W, 0, IW, I3, INFO )
        // CALL Chkxer( 'Rsbevx_2stage', INFOT, NOUT, LERR, OK )
        infot = 11;
        Rsbevx_2stage("N", "V", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        infot = 12;
        Rsbevx_2stage("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        infot = 12;
        Rsbevx_2stage("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        infot = 13;
        Rsbevx_2stage("N", "I", "U", 2, 0, a, 1, q, 1, 0.0, 0.0, 2, 1, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        infot = 13;
        Rsbevx_2stage("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        // INFOT = 18
        // CALL Rsbevx_2stage( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0D0,
        // $          0.0D0, 0, 0, 0.0D0, M, X, Z, 1, W, 0, IW, I3, INFO )
        // CALL Chkxer( 'Rsbevx_2stage', INFOT, NOUT, LERR, OK )
        infot = 20;
        Rsbevx_2stage("N", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, iw, i3, info);
        Chkxer("DSBEVX_2STAGE", infot, nout, lerr, ok);
        // NT = NT + 15
        nt += 13;
    }
    //
    // Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a3,' routines passed the tests of the error exits',' (',i3,"
                    "' tests done)')"),
            path, nt;
    } else {
        write(nout, "(' *** ',a3,' routines failed the tests of the error ','exits ***')"), path;
    }
    //
    // End of Rerrst
    //
}
