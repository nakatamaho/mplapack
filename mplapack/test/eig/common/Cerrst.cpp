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

void Cerrst(const char *path, INTEGER const nunit) {
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
    const INTEGER nmax = 3;
    INTEGER i = 0;
    arr_2d<nmax, nmax, COMPLEX> a;
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / (i + j).real();
        }
    }
    arr_1d<nmax, REAL> d;
    arr_1d<nmax, REAL> e;
    arr_1d<nmax, int> i1;
    arr_1d<nmax, int> i2;
    arr_1d<nmax, COMPLEX> tau;
    for (j = 1; j <= nmax; j = j + 1) {
        d[j - 1] = j.real();
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
    arr_1d<lw, COMPLEX> w;
    INTEGER info = 0;
    arr_2d<nmax, nmax, COMPLEX> c;
    arr_2d<nmax, nmax, COMPLEX> z;
    arr_1d<lw, REAL> rw;
    arr_1d<nmax, REAL> x;
    const INTEGER liw = 12 * nmax;
    arr_1d<liw, int> iw;
    arr_1d<nmax, int> i3;
    INTEGER m = 0;
    INTEGER n = 0;
    arr_1d<lw, REAL> r;
    arr_2d<nmax, nmax, COMPLEX> q;
    if (Mlsamen(2, c2, "ST")) {
        //
        //        Chetrd
        //
        srnamt = "Chetrd";
        infot = 1;
        Chetrd("/", 0, a, 1, d, e, tau, w, 1, info);
        chkxer("Chetrd", infot, nout, lerr, ok);
        infot = 2;
        Chetrd("U", -1, a, 1, d, e, tau, w, 1, info);
        chkxer("Chetrd", infot, nout, lerr, ok);
        infot = 4;
        Chetrd("U", 2, a, 1, d, e, tau, w, 1, info);
        chkxer("Chetrd", infot, nout, lerr, ok);
        infot = 9;
        Chetrd("U", 0, a, 1, d, e, tau, w, 0, info);
        chkxer("Chetrd", infot, nout, lerr, ok);
        nt += 4;
        //
        //        Chetrd_2stage
        //
        srnamt = "Chetrd_2stage";
        infot = 1;
        Chetrd_2stage("/", "U", 0, a, 1, d, e, tau, c, 1, w, 1, info);
        chkxer("Chetrd_2stage", infot, nout, lerr, ok);
        infot = 1;
        Chetrd_2stage("H", "U", 0, a, 1, d, e, tau, c, 1, w, 1, info);
        chkxer("Chetrd_2stage", infot, nout, lerr, ok);
        infot = 2;
        Chetrd_2stage("N", "/", 0, a, 1, d, e, tau, c, 1, w, 1, info);
        chkxer("Chetrd_2stage", infot, nout, lerr, ok);
        infot = 3;
        Chetrd_2stage("N", "U", -1, a, 1, d, e, tau, c, 1, w, 1, info);
        chkxer("Chetrd_2stage", infot, nout, lerr, ok);
        infot = 5;
        Chetrd_2stage("N", "U", 2, a, 1, d, e, tau, c, 1, w, 1, info);
        chkxer("Chetrd_2stage", infot, nout, lerr, ok);
        infot = 10;
        Chetrd_2stage("N", "U", 0, a, 1, d, e, tau, c, 0, w, 1, info);
        chkxer("Chetrd_2stage", infot, nout, lerr, ok);
        infot = 12;
        Chetrd_2stage("N", "U", 0, a, 1, d, e, tau, c, 1, w, 0, info);
        chkxer("Chetrd_2stage", infot, nout, lerr, ok);
        nt += 7;
        //
        //        Chetrd_he2hb
        //
        srnamt = "Chetrd_he2hb";
        infot = 1;
        Chetrd_he2hb("/", 0, 0, a, 1, c, 1, tau, w, 1, info);
        chkxer("Chetrd_he2hb", infot, nout, lerr, ok);
        infot = 2;
        Chetrd_he2hb("U", -1, 0, a, 1, c, 1, tau, w, 1, info);
        chkxer("Chetrd_he2hb", infot, nout, lerr, ok);
        infot = 3;
        Chetrd_he2hb("U", 0, -1, a, 1, c, 1, tau, w, 1, info);
        chkxer("Chetrd_he2hb", infot, nout, lerr, ok);
        infot = 5;
        Chetrd_he2hb("U", 2, 0, a, 1, c, 1, tau, w, 1, info);
        chkxer("Chetrd_he2hb", infot, nout, lerr, ok);
        infot = 7;
        Chetrd_he2hb("U", 0, 2, a, 1, c, 1, tau, w, 1, info);
        chkxer("Chetrd_he2hb", infot, nout, lerr, ok);
        infot = 10;
        Chetrd_he2hb("U", 0, 0, a, 1, c, 1, tau, w, 0, info);
        chkxer("Chetrd_he2hb", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Chetrd_HB2ST
        //
        srnamt = "Chetrd_HB2ST";
        infot = 1;
        Chetrd_hb2st("/", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 2;
        Chetrd_hb2st("Y", "/", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 2;
        Chetrd_hb2st("Y", "H", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 3;
        Chetrd_hb2st("Y", "N", "/", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 4;
        Chetrd_hb2st("Y", "N", "U", -1, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 5;
        Chetrd_hb2st("Y", "N", "U", 0, -1, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 7;
        Chetrd_hb2st("Y", "N", "U", 0, 1, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 11;
        Chetrd_hb2st("Y", "N", "U", 0, 0, a, 1, d, e, c, 0, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 13;
        Chetrd_hb2st("Y", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 0, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Cungtr
        //
        srnamt = "Cungtr";
        infot = 1;
        Cungtr("/", 0, a, 1, tau, w, 1, info);
        chkxer("Cungtr", infot, nout, lerr, ok);
        infot = 2;
        Cungtr("U", -1, a, 1, tau, w, 1, info);
        chkxer("Cungtr", infot, nout, lerr, ok);
        infot = 4;
        Cungtr("U", 2, a, 1, tau, w, 1, info);
        chkxer("Cungtr", infot, nout, lerr, ok);
        infot = 7;
        Cungtr("U", 3, a, 3, tau, w, 1, info);
        chkxer("Cungtr", infot, nout, lerr, ok);
        nt += 4;
        //
        //        Cunmtr
        //
        srnamt = "Cunmtr";
        infot = 1;
        Cunmtr("/", "U", "N", 0, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmtr", infot, nout, lerr, ok);
        infot = 2;
        Cunmtr("L", "/", "N", 0, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmtr", infot, nout, lerr, ok);
        infot = 3;
        Cunmtr("L", "U", "/", 0, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmtr", infot, nout, lerr, ok);
        infot = 4;
        Cunmtr("L", "U", "N", -1, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmtr", infot, nout, lerr, ok);
        infot = 5;
        Cunmtr("L", "U", "N", 0, -1, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmtr", infot, nout, lerr, ok);
        infot = 7;
        Cunmtr("L", "U", "N", 2, 0, a, 1, tau, c, 2, w, 1, info);
        chkxer("Cunmtr", infot, nout, lerr, ok);
        infot = 7;
        Cunmtr("R", "U", "N", 0, 2, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmtr", infot, nout, lerr, ok);
        infot = 10;
        Cunmtr("L", "U", "N", 2, 0, a, 2, tau, c, 1, w, 1, info);
        chkxer("Cunmtr", infot, nout, lerr, ok);
        infot = 12;
        Cunmtr("L", "U", "N", 0, 2, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmtr", infot, nout, lerr, ok);
        infot = 12;
        Cunmtr("R", "U", "N", 2, 0, a, 1, tau, c, 2, w, 1, info);
        chkxer("Cunmtr", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Chptrd
        //
        srnamt = "Chptrd";
        infot = 1;
        Chptrd("/", 0, a, d, e, tau, info);
        chkxer("Chptrd", infot, nout, lerr, ok);
        infot = 2;
        Chptrd("U", -1, a, d, e, tau, info);
        chkxer("Chptrd", infot, nout, lerr, ok);
        nt += 2;
        //
        //        Cupgtr
        //
        srnamt = "Cupgtr";
        infot = 1;
        Cupgtr("/", 0, a, tau, z, 1, w, info);
        chkxer("Cupgtr", infot, nout, lerr, ok);
        infot = 2;
        Cupgtr("U", -1, a, tau, z, 1, w, info);
        chkxer("Cupgtr", infot, nout, lerr, ok);
        infot = 6;
        Cupgtr("U", 2, a, tau, z, 1, w, info);
        chkxer("Cupgtr", infot, nout, lerr, ok);
        nt += 3;
        //
        //        Cupmtr
        //
        srnamt = "Cupmtr";
        infot = 1;
        Cupmtr("/", "U", "N", 0, 0, a, tau, c, 1, w, info);
        chkxer("Cupmtr", infot, nout, lerr, ok);
        infot = 2;
        Cupmtr("L", "/", "N", 0, 0, a, tau, c, 1, w, info);
        chkxer("Cupmtr", infot, nout, lerr, ok);
        infot = 3;
        Cupmtr("L", "U", "/", 0, 0, a, tau, c, 1, w, info);
        chkxer("Cupmtr", infot, nout, lerr, ok);
        infot = 4;
        Cupmtr("L", "U", "N", -1, 0, a, tau, c, 1, w, info);
        chkxer("Cupmtr", infot, nout, lerr, ok);
        infot = 5;
        Cupmtr("L", "U", "N", 0, -1, a, tau, c, 1, w, info);
        chkxer("Cupmtr", infot, nout, lerr, ok);
        infot = 9;
        Cupmtr("L", "U", "N", 2, 0, a, tau, c, 1, w, info);
        chkxer("Cupmtr", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Cpteqr
        //
        srnamt = "Cpteqr";
        infot = 1;
        Cpteqr("/", 0, d, e, z, 1, rw, info);
        chkxer("Cpteqr", infot, nout, lerr, ok);
        infot = 2;
        Cpteqr("N", -1, d, e, z, 1, rw, info);
        chkxer("Cpteqr", infot, nout, lerr, ok);
        infot = 6;
        Cpteqr("V", 2, d, e, z, 1, rw, info);
        chkxer("Cpteqr", infot, nout, lerr, ok);
        nt += 3;
        //
        //        Cstein
        //
        srnamt = "Cstein";
        infot = 1;
        Cstein(-1, d, e, 0, x, i1, i2, z, 1, rw, iw, i3, info);
        chkxer("Cstein", infot, nout, lerr, ok);
        infot = 4;
        Cstein(0, d, e, -1, x, i1, i2, z, 1, rw, iw, i3, info);
        chkxer("Cstein", infot, nout, lerr, ok);
        infot = 4;
        Cstein(0, d, e, 1, x, i1, i2, z, 1, rw, iw, i3, info);
        chkxer("Cstein", infot, nout, lerr, ok);
        infot = 9;
        Cstein(2, d, e, 0, x, i1, i2, z, 1, rw, iw, i3, info);
        chkxer("Cstein", infot, nout, lerr, ok);
        nt += 4;
        //
        //        Csteqr
        //
        srnamt = "Csteqr";
        infot = 1;
        Csteqr("/", 0, d, e, z, 1, rw, info);
        chkxer("Csteqr", infot, nout, lerr, ok);
        infot = 2;
        Csteqr("N", -1, d, e, z, 1, rw, info);
        chkxer("Csteqr", infot, nout, lerr, ok);
        infot = 6;
        Csteqr("V", 2, d, e, z, 1, rw, info);
        chkxer("Csteqr", infot, nout, lerr, ok);
        nt += 3;
        //
        //        Cstedc
        //
        srnamt = "Cstedc";
        infot = 1;
        Cstedc("/", 0, d, e, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        infot = 2;
        Cstedc("N", -1, d, e, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        infot = 6;
        Cstedc("V", 2, d, e, z, 1, w, 4, rw, 23, iw, 28, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        infot = 8;
        Cstedc("N", 2, d, e, z, 1, w, 0, rw, 1, iw, 1, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        infot = 8;
        Cstedc("V", 2, d, e, z, 2, w, 0, rw, 23, iw, 28, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        infot = 10;
        Cstedc("N", 2, d, e, z, 1, w, 1, rw, 0, iw, 1, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        infot = 10;
        Cstedc("I", 2, d, e, z, 2, w, 1, rw, 1, iw, 12, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        infot = 10;
        Cstedc("V", 2, d, e, z, 2, w, 4, rw, 1, iw, 28, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        infot = 12;
        Cstedc("N", 2, d, e, z, 1, w, 1, rw, 1, iw, 0, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        infot = 12;
        Cstedc("I", 2, d, e, z, 2, w, 1, rw, 23, iw, 0, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        infot = 12;
        Cstedc("V", 2, d, e, z, 2, w, 4, rw, 23, iw, 0, info);
        chkxer("Cstedc", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Cheevd
        //
        srnamt = "Cheevd";
        infot = 1;
        Cheevd("/", "U", 0, a, 1, x, w, 1, rw, 1, iw, 1, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 2;
        Cheevd("N", "/", 0, a, 1, x, w, 1, rw, 1, iw, 1, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 3;
        Cheevd("N", "U", -1, a, 1, x, w, 1, rw, 1, iw, 1, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 5;
        Cheevd("N", "U", 2, a, 1, x, w, 3, rw, 2, iw, 1, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 8;
        Cheevd("N", "U", 1, a, 1, x, w, 0, rw, 1, iw, 1, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 8;
        Cheevd("N", "U", 2, a, 2, x, w, 2, rw, 2, iw, 1, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 8;
        Cheevd("V", "U", 2, a, 2, x, w, 3, rw, 25, iw, 12, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 10;
        Cheevd("N", "U", 1, a, 1, x, w, 1, rw, 0, iw, 1, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 10;
        Cheevd("N", "U", 2, a, 2, x, w, 3, rw, 1, iw, 1, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 10;
        Cheevd("V", "U", 2, a, 2, x, w, 8, rw, 18, iw, 12, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 12;
        Cheevd("N", "U", 1, a, 1, x, w, 1, rw, 1, iw, 0, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        infot = 12;
        Cheevd("V", "U", 2, a, 2, x, w, 8, rw, 25, iw, 11, info);
        chkxer("Cheevd", infot, nout, lerr, ok);
        nt += 12;
        //
        //        Cheevd_2stage
        //
        srnamt = "Cheevd_2stage";
        infot = 1;
        Cheevd_2stage("/", "U", 0, a, 1, x, w, 1, rw, 1, iw, 1, info);
        chkxer("Cheevd_2stage", infot, nout, lerr, ok);
        infot = 1;
        Cheevd_2stage("V", "U", 0, a, 1, x, w, 1, rw, 1, iw, 1, info);
        chkxer("Cheevd_2stage", infot, nout, lerr, ok);
        infot = 2;
        Cheevd_2stage("N", "/", 0, a, 1, x, w, 1, rw, 1, iw, 1, info);
        chkxer("Cheevd_2stage", infot, nout, lerr, ok);
        infot = 3;
        Cheevd_2stage("N", "U", -1, a, 1, x, w, 1, rw, 1, iw, 1, info);
        chkxer("Cheevd_2stage", infot, nout, lerr, ok);
        infot = 5;
        Cheevd_2stage("N", "U", 2, a, 1, x, w, 3, rw, 2, iw, 1, info);
        chkxer("Cheevd_2stage", infot, nout, lerr, ok);
        infot = 8;
        Cheevd_2stage("N", "U", 1, a, 1, x, w, 0, rw, 1, iw, 1, info);
        chkxer("Cheevd_2stage", infot, nout, lerr, ok);
        infot = 8;
        Cheevd_2stage("N", "U", 2, a, 2, x, w, 2, rw, 2, iw, 1, info);
        chkxer("Cheevd_2stage", infot, nout, lerr, ok);
        //         INFOT = 8
        //         CALL Cheevd_2stage( 'V', 'U', 2, A, 2, X, W, 3,
        //     $                            RW, 25, IW, 12, INFO )
        //         CALL CHKXER( 'Cheevd_2stage', INFOT, NOUT, LERR, OK )
        infot = 10;
        Cheevd_2stage("N", "U", 1, a, 1, x, w, 1, rw, 0, iw, 1, info);
        chkxer("Cheevd_2stage", infot, nout, lerr, ok);
        infot = 10;
        Cheevd_2stage("N", "U", 2, a, 2, x, w, 25, rw, 1, iw, 1, info);
        chkxer("Cheevd_2stage", infot, nout, lerr, ok);
        //         INFOT = 10
        //         CALL Cheevd_2stage( 'V', 'U', 2, A, 2, X, W, 8,
        //     $                            RW, 18, IW, 12, INFO )
        //         CALL CHKXER( 'Cheevd_2stage', INFOT, NOUT, LERR, OK )
        infot = 12;
        Cheevd_2stage("N", "U", 1, a, 1, x, w, 1, rw, 1, iw, 0, info);
        chkxer("Cheevd_2stage", infot, nout, lerr, ok);
        infot = 12;
        //         CALL Cheevd_2stage( 'V', 'U', 2, A, 2, X, W, 8,
        //     $                            RW, 25, IW, 11, INFO )
        //         CALL CHKXER( 'Cheevd_2stage', INFOT, NOUT, LERR, OK )
        nt += 10;
        //
        //        Cheev
        //
        srnamt = "Cheev ";
        infot = 1;
        Cheev("/", "U", 0, a, 1, x, w, 1, rw, info);
        chkxer("Cheev ", infot, nout, lerr, ok);
        infot = 2;
        Cheev("N", "/", 0, a, 1, x, w, 1, rw, info);
        chkxer("Cheev ", infot, nout, lerr, ok);
        infot = 3;
        Cheev("N", "U", -1, a, 1, x, w, 1, rw, info);
        chkxer("Cheev ", infot, nout, lerr, ok);
        infot = 5;
        Cheev("N", "U", 2, a, 1, x, w, 3, rw, info);
        chkxer("Cheev ", infot, nout, lerr, ok);
        infot = 8;
        Cheev("N", "U", 2, a, 2, x, w, 2, rw, info);
        chkxer("Cheev ", infot, nout, lerr, ok);
        nt += 5;
        //
        //        Cheev_2stage
        //
        srnamt = "Cheev_2stage ";
        infot = 1;
        Cheev_2stage("/", "U", 0, a, 1, x, w, 1, rw, info);
        chkxer("Cheev_2stage ", infot, nout, lerr, ok);
        infot = 1;
        Cheev_2stage("V", "U", 0, a, 1, x, w, 1, rw, info);
        chkxer("Cheev_2stage ", infot, nout, lerr, ok);
        infot = 2;
        Cheev_2stage("N", "/", 0, a, 1, x, w, 1, rw, info);
        chkxer("Cheev_2stage ", infot, nout, lerr, ok);
        infot = 3;
        Cheev_2stage("N", "U", -1, a, 1, x, w, 1, rw, info);
        chkxer("Cheev_2stage ", infot, nout, lerr, ok);
        infot = 5;
        Cheev_2stage("N", "U", 2, a, 1, x, w, 3, rw, info);
        chkxer("Cheev_2stage ", infot, nout, lerr, ok);
        infot = 8;
        Cheev_2stage("N", "U", 2, a, 2, x, w, 2, rw, info);
        chkxer("Cheev_2stage ", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Cheevx
        //
        srnamt = "Cheevx";
        infot = 1;
        Cheevx("/", "A", "U", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx", infot, nout, lerr, ok);
        infot = 2;
        Cheevx("V", "/", "U", 0, a, 1, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx", infot, nout, lerr, ok);
        infot = 3;
        Cheevx("V", "A", "/", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        infot = 4;
        Cheevx("V", "A", "U", -1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx", infot, nout, lerr, ok);
        infot = 6;
        Cheevx("V", "A", "U", 2, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 2, w, 3, rw, iw, i3, info);
        chkxer("Cheevx", infot, nout, lerr, ok);
        infot = 8;
        Cheevx("V", "V", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx", infot, nout, lerr, ok);
        infot = 9;
        Cheevx("V", "I", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx", infot, nout, lerr, ok);
        infot = 10;
        Cheevx("V", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, x, z, 2, w, 3, rw, iw, i3, info);
        chkxer("Cheevx", infot, nout, lerr, ok);
        infot = 15;
        Cheevx("V", "A", "U", 2, a, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 3, rw, iw, i3, info);
        chkxer("Cheevx", infot, nout, lerr, ok);
        infot = 17;
        Cheevx("V", "A", "U", 2, a, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 2, w, 2, rw, iw, i1, info);
        chkxer("Cheevx", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Cheevx_2stage
        //
        srnamt = "Cheevx_2stage";
        infot = 1;
        Cheevx_2stage("/", "A", "U", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx_2stage", infot, nout, lerr, ok);
        infot = 1;
        Cheevx_2stage("V", "A", "U", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx_2stage", infot, nout, lerr, ok);
        infot = 2;
        Cheevx_2stage("N", "/", "U", 0, a, 1, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx_2stage", infot, nout, lerr, ok);
        infot = 3;
        Cheevx_2stage("N", "A", "/", 0, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        infot = 4;
        Cheevx_2stage("N", "A", "U", -1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx_2stage", infot, nout, lerr, ok);
        infot = 6;
        Cheevx_2stage("N", "A", "U", 2, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 2, w, 3, rw, iw, i3, info);
        chkxer("Cheevx_2stage", infot, nout, lerr, ok);
        infot = 8;
        Cheevx_2stage("N", "V", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx_2stage", infot, nout, lerr, ok);
        infot = 9;
        Cheevx_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 1, rw, iw, i3, info);
        chkxer("Cheevx_2stage", infot, nout, lerr, ok);
        infot = 10;
        Cheevx_2stage("N", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, x, z, 2, w, 3, rw, iw, i3, info);
        chkxer("Cheevx_2stage", infot, nout, lerr, ok);
        infot = 15;
        Cheevx_2stage("N", "A", "U", 2, a, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 0, w, 3, rw, iw, i3, info);
        chkxer("Cheevx_2stage", infot, nout, lerr, ok);
        infot = 17;
        Cheevx_2stage("N", "A", "U", 2, a, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 2, w, 0, rw, iw, i1, info);
        chkxer("Cheevx_2stage", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Cheevr
        //
        srnamt = "Cheevr";
        n = 1;
        infot = 1;
        Cheevr("/", "A", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 2;
        Cheevr("V", "/", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 3;
        Cheevr("V", "A", "/", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 4;
        Cheevr("V", "A", "U", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 6;
        Cheevr("V", "A", "U", 2, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 8;
        Cheevr("V", "V", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 9;
        Cheevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 0, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 10;
        //
        Cheevr("V", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 15;
        Cheevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 0, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 18;
        Cheevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n - 1, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 20;
        Cheevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n - 1, iw[(2 * n - 1) - 1], 10 * n, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        infot = 22;
        Cheevr("V", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw, 10 * n - 1, info);
        chkxer("Cheevr", infot, nout, lerr, ok);
        nt += 12;
        //
        //        Cheevr_2stage
        //
        srnamt = "Cheevr_2stage";
        n = 1;
        infot = 1;
        Cheevr_2stage("/", "A", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 1;
        Cheevr_2stage("V", "A", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 2;
        Cheevr_2stage("N", "/", "U", 0, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 3;
        Cheevr_2stage("N", "A", "/", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 4;
        Cheevr_2stage("N", "A", "U", -1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 6;
        Cheevr_2stage("N", "A", "U", 2, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 8;
        Cheevr_2stage("N", "V", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 9;
        Cheevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 0, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 10;
        Cheevr_2stage("N", "I", "U", 2, a, 2, 0.0, 0.0, 2, 1, 0.0, m, r, z, 1, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 15;
        Cheevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 0, iw, q, 2 * n, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 18;
        Cheevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 2 * n - 1, rw, 24 * n, iw[(2 * n + 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 20;
        Cheevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, rw, 24 * n - 1, iw[(2 * n - 1) - 1], 10 * n, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        infot = 22;
        Cheevr_2stage("N", "I", "U", 1, a, 1, 0.0, 0.0, 1, 1, 0.0, m, r, z, 1, iw, q, 26 * n, rw, 24 * n, iw, 10 * n - 1, info);
        chkxer("Cheevr_2stage", infot, nout, lerr, ok);
        nt += 13;
        //
        //        Chpevd
        //
        srnamt = "Chpevd";
        infot = 1;
        Chpevd("/", "U", 0, a, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 2;
        Chpevd("N", "/", 0, a, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 3;
        Chpevd("N", "U", -1, a, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 7;
        Chpevd("V", "U", 2, a, x, z, 1, w, 4, rw, 25, iw, 12, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 9;
        Chpevd("N", "U", 1, a, x, z, 1, w, 0, rw, 1, iw, 1, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 9;
        Chpevd("N", "U", 2, a, x, z, 2, w, 1, rw, 2, iw, 1, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 9;
        Chpevd("V", "U", 2, a, x, z, 2, w, 2, rw, 25, iw, 12, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 11;
        Chpevd("N", "U", 1, a, x, z, 1, w, 1, rw, 0, iw, 1, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 11;
        Chpevd("N", "U", 2, a, x, z, 2, w, 2, rw, 1, iw, 1, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 11;
        Chpevd("V", "U", 2, a, x, z, 2, w, 4, rw, 18, iw, 12, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 13;
        Chpevd("N", "U", 1, a, x, z, 1, w, 1, rw, 1, iw, 0, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 13;
        Chpevd("N", "U", 2, a, x, z, 2, w, 2, rw, 2, iw, 0, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        infot = 13;
        Chpevd("V", "U", 2, a, x, z, 2, w, 4, rw, 25, iw, 2, info);
        chkxer("Chpevd", infot, nout, lerr, ok);
        nt += 13;
        //
        //        Chpev
        //
        srnamt = "Chpev ";
        infot = 1;
        Chpev("/", "U", 0, a, x, z, 1, w, rw, info);
        chkxer("Chpev ", infot, nout, lerr, ok);
        infot = 2;
        Chpev("N", "/", 0, a, x, z, 1, w, rw, info);
        chkxer("Chpev ", infot, nout, lerr, ok);
        infot = 3;
        Chpev("N", "U", -1, a, x, z, 1, w, rw, info);
        chkxer("Chpev ", infot, nout, lerr, ok);
        infot = 7;
        Chpev("V", "U", 2, a, x, z, 1, w, rw, info);
        chkxer("Chpev ", infot, nout, lerr, ok);
        nt += 4;
        //
        //        Chpevx
        //
        srnamt = "Chpevx";
        infot = 1;
        Chpevx("/", "A", "U", 0, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chpevx", infot, nout, lerr, ok);
        infot = 2;
        Chpevx("V", "/", "U", 0, a, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chpevx", infot, nout, lerr, ok);
        infot = 3;
        Chpevx("V", "A", "/", 0, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chpevx", infot, nout, lerr, ok);
        infot = 4;
        Chpevx("V", "A", "U", -1, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chpevx", infot, nout, lerr, ok);
        infot = 7;
        Chpevx("V", "V", "U", 1, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chpevx", infot, nout, lerr, ok);
        infot = 8;
        Chpevx("V", "I", "U", 1, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chpevx", infot, nout, lerr, ok);
        infot = 9;
        Chpevx("V", "I", "U", 2, a, 0.0, 0.0, 2, 1, 0.0, m, x, z, 2, w, rw, iw, i3, info);
        chkxer("Chpevx", infot, nout, lerr, ok);
        infot = 14;
        Chpevx("V", "A", "U", 2, a, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chpevx", infot, nout, lerr, ok);
        nt += 8;
        //
        //     Test error exits for the HB path.
        //
    } else if (Mlsamen(2, c2, "HB")) {
        //
        //        Chbtrd
        //
        srnamt = "Chbtrd";
        infot = 1;
        Chbtrd("/", "U", 0, 0, a, 1, d, e, z, 1, w, info);
        chkxer("Chbtrd", infot, nout, lerr, ok);
        infot = 2;
        Chbtrd("N", "/", 0, 0, a, 1, d, e, z, 1, w, info);
        chkxer("Chbtrd", infot, nout, lerr, ok);
        infot = 3;
        Chbtrd("N", "U", -1, 0, a, 1, d, e, z, 1, w, info);
        chkxer("Chbtrd", infot, nout, lerr, ok);
        infot = 4;
        Chbtrd("N", "U", 0, -1, a, 1, d, e, z, 1, w, info);
        chkxer("Chbtrd", infot, nout, lerr, ok);
        infot = 6;
        Chbtrd("N", "U", 1, 1, a, 1, d, e, z, 1, w, info);
        chkxer("Chbtrd", infot, nout, lerr, ok);
        infot = 10;
        Chbtrd("V", "U", 2, 0, a, 1, d, e, z, 1, w, info);
        chkxer("Chbtrd", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Chetrd_HB2ST
        //
        srnamt = "Chetrd_HB2ST";
        infot = 1;
        Chetrd_hb2st("/", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 2;
        Chetrd_hb2st("N", "/", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 2;
        Chetrd_hb2st("N", "H", "U", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 3;
        Chetrd_hb2st("N", "N", "/", 0, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 4;
        Chetrd_hb2st("N", "N", "U", -1, 0, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 5;
        Chetrd_hb2st("N", "N", "U", 0, -1, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 7;
        Chetrd_hb2st("N", "N", "U", 0, 1, a, 1, d, e, c, 1, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 11;
        Chetrd_hb2st("N", "N", "U", 0, 0, a, 1, d, e, c, 0, w, 1, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        infot = 13;
        Chetrd_hb2st("N", "N", "U", 0, 0, a, 1, d, e, c, 1, w, 0, info);
        chkxer("Chetrd_HB2ST", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Chbevd
        //
        srnamt = "Chbevd";
        infot = 1;
        Chbevd("/", "U", 0, 0, a, 1, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 2;
        Chbevd("N", "/", 0, 0, a, 1, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 3;
        Chbevd("N", "U", -1, 0, a, 1, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 4;
        Chbevd("N", "U", 0, -1, a, 1, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 6;
        Chbevd("N", "U", 2, 1, a, 1, x, z, 1, w, 2, rw, 2, iw, 1, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 9;
        Chbevd("V", "U", 2, 1, a, 2, x, z, 1, w, 8, rw, 25, iw, 12, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 11;
        Chbevd("N", "U", 1, 0, a, 1, x, z, 1, w, 0, rw, 1, iw, 1, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 11;
        Chbevd("N", "U", 2, 1, a, 2, x, z, 2, w, 1, rw, 2, iw, 1, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 11;
        Chbevd("V", "U", 2, 1, a, 2, x, z, 2, w, 2, rw, 25, iw, 12, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 13;
        Chbevd("N", "U", 1, 0, a, 1, x, z, 1, w, 1, rw, 0, iw, 1, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 13;
        Chbevd("N", "U", 2, 1, a, 2, x, z, 2, w, 2, rw, 1, iw, 1, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 13;
        Chbevd("V", "U", 2, 1, a, 2, x, z, 2, w, 8, rw, 2, iw, 12, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 15;
        Chbevd("N", "U", 1, 0, a, 1, x, z, 1, w, 1, rw, 1, iw, 0, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 15;
        Chbevd("N", "U", 2, 1, a, 2, x, z, 2, w, 2, rw, 2, iw, 0, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        infot = 15;
        Chbevd("V", "U", 2, 1, a, 2, x, z, 2, w, 8, rw, 25, iw, 2, info);
        chkxer("Chbevd", infot, nout, lerr, ok);
        nt += 15;
        //
        //        Chbevd_2stage
        //
        srnamt = "Chbevd_2stage";
        infot = 1;
        Chbevd_2stage("/", "U", 0, 0, a, 1, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        infot = 1;
        Chbevd_2stage("V", "U", 0, 0, a, 1, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        infot = 2;
        Chbevd_2stage("N", "/", 0, 0, a, 1, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        infot = 3;
        Chbevd_2stage("N", "U", -1, 0, a, 1, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        infot = 4;
        Chbevd_2stage("N", "U", 0, -1, a, 1, x, z, 1, w, 1, rw, 1, iw, 1, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        infot = 6;
        Chbevd_2stage("N", "U", 2, 1, a, 1, x, z, 1, w, 2, rw, 2, iw, 1, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        infot = 9;
        Chbevd_2stage("N", "U", 2, 1, a, 2, x, z, 0, w, 8, rw, 25, iw, 12, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        infot = 11;
        Chbevd_2stage("N", "U", 1, 0, a, 1, x, z, 1, w, 0, rw, 1, iw, 1, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        infot = 11;
        Chbevd_2stage("N", "U", 2, 1, a, 2, x, z, 2, w, 1, rw, 2, iw, 1, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        //         INFOT = 11
        //         CALL Chbevd_2stage( 'V', 'U', 2, 1, A, 2, X, Z, 2,
        //     $                         W, 2, RW, 25, IW, 12, INFO )
        //         CALL CHKXER( 'Chbevd_2stage', INFOT, NOUT, LERR, OK )
        infot = 13;
        Chbevd_2stage("N", "U", 1, 0, a, 1, x, z, 1, w, 1, rw, 0, iw, 1, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        infot = 13;
        Chbevd_2stage("N", "U", 2, 1, a, 2, x, z, 2, w, 25, rw, 1, iw, 1, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        //         INFOT = 13
        //         CALL Chbevd_2stage( 'V', 'U', 2, 1, A, 2, X, Z, 2,
        //     $                          W, 25, RW, 2, IW, 12, INFO )
        //         CALL CHKXER( 'Chbevd_2stage', INFOT, NOUT, LERR, OK )
        infot = 15;
        Chbevd_2stage("N", "U", 1, 0, a, 1, x, z, 1, w, 1, rw, 1, iw, 0, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        infot = 15;
        Chbevd_2stage("N", "U", 2, 1, a, 2, x, z, 2, w, 25, rw, 2, iw, 0, info);
        chkxer("Chbevd_2stage", infot, nout, lerr, ok);
        //         INFOT = 15
        //         CALL Chbevd_2stage( 'V', 'U', 2, 1, A, 2, X, Z, 2,
        //     $                          W, 25, RW, 25, IW, 2, INFO )
        //         CALL CHKXER( 'Chbevd_2stage', INFOT, NOUT, LERR, OK )
        nt += 13;
        //
        //        Chbev
        //
        srnamt = "Chbev ";
        infot = 1;
        Chbev("/", "U", 0, 0, a, 1, x, z, 1, w, rw, info);
        chkxer("Chbev ", infot, nout, lerr, ok);
        infot = 2;
        Chbev("N", "/", 0, 0, a, 1, x, z, 1, w, rw, info);
        chkxer("Chbev ", infot, nout, lerr, ok);
        infot = 3;
        Chbev("N", "U", -1, 0, a, 1, x, z, 1, w, rw, info);
        chkxer("Chbev ", infot, nout, lerr, ok);
        infot = 4;
        Chbev("N", "U", 0, -1, a, 1, x, z, 1, w, rw, info);
        chkxer("Chbev ", infot, nout, lerr, ok);
        infot = 6;
        Chbev("N", "U", 2, 1, a, 1, x, z, 1, w, rw, info);
        chkxer("Chbev ", infot, nout, lerr, ok);
        infot = 9;
        Chbev("V", "U", 2, 0, a, 1, x, z, 1, w, rw, info);
        chkxer("Chbev ", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Chbev_2stage
        //
        srnamt = "Chbev_2stage ";
        infot = 1;
        Chbev_2stage("/", "U", 0, 0, a, 1, x, z, 1, w, 0, rw, info);
        chkxer("Chbev_2stage ", infot, nout, lerr, ok);
        infot = 1;
        Chbev_2stage("V", "U", 0, 0, a, 1, x, z, 1, w, 0, rw, info);
        chkxer("Chbev_2stage ", infot, nout, lerr, ok);
        infot = 2;
        Chbev_2stage("N", "/", 0, 0, a, 1, x, z, 1, w, 0, rw, info);
        chkxer("Chbev_2stage ", infot, nout, lerr, ok);
        infot = 3;
        Chbev_2stage("N", "U", -1, 0, a, 1, x, z, 1, w, 0, rw, info);
        chkxer("Chbev_2stage ", infot, nout, lerr, ok);
        infot = 4;
        Chbev_2stage("N", "U", 0, -1, a, 1, x, z, 1, w, 0, rw, info);
        chkxer("Chbev_2stage ", infot, nout, lerr, ok);
        infot = 6;
        Chbev_2stage("N", "U", 2, 1, a, 1, x, z, 1, w, 0, rw, info);
        chkxer("Chbev_2stage ", infot, nout, lerr, ok);
        infot = 9;
        Chbev_2stage("N", "U", 2, 0, a, 1, x, z, 0, w, 0, rw, info);
        chkxer("Chbev_2stage ", infot, nout, lerr, ok);
        infot = 11;
        Chbev_2stage("N", "U", 2, 0, a, 1, x, z, 1, w, 0, rw, info);
        chkxer("Chbev_2stage ", infot, nout, lerr, ok);
        nt += 8;
        //
        //        Chbevx
        //
        srnamt = "Chbevx";
        infot = 1;
        Chbevx("/", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chbevx", infot, nout, lerr, ok);
        infot = 2;
        Chbevx("V", "/", "U", 0, 0, a, 1, q, 1, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chbevx", infot, nout, lerr, ok);
        infot = 3;
        Chbevx("V", "A", "/", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        infot = 4;
        Chbevx("V", "A", "U", -1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chbevx", infot, nout, lerr, ok);
        infot = 5;
        Chbevx("V", "A", "U", 0, -1, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chbevx", infot, nout, lerr, ok);
        infot = 7;
        Chbevx("V", "A", "U", 2, 1, a, 1, q, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 2, w, rw, iw, i3, info);
        chkxer("Chbevx", infot, nout, lerr, ok);
        infot = 9;
        Chbevx("V", "A", "U", 2, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 2, w, rw, iw, i3, info);
        chkxer("Chbevx", infot, nout, lerr, ok);
        infot = 11;
        Chbevx("V", "V", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chbevx", infot, nout, lerr, ok);
        infot = 12;
        Chbevx("V", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chbevx", infot, nout, lerr, ok);
        infot = 13;
        Chbevx("V", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chbevx", infot, nout, lerr, ok);
        infot = 18;
        Chbevx("V", "A", "U", 2, 0, a, 1, q, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, rw, iw, i3, info);
        chkxer("Chbevx", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Chbevx_2stage
        //
        srnamt = "Chbevx_2stage";
        infot = 1;
        Chbevx_2stage("/", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, rw, iw, i3, info);
        infot = 1;
        Chbevx_2stage("V", "A", "U", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, rw, iw, i3, info);
        chkxer("Chbevx_2stage", infot, nout, lerr, ok);
        infot = 2;
        Chbevx_2stage("N", "/", "U", 0, 0, a, 1, q, 1, 0.0, 1.0, 1, 0, 0.0, m, x, z, 1, w, 0, rw, iw, i3, info);
        chkxer("Chbevx_2stage", infot, nout, lerr, ok);
        infot = 3;
        Chbevx_2stage("N", "A", "/", 0, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, rw, iw, i3, info);
        infot = 4;
        Chbevx_2stage("N", "A", "U", -1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, rw, iw, i3, info);
        chkxer("Chbevx_2stage", infot, nout, lerr, ok);
        infot = 5;
        Chbevx_2stage("N", "A", "U", 0, -1, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, rw, iw, i3, info);
        chkxer("Chbevx_2stage", infot, nout, lerr, ok);
        infot = 7;
        Chbevx_2stage("N", "A", "U", 2, 1, a, 1, q, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 2, w, 0, rw, iw, i3, info);
        chkxer("Chbevx_2stage", infot, nout, lerr, ok);
        //         INFOT = 9
        //         CALL Chbevx_2stage( 'V', 'A', 'U', 2, 0, A, 1, Q, 1,
        //     $                       0.0D0, 0.0D0, 0, 0, 0.0D0,
        //     $                       M, X, Z, 2, W, 0, RW, IW, I3, INFO )
        //         CALL CHKXER( 'Chbevx_2stage', INFOT, NOUT, LERR, OK )
        infot = 11;
        Chbevx_2stage("N", "V", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, rw, iw, i3, info);
        chkxer("Chbevx_2stage", infot, nout, lerr, ok);
        infot = 12;
        Chbevx_2stage("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, rw, iw, i3, info);
        chkxer("Chbevx_2stage", infot, nout, lerr, ok);
        infot = 13;
        Chbevx_2stage("N", "I", "U", 1, 0, a, 1, q, 1, 0.0, 0.0, 1, 2, 0.0, m, x, z, 1, w, 0, rw, iw, i3, info);
        chkxer("Chbevx_2stage", infot, nout, lerr, ok);
        infot = 18;
        Chbevx_2stage("N", "A", "U", 2, 0, a, 1, q, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 0, w, 0, rw, iw, i3, info);
        chkxer("Chbevx_2stage", infot, nout, lerr, ok);
        infot = 20;
        Chbevx_2stage("N", "A", "U", 2, 0, a, 1, q, 2, 0.0, 0.0, 0, 0, 0.0, m, x, z, 1, w, 0, rw, iw, i3, info);
        chkxer("Chbevx_2stage", infot, nout, lerr, ok);
        nt += 12;
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
    //     End of Cerrst
    //
}
