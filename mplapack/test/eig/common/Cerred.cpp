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

void Cerred(const char *path, INTEGER const nunit) {
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
    COMPLEX a[nmax * nmax];
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
    COMPLEX x[nmax];
    COMPLEX vl[nmax * nmax];
    COMPLEX vr[nmax * nmax];
    COMPLEX w[10 * nmax];
    const INTEGER lw = 5 * nmax;
    REAL rw[lw];
    INTEGER info = 0;
    INTEGER sdim = 0;
    bool b[nmax];
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL s[nmax];
    REAL abnrm = 0.0;
    REAL r1[nmax];
    REAL r2[nmax];
    COMPLEX u[nmax * nmax];
    COMPLEX vt[nmax * nmax];
    INTEGER iw[4 * nmax];
    INTEGER ns = 0;
    if (Mlsamen(2, c2, "EV")) {
        //
        //        Test Cgeev
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Cgeev ", 16);
        infot = 1;
        Cgeev("X", "N", 0, a, 1, x, vl, 1, vr, 1, w, 1, rw, info);
        chkxer("Cgeev ", infot, nout, lerr, ok);
        infot = 2;
        Cgeev("N", "X", 0, a, 1, x, vl, 1, vr, 1, w, 1, rw, info);
        chkxer("Cgeev ", infot, nout, lerr, ok);
        infot = 3;
        Cgeev("N", "N", -1, a, 1, x, vl, 1, vr, 1, w, 1, rw, info);
        chkxer("Cgeev ", infot, nout, lerr, ok);
        infot = 5;
        Cgeev("N", "N", 2, a, 1, x, vl, 1, vr, 1, w, 4, rw, info);
        chkxer("Cgeev ", infot, nout, lerr, ok);
        infot = 8;
        Cgeev("V", "N", 2, a, 2, x, vl, 1, vr, 1, w, 4, rw, info);
        chkxer("Cgeev ", infot, nout, lerr, ok);
        infot = 10;
        Cgeev("N", "V", 2, a, 2, x, vl, 1, vr, 1, w, 4, rw, info);
        chkxer("Cgeev ", infot, nout, lerr, ok);
        infot = 12;
        Cgeev("V", "V", 1, a, 1, x, vl, 1, vr, 1, w, 1, rw, info);
        chkxer("Cgeev ", infot, nout, lerr, ok);
        nt += 7;
        //
    } else if (Mlsamen(2, c2, "ES")) {
        //
        //        Test Cgees
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Cgees ", 16);
        infot = 1;
        Cgees("X", "N", Cslect, 0, a, 1, sdim, x, vl, 1, w, 1, rw, b, info);
        chkxer("Cgees ", infot, nout, lerr, ok);
        infot = 2;
        Cgees("N", "X", Cslect, 0, a, 1, sdim, x, vl, 1, w, 1, rw, b, info);
        chkxer("Cgees ", infot, nout, lerr, ok);
        infot = 4;
        Cgees("N", "S", Cslect, -1, a, 1, sdim, x, vl, 1, w, 1, rw, b, info);
        chkxer("Cgees ", infot, nout, lerr, ok);
        infot = 6;
        Cgees("N", "S", Cslect, 2, a, 1, sdim, x, vl, 1, w, 4, rw, b, info);
        chkxer("Cgees ", infot, nout, lerr, ok);
        infot = 10;
        Cgees("V", "S", Cslect, 2, a, 2, sdim, x, vl, 1, w, 4, rw, b, info);
        chkxer("Cgees ", infot, nout, lerr, ok);
        infot = 12;
        Cgees("N", "S", Cslect, 1, a, 1, sdim, x, vl, 1, w, 1, rw, b, info);
        chkxer("Cgees ", infot, nout, lerr, ok);
        nt += 6;
        //
    } else if (Mlsamen(2, c2, "VX")) {
        //
        //        Test Cgeevx
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Cgeevx", 16);
        infot = 1;
        Cgeevx("X", "N", "N", "N", 0, a, 1, x, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, rw, info);
        chkxer("Cgeevx", infot, nout, lerr, ok);
        infot = 2;
        Cgeevx("N", "X", "N", "N", 0, a, 1, x, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, rw, info);
        chkxer("Cgeevx", infot, nout, lerr, ok);
        infot = 3;
        Cgeevx("N", "N", "X", "N", 0, a, 1, x, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, rw, info);
        chkxer("Cgeevx", infot, nout, lerr, ok);
        infot = 4;
        Cgeevx("N", "N", "N", "X", 0, a, 1, x, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, rw, info);
        chkxer("Cgeevx", infot, nout, lerr, ok);
        infot = 5;
        Cgeevx("N", "N", "N", "N", -1, a, 1, x, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, rw, info);
        chkxer("Cgeevx", infot, nout, lerr, ok);
        infot = 7;
        Cgeevx("N", "N", "N", "N", 2, a, 1, x, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 4, rw, info);
        chkxer("Cgeevx", infot, nout, lerr, ok);
        infot = 10;
        Cgeevx("N", "V", "N", "N", 2, a, 2, x, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 4, rw, info);
        chkxer("Cgeevx", infot, nout, lerr, ok);
        infot = 12;
        Cgeevx("N", "N", "V", "N", 2, a, 2, x, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 4, rw, info);
        chkxer("Cgeevx", infot, nout, lerr, ok);
        infot = 20;
        Cgeevx("N", "N", "N", "N", 1, a, 1, x, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 1, rw, info);
        chkxer("Cgeevx", infot, nout, lerr, ok);
        infot = 20;
        Cgeevx("N", "N", "V", "V", 1, a, 1, x, vl, 1, vr, 1, ilo, ihi, s, abnrm, r1, r2, w, 2, rw, info);
        chkxer("Cgeevx", infot, nout, lerr, ok);
        nt += 10;
        //
    } else if (Mlsamen(2, c2, "SX")) {
        //
        //        Test Cgeesx
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Cgeesx", 16);
        infot = 1;
        Cgeesx("X", "N", Cslect, "N", 0, a, 1, sdim, x, vl, 1, r1[1 - 1], r2[1 - 1], w, 1, rw, b, info);
        chkxer("Cgeesx", infot, nout, lerr, ok);
        infot = 2;
        Cgeesx("N", "X", Cslect, "N", 0, a, 1, sdim, x, vl, 1, r1[1 - 1], r2[1 - 1], w, 1, rw, b, info);
        chkxer("Cgeesx", infot, nout, lerr, ok);
        infot = 4;
        Cgeesx("N", "N", Cslect, "X", 0, a, 1, sdim, x, vl, 1, r1[1 - 1], r2[1 - 1], w, 1, rw, b, info);
        chkxer("Cgeesx", infot, nout, lerr, ok);
        infot = 5;
        Cgeesx("N", "N", Cslect, "N", -1, a, 1, sdim, x, vl, 1, r1[1 - 1], r2[1 - 1], w, 1, rw, b, info);
        chkxer("Cgeesx", infot, nout, lerr, ok);
        infot = 7;
        Cgeesx("N", "N", Cslect, "N", 2, a, 1, sdim, x, vl, 1, r1[1 - 1], r2[1 - 1], w, 4, rw, b, info);
        chkxer("Cgeesx", infot, nout, lerr, ok);
        infot = 11;
        Cgeesx("V", "N", Cslect, "N", 2, a, 2, sdim, x, vl, 1, r1[1 - 1], r2[1 - 1], w, 4, rw, b, info);
        chkxer("Cgeesx", infot, nout, lerr, ok);
        infot = 15;
        Cgeesx("N", "N", Cslect, "N", 1, a, 1, sdim, x, vl, 1, r1[1 - 1], r2[1 - 1], w, 1, rw, b, info);
        chkxer("Cgeesx", infot, nout, lerr, ok);
        nt += 7;
        //
    } else if (Mlsamen(2, c2, "BD")) {
        //
        //        Test Cgesvd
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Cgesvd", 16);
        infot = 1;
        Cgesvd("X", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, info);
        chkxer("Cgesvd", infot, nout, lerr, ok);
        infot = 2;
        Cgesvd("N", "X", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, info);
        chkxer("Cgesvd", infot, nout, lerr, ok);
        infot = 2;
        Cgesvd("O", "O", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, info);
        chkxer("Cgesvd", infot, nout, lerr, ok);
        infot = 3;
        Cgesvd("N", "N", -1, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, info);
        chkxer("Cgesvd", infot, nout, lerr, ok);
        infot = 4;
        Cgesvd("N", "N", 0, -1, a, 1, s, u, 1, vt, 1, w, 1, rw, info);
        chkxer("Cgesvd", infot, nout, lerr, ok);
        infot = 6;
        Cgesvd("N", "N", 2, 1, a, 1, s, u, 1, vt, 1, w, 5, rw, info);
        chkxer("Cgesvd", infot, nout, lerr, ok);
        infot = 9;
        Cgesvd("A", "N", 2, 1, a, 2, s, u, 1, vt, 1, w, 5, rw, info);
        chkxer("Cgesvd", infot, nout, lerr, ok);
        infot = 11;
        Cgesvd("N", "A", 1, 2, a, 1, s, u, 1, vt, 1, w, 5, rw, info);
        chkxer("Cgesvd", infot, nout, lerr, ok);
        nt += 8;
        if (ok) {
            write(nout, format_9999), srnamt, nt;
        } else {
            write(nout, format_9998);
        }
        //
        //        Test Cgesdd
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Cgesdd", 16);
        infot = 1;
        Cgesdd("X", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesdd", infot, nout, lerr, ok);
        infot = 2;
        Cgesdd("N", -1, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesdd", infot, nout, lerr, ok);
        infot = 3;
        Cgesdd("N", 0, -1, a, 1, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesdd", infot, nout, lerr, ok);
        infot = 5;
        Cgesdd("N", 2, 1, a, 1, s, u, 1, vt, 1, w, 5, rw, iw, info);
        chkxer("Cgesdd", infot, nout, lerr, ok);
        infot = 8;
        Cgesdd("A", 2, 1, a, 2, s, u, 1, vt, 1, w, 5, rw, iw, info);
        chkxer("Cgesdd", infot, nout, lerr, ok);
        infot = 10;
        Cgesdd("A", 1, 2, a, 1, s, u, 1, vt, 1, w, 5, rw, iw, info);
        chkxer("Cgesdd", infot, nout, lerr, ok);
        nt = nt - 2;
        if (ok) {
            write(nout, format_9999), srnamt, nt;
        } else {
            write(nout, format_9998);
        }
        //
        //        Test Cgejsv
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Cgejsv", 16);
        infot = 1;
        Cgejsv("X", "U", "V", "R", "N", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        infot = 2;
        Cgejsv("G", "X", "V", "R", "N", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        infot = 3;
        Cgejsv("G", "U", "X", "R", "N", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        infot = 4;
        Cgejsv("G", "U", "V", "X", "N", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        infot = 5;
        Cgejsv("G", "U", "V", "R", "X", "N", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        infot = 6;
        Cgejsv("G", "U", "V", "R", "N", "X", 0, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        infot = 7;
        Cgejsv("G", "U", "V", "R", "N", "N", -1, 0, a, 1, s, u, 1, vt, 1, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        infot = 8;
        Cgejsv("G", "U", "V", "R", "N", "N", 0, -1, a, 1, s, u, 1, vt, 1, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        infot = 10;
        Cgejsv("G", "U", "V", "R", "N", "N", 2, 1, a, 1, s, u, 1, vt, 1, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        infot = 13;
        Cgejsv("G", "U", "V", "R", "N", "N", 2, 2, a, 2, s, u, 1, vt, 2, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        infot = 15;
        Cgejsv("G", "U", "V", "R", "N", "N", 2, 2, a, 2, s, u, 2, vt, 1, w, 1, rw, 1, iw, info);
        chkxer("Cgejsv", infot, nout, lerr, ok);
        nt = 11;
        if (ok) {
            write(nout, format_9999), srnamt, nt;
        } else {
            write(nout, format_9998);
        }
        //
        //        Test Cgesvdx
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Cgesvdx", 16);
        infot = 1;
        Cgesvdx("X", "N", "A", 0, 0, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 2;
        Cgesvdx("N", "X", "A", 0, 0, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 3;
        Cgesvdx("N", "N", "X", 0, 0, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 4;
        Cgesvdx("N", "N", "A", -1, 0, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 5;
        Cgesvdx("N", "N", "A", 0, -1, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 7;
        Cgesvdx("N", "N", "A", 2, 1, a, 1, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 8;
        Cgesvdx("N", "N", "V", 2, 1, a, 2, -one, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 9;
        Cgesvdx("N", "N", "V", 2, 1, a, 2, one, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 10;
        Cgesvdx("N", "N", "I", 2, 2, a, 2, zero, zero, 0, 1, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 11;
        Cgesvdx("V", "N", "I", 2, 2, a, 2, zero, zero, 1, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 15;
        Cgesvdx("V", "N", "A", 2, 2, a, 2, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        infot = 17;
        Cgesvdx("N", "V", "A", 2, 2, a, 2, zero, zero, 0, 0, ns, s, u, 1, vt, 1, w, 1, rw, iw, info);
        chkxer("Cgesvdx", infot, nout, lerr, ok);
        nt = 12;
        if (ok) {
            write(nout, format_9999), srnamt, nt;
        } else {
            write(nout, format_9998);
        }
        //
        //        Test Cgesvdq
        //
        memset(srnamt, 0, sizeof(srnamt));
        strncpy(srnamt, "Cgesvdq", 16);
        infot = 1;
        Cgesvdq("X", "P", "T", "A", "A", 0, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
        infot = 2;
        Cgesvdq("A", "X", "T", "A", "A", 0, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
        infot = 3;
        Cgesvdq("A", "P", "X", "A", "A", 0, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
        infot = 4;
        Cgesvdq("A", "P", "T", "X", "A", 0, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
        infot = 5;
        Cgesvdq("A", "P", "T", "A", "X", 0, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
        infot = 6;
        Cgesvdq("A", "P", "T", "A", "A", -1, 0, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
        infot = 7;
        Cgesvdq("A", "P", "T", "A", "A", 0, 1, a, 1, s, u, 0, vt, 0, ns, iw, 1, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
        infot = 9;
        Cgesvdq("A", "P", "T", "A", "A", 1, 1, a, 0, s, u, 0, vt, 0, ns, iw, 1, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
        infot = 12;
        Cgesvdq("A", "P", "T", "A", "A", 1, 1, a, 1, s, u, -1, vt, 0, ns, iw, 1, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
        infot = 14;
        Cgesvdq("A", "P", "T", "A", "A", 1, 1, a, 1, s, u, 1, vt, -1, ns, iw, 1, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
        infot = 17;
        Cgesvdq("A", "P", "T", "A", "A", 1, 1, a, 1, s, u, 1, vt, 1, ns, iw, -5, w, 1, rw, 1, info);
        chkxer("Cgesvdq", infot, nout, lerr, ok);
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
    //     End of Cerred
    //
}
