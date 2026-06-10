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

// Derived from LAPACK routine DERRHS.
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

void Rerrhs(fem::str_cref path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    static const char *format_9999 = "(1x,a3,' routines passed the tests of the error exits',' (',i3,"
                                     "' tests done)')";
    static const char *format_9998 = "(' *** ',a3,' routines failed the tests of the error ','exits ***')";
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
    REAL wi[nmax];
    bool sel[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * nmax] = 1.0 / castREAL(i + j);
        }
        wi[j - 1] = castREAL(j);
        sel[j - 1] = true;
    }
    ok = true;
    INTEGER nt = 0;
    //
    // Test error exits of the nonsymmetric eigenvalue routines.
    //
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL s[nmax];
    INTEGER info = 0;
    REAL tau[nmax];
    const INTEGER lw = (nmax + 2) * (nmax + 2) + nmax;
    REAL w[lw];
    REAL c[nmax * nmax];
    REAL wr[nmax];
    REAL vl[nmax * nmax];
    REAL vr[nmax * nmax];
    INTEGER m = 0;
    INTEGER ifaill[nmax];
    INTEGER ifailr[nmax];
    if (Mlsamen(2, c2.elems, "HS")) {
        //
        // Rgebal
        //
        srnamt = "Rgebal";
        infot = 1;
        Rgebal("/", 0, a, 1, ilo, ihi, s, info);
        Chkxer("Rgebal", infot, nout, lerr, ok);
        infot = 2;
        Rgebal("N", -1, a, 1, ilo, ihi, s, info);
        Chkxer("Rgebal", infot, nout, lerr, ok);
        infot = 4;
        Rgebal("N", 2, a, 1, ilo, ihi, s, info);
        Chkxer("Rgebal", infot, nout, lerr, ok);
        nt += 3;
        //
        // Rgebak
        //
        srnamt = "Rgebak";
        infot = 1;
        Rgebak("/", "R", 0, 1, 0, s, 0, a, 1, info);
        Chkxer("Rgebak", infot, nout, lerr, ok);
        infot = 2;
        Rgebak("N", "/", 0, 1, 0, s, 0, a, 1, info);
        Chkxer("Rgebak", infot, nout, lerr, ok);
        infot = 3;
        Rgebak("N", "R", -1, 1, 0, s, 0, a, 1, info);
        Chkxer("Rgebak", infot, nout, lerr, ok);
        infot = 4;
        Rgebak("N", "R", 0, 0, 0, s, 0, a, 1, info);
        Chkxer("Rgebak", infot, nout, lerr, ok);
        infot = 4;
        Rgebak("N", "R", 0, 2, 0, s, 0, a, 1, info);
        Chkxer("Rgebak", infot, nout, lerr, ok);
        infot = 5;
        Rgebak("N", "R", 2, 2, 1, s, 0, a, 2, info);
        Chkxer("Rgebak", infot, nout, lerr, ok);
        infot = 5;
        Rgebak("N", "R", 0, 1, 1, s, 0, a, 1, info);
        Chkxer("Rgebak", infot, nout, lerr, ok);
        infot = 7;
        Rgebak("N", "R", 0, 1, 0, s, -1, a, 1, info);
        Chkxer("Rgebak", infot, nout, lerr, ok);
        infot = 9;
        Rgebak("N", "R", 2, 1, 2, s, 0, a, 1, info);
        Chkxer("Rgebak", infot, nout, lerr, ok);
        nt += 9;
        //
        // Rgehrd
        //
        srnamt = "Rgehrd";
        infot = 1;
        Rgehrd(-1, 1, 1, a, 1, tau, w, 1, info);
        Chkxer("Rgehrd", infot, nout, lerr, ok);
        infot = 2;
        Rgehrd(0, 0, 0, a, 1, tau, w, 1, info);
        Chkxer("Rgehrd", infot, nout, lerr, ok);
        infot = 2;
        Rgehrd(0, 2, 0, a, 1, tau, w, 1, info);
        Chkxer("Rgehrd", infot, nout, lerr, ok);
        infot = 3;
        Rgehrd(1, 1, 0, a, 1, tau, w, 1, info);
        Chkxer("Rgehrd", infot, nout, lerr, ok);
        infot = 3;
        Rgehrd(0, 1, 1, a, 1, tau, w, 1, info);
        Chkxer("Rgehrd", infot, nout, lerr, ok);
        infot = 5;
        Rgehrd(2, 1, 1, a, 1, tau, w, 2, info);
        Chkxer("Rgehrd", infot, nout, lerr, ok);
        infot = 8;
        Rgehrd(2, 1, 2, a, 2, tau, w, 1, info);
        Chkxer("Rgehrd", infot, nout, lerr, ok);
        nt += 7;
        //
        // Rgehd2
        //
        srnamt = "Rgehd2";
        infot = 1;
        Rgehd2(-1, 1, 1, a, 1, tau, w, info);
        Chkxer("Rgehd2", infot, nout, lerr, ok);
        infot = 2;
        Rgehd2(0, 0, 0, a, 1, tau, w, info);
        Chkxer("Rgehd2", infot, nout, lerr, ok);
        infot = 2;
        Rgehd2(0, 2, 0, a, 1, tau, w, info);
        Chkxer("Rgehd2", infot, nout, lerr, ok);
        infot = 3;
        Rgehd2(1, 1, 0, a, 1, tau, w, info);
        Chkxer("Rgehd2", infot, nout, lerr, ok);
        infot = 3;
        Rgehd2(0, 1, 1, a, 1, tau, w, info);
        Chkxer("Rgehd2", infot, nout, lerr, ok);
        infot = 5;
        Rgehd2(2, 1, 1, a, 1, tau, w, info);
        Chkxer("Rgehd2", infot, nout, lerr, ok);
        nt += 6;
        //
        // Rorghr
        //
        srnamt = "Rorghr";
        infot = 1;
        Rorghr(-1, 1, 1, a, 1, tau, w, 1, info);
        Chkxer("Rorghr", infot, nout, lerr, ok);
        infot = 2;
        Rorghr(0, 0, 0, a, 1, tau, w, 1, info);
        Chkxer("Rorghr", infot, nout, lerr, ok);
        infot = 2;
        Rorghr(0, 2, 0, a, 1, tau, w, 1, info);
        Chkxer("Rorghr", infot, nout, lerr, ok);
        infot = 3;
        Rorghr(1, 1, 0, a, 1, tau, w, 1, info);
        Chkxer("Rorghr", infot, nout, lerr, ok);
        infot = 3;
        Rorghr(0, 1, 1, a, 1, tau, w, 1, info);
        Chkxer("Rorghr", infot, nout, lerr, ok);
        infot = 5;
        Rorghr(2, 1, 1, a, 1, tau, w, 1, info);
        Chkxer("Rorghr", infot, nout, lerr, ok);
        infot = 8;
        Rorghr(3, 1, 3, a, 3, tau, w, 1, info);
        Chkxer("Rorghr", infot, nout, lerr, ok);
        nt += 7;
        //
        // Rormhr
        //
        srnamt = "Rormhr";
        infot = 1;
        Rormhr("/", "N", 0, 0, 1, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 2;
        Rormhr("L", "/", 0, 0, 1, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 3;
        Rormhr("L", "N", -1, 0, 1, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 4;
        Rormhr("L", "N", 0, -1, 1, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 5;
        Rormhr("L", "N", 0, 0, 0, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 5;
        Rormhr("L", "N", 0, 0, 2, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 5;
        Rormhr("L", "N", 1, 2, 2, 1, a, 1, tau, c, 1, w, 2, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 5;
        Rormhr("R", "N", 2, 1, 2, 1, a, 1, tau, c, 2, w, 2, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 6;
        Rormhr("L", "N", 1, 1, 1, 0, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 6;
        Rormhr("L", "N", 0, 1, 1, 1, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 6;
        Rormhr("R", "N", 1, 0, 1, 1, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 8;
        Rormhr("L", "N", 2, 1, 1, 1, a, 1, tau, c, 2, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 8;
        Rormhr("R", "N", 1, 2, 1, 1, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 11;
        Rormhr("L", "N", 2, 1, 1, 1, a, 2, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 13;
        Rormhr("L", "N", 1, 2, 1, 1, a, 1, tau, c, 1, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        infot = 13;
        Rormhr("R", "N", 2, 1, 1, 1, a, 1, tau, c, 2, w, 1, info);
        Chkxer("Rormhr", infot, nout, lerr, ok);
        nt += 16;
        //
        // Rhseqr
        //
        srnamt = "Rhseqr";
        infot = 1;
        Rhseqr("/", "N", 0, 1, 0, a, 1, wr, wi, c, 1, w, 1, info);
        Chkxer("Rhseqr", infot, nout, lerr, ok);
        infot = 2;
        Rhseqr("E", "/", 0, 1, 0, a, 1, wr, wi, c, 1, w, 1, info);
        Chkxer("Rhseqr", infot, nout, lerr, ok);
        infot = 3;
        Rhseqr("E", "N", -1, 1, 0, a, 1, wr, wi, c, 1, w, 1, info);
        Chkxer("Rhseqr", infot, nout, lerr, ok);
        infot = 4;
        Rhseqr("E", "N", 0, 0, 0, a, 1, wr, wi, c, 1, w, 1, info);
        Chkxer("Rhseqr", infot, nout, lerr, ok);
        infot = 4;
        Rhseqr("E", "N", 0, 2, 0, a, 1, wr, wi, c, 1, w, 1, info);
        Chkxer("Rhseqr", infot, nout, lerr, ok);
        infot = 5;
        Rhseqr("E", "N", 1, 1, 0, a, 1, wr, wi, c, 1, w, 1, info);
        Chkxer("Rhseqr", infot, nout, lerr, ok);
        infot = 5;
        Rhseqr("E", "N", 1, 1, 2, a, 1, wr, wi, c, 1, w, 1, info);
        Chkxer("Rhseqr", infot, nout, lerr, ok);
        infot = 7;
        Rhseqr("E", "N", 2, 1, 2, a, 1, wr, wi, c, 2, w, 1, info);
        Chkxer("Rhseqr", infot, nout, lerr, ok);
        infot = 11;
        Rhseqr("E", "V", 2, 1, 2, a, 2, wr, wi, c, 1, w, 1, info);
        Chkxer("Rhseqr", infot, nout, lerr, ok);
        infot = 13;
        Rhseqr("E", "N", 2, 1, 2, a, 2, wr, wi, c, 1, w, 1, info);
        Chkxer("Rhseqr", infot, nout, lerr, ok);
        nt += 10;
        //
        // Rhsein
        //
        srnamt = "Rhsein";
        infot = 1;
        Rhsein("/", "N", "N", sel, 0, a, 1, wr, wi, vl, 1, vr, 1, 0, m, w, ifaill, ifailr, info);
        Chkxer("Rhsein", infot, nout, lerr, ok);
        infot = 2;
        Rhsein("R", "/", "N", sel, 0, a, 1, wr, wi, vl, 1, vr, 1, 0, m, w, ifaill, ifailr, info);
        Chkxer("Rhsein", infot, nout, lerr, ok);
        infot = 3;
        Rhsein("R", "N", "/", sel, 0, a, 1, wr, wi, vl, 1, vr, 1, 0, m, w, ifaill, ifailr, info);
        Chkxer("Rhsein", infot, nout, lerr, ok);
        infot = 5;
        Rhsein("R", "N", "N", sel, -1, a, 1, wr, wi, vl, 1, vr, 1, 0, m, w, ifaill, ifailr, info);
        Chkxer("Rhsein", infot, nout, lerr, ok);
        infot = 7;
        Rhsein("R", "N", "N", sel, 2, a, 1, wr, wi, vl, 1, vr, 2, 4, m, w, ifaill, ifailr, info);
        Chkxer("Rhsein", infot, nout, lerr, ok);
        infot = 11;
        Rhsein("L", "N", "N", sel, 2, a, 2, wr, wi, vl, 1, vr, 1, 4, m, w, ifaill, ifailr, info);
        Chkxer("Rhsein", infot, nout, lerr, ok);
        infot = 13;
        Rhsein("R", "N", "N", sel, 2, a, 2, wr, wi, vl, 1, vr, 1, 4, m, w, ifaill, ifailr, info);
        Chkxer("Rhsein", infot, nout, lerr, ok);
        infot = 14;
        Rhsein("R", "N", "N", sel, 2, a, 2, wr, wi, vl, 1, vr, 2, 1, m, w, ifaill, ifailr, info);
        Chkxer("Rhsein", infot, nout, lerr, ok);
        nt += 8;
        //
        // Rtrevc
        //
        srnamt = "Rtrevc";
        infot = 1;
        Rtrevc("/", "A", sel, 0, a, 1, vl, 1, vr, 1, 0, m, w, info);
        Chkxer("Rtrevc", infot, nout, lerr, ok);
        infot = 2;
        Rtrevc("L", "/", sel, 0, a, 1, vl, 1, vr, 1, 0, m, w, info);
        Chkxer("Rtrevc", infot, nout, lerr, ok);
        infot = 4;
        Rtrevc("L", "A", sel, -1, a, 1, vl, 1, vr, 1, 0, m, w, info);
        Chkxer("Rtrevc", infot, nout, lerr, ok);
        infot = 6;
        Rtrevc("L", "A", sel, 2, a, 1, vl, 2, vr, 1, 4, m, w, info);
        Chkxer("Rtrevc", infot, nout, lerr, ok);
        infot = 8;
        Rtrevc("L", "A", sel, 2, a, 2, vl, 1, vr, 1, 4, m, w, info);
        Chkxer("Rtrevc", infot, nout, lerr, ok);
        infot = 10;
        Rtrevc("R", "A", sel, 2, a, 2, vl, 1, vr, 1, 4, m, w, info);
        Chkxer("Rtrevc", infot, nout, lerr, ok);
        infot = 11;
        Rtrevc("L", "A", sel, 2, a, 2, vl, 2, vr, 1, 1, m, w, info);
        Chkxer("Rtrevc", infot, nout, lerr, ok);
        nt += 7;
        //
        // Rtrevc3
        //
        srnamt = "Rtrevc3";
        infot = 1;
        Rtrevc3("/", "A", sel, 0, a, 1, vl, 1, vr, 1, 0, m, w, lw, info);
        Chkxer("Rtrevc3", infot, nout, lerr, ok);
        infot = 2;
        Rtrevc3("L", "/", sel, 0, a, 1, vl, 1, vr, 1, 0, m, w, lw, info);
        Chkxer("Rtrevc3", infot, nout, lerr, ok);
        infot = 4;
        Rtrevc3("L", "A", sel, -1, a, 1, vl, 1, vr, 1, 0, m, w, lw, info);
        Chkxer("Rtrevc3", infot, nout, lerr, ok);
        infot = 6;
        Rtrevc3("L", "A", sel, 2, a, 1, vl, 2, vr, 1, 4, m, w, lw, info);
        Chkxer("Rtrevc3", infot, nout, lerr, ok);
        infot = 8;
        Rtrevc3("L", "A", sel, 2, a, 2, vl, 1, vr, 1, 4, m, w, lw, info);
        Chkxer("Rtrevc3", infot, nout, lerr, ok);
        infot = 10;
        Rtrevc3("R", "A", sel, 2, a, 2, vl, 1, vr, 1, 4, m, w, lw, info);
        Chkxer("Rtrevc3", infot, nout, lerr, ok);
        infot = 11;
        Rtrevc3("L", "A", sel, 2, a, 2, vl, 2, vr, 1, 1, m, w, lw, info);
        Chkxer("Rtrevc3", infot, nout, lerr, ok);
        infot = 14;
        Rtrevc3("L", "A", sel, 2, a, 2, vl, 2, vr, 1, 2, m, w, 2, info);
        Chkxer("Rtrevc3", infot, nout, lerr, ok);
        nt += 8;
    }
    //
    // Print a summary line.
    //
    if (ok) {
        write(nout, format_9999), path, nt;
    } else {
        write(nout, format_9998), path;
    }
    //
    // End of Rerrhs
    //
}
