/*
 * Copyright (c) 2021-2022
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

void Cerrhs(const char *path, INTEGER const nunit) {
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
    COMPLEX a[nmax * nmax];
    INTEGER lda = nmax;
    bool sel[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / castREAL(i + j);
        }
        sel[j - 1] = true;
    }
    ok = true;
    INTEGER nt = 0;
    //
    //     Test error exits of the nonsymmetric eigenvalue routines.
    //
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL s[nmax];
    INTEGER info = 0;
    COMPLEX tau[nmax];
    const INTEGER lw = nmax * nmax;
    COMPLEX w[lw];
    COMPLEX c[nmax * nmax];
    COMPLEX x[nmax];
    COMPLEX vl[nmax * nmax];
    COMPLEX vr[nmax * nmax];
    INTEGER m = 0;
    REAL rw[nmax];
    INTEGER ifaill[nmax];
    INTEGER ifailr[nmax];
    if (Mlsamen(2, c2, "HS")) {
        //
        //        Cgebal
        //
        strncpy(srnamt, "Cgebal", srnamt_len);
        infot = 1;
        Cgebal("/", 0, a, 1, ilo, ihi, s, info);
        chkxer("Cgebal", infot, nout, lerr, ok);
        infot = 2;
        Cgebal("N", -1, a, 1, ilo, ihi, s, info);
        chkxer("Cgebal", infot, nout, lerr, ok);
        infot = 4;
        Cgebal("N", 2, a, 1, ilo, ihi, s, info);
        chkxer("Cgebal", infot, nout, lerr, ok);
        nt += 3;
        //
        //        Cgebak
        //
        strncpy(srnamt, "Cgebak", srnamt_len);
        infot = 1;
        Cgebak("/", "R", 0, 1, 0, s, 0, a, 1, info);
        chkxer("Cgebak", infot, nout, lerr, ok);
        infot = 2;
        Cgebak("N", "/", 0, 1, 0, s, 0, a, 1, info);
        chkxer("Cgebak", infot, nout, lerr, ok);
        infot = 3;
        Cgebak("N", "R", -1, 1, 0, s, 0, a, 1, info);
        chkxer("Cgebak", infot, nout, lerr, ok);
        infot = 4;
        Cgebak("N", "R", 0, 0, 0, s, 0, a, 1, info);
        chkxer("Cgebak", infot, nout, lerr, ok);
        infot = 4;
        Cgebak("N", "R", 0, 2, 0, s, 0, a, 1, info);
        chkxer("Cgebak", infot, nout, lerr, ok);
        infot = 5;
        Cgebak("N", "R", 2, 2, 1, s, 0, a, 2, info);
        chkxer("Cgebak", infot, nout, lerr, ok);
        infot = 5;
        Cgebak("N", "R", 0, 1, 1, s, 0, a, 1, info);
        chkxer("Cgebak", infot, nout, lerr, ok);
        infot = 7;
        Cgebak("N", "R", 0, 1, 0, s, -1, a, 1, info);
        chkxer("Cgebak", infot, nout, lerr, ok);
        infot = 9;
        Cgebak("N", "R", 2, 1, 2, s, 0, a, 1, info);
        chkxer("Cgebak", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Cgehrd
        //
        strncpy(srnamt, "Cgehrd", srnamt_len);
        infot = 1;
        Cgehrd(-1, 1, 1, a, 1, tau, w, 1, info);
        chkxer("Cgehrd", infot, nout, lerr, ok);
        infot = 2;
        Cgehrd(0, 0, 0, a, 1, tau, w, 1, info);
        chkxer("Cgehrd", infot, nout, lerr, ok);
        infot = 2;
        Cgehrd(0, 2, 0, a, 1, tau, w, 1, info);
        chkxer("Cgehrd", infot, nout, lerr, ok);
        infot = 3;
        Cgehrd(1, 1, 0, a, 1, tau, w, 1, info);
        chkxer("Cgehrd", infot, nout, lerr, ok);
        infot = 3;
        Cgehrd(0, 1, 1, a, 1, tau, w, 1, info);
        chkxer("Cgehrd", infot, nout, lerr, ok);
        infot = 5;
        Cgehrd(2, 1, 1, a, 1, tau, w, 2, info);
        chkxer("Cgehrd", infot, nout, lerr, ok);
        infot = 8;
        Cgehrd(2, 1, 2, a, 2, tau, w, 1, info);
        chkxer("Cgehrd", infot, nout, lerr, ok);
        nt += 7;
        //
        //        Cunghr
        //
        strncpy(srnamt, "Cunghr", srnamt_len);
        infot = 1;
        Cunghr(-1, 1, 1, a, 1, tau, w, 1, info);
        chkxer("Cunghr", infot, nout, lerr, ok);
        infot = 2;
        Cunghr(0, 0, 0, a, 1, tau, w, 1, info);
        chkxer("Cunghr", infot, nout, lerr, ok);
        infot = 2;
        Cunghr(0, 2, 0, a, 1, tau, w, 1, info);
        chkxer("Cunghr", infot, nout, lerr, ok);
        infot = 3;
        Cunghr(1, 1, 0, a, 1, tau, w, 1, info);
        chkxer("Cunghr", infot, nout, lerr, ok);
        infot = 3;
        Cunghr(0, 1, 1, a, 1, tau, w, 1, info);
        chkxer("Cunghr", infot, nout, lerr, ok);
        infot = 5;
        Cunghr(2, 1, 1, a, 1, tau, w, 1, info);
        chkxer("Cunghr", infot, nout, lerr, ok);
        infot = 8;
        Cunghr(3, 1, 3, a, 3, tau, w, 1, info);
        chkxer("Cunghr", infot, nout, lerr, ok);
        nt += 7;
        //
        //        Cunmhr
        //
        strncpy(srnamt, "Cunmhr", srnamt_len);
        infot = 1;
        Cunmhr("/", "N", 0, 0, 1, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 2;
        Cunmhr("L", "/", 0, 0, 1, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 3;
        Cunmhr("L", "N", -1, 0, 1, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 4;
        Cunmhr("L", "N", 0, -1, 1, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 5;
        Cunmhr("L", "N", 0, 0, 0, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 5;
        Cunmhr("L", "N", 0, 0, 2, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 5;
        Cunmhr("L", "N", 1, 2, 2, 1, a, 1, tau, c, 1, w, 2, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 5;
        Cunmhr("R", "N", 2, 1, 2, 1, a, 1, tau, c, 2, w, 2, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 6;
        Cunmhr("L", "N", 1, 1, 1, 0, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 6;
        Cunmhr("L", "N", 0, 1, 1, 1, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 6;
        Cunmhr("R", "N", 1, 0, 1, 1, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 8;
        Cunmhr("L", "N", 2, 1, 1, 1, a, 1, tau, c, 2, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 8;
        Cunmhr("R", "N", 1, 2, 1, 1, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 11;
        Cunmhr("L", "N", 2, 1, 1, 1, a, 2, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 13;
        Cunmhr("L", "N", 1, 2, 1, 1, a, 1, tau, c, 1, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        infot = 13;
        Cunmhr("R", "N", 2, 1, 1, 1, a, 1, tau, c, 2, w, 1, info);
        chkxer("Cunmhr", infot, nout, lerr, ok);
        nt += 16;
        //
        //        Chseqr
        //
        strncpy(srnamt, "Chseqr", srnamt_len);
        infot = 1;
        Chseqr("/", "N", 0, 1, 0, a, 1, x, c, 1, w, 1, info);
        chkxer("Chseqr", infot, nout, lerr, ok);
        infot = 2;
        Chseqr("E", "/", 0, 1, 0, a, 1, x, c, 1, w, 1, info);
        chkxer("Chseqr", infot, nout, lerr, ok);
        infot = 3;
        Chseqr("E", "N", -1, 1, 0, a, 1, x, c, 1, w, 1, info);
        chkxer("Chseqr", infot, nout, lerr, ok);
        infot = 4;
        Chseqr("E", "N", 0, 0, 0, a, 1, x, c, 1, w, 1, info);
        chkxer("Chseqr", infot, nout, lerr, ok);
        infot = 4;
        Chseqr("E", "N", 0, 2, 0, a, 1, x, c, 1, w, 1, info);
        chkxer("Chseqr", infot, nout, lerr, ok);
        infot = 5;
        Chseqr("E", "N", 1, 1, 0, a, 1, x, c, 1, w, 1, info);
        chkxer("Chseqr", infot, nout, lerr, ok);
        infot = 5;
        Chseqr("E", "N", 1, 1, 2, a, 1, x, c, 1, w, 1, info);
        chkxer("Chseqr", infot, nout, lerr, ok);
        infot = 7;
        Chseqr("E", "N", 2, 1, 2, a, 1, x, c, 2, w, 1, info);
        chkxer("Chseqr", infot, nout, lerr, ok);
        infot = 10;
        Chseqr("E", "V", 2, 1, 2, a, 2, x, c, 1, w, 1, info);
        chkxer("Chseqr", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Chsein
        //
        strncpy(srnamt, "Chsein", srnamt_len);
        infot = 1;
        Chsein("/", "N", "N", sel, 0, a, 1, x, vl, 1, vr, 1, 0, m, w, rw, ifaill, ifailr, info);
        chkxer("Chsein", infot, nout, lerr, ok);
        infot = 2;
        Chsein("R", "/", "N", sel, 0, a, 1, x, vl, 1, vr, 1, 0, m, w, rw, ifaill, ifailr, info);
        chkxer("Chsein", infot, nout, lerr, ok);
        infot = 3;
        Chsein("R", "N", "/", sel, 0, a, 1, x, vl, 1, vr, 1, 0, m, w, rw, ifaill, ifailr, info);
        chkxer("Chsein", infot, nout, lerr, ok);
        infot = 5;
        Chsein("R", "N", "N", sel, -1, a, 1, x, vl, 1, vr, 1, 0, m, w, rw, ifaill, ifailr, info);
        chkxer("Chsein", infot, nout, lerr, ok);
        infot = 7;
        Chsein("R", "N", "N", sel, 2, a, 1, x, vl, 1, vr, 2, 4, m, w, rw, ifaill, ifailr, info);
        chkxer("Chsein", infot, nout, lerr, ok);
        infot = 10;
        Chsein("L", "N", "N", sel, 2, a, 2, x, vl, 1, vr, 1, 4, m, w, rw, ifaill, ifailr, info);
        chkxer("Chsein", infot, nout, lerr, ok);
        infot = 12;
        Chsein("R", "N", "N", sel, 2, a, 2, x, vl, 1, vr, 1, 4, m, w, rw, ifaill, ifailr, info);
        chkxer("Chsein", infot, nout, lerr, ok);
        infot = 13;
        Chsein("R", "N", "N", sel, 2, a, 2, x, vl, 1, vr, 2, 1, m, w, rw, ifaill, ifailr, info);
        chkxer("Chsein", infot, nout, lerr, ok);
        nt += 8;
        //
        //        Ctrevc
        //
        strncpy(srnamt, "Ctrevc", srnamt_len);
        infot = 1;
        Ctrevc("/", "A", sel, 0, a, 1, vl, 1, vr, 1, 0, m, w, rw, info);
        chkxer("Ctrevc", infot, nout, lerr, ok);
        infot = 2;
        Ctrevc("L", "/", sel, 0, a, 1, vl, 1, vr, 1, 0, m, w, rw, info);
        chkxer("Ctrevc", infot, nout, lerr, ok);
        infot = 4;
        Ctrevc("L", "A", sel, -1, a, 1, vl, 1, vr, 1, 0, m, w, rw, info);
        chkxer("Ctrevc", infot, nout, lerr, ok);
        infot = 6;
        Ctrevc("L", "A", sel, 2, a, 1, vl, 2, vr, 1, 4, m, w, rw, info);
        chkxer("Ctrevc", infot, nout, lerr, ok);
        infot = 8;
        Ctrevc("L", "A", sel, 2, a, 2, vl, 1, vr, 1, 4, m, w, rw, info);
        chkxer("Ctrevc", infot, nout, lerr, ok);
        infot = 10;
        Ctrevc("R", "A", sel, 2, a, 2, vl, 1, vr, 1, 4, m, w, rw, info);
        chkxer("Ctrevc", infot, nout, lerr, ok);
        infot = 11;
        Ctrevc("L", "A", sel, 2, a, 2, vl, 2, vr, 1, 1, m, w, rw, info);
        chkxer("Ctrevc", infot, nout, lerr, ok);
        nt += 7;
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
    //     End of Cerrhs
    //
}
