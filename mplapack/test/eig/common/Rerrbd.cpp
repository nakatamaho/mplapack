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

void Rerrbd(const char *path, INTEGER const nunit) {
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
    const INTEGER nmax = 4;
    INTEGER i = 0;
    REAL a[nmax * nmax];
    INTEGER lda = nmax;
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / castREAL(i + j);
        }
    }
    ok = true;
    INTEGER nt = 0;
    //
    //     Test error exits of the SVD routines.
    //
    REAL d[nmax];
    REAL e[nmax];
    REAL tq[nmax];
    REAL tp[nmax];
    const INTEGER lw = nmax;
    REAL w[lw];
    INTEGER info = 0;
    REAL u[nmax * nmax];
    REAL v[nmax * nmax];
    REAL q[nmax * nmax];
    INTEGER iq[nmax * nmax];
    INTEGER iw[nmax];
    const REAL zero = 0.0;
    const REAL one = 1.0;
    INTEGER ns = 0;
    REAL s[nmax];
    if (Mlsamen(2, c2, "BD")) {
        //
        //        Rgebrd
        //
        infot = 1;
        strncpy(srnamt, "Rgebrd", srnamt_len);
        Rgebrd(-1, 0, a, 1, d, e, tq, tp, w, 1, info);
        chkxer("Rgebrd", infot, nout, lerr, ok);
        infot = 2;
        Rgebrd(0, -1, a, 1, d, e, tq, tp, w, 1, info);
        chkxer("Rgebrd", infot, nout, lerr, ok);
        infot = 4;
        Rgebrd(2, 1, a, 1, d, e, tq, tp, w, 2, info);
        chkxer("Rgebrd", infot, nout, lerr, ok);
        infot = 10;
        Rgebrd(2, 1, a, 2, d, e, tq, tp, w, 1, info);
        chkxer("Rgebrd", infot, nout, lerr, ok);
        nt += 4;
        //
        //        Rgebd2
        //
        infot = 1;
        Rgebd2(-1, 0, a, 1, d, e, tq, tp, w, info);
        chkxer("Rgebd2", infot, nout, lerr, ok);
        infot = 2;
        Rgebd2(0, -1, a, 1, d, e, tq, tp, w, info);
        chkxer("Rgebd2", infot, nout, lerr, ok);
        infot = 4;
        Rgebd2(2, 1, a, 1, d, e, tq, tp, w, info);
        chkxer("Rgebd2", infot, nout, lerr, ok);
        nt += 3;
        //
        //        Rorgbr
        //
        infot = 1;
        Rorgbr("/", 0, 0, 0, a, 1, tq, w, 1, info);
        chkxer("Rorgbr", infot, nout, lerr, ok);
        infot = 2;
        Rorgbr("Q", -1, 0, 0, a, 1, tq, w, 1, info);
        chkxer("Rorgbr", infot, nout, lerr, ok);
        infot = 3;
        Rorgbr("Q", 0, -1, 0, a, 1, tq, w, 1, info);
        chkxer("Rorgbr", infot, nout, lerr, ok);
        infot = 3;
        Rorgbr("Q", 0, 1, 0, a, 1, tq, w, 1, info);
        chkxer("Rorgbr", infot, nout, lerr, ok);
        infot = 3;
        Rorgbr("Q", 1, 0, 1, a, 1, tq, w, 1, info);
        chkxer("Rorgbr", infot, nout, lerr, ok);
        infot = 3;
        Rorgbr("P", 1, 0, 0, a, 1, tq, w, 1, info);
        chkxer("Rorgbr", infot, nout, lerr, ok);
        infot = 3;
        Rorgbr("P", 0, 1, 1, a, 1, tq, w, 1, info);
        chkxer("Rorgbr", infot, nout, lerr, ok);
        infot = 4;
        Rorgbr("Q", 0, 0, -1, a, 1, tq, w, 1, info);
        chkxer("Rorgbr", infot, nout, lerr, ok);
        infot = 6;
        Rorgbr("Q", 2, 1, 1, a, 1, tq, w, 1, info);
        chkxer("Rorgbr", infot, nout, lerr, ok);
        infot = 9;
        Rorgbr("Q", 2, 2, 1, a, 2, tq, w, 1, info);
        chkxer("Rorgbr", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Rormbr
        //
        infot = 1;
        Rormbr("/", "L", "T", 0, 0, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 2;
        Rormbr("Q", "/", "T", 0, 0, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 3;
        Rormbr("Q", "L", "/", 0, 0, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 4;
        Rormbr("Q", "L", "T", -1, 0, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 5;
        Rormbr("Q", "L", "T", 0, -1, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 6;
        Rormbr("Q", "L", "T", 0, 0, -1, a, 1, tq, u, 1, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 8;
        Rormbr("Q", "L", "T", 2, 0, 0, a, 1, tq, u, 2, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 8;
        Rormbr("Q", "R", "T", 0, 2, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 8;
        Rormbr("P", "L", "T", 2, 0, 2, a, 1, tq, u, 2, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 8;
        Rormbr("P", "R", "T", 0, 2, 2, a, 1, tq, u, 1, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 11;
        Rormbr("Q", "R", "T", 2, 0, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 13;
        Rormbr("Q", "L", "T", 0, 2, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        infot = 13;
        Rormbr("Q", "R", "T", 2, 0, 0, a, 1, tq, u, 2, w, 1, info);
        chkxer("Rormbr", infot, nout, lerr, ok);
        nt += 13;
        //
        //        Rbdsqr
        //
        infot = 1;
        Rbdsqr("/", 0, 0, 0, 0, d, e, v, 1, u, 1, a, 1, w, info);
        chkxer("Rbdsqr", infot, nout, lerr, ok);
        infot = 2;
        Rbdsqr("U", -1, 0, 0, 0, d, e, v, 1, u, 1, a, 1, w, info);
        chkxer("Rbdsqr", infot, nout, lerr, ok);
        infot = 3;
        Rbdsqr("U", 0, -1, 0, 0, d, e, v, 1, u, 1, a, 1, w, info);
        chkxer("Rbdsqr", infot, nout, lerr, ok);
        infot = 4;
        Rbdsqr("U", 0, 0, -1, 0, d, e, v, 1, u, 1, a, 1, w, info);
        chkxer("Rbdsqr", infot, nout, lerr, ok);
        infot = 5;
        Rbdsqr("U", 0, 0, 0, -1, d, e, v, 1, u, 1, a, 1, w, info);
        chkxer("Rbdsqr", infot, nout, lerr, ok);
        infot = 9;
        Rbdsqr("U", 2, 1, 0, 0, d, e, v, 1, u, 1, a, 1, w, info);
        chkxer("Rbdsqr", infot, nout, lerr, ok);
        infot = 11;
        Rbdsqr("U", 0, 0, 2, 0, d, e, v, 1, u, 1, a, 1, w, info);
        chkxer("Rbdsqr", infot, nout, lerr, ok);
        infot = 13;
        Rbdsqr("U", 2, 0, 0, 1, d, e, v, 1, u, 1, a, 1, w, info);
        chkxer("Rbdsqr", infot, nout, lerr, ok);
        nt += 8;
        //
        //        Rbdsdc
        //
        infot = 1;
        Rbdsdc("/", "N", 0, d, e, u, 1, v, 1, q, iq, w, iw, info);
        chkxer("Rbdsdc", infot, nout, lerr, ok);
        infot = 2;
        Rbdsdc("U", "/", 0, d, e, u, 1, v, 1, q, iq, w, iw, info);
        chkxer("Rbdsdc", infot, nout, lerr, ok);
        infot = 3;
        Rbdsdc("U", "N", -1, d, e, u, 1, v, 1, q, iq, w, iw, info);
        chkxer("Rbdsdc", infot, nout, lerr, ok);
        infot = 7;
        Rbdsdc("U", "I", 2, d, e, u, 1, v, 1, q, iq, w, iw, info);
        chkxer("Rbdsdc", infot, nout, lerr, ok);
        infot = 9;
        Rbdsdc("U", "I", 2, d, e, u, 2, v, 1, q, iq, w, iw, info);
        chkxer("Rbdsdc", infot, nout, lerr, ok);
        nt += 5;
        //
        //        Rbdsvdx
        //
        infot = 1;
        Rbdsvdx("X", "N", "A", 1, d, e, zero, one, 0, 0, ns, s, q, 1, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 2;
        Rbdsvdx("U", "X", "A", 1, d, e, zero, one, 0, 0, ns, s, q, 1, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 3;
        Rbdsvdx("U", "V", "X", 1, d, e, zero, one, 0, 0, ns, s, q, 1, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 4;
        Rbdsvdx("U", "V", "A", -1, d, e, zero, one, 0, 0, ns, s, q, 1, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 7;
        Rbdsvdx("U", "V", "V", 2, d, e, -one, zero, 0, 0, ns, s, q, 1, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 8;
        Rbdsvdx("U", "V", "V", 2, d, e, one, zero, 0, 0, ns, s, q, 1, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 9;
        Rbdsvdx("L", "V", "I", 2, d, e, zero, zero, 0, 2, ns, s, q, 1, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 9;
        Rbdsvdx("L", "V", "I", 4, d, e, zero, zero, 5, 2, ns, s, q, 1, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 10;
        Rbdsvdx("L", "V", "I", 4, d, e, zero, zero, 3, 2, ns, s, q, 1, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 10;
        Rbdsvdx("L", "V", "I", 4, d, e, zero, zero, 3, 5, ns, s, q, 1, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 14;
        Rbdsvdx("L", "V", "A", 4, d, e, zero, zero, 0, 0, ns, s, q, 0, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
        infot = 14;
        Rbdsvdx("L", "V", "A", 4, d, e, zero, zero, 0, 0, ns, s, q, 2, w, iw, info);
        chkxer("Rbdsvdx", infot, nout, lerr, ok);
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
    //     End of Rerrbd
    //
}
