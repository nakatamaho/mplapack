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

void Cerrbd(const char *path, INTEGER const nunit) {
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
    const INTEGER nmax = 4;
    INTEGER i = 0;
    arr_2d<nmax, nmax, COMPLEX> a;
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / (i + j).real();
        }
    }
    ok = true;
    INTEGER nt = 0;
    //
    //     Test error exits of the SVD routines.
    //
    arr_1d<nmax, REAL> d;
    arr_1d<nmax, REAL> e;
    arr_1d<nmax, COMPLEX> tq;
    arr_1d<nmax, COMPLEX> tp;
    const INTEGER lw = nmax;
    arr_1d<lw, COMPLEX> w;
    INTEGER info = 0;
    arr_2d<nmax, nmax, COMPLEX> u;
    arr_2d<nmax, nmax, COMPLEX> v;
    arr_1d<4 * nmax, REAL> rw;
    if (Mlsamen(2, c2, "BD")) {
        //
        //        Cgebrd
        //
        srnamt = "Cgebrd";
        infot = 1;
        Cgebrd(-1, 0, a, 1, d, e, tq, tp, w, 1, info);
        chkxer("Cgebrd", infot, nout, lerr, ok);
        infot = 2;
        Cgebrd(0, -1, a, 1, d, e, tq, tp, w, 1, info);
        chkxer("Cgebrd", infot, nout, lerr, ok);
        infot = 4;
        Cgebrd(2, 1, a, 1, d, e, tq, tp, w, 2, info);
        chkxer("Cgebrd", infot, nout, lerr, ok);
        infot = 10;
        Cgebrd(2, 1, a, 2, d, e, tq, tp, w, 1, info);
        chkxer("Cgebrd", infot, nout, lerr, ok);
        nt += 4;
        //
        //        Cungbr
        //
        srnamt = "Cungbr";
        infot = 1;
        Cungbr("/", 0, 0, 0, a, 1, tq, w, 1, info);
        chkxer("Cungbr", infot, nout, lerr, ok);
        infot = 2;
        Cungbr("Q", -1, 0, 0, a, 1, tq, w, 1, info);
        chkxer("Cungbr", infot, nout, lerr, ok);
        infot = 3;
        Cungbr("Q", 0, -1, 0, a, 1, tq, w, 1, info);
        chkxer("Cungbr", infot, nout, lerr, ok);
        infot = 3;
        Cungbr("Q", 0, 1, 0, a, 1, tq, w, 1, info);
        chkxer("Cungbr", infot, nout, lerr, ok);
        infot = 3;
        Cungbr("Q", 1, 0, 1, a, 1, tq, w, 1, info);
        chkxer("Cungbr", infot, nout, lerr, ok);
        infot = 3;
        Cungbr("P", 1, 0, 0, a, 1, tq, w, 1, info);
        chkxer("Cungbr", infot, nout, lerr, ok);
        infot = 3;
        Cungbr("P", 0, 1, 1, a, 1, tq, w, 1, info);
        chkxer("Cungbr", infot, nout, lerr, ok);
        infot = 4;
        Cungbr("Q", 0, 0, -1, a, 1, tq, w, 1, info);
        chkxer("Cungbr", infot, nout, lerr, ok);
        infot = 6;
        Cungbr("Q", 2, 1, 1, a, 1, tq, w, 1, info);
        chkxer("Cungbr", infot, nout, lerr, ok);
        infot = 9;
        Cungbr("Q", 2, 2, 1, a, 2, tq, w, 1, info);
        chkxer("Cungbr", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Cunmbr
        //
        srnamt = "Cunmbr";
        infot = 1;
        Cunmbr("/", "L", "T", 0, 0, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 2;
        Cunmbr("Q", "/", "T", 0, 0, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 3;
        Cunmbr("Q", "L", "/", 0, 0, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 4;
        Cunmbr("Q", "L", "C", -1, 0, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 5;
        Cunmbr("Q", "L", "C", 0, -1, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 6;
        Cunmbr("Q", "L", "C", 0, 0, -1, a, 1, tq, u, 1, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 8;
        Cunmbr("Q", "L", "C", 2, 0, 0, a, 1, tq, u, 2, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 8;
        Cunmbr("Q", "R", "C", 0, 2, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 8;
        Cunmbr("P", "L", "C", 2, 0, 2, a, 1, tq, u, 2, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 8;
        Cunmbr("P", "R", "C", 0, 2, 2, a, 1, tq, u, 1, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 11;
        Cunmbr("Q", "R", "C", 2, 0, 0, a, 1, tq, u, 1, w, 1, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 13;
        Cunmbr("Q", "L", "C", 0, 2, 0, a, 1, tq, u, 1, w, 0, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        infot = 13;
        Cunmbr("Q", "R", "C", 2, 0, 0, a, 1, tq, u, 2, w, 0, info);
        chkxer("Cunmbr", infot, nout, lerr, ok);
        nt += 13;
        //
        //        Cbdsqr
        //
        srnamt = "Cbdsqr";
        infot = 1;
        Cbdsqr("/", 0, 0, 0, 0, d, e, v, 1, u, 1, a, 1, rw, info);
        chkxer("Cbdsqr", infot, nout, lerr, ok);
        infot = 2;
        Cbdsqr("U", -1, 0, 0, 0, d, e, v, 1, u, 1, a, 1, rw, info);
        chkxer("Cbdsqr", infot, nout, lerr, ok);
        infot = 3;
        Cbdsqr("U", 0, -1, 0, 0, d, e, v, 1, u, 1, a, 1, rw, info);
        chkxer("Cbdsqr", infot, nout, lerr, ok);
        infot = 4;
        Cbdsqr("U", 0, 0, -1, 0, d, e, v, 1, u, 1, a, 1, rw, info);
        chkxer("Cbdsqr", infot, nout, lerr, ok);
        infot = 5;
        Cbdsqr("U", 0, 0, 0, -1, d, e, v, 1, u, 1, a, 1, rw, info);
        chkxer("Cbdsqr", infot, nout, lerr, ok);
        infot = 9;
        Cbdsqr("U", 2, 1, 0, 0, d, e, v, 1, u, 1, a, 1, rw, info);
        chkxer("Cbdsqr", infot, nout, lerr, ok);
        infot = 11;
        Cbdsqr("U", 0, 0, 2, 0, d, e, v, 1, u, 1, a, 1, rw, info);
        chkxer("Cbdsqr", infot, nout, lerr, ok);
        infot = 13;
        Cbdsqr("U", 2, 0, 0, 1, d, e, v, 1, u, 1, a, 1, rw, info);
        chkxer("Cbdsqr", infot, nout, lerr, ok);
        nt += 8;
    }
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a3,' routines passed the tests of the error exits (',i3,"
                    "' tests done)')"),
            path, nt;
    } else {
        write(nout, "(' *** ',a3,' routines failed the tests of the error ','exits ***')"), path;
    }
    //
    //     End of Cerrbd
    //
}
