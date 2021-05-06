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

void Cerrlqtp(const char *path, INTEGER const nunit) {
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
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 2;
    INTEGER i = 0;
    COMPLEX a[nmax * nmax];
    COMPLEX c[nmax * nmax];
    COMPLEX t[nmax * nmax];
    COMPLEX w[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / COMPLEX(castREAL(i + j), 0.0);
            c[(i - 1) + (j - 1) * ldc] = 1.0 / COMPLEX(castREAL(i + j), 0.0);
            t[(i - 1) + (j - 1) * ldt] = 1.0 / COMPLEX(castREAL(i + j), 0.0);
        }
        w[j - 1] = 0.0f;
    }
    ok = true;
    //
    //     Error exits for TPLQT factorization
    //
    //     Ctplqt
    //
    infot = 1;
    COMPLEX b[nmax * nmax];
    INTEGER info = 0;
    Ctplqt(-1, 1, 0, 1, a, 1, b, 1, t, 1, w, info);
    chkxer("Ctplqt", infot, nout, lerr, ok);
    infot = 2;
    Ctplqt(1, -1, 0, 1, a, 1, b, 1, t, 1, w, info);
    chkxer("Ctplqt", infot, nout, lerr, ok);
    infot = 3;
    Ctplqt(0, 1, -1, 1, a, 1, b, 1, t, 1, w, info);
    chkxer("Ctplqt", infot, nout, lerr, ok);
    infot = 3;
    Ctplqt(0, 1, 1, 1, a, 1, b, 1, t, 1, w, info);
    chkxer("Ctplqt", infot, nout, lerr, ok);
    infot = 4;
    Ctplqt(0, 1, 0, 0, a, 1, b, 1, t, 1, w, info);
    chkxer("Ctplqt", infot, nout, lerr, ok);
    infot = 4;
    Ctplqt(1, 1, 0, 2, a, 1, b, 1, t, 1, w, info);
    chkxer("Ctplqt", infot, nout, lerr, ok);
    infot = 6;
    Ctplqt(2, 1, 0, 2, a, 1, b, 1, t, 1, w, info);
    chkxer("Ctplqt", infot, nout, lerr, ok);
    infot = 8;
    Ctplqt(2, 1, 0, 1, a, 2, b, 1, t, 1, w, info);
    chkxer("Ctplqt", infot, nout, lerr, ok);
    infot = 10;
    Ctplqt(2, 2, 1, 2, a, 2, b, 2, t, 1, w, info);
    chkxer("Ctplqt", infot, nout, lerr, ok);
    //
    //     Ctplqt2
    //
    infot = 1;
    Ctplqt2(-1, 0, 0, a, 1, b, 1, t, 1, info);
    chkxer("Ctplqt2", infot, nout, lerr, ok);
    infot = 2;
    Ctplqt2(0, -1, 0, a, 1, b, 1, t, 1, info);
    chkxer("Ctplqt2", infot, nout, lerr, ok);
    infot = 3;
    Ctplqt2(0, 0, -1, a, 1, b, 1, t, 1, info);
    chkxer("Ctplqt2", infot, nout, lerr, ok);
    infot = 5;
    Ctplqt2(2, 2, 0, a, 1, b, 2, t, 2, info);
    chkxer("Ctplqt2", infot, nout, lerr, ok);
    infot = 7;
    Ctplqt2(2, 2, 0, a, 2, b, 1, t, 2, info);
    chkxer("Ctplqt2", infot, nout, lerr, ok);
    infot = 9;
    Ctplqt2(2, 2, 0, a, 2, b, 2, t, 1, info);
    chkxer("Ctplqt2", infot, nout, lerr, ok);
    //
    //     Ctpmlqt
    //
    infot = 1;
    Ctpmlqt("/", "N", 0, 0, 0, 0, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Ctpmlqt", infot, nout, lerr, ok);
    infot = 2;
    Ctpmlqt("L", "/", 0, 0, 0, 0, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Ctpmlqt", infot, nout, lerr, ok);
    infot = 3;
    Ctpmlqt("L", "N", -1, 0, 0, 0, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Ctpmlqt", infot, nout, lerr, ok);
    infot = 4;
    Ctpmlqt("L", "N", 0, -1, 0, 0, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Ctpmlqt", infot, nout, lerr, ok);
    infot = 5;
    Ctpmlqt("L", "N", 0, 0, -1, 0, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    infot = 6;
    Ctpmlqt("L", "N", 0, 0, 0, -1, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Ctpmlqt", infot, nout, lerr, ok);
    infot = 7;
    Ctpmlqt("L", "N", 0, 0, 0, 0, 0, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Ctpmlqt", infot, nout, lerr, ok);
    infot = 9;
    Ctpmlqt("R", "N", 2, 2, 2, 1, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Ctpmlqt", infot, nout, lerr, ok);
    infot = 11;
    Ctpmlqt("R", "N", 1, 1, 1, 1, 1, a, 1, t, 0, b, 1, c, 1, w, info);
    chkxer("Ctpmlqt", infot, nout, lerr, ok);
    infot = 13;
    Ctpmlqt("L", "N", 1, 1, 1, 1, 1, a, 1, t, 1, b, 0, c, 1, w, info);
    chkxer("Ctpmlqt", infot, nout, lerr, ok);
    infot = 15;
    Ctpmlqt("L", "N", 1, 1, 1, 1, 1, a, 1, t, 1, b, 1, c, 0, w, info);
    chkxer("Ctpmlqt", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrlqt
    //
}
