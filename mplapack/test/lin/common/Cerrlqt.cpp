/*
 * Copyright (c) 2008-2021
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

void Cerrlqt(common &cmn, const char *path, INTEGER const nunit) {
    common_write write(cmn);
    // COMMON infoc
    INTEGER &infot = cmn.infot;
    INTEGER &nout = cmn.nout;
    bool &ok = cmn.ok;
    bool &lerr = cmn.lerr;
    // COMMON srnamc
    str<32> &srnamt = cmn.srnamt;
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
    arr_2d<nmax, nmax, COMPLEX> a(fill0);
    arr_2d<nmax, nmax, COMPLEX> c(fill0);
    arr_2d<nmax, nmax, COMPLEX> t(fill0);
    arr_1d<nmax, COMPLEX> w(fill0);
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / COMPLEX((i + j).real(), 0.0);
            c[(i - 1) + (j - 1) * ldc] = 1.0 / COMPLEX((i + j).real(), 0.0);
            t[(i - 1) + (j - 1) * ldt] = 1.0 / COMPLEX((i + j).real(), 0.0);
        }
        w[j - 1] = 0.0;
    }
    ok = true;
    //
    //     Error exits for LQT factorization
    //
    //     ZGELQT
    //
    srnamt = "ZGELQT";
    infot = 1;
    INTEGER info = 0;
    zgelqt(-1, 0, 1, a, 1, t, 1, w, info);
    chkxer("ZGELQT", infot, nout, lerr, ok);
    infot = 2;
    zgelqt(0, -1, 1, a, 1, t, 1, w, info);
    chkxer("ZGELQT", infot, nout, lerr, ok);
    infot = 3;
    zgelqt(0, 0, 0, a, 1, t, 1, w, info);
    chkxer("ZGELQT", infot, nout, lerr, ok);
    infot = 5;
    zgelqt(2, 1, 1, a, 1, t, 1, w, info);
    chkxer("ZGELQT", infot, nout, lerr, ok);
    infot = 7;
    zgelqt(2, 2, 2, a, 2, t, 1, w, info);
    chkxer("ZGELQT", infot, nout, lerr, ok);
    //
    //     ZGELQT3
    //
    srnamt = "ZGELQT3";
    infot = 1;
    zgelqt3(-1, 0, a, 1, t, 1, info);
    chkxer("ZGELQT3", infot, nout, lerr, ok);
    infot = 2;
    zgelqt3(0, -1, a, 1, t, 1, info);
    chkxer("ZGELQT3", infot, nout, lerr, ok);
    infot = 4;
    zgelqt3(2, 2, a, 1, t, 1, info);
    chkxer("ZGELQT3", infot, nout, lerr, ok);
    infot = 6;
    zgelqt3(2, 2, a, 2, t, 1, info);
    chkxer("ZGELQT3", infot, nout, lerr, ok);
    //
    //     ZGEMLQT
    //
    srnamt = "ZGEMLQT";
    infot = 1;
    zgemlqt("/", "N", 0, 0, 0, 1, a, 1, t, 1, c, 1, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    infot = 2;
    zgemlqt("L", "/", 0, 0, 0, 1, a, 1, t, 1, c, 1, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    infot = 3;
    zgemlqt("L", "N", -1, 0, 0, 1, a, 1, t, 1, c, 1, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    infot = 4;
    zgemlqt("L", "N", 0, -1, 0, 1, a, 1, t, 1, c, 1, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    infot = 5;
    zgemlqt("L", "N", 0, 0, -1, 1, a, 1, t, 1, c, 1, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    infot = 5;
    zgemlqt("R", "N", 0, 0, -1, 1, a, 1, t, 1, c, 1, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    infot = 6;
    zgemlqt("L", "N", 0, 0, 0, 0, a, 1, t, 1, c, 1, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    infot = 8;
    zgemlqt("R", "N", 2, 2, 2, 1, a, 1, t, 1, c, 1, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    infot = 8;
    zgemlqt("L", "N", 2, 2, 2, 1, a, 1, t, 1, c, 1, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    infot = 10;
    zgemlqt("R", "N", 1, 1, 1, 1, a, 1, t, 0, c, 1, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    infot = 12;
    zgemlqt("L", "N", 1, 1, 1, 1, a, 1, t, 1, c, 0, w, info);
    chkxer("ZGEMLQT", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrlqt
    //
}
