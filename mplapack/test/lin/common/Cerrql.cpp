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

void Cerrql(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_2d<nmax, nmax, COMPLEX> af(fill0);
    arr_1d<nmax, COMPLEX> b(fill0);
    arr_1d<nmax, COMPLEX> w(fill0);
    arr_1d<nmax, COMPLEX> x(fill0);
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = COMPLEX(1.0 / (i + j).real(), -1.0 / (i + j).real());
            af[(i - 1) + (j - 1) * ldaf] = COMPLEX(1.0 / (i + j).real(), -1.0 / (i + j).real());
        }
        b[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
    }
    ok = true;
    //
    //     Error exits for QL factorization
    //
    //     ZGEQLF
    //
    srnamt = "ZGEQLF";
    infot = 1;
    INTEGER info = 0;
    zgeqlf(-1, 0, a, 1, b, w, 1, info);
    chkxer("ZGEQLF", infot, nout, lerr, ok);
    infot = 2;
    zgeqlf(0, -1, a, 1, b, w, 1, info);
    chkxer("ZGEQLF", infot, nout, lerr, ok);
    infot = 4;
    zgeqlf(2, 1, a, 1, b, w, 1, info);
    chkxer("ZGEQLF", infot, nout, lerr, ok);
    infot = 7;
    zgeqlf(1, 2, a, 1, b, w, 1, info);
    chkxer("ZGEQLF", infot, nout, lerr, ok);
    //
    //     ZGEQL2
    //
    srnamt = "ZGEQL2";
    infot = 1;
    zgeql2(-1, 0, a, 1, b, w, info);
    chkxer("ZGEQL2", infot, nout, lerr, ok);
    infot = 2;
    zgeql2(0, -1, a, 1, b, w, info);
    chkxer("ZGEQL2", infot, nout, lerr, ok);
    infot = 4;
    zgeql2(2, 1, a, 1, b, w, info);
    chkxer("ZGEQL2", infot, nout, lerr, ok);
    //
    //     Cgeqls
    //
    srnamt = "Cgeqls";
    infot = 1;
    Cgeqls(-1, 0, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgeqls", infot, nout, lerr, ok);
    infot = 2;
    Cgeqls(0, -1, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgeqls", infot, nout, lerr, ok);
    infot = 2;
    Cgeqls(1, 2, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgeqls", infot, nout, lerr, ok);
    infot = 3;
    Cgeqls(0, 0, -1, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgeqls", infot, nout, lerr, ok);
    infot = 5;
    Cgeqls(2, 1, 0, a, 1, x, b, 2, w, 1, info);
    chkxer("Cgeqls", infot, nout, lerr, ok);
    infot = 8;
    Cgeqls(2, 1, 0, a, 2, x, b, 1, w, 1, info);
    chkxer("Cgeqls", infot, nout, lerr, ok);
    infot = 10;
    Cgeqls(1, 1, 2, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgeqls", infot, nout, lerr, ok);
    //
    //     ZUNGQL
    //
    srnamt = "ZUNGQL";
    infot = 1;
    zungql(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("ZUNGQL", infot, nout, lerr, ok);
    infot = 2;
    zungql(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("ZUNGQL", infot, nout, lerr, ok);
    infot = 2;
    zungql(1, 2, 0, a, 1, x, w, 2, info);
    chkxer("ZUNGQL", infot, nout, lerr, ok);
    infot = 3;
    zungql(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("ZUNGQL", infot, nout, lerr, ok);
    infot = 3;
    zungql(1, 1, 2, a, 1, x, w, 1, info);
    chkxer("ZUNGQL", infot, nout, lerr, ok);
    infot = 5;
    zungql(2, 1, 0, a, 1, x, w, 1, info);
    chkxer("ZUNGQL", infot, nout, lerr, ok);
    infot = 8;
    zungql(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("ZUNGQL", infot, nout, lerr, ok);
    //
    //     ZUNG2L
    //
    srnamt = "ZUNG2L";
    infot = 1;
    zung2l(-1, 0, 0, a, 1, x, w, info);
    chkxer("ZUNG2L", infot, nout, lerr, ok);
    infot = 2;
    zung2l(0, -1, 0, a, 1, x, w, info);
    chkxer("ZUNG2L", infot, nout, lerr, ok);
    infot = 2;
    zung2l(1, 2, 0, a, 1, x, w, info);
    chkxer("ZUNG2L", infot, nout, lerr, ok);
    infot = 3;
    zung2l(0, 0, -1, a, 1, x, w, info);
    chkxer("ZUNG2L", infot, nout, lerr, ok);
    infot = 3;
    zung2l(2, 1, 2, a, 2, x, w, info);
    chkxer("ZUNG2L", infot, nout, lerr, ok);
    infot = 5;
    zung2l(2, 1, 0, a, 1, x, w, info);
    chkxer("ZUNG2L", infot, nout, lerr, ok);
    //
    //     ZUNMQL
    //
    srnamt = "ZUNMQL";
    infot = 1;
    zunmql("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 2;
    zunmql("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 3;
    zunmql("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 4;
    zunmql("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 5;
    zunmql("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 5;
    zunmql("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 5;
    zunmql("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 7;
    zunmql("L", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 7;
    zunmql("R", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 10;
    zunmql("L", "N", 2, 1, 0, a, 2, x, af, 1, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 12;
    zunmql("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    infot = 12;
    zunmql("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("ZUNMQL", infot, nout, lerr, ok);
    //
    //     ZUNM2L
    //
    srnamt = "ZUNM2L";
    infot = 1;
    zunm2l("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2L", infot, nout, lerr, ok);
    infot = 2;
    zunm2l("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2L", infot, nout, lerr, ok);
    infot = 3;
    zunm2l("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2L", infot, nout, lerr, ok);
    infot = 4;
    zunm2l("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2L", infot, nout, lerr, ok);
    infot = 5;
    zunm2l("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2L", infot, nout, lerr, ok);
    infot = 5;
    zunm2l("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2L", infot, nout, lerr, ok);
    infot = 5;
    zunm2l("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2L", infot, nout, lerr, ok);
    infot = 7;
    zunm2l("L", "N", 2, 1, 0, a, 1, x, af, 2, w, info);
    chkxer("ZUNM2L", infot, nout, lerr, ok);
    infot = 7;
    zunm2l("R", "N", 1, 2, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2L", infot, nout, lerr, ok);
    infot = 10;
    zunm2l("L", "N", 2, 1, 0, a, 2, x, af, 1, w, info);
    chkxer("ZUNM2L", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrql
    //
}
