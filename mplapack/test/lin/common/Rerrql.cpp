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

void Rerrql(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_2d<nmax, nmax, REAL> a(fill0);
    arr_2d<nmax, nmax, REAL> af(fill0);
    arr_1d<nmax, REAL> b(fill0);
    arr_1d<nmax, REAL> w(fill0);
    arr_1d<nmax, REAL> x(fill0);
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / (i + j).real();
            af[(i - 1) + (j - 1) * ldaf] = 1.0 / (i + j).real();
        }
        b[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
    }
    ok = true;
    //
    //     Error exits for QL factorization
    //
    //     DGEQLF
    //
    srnamt = "DGEQLF";
    infot = 1;
    INTEGER info = 0;
    dgeqlf(-1, 0, a, 1, b, w, 1, info);
    chkxer("DGEQLF", infot, nout, lerr, ok);
    infot = 2;
    dgeqlf(0, -1, a, 1, b, w, 1, info);
    chkxer("DGEQLF", infot, nout, lerr, ok);
    infot = 4;
    dgeqlf(2, 1, a, 1, b, w, 1, info);
    chkxer("DGEQLF", infot, nout, lerr, ok);
    infot = 7;
    dgeqlf(1, 2, a, 1, b, w, 1, info);
    chkxer("DGEQLF", infot, nout, lerr, ok);
    //
    //     DGEQL2
    //
    srnamt = "DGEQL2";
    infot = 1;
    dgeql2(-1, 0, a, 1, b, w, info);
    chkxer("DGEQL2", infot, nout, lerr, ok);
    infot = 2;
    dgeql2(0, -1, a, 1, b, w, info);
    chkxer("DGEQL2", infot, nout, lerr, ok);
    infot = 4;
    dgeql2(2, 1, a, 1, b, w, info);
    chkxer("DGEQL2", infot, nout, lerr, ok);
    //
    //     Rgeqls
    //
    srnamt = "Rgeqls";
    infot = 1;
    Rgeqls(-1, 0, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgeqls", infot, nout, lerr, ok);
    infot = 2;
    Rgeqls(0, -1, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgeqls", infot, nout, lerr, ok);
    infot = 2;
    Rgeqls(1, 2, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgeqls", infot, nout, lerr, ok);
    infot = 3;
    Rgeqls(0, 0, -1, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgeqls", infot, nout, lerr, ok);
    infot = 5;
    Rgeqls(2, 1, 0, a, 1, x, b, 2, w, 1, info);
    chkxer("Rgeqls", infot, nout, lerr, ok);
    infot = 8;
    Rgeqls(2, 1, 0, a, 2, x, b, 1, w, 1, info);
    chkxer("Rgeqls", infot, nout, lerr, ok);
    infot = 10;
    Rgeqls(1, 1, 2, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgeqls", infot, nout, lerr, ok);
    //
    //     DORGQL
    //
    srnamt = "DORGQL";
    infot = 1;
    dorgql(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("DORGQL", infot, nout, lerr, ok);
    infot = 2;
    dorgql(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("DORGQL", infot, nout, lerr, ok);
    infot = 2;
    dorgql(1, 2, 0, a, 1, x, w, 2, info);
    chkxer("DORGQL", infot, nout, lerr, ok);
    infot = 3;
    dorgql(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("DORGQL", infot, nout, lerr, ok);
    infot = 3;
    dorgql(1, 1, 2, a, 1, x, w, 1, info);
    chkxer("DORGQL", infot, nout, lerr, ok);
    infot = 5;
    dorgql(2, 1, 0, a, 1, x, w, 1, info);
    chkxer("DORGQL", infot, nout, lerr, ok);
    infot = 8;
    dorgql(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("DORGQL", infot, nout, lerr, ok);
    //
    //     DORG2L
    //
    srnamt = "DORG2L";
    infot = 1;
    dorg2l(-1, 0, 0, a, 1, x, w, info);
    chkxer("DORG2L", infot, nout, lerr, ok);
    infot = 2;
    dorg2l(0, -1, 0, a, 1, x, w, info);
    chkxer("DORG2L", infot, nout, lerr, ok);
    infot = 2;
    dorg2l(1, 2, 0, a, 1, x, w, info);
    chkxer("DORG2L", infot, nout, lerr, ok);
    infot = 3;
    dorg2l(0, 0, -1, a, 1, x, w, info);
    chkxer("DORG2L", infot, nout, lerr, ok);
    infot = 3;
    dorg2l(2, 1, 2, a, 2, x, w, info);
    chkxer("DORG2L", infot, nout, lerr, ok);
    infot = 5;
    dorg2l(2, 1, 0, a, 1, x, w, info);
    chkxer("DORG2L", infot, nout, lerr, ok);
    //
    //     DORMQL
    //
    srnamt = "DORMQL";
    infot = 1;
    dormql("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 2;
    dormql("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 3;
    dormql("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 4;
    dormql("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 5;
    dormql("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 5;
    dormql("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 5;
    dormql("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 7;
    dormql("L", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 7;
    dormql("R", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 10;
    dormql("L", "N", 2, 1, 0, a, 2, x, af, 1, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 12;
    dormql("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    infot = 12;
    dormql("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("DORMQL", infot, nout, lerr, ok);
    //
    //     DORM2L
    //
    srnamt = "DORM2L";
    infot = 1;
    dorm2l("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("DORM2L", infot, nout, lerr, ok);
    infot = 2;
    dorm2l("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("DORM2L", infot, nout, lerr, ok);
    infot = 3;
    dorm2l("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("DORM2L", infot, nout, lerr, ok);
    infot = 4;
    dorm2l("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("DORM2L", infot, nout, lerr, ok);
    infot = 5;
    dorm2l("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("DORM2L", infot, nout, lerr, ok);
    infot = 5;
    dorm2l("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("DORM2L", infot, nout, lerr, ok);
    infot = 5;
    dorm2l("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("DORM2L", infot, nout, lerr, ok);
    infot = 7;
    dorm2l("L", "N", 2, 1, 0, a, 1, x, af, 2, w, info);
    chkxer("DORM2L", infot, nout, lerr, ok);
    infot = 7;
    dorm2l("R", "N", 1, 2, 0, a, 1, x, af, 1, w, info);
    chkxer("DORM2L", infot, nout, lerr, ok);
    infot = 10;
    dorm2l("L", "N", 2, 1, 0, a, 2, x, af, 1, w, info);
    chkxer("DORM2L", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrql
    //
}
