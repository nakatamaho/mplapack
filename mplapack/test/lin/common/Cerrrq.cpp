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

void Cerrrq(common &cmn, const char *path, INTEGER const nunit) {
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
    //     Error exits for RQ factorization
    //
    //     ZGERQF
    //
    srnamt = "ZGERQF";
    infot = 1;
    INTEGER info = 0;
    zgerqf(-1, 0, a, 1, b, w, 1, info);
    chkxer("ZGERQF", infot, nout, lerr, ok);
    infot = 2;
    zgerqf(0, -1, a, 1, b, w, 1, info);
    chkxer("ZGERQF", infot, nout, lerr, ok);
    infot = 4;
    zgerqf(2, 1, a, 1, b, w, 2, info);
    chkxer("ZGERQF", infot, nout, lerr, ok);
    infot = 7;
    zgerqf(2, 1, a, 2, b, w, 1, info);
    chkxer("ZGERQF", infot, nout, lerr, ok);
    //
    //     ZGERQ2
    //
    srnamt = "ZGERQ2";
    infot = 1;
    zgerq2(-1, 0, a, 1, b, w, info);
    chkxer("ZGERQ2", infot, nout, lerr, ok);
    infot = 2;
    zgerq2(0, -1, a, 1, b, w, info);
    chkxer("ZGERQ2", infot, nout, lerr, ok);
    infot = 4;
    zgerq2(2, 1, a, 1, b, w, info);
    chkxer("ZGERQ2", infot, nout, lerr, ok);
    //
    //     Cgerqs
    //
    srnamt = "Cgerqs";
    infot = 1;
    Cgerqs(-1, 0, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgerqs", infot, nout, lerr, ok);
    infot = 2;
    Cgerqs(0, -1, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgerqs", infot, nout, lerr, ok);
    infot = 2;
    Cgerqs(2, 1, 0, a, 2, x, b, 1, w, 1, info);
    chkxer("Cgerqs", infot, nout, lerr, ok);
    infot = 3;
    Cgerqs(0, 0, -1, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgerqs", infot, nout, lerr, ok);
    infot = 5;
    Cgerqs(2, 2, 0, a, 1, x, b, 2, w, 1, info);
    chkxer("Cgerqs", infot, nout, lerr, ok);
    infot = 8;
    Cgerqs(2, 2, 0, a, 2, x, b, 1, w, 1, info);
    chkxer("Cgerqs", infot, nout, lerr, ok);
    infot = 10;
    Cgerqs(1, 1, 2, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgerqs", infot, nout, lerr, ok);
    //
    //     ZUNGRQ
    //
    srnamt = "ZUNGRQ";
    infot = 1;
    zungrq(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("ZUNGRQ", infot, nout, lerr, ok);
    infot = 2;
    zungrq(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("ZUNGRQ", infot, nout, lerr, ok);
    infot = 2;
    zungrq(2, 1, 0, a, 2, x, w, 2, info);
    chkxer("ZUNGRQ", infot, nout, lerr, ok);
    infot = 3;
    zungrq(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("ZUNGRQ", infot, nout, lerr, ok);
    infot = 3;
    zungrq(1, 2, 2, a, 1, x, w, 1, info);
    chkxer("ZUNGRQ", infot, nout, lerr, ok);
    infot = 5;
    zungrq(2, 2, 0, a, 1, x, w, 2, info);
    chkxer("ZUNGRQ", infot, nout, lerr, ok);
    infot = 8;
    zungrq(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("ZUNGRQ", infot, nout, lerr, ok);
    //
    //     ZUNGR2
    //
    srnamt = "ZUNGR2";
    infot = 1;
    zungr2(-1, 0, 0, a, 1, x, w, info);
    chkxer("ZUNGR2", infot, nout, lerr, ok);
    infot = 2;
    zungr2(0, -1, 0, a, 1, x, w, info);
    chkxer("ZUNGR2", infot, nout, lerr, ok);
    infot = 2;
    zungr2(2, 1, 0, a, 2, x, w, info);
    chkxer("ZUNGR2", infot, nout, lerr, ok);
    infot = 3;
    zungr2(0, 0, -1, a, 1, x, w, info);
    chkxer("ZUNGR2", infot, nout, lerr, ok);
    infot = 3;
    zungr2(1, 2, 2, a, 2, x, w, info);
    chkxer("ZUNGR2", infot, nout, lerr, ok);
    infot = 5;
    zungr2(2, 2, 0, a, 1, x, w, info);
    chkxer("ZUNGR2", infot, nout, lerr, ok);
    //
    //     ZUNMRQ
    //
    srnamt = "ZUNMRQ";
    infot = 1;
    zunmrq("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 2;
    zunmrq("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 3;
    zunmrq("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 4;
    zunmrq("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 5;
    zunmrq("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 5;
    zunmrq("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 5;
    zunmrq("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 7;
    zunmrq("L", "N", 2, 1, 2, a, 1, x, af, 2, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 7;
    zunmrq("R", "N", 1, 2, 2, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 10;
    zunmrq("L", "N", 2, 1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 12;
    zunmrq("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    infot = 12;
    zunmrq("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("ZUNMRQ", infot, nout, lerr, ok);
    //
    //     ZUNMR2
    //
    srnamt = "ZUNMR2";
    infot = 1;
    zunmr2("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNMR2", infot, nout, lerr, ok);
    infot = 2;
    zunmr2("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNMR2", infot, nout, lerr, ok);
    infot = 3;
    zunmr2("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNMR2", infot, nout, lerr, ok);
    infot = 4;
    zunmr2("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNMR2", infot, nout, lerr, ok);
    infot = 5;
    zunmr2("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("ZUNMR2", infot, nout, lerr, ok);
    infot = 5;
    zunmr2("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("ZUNMR2", infot, nout, lerr, ok);
    infot = 5;
    zunmr2("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("ZUNMR2", infot, nout, lerr, ok);
    infot = 7;
    zunmr2("L", "N", 2, 1, 2, a, 1, x, af, 2, w, info);
    chkxer("ZUNMR2", infot, nout, lerr, ok);
    infot = 7;
    zunmr2("R", "N", 1, 2, 2, a, 1, x, af, 1, w, info);
    chkxer("ZUNMR2", infot, nout, lerr, ok);
    infot = 10;
    zunmr2("L", "N", 2, 1, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNMR2", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrrq
    //
}
