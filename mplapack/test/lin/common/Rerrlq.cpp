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

void Rerrlq(common &cmn, const char *path, INTEGER const nunit) {
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
    //     Error exits for LQ factorization
    //
    //     DGELQF
    //
    srnamt = "DGELQF";
    infot = 1;
    INTEGER info = 0;
    dgelqf(-1, 0, a, 1, b, w, 1, info);
    chkxer("DGELQF", infot, nout, lerr, ok);
    infot = 2;
    dgelqf(0, -1, a, 1, b, w, 1, info);
    chkxer("DGELQF", infot, nout, lerr, ok);
    infot = 4;
    dgelqf(2, 1, a, 1, b, w, 2, info);
    chkxer("DGELQF", infot, nout, lerr, ok);
    infot = 7;
    dgelqf(2, 1, a, 2, b, w, 1, info);
    chkxer("DGELQF", infot, nout, lerr, ok);
    //
    //     DGELQ2
    //
    srnamt = "DGELQ2";
    infot = 1;
    dgelq2(-1, 0, a, 1, b, w, info);
    chkxer("DGELQ2", infot, nout, lerr, ok);
    infot = 2;
    dgelq2(0, -1, a, 1, b, w, info);
    chkxer("DGELQ2", infot, nout, lerr, ok);
    infot = 4;
    dgelq2(2, 1, a, 1, b, w, info);
    chkxer("DGELQ2", infot, nout, lerr, ok);
    //
    //     Rgelqs
    //
    srnamt = "Rgelqs";
    infot = 1;
    Rgelqs(-1, 0, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgelqs", infot, nout, lerr, ok);
    infot = 2;
    Rgelqs(0, -1, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgelqs", infot, nout, lerr, ok);
    infot = 2;
    Rgelqs(2, 1, 0, a, 2, x, b, 1, w, 1, info);
    chkxer("Rgelqs", infot, nout, lerr, ok);
    infot = 3;
    Rgelqs(0, 0, -1, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgelqs", infot, nout, lerr, ok);
    infot = 5;
    Rgelqs(2, 2, 0, a, 1, x, b, 2, w, 1, info);
    chkxer("Rgelqs", infot, nout, lerr, ok);
    infot = 8;
    Rgelqs(1, 2, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgelqs", infot, nout, lerr, ok);
    infot = 10;
    Rgelqs(1, 1, 2, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgelqs", infot, nout, lerr, ok);
    //
    //     DORGLQ
    //
    srnamt = "DORGLQ";
    infot = 1;
    dorglq(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("DORGLQ", infot, nout, lerr, ok);
    infot = 2;
    dorglq(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("DORGLQ", infot, nout, lerr, ok);
    infot = 2;
    dorglq(2, 1, 0, a, 2, x, w, 2, info);
    chkxer("DORGLQ", infot, nout, lerr, ok);
    infot = 3;
    dorglq(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("DORGLQ", infot, nout, lerr, ok);
    infot = 3;
    dorglq(1, 1, 2, a, 1, x, w, 1, info);
    chkxer("DORGLQ", infot, nout, lerr, ok);
    infot = 5;
    dorglq(2, 2, 0, a, 1, x, w, 2, info);
    chkxer("DORGLQ", infot, nout, lerr, ok);
    infot = 8;
    dorglq(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("DORGLQ", infot, nout, lerr, ok);
    //
    //     DORGL2
    //
    srnamt = "DORGL2";
    infot = 1;
    dorgl2(-1, 0, 0, a, 1, x, w, info);
    chkxer("DORGL2", infot, nout, lerr, ok);
    infot = 2;
    dorgl2(0, -1, 0, a, 1, x, w, info);
    chkxer("DORGL2", infot, nout, lerr, ok);
    infot = 2;
    dorgl2(2, 1, 0, a, 2, x, w, info);
    chkxer("DORGL2", infot, nout, lerr, ok);
    infot = 3;
    dorgl2(0, 0, -1, a, 1, x, w, info);
    chkxer("DORGL2", infot, nout, lerr, ok);
    infot = 3;
    dorgl2(1, 1, 2, a, 1, x, w, info);
    chkxer("DORGL2", infot, nout, lerr, ok);
    infot = 5;
    dorgl2(2, 2, 0, a, 1, x, w, info);
    chkxer("DORGL2", infot, nout, lerr, ok);
    //
    //     DORMLQ
    //
    srnamt = "DORMLQ";
    infot = 1;
    dormlq("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 2;
    dormlq("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 3;
    dormlq("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 4;
    dormlq("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 5;
    dormlq("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 5;
    dormlq("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 5;
    dormlq("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 7;
    dormlq("L", "N", 2, 0, 2, a, 1, x, af, 2, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 7;
    dormlq("R", "N", 0, 2, 2, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 10;
    dormlq("L", "N", 2, 1, 0, a, 2, x, af, 1, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 12;
    dormlq("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    infot = 12;
    dormlq("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("DORMLQ", infot, nout, lerr, ok);
    //
    //     DORML2
    //
    srnamt = "DORML2";
    infot = 1;
    dorml2("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("DORML2", infot, nout, lerr, ok);
    infot = 2;
    dorml2("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("DORML2", infot, nout, lerr, ok);
    infot = 3;
    dorml2("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("DORML2", infot, nout, lerr, ok);
    infot = 4;
    dorml2("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("DORML2", infot, nout, lerr, ok);
    infot = 5;
    dorml2("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("DORML2", infot, nout, lerr, ok);
    infot = 5;
    dorml2("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("DORML2", infot, nout, lerr, ok);
    infot = 5;
    dorml2("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("DORML2", infot, nout, lerr, ok);
    infot = 7;
    dorml2("L", "N", 2, 1, 2, a, 1, x, af, 2, w, info);
    chkxer("DORML2", infot, nout, lerr, ok);
    infot = 7;
    dorml2("R", "N", 1, 2, 2, a, 1, x, af, 1, w, info);
    chkxer("DORML2", infot, nout, lerr, ok);
    infot = 10;
    dorml2("L", "N", 2, 1, 0, a, 2, x, af, 1, w, info);
    chkxer("DORML2", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrlq
    //
}
