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

void Cerrqr(common &cmn, const char *path, INTEGER const nunit) {
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
    //  -- LAPACK test routine (--
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
    //     Error exits for QR factorization
    //
    //     ZGEQRF
    //
    srnamt = "ZGEQRF";
    infot = 1;
    INTEGER info = 0;
    zgeqrf(-1, 0, a, 1, b, w, 1, info);
    chkxer("ZGEQRF", infot, nout, lerr, ok);
    infot = 2;
    zgeqrf(0, -1, a, 1, b, w, 1, info);
    chkxer("ZGEQRF", infot, nout, lerr, ok);
    infot = 4;
    zgeqrf(2, 1, a, 1, b, w, 1, info);
    chkxer("ZGEQRF", infot, nout, lerr, ok);
    infot = 7;
    zgeqrf(1, 2, a, 1, b, w, 1, info);
    chkxer("ZGEQRF", infot, nout, lerr, ok);
    //
    //     ZGEQRFP
    //
    srnamt = "ZGEQRFP";
    infot = 1;
    zgeqrfp(-1, 0, a, 1, b, w, 1, info);
    chkxer("ZGEQRFP", infot, nout, lerr, ok);
    infot = 2;
    zgeqrfp(0, -1, a, 1, b, w, 1, info);
    chkxer("ZGEQRFP", infot, nout, lerr, ok);
    infot = 4;
    zgeqrfp(2, 1, a, 1, b, w, 1, info);
    chkxer("ZGEQRFP", infot, nout, lerr, ok);
    infot = 7;
    zgeqrfp(1, 2, a, 1, b, w, 1, info);
    chkxer("ZGEQRFP", infot, nout, lerr, ok);
    //
    //     ZGEQR2
    //
    srnamt = "ZGEQR2";
    infot = 1;
    zgeqr2(-1, 0, a, 1, b, w, info);
    chkxer("ZGEQR2", infot, nout, lerr, ok);
    infot = 2;
    zgeqr2(0, -1, a, 1, b, w, info);
    chkxer("ZGEQR2", infot, nout, lerr, ok);
    infot = 4;
    zgeqr2(2, 1, a, 1, b, w, info);
    chkxer("ZGEQR2", infot, nout, lerr, ok);
    //
    //     ZGEQR2P
    //
    srnamt = "ZGEQR2P";
    infot = 1;
    zgeqr2p(-1, 0, a, 1, b, w, info);
    chkxer("ZGEQR2P", infot, nout, lerr, ok);
    infot = 2;
    zgeqr2p(0, -1, a, 1, b, w, info);
    chkxer("ZGEQR2P", infot, nout, lerr, ok);
    infot = 4;
    zgeqr2p(2, 1, a, 1, b, w, info);
    chkxer("ZGEQR2P", infot, nout, lerr, ok);
    //
    //     Cgeqrs
    //
    srnamt = "Cgeqrs";
    infot = 1;
    Cgeqrs(-1, 0, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgeqrs", infot, nout, lerr, ok);
    infot = 2;
    Cgeqrs(0, -1, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgeqrs", infot, nout, lerr, ok);
    infot = 2;
    Cgeqrs(1, 2, 0, a, 2, x, b, 2, w, 1, info);
    chkxer("Cgeqrs", infot, nout, lerr, ok);
    infot = 3;
    Cgeqrs(0, 0, -1, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgeqrs", infot, nout, lerr, ok);
    infot = 5;
    Cgeqrs(2, 1, 0, a, 1, x, b, 2, w, 1, info);
    chkxer("Cgeqrs", infot, nout, lerr, ok);
    infot = 8;
    Cgeqrs(2, 1, 0, a, 2, x, b, 1, w, 1, info);
    chkxer("Cgeqrs", infot, nout, lerr, ok);
    infot = 10;
    Cgeqrs(1, 1, 2, a, 1, x, b, 1, w, 1, info);
    chkxer("Cgeqrs", infot, nout, lerr, ok);
    //
    //     ZUNGQR
    //
    srnamt = "ZUNGQR";
    infot = 1;
    zungqr(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("ZUNGQR", infot, nout, lerr, ok);
    infot = 2;
    zungqr(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("ZUNGQR", infot, nout, lerr, ok);
    infot = 2;
    zungqr(1, 2, 0, a, 1, x, w, 2, info);
    chkxer("ZUNGQR", infot, nout, lerr, ok);
    infot = 3;
    zungqr(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("ZUNGQR", infot, nout, lerr, ok);
    infot = 3;
    zungqr(1, 1, 2, a, 1, x, w, 1, info);
    chkxer("ZUNGQR", infot, nout, lerr, ok);
    infot = 5;
    zungqr(2, 2, 0, a, 1, x, w, 2, info);
    chkxer("ZUNGQR", infot, nout, lerr, ok);
    infot = 8;
    zungqr(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("ZUNGQR", infot, nout, lerr, ok);
    //
    //     ZUNG2R
    //
    srnamt = "ZUNG2R";
    infot = 1;
    zung2r(-1, 0, 0, a, 1, x, w, info);
    chkxer("ZUNG2R", infot, nout, lerr, ok);
    infot = 2;
    zung2r(0, -1, 0, a, 1, x, w, info);
    chkxer("ZUNG2R", infot, nout, lerr, ok);
    infot = 2;
    zung2r(1, 2, 0, a, 1, x, w, info);
    chkxer("ZUNG2R", infot, nout, lerr, ok);
    infot = 3;
    zung2r(0, 0, -1, a, 1, x, w, info);
    chkxer("ZUNG2R", infot, nout, lerr, ok);
    infot = 3;
    zung2r(2, 1, 2, a, 2, x, w, info);
    chkxer("ZUNG2R", infot, nout, lerr, ok);
    infot = 5;
    zung2r(2, 1, 0, a, 1, x, w, info);
    chkxer("ZUNG2R", infot, nout, lerr, ok);
    //
    //     ZUNMQR
    //
    srnamt = "ZUNMQR";
    infot = 1;
    zunmqr("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 2;
    zunmqr("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 3;
    zunmqr("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 4;
    zunmqr("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 5;
    zunmqr("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 5;
    zunmqr("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 5;
    zunmqr("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 7;
    zunmqr("L", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 7;
    zunmqr("R", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 10;
    zunmqr("L", "N", 2, 1, 0, a, 2, x, af, 1, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 12;
    zunmqr("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    infot = 12;
    zunmqr("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("ZUNMQR", infot, nout, lerr, ok);
    //
    //     ZUNM2R
    //
    srnamt = "ZUNM2R";
    infot = 1;
    zunm2r("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2R", infot, nout, lerr, ok);
    infot = 2;
    zunm2r("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2R", infot, nout, lerr, ok);
    infot = 3;
    zunm2r("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2R", infot, nout, lerr, ok);
    infot = 4;
    zunm2r("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2R", infot, nout, lerr, ok);
    infot = 5;
    zunm2r("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2R", infot, nout, lerr, ok);
    infot = 5;
    zunm2r("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2R", infot, nout, lerr, ok);
    infot = 5;
    zunm2r("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2R", infot, nout, lerr, ok);
    infot = 7;
    zunm2r("L", "N", 2, 1, 0, a, 1, x, af, 2, w, info);
    chkxer("ZUNM2R", infot, nout, lerr, ok);
    infot = 7;
    zunm2r("R", "N", 1, 2, 0, a, 1, x, af, 1, w, info);
    chkxer("ZUNM2R", infot, nout, lerr, ok);
    infot = 10;
    zunm2r("L", "N", 2, 1, 0, a, 2, x, af, 1, w, info);
    chkxer("ZUNM2R", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrqr
    //
}
