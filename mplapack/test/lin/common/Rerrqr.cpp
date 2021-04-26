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

void Rerrqr(common &cmn, const char *path, INTEGER const nunit) {
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
    //     Error exits for QR factorization
    //
    //     DGEQRF
    //
    srnamt = "DGEQRF";
    infot = 1;
    INTEGER info = 0;
    dgeqrf(-1, 0, a, 1, b, w, 1, info);
    chkxer("DGEQRF", infot, nout, lerr, ok);
    infot = 2;
    dgeqrf(0, -1, a, 1, b, w, 1, info);
    chkxer("DGEQRF", infot, nout, lerr, ok);
    infot = 4;
    dgeqrf(2, 1, a, 1, b, w, 1, info);
    chkxer("DGEQRF", infot, nout, lerr, ok);
    infot = 7;
    dgeqrf(1, 2, a, 1, b, w, 1, info);
    chkxer("DGEQRF", infot, nout, lerr, ok);
    //
    //     DGEQRFP
    //
    srnamt = "DGEQRFP";
    infot = 1;
    dgeqrfp(-1, 0, a, 1, b, w, 1, info);
    chkxer("DGEQRFP", infot, nout, lerr, ok);
    infot = 2;
    dgeqrfp(0, -1, a, 1, b, w, 1, info);
    chkxer("DGEQRFP", infot, nout, lerr, ok);
    infot = 4;
    dgeqrfp(2, 1, a, 1, b, w, 1, info);
    chkxer("DGEQRFP", infot, nout, lerr, ok);
    infot = 7;
    dgeqrfp(1, 2, a, 1, b, w, 1, info);
    chkxer("DGEQRFP", infot, nout, lerr, ok);
    //
    //     DGEQR2
    //
    srnamt = "DGEQR2";
    infot = 1;
    dgeqr2(-1, 0, a, 1, b, w, info);
    chkxer("DGEQR2", infot, nout, lerr, ok);
    infot = 2;
    dgeqr2(0, -1, a, 1, b, w, info);
    chkxer("DGEQR2", infot, nout, lerr, ok);
    infot = 4;
    dgeqr2(2, 1, a, 1, b, w, info);
    chkxer("DGEQR2", infot, nout, lerr, ok);
    //
    //     DGEQR2P
    //
    srnamt = "DGEQR2P";
    infot = 1;
    dgeqr2p(-1, 0, a, 1, b, w, info);
    chkxer("DGEQR2P", infot, nout, lerr, ok);
    infot = 2;
    dgeqr2p(0, -1, a, 1, b, w, info);
    chkxer("DGEQR2P", infot, nout, lerr, ok);
    infot = 4;
    dgeqr2p(2, 1, a, 1, b, w, info);
    chkxer("DGEQR2P", infot, nout, lerr, ok);
    //
    //     Rgeqrs
    //
    srnamt = "Rgeqrs";
    infot = 1;
    Rgeqrs(-1, 0, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgeqrs", infot, nout, lerr, ok);
    infot = 2;
    Rgeqrs(0, -1, 0, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgeqrs", infot, nout, lerr, ok);
    infot = 2;
    Rgeqrs(1, 2, 0, a, 2, x, b, 2, w, 1, info);
    chkxer("Rgeqrs", infot, nout, lerr, ok);
    infot = 3;
    Rgeqrs(0, 0, -1, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgeqrs", infot, nout, lerr, ok);
    infot = 5;
    Rgeqrs(2, 1, 0, a, 1, x, b, 2, w, 1, info);
    chkxer("Rgeqrs", infot, nout, lerr, ok);
    infot = 8;
    Rgeqrs(2, 1, 0, a, 2, x, b, 1, w, 1, info);
    chkxer("Rgeqrs", infot, nout, lerr, ok);
    infot = 10;
    Rgeqrs(1, 1, 2, a, 1, x, b, 1, w, 1, info);
    chkxer("Rgeqrs", infot, nout, lerr, ok);
    //
    //     DORGQR
    //
    srnamt = "DORGQR";
    infot = 1;
    dorgqr(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("DORGQR", infot, nout, lerr, ok);
    infot = 2;
    dorgqr(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("DORGQR", infot, nout, lerr, ok);
    infot = 2;
    dorgqr(1, 2, 0, a, 1, x, w, 2, info);
    chkxer("DORGQR", infot, nout, lerr, ok);
    infot = 3;
    dorgqr(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("DORGQR", infot, nout, lerr, ok);
    infot = 3;
    dorgqr(1, 1, 2, a, 1, x, w, 1, info);
    chkxer("DORGQR", infot, nout, lerr, ok);
    infot = 5;
    dorgqr(2, 2, 0, a, 1, x, w, 2, info);
    chkxer("DORGQR", infot, nout, lerr, ok);
    infot = 8;
    dorgqr(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("DORGQR", infot, nout, lerr, ok);
    //
    //     DORG2R
    //
    srnamt = "DORG2R";
    infot = 1;
    dorg2r(-1, 0, 0, a, 1, x, w, info);
    chkxer("DORG2R", infot, nout, lerr, ok);
    infot = 2;
    dorg2r(0, -1, 0, a, 1, x, w, info);
    chkxer("DORG2R", infot, nout, lerr, ok);
    infot = 2;
    dorg2r(1, 2, 0, a, 1, x, w, info);
    chkxer("DORG2R", infot, nout, lerr, ok);
    infot = 3;
    dorg2r(0, 0, -1, a, 1, x, w, info);
    chkxer("DORG2R", infot, nout, lerr, ok);
    infot = 3;
    dorg2r(2, 1, 2, a, 2, x, w, info);
    chkxer("DORG2R", infot, nout, lerr, ok);
    infot = 5;
    dorg2r(2, 1, 0, a, 1, x, w, info);
    chkxer("DORG2R", infot, nout, lerr, ok);
    //
    //     DORMQR
    //
    srnamt = "DORMQR";
    infot = 1;
    dormqr("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 2;
    dormqr("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 3;
    dormqr("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 4;
    dormqr("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 5;
    dormqr("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 5;
    dormqr("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 5;
    dormqr("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 7;
    dormqr("L", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 7;
    dormqr("R", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 10;
    dormqr("L", "N", 2, 1, 0, a, 2, x, af, 1, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 12;
    dormqr("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    infot = 12;
    dormqr("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("DORMQR", infot, nout, lerr, ok);
    //
    //     DORM2R
    //
    srnamt = "DORM2R";
    infot = 1;
    dorm2r("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("DORM2R", infot, nout, lerr, ok);
    infot = 2;
    dorm2r("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("DORM2R", infot, nout, lerr, ok);
    infot = 3;
    dorm2r("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("DORM2R", infot, nout, lerr, ok);
    infot = 4;
    dorm2r("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("DORM2R", infot, nout, lerr, ok);
    infot = 5;
    dorm2r("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("DORM2R", infot, nout, lerr, ok);
    infot = 5;
    dorm2r("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("DORM2R", infot, nout, lerr, ok);
    infot = 5;
    dorm2r("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("DORM2R", infot, nout, lerr, ok);
    infot = 7;
    dorm2r("L", "N", 2, 1, 0, a, 1, x, af, 2, w, info);
    chkxer("DORM2R", infot, nout, lerr, ok);
    infot = 7;
    dorm2r("R", "N", 1, 2, 0, a, 1, x, af, 1, w, info);
    chkxer("DORM2R", infot, nout, lerr, ok);
    infot = 10;
    dorm2r("L", "N", 2, 1, 0, a, 2, x, af, 1, w, info);
    chkxer("DORM2R", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrqr
    //
}
