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

void Rerrqr(const char *path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);

    INTEGER infot;
    INTEGER nout;
    bool ok;
    bool lerr;
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
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 2;
    INTEGER i = 0;
    REAL a[nmax * nmax];
    REAL af[nmax * nmax];
    REAL b[nmax];
    REAL w[nmax];
    REAL x[nmax];
    INTEGER lda = nmax;
    INTEGER ldaf = nmax;
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / castREAL(i + j);
            af[(i - 1) + (j - 1) * ldaf] = 1.0 / castREAL(i + j);
        }
        b[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
    }
    ok = true;
    //
    //     Error exits for QR factorization
    //
    //     Rgeqrf
    //
    infot = 1;
    INTEGER info = 0;
    Rgeqrf(-1, 0, a, 1, b, w, 1, info);
    chkxer("Rgeqrf", infot, nout, lerr, ok);
    infot = 2;
    Rgeqrf(0, -1, a, 1, b, w, 1, info);
    chkxer("Rgeqrf", infot, nout, lerr, ok);
    infot = 4;
    Rgeqrf(2, 1, a, 1, b, w, 1, info);
    chkxer("Rgeqrf", infot, nout, lerr, ok);
    infot = 7;
    Rgeqrf(1, 2, a, 1, b, w, 1, info);
    chkxer("Rgeqrf", infot, nout, lerr, ok);
    //
    //     Rgeqrfp
    //
    infot = 1;
    Rgeqrfp(-1, 0, a, 1, b, w, 1, info);
    chkxer("Rgeqrfp", infot, nout, lerr, ok);
    infot = 2;
    Rgeqrfp(0, -1, a, 1, b, w, 1, info);
    chkxer("Rgeqrfp", infot, nout, lerr, ok);
    infot = 4;
    Rgeqrfp(2, 1, a, 1, b, w, 1, info);
    chkxer("Rgeqrfp", infot, nout, lerr, ok);
    infot = 7;
    Rgeqrfp(1, 2, a, 1, b, w, 1, info);
    chkxer("Rgeqrfp", infot, nout, lerr, ok);
    //
    //     Rgeqr2
    //
    infot = 1;
    Rgeqr2(-1, 0, a, 1, b, w, info);
    chkxer("Rgeqr2", infot, nout, lerr, ok);
    infot = 2;
    Rgeqr2(0, -1, a, 1, b, w, info);
    chkxer("Rgeqr2", infot, nout, lerr, ok);
    infot = 4;
    Rgeqr2(2, 1, a, 1, b, w, info);
    chkxer("Rgeqr2", infot, nout, lerr, ok);
    //
    //     Rgeqr2p
    //
    infot = 1;
    Rgeqr2p(-1, 0, a, 1, b, w, info);
    chkxer("Rgeqr2p", infot, nout, lerr, ok);
    infot = 2;
    Rgeqr2p(0, -1, a, 1, b, w, info);
    chkxer("Rgeqr2p", infot, nout, lerr, ok);
    infot = 4;
    Rgeqr2p(2, 1, a, 1, b, w, info);
    chkxer("Rgeqr2p", infot, nout, lerr, ok);
    //
    //     Rgeqrs
    //
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
    //     Rorgqr
    //
    infot = 1;
    Rorgqr(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("Rorgqr", infot, nout, lerr, ok);
    infot = 2;
    Rorgqr(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("Rorgqr", infot, nout, lerr, ok);
    infot = 2;
    Rorgqr(1, 2, 0, a, 1, x, w, 2, info);
    chkxer("Rorgqr", infot, nout, lerr, ok);
    infot = 3;
    Rorgqr(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("Rorgqr", infot, nout, lerr, ok);
    infot = 3;
    Rorgqr(1, 1, 2, a, 1, x, w, 1, info);
    chkxer("Rorgqr", infot, nout, lerr, ok);
    infot = 5;
    Rorgqr(2, 2, 0, a, 1, x, w, 2, info);
    chkxer("Rorgqr", infot, nout, lerr, ok);
    infot = 8;
    Rorgqr(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("Rorgqr", infot, nout, lerr, ok);
    //
    //     Rorg2r
    //
    infot = 1;
    Rorg2r(-1, 0, 0, a, 1, x, w, info);
    chkxer("Rorg2r", infot, nout, lerr, ok);
    infot = 2;
    Rorg2r(0, -1, 0, a, 1, x, w, info);
    chkxer("Rorg2r", infot, nout, lerr, ok);
    infot = 2;
    Rorg2r(1, 2, 0, a, 1, x, w, info);
    chkxer("Rorg2r", infot, nout, lerr, ok);
    infot = 3;
    Rorg2r(0, 0, -1, a, 1, x, w, info);
    chkxer("Rorg2r", infot, nout, lerr, ok);
    infot = 3;
    Rorg2r(2, 1, 2, a, 2, x, w, info);
    chkxer("Rorg2r", infot, nout, lerr, ok);
    infot = 5;
    Rorg2r(2, 1, 0, a, 1, x, w, info);
    chkxer("Rorg2r", infot, nout, lerr, ok);
    //
    //     Rormqr
    //
    infot = 1;
    Rormqr("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 2;
    Rormqr("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 3;
    Rormqr("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 4;
    Rormqr("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 5;
    Rormqr("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 5;
    Rormqr("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 5;
    Rormqr("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 7;
    Rormqr("L", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 7;
    Rormqr("R", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 10;
    Rormqr("L", "N", 2, 1, 0, a, 2, x, af, 1, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 12;
    Rormqr("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    infot = 12;
    Rormqr("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("Rormqr", infot, nout, lerr, ok);
    //
    //     Rorm2r
    //
    infot = 1;
    Rorm2r("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorm2r", infot, nout, lerr, ok);
    infot = 2;
    Rorm2r("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorm2r", infot, nout, lerr, ok);
    infot = 3;
    Rorm2r("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorm2r", infot, nout, lerr, ok);
    infot = 4;
    Rorm2r("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorm2r", infot, nout, lerr, ok);
    infot = 5;
    Rorm2r("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("Rorm2r", infot, nout, lerr, ok);
    infot = 5;
    Rorm2r("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("Rorm2r", infot, nout, lerr, ok);
    infot = 5;
    Rorm2r("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("Rorm2r", infot, nout, lerr, ok);
    infot = 7;
    Rorm2r("L", "N", 2, 1, 0, a, 1, x, af, 2, w, info);
    chkxer("Rorm2r", infot, nout, lerr, ok);
    infot = 7;
    Rorm2r("R", "N", 1, 2, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorm2r", infot, nout, lerr, ok);
    infot = 10;
    Rorm2r("L", "N", 2, 1, 0, a, 2, x, af, 1, w, info);
    chkxer("Rorm2r", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrqr
    //
}
