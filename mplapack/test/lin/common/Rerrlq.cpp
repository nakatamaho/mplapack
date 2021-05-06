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

void Rerrlq(const char *path, INTEGER const nunit) {
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
    REAL a[nmax * nmax];
    REAL af[nmax * nmax];
    REAL b[nmax];
    REAL w[nmax];
    REAL x[nmax];
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
    //     Error exits for LQ factorization
    //
    //     Rgelqf
    //
    infot = 1;
    INTEGER info = 0;
    Rgelqf(-1, 0, a, 1, b, w, 1, info);
    chkxer("Rgelqf", infot, nout, lerr, ok);
    infot = 2;
    Rgelqf(0, -1, a, 1, b, w, 1, info);
    chkxer("Rgelqf", infot, nout, lerr, ok);
    infot = 4;
    Rgelqf(2, 1, a, 1, b, w, 2, info);
    chkxer("Rgelqf", infot, nout, lerr, ok);
    infot = 7;
    Rgelqf(2, 1, a, 2, b, w, 1, info);
    chkxer("Rgelqf", infot, nout, lerr, ok);
    //
    //     Rgelq2
    //
    infot = 1;
    Rgelq2(-1, 0, a, 1, b, w, info);
    chkxer("Rgelq2", infot, nout, lerr, ok);
    infot = 2;
    Rgelq2(0, -1, a, 1, b, w, info);
    chkxer("Rgelq2", infot, nout, lerr, ok);
    infot = 4;
    Rgelq2(2, 1, a, 1, b, w, info);
    chkxer("Rgelq2", infot, nout, lerr, ok);
    //
    //     Rgelqs
    //
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
    //     Rorglq
    //
    infot = 1;
    Rorglq(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("Rorglq", infot, nout, lerr, ok);
    infot = 2;
    Rorglq(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("Rorglq", infot, nout, lerr, ok);
    infot = 2;
    Rorglq(2, 1, 0, a, 2, x, w, 2, info);
    chkxer("Rorglq", infot, nout, lerr, ok);
    infot = 3;
    Rorglq(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("Rorglq", infot, nout, lerr, ok);
    infot = 3;
    Rorglq(1, 1, 2, a, 1, x, w, 1, info);
    chkxer("Rorglq", infot, nout, lerr, ok);
    infot = 5;
    Rorglq(2, 2, 0, a, 1, x, w, 2, info);
    chkxer("Rorglq", infot, nout, lerr, ok);
    infot = 8;
    Rorglq(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("Rorglq", infot, nout, lerr, ok);
    //
    //     Rorgl2
    //
    infot = 1;
    Rorgl2(-1, 0, 0, a, 1, x, w, info);
    chkxer("Rorgl2", infot, nout, lerr, ok);
    infot = 2;
    Rorgl2(0, -1, 0, a, 1, x, w, info);
    chkxer("Rorgl2", infot, nout, lerr, ok);
    infot = 2;
    Rorgl2(2, 1, 0, a, 2, x, w, info);
    chkxer("Rorgl2", infot, nout, lerr, ok);
    infot = 3;
    Rorgl2(0, 0, -1, a, 1, x, w, info);
    chkxer("Rorgl2", infot, nout, lerr, ok);
    infot = 3;
    Rorgl2(1, 1, 2, a, 1, x, w, info);
    chkxer("Rorgl2", infot, nout, lerr, ok);
    infot = 5;
    Rorgl2(2, 2, 0, a, 1, x, w, info);
    chkxer("Rorgl2", infot, nout, lerr, ok);
    //
    //     Rormlq
    //
    infot = 1;
    Rormlq("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 2;
    Rormlq("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 3;
    Rormlq("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 4;
    Rormlq("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 5;
    Rormlq("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 5;
    Rormlq("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 5;
    Rormlq("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 7;
    Rormlq("L", "N", 2, 0, 2, a, 1, x, af, 2, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 7;
    Rormlq("R", "N", 0, 2, 2, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 10;
    Rormlq("L", "N", 2, 1, 0, a, 2, x, af, 1, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 12;
    Rormlq("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    infot = 12;
    Rormlq("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("Rormlq", infot, nout, lerr, ok);
    //
    //     Rorml2
    //
    infot = 1;
    Rorml2("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorml2", infot, nout, lerr, ok);
    infot = 2;
    Rorml2("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorml2", infot, nout, lerr, ok);
    infot = 3;
    Rorml2("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorml2", infot, nout, lerr, ok);
    infot = 4;
    Rorml2("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorml2", infot, nout, lerr, ok);
    infot = 5;
    Rorml2("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("Rorml2", infot, nout, lerr, ok);
    infot = 5;
    Rorml2("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("Rorml2", infot, nout, lerr, ok);
    infot = 5;
    Rorml2("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("Rorml2", infot, nout, lerr, ok);
    infot = 7;
    Rorml2("L", "N", 2, 1, 2, a, 1, x, af, 2, w, info);
    chkxer("Rorml2", infot, nout, lerr, ok);
    infot = 7;
    Rorml2("R", "N", 1, 2, 2, a, 1, x, af, 1, w, info);
    chkxer("Rorml2", infot, nout, lerr, ok);
    infot = 10;
    Rorml2("L", "N", 2, 1, 0, a, 2, x, af, 1, w, info);
    chkxer("Rorml2", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrlq
    //
}
