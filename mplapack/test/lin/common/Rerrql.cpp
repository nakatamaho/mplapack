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

#include <mplapack_debug.h>

void Rerrql(const char *path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
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
    INTEGER lda = nmax;
    INTEGER ldaf = nmax;
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
    //     Error exits for QL factorization
    //
    //     Rgeqlf
    //
    infot = 1;
    INTEGER info = 0;
    strncpy(srnamt, "Rgeqlf", srnamt_len);
    Rgeqlf(-1, 0, a, 1, b, w, 1, info);
    chkxer("Rgeqlf", infot, nout, lerr, ok);
    infot = 2;
    Rgeqlf(0, -1, a, 1, b, w, 1, info);
    chkxer("Rgeqlf", infot, nout, lerr, ok);
    infot = 4;
    Rgeqlf(2, 1, a, 1, b, w, 1, info);
    chkxer("Rgeqlf", infot, nout, lerr, ok);
    infot = 7;
    Rgeqlf(1, 2, a, 1, b, w, 1, info);
    chkxer("Rgeqlf", infot, nout, lerr, ok);
    //
    //     Rgeql2
    //
    infot = 1;
    strncpy(srnamt, "Rgeql2", srnamt_len);
    Rgeql2(-1, 0, a, 1, b, w, info);
    chkxer("Rgeql2", infot, nout, lerr, ok);
    infot = 2;
    Rgeql2(0, -1, a, 1, b, w, info);
    chkxer("Rgeql2", infot, nout, lerr, ok);
    infot = 4;
    Rgeql2(2, 1, a, 1, b, w, info);
    chkxer("Rgeql2", infot, nout, lerr, ok);
    //
    //     Rgeqls
    //
    infot = 1;
    strncpy(srnamt, "Rgeqls", srnamt_len);
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
    //     Rorgql
    //
    infot = 1;
    strncpy(srnamt, "Rorgql", srnamt_len);
    Rorgql(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("Rorgql", infot, nout, lerr, ok);
    infot = 2;
    Rorgql(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("Rorgql", infot, nout, lerr, ok);
    infot = 2;
    Rorgql(1, 2, 0, a, 1, x, w, 2, info);
    chkxer("Rorgql", infot, nout, lerr, ok);
    infot = 3;
    Rorgql(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("Rorgql", infot, nout, lerr, ok);
    infot = 3;
    Rorgql(1, 1, 2, a, 1, x, w, 1, info);
    chkxer("Rorgql", infot, nout, lerr, ok);
    infot = 5;
    Rorgql(2, 1, 0, a, 1, x, w, 1, info);
    chkxer("Rorgql", infot, nout, lerr, ok);
    infot = 8;
    Rorgql(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("Rorgql", infot, nout, lerr, ok);
    //
    //     Rorg2l
    //
    infot = 1;
    strncpy(srnamt, "Rorg2l", srnamt_len);
    Rorg2l(-1, 0, 0, a, 1, x, w, info);
    chkxer("Rorg2l", infot, nout, lerr, ok);
    infot = 2;
    Rorg2l(0, -1, 0, a, 1, x, w, info);
    chkxer("Rorg2l", infot, nout, lerr, ok);
    infot = 2;
    Rorg2l(1, 2, 0, a, 1, x, w, info);
    chkxer("Rorg2l", infot, nout, lerr, ok);
    infot = 3;
    Rorg2l(0, 0, -1, a, 1, x, w, info);
    chkxer("Rorg2l", infot, nout, lerr, ok);
    infot = 3;
    Rorg2l(2, 1, 2, a, 2, x, w, info);
    chkxer("Rorg2l", infot, nout, lerr, ok);
    infot = 5;
    Rorg2l(2, 1, 0, a, 1, x, w, info);
    chkxer("Rorg2l", infot, nout, lerr, ok);
    //
    //     Rormql
    //
    infot = 1;
    strncpy(srnamt, "Rormql", srnamt_len);
    Rormql("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 2;
    Rormql("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 3;
    Rormql("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 4;
    Rormql("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 5;
    Rormql("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 5;
    Rormql("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 5;
    Rormql("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 7;
    Rormql("L", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 7;
    Rormql("R", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 10;
    Rormql("L", "N", 2, 1, 0, a, 2, x, af, 1, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 12;
    Rormql("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    infot = 12;
    Rormql("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("Rormql", infot, nout, lerr, ok);
    //
    //     Rorm2l
    //
    infot = 1;
    strncpy(srnamt, "Rorm2l", srnamt_len);
    Rorm2l("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorm2l", infot, nout, lerr, ok);
    infot = 2;
    Rorm2l("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorm2l", infot, nout, lerr, ok);
    infot = 3;
    Rorm2l("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorm2l", infot, nout, lerr, ok);
    infot = 4;
    Rorm2l("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorm2l", infot, nout, lerr, ok);
    infot = 5;
    Rorm2l("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("Rorm2l", infot, nout, lerr, ok);
    infot = 5;
    Rorm2l("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("Rorm2l", infot, nout, lerr, ok);
    infot = 5;
    Rorm2l("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("Rorm2l", infot, nout, lerr, ok);
    infot = 7;
    Rorm2l("L", "N", 2, 1, 0, a, 1, x, af, 2, w, info);
    chkxer("Rorm2l", infot, nout, lerr, ok);
    infot = 7;
    Rorm2l("R", "N", 1, 2, 0, a, 1, x, af, 1, w, info);
    chkxer("Rorm2l", infot, nout, lerr, ok);
    infot = 10;
    Rorm2l("L", "N", 2, 1, 0, a, 2, x, af, 1, w, info);
    chkxer("Rorm2l", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrql
    //
}
