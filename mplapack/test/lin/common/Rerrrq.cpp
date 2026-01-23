/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine DERRRQ.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

void Rerrrq(fem::str_cref path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    nout = nunit;
    write(nout, star);
    //
    // Set the variables to innocuous values.
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
            a[(i - 1) + (j - 1) * nmax] = 1.0 / castREAL(i + j);
            af[(i - 1) + (j - 1) * nmax] = 1.0 / castREAL(i + j);
        }
        b[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
    }
    ok = true;
    //
    // Error exits for RQ factorization
    //
    // Rgerqf
    //
    srnamt = "Rgerqf";
    infot = 1;
    INTEGER info = 0;
    Rgerqf(-1, 0, a, 1, b, w, 1, info);
    Chkxer("Rgerqf", infot, nout, lerr, ok);
    infot = 2;
    Rgerqf(0, -1, a, 1, b, w, 1, info);
    Chkxer("Rgerqf", infot, nout, lerr, ok);
    infot = 4;
    Rgerqf(2, 1, a, 1, b, w, 2, info);
    Chkxer("Rgerqf", infot, nout, lerr, ok);
    infot = 7;
    Rgerqf(2, 1, a, 2, b, w, 1, info);
    Chkxer("Rgerqf", infot, nout, lerr, ok);
    //
    // Rgerq2
    //
    srnamt = "Rgerq2";
    infot = 1;
    Rgerq2(-1, 0, a, 1, b, w, info);
    Chkxer("Rgerq2", infot, nout, lerr, ok);
    infot = 2;
    Rgerq2(0, -1, a, 1, b, w, info);
    Chkxer("Rgerq2", infot, nout, lerr, ok);
    infot = 4;
    Rgerq2(2, 1, a, 1, b, w, info);
    Chkxer("Rgerq2", infot, nout, lerr, ok);
    //
    // Rgerqs
    //
    srnamt = "Rgerqs";
    infot = 1;
    Rgerqs(-1, 0, 0, a, 1, x, b, 1, w, 1, info);
    Chkxer("Rgerqs", infot, nout, lerr, ok);
    infot = 2;
    Rgerqs(0, -1, 0, a, 1, x, b, 1, w, 1, info);
    Chkxer("Rgerqs", infot, nout, lerr, ok);
    infot = 2;
    Rgerqs(2, 1, 0, a, 2, x, b, 1, w, 1, info);
    Chkxer("Rgerqs", infot, nout, lerr, ok);
    infot = 3;
    Rgerqs(0, 0, -1, a, 1, x, b, 1, w, 1, info);
    Chkxer("Rgerqs", infot, nout, lerr, ok);
    infot = 5;
    Rgerqs(2, 2, 0, a, 1, x, b, 2, w, 1, info);
    Chkxer("Rgerqs", infot, nout, lerr, ok);
    infot = 8;
    Rgerqs(2, 2, 0, a, 2, x, b, 1, w, 1, info);
    Chkxer("Rgerqs", infot, nout, lerr, ok);
    infot = 10;
    Rgerqs(1, 1, 2, a, 1, x, b, 1, w, 1, info);
    Chkxer("Rgerqs", infot, nout, lerr, ok);
    //
    // Rorgrq
    //
    srnamt = "Rorgrq";
    infot = 1;
    Rorgrq(-1, 0, 0, a, 1, x, w, 1, info);
    Chkxer("Rorgrq", infot, nout, lerr, ok);
    infot = 2;
    Rorgrq(0, -1, 0, a, 1, x, w, 1, info);
    Chkxer("Rorgrq", infot, nout, lerr, ok);
    infot = 2;
    Rorgrq(2, 1, 0, a, 2, x, w, 2, info);
    Chkxer("Rorgrq", infot, nout, lerr, ok);
    infot = 3;
    Rorgrq(0, 0, -1, a, 1, x, w, 1, info);
    Chkxer("Rorgrq", infot, nout, lerr, ok);
    infot = 3;
    Rorgrq(1, 2, 2, a, 1, x, w, 1, info);
    Chkxer("Rorgrq", infot, nout, lerr, ok);
    infot = 5;
    Rorgrq(2, 2, 0, a, 1, x, w, 2, info);
    Chkxer("Rorgrq", infot, nout, lerr, ok);
    infot = 8;
    Rorgrq(2, 2, 0, a, 2, x, w, 1, info);
    Chkxer("Rorgrq", infot, nout, lerr, ok);
    //
    // Rorgr2
    //
    srnamt = "Rorgr2";
    infot = 1;
    Rorgr2(-1, 0, 0, a, 1, x, w, info);
    Chkxer("Rorgr2", infot, nout, lerr, ok);
    infot = 2;
    Rorgr2(0, -1, 0, a, 1, x, w, info);
    Chkxer("Rorgr2", infot, nout, lerr, ok);
    infot = 2;
    Rorgr2(2, 1, 0, a, 2, x, w, info);
    Chkxer("Rorgr2", infot, nout, lerr, ok);
    infot = 3;
    Rorgr2(0, 0, -1, a, 1, x, w, info);
    Chkxer("Rorgr2", infot, nout, lerr, ok);
    infot = 3;
    Rorgr2(1, 2, 2, a, 2, x, w, info);
    Chkxer("Rorgr2", infot, nout, lerr, ok);
    infot = 5;
    Rorgr2(2, 2, 0, a, 1, x, w, info);
    Chkxer("Rorgr2", infot, nout, lerr, ok);
    //
    // Rormrq
    //
    srnamt = "Rormrq";
    infot = 1;
    Rormrq("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 2;
    Rormrq("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 3;
    Rormrq("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 4;
    Rormrq("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 5;
    Rormrq("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 5;
    Rormrq("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 5;
    Rormrq("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 7;
    Rormrq("L", "N", 2, 1, 2, a, 1, x, af, 2, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 7;
    Rormrq("R", "N", 1, 2, 2, a, 1, x, af, 1, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 10;
    Rormrq("L", "N", 2, 1, 0, a, 1, x, af, 1, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 12;
    Rormrq("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    infot = 12;
    Rormrq("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    Chkxer("Rormrq", infot, nout, lerr, ok);
    //
    // Rormr2
    //
    srnamt = "Rormr2";
    infot = 1;
    Rormr2("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    Chkxer("Rormr2", infot, nout, lerr, ok);
    infot = 2;
    Rormr2("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    Chkxer("Rormr2", infot, nout, lerr, ok);
    infot = 3;
    Rormr2("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    Chkxer("Rormr2", infot, nout, lerr, ok);
    infot = 4;
    Rormr2("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    Chkxer("Rormr2", infot, nout, lerr, ok);
    infot = 5;
    Rormr2("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    Chkxer("Rormr2", infot, nout, lerr, ok);
    infot = 5;
    Rormr2("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    Chkxer("Rormr2", infot, nout, lerr, ok);
    infot = 5;
    Rormr2("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    Chkxer("Rormr2", infot, nout, lerr, ok);
    infot = 7;
    Rormr2("L", "N", 2, 1, 2, a, 1, x, af, 2, w, info);
    Chkxer("Rormr2", infot, nout, lerr, ok);
    infot = 7;
    Rormr2("R", "N", 1, 2, 2, a, 1, x, af, 1, w, info);
    Chkxer("Rormr2", infot, nout, lerr, ok);
    infot = 10;
    Rormr2("L", "N", 2, 1, 0, a, 1, x, af, 1, w, info);
    Chkxer("Rormr2", infot, nout, lerr, ok);
    //
    // Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    // End of Rerrrq
    //
}
