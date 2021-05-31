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

void Cerrql(const char *path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    nout = nunit;
    write(nout, star);
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 2;
    INTEGER i = 0;
    COMPLEX a[nmax * nmax];
    COMPLEX af[nmax * nmax];
    INTEGER lda = nmax;
    INTEGER ldaf = nmax;
    COMPLEX b[nmax];
    COMPLEX w[nmax];
    COMPLEX x[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = COMPLEX(1.0 / castREAL(i + j), -1.0 / castREAL(i + j));
            af[(i - 1) + (j - 1) * ldaf] = COMPLEX(1.0 / castREAL(i + j), -1.0 / castREAL(i + j));
        }
        b[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
    }
    ok = true;
    //
    //     Error exits for QL factorization
    //
    //     Cgeqlf
    //
    strncpy(srnamt, "Cgeqlf", srnamt_len);
    infot = 1;
    INTEGER info = 0;
    Cgeqlf(-1, 0, a, 1, b, w, 1, info);
    chkxer("Cgeqlf", infot, nout, lerr, ok);
    infot = 2;
    Cgeqlf(0, -1, a, 1, b, w, 1, info);
    chkxer("Cgeqlf", infot, nout, lerr, ok);
    infot = 4;
    Cgeqlf(2, 1, a, 1, b, w, 1, info);
    chkxer("Cgeqlf", infot, nout, lerr, ok);
    infot = 7;
    Cgeqlf(1, 2, a, 1, b, w, 1, info);
    chkxer("Cgeqlf", infot, nout, lerr, ok);
    //
    //     Cgeql2
    //
    strncpy(srnamt, "Cgeql2", srnamt_len);
    infot = 1;
    Cgeql2(-1, 0, a, 1, b, w, info);
    chkxer("Cgeql2", infot, nout, lerr, ok);
    infot = 2;
    Cgeql2(0, -1, a, 1, b, w, info);
    chkxer("Cgeql2", infot, nout, lerr, ok);
    infot = 4;
    Cgeql2(2, 1, a, 1, b, w, info);
    chkxer("Cgeql2", infot, nout, lerr, ok);
    //
    //     Cgeqls
    //
    strncpy(srnamt, "Cgeqls", srnamt_len);
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
    //     Cungql
    //
    strncpy(srnamt, "Cungql", srnamt_len);
    infot = 1;
    Cungql(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("Cungql", infot, nout, lerr, ok);
    infot = 2;
    Cungql(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("Cungql", infot, nout, lerr, ok);
    infot = 2;
    Cungql(1, 2, 0, a, 1, x, w, 2, info);
    chkxer("Cungql", infot, nout, lerr, ok);
    infot = 3;
    Cungql(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("Cungql", infot, nout, lerr, ok);
    infot = 3;
    Cungql(1, 1, 2, a, 1, x, w, 1, info);
    chkxer("Cungql", infot, nout, lerr, ok);
    infot = 5;
    Cungql(2, 1, 0, a, 1, x, w, 1, info);
    chkxer("Cungql", infot, nout, lerr, ok);
    infot = 8;
    Cungql(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("Cungql", infot, nout, lerr, ok);
    //
    //     Cung2l
    //
    strncpy(srnamt, "Cung2l", srnamt_len);
    infot = 1;
    Cung2l(-1, 0, 0, a, 1, x, w, info);
    chkxer("Cung2l", infot, nout, lerr, ok);
    infot = 2;
    Cung2l(0, -1, 0, a, 1, x, w, info);
    chkxer("Cung2l", infot, nout, lerr, ok);
    infot = 2;
    Cung2l(1, 2, 0, a, 1, x, w, info);
    chkxer("Cung2l", infot, nout, lerr, ok);
    infot = 3;
    Cung2l(0, 0, -1, a, 1, x, w, info);
    chkxer("Cung2l", infot, nout, lerr, ok);
    infot = 3;
    Cung2l(2, 1, 2, a, 2, x, w, info);
    chkxer("Cung2l", infot, nout, lerr, ok);
    infot = 5;
    Cung2l(2, 1, 0, a, 1, x, w, info);
    chkxer("Cung2l", infot, nout, lerr, ok);
    //
    //     Cunmql
    //
    strncpy(srnamt, "Cunmql", srnamt_len);
    infot = 1;
    Cunmql("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 2;
    Cunmql("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 3;
    Cunmql("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 4;
    Cunmql("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 5;
    Cunmql("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 5;
    Cunmql("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 5;
    Cunmql("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 7;
    Cunmql("L", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 7;
    Cunmql("R", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 10;
    Cunmql("L", "N", 2, 1, 0, a, 2, x, af, 1, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 12;
    Cunmql("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    infot = 12;
    Cunmql("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("Cunmql", infot, nout, lerr, ok);
    //
    //     Cunm2l
    //
    strncpy(srnamt, "Cunm2l", srnamt_len);
    infot = 1;
    Cunm2l("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Cunm2l", infot, nout, lerr, ok);
    infot = 2;
    Cunm2l("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Cunm2l", infot, nout, lerr, ok);
    infot = 3;
    Cunm2l("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Cunm2l", infot, nout, lerr, ok);
    infot = 4;
    Cunm2l("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("Cunm2l", infot, nout, lerr, ok);
    infot = 5;
    Cunm2l("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("Cunm2l", infot, nout, lerr, ok);
    infot = 5;
    Cunm2l("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("Cunm2l", infot, nout, lerr, ok);
    infot = 5;
    Cunm2l("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("Cunm2l", infot, nout, lerr, ok);
    infot = 7;
    Cunm2l("L", "N", 2, 1, 0, a, 1, x, af, 2, w, info);
    chkxer("Cunm2l", infot, nout, lerr, ok);
    infot = 7;
    Cunm2l("R", "N", 1, 2, 0, a, 1, x, af, 1, w, info);
    chkxer("Cunm2l", infot, nout, lerr, ok);
    infot = 10;
    Cunm2l("L", "N", 2, 1, 0, a, 2, x, af, 1, w, info);
    chkxer("Cunm2l", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrql
    //
}
