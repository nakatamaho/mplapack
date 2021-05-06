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

void Cerrrq(const char *path, INTEGER const nunit) {
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
    COMPLEX a[nmax * nmax];
    COMPLEX af[nmax * nmax];
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
    //     Error exits for RQ factorization
    //
    //     Cgerqf
    //
    infot = 1;
    INTEGER info = 0;
    Cgerqf(-1, 0, a, 1, b, w, 1, info);
    chkxer("Cgerqf", infot, nout, lerr, ok);
    infot = 2;
    Cgerqf(0, -1, a, 1, b, w, 1, info);
    chkxer("Cgerqf", infot, nout, lerr, ok);
    infot = 4;
    Cgerqf(2, 1, a, 1, b, w, 2, info);
    chkxer("Cgerqf", infot, nout, lerr, ok);
    infot = 7;
    Cgerqf(2, 1, a, 2, b, w, 1, info);
    chkxer("Cgerqf", infot, nout, lerr, ok);
    //
    //     Cgerq2
    //
    infot = 1;
    Cgerq2(-1, 0, a, 1, b, w, info);
    chkxer("Cgerq2", infot, nout, lerr, ok);
    infot = 2;
    Cgerq2(0, -1, a, 1, b, w, info);
    chkxer("Cgerq2", infot, nout, lerr, ok);
    infot = 4;
    Cgerq2(2, 1, a, 1, b, w, info);
    chkxer("Cgerq2", infot, nout, lerr, ok);
    //
    //     Cgerqs
    //
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
    //     Cungrq
    //
    infot = 1;
    Cungrq(-1, 0, 0, a, 1, x, w, 1, info);
    chkxer("Cungrq", infot, nout, lerr, ok);
    infot = 2;
    Cungrq(0, -1, 0, a, 1, x, w, 1, info);
    chkxer("Cungrq", infot, nout, lerr, ok);
    infot = 2;
    Cungrq(2, 1, 0, a, 2, x, w, 2, info);
    chkxer("Cungrq", infot, nout, lerr, ok);
    infot = 3;
    Cungrq(0, 0, -1, a, 1, x, w, 1, info);
    chkxer("Cungrq", infot, nout, lerr, ok);
    infot = 3;
    Cungrq(1, 2, 2, a, 1, x, w, 1, info);
    chkxer("Cungrq", infot, nout, lerr, ok);
    infot = 5;
    Cungrq(2, 2, 0, a, 1, x, w, 2, info);
    chkxer("Cungrq", infot, nout, lerr, ok);
    infot = 8;
    Cungrq(2, 2, 0, a, 2, x, w, 1, info);
    chkxer("Cungrq", infot, nout, lerr, ok);
    //
    //     Cungr2
    //
    infot = 1;
    Cungr2(-1, 0, 0, a, 1, x, w, info);
    chkxer("Cungr2", infot, nout, lerr, ok);
    infot = 2;
    Cungr2(0, -1, 0, a, 1, x, w, info);
    chkxer("Cungr2", infot, nout, lerr, ok);
    infot = 2;
    Cungr2(2, 1, 0, a, 2, x, w, info);
    chkxer("Cungr2", infot, nout, lerr, ok);
    infot = 3;
    Cungr2(0, 0, -1, a, 1, x, w, info);
    chkxer("Cungr2", infot, nout, lerr, ok);
    infot = 3;
    Cungr2(1, 2, 2, a, 2, x, w, info);
    chkxer("Cungr2", infot, nout, lerr, ok);
    infot = 5;
    Cungr2(2, 2, 0, a, 1, x, w, info);
    chkxer("Cungr2", infot, nout, lerr, ok);
    //
    //     Cunmrq
    //
    infot = 1;
    Cunmrq("/", "N", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 2;
    Cunmrq("L", "/", 0, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 3;
    Cunmrq("L", "N", -1, 0, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 4;
    Cunmrq("L", "N", 0, -1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 5;
    Cunmrq("L", "N", 0, 0, -1, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 5;
    Cunmrq("L", "N", 0, 1, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 5;
    Cunmrq("R", "N", 1, 0, 1, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 7;
    Cunmrq("L", "N", 2, 1, 2, a, 1, x, af, 2, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 7;
    Cunmrq("R", "N", 1, 2, 2, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 10;
    Cunmrq("L", "N", 2, 1, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 12;
    Cunmrq("L", "N", 1, 2, 0, a, 1, x, af, 1, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    infot = 12;
    Cunmrq("R", "N", 2, 1, 0, a, 1, x, af, 2, w, 1, info);
    chkxer("Cunmrq", infot, nout, lerr, ok);
    //
    //     Cunmr2
    //
    infot = 1;
    Cunmr2("/", "N", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Cunmr2", infot, nout, lerr, ok);
    infot = 2;
    Cunmr2("L", "/", 0, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Cunmr2", infot, nout, lerr, ok);
    infot = 3;
    Cunmr2("L", "N", -1, 0, 0, a, 1, x, af, 1, w, info);
    chkxer("Cunmr2", infot, nout, lerr, ok);
    infot = 4;
    Cunmr2("L", "N", 0, -1, 0, a, 1, x, af, 1, w, info);
    chkxer("Cunmr2", infot, nout, lerr, ok);
    infot = 5;
    Cunmr2("L", "N", 0, 0, -1, a, 1, x, af, 1, w, info);
    chkxer("Cunmr2", infot, nout, lerr, ok);
    infot = 5;
    Cunmr2("L", "N", 0, 1, 1, a, 1, x, af, 1, w, info);
    chkxer("Cunmr2", infot, nout, lerr, ok);
    infot = 5;
    Cunmr2("R", "N", 1, 0, 1, a, 1, x, af, 1, w, info);
    chkxer("Cunmr2", infot, nout, lerr, ok);
    infot = 7;
    Cunmr2("L", "N", 2, 1, 2, a, 1, x, af, 2, w, info);
    chkxer("Cunmr2", infot, nout, lerr, ok);
    infot = 7;
    Cunmr2("R", "N", 1, 2, 2, a, 1, x, af, 1, w, info);
    chkxer("Cunmr2", infot, nout, lerr, ok);
    infot = 10;
    Cunmr2("L", "N", 2, 1, 0, a, 1, x, af, 1, w, info);
    chkxer("Cunmr2", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrrq
    //
}
