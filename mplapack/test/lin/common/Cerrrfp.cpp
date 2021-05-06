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

void Cerrrfp(INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    INTEGER infot;
    INTEGER nout;
    bool ok;
    bool lerr;
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
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Executable Statements ..
    //
    nout = nunit;
    ok = true;
    COMPLEX a[1 * 1];
    COMPLEX b[1 * 1];
    a[0] = COMPLEX(1.0, 1.0);
    b[0] = COMPLEX(1.0, 1.0);
    REAL alpha = 1.0;
    COMPLEX calpha = COMPLEX(1.0, 1.0);
    REAL beta = 1.0;
    //
    infot = 1;
    INTEGER info = 0;
    Cpftrf("/", "U", 0, a, info);
    chkxer("Cpftrf", infot, nout, lerr, ok);
    infot = 2;
    Cpftrf("N", "/", 0, a, info);
    chkxer("Cpftrf", infot, nout, lerr, ok);
    infot = 3;
    Cpftrf("N", "U", -1, a, info);
    chkxer("Cpftrf", infot, nout, lerr, ok);
    //
    infot = 1;
    Cpftrs("/", "U", 0, 0, a, b, 1, info);
    chkxer("Cpftrs", infot, nout, lerr, ok);
    infot = 2;
    Cpftrs("N", "/", 0, 0, a, b, 1, info);
    chkxer("Cpftrs", infot, nout, lerr, ok);
    infot = 3;
    Cpftrs("N", "U", -1, 0, a, b, 1, info);
    chkxer("Cpftrs", infot, nout, lerr, ok);
    infot = 4;
    Cpftrs("N", "U", 0, -1, a, b, 1, info);
    chkxer("Cpftrs", infot, nout, lerr, ok);
    infot = 7;
    Cpftrs("N", "U", 0, 0, a, b, 0, info);
    chkxer("Cpftrs", infot, nout, lerr, ok);
    //
    infot = 1;
    Cpftri("/", "U", 0, a, info);
    chkxer("Cpftri", infot, nout, lerr, ok);
    infot = 2;
    Cpftri("N", "/", 0, a, info);
    chkxer("Cpftri", infot, nout, lerr, ok);
    infot = 3;
    Cpftri("N", "U", -1, a, info);
    chkxer("Cpftri", infot, nout, lerr, ok);
    //
    infot = 1;
    Ctfsm("/", "L", "U", "C", "U", 0, 0, calpha, a, b, 1);
    chkxer("Ctfsm ", infot, nout, lerr, ok);
    infot = 2;
    Ctfsm("N", "/", "U", "C", "U", 0, 0, calpha, a, b, 1);
    chkxer("Ctfsm ", infot, nout, lerr, ok);
    infot = 3;
    Ctfsm("N", "L", "/", "C", "U", 0, 0, calpha, a, b, 1);
    chkxer("Ctfsm ", infot, nout, lerr, ok);
    infot = 4;
    Ctfsm("N", "L", "U", "/", "U", 0, 0, calpha, a, b, 1);
    chkxer("Ctfsm ", infot, nout, lerr, ok);
    infot = 5;
    Ctfsm("N", "L", "U", "C", "/", 0, 0, calpha, a, b, 1);
    chkxer("Ctfsm ", infot, nout, lerr, ok);
    infot = 6;
    Ctfsm("N", "L", "U", "C", "U", -1, 0, calpha, a, b, 1);
    chkxer("Ctfsm ", infot, nout, lerr, ok);
    infot = 7;
    Ctfsm("N", "L", "U", "C", "U", 0, -1, calpha, a, b, 1);
    chkxer("Ctfsm ", infot, nout, lerr, ok);
    infot = 11;
    Ctfsm("N", "L", "U", "C", "U", 0, 0, calpha, a, b, 0);
    chkxer("Ctfsm ", infot, nout, lerr, ok);
    //
    infot = 1;
    Ctftri("/", "L", "N", 0, a, info);
    chkxer("Ctftri", infot, nout, lerr, ok);
    infot = 2;
    Ctftri("N", "/", "N", 0, a, info);
    chkxer("Ctftri", infot, nout, lerr, ok);
    infot = 3;
    Ctftri("N", "L", "/", 0, a, info);
    chkxer("Ctftri", infot, nout, lerr, ok);
    infot = 4;
    Ctftri("N", "L", "N", -1, a, info);
    chkxer("Ctftri", infot, nout, lerr, ok);
    //
    infot = 1;
    Ctfttr("/", "U", 0, a, b, 1, info);
    chkxer("Ctfttr", infot, nout, lerr, ok);
    infot = 2;
    Ctfttr("N", "/", 0, a, b, 1, info);
    chkxer("Ctfttr", infot, nout, lerr, ok);
    infot = 3;
    Ctfttr("N", "U", -1, a, b, 1, info);
    chkxer("Ctfttr", infot, nout, lerr, ok);
    infot = 6;
    Ctfttr("N", "U", 0, a, b, 0, info);
    chkxer("Ctfttr", infot, nout, lerr, ok);
    //
    infot = 1;
    Ctrttf("/", "U", 0, a, 1, b, info);
    chkxer("Ctrttf", infot, nout, lerr, ok);
    infot = 2;
    Ctrttf("N", "/", 0, a, 1, b, info);
    chkxer("Ctrttf", infot, nout, lerr, ok);
    infot = 3;
    Ctrttf("N", "U", -1, a, 1, b, info);
    chkxer("Ctrttf", infot, nout, lerr, ok);
    infot = 5;
    Ctrttf("N", "U", 0, a, 0, b, info);
    chkxer("Ctrttf", infot, nout, lerr, ok);
    //
    infot = 1;
    Ctfttp("/", "U", 0, a, b, info);
    chkxer("Ctfttp", infot, nout, lerr, ok);
    infot = 2;
    Ctfttp("N", "/", 0, a, b, info);
    chkxer("Ctfttp", infot, nout, lerr, ok);
    infot = 3;
    Ctfttp("N", "U", -1, a, b, info);
    chkxer("Ctfttp", infot, nout, lerr, ok);
    //
    infot = 1;
    Ctpttf("/", "U", 0, a, b, info);
    chkxer("Ctpttf", infot, nout, lerr, ok);
    infot = 2;
    Ctpttf("N", "/", 0, a, b, info);
    chkxer("Ctpttf", infot, nout, lerr, ok);
    infot = 3;
    Ctpttf("N", "U", -1, a, b, info);
    chkxer("Ctpttf", infot, nout, lerr, ok);
    //
    infot = 1;
    Ctrttp("/", 0, a, 1, b, info);
    chkxer("Ctrttp", infot, nout, lerr, ok);
    infot = 2;
    Ctrttp("U", -1, a, 1, b, info);
    chkxer("Ctrttp", infot, nout, lerr, ok);
    infot = 4;
    Ctrttp("U", 0, a, 0, b, info);
    chkxer("Ctrttp", infot, nout, lerr, ok);
    //
    infot = 1;
    Ctpttr("/", 0, a, b, 1, info);
    chkxer("Ctpttr", infot, nout, lerr, ok);
    infot = 2;
    Ctpttr("U", -1, a, b, 1, info);
    chkxer("Ctpttr", infot, nout, lerr, ok);
    infot = 5;
    Ctpttr("U", 0, a, b, 0, info);
    chkxer("Ctpttr", infot, nout, lerr, ok);
    //
    infot = 1;
    Chfrk("/", "U", "N", 0, 0, alpha, a, 1, beta, b);
    chkxer("Chfrk ", infot, nout, lerr, ok);
    infot = 2;
    Chfrk("N", "/", "N", 0, 0, alpha, a, 1, beta, b);
    chkxer("Chfrk ", infot, nout, lerr, ok);
    infot = 3;
    Chfrk("N", "U", "/", 0, 0, alpha, a, 1, beta, b);
    chkxer("Chfrk ", infot, nout, lerr, ok);
    infot = 4;
    Chfrk("N", "U", "N", -1, 0, alpha, a, 1, beta, b);
    chkxer("Chfrk ", infot, nout, lerr, ok);
    infot = 5;
    Chfrk("N", "U", "N", 0, -1, alpha, a, 1, beta, b);
    chkxer("Chfrk ", infot, nout, lerr, ok);
    infot = 8;
    Chfrk("N", "U", "N", 0, 0, alpha, a, 0, beta, b);
    chkxer("Chfrk ", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,'COMPLEX*16 RFP routines passed the tests of the ','error exits')");
    } else {
        write(nout, "(' *** RFP routines failed the tests of the error ','exits ***')");
    }
    //
    //     End of Cerrrfp
    //
}
