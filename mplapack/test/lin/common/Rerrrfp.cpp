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

void Rerrrfp(INTEGER const nunit) {
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
    //     .. Executable Statements ..
    //
    nout = nunit;
    ok = true;
    INTEGER &a[(1 - 1) + (1 - 1) * lda] = 1.0;
    INTEGER &b[(1 - 1) + (1 - 1) * ldb] = 1.0;
    REAL alpha = 1.0;
    REAL beta = 1.0;
    //
    infot = 1;
    REAL a[1 * 1];
    INTEGER info = 0;
    Rpftrf("/", "U", 0, a, info);
    chkxer("Rpftrf", infot, nout, lerr, ok);
    infot = 2;
    Rpftrf("N", "/", 0, a, info);
    chkxer("Rpftrf", infot, nout, lerr, ok);
    infot = 3;
    Rpftrf("N", "U", -1, a, info);
    chkxer("Rpftrf", infot, nout, lerr, ok);
    //
    infot = 1;
    REAL b[1 * 1];
    Rpftrs("/", "U", 0, 0, a, b, 1, info);
    chkxer("Rpftrs", infot, nout, lerr, ok);
    infot = 2;
    Rpftrs("N", "/", 0, 0, a, b, 1, info);
    chkxer("Rpftrs", infot, nout, lerr, ok);
    infot = 3;
    Rpftrs("N", "U", -1, 0, a, b, 1, info);
    chkxer("Rpftrs", infot, nout, lerr, ok);
    infot = 4;
    Rpftrs("N", "U", 0, -1, a, b, 1, info);
    chkxer("Rpftrs", infot, nout, lerr, ok);
    infot = 7;
    Rpftrs("N", "U", 0, 0, a, b, 0, info);
    chkxer("Rpftrs", infot, nout, lerr, ok);
    //
    infot = 1;
    Rpftri("/", "U", 0, a, info);
    chkxer("Rpftri", infot, nout, lerr, ok);
    infot = 2;
    Rpftri("N", "/", 0, a, info);
    chkxer("Rpftri", infot, nout, lerr, ok);
    infot = 3;
    Rpftri("N", "U", -1, a, info);
    chkxer("Rpftri", infot, nout, lerr, ok);
    //
    infot = 1;
    Rtfsm("/", "L", "U", "T", "U", 0, 0, alpha, a, b, 1);
    chkxer("Rtfsm ", infot, nout, lerr, ok);
    infot = 2;
    Rtfsm("N", "/", "U", "T", "U", 0, 0, alpha, a, b, 1);
    chkxer("Rtfsm ", infot, nout, lerr, ok);
    infot = 3;
    Rtfsm("N", "L", "/", "T", "U", 0, 0, alpha, a, b, 1);
    chkxer("Rtfsm ", infot, nout, lerr, ok);
    infot = 4;
    Rtfsm("N", "L", "U", "/", "U", 0, 0, alpha, a, b, 1);
    chkxer("Rtfsm ", infot, nout, lerr, ok);
    infot = 5;
    Rtfsm("N", "L", "U", "T", "/", 0, 0, alpha, a, b, 1);
    chkxer("Rtfsm ", infot, nout, lerr, ok);
    infot = 6;
    Rtfsm("N", "L", "U", "T", "U", -1, 0, alpha, a, b, 1);
    chkxer("Rtfsm ", infot, nout, lerr, ok);
    infot = 7;
    Rtfsm("N", "L", "U", "T", "U", 0, -1, alpha, a, b, 1);
    chkxer("Rtfsm ", infot, nout, lerr, ok);
    infot = 11;
    Rtfsm("N", "L", "U", "T", "U", 0, 0, alpha, a, b, 0);
    chkxer("Rtfsm ", infot, nout, lerr, ok);
    //
    infot = 1;
    Rtftri("/", "L", "N", 0, a, info);
    chkxer("Rtftri", infot, nout, lerr, ok);
    infot = 2;
    Rtftri("N", "/", "N", 0, a, info);
    chkxer("Rtftri", infot, nout, lerr, ok);
    infot = 3;
    Rtftri("N", "L", "/", 0, a, info);
    chkxer("Rtftri", infot, nout, lerr, ok);
    infot = 4;
    Rtftri("N", "L", "N", -1, a, info);
    chkxer("Rtftri", infot, nout, lerr, ok);
    //
    infot = 1;
    Rtfttr("/", "U", 0, a, b, 1, info);
    chkxer("Rtfttr", infot, nout, lerr, ok);
    infot = 2;
    Rtfttr("N", "/", 0, a, b, 1, info);
    chkxer("Rtfttr", infot, nout, lerr, ok);
    infot = 3;
    Rtfttr("N", "U", -1, a, b, 1, info);
    chkxer("Rtfttr", infot, nout, lerr, ok);
    infot = 6;
    Rtfttr("N", "U", 0, a, b, 0, info);
    chkxer("Rtfttr", infot, nout, lerr, ok);
    //
    infot = 1;
    Rtrttf("/", "U", 0, a, 1, b, info);
    chkxer("Rtrttf", infot, nout, lerr, ok);
    infot = 2;
    Rtrttf("N", "/", 0, a, 1, b, info);
    chkxer("Rtrttf", infot, nout, lerr, ok);
    infot = 3;
    Rtrttf("N", "U", -1, a, 1, b, info);
    chkxer("Rtrttf", infot, nout, lerr, ok);
    infot = 5;
    Rtrttf("N", "U", 0, a, 0, b, info);
    chkxer("Rtrttf", infot, nout, lerr, ok);
    //
    infot = 1;
    Rtfttp("/", "U", 0, a, b, info);
    chkxer("Rtfttp", infot, nout, lerr, ok);
    infot = 2;
    Rtfttp("N", "/", 0, a, b, info);
    chkxer("Rtfttp", infot, nout, lerr, ok);
    infot = 3;
    Rtfttp("N", "U", -1, a, b, info);
    chkxer("Rtfttp", infot, nout, lerr, ok);
    //
    infot = 1;
    Rtpttf("/", "U", 0, a, b, info);
    chkxer("Rtpttf", infot, nout, lerr, ok);
    infot = 2;
    Rtpttf("N", "/", 0, a, b, info);
    chkxer("Rtpttf", infot, nout, lerr, ok);
    infot = 3;
    Rtpttf("N", "U", -1, a, b, info);
    chkxer("Rtpttf", infot, nout, lerr, ok);
    //
    infot = 1;
    Rtrttp("/", 0, a, 1, b, info);
    chkxer("Rtrttp", infot, nout, lerr, ok);
    infot = 2;
    Rtrttp("U", -1, a, 1, b, info);
    chkxer("Rtrttp", infot, nout, lerr, ok);
    infot = 4;
    Rtrttp("U", 0, a, 0, b, info);
    chkxer("Rtrttp", infot, nout, lerr, ok);
    //
    infot = 1;
    Rtpttr("/", 0, a, b, 1, info);
    chkxer("Rtpttr", infot, nout, lerr, ok);
    infot = 2;
    Rtpttr("U", -1, a, b, 1, info);
    chkxer("Rtpttr", infot, nout, lerr, ok);
    infot = 5;
    Rtpttr("U", 0, a, b, 0, info);
    chkxer("Rtpttr", infot, nout, lerr, ok);
    //
    infot = 1;
    Rsfrk("/", "U", "N", 0, 0, alpha, a, 1, beta, b);
    chkxer("Rsfrk ", infot, nout, lerr, ok);
    infot = 2;
    Rsfrk("N", "/", "N", 0, 0, alpha, a, 1, beta, b);
    chkxer("Rsfrk ", infot, nout, lerr, ok);
    infot = 3;
    Rsfrk("N", "U", "/", 0, 0, alpha, a, 1, beta, b);
    chkxer("Rsfrk ", infot, nout, lerr, ok);
    infot = 4;
    Rsfrk("N", "U", "N", -1, 0, alpha, a, 1, beta, b);
    chkxer("Rsfrk ", infot, nout, lerr, ok);
    infot = 5;
    Rsfrk("N", "U", "N", 0, -1, alpha, a, 1, beta, b);
    chkxer("Rsfrk ", infot, nout, lerr, ok);
    infot = 8;
    Rsfrk("N", "U", "N", 0, 0, alpha, a, 0, beta, b);
    chkxer("Rsfrk ", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,'DOUBLE PRECISION RFP routines passed the tests of ',"
                    "'the error exits')");
    } else {
        write(nout, "(' *** RFP routines failed the tests of the error ','exits ***')");
    }
    //
    //     End of Rerrrfp
    //
}
