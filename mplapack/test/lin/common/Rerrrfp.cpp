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

void Rerrrfp(common &cmn, INTEGER const nunit) {
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
    arr_2d<1, 1, REAL> a(fill0);
    a[(1 - 1)] = 1.0;
    arr_2d<1, 1, REAL> b(fill0);
    b[(1 - 1)] = 1.0;
    REAL alpha = 1.0;
    REAL beta = 1.0;
    //
    srnamt = "DPFTRF";
    infot = 1;
    INTEGER info = 0;
    dpftrf("/", "U", 0, a, info);
    chkxer("DPFTRF", infot, nout, lerr, ok);
    infot = 2;
    dpftrf("N", "/", 0, a, info);
    chkxer("DPFTRF", infot, nout, lerr, ok);
    infot = 3;
    dpftrf("N", "U", -1, a, info);
    chkxer("DPFTRF", infot, nout, lerr, ok);
    //
    srnamt = "DPFTRS";
    infot = 1;
    dpftrs("/", "U", 0, 0, a, b, 1, info);
    chkxer("DPFTRS", infot, nout, lerr, ok);
    infot = 2;
    dpftrs("N", "/", 0, 0, a, b, 1, info);
    chkxer("DPFTRS", infot, nout, lerr, ok);
    infot = 3;
    dpftrs("N", "U", -1, 0, a, b, 1, info);
    chkxer("DPFTRS", infot, nout, lerr, ok);
    infot = 4;
    dpftrs("N", "U", 0, -1, a, b, 1, info);
    chkxer("DPFTRS", infot, nout, lerr, ok);
    infot = 7;
    dpftrs("N", "U", 0, 0, a, b, 0, info);
    chkxer("DPFTRS", infot, nout, lerr, ok);
    //
    srnamt = "DPFTRI";
    infot = 1;
    dpftri("/", "U", 0, a, info);
    chkxer("DPFTRI", infot, nout, lerr, ok);
    infot = 2;
    dpftri("N", "/", 0, a, info);
    chkxer("DPFTRI", infot, nout, lerr, ok);
    infot = 3;
    dpftri("N", "U", -1, a, info);
    chkxer("DPFTRI", infot, nout, lerr, ok);
    //
    srnamt = "DTFSM ";
    infot = 1;
    dtfsm("/", "L", "U", "T", "U", 0, 0, alpha, a, b, 1);
    chkxer("DTFSM ", infot, nout, lerr, ok);
    infot = 2;
    dtfsm("N", "/", "U", "T", "U", 0, 0, alpha, a, b, 1);
    chkxer("DTFSM ", infot, nout, lerr, ok);
    infot = 3;
    dtfsm("N", "L", "/", "T", "U", 0, 0, alpha, a, b, 1);
    chkxer("DTFSM ", infot, nout, lerr, ok);
    infot = 4;
    dtfsm("N", "L", "U", "/", "U", 0, 0, alpha, a, b, 1);
    chkxer("DTFSM ", infot, nout, lerr, ok);
    infot = 5;
    dtfsm("N", "L", "U", "T", "/", 0, 0, alpha, a, b, 1);
    chkxer("DTFSM ", infot, nout, lerr, ok);
    infot = 6;
    dtfsm("N", "L", "U", "T", "U", -1, 0, alpha, a, b, 1);
    chkxer("DTFSM ", infot, nout, lerr, ok);
    infot = 7;
    dtfsm("N", "L", "U", "T", "U", 0, -1, alpha, a, b, 1);
    chkxer("DTFSM ", infot, nout, lerr, ok);
    infot = 11;
    dtfsm("N", "L", "U", "T", "U", 0, 0, alpha, a, b, 0);
    chkxer("DTFSM ", infot, nout, lerr, ok);
    //
    srnamt = "DTFTRI";
    infot = 1;
    dtftri("/", "L", "N", 0, a, info);
    chkxer("DTFTRI", infot, nout, lerr, ok);
    infot = 2;
    dtftri("N", "/", "N", 0, a, info);
    chkxer("DTFTRI", infot, nout, lerr, ok);
    infot = 3;
    dtftri("N", "L", "/", 0, a, info);
    chkxer("DTFTRI", infot, nout, lerr, ok);
    infot = 4;
    dtftri("N", "L", "N", -1, a, info);
    chkxer("DTFTRI", infot, nout, lerr, ok);
    //
    srnamt = "DTFTTR";
    infot = 1;
    dtfttr("/", "U", 0, a, b, 1, info);
    chkxer("DTFTTR", infot, nout, lerr, ok);
    infot = 2;
    dtfttr("N", "/", 0, a, b, 1, info);
    chkxer("DTFTTR", infot, nout, lerr, ok);
    infot = 3;
    dtfttr("N", "U", -1, a, b, 1, info);
    chkxer("DTFTTR", infot, nout, lerr, ok);
    infot = 6;
    dtfttr("N", "U", 0, a, b, 0, info);
    chkxer("DTFTTR", infot, nout, lerr, ok);
    //
    srnamt = "DTRTTF";
    infot = 1;
    dtrttf("/", "U", 0, a, 1, b, info);
    chkxer("DTRTTF", infot, nout, lerr, ok);
    infot = 2;
    dtrttf("N", "/", 0, a, 1, b, info);
    chkxer("DTRTTF", infot, nout, lerr, ok);
    infot = 3;
    dtrttf("N", "U", -1, a, 1, b, info);
    chkxer("DTRTTF", infot, nout, lerr, ok);
    infot = 5;
    dtrttf("N", "U", 0, a, 0, b, info);
    chkxer("DTRTTF", infot, nout, lerr, ok);
    //
    srnamt = "DTFTTP";
    infot = 1;
    dtfttp("/", "U", 0, a, b, info);
    chkxer("DTFTTP", infot, nout, lerr, ok);
    infot = 2;
    dtfttp("N", "/", 0, a, b, info);
    chkxer("DTFTTP", infot, nout, lerr, ok);
    infot = 3;
    dtfttp("N", "U", -1, a, b, info);
    chkxer("DTFTTP", infot, nout, lerr, ok);
    //
    srnamt = "DTPTTF";
    infot = 1;
    dtpttf("/", "U", 0, a, b, info);
    chkxer("DTPTTF", infot, nout, lerr, ok);
    infot = 2;
    dtpttf("N", "/", 0, a, b, info);
    chkxer("DTPTTF", infot, nout, lerr, ok);
    infot = 3;
    dtpttf("N", "U", -1, a, b, info);
    chkxer("DTPTTF", infot, nout, lerr, ok);
    //
    srnamt = "DTRTTP";
    infot = 1;
    dtrttp("/", 0, a, 1, b, info);
    chkxer("DTRTTP", infot, nout, lerr, ok);
    infot = 2;
    dtrttp("U", -1, a, 1, b, info);
    chkxer("DTRTTP", infot, nout, lerr, ok);
    infot = 4;
    dtrttp("U", 0, a, 0, b, info);
    chkxer("DTRTTP", infot, nout, lerr, ok);
    //
    srnamt = "DTPTTR";
    infot = 1;
    dtpttr("/", 0, a, b, 1, info);
    chkxer("DTPTTR", infot, nout, lerr, ok);
    infot = 2;
    dtpttr("U", -1, a, b, 1, info);
    chkxer("DTPTTR", infot, nout, lerr, ok);
    infot = 5;
    dtpttr("U", 0, a, b, 0, info);
    chkxer("DTPTTR", infot, nout, lerr, ok);
    //
    srnamt = "DSFRK ";
    infot = 1;
    dsfrk("/", "U", "N", 0, 0, alpha, a, 1, beta, b);
    chkxer("DSFRK ", infot, nout, lerr, ok);
    infot = 2;
    dsfrk("N", "/", "N", 0, 0, alpha, a, 1, beta, b);
    chkxer("DSFRK ", infot, nout, lerr, ok);
    infot = 3;
    dsfrk("N", "U", "/", 0, 0, alpha, a, 1, beta, b);
    chkxer("DSFRK ", infot, nout, lerr, ok);
    infot = 4;
    dsfrk("N", "U", "N", -1, 0, alpha, a, 1, beta, b);
    chkxer("DSFRK ", infot, nout, lerr, ok);
    infot = 5;
    dsfrk("N", "U", "N", 0, -1, alpha, a, 1, beta, b);
    chkxer("DSFRK ", infot, nout, lerr, ok);
    infot = 8;
    dsfrk("N", "U", "N", 0, 0, alpha, a, 0, beta, b);
    chkxer("DSFRK ", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
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
