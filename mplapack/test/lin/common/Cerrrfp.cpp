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

void Cerrrfp(common &cmn, INTEGER const nunit) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Executable Statements ..
    //
    nout = nunit;
    ok = true;
    arr_2d<1, 1, COMPLEX> a(fill0);
    a[(1 - 1)] = COMPLEX(1.0, 1.0);
    arr_2d<1, 1, COMPLEX> b(fill0);
    b[(1 - 1)] = COMPLEX(1.0, 1.0);
    REAL alpha = 1.0;
    COMPLEX calpha = COMPLEX(1.0, 1.0);
    REAL beta = 1.0;
    //
    srnamt = "ZPFTRF";
    infot = 1;
    INTEGER info = 0;
    zpftrf("/", "U", 0, a, info);
    chkxer("ZPFTRF", infot, nout, lerr, ok);
    infot = 2;
    zpftrf("N", "/", 0, a, info);
    chkxer("ZPFTRF", infot, nout, lerr, ok);
    infot = 3;
    zpftrf("N", "U", -1, a, info);
    chkxer("ZPFTRF", infot, nout, lerr, ok);
    //
    srnamt = "ZPFTRS";
    infot = 1;
    zpftrs("/", "U", 0, 0, a, b, 1, info);
    chkxer("ZPFTRS", infot, nout, lerr, ok);
    infot = 2;
    zpftrs("N", "/", 0, 0, a, b, 1, info);
    chkxer("ZPFTRS", infot, nout, lerr, ok);
    infot = 3;
    zpftrs("N", "U", -1, 0, a, b, 1, info);
    chkxer("ZPFTRS", infot, nout, lerr, ok);
    infot = 4;
    zpftrs("N", "U", 0, -1, a, b, 1, info);
    chkxer("ZPFTRS", infot, nout, lerr, ok);
    infot = 7;
    zpftrs("N", "U", 0, 0, a, b, 0, info);
    chkxer("ZPFTRS", infot, nout, lerr, ok);
    //
    srnamt = "ZPFTRI";
    infot = 1;
    zpftri("/", "U", 0, a, info);
    chkxer("ZPFTRI", infot, nout, lerr, ok);
    infot = 2;
    zpftri("N", "/", 0, a, info);
    chkxer("ZPFTRI", infot, nout, lerr, ok);
    infot = 3;
    zpftri("N", "U", -1, a, info);
    chkxer("ZPFTRI", infot, nout, lerr, ok);
    //
    srnamt = "ZTFSM ";
    infot = 1;
    ztfsm("/", "L", "U", "C", "U", 0, 0, calpha, a, b, 1);
    chkxer("ZTFSM ", infot, nout, lerr, ok);
    infot = 2;
    ztfsm("N", "/", "U", "C", "U", 0, 0, calpha, a, b, 1);
    chkxer("ZTFSM ", infot, nout, lerr, ok);
    infot = 3;
    ztfsm("N", "L", "/", "C", "U", 0, 0, calpha, a, b, 1);
    chkxer("ZTFSM ", infot, nout, lerr, ok);
    infot = 4;
    ztfsm("N", "L", "U", "/", "U", 0, 0, calpha, a, b, 1);
    chkxer("ZTFSM ", infot, nout, lerr, ok);
    infot = 5;
    ztfsm("N", "L", "U", "C", "/", 0, 0, calpha, a, b, 1);
    chkxer("ZTFSM ", infot, nout, lerr, ok);
    infot = 6;
    ztfsm("N", "L", "U", "C", "U", -1, 0, calpha, a, b, 1);
    chkxer("ZTFSM ", infot, nout, lerr, ok);
    infot = 7;
    ztfsm("N", "L", "U", "C", "U", 0, -1, calpha, a, b, 1);
    chkxer("ZTFSM ", infot, nout, lerr, ok);
    infot = 11;
    ztfsm("N", "L", "U", "C", "U", 0, 0, calpha, a, b, 0);
    chkxer("ZTFSM ", infot, nout, lerr, ok);
    //
    srnamt = "ZTFTRI";
    infot = 1;
    ztftri("/", "L", "N", 0, a, info);
    chkxer("ZTFTRI", infot, nout, lerr, ok);
    infot = 2;
    ztftri("N", "/", "N", 0, a, info);
    chkxer("ZTFTRI", infot, nout, lerr, ok);
    infot = 3;
    ztftri("N", "L", "/", 0, a, info);
    chkxer("ZTFTRI", infot, nout, lerr, ok);
    infot = 4;
    ztftri("N", "L", "N", -1, a, info);
    chkxer("ZTFTRI", infot, nout, lerr, ok);
    //
    srnamt = "ZTFTTR";
    infot = 1;
    ztfttr("/", "U", 0, a, b, 1, info);
    chkxer("ZTFTTR", infot, nout, lerr, ok);
    infot = 2;
    ztfttr("N", "/", 0, a, b, 1, info);
    chkxer("ZTFTTR", infot, nout, lerr, ok);
    infot = 3;
    ztfttr("N", "U", -1, a, b, 1, info);
    chkxer("ZTFTTR", infot, nout, lerr, ok);
    infot = 6;
    ztfttr("N", "U", 0, a, b, 0, info);
    chkxer("ZTFTTR", infot, nout, lerr, ok);
    //
    srnamt = "ZTRTTF";
    infot = 1;
    ztrttf("/", "U", 0, a, 1, b, info);
    chkxer("ZTRTTF", infot, nout, lerr, ok);
    infot = 2;
    ztrttf("N", "/", 0, a, 1, b, info);
    chkxer("ZTRTTF", infot, nout, lerr, ok);
    infot = 3;
    ztrttf("N", "U", -1, a, 1, b, info);
    chkxer("ZTRTTF", infot, nout, lerr, ok);
    infot = 5;
    ztrttf("N", "U", 0, a, 0, b, info);
    chkxer("ZTRTTF", infot, nout, lerr, ok);
    //
    srnamt = "ZTFTTP";
    infot = 1;
    ztfttp("/", "U", 0, a, b, info);
    chkxer("ZTFTTP", infot, nout, lerr, ok);
    infot = 2;
    ztfttp("N", "/", 0, a, b, info);
    chkxer("ZTFTTP", infot, nout, lerr, ok);
    infot = 3;
    ztfttp("N", "U", -1, a, b, info);
    chkxer("ZTFTTP", infot, nout, lerr, ok);
    //
    srnamt = "ZTPTTF";
    infot = 1;
    ztpttf("/", "U", 0, a, b, info);
    chkxer("ZTPTTF", infot, nout, lerr, ok);
    infot = 2;
    ztpttf("N", "/", 0, a, b, info);
    chkxer("ZTPTTF", infot, nout, lerr, ok);
    infot = 3;
    ztpttf("N", "U", -1, a, b, info);
    chkxer("ZTPTTF", infot, nout, lerr, ok);
    //
    srnamt = "ZTRTTP";
    infot = 1;
    ztrttp("/", 0, a, 1, b, info);
    chkxer("ZTRTTP", infot, nout, lerr, ok);
    infot = 2;
    ztrttp("U", -1, a, 1, b, info);
    chkxer("ZTRTTP", infot, nout, lerr, ok);
    infot = 4;
    ztrttp("U", 0, a, 0, b, info);
    chkxer("ZTRTTP", infot, nout, lerr, ok);
    //
    srnamt = "ZTPTTR";
    infot = 1;
    ztpttr("/", 0, a, b, 1, info);
    chkxer("ZTPTTR", infot, nout, lerr, ok);
    infot = 2;
    ztpttr("U", -1, a, b, 1, info);
    chkxer("ZTPTTR", infot, nout, lerr, ok);
    infot = 5;
    ztpttr("U", 0, a, b, 0, info);
    chkxer("ZTPTTR", infot, nout, lerr, ok);
    //
    srnamt = "ZHFRK ";
    infot = 1;
    zhfrk("/", "U", "N", 0, 0, alpha, a, 1, beta, b);
    chkxer("ZHFRK ", infot, nout, lerr, ok);
    infot = 2;
    zhfrk("N", "/", "N", 0, 0, alpha, a, 1, beta, b);
    chkxer("ZHFRK ", infot, nout, lerr, ok);
    infot = 3;
    zhfrk("N", "U", "/", 0, 0, alpha, a, 1, beta, b);
    chkxer("ZHFRK ", infot, nout, lerr, ok);
    infot = 4;
    zhfrk("N", "U", "N", -1, 0, alpha, a, 1, beta, b);
    chkxer("ZHFRK ", infot, nout, lerr, ok);
    infot = 5;
    zhfrk("N", "U", "N", 0, -1, alpha, a, 1, beta, b);
    chkxer("ZHFRK ", infot, nout, lerr, ok);
    infot = 8;
    zhfrk("N", "U", "N", 0, 0, alpha, a, 0, beta, b);
    chkxer("ZHFRK ", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
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
