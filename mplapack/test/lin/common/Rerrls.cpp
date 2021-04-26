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

void Rerrls(common &cmn, const char *path, INTEGER const nunit) {
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
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
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
    write(nout, star);
    str<2> c2 = path[(2 - 1) + (3 - 1) * ldpath];
    const INTEGER nmax = 2;
    arr_2d<nmax, nmax, REAL> a(fill0);
    a[(1 - 1)] = 1.0;
    a[(2 - 1) * lda] = 2.0e+0;
    a[(2 - 1) + (2 - 1) * lda] = 3.0e+0;
    a[(2 - 1)] = 4.0e+0;
    ok = true;
    //
    arr_2d<nmax, nmax, REAL> b(fill0);
    arr_1d<nmax, REAL> w(fill0);
    INTEGER info = 0;
    arr_1d<nmax, REAL> s(fill0);
    REAL rcond = 0.0;
    INTEGER irnk = 0;
    arr_1d<nmax, int> ip(fill0);
    if (Mlsamen2, c2, "LS") {
        //
        //        Test error exits for the least squares driver routines.
        //
        //        DGELS
        //
        srnamt = "DGELS ";
        infot = 1;
        dgels("/", 0, 0, 0, a, 1, b, 1, w, 1, info);
        chkxer("DGELS ", infot, nout, lerr, ok);
        infot = 2;
        dgels("N", -1, 0, 0, a, 1, b, 1, w, 1, info);
        chkxer("DGELS ", infot, nout, lerr, ok);
        infot = 3;
        dgels("N", 0, -1, 0, a, 1, b, 1, w, 1, info);
        chkxer("DGELS ", infot, nout, lerr, ok);
        infot = 4;
        dgels("N", 0, 0, -1, a, 1, b, 1, w, 1, info);
        chkxer("DGELS ", infot, nout, lerr, ok);
        infot = 6;
        dgels("N", 2, 0, 0, a, 1, b, 2, w, 2, info);
        chkxer("DGELS ", infot, nout, lerr, ok);
        infot = 8;
        dgels("N", 2, 0, 0, a, 2, b, 1, w, 2, info);
        chkxer("DGELS ", infot, nout, lerr, ok);
        infot = 10;
        dgels("N", 1, 1, 0, a, 1, b, 1, w, 1, info);
        chkxer("DGELS ", infot, nout, lerr, ok);
        //
        //        DGELSS
        //
        srnamt = "DGELSS";
        infot = 1;
        dgelss(-1, 0, 0, a, 1, b, 1, s, rcond, irnk, w, 1, info);
        chkxer("DGELSS", infot, nout, lerr, ok);
        infot = 2;
        dgelss(0, -1, 0, a, 1, b, 1, s, rcond, irnk, w, 1, info);
        chkxer("DGELSS", infot, nout, lerr, ok);
        infot = 3;
        dgelss(0, 0, -1, a, 1, b, 1, s, rcond, irnk, w, 1, info);
        chkxer("DGELSS", infot, nout, lerr, ok);
        infot = 5;
        dgelss(2, 0, 0, a, 1, b, 2, s, rcond, irnk, w, 2, info);
        chkxer("DGELSS", infot, nout, lerr, ok);
        infot = 7;
        dgelss(2, 0, 0, a, 2, b, 1, s, rcond, irnk, w, 2, info);
        chkxer("DGELSS", infot, nout, lerr, ok);
        //
        //        DGELSY
        //
        srnamt = "DGELSY";
        infot = 1;
        dgelsy(-1, 0, 0, a, 1, b, 1, ip, rcond, irnk, w, 10, info);
        chkxer("DGELSY", infot, nout, lerr, ok);
        infot = 2;
        dgelsy(0, -1, 0, a, 1, b, 1, ip, rcond, irnk, w, 10, info);
        chkxer("DGELSY", infot, nout, lerr, ok);
        infot = 3;
        dgelsy(0, 0, -1, a, 1, b, 1, ip, rcond, irnk, w, 10, info);
        chkxer("DGELSY", infot, nout, lerr, ok);
        infot = 5;
        dgelsy(2, 0, 0, a, 1, b, 2, ip, rcond, irnk, w, 10, info);
        chkxer("DGELSY", infot, nout, lerr, ok);
        infot = 7;
        dgelsy(2, 0, 0, a, 2, b, 1, ip, rcond, irnk, w, 10, info);
        chkxer("DGELSY", infot, nout, lerr, ok);
        infot = 12;
        dgelsy(2, 2, 1, a, 2, b, 2, ip, rcond, irnk, w, 1, info);
        chkxer("DGELSY", infot, nout, lerr, ok);
        //
        //        DGELSD
        //
        srnamt = "DGELSD";
        infot = 1;
        dgelsd(-1, 0, 0, a, 1, b, 1, s, rcond, irnk, w, 10, ip, info);
        chkxer("DGELSD", infot, nout, lerr, ok);
        infot = 2;
        dgelsd(0, -1, 0, a, 1, b, 1, s, rcond, irnk, w, 10, ip, info);
        chkxer("DGELSD", infot, nout, lerr, ok);
        infot = 3;
        dgelsd(0, 0, -1, a, 1, b, 1, s, rcond, irnk, w, 10, ip, info);
        chkxer("DGELSD", infot, nout, lerr, ok);
        infot = 5;
        dgelsd(2, 0, 0, a, 1, b, 2, s, rcond, irnk, w, 10, ip, info);
        chkxer("DGELSD", infot, nout, lerr, ok);
        infot = 7;
        dgelsd(2, 0, 0, a, 2, b, 1, s, rcond, irnk, w, 10, ip, info);
        chkxer("DGELSD", infot, nout, lerr, ok);
        infot = 12;
        dgelsd(2, 2, 1, a, 2, b, 2, s, rcond, irnk, w, 1, ip, info);
        chkxer("DGELSD", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrls
    //
}
