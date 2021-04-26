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

void Rerrgt(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_1d<nmax, REAL> d(fill0);
    d[1 - 1] = 1.0;
    d[2 - 1] = 2.0;
    arr_1d<nmax, REAL> df(fill0);
    df[1 - 1] = 1.0;
    df[2 - 1] = 2.0;
    arr_1d<nmax, REAL> e(fill0);
    e[1 - 1] = 3.e0;
    e[2 - 1] = 4.e0;
    arr_1d<nmax, REAL> ef(fill0);
    ef[1 - 1] = 3.e0;
    ef[2 - 1] = 4.e0;
    REAL anorm = 1.0;
    ok = true;
    //
    arr_1d<nmax, REAL> c(fill0);
    arr_1d<nmax, REAL> f(fill0);
    arr_1d<nmax, int> ip(fill0);
    INTEGER info = 0;
    arr_1d<nmax, REAL> x(fill0);
    arr_1d<nmax, REAL> cf(fill0);
    arr_1d<nmax, REAL> b(fill0);
    arr_1d<nmax, REAL> r1(fill0);
    arr_1d<nmax, REAL> r2(fill0);
    arr_1d<nmax, REAL> w(fill0);
    arr_1d<nmax, int> iw(fill0);
    REAL rcond = 0.0;
    if (Mlsamen2, c2, "GT") {
        //
        //        Test error exits for the general tridiagonal routines.
        //
        //        DGTTRF
        //
        srnamt = "DGTTRF";
        infot = 1;
        dgttrf(-1, c, d, e, f, ip, info);
        chkxer("DGTTRF", infot, nout, lerr, ok);
        //
        //        DGTTRS
        //
        srnamt = "DGTTRS";
        infot = 1;
        dgttrs("/", 0, 0, c, d, e, f, ip, x, 1, info);
        chkxer("DGTTRS", infot, nout, lerr, ok);
        infot = 2;
        dgttrs("N", -1, 0, c, d, e, f, ip, x, 1, info);
        chkxer("DGTTRS", infot, nout, lerr, ok);
        infot = 3;
        dgttrs("N", 0, -1, c, d, e, f, ip, x, 1, info);
        chkxer("DGTTRS", infot, nout, lerr, ok);
        infot = 10;
        dgttrs("N", 2, 1, c, d, e, f, ip, x, 1, info);
        chkxer("DGTTRS", infot, nout, lerr, ok);
        //
        //        DGTRFS
        //
        srnamt = "DGTRFS";
        infot = 1;
        dgtrfs("/", 0, 0, c, d, e, cf, df, ef, f, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DGTRFS", infot, nout, lerr, ok);
        infot = 2;
        dgtrfs("N", -1, 0, c, d, e, cf, df, ef, f, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DGTRFS", infot, nout, lerr, ok);
        infot = 3;
        dgtrfs("N", 0, -1, c, d, e, cf, df, ef, f, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DGTRFS", infot, nout, lerr, ok);
        infot = 13;
        dgtrfs("N", 2, 1, c, d, e, cf, df, ef, f, ip, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("DGTRFS", infot, nout, lerr, ok);
        infot = 15;
        dgtrfs("N", 2, 1, c, d, e, cf, df, ef, f, ip, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("DGTRFS", infot, nout, lerr, ok);
        //
        //        DGTCON
        //
        srnamt = "DGTCON";
        infot = 1;
        dgtcon("/", 0, c, d, e, f, ip, anorm, rcond, w, iw, info);
        chkxer("DGTCON", infot, nout, lerr, ok);
        infot = 2;
        dgtcon("I", -1, c, d, e, f, ip, anorm, rcond, w, iw, info);
        chkxer("DGTCON", infot, nout, lerr, ok);
        infot = 8;
        dgtcon("I", 0, c, d, e, f, ip, -anorm, rcond, w, iw, info);
        chkxer("DGTCON", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PT") {
        //
        //        Test error exits for the positive definite tridiagonal
        //        routines.
        //
        //        DPTTRF
        //
        srnamt = "DPTTRF";
        infot = 1;
        dpttrf(-1, d, e, info);
        chkxer("DPTTRF", infot, nout, lerr, ok);
        //
        //        DPTTRS
        //
        srnamt = "DPTTRS";
        infot = 1;
        dpttrs(-1, 0, d, e, x, 1, info);
        chkxer("DPTTRS", infot, nout, lerr, ok);
        infot = 2;
        dpttrs(0, -1, d, e, x, 1, info);
        chkxer("DPTTRS", infot, nout, lerr, ok);
        infot = 6;
        dpttrs(2, 1, d, e, x, 1, info);
        chkxer("DPTTRS", infot, nout, lerr, ok);
        //
        //        DPTRFS
        //
        srnamt = "DPTRFS";
        infot = 1;
        dptrfs(-1, 0, d, e, df, ef, b, 1, x, 1, r1, r2, w, info);
        chkxer("DPTRFS", infot, nout, lerr, ok);
        infot = 2;
        dptrfs(0, -1, d, e, df, ef, b, 1, x, 1, r1, r2, w, info);
        chkxer("DPTRFS", infot, nout, lerr, ok);
        infot = 8;
        dptrfs(2, 1, d, e, df, ef, b, 1, x, 2, r1, r2, w, info);
        chkxer("DPTRFS", infot, nout, lerr, ok);
        infot = 10;
        dptrfs(2, 1, d, e, df, ef, b, 2, x, 1, r1, r2, w, info);
        chkxer("DPTRFS", infot, nout, lerr, ok);
        //
        //        DPTCON
        //
        srnamt = "DPTCON";
        infot = 1;
        dptcon(-1, d, e, anorm, rcond, w, info);
        chkxer("DPTCON", infot, nout, lerr, ok);
        infot = 4;
        dptcon(0, d, e, -anorm, rcond, w, info);
        chkxer("DPTCON", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrgt
    //
}
