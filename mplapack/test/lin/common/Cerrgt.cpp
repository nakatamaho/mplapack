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

void Cerrgt(common &cmn, const char *path, INTEGER const nunit) {
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
    INTEGER i = 0;
    const INTEGER nmax = 2;
    arr_1d<nmax, REAL> d(fill0);
    arr_1d<nmax, COMPLEX> e(fill0);
    arr_1d<nmax, COMPLEX> dl(fill0);
    arr_1d<nmax, COMPLEX> du(fill0);
    for (i = 1; i <= nmax; i = i + 1) {
        d[i - 1] = 1.0;
        e[i - 1] = 2.0;
        dl[i - 1] = 3.e0;
        du[i - 1] = 4.e0;
    }
    REAL anorm = 1.0;
    ok = true;
    //
    arr_1d<nmax, COMPLEX> du2(fill0);
    arr_1d<nmax, int> ip(fill0);
    INTEGER info = 0;
    arr_1d<nmax, COMPLEX> x(fill0);
    arr_1d<nmax, COMPLEX> dlf(fill0);
    arr_1d<nmax, COMPLEX> ef(fill0);
    arr_1d<nmax, COMPLEX> duf(fill0);
    arr_1d<nmax, COMPLEX> b(fill0);
    arr_1d<nmax, REAL> r1(fill0);
    arr_1d<nmax, REAL> r2(fill0);
    arr_1d<nmax, COMPLEX> w(fill0);
    arr_1d<nmax, REAL> rw(fill0);
    REAL rcond = 0.0;
    arr_1d<nmax, REAL> df(fill0);
    if (Mlsamen2, c2, "GT") {
        //
        //        Test error exits for the general tridiagonal routines.
        //
        //        ZGTTRF
        //
        srnamt = "ZGTTRF";
        infot = 1;
        zgttrf(-1, dl, e, du, du2, ip, info);
        chkxer("ZGTTRF", infot, nout, lerr, ok);
        //
        //        ZGTTRS
        //
        srnamt = "ZGTTRS";
        infot = 1;
        zgttrs("/", 0, 0, dl, e, du, du2, ip, x, 1, info);
        chkxer("ZGTTRS", infot, nout, lerr, ok);
        infot = 2;
        zgttrs("N", -1, 0, dl, e, du, du2, ip, x, 1, info);
        chkxer("ZGTTRS", infot, nout, lerr, ok);
        infot = 3;
        zgttrs("N", 0, -1, dl, e, du, du2, ip, x, 1, info);
        chkxer("ZGTTRS", infot, nout, lerr, ok);
        infot = 10;
        zgttrs("N", 2, 1, dl, e, du, du2, ip, x, 1, info);
        chkxer("ZGTTRS", infot, nout, lerr, ok);
        //
        //        ZGTRFS
        //
        srnamt = "ZGTRFS";
        infot = 1;
        zgtrfs("/", 0, 0, dl, e, du, dlf, ef, duf, du2, ip, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZGTRFS", infot, nout, lerr, ok);
        infot = 2;
        zgtrfs("N", -1, 0, dl, e, du, dlf, ef, duf, du2, ip, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZGTRFS", infot, nout, lerr, ok);
        infot = 3;
        zgtrfs("N", 0, -1, dl, e, du, dlf, ef, duf, du2, ip, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZGTRFS", infot, nout, lerr, ok);
        infot = 13;
        zgtrfs("N", 2, 1, dl, e, du, dlf, ef, duf, du2, ip, b, 1, x, 2, r1, r2, w, rw, info);
        chkxer("ZGTRFS", infot, nout, lerr, ok);
        infot = 15;
        zgtrfs("N", 2, 1, dl, e, du, dlf, ef, duf, du2, ip, b, 2, x, 1, r1, r2, w, rw, info);
        chkxer("ZGTRFS", infot, nout, lerr, ok);
        //
        //        ZGTCON
        //
        srnamt = "ZGTCON";
        infot = 1;
        zgtcon("/", 0, dl, e, du, du2, ip, anorm, rcond, w, info);
        chkxer("ZGTCON", infot, nout, lerr, ok);
        infot = 2;
        zgtcon("I", -1, dl, e, du, du2, ip, anorm, rcond, w, info);
        chkxer("ZGTCON", infot, nout, lerr, ok);
        infot = 8;
        zgtcon("I", 0, dl, e, du, du2, ip, -anorm, rcond, w, info);
        chkxer("ZGTCON", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PT") {
        //
        //        Test error exits for the positive definite tridiagonal
        //        routines.
        //
        //        ZPTTRF
        //
        srnamt = "ZPTTRF";
        infot = 1;
        zpttrf(-1, d, e, info);
        chkxer("ZPTTRF", infot, nout, lerr, ok);
        //
        //        ZPTTRS
        //
        srnamt = "ZPTTRS";
        infot = 1;
        zpttrs("/", 1, 0, d, e, x, 1, info);
        chkxer("ZPTTRS", infot, nout, lerr, ok);
        infot = 2;
        zpttrs("U", -1, 0, d, e, x, 1, info);
        chkxer("ZPTTRS", infot, nout, lerr, ok);
        infot = 3;
        zpttrs("U", 0, -1, d, e, x, 1, info);
        chkxer("ZPTTRS", infot, nout, lerr, ok);
        infot = 7;
        zpttrs("U", 2, 1, d, e, x, 1, info);
        chkxer("ZPTTRS", infot, nout, lerr, ok);
        //
        //        ZPTRFS
        //
        srnamt = "ZPTRFS";
        infot = 1;
        zptrfs("/", 1, 0, d, e, df, ef, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZPTRFS", infot, nout, lerr, ok);
        infot = 2;
        zptrfs("U", -1, 0, d, e, df, ef, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZPTRFS", infot, nout, lerr, ok);
        infot = 3;
        zptrfs("U", 0, -1, d, e, df, ef, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZPTRFS", infot, nout, lerr, ok);
        infot = 9;
        zptrfs("U", 2, 1, d, e, df, ef, b, 1, x, 2, r1, r2, w, rw, info);
        chkxer("ZPTRFS", infot, nout, lerr, ok);
        infot = 11;
        zptrfs("U", 2, 1, d, e, df, ef, b, 2, x, 1, r1, r2, w, rw, info);
        chkxer("ZPTRFS", infot, nout, lerr, ok);
        //
        //        ZPTCON
        //
        srnamt = "ZPTCON";
        infot = 1;
        zptcon(-1, d, e, anorm, rcond, rw, info);
        chkxer("ZPTCON", infot, nout, lerr, ok);
        infot = 4;
        zptcon(0, d, e, -anorm, rcond, rw, info);
        chkxer("ZPTCON", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrgt
    //
}
