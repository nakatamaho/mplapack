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
    if (Mlsamen(2, c2, "GT")) {
        //
        //        Test error exits for the general tridiagonal routines.
        //
        //        Cgttrf
        //
        srnamt = "Cgttrf";
        infot = 1;
        Cgttrf(-1, dl, e, du, du2, ip, info);
        chkxer("Cgttrf", infot, nout, lerr, ok);
        //
        //        Cgttrs
        //
        srnamt = "Cgttrs";
        infot = 1;
        Cgttrs("/", 0, 0, dl, e, du, du2, ip, x, 1, info);
        chkxer("Cgttrs", infot, nout, lerr, ok);
        infot = 2;
        Cgttrs("N", -1, 0, dl, e, du, du2, ip, x, 1, info);
        chkxer("Cgttrs", infot, nout, lerr, ok);
        infot = 3;
        Cgttrs("N", 0, -1, dl, e, du, du2, ip, x, 1, info);
        chkxer("Cgttrs", infot, nout, lerr, ok);
        infot = 10;
        Cgttrs("N", 2, 1, dl, e, du, du2, ip, x, 1, info);
        chkxer("Cgttrs", infot, nout, lerr, ok);
        //
        //        Cgtrfs
        //
        srnamt = "Cgtrfs";
        infot = 1;
        Cgtrfs("/", 0, 0, dl, e, du, dlf, ef, duf, du2, ip, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Cgtrfs", infot, nout, lerr, ok);
        infot = 2;
        Cgtrfs("N", -1, 0, dl, e, du, dlf, ef, duf, du2, ip, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Cgtrfs", infot, nout, lerr, ok);
        infot = 3;
        Cgtrfs("N", 0, -1, dl, e, du, dlf, ef, duf, du2, ip, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Cgtrfs", infot, nout, lerr, ok);
        infot = 13;
        Cgtrfs("N", 2, 1, dl, e, du, dlf, ef, duf, du2, ip, b, 1, x, 2, r1, r2, w, rw, info);
        chkxer("Cgtrfs", infot, nout, lerr, ok);
        infot = 15;
        Cgtrfs("N", 2, 1, dl, e, du, dlf, ef, duf, du2, ip, b, 2, x, 1, r1, r2, w, rw, info);
        chkxer("Cgtrfs", infot, nout, lerr, ok);
        //
        //        Cgtcon
        //
        srnamt = "Cgtcon";
        infot = 1;
        Cgtcon("/", 0, dl, e, du, du2, ip, anorm, rcond, w, info);
        chkxer("Cgtcon", infot, nout, lerr, ok);
        infot = 2;
        Cgtcon("I", -1, dl, e, du, du2, ip, anorm, rcond, w, info);
        chkxer("Cgtcon", infot, nout, lerr, ok);
        infot = 8;
        Cgtcon("I", 0, dl, e, du, du2, ip, -anorm, rcond, w, info);
        chkxer("Cgtcon", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PT")) {
        //
        //        Test error exits for the positive definite tridiagonal
        //        routines.
        //
        //        Cpttrf
        //
        srnamt = "Cpttrf";
        infot = 1;
        Cpttrf(-1, d, e, info);
        chkxer("Cpttrf", infot, nout, lerr, ok);
        //
        //        Cpttrs
        //
        srnamt = "Cpttrs";
        infot = 1;
        Cpttrs("/", 1, 0, d, e, x, 1, info);
        chkxer("Cpttrs", infot, nout, lerr, ok);
        infot = 2;
        Cpttrs("U", -1, 0, d, e, x, 1, info);
        chkxer("Cpttrs", infot, nout, lerr, ok);
        infot = 3;
        Cpttrs("U", 0, -1, d, e, x, 1, info);
        chkxer("Cpttrs", infot, nout, lerr, ok);
        infot = 7;
        Cpttrs("U", 2, 1, d, e, x, 1, info);
        chkxer("Cpttrs", infot, nout, lerr, ok);
        //
        //        Cptrfs
        //
        srnamt = "Cptrfs";
        infot = 1;
        Cptrfs("/", 1, 0, d, e, df, ef, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Cptrfs", infot, nout, lerr, ok);
        infot = 2;
        Cptrfs("U", -1, 0, d, e, df, ef, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Cptrfs", infot, nout, lerr, ok);
        infot = 3;
        Cptrfs("U", 0, -1, d, e, df, ef, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Cptrfs", infot, nout, lerr, ok);
        infot = 9;
        Cptrfs("U", 2, 1, d, e, df, ef, b, 1, x, 2, r1, r2, w, rw, info);
        chkxer("Cptrfs", infot, nout, lerr, ok);
        infot = 11;
        Cptrfs("U", 2, 1, d, e, df, ef, b, 2, x, 1, r1, r2, w, rw, info);
        chkxer("Cptrfs", infot, nout, lerr, ok);
        //
        //        Cptcon
        //
        srnamt = "Cptcon";
        infot = 1;
        Cptcon(-1, d, e, anorm, rcond, rw, info);
        chkxer("Cptcon", infot, nout, lerr, ok);
        infot = 4;
        Cptcon(0, d, e, -anorm, rcond, rw, info);
        chkxer("Cptcon", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrgt
    //
}
