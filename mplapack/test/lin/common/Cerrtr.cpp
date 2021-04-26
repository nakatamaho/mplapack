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

void Cerrtr(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_2d<nmax, nmax, COMPLEX> a(fill0);
    a[(1 - 1)] = 1.0;
    a[(2 - 1) * lda] = 2.0;
    a[(2 - 1) + (2 - 1) * lda] = 3.e0;
    a[(2 - 1)] = 4.e0;
    ok = true;
    //
    //     Test error exits for the general triangular routines.
    //
    INTEGER info = 0;
    arr_1d<nmax, COMPLEX> x(fill0);
    arr_1d<nmax, COMPLEX> b(fill0);
    arr_1d<nmax, REAL> r1(fill0);
    arr_1d<nmax, REAL> r2(fill0);
    arr_1d<nmax, COMPLEX> w(fill0);
    arr_1d<nmax, REAL> rw(fill0);
    REAL rcond = 0.0;
    REAL scale = 0.0;
    if (Mlsamen2, c2, "TR") {
        //
        //        ZTRTRI
        //
        srnamt = "ZTRTRI";
        infot = 1;
        ztrtri("/", "N", 0, a, 1, info);
        chkxer("ZTRTRI", infot, nout, lerr, ok);
        infot = 2;
        ztrtri("U", "/", 0, a, 1, info);
        chkxer("ZTRTRI", infot, nout, lerr, ok);
        infot = 3;
        ztrtri("U", "N", -1, a, 1, info);
        chkxer("ZTRTRI", infot, nout, lerr, ok);
        infot = 5;
        ztrtri("U", "N", 2, a, 1, info);
        chkxer("ZTRTRI", infot, nout, lerr, ok);
        //
        //        ZTRTI2
        //
        srnamt = "ZTRTI2";
        infot = 1;
        ztrti2("/", "N", 0, a, 1, info);
        chkxer("ZTRTI2", infot, nout, lerr, ok);
        infot = 2;
        ztrti2("U", "/", 0, a, 1, info);
        chkxer("ZTRTI2", infot, nout, lerr, ok);
        infot = 3;
        ztrti2("U", "N", -1, a, 1, info);
        chkxer("ZTRTI2", infot, nout, lerr, ok);
        infot = 5;
        ztrti2("U", "N", 2, a, 1, info);
        chkxer("ZTRTI2", infot, nout, lerr, ok);
        //
        //        ZTRTRS
        //
        srnamt = "ZTRTRS";
        infot = 1;
        ztrtrs("/", "N", "N", 0, 0, a, 1, x, 1, info);
        chkxer("ZTRTRS", infot, nout, lerr, ok);
        infot = 2;
        ztrtrs("U", "/", "N", 0, 0, a, 1, x, 1, info);
        chkxer("ZTRTRS", infot, nout, lerr, ok);
        infot = 3;
        ztrtrs("U", "N", "/", 0, 0, a, 1, x, 1, info);
        chkxer("ZTRTRS", infot, nout, lerr, ok);
        infot = 4;
        ztrtrs("U", "N", "N", -1, 0, a, 1, x, 1, info);
        chkxer("ZTRTRS", infot, nout, lerr, ok);
        infot = 5;
        ztrtrs("U", "N", "N", 0, -1, a, 1, x, 1, info);
        chkxer("ZTRTRS", infot, nout, lerr, ok);
        infot = 7;
        //
        //        ZTRRFS
        //
        srnamt = "ZTRRFS";
        infot = 1;
        ztrrfs("/", "N", "N", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTRRFS", infot, nout, lerr, ok);
        infot = 2;
        ztrrfs("U", "/", "N", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTRRFS", infot, nout, lerr, ok);
        infot = 3;
        ztrrfs("U", "N", "/", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTRRFS", infot, nout, lerr, ok);
        infot = 4;
        ztrrfs("U", "N", "N", -1, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTRRFS", infot, nout, lerr, ok);
        infot = 5;
        ztrrfs("U", "N", "N", 0, -1, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTRRFS", infot, nout, lerr, ok);
        infot = 7;
        ztrrfs("U", "N", "N", 2, 1, a, 1, b, 2, x, 2, r1, r2, w, rw, info);
        chkxer("ZTRRFS", infot, nout, lerr, ok);
        infot = 9;
        ztrrfs("U", "N", "N", 2, 1, a, 2, b, 1, x, 2, r1, r2, w, rw, info);
        chkxer("ZTRRFS", infot, nout, lerr, ok);
        infot = 11;
        ztrrfs("U", "N", "N", 2, 1, a, 2, b, 2, x, 1, r1, r2, w, rw, info);
        chkxer("ZTRRFS", infot, nout, lerr, ok);
        //
        //        ZTRCON
        //
        srnamt = "ZTRCON";
        infot = 1;
        ztrcon("/", "U", "N", 0, a, 1, rcond, w, rw, info);
        chkxer("ZTRCON", infot, nout, lerr, ok);
        infot = 2;
        ztrcon("1", "/", "N", 0, a, 1, rcond, w, rw, info);
        chkxer("ZTRCON", infot, nout, lerr, ok);
        infot = 3;
        ztrcon("1", "U", "/", 0, a, 1, rcond, w, rw, info);
        chkxer("ZTRCON", infot, nout, lerr, ok);
        infot = 4;
        ztrcon("1", "U", "N", -1, a, 1, rcond, w, rw, info);
        chkxer("ZTRCON", infot, nout, lerr, ok);
        infot = 6;
        ztrcon("1", "U", "N", 2, a, 1, rcond, w, rw, info);
        chkxer("ZTRCON", infot, nout, lerr, ok);
        //
        //        ZLATRS
        //
        srnamt = "ZLATRS";
        infot = 1;
        zlatrs("/", "N", "N", "N", 0, a, 1, x, scale, rw, info);
        chkxer("ZLATRS", infot, nout, lerr, ok);
        infot = 2;
        zlatrs("U", "/", "N", "N", 0, a, 1, x, scale, rw, info);
        chkxer("ZLATRS", infot, nout, lerr, ok);
        infot = 3;
        zlatrs("U", "N", "/", "N", 0, a, 1, x, scale, rw, info);
        chkxer("ZLATRS", infot, nout, lerr, ok);
        infot = 4;
        zlatrs("U", "N", "N", "/", 0, a, 1, x, scale, rw, info);
        chkxer("ZLATRS", infot, nout, lerr, ok);
        infot = 5;
        zlatrs("U", "N", "N", "N", -1, a, 1, x, scale, rw, info);
        chkxer("ZLATRS", infot, nout, lerr, ok);
        infot = 7;
        zlatrs("U", "N", "N", "N", 2, a, 1, x, scale, rw, info);
        chkxer("ZLATRS", infot, nout, lerr, ok);
        //
        //     Test error exits for the packed triangular routines.
        //
    } else if (Mlsamen2, c2, "TP") {
        //
        //        ZTPTRI
        //
        srnamt = "ZTPTRI";
        infot = 1;
        ztptri("/", "N", 0, a, info);
        chkxer("ZTPTRI", infot, nout, lerr, ok);
        infot = 2;
        ztptri("U", "/", 0, a, info);
        chkxer("ZTPTRI", infot, nout, lerr, ok);
        infot = 3;
        ztptri("U", "N", -1, a, info);
        chkxer("ZTPTRI", infot, nout, lerr, ok);
        //
        //        ZTPTRS
        //
        srnamt = "ZTPTRS";
        infot = 1;
        ztptrs("/", "N", "N", 0, 0, a, x, 1, info);
        chkxer("ZTPTRS", infot, nout, lerr, ok);
        infot = 2;
        ztptrs("U", "/", "N", 0, 0, a, x, 1, info);
        chkxer("ZTPTRS", infot, nout, lerr, ok);
        infot = 3;
        ztptrs("U", "N", "/", 0, 0, a, x, 1, info);
        chkxer("ZTPTRS", infot, nout, lerr, ok);
        infot = 4;
        ztptrs("U", "N", "N", -1, 0, a, x, 1, info);
        chkxer("ZTPTRS", infot, nout, lerr, ok);
        infot = 5;
        ztptrs("U", "N", "N", 0, -1, a, x, 1, info);
        chkxer("ZTPTRS", infot, nout, lerr, ok);
        infot = 8;
        ztptrs("U", "N", "N", 2, 1, a, x, 1, info);
        chkxer("ZTPTRS", infot, nout, lerr, ok);
        //
        //        ZTPRFS
        //
        srnamt = "ZTPRFS";
        infot = 1;
        ztprfs("/", "N", "N", 0, 0, a, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTPRFS", infot, nout, lerr, ok);
        infot = 2;
        ztprfs("U", "/", "N", 0, 0, a, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTPRFS", infot, nout, lerr, ok);
        infot = 3;
        ztprfs("U", "N", "/", 0, 0, a, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTPRFS", infot, nout, lerr, ok);
        infot = 4;
        ztprfs("U", "N", "N", -1, 0, a, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTPRFS", infot, nout, lerr, ok);
        infot = 5;
        ztprfs("U", "N", "N", 0, -1, a, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTPRFS", infot, nout, lerr, ok);
        infot = 8;
        ztprfs("U", "N", "N", 2, 1, a, b, 1, x, 2, r1, r2, w, rw, info);
        chkxer("ZTPRFS", infot, nout, lerr, ok);
        infot = 10;
        ztprfs("U", "N", "N", 2, 1, a, b, 2, x, 1, r1, r2, w, rw, info);
        chkxer("ZTPRFS", infot, nout, lerr, ok);
        //
        //        ZTPCON
        //
        srnamt = "ZTPCON";
        infot = 1;
        ztpcon("/", "U", "N", 0, a, rcond, w, rw, info);
        chkxer("ZTPCON", infot, nout, lerr, ok);
        infot = 2;
        ztpcon("1", "/", "N", 0, a, rcond, w, rw, info);
        chkxer("ZTPCON", infot, nout, lerr, ok);
        infot = 3;
        ztpcon("1", "U", "/", 0, a, rcond, w, rw, info);
        chkxer("ZTPCON", infot, nout, lerr, ok);
        infot = 4;
        ztpcon("1", "U", "N", -1, a, rcond, w, rw, info);
        chkxer("ZTPCON", infot, nout, lerr, ok);
        //
        //        ZLATPS
        //
        srnamt = "ZLATPS";
        infot = 1;
        zlatps("/", "N", "N", "N", 0, a, x, scale, rw, info);
        chkxer("ZLATPS", infot, nout, lerr, ok);
        infot = 2;
        zlatps("U", "/", "N", "N", 0, a, x, scale, rw, info);
        chkxer("ZLATPS", infot, nout, lerr, ok);
        infot = 3;
        zlatps("U", "N", "/", "N", 0, a, x, scale, rw, info);
        chkxer("ZLATPS", infot, nout, lerr, ok);
        infot = 4;
        zlatps("U", "N", "N", "/", 0, a, x, scale, rw, info);
        chkxer("ZLATPS", infot, nout, lerr, ok);
        infot = 5;
        zlatps("U", "N", "N", "N", -1, a, x, scale, rw, info);
        chkxer("ZLATPS", infot, nout, lerr, ok);
        //
        //     Test error exits for the banded triangular routines.
        //
    } else if (Mlsamen2, c2, "TB") {
        //
        //        ZTBTRS
        //
        srnamt = "ZTBTRS";
        infot = 1;
        ztbtrs("/", "N", "N", 0, 0, 0, a, 1, x, 1, info);
        chkxer("ZTBTRS", infot, nout, lerr, ok);
        infot = 2;
        ztbtrs("U", "/", "N", 0, 0, 0, a, 1, x, 1, info);
        chkxer("ZTBTRS", infot, nout, lerr, ok);
        infot = 3;
        ztbtrs("U", "N", "/", 0, 0, 0, a, 1, x, 1, info);
        chkxer("ZTBTRS", infot, nout, lerr, ok);
        infot = 4;
        ztbtrs("U", "N", "N", -1, 0, 0, a, 1, x, 1, info);
        chkxer("ZTBTRS", infot, nout, lerr, ok);
        infot = 5;
        ztbtrs("U", "N", "N", 0, -1, 0, a, 1, x, 1, info);
        chkxer("ZTBTRS", infot, nout, lerr, ok);
        infot = 6;
        ztbtrs("U", "N", "N", 0, 0, -1, a, 1, x, 1, info);
        chkxer("ZTBTRS", infot, nout, lerr, ok);
        infot = 8;
        ztbtrs("U", "N", "N", 2, 1, 1, a, 1, x, 2, info);
        chkxer("ZTBTRS", infot, nout, lerr, ok);
        infot = 10;
        ztbtrs("U", "N", "N", 2, 0, 1, a, 1, x, 1, info);
        chkxer("ZTBTRS", infot, nout, lerr, ok);
        //
        //        ZTBRFS
        //
        srnamt = "ZTBRFS";
        infot = 1;
        ztbrfs("/", "N", "N", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTBRFS", infot, nout, lerr, ok);
        infot = 2;
        ztbrfs("U", "/", "N", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTBRFS", infot, nout, lerr, ok);
        infot = 3;
        ztbrfs("U", "N", "/", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTBRFS", infot, nout, lerr, ok);
        infot = 4;
        ztbrfs("U", "N", "N", -1, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTBRFS", infot, nout, lerr, ok);
        infot = 5;
        ztbrfs("U", "N", "N", 0, -1, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTBRFS", infot, nout, lerr, ok);
        infot = 6;
        ztbrfs("U", "N", "N", 0, 0, -1, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("ZTBRFS", infot, nout, lerr, ok);
        infot = 8;
        ztbrfs("U", "N", "N", 2, 1, 1, a, 1, b, 2, x, 2, r1, r2, w, rw, info);
        chkxer("ZTBRFS", infot, nout, lerr, ok);
        infot = 10;
        ztbrfs("U", "N", "N", 2, 1, 1, a, 2, b, 1, x, 2, r1, r2, w, rw, info);
        chkxer("ZTBRFS", infot, nout, lerr, ok);
        infot = 12;
        ztbrfs("U", "N", "N", 2, 1, 1, a, 2, b, 2, x, 1, r1, r2, w, rw, info);
        chkxer("ZTBRFS", infot, nout, lerr, ok);
        //
        //        ZTBCON
        //
        srnamt = "ZTBCON";
        infot = 1;
        ztbcon("/", "U", "N", 0, 0, a, 1, rcond, w, rw, info);
        chkxer("ZTBCON", infot, nout, lerr, ok);
        infot = 2;
        ztbcon("1", "/", "N", 0, 0, a, 1, rcond, w, rw, info);
        chkxer("ZTBCON", infot, nout, lerr, ok);
        infot = 3;
        ztbcon("1", "U", "/", 0, 0, a, 1, rcond, w, rw, info);
        chkxer("ZTBCON", infot, nout, lerr, ok);
        infot = 4;
        ztbcon("1", "U", "N", -1, 0, a, 1, rcond, w, rw, info);
        chkxer("ZTBCON", infot, nout, lerr, ok);
        infot = 5;
        ztbcon("1", "U", "N", 0, -1, a, 1, rcond, w, rw, info);
        chkxer("ZTBCON", infot, nout, lerr, ok);
        infot = 7;
        ztbcon("1", "U", "N", 2, 1, a, 1, rcond, w, rw, info);
        chkxer("ZTBCON", infot, nout, lerr, ok);
        //
        //        ZLATBS
        //
        srnamt = "ZLATBS";
        infot = 1;
        zlatbs("/", "N", "N", "N", 0, 0, a, 1, x, scale, rw, info);
        chkxer("ZLATBS", infot, nout, lerr, ok);
        infot = 2;
        zlatbs("U", "/", "N", "N", 0, 0, a, 1, x, scale, rw, info);
        chkxer("ZLATBS", infot, nout, lerr, ok);
        infot = 3;
        zlatbs("U", "N", "/", "N", 0, 0, a, 1, x, scale, rw, info);
        chkxer("ZLATBS", infot, nout, lerr, ok);
        infot = 4;
        zlatbs("U", "N", "N", "/", 0, 0, a, 1, x, scale, rw, info);
        chkxer("ZLATBS", infot, nout, lerr, ok);
        infot = 5;
        zlatbs("U", "N", "N", "N", -1, 0, a, 1, x, scale, rw, info);
        chkxer("ZLATBS", infot, nout, lerr, ok);
        infot = 6;
        zlatbs("U", "N", "N", "N", 1, -1, a, 1, x, scale, rw, info);
        chkxer("ZLATBS", infot, nout, lerr, ok);
        infot = 8;
        zlatbs("U", "N", "N", "N", 2, 1, a, 1, x, scale, rw, info);
        chkxer("ZLATBS", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrtr
    //
}
