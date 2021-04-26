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

void Rerrtr(common &cmn, const char *path, INTEGER const nunit) {
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
    a[(2 - 1) * lda] = 2.0;
    a[(2 - 1) + (2 - 1) * lda] = 3.e0;
    a[(2 - 1)] = 4.e0;
    ok = true;
    //
    INTEGER info = 0;
    arr_1d<nmax, REAL> x(fill0);
    arr_1d<nmax, REAL> b(fill0);
    arr_1d<nmax, REAL> r1(fill0);
    arr_1d<nmax, REAL> r2(fill0);
    arr_1d<nmax, REAL> w(fill0);
    arr_1d<nmax, int> iw(fill0);
    REAL rcond = 0.0;
    REAL scale = 0.0;
    if (Mlsamen2, c2, "TR") {
        //
        //        Test error exits for the general triangular routines.
        //
        //        DTRTRI
        //
        srnamt = "DTRTRI";
        infot = 1;
        dtrtri("/", "N", 0, a, 1, info);
        chkxer("DTRTRI", infot, nout, lerr, ok);
        infot = 2;
        dtrtri("U", "/", 0, a, 1, info);
        chkxer("DTRTRI", infot, nout, lerr, ok);
        infot = 3;
        dtrtri("U", "N", -1, a, 1, info);
        chkxer("DTRTRI", infot, nout, lerr, ok);
        infot = 5;
        dtrtri("U", "N", 2, a, 1, info);
        chkxer("DTRTRI", infot, nout, lerr, ok);
        //
        //        DTRTI2
        //
        srnamt = "DTRTI2";
        infot = 1;
        dtrti2("/", "N", 0, a, 1, info);
        chkxer("DTRTI2", infot, nout, lerr, ok);
        infot = 2;
        dtrti2("U", "/", 0, a, 1, info);
        chkxer("DTRTI2", infot, nout, lerr, ok);
        infot = 3;
        dtrti2("U", "N", -1, a, 1, info);
        chkxer("DTRTI2", infot, nout, lerr, ok);
        infot = 5;
        dtrti2("U", "N", 2, a, 1, info);
        chkxer("DTRTI2", infot, nout, lerr, ok);
        //
        //        DTRTRS
        //
        srnamt = "DTRTRS";
        infot = 1;
        dtrtrs("/", "N", "N", 0, 0, a, 1, x, 1, info);
        chkxer("DTRTRS", infot, nout, lerr, ok);
        infot = 2;
        dtrtrs("U", "/", "N", 0, 0, a, 1, x, 1, info);
        chkxer("DTRTRS", infot, nout, lerr, ok);
        infot = 3;
        dtrtrs("U", "N", "/", 0, 0, a, 1, x, 1, info);
        chkxer("DTRTRS", infot, nout, lerr, ok);
        infot = 4;
        dtrtrs("U", "N", "N", -1, 0, a, 1, x, 1, info);
        chkxer("DTRTRS", infot, nout, lerr, ok);
        infot = 5;
        dtrtrs("U", "N", "N", 0, -1, a, 1, x, 1, info);
        chkxer("DTRTRS", infot, nout, lerr, ok);
        infot = 7;
        dtrtrs("U", "N", "N", 2, 1, a, 1, x, 2, info);
        chkxer("DTRTRS", infot, nout, lerr, ok);
        infot = 9;
        dtrtrs("U", "N", "N", 2, 1, a, 2, x, 1, info);
        chkxer("DTRTRS", infot, nout, lerr, ok);
        //
        //        DTRRFS
        //
        srnamt = "DTRRFS";
        infot = 1;
        dtrrfs("/", "N", "N", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTRRFS", infot, nout, lerr, ok);
        infot = 2;
        dtrrfs("U", "/", "N", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTRRFS", infot, nout, lerr, ok);
        infot = 3;
        dtrrfs("U", "N", "/", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTRRFS", infot, nout, lerr, ok);
        infot = 4;
        dtrrfs("U", "N", "N", -1, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTRRFS", infot, nout, lerr, ok);
        infot = 5;
        dtrrfs("U", "N", "N", 0, -1, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTRRFS", infot, nout, lerr, ok);
        infot = 7;
        dtrrfs("U", "N", "N", 2, 1, a, 1, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("DTRRFS", infot, nout, lerr, ok);
        infot = 9;
        dtrrfs("U", "N", "N", 2, 1, a, 2, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("DTRRFS", infot, nout, lerr, ok);
        infot = 11;
        dtrrfs("U", "N", "N", 2, 1, a, 2, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("DTRRFS", infot, nout, lerr, ok);
        //
        //        DTRCON
        //
        srnamt = "DTRCON";
        infot = 1;
        dtrcon("/", "U", "N", 0, a, 1, rcond, w, iw, info);
        chkxer("DTRCON", infot, nout, lerr, ok);
        infot = 2;
        dtrcon("1", "/", "N", 0, a, 1, rcond, w, iw, info);
        chkxer("DTRCON", infot, nout, lerr, ok);
        infot = 3;
        dtrcon("1", "U", "/", 0, a, 1, rcond, w, iw, info);
        chkxer("DTRCON", infot, nout, lerr, ok);
        infot = 4;
        dtrcon("1", "U", "N", -1, a, 1, rcond, w, iw, info);
        chkxer("DTRCON", infot, nout, lerr, ok);
        infot = 6;
        dtrcon("1", "U", "N", 2, a, 1, rcond, w, iw, info);
        chkxer("DTRCON", infot, nout, lerr, ok);
        //
        //        DLATRS
        //
        srnamt = "DLATRS";
        infot = 1;
        dlatrs("/", "N", "N", "N", 0, a, 1, x, scale, w, info);
        chkxer("DLATRS", infot, nout, lerr, ok);
        infot = 2;
        dlatrs("U", "/", "N", "N", 0, a, 1, x, scale, w, info);
        chkxer("DLATRS", infot, nout, lerr, ok);
        infot = 3;
        dlatrs("U", "N", "/", "N", 0, a, 1, x, scale, w, info);
        chkxer("DLATRS", infot, nout, lerr, ok);
        infot = 4;
        dlatrs("U", "N", "N", "/", 0, a, 1, x, scale, w, info);
        chkxer("DLATRS", infot, nout, lerr, ok);
        infot = 5;
        dlatrs("U", "N", "N", "N", -1, a, 1, x, scale, w, info);
        chkxer("DLATRS", infot, nout, lerr, ok);
        infot = 7;
        dlatrs("U", "N", "N", "N", 2, a, 1, x, scale, w, info);
        chkxer("DLATRS", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "TP") {
        //
        //        Test error exits for the packed triangular routines.
        //
        //        DTPTRI
        //
        srnamt = "DTPTRI";
        infot = 1;
        dtptri("/", "N", 0, a, info);
        chkxer("DTPTRI", infot, nout, lerr, ok);
        infot = 2;
        dtptri("U", "/", 0, a, info);
        chkxer("DTPTRI", infot, nout, lerr, ok);
        infot = 3;
        dtptri("U", "N", -1, a, info);
        chkxer("DTPTRI", infot, nout, lerr, ok);
        //
        //        DTPTRS
        //
        srnamt = "DTPTRS";
        infot = 1;
        dtptrs("/", "N", "N", 0, 0, a, x, 1, info);
        chkxer("DTPTRS", infot, nout, lerr, ok);
        infot = 2;
        dtptrs("U", "/", "N", 0, 0, a, x, 1, info);
        chkxer("DTPTRS", infot, nout, lerr, ok);
        infot = 3;
        dtptrs("U", "N", "/", 0, 0, a, x, 1, info);
        chkxer("DTPTRS", infot, nout, lerr, ok);
        infot = 4;
        dtptrs("U", "N", "N", -1, 0, a, x, 1, info);
        chkxer("DTPTRS", infot, nout, lerr, ok);
        infot = 5;
        dtptrs("U", "N", "N", 0, -1, a, x, 1, info);
        chkxer("DTPTRS", infot, nout, lerr, ok);
        infot = 8;
        dtptrs("U", "N", "N", 2, 1, a, x, 1, info);
        chkxer("DTPTRS", infot, nout, lerr, ok);
        //
        //        DTPRFS
        //
        srnamt = "DTPRFS";
        infot = 1;
        dtprfs("/", "N", "N", 0, 0, a, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTPRFS", infot, nout, lerr, ok);
        infot = 2;
        dtprfs("U", "/", "N", 0, 0, a, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTPRFS", infot, nout, lerr, ok);
        infot = 3;
        dtprfs("U", "N", "/", 0, 0, a, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTPRFS", infot, nout, lerr, ok);
        infot = 4;
        dtprfs("U", "N", "N", -1, 0, a, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTPRFS", infot, nout, lerr, ok);
        infot = 5;
        dtprfs("U", "N", "N", 0, -1, a, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTPRFS", infot, nout, lerr, ok);
        infot = 8;
        dtprfs("U", "N", "N", 2, 1, a, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("DTPRFS", infot, nout, lerr, ok);
        infot = 10;
        dtprfs("U", "N", "N", 2, 1, a, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("DTPRFS", infot, nout, lerr, ok);
        //
        //        DTPCON
        //
        srnamt = "DTPCON";
        infot = 1;
        dtpcon("/", "U", "N", 0, a, rcond, w, iw, info);
        chkxer("DTPCON", infot, nout, lerr, ok);
        infot = 2;
        dtpcon("1", "/", "N", 0, a, rcond, w, iw, info);
        chkxer("DTPCON", infot, nout, lerr, ok);
        infot = 3;
        dtpcon("1", "U", "/", 0, a, rcond, w, iw, info);
        chkxer("DTPCON", infot, nout, lerr, ok);
        infot = 4;
        dtpcon("1", "U", "N", -1, a, rcond, w, iw, info);
        chkxer("DTPCON", infot, nout, lerr, ok);
        //
        //        DLATPS
        //
        srnamt = "DLATPS";
        infot = 1;
        dlatps("/", "N", "N", "N", 0, a, x, scale, w, info);
        chkxer("DLATPS", infot, nout, lerr, ok);
        infot = 2;
        dlatps("U", "/", "N", "N", 0, a, x, scale, w, info);
        chkxer("DLATPS", infot, nout, lerr, ok);
        infot = 3;
        dlatps("U", "N", "/", "N", 0, a, x, scale, w, info);
        chkxer("DLATPS", infot, nout, lerr, ok);
        infot = 4;
        dlatps("U", "N", "N", "/", 0, a, x, scale, w, info);
        chkxer("DLATPS", infot, nout, lerr, ok);
        infot = 5;
        dlatps("U", "N", "N", "N", -1, a, x, scale, w, info);
        chkxer("DLATPS", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "TB") {
        //
        //        Test error exits for the banded triangular routines.
        //
        //        DTBTRS
        //
        srnamt = "DTBTRS";
        infot = 1;
        dtbtrs("/", "N", "N", 0, 0, 0, a, 1, x, 1, info);
        chkxer("DTBTRS", infot, nout, lerr, ok);
        infot = 2;
        dtbtrs("U", "/", "N", 0, 0, 0, a, 1, x, 1, info);
        chkxer("DTBTRS", infot, nout, lerr, ok);
        infot = 3;
        dtbtrs("U", "N", "/", 0, 0, 0, a, 1, x, 1, info);
        chkxer("DTBTRS", infot, nout, lerr, ok);
        infot = 4;
        dtbtrs("U", "N", "N", -1, 0, 0, a, 1, x, 1, info);
        chkxer("DTBTRS", infot, nout, lerr, ok);
        infot = 5;
        dtbtrs("U", "N", "N", 0, -1, 0, a, 1, x, 1, info);
        chkxer("DTBTRS", infot, nout, lerr, ok);
        infot = 6;
        dtbtrs("U", "N", "N", 0, 0, -1, a, 1, x, 1, info);
        chkxer("DTBTRS", infot, nout, lerr, ok);
        infot = 8;
        dtbtrs("U", "N", "N", 2, 1, 1, a, 1, x, 2, info);
        chkxer("DTBTRS", infot, nout, lerr, ok);
        infot = 10;
        dtbtrs("U", "N", "N", 2, 0, 1, a, 1, x, 1, info);
        chkxer("DTBTRS", infot, nout, lerr, ok);
        //
        //        DTBRFS
        //
        srnamt = "DTBRFS";
        infot = 1;
        dtbrfs("/", "N", "N", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTBRFS", infot, nout, lerr, ok);
        infot = 2;
        dtbrfs("U", "/", "N", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTBRFS", infot, nout, lerr, ok);
        infot = 3;
        dtbrfs("U", "N", "/", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTBRFS", infot, nout, lerr, ok);
        infot = 4;
        dtbrfs("U", "N", "N", -1, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTBRFS", infot, nout, lerr, ok);
        infot = 5;
        dtbrfs("U", "N", "N", 0, -1, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTBRFS", infot, nout, lerr, ok);
        infot = 6;
        dtbrfs("U", "N", "N", 0, 0, -1, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DTBRFS", infot, nout, lerr, ok);
        infot = 8;
        dtbrfs("U", "N", "N", 2, 1, 1, a, 1, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("DTBRFS", infot, nout, lerr, ok);
        infot = 10;
        dtbrfs("U", "N", "N", 2, 1, 1, a, 2, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("DTBRFS", infot, nout, lerr, ok);
        infot = 12;
        dtbrfs("U", "N", "N", 2, 1, 1, a, 2, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("DTBRFS", infot, nout, lerr, ok);
        //
        //        DTBCON
        //
        srnamt = "DTBCON";
        infot = 1;
        dtbcon("/", "U", "N", 0, 0, a, 1, rcond, w, iw, info);
        chkxer("DTBCON", infot, nout, lerr, ok);
        infot = 2;
        dtbcon("1", "/", "N", 0, 0, a, 1, rcond, w, iw, info);
        chkxer("DTBCON", infot, nout, lerr, ok);
        infot = 3;
        dtbcon("1", "U", "/", 0, 0, a, 1, rcond, w, iw, info);
        chkxer("DTBCON", infot, nout, lerr, ok);
        infot = 4;
        dtbcon("1", "U", "N", -1, 0, a, 1, rcond, w, iw, info);
        chkxer("DTBCON", infot, nout, lerr, ok);
        infot = 5;
        dtbcon("1", "U", "N", 0, -1, a, 1, rcond, w, iw, info);
        chkxer("DTBCON", infot, nout, lerr, ok);
        infot = 7;
        dtbcon("1", "U", "N", 2, 1, a, 1, rcond, w, iw, info);
        chkxer("DTBCON", infot, nout, lerr, ok);
        //
        //        DLATBS
        //
        srnamt = "DLATBS";
        infot = 1;
        dlatbs("/", "N", "N", "N", 0, 0, a, 1, x, scale, w, info);
        chkxer("DLATBS", infot, nout, lerr, ok);
        infot = 2;
        dlatbs("U", "/", "N", "N", 0, 0, a, 1, x, scale, w, info);
        chkxer("DLATBS", infot, nout, lerr, ok);
        infot = 3;
        dlatbs("U", "N", "/", "N", 0, 0, a, 1, x, scale, w, info);
        chkxer("DLATBS", infot, nout, lerr, ok);
        infot = 4;
        dlatbs("U", "N", "N", "/", 0, 0, a, 1, x, scale, w, info);
        chkxer("DLATBS", infot, nout, lerr, ok);
        infot = 5;
        dlatbs("U", "N", "N", "N", -1, 0, a, 1, x, scale, w, info);
        chkxer("DLATBS", infot, nout, lerr, ok);
        infot = 6;
        dlatbs("U", "N", "N", "N", 1, -1, a, 1, x, scale, w, info);
        chkxer("DLATBS", infot, nout, lerr, ok);
        infot = 8;
        dlatbs("U", "N", "N", "N", 2, 1, a, 1, x, scale, w, info);
        chkxer("DLATBS", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrtr
    //
}
