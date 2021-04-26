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

void Cerrpo(common &cmn, const char *path, INTEGER const nunit) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    nout = nunit;
    write(nout, star);
    str<2> c2 = path[(2 - 1) + (3 - 1) * ldpath];
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 4;
    INTEGER i = 0;
    arr_2d<nmax, nmax, COMPLEX> a(fill0);
    arr_2d<nmax, nmax, COMPLEX> af(fill0);
    arr_1d<nmax, COMPLEX> b(fill0);
    arr_1d<nmax, REAL> r1(fill0);
    arr_1d<nmax, REAL> r2(fill0);
    arr_1d<2 * nmax, COMPLEX> w(fill0);
    arr_1d<nmax, COMPLEX> x(fill0);
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = COMPLEX(1.0 / (i + j).real(), -1.0 / (i + j).real());
            af[(i - 1) + (j - 1) * ldaf] = COMPLEX(1.0 / (i + j).real(), -1.0 / (i + j).real());
        }
        b[j - 1] = 0.0;
        r1[j - 1] = 0.0;
        r2[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
    }
    REAL anrm = 1.0;
    ok = true;
    //
    //     Test error exits of the routines that use the Cholesky
    //     decomposition of a Hermitian positive definite matrix.
    //
    INTEGER info = 0;
    arr_1d<nmax, REAL> r(fill0);
    REAL rcond = 0.0;
    if (Mlsamen2, c2, "PO") {
        //
        //        ZPOTRF
        //
        srnamt = "ZPOTRF";
        infot = 1;
        zpotrf("/", 0, a, 1, info);
        chkxer("ZPOTRF", infot, nout, lerr, ok);
        infot = 2;
        zpotrf("U", -1, a, 1, info);
        chkxer("ZPOTRF", infot, nout, lerr, ok);
        infot = 4;
        zpotrf("U", 2, a, 1, info);
        chkxer("ZPOTRF", infot, nout, lerr, ok);
        //
        //        ZPOTF2
        //
        srnamt = "ZPOTF2";
        infot = 1;
        zpotf2("/", 0, a, 1, info);
        chkxer("ZPOTF2", infot, nout, lerr, ok);
        infot = 2;
        zpotf2("U", -1, a, 1, info);
        chkxer("ZPOTF2", infot, nout, lerr, ok);
        infot = 4;
        zpotf2("U", 2, a, 1, info);
        chkxer("ZPOTF2", infot, nout, lerr, ok);
        //
        //        ZPOTRI
        //
        srnamt = "ZPOTRI";
        infot = 1;
        zpotri("/", 0, a, 1, info);
        chkxer("ZPOTRI", infot, nout, lerr, ok);
        infot = 2;
        zpotri("U", -1, a, 1, info);
        chkxer("ZPOTRI", infot, nout, lerr, ok);
        infot = 4;
        zpotri("U", 2, a, 1, info);
        chkxer("ZPOTRI", infot, nout, lerr, ok);
        //
        //        ZPOTRS
        //
        srnamt = "ZPOTRS";
        infot = 1;
        zpotrs("/", 0, 0, a, 1, b, 1, info);
        chkxer("ZPOTRS", infot, nout, lerr, ok);
        infot = 2;
        zpotrs("U", -1, 0, a, 1, b, 1, info);
        chkxer("ZPOTRS", infot, nout, lerr, ok);
        infot = 3;
        zpotrs("U", 0, -1, a, 1, b, 1, info);
        chkxer("ZPOTRS", infot, nout, lerr, ok);
        infot = 5;
        zpotrs("U", 2, 1, a, 1, b, 2, info);
        chkxer("ZPOTRS", infot, nout, lerr, ok);
        infot = 7;
        zpotrs("U", 2, 1, a, 2, b, 1, info);
        chkxer("ZPOTRS", infot, nout, lerr, ok);
        //
        //        ZPORFS
        //
        srnamt = "ZPORFS";
        infot = 1;
        zporfs("/", 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZPORFS", infot, nout, lerr, ok);
        infot = 2;
        zporfs("U", -1, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZPORFS", infot, nout, lerr, ok);
        infot = 3;
        zporfs("U", 0, -1, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZPORFS", infot, nout, lerr, ok);
        infot = 5;
        zporfs("U", 2, 1, a, 1, af, 2, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("ZPORFS", infot, nout, lerr, ok);
        infot = 7;
        zporfs("U", 2, 1, a, 2, af, 1, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("ZPORFS", infot, nout, lerr, ok);
        infot = 9;
        zporfs("U", 2, 1, a, 2, af, 2, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("ZPORFS", infot, nout, lerr, ok);
        infot = 11;
        zporfs("U", 2, 1, a, 2, af, 2, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("ZPORFS", infot, nout, lerr, ok);
        //
        //        ZPOCON
        //
        srnamt = "ZPOCON";
        infot = 1;
        zpocon("/", 0, a, 1, anrm, rcond, w, r, info);
        chkxer("ZPOCON", infot, nout, lerr, ok);
        infot = 2;
        zpocon("U", -1, a, 1, anrm, rcond, w, r, info);
        chkxer("ZPOCON", infot, nout, lerr, ok);
        infot = 4;
        zpocon("U", 2, a, 1, anrm, rcond, w, r, info);
        chkxer("ZPOCON", infot, nout, lerr, ok);
        infot = 5;
        zpocon("U", 1, a, 1, -anrm, rcond, w, r, info);
        chkxer("ZPOCON", infot, nout, lerr, ok);
        //
        //        ZPOEQU
        //
        srnamt = "ZPOEQU";
        infot = 1;
        zpoequ(-1, a, 1, r1, rcond, anrm, info);
        chkxer("ZPOEQU", infot, nout, lerr, ok);
        infot = 3;
        zpoequ(2, a, 1, r1, rcond, anrm, info);
        chkxer("ZPOEQU", infot, nout, lerr, ok);
        //
        //     Test error exits of the routines that use the Cholesky
        //     decomposition of a Hermitian positive definite packed matrix.
        //
    } else if (Mlsamen2, c2, "PP") {
        //
        //        ZPPTRF
        //
        srnamt = "ZPPTRF";
        infot = 1;
        zpptrf("/", 0, a, info);
        chkxer("ZPPTRF", infot, nout, lerr, ok);
        infot = 2;
        zpptrf("U", -1, a, info);
        chkxer("ZPPTRF", infot, nout, lerr, ok);
        //
        //        ZPPTRI
        //
        srnamt = "ZPPTRI";
        infot = 1;
        zpptri("/", 0, a, info);
        chkxer("ZPPTRI", infot, nout, lerr, ok);
        infot = 2;
        zpptri("U", -1, a, info);
        chkxer("ZPPTRI", infot, nout, lerr, ok);
        //
        //        ZPPTRS
        //
        srnamt = "ZPPTRS";
        infot = 1;
        zpptrs("/", 0, 0, a, b, 1, info);
        chkxer("ZPPTRS", infot, nout, lerr, ok);
        infot = 2;
        zpptrs("U", -1, 0, a, b, 1, info);
        chkxer("ZPPTRS", infot, nout, lerr, ok);
        infot = 3;
        zpptrs("U", 0, -1, a, b, 1, info);
        chkxer("ZPPTRS", infot, nout, lerr, ok);
        infot = 6;
        zpptrs("U", 2, 1, a, b, 1, info);
        chkxer("ZPPTRS", infot, nout, lerr, ok);
        //
        //        ZPPRFS
        //
        srnamt = "ZPPRFS";
        infot = 1;
        zpprfs("/", 0, 0, a, af, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZPPRFS", infot, nout, lerr, ok);
        infot = 2;
        zpprfs("U", -1, 0, a, af, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZPPRFS", infot, nout, lerr, ok);
        infot = 3;
        zpprfs("U", 0, -1, a, af, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZPPRFS", infot, nout, lerr, ok);
        infot = 7;
        zpprfs("U", 2, 1, a, af, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("ZPPRFS", infot, nout, lerr, ok);
        infot = 9;
        zpprfs("U", 2, 1, a, af, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("ZPPRFS", infot, nout, lerr, ok);
        //
        //        ZPPCON
        //
        srnamt = "ZPPCON";
        infot = 1;
        zppcon("/", 0, a, anrm, rcond, w, r, info);
        chkxer("ZPPCON", infot, nout, lerr, ok);
        infot = 2;
        zppcon("U", -1, a, anrm, rcond, w, r, info);
        chkxer("ZPPCON", infot, nout, lerr, ok);
        infot = 4;
        zppcon("U", 1, a, -anrm, rcond, w, r, info);
        chkxer("ZPPCON", infot, nout, lerr, ok);
        //
        //        ZPPEQU
        //
        srnamt = "ZPPEQU";
        infot = 1;
        zppequ("/", 0, a, r1, rcond, anrm, info);
        chkxer("ZPPEQU", infot, nout, lerr, ok);
        infot = 2;
        zppequ("U", -1, a, r1, rcond, anrm, info);
        chkxer("ZPPEQU", infot, nout, lerr, ok);
        //
        //     Test error exits of the routines that use the Cholesky
        //     decomposition of a Hermitian positive definite band matrix.
        //
    } else if (Mlsamen2, c2, "PB") {
        //
        //        ZPBTRF
        //
        srnamt = "ZPBTRF";
        infot = 1;
        zpbtrf("/", 0, 0, a, 1, info);
        chkxer("ZPBTRF", infot, nout, lerr, ok);
        infot = 2;
        zpbtrf("U", -1, 0, a, 1, info);
        chkxer("ZPBTRF", infot, nout, lerr, ok);
        infot = 3;
        zpbtrf("U", 1, -1, a, 1, info);
        chkxer("ZPBTRF", infot, nout, lerr, ok);
        infot = 5;
        zpbtrf("U", 2, 1, a, 1, info);
        chkxer("ZPBTRF", infot, nout, lerr, ok);
        //
        //        ZPBTF2
        //
        srnamt = "ZPBTF2";
        infot = 1;
        zpbtf2("/", 0, 0, a, 1, info);
        chkxer("ZPBTF2", infot, nout, lerr, ok);
        infot = 2;
        zpbtf2("U", -1, 0, a, 1, info);
        chkxer("ZPBTF2", infot, nout, lerr, ok);
        infot = 3;
        zpbtf2("U", 1, -1, a, 1, info);
        chkxer("ZPBTF2", infot, nout, lerr, ok);
        infot = 5;
        zpbtf2("U", 2, 1, a, 1, info);
        chkxer("ZPBTF2", infot, nout, lerr, ok);
        //
        //        ZPBTRS
        //
        srnamt = "ZPBTRS";
        infot = 1;
        zpbtrs("/", 0, 0, 0, a, 1, b, 1, info);
        chkxer("ZPBTRS", infot, nout, lerr, ok);
        infot = 2;
        zpbtrs("U", -1, 0, 0, a, 1, b, 1, info);
        chkxer("ZPBTRS", infot, nout, lerr, ok);
        infot = 3;
        zpbtrs("U", 1, -1, 0, a, 1, b, 1, info);
        chkxer("ZPBTRS", infot, nout, lerr, ok);
        infot = 4;
        zpbtrs("U", 0, 0, -1, a, 1, b, 1, info);
        chkxer("ZPBTRS", infot, nout, lerr, ok);
        infot = 6;
        zpbtrs("U", 2, 1, 1, a, 1, b, 1, info);
        chkxer("ZPBTRS", infot, nout, lerr, ok);
        infot = 8;
        zpbtrs("U", 2, 0, 1, a, 1, b, 1, info);
        chkxer("ZPBTRS", infot, nout, lerr, ok);
        //
        //        ZPBRFS
        //
        srnamt = "ZPBRFS";
        infot = 1;
        zpbrfs("/", 0, 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZPBRFS", infot, nout, lerr, ok);
        infot = 2;
        zpbrfs("U", -1, 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZPBRFS", infot, nout, lerr, ok);
        infot = 3;
        zpbrfs("U", 1, -1, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZPBRFS", infot, nout, lerr, ok);
        infot = 4;
        zpbrfs("U", 0, 0, -1, a, 1, af, 1, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZPBRFS", infot, nout, lerr, ok);
        infot = 6;
        zpbrfs("U", 2, 1, 1, a, 1, af, 2, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("ZPBRFS", infot, nout, lerr, ok);
        infot = 8;
        zpbrfs("U", 2, 1, 1, a, 2, af, 1, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("ZPBRFS", infot, nout, lerr, ok);
        infot = 10;
        zpbrfs("U", 2, 0, 1, a, 1, af, 1, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("ZPBRFS", infot, nout, lerr, ok);
        infot = 12;
        zpbrfs("U", 2, 0, 1, a, 1, af, 1, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("ZPBRFS", infot, nout, lerr, ok);
        //
        //        ZPBCON
        //
        srnamt = "ZPBCON";
        infot = 1;
        zpbcon("/", 0, 0, a, 1, anrm, rcond, w, r, info);
        chkxer("ZPBCON", infot, nout, lerr, ok);
        infot = 2;
        zpbcon("U", -1, 0, a, 1, anrm, rcond, w, r, info);
        chkxer("ZPBCON", infot, nout, lerr, ok);
        infot = 3;
        zpbcon("U", 1, -1, a, 1, anrm, rcond, w, r, info);
        chkxer("ZPBCON", infot, nout, lerr, ok);
        infot = 5;
        zpbcon("U", 2, 1, a, 1, anrm, rcond, w, r, info);
        chkxer("ZPBCON", infot, nout, lerr, ok);
        infot = 6;
        zpbcon("U", 1, 0, a, 1, -anrm, rcond, w, r, info);
        chkxer("ZPBCON", infot, nout, lerr, ok);
        //
        //        ZPBEQU
        //
        srnamt = "ZPBEQU";
        infot = 1;
        zpbequ("/", 0, 0, a, 1, r1, rcond, anrm, info);
        chkxer("ZPBEQU", infot, nout, lerr, ok);
        infot = 2;
        zpbequ("U", -1, 0, a, 1, r1, rcond, anrm, info);
        chkxer("ZPBEQU", infot, nout, lerr, ok);
        infot = 3;
        zpbequ("U", 1, -1, a, 1, r1, rcond, anrm, info);
        chkxer("ZPBEQU", infot, nout, lerr, ok);
        infot = 5;
        zpbequ("U", 2, 1, a, 1, r1, rcond, anrm, info);
        chkxer("ZPBEQU", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrpo
    //
}
