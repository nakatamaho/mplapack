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

void Rerrpo(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_2d<nmax, nmax, REAL> a(fill0);
    arr_2d<nmax, nmax, REAL> af(fill0);
    arr_1d<nmax, REAL> b(fill0);
    arr_1d<nmax, REAL> r1(fill0);
    arr_1d<nmax, REAL> r2(fill0);
    arr_1d<3 * nmax, REAL> w(fill0);
    arr_1d<nmax, REAL> x(fill0);
    arr_1d<nmax, int> iw(fill0);
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / (i + j).real();
            af[(i - 1) + (j - 1) * ldaf] = 1.0 / (i + j).real();
        }
        b[j - 1] = 0.0;
        r1[j - 1] = 0.0;
        r2[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
        iw[j - 1] = j;
    }
    ok = true;
    //
    INTEGER info = 0;
    REAL anrm = 0.0;
    REAL rcond = 0.0;
    if (Mlsamen2, c2, "PO") {
        //
        //        Test error exits of the routines that use the Cholesky
        //        decomposition of a symmetric positive definite matrix.
        //
        //        DPOTRF
        //
        srnamt = "DPOTRF";
        infot = 1;
        dpotrf("/", 0, a, 1, info);
        chkxer("DPOTRF", infot, nout, lerr, ok);
        infot = 2;
        dpotrf("U", -1, a, 1, info);
        chkxer("DPOTRF", infot, nout, lerr, ok);
        infot = 4;
        dpotrf("U", 2, a, 1, info);
        chkxer("DPOTRF", infot, nout, lerr, ok);
        //
        //        DPOTF2
        //
        srnamt = "DPOTF2";
        infot = 1;
        dpotf2("/", 0, a, 1, info);
        chkxer("DPOTF2", infot, nout, lerr, ok);
        infot = 2;
        dpotf2("U", -1, a, 1, info);
        chkxer("DPOTF2", infot, nout, lerr, ok);
        infot = 4;
        dpotf2("U", 2, a, 1, info);
        chkxer("DPOTF2", infot, nout, lerr, ok);
        //
        //        DPOTRI
        //
        srnamt = "DPOTRI";
        infot = 1;
        dpotri("/", 0, a, 1, info);
        chkxer("DPOTRI", infot, nout, lerr, ok);
        infot = 2;
        dpotri("U", -1, a, 1, info);
        chkxer("DPOTRI", infot, nout, lerr, ok);
        infot = 4;
        dpotri("U", 2, a, 1, info);
        chkxer("DPOTRI", infot, nout, lerr, ok);
        //
        //        DPOTRS
        //
        srnamt = "DPOTRS";
        infot = 1;
        dpotrs("/", 0, 0, a, 1, b, 1, info);
        chkxer("DPOTRS", infot, nout, lerr, ok);
        infot = 2;
        dpotrs("U", -1, 0, a, 1, b, 1, info);
        chkxer("DPOTRS", infot, nout, lerr, ok);
        infot = 3;
        dpotrs("U", 0, -1, a, 1, b, 1, info);
        chkxer("DPOTRS", infot, nout, lerr, ok);
        infot = 5;
        dpotrs("U", 2, 1, a, 1, b, 2, info);
        chkxer("DPOTRS", infot, nout, lerr, ok);
        infot = 7;
        dpotrs("U", 2, 1, a, 2, b, 1, info);
        chkxer("DPOTRS", infot, nout, lerr, ok);
        //
        //        DPORFS
        //
        srnamt = "DPORFS";
        infot = 1;
        dporfs("/", 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DPORFS", infot, nout, lerr, ok);
        infot = 2;
        dporfs("U", -1, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DPORFS", infot, nout, lerr, ok);
        infot = 3;
        dporfs("U", 0, -1, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DPORFS", infot, nout, lerr, ok);
        infot = 5;
        dporfs("U", 2, 1, a, 1, af, 2, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("DPORFS", infot, nout, lerr, ok);
        infot = 7;
        dporfs("U", 2, 1, a, 2, af, 1, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("DPORFS", infot, nout, lerr, ok);
        infot = 9;
        dporfs("U", 2, 1, a, 2, af, 2, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("DPORFS", infot, nout, lerr, ok);
        infot = 11;
        dporfs("U", 2, 1, a, 2, af, 2, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("DPORFS", infot, nout, lerr, ok);
        //
        //        DPOCON
        //
        srnamt = "DPOCON";
        infot = 1;
        dpocon("/", 0, a, 1, anrm, rcond, w, iw, info);
        chkxer("DPOCON", infot, nout, lerr, ok);
        infot = 2;
        dpocon("U", -1, a, 1, anrm, rcond, w, iw, info);
        chkxer("DPOCON", infot, nout, lerr, ok);
        infot = 4;
        dpocon("U", 2, a, 1, anrm, rcond, w, iw, info);
        chkxer("DPOCON", infot, nout, lerr, ok);
        //
        //        DPOEQU
        //
        srnamt = "DPOEQU";
        infot = 1;
        dpoequ(-1, a, 1, r1, rcond, anrm, info);
        chkxer("DPOEQU", infot, nout, lerr, ok);
        infot = 3;
        dpoequ(2, a, 1, r1, rcond, anrm, info);
        chkxer("DPOEQU", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PP") {
        //
        //        Test error exits of the routines that use the Cholesky
        //        decomposition of a symmetric positive definite packed matrix.
        //
        //        DPPTRF
        //
        srnamt = "DPPTRF";
        infot = 1;
        dpptrf("/", 0, a, info);
        chkxer("DPPTRF", infot, nout, lerr, ok);
        infot = 2;
        dpptrf("U", -1, a, info);
        chkxer("DPPTRF", infot, nout, lerr, ok);
        //
        //        DPPTRI
        //
        srnamt = "DPPTRI";
        infot = 1;
        dpptri("/", 0, a, info);
        chkxer("DPPTRI", infot, nout, lerr, ok);
        infot = 2;
        dpptri("U", -1, a, info);
        chkxer("DPPTRI", infot, nout, lerr, ok);
        //
        //        DPPTRS
        //
        srnamt = "DPPTRS";
        infot = 1;
        dpptrs("/", 0, 0, a, b, 1, info);
        chkxer("DPPTRS", infot, nout, lerr, ok);
        infot = 2;
        dpptrs("U", -1, 0, a, b, 1, info);
        chkxer("DPPTRS", infot, nout, lerr, ok);
        infot = 3;
        dpptrs("U", 0, -1, a, b, 1, info);
        chkxer("DPPTRS", infot, nout, lerr, ok);
        infot = 6;
        dpptrs("U", 2, 1, a, b, 1, info);
        chkxer("DPPTRS", infot, nout, lerr, ok);
        //
        //        DPPRFS
        //
        srnamt = "DPPRFS";
        infot = 1;
        dpprfs("/", 0, 0, a, af, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DPPRFS", infot, nout, lerr, ok);
        infot = 2;
        dpprfs("U", -1, 0, a, af, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DPPRFS", infot, nout, lerr, ok);
        infot = 3;
        dpprfs("U", 0, -1, a, af, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DPPRFS", infot, nout, lerr, ok);
        infot = 7;
        dpprfs("U", 2, 1, a, af, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("DPPRFS", infot, nout, lerr, ok);
        infot = 9;
        dpprfs("U", 2, 1, a, af, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("DPPRFS", infot, nout, lerr, ok);
        //
        //        DPPCON
        //
        srnamt = "DPPCON";
        infot = 1;
        dppcon("/", 0, a, anrm, rcond, w, iw, info);
        chkxer("DPPCON", infot, nout, lerr, ok);
        infot = 2;
        dppcon("U", -1, a, anrm, rcond, w, iw, info);
        chkxer("DPPCON", infot, nout, lerr, ok);
        //
        //        DPPEQU
        //
        srnamt = "DPPEQU";
        infot = 1;
        dppequ("/", 0, a, r1, rcond, anrm, info);
        chkxer("DPPEQU", infot, nout, lerr, ok);
        infot = 2;
        dppequ("U", -1, a, r1, rcond, anrm, info);
        chkxer("DPPEQU", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PB") {
        //
        //        Test error exits of the routines that use the Cholesky
        //        decomposition of a symmetric positive definite band matrix.
        //
        //        DPBTRF
        //
        srnamt = "DPBTRF";
        infot = 1;
        dpbtrf("/", 0, 0, a, 1, info);
        chkxer("DPBTRF", infot, nout, lerr, ok);
        infot = 2;
        dpbtrf("U", -1, 0, a, 1, info);
        chkxer("DPBTRF", infot, nout, lerr, ok);
        infot = 3;
        dpbtrf("U", 1, -1, a, 1, info);
        chkxer("DPBTRF", infot, nout, lerr, ok);
        infot = 5;
        dpbtrf("U", 2, 1, a, 1, info);
        chkxer("DPBTRF", infot, nout, lerr, ok);
        //
        //        DPBTF2
        //
        srnamt = "DPBTF2";
        infot = 1;
        dpbtf2("/", 0, 0, a, 1, info);
        chkxer("DPBTF2", infot, nout, lerr, ok);
        infot = 2;
        dpbtf2("U", -1, 0, a, 1, info);
        chkxer("DPBTF2", infot, nout, lerr, ok);
        infot = 3;
        dpbtf2("U", 1, -1, a, 1, info);
        chkxer("DPBTF2", infot, nout, lerr, ok);
        infot = 5;
        dpbtf2("U", 2, 1, a, 1, info);
        chkxer("DPBTF2", infot, nout, lerr, ok);
        //
        //        DPBTRS
        //
        srnamt = "DPBTRS";
        infot = 1;
        dpbtrs("/", 0, 0, 0, a, 1, b, 1, info);
        chkxer("DPBTRS", infot, nout, lerr, ok);
        infot = 2;
        dpbtrs("U", -1, 0, 0, a, 1, b, 1, info);
        chkxer("DPBTRS", infot, nout, lerr, ok);
        infot = 3;
        dpbtrs("U", 1, -1, 0, a, 1, b, 1, info);
        chkxer("DPBTRS", infot, nout, lerr, ok);
        infot = 4;
        dpbtrs("U", 0, 0, -1, a, 1, b, 1, info);
        chkxer("DPBTRS", infot, nout, lerr, ok);
        infot = 6;
        dpbtrs("U", 2, 1, 1, a, 1, b, 1, info);
        chkxer("DPBTRS", infot, nout, lerr, ok);
        infot = 8;
        dpbtrs("U", 2, 0, 1, a, 1, b, 1, info);
        chkxer("DPBTRS", infot, nout, lerr, ok);
        //
        //        DPBRFS
        //
        srnamt = "DPBRFS";
        infot = 1;
        dpbrfs("/", 0, 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DPBRFS", infot, nout, lerr, ok);
        infot = 2;
        dpbrfs("U", -1, 0, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DPBRFS", infot, nout, lerr, ok);
        infot = 3;
        dpbrfs("U", 1, -1, 0, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DPBRFS", infot, nout, lerr, ok);
        infot = 4;
        dpbrfs("U", 0, 0, -1, a, 1, af, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DPBRFS", infot, nout, lerr, ok);
        infot = 6;
        dpbrfs("U", 2, 1, 1, a, 1, af, 2, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("DPBRFS", infot, nout, lerr, ok);
        infot = 8;
        dpbrfs("U", 2, 1, 1, a, 2, af, 1, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("DPBRFS", infot, nout, lerr, ok);
        infot = 10;
        dpbrfs("U", 2, 0, 1, a, 1, af, 1, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("DPBRFS", infot, nout, lerr, ok);
        infot = 12;
        dpbrfs("U", 2, 0, 1, a, 1, af, 1, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("DPBRFS", infot, nout, lerr, ok);
        //
        //        DPBCON
        //
        srnamt = "DPBCON";
        infot = 1;
        dpbcon("/", 0, 0, a, 1, anrm, rcond, w, iw, info);
        chkxer("DPBCON", infot, nout, lerr, ok);
        infot = 2;
        dpbcon("U", -1, 0, a, 1, anrm, rcond, w, iw, info);
        chkxer("DPBCON", infot, nout, lerr, ok);
        infot = 3;
        dpbcon("U", 1, -1, a, 1, anrm, rcond, w, iw, info);
        chkxer("DPBCON", infot, nout, lerr, ok);
        infot = 5;
        dpbcon("U", 2, 1, a, 1, anrm, rcond, w, iw, info);
        chkxer("DPBCON", infot, nout, lerr, ok);
        //
        //        DPBEQU
        //
        srnamt = "DPBEQU";
        infot = 1;
        dpbequ("/", 0, 0, a, 1, r1, rcond, anrm, info);
        chkxer("DPBEQU", infot, nout, lerr, ok);
        infot = 2;
        dpbequ("U", -1, 0, a, 1, r1, rcond, anrm, info);
        chkxer("DPBEQU", infot, nout, lerr, ok);
        infot = 3;
        dpbequ("U", 1, -1, a, 1, r1, rcond, anrm, info);
        chkxer("DPBEQU", infot, nout, lerr, ok);
        infot = 5;
        dpbequ("U", 2, 1, a, 1, r1, rcond, anrm, info);
        chkxer("DPBEQU", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrpo
    //
}
