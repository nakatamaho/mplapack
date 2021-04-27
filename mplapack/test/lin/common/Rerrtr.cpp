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
    if (Mlsamen(2, c2, "TR")) {
        //
        //        Test error exits for the general triangular routines.
        //
        //        Rtrtri
        //
        srnamt = "Rtrtri";
        infot = 1;
        Rtrtri("/", "N", 0, a, 1, info);
        chkxer("Rtrtri", infot, nout, lerr, ok);
        infot = 2;
        Rtrtri("U", "/", 0, a, 1, info);
        chkxer("Rtrtri", infot, nout, lerr, ok);
        infot = 3;
        Rtrtri("U", "N", -1, a, 1, info);
        chkxer("Rtrtri", infot, nout, lerr, ok);
        infot = 5;
        Rtrtri("U", "N", 2, a, 1, info);
        chkxer("Rtrtri", infot, nout, lerr, ok);
        //
        //        Rtrti2
        //
        srnamt = "Rtrti2";
        infot = 1;
        Rtrti2("/", "N", 0, a, 1, info);
        chkxer("Rtrti2", infot, nout, lerr, ok);
        infot = 2;
        Rtrti2("U", "/", 0, a, 1, info);
        chkxer("Rtrti2", infot, nout, lerr, ok);
        infot = 3;
        Rtrti2("U", "N", -1, a, 1, info);
        chkxer("Rtrti2", infot, nout, lerr, ok);
        infot = 5;
        Rtrti2("U", "N", 2, a, 1, info);
        chkxer("Rtrti2", infot, nout, lerr, ok);
        //
        //        Rtrtrs
        //
        srnamt = "Rtrtrs";
        infot = 1;
        Rtrtrs("/", "N", "N", 0, 0, a, 1, x, 1, info);
        chkxer("Rtrtrs", infot, nout, lerr, ok);
        infot = 2;
        Rtrtrs("U", "/", "N", 0, 0, a, 1, x, 1, info);
        chkxer("Rtrtrs", infot, nout, lerr, ok);
        infot = 3;
        Rtrtrs("U", "N", "/", 0, 0, a, 1, x, 1, info);
        chkxer("Rtrtrs", infot, nout, lerr, ok);
        infot = 4;
        Rtrtrs("U", "N", "N", -1, 0, a, 1, x, 1, info);
        chkxer("Rtrtrs", infot, nout, lerr, ok);
        infot = 5;
        Rtrtrs("U", "N", "N", 0, -1, a, 1, x, 1, info);
        chkxer("Rtrtrs", infot, nout, lerr, ok);
        infot = 7;
        Rtrtrs("U", "N", "N", 2, 1, a, 1, x, 2, info);
        chkxer("Rtrtrs", infot, nout, lerr, ok);
        infot = 9;
        Rtrtrs("U", "N", "N", 2, 1, a, 2, x, 1, info);
        chkxer("Rtrtrs", infot, nout, lerr, ok);
        //
        //        Rtrrfs
        //
        srnamt = "Rtrrfs";
        infot = 1;
        Rtrrfs("/", "N", "N", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtrrfs", infot, nout, lerr, ok);
        infot = 2;
        Rtrrfs("U", "/", "N", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtrrfs", infot, nout, lerr, ok);
        infot = 3;
        Rtrrfs("U", "N", "/", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtrrfs", infot, nout, lerr, ok);
        infot = 4;
        Rtrrfs("U", "N", "N", -1, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtrrfs", infot, nout, lerr, ok);
        infot = 5;
        Rtrrfs("U", "N", "N", 0, -1, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtrrfs", infot, nout, lerr, ok);
        infot = 7;
        Rtrrfs("U", "N", "N", 2, 1, a, 1, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("Rtrrfs", infot, nout, lerr, ok);
        infot = 9;
        Rtrrfs("U", "N", "N", 2, 1, a, 2, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("Rtrrfs", infot, nout, lerr, ok);
        infot = 11;
        Rtrrfs("U", "N", "N", 2, 1, a, 2, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("Rtrrfs", infot, nout, lerr, ok);
        //
        //        Rtrcon
        //
        srnamt = "Rtrcon";
        infot = 1;
        Rtrcon("/", "U", "N", 0, a, 1, rcond, w, iw, info);
        chkxer("Rtrcon", infot, nout, lerr, ok);
        infot = 2;
        Rtrcon("1", "/", "N", 0, a, 1, rcond, w, iw, info);
        chkxer("Rtrcon", infot, nout, lerr, ok);
        infot = 3;
        Rtrcon("1", "U", "/", 0, a, 1, rcond, w, iw, info);
        chkxer("Rtrcon", infot, nout, lerr, ok);
        infot = 4;
        Rtrcon("1", "U", "N", -1, a, 1, rcond, w, iw, info);
        chkxer("Rtrcon", infot, nout, lerr, ok);
        infot = 6;
        Rtrcon("1", "U", "N", 2, a, 1, rcond, w, iw, info);
        chkxer("Rtrcon", infot, nout, lerr, ok);
        //
        //        Rlatrs
        //
        srnamt = "Rlatrs";
        infot = 1;
        Rlatrs("/", "N", "N", "N", 0, a, 1, x, scale, w, info);
        chkxer("Rlatrs", infot, nout, lerr, ok);
        infot = 2;
        Rlatrs("U", "/", "N", "N", 0, a, 1, x, scale, w, info);
        chkxer("Rlatrs", infot, nout, lerr, ok);
        infot = 3;
        Rlatrs("U", "N", "/", "N", 0, a, 1, x, scale, w, info);
        chkxer("Rlatrs", infot, nout, lerr, ok);
        infot = 4;
        Rlatrs("U", "N", "N", "/", 0, a, 1, x, scale, w, info);
        chkxer("Rlatrs", infot, nout, lerr, ok);
        infot = 5;
        Rlatrs("U", "N", "N", "N", -1, a, 1, x, scale, w, info);
        chkxer("Rlatrs", infot, nout, lerr, ok);
        infot = 7;
        Rlatrs("U", "N", "N", "N", 2, a, 1, x, scale, w, info);
        chkxer("Rlatrs", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "TP")) {
        //
        //        Test error exits for the packed triangular routines.
        //
        //        Rtptri
        //
        srnamt = "Rtptri";
        infot = 1;
        Rtptri("/", "N", 0, a, info);
        chkxer("Rtptri", infot, nout, lerr, ok);
        infot = 2;
        Rtptri("U", "/", 0, a, info);
        chkxer("Rtptri", infot, nout, lerr, ok);
        infot = 3;
        Rtptri("U", "N", -1, a, info);
        chkxer("Rtptri", infot, nout, lerr, ok);
        //
        //        Rtptrs
        //
        srnamt = "Rtptrs";
        infot = 1;
        Rtptrs("/", "N", "N", 0, 0, a, x, 1, info);
        chkxer("Rtptrs", infot, nout, lerr, ok);
        infot = 2;
        Rtptrs("U", "/", "N", 0, 0, a, x, 1, info);
        chkxer("Rtptrs", infot, nout, lerr, ok);
        infot = 3;
        Rtptrs("U", "N", "/", 0, 0, a, x, 1, info);
        chkxer("Rtptrs", infot, nout, lerr, ok);
        infot = 4;
        Rtptrs("U", "N", "N", -1, 0, a, x, 1, info);
        chkxer("Rtptrs", infot, nout, lerr, ok);
        infot = 5;
        Rtptrs("U", "N", "N", 0, -1, a, x, 1, info);
        chkxer("Rtptrs", infot, nout, lerr, ok);
        infot = 8;
        Rtptrs("U", "N", "N", 2, 1, a, x, 1, info);
        chkxer("Rtptrs", infot, nout, lerr, ok);
        //
        //        Rtprfs
        //
        srnamt = "Rtprfs";
        infot = 1;
        Rtprfs("/", "N", "N", 0, 0, a, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtprfs", infot, nout, lerr, ok);
        infot = 2;
        Rtprfs("U", "/", "N", 0, 0, a, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtprfs", infot, nout, lerr, ok);
        infot = 3;
        Rtprfs("U", "N", "/", 0, 0, a, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtprfs", infot, nout, lerr, ok);
        infot = 4;
        Rtprfs("U", "N", "N", -1, 0, a, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtprfs", infot, nout, lerr, ok);
        infot = 5;
        Rtprfs("U", "N", "N", 0, -1, a, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtprfs", infot, nout, lerr, ok);
        infot = 8;
        Rtprfs("U", "N", "N", 2, 1, a, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("Rtprfs", infot, nout, lerr, ok);
        infot = 10;
        Rtprfs("U", "N", "N", 2, 1, a, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("Rtprfs", infot, nout, lerr, ok);
        //
        //        Rtpcon
        //
        srnamt = "Rtpcon";
        infot = 1;
        Rtpcon("/", "U", "N", 0, a, rcond, w, iw, info);
        chkxer("Rtpcon", infot, nout, lerr, ok);
        infot = 2;
        Rtpcon("1", "/", "N", 0, a, rcond, w, iw, info);
        chkxer("Rtpcon", infot, nout, lerr, ok);
        infot = 3;
        Rtpcon("1", "U", "/", 0, a, rcond, w, iw, info);
        chkxer("Rtpcon", infot, nout, lerr, ok);
        infot = 4;
        Rtpcon("1", "U", "N", -1, a, rcond, w, iw, info);
        chkxer("Rtpcon", infot, nout, lerr, ok);
        //
        //        Rlatps
        //
        srnamt = "Rlatps";
        infot = 1;
        Rlatps("/", "N", "N", "N", 0, a, x, scale, w, info);
        chkxer("Rlatps", infot, nout, lerr, ok);
        infot = 2;
        Rlatps("U", "/", "N", "N", 0, a, x, scale, w, info);
        chkxer("Rlatps", infot, nout, lerr, ok);
        infot = 3;
        Rlatps("U", "N", "/", "N", 0, a, x, scale, w, info);
        chkxer("Rlatps", infot, nout, lerr, ok);
        infot = 4;
        Rlatps("U", "N", "N", "/", 0, a, x, scale, w, info);
        chkxer("Rlatps", infot, nout, lerr, ok);
        infot = 5;
        Rlatps("U", "N", "N", "N", -1, a, x, scale, w, info);
        chkxer("Rlatps", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "TB")) {
        //
        //        Test error exits for the banded triangular routines.
        //
        //        Rtbtrs
        //
        srnamt = "Rtbtrs";
        infot = 1;
        Rtbtrs("/", "N", "N", 0, 0, 0, a, 1, x, 1, info);
        chkxer("Rtbtrs", infot, nout, lerr, ok);
        infot = 2;
        Rtbtrs("U", "/", "N", 0, 0, 0, a, 1, x, 1, info);
        chkxer("Rtbtrs", infot, nout, lerr, ok);
        infot = 3;
        Rtbtrs("U", "N", "/", 0, 0, 0, a, 1, x, 1, info);
        chkxer("Rtbtrs", infot, nout, lerr, ok);
        infot = 4;
        Rtbtrs("U", "N", "N", -1, 0, 0, a, 1, x, 1, info);
        chkxer("Rtbtrs", infot, nout, lerr, ok);
        infot = 5;
        Rtbtrs("U", "N", "N", 0, -1, 0, a, 1, x, 1, info);
        chkxer("Rtbtrs", infot, nout, lerr, ok);
        infot = 6;
        Rtbtrs("U", "N", "N", 0, 0, -1, a, 1, x, 1, info);
        chkxer("Rtbtrs", infot, nout, lerr, ok);
        infot = 8;
        Rtbtrs("U", "N", "N", 2, 1, 1, a, 1, x, 2, info);
        chkxer("Rtbtrs", infot, nout, lerr, ok);
        infot = 10;
        Rtbtrs("U", "N", "N", 2, 0, 1, a, 1, x, 1, info);
        chkxer("Rtbtrs", infot, nout, lerr, ok);
        //
        //        Rtbrfs
        //
        srnamt = "Rtbrfs";
        infot = 1;
        Rtbrfs("/", "N", "N", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtbrfs", infot, nout, lerr, ok);
        infot = 2;
        Rtbrfs("U", "/", "N", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtbrfs", infot, nout, lerr, ok);
        infot = 3;
        Rtbrfs("U", "N", "/", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtbrfs", infot, nout, lerr, ok);
        infot = 4;
        Rtbrfs("U", "N", "N", -1, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtbrfs", infot, nout, lerr, ok);
        infot = 5;
        Rtbrfs("U", "N", "N", 0, -1, 0, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtbrfs", infot, nout, lerr, ok);
        infot = 6;
        Rtbrfs("U", "N", "N", 0, 0, -1, a, 1, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rtbrfs", infot, nout, lerr, ok);
        infot = 8;
        Rtbrfs("U", "N", "N", 2, 1, 1, a, 1, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("Rtbrfs", infot, nout, lerr, ok);
        infot = 10;
        Rtbrfs("U", "N", "N", 2, 1, 1, a, 2, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("Rtbrfs", infot, nout, lerr, ok);
        infot = 12;
        Rtbrfs("U", "N", "N", 2, 1, 1, a, 2, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("Rtbrfs", infot, nout, lerr, ok);
        //
        //        Rtbcon
        //
        srnamt = "Rtbcon";
        infot = 1;
        Rtbcon("/", "U", "N", 0, 0, a, 1, rcond, w, iw, info);
        chkxer("Rtbcon", infot, nout, lerr, ok);
        infot = 2;
        Rtbcon("1", "/", "N", 0, 0, a, 1, rcond, w, iw, info);
        chkxer("Rtbcon", infot, nout, lerr, ok);
        infot = 3;
        Rtbcon("1", "U", "/", 0, 0, a, 1, rcond, w, iw, info);
        chkxer("Rtbcon", infot, nout, lerr, ok);
        infot = 4;
        Rtbcon("1", "U", "N", -1, 0, a, 1, rcond, w, iw, info);
        chkxer("Rtbcon", infot, nout, lerr, ok);
        infot = 5;
        Rtbcon("1", "U", "N", 0, -1, a, 1, rcond, w, iw, info);
        chkxer("Rtbcon", infot, nout, lerr, ok);
        infot = 7;
        Rtbcon("1", "U", "N", 2, 1, a, 1, rcond, w, iw, info);
        chkxer("Rtbcon", infot, nout, lerr, ok);
        //
        //        Rlatbs
        //
        srnamt = "Rlatbs";
        infot = 1;
        Rlatbs("/", "N", "N", "N", 0, 0, a, 1, x, scale, w, info);
        chkxer("Rlatbs", infot, nout, lerr, ok);
        infot = 2;
        Rlatbs("U", "/", "N", "N", 0, 0, a, 1, x, scale, w, info);
        chkxer("Rlatbs", infot, nout, lerr, ok);
        infot = 3;
        Rlatbs("U", "N", "/", "N", 0, 0, a, 1, x, scale, w, info);
        chkxer("Rlatbs", infot, nout, lerr, ok);
        infot = 4;
        Rlatbs("U", "N", "N", "/", 0, 0, a, 1, x, scale, w, info);
        chkxer("Rlatbs", infot, nout, lerr, ok);
        infot = 5;
        Rlatbs("U", "N", "N", "N", -1, 0, a, 1, x, scale, w, info);
        chkxer("Rlatbs", infot, nout, lerr, ok);
        infot = 6;
        Rlatbs("U", "N", "N", "N", 1, -1, a, 1, x, scale, w, info);
        chkxer("Rlatbs", infot, nout, lerr, ok);
        infot = 8;
        Rlatbs("U", "N", "N", "N", 2, 1, a, 1, x, scale, w, info);
        chkxer("Rlatbs", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrtr
    //
}
