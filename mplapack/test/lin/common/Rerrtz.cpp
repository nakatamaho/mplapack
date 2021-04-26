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

void Rerrtz(common &cmn, const char *path, INTEGER const nunit) {
    common_write write(cmn);
    // COMMON infoc
    INTEGER &infot = cmn.infot;
    INTEGER &nout = cmn.nout;
    bool &ok = cmn.ok;
    bool &lerr = cmn.lerr;
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
    a[(1 - 1)] = 1.e+0;
    a[(2 - 1) * lda] = 2.e+0;
    a[(2 - 1) + (2 - 1) * lda] = 3.e+0;
    a[(2 - 1)] = 4.e+0;
    arr_1d<nmax, REAL> w(fill0);
    w[1 - 1] = 0.0;
    w[2 - 1] = 0.0;
    ok = true;
    //
    arr_1d<nmax, REAL> tau(fill0);
    INTEGER info = 0;
    if (Mlsamen2, c2, "TZ") {
        //
        //        Test error exits for the trapezoidal routines.
        //
        //        DTZRZF
        //
        cmn.srnamt = "DTZRZF";
        infot = 1;
        dtzrzf(-1, 0, a, 1, tau, w, 1, info);
        chkxer("DTZRZF", infot, nout, lerr, ok);
        infot = 2;
        dtzrzf(1, 0, a, 1, tau, w, 1, info);
        chkxer("DTZRZF", infot, nout, lerr, ok);
        infot = 4;
        dtzrzf(2, 2, a, 1, tau, w, 1, info);
        chkxer("DTZRZF", infot, nout, lerr, ok);
        infot = 7;
        dtzrzf(2, 2, a, 2, tau, w, 0, info);
        chkxer("DTZRZF", infot, nout, lerr, ok);
        infot = 7;
        dtzrzf(2, 3, a, 2, tau, w, 1, info);
        chkxer("DTZRZF", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrtz
    //
}
