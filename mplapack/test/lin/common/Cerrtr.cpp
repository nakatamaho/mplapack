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
#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;
#include <mplapack_lin.h>
#include <mplapack.h>

void Cerrtr(const char *path, INTEGER const nunit) {
    common_write write(cmn);
    // COMMON infoc
    INTEGER &infot = cmn.infot;
    INTEGER &nout = cmn.nout;
    bool &ok = cmn.ok;
    bool &lerr = cmn.lerr;
    // COMMON srnamc
    char[32] &srnamt = cmn.srnamt;
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
    char[2] c2 = path[(2 - 1) + (3 - 1) * ldpath];
    const INTEGER nmax = 2;
    COMPLEX a[nmax * nmax];
    a[(1 - 1)] = 1.0;
    a[(2 - 1) * lda] = 2.0;
    a[(2 - 1) + (2 - 1) * lda] = 3.e0;
    a[(2 - 1)] = 4.e0;
    ok = true;
    //
    //     Test error exits for the general triangular routines.
    //
    INTEGER info = 0;
    COMPLEX x[nmax];
    COMPLEX b[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    COMPLEX w[nmax];
    REAL rw[nmax];
    REAL rcond = 0.0;
    REAL scale = 0.0;
    if (Mlsamen(2, c2, "TR")) {
        //
        //        Ctrtri
        //
        srnamt = "Ctrtri";
        infot = 1;
        Ctrtri("/", "N", 0, a, 1, info);
        chkxer("Ctrtri", infot, nout, lerr, ok);
        infot = 2;
        Ctrtri("U", "/", 0, a, 1, info);
        chkxer("Ctrtri", infot, nout, lerr, ok);
        infot = 3;
        Ctrtri("U", "N", -1, a, 1, info);
        chkxer("Ctrtri", infot, nout, lerr, ok);
        infot = 5;
        Ctrtri("U", "N", 2, a, 1, info);
        chkxer("Ctrtri", infot, nout, lerr, ok);
        //
        //        Ctrti2
        //
        srnamt = "Ctrti2";
        infot = 1;
        Ctrti2("/", "N", 0, a, 1, info);
        chkxer("Ctrti2", infot, nout, lerr, ok);
        infot = 2;
        Ctrti2("U", "/", 0, a, 1, info);
        chkxer("Ctrti2", infot, nout, lerr, ok);
        infot = 3;
        Ctrti2("U", "N", -1, a, 1, info);
        chkxer("Ctrti2", infot, nout, lerr, ok);
        infot = 5;
        Ctrti2("U", "N", 2, a, 1, info);
        chkxer("Ctrti2", infot, nout, lerr, ok);
        //
        //        Ctrtrs
        //
        srnamt = "Ctrtrs";
        infot = 1;
        Ctrtrs("/", "N", "N", 0, 0, a, 1, x, 1, info);
        chkxer("Ctrtrs", infot, nout, lerr, ok);
        infot = 2;
        Ctrtrs("U", "/", "N", 0, 0, a, 1, x, 1, info);
        chkxer("Ctrtrs", infot, nout, lerr, ok);
        infot = 3;
        Ctrtrs("U", "N", "/", 0, 0, a, 1, x, 1, info);
        chkxer("Ctrtrs", infot, nout, lerr, ok);
        infot = 4;
        Ctrtrs("U", "N", "N", -1, 0, a, 1, x, 1, info);
        chkxer("Ctrtrs", infot, nout, lerr, ok);
        infot = 5;
        Ctrtrs("U", "N", "N", 0, -1, a, 1, x, 1, info);
        chkxer("Ctrtrs", infot, nout, lerr, ok);
        infot = 7;
        //
        //        Ctrrfs
        //
        srnamt = "Ctrrfs";
        infot = 1;
        Ctrrfs("/", "N", "N", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctrrfs", infot, nout, lerr, ok);
        infot = 2;
        Ctrrfs("U", "/", "N", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctrrfs", infot, nout, lerr, ok);
        infot = 3;
        Ctrrfs("U", "N", "/", 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctrrfs", infot, nout, lerr, ok);
        infot = 4;
        Ctrrfs("U", "N", "N", -1, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctrrfs", infot, nout, lerr, ok);
        infot = 5;
        Ctrrfs("U", "N", "N", 0, -1, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctrrfs", infot, nout, lerr, ok);
        infot = 7;
        Ctrrfs("U", "N", "N", 2, 1, a, 1, b, 2, x, 2, r1, r2, w, rw, info);
        chkxer("Ctrrfs", infot, nout, lerr, ok);
        infot = 9;
        Ctrrfs("U", "N", "N", 2, 1, a, 2, b, 1, x, 2, r1, r2, w, rw, info);
        chkxer("Ctrrfs", infot, nout, lerr, ok);
        infot = 11;
        Ctrrfs("U", "N", "N", 2, 1, a, 2, b, 2, x, 1, r1, r2, w, rw, info);
        chkxer("Ctrrfs", infot, nout, lerr, ok);
        //
        //        Ctrcon
        //
        srnamt = "Ctrcon";
        infot = 1;
        Ctrcon("/", "U", "N", 0, a, 1, rcond, w, rw, info);
        chkxer("Ctrcon", infot, nout, lerr, ok);
        infot = 2;
        Ctrcon("1", "/", "N", 0, a, 1, rcond, w, rw, info);
        chkxer("Ctrcon", infot, nout, lerr, ok);
        infot = 3;
        Ctrcon("1", "U", "/", 0, a, 1, rcond, w, rw, info);
        chkxer("Ctrcon", infot, nout, lerr, ok);
        infot = 4;
        Ctrcon("1", "U", "N", -1, a, 1, rcond, w, rw, info);
        chkxer("Ctrcon", infot, nout, lerr, ok);
        infot = 6;
        Ctrcon("1", "U", "N", 2, a, 1, rcond, w, rw, info);
        chkxer("Ctrcon", infot, nout, lerr, ok);
        //
        //        Clatrs
        //
        srnamt = "Clatrs";
        infot = 1;
        Clatrs("/", "N", "N", "N", 0, a, 1, x, scale, rw, info);
        chkxer("Clatrs", infot, nout, lerr, ok);
        infot = 2;
        Clatrs("U", "/", "N", "N", 0, a, 1, x, scale, rw, info);
        chkxer("Clatrs", infot, nout, lerr, ok);
        infot = 3;
        Clatrs("U", "N", "/", "N", 0, a, 1, x, scale, rw, info);
        chkxer("Clatrs", infot, nout, lerr, ok);
        infot = 4;
        Clatrs("U", "N", "N", "/", 0, a, 1, x, scale, rw, info);
        chkxer("Clatrs", infot, nout, lerr, ok);
        infot = 5;
        Clatrs("U", "N", "N", "N", -1, a, 1, x, scale, rw, info);
        chkxer("Clatrs", infot, nout, lerr, ok);
        infot = 7;
        Clatrs("U", "N", "N", "N", 2, a, 1, x, scale, rw, info);
        chkxer("Clatrs", infot, nout, lerr, ok);
        //
        //     Test error exits for the packed triangular routines.
        //
    } else if (Mlsamen(2, c2, "TP")) {
        //
        //        Ctptri
        //
        srnamt = "Ctptri";
        infot = 1;
        Ctptri("/", "N", 0, a, info);
        chkxer("Ctptri", infot, nout, lerr, ok);
        infot = 2;
        Ctptri("U", "/", 0, a, info);
        chkxer("Ctptri", infot, nout, lerr, ok);
        infot = 3;
        Ctptri("U", "N", -1, a, info);
        chkxer("Ctptri", infot, nout, lerr, ok);
        //
        //        Ctptrs
        //
        srnamt = "Ctptrs";
        infot = 1;
        Ctptrs("/", "N", "N", 0, 0, a, x, 1, info);
        chkxer("Ctptrs", infot, nout, lerr, ok);
        infot = 2;
        Ctptrs("U", "/", "N", 0, 0, a, x, 1, info);
        chkxer("Ctptrs", infot, nout, lerr, ok);
        infot = 3;
        Ctptrs("U", "N", "/", 0, 0, a, x, 1, info);
        chkxer("Ctptrs", infot, nout, lerr, ok);
        infot = 4;
        Ctptrs("U", "N", "N", -1, 0, a, x, 1, info);
        chkxer("Ctptrs", infot, nout, lerr, ok);
        infot = 5;
        Ctptrs("U", "N", "N", 0, -1, a, x, 1, info);
        chkxer("Ctptrs", infot, nout, lerr, ok);
        infot = 8;
        Ctptrs("U", "N", "N", 2, 1, a, x, 1, info);
        chkxer("Ctptrs", infot, nout, lerr, ok);
        //
        //        Ctprfs
        //
        srnamt = "Ctprfs";
        infot = 1;
        Ctprfs("/", "N", "N", 0, 0, a, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctprfs", infot, nout, lerr, ok);
        infot = 2;
        Ctprfs("U", "/", "N", 0, 0, a, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctprfs", infot, nout, lerr, ok);
        infot = 3;
        Ctprfs("U", "N", "/", 0, 0, a, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctprfs", infot, nout, lerr, ok);
        infot = 4;
        Ctprfs("U", "N", "N", -1, 0, a, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctprfs", infot, nout, lerr, ok);
        infot = 5;
        Ctprfs("U", "N", "N", 0, -1, a, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctprfs", infot, nout, lerr, ok);
        infot = 8;
        Ctprfs("U", "N", "N", 2, 1, a, b, 1, x, 2, r1, r2, w, rw, info);
        chkxer("Ctprfs", infot, nout, lerr, ok);
        infot = 10;
        Ctprfs("U", "N", "N", 2, 1, a, b, 2, x, 1, r1, r2, w, rw, info);
        chkxer("Ctprfs", infot, nout, lerr, ok);
        //
        //        Ctpcon
        //
        srnamt = "Ctpcon";
        infot = 1;
        Ctpcon("/", "U", "N", 0, a, rcond, w, rw, info);
        chkxer("Ctpcon", infot, nout, lerr, ok);
        infot = 2;
        Ctpcon("1", "/", "N", 0, a, rcond, w, rw, info);
        chkxer("Ctpcon", infot, nout, lerr, ok);
        infot = 3;
        Ctpcon("1", "U", "/", 0, a, rcond, w, rw, info);
        chkxer("Ctpcon", infot, nout, lerr, ok);
        infot = 4;
        Ctpcon("1", "U", "N", -1, a, rcond, w, rw, info);
        chkxer("Ctpcon", infot, nout, lerr, ok);
        //
        //        Clatps
        //
        srnamt = "Clatps";
        infot = 1;
        Clatps("/", "N", "N", "N", 0, a, x, scale, rw, info);
        chkxer("Clatps", infot, nout, lerr, ok);
        infot = 2;
        Clatps("U", "/", "N", "N", 0, a, x, scale, rw, info);
        chkxer("Clatps", infot, nout, lerr, ok);
        infot = 3;
        Clatps("U", "N", "/", "N", 0, a, x, scale, rw, info);
        chkxer("Clatps", infot, nout, lerr, ok);
        infot = 4;
        Clatps("U", "N", "N", "/", 0, a, x, scale, rw, info);
        chkxer("Clatps", infot, nout, lerr, ok);
        infot = 5;
        Clatps("U", "N", "N", "N", -1, a, x, scale, rw, info);
        chkxer("Clatps", infot, nout, lerr, ok);
        //
        //     Test error exits for the banded triangular routines.
        //
    } else if (Mlsamen(2, c2, "TB")) {
        //
        //        Ctbtrs
        //
        srnamt = "Ctbtrs";
        infot = 1;
        Ctbtrs("/", "N", "N", 0, 0, 0, a, 1, x, 1, info);
        chkxer("Ctbtrs", infot, nout, lerr, ok);
        infot = 2;
        Ctbtrs("U", "/", "N", 0, 0, 0, a, 1, x, 1, info);
        chkxer("Ctbtrs", infot, nout, lerr, ok);
        infot = 3;
        Ctbtrs("U", "N", "/", 0, 0, 0, a, 1, x, 1, info);
        chkxer("Ctbtrs", infot, nout, lerr, ok);
        infot = 4;
        Ctbtrs("U", "N", "N", -1, 0, 0, a, 1, x, 1, info);
        chkxer("Ctbtrs", infot, nout, lerr, ok);
        infot = 5;
        Ctbtrs("U", "N", "N", 0, -1, 0, a, 1, x, 1, info);
        chkxer("Ctbtrs", infot, nout, lerr, ok);
        infot = 6;
        Ctbtrs("U", "N", "N", 0, 0, -1, a, 1, x, 1, info);
        chkxer("Ctbtrs", infot, nout, lerr, ok);
        infot = 8;
        Ctbtrs("U", "N", "N", 2, 1, 1, a, 1, x, 2, info);
        chkxer("Ctbtrs", infot, nout, lerr, ok);
        infot = 10;
        Ctbtrs("U", "N", "N", 2, 0, 1, a, 1, x, 1, info);
        chkxer("Ctbtrs", infot, nout, lerr, ok);
        //
        //        Ctbrfs
        //
        srnamt = "Ctbrfs";
        infot = 1;
        Ctbrfs("/", "N", "N", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctbrfs", infot, nout, lerr, ok);
        infot = 2;
        Ctbrfs("U", "/", "N", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctbrfs", infot, nout, lerr, ok);
        infot = 3;
        Ctbrfs("U", "N", "/", 0, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctbrfs", infot, nout, lerr, ok);
        infot = 4;
        Ctbrfs("U", "N", "N", -1, 0, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctbrfs", infot, nout, lerr, ok);
        infot = 5;
        Ctbrfs("U", "N", "N", 0, -1, 0, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctbrfs", infot, nout, lerr, ok);
        infot = 6;
        Ctbrfs("U", "N", "N", 0, 0, -1, a, 1, b, 1, x, 1, r1, r2, w, rw, info);
        chkxer("Ctbrfs", infot, nout, lerr, ok);
        infot = 8;
        Ctbrfs("U", "N", "N", 2, 1, 1, a, 1, b, 2, x, 2, r1, r2, w, rw, info);
        chkxer("Ctbrfs", infot, nout, lerr, ok);
        infot = 10;
        Ctbrfs("U", "N", "N", 2, 1, 1, a, 2, b, 1, x, 2, r1, r2, w, rw, info);
        chkxer("Ctbrfs", infot, nout, lerr, ok);
        infot = 12;
        Ctbrfs("U", "N", "N", 2, 1, 1, a, 2, b, 2, x, 1, r1, r2, w, rw, info);
        chkxer("Ctbrfs", infot, nout, lerr, ok);
        //
        //        Ctbcon
        //
        srnamt = "Ctbcon";
        infot = 1;
        Ctbcon("/", "U", "N", 0, 0, a, 1, rcond, w, rw, info);
        chkxer("Ctbcon", infot, nout, lerr, ok);
        infot = 2;
        Ctbcon("1", "/", "N", 0, 0, a, 1, rcond, w, rw, info);
        chkxer("Ctbcon", infot, nout, lerr, ok);
        infot = 3;
        Ctbcon("1", "U", "/", 0, 0, a, 1, rcond, w, rw, info);
        chkxer("Ctbcon", infot, nout, lerr, ok);
        infot = 4;
        Ctbcon("1", "U", "N", -1, 0, a, 1, rcond, w, rw, info);
        chkxer("Ctbcon", infot, nout, lerr, ok);
        infot = 5;
        Ctbcon("1", "U", "N", 0, -1, a, 1, rcond, w, rw, info);
        chkxer("Ctbcon", infot, nout, lerr, ok);
        infot = 7;
        Ctbcon("1", "U", "N", 2, 1, a, 1, rcond, w, rw, info);
        chkxer("Ctbcon", infot, nout, lerr, ok);
        //
        //        Clatbs
        //
        srnamt = "Clatbs";
        infot = 1;
        Clatbs("/", "N", "N", "N", 0, 0, a, 1, x, scale, rw, info);
        chkxer("Clatbs", infot, nout, lerr, ok);
        infot = 2;
        Clatbs("U", "/", "N", "N", 0, 0, a, 1, x, scale, rw, info);
        chkxer("Clatbs", infot, nout, lerr, ok);
        infot = 3;
        Clatbs("U", "N", "/", "N", 0, 0, a, 1, x, scale, rw, info);
        chkxer("Clatbs", infot, nout, lerr, ok);
        infot = 4;
        Clatbs("U", "N", "N", "/", 0, 0, a, 1, x, scale, rw, info);
        chkxer("Clatbs", infot, nout, lerr, ok);
        infot = 5;
        Clatbs("U", "N", "N", "N", -1, 0, a, 1, x, scale, rw, info);
        chkxer("Clatbs", infot, nout, lerr, ok);
        infot = 6;
        Clatbs("U", "N", "N", "N", 1, -1, a, 1, x, scale, rw, info);
        chkxer("Clatbs", infot, nout, lerr, ok);
        infot = 8;
        Clatbs("U", "N", "N", "N", 2, 1, a, 1, x, scale, rw, info);
        chkxer("Clatbs", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrtr
    //
}
