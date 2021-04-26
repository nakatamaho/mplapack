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

void Rerrge(common &cmn, const char *path, INTEGER const nunit) {
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
    const INTEGER lw = 3 * nmax;
    arr_1d<lw, REAL> w(fill0);
    arr_1d<nmax, REAL> x(fill0);
    arr_1d<nmax, int> ip(fill0);
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
        ip[j - 1] = j;
        iw[j - 1] = j;
    }
    ok = true;
    //
    INTEGER info = 0;
    REAL anrm = 0.0;
    REAL rcond = 0.0;
    REAL ccond = 0.0;
    if (Mlsamen2, c2, "GE") {
        //
        //        Test error exits of the routines that use the LU decomposition
        //        of a general matrix.
        //
        //        DGETRF
        //
        srnamt = "DGETRF";
        infot = 1;
        dgetrf(-1, 0, a, 1, ip, info);
        chkxer("DGETRF", infot, nout, lerr, ok);
        infot = 2;
        dgetrf(0, -1, a, 1, ip, info);
        chkxer("DGETRF", infot, nout, lerr, ok);
        infot = 4;
        dgetrf(2, 1, a, 1, ip, info);
        chkxer("DGETRF", infot, nout, lerr, ok);
        //
        //        DGETF2
        //
        srnamt = "DGETF2";
        infot = 1;
        dgetf2(-1, 0, a, 1, ip, info);
        chkxer("DGETF2", infot, nout, lerr, ok);
        infot = 2;
        dgetf2(0, -1, a, 1, ip, info);
        chkxer("DGETF2", infot, nout, lerr, ok);
        infot = 4;
        dgetf2(2, 1, a, 1, ip, info);
        chkxer("DGETF2", infot, nout, lerr, ok);
        //
        //        DGETRI
        //
        srnamt = "DGETRI";
        infot = 1;
        dgetri(-1, a, 1, ip, w, lw, info);
        chkxer("DGETRI", infot, nout, lerr, ok);
        infot = 3;
        dgetri(2, a, 1, ip, w, lw, info);
        chkxer("DGETRI", infot, nout, lerr, ok);
        //
        //        DGETRS
        //
        srnamt = "DGETRS";
        infot = 1;
        dgetrs("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("DGETRS", infot, nout, lerr, ok);
        infot = 2;
        dgetrs("N", -1, 0, a, 1, ip, b, 1, info);
        chkxer("DGETRS", infot, nout, lerr, ok);
        infot = 3;
        dgetrs("N", 0, -1, a, 1, ip, b, 1, info);
        chkxer("DGETRS", infot, nout, lerr, ok);
        infot = 5;
        dgetrs("N", 2, 1, a, 1, ip, b, 2, info);
        chkxer("DGETRS", infot, nout, lerr, ok);
        infot = 8;
        dgetrs("N", 2, 1, a, 2, ip, b, 1, info);
        chkxer("DGETRS", infot, nout, lerr, ok);
        //
        //        RgerFS
        //
        srnamt = "RgerFS";
        infot = 1;
        Rgerfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("RgerFS", infot, nout, lerr, ok);
        infot = 2;
        Rgerfs("N", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("RgerFS", infot, nout, lerr, ok);
        infot = 3;
        Rgerfs("N", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("RgerFS", infot, nout, lerr, ok);
        infot = 5;
        Rgerfs("N", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("RgerFS", infot, nout, lerr, ok);
        infot = 7;
        Rgerfs("N", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("RgerFS", infot, nout, lerr, ok);
        infot = 10;
        Rgerfs("N", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("RgerFS", infot, nout, lerr, ok);
        infot = 12;
        Rgerfs("N", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("RgerFS", infot, nout, lerr, ok);
        //
        //        DGECON
        //
        srnamt = "DGECON";
        infot = 1;
        dgecon("/", 0, a, 1, anrm, rcond, w, iw, info);
        chkxer("DGECON", infot, nout, lerr, ok);
        infot = 2;
        dgecon("1", -1, a, 1, anrm, rcond, w, iw, info);
        chkxer("DGECON", infot, nout, lerr, ok);
        infot = 4;
        dgecon("1", 2, a, 1, anrm, rcond, w, iw, info);
        chkxer("DGECON", infot, nout, lerr, ok);
        //
        //        DGEEQU
        //
        srnamt = "DGEEQU";
        infot = 1;
        dgeequ(-1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("DGEEQU", infot, nout, lerr, ok);
        infot = 2;
        dgeequ(0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("DGEEQU", infot, nout, lerr, ok);
        infot = 4;
        dgeequ(2, 2, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("DGEEQU", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "GB") {
        //
        //        Test error exits of the routines that use the LU decomposition
        //        of a general band matrix.
        //
        //        DGBTRF
        //
        srnamt = "DGBTRF";
        infot = 1;
        dgbtrf(-1, 0, 0, 0, a, 1, ip, info);
        chkxer("DGBTRF", infot, nout, lerr, ok);
        infot = 2;
        dgbtrf(0, -1, 0, 0, a, 1, ip, info);
        chkxer("DGBTRF", infot, nout, lerr, ok);
        infot = 3;
        dgbtrf(1, 1, -1, 0, a, 1, ip, info);
        chkxer("DGBTRF", infot, nout, lerr, ok);
        infot = 4;
        dgbtrf(1, 1, 0, -1, a, 1, ip, info);
        chkxer("DGBTRF", infot, nout, lerr, ok);
        infot = 6;
        dgbtrf(2, 2, 1, 1, a, 3, ip, info);
        chkxer("DGBTRF", infot, nout, lerr, ok);
        //
        //        DGBTF2
        //
        srnamt = "DGBTF2";
        infot = 1;
        dgbtf2(-1, 0, 0, 0, a, 1, ip, info);
        chkxer("DGBTF2", infot, nout, lerr, ok);
        infot = 2;
        dgbtf2(0, -1, 0, 0, a, 1, ip, info);
        chkxer("DGBTF2", infot, nout, lerr, ok);
        infot = 3;
        dgbtf2(1, 1, -1, 0, a, 1, ip, info);
        chkxer("DGBTF2", infot, nout, lerr, ok);
        infot = 4;
        dgbtf2(1, 1, 0, -1, a, 1, ip, info);
        chkxer("DGBTF2", infot, nout, lerr, ok);
        infot = 6;
        dgbtf2(2, 2, 1, 1, a, 3, ip, info);
        chkxer("DGBTF2", infot, nout, lerr, ok);
        //
        //        DGBTRS
        //
        srnamt = "DGBTRS";
        infot = 1;
        dgbtrs("/", 0, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 2;
        dgbtrs("N", -1, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 3;
        dgbtrs("N", 1, -1, 0, 1, a, 1, ip, b, 1, info);
        chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 4;
        dgbtrs("N", 1, 0, -1, 1, a, 1, ip, b, 1, info);
        chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 5;
        dgbtrs("N", 1, 0, 0, -1, a, 1, ip, b, 1, info);
        chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 7;
        dgbtrs("N", 2, 1, 1, 1, a, 3, ip, b, 2, info);
        chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 10;
        dgbtrs("N", 2, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("DGBTRS", infot, nout, lerr, ok);
        //
        //        DGBRFS
        //
        srnamt = "DGBRFS";
        infot = 1;
        dgbrfs("/", 0, 0, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 2;
        dgbrfs("N", -1, 0, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 3;
        dgbrfs("N", 1, -1, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 4;
        dgbrfs("N", 1, 0, -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 5;
        dgbrfs("N", 1, 0, 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 7;
        dgbrfs("N", 2, 1, 1, 1, a, 2, af, 4, ip, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 9;
        dgbrfs("N", 2, 1, 1, 1, a, 3, af, 3, ip, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 12;
        dgbrfs("N", 2, 0, 0, 1, a, 1, af, 1, ip, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 14;
        dgbrfs("N", 2, 0, 0, 1, a, 1, af, 1, ip, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("DGBRFS", infot, nout, lerr, ok);
        //
        //        DGBCON
        //
        srnamt = "DGBCON";
        infot = 1;
        dgbcon("/", 0, 0, 0, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("DGBCON", infot, nout, lerr, ok);
        infot = 2;
        dgbcon("1", -1, 0, 0, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("DGBCON", infot, nout, lerr, ok);
        infot = 3;
        dgbcon("1", 1, -1, 0, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("DGBCON", infot, nout, lerr, ok);
        infot = 4;
        dgbcon("1", 1, 0, -1, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("DGBCON", infot, nout, lerr, ok);
        infot = 6;
        dgbcon("1", 2, 1, 1, a, 3, ip, anrm, rcond, w, iw, info);
        chkxer("DGBCON", infot, nout, lerr, ok);
        //
        //        DGBEQU
        //
        srnamt = "DGBEQU";
        infot = 1;
        dgbequ(-1, 0, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("DGBEQU", infot, nout, lerr, ok);
        infot = 2;
        dgbequ(0, -1, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("DGBEQU", infot, nout, lerr, ok);
        infot = 3;
        dgbequ(1, 1, -1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("DGBEQU", infot, nout, lerr, ok);
        infot = 4;
        dgbequ(1, 1, 0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("DGBEQU", infot, nout, lerr, ok);
        infot = 6;
        dgbequ(2, 2, 1, 1, a, 2, r1, r2, rcond, ccond, anrm, info);
        chkxer("DGBEQU", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrge
    //
}
