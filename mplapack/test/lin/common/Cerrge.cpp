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

void Cerrge(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_1d<nmax, int> ip(fill0);
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
        ip[j - 1] = j;
    }
    ok = true;
    //
    //     Test error exits of the routines that use the LU decomposition
    //     of a general matrix.
    //
    INTEGER info = 0;
    arr_1d<nmax, REAL> r(fill0);
    REAL anrm = 0.0;
    REAL rcond = 0.0;
    REAL ccond = 0.0;
    if (Mlsamen2, c2, "GE") {
        //
        //        ZGETRF
        //
        srnamt = "ZGETRF";
        infot = 1;
        zgetrf(-1, 0, a, 1, ip, info);
        chkxer("ZGETRF", infot, nout, lerr, ok);
        infot = 2;
        zgetrf(0, -1, a, 1, ip, info);
        chkxer("ZGETRF", infot, nout, lerr, ok);
        infot = 4;
        zgetrf(2, 1, a, 1, ip, info);
        chkxer("ZGETRF", infot, nout, lerr, ok);
        //
        //        ZGETF2
        //
        srnamt = "ZGETF2";
        infot = 1;
        zgetf2(-1, 0, a, 1, ip, info);
        chkxer("ZGETF2", infot, nout, lerr, ok);
        infot = 2;
        zgetf2(0, -1, a, 1, ip, info);
        chkxer("ZGETF2", infot, nout, lerr, ok);
        infot = 4;
        zgetf2(2, 1, a, 1, ip, info);
        chkxer("ZGETF2", infot, nout, lerr, ok);
        //
        //        ZGETRI
        //
        srnamt = "ZGETRI";
        infot = 1;
        zgetri(-1, a, 1, ip, w, 1, info);
        chkxer("ZGETRI", infot, nout, lerr, ok);
        infot = 3;
        zgetri(2, a, 1, ip, w, 2, info);
        chkxer("ZGETRI", infot, nout, lerr, ok);
        infot = 6;
        zgetri(2, a, 2, ip, w, 1, info);
        chkxer("ZGETRI", infot, nout, lerr, ok);
        //
        //        ZGETRS
        //
        srnamt = "ZGETRS";
        infot = 1;
        zgetrs("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("ZGETRS", infot, nout, lerr, ok);
        infot = 2;
        zgetrs("N", -1, 0, a, 1, ip, b, 1, info);
        chkxer("ZGETRS", infot, nout, lerr, ok);
        infot = 3;
        zgetrs("N", 0, -1, a, 1, ip, b, 1, info);
        chkxer("ZGETRS", infot, nout, lerr, ok);
        infot = 5;
        zgetrs("N", 2, 1, a, 1, ip, b, 2, info);
        chkxer("ZGETRS", infot, nout, lerr, ok);
        infot = 8;
        zgetrs("N", 2, 1, a, 2, ip, b, 1, info);
        chkxer("ZGETRS", infot, nout, lerr, ok);
        //
        //        ZGERFS
        //
        srnamt = "ZGERFS";
        infot = 1;
        zgerfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZGERFS", infot, nout, lerr, ok);
        infot = 2;
        zgerfs("N", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZGERFS", infot, nout, lerr, ok);
        infot = 3;
        zgerfs("N", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZGERFS", infot, nout, lerr, ok);
        infot = 5;
        zgerfs("N", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("ZGERFS", infot, nout, lerr, ok);
        infot = 7;
        zgerfs("N", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("ZGERFS", infot, nout, lerr, ok);
        infot = 10;
        zgerfs("N", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("ZGERFS", infot, nout, lerr, ok);
        infot = 12;
        zgerfs("N", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("ZGERFS", infot, nout, lerr, ok);
        //
        //        ZGECON
        //
        srnamt = "ZGECON";
        infot = 1;
        zgecon("/", 0, a, 1, anrm, rcond, w, r, info);
        chkxer("ZGECON", infot, nout, lerr, ok);
        infot = 2;
        zgecon("1", -1, a, 1, anrm, rcond, w, r, info);
        chkxer("ZGECON", infot, nout, lerr, ok);
        infot = 4;
        zgecon("1", 2, a, 1, anrm, rcond, w, r, info);
        chkxer("ZGECON", infot, nout, lerr, ok);
        //
        //        ZGEEQU
        //
        srnamt = "ZGEEQU";
        infot = 1;
        zgeequ(-1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("ZGEEQU", infot, nout, lerr, ok);
        infot = 2;
        zgeequ(0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("ZGEEQU", infot, nout, lerr, ok);
        infot = 4;
        zgeequ(2, 2, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("ZGEEQU", infot, nout, lerr, ok);
        //
        //     Test error exits of the routines that use the LU decomposition
        //     of a general band matrix.
        //
    } else if (Mlsamen2, c2, "GB") {
        //
        //        ZGBTRF
        //
        srnamt = "ZGBTRF";
        infot = 1;
        zgbtrf(-1, 0, 0, 0, a, 1, ip, info);
        chkxer("ZGBTRF", infot, nout, lerr, ok);
        infot = 2;
        zgbtrf(0, -1, 0, 0, a, 1, ip, info);
        chkxer("ZGBTRF", infot, nout, lerr, ok);
        infot = 3;
        zgbtrf(1, 1, -1, 0, a, 1, ip, info);
        chkxer("ZGBTRF", infot, nout, lerr, ok);
        infot = 4;
        zgbtrf(1, 1, 0, -1, a, 1, ip, info);
        chkxer("ZGBTRF", infot, nout, lerr, ok);
        infot = 6;
        zgbtrf(2, 2, 1, 1, a, 3, ip, info);
        chkxer("ZGBTRF", infot, nout, lerr, ok);
        //
        //        ZGBTF2
        //
        srnamt = "ZGBTF2";
        infot = 1;
        zgbtf2(-1, 0, 0, 0, a, 1, ip, info);
        chkxer("ZGBTF2", infot, nout, lerr, ok);
        infot = 2;
        zgbtf2(0, -1, 0, 0, a, 1, ip, info);
        chkxer("ZGBTF2", infot, nout, lerr, ok);
        infot = 3;
        zgbtf2(1, 1, -1, 0, a, 1, ip, info);
        chkxer("ZGBTF2", infot, nout, lerr, ok);
        infot = 4;
        zgbtf2(1, 1, 0, -1, a, 1, ip, info);
        chkxer("ZGBTF2", infot, nout, lerr, ok);
        infot = 6;
        zgbtf2(2, 2, 1, 1, a, 3, ip, info);
        chkxer("ZGBTF2", infot, nout, lerr, ok);
        //
        //        ZGBTRS
        //
        srnamt = "ZGBTRS";
        infot = 1;
        zgbtrs("/", 0, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("ZGBTRS", infot, nout, lerr, ok);
        infot = 2;
        zgbtrs("N", -1, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("ZGBTRS", infot, nout, lerr, ok);
        infot = 3;
        zgbtrs("N", 1, -1, 0, 1, a, 1, ip, b, 1, info);
        chkxer("ZGBTRS", infot, nout, lerr, ok);
        infot = 4;
        zgbtrs("N", 1, 0, -1, 1, a, 1, ip, b, 1, info);
        chkxer("ZGBTRS", infot, nout, lerr, ok);
        infot = 5;
        zgbtrs("N", 1, 0, 0, -1, a, 1, ip, b, 1, info);
        chkxer("ZGBTRS", infot, nout, lerr, ok);
        infot = 7;
        zgbtrs("N", 2, 1, 1, 1, a, 3, ip, b, 2, info);
        chkxer("ZGBTRS", infot, nout, lerr, ok);
        infot = 10;
        zgbtrs("N", 2, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("ZGBTRS", infot, nout, lerr, ok);
        //
        //        ZGBRFS
        //
        srnamt = "ZGBRFS";
        infot = 1;
        zgbrfs("/", 0, 0, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZGBRFS", infot, nout, lerr, ok);
        infot = 2;
        zgbrfs("N", -1, 0, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZGBRFS", infot, nout, lerr, ok);
        infot = 3;
        zgbrfs("N", 1, -1, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZGBRFS", infot, nout, lerr, ok);
        infot = 4;
        zgbrfs("N", 1, 0, -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZGBRFS", infot, nout, lerr, ok);
        infot = 5;
        zgbrfs("N", 1, 0, 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZGBRFS", infot, nout, lerr, ok);
        infot = 7;
        zgbrfs("N", 2, 1, 1, 1, a, 2, af, 4, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("ZGBRFS", infot, nout, lerr, ok);
        infot = 9;
        zgbrfs("N", 2, 1, 1, 1, a, 3, af, 3, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("ZGBRFS", infot, nout, lerr, ok);
        infot = 12;
        zgbrfs("N", 2, 0, 0, 1, a, 1, af, 1, ip, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("ZGBRFS", infot, nout, lerr, ok);
        infot = 14;
        zgbrfs("N", 2, 0, 0, 1, a, 1, af, 1, ip, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("ZGBRFS", infot, nout, lerr, ok);
        //
        //        ZGBCON
        //
        srnamt = "ZGBCON";
        infot = 1;
        zgbcon("/", 0, 0, 0, a, 1, ip, anrm, rcond, w, r, info);
        chkxer("ZGBCON", infot, nout, lerr, ok);
        infot = 2;
        zgbcon("1", -1, 0, 0, a, 1, ip, anrm, rcond, w, r, info);
        chkxer("ZGBCON", infot, nout, lerr, ok);
        infot = 3;
        zgbcon("1", 1, -1, 0, a, 1, ip, anrm, rcond, w, r, info);
        chkxer("ZGBCON", infot, nout, lerr, ok);
        infot = 4;
        zgbcon("1", 1, 0, -1, a, 1, ip, anrm, rcond, w, r, info);
        chkxer("ZGBCON", infot, nout, lerr, ok);
        infot = 6;
        zgbcon("1", 2, 1, 1, a, 3, ip, anrm, rcond, w, r, info);
        chkxer("ZGBCON", infot, nout, lerr, ok);
        //
        //        ZGBEQU
        //
        srnamt = "ZGBEQU";
        infot = 1;
        zgbequ(-1, 0, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("ZGBEQU", infot, nout, lerr, ok);
        infot = 2;
        zgbequ(0, -1, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("ZGBEQU", infot, nout, lerr, ok);
        infot = 3;
        zgbequ(1, 1, -1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("ZGBEQU", infot, nout, lerr, ok);
        infot = 4;
        zgbequ(1, 1, 0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("ZGBEQU", infot, nout, lerr, ok);
        infot = 6;
        zgbequ(2, 2, 1, 1, a, 2, r1, r2, rcond, ccond, anrm, info);
        chkxer("ZGBEQU", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrge
    //
}
