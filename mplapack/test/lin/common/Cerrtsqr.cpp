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

void Cerrtsqr(common &cmn, const char *path, INTEGER const nunit) {
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
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 2;
    INTEGER i = 0;
    arr_2d<nmax, nmax, COMPLEX> a(fill0);
    arr_2d<nmax, nmax, COMPLEX> c(fill0);
    arr_2d<nmax, nmax, COMPLEX> t(fill0);
    arr_1d<nmax, COMPLEX> w(fill0);
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / (i + j).real();
            c[(i - 1) + (j - 1) * ldc] = 1.0 / (i + j).real();
            t[(i - 1) + (j - 1) * ldt] = 1.0 / (i + j).real();
        }
        w[j - 1] = 0.0;
    }
    ok = true;
    //
    //     Error exits for TS factorization
    //
    //     ZGEQR
    //
    srnamt = "ZGEQR";
    infot = 1;
    arr_1d<nmax, COMPLEX> tau(fill0);
    INTEGER info = 0;
    zgeqr(-1, 0, a, 1, tau, 1, w, 1, info);
    chkxer("ZGEQR", infot, nout, lerr, ok);
    infot = 2;
    zgeqr(0, -1, a, 1, tau, 1, w, 1, info);
    chkxer("ZGEQR", infot, nout, lerr, ok);
    infot = 4;
    zgeqr(1, 1, a, 0, tau, 1, w, 1, info);
    chkxer("ZGEQR", infot, nout, lerr, ok);
    infot = 6;
    zgeqr(3, 2, a, 3, tau, 1, w, 1, info);
    chkxer("ZGEQR", infot, nout, lerr, ok);
    infot = 8;
    zgeqr(3, 2, a, 3, tau, 8, w, 0, info);
    chkxer("ZGEQR", infot, nout, lerr, ok);
    //
    //     ZGEMQR
    //
    tau[1 - 1] = 1;
    tau[2 - 1] = 1;
    srnamt = "ZGEMQR";
    INTEGER nb = 1;
    infot = 1;
    zgemqr("/", "N", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    infot = 2;
    zgemqr("L", "/", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    infot = 3;
    zgemqr("L", "N", -1, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    infot = 4;
    zgemqr("L", "N", 0, -1, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    infot = 5;
    zgemqr("L", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    infot = 5;
    zgemqr("R", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    infot = 7;
    zgemqr("L", "N", 2, 1, 0, a, 0, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    infot = 9;
    zgemqr("R", "N", 2, 2, 1, a, 2, tau, 0, c, 1, w, 1, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    infot = 9;
    zgemqr("L", "N", 2, 2, 1, a, 2, tau, 0, c, 1, w, 1, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    infot = 11;
    zgemqr("L", "N", 2, 1, 1, a, 2, tau, 6, c, 0, w, 1, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    infot = 13;
    zgemqr("L", "N", 2, 2, 1, a, 2, tau, 6, c, 2, w, 0, info);
    chkxer("ZGEMQR", infot, nout, lerr, ok);
    //
    //     ZGELQ
    //
    srnamt = "ZGELQ";
    infot = 1;
    zgelq(-1, 0, a, 1, tau, 1, w, 1, info);
    chkxer("ZGELQ", infot, nout, lerr, ok);
    infot = 2;
    zgelq(0, -1, a, 1, tau, 1, w, 1, info);
    chkxer("ZGELQ", infot, nout, lerr, ok);
    infot = 4;
    zgelq(1, 1, a, 0, tau, 1, w, 1, info);
    chkxer("ZGELQ", infot, nout, lerr, ok);
    infot = 6;
    zgelq(2, 3, a, 3, tau, 1, w, 1, info);
    chkxer("ZGELQ", infot, nout, lerr, ok);
    infot = 8;
    zgelq(2, 3, a, 3, tau, 8, w, 0, info);
    chkxer("ZGELQ", infot, nout, lerr, ok);
    //
    //     ZGEMLQ
    //
    tau[1 - 1] = 1;
    tau[2 - 1] = 1;
    srnamt = "ZGEMLQ";
    nb = 1;
    infot = 1;
    zgemlq("/", "N", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    infot = 2;
    zgemlq("L", "/", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    infot = 3;
    zgemlq("L", "N", -1, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    infot = 4;
    zgemlq("L", "N", 0, -1, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    infot = 5;
    zgemlq("L", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    infot = 5;
    zgemlq("R", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    infot = 7;
    zgemlq("L", "N", 1, 2, 0, a, 0, tau, 1, c, 1, w, 1, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    infot = 9;
    zgemlq("R", "N", 2, 2, 1, a, 1, tau, 0, c, 1, w, 1, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    infot = 9;
    zgemlq("L", "N", 2, 2, 1, a, 1, tau, 0, c, 1, w, 1, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    infot = 11;
    zgemlq("L", "N", 1, 2, 1, a, 1, tau, 6, c, 0, w, 1, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    infot = 13;
    zgemlq("L", "N", 2, 2, 1, a, 2, tau, 6, c, 2, w, 0, info);
    chkxer("ZGEMLQ", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrtsqr
    //
}
