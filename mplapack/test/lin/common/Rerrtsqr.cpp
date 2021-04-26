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

void Rerrtsqr(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_2d<nmax, nmax, REAL> a(fill0);
    arr_2d<nmax, nmax, REAL> c(fill0);
    arr_2d<nmax, nmax, REAL> t(fill0);
    arr_1d<nmax, REAL> w(fill0);
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
    //     DGEQR
    //
    srnamt = "DGEQR";
    infot = 1;
    arr_1d<nmax * 2, REAL> tau(fill0);
    INTEGER info = 0;
    dgeqr(-1, 0, a, 1, tau, 1, w, 1, info);
    chkxer("DGEQR", infot, nout, lerr, ok);
    infot = 2;
    dgeqr(0, -1, a, 1, tau, 1, w, 1, info);
    chkxer("DGEQR", infot, nout, lerr, ok);
    infot = 4;
    dgeqr(1, 1, a, 0, tau, 1, w, 1, info);
    chkxer("DGEQR", infot, nout, lerr, ok);
    infot = 6;
    dgeqr(3, 2, a, 3, tau, 1, w, 1, info);
    chkxer("DGEQR", infot, nout, lerr, ok);
    infot = 8;
    dgeqr(3, 2, a, 3, tau, 7, w, 0, info);
    chkxer("DGEQR", infot, nout, lerr, ok);
    //
    //     DGEMQR
    //
    tau[1 - 1] = 1;
    tau[2 - 1] = 1;
    tau[3 - 1] = 1;
    tau[4 - 1] = 1;
    srnamt = "DGEMQR";
    INTEGER nb = 1;
    infot = 1;
    dgemqr("/", "N", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    infot = 2;
    dgemqr("L", "/", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    infot = 3;
    dgemqr("L", "N", -1, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    infot = 4;
    dgemqr("L", "N", 0, -1, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    infot = 5;
    dgemqr("L", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    infot = 5;
    dgemqr("R", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    infot = 7;
    dgemqr("L", "N", 2, 1, 0, a, 0, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    infot = 9;
    dgemqr("R", "N", 2, 2, 1, a, 2, tau, 0, c, 1, w, 1, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    infot = 9;
    dgemqr("L", "N", 2, 2, 1, a, 2, tau, 0, c, 1, w, 1, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    infot = 11;
    dgemqr("L", "N", 2, 1, 1, a, 2, tau, 6, c, 0, w, 1, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    infot = 13;
    dgemqr("L", "N", 2, 2, 1, a, 2, tau, 6, c, 2, w, 0, info);
    chkxer("DGEMQR", infot, nout, lerr, ok);
    //
    //     DGELQ
    //
    srnamt = "DGELQ";
    infot = 1;
    dgelq(-1, 0, a, 1, tau, 1, w, 1, info);
    chkxer("DGELQ", infot, nout, lerr, ok);
    infot = 2;
    dgelq(0, -1, a, 1, tau, 1, w, 1, info);
    chkxer("DGELQ", infot, nout, lerr, ok);
    infot = 4;
    dgelq(1, 1, a, 0, tau, 1, w, 1, info);
    chkxer("DGELQ", infot, nout, lerr, ok);
    infot = 6;
    dgelq(2, 3, a, 3, tau, 1, w, 1, info);
    chkxer("DGELQ", infot, nout, lerr, ok);
    infot = 8;
    dgelq(2, 3, a, 3, tau, 7, w, 0, info);
    chkxer("DGELQ", infot, nout, lerr, ok);
    //
    //     DGEMLQ
    //
    tau[1 - 1] = 1;
    tau[2 - 1] = 1;
    srnamt = "DGEMLQ";
    nb = 1;
    infot = 1;
    dgemlq("/", "N", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    infot = 2;
    dgemlq("L", "/", 0, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    infot = 3;
    dgemlq("L", "N", -1, 0, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    infot = 4;
    dgemlq("L", "N", 0, -1, 0, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    infot = 5;
    dgemlq("L", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    infot = 5;
    dgemlq("R", "N", 0, 0, -1, a, 1, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    infot = 7;
    dgemlq("L", "N", 1, 2, 0, a, 0, tau, 1, c, 1, w, 1, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    infot = 9;
    dgemlq("R", "N", 2, 2, 1, a, 1, tau, 0, c, 1, w, 1, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    infot = 9;
    dgemlq("L", "N", 2, 2, 1, a, 1, tau, 0, c, 1, w, 1, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    infot = 11;
    dgemlq("L", "N", 1, 2, 1, a, 1, tau, 6, c, 0, w, 1, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    infot = 13;
    dgemlq("L", "N", 2, 2, 1, a, 2, tau, 6, c, 2, w, 0, info);
    chkxer("DGEMLQ", infot, nout, lerr, ok);
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrtsqr
    //
}
