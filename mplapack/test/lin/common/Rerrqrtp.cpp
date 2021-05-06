/*
 * Copyright (c) 2021
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

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

void Rerrqrtp(const char *path, INTEGER const nunit) {
    common_write write(cmn);
    // COMMON infoc
    INTEGER &infot = cmn.infot;
    INTEGER &nout = cmn.nout;
    bool &ok = cmn.ok;
    bool &lerr = cmn.lerr;
    // COMMON srnamc
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
    REAL a[nmax * nmax];
    REAL c[nmax * nmax];
    REAL t[nmax * nmax];
    REAL w[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / castREAL(i + j);
            c[(i - 1) + (j - 1) * ldc] = 1.0 / castREAL(i + j);
            t[(i - 1) + (j - 1) * ldt] = 1.0 / castREAL(i + j);
        }
        w[j - 1] = 0.0f;
    }
    ok = true;
    //
    //     Error exits for TPQRT factorization
    //
    //     Rtpqrt
    //
    infot = 1;
    REAL b[nmax * nmax];
    INTEGER info = 0;
    Rtpqrt(-1, 1, 0, 1, a, 1, b, 1, t, 1, w, info);
    chkxer("Rtpqrt", infot, nout, lerr, ok);
    infot = 2;
    Rtpqrt(1, -1, 0, 1, a, 1, b, 1, t, 1, w, info);
    chkxer("Rtpqrt", infot, nout, lerr, ok);
    infot = 3;
    Rtpqrt(0, 1, -1, 1, a, 1, b, 1, t, 1, w, info);
    chkxer("Rtpqrt", infot, nout, lerr, ok);
    infot = 3;
    Rtpqrt(0, 1, 1, 1, a, 1, b, 1, t, 1, w, info);
    chkxer("Rtpqrt", infot, nout, lerr, ok);
    infot = 4;
    Rtpqrt(0, 1, 0, 0, a, 1, b, 1, t, 1, w, info);
    chkxer("Rtpqrt", infot, nout, lerr, ok);
    infot = 4;
    Rtpqrt(0, 1, 0, 2, a, 1, b, 1, t, 1, w, info);
    chkxer("Rtpqrt", infot, nout, lerr, ok);
    infot = 6;
    Rtpqrt(1, 2, 0, 2, a, 1, b, 1, t, 1, w, info);
    chkxer("Rtpqrt", infot, nout, lerr, ok);
    infot = 8;
    Rtpqrt(2, 1, 0, 1, a, 1, b, 1, t, 1, w, info);
    chkxer("Rtpqrt", infot, nout, lerr, ok);
    infot = 10;
    Rtpqrt(2, 2, 1, 2, a, 2, b, 2, t, 1, w, info);
    chkxer("Rtpqrt", infot, nout, lerr, ok);
    //
    //     Rtpqrt2
    //
    infot = 1;
    Rtpqrt2(-1, 0, 0, a, 1, b, 1, t, 1, info);
    chkxer("Rtpqrt2", infot, nout, lerr, ok);
    infot = 2;
    Rtpqrt2(0, -1, 0, a, 1, b, 1, t, 1, info);
    chkxer("Rtpqrt2", infot, nout, lerr, ok);
    infot = 3;
    Rtpqrt2(0, 0, -1, a, 1, b, 1, t, 1, info);
    chkxer("Rtpqrt2", infot, nout, lerr, ok);
    infot = 5;
    Rtpqrt2(2, 2, 0, a, 1, b, 2, t, 2, info);
    chkxer("Rtpqrt2", infot, nout, lerr, ok);
    infot = 7;
    Rtpqrt2(2, 2, 0, a, 2, b, 1, t, 2, info);
    chkxer("Rtpqrt2", infot, nout, lerr, ok);
    infot = 9;
    Rtpqrt2(2, 2, 0, a, 2, b, 2, t, 1, info);
    chkxer("Rtpqrt2", infot, nout, lerr, ok);
    //
    //     Rtpmqrt
    //
    infot = 1;
    Rtpmqrt("/", "N", 0, 0, 0, 0, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    infot = 2;
    Rtpmqrt("L", "/", 0, 0, 0, 0, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    infot = 3;
    Rtpmqrt("L", "N", -1, 0, 0, 0, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    infot = 4;
    Rtpmqrt("L", "N", 0, -1, 0, 0, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    infot = 5;
    Rtpmqrt("L", "N", 0, 0, -1, 0, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    infot = 6;
    Rtpmqrt("L", "N", 0, 0, 0, -1, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    infot = 7;
    Rtpmqrt("L", "N", 0, 0, 0, 0, 0, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    infot = 9;
    Rtpmqrt("R", "N", 1, 2, 1, 1, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    infot = 9;
    Rtpmqrt("L", "N", 2, 1, 1, 1, 1, a, 1, t, 1, b, 1, c, 1, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    infot = 11;
    Rtpmqrt("R", "N", 1, 1, 1, 1, 1, a, 1, t, 0, b, 1, c, 1, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    infot = 13;
    Rtpmqrt("L", "N", 1, 1, 1, 1, 1, a, 1, t, 1, b, 0, c, 1, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    infot = 15;
    Rtpmqrt("L", "N", 1, 1, 1, 1, 1, a, 1, t, 1, b, 1, c, 0, w, info);
    chkxer("Rtpmqrt", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrqrt
    //
}
