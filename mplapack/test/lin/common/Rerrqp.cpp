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

void Rerrqp(const char *path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    bool ok;
    bool lerr;
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
    INTEGER nout = nunit;
    char c2[2];
    c2[0] = path[1];
    c2[1] = path[2];
    const INTEGER nmax = 3;
    INTEGER lw = 3 * nmax + 1;
    REAL a[nmax * nmax];
    INTEGER lda = nmax;
    a[(1 - 1) + (1 - 1) * lda] = 1.0;
    a[(1 - 1) + (2 - 1) * lda] = 2.0e+0;
    a[(2 - 1) + (2 - 1) * lda] = 3.0e+0;
    a[(2 - 1)] = 4.0e+0;
    ok = true;
    //
    INTEGER ip[nmax];
    REAL tau[nmax];
    REAL w[3 * nmax + 1];
    INTEGER info = 0;
    bool infot = 0;
    if (Mlsamen(2, c2, "QP")) {
        //
        //        Test error exits for QR factorization with pivoting
        //
        //        Rgeqp3
        //
        infot = 1;
        Rgeqp3(-1, 0, a, 1, ip, tau, w, lw, info);
        chkxer("Rgeqp3", infot, nout, lerr, ok);
        infot = 2;
        Rgeqp3(1, -1, a, 1, ip, tau, w, lw, info);
        chkxer("Rgeqp3", infot, nout, lerr, ok);
        infot = 4;
        Rgeqp3(2, 3, a, 1, ip, tau, w, lw, info);
        chkxer("Rgeqp3", infot, nout, lerr, ok);
        infot = 8;
        Rgeqp3(2, 2, a, 2, ip, tau, w, lw - 10, info);
        chkxer("Rgeqp3", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrqp
    //
}
