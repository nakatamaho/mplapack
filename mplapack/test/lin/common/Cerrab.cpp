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

void Cerrab(INTEGER const nunit) {
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
    const INTEGER nmax = 4;
    INTEGER i = 0;
    COMPLEX a[nmax * nmax];
    COMPLEX af[nmax * nmax];
    COMPLEX b[nmax];
    COMPLEX r1[nmax];
    COMPLEX r2[nmax];
    COMPLEX w[2 * nmax];
    COMPLEX x[nmax];
    COMPLEX c[nmax];
    COMPLEX r[nmax];
    INTEGER ip[nmax];
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
        c[j - 1] = 0.0;
        r[j - 1] = 0.0;
        ip[j - 1] = j;
    }
    ok = true;
    //
    infot = 1;
    COMPLEX work[1];
    std::complex<float> swork[1];
    REAL rwork[1];
    INTEGER iter = 0;
    INTEGER info = 0;
    Ccgesv(-1, 0, a, 1, ip, b, 1, x, 1, work, swork, rwork, iter, info);
    chkxer("Ccgesv", infot, nout, lerr, ok);
    infot = 2;
    Ccgesv(0, -1, a, 1, ip, b, 1, x, 1, work, swork, rwork, iter, info);
    chkxer("Ccgesv", infot, nout, lerr, ok);
    infot = 4;
    Ccgesv(2, 1, a, 1, ip, b, 2, x, 2, work, swork, rwork, iter, info);
    chkxer("Ccgesv", infot, nout, lerr, ok);
    infot = 7;
    Ccgesv(2, 1, a, 2, ip, b, 1, x, 2, work, swork, rwork, iter, info);
    chkxer("Ccgesv", infot, nout, lerr, ok);
    infot = 9;
    Ccgesv(2, 1, a, 2, ip, b, 2, x, 1, work, swork, rwork, iter, info);
    chkxer("Ccgesv", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a6,' drivers passed the tests of the error exits')"), "Ccgesv";
    } else {
        write(nout, "(' *** ',a6,' drivers failed the tests of the error ','exits ***')"), "Ccgesv";
    }
    //
    //     End of Cerrab
    //
}
