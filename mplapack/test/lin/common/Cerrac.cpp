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

void Cerrac(INTEGER const nunit) {
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
    arr_2d<nmax, nmax, COMPLEX> a;
    arr_2d<nmax, nmax, COMPLEX> af;
    arr_1d<nmax, COMPLEX> b;
    arr_1d<nmax, COMPLEX> r1;
    arr_1d<nmax, COMPLEX> r2;
    arr_1d<2 * nmax, COMPLEX> w;
    arr_1d<nmax, COMPLEX> x;
    arr_1d<nmax, COMPLEX> c;
    arr_1d<nmax, COMPLEX> r;
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
    }
    ok = true;
    //
    cmn.srnamt = "Ccposv";
    infot = 1;
    arr_1d<nmax * nmax, COMPLEX> work;
    arr_1d<nmax * nmax, std::complex<float>> swork;
    arr_1d<nmax, REAL> rwork;
    INTEGER iter = 0;
    INTEGER info = 0;
    Ccposv("/", 0, 0, a, 1, b, 1, x, 1, work, swork, rwork, iter, info);
    chkxer("Ccposv", infot, nout, lerr, ok);
    infot = 2;
    Ccposv("U", -1, 0, a, 1, b, 1, x, 1, work, swork, rwork, iter, info);
    chkxer("Ccposv", infot, nout, lerr, ok);
    infot = 3;
    Ccposv("U", 0, -1, a, 1, b, 1, x, 1, work, swork, rwork, iter, info);
    chkxer("Ccposv", infot, nout, lerr, ok);
    infot = 5;
    Ccposv("U", 2, 1, a, 1, b, 2, x, 2, work, swork, rwork, iter, info);
    chkxer("Ccposv", infot, nout, lerr, ok);
    infot = 7;
    Ccposv("U", 2, 1, a, 2, b, 1, x, 2, work, swork, rwork, iter, info);
    chkxer("Ccposv", infot, nout, lerr, ok);
    infot = 9;
    Ccposv("U", 2, 1, a, 2, b, 2, x, 1, work, swork, rwork, iter, info);
    chkxer("Ccposv", infot, nout, lerr, ok);
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a6,' drivers passed the tests of the error exits')"), "Ccposv";
    } else {
        write(nout, "(' *** ',a6,' drivers failed the tests of the error ','exits ***')"), "Ccposv";
    }
    //
    //     End of Cerrac
    //
}
