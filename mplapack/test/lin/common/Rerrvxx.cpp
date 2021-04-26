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

void Rerrvx(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_1d<nmax, REAL> e(fill0);
    arr_1d<nmax, REAL> r1(fill0);
    arr_1d<nmax, REAL> r2(fill0);
    arr_1d<2 * nmax, REAL> w(fill0);
    arr_1d<nmax, REAL> x(fill0);
    arr_1d<nmax, REAL> c(fill0);
    arr_1d<nmax, REAL> r(fill0);
    arr_1d<nmax, int> ip(fill0);
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / (i + j).real();
            af[(i - 1) + (j - 1) * ldaf] = 1.0 / (i + j).real();
        }
        b[j - 1] = 0.e+0;
        e[j - 1] = 0.e+0;
        r1[j - 1] = 0.e+0;
        r2[j - 1] = 0.e+0;
        w[j - 1] = 0.e+0;
        x[j - 1] = 0.e+0;
        c[j - 1] = 0.e+0;
        r[j - 1] = 0.e+0;
        ip[j - 1] = j;
    }
    char eq = " ";
    ok = true;
    //
    INTEGER info = 0;
    REAL rcond = 0.0;
    arr_1d<nmax, int> iw(fill0);
    INTEGER n_err_bnds = 0;
    INTEGER nparams = 0;
    REAL rpvgrw = 0.0;
    REAL berr = 0.0;
    arr_2d<nmax, 3, REAL> err_bnds_n(fill0);
    arr_2d<nmax, 3, REAL> err_bnds_c(fill0);
    arr_1d<1, REAL> params(fill0);
    const float one = 1.0;
    if (Mlsamen2, c2, "GE") {
        //
        //        DGESV
        //
        srnamt = "DGESV ";
        infot = 1;
        dgesv(-1, 0, a, 1, ip, b, 1, info);
        chkxer("DGESV ", infot, nout, lerr, ok);
        infot = 2;
        dgesv(0, -1, a, 1, ip, b, 1, info);
        chkxer("DGESV ", infot, nout, lerr, ok);
        infot = 4;
        dgesv(2, 1, a, 1, ip, b, 2, info);
        chkxer("DGESV ", infot, nout, lerr, ok);
        infot = 7;
        dgesv(2, 1, a, 2, ip, b, 1, info);
        chkxer("DGESV ", infot, nout, lerr, ok);
        //
        //        DGESVX
        //
        srnamt = "DGESVX";
        infot = 1;
        dgesvx("/", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        infot = 2;
        dgesvx("N", "/", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        infot = 3;
        dgesvx("N", "N", -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        infot = 4;
        dgesvx("N", "N", 0, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        infot = 6;
        dgesvx("N", "N", 2, 1, a, 1, af, 2, ip, eq, r, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        infot = 8;
        dgesvx("N", "N", 2, 1, a, 2, af, 1, ip, eq, r, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        infot = 10;
        eq = "/";
        dgesvx("F", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        infot = 11;
        eq = "R";
        dgesvx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        infot = 12;
        eq = "C";
        dgesvx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        infot = 14;
        dgesvx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        infot = 16;
        dgesvx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGESVX", infot, nout, lerr, ok);
        //
        //        DGESVXX
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "DGESVXX";
        infot = 1;
        dgesvxx("/", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        infot = 2;
        dgesvxx("N", "/", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        infot = 3;
        dgesvxx("N", "N", -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        infot = 4;
        dgesvxx("N", "N", 0, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        infot = 6;
        dgesvxx("N", "N", 2, 1, a, 1, af, 2, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        infot = 8;
        dgesvxx("N", "N", 2, 1, a, 2, af, 1, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        infot = 10;
        eq = "/";
        dgesvxx("F", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        infot = 11;
        eq = "R";
        dgesvxx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        infot = 12;
        eq = "C";
        dgesvxx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        infot = 14;
        dgesvxx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        infot = 16;
        dgesvxx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGESVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "GB") {
        //
        //        DGBSV
        //
        srnamt = "DGBSV ";
        infot = 1;
        dgbsv(-1, 0, 0, 0, a, 1, ip, b, 1, info);
        chkxer("DGBSV ", infot, nout, lerr, ok);
        infot = 2;
        dgbsv(1, -1, 0, 0, a, 1, ip, b, 1, info);
        chkxer("DGBSV ", infot, nout, lerr, ok);
        infot = 3;
        dgbsv(1, 0, -1, 0, a, 1, ip, b, 1, info);
        chkxer("DGBSV ", infot, nout, lerr, ok);
        infot = 4;
        dgbsv(0, 0, 0, -1, a, 1, ip, b, 1, info);
        chkxer("DGBSV ", infot, nout, lerr, ok);
        infot = 6;
        dgbsv(1, 1, 1, 0, a, 3, ip, b, 1, info);
        chkxer("DGBSV ", infot, nout, lerr, ok);
        infot = 9;
        dgbsv(2, 0, 0, 0, a, 1, ip, b, 1, info);
        chkxer("DGBSV ", infot, nout, lerr, ok);
        //
        //        DGBSVX
        //
        srnamt = "DGBSVX";
        infot = 1;
        dgbsvx("/", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 2;
        dgbsvx("N", "/", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 3;
        dgbsvx("N", "N", -1, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 4;
        dgbsvx("N", "N", 1, -1, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 5;
        dgbsvx("N", "N", 1, 0, -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 6;
        dgbsvx("N", "N", 0, 0, 0, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 8;
        dgbsvx("N", "N", 1, 1, 1, 0, a, 2, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 10;
        dgbsvx("N", "N", 1, 1, 1, 0, a, 3, af, 3, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 12;
        eq = "/";
        dgbsvx("F", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 13;
        eq = "R";
        dgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 14;
        eq = "C";
        dgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 16;
        dgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        infot = 18;
        dgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGBSVX", infot, nout, lerr, ok);
        //
        //        DGBSVXX
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "DGBSVXX";
        infot = 1;
        dgbsvxx("/", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 2;
        dgbsvxx("N", "/", 0, 1, 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 3;
        dgbsvxx("N", "N", -1, 1, 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 4;
        dgbsvxx("N", "N", 2, -1, 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 5;
        dgbsvxx("N", "N", 2, 1, -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 6;
        dgbsvxx("N", "N", 0, 1, 1, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 8;
        dgbsvxx("N", "N", 2, 1, 1, 1, a, 2, af, 2, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 10;
        dgbsvxx("N", "N", 2, 1, 1, 1, a, 3, af, 3, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 12;
        eq = "/";
        dgbsvxx("F", "N", 0, 1, 1, 0, a, 3, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 13;
        eq = "R";
        dgbsvxx("F", "N", 1, 1, 1, 0, a, 3, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 14;
        eq = "C";
        dgbsvxx("F", "N", 1, 1, 1, 0, a, 3, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 15;
        dgbsvxx("N", "N", 2, 1, 1, 1, a, 3, af, 4, ip, eq, r, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        infot = 16;
        dgbsvxx("N", "N", 2, 1, 1, 1, a, 3, af, 4, ip, eq, r, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DGBSVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "GT") {
        //
        //        DGTSV
        //
        srnamt = "DGTSV ";
        infot = 1;
        dgtsv(-1, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], b, 1, info);
        chkxer("DGTSV ", infot, nout, lerr, ok);
        infot = 2;
        dgtsv(0, -1, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], b, 1, info);
        chkxer("DGTSV ", infot, nout, lerr, ok);
        infot = 7;
        dgtsv(2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], b, 1, info);
        chkxer("DGTSV ", infot, nout, lerr, ok);
        //
        //        DGTSVX
        //
        srnamt = "DGTSVX";
        infot = 1;
        dgtsvx("/", "N", 0, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGTSVX", infot, nout, lerr, ok);
        infot = 2;
        dgtsvx("N", "/", 0, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGTSVX", infot, nout, lerr, ok);
        infot = 3;
        dgtsvx("N", "N", -1, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGTSVX", infot, nout, lerr, ok);
        infot = 4;
        dgtsvx("N", "N", 0, -1, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGTSVX", infot, nout, lerr, ok);
        infot = 14;
        dgtsvx("N", "N", 2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DGTSVX", infot, nout, lerr, ok);
        infot = 16;
        dgtsvx("N", "N", 2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DGTSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PO") {
        //
        //        DPOSV
        //
        srnamt = "DPOSV ";
        infot = 1;
        dposv("/", 0, 0, a, 1, b, 1, info);
        chkxer("DPOSV ", infot, nout, lerr, ok);
        infot = 2;
        dposv("U", -1, 0, a, 1, b, 1, info);
        chkxer("DPOSV ", infot, nout, lerr, ok);
        infot = 3;
        dposv("U", 0, -1, a, 1, b, 1, info);
        chkxer("DPOSV ", infot, nout, lerr, ok);
        infot = 5;
        dposv("U", 2, 0, a, 1, b, 2, info);
        chkxer("DPOSV ", infot, nout, lerr, ok);
        infot = 7;
        dposv("U", 2, 0, a, 2, b, 1, info);
        chkxer("DPOSV ", infot, nout, lerr, ok);
        //
        //        DPOSVX
        //
        srnamt = "DPOSVX";
        infot = 1;
        dposvx("/", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPOSVX", infot, nout, lerr, ok);
        infot = 2;
        dposvx("N", "/", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPOSVX", infot, nout, lerr, ok);
        infot = 3;
        dposvx("N", "U", -1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPOSVX", infot, nout, lerr, ok);
        infot = 4;
        dposvx("N", "U", 0, -1, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPOSVX", infot, nout, lerr, ok);
        infot = 6;
        dposvx("N", "U", 2, 0, a, 1, af, 2, eq, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DPOSVX", infot, nout, lerr, ok);
        infot = 8;
        dposvx("N", "U", 2, 0, a, 2, af, 1, eq, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DPOSVX", infot, nout, lerr, ok);
        infot = 9;
        eq = "/";
        dposvx("F", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPOSVX", infot, nout, lerr, ok);
        infot = 10;
        eq = "Y";
        dposvx("F", "U", 1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPOSVX", infot, nout, lerr, ok);
        infot = 12;
        dposvx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DPOSVX", infot, nout, lerr, ok);
        infot = 14;
        dposvx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPOSVX", infot, nout, lerr, ok);
        //
        //        DPOSVXX
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "DPOSVXX";
        infot = 1;
        dposvxx("/", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DPOSVXX", infot, nout, lerr, ok);
        infot = 2;
        dposvxx("N", "/", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DPOSVXX", infot, nout, lerr, ok);
        infot = 3;
        dposvxx("N", "U", -1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DPOSVXX", infot, nout, lerr, ok);
        infot = 4;
        dposvxx("N", "U", 0, -1, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DPOSVXX", infot, nout, lerr, ok);
        infot = 6;
        dposvxx("N", "U", 2, 0, a, 1, af, 2, eq, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DPOSVXX", infot, nout, lerr, ok);
        infot = 8;
        dposvxx("N", "U", 2, 0, a, 2, af, 1, eq, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DPOSVXX", infot, nout, lerr, ok);
        infot = 9;
        eq = "/";
        dposvxx("F", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DPOSVXX", infot, nout, lerr, ok);
        infot = 10;
        eq = "Y";
        dposvxx("F", "U", 1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DPOSVXX", infot, nout, lerr, ok);
        infot = 12;
        dposvxx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DPOSVXX", infot, nout, lerr, ok);
        infot = 14;
        dposvxx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DPOSVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PP") {
        //
        //        DPPSV
        //
        srnamt = "DPPSV ";
        infot = 1;
        dppsv("/", 0, 0, a, b, 1, info);
        chkxer("DPPSV ", infot, nout, lerr, ok);
        infot = 2;
        dppsv("U", -1, 0, a, b, 1, info);
        chkxer("DPPSV ", infot, nout, lerr, ok);
        infot = 3;
        dppsv("U", 0, -1, a, b, 1, info);
        chkxer("DPPSV ", infot, nout, lerr, ok);
        infot = 6;
        dppsv("U", 2, 0, a, b, 1, info);
        chkxer("DPPSV ", infot, nout, lerr, ok);
        //
        //        DPPSVX
        //
        srnamt = "DPPSVX";
        infot = 1;
        dppsvx("/", "U", 0, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPPSVX", infot, nout, lerr, ok);
        infot = 2;
        dppsvx("N", "/", 0, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPPSVX", infot, nout, lerr, ok);
        infot = 3;
        dppsvx("N", "U", -1, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPPSVX", infot, nout, lerr, ok);
        infot = 4;
        dppsvx("N", "U", 0, -1, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPPSVX", infot, nout, lerr, ok);
        infot = 7;
        eq = "/";
        dppsvx("F", "U", 0, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPPSVX", infot, nout, lerr, ok);
        infot = 8;
        eq = "Y";
        dppsvx("F", "U", 1, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPPSVX", infot, nout, lerr, ok);
        infot = 10;
        dppsvx("N", "U", 2, 0, a, af, eq, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DPPSVX", infot, nout, lerr, ok);
        infot = 12;
        dppsvx("N", "U", 2, 0, a, af, eq, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPPSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PB") {
        //
        //        DPBSV
        //
        srnamt = "DPBSV ";
        infot = 1;
        dpbsv("/", 0, 0, 0, a, 1, b, 1, info);
        chkxer("DPBSV ", infot, nout, lerr, ok);
        infot = 2;
        dpbsv("U", -1, 0, 0, a, 1, b, 1, info);
        chkxer("DPBSV ", infot, nout, lerr, ok);
        infot = 3;
        dpbsv("U", 1, -1, 0, a, 1, b, 1, info);
        chkxer("DPBSV ", infot, nout, lerr, ok);
        infot = 4;
        dpbsv("U", 0, 0, -1, a, 1, b, 1, info);
        chkxer("DPBSV ", infot, nout, lerr, ok);
        infot = 6;
        dpbsv("U", 1, 1, 0, a, 1, b, 2, info);
        chkxer("DPBSV ", infot, nout, lerr, ok);
        infot = 8;
        dpbsv("U", 2, 0, 0, a, 1, b, 1, info);
        chkxer("DPBSV ", infot, nout, lerr, ok);
        //
        //        DPBSVX
        //
        srnamt = "DPBSVX";
        infot = 1;
        dpbsvx("/", "U", 0, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        infot = 2;
        dpbsvx("N", "/", 0, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        infot = 3;
        dpbsvx("N", "U", -1, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        infot = 4;
        dpbsvx("N", "U", 1, -1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        infot = 5;
        dpbsvx("N", "U", 0, 0, -1, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        infot = 7;
        dpbsvx("N", "U", 1, 1, 0, a, 1, af, 2, eq, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        infot = 9;
        dpbsvx("N", "U", 1, 1, 0, a, 2, af, 1, eq, c, b, 2, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        infot = 10;
        eq = "/";
        dpbsvx("F", "U", 0, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        infot = 11;
        eq = "Y";
        dpbsvx("F", "U", 1, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        infot = 13;
        dpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        infot = 15;
        dpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, eq, c, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DPBSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PT") {
        //
        //        DPTSV
        //
        srnamt = "DPTSV ";
        infot = 1;
        dptsv(-1, 0, &a[(1 - 1)], &a[(2 - 1) * lda], b, 1, info);
        chkxer("DPTSV ", infot, nout, lerr, ok);
        infot = 2;
        dptsv(0, -1, &a[(1 - 1)], &a[(2 - 1) * lda], b, 1, info);
        chkxer("DPTSV ", infot, nout, lerr, ok);
        infot = 6;
        dptsv(2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], b, 1, info);
        chkxer("DPTSV ", infot, nout, lerr, ok);
        //
        //        DPTSVX
        //
        srnamt = "DPTSVX";
        infot = 1;
        dptsvx("/", 0, 0, &a[(1 - 1)], &a[(2 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], b, 1, x, 1, rcond, r1, r2, w, info);
        chkxer("DPTSVX", infot, nout, lerr, ok);
        infot = 2;
        dptsvx("N", -1, 0, &a[(1 - 1)], &a[(2 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], b, 1, x, 1, rcond, r1, r2, w, info);
        chkxer("DPTSVX", infot, nout, lerr, ok);
        infot = 3;
        dptsvx("N", 0, -1, &a[(1 - 1)], &a[(2 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], b, 1, x, 1, rcond, r1, r2, w, info);
        chkxer("DPTSVX", infot, nout, lerr, ok);
        infot = 9;
        dptsvx("N", 2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], b, 1, x, 2, rcond, r1, r2, w, info);
        chkxer("DPTSVX", infot, nout, lerr, ok);
        infot = 11;
        dptsvx("N", 2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], b, 2, x, 1, rcond, r1, r2, w, info);
        chkxer("DPTSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SY") {
        //
        //        DSYSV
        //
        srnamt = "DSYSV ";
        infot = 1;
        dsysv("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("DSYSV ", infot, nout, lerr, ok);
        infot = 2;
        dsysv("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("DSYSV ", infot, nout, lerr, ok);
        infot = 3;
        dsysv("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("DSYSV ", infot, nout, lerr, ok);
        infot = 5;
        dsysv("U", 2, 0, a, 1, ip, b, 2, w, 1, info);
        chkxer("DSYSV_ROOK", infot, nout, lerr, ok);
        infot = 8;
        dsysv("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("DSYSV ", infot, nout, lerr, ok);
        infot = 10;
        dsysv("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("DSYSV ", infot, nout, lerr, ok);
        infot = 10;
        dsysv("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        chkxer("DSYSV ", infot, nout, lerr, ok);
        //
        //        DSYSVX
        //
        srnamt = "DSYSVX";
        infot = 1;
        dsysvx("/", "U", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        chkxer("DSYSVX", infot, nout, lerr, ok);
        infot = 2;
        dsysvx("N", "/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        chkxer("DSYSVX", infot, nout, lerr, ok);
        infot = 3;
        dsysvx("N", "U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        chkxer("DSYSVX", infot, nout, lerr, ok);
        infot = 4;
        dsysvx("N", "U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, iw, info);
        chkxer("DSYSVX", infot, nout, lerr, ok);
        infot = 6;
        dsysvx("N", "U", 2, 0, a, 1, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 4, iw, info);
        chkxer("DSYSVX", infot, nout, lerr, ok);
        infot = 8;
        dsysvx("N", "U", 2, 0, a, 2, af, 1, ip, b, 2, x, 2, rcond, r1, r2, w, 4, iw, info);
        chkxer("DSYSVX", infot, nout, lerr, ok);
        infot = 11;
        dsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 1, x, 2, rcond, r1, r2, w, 4, iw, info);
        chkxer("DSYSVX", infot, nout, lerr, ok);
        infot = 13;
        dsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 1, rcond, r1, r2, w, 4, iw, info);
        chkxer("DSYSVX", infot, nout, lerr, ok);
        infot = 18;
        dsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 3, iw, info);
        chkxer("DSYSVX", infot, nout, lerr, ok);
        //
        //        DSYSVXX
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "DSYSVXX";
        infot = 1;
        eq = "N";
        dsysvxx("/", "U", 0, 0, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        infot = 2;
        dsysvxx("N", "/", 0, 0, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        infot = 3;
        dsysvxx("N", "U", -1, 0, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        infot = 4;
        eq = "/";
        dsysvxx("N", "U", 0, -1, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        eq = "Y";
        infot = 6;
        dsysvxx("N", "U", 2, 0, a, 1, af, 2, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        infot = 8;
        dsysvxx("N", "U", 2, 0, a, 2, af, 1, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        infot = 10;
        dsysvxx("F", "U", 2, 0, a, 2, af, 2, ip, "A", r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        infot = 11;
        eq = "Y";
        dsysvxx("F", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        infot = 11;
        eq = "Y";
        r[1 - 1] = -one;
        dsysvxx("F", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        infot = 13;
        eq = "N";
        dsysvxx("N", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        infot = 15;
        dsysvxx("N", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("DSYSVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SR") {
        //
        //        DSYSV_ROOK
        //
        srnamt = "DSYSV_ROOK";
        infot = 1;
        dsysv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("DSYSV_ROOK", infot, nout, lerr, ok);
        infot = 2;
        dsysv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("DSYSV_ROOK", infot, nout, lerr, ok);
        infot = 3;
        dsysv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("DSYSV_ROOK", infot, nout, lerr, ok);
        infot = 5;
        dsysv_rook("U", 2, 0, a, 1, ip, b, 2, w, 1, info);
        chkxer("DSYSV_ROOK", infot, nout, lerr, ok);
        infot = 8;
        dsysv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("DSYSV_ROOK", infot, nout, lerr, ok);
        infot = 10;
        dsysv_rook("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("DSYSV_ROOK", infot, nout, lerr, ok);
        infot = 10;
        dsysv_rook("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        chkxer("DSYSV_ROOK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SK") {
        //
        //        DSYSV_RK
        //
        //        Test error exits of the driver that uses factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting with the new storage
        //        format for factors L ( or U) and D.
        //
        //        L (or U) is stored in A, diagonal of D is stored on the
        //        diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        srnamt = "DSYSV_RK";
        infot = 1;
        dsysv_rk("/", 0, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("DSYSV_RK", infot, nout, lerr, ok);
        infot = 2;
        dsysv_rk("U", -1, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("DSYSV_RK", infot, nout, lerr, ok);
        infot = 3;
        dsysv_rk("U", 0, -1, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("DSYSV_RK", infot, nout, lerr, ok);
        infot = 5;
        dsysv_rk("U", 2, 0, a, 1, e, ip, b, 2, w, 1, info);
        chkxer("DSYSV_RK", infot, nout, lerr, ok);
        infot = 9;
        dsysv_rk("U", 2, 0, a, 2, e, ip, b, 1, w, 1, info);
        chkxer("DSYSV_RK", infot, nout, lerr, ok);
        infot = 11;
        dsysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, 0, info);
        chkxer("DSYSV_RK", infot, nout, lerr, ok);
        infot = 11;
        dsysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, -2, info);
        chkxer("DSYSV_RK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SP") {
        //
        //        DSPSV
        //
        srnamt = "DSPSV ";
        infot = 1;
        dspsv("/", 0, 0, a, ip, b, 1, info);
        chkxer("DSPSV ", infot, nout, lerr, ok);
        infot = 2;
        dspsv("U", -1, 0, a, ip, b, 1, info);
        chkxer("DSPSV ", infot, nout, lerr, ok);
        infot = 3;
        dspsv("U", 0, -1, a, ip, b, 1, info);
        chkxer("DSPSV ", infot, nout, lerr, ok);
        infot = 7;
        dspsv("U", 2, 0, a, ip, b, 1, info);
        chkxer("DSPSV ", infot, nout, lerr, ok);
        //
        //        DSPSVX
        //
        srnamt = "DSPSVX";
        infot = 1;
        dspsvx("/", "U", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DSPSVX", infot, nout, lerr, ok);
        infot = 2;
        dspsvx("N", "/", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DSPSVX", infot, nout, lerr, ok);
        infot = 3;
        dspsvx("N", "U", -1, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DSPSVX", infot, nout, lerr, ok);
        infot = 4;
        dspsvx("N", "U", 0, -1, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DSPSVX", infot, nout, lerr, ok);
        infot = 9;
        dspsvx("N", "U", 2, 0, a, af, ip, b, 1, x, 2, rcond, r1, r2, w, iw, info);
        chkxer("DSPSVX", infot, nout, lerr, ok);
        infot = 11;
        dspsvx("N", "U", 2, 0, a, af, ip, b, 2, x, 1, rcond, r1, r2, w, iw, info);
        chkxer("DSPSVX", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    if (ok) {
        write(nout, "(1x,a3,' drivers passed the tests of the error exits')"), path;
    } else {
        write(nout, "(' *** ',a3,' drivers failed the tests of the error ','exits ***')"), path;
    }
    //
    //     End of Rerrvx
    //
}
