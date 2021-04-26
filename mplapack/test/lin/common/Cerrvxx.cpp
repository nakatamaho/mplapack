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

void Cerrvx(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_1d<nmax, COMPLEX> e(fill0);
    arr_1d<nmax, REAL> r1(fill0);
    arr_1d<nmax, REAL> r2(fill0);
    arr_1d<2 * nmax, COMPLEX> w(fill0);
    arr_1d<nmax, COMPLEX> x(fill0);
    arr_1d<nmax, REAL> c(fill0);
    arr_1d<nmax, REAL> r(fill0);
    arr_1d<nmax, int> ip(fill0);
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = COMPLEX(1.0 / (i + j).real(), -1.0 / (i + j).real());
            af[(i - 1) + (j - 1) * ldaf] = COMPLEX(1.0 / (i + j).real(), -1.0 / (i + j).real());
        }
        b[j - 1] = 0.0;
        e[j - 1] = 0.0;
        r1[j - 1] = 0.0;
        r2[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
        c[j - 1] = 0.0;
        r[j - 1] = 0.0;
        ip[j - 1] = j;
    }
    char eq = " ";
    ok = true;
    //
    INTEGER info = 0;
    REAL rcond = 0.0;
    arr_1d<nmax, REAL> rw(fill0);
    INTEGER n_err_bnds = 0;
    INTEGER nparams = 0;
    REAL rpvgrw = 0.0;
    REAL berr = 0.0;
    arr_2d<nmax, 3, REAL> err_bnds_n(fill0);
    arr_2d<nmax, 3, REAL> err_bnds_c(fill0);
    arr_1d<1, REAL> params(fill0);
    arr_1d<nmax, REAL> rf(fill0);
    const float one = 1.0;
    if (Mlsamen2, c2, "GE") {
        //
        //        ZGESV
        //
        srnamt = "ZGESV ";
        infot = 1;
        zgesv(-1, 0, a, 1, ip, b, 1, info);
        chkxer("ZGESV ", infot, nout, lerr, ok);
        infot = 2;
        zgesv(0, -1, a, 1, ip, b, 1, info);
        chkxer("ZGESV ", infot, nout, lerr, ok);
        infot = 4;
        zgesv(2, 1, a, 1, ip, b, 2, info);
        chkxer("ZGESV ", infot, nout, lerr, ok);
        infot = 7;
        zgesv(2, 1, a, 2, ip, b, 1, info);
        chkxer("ZGESV ", infot, nout, lerr, ok);
        //
        //        ZGESVX
        //
        srnamt = "ZGESVX";
        infot = 1;
        zgesvx("/", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 2;
        zgesvx("N", "/", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 3;
        zgesvx("N", "N", -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 4;
        zgesvx("N", "N", 0, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 6;
        zgesvx("N", "N", 2, 1, a, 1, af, 2, ip, eq, r, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 8;
        zgesvx("N", "N", 2, 1, a, 2, af, 1, ip, eq, r, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 10;
        eq = "/";
        zgesvx("F", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 11;
        eq = "R";
        zgesvx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 12;
        eq = "C";
        zgesvx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 14;
        zgesvx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 16;
        zgesvx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGESVX", infot, nout, lerr, ok);
        //
        //        ZGESVXX
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "ZGESVXX";
        infot = 1;
        zgesvxx("/", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 2;
        zgesvxx("N", "/", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 3;
        zgesvxx("N", "N", -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 4;
        zgesvxx("N", "N", 0, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 6;
        zgesvxx("N", "N", 2, 1, a, 1, af, 2, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 8;
        zgesvxx("N", "N", 2, 1, a, 2, af, 1, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 10;
        eq = "/";
        zgesvxx("F", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 11;
        eq = "R";
        zgesvxx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 12;
        eq = "C";
        zgesvxx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 14;
        zgesvxx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 16;
        zgesvxx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGESVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "GB") {
        //
        //        ZGBSV
        //
        srnamt = "ZGBSV ";
        infot = 1;
        zgbsv(-1, 0, 0, 0, a, 1, ip, b, 1, info);
        chkxer("ZGBSV ", infot, nout, lerr, ok);
        infot = 2;
        zgbsv(1, -1, 0, 0, a, 1, ip, b, 1, info);
        chkxer("ZGBSV ", infot, nout, lerr, ok);
        infot = 3;
        zgbsv(1, 0, -1, 0, a, 1, ip, b, 1, info);
        chkxer("ZGBSV ", infot, nout, lerr, ok);
        infot = 4;
        zgbsv(0, 0, 0, -1, a, 1, ip, b, 1, info);
        chkxer("ZGBSV ", infot, nout, lerr, ok);
        infot = 6;
        zgbsv(1, 1, 1, 0, a, 3, ip, b, 1, info);
        chkxer("ZGBSV ", infot, nout, lerr, ok);
        infot = 9;
        zgbsv(2, 0, 0, 0, a, 1, ip, b, 1, info);
        chkxer("ZGBSV ", infot, nout, lerr, ok);
        //
        //        ZGBSVX
        //
        srnamt = "ZGBSVX";
        infot = 1;
        zgbsvx("/", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 2;
        zgbsvx("N", "/", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 3;
        zgbsvx("N", "N", -1, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 4;
        zgbsvx("N", "N", 1, -1, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 5;
        zgbsvx("N", "N", 1, 0, -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 6;
        zgbsvx("N", "N", 0, 0, 0, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 8;
        zgbsvx("N", "N", 1, 1, 1, 0, a, 2, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 10;
        zgbsvx("N", "N", 1, 1, 1, 0, a, 3, af, 3, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 12;
        eq = "/";
        zgbsvx("F", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 13;
        eq = "R";
        zgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 14;
        eq = "C";
        zgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 16;
        zgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 18;
        zgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGBSVX", infot, nout, lerr, ok);
        //
        //        ZGBSVXX
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "ZGBSVXX";
        infot = 1;
        zgbsvxx("/", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 2;
        zgbsvxx("N", "/", 0, 1, 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 3;
        zgbsvxx("N", "N", -1, 1, 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 4;
        zgbsvxx("N", "N", 2, -1, 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 5;
        zgbsvxx("N", "N", 2, 1, -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 6;
        zgbsvxx("N", "N", 0, 1, 1, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 8;
        zgbsvxx("N", "N", 2, 1, 1, 1, a, 2, af, 2, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 10;
        zgbsvxx("N", "N", 2, 1, 1, 1, a, 3, af, 3, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 12;
        eq = "/";
        zgbsvxx("F", "N", 0, 1, 1, 0, a, 3, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 13;
        eq = "R";
        zgbsvxx("F", "N", 1, 1, 1, 0, a, 3, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 14;
        eq = "C";
        zgbsvxx("F", "N", 1, 1, 1, 0, a, 3, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 15;
        zgbsvxx("N", "N", 2, 1, 1, 1, a, 3, af, 4, ip, eq, r, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 16;
        zgbsvxx("N", "N", 2, 1, 1, 1, a, 3, af, 4, ip, eq, r, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZGBSVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "GT") {
        //
        //        ZGTSV
        //
        srnamt = "ZGTSV ";
        infot = 1;
        zgtsv(-1, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], b, 1, info);
        chkxer("ZGTSV ", infot, nout, lerr, ok);
        infot = 2;
        zgtsv(0, -1, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], b, 1, info);
        chkxer("ZGTSV ", infot, nout, lerr, ok);
        infot = 7;
        zgtsv(2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], b, 1, info);
        chkxer("ZGTSV ", infot, nout, lerr, ok);
        //
        //        ZGTSVX
        //
        srnamt = "ZGTSVX";
        infot = 1;
        zgtsvx("/", "N", 0, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGTSVX", infot, nout, lerr, ok);
        infot = 2;
        zgtsvx("N", "/", 0, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGTSVX", infot, nout, lerr, ok);
        infot = 3;
        zgtsvx("N", "N", -1, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGTSVX", infot, nout, lerr, ok);
        infot = 4;
        zgtsvx("N", "N", 0, -1, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGTSVX", infot, nout, lerr, ok);
        infot = 14;
        zgtsvx("N", "N", 2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZGTSVX", infot, nout, lerr, ok);
        infot = 16;
        zgtsvx("N", "N", 2, 0, &a[(1 - 1)], &a[(2 - 1) * lda], &a[(3 - 1) * lda], af[(1 - 1)], af[(2 - 1) * ldaf], af[(3 - 1) * ldaf], af[(4 - 1) * ldaf], ip, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZGTSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "HR") {
        //
        //        ZHESV_ROOK
        //
        srnamt = "ZHESV_ROOK";
        infot = 1;
        zhesv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zhesv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 3;
        zhesv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 8;
        zhesv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PO") {
        //
        //        ZPOSV
        //
        srnamt = "ZPOSV ";
        infot = 1;
        zposv("/", 0, 0, a, 1, b, 1, info);
        chkxer("ZPOSV ", infot, nout, lerr, ok);
        infot = 2;
        zposv("U", -1, 0, a, 1, b, 1, info);
        chkxer("ZPOSV ", infot, nout, lerr, ok);
        infot = 3;
        zposv("U", 0, -1, a, 1, b, 1, info);
        chkxer("ZPOSV ", infot, nout, lerr, ok);
        infot = 5;
        zposv("U", 2, 0, a, 1, b, 2, info);
        chkxer("ZPOSV ", infot, nout, lerr, ok);
        infot = 7;
        zposv("U", 2, 0, a, 2, b, 1, info);
        chkxer("ZPOSV ", infot, nout, lerr, ok);
        //
        //        ZPOSVX
        //
        srnamt = "ZPOSVX";
        infot = 1;
        zposvx("/", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 2;
        zposvx("N", "/", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 3;
        zposvx("N", "U", -1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 4;
        zposvx("N", "U", 0, -1, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 6;
        zposvx("N", "U", 2, 0, a, 1, af, 2, eq, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 8;
        zposvx("N", "U", 2, 0, a, 2, af, 1, eq, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 9;
        eq = "/";
        zposvx("F", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 10;
        eq = "Y";
        zposvx("F", "U", 1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 12;
        zposvx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 14;
        zposvx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPOSVX", infot, nout, lerr, ok);
        //
        //        ZPOSVXX
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "ZPOSVXX";
        infot = 1;
        zposvxx("/", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 2;
        zposvxx("N", "/", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 3;
        zposvxx("N", "U", -1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 4;
        zposvxx("N", "U", 0, -1, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 6;
        zposvxx("N", "U", 2, 0, a, 1, af, 2, eq, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 8;
        zposvxx("N", "U", 2, 0, a, 2, af, 1, eq, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 9;
        eq = "/";
        zposvxx("F", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 10;
        eq = "Y";
        zposvxx("F", "U", 1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 12;
        zposvxx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 14;
        zposvxx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZPOSVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PP") {
        //
        //        ZPPSV
        //
        srnamt = "ZPPSV ";
        infot = 1;
        zppsv("/", 0, 0, a, b, 1, info);
        chkxer("ZPPSV ", infot, nout, lerr, ok);
        infot = 2;
        zppsv("U", -1, 0, a, b, 1, info);
        chkxer("ZPPSV ", infot, nout, lerr, ok);
        infot = 3;
        zppsv("U", 0, -1, a, b, 1, info);
        chkxer("ZPPSV ", infot, nout, lerr, ok);
        infot = 6;
        zppsv("U", 2, 0, a, b, 1, info);
        chkxer("ZPPSV ", infot, nout, lerr, ok);
        //
        //        ZPPSVX
        //
        srnamt = "ZPPSVX";
        infot = 1;
        zppsvx("/", "U", 0, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 2;
        zppsvx("N", "/", 0, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 3;
        zppsvx("N", "U", -1, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 4;
        zppsvx("N", "U", 0, -1, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 7;
        eq = "/";
        zppsvx("F", "U", 0, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 8;
        eq = "Y";
        zppsvx("F", "U", 1, 0, a, af, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 10;
        zppsvx("N", "U", 2, 0, a, af, eq, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 12;
        zppsvx("N", "U", 2, 0, a, af, eq, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPPSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PB") {
        //
        //        ZPBSV
        //
        srnamt = "ZPBSV ";
        infot = 1;
        zpbsv("/", 0, 0, 0, a, 1, b, 1, info);
        chkxer("ZPBSV ", infot, nout, lerr, ok);
        infot = 2;
        zpbsv("U", -1, 0, 0, a, 1, b, 1, info);
        chkxer("ZPBSV ", infot, nout, lerr, ok);
        infot = 3;
        zpbsv("U", 1, -1, 0, a, 1, b, 1, info);
        chkxer("ZPBSV ", infot, nout, lerr, ok);
        infot = 4;
        zpbsv("U", 0, 0, -1, a, 1, b, 1, info);
        chkxer("ZPBSV ", infot, nout, lerr, ok);
        infot = 6;
        zpbsv("U", 1, 1, 0, a, 1, b, 2, info);
        chkxer("ZPBSV ", infot, nout, lerr, ok);
        infot = 8;
        zpbsv("U", 2, 0, 0, a, 1, b, 1, info);
        chkxer("ZPBSV ", infot, nout, lerr, ok);
        //
        //        ZPBSVX
        //
        srnamt = "ZPBSVX";
        infot = 1;
        zpbsvx("/", "U", 0, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 2;
        zpbsvx("N", "/", 0, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 3;
        zpbsvx("N", "U", -1, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 4;
        zpbsvx("N", "U", 1, -1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 5;
        zpbsvx("N", "U", 0, 0, -1, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 7;
        zpbsvx("N", "U", 1, 1, 0, a, 1, af, 2, eq, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 9;
        zpbsvx("N", "U", 1, 1, 0, a, 2, af, 1, eq, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 10;
        eq = "/";
        zpbsvx("F", "U", 0, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 11;
        eq = "Y";
        zpbsvx("F", "U", 1, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 13;
        zpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, eq, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 15;
        zpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, eq, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPBSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "PT") {
        //
        //        ZPTSV
        //
        srnamt = "ZPTSV ";
        infot = 1;
        zptsv(-1, 0, r, &a[(1 - 1)], b, 1, info);
        chkxer("ZPTSV ", infot, nout, lerr, ok);
        infot = 2;
        zptsv(0, -1, r, &a[(1 - 1)], b, 1, info);
        chkxer("ZPTSV ", infot, nout, lerr, ok);
        infot = 6;
        zptsv(2, 0, r, &a[(1 - 1)], b, 1, info);
        chkxer("ZPTSV ", infot, nout, lerr, ok);
        //
        //        ZPTSVX
        //
        srnamt = "ZPTSVX";
        infot = 1;
        zptsvx("/", 0, 0, r, &a[(1 - 1)], rf, af[(1 - 1)], b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPTSVX", infot, nout, lerr, ok);
        infot = 2;
        zptsvx("N", -1, 0, r, &a[(1 - 1)], rf, af[(1 - 1)], b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPTSVX", infot, nout, lerr, ok);
        infot = 3;
        zptsvx("N", 0, -1, r, &a[(1 - 1)], rf, af[(1 - 1)], b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPTSVX", infot, nout, lerr, ok);
        infot = 9;
        zptsvx("N", 2, 0, r, &a[(1 - 1)], rf, af[(1 - 1)], b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZPTSVX", infot, nout, lerr, ok);
        infot = 11;
        zptsvx("N", 2, 0, r, &a[(1 - 1)], rf, af[(1 - 1)], b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZPTSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "HE") {
        //
        //        ZHESV
        //
        srnamt = "ZHESV ";
        infot = 1;
        zhesv("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 2;
        zhesv("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 3;
        zhesv("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 5;
        zhesv("U", 2, 0, a, 1, ip, b, 2, w, 1, info);
        chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 8;
        zhesv("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 10;
        zhesv("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 10;
        zhesv("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        chkxer("ZHESV ", infot, nout, lerr, ok);
        //
        //        ZHESVX
        //
        srnamt = "ZHESVX";
        infot = 1;
        zhesvx("/", "U", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 2;
        zhesvx("N", "/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 3;
        zhesvx("N", "U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 4;
        zhesvx("N", "U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 6;
        zhesvx("N", "U", 2, 0, a, 1, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 8;
        zhesvx("N", "U", 2, 0, a, 2, af, 1, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 11;
        zhesvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 1, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 13;
        zhesvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 1, rcond, r1, r2, w, 4, rw, info);
        chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 18;
        zhesvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 3, rw, info);
        chkxer("ZHESVX", infot, nout, lerr, ok);
        //
        //        ZHESVXX
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "ZHESVXX";
        infot = 1;
        zhesvxx("/", "U", 0, 0, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 2;
        zhesvxx("N", "/", 0, 0, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 3;
        zhesvxx("N", "U", -1, 0, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 4;
        zhesvxx("N", "U", 0, -1, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 6;
        zhesvxx("N", "U", 2, 0, a, 1, af, 2, ip, eq, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 8;
        zhesvxx("N", "U", 2, 0, a, 2, af, 1, ip, eq, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 9;
        eq = "/";
        zhesvxx("F", "U", 0, 0, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 10;
        eq = "Y";
        zhesvxx("F", "U", 1, 0, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 12;
        zhesvxx("N", "U", 2, 0, a, 2, af, 2, ip, eq, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 14;
        zhesvxx("N", "U", 2, 0, a, 2, af, 2, ip, eq, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZHESVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "HR") {
        //
        //        ZHESV_ROOK
        //
        srnamt = "ZHESV_ROOK";
        infot = 1;
        zhesv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zhesv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 3;
        zhesv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 8;
        zhesv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 10;
        zhesv_rook("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 10;
        zhesv_rook("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "HK") {
        //
        //        ZSYSV_RK
        //
        //        Test error exits of the driver that uses factorization
        //        of a Hermitian indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting with the new storage
        //        format for factors L ( or U) and D.
        //
        //        L (or U) is stored in A, diagonal of D is stored on the
        //        diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        srnamt = "ZHESV_RK";
        infot = 1;
        zhesv_rk("/", 0, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 2;
        zhesv_rk("U", -1, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 3;
        zhesv_rk("U", 0, -1, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 5;
        zhesv_rk("U", 2, 0, a, 1, e, ip, b, 2, w, 1, info);
        chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 9;
        zhesv_rk("U", 2, 0, a, 2, e, ip, b, 1, w, 1, info);
        chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 11;
        zhesv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, 0, info);
        chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 11;
        zhesv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, -2, info);
        chkxer("ZHESV_RK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "HP") {
        //
        //        ZHPSV
        //
        srnamt = "ZHPSV ";
        infot = 1;
        zhpsv("/", 0, 0, a, ip, b, 1, info);
        chkxer("ZHPSV ", infot, nout, lerr, ok);
        infot = 2;
        zhpsv("U", -1, 0, a, ip, b, 1, info);
        chkxer("ZHPSV ", infot, nout, lerr, ok);
        infot = 3;
        zhpsv("U", 0, -1, a, ip, b, 1, info);
        chkxer("ZHPSV ", infot, nout, lerr, ok);
        infot = 7;
        zhpsv("U", 2, 0, a, ip, b, 1, info);
        chkxer("ZHPSV ", infot, nout, lerr, ok);
        //
        //        ZHPSVX
        //
        srnamt = "ZHPSVX";
        infot = 1;
        zhpsvx("/", "U", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZHPSVX", infot, nout, lerr, ok);
        infot = 2;
        zhpsvx("N", "/", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZHPSVX", infot, nout, lerr, ok);
        infot = 3;
        zhpsvx("N", "U", -1, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZHPSVX", infot, nout, lerr, ok);
        infot = 4;
        zhpsvx("N", "U", 0, -1, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZHPSVX", infot, nout, lerr, ok);
        infot = 9;
        zhpsvx("N", "U", 2, 0, a, af, ip, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZHPSVX", infot, nout, lerr, ok);
        infot = 11;
        zhpsvx("N", "U", 2, 0, a, af, ip, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZHPSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SY") {
        //
        //        ZSYSV
        //
        srnamt = "ZSYSV ";
        infot = 1;
        zsysv("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZSYSV ", infot, nout, lerr, ok);
        infot = 2;
        zsysv("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZSYSV ", infot, nout, lerr, ok);
        infot = 3;
        zsysv("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZSYSV ", infot, nout, lerr, ok);
        infot = 8;
        zsysv("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("ZSYSV ", infot, nout, lerr, ok);
        infot = 10;
        zsysv("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("ZSYSV ", infot, nout, lerr, ok);
        infot = 10;
        zsysv("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        chkxer("ZSYSV ", infot, nout, lerr, ok);
        //
        //        ZSYSVX
        //
        srnamt = "ZSYSVX";
        infot = 1;
        zsysvx("/", "U", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 2;
        zsysvx("N", "/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 3;
        zsysvx("N", "U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 4;
        zsysvx("N", "U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 6;
        zsysvx("N", "U", 2, 0, a, 1, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 8;
        zsysvx("N", "U", 2, 0, a, 2, af, 1, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 11;
        zsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 1, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 13;
        zsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 1, rcond, r1, r2, w, 4, rw, info);
        chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 18;
        zsysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 3, rw, info);
        chkxer("ZSYSVX", infot, nout, lerr, ok);
        //
        //        ZSYSVXX
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "ZSYSVXX";
        infot = 1;
        eq = "N";
        zsysvxx("/", "U", 0, 0, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 2;
        zsysvxx("N", "/", 0, 0, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 3;
        zsysvxx("N", "U", -1, 0, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 4;
        eq = "/";
        zsysvxx("N", "U", 0, -1, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        eq = "Y";
        infot = 6;
        zsysvxx("N", "U", 2, 0, a, 1, af, 2, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 8;
        zsysvxx("N", "U", 2, 0, a, 2, af, 1, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 10;
        zsysvxx("F", "U", 2, 0, a, 2, af, 2, ip, "A", r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 11;
        eq = "Y";
        zsysvxx("F", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 11;
        eq = "Y";
        r[1 - 1] = -one;
        zsysvxx("F", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 13;
        eq = "N";
        zsysvxx("N", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 15;
        zsysvxx("N", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        chkxer("ZSYSVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SR") {
        //
        //        ZSYSV_ROOK
        //
        srnamt = "ZSYSV_ROOK";
        infot = 1;
        zsysv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZSYSV_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zsysv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZSYSV_ROOK", infot, nout, lerr, ok);
        infot = 3;
        zsysv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZSYSV_ROOK", infot, nout, lerr, ok);
        infot = 8;
        zsysv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("ZSYSV_ROOK", infot, nout, lerr, ok);
        infot = 10;
        zsysv_rook("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("ZSYSV_ROOK", infot, nout, lerr, ok);
        infot = 10;
        zsysv_rook("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        //
    } else if (Mlsamen2, c2, "SK") {
        //
        //        ZSYSV_RK
        //
        //        Test error exits of the driver that uses factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting with the new storage
        //        format for factors L ( or U) and D.
        //
        //        L (or U) is stored in A, diagonal of D is stored on the
        //        diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        srnamt = "ZSYSV_RK";
        infot = 1;
        zsysv_rk("/", 0, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 2;
        zsysv_rk("U", -1, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 3;
        zsysv_rk("U", 0, -1, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 5;
        zsysv_rk("U", 2, 0, a, 1, e, ip, b, 2, w, 1, info);
        chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 9;
        zsysv_rk("U", 2, 0, a, 2, e, ip, b, 1, w, 1, info);
        chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 11;
        zsysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, 0, info);
        chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 11;
        zsysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, -2, info);
        chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SP") {
        //
        //        ZSPSV
        //
        srnamt = "ZSPSV ";
        infot = 1;
        zspsv("/", 0, 0, a, ip, b, 1, info);
        chkxer("ZSPSV ", infot, nout, lerr, ok);
        infot = 2;
        zspsv("U", -1, 0, a, ip, b, 1, info);
        chkxer("ZSPSV ", infot, nout, lerr, ok);
        infot = 3;
        zspsv("U", 0, -1, a, ip, b, 1, info);
        chkxer("ZSPSV ", infot, nout, lerr, ok);
        infot = 7;
        zspsv("U", 2, 0, a, ip, b, 1, info);
        chkxer("ZSPSV ", infot, nout, lerr, ok);
        //
        //        ZSPSVX
        //
        srnamt = "ZSPSVX";
        infot = 1;
        zspsvx("/", "U", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZSPSVX", infot, nout, lerr, ok);
        infot = 2;
        zspsvx("N", "/", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZSPSVX", infot, nout, lerr, ok);
        infot = 3;
        zspsvx("N", "U", -1, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZSPSVX", infot, nout, lerr, ok);
        infot = 4;
        zspsvx("N", "U", 0, -1, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZSPSVX", infot, nout, lerr, ok);
        infot = 9;
        zspsvx("N", "U", 2, 0, a, af, ip, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("ZSPSVX", infot, nout, lerr, ok);
        infot = 11;
        zspsvx("N", "U", 2, 0, a, af, ip, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("ZSPSVX", infot, nout, lerr, ok);
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
    //     End of Cerrvx
    //
}
