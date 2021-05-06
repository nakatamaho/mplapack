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

void Rerrge(const char *path, INTEGER const nunit) {
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
    char c2[2] = path[(2 - 1) + (3 - 1) * ldpath];
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 4;
    INTEGER i = 0;
    REAL a[nmax * nmax];
    REAL af[nmax * nmax];
    REAL b[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    const INTEGER lw = 3 * nmax;
    REAL w[lw];
    REAL x[nmax];
    REAL c[nmax];
    REAL r[nmax];
    INTEGER ip[nmax];
    INTEGER iw[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / castREAL(i + j);
            af[(i - 1) + (j - 1) * ldaf] = 1.0 / castREAL(i + j);
        }
        b[j - 1] = 0.0;
        r1[j - 1] = 0.0;
        r2[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
        c[j - 1] = 0.0;
        r[j - 1] = 0.0;
        ip[j - 1] = j;
        iw[j - 1] = j;
    }
    ok = true;
    //
    INTEGER info = 0;
    INTEGER n_err_bnds = 0;
    INTEGER nparams = 0;
    char eq;
    REAL rcond = 0.0;
    REAL berr = 0.0;
    REAL err_bnds_n[nmax * 3];
    REAL err_bnds_c[nmax * 3];
    REAL params[1];
    REAL anrm = 0.0;
    REAL ccond = 0.0;
    if (Mlsamen(2, c2, "GE")) {
        //
        //        Test error exits of the routines that use the LU decomposition
        //        of a general matrix.
        //
        //        Rgetrf
        //
        infot = 1;
        Rgetrf(-1, 0, a, 1, ip, info);
        chkxer("Rgetrf", infot, nout, lerr, ok);
        infot = 2;
        Rgetrf(0, -1, a, 1, ip, info);
        chkxer("Rgetrf", infot, nout, lerr, ok);
        infot = 4;
        Rgetrf(2, 1, a, 1, ip, info);
        chkxer("Rgetrf", infot, nout, lerr, ok);
        //
        //        Rgetf2
        //
        infot = 1;
        Rgetf2(-1, 0, a, 1, ip, info);
        chkxer("Rgetf2", infot, nout, lerr, ok);
        infot = 2;
        Rgetf2(0, -1, a, 1, ip, info);
        chkxer("Rgetf2", infot, nout, lerr, ok);
        infot = 4;
        Rgetf2(2, 1, a, 1, ip, info);
        chkxer("Rgetf2", infot, nout, lerr, ok);
        //
        //        Rgetri
        //
        infot = 1;
        Rgetri(-1, a, 1, ip, w, lw, info);
        chkxer("Rgetri", infot, nout, lerr, ok);
        infot = 3;
        Rgetri(2, a, 1, ip, w, lw, info);
        chkxer("Rgetri", infot, nout, lerr, ok);
        //
        //        Rgetrs
        //
        infot = 1;
        Rgetrs("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("Rgetrs", infot, nout, lerr, ok);
        infot = 2;
        Rgetrs("N", -1, 0, a, 1, ip, b, 1, info);
        chkxer("Rgetrs", infot, nout, lerr, ok);
        infot = 3;
        Rgetrs("N", 0, -1, a, 1, ip, b, 1, info);
        chkxer("Rgetrs", infot, nout, lerr, ok);
        infot = 5;
        Rgetrs("N", 2, 1, a, 1, ip, b, 2, info);
        chkxer("Rgetrs", infot, nout, lerr, ok);
        infot = 8;
        Rgetrs("N", 2, 1, a, 2, ip, b, 1, info);
        chkxer("Rgetrs", infot, nout, lerr, ok);
        //
        //        RgerFS
        //
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
        //        RgerFSX
        //
        n_err_bnds = 3;
        nparams = 0;
        infot = 1;
        Rgerfsx("/", eq, 0, 0, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RgerFSX", infot, nout, lerr, ok);
        infot = 2;
        eq = "/";
        Rgerfsx("N", eq, 2, 1, a, 1, af, 2, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RgerFSX", infot, nout, lerr, ok);
        infot = 3;
        eq = "R";
        Rgerfsx("N", eq, -1, 0, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RgerFSX", infot, nout, lerr, ok);
        infot = 4;
        Rgerfsx("N", eq, 0, -1, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RgerFSX", infot, nout, lerr, ok);
        infot = 6;
        Rgerfsx("N", eq, 2, 1, a, 1, af, 2, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RgerFSX", infot, nout, lerr, ok);
        infot = 8;
        Rgerfsx("N", eq, 2, 1, a, 2, af, 1, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RgerFSX", infot, nout, lerr, ok);
        infot = 13;
        eq = "C";
        Rgerfsx("N", eq, 2, 1, a, 2, af, 2, ip, r, c, b, 1, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RgerFSX", infot, nout, lerr, ok);
        infot = 15;
        Rgerfsx("N", eq, 2, 1, a, 2, af, 2, ip, r, c, b, 2, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RgerFSX", infot, nout, lerr, ok);
        //
        //        Rgecon
        //
        infot = 1;
        Rgecon("/", 0, a, 1, anrm, rcond, w, iw, info);
        chkxer("Rgecon", infot, nout, lerr, ok);
        infot = 2;
        Rgecon("1", -1, a, 1, anrm, rcond, w, iw, info);
        chkxer("Rgecon", infot, nout, lerr, ok);
        infot = 4;
        Rgecon("1", 2, a, 1, anrm, rcond, w, iw, info);
        chkxer("Rgecon", infot, nout, lerr, ok);
        //
        //        Rgeequ
        //
        infot = 1;
        Rgeequ(-1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgeequ", infot, nout, lerr, ok);
        infot = 2;
        Rgeequ(0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgeequ", infot, nout, lerr, ok);
        infot = 4;
        Rgeequ(2, 2, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgeequ", infot, nout, lerr, ok);
        //
        //        Rgeequb
        //
        infot = 1;
        Rgeequb(-1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgeequb", infot, nout, lerr, ok);
        infot = 2;
        Rgeequb(0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgeequb", infot, nout, lerr, ok);
        infot = 4;
        Rgeequb(2, 2, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgeequb", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "GB")) {
        //
        //        Test error exits of the routines that use the LU decomposition
        //        of a general band matrix.
        //
        //        Rgbtrf
        //
        infot = 1;
        Rgbtrf(-1, 0, 0, 0, a, 1, ip, info);
        chkxer("Rgbtrf", infot, nout, lerr, ok);
        infot = 2;
        Rgbtrf(0, -1, 0, 0, a, 1, ip, info);
        chkxer("Rgbtrf", infot, nout, lerr, ok);
        infot = 3;
        Rgbtrf(1, 1, -1, 0, a, 1, ip, info);
        chkxer("Rgbtrf", infot, nout, lerr, ok);
        infot = 4;
        Rgbtrf(1, 1, 0, -1, a, 1, ip, info);
        chkxer("Rgbtrf", infot, nout, lerr, ok);
        infot = 6;
        Rgbtrf(2, 2, 1, 1, a, 3, ip, info);
        chkxer("Rgbtrf", infot, nout, lerr, ok);
        //
        //        Rgbtf2
        //
        infot = 1;
        Rgbtf2(-1, 0, 0, 0, a, 1, ip, info);
        chkxer("Rgbtf2", infot, nout, lerr, ok);
        infot = 2;
        Rgbtf2(0, -1, 0, 0, a, 1, ip, info);
        chkxer("Rgbtf2", infot, nout, lerr, ok);
        infot = 3;
        Rgbtf2(1, 1, -1, 0, a, 1, ip, info);
        chkxer("Rgbtf2", infot, nout, lerr, ok);
        infot = 4;
        Rgbtf2(1, 1, 0, -1, a, 1, ip, info);
        chkxer("Rgbtf2", infot, nout, lerr, ok);
        infot = 6;
        Rgbtf2(2, 2, 1, 1, a, 3, ip, info);
        chkxer("Rgbtf2", infot, nout, lerr, ok);
        //
        //        Rgbtrs
        //
        infot = 1;
        Rgbtrs("/", 0, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("Rgbtrs", infot, nout, lerr, ok);
        infot = 2;
        Rgbtrs("N", -1, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("Rgbtrs", infot, nout, lerr, ok);
        infot = 3;
        Rgbtrs("N", 1, -1, 0, 1, a, 1, ip, b, 1, info);
        chkxer("Rgbtrs", infot, nout, lerr, ok);
        infot = 4;
        Rgbtrs("N", 1, 0, -1, 1, a, 1, ip, b, 1, info);
        chkxer("Rgbtrs", infot, nout, lerr, ok);
        infot = 5;
        Rgbtrs("N", 1, 0, 0, -1, a, 1, ip, b, 1, info);
        chkxer("Rgbtrs", infot, nout, lerr, ok);
        infot = 7;
        Rgbtrs("N", 2, 1, 1, 1, a, 3, ip, b, 2, info);
        chkxer("Rgbtrs", infot, nout, lerr, ok);
        infot = 10;
        Rgbtrs("N", 2, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("Rgbtrs", infot, nout, lerr, ok);
        //
        //        Rgbrfs
        //
        infot = 1;
        Rgbrfs("/", 0, 0, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rgbrfs", infot, nout, lerr, ok);
        infot = 2;
        Rgbrfs("N", -1, 0, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rgbrfs", infot, nout, lerr, ok);
        infot = 3;
        Rgbrfs("N", 1, -1, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rgbrfs", infot, nout, lerr, ok);
        infot = 4;
        Rgbrfs("N", 1, 0, -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rgbrfs", infot, nout, lerr, ok);
        infot = 5;
        Rgbrfs("N", 1, 0, 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rgbrfs", infot, nout, lerr, ok);
        infot = 7;
        Rgbrfs("N", 2, 1, 1, 1, a, 2, af, 4, ip, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("Rgbrfs", infot, nout, lerr, ok);
        infot = 9;
        Rgbrfs("N", 2, 1, 1, 1, a, 3, af, 3, ip, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("Rgbrfs", infot, nout, lerr, ok);
        infot = 12;
        Rgbrfs("N", 2, 0, 0, 1, a, 1, af, 1, ip, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("Rgbrfs", infot, nout, lerr, ok);
        infot = 14;
        Rgbrfs("N", 2, 0, 0, 1, a, 1, af, 1, ip, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("Rgbrfs", infot, nout, lerr, ok);
        //
        //        Rgbrfsx
        //
        n_err_bnds = 3;
        nparams = 0;
        infot = 1;
        Rgbrfsx("/", eq, 0, 0, 0, 0, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("Rgbrfsx", infot, nout, lerr, ok);
        infot = 2;
        eq = "/";
        Rgbrfsx("N", eq, 2, 1, 1, 1, a, 1, af, 2, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("Rgbrfsx", infot, nout, lerr, ok);
        infot = 3;
        eq = "R";
        Rgbrfsx("N", eq, -1, 1, 1, 0, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("Rgbrfsx", infot, nout, lerr, ok);
        infot = 4;
        eq = "R";
        Rgbrfsx("N", eq, 2, -1, 1, 1, a, 3, af, 4, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("Rgbrfsx", infot, nout, lerr, ok);
        infot = 5;
        eq = "R";
        Rgbrfsx("N", eq, 2, 1, -1, 1, a, 3, af, 4, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("Rgbrfsx", infot, nout, lerr, ok);
        infot = 6;
        Rgbrfsx("N", eq, 0, 0, 0, -1, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("Rgbrfsx", infot, nout, lerr, ok);
        infot = 8;
        Rgbrfsx("N", eq, 2, 1, 1, 1, a, 1, af, 2, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("Rgbrfsx", infot, nout, lerr, ok);
        infot = 10;
        Rgbrfsx("N", eq, 2, 1, 1, 1, a, 3, af, 3, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("Rgbrfsx", infot, nout, lerr, ok);
        infot = 13;
        eq = "C";
        Rgbrfsx("N", eq, 2, 1, 1, 1, a, 3, af, 5, ip, r, c, b, 1, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("Rgbrfsx", infot, nout, lerr, ok);
        infot = 15;
        Rgbrfsx("N", eq, 2, 1, 1, 1, a, 3, af, 5, ip, r, c, b, 2, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("Rgbrfsx", infot, nout, lerr, ok);
        //
        //        Rgbcon
        //
        infot = 1;
        Rgbcon("/", 0, 0, 0, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("Rgbcon", infot, nout, lerr, ok);
        infot = 2;
        Rgbcon("1", -1, 0, 0, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("Rgbcon", infot, nout, lerr, ok);
        infot = 3;
        Rgbcon("1", 1, -1, 0, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("Rgbcon", infot, nout, lerr, ok);
        infot = 4;
        Rgbcon("1", 1, 0, -1, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("Rgbcon", infot, nout, lerr, ok);
        infot = 6;
        Rgbcon("1", 2, 1, 1, a, 3, ip, anrm, rcond, w, iw, info);
        chkxer("Rgbcon", infot, nout, lerr, ok);
        //
        //        Rgbequ
        //
        infot = 1;
        Rgbequ(-1, 0, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgbequ", infot, nout, lerr, ok);
        infot = 2;
        Rgbequ(0, -1, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgbequ", infot, nout, lerr, ok);
        infot = 3;
        Rgbequ(1, 1, -1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgbequ", infot, nout, lerr, ok);
        infot = 4;
        Rgbequ(1, 1, 0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgbequ", infot, nout, lerr, ok);
        infot = 6;
        Rgbequ(2, 2, 1, 1, a, 2, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgbequ", infot, nout, lerr, ok);
        //
        //        Rgbequb
        //
        infot = 1;
        Rgbequb(-1, 0, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgbequb", infot, nout, lerr, ok);
        infot = 2;
        Rgbequb(0, -1, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgbequb", infot, nout, lerr, ok);
        infot = 3;
        Rgbequb(1, 1, -1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgbequb", infot, nout, lerr, ok);
        infot = 4;
        Rgbequb(1, 1, 0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgbequb", infot, nout, lerr, ok);
        infot = 6;
        Rgbequb(2, 2, 1, 1, a, 2, r1, r2, rcond, ccond, anrm, info);
        chkxer("Rgbequb", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrge
    //
}
