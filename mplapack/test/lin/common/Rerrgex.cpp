/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine DERRGEX.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

void Rerrge(fem::str_cref path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    nout = nunit;
    write(nout, star);
    fem::str<2> c2 = path(2, 3);
    //
    // Set the variables to innocuous values.
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
            a[(i - 1) + (j - 1) * nmax] = 1.0 / castREAL(i + j);
            af[(i - 1) + (j - 1) * nmax] = 1.0 / castREAL(i + j);
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
    fem::str<1> eq;
    REAL rcond = 0.0;
    REAL berr = 0.0;
    REAL err_bnds_n[nmax * 3];
    REAL err_bnds_c[nmax * 3];
    REAL params[1];
    REAL anrm = 0.0;
    REAL ccond = 0.0;
    if (Mlsamen(2, c2.elems, "GE")) {
        //
        // Test error exits of the routines that use the LU decomposition
        // of a general matrix.
        //
        // Rgetrf
        //
        srnamt = "DGETRF";
        infot = 1;
        Rgetrf(-1, 0, a, 1, ip, info);
        Chkxer("DGETRF", infot, nout, lerr, ok);
        infot = 2;
        Rgetrf(0, -1, a, 1, ip, info);
        Chkxer("DGETRF", infot, nout, lerr, ok);
        infot = 4;
        Rgetrf(2, 1, a, 1, ip, info);
        Chkxer("DGETRF", infot, nout, lerr, ok);
        //
        // Rgetf2
        //
        srnamt = "DGETF2";
        infot = 1;
        Rgetf2(-1, 0, a, 1, ip, info);
        Chkxer("DGETF2", infot, nout, lerr, ok);
        infot = 2;
        Rgetf2(0, -1, a, 1, ip, info);
        Chkxer("DGETF2", infot, nout, lerr, ok);
        infot = 4;
        Rgetf2(2, 1, a, 1, ip, info);
        Chkxer("DGETF2", infot, nout, lerr, ok);
        //
        // Rgetri
        //
        srnamt = "DGETRI";
        infot = 1;
        Rgetri(-1, a, 1, ip, w, lw, info);
        Chkxer("DGETRI", infot, nout, lerr, ok);
        infot = 3;
        Rgetri(2, a, 1, ip, w, lw, info);
        Chkxer("DGETRI", infot, nout, lerr, ok);
        //
        // Rgetrs
        //
        srnamt = "DGETRS";
        infot = 1;
        Rgetrs("/", 0, 0, a, 1, ip, b, 1, info);
        Chkxer("DGETRS", infot, nout, lerr, ok);
        infot = 2;
        Rgetrs("N", -1, 0, a, 1, ip, b, 1, info);
        Chkxer("DGETRS", infot, nout, lerr, ok);
        infot = 3;
        Rgetrs("N", 0, -1, a, 1, ip, b, 1, info);
        Chkxer("DGETRS", infot, nout, lerr, ok);
        infot = 5;
        Rgetrs("N", 2, 1, a, 1, ip, b, 2, info);
        Chkxer("DGETRS", infot, nout, lerr, ok);
        infot = 8;
        Rgetrs("N", 2, 1, a, 2, ip, b, 1, info);
        Chkxer("DGETRS", infot, nout, lerr, ok);
        //
        // Rgerfs
        //
        srnamt = "DGERFS";
        infot = 1;
        Rgerfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DGERFS", infot, nout, lerr, ok);
        infot = 2;
        Rgerfs("N", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DGERFS", infot, nout, lerr, ok);
        infot = 3;
        Rgerfs("N", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DGERFS", infot, nout, lerr, ok);
        infot = 5;
        Rgerfs("N", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, iw, info);
        Chkxer("DGERFS", infot, nout, lerr, ok);
        infot = 7;
        Rgerfs("N", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, iw, info);
        Chkxer("DGERFS", infot, nout, lerr, ok);
        infot = 10;
        Rgerfs("N", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, iw, info);
        Chkxer("DGERFS", infot, nout, lerr, ok);
        infot = 12;
        Rgerfs("N", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, iw, info);
        Chkxer("DGERFS", infot, nout, lerr, ok);
        //
        // Rgerfsx
        //
        n_err_bnds = 3;
        nparams = 0;
        srnamt = "DGERFSX";
        infot = 1;
        Rgerfsx("/", eq, 0, 0, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGERFSX", infot, nout, lerr, ok);
        infot = 2;
        eq = "/";
        Rgerfsx("N", eq, 2, 1, a, 1, af, 2, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGERFSX", infot, nout, lerr, ok);
        infot = 3;
        eq = "R";
        Rgerfsx("N", eq, -1, 0, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGERFSX", infot, nout, lerr, ok);
        infot = 4;
        Rgerfsx("N", eq, 0, -1, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGERFSX", infot, nout, lerr, ok);
        infot = 6;
        Rgerfsx("N", eq, 2, 1, a, 1, af, 2, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGERFSX", infot, nout, lerr, ok);
        infot = 8;
        Rgerfsx("N", eq, 2, 1, a, 2, af, 1, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGERFSX", infot, nout, lerr, ok);
        infot = 13;
        eq = "C";
        Rgerfsx("N", eq, 2, 1, a, 2, af, 2, ip, r, c, b, 1, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGERFSX", infot, nout, lerr, ok);
        infot = 15;
        Rgerfsx("N", eq, 2, 1, a, 2, af, 2, ip, r, c, b, 2, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGERFSX", infot, nout, lerr, ok);
        //
        // Rgecon
        //
        srnamt = "DGECON";
        infot = 1;
        Rgecon("/", 0, a, 1, anrm, rcond, w, iw, info);
        Chkxer("DGECON", infot, nout, lerr, ok);
        infot = 2;
        Rgecon("1", -1, a, 1, anrm, rcond, w, iw, info);
        Chkxer("DGECON", infot, nout, lerr, ok);
        infot = 4;
        Rgecon("1", 2, a, 1, anrm, rcond, w, iw, info);
        Chkxer("DGECON", infot, nout, lerr, ok);
        //
        // Rgeequ
        //
        srnamt = "DGEEQU";
        infot = 1;
        Rgeequ(-1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGEEQU", infot, nout, lerr, ok);
        infot = 2;
        Rgeequ(0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGEEQU", infot, nout, lerr, ok);
        infot = 4;
        Rgeequ(2, 2, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGEEQU", infot, nout, lerr, ok);
        //
        // Rgeequb
        //
        srnamt = "DGEEQUB";
        infot = 1;
        Rgeequb(-1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGEEQUB", infot, nout, lerr, ok);
        infot = 2;
        Rgeequb(0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGEEQUB", infot, nout, lerr, ok);
        infot = 4;
        Rgeequb(2, 2, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGEEQUB", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "GB")) {
        //
        // Test error exits of the routines that use the LU decomposition
        // of a general band matrix.
        //
        // Rgbtrf
        //
        srnamt = "DGBTRF";
        infot = 1;
        Rgbtrf(-1, 0, 0, 0, a, 1, ip, info);
        Chkxer("DGBTRF", infot, nout, lerr, ok);
        infot = 2;
        Rgbtrf(0, -1, 0, 0, a, 1, ip, info);
        Chkxer("DGBTRF", infot, nout, lerr, ok);
        infot = 3;
        Rgbtrf(1, 1, -1, 0, a, 1, ip, info);
        Chkxer("DGBTRF", infot, nout, lerr, ok);
        infot = 4;
        Rgbtrf(1, 1, 0, -1, a, 1, ip, info);
        Chkxer("DGBTRF", infot, nout, lerr, ok);
        infot = 6;
        Rgbtrf(2, 2, 1, 1, a, 3, ip, info);
        Chkxer("DGBTRF", infot, nout, lerr, ok);
        //
        // Rgbtf2
        //
        srnamt = "DGBTF2";
        infot = 1;
        Rgbtf2(-1, 0, 0, 0, a, 1, ip, info);
        Chkxer("DGBTF2", infot, nout, lerr, ok);
        infot = 2;
        Rgbtf2(0, -1, 0, 0, a, 1, ip, info);
        Chkxer("DGBTF2", infot, nout, lerr, ok);
        infot = 3;
        Rgbtf2(1, 1, -1, 0, a, 1, ip, info);
        Chkxer("DGBTF2", infot, nout, lerr, ok);
        infot = 4;
        Rgbtf2(1, 1, 0, -1, a, 1, ip, info);
        Chkxer("DGBTF2", infot, nout, lerr, ok);
        infot = 6;
        Rgbtf2(2, 2, 1, 1, a, 3, ip, info);
        Chkxer("DGBTF2", infot, nout, lerr, ok);
        //
        // Rgbtrs
        //
        srnamt = "DGBTRS";
        infot = 1;
        Rgbtrs("/", 0, 0, 0, 1, a, 1, ip, b, 1, info);
        Chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 2;
        Rgbtrs("N", -1, 0, 0, 1, a, 1, ip, b, 1, info);
        Chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 3;
        Rgbtrs("N", 1, -1, 0, 1, a, 1, ip, b, 1, info);
        Chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 4;
        Rgbtrs("N", 1, 0, -1, 1, a, 1, ip, b, 1, info);
        Chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 5;
        Rgbtrs("N", 1, 0, 0, -1, a, 1, ip, b, 1, info);
        Chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 7;
        Rgbtrs("N", 2, 1, 1, 1, a, 3, ip, b, 2, info);
        Chkxer("DGBTRS", infot, nout, lerr, ok);
        infot = 10;
        Rgbtrs("N", 2, 0, 0, 1, a, 1, ip, b, 1, info);
        Chkxer("DGBTRS", infot, nout, lerr, ok);
        //
        // Rgbrfs
        //
        srnamt = "DGBRFS";
        infot = 1;
        Rgbrfs("/", 0, 0, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 2;
        Rgbrfs("N", -1, 0, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 3;
        Rgbrfs("N", 1, -1, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 4;
        Rgbrfs("N", 1, 0, -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 5;
        Rgbrfs("N", 1, 0, 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 7;
        Rgbrfs("N", 2, 1, 1, 1, a, 2, af, 4, ip, b, 2, x, 2, r1, r2, w, iw, info);
        Chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 9;
        Rgbrfs("N", 2, 1, 1, 1, a, 3, af, 3, ip, b, 2, x, 2, r1, r2, w, iw, info);
        Chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 12;
        Rgbrfs("N", 2, 0, 0, 1, a, 1, af, 1, ip, b, 1, x, 2, r1, r2, w, iw, info);
        Chkxer("DGBRFS", infot, nout, lerr, ok);
        infot = 14;
        Rgbrfs("N", 2, 0, 0, 1, a, 1, af, 1, ip, b, 2, x, 1, r1, r2, w, iw, info);
        Chkxer("DGBRFS", infot, nout, lerr, ok);
        //
        // Rgbrfsx
        //
        n_err_bnds = 3;
        nparams = 0;
        srnamt = "DGBRFSX";
        infot = 1;
        Rgbrfsx("/", eq, 0, 0, 0, 0, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGBRFSX", infot, nout, lerr, ok);
        infot = 2;
        eq = "/";
        Rgbrfsx("N", eq, 2, 1, 1, 1, a, 1, af, 2, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGBRFSX", infot, nout, lerr, ok);
        infot = 3;
        eq = "R";
        Rgbrfsx("N", eq, -1, 1, 1, 0, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGBRFSX", infot, nout, lerr, ok);
        infot = 4;
        eq = "R";
        Rgbrfsx("N", eq, 2, -1, 1, 1, a, 3, af, 4, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGBRFSX", infot, nout, lerr, ok);
        infot = 5;
        eq = "R";
        Rgbrfsx("N", eq, 2, 1, -1, 1, a, 3, af, 4, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGBRFSX", infot, nout, lerr, ok);
        infot = 6;
        Rgbrfsx("N", eq, 0, 0, 0, -1, a, 1, af, 1, ip, r, c, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGBRFSX", infot, nout, lerr, ok);
        infot = 8;
        Rgbrfsx("N", eq, 2, 1, 1, 1, a, 1, af, 2, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGBRFSX", infot, nout, lerr, ok);
        infot = 10;
        Rgbrfsx("N", eq, 2, 1, 1, 1, a, 3, af, 3, ip, r, c, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGBRFSX", infot, nout, lerr, ok);
        infot = 13;
        eq = "C";
        Rgbrfsx("N", eq, 2, 1, 1, 1, a, 3, af, 5, ip, r, c, b, 1, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGBRFSX", infot, nout, lerr, ok);
        infot = 15;
        Rgbrfsx("N", eq, 2, 1, 1, 1, a, 3, af, 5, ip, r, c, b, 2, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        Chkxer("DGBRFSX", infot, nout, lerr, ok);
        //
        // Rgbcon
        //
        srnamt = "DGBCON";
        infot = 1;
        Rgbcon("/", 0, 0, 0, a, 1, ip, anrm, rcond, w, iw, info);
        Chkxer("DGBCON", infot, nout, lerr, ok);
        infot = 2;
        Rgbcon("1", -1, 0, 0, a, 1, ip, anrm, rcond, w, iw, info);
        Chkxer("DGBCON", infot, nout, lerr, ok);
        infot = 3;
        Rgbcon("1", 1, -1, 0, a, 1, ip, anrm, rcond, w, iw, info);
        Chkxer("DGBCON", infot, nout, lerr, ok);
        infot = 4;
        Rgbcon("1", 1, 0, -1, a, 1, ip, anrm, rcond, w, iw, info);
        Chkxer("DGBCON", infot, nout, lerr, ok);
        infot = 6;
        Rgbcon("1", 2, 1, 1, a, 3, ip, anrm, rcond, w, iw, info);
        Chkxer("DGBCON", infot, nout, lerr, ok);
        //
        // Rgbequ
        //
        srnamt = "DGBEQU";
        infot = 1;
        Rgbequ(-1, 0, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGBEQU", infot, nout, lerr, ok);
        infot = 2;
        Rgbequ(0, -1, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGBEQU", infot, nout, lerr, ok);
        infot = 3;
        Rgbequ(1, 1, -1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGBEQU", infot, nout, lerr, ok);
        infot = 4;
        Rgbequ(1, 1, 0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGBEQU", infot, nout, lerr, ok);
        infot = 6;
        Rgbequ(2, 2, 1, 1, a, 2, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGBEQU", infot, nout, lerr, ok);
        //
        // Rgbequb
        //
        srnamt = "DGBEQUB";
        infot = 1;
        Rgbequb(-1, 0, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGBEQUB", infot, nout, lerr, ok);
        infot = 2;
        Rgbequb(0, -1, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGBEQUB", infot, nout, lerr, ok);
        infot = 3;
        Rgbequb(1, 1, -1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGBEQUB", infot, nout, lerr, ok);
        infot = 4;
        Rgbequb(1, 1, 0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGBEQUB", infot, nout, lerr, ok);
        infot = 6;
        Rgbequb(2, 2, 1, 1, a, 2, r1, r2, rcond, ccond, anrm, info);
        Chkxer("DGBEQUB", infot, nout, lerr, ok);
    }
    //
    // Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    // End of Rerrge
    //
}
