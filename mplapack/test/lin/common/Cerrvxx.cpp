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

// Derived from LAPACK routine ZERRVXX.
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

void Cerrvx(fem::str_cref path, INTEGER const nunit) {
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
    COMPLEX a[nmax * nmax];
    COMPLEX af[nmax * nmax];
    COMPLEX b[nmax];
    COMPLEX e[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    COMPLEX w[2 * nmax];
    COMPLEX x[nmax];
    REAL c[nmax];
    REAL r[nmax];
    INTEGER ip[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * nmax] = COMPLEX(1.0 / castREAL(i + j), -1.0 / castREAL(i + j));
            af[(i - 1) + (j - 1) * nmax] = COMPLEX(1.0 / castREAL(i + j), -1.0 / castREAL(i + j));
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
    fem::str<1> eq = ' ';
    ok = true;
    //
    INTEGER info = 0;
    REAL rcond = 0.0;
    REAL rw[nmax];
    INTEGER n_err_bnds = 0;
    INTEGER nparams = 0;
    REAL rpvgrw = 0.0;
    REAL berr = 0.0;
    REAL err_bnds_n[nmax * 3];
    REAL err_bnds_c[nmax * 3];
    REAL params[1];
    REAL rf[nmax];
    const REAL one = 1.0;
    if (Mlsamen(2, c2.elems, "GE")) {
        //
        // Cgesv
        //
        srnamt = "ZGESV ";
        infot = 1;
        Cgesv(-1, 0, a, 1, ip, b, 1, info);
        Chkxer("ZGESV ", infot, nout, lerr, ok);
        infot = 2;
        Cgesv(0, -1, a, 1, ip, b, 1, info);
        Chkxer("ZGESV ", infot, nout, lerr, ok);
        infot = 4;
        Cgesv(2, 1, a, 1, ip, b, 2, info);
        Chkxer("ZGESV ", infot, nout, lerr, ok);
        infot = 7;
        Cgesv(2, 1, a, 2, ip, b, 1, info);
        Chkxer("ZGESV ", infot, nout, lerr, ok);
        //
        // Cgesvx
        //
        srnamt = "ZGESVX";
        infot = 1;
        Cgesvx("/", "N", 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 2;
        Cgesvx("N", "/", 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 3;
        Cgesvx("N", "N", -1, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 4;
        Cgesvx("N", "N", 0, -1, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 6;
        Cgesvx("N", "N", 2, 1, a, 1, af, 2, ip, eq.elems, r, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 8;
        Cgesvx("N", "N", 2, 1, a, 2, af, 1, ip, eq.elems, r, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 10;
        eq = '/';
        Cgesvx("F", "N", 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 11;
        eq = 'R';
        Cgesvx("F", "N", 1, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 12;
        eq = 'C';
        Cgesvx("F", "N", 1, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 14;
        Cgesvx("N", "N", 2, 1, a, 2, af, 2, ip, eq.elems, r, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        infot = 16;
        Cgesvx("N", "N", 2, 1, a, 2, af, 2, ip, eq.elems, r, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGESVX", infot, nout, lerr, ok);
        //
        // Cgesvxx
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "ZGESVXX";
        infot = 1;
        Cgesvxx("/", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 2;
        Cgesvxx("N", "/", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 3;
        Cgesvxx("N", "N", -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 4;
        Cgesvxx("N", "N", 0, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 6;
        Cgesvxx("N", "N", 2, 1, a, 1, af, 2, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 8;
        Cgesvxx("N", "N", 2, 1, a, 2, af, 1, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 10;
        eq = '/';
        Cgesvxx("F", "N", 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 11;
        eq = 'R';
        Cgesvxx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 12;
        eq = 'C';
        Cgesvxx("F", "N", 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 14;
        Cgesvxx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        infot = 16;
        Cgesvxx("N", "N", 2, 1, a, 2, af, 2, ip, eq, r, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGESVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "GB")) {
        //
        // Cgbsv
        //
        srnamt = "ZGBSV ";
        infot = 1;
        Cgbsv(-1, 0, 0, 0, a, 1, ip, b, 1, info);
        Chkxer("ZGBSV ", infot, nout, lerr, ok);
        infot = 2;
        Cgbsv(1, -1, 0, 0, a, 1, ip, b, 1, info);
        Chkxer("ZGBSV ", infot, nout, lerr, ok);
        infot = 3;
        Cgbsv(1, 0, -1, 0, a, 1, ip, b, 1, info);
        Chkxer("ZGBSV ", infot, nout, lerr, ok);
        infot = 4;
        Cgbsv(0, 0, 0, -1, a, 1, ip, b, 1, info);
        Chkxer("ZGBSV ", infot, nout, lerr, ok);
        infot = 6;
        Cgbsv(1, 1, 1, 0, a, 3, ip, b, 1, info);
        Chkxer("ZGBSV ", infot, nout, lerr, ok);
        infot = 9;
        Cgbsv(2, 0, 0, 0, a, 1, ip, b, 1, info);
        Chkxer("ZGBSV ", infot, nout, lerr, ok);
        //
        // Cgbsvx
        //
        srnamt = "ZGBSVX";
        infot = 1;
        Cgbsvx("/", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 2;
        Cgbsvx("N", "/", 0, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 3;
        Cgbsvx("N", "N", -1, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 4;
        Cgbsvx("N", "N", 1, -1, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 5;
        Cgbsvx("N", "N", 1, 0, -1, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 6;
        Cgbsvx("N", "N", 0, 0, 0, -1, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 8;
        Cgbsvx("N", "N", 1, 1, 1, 0, a, 2, af, 4, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 10;
        Cgbsvx("N", "N", 1, 1, 1, 0, a, 3, af, 3, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 12;
        eq = '/';
        Cgbsvx("F", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 13;
        eq = 'R';
        Cgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 14;
        eq = 'C';
        Cgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 16;
        Cgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        infot = 18;
        Cgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, eq.elems, r, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGBSVX", infot, nout, lerr, ok);
        //
        // Cgbsvxx
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "ZGBSVXX";
        infot = 1;
        Cgbsvxx("/", "N", 0, 0, 0, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 2;
        Cgbsvxx("N", "/", 0, 1, 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 3;
        Cgbsvxx("N", "N", -1, 1, 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 4;
        Cgbsvxx("N", "N", 2, -1, 1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 5;
        Cgbsvxx("N", "N", 2, 1, -1, 0, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 6;
        Cgbsvxx("N", "N", 0, 1, 1, -1, a, 1, af, 1, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 8;
        Cgbsvxx("N", "N", 2, 1, 1, 1, a, 2, af, 2, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 10;
        Cgbsvxx("N", "N", 2, 1, 1, 1, a, 3, af, 3, ip, eq, r, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 12;
        eq = '/';
        Cgbsvxx("F", "N", 0, 1, 1, 0, a, 3, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 13;
        eq = 'R';
        Cgbsvxx("F", "N", 1, 1, 1, 0, a, 3, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 14;
        eq = 'C';
        Cgbsvxx("F", "N", 1, 1, 1, 0, a, 3, af, 4, ip, eq, r, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 15;
        Cgbsvxx("N", "N", 2, 1, 1, 1, a, 3, af, 4, ip, eq, r, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        infot = 16;
        Cgbsvxx("N", "N", 2, 1, 1, 1, a, 3, af, 4, ip, eq, r, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZGBSVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "GT")) {
        //
        // Cgtsv
        //
        srnamt = "ZGTSV ";
        infot = 1;
        Cgtsv(-1, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], b, 1, info);
        Chkxer("ZGTSV ", infot, nout, lerr, ok);
        infot = 2;
        Cgtsv(0, -1, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], b, 1, info);
        Chkxer("ZGTSV ", infot, nout, lerr, ok);
        infot = 7;
        Cgtsv(2, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], b, 1, info);
        Chkxer("ZGTSV ", infot, nout, lerr, ok);
        //
        // Cgtsvx
        //
        srnamt = "ZGTSVX";
        infot = 1;
        Cgtsvx("/", "N", 0, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGTSVX", infot, nout, lerr, ok);
        infot = 2;
        Cgtsvx("N", "/", 0, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGTSVX", infot, nout, lerr, ok);
        infot = 3;
        Cgtsvx("N", "N", -1, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGTSVX", infot, nout, lerr, ok);
        infot = 4;
        Cgtsvx("N", "N", 0, -1, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGTSVX", infot, nout, lerr, ok);
        infot = 14;
        Cgtsvx("N", "N", 2, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZGTSVX", infot, nout, lerr, ok);
        infot = 16;
        Cgtsvx("N", "N", 2, 0, &a[0], &a[(2 - 1) * nmax], &a[(3 - 1) * nmax], &af[0], &af[(2 - 1) * nmax], &af[(3 - 1) * nmax], &af[(4 - 1) * nmax], ip, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZGTSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "HR")) {
        //
        // Chesv_rook
        //
        srnamt = "ZHESV_ROOK";
        infot = 1;
        Chesv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Chesv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 3;
        Chesv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 8;
        Chesv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "PO")) {
        //
        // Cposv
        //
        srnamt = "ZPOSV ";
        infot = 1;
        Cposv("/", 0, 0, a, 1, b, 1, info);
        Chkxer("ZPOSV ", infot, nout, lerr, ok);
        infot = 2;
        Cposv("U", -1, 0, a, 1, b, 1, info);
        Chkxer("ZPOSV ", infot, nout, lerr, ok);
        infot = 3;
        Cposv("U", 0, -1, a, 1, b, 1, info);
        Chkxer("ZPOSV ", infot, nout, lerr, ok);
        infot = 5;
        Cposv("U", 2, 0, a, 1, b, 2, info);
        Chkxer("ZPOSV ", infot, nout, lerr, ok);
        infot = 7;
        Cposv("U", 2, 0, a, 2, b, 1, info);
        Chkxer("ZPOSV ", infot, nout, lerr, ok);
        //
        // Cposvx
        //
        srnamt = "ZPOSVX";
        infot = 1;
        Cposvx("/", "U", 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 2;
        Cposvx("N", "/", 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 3;
        Cposvx("N", "U", -1, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 4;
        Cposvx("N", "U", 0, -1, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 6;
        Cposvx("N", "U", 2, 0, a, 1, af, 2, eq.elems, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 8;
        Cposvx("N", "U", 2, 0, a, 2, af, 1, eq.elems, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 9;
        eq = '/';
        Cposvx("F", "U", 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 10;
        eq = 'Y';
        Cposvx("F", "U", 1, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 12;
        Cposvx("N", "U", 2, 0, a, 2, af, 2, eq.elems, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZPOSVX", infot, nout, lerr, ok);
        infot = 14;
        Cposvx("N", "U", 2, 0, a, 2, af, 2, eq.elems, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPOSVX", infot, nout, lerr, ok);
        //
        // Cposvxx
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "ZPOSVXX";
        infot = 1;
        Cposvxx("/", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 2;
        Cposvxx("N", "/", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 3;
        Cposvxx("N", "U", -1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 4;
        Cposvxx("N", "U", 0, -1, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 6;
        Cposvxx("N", "U", 2, 0, a, 1, af, 2, eq, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 8;
        Cposvxx("N", "U", 2, 0, a, 2, af, 1, eq, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 9;
        eq = '/';
        Cposvxx("F", "U", 0, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 10;
        eq = 'Y';
        Cposvxx("F", "U", 1, 0, a, 1, af, 1, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 12;
        Cposvxx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZPOSVXX", infot, nout, lerr, ok);
        infot = 14;
        Cposvxx("N", "U", 2, 0, a, 2, af, 2, eq, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZPOSVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "PP")) {
        //
        // Cppsv
        //
        srnamt = "ZPPSV ";
        infot = 1;
        Cppsv("/", 0, 0, a, b, 1, info);
        Chkxer("ZPPSV ", infot, nout, lerr, ok);
        infot = 2;
        Cppsv("U", -1, 0, a, b, 1, info);
        Chkxer("ZPPSV ", infot, nout, lerr, ok);
        infot = 3;
        Cppsv("U", 0, -1, a, b, 1, info);
        Chkxer("ZPPSV ", infot, nout, lerr, ok);
        infot = 6;
        Cppsv("U", 2, 0, a, b, 1, info);
        Chkxer("ZPPSV ", infot, nout, lerr, ok);
        //
        // Cppsvx
        //
        srnamt = "ZPPSVX";
        infot = 1;
        Cppsvx("/", "U", 0, 0, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 2;
        Cppsvx("N", "/", 0, 0, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 3;
        Cppsvx("N", "U", -1, 0, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 4;
        Cppsvx("N", "U", 0, -1, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 7;
        eq = '/';
        Cppsvx("F", "U", 0, 0, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 8;
        eq = 'Y';
        Cppsvx("F", "U", 1, 0, a, af, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 10;
        Cppsvx("N", "U", 2, 0, a, af, eq.elems, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZPPSVX", infot, nout, lerr, ok);
        infot = 12;
        Cppsvx("N", "U", 2, 0, a, af, eq.elems, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPPSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "PB")) {
        //
        // Cpbsv
        //
        srnamt = "ZPBSV ";
        infot = 1;
        Cpbsv("/", 0, 0, 0, a, 1, b, 1, info);
        Chkxer("ZPBSV ", infot, nout, lerr, ok);
        infot = 2;
        Cpbsv("U", -1, 0, 0, a, 1, b, 1, info);
        Chkxer("ZPBSV ", infot, nout, lerr, ok);
        infot = 3;
        Cpbsv("U", 1, -1, 0, a, 1, b, 1, info);
        Chkxer("ZPBSV ", infot, nout, lerr, ok);
        infot = 4;
        Cpbsv("U", 0, 0, -1, a, 1, b, 1, info);
        Chkxer("ZPBSV ", infot, nout, lerr, ok);
        infot = 6;
        Cpbsv("U", 1, 1, 0, a, 1, b, 2, info);
        Chkxer("ZPBSV ", infot, nout, lerr, ok);
        infot = 8;
        Cpbsv("U", 2, 0, 0, a, 1, b, 1, info);
        Chkxer("ZPBSV ", infot, nout, lerr, ok);
        //
        // Cpbsvx
        //
        srnamt = "ZPBSVX";
        infot = 1;
        Cpbsvx("/", "U", 0, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 2;
        Cpbsvx("N", "/", 0, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 3;
        Cpbsvx("N", "U", -1, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 4;
        Cpbsvx("N", "U", 1, -1, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 5;
        Cpbsvx("N", "U", 0, 0, -1, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 7;
        Cpbsvx("N", "U", 1, 1, 0, a, 1, af, 2, eq.elems, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 9;
        Cpbsvx("N", "U", 1, 1, 0, a, 2, af, 1, eq.elems, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 10;
        eq = '/';
        Cpbsvx("F", "U", 0, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 11;
        eq = 'Y';
        Cpbsvx("F", "U", 1, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 13;
        Cpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, eq.elems, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        infot = 15;
        Cpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, eq.elems, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPBSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "PT")) {
        //
        // Cptsv
        //
        srnamt = "ZPTSV ";
        infot = 1;
        Cptsv(-1, 0, r, &a[0], b, 1, info);
        Chkxer("ZPTSV ", infot, nout, lerr, ok);
        infot = 2;
        Cptsv(0, -1, r, &a[0], b, 1, info);
        Chkxer("ZPTSV ", infot, nout, lerr, ok);
        infot = 6;
        Cptsv(2, 0, r, &a[0], b, 1, info);
        Chkxer("ZPTSV ", infot, nout, lerr, ok);
        //
        // Cptsvx
        //
        srnamt = "ZPTSVX";
        infot = 1;
        Cptsvx("/", 0, 0, r, &a[0], rf, &af[0], b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPTSVX", infot, nout, lerr, ok);
        infot = 2;
        Cptsvx("N", -1, 0, r, &a[0], rf, &af[0], b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPTSVX", infot, nout, lerr, ok);
        infot = 3;
        Cptsvx("N", 0, -1, r, &a[0], rf, &af[0], b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPTSVX", infot, nout, lerr, ok);
        infot = 9;
        Cptsvx("N", 2, 0, r, &a[0], rf, &af[0], b, 1, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZPTSVX", infot, nout, lerr, ok);
        infot = 11;
        Cptsvx("N", 2, 0, r, &a[0], rf, &af[0], b, 2, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZPTSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "HE")) {
        //
        // Chesv
        //
        srnamt = "ZHESV ";
        infot = 1;
        Chesv("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 2;
        Chesv("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 3;
        Chesv("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 5;
        Chesv("U", 2, 0, a, 1, ip, b, 2, w, 1, info);
        Chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 8;
        Chesv("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        Chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 10;
        Chesv("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        Chkxer("ZHESV ", infot, nout, lerr, ok);
        infot = 10;
        Chesv("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        Chkxer("ZHESV ", infot, nout, lerr, ok);
        //
        // Chesvx
        //
        srnamt = "ZHESVX";
        infot = 1;
        Chesvx("/", "U", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        Chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 2;
        Chesvx("N", "/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        Chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 3;
        Chesvx("N", "U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        Chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 4;
        Chesvx("N", "U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        Chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 6;
        Chesvx("N", "U", 2, 0, a, 1, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        Chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 8;
        Chesvx("N", "U", 2, 0, a, 2, af, 1, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        Chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 11;
        Chesvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 1, x, 2, rcond, r1, r2, w, 4, rw, info);
        Chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 13;
        Chesvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 1, rcond, r1, r2, w, 4, rw, info);
        Chkxer("ZHESVX", infot, nout, lerr, ok);
        infot = 18;
        Chesvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 3, rw, info);
        Chkxer("ZHESVX", infot, nout, lerr, ok);
        //
        // Chesvxx
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "ZHESVXX";
        infot = 1;
        Chesvxx("/", "U", 0, 0, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 2;
        Chesvxx("N", "/", 0, 0, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 3;
        Chesvxx("N", "U", -1, 0, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 4;
        Chesvxx("N", "U", 0, -1, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 6;
        Chesvxx("N", "U", 2, 0, a, 1, af, 2, ip, eq, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 8;
        Chesvxx("N", "U", 2, 0, a, 2, af, 1, ip, eq, c, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 9;
        eq = '/';
        Chesvxx("F", "U", 0, 0, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 10;
        eq = 'Y';
        Chesvxx("F", "U", 1, 0, a, 1, af, 1, ip, eq, c, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 12;
        Chesvxx("N", "U", 2, 0, a, 2, af, 2, ip, eq, c, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZHESVXX", infot, nout, lerr, ok);
        infot = 14;
        Chesvxx("N", "U", 2, 0, a, 2, af, 2, ip, eq, c, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZHESVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "HR")) {
        //
        // Chesv_rook
        //
        srnamt = "ZHESV_ROOK";
        infot = 1;
        Chesv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Chesv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 3;
        Chesv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 8;
        Chesv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 10;
        Chesv_rook("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        Chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        infot = 10;
        Chesv_rook("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        Chkxer("ZHESV_ROOK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "HK")) {
        //
        // Csysv_rk
        //
        // Test error exits of the driver that uses factorization
        // of a Hermitian indefinite matrix with rook
        // (bounded Bunch-Kaufman) pivoting with the new storage
        // format for factors L ( or U) and D.
        //
        // L (or U) is stored in A, diagonal of D is stored on the
        // diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        srnamt = "ZHESV_RK";
        infot = 1;
        Chesv_rk("/", 0, 0, a, 1, e, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 2;
        Chesv_rk("U", -1, 0, a, 1, e, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 3;
        Chesv_rk("U", 0, -1, a, 1, e, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 5;
        Chesv_rk("U", 2, 0, a, 1, e, ip, b, 2, w, 1, info);
        Chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 9;
        Chesv_rk("U", 2, 0, a, 2, e, ip, b, 1, w, 1, info);
        Chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 11;
        Chesv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, 0, info);
        Chkxer("ZHESV_RK", infot, nout, lerr, ok);
        infot = 11;
        Chesv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, -2, info);
        Chkxer("ZHESV_RK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "HP")) {
        //
        // Chpsv
        //
        srnamt = "ZHPSV ";
        infot = 1;
        Chpsv("/", 0, 0, a, ip, b, 1, info);
        Chkxer("ZHPSV ", infot, nout, lerr, ok);
        infot = 2;
        Chpsv("U", -1, 0, a, ip, b, 1, info);
        Chkxer("ZHPSV ", infot, nout, lerr, ok);
        infot = 3;
        Chpsv("U", 0, -1, a, ip, b, 1, info);
        Chkxer("ZHPSV ", infot, nout, lerr, ok);
        infot = 7;
        Chpsv("U", 2, 0, a, ip, b, 1, info);
        Chkxer("ZHPSV ", infot, nout, lerr, ok);
        //
        // Chpsvx
        //
        srnamt = "ZHPSVX";
        infot = 1;
        Chpsvx("/", "U", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZHPSVX", infot, nout, lerr, ok);
        infot = 2;
        Chpsvx("N", "/", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZHPSVX", infot, nout, lerr, ok);
        infot = 3;
        Chpsvx("N", "U", -1, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZHPSVX", infot, nout, lerr, ok);
        infot = 4;
        Chpsvx("N", "U", 0, -1, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZHPSVX", infot, nout, lerr, ok);
        infot = 9;
        Chpsvx("N", "U", 2, 0, a, af, ip, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZHPSVX", infot, nout, lerr, ok);
        infot = 11;
        Chpsvx("N", "U", 2, 0, a, af, ip, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZHPSVX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SY")) {
        //
        // Csysv
        //
        srnamt = "ZSYSV ";
        infot = 1;
        Csysv("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV ", infot, nout, lerr, ok);
        infot = 2;
        Csysv("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV ", infot, nout, lerr, ok);
        infot = 3;
        Csysv("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV ", infot, nout, lerr, ok);
        infot = 8;
        Csysv("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV ", infot, nout, lerr, ok);
        infot = 10;
        Csysv("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        Chkxer("ZSYSV ", infot, nout, lerr, ok);
        infot = 10;
        Csysv("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        Chkxer("ZSYSV ", infot, nout, lerr, ok);
        //
        // Csysvx
        //
        srnamt = "ZSYSVX";
        infot = 1;
        Csysvx("/", "U", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        Chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 2;
        Csysvx("N", "/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        Chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 3;
        Csysvx("N", "U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        Chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 4;
        Csysvx("N", "U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        Chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 6;
        Csysvx("N", "U", 2, 0, a, 1, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        Chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 8;
        Csysvx("N", "U", 2, 0, a, 2, af, 1, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        Chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 11;
        Csysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 1, x, 2, rcond, r1, r2, w, 4, rw, info);
        Chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 13;
        Csysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 1, rcond, r1, r2, w, 4, rw, info);
        Chkxer("ZSYSVX", infot, nout, lerr, ok);
        infot = 18;
        Csysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 3, rw, info);
        Chkxer("ZSYSVX", infot, nout, lerr, ok);
        //
        // Csysvxx
        //
        n_err_bnds = 3;
        nparams = 1;
        srnamt = "ZSYSVXX";
        infot = 1;
        eq = 'N';
        Csysvxx("/", "U", 0, 0, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 2;
        Csysvxx("N", "/", 0, 0, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 3;
        Csysvxx("N", "U", -1, 0, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 4;
        eq = '/';
        Csysvxx("N", "U", 0, -1, a, 1, af, 1, ip, eq, r, b, 1, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        eq = 'Y';
        infot = 6;
        Csysvxx("N", "U", 2, 0, a, 1, af, 2, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 8;
        Csysvxx("N", "U", 2, 0, a, 2, af, 1, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 10;
        Csysvxx("F", "U", 2, 0, a, 2, af, 2, ip, "A", r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 11;
        eq = 'Y';
        Csysvxx("F", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 11;
        eq = 'Y';
        r[1 - 1] = -one;
        Csysvxx("F", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 2, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 13;
        eq = 'N';
        Csysvxx("N", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 1, x, 2, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        infot = 15;
        Csysvxx("N", "U", 2, 0, a, 2, af, 2, ip, eq, r, b, 2, x, 1, rcond, rpvgrw, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, rw, info);
        Chkxer("ZSYSVXX", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SR")) {
        //
        // Csysv_rook
        //
        srnamt = "ZSYSV_ROOK";
        infot = 1;
        Csysv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Csysv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV_ROOK", infot, nout, lerr, ok);
        infot = 3;
        Csysv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV_ROOK", infot, nout, lerr, ok);
        infot = 8;
        Csysv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV_ROOK", infot, nout, lerr, ok);
        infot = 10;
        Csysv_rook("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        Chkxer("ZSYSV_ROOK", infot, nout, lerr, ok);
        infot = 10;
        Csysv_rook("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        //
    } else if (Mlsamen(2, c2.elems, "SK")) {
        //
        // Csysv_rk
        //
        // Test error exits of the driver that uses factorization
        // of a symmetric indefinite matrix with rook
        // (bounded Bunch-Kaufman) pivoting with the new storage
        // format for factors L ( or U) and D.
        //
        // L (or U) is stored in A, diagonal of D is stored on the
        // diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        srnamt = "ZSYSV_RK";
        infot = 1;
        Csysv_rk("/", 0, 0, a, 1, e, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 2;
        Csysv_rk("U", -1, 0, a, 1, e, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 3;
        Csysv_rk("U", 0, -1, a, 1, e, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 5;
        Csysv_rk("U", 2, 0, a, 1, e, ip, b, 2, w, 1, info);
        Chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 9;
        Csysv_rk("U", 2, 0, a, 2, e, ip, b, 1, w, 1, info);
        Chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 11;
        Csysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, 0, info);
        Chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        infot = 11;
        Csysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, -2, info);
        Chkxer("ZSYSV_RK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SP")) {
        //
        // Cspsv
        //
        srnamt = "ZSPSV ";
        infot = 1;
        Cspsv("/", 0, 0, a, ip, b, 1, info);
        Chkxer("ZSPSV ", infot, nout, lerr, ok);
        infot = 2;
        Cspsv("U", -1, 0, a, ip, b, 1, info);
        Chkxer("ZSPSV ", infot, nout, lerr, ok);
        infot = 3;
        Cspsv("U", 0, -1, a, ip, b, 1, info);
        Chkxer("ZSPSV ", infot, nout, lerr, ok);
        infot = 7;
        Cspsv("U", 2, 0, a, ip, b, 1, info);
        Chkxer("ZSPSV ", infot, nout, lerr, ok);
        //
        // Cspsvx
        //
        srnamt = "ZSPSVX";
        infot = 1;
        Cspsvx("/", "U", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZSPSVX", infot, nout, lerr, ok);
        infot = 2;
        Cspsvx("N", "/", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZSPSVX", infot, nout, lerr, ok);
        infot = 3;
        Cspsvx("N", "U", -1, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZSPSVX", infot, nout, lerr, ok);
        infot = 4;
        Cspsvx("N", "U", 0, -1, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZSPSVX", infot, nout, lerr, ok);
        infot = 9;
        Cspsvx("N", "U", 2, 0, a, af, ip, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        Chkxer("ZSPSVX", infot, nout, lerr, ok);
        infot = 11;
        Cspsvx("N", "U", 2, 0, a, af, ip, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        Chkxer("ZSPSVX", infot, nout, lerr, ok);
    }
    //
    // Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a3,' drivers passed the tests of the error exits')"), path;
    } else {
        write(nout, "(' *** ',a3,' drivers failed the tests of the error ','exits ***')"), path;
    }
    //
    // End of Cerrvx
    //
}
