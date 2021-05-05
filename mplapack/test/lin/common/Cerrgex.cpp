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

void Cerrge(const char *path, INTEGER const nunit) {
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
    COMPLEX a[nmax * nmax];
    COMPLEX af[nmax * nmax];
    COMPLEX b[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    COMPLEX w[2 * nmax];
    COMPLEX x[nmax];
    REAL cs[nmax];
    REAL rs[nmax];
    INTEGER ip[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = COMPLEX(1.0 / (i + j).real(), -1.0 / (i + j).real());
            af[(i - 1) + (j - 1) * ldaf] = COMPLEX(1.0 / (i + j).real(), -1.0 / (i + j).real());
        }
        b[j - 1] = 0.0;
        r1[j - 1] = 0.0;
        r2[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
        cs[j - 1] = 0.0;
        rs[j - 1] = 0.0;
        ip[j - 1] = j;
    }
    ok = true;
    //
    //     Test error exits of the routines that use the LU decomposition
    //     of a general matrix.
    //
    INTEGER info = 0;
    REAL r[nmax];
    INTEGER n_err_bnds = 0;
    INTEGER nparams = 0;
    char eq;
    REAL rcond = 0.0;
    REAL berr = 0.0;
    COMPLEX err_bnds_n[nmax * 3];
    COMPLEX err_bnds_c[nmax * 3];
    COMPLEX params = 0.0;
    REAL anrm = 0.0;
    REAL ccond = 0.0;
    if (Mlsamen(2, c2, "GE")) {
        //
        //        Cgetrf
        //
        infot = 1;
        Cgetrf(-1, 0, a, 1, ip, info);
        chkxer("Cgetrf", infot, nout, lerr, ok);
        infot = 2;
        Cgetrf(0, -1, a, 1, ip, info);
        chkxer("Cgetrf", infot, nout, lerr, ok);
        infot = 4;
        Cgetrf(2, 1, a, 1, ip, info);
        chkxer("Cgetrf", infot, nout, lerr, ok);
        //
        //        Cgetf2
        //
        infot = 1;
        Cgetf2(-1, 0, a, 1, ip, info);
        chkxer("Cgetf2", infot, nout, lerr, ok);
        infot = 2;
        Cgetf2(0, -1, a, 1, ip, info);
        chkxer("Cgetf2", infot, nout, lerr, ok);
        infot = 4;
        Cgetf2(2, 1, a, 1, ip, info);
        chkxer("Cgetf2", infot, nout, lerr, ok);
        //
        //        Cgetri
        //
        infot = 1;
        Cgetri(-1, a, 1, ip, w, 1, info);
        chkxer("Cgetri", infot, nout, lerr, ok);
        infot = 3;
        Cgetri(2, a, 1, ip, w, 2, info);
        chkxer("Cgetri", infot, nout, lerr, ok);
        infot = 6;
        Cgetri(2, a, 2, ip, w, 1, info);
        chkxer("Cgetri", infot, nout, lerr, ok);
        //
        //        Cgetrs
        //
        infot = 1;
        Cgetrs("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("Cgetrs", infot, nout, lerr, ok);
        infot = 2;
        Cgetrs("N", -1, 0, a, 1, ip, b, 1, info);
        chkxer("Cgetrs", infot, nout, lerr, ok);
        infot = 3;
        Cgetrs("N", 0, -1, a, 1, ip, b, 1, info);
        chkxer("Cgetrs", infot, nout, lerr, ok);
        infot = 5;
        Cgetrs("N", 2, 1, a, 1, ip, b, 2, info);
        chkxer("Cgetrs", infot, nout, lerr, ok);
        infot = 8;
        Cgetrs("N", 2, 1, a, 2, ip, b, 1, info);
        chkxer("Cgetrs", infot, nout, lerr, ok);
        //
        //        Cgerfs
        //
        infot = 1;
        Cgerfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cgerfs", infot, nout, lerr, ok);
        infot = 2;
        Cgerfs("N", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cgerfs", infot, nout, lerr, ok);
        infot = 3;
        Cgerfs("N", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cgerfs", infot, nout, lerr, ok);
        infot = 5;
        Cgerfs("N", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("Cgerfs", infot, nout, lerr, ok);
        infot = 7;
        Cgerfs("N", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("Cgerfs", infot, nout, lerr, ok);
        infot = 10;
        Cgerfs("N", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("Cgerfs", infot, nout, lerr, ok);
        infot = 12;
        Cgerfs("N", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("Cgerfs", infot, nout, lerr, ok);
        //
        //        Cgerfsx
        //
        n_err_bnds = 3;
        nparams = 0;
        infot = 1;
        Cgerfsx("/", eq, 0, 0, a, 1, af, 1, ip, rs, cs, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgerfsx", infot, nout, lerr, ok);
        infot = 2;
        eq = "/";
        Cgerfsx("N", eq, 2, 1, a, 1, af, 2, ip, rs, cs, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgerfsx", infot, nout, lerr, ok);
        infot = 3;
        eq = "R";
        Cgerfsx("N", eq, -1, 0, a, 1, af, 1, ip, rs, cs, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgerfsx", infot, nout, lerr, ok);
        infot = 4;
        Cgerfsx("N", eq, 0, -1, a, 1, af, 1, ip, rs, cs, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgerfsx", infot, nout, lerr, ok);
        infot = 6;
        Cgerfsx("N", eq, 2, 1, a, 1, af, 2, ip, rs, cs, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgerfsx", infot, nout, lerr, ok);
        infot = 8;
        Cgerfsx("N", eq, 2, 1, a, 2, af, 1, ip, rs, cs, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgerfsx", infot, nout, lerr, ok);
        infot = 13;
        eq = "C";
        Cgerfsx("N", eq, 2, 1, a, 2, af, 2, ip, rs, cs, b, 1, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgerfsx", infot, nout, lerr, ok);
        infot = 15;
        Cgerfsx("N", eq, 2, 1, a, 2, af, 2, ip, rs, cs, b, 2, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgerfsx", infot, nout, lerr, ok);
        //
        //        Cgecon
        //
        infot = 1;
        Cgecon("/", 0, a, 1, anrm, rcond, w, r, info);
        chkxer("Cgecon", infot, nout, lerr, ok);
        infot = 2;
        Cgecon("1", -1, a, 1, anrm, rcond, w, r, info);
        chkxer("Cgecon", infot, nout, lerr, ok);
        infot = 4;
        Cgecon("1", 2, a, 1, anrm, rcond, w, r, info);
        chkxer("Cgecon", infot, nout, lerr, ok);
        //
        //        Cgeequ
        //
        infot = 1;
        Cgeequ(-1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgeequ", infot, nout, lerr, ok);
        infot = 2;
        Cgeequ(0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgeequ", infot, nout, lerr, ok);
        infot = 4;
        Cgeequ(2, 2, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgeequ", infot, nout, lerr, ok);
        //
        //        Cgeequb
        //
        infot = 1;
        Cgeequb(-1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgeequb", infot, nout, lerr, ok);
        infot = 2;
        Cgeequb(0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgeequb", infot, nout, lerr, ok);
        infot = 4;
        Cgeequb(2, 2, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgeequb", infot, nout, lerr, ok);
        //
        //     Test error exits of the routines that use the LU decomposition
        //     of a general band matrix.
        //
    } else if (Mlsamen(2, c2, "GB")) {
        //
        //        Cgbtrf
        //
        infot = 1;
        Cgbtrf(-1, 0, 0, 0, a, 1, ip, info);
        chkxer("Cgbtrf", infot, nout, lerr, ok);
        infot = 2;
        Cgbtrf(0, -1, 0, 0, a, 1, ip, info);
        chkxer("Cgbtrf", infot, nout, lerr, ok);
        infot = 3;
        Cgbtrf(1, 1, -1, 0, a, 1, ip, info);
        chkxer("Cgbtrf", infot, nout, lerr, ok);
        infot = 4;
        Cgbtrf(1, 1, 0, -1, a, 1, ip, info);
        chkxer("Cgbtrf", infot, nout, lerr, ok);
        infot = 6;
        Cgbtrf(2, 2, 1, 1, a, 3, ip, info);
        chkxer("Cgbtrf", infot, nout, lerr, ok);
        //
        //        Cgbtf2
        //
        infot = 1;
        Cgbtf2(-1, 0, 0, 0, a, 1, ip, info);
        chkxer("Cgbtf2", infot, nout, lerr, ok);
        infot = 2;
        Cgbtf2(0, -1, 0, 0, a, 1, ip, info);
        chkxer("Cgbtf2", infot, nout, lerr, ok);
        infot = 3;
        Cgbtf2(1, 1, -1, 0, a, 1, ip, info);
        chkxer("Cgbtf2", infot, nout, lerr, ok);
        infot = 4;
        Cgbtf2(1, 1, 0, -1, a, 1, ip, info);
        chkxer("Cgbtf2", infot, nout, lerr, ok);
        infot = 6;
        Cgbtf2(2, 2, 1, 1, a, 3, ip, info);
        chkxer("Cgbtf2", infot, nout, lerr, ok);
        //
        //        Cgbtrs
        //
        infot = 1;
        Cgbtrs("/", 0, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("Cgbtrs", infot, nout, lerr, ok);
        infot = 2;
        Cgbtrs("N", -1, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("Cgbtrs", infot, nout, lerr, ok);
        infot = 3;
        Cgbtrs("N", 1, -1, 0, 1, a, 1, ip, b, 1, info);
        chkxer("Cgbtrs", infot, nout, lerr, ok);
        infot = 4;
        Cgbtrs("N", 1, 0, -1, 1, a, 1, ip, b, 1, info);
        chkxer("Cgbtrs", infot, nout, lerr, ok);
        infot = 5;
        Cgbtrs("N", 1, 0, 0, -1, a, 1, ip, b, 1, info);
        chkxer("Cgbtrs", infot, nout, lerr, ok);
        infot = 7;
        Cgbtrs("N", 2, 1, 1, 1, a, 3, ip, b, 2, info);
        chkxer("Cgbtrs", infot, nout, lerr, ok);
        infot = 10;
        Cgbtrs("N", 2, 0, 0, 1, a, 1, ip, b, 1, info);
        chkxer("Cgbtrs", infot, nout, lerr, ok);
        //
        //        Cgbrfs
        //
        infot = 1;
        Cgbrfs("/", 0, 0, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cgbrfs", infot, nout, lerr, ok);
        infot = 2;
        Cgbrfs("N", -1, 0, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cgbrfs", infot, nout, lerr, ok);
        infot = 3;
        Cgbrfs("N", 1, -1, 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cgbrfs", infot, nout, lerr, ok);
        infot = 4;
        Cgbrfs("N", 1, 0, -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cgbrfs", infot, nout, lerr, ok);
        infot = 5;
        Cgbrfs("N", 1, 0, 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Cgbrfs", infot, nout, lerr, ok);
        infot = 7;
        Cgbrfs("N", 2, 1, 1, 1, a, 2, af, 4, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("Cgbrfs", infot, nout, lerr, ok);
        infot = 9;
        Cgbrfs("N", 2, 1, 1, 1, a, 3, af, 3, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("Cgbrfs", infot, nout, lerr, ok);
        infot = 12;
        Cgbrfs("N", 2, 0, 0, 1, a, 1, af, 1, ip, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("Cgbrfs", infot, nout, lerr, ok);
        infot = 14;
        Cgbrfs("N", 2, 0, 0, 1, a, 1, af, 1, ip, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("Cgbrfs", infot, nout, lerr, ok);
        //
        //        Cgbrfsx
        //
        n_err_bnds = 3;
        nparams = 0;
        infot = 1;
        Cgbrfsx("/", eq, 0, 0, 0, 0, a, 1, af, 1, ip, rs, cs, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgbrfsx", infot, nout, lerr, ok);
        infot = 2;
        eq = "/";
        Cgbrfsx("N", eq, 2, 1, 1, 1, a, 1, af, 2, ip, rs, cs, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgbrfsx", infot, nout, lerr, ok);
        infot = 3;
        eq = "R";
        Cgbrfsx("N", eq, -1, 1, 1, 0, a, 1, af, 1, ip, rs, cs, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgbrfsx", infot, nout, lerr, ok);
        infot = 4;
        eq = "R";
        Cgbrfsx("N", eq, 2, -1, 1, 1, a, 3, af, 4, ip, rs, cs, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgbrfsx", infot, nout, lerr, ok);
        infot = 5;
        eq = "R";
        Cgbrfsx("N", eq, 2, 1, -1, 1, a, 3, af, 4, ip, rs, cs, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgbrfsx", infot, nout, lerr, ok);
        infot = 6;
        Cgbrfsx("N", eq, 0, 0, 0, -1, a, 1, af, 1, ip, rs, cs, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgbrfsx", infot, nout, lerr, ok);
        infot = 8;
        Cgbrfsx("N", eq, 2, 1, 1, 1, a, 1, af, 2, ip, rs, cs, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgbrfsx", infot, nout, lerr, ok);
        infot = 10;
        Cgbrfsx("N", eq, 2, 1, 1, 1, a, 3, af, 3, ip, rs, cs, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgbrfsx", infot, nout, lerr, ok);
        infot = 13;
        eq = "C";
        Cgbrfsx("N", eq, 2, 1, 1, 1, a, 3, af, 5, ip, rs, cs, b, 1, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgbrfsx", infot, nout, lerr, ok);
        infot = 15;
        Cgbrfsx("N", eq, 2, 1, 1, 1, a, 3, af, 5, ip, rs, cs, b, 2, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("Cgbrfsx", infot, nout, lerr, ok);
        //
        //        Cgbcon
        //
        infot = 1;
        Cgbcon("/", 0, 0, 0, a, 1, ip, anrm, rcond, w, r, info);
        chkxer("Cgbcon", infot, nout, lerr, ok);
        infot = 2;
        Cgbcon("1", -1, 0, 0, a, 1, ip, anrm, rcond, w, r, info);
        chkxer("Cgbcon", infot, nout, lerr, ok);
        infot = 3;
        Cgbcon("1", 1, -1, 0, a, 1, ip, anrm, rcond, w, r, info);
        chkxer("Cgbcon", infot, nout, lerr, ok);
        infot = 4;
        Cgbcon("1", 1, 0, -1, a, 1, ip, anrm, rcond, w, r, info);
        chkxer("Cgbcon", infot, nout, lerr, ok);
        infot = 6;
        Cgbcon("1", 2, 1, 1, a, 3, ip, anrm, rcond, w, r, info);
        chkxer("Cgbcon", infot, nout, lerr, ok);
        //
        //        Cgbequ
        //
        infot = 1;
        Cgbequ(-1, 0, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgbequ", infot, nout, lerr, ok);
        infot = 2;
        Cgbequ(0, -1, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgbequ", infot, nout, lerr, ok);
        infot = 3;
        Cgbequ(1, 1, -1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgbequ", infot, nout, lerr, ok);
        infot = 4;
        Cgbequ(1, 1, 0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgbequ", infot, nout, lerr, ok);
        infot = 6;
        Cgbequ(2, 2, 1, 1, a, 2, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgbequ", infot, nout, lerr, ok);
        //
        //        Cgbequb
        //
        infot = 1;
        Cgbequb(-1, 0, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgbequb", infot, nout, lerr, ok);
        infot = 2;
        Cgbequb(0, -1, 0, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgbequb", infot, nout, lerr, ok);
        infot = 3;
        Cgbequb(1, 1, -1, 0, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgbequb", infot, nout, lerr, ok);
        infot = 4;
        Cgbequb(1, 1, 0, -1, a, 1, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgbequb", infot, nout, lerr, ok);
        infot = 6;
        Cgbequb(2, 2, 1, 1, a, 2, r1, r2, rcond, ccond, anrm, info);
        chkxer("Cgbequb", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrge
    //
}
