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

// Derived from LAPACK routine ZERRHEX.
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

void Cerrhe(fem::str_cref path, INTEGER const nunit) {
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
    REAL s[nmax];
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
        s[j - 1] = 0.0;
        ip[j - 1] = j;
    }
    REAL anrm = 1.0;
    ok = true;
    //
    // Test error exits of the routines that use factorization
    // of a Hermitian indefinite matrix with patrial
    // (Bunch-Kaufman) diagonal pivoting method.
    //
    INTEGER info = 0;
    REAL r[nmax];
    INTEGER n_err_bnds = 0;
    INTEGER nparams = 0;
    fem::str<1> eq;
    REAL rcond = 0.0;
    REAL berr = 0.0;
    REAL err_bnds_n[nmax * 3];
    REAL err_bnds_c[nmax * 3];
    REAL params[1];
    if (Mlsamen(2, c2.elems, "HE")) {
        //
        // Chetrf
        //
        srnamt = "ZHETRF";
        infot = 1;
        Chetrf("/", 0, a, 1, ip, w, 1, info);
        Chkxer("ZHETRF", infot, nout, lerr, ok);
        infot = 2;
        Chetrf("U", -1, a, 1, ip, w, 1, info);
        Chkxer("ZHETRF", infot, nout, lerr, ok);
        infot = 4;
        Chetrf("U", 2, a, 1, ip, w, 4, info);
        Chkxer("ZHETRF", infot, nout, lerr, ok);
        infot = 7;
        Chetrf("U", 0, a, 1, ip, w, 0, info);
        Chkxer("ZHETRF", infot, nout, lerr, ok);
        infot = 7;
        Chetrf("U", 0, a, 1, ip, w, -2, info);
        Chkxer("ZHETRF", infot, nout, lerr, ok);
        //
        // Chetf2
        //
        srnamt = "ZHETF2";
        infot = 1;
        Chetf2("/", 0, a, 1, ip, info);
        Chkxer("ZHETF2", infot, nout, lerr, ok);
        infot = 2;
        Chetf2("U", -1, a, 1, ip, info);
        Chkxer("ZHETF2", infot, nout, lerr, ok);
        infot = 4;
        Chetf2("U", 2, a, 1, ip, info);
        Chkxer("ZHETF2", infot, nout, lerr, ok);
        //
        // Chetri
        //
        srnamt = "ZHETRI";
        infot = 1;
        Chetri("/", 0, a, 1, ip, w, info);
        Chkxer("ZHETRI", infot, nout, lerr, ok);
        infot = 2;
        Chetri("U", -1, a, 1, ip, w, info);
        Chkxer("ZHETRI", infot, nout, lerr, ok);
        infot = 4;
        Chetri("U", 2, a, 1, ip, w, info);
        Chkxer("ZHETRI", infot, nout, lerr, ok);
        //
        // Chetri2
        //
        srnamt = "ZHETRI2";
        infot = 1;
        Chetri2("/", 0, a, 1, ip, w, 1, info);
        Chkxer("ZHETRI2", infot, nout, lerr, ok);
        infot = 2;
        Chetri2("U", -1, a, 1, ip, w, 1, info);
        Chkxer("ZHETRI2", infot, nout, lerr, ok);
        infot = 4;
        Chetri2("U", 2, a, 1, ip, w, 1, info);
        Chkxer("ZHETRI2", infot, nout, lerr, ok);
        //
        // Chetri2x
        //
        srnamt = "ZHETRI2X";
        infot = 1;
        Chetri2x("/", 0, a, 1, ip, w, 1, info);
        Chkxer("ZHETRI2X", infot, nout, lerr, ok);
        infot = 2;
        Chetri2x("U", -1, a, 1, ip, w, 1, info);
        Chkxer("ZHETRI2X", infot, nout, lerr, ok);
        infot = 4;
        Chetri2x("U", 2, a, 1, ip, w, 1, info);
        Chkxer("ZHETRI2X", infot, nout, lerr, ok);
        //
        // Chetrs
        //
        srnamt = "ZHETRS";
        infot = 1;
        Chetrs("/", 0, 0, a, 1, ip, b, 1, info);
        Chkxer("ZHETRS", infot, nout, lerr, ok);
        infot = 2;
        Chetrs("U", -1, 0, a, 1, ip, b, 1, info);
        Chkxer("ZHETRS", infot, nout, lerr, ok);
        infot = 3;
        Chetrs("U", 0, -1, a, 1, ip, b, 1, info);
        Chkxer("ZHETRS", infot, nout, lerr, ok);
        infot = 5;
        Chetrs("U", 2, 1, a, 1, ip, b, 2, info);
        Chkxer("ZHETRS", infot, nout, lerr, ok);
        infot = 8;
        Chetrs("U", 2, 1, a, 2, ip, b, 1, info);
        Chkxer("ZHETRS", infot, nout, lerr, ok);
        //
        // Cherfs
        //
        srnamt = "ZHERFS";
        infot = 1;
        Cherfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("ZHERFS", infot, nout, lerr, ok);
        infot = 2;
        Cherfs("U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("ZHERFS", infot, nout, lerr, ok);
        infot = 3;
        Cherfs("U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("ZHERFS", infot, nout, lerr, ok);
        infot = 5;
        Cherfs("U", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, r, info);
        Chkxer("ZHERFS", infot, nout, lerr, ok);
        infot = 7;
        Cherfs("U", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, r, info);
        Chkxer("ZHERFS", infot, nout, lerr, ok);
        infot = 10;
        Cherfs("U", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, r, info);
        Chkxer("ZHERFS", infot, nout, lerr, ok);
        infot = 12;
        Cherfs("U", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, r, info);
        Chkxer("ZHERFS", infot, nout, lerr, ok);
        //
        // Cherfsx
        //
        n_err_bnds = 3;
        nparams = 0;
        srnamt = "ZHERFSX";
        infot = 1;
        Cherfsx("/", eq, 0, 0, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        Chkxer("ZHERFSX", infot, nout, lerr, ok);
        infot = 2;
        Cherfsx("U", eq, -1, 0, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        Chkxer("ZHERFSX", infot, nout, lerr, ok);
        eq = 'N';
        infot = 3;
        Cherfsx("U", eq, -1, 0, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        Chkxer("ZHERFSX", infot, nout, lerr, ok);
        infot = 4;
        Cherfsx("U", eq, 0, -1, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        Chkxer("ZHERFSX", infot, nout, lerr, ok);
        infot = 6;
        Cherfsx("U", eq, 2, 1, a, 1, af, 2, ip, s, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        Chkxer("ZHERFSX", infot, nout, lerr, ok);
        infot = 8;
        Cherfsx("U", eq, 2, 1, a, 2, af, 1, ip, s, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        Chkxer("ZHERFSX", infot, nout, lerr, ok);
        infot = 12;
        Cherfsx("U", eq, 2, 1, a, 2, af, 2, ip, s, b, 1, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        Chkxer("ZHERFSX", infot, nout, lerr, ok);
        infot = 14;
        Cherfsx("U", eq, 2, 1, a, 2, af, 2, ip, s, b, 2, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        Chkxer("ZHERFSX", infot, nout, lerr, ok);
        //
        // Checon
        //
        srnamt = "ZHECON";
        infot = 1;
        Checon("/", 0, a, 1, ip, anrm, rcond, w, info);
        Chkxer("ZHECON", infot, nout, lerr, ok);
        infot = 2;
        Checon("U", -1, a, 1, ip, anrm, rcond, w, info);
        Chkxer("ZHECON", infot, nout, lerr, ok);
        infot = 4;
        Checon("U", 2, a, 1, ip, anrm, rcond, w, info);
        Chkxer("ZHECON", infot, nout, lerr, ok);
        infot = 6;
        Checon("U", 1, a, 1, ip, -anrm, rcond, w, info);
        Chkxer("ZHECON", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "HR")) {
        //
        // Test error exits of the routines that use factorization
        // of a Hermitian indefinite matrix with rook
        // (bounded Bunch-Kaufman) diagonal pivoting method.
        //
        // Chetrf_rook
        //
        srnamt = "ZHETRF_ROOK";
        infot = 1;
        Chetrf_rook("/", 0, a, 1, ip, w, 1, info);
        Chkxer("ZHETRF_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Chetrf_rook("U", -1, a, 1, ip, w, 1, info);
        Chkxer("ZHETRF_ROOK", infot, nout, lerr, ok);
        infot = 4;
        Chetrf_rook("U", 2, a, 1, ip, w, 4, info);
        Chkxer("ZHETRF_ROOK", infot, nout, lerr, ok);
        infot = 7;
        Chetrf_rook("U", 0, a, 1, ip, w, 0, info);
        Chkxer("ZHETRF_ROOK", infot, nout, lerr, ok);
        infot = 7;
        Chetrf_rook("U", 0, a, 1, ip, w, -2, info);
        Chkxer("ZHETRF_ROOK", infot, nout, lerr, ok);
        //
        // Chetf2_rook
        //
        srnamt = "ZHETF2_ROOK";
        infot = 1;
        Chetf2_rook("/", 0, a, 1, ip, info);
        Chkxer("ZHETF2_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Chetf2_rook("U", -1, a, 1, ip, info);
        Chkxer("ZHETF2_ROOK", infot, nout, lerr, ok);
        infot = 4;
        Chetf2_rook("U", 2, a, 1, ip, info);
        Chkxer("ZHETF2_ROOK", infot, nout, lerr, ok);
        //
        // Chetri_rook
        //
        srnamt = "ZHETRI_ROOK";
        infot = 1;
        Chetri_rook("/", 0, a, 1, ip, w, info);
        Chkxer("ZHETRI_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Chetri_rook("U", -1, a, 1, ip, w, info);
        Chkxer("ZHETRI_ROOK", infot, nout, lerr, ok);
        infot = 4;
        Chetri_rook("U", 2, a, 1, ip, w, info);
        Chkxer("ZHETRI_ROOK", infot, nout, lerr, ok);
        //
        // Chetrs_rook
        //
        srnamt = "ZHETRS_ROOK";
        infot = 1;
        Chetrs_rook("/", 0, 0, a, 1, ip, b, 1, info);
        Chkxer("ZHETRS_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Chetrs_rook("U", -1, 0, a, 1, ip, b, 1, info);
        Chkxer("ZHETRS_ROOK", infot, nout, lerr, ok);
        infot = 3;
        Chetrs_rook("U", 0, -1, a, 1, ip, b, 1, info);
        Chkxer("ZHETRS_ROOK", infot, nout, lerr, ok);
        infot = 5;
        Chetrs_rook("U", 2, 1, a, 1, ip, b, 2, info);
        Chkxer("ZHETRS_ROOK", infot, nout, lerr, ok);
        infot = 8;
        Chetrs_rook("U", 2, 1, a, 2, ip, b, 1, info);
        Chkxer("ZHETRS_ROOK", infot, nout, lerr, ok);
        //
        // Checon_rook
        //
        srnamt = "ZHECON_ROOK";
        infot = 1;
        Checon_rook("/", 0, a, 1, ip, anrm, rcond, w, info);
        Chkxer("ZHECON_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Checon_rook("U", -1, a, 1, ip, anrm, rcond, w, info);
        Chkxer("ZHECON_ROOK", infot, nout, lerr, ok);
        infot = 4;
        Checon_rook("U", 2, a, 1, ip, anrm, rcond, w, info);
        Chkxer("ZHECON_ROOK", infot, nout, lerr, ok);
        infot = 6;
        Checon_rook("U", 1, a, 1, ip, -anrm, rcond, w, info);
        Chkxer("ZHECON_ROOK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "HK")) {
        //
        // Test error exits of the routines that use factorization
        // of a symmetric indefinite matrix with rook
        // (bounded Bunch-Kaufman) pivoting with the new storage
        // format for factors L ( or U) and D.
        //
        // L (or U) is stored in A, diagonal of D is stored on the
        // diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        // Chetrf_rk
        //
        srnamt = "ZHETRF_RK";
        infot = 1;
        Chetrf_rk("/", 0, a, 1, e, ip, w, 1, info);
        Chkxer("ZHETRF_RK", infot, nout, lerr, ok);
        infot = 2;
        Chetrf_rk("U", -1, a, 1, e, ip, w, 1, info);
        Chkxer("ZHETRF_RK", infot, nout, lerr, ok);
        infot = 4;
        Chetrf_rk("U", 2, a, 1, e, ip, w, 4, info);
        Chkxer("ZHETRF_RK", infot, nout, lerr, ok);
        infot = 8;
        Chetrf_rk("U", 0, a, 1, e, ip, w, 0, info);
        Chkxer("ZHETRF_RK", infot, nout, lerr, ok);
        infot = 8;
        Chetrf_rk("U", 0, a, 1, e, ip, w, -2, info);
        Chkxer("ZHETRF_RK", infot, nout, lerr, ok);
        //
        // Chetf2_rk
        //
        srnamt = "ZHETF2_RK";
        infot = 1;
        Chetf2_rk("/", 0, a, 1, e, ip, info);
        Chkxer("ZHETF2_RK", infot, nout, lerr, ok);
        infot = 2;
        Chetf2_rk("U", -1, a, 1, e, ip, info);
        Chkxer("ZHETF2_RK", infot, nout, lerr, ok);
        infot = 4;
        Chetf2_rk("U", 2, a, 1, e, ip, info);
        Chkxer("ZHETF2_RK", infot, nout, lerr, ok);
        //
        // Chetri_3
        //
        srnamt = "ZHETRI_3";
        infot = 1;
        Chetri_3("/", 0, a, 1, e, ip, w, 1, info);
        Chkxer("ZHETRI_3", infot, nout, lerr, ok);
        infot = 2;
        Chetri_3("U", -1, a, 1, e, ip, w, 1, info);
        Chkxer("ZHETRI_3", infot, nout, lerr, ok);
        infot = 4;
        Chetri_3("U", 2, a, 1, e, ip, w, 1, info);
        Chkxer("ZHETRI_3", infot, nout, lerr, ok);
        infot = 8;
        Chetri_3("U", 0, a, 1, e, ip, w, 0, info);
        Chkxer("ZHETRI_3", infot, nout, lerr, ok);
        infot = 8;
        Chetri_3("U", 0, a, 1, e, ip, w, -2, info);
        Chkxer("ZHETRI_3", infot, nout, lerr, ok);
        //
        // Chetri_3x
        //
        srnamt = "ZHETRI_3X";
        infot = 1;
        Chetri_3x("/", 0, a, 1, e, ip, w, 1, info);
        Chkxer("ZHETRI_3X", infot, nout, lerr, ok);
        infot = 2;
        Chetri_3x("U", -1, a, 1, e, ip, w, 1, info);
        Chkxer("ZHETRI_3X", infot, nout, lerr, ok);
        infot = 4;
        Chetri_3x("U", 2, a, 1, e, ip, w, 1, info);
        Chkxer("ZHETRI_3X", infot, nout, lerr, ok);
        //
        // Chetrs_3
        //
        srnamt = "ZHETRS_3";
        infot = 1;
        Chetrs_3("/", 0, 0, a, 1, e, ip, b, 1, info);
        Chkxer("ZHETRS_3", infot, nout, lerr, ok);
        infot = 2;
        Chetrs_3("U", -1, 0, a, 1, e, ip, b, 1, info);
        Chkxer("ZHETRS_3", infot, nout, lerr, ok);
        infot = 3;
        Chetrs_3("U", 0, -1, a, 1, e, ip, b, 1, info);
        Chkxer("ZHETRS_3", infot, nout, lerr, ok);
        infot = 5;
        Chetrs_3("U", 2, 1, a, 1, e, ip, b, 2, info);
        Chkxer("ZHETRS_3", infot, nout, lerr, ok);
        infot = 9;
        Chetrs_3("U", 2, 1, a, 2, e, ip, b, 1, info);
        Chkxer("ZHETRS_3", infot, nout, lerr, ok);
        //
        // Checon_3
        //
        srnamt = "ZHECON_3";
        infot = 1;
        Checon_3("/", 0, a, 1, e, ip, anrm, rcond, w, info);
        Chkxer("ZHECON_3", infot, nout, lerr, ok);
        infot = 2;
        Checon_3("U", -1, a, 1, e, ip, anrm, rcond, w, info);
        Chkxer("ZHECON_3", infot, nout, lerr, ok);
        infot = 4;
        Checon_3("U", 2, a, 1, e, ip, anrm, rcond, w, info);
        Chkxer("ZHECON_3", infot, nout, lerr, ok);
        infot = 7;
        Checon_3("U", 1, a, 1, e, ip, -1.0, rcond, w, info);
        Chkxer("ZHECON_3", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "HP")) {
        //
        // Test error exits of the routines that use factorization
        // of a Hermitian indefinite packed matrix with patrial
        // (Bunch-Kaufman) diagonal pivoting method.
        //
        // Chptrf
        //
        srnamt = "ZHPTRF";
        infot = 1;
        Chptrf("/", 0, a, ip, info);
        Chkxer("ZHPTRF", infot, nout, lerr, ok);
        infot = 2;
        Chptrf("U", -1, a, ip, info);
        Chkxer("ZHPTRF", infot, nout, lerr, ok);
        //
        // Chptri
        //
        srnamt = "ZHPTRI";
        infot = 1;
        Chptri("/", 0, a, ip, w, info);
        Chkxer("ZHPTRI", infot, nout, lerr, ok);
        infot = 2;
        Chptri("U", -1, a, ip, w, info);
        Chkxer("ZHPTRI", infot, nout, lerr, ok);
        //
        // Chptrs
        //
        srnamt = "ZHPTRS";
        infot = 1;
        Chptrs("/", 0, 0, a, ip, b, 1, info);
        Chkxer("ZHPTRS", infot, nout, lerr, ok);
        infot = 2;
        Chptrs("U", -1, 0, a, ip, b, 1, info);
        Chkxer("ZHPTRS", infot, nout, lerr, ok);
        infot = 3;
        Chptrs("U", 0, -1, a, ip, b, 1, info);
        Chkxer("ZHPTRS", infot, nout, lerr, ok);
        infot = 7;
        Chptrs("U", 2, 1, a, ip, b, 1, info);
        Chkxer("ZHPTRS", infot, nout, lerr, ok);
        //
        // Chprfs
        //
        srnamt = "ZHPRFS";
        infot = 1;
        Chprfs("/", 0, 0, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("ZHPRFS", infot, nout, lerr, ok);
        infot = 2;
        Chprfs("U", -1, 0, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("ZHPRFS", infot, nout, lerr, ok);
        infot = 3;
        Chprfs("U", 0, -1, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        Chkxer("ZHPRFS", infot, nout, lerr, ok);
        infot = 8;
        Chprfs("U", 2, 1, a, af, ip, b, 1, x, 2, r1, r2, w, r, info);
        Chkxer("ZHPRFS", infot, nout, lerr, ok);
        infot = 10;
        Chprfs("U", 2, 1, a, af, ip, b, 2, x, 1, r1, r2, w, r, info);
        Chkxer("ZHPRFS", infot, nout, lerr, ok);
        //
        // Chpcon
        //
        srnamt = "ZHPCON";
        infot = 1;
        Chpcon("/", 0, a, ip, anrm, rcond, w, info);
        Chkxer("ZHPCON", infot, nout, lerr, ok);
        infot = 2;
        Chpcon("U", -1, a, ip, anrm, rcond, w, info);
        Chkxer("ZHPCON", infot, nout, lerr, ok);
        infot = 5;
        Chpcon("U", 1, a, ip, -anrm, rcond, w, info);
        Chkxer("ZHPCON", infot, nout, lerr, ok);
    }
    //
    // Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    // End of Cerrhe
    //
}
