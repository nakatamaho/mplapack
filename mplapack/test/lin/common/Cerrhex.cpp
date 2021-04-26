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

void Cerrhe(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_1d<nmax, REAL> s(fill0);
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
        s[j - 1] = 0.0;
        ip[j - 1] = j;
    }
    REAL anrm = 1.0;
    ok = true;
    //
    //     Test error exits of the routines that use factorization
    //     of a Hermitian indefinite matrix with patrial
    //     (Bunch-Kaufman) diagonal pivoting method.
    //
    INTEGER info = 0;
    arr_1d<nmax, REAL> r(fill0);
    INTEGER n_err_bnds = 0;
    INTEGER nparams = 0;
    char eq = char0;
    REAL rcond = 0.0;
    REAL berr = 0.0;
    arr_2d<nmax, 3, REAL> err_bnds_n(fill0);
    arr_2d<nmax, 3, REAL> err_bnds_c(fill0);
    arr_1d<1, REAL> params(fill0);
    if (Mlsamen2, c2, "HE") {
        //
        //        ZHETRF
        //
        srnamt = "ZHETRF";
        infot = 1;
        zhetrf("/", 0, a, 1, ip, w, 1, info);
        chkxer("ZHETRF", infot, nout, lerr, ok);
        infot = 2;
        zhetrf("U", -1, a, 1, ip, w, 1, info);
        chkxer("ZHETRF", infot, nout, lerr, ok);
        infot = 4;
        zhetrf("U", 2, a, 1, ip, w, 4, info);
        chkxer("ZHETRF", infot, nout, lerr, ok);
        infot = 7;
        zhetrf("U", 0, a, 1, ip, w, 0, info);
        chkxer("ZHETRF", infot, nout, lerr, ok);
        infot = 7;
        zhetrf("U", 0, a, 1, ip, w, -2, info);
        chkxer("ZHETRF", infot, nout, lerr, ok);
        //
        //        ZHETF2
        //
        srnamt = "ZHETF2";
        infot = 1;
        zhetf2("/", 0, a, 1, ip, info);
        chkxer("ZHETF2", infot, nout, lerr, ok);
        infot = 2;
        zhetf2("U", -1, a, 1, ip, info);
        chkxer("ZHETF2", infot, nout, lerr, ok);
        infot = 4;
        zhetf2("U", 2, a, 1, ip, info);
        chkxer("ZHETF2", infot, nout, lerr, ok);
        //
        //        ZHETRI
        //
        srnamt = "ZHETRI";
        infot = 1;
        zhetri("/", 0, a, 1, ip, w, info);
        chkxer("ZHETRI", infot, nout, lerr, ok);
        infot = 2;
        zhetri("U", -1, a, 1, ip, w, info);
        chkxer("ZHETRI", infot, nout, lerr, ok);
        infot = 4;
        zhetri("U", 2, a, 1, ip, w, info);
        chkxer("ZHETRI", infot, nout, lerr, ok);
        //
        //        ZHETRI2
        //
        srnamt = "ZHETRI2";
        infot = 1;
        zhetri2("/", 0, a, 1, ip, w, 1, info);
        chkxer("ZHETRI2", infot, nout, lerr, ok);
        infot = 2;
        zhetri2("U", -1, a, 1, ip, w, 1, info);
        chkxer("ZHETRI2", infot, nout, lerr, ok);
        infot = 4;
        zhetri2("U", 2, a, 1, ip, w, 1, info);
        chkxer("ZHETRI2", infot, nout, lerr, ok);
        //
        //        ZHETRI2X
        //
        srnamt = "ZHETRI2X";
        infot = 1;
        zhetri2x("/", 0, a, 1, ip, w, 1, info);
        chkxer("ZHETRI2X", infot, nout, lerr, ok);
        infot = 2;
        zhetri2x("U", -1, a, 1, ip, w, 1, info);
        chkxer("ZHETRI2X", infot, nout, lerr, ok);
        infot = 4;
        zhetri2x("U", 2, a, 1, ip, w, 1, info);
        chkxer("ZHETRI2X", infot, nout, lerr, ok);
        //
        //        ZHETRS
        //
        srnamt = "ZHETRS";
        infot = 1;
        zhetrs("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("ZHETRS", infot, nout, lerr, ok);
        infot = 2;
        zhetrs("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("ZHETRS", infot, nout, lerr, ok);
        infot = 3;
        zhetrs("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("ZHETRS", infot, nout, lerr, ok);
        infot = 5;
        zhetrs("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("ZHETRS", infot, nout, lerr, ok);
        infot = 8;
        zhetrs("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("ZHETRS", infot, nout, lerr, ok);
        //
        //        CherFS
        //
        srnamt = "CherFS";
        infot = 1;
        Cherfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("CherFS", infot, nout, lerr, ok);
        infot = 2;
        Cherfs("U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("CherFS", infot, nout, lerr, ok);
        infot = 3;
        Cherfs("U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("CherFS", infot, nout, lerr, ok);
        infot = 5;
        Cherfs("U", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("CherFS", infot, nout, lerr, ok);
        infot = 7;
        Cherfs("U", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("CherFS", infot, nout, lerr, ok);
        infot = 10;
        Cherfs("U", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("CherFS", infot, nout, lerr, ok);
        infot = 12;
        Cherfs("U", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("CherFS", infot, nout, lerr, ok);
        //
        //        CherFSX
        //
        n_err_bnds = 3;
        nparams = 0;
        srnamt = "CherFSX";
        infot = 1;
        Cherfsx("/", eq, 0, 0, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("CherFSX", infot, nout, lerr, ok);
        infot = 2;
        Cherfsx("U", eq, -1, 0, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("CherFSX", infot, nout, lerr, ok);
        eq = "N";
        infot = 3;
        Cherfsx("U", eq, -1, 0, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("CherFSX", infot, nout, lerr, ok);
        infot = 4;
        Cherfsx("U", eq, 0, -1, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("CherFSX", infot, nout, lerr, ok);
        infot = 6;
        Cherfsx("U", eq, 2, 1, a, 1, af, 2, ip, s, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("CherFSX", infot, nout, lerr, ok);
        infot = 8;
        Cherfsx("U", eq, 2, 1, a, 2, af, 1, ip, s, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("CherFSX", infot, nout, lerr, ok);
        infot = 12;
        Cherfsx("U", eq, 2, 1, a, 2, af, 2, ip, s, b, 1, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("CherFSX", infot, nout, lerr, ok);
        infot = 14;
        Cherfsx("U", eq, 2, 1, a, 2, af, 2, ip, s, b, 2, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, r, info);
        chkxer("CherFSX", infot, nout, lerr, ok);
        //
        //        ZHECON
        //
        srnamt = "ZHECON";
        infot = 1;
        zhecon("/", 0, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZHECON", infot, nout, lerr, ok);
        infot = 2;
        zhecon("U", -1, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZHECON", infot, nout, lerr, ok);
        infot = 4;
        zhecon("U", 2, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZHECON", infot, nout, lerr, ok);
        infot = 6;
        zhecon("U", 1, a, 1, ip, -anrm, rcond, w, info);
        chkxer("ZHECON", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "HR") {
        //
        //        Test error exits of the routines that use factorization
        //        of a Hermitian indefinite matrix with rook
        //        (bounded Bunch-Kaufman) diagonal pivoting method.
        //
        //        ZHETRF_ROOK
        //
        srnamt = "ZHETRF_ROOK";
        infot = 1;
        zhetrf_rook("/", 0, a, 1, ip, w, 1, info);
        chkxer("ZHETRF_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zhetrf_rook("U", -1, a, 1, ip, w, 1, info);
        chkxer("ZHETRF_ROOK", infot, nout, lerr, ok);
        infot = 4;
        zhetrf_rook("U", 2, a, 1, ip, w, 4, info);
        chkxer("ZHETRF_ROOK", infot, nout, lerr, ok);
        infot = 7;
        zhetrf_rook("U", 0, a, 1, ip, w, 0, info);
        chkxer("ZHETRF_ROOK", infot, nout, lerr, ok);
        infot = 7;
        zhetrf_rook("U", 0, a, 1, ip, w, -2, info);
        chkxer("ZHETRF_ROOK", infot, nout, lerr, ok);
        //
        //        ZHETF2_ROOK
        //
        srnamt = "ZHETF2_ROOK";
        infot = 1;
        zhetf2_rook("/", 0, a, 1, ip, info);
        chkxer("ZHETF2_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zhetf2_rook("U", -1, a, 1, ip, info);
        chkxer("ZHETF2_ROOK", infot, nout, lerr, ok);
        infot = 4;
        zhetf2_rook("U", 2, a, 1, ip, info);
        chkxer("ZHETF2_ROOK", infot, nout, lerr, ok);
        //
        //        ZHETRI_ROOK
        //
        srnamt = "ZHETRI_ROOK";
        infot = 1;
        zhetri_rook("/", 0, a, 1, ip, w, info);
        chkxer("ZHETRI_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zhetri_rook("U", -1, a, 1, ip, w, info);
        chkxer("ZHETRI_ROOK", infot, nout, lerr, ok);
        infot = 4;
        zhetri_rook("U", 2, a, 1, ip, w, info);
        chkxer("ZHETRI_ROOK", infot, nout, lerr, ok);
        //
        //        ZHETRS_ROOK
        //
        srnamt = "ZHETRS_ROOK";
        infot = 1;
        zhetrs_rook("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("ZHETRS_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zhetrs_rook("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("ZHETRS_ROOK", infot, nout, lerr, ok);
        infot = 3;
        zhetrs_rook("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("ZHETRS_ROOK", infot, nout, lerr, ok);
        infot = 5;
        zhetrs_rook("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("ZHETRS_ROOK", infot, nout, lerr, ok);
        infot = 8;
        zhetrs_rook("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("ZHETRS_ROOK", infot, nout, lerr, ok);
        //
        //        ZHECON_ROOK
        //
        srnamt = "ZHECON_ROOK";
        infot = 1;
        zhecon_rook("/", 0, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZHECON_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zhecon_rook("U", -1, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZHECON_ROOK", infot, nout, lerr, ok);
        infot = 4;
        zhecon_rook("U", 2, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZHECON_ROOK", infot, nout, lerr, ok);
        infot = 6;
        zhecon_rook("U", 1, a, 1, ip, -anrm, rcond, w, info);
        chkxer("ZHECON_ROOK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "HK") {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting with the new storage
        //        format for factors L ( or U) and D.
        //
        //        L (or U) is stored in A, diagonal of D is stored on the
        //        diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        //        ZHETRF_RK
        //
        srnamt = "ZHETRF_RK";
        infot = 1;
        zhetrf_rk("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("ZHETRF_RK", infot, nout, lerr, ok);
        infot = 2;
        zhetrf_rk("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("ZHETRF_RK", infot, nout, lerr, ok);
        infot = 4;
        zhetrf_rk("U", 2, a, 1, e, ip, w, 4, info);
        chkxer("ZHETRF_RK", infot, nout, lerr, ok);
        infot = 8;
        zhetrf_rk("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("ZHETRF_RK", infot, nout, lerr, ok);
        infot = 8;
        zhetrf_rk("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("ZHETRF_RK", infot, nout, lerr, ok);
        //
        //        ZHETF2_RK
        //
        srnamt = "ZHETF2_RK";
        infot = 1;
        zhetf2_rk("/", 0, a, 1, e, ip, info);
        chkxer("ZHETF2_RK", infot, nout, lerr, ok);
        infot = 2;
        zhetf2_rk("U", -1, a, 1, e, ip, info);
        chkxer("ZHETF2_RK", infot, nout, lerr, ok);
        infot = 4;
        zhetf2_rk("U", 2, a, 1, e, ip, info);
        chkxer("ZHETF2_RK", infot, nout, lerr, ok);
        //
        //        ZHETRI_3
        //
        srnamt = "ZHETRI_3";
        infot = 1;
        zhetri_3("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("ZHETRI_3", infot, nout, lerr, ok);
        infot = 2;
        zhetri_3("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("ZHETRI_3", infot, nout, lerr, ok);
        infot = 4;
        zhetri_3("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("ZHETRI_3", infot, nout, lerr, ok);
        infot = 8;
        zhetri_3("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("ZHETRI_3", infot, nout, lerr, ok);
        infot = 8;
        zhetri_3("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("ZHETRI_3", infot, nout, lerr, ok);
        //
        //        ZHETRI_3X
        //
        srnamt = "ZHETRI_3X";
        infot = 1;
        zhetri_3x("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("ZHETRI_3X", infot, nout, lerr, ok);
        infot = 2;
        zhetri_3x("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("ZHETRI_3X", infot, nout, lerr, ok);
        infot = 4;
        zhetri_3x("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("ZHETRI_3X", infot, nout, lerr, ok);
        //
        //        ZHETRS_3
        //
        srnamt = "ZHETRS_3";
        infot = 1;
        zhetrs_3("/", 0, 0, a, 1, e, ip, b, 1, info);
        chkxer("ZHETRS_3", infot, nout, lerr, ok);
        infot = 2;
        zhetrs_3("U", -1, 0, a, 1, e, ip, b, 1, info);
        chkxer("ZHETRS_3", infot, nout, lerr, ok);
        infot = 3;
        zhetrs_3("U", 0, -1, a, 1, e, ip, b, 1, info);
        chkxer("ZHETRS_3", infot, nout, lerr, ok);
        infot = 5;
        zhetrs_3("U", 2, 1, a, 1, e, ip, b, 2, info);
        chkxer("ZHETRS_3", infot, nout, lerr, ok);
        infot = 9;
        zhetrs_3("U", 2, 1, a, 2, e, ip, b, 1, info);
        chkxer("ZHETRS_3", infot, nout, lerr, ok);
        //
        //        ZHECON_3
        //
        srnamt = "ZHECON_3";
        infot = 1;
        zhecon_3("/", 0, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("ZHECON_3", infot, nout, lerr, ok);
        infot = 2;
        zhecon_3("U", -1, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("ZHECON_3", infot, nout, lerr, ok);
        infot = 4;
        zhecon_3("U", 2, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("ZHECON_3", infot, nout, lerr, ok);
        infot = 7;
        zhecon_3("U", 1, a, 1, e, ip, -1.0, rcond, w, info);
        chkxer("ZHECON_3", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "HP") {
        //
        //        Test error exits of the routines that use factorization
        //        of a Hermitian indefinite packed matrix with patrial
        //        (Bunch-Kaufman) diagonal pivoting method.
        //
        //        ZHPTRF
        //
        srnamt = "ZHPTRF";
        infot = 1;
        zhptrf("/", 0, a, ip, info);
        chkxer("ZHPTRF", infot, nout, lerr, ok);
        infot = 2;
        zhptrf("U", -1, a, ip, info);
        chkxer("ZHPTRF", infot, nout, lerr, ok);
        //
        //        ZHPTRI
        //
        srnamt = "ZHPTRI";
        infot = 1;
        zhptri("/", 0, a, ip, w, info);
        chkxer("ZHPTRI", infot, nout, lerr, ok);
        infot = 2;
        zhptri("U", -1, a, ip, w, info);
        chkxer("ZHPTRI", infot, nout, lerr, ok);
        //
        //        ZHPTRS
        //
        srnamt = "ZHPTRS";
        infot = 1;
        zhptrs("/", 0, 0, a, ip, b, 1, info);
        chkxer("ZHPTRS", infot, nout, lerr, ok);
        infot = 2;
        zhptrs("U", -1, 0, a, ip, b, 1, info);
        chkxer("ZHPTRS", infot, nout, lerr, ok);
        infot = 3;
        zhptrs("U", 0, -1, a, ip, b, 1, info);
        chkxer("ZHPTRS", infot, nout, lerr, ok);
        infot = 7;
        zhptrs("U", 2, 1, a, ip, b, 1, info);
        chkxer("ZHPTRS", infot, nout, lerr, ok);
        //
        //        ChprFS
        //
        srnamt = "ChprFS";
        infot = 1;
        Chprfs("/", 0, 0, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ChprFS", infot, nout, lerr, ok);
        infot = 2;
        Chprfs("U", -1, 0, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ChprFS", infot, nout, lerr, ok);
        infot = 3;
        Chprfs("U", 0, -1, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ChprFS", infot, nout, lerr, ok);
        infot = 8;
        Chprfs("U", 2, 1, a, af, ip, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("ChprFS", infot, nout, lerr, ok);
        infot = 10;
        Chprfs("U", 2, 1, a, af, ip, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("ChprFS", infot, nout, lerr, ok);
        //
        //        ZHPCON
        //
        srnamt = "ZHPCON";
        infot = 1;
        zhpcon("/", 0, a, ip, anrm, rcond, w, info);
        chkxer("ZHPCON", infot, nout, lerr, ok);
        infot = 2;
        zhpcon("U", -1, a, ip, anrm, rcond, w, info);
        chkxer("ZHPCON", infot, nout, lerr, ok);
        infot = 5;
        zhpcon("U", 1, a, ip, -anrm, rcond, w, info);
        chkxer("ZHPCON", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrhe
    //
}
