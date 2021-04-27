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
    if (Mlsamen(2, c2, "HE")) {
        //
        //        Chetrf
        //
        srnamt = "Chetrf";
        infot = 1;
        Chetrf("/", 0, a, 1, ip, w, 1, info);
        chkxer("Chetrf", infot, nout, lerr, ok);
        infot = 2;
        Chetrf("U", -1, a, 1, ip, w, 1, info);
        chkxer("Chetrf", infot, nout, lerr, ok);
        infot = 4;
        Chetrf("U", 2, a, 1, ip, w, 4, info);
        chkxer("Chetrf", infot, nout, lerr, ok);
        infot = 7;
        Chetrf("U", 0, a, 1, ip, w, 0, info);
        chkxer("Chetrf", infot, nout, lerr, ok);
        infot = 7;
        Chetrf("U", 0, a, 1, ip, w, -2, info);
        chkxer("Chetrf", infot, nout, lerr, ok);
        //
        //        Chetf2
        //
        srnamt = "Chetf2";
        infot = 1;
        Chetf2("/", 0, a, 1, ip, info);
        chkxer("Chetf2", infot, nout, lerr, ok);
        infot = 2;
        Chetf2("U", -1, a, 1, ip, info);
        chkxer("Chetf2", infot, nout, lerr, ok);
        infot = 4;
        Chetf2("U", 2, a, 1, ip, info);
        chkxer("Chetf2", infot, nout, lerr, ok);
        //
        //        Chetri
        //
        srnamt = "Chetri";
        infot = 1;
        Chetri("/", 0, a, 1, ip, w, info);
        chkxer("Chetri", infot, nout, lerr, ok);
        infot = 2;
        Chetri("U", -1, a, 1, ip, w, info);
        chkxer("Chetri", infot, nout, lerr, ok);
        infot = 4;
        Chetri("U", 2, a, 1, ip, w, info);
        chkxer("Chetri", infot, nout, lerr, ok);
        //
        //        Chetri2
        //
        srnamt = "Chetri2";
        infot = 1;
        Chetri2("/", 0, a, 1, ip, w, 1, info);
        chkxer("Chetri2", infot, nout, lerr, ok);
        infot = 2;
        Chetri2("U", -1, a, 1, ip, w, 1, info);
        chkxer("Chetri2", infot, nout, lerr, ok);
        infot = 4;
        Chetri2("U", 2, a, 1, ip, w, 1, info);
        chkxer("Chetri2", infot, nout, lerr, ok);
        //
        //        Chetri2x
        //
        srnamt = "Chetri2x";
        infot = 1;
        Chetri2x("/", 0, a, 1, ip, w, 1, info);
        chkxer("Chetri2x", infot, nout, lerr, ok);
        infot = 2;
        Chetri2x("U", -1, a, 1, ip, w, 1, info);
        chkxer("Chetri2x", infot, nout, lerr, ok);
        infot = 4;
        Chetri2x("U", 2, a, 1, ip, w, 1, info);
        chkxer("Chetri2x", infot, nout, lerr, ok);
        //
        //        Chetrs
        //
        srnamt = "Chetrs";
        infot = 1;
        Chetrs("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("Chetrs", infot, nout, lerr, ok);
        infot = 2;
        Chetrs("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("Chetrs", infot, nout, lerr, ok);
        infot = 3;
        Chetrs("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("Chetrs", infot, nout, lerr, ok);
        infot = 5;
        Chetrs("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("Chetrs", infot, nout, lerr, ok);
        infot = 8;
        Chetrs("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("Chetrs", infot, nout, lerr, ok);
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
        //        Checon
        //
        srnamt = "Checon";
        infot = 1;
        Checon("/", 0, a, 1, ip, anrm, rcond, w, info);
        chkxer("Checon", infot, nout, lerr, ok);
        infot = 2;
        Checon("U", -1, a, 1, ip, anrm, rcond, w, info);
        chkxer("Checon", infot, nout, lerr, ok);
        infot = 4;
        Checon("U", 2, a, 1, ip, anrm, rcond, w, info);
        chkxer("Checon", infot, nout, lerr, ok);
        infot = 6;
        Checon("U", 1, a, 1, ip, -anrm, rcond, w, info);
        chkxer("Checon", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "HR")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a Hermitian indefinite matrix with rook
        //        (bounded Bunch-Kaufman) diagonal pivoting method.
        //
        //        Chetrf_rook
        //
        srnamt = "Chetrf_rook";
        infot = 1;
        Chetrf_rook("/", 0, a, 1, ip, w, 1, info);
        chkxer("Chetrf_rook", infot, nout, lerr, ok);
        infot = 2;
        Chetrf_rook("U", -1, a, 1, ip, w, 1, info);
        chkxer("Chetrf_rook", infot, nout, lerr, ok);
        infot = 4;
        Chetrf_rook("U", 2, a, 1, ip, w, 4, info);
        chkxer("Chetrf_rook", infot, nout, lerr, ok);
        infot = 7;
        Chetrf_rook("U", 0, a, 1, ip, w, 0, info);
        chkxer("Chetrf_rook", infot, nout, lerr, ok);
        infot = 7;
        Chetrf_rook("U", 0, a, 1, ip, w, -2, info);
        chkxer("Chetrf_rook", infot, nout, lerr, ok);
        //
        //        Chetf2_rook
        //
        srnamt = "Chetf2_rook";
        infot = 1;
        Chetf2_rook("/", 0, a, 1, ip, info);
        chkxer("Chetf2_rook", infot, nout, lerr, ok);
        infot = 2;
        Chetf2_rook("U", -1, a, 1, ip, info);
        chkxer("Chetf2_rook", infot, nout, lerr, ok);
        infot = 4;
        Chetf2_rook("U", 2, a, 1, ip, info);
        chkxer("Chetf2_rook", infot, nout, lerr, ok);
        //
        //        Chetri_rook
        //
        srnamt = "Chetri_rook";
        infot = 1;
        Chetri_rook("/", 0, a, 1, ip, w, info);
        chkxer("Chetri_rook", infot, nout, lerr, ok);
        infot = 2;
        Chetri_rook("U", -1, a, 1, ip, w, info);
        chkxer("Chetri_rook", infot, nout, lerr, ok);
        infot = 4;
        Chetri_rook("U", 2, a, 1, ip, w, info);
        chkxer("Chetri_rook", infot, nout, lerr, ok);
        //
        //        Chetrs_rook
        //
        srnamt = "Chetrs_rook";
        infot = 1;
        Chetrs_rook("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("Chetrs_rook", infot, nout, lerr, ok);
        infot = 2;
        Chetrs_rook("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("Chetrs_rook", infot, nout, lerr, ok);
        infot = 3;
        Chetrs_rook("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("Chetrs_rook", infot, nout, lerr, ok);
        infot = 5;
        Chetrs_rook("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("Chetrs_rook", infot, nout, lerr, ok);
        infot = 8;
        Chetrs_rook("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("Chetrs_rook", infot, nout, lerr, ok);
        //
        //        Checon_rook
        //
        srnamt = "Checon_rook";
        infot = 1;
        Checon_rook("/", 0, a, 1, ip, anrm, rcond, w, info);
        chkxer("Checon_rook", infot, nout, lerr, ok);
        infot = 2;
        Checon_rook("U", -1, a, 1, ip, anrm, rcond, w, info);
        chkxer("Checon_rook", infot, nout, lerr, ok);
        infot = 4;
        Checon_rook("U", 2, a, 1, ip, anrm, rcond, w, info);
        chkxer("Checon_rook", infot, nout, lerr, ok);
        infot = 6;
        Checon_rook("U", 1, a, 1, ip, -anrm, rcond, w, info);
        chkxer("Checon_rook", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "HK")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting with the new storage
        //        format for factors L ( or U) and D.
        //
        //        L (or U) is stored in A, diagonal of D is stored on the
        //        diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        //        Chetrf_rk
        //
        srnamt = "Chetrf_rk";
        infot = 1;
        Chetrf_rk("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("Chetrf_rk", infot, nout, lerr, ok);
        infot = 2;
        Chetrf_rk("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("Chetrf_rk", infot, nout, lerr, ok);
        infot = 4;
        Chetrf_rk("U", 2, a, 1, e, ip, w, 4, info);
        chkxer("Chetrf_rk", infot, nout, lerr, ok);
        infot = 8;
        Chetrf_rk("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("Chetrf_rk", infot, nout, lerr, ok);
        infot = 8;
        Chetrf_rk("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("Chetrf_rk", infot, nout, lerr, ok);
        //
        //        Chetf2_rk
        //
        srnamt = "Chetf2_rk";
        infot = 1;
        Chetf2_rk("/", 0, a, 1, e, ip, info);
        chkxer("Chetf2_rk", infot, nout, lerr, ok);
        infot = 2;
        Chetf2_rk("U", -1, a, 1, e, ip, info);
        chkxer("Chetf2_rk", infot, nout, lerr, ok);
        infot = 4;
        Chetf2_rk("U", 2, a, 1, e, ip, info);
        chkxer("Chetf2_rk", infot, nout, lerr, ok);
        //
        //        Chetri_3
        //
        srnamt = "Chetri_3";
        infot = 1;
        Chetri_3("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("Chetri_3", infot, nout, lerr, ok);
        infot = 2;
        Chetri_3("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("Chetri_3", infot, nout, lerr, ok);
        infot = 4;
        Chetri_3("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("Chetri_3", infot, nout, lerr, ok);
        infot = 8;
        Chetri_3("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("Chetri_3", infot, nout, lerr, ok);
        infot = 8;
        Chetri_3("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("Chetri_3", infot, nout, lerr, ok);
        //
        //        Chetri_3x
        //
        srnamt = "Chetri_3x";
        infot = 1;
        Chetri_3x("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("Chetri_3x", infot, nout, lerr, ok);
        infot = 2;
        Chetri_3x("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("Chetri_3x", infot, nout, lerr, ok);
        infot = 4;
        Chetri_3x("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("Chetri_3x", infot, nout, lerr, ok);
        //
        //        Chetrs_3
        //
        srnamt = "Chetrs_3";
        infot = 1;
        Chetrs_3("/", 0, 0, a, 1, e, ip, b, 1, info);
        chkxer("Chetrs_3", infot, nout, lerr, ok);
        infot = 2;
        Chetrs_3("U", -1, 0, a, 1, e, ip, b, 1, info);
        chkxer("Chetrs_3", infot, nout, lerr, ok);
        infot = 3;
        Chetrs_3("U", 0, -1, a, 1, e, ip, b, 1, info);
        chkxer("Chetrs_3", infot, nout, lerr, ok);
        infot = 5;
        Chetrs_3("U", 2, 1, a, 1, e, ip, b, 2, info);
        chkxer("Chetrs_3", infot, nout, lerr, ok);
        infot = 9;
        Chetrs_3("U", 2, 1, a, 2, e, ip, b, 1, info);
        chkxer("Chetrs_3", infot, nout, lerr, ok);
        //
        //        Checon_3
        //
        srnamt = "Checon_3";
        infot = 1;
        Checon_3("/", 0, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("Checon_3", infot, nout, lerr, ok);
        infot = 2;
        Checon_3("U", -1, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("Checon_3", infot, nout, lerr, ok);
        infot = 4;
        Checon_3("U", 2, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("Checon_3", infot, nout, lerr, ok);
        infot = 7;
        Checon_3("U", 1, a, 1, e, ip, -1.0, rcond, w, info);
        chkxer("Checon_3", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "HP")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a Hermitian indefinite packed matrix with patrial
        //        (Bunch-Kaufman) diagonal pivoting method.
        //
        //        Chptrf
        //
        srnamt = "Chptrf";
        infot = 1;
        Chptrf("/", 0, a, ip, info);
        chkxer("Chptrf", infot, nout, lerr, ok);
        infot = 2;
        Chptrf("U", -1, a, ip, info);
        chkxer("Chptrf", infot, nout, lerr, ok);
        //
        //        Chptri
        //
        srnamt = "Chptri";
        infot = 1;
        Chptri("/", 0, a, ip, w, info);
        chkxer("Chptri", infot, nout, lerr, ok);
        infot = 2;
        Chptri("U", -1, a, ip, w, info);
        chkxer("Chptri", infot, nout, lerr, ok);
        //
        //        Chptrs
        //
        srnamt = "Chptrs";
        infot = 1;
        Chptrs("/", 0, 0, a, ip, b, 1, info);
        chkxer("Chptrs", infot, nout, lerr, ok);
        infot = 2;
        Chptrs("U", -1, 0, a, ip, b, 1, info);
        chkxer("Chptrs", infot, nout, lerr, ok);
        infot = 3;
        Chptrs("U", 0, -1, a, ip, b, 1, info);
        chkxer("Chptrs", infot, nout, lerr, ok);
        infot = 7;
        Chptrs("U", 2, 1, a, ip, b, 1, info);
        chkxer("Chptrs", infot, nout, lerr, ok);
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
        //        Chpcon
        //
        srnamt = "Chpcon";
        infot = 1;
        Chpcon("/", 0, a, ip, anrm, rcond, w, info);
        chkxer("Chpcon", infot, nout, lerr, ok);
        infot = 2;
        Chpcon("U", -1, a, ip, anrm, rcond, w, info);
        chkxer("Chpcon", infot, nout, lerr, ok);
        infot = 5;
        Chpcon("U", 1, a, ip, -anrm, rcond, w, info);
        chkxer("Chpcon", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrhe
    //
}
