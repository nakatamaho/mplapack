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

void Rerrsy(common &cmn, const char *path, INTEGER const nunit) {
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
    arr_1d<3 * nmax, REAL> w(fill0);
    arr_1d<nmax, REAL> x(fill0);
    arr_1d<nmax, REAL> s(fill0);
    arr_1d<nmax, int> ip(fill0);
    arr_1d<nmax, int> iw(fill0);
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / (i + j).real();
            af[(i - 1) + (j - 1) * ldaf] = 1.0 / (i + j).real();
        }
        b[j - 1] = 0.0;
        e[j - 1] = 0.0;
        r1[j - 1] = 0.0;
        r2[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
        s[j - 1] = 0.0;
        ip[j - 1] = j;
        iw[j - 1] = j;
    }
    REAL anrm = 1.0;
    REAL rcond = 1.0;
    ok = true;
    //
    INTEGER info = 0;
    INTEGER n_err_bnds = 0;
    INTEGER nparams = 0;
    char eq = char0;
    REAL berr = 0.0;
    arr_2d<nmax, 3, REAL> err_bnds_n(fill0);
    arr_2d<nmax, 3, REAL> err_bnds_c(fill0);
    arr_1d<1, REAL> params(fill0);
    if (Mlsamen2, c2, "SY") {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with patrial
        //        (Bunch-Kaufman) pivoting.
        //
        //        DSYTRF
        //
        srnamt = "DSYTRF";
        infot = 1;
        dsytrf("/", 0, a, 1, ip, w, 1, info);
        chkxer("DSYTRF", infot, nout, lerr, ok);
        infot = 2;
        dsytrf("U", -1, a, 1, ip, w, 1, info);
        chkxer("DSYTRF", infot, nout, lerr, ok);
        infot = 4;
        dsytrf("U", 2, a, 1, ip, w, 4, info);
        chkxer("DSYTRF", infot, nout, lerr, ok);
        infot = 7;
        dsytrf("U", 0, a, 1, ip, w, 0, info);
        chkxer("DSYTRF", infot, nout, lerr, ok);
        infot = 7;
        dsytrf("U", 0, a, 1, ip, w, -2, info);
        chkxer("DSYTRF", infot, nout, lerr, ok);
        //
        //        DSYTF2
        //
        srnamt = "DSYTF2";
        infot = 1;
        dsytf2("/", 0, a, 1, ip, info);
        chkxer("DSYTF2", infot, nout, lerr, ok);
        infot = 2;
        dsytf2("U", -1, a, 1, ip, info);
        chkxer("DSYTF2", infot, nout, lerr, ok);
        infot = 4;
        dsytf2("U", 2, a, 1, ip, info);
        chkxer("DSYTF2", infot, nout, lerr, ok);
        //
        //        DSYTRI
        //
        srnamt = "DSYTRI";
        infot = 1;
        dsytri("/", 0, a, 1, ip, w, info);
        chkxer("DSYTRI", infot, nout, lerr, ok);
        infot = 2;
        dsytri("U", -1, a, 1, ip, w, info);
        chkxer("DSYTRI", infot, nout, lerr, ok);
        infot = 4;
        dsytri("U", 2, a, 1, ip, w, info);
        chkxer("DSYTRI", infot, nout, lerr, ok);
        //
        //        DSYTRI2
        //
        srnamt = "DSYTRI2";
        infot = 1;
        dsytri2("/", 0, a, 1, ip, w, iw, info);
        chkxer("DSYTRI2", infot, nout, lerr, ok);
        infot = 2;
        dsytri2("U", -1, a, 1, ip, w, iw, info);
        chkxer("DSYTRI2", infot, nout, lerr, ok);
        infot = 4;
        dsytri2("U", 2, a, 1, ip, w, iw, info);
        chkxer("DSYTRI2", infot, nout, lerr, ok);
        //
        //        DSYTRI2X
        //
        srnamt = "DSYTRI2X";
        infot = 1;
        dsytri2x("/", 0, a, 1, ip, w, 1, info);
        chkxer("DSYTRI2X", infot, nout, lerr, ok);
        infot = 2;
        dsytri2x("U", -1, a, 1, ip, w, 1, info);
        chkxer("DSYTRI2X", infot, nout, lerr, ok);
        infot = 4;
        dsytri2x("U", 2, a, 1, ip, w, 1, info);
        chkxer("DSYTRI2X", infot, nout, lerr, ok);
        //
        //        DSYTRS
        //
        srnamt = "DSYTRS";
        infot = 1;
        dsytrs("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("DSYTRS", infot, nout, lerr, ok);
        infot = 2;
        dsytrs("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("DSYTRS", infot, nout, lerr, ok);
        infot = 3;
        dsytrs("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("DSYTRS", infot, nout, lerr, ok);
        infot = 5;
        dsytrs("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("DSYTRS", infot, nout, lerr, ok);
        infot = 8;
        dsytrs("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("DSYTRS", infot, nout, lerr, ok);
        //
        //        RsyrFS
        //
        srnamt = "RsyrFS";
        infot = 1;
        Rsyrfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("RsyrFS", infot, nout, lerr, ok);
        infot = 2;
        Rsyrfs("U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("RsyrFS", infot, nout, lerr, ok);
        infot = 3;
        Rsyrfs("U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("RsyrFS", infot, nout, lerr, ok);
        infot = 5;
        Rsyrfs("U", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("RsyrFS", infot, nout, lerr, ok);
        infot = 7;
        Rsyrfs("U", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("RsyrFS", infot, nout, lerr, ok);
        infot = 10;
        Rsyrfs("U", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("RsyrFS", infot, nout, lerr, ok);
        infot = 12;
        Rsyrfs("U", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("RsyrFS", infot, nout, lerr, ok);
        //
        //        RsyrFSX
        //
        n_err_bnds = 3;
        nparams = 0;
        srnamt = "RsyrFSX";
        infot = 1;
        Rsyrfsx("/", eq, 0, 0, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RsyrFSX", infot, nout, lerr, ok);
        infot = 2;
        Rsyrfsx("U", eq, -1, 0, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RsyrFSX", infot, nout, lerr, ok);
        eq = "N";
        infot = 3;
        Rsyrfsx("U", eq, -1, 0, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RsyrFSX", infot, nout, lerr, ok);
        infot = 4;
        Rsyrfsx("U", eq, 0, -1, a, 1, af, 1, ip, s, b, 1, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RsyrFSX", infot, nout, lerr, ok);
        infot = 6;
        Rsyrfsx("U", eq, 2, 1, a, 1, af, 2, ip, s, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RsyrFSX", infot, nout, lerr, ok);
        infot = 8;
        Rsyrfsx("U", eq, 2, 1, a, 2, af, 1, ip, s, b, 2, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RsyrFSX", infot, nout, lerr, ok);
        infot = 12;
        Rsyrfsx("U", eq, 2, 1, a, 2, af, 2, ip, s, b, 1, x, 2, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RsyrFSX", infot, nout, lerr, ok);
        infot = 14;
        Rsyrfsx("U", eq, 2, 1, a, 2, af, 2, ip, s, b, 2, x, 1, rcond, berr, n_err_bnds, err_bnds_n, err_bnds_c, nparams, params, w, iw, info);
        chkxer("RsyrFSX", infot, nout, lerr, ok);
        //
        //        DSYCON
        //
        srnamt = "DSYCON";
        infot = 1;
        dsycon("/", 0, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("DSYCON", infot, nout, lerr, ok);
        infot = 2;
        dsycon("U", -1, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("DSYCON", infot, nout, lerr, ok);
        infot = 4;
        dsycon("U", 2, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("DSYCON", infot, nout, lerr, ok);
        infot = 6;
        dsycon("U", 1, a, 1, ip, -1.0, rcond, w, iw, info);
        chkxer("DSYCON", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SR") {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting.
        //
        //        DSYTRF_ROOK
        //
        srnamt = "DSYTRF_ROOK";
        infot = 1;
        dsytrf_rook("/", 0, a, 1, ip, w, 1, info);
        chkxer("DSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 2;
        dsytrf_rook("U", -1, a, 1, ip, w, 1, info);
        chkxer("DSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 4;
        dsytrf_rook("U", 2, a, 1, ip, w, 4, info);
        chkxer("DSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 7;
        dsytrf_rook("U", 0, a, 1, ip, w, 0, info);
        chkxer("DSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 7;
        dsytrf_rook("U", 0, a, 1, ip, w, -2, info);
        chkxer("DSYTRF_ROOK", infot, nout, lerr, ok);
        //
        //        DSYTF2_ROOK
        //
        srnamt = "DSYTF2_ROOK";
        infot = 1;
        dsytf2_rook("/", 0, a, 1, ip, info);
        chkxer("DSYTF2_ROOK", infot, nout, lerr, ok);
        infot = 2;
        dsytf2_rook("U", -1, a, 1, ip, info);
        chkxer("DSYTF2_ROOK", infot, nout, lerr, ok);
        infot = 4;
        dsytf2_rook("U", 2, a, 1, ip, info);
        chkxer("DSYTF2_ROOK", infot, nout, lerr, ok);
        //
        //        DSYTRI_ROOK
        //
        srnamt = "DSYTRI_ROOK";
        infot = 1;
        dsytri_rook("/", 0, a, 1, ip, w, info);
        chkxer("DSYTRI_ROOK", infot, nout, lerr, ok);
        infot = 2;
        dsytri_rook("U", -1, a, 1, ip, w, info);
        chkxer("DSYTRI_ROOK", infot, nout, lerr, ok);
        infot = 4;
        dsytri_rook("U", 2, a, 1, ip, w, info);
        chkxer("DSYTRI_ROOK", infot, nout, lerr, ok);
        //
        //        DSYTRS_ROOK
        //
        srnamt = "DSYTRS_ROOK";
        infot = 1;
        dsytrs_rook("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("DSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 2;
        dsytrs_rook("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("DSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 3;
        dsytrs_rook("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("DSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 5;
        dsytrs_rook("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("DSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 8;
        dsytrs_rook("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("DSYTRS_ROOK", infot, nout, lerr, ok);
        //
        //        DSYCON_ROOK
        //
        srnamt = "DSYCON_ROOK";
        infot = 1;
        dsycon_rook("/", 0, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("DSYCON_ROOK", infot, nout, lerr, ok);
        infot = 2;
        dsycon_rook("U", -1, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("DSYCON_ROOK", infot, nout, lerr, ok);
        infot = 4;
        dsycon_rook("U", 2, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("DSYCON_ROOK", infot, nout, lerr, ok);
        infot = 6;
        dsycon_rook("U", 1, a, 1, ip, -1.0, rcond, w, iw, info);
        chkxer("DSYCON_ROOK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SK") {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting with the new storage
        //        format for factors L ( or U) and D.
        //
        //        L (or U) is stored in A, diagonal of D is stored on the
        //        diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        //        DSYTRF_RK
        //
        srnamt = "DSYTRF_RK";
        infot = 1;
        dsytrf_rk("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("DSYTRF_RK", infot, nout, lerr, ok);
        infot = 2;
        dsytrf_rk("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("DSYTRF_RK", infot, nout, lerr, ok);
        infot = 4;
        dsytrf_rk("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("DSYTRF_RK", infot, nout, lerr, ok);
        infot = 8;
        dsytrf_rk("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("DSYTRF_RK", infot, nout, lerr, ok);
        infot = 8;
        dsytrf_rk("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("DSYTRF_RK", infot, nout, lerr, ok);
        //
        //        DSYTF2_RK
        //
        srnamt = "DSYTF2_RK";
        infot = 1;
        dsytf2_rk("/", 0, a, 1, e, ip, info);
        chkxer("DSYTF2_RK", infot, nout, lerr, ok);
        infot = 2;
        dsytf2_rk("U", -1, a, 1, e, ip, info);
        chkxer("DSYTF2_RK", infot, nout, lerr, ok);
        infot = 4;
        dsytf2_rk("U", 2, a, 1, e, ip, info);
        chkxer("DSYTF2_RK", infot, nout, lerr, ok);
        //
        //        DSYTRI_3
        //
        srnamt = "DSYTRI_3";
        infot = 1;
        dsytri_3("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("DSYTRI_3", infot, nout, lerr, ok);
        infot = 2;
        dsytri_3("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("DSYTRI_3", infot, nout, lerr, ok);
        infot = 4;
        dsytri_3("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("DSYTRI_3", infot, nout, lerr, ok);
        infot = 8;
        dsytri_3("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("DSYTRI_3", infot, nout, lerr, ok);
        infot = 8;
        dsytri_3("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("DSYTRI_3", infot, nout, lerr, ok);
        //
        //        DSYTRI_3X
        //
        srnamt = "DSYTRI_3X";
        infot = 1;
        dsytri_3x("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("DSYTRI_3X", infot, nout, lerr, ok);
        infot = 2;
        dsytri_3x("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("DSYTRI_3X", infot, nout, lerr, ok);
        infot = 4;
        dsytri_3x("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("DSYTRI_3X", infot, nout, lerr, ok);
        //
        //        DSYTRS_3
        //
        srnamt = "DSYTRS_3";
        infot = 1;
        dsytrs_3("/", 0, 0, a, 1, e, ip, b, 1, info);
        chkxer("DSYTRS_3", infot, nout, lerr, ok);
        infot = 2;
        dsytrs_3("U", -1, 0, a, 1, e, ip, b, 1, info);
        chkxer("DSYTRS_3", infot, nout, lerr, ok);
        infot = 3;
        dsytrs_3("U", 0, -1, a, 1, e, ip, b, 1, info);
        chkxer("DSYTRS_3", infot, nout, lerr, ok);
        infot = 5;
        dsytrs_3("U", 2, 1, a, 1, e, ip, b, 2, info);
        chkxer("DSYTRS_3", infot, nout, lerr, ok);
        infot = 9;
        dsytrs_3("U", 2, 1, a, 2, e, ip, b, 1, info);
        chkxer("DSYTRS_3", infot, nout, lerr, ok);
        //
        //        DSYCON_3
        //
        srnamt = "DSYCON_3";
        infot = 1;
        dsycon_3("/", 0, a, 1, e, ip, anrm, rcond, w, iw, info);
        chkxer("DSYCON_3", infot, nout, lerr, ok);
        infot = 2;
        dsycon_3("U", -1, a, 1, e, ip, anrm, rcond, w, iw, info);
        chkxer("DSYCON_3", infot, nout, lerr, ok);
        infot = 4;
        dsycon_3("U", 2, a, 1, e, ip, anrm, rcond, w, iw, info);
        chkxer("DSYCON_3", infot, nout, lerr, ok);
        infot = 7;
        dsycon_3("U", 1, a, 1, e, ip, -1.0, rcond, w, iw, info);
        chkxer("DSYCON_3", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SP") {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite packed matrix with patrial
        //        (Bunch-Kaufman) pivoting.
        //
        //        DSPTRF
        //
        srnamt = "DSPTRF";
        infot = 1;
        dsptrf("/", 0, a, ip, info);
        chkxer("DSPTRF", infot, nout, lerr, ok);
        infot = 2;
        dsptrf("U", -1, a, ip, info);
        chkxer("DSPTRF", infot, nout, lerr, ok);
        //
        //        DSPTRI
        //
        srnamt = "DSPTRI";
        infot = 1;
        dsptri("/", 0, a, ip, w, info);
        chkxer("DSPTRI", infot, nout, lerr, ok);
        infot = 2;
        dsptri("U", -1, a, ip, w, info);
        chkxer("DSPTRI", infot, nout, lerr, ok);
        //
        //        DSPTRS
        //
        srnamt = "DSPTRS";
        infot = 1;
        dsptrs("/", 0, 0, a, ip, b, 1, info);
        chkxer("DSPTRS", infot, nout, lerr, ok);
        infot = 2;
        dsptrs("U", -1, 0, a, ip, b, 1, info);
        chkxer("DSPTRS", infot, nout, lerr, ok);
        infot = 3;
        dsptrs("U", 0, -1, a, ip, b, 1, info);
        chkxer("DSPTRS", infot, nout, lerr, ok);
        infot = 7;
        dsptrs("U", 2, 1, a, ip, b, 1, info);
        chkxer("DSPTRS", infot, nout, lerr, ok);
        //
        //        RsprFS
        //
        srnamt = "RsprFS";
        infot = 1;
        Rsprfs("/", 0, 0, a, af, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("RsprFS", infot, nout, lerr, ok);
        infot = 2;
        Rsprfs("U", -1, 0, a, af, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("RsprFS", infot, nout, lerr, ok);
        infot = 3;
        Rsprfs("U", 0, -1, a, af, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("RsprFS", infot, nout, lerr, ok);
        infot = 8;
        Rsprfs("U", 2, 1, a, af, ip, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("RsprFS", infot, nout, lerr, ok);
        infot = 10;
        Rsprfs("U", 2, 1, a, af, ip, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("RsprFS", infot, nout, lerr, ok);
        //
        //        DSPCON
        //
        srnamt = "DSPCON";
        infot = 1;
        dspcon("/", 0, a, ip, anrm, rcond, w, iw, info);
        chkxer("DSPCON", infot, nout, lerr, ok);
        infot = 2;
        dspcon("U", -1, a, ip, anrm, rcond, w, iw, info);
        chkxer("DSPCON", infot, nout, lerr, ok);
        infot = 5;
        dspcon("U", 1, a, ip, -1.0, rcond, w, iw, info);
        chkxer("DSPCON", infot, nout, lerr, ok);
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrsy
    //
}
