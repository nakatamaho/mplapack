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
        ip[j - 1] = j;
        iw[j - 1] = j;
    }
    REAL anrm = 1.0;
    REAL rcond = 1.0;
    ok = true;
    //
    INTEGER info = 0;
    if (Mlsamen(2, c2, "SY")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with patrial
        //        (Bunch-Kaufman) pivoting.
        //
        //        Rsytrf
        //
        srnamt = "Rsytrf";
        infot = 1;
        Rsytrf("/", 0, a, 1, ip, w, 1, info);
        chkxer("Rsytrf", infot, nout, lerr, ok);
        infot = 2;
        Rsytrf("U", -1, a, 1, ip, w, 1, info);
        chkxer("Rsytrf", infot, nout, lerr, ok);
        infot = 4;
        Rsytrf("U", 2, a, 1, ip, w, 4, info);
        chkxer("Rsytrf", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf("U", 0, a, 1, ip, w, 0, info);
        chkxer("Rsytrf", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf("U", 0, a, 1, ip, w, -2, info);
        chkxer("Rsytrf", infot, nout, lerr, ok);
        //
        //        Rsytf2
        //
        srnamt = "Rsytf2";
        infot = 1;
        Rsytf2("/", 0, a, 1, ip, info);
        chkxer("Rsytf2", infot, nout, lerr, ok);
        infot = 2;
        Rsytf2("U", -1, a, 1, ip, info);
        chkxer("Rsytf2", infot, nout, lerr, ok);
        infot = 4;
        Rsytf2("U", 2, a, 1, ip, info);
        chkxer("Rsytf2", infot, nout, lerr, ok);
        //
        //        Rsytri
        //
        srnamt = "Rsytri";
        infot = 1;
        Rsytri("/", 0, a, 1, ip, w, info);
        chkxer("Rsytri", infot, nout, lerr, ok);
        infot = 2;
        Rsytri("U", -1, a, 1, ip, w, info);
        chkxer("Rsytri", infot, nout, lerr, ok);
        infot = 4;
        Rsytri("U", 2, a, 1, ip, w, info);
        chkxer("Rsytri", infot, nout, lerr, ok);
        //
        //        Rsytri2
        //
        srnamt = "Rsytri2";
        infot = 1;
        Rsytri2("/", 0, a, 1, ip, w, iw[1 - 1], info);
        chkxer("Rsytri2", infot, nout, lerr, ok);
        infot = 2;
        Rsytri2("U", -1, a, 1, ip, w, iw[1 - 1], info);
        chkxer("Rsytri2", infot, nout, lerr, ok);
        infot = 4;
        Rsytri2("U", 2, a, 1, ip, w, iw[1 - 1], info);
        chkxer("Rsytri2", infot, nout, lerr, ok);
        //
        //        Rsytri2x
        //
        srnamt = "Rsytri2x";
        infot = 1;
        Rsytri2x("/", 0, a, 1, ip, w, 1, info);
        chkxer("Rsytri2x", infot, nout, lerr, ok);
        infot = 2;
        Rsytri2x("U", -1, a, 1, ip, w, 1, info);
        chkxer("Rsytri2x", infot, nout, lerr, ok);
        infot = 4;
        Rsytri2x("U", 2, a, 1, ip, w, 1, info);
        chkxer("Rsytri2x", infot, nout, lerr, ok);
        //
        //        Rsytrs
        //
        srnamt = "Rsytrs";
        infot = 1;
        Rsytrs("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("Rsytrs", infot, nout, lerr, ok);
        infot = 2;
        Rsytrs("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("Rsytrs", infot, nout, lerr, ok);
        infot = 3;
        Rsytrs("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("Rsytrs", infot, nout, lerr, ok);
        infot = 5;
        Rsytrs("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("Rsytrs", infot, nout, lerr, ok);
        infot = 8;
        Rsytrs("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("Rsytrs", infot, nout, lerr, ok);
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
        //        Rsycon
        //
        srnamt = "Rsycon";
        infot = 1;
        Rsycon("/", 0, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("Rsycon", infot, nout, lerr, ok);
        infot = 2;
        Rsycon("U", -1, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("Rsycon", infot, nout, lerr, ok);
        infot = 4;
        Rsycon("U", 2, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("Rsycon", infot, nout, lerr, ok);
        infot = 6;
        Rsycon("U", 1, a, 1, ip, -1.0, rcond, w, iw, info);
        chkxer("Rsycon", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SR")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting.
        //
        //        Rsytrf_rook
        //
        srnamt = "Rsytrf_rook";
        infot = 1;
        Rsytrf_rook("/", 0, a, 1, ip, w, 1, info);
        chkxer("Rsytrf_rook", infot, nout, lerr, ok);
        infot = 2;
        Rsytrf_rook("U", -1, a, 1, ip, w, 1, info);
        chkxer("Rsytrf_rook", infot, nout, lerr, ok);
        infot = 4;
        Rsytrf_rook("U", 2, a, 1, ip, w, 4, info);
        chkxer("Rsytrf_rook", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf_rook("U", 0, a, 1, ip, w, 0, info);
        chkxer("Rsytrf_rook", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf_rook("U", 0, a, 1, ip, w, -2, info);
        chkxer("Rsytrf_rook", infot, nout, lerr, ok);
        //
        //        Rsytf2_rook
        //
        srnamt = "Rsytf2_rook";
        infot = 1;
        Rsytf2_rook("/", 0, a, 1, ip, info);
        chkxer("Rsytf2_rook", infot, nout, lerr, ok);
        infot = 2;
        Rsytf2_rook("U", -1, a, 1, ip, info);
        chkxer("Rsytf2_rook", infot, nout, lerr, ok);
        infot = 4;
        Rsytf2_rook("U", 2, a, 1, ip, info);
        chkxer("Rsytf2_rook", infot, nout, lerr, ok);
        //
        //        Rsytri_rook
        //
        srnamt = "Rsytri_rook";
        infot = 1;
        Rsytri_rook("/", 0, a, 1, ip, w, info);
        chkxer("Rsytri_rook", infot, nout, lerr, ok);
        infot = 2;
        Rsytri_rook("U", -1, a, 1, ip, w, info);
        chkxer("Rsytri_rook", infot, nout, lerr, ok);
        infot = 4;
        Rsytri_rook("U", 2, a, 1, ip, w, info);
        chkxer("Rsytri_rook", infot, nout, lerr, ok);
        //
        //        Rsytrs_rook
        //
        srnamt = "Rsytrs_rook";
        infot = 1;
        Rsytrs_rook("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("Rsytrs_rook", infot, nout, lerr, ok);
        infot = 2;
        Rsytrs_rook("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("Rsytrs_rook", infot, nout, lerr, ok);
        infot = 3;
        Rsytrs_rook("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("Rsytrs_rook", infot, nout, lerr, ok);
        infot = 5;
        Rsytrs_rook("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("Rsytrs_rook", infot, nout, lerr, ok);
        infot = 8;
        Rsytrs_rook("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("Rsytrs_rook", infot, nout, lerr, ok);
        //
        //        Rsycon_rook
        //
        srnamt = "Rsycon_rook";
        infot = 1;
        Rsycon_rook("/", 0, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("Rsycon_rook", infot, nout, lerr, ok);
        infot = 2;
        Rsycon_rook("U", -1, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("Rsycon_rook", infot, nout, lerr, ok);
        infot = 4;
        Rsycon_rook("U", 2, a, 1, ip, anrm, rcond, w, iw, info);
        chkxer("Rsycon_rook", infot, nout, lerr, ok);
        infot = 6;
        Rsycon_rook("U", 1, a, 1, ip, -1.0, rcond, w, iw, info);
        chkxer("Rsycon_rook", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SK")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting with the new storage
        //        format for factors L ( or U) and D.
        //
        //        L (or U) is stored in A, diagonal of D is stored on the
        //        diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        //        Rsytrf_rk
        //
        srnamt = "Rsytrf_rk";
        infot = 1;
        Rsytrf_rk("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("Rsytrf_rk", infot, nout, lerr, ok);
        infot = 2;
        Rsytrf_rk("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("Rsytrf_rk", infot, nout, lerr, ok);
        infot = 4;
        Rsytrf_rk("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("Rsytrf_rk", infot, nout, lerr, ok);
        infot = 8;
        Rsytrf_rk("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("Rsytrf_rk", infot, nout, lerr, ok);
        infot = 8;
        Rsytrf_rk("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("Rsytrf_rk", infot, nout, lerr, ok);
        //
        //        Rsytf2_rk
        //
        srnamt = "Rsytf2_rk";
        infot = 1;
        Rsytf2_rk("/", 0, a, 1, e, ip, info);
        chkxer("Rsytf2_rk", infot, nout, lerr, ok);
        infot = 2;
        Rsytf2_rk("U", -1, a, 1, e, ip, info);
        chkxer("Rsytf2_rk", infot, nout, lerr, ok);
        infot = 4;
        Rsytf2_rk("U", 2, a, 1, e, ip, info);
        chkxer("Rsytf2_rk", infot, nout, lerr, ok);
        //
        //        Rsytri_3
        //
        srnamt = "Rsytri_3";
        infot = 1;
        Rsytri_3("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("Rsytri_3", infot, nout, lerr, ok);
        infot = 2;
        Rsytri_3("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("Rsytri_3", infot, nout, lerr, ok);
        infot = 4;
        Rsytri_3("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("Rsytri_3", infot, nout, lerr, ok);
        infot = 8;
        Rsytri_3("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("Rsytri_3", infot, nout, lerr, ok);
        infot = 8;
        Rsytri_3("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("Rsytri_3", infot, nout, lerr, ok);
        //
        //        Rsytri_3x
        //
        srnamt = "Rsytri_3x";
        infot = 1;
        Rsytri_3x("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("Rsytri_3x", infot, nout, lerr, ok);
        infot = 2;
        Rsytri_3x("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("Rsytri_3x", infot, nout, lerr, ok);
        infot = 4;
        Rsytri_3x("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("Rsytri_3x", infot, nout, lerr, ok);
        //
        //        Rsytrs_3
        //
        srnamt = "Rsytrs_3";
        infot = 1;
        Rsytrs_3("/", 0, 0, a, 1, e, ip, b, 1, info);
        chkxer("Rsytrs_3", infot, nout, lerr, ok);
        infot = 2;
        Rsytrs_3("U", -1, 0, a, 1, e, ip, b, 1, info);
        chkxer("Rsytrs_3", infot, nout, lerr, ok);
        infot = 3;
        Rsytrs_3("U", 0, -1, a, 1, e, ip, b, 1, info);
        chkxer("Rsytrs_3", infot, nout, lerr, ok);
        infot = 5;
        Rsytrs_3("U", 2, 1, a, 1, e, ip, b, 2, info);
        chkxer("Rsytrs_3", infot, nout, lerr, ok);
        infot = 9;
        Rsytrs_3("U", 2, 1, a, 2, e, ip, b, 1, info);
        chkxer("Rsytrs_3", infot, nout, lerr, ok);
        //
        //        Rsycon_3
        //
        srnamt = "Rsycon_3";
        infot = 1;
        Rsycon_3("/", 0, a, 1, e, ip, anrm, rcond, w, iw, info);
        chkxer("Rsycon_3", infot, nout, lerr, ok);
        infot = 2;
        Rsycon_3("U", -1, a, 1, e, ip, anrm, rcond, w, iw, info);
        chkxer("Rsycon_3", infot, nout, lerr, ok);
        infot = 4;
        Rsycon_3("U", 2, a, 1, e, ip, anrm, rcond, w, iw, info);
        chkxer("Rsycon_3", infot, nout, lerr, ok);
        infot = 7;
        Rsycon_3("U", 1, a, 1, e, ip, -1.0, rcond, w, iw, info);
        chkxer("Rsycon_3", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SA")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with Aasen's algorithm.
        //
        //        Rsytrf_aa
        //
        srnamt = "Rsytrf_aa";
        infot = 1;
        Rsytrf_aa("/", 0, a, 1, ip, w, 1, info);
        chkxer("Rsytrf_aa", infot, nout, lerr, ok);
        infot = 2;
        Rsytrf_aa("U", -1, a, 1, ip, w, 1, info);
        chkxer("Rsytrf_aa", infot, nout, lerr, ok);
        infot = 4;
        Rsytrf_aa("U", 2, a, 1, ip, w, 4, info);
        chkxer("Rsytrf_aa", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf_aa("U", 0, a, 1, ip, w, 0, info);
        chkxer("Rsytrf_aa", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf_aa("U", 0, a, 1, ip, w, -2, info);
        chkxer("Rsytrf_aa", infot, nout, lerr, ok);
        //
        //        Rsytrs_aa
        //
        srnamt = "Rsytrs_aa";
        infot = 1;
        Rsytrs_aa("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsytrs_aa", infot, nout, lerr, ok);
        infot = 2;
        Rsytrs_aa("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsytrs_aa", infot, nout, lerr, ok);
        infot = 3;
        Rsytrs_aa("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("Rsytrs_aa", infot, nout, lerr, ok);
        infot = 5;
        Rsytrs_aa("U", 2, 1, a, 1, ip, b, 2, w, 1, info);
        chkxer("Rsytrs_aa", infot, nout, lerr, ok);
        infot = 8;
        Rsytrs_aa("U", 2, 1, a, 2, ip, b, 1, w, 1, info);
        chkxer("Rsytrs_aa", infot, nout, lerr, ok);
        infot = 10;
        Rsytrs_aa("U", 0, 1, a, 2, ip, b, 1, w, 0, info);
        chkxer("Rsytrs_aa", infot, nout, lerr, ok);
        infot = 10;
        Rsytrs_aa("U", 0, 1, a, 2, ip, b, 1, w, -2, info);
        chkxer("Rsytrs_aa", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "S2")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with Aasen's algorithm.
        //
        //        Rsytrf_aa_2stage
        //
        srnamt = "Rsytrf_aa_2stage";
        infot = 1;
        Rsytrf_aa_2stage("/", 0, a, 1, a, 1, ip, ip, w, 1, info);
        chkxer("Rsytrf_aa_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsytrf_aa_2stage("U", -1, a, 1, a, 1, ip, ip, w, 1, info);
        chkxer("Rsytrf_aa_2stage", infot, nout, lerr, ok);
        infot = 4;
        Rsytrf_aa_2stage("U", 2, a, 1, a, 2, ip, ip, w, 1, info);
        chkxer("Rsytrf_aa_2stage", infot, nout, lerr, ok);
        infot = 6;
        Rsytrf_aa_2stage("U", 2, a, 2, a, 1, ip, ip, w, 1, info);
        chkxer("Rsytrf_aa_2stage", infot, nout, lerr, ok);
        infot = 10;
        Rsytrf_aa_2stage("U", 2, a, 2, a, 8, ip, ip, w, 0, info);
        chkxer("Rsytrf_aa_2stage", infot, nout, lerr, ok);
        //
        //        Rsytrs_aa_2stage
        //
        srnamt = "Rsytrs_aa_2stage";
        infot = 1;
        Rsytrs_aa_2stage("/", 0, 0, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("Rsytrs_aa_2stage", infot, nout, lerr, ok);
        infot = 2;
        Rsytrs_aa_2stage("U", -1, 0, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("Rsytrs_aa_2stage", infot, nout, lerr, ok);
        infot = 3;
        Rsytrs_aa_2stage("U", 0, -1, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("Rsytrs_aa_2stage", infot, nout, lerr, ok);
        infot = 5;
        Rsytrs_aa_2stage("U", 2, 1, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("Rsytrs_aa_2stage", infot, nout, lerr, ok);
        infot = 7;
        Rsytrs_aa_2stage("U", 2, 1, a, 2, a, 1, ip, ip, b, 1, info);
        chkxer("Rsytrs_aa_2stage", infot, nout, lerr, ok);
        infot = 11;
        Rsytrs_aa_2stage("U", 2, 1, a, 2, a, 8, ip, ip, b, 1, info);
        chkxer("Rsytrs_aa_STAGE", infot, nout, lerr, ok);
    } else if (Mlsamen(2, c2, "SP")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite packed matrix with patrial
        //        (Bunch-Kaufman) pivoting.
        //
        //        Rsptrf
        //
        srnamt = "Rsptrf";
        infot = 1;
        Rsptrf("/", 0, a, ip, info);
        chkxer("Rsptrf", infot, nout, lerr, ok);
        infot = 2;
        Rsptrf("U", -1, a, ip, info);
        chkxer("Rsptrf", infot, nout, lerr, ok);
        //
        //        Rsptri
        //
        srnamt = "Rsptri";
        infot = 1;
        Rsptri("/", 0, a, ip, w, info);
        chkxer("Rsptri", infot, nout, lerr, ok);
        infot = 2;
        Rsptri("U", -1, a, ip, w, info);
        chkxer("Rsptri", infot, nout, lerr, ok);
        //
        //        Rsptrs
        //
        srnamt = "Rsptrs";
        infot = 1;
        Rsptrs("/", 0, 0, a, ip, b, 1, info);
        chkxer("Rsptrs", infot, nout, lerr, ok);
        infot = 2;
        Rsptrs("U", -1, 0, a, ip, b, 1, info);
        chkxer("Rsptrs", infot, nout, lerr, ok);
        infot = 3;
        Rsptrs("U", 0, -1, a, ip, b, 1, info);
        chkxer("Rsptrs", infot, nout, lerr, ok);
        infot = 7;
        Rsptrs("U", 2, 1, a, ip, b, 1, info);
        chkxer("Rsptrs", infot, nout, lerr, ok);
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
        //        Rspcon
        //
        srnamt = "Rspcon";
        infot = 1;
        Rspcon("/", 0, a, ip, anrm, rcond, w, iw, info);
        chkxer("Rspcon", infot, nout, lerr, ok);
        infot = 2;
        Rspcon("U", -1, a, ip, anrm, rcond, w, iw, info);
        chkxer("Rspcon", infot, nout, lerr, ok);
        infot = 5;
        Rspcon("U", 1, a, ip, -1.0, rcond, w, iw, info);
        chkxer("Rspcon", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Rerrsy
    //
}
