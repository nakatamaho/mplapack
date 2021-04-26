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

void Cerrsy(common &cmn, const char *path, INTEGER const nunit) {
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
        ip[j - 1] = j;
    }
    REAL anrm = 1.0;
    ok = true;
    //
    INTEGER info = 0;
    arr_1d<nmax, REAL> r(fill0);
    REAL rcond = 0.0;
    if (Mlsamen2, c2, "SY") {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with patrial
        //        (Bunch-Kaufman) diagonal pivoting method.
        //
        //        ZSYTRF
        //
        srnamt = "ZSYTRF";
        infot = 1;
        zsytrf("/", 0, a, 1, ip, w, 1, info);
        chkxer("ZSYTRF", infot, nout, lerr, ok);
        infot = 2;
        zsytrf("U", -1, a, 1, ip, w, 1, info);
        chkxer("ZSYTRF", infot, nout, lerr, ok);
        infot = 4;
        zsytrf("U", 2, a, 1, ip, w, 4, info);
        chkxer("ZSYTRF", infot, nout, lerr, ok);
        infot = 7;
        zsytrf("U", 0, a, 1, ip, w, 0, info);
        chkxer("ZSYTRF", infot, nout, lerr, ok);
        infot = 7;
        zsytrf("U", 0, a, 1, ip, w, -2, info);
        chkxer("ZSYTRF", infot, nout, lerr, ok);
        //
        //        ZSYTF2
        //
        srnamt = "ZSYTF2";
        infot = 1;
        zsytf2("/", 0, a, 1, ip, info);
        chkxer("ZSYTF2", infot, nout, lerr, ok);
        infot = 2;
        zsytf2("U", -1, a, 1, ip, info);
        chkxer("ZSYTF2", infot, nout, lerr, ok);
        infot = 4;
        zsytf2("U", 2, a, 1, ip, info);
        chkxer("ZSYTF2", infot, nout, lerr, ok);
        //
        //        ZSYTRI
        //
        srnamt = "ZSYTRI";
        infot = 1;
        zsytri("/", 0, a, 1, ip, w, info);
        chkxer("ZSYTRI", infot, nout, lerr, ok);
        infot = 2;
        zsytri("U", -1, a, 1, ip, w, info);
        chkxer("ZSYTRI", infot, nout, lerr, ok);
        infot = 4;
        zsytri("U", 2, a, 1, ip, w, info);
        chkxer("ZSYTRI", infot, nout, lerr, ok);
        //
        //        ZSYTRI2
        //
        srnamt = "ZSYTRI2";
        infot = 1;
        zsytri2("/", 0, a, 1, ip, w, 1, info);
        chkxer("ZSYTRI2", infot, nout, lerr, ok);
        infot = 2;
        zsytri2("U", -1, a, 1, ip, w, 1, info);
        chkxer("ZSYTRI2", infot, nout, lerr, ok);
        infot = 4;
        zsytri2("U", 2, a, 1, ip, w, 1, info);
        chkxer("ZSYTRI2", infot, nout, lerr, ok);
        //
        //        ZSYTRI2X
        //
        srnamt = "ZSYTRI2X";
        infot = 1;
        zsytri2x("/", 0, a, 1, ip, w, 1, info);
        chkxer("ZSYTRI2X", infot, nout, lerr, ok);
        infot = 2;
        zsytri2x("U", -1, a, 1, ip, w, 1, info);
        chkxer("ZSYTRI2X", infot, nout, lerr, ok);
        infot = 4;
        zsytri2x("U", 2, a, 1, ip, w, 1, info);
        chkxer("ZSYTRI2X", infot, nout, lerr, ok);
        //
        //        ZSYTRS
        //
        srnamt = "ZSYTRS";
        infot = 1;
        zsytrs("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("ZSYTRS", infot, nout, lerr, ok);
        infot = 2;
        zsytrs("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("ZSYTRS", infot, nout, lerr, ok);
        infot = 3;
        zsytrs("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("ZSYTRS", infot, nout, lerr, ok);
        infot = 5;
        zsytrs("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("ZSYTRS", infot, nout, lerr, ok);
        infot = 8;
        zsytrs("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("ZSYTRS", infot, nout, lerr, ok);
        //
        //        ZSYRFS
        //
        srnamt = "ZSYRFS";
        infot = 1;
        zsyrfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZSYRFS", infot, nout, lerr, ok);
        infot = 2;
        zsyrfs("U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZSYRFS", infot, nout, lerr, ok);
        infot = 3;
        zsyrfs("U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZSYRFS", infot, nout, lerr, ok);
        infot = 5;
        zsyrfs("U", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("ZSYRFS", infot, nout, lerr, ok);
        infot = 7;
        zsyrfs("U", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("ZSYRFS", infot, nout, lerr, ok);
        infot = 10;
        zsyrfs("U", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("ZSYRFS", infot, nout, lerr, ok);
        infot = 12;
        zsyrfs("U", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("ZSYRFS", infot, nout, lerr, ok);
        //
        //        ZSYCON
        //
        srnamt = "ZSYCON";
        infot = 1;
        zsycon("/", 0, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZSYCON", infot, nout, lerr, ok);
        infot = 2;
        zsycon("U", -1, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZSYCON", infot, nout, lerr, ok);
        infot = 4;
        zsycon("U", 2, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZSYCON", infot, nout, lerr, ok);
        infot = 6;
        zsycon("U", 1, a, 1, ip, -anrm, rcond, w, info);
        chkxer("ZSYCON", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SR") {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) diagonal pivoting method.
        //
        //        ZSYTRF_ROOK
        //
        srnamt = "ZSYTRF_ROOK";
        infot = 1;
        zsytrf_rook("/", 0, a, 1, ip, w, 1, info);
        chkxer("ZSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zsytrf_rook("U", -1, a, 1, ip, w, 1, info);
        chkxer("ZSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 4;
        zsytrf_rook("U", 2, a, 1, ip, w, 4, info);
        chkxer("ZSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 7;
        zsytrf_rook("U", 0, a, 1, ip, w, 0, info);
        chkxer("ZSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 7;
        zsytrf_rook("U", 0, a, 1, ip, w, -2, info);
        chkxer("ZSYTRF_ROOK", infot, nout, lerr, ok);
        //
        //        ZSYTF2_ROOK
        //
        srnamt = "ZSYTF2_ROOK";
        infot = 1;
        zsytf2_rook("/", 0, a, 1, ip, info);
        chkxer("ZSYTF2_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zsytf2_rook("U", -1, a, 1, ip, info);
        chkxer("ZSYTF2_ROOK", infot, nout, lerr, ok);
        infot = 4;
        zsytf2_rook("U", 2, a, 1, ip, info);
        chkxer("ZSYTF2_ROOK", infot, nout, lerr, ok);
        //
        //        ZSYTRI_ROOK
        //
        srnamt = "ZSYTRI_ROOK";
        infot = 1;
        zsytri_rook("/", 0, a, 1, ip, w, info);
        chkxer("ZSYTRI_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zsytri_rook("U", -1, a, 1, ip, w, info);
        chkxer("ZSYTRI_ROOK", infot, nout, lerr, ok);
        infot = 4;
        zsytri_rook("U", 2, a, 1, ip, w, info);
        chkxer("ZSYTRI_ROOK", infot, nout, lerr, ok);
        //
        //        ZSYTRS_ROOK
        //
        srnamt = "ZSYTRS_ROOK";
        infot = 1;
        zsytrs_rook("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("ZSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zsytrs_rook("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("ZSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 3;
        zsytrs_rook("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("ZSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 5;
        zsytrs_rook("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("ZSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 8;
        zsytrs_rook("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("ZSYTRS_ROOK", infot, nout, lerr, ok);
        //
        //        ZSYCON_ROOK
        //
        srnamt = "ZSYCON_ROOK";
        infot = 1;
        zsycon_rook("/", 0, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZSYCON_ROOK", infot, nout, lerr, ok);
        infot = 2;
        zsycon_rook("U", -1, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZSYCON_ROOK", infot, nout, lerr, ok);
        infot = 4;
        zsycon_rook("U", 2, a, 1, ip, anrm, rcond, w, info);
        chkxer("ZSYCON_ROOK", infot, nout, lerr, ok);
        infot = 6;
        zsycon_rook("U", 1, a, 1, ip, -anrm, rcond, w, info);
        chkxer("ZSYCON_ROOK", infot, nout, lerr, ok);
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
        //        ZSYTRF_RK
        //
        srnamt = "ZSYTRF_RK";
        infot = 1;
        zsytrf_rk("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("ZSYTRF_RK", infot, nout, lerr, ok);
        infot = 2;
        zsytrf_rk("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("ZSYTRF_RK", infot, nout, lerr, ok);
        infot = 4;
        zsytrf_rk("U", 2, a, 1, e, ip, w, 4, info);
        chkxer("ZSYTRF_RK", infot, nout, lerr, ok);
        infot = 8;
        zsytrf_rk("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("ZSYTRF_RK", infot, nout, lerr, ok);
        infot = 8;
        zsytrf_rk("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("ZSYTRF_RK", infot, nout, lerr, ok);
        //
        //        ZSYTF2_RK
        //
        srnamt = "ZSYTF2_RK";
        infot = 1;
        zsytf2_rk("/", 0, a, 1, e, ip, info);
        chkxer("ZSYTF2_RK", infot, nout, lerr, ok);
        infot = 2;
        zsytf2_rk("U", -1, a, 1, e, ip, info);
        chkxer("ZSYTF2_RK", infot, nout, lerr, ok);
        infot = 4;
        zsytf2_rk("U", 2, a, 1, e, ip, info);
        chkxer("ZSYTF2_RK", infot, nout, lerr, ok);
        //
        //        ZSYTRI_3
        //
        srnamt = "ZSYTRI_3";
        infot = 1;
        zsytri_3("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("ZSYTRI_3", infot, nout, lerr, ok);
        infot = 2;
        zsytri_3("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("ZSYTRI_3", infot, nout, lerr, ok);
        infot = 4;
        zsytri_3("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("ZSYTRI_3", infot, nout, lerr, ok);
        infot = 8;
        zsytri_3("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("ZSYTRI_3", infot, nout, lerr, ok);
        infot = 8;
        zsytri_3("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("ZSYTRI_3", infot, nout, lerr, ok);
        //
        //        ZSYTRI_3X
        //
        srnamt = "ZSYTRI_3X";
        infot = 1;
        zsytri_3x("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("ZSYTRI_3X", infot, nout, lerr, ok);
        infot = 2;
        zsytri_3x("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("ZSYTRI_3X", infot, nout, lerr, ok);
        infot = 4;
        zsytri_3x("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("ZSYTRI_3X", infot, nout, lerr, ok);
        //
        //        ZSYTRS_3
        //
        srnamt = "ZSYTRS_3";
        infot = 1;
        zsytrs_3("/", 0, 0, a, 1, e, ip, b, 1, info);
        chkxer("ZSYTRS_3", infot, nout, lerr, ok);
        infot = 2;
        zsytrs_3("U", -1, 0, a, 1, e, ip, b, 1, info);
        chkxer("ZSYTRS_3", infot, nout, lerr, ok);
        infot = 3;
        zsytrs_3("U", 0, -1, a, 1, e, ip, b, 1, info);
        chkxer("ZSYTRS_3", infot, nout, lerr, ok);
        infot = 5;
        zsytrs_3("U", 2, 1, a, 1, e, ip, b, 2, info);
        chkxer("ZSYTRS_3", infot, nout, lerr, ok);
        infot = 9;
        zsytrs_3("U", 2, 1, a, 2, e, ip, b, 1, info);
        chkxer("ZSYTRS_3", infot, nout, lerr, ok);
        //
        //        ZSYCON_3
        //
        srnamt = "ZSYCON_3";
        infot = 1;
        zsycon_3("/", 0, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("ZSYCON_3", infot, nout, lerr, ok);
        infot = 2;
        zsycon_3("U", -1, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("ZSYCON_3", infot, nout, lerr, ok);
        infot = 4;
        zsycon_3("U", 2, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("ZSYCON_3", infot, nout, lerr, ok);
        infot = 7;
        zsycon_3("U", 1, a, 1, e, ip, -1.0, rcond, w, info);
        chkxer("ZSYCON_3", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SP") {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite packed matrix with patrial
        //        (Bunch-Kaufman) pivoting.
        //
        //        ZSPTRF
        //
        srnamt = "ZSPTRF";
        infot = 1;
        zsptrf("/", 0, a, ip, info);
        chkxer("ZSPTRF", infot, nout, lerr, ok);
        infot = 2;
        zsptrf("U", -1, a, ip, info);
        chkxer("ZSPTRF", infot, nout, lerr, ok);
        //
        //        ZSPTRI
        //
        srnamt = "ZSPTRI";
        infot = 1;
        zsptri("/", 0, a, ip, w, info);
        chkxer("ZSPTRI", infot, nout, lerr, ok);
        infot = 2;
        zsptri("U", -1, a, ip, w, info);
        chkxer("ZSPTRI", infot, nout, lerr, ok);
        //
        //        ZSPTRS
        //
        srnamt = "ZSPTRS";
        infot = 1;
        zsptrs("/", 0, 0, a, ip, b, 1, info);
        chkxer("ZSPTRS", infot, nout, lerr, ok);
        infot = 2;
        zsptrs("U", -1, 0, a, ip, b, 1, info);
        chkxer("ZSPTRS", infot, nout, lerr, ok);
        infot = 3;
        zsptrs("U", 0, -1, a, ip, b, 1, info);
        chkxer("ZSPTRS", infot, nout, lerr, ok);
        infot = 7;
        zsptrs("U", 2, 1, a, ip, b, 1, info);
        chkxer("ZSPTRS", infot, nout, lerr, ok);
        //
        //        ZSPRFS
        //
        srnamt = "ZSPRFS";
        infot = 1;
        zsprfs("/", 0, 0, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZSPRFS", infot, nout, lerr, ok);
        infot = 2;
        zsprfs("U", -1, 0, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZSPRFS", infot, nout, lerr, ok);
        infot = 3;
        zsprfs("U", 0, -1, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("ZSPRFS", infot, nout, lerr, ok);
        infot = 8;
        zsprfs("U", 2, 1, a, af, ip, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("ZSPRFS", infot, nout, lerr, ok);
        infot = 10;
        zsprfs("U", 2, 1, a, af, ip, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("ZSPRFS", infot, nout, lerr, ok);
        //
        //        ZSPCON
        //
        srnamt = "ZSPCON";
        infot = 1;
        zspcon("/", 0, a, ip, anrm, rcond, w, info);
        chkxer("ZSPCON", infot, nout, lerr, ok);
        infot = 2;
        zspcon("U", -1, a, ip, anrm, rcond, w, info);
        chkxer("ZSPCON", infot, nout, lerr, ok);
        infot = 5;
        zspcon("U", 1, a, ip, -anrm, rcond, w, info);
        chkxer("ZSPCON", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "SA") {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with Aasen's algorithm.
        //
        //        ZSYTRF_AA
        //
        srnamt = "ZSYTRF_AA";
        infot = 1;
        zsytrf_aa("/", 0, a, 1, ip, w, 1, info);
        chkxer("ZSYTRF_AA", infot, nout, lerr, ok);
        infot = 2;
        zsytrf_aa("U", -1, a, 1, ip, w, 1, info);
        chkxer("ZSYTRF_AA", infot, nout, lerr, ok);
        infot = 4;
        zsytrf_aa("U", 2, a, 1, ip, w, 4, info);
        chkxer("ZSYTRF_AA", infot, nout, lerr, ok);
        infot = 7;
        zsytrf_aa("U", 0, a, 1, ip, w, 0, info);
        chkxer("ZSYTRF_AA", infot, nout, lerr, ok);
        infot = 7;
        zsytrf_aa("U", 0, a, 1, ip, w, -2, info);
        chkxer("ZSYTRF_AA", infot, nout, lerr, ok);
        //
        //        ZSYTRS_AA
        //
        srnamt = "ZSYTRS_AA";
        infot = 1;
        zsytrs_aa("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZSYTRS_AA", infot, nout, lerr, ok);
        infot = 2;
        zsytrs_aa("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZSYTRS_AA", infot, nout, lerr, ok);
        infot = 3;
        zsytrs_aa("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("ZSYTRS_AA", infot, nout, lerr, ok);
        infot = 5;
        zsytrs_aa("U", 2, 1, a, 1, ip, b, 2, w, 1, info);
        chkxer("ZSYTRS_AA", infot, nout, lerr, ok);
        infot = 8;
        zsytrs_aa("U", 2, 1, a, 2, ip, b, 1, w, 1, info);
        chkxer("ZSYTRS_AA", infot, nout, lerr, ok);
        //
    } else if (Mlsamen2, c2, "S2") {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with Aasen's algorithm.
        //
        //        ZSYTRF_AA_2STAGE
        //
        srnamt = "ZSYTRF_AA_2STAGE";
        infot = 1;
        zsytrf_aa_2stage("/", 0, a, 1, a, 1, ip, ip, w, 1, info);
        chkxer("ZSYTRF_AA_2STAGE", infot, nout, lerr, ok);
        infot = 2;
        zsytrf_aa_2stage("U", -1, a, 1, a, 1, ip, ip, w, 1, info);
        chkxer("ZSYTRF_AA_2STAGE", infot, nout, lerr, ok);
        infot = 4;
        zsytrf_aa_2stage("U", 2, a, 1, a, 2, ip, ip, w, 1, info);
        chkxer("ZSYTRF_AA_2STAGE", infot, nout, lerr, ok);
        infot = 6;
        zsytrf_aa_2stage("U", 2, a, 2, a, 1, ip, ip, w, 1, info);
        chkxer("ZSYTRF_AA_2STAGE", infot, nout, lerr, ok);
        infot = 10;
        zsytrf_aa_2stage("U", 2, a, 2, a, 8, ip, ip, w, 0, info);
        chkxer("ZSYTRF_AA_2STAGE", infot, nout, lerr, ok);
        //
        //        CHETRS_AA_2STAGE
        //
        srnamt = "ZSYTRS_AA_2STAGE";
        infot = 1;
        zsytrs_aa_2stage("/", 0, 0, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("ZSYTRS_AA_2STAGE", infot, nout, lerr, ok);
        infot = 2;
        zsytrs_aa_2stage("U", -1, 0, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("ZSYTRS_AA_2STAGE", infot, nout, lerr, ok);
        infot = 3;
        zsytrs_aa_2stage("U", 0, -1, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("ZSYTRS_AA_2STAGE", infot, nout, lerr, ok);
        infot = 5;
        zsytrs_aa_2stage("U", 2, 1, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("ZSYTRS_AA_2STAGE", infot, nout, lerr, ok);
        infot = 7;
        zsytrs_aa_2stage("U", 2, 1, a, 2, a, 1, ip, ip, b, 1, info);
        chkxer("ZSYTRS_AA_2STAGE", infot, nout, lerr, ok);
        infot = 11;
        zsytrs_aa_2stage("U", 2, 1, a, 2, a, 8, ip, ip, b, 1, info);
        chkxer("ZSYTRS_AA_STAGE", infot, nout, lerr, ok);
        //
    }
    //
    //     PrINTEGER a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrsy
    //
}
