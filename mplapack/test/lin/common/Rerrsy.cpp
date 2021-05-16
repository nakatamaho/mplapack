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
#include <mplapack_debug.h>

void Rerrsy(const char *path, INTEGER const nunit) {
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
    char c2[2];
    c2[0] = path[1];
    c2[1] = path[2];
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 4;
    INTEGER i = 0;
    REAL a[nmax * nmax];
    REAL af[nmax * nmax];
    INTEGER lda = nmax;
    INTEGER ldaf = nmax;
    REAL b[nmax];
    REAL e[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    REAL w[3 * nmax];
    REAL x[nmax];
    INTEGER ip[nmax];
    INTEGER iw[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = 1.0 / castREAL(i + j);
            af[(i - 1) + (j - 1) * ldaf] = 1.0 / castREAL(i + j);
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
    bool infot = 0;
    if (Mlsamen(2, c2, "SY")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with patrial
        //        (Bunch-Kaufman) pivoting.
        //
        //        Rsytrf
        //
        strncpy(srnamt, "Rsytrf", srnamt_len);
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
        strncpy(srnamt, "Rsytf2", srnamt_len);
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
        strncpy(srnamt, "Rsytri", srnamt_len);
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
        strncpy(srnamt, "Rsytri2", srnamt_len);
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
        strncpy(srnamt, "Rsytri2x", srnamt_len);
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
        strncpy(srnamt, "Rsytrs", srnamt_len);
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
        //        Rsyrfs
        //
        strncpy(srnamt, "Rsyrfs", srnamt_len);
        infot = 1;
        Rsyrfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rsyrfs", infot, nout, lerr, ok);
        infot = 2;
        Rsyrfs("U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rsyrfs", infot, nout, lerr, ok);
        infot = 3;
        Rsyrfs("U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rsyrfs", infot, nout, lerr, ok);
        infot = 5;
        Rsyrfs("U", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("Rsyrfs", infot, nout, lerr, ok);
        infot = 7;
        Rsyrfs("U", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, iw, info);
        chkxer("Rsyrfs", infot, nout, lerr, ok);
        infot = 10;
        Rsyrfs("U", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("Rsyrfs", infot, nout, lerr, ok);
        infot = 12;
        Rsyrfs("U", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("Rsyrfs", infot, nout, lerr, ok);
        //
        //        Rsycon
        //
        strncpy(srnamt, "Rsycon", srnamt_len);
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
        strncpy(srnamt, "Rsytrf_rook", srnamt_len);
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
        strncpy(srnamt, "Rsytf2_rook", srnamt_len);
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
        strncpy(srnamt, "Rsytri_rook", srnamt_len);
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
        strncpy(srnamt, "Rsytrs_rook", srnamt_len);
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
        strncpy(srnamt, "Rsycon_rook", srnamt_len);
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
        strncpy(srnamt, "Rsytrf_rk", srnamt_len);
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
        strncpy(srnamt, "Rsytf2_rk", srnamt_len);
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
        strncpy(srnamt, "Rsytri_3", srnamt_len);
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
        strncpy(srnamt, "Rsytri_3x", srnamt_len);
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
        strncpy(srnamt, "Rsytrs_3", srnamt_len);
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
        strncpy(srnamt, "Rsycon_3", srnamt_len);
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
        strncpy(srnamt, "Rsytrf_aa", srnamt_len);
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
        strncpy(srnamt, "Rsytrs_aa", srnamt_len);
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
        strncpy(srnamt, "Rsytrf_aa_2stage", srnamt_len);
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
        strncpy(srnamt, "Rsytrs_aa_2stage", srnamt_len);
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
        strncpy(srnamt, "Rsptrf", srnamt_len);
        infot = 1;
        Rsptrf("/", 0, a, ip, info);
        chkxer("Rsptrf", infot, nout, lerr, ok);
        infot = 2;
        Rsptrf("U", -1, a, ip, info);
        chkxer("Rsptrf", infot, nout, lerr, ok);
        //
        //        Rsptri
        //
        strncpy(srnamt, "Rsptri", srnamt_len);
        infot = 1;
        Rsptri("/", 0, a, ip, w, info);
        chkxer("Rsptri", infot, nout, lerr, ok);
        infot = 2;
        Rsptri("U", -1, a, ip, w, info);
        chkxer("Rsptri", infot, nout, lerr, ok);
        //
        //        Rsptrs
        //
        strncpy(srnamt, "Rsptrs", srnamt_len);
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
        //        Rsprfs
        //
        strncpy(srnamt, "Rsprfs", srnamt_len);
        infot = 1;
        Rsprfs("/", 0, 0, a, af, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rsprfs", infot, nout, lerr, ok);
        infot = 2;
        Rsprfs("U", -1, 0, a, af, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rsprfs", infot, nout, lerr, ok);
        infot = 3;
        Rsprfs("U", 0, -1, a, af, ip, b, 1, x, 1, r1, r2, w, iw, info);
        chkxer("Rsprfs", infot, nout, lerr, ok);
        infot = 8;
        Rsprfs("U", 2, 1, a, af, ip, b, 1, x, 2, r1, r2, w, iw, info);
        chkxer("Rsprfs", infot, nout, lerr, ok);
        infot = 10;
        Rsprfs("U", 2, 1, a, af, ip, b, 2, x, 1, r1, r2, w, iw, info);
        chkxer("Rsprfs", infot, nout, lerr, ok);
        //
        //        Rspcon
        //
        infot = 1;
        strncpy(srnamt, "Rspcon", srnamt_len);
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
