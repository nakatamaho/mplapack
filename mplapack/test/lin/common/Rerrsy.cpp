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

// Derived from LAPACK routine DERRSY.
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

void Rerrsy(fem::str_cref path, INTEGER const nunit) {
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
    REAL e[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    REAL w[3 * nmax];
    REAL x[nmax];
    INTEGER ip[nmax];
    INTEGER iw[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * nmax] = 1.0 / castREAL(i + j);
            af[(i - 1) + (j - 1) * nmax] = 1.0 / castREAL(i + j);
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
    if (Mlsamen(2, c2.elems, "SY")) {
        //
        // Test error exits of the routines that use factorization
        // of a symmetric indefinite matrix with patrial
        // (Bunch-Kaufman) pivoting.
        //
        // Rsytrf
        //
        srnamt = "DSYTRF";
        infot = 1;
        Rsytrf("/", 0, a, 1, ip, w, 1, info);
        Chkxer("DSYTRF", infot, nout, lerr, ok);
        infot = 2;
        Rsytrf("U", -1, a, 1, ip, w, 1, info);
        Chkxer("DSYTRF", infot, nout, lerr, ok);
        infot = 4;
        Rsytrf("U", 2, a, 1, ip, w, 4, info);
        Chkxer("DSYTRF", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf("U", 0, a, 1, ip, w, 0, info);
        Chkxer("DSYTRF", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf("U", 0, a, 1, ip, w, -2, info);
        Chkxer("DSYTRF", infot, nout, lerr, ok);
        //
        // Rsytf2
        //
        srnamt = "DSYTF2";
        infot = 1;
        Rsytf2("/", 0, a, 1, ip, info);
        Chkxer("DSYTF2", infot, nout, lerr, ok);
        infot = 2;
        Rsytf2("U", -1, a, 1, ip, info);
        Chkxer("DSYTF2", infot, nout, lerr, ok);
        infot = 4;
        Rsytf2("U", 2, a, 1, ip, info);
        Chkxer("DSYTF2", infot, nout, lerr, ok);
        //
        // Rsytri
        //
        srnamt = "DSYTRI";
        infot = 1;
        Rsytri("/", 0, a, 1, ip, w, info);
        Chkxer("DSYTRI", infot, nout, lerr, ok);
        infot = 2;
        Rsytri("U", -1, a, 1, ip, w, info);
        Chkxer("DSYTRI", infot, nout, lerr, ok);
        infot = 4;
        Rsytri("U", 2, a, 1, ip, w, info);
        Chkxer("DSYTRI", infot, nout, lerr, ok);
        //
        // Rsytri2
        //
        srnamt = "DSYTRI2";
        infot = 1;
        Rsytri2("/", 0, a, 1, ip, w, iw[1 - 1], info);
        Chkxer("DSYTRI2", infot, nout, lerr, ok);
        infot = 2;
        Rsytri2("U", -1, a, 1, ip, w, iw[1 - 1], info);
        Chkxer("DSYTRI2", infot, nout, lerr, ok);
        infot = 4;
        Rsytri2("U", 2, a, 1, ip, w, iw[1 - 1], info);
        Chkxer("DSYTRI2", infot, nout, lerr, ok);
        //
        // Rsytri2x
        //
        srnamt = "DSYTRI2X";
        infot = 1;
        Rsytri2x("/", 0, a, 1, ip, w, 1, info);
        Chkxer("DSYTRI2X", infot, nout, lerr, ok);
        infot = 2;
        Rsytri2x("U", -1, a, 1, ip, w, 1, info);
        Chkxer("DSYTRI2X", infot, nout, lerr, ok);
        infot = 4;
        Rsytri2x("U", 2, a, 1, ip, w, 1, info);
        Chkxer("DSYTRI2X", infot, nout, lerr, ok);
        //
        // Rsytrs
        //
        srnamt = "DSYTRS";
        infot = 1;
        Rsytrs("/", 0, 0, a, 1, ip, b, 1, info);
        Chkxer("DSYTRS", infot, nout, lerr, ok);
        infot = 2;
        Rsytrs("U", -1, 0, a, 1, ip, b, 1, info);
        Chkxer("DSYTRS", infot, nout, lerr, ok);
        infot = 3;
        Rsytrs("U", 0, -1, a, 1, ip, b, 1, info);
        Chkxer("DSYTRS", infot, nout, lerr, ok);
        infot = 5;
        Rsytrs("U", 2, 1, a, 1, ip, b, 2, info);
        Chkxer("DSYTRS", infot, nout, lerr, ok);
        infot = 8;
        Rsytrs("U", 2, 1, a, 2, ip, b, 1, info);
        Chkxer("DSYTRS", infot, nout, lerr, ok);
        //
        // Rsyrfs
        //
        srnamt = "DSYRFS";
        infot = 1;
        Rsyrfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DSYRFS", infot, nout, lerr, ok);
        infot = 2;
        Rsyrfs("U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DSYRFS", infot, nout, lerr, ok);
        infot = 3;
        Rsyrfs("U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DSYRFS", infot, nout, lerr, ok);
        infot = 5;
        Rsyrfs("U", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, iw, info);
        Chkxer("DSYRFS", infot, nout, lerr, ok);
        infot = 7;
        Rsyrfs("U", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, iw, info);
        Chkxer("DSYRFS", infot, nout, lerr, ok);
        infot = 10;
        Rsyrfs("U", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, iw, info);
        Chkxer("DSYRFS", infot, nout, lerr, ok);
        infot = 12;
        Rsyrfs("U", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, iw, info);
        Chkxer("DSYRFS", infot, nout, lerr, ok);
        //
        // Rsycon
        //
        srnamt = "DSYCON";
        infot = 1;
        Rsycon("/", 0, a, 1, ip, anrm, rcond, w, iw, info);
        Chkxer("DSYCON", infot, nout, lerr, ok);
        infot = 2;
        Rsycon("U", -1, a, 1, ip, anrm, rcond, w, iw, info);
        Chkxer("DSYCON", infot, nout, lerr, ok);
        infot = 4;
        Rsycon("U", 2, a, 1, ip, anrm, rcond, w, iw, info);
        Chkxer("DSYCON", infot, nout, lerr, ok);
        infot = 6;
        Rsycon("U", 1, a, 1, ip, -1.0, rcond, w, iw, info);
        Chkxer("DSYCON", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SR")) {
        //
        // Test error exits of the routines that use factorization
        // of a symmetric indefinite matrix with rook
        // (bounded Bunch-Kaufman) pivoting.
        //
        // Rsytrf_rook
        //
        srnamt = "DSYTRF_ROOK";
        infot = 1;
        Rsytrf_rook("/", 0, a, 1, ip, w, 1, info);
        Chkxer("DSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Rsytrf_rook("U", -1, a, 1, ip, w, 1, info);
        Chkxer("DSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 4;
        Rsytrf_rook("U", 2, a, 1, ip, w, 4, info);
        Chkxer("DSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf_rook("U", 0, a, 1, ip, w, 0, info);
        Chkxer("DSYTRF_ROOK", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf_rook("U", 0, a, 1, ip, w, -2, info);
        Chkxer("DSYTRF_ROOK", infot, nout, lerr, ok);
        //
        // Rsytf2_rook
        //
        srnamt = "DSYTF2_ROOK";
        infot = 1;
        Rsytf2_rook("/", 0, a, 1, ip, info);
        Chkxer("DSYTF2_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Rsytf2_rook("U", -1, a, 1, ip, info);
        Chkxer("DSYTF2_ROOK", infot, nout, lerr, ok);
        infot = 4;
        Rsytf2_rook("U", 2, a, 1, ip, info);
        Chkxer("DSYTF2_ROOK", infot, nout, lerr, ok);
        //
        // Rsytri_rook
        //
        srnamt = "DSYTRI_ROOK";
        infot = 1;
        Rsytri_rook("/", 0, a, 1, ip, w, info);
        Chkxer("DSYTRI_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Rsytri_rook("U", -1, a, 1, ip, w, info);
        Chkxer("DSYTRI_ROOK", infot, nout, lerr, ok);
        infot = 4;
        Rsytri_rook("U", 2, a, 1, ip, w, info);
        Chkxer("DSYTRI_ROOK", infot, nout, lerr, ok);
        //
        // Rsytrs_rook
        //
        srnamt = "DSYTRS_ROOK";
        infot = 1;
        Rsytrs_rook("/", 0, 0, a, 1, ip, b, 1, info);
        Chkxer("DSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Rsytrs_rook("U", -1, 0, a, 1, ip, b, 1, info);
        Chkxer("DSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 3;
        Rsytrs_rook("U", 0, -1, a, 1, ip, b, 1, info);
        Chkxer("DSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 5;
        Rsytrs_rook("U", 2, 1, a, 1, ip, b, 2, info);
        Chkxer("DSYTRS_ROOK", infot, nout, lerr, ok);
        infot = 8;
        Rsytrs_rook("U", 2, 1, a, 2, ip, b, 1, info);
        Chkxer("DSYTRS_ROOK", infot, nout, lerr, ok);
        //
        // Rsycon_rook
        //
        srnamt = "DSYCON_ROOK";
        infot = 1;
        Rsycon_rook("/", 0, a, 1, ip, anrm, rcond, w, iw, info);
        Chkxer("DSYCON_ROOK", infot, nout, lerr, ok);
        infot = 2;
        Rsycon_rook("U", -1, a, 1, ip, anrm, rcond, w, iw, info);
        Chkxer("DSYCON_ROOK", infot, nout, lerr, ok);
        infot = 4;
        Rsycon_rook("U", 2, a, 1, ip, anrm, rcond, w, iw, info);
        Chkxer("DSYCON_ROOK", infot, nout, lerr, ok);
        infot = 6;
        Rsycon_rook("U", 1, a, 1, ip, -1.0, rcond, w, iw, info);
        Chkxer("DSYCON_ROOK", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SK")) {
        //
        // Test error exits of the routines that use factorization
        // of a symmetric indefinite matrix with rook
        // (bounded Bunch-Kaufman) pivoting with the new storage
        // format for factors L ( or U) and D.
        //
        // L (or U) is stored in A, diagonal of D is stored on the
        // diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        // Rsytrf_rk
        //
        srnamt = "DSYTRF_RK";
        infot = 1;
        Rsytrf_rk("/", 0, a, 1, e, ip, w, 1, info);
        Chkxer("DSYTRF_RK", infot, nout, lerr, ok);
        infot = 2;
        Rsytrf_rk("U", -1, a, 1, e, ip, w, 1, info);
        Chkxer("DSYTRF_RK", infot, nout, lerr, ok);
        infot = 4;
        Rsytrf_rk("U", 2, a, 1, e, ip, w, 1, info);
        Chkxer("DSYTRF_RK", infot, nout, lerr, ok);
        infot = 8;
        Rsytrf_rk("U", 0, a, 1, e, ip, w, 0, info);
        Chkxer("DSYTRF_RK", infot, nout, lerr, ok);
        infot = 8;
        Rsytrf_rk("U", 0, a, 1, e, ip, w, -2, info);
        Chkxer("DSYTRF_RK", infot, nout, lerr, ok);
        //
        // Rsytf2_rk
        //
        srnamt = "DSYTF2_RK";
        infot = 1;
        Rsytf2_rk("/", 0, a, 1, e, ip, info);
        Chkxer("DSYTF2_RK", infot, nout, lerr, ok);
        infot = 2;
        Rsytf2_rk("U", -1, a, 1, e, ip, info);
        Chkxer("DSYTF2_RK", infot, nout, lerr, ok);
        infot = 4;
        Rsytf2_rk("U", 2, a, 1, e, ip, info);
        Chkxer("DSYTF2_RK", infot, nout, lerr, ok);
        //
        // Rsytri_3
        //
        srnamt = "DSYTRI_3";
        infot = 1;
        Rsytri_3("/", 0, a, 1, e, ip, w, 1, info);
        Chkxer("DSYTRI_3", infot, nout, lerr, ok);
        infot = 2;
        Rsytri_3("U", -1, a, 1, e, ip, w, 1, info);
        Chkxer("DSYTRI_3", infot, nout, lerr, ok);
        infot = 4;
        Rsytri_3("U", 2, a, 1, e, ip, w, 1, info);
        Chkxer("DSYTRI_3", infot, nout, lerr, ok);
        infot = 8;
        Rsytri_3("U", 0, a, 1, e, ip, w, 0, info);
        Chkxer("DSYTRI_3", infot, nout, lerr, ok);
        infot = 8;
        Rsytri_3("U", 0, a, 1, e, ip, w, -2, info);
        Chkxer("DSYTRI_3", infot, nout, lerr, ok);
        //
        // Rsytri_3x
        //
        srnamt = "DSYTRI_3X";
        infot = 1;
        Rsytri_3x("/", 0, a, 1, e, ip, w, 1, info);
        Chkxer("DSYTRI_3X", infot, nout, lerr, ok);
        infot = 2;
        Rsytri_3x("U", -1, a, 1, e, ip, w, 1, info);
        Chkxer("DSYTRI_3X", infot, nout, lerr, ok);
        infot = 4;
        Rsytri_3x("U", 2, a, 1, e, ip, w, 1, info);
        Chkxer("DSYTRI_3X", infot, nout, lerr, ok);
        //
        // Rsytrs_3
        //
        srnamt = "DSYTRS_3";
        infot = 1;
        Rsytrs_3("/", 0, 0, a, 1, e, ip, b, 1, info);
        Chkxer("DSYTRS_3", infot, nout, lerr, ok);
        infot = 2;
        Rsytrs_3("U", -1, 0, a, 1, e, ip, b, 1, info);
        Chkxer("DSYTRS_3", infot, nout, lerr, ok);
        infot = 3;
        Rsytrs_3("U", 0, -1, a, 1, e, ip, b, 1, info);
        Chkxer("DSYTRS_3", infot, nout, lerr, ok);
        infot = 5;
        Rsytrs_3("U", 2, 1, a, 1, e, ip, b, 2, info);
        Chkxer("DSYTRS_3", infot, nout, lerr, ok);
        infot = 9;
        Rsytrs_3("U", 2, 1, a, 2, e, ip, b, 1, info);
        Chkxer("DSYTRS_3", infot, nout, lerr, ok);
        //
        // Rsycon_3
        //
        srnamt = "DSYCON_3";
        infot = 1;
        Rsycon_3("/", 0, a, 1, e, ip, anrm, rcond, w, iw, info);
        Chkxer("DSYCON_3", infot, nout, lerr, ok);
        infot = 2;
        Rsycon_3("U", -1, a, 1, e, ip, anrm, rcond, w, iw, info);
        Chkxer("DSYCON_3", infot, nout, lerr, ok);
        infot = 4;
        Rsycon_3("U", 2, a, 1, e, ip, anrm, rcond, w, iw, info);
        Chkxer("DSYCON_3", infot, nout, lerr, ok);
        infot = 7;
        Rsycon_3("U", 1, a, 1, e, ip, -1.0, rcond, w, iw, info);
        Chkxer("DSYCON_3", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "SA")) {
        //
        // Test error exits of the routines that use factorization
        // of a symmetric indefinite matrix with Aasen's algorithm.
        //
        // Rsytrf_aa
        //
        srnamt = "DSYTRF_AA";
        infot = 1;
        Rsytrf_aa("/", 0, a, 1, ip, w, 1, info);
        Chkxer("DSYTRF_AA", infot, nout, lerr, ok);
        infot = 2;
        Rsytrf_aa("U", -1, a, 1, ip, w, 1, info);
        Chkxer("DSYTRF_AA", infot, nout, lerr, ok);
        infot = 4;
        Rsytrf_aa("U", 2, a, 1, ip, w, 4, info);
        Chkxer("DSYTRF_AA", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf_aa("U", 0, a, 1, ip, w, 0, info);
        Chkxer("DSYTRF_AA", infot, nout, lerr, ok);
        infot = 7;
        Rsytrf_aa("U", 0, a, 1, ip, w, -2, info);
        Chkxer("DSYTRF_AA", infot, nout, lerr, ok);
        //
        // Rsytrs_aa
        //
        srnamt = "DSYTRS_AA";
        infot = 1;
        Rsytrs_aa("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("DSYTRS_AA", infot, nout, lerr, ok);
        infot = 2;
        Rsytrs_aa("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        Chkxer("DSYTRS_AA", infot, nout, lerr, ok);
        infot = 3;
        Rsytrs_aa("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        Chkxer("DSYTRS_AA", infot, nout, lerr, ok);
        infot = 5;
        Rsytrs_aa("U", 2, 1, a, 1, ip, b, 2, w, 1, info);
        Chkxer("DSYTRS_AA", infot, nout, lerr, ok);
        infot = 8;
        Rsytrs_aa("U", 2, 1, a, 2, ip, b, 1, w, 1, info);
        Chkxer("DSYTRS_AA", infot, nout, lerr, ok);
        infot = 10;
        Rsytrs_aa("U", 0, 1, a, 2, ip, b, 1, w, 0, info);
        Chkxer("DSYTRS_AA", infot, nout, lerr, ok);
        infot = 10;
        Rsytrs_aa("U", 0, 1, a, 2, ip, b, 1, w, -2, info);
        Chkxer("DSYTRS_AA", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2.elems, "S2")) {
        //
        // Test error exits of the routines that use factorization
        // of a symmetric indefinite matrix with Aasen's algorithm.
        //
        // Rsytrf_aa_2stage
        //
        srnamt = "DSYTRF_AA_2STAGE";
        infot = 1;
        Rsytrf_aa_2stage("/", 0, a, 1, a, 1, ip, ip, w, 1, info);
        Chkxer("DSYTRF_AA_2STAGE", infot, nout, lerr, ok);
        infot = 2;
        Rsytrf_aa_2stage("U", -1, a, 1, a, 1, ip, ip, w, 1, info);
        Chkxer("DSYTRF_AA_2STAGE", infot, nout, lerr, ok);
        infot = 4;
        Rsytrf_aa_2stage("U", 2, a, 1, a, 2, ip, ip, w, 1, info);
        Chkxer("DSYTRF_AA_2STAGE", infot, nout, lerr, ok);
        infot = 6;
        Rsytrf_aa_2stage("U", 2, a, 2, a, 1, ip, ip, w, 1, info);
        Chkxer("DSYTRF_AA_2STAGE", infot, nout, lerr, ok);
        infot = 10;
        Rsytrf_aa_2stage("U", 2, a, 2, a, 8, ip, ip, w, 0, info);
        Chkxer("DSYTRF_AA_2STAGE", infot, nout, lerr, ok);
        //
        // Rsytrs_aa_2stage
        //
        srnamt = "DSYTRS_AA_2STAGE";
        infot = 1;
        Rsytrs_aa_2stage("/", 0, 0, a, 1, a, 1, ip, ip, b, 1, info);
        Chkxer("DSYTRS_AA_2STAGE", infot, nout, lerr, ok);
        infot = 2;
        Rsytrs_aa_2stage("U", -1, 0, a, 1, a, 1, ip, ip, b, 1, info);
        Chkxer("DSYTRS_AA_2STAGE", infot, nout, lerr, ok);
        infot = 3;
        Rsytrs_aa_2stage("U", 0, -1, a, 1, a, 1, ip, ip, b, 1, info);
        Chkxer("DSYTRS_AA_2STAGE", infot, nout, lerr, ok);
        infot = 5;
        Rsytrs_aa_2stage("U", 2, 1, a, 1, a, 1, ip, ip, b, 1, info);
        Chkxer("DSYTRS_AA_2STAGE", infot, nout, lerr, ok);
        infot = 7;
        Rsytrs_aa_2stage("U", 2, 1, a, 2, a, 1, ip, ip, b, 1, info);
        Chkxer("DSYTRS_AA_2STAGE", infot, nout, lerr, ok);
        infot = 11;
        Rsytrs_aa_2stage("U", 2, 1, a, 2, a, 8, ip, ip, b, 1, info);
        Chkxer("DSYTRS_AA_STAGE", infot, nout, lerr, ok);
    } else if (Mlsamen(2, c2.elems, "SP")) {
        //
        // Test error exits of the routines that use factorization
        // of a symmetric indefinite packed matrix with patrial
        // (Bunch-Kaufman) pivoting.
        //
        // Rsptrf
        //
        srnamt = "DSPTRF";
        infot = 1;
        Rsptrf("/", 0, a, ip, info);
        Chkxer("DSPTRF", infot, nout, lerr, ok);
        infot = 2;
        Rsptrf("U", -1, a, ip, info);
        Chkxer("DSPTRF", infot, nout, lerr, ok);
        //
        // Rsptri
        //
        srnamt = "DSPTRI";
        infot = 1;
        Rsptri("/", 0, a, ip, w, info);
        Chkxer("DSPTRI", infot, nout, lerr, ok);
        infot = 2;
        Rsptri("U", -1, a, ip, w, info);
        Chkxer("DSPTRI", infot, nout, lerr, ok);
        //
        // Rsptrs
        //
        srnamt = "DSPTRS";
        infot = 1;
        Rsptrs("/", 0, 0, a, ip, b, 1, info);
        Chkxer("DSPTRS", infot, nout, lerr, ok);
        infot = 2;
        Rsptrs("U", -1, 0, a, ip, b, 1, info);
        Chkxer("DSPTRS", infot, nout, lerr, ok);
        infot = 3;
        Rsptrs("U", 0, -1, a, ip, b, 1, info);
        Chkxer("DSPTRS", infot, nout, lerr, ok);
        infot = 7;
        Rsptrs("U", 2, 1, a, ip, b, 1, info);
        Chkxer("DSPTRS", infot, nout, lerr, ok);
        //
        // Rsprfs
        //
        srnamt = "DSPRFS";
        infot = 1;
        Rsprfs("/", 0, 0, a, af, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DSPRFS", infot, nout, lerr, ok);
        infot = 2;
        Rsprfs("U", -1, 0, a, af, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DSPRFS", infot, nout, lerr, ok);
        infot = 3;
        Rsprfs("U", 0, -1, a, af, ip, b, 1, x, 1, r1, r2, w, iw, info);
        Chkxer("DSPRFS", infot, nout, lerr, ok);
        infot = 8;
        Rsprfs("U", 2, 1, a, af, ip, b, 1, x, 2, r1, r2, w, iw, info);
        Chkxer("DSPRFS", infot, nout, lerr, ok);
        infot = 10;
        Rsprfs("U", 2, 1, a, af, ip, b, 2, x, 1, r1, r2, w, iw, info);
        Chkxer("DSPRFS", infot, nout, lerr, ok);
        //
        // Rspcon
        //
        srnamt = "DSPCON";
        infot = 1;
        Rspcon("/", 0, a, ip, anrm, rcond, w, iw, info);
        Chkxer("DSPCON", infot, nout, lerr, ok);
        infot = 2;
        Rspcon("U", -1, a, ip, anrm, rcond, w, iw, info);
        Chkxer("DSPCON", infot, nout, lerr, ok);
        infot = 5;
        Rspcon("U", 1, a, ip, -1.0, rcond, w, iw, info);
        Chkxer("DSPCON", infot, nout, lerr, ok);
    }
    //
    // Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    // End of Rerrsy
    //
}
