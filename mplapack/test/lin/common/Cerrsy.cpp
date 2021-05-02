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
#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;
#include <mplapack_lin.h>
#include <mplapack.h>

void Cerrsy(const char *path, INTEGER const nunit) {
    common_write write(cmn);
    // COMMON infoc
    INTEGER &infot = cmn.infot;
    INTEGER &nout = cmn.nout;
    bool &ok = cmn.ok;
    bool &lerr = cmn.lerr;
    // COMMON srnamc
    char[32] &srnamt = cmn.srnamt;
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
    char[2] c2 = path[(2 - 1) + (3 - 1) * ldpath];
    //
    //     Set the variables to innocuous values.
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
    INTEGER ip[nmax];
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
    REAL r[nmax];
    REAL rcond = 0.0;
    if (Mlsamen(2, c2, "SY")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with patrial
        //        (Bunch-Kaufman) diagonal pivoting method.
        //
        //        Csytrf
        //
        srnamt = "Csytrf";
        infot = 1;
        Csytrf("/", 0, a, 1, ip, w, 1, info);
        chkxer("Csytrf", infot, nout, lerr, ok);
        infot = 2;
        Csytrf("U", -1, a, 1, ip, w, 1, info);
        chkxer("Csytrf", infot, nout, lerr, ok);
        infot = 4;
        Csytrf("U", 2, a, 1, ip, w, 4, info);
        chkxer("Csytrf", infot, nout, lerr, ok);
        infot = 7;
        Csytrf("U", 0, a, 1, ip, w, 0, info);
        chkxer("Csytrf", infot, nout, lerr, ok);
        infot = 7;
        Csytrf("U", 0, a, 1, ip, w, -2, info);
        chkxer("Csytrf", infot, nout, lerr, ok);
        //
        //        Csytf2
        //
        srnamt = "Csytf2";
        infot = 1;
        Csytf2("/", 0, a, 1, ip, info);
        chkxer("Csytf2", infot, nout, lerr, ok);
        infot = 2;
        Csytf2("U", -1, a, 1, ip, info);
        chkxer("Csytf2", infot, nout, lerr, ok);
        infot = 4;
        Csytf2("U", 2, a, 1, ip, info);
        chkxer("Csytf2", infot, nout, lerr, ok);
        //
        //        Csytri
        //
        srnamt = "Csytri";
        infot = 1;
        Csytri("/", 0, a, 1, ip, w, info);
        chkxer("Csytri", infot, nout, lerr, ok);
        infot = 2;
        Csytri("U", -1, a, 1, ip, w, info);
        chkxer("Csytri", infot, nout, lerr, ok);
        infot = 4;
        Csytri("U", 2, a, 1, ip, w, info);
        chkxer("Csytri", infot, nout, lerr, ok);
        //
        //        Csytri2
        //
        srnamt = "Csytri2";
        infot = 1;
        Csytri2("/", 0, a, 1, ip, w, 1, info);
        chkxer("Csytri2", infot, nout, lerr, ok);
        infot = 2;
        Csytri2("U", -1, a, 1, ip, w, 1, info);
        chkxer("Csytri2", infot, nout, lerr, ok);
        infot = 4;
        Csytri2("U", 2, a, 1, ip, w, 1, info);
        chkxer("Csytri2", infot, nout, lerr, ok);
        //
        //        Csytri2x
        //
        srnamt = "Csytri2x";
        infot = 1;
        Csytri2x("/", 0, a, 1, ip, w, 1, info);
        chkxer("Csytri2x", infot, nout, lerr, ok);
        infot = 2;
        Csytri2x("U", -1, a, 1, ip, w, 1, info);
        chkxer("Csytri2x", infot, nout, lerr, ok);
        infot = 4;
        Csytri2x("U", 2, a, 1, ip, w, 1, info);
        chkxer("Csytri2x", infot, nout, lerr, ok);
        //
        //        Csytrs
        //
        srnamt = "Csytrs";
        infot = 1;
        Csytrs("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("Csytrs", infot, nout, lerr, ok);
        infot = 2;
        Csytrs("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("Csytrs", infot, nout, lerr, ok);
        infot = 3;
        Csytrs("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("Csytrs", infot, nout, lerr, ok);
        infot = 5;
        Csytrs("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("Csytrs", infot, nout, lerr, ok);
        infot = 8;
        Csytrs("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("Csytrs", infot, nout, lerr, ok);
        //
        //        Csyrfs
        //
        srnamt = "Csyrfs";
        infot = 1;
        Csyrfs("/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Csyrfs", infot, nout, lerr, ok);
        infot = 2;
        Csyrfs("U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Csyrfs", infot, nout, lerr, ok);
        infot = 3;
        Csyrfs("U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Csyrfs", infot, nout, lerr, ok);
        infot = 5;
        Csyrfs("U", 2, 1, a, 1, af, 2, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("Csyrfs", infot, nout, lerr, ok);
        infot = 7;
        Csyrfs("U", 2, 1, a, 2, af, 1, ip, b, 2, x, 2, r1, r2, w, r, info);
        chkxer("Csyrfs", infot, nout, lerr, ok);
        infot = 10;
        Csyrfs("U", 2, 1, a, 2, af, 2, ip, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("Csyrfs", infot, nout, lerr, ok);
        infot = 12;
        Csyrfs("U", 2, 1, a, 2, af, 2, ip, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("Csyrfs", infot, nout, lerr, ok);
        //
        //        Csycon
        //
        srnamt = "Csycon";
        infot = 1;
        Csycon("/", 0, a, 1, ip, anrm, rcond, w, info);
        chkxer("Csycon", infot, nout, lerr, ok);
        infot = 2;
        Csycon("U", -1, a, 1, ip, anrm, rcond, w, info);
        chkxer("Csycon", infot, nout, lerr, ok);
        infot = 4;
        Csycon("U", 2, a, 1, ip, anrm, rcond, w, info);
        chkxer("Csycon", infot, nout, lerr, ok);
        infot = 6;
        Csycon("U", 1, a, 1, ip, -anrm, rcond, w, info);
        chkxer("Csycon", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SR")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) diagonal pivoting method.
        //
        //        Csytrf_rook
        //
        srnamt = "Csytrf_rook";
        infot = 1;
        Csytrf_rook("/", 0, a, 1, ip, w, 1, info);
        chkxer("Csytrf_rook", infot, nout, lerr, ok);
        infot = 2;
        Csytrf_rook("U", -1, a, 1, ip, w, 1, info);
        chkxer("Csytrf_rook", infot, nout, lerr, ok);
        infot = 4;
        Csytrf_rook("U", 2, a, 1, ip, w, 4, info);
        chkxer("Csytrf_rook", infot, nout, lerr, ok);
        infot = 7;
        Csytrf_rook("U", 0, a, 1, ip, w, 0, info);
        chkxer("Csytrf_rook", infot, nout, lerr, ok);
        infot = 7;
        Csytrf_rook("U", 0, a, 1, ip, w, -2, info);
        chkxer("Csytrf_rook", infot, nout, lerr, ok);
        //
        //        Csytf2_rook
        //
        srnamt = "Csytf2_rook";
        infot = 1;
        Csytf2_rook("/", 0, a, 1, ip, info);
        chkxer("Csytf2_rook", infot, nout, lerr, ok);
        infot = 2;
        Csytf2_rook("U", -1, a, 1, ip, info);
        chkxer("Csytf2_rook", infot, nout, lerr, ok);
        infot = 4;
        Csytf2_rook("U", 2, a, 1, ip, info);
        chkxer("Csytf2_rook", infot, nout, lerr, ok);
        //
        //        Csytri_rook
        //
        srnamt = "Csytri_rook";
        infot = 1;
        Csytri_rook("/", 0, a, 1, ip, w, info);
        chkxer("Csytri_rook", infot, nout, lerr, ok);
        infot = 2;
        Csytri_rook("U", -1, a, 1, ip, w, info);
        chkxer("Csytri_rook", infot, nout, lerr, ok);
        infot = 4;
        Csytri_rook("U", 2, a, 1, ip, w, info);
        chkxer("Csytri_rook", infot, nout, lerr, ok);
        //
        //        Csytrs_rook
        //
        srnamt = "Csytrs_rook";
        infot = 1;
        Csytrs_rook("/", 0, 0, a, 1, ip, b, 1, info);
        chkxer("Csytrs_rook", infot, nout, lerr, ok);
        infot = 2;
        Csytrs_rook("U", -1, 0, a, 1, ip, b, 1, info);
        chkxer("Csytrs_rook", infot, nout, lerr, ok);
        infot = 3;
        Csytrs_rook("U", 0, -1, a, 1, ip, b, 1, info);
        chkxer("Csytrs_rook", infot, nout, lerr, ok);
        infot = 5;
        Csytrs_rook("U", 2, 1, a, 1, ip, b, 2, info);
        chkxer("Csytrs_rook", infot, nout, lerr, ok);
        infot = 8;
        Csytrs_rook("U", 2, 1, a, 2, ip, b, 1, info);
        chkxer("Csytrs_rook", infot, nout, lerr, ok);
        //
        //        Csycon_rook
        //
        srnamt = "Csycon_rook";
        infot = 1;
        Csycon_rook("/", 0, a, 1, ip, anrm, rcond, w, info);
        chkxer("Csycon_rook", infot, nout, lerr, ok);
        infot = 2;
        Csycon_rook("U", -1, a, 1, ip, anrm, rcond, w, info);
        chkxer("Csycon_rook", infot, nout, lerr, ok);
        infot = 4;
        Csycon_rook("U", 2, a, 1, ip, anrm, rcond, w, info);
        chkxer("Csycon_rook", infot, nout, lerr, ok);
        infot = 6;
        Csycon_rook("U", 1, a, 1, ip, -anrm, rcond, w, info);
        chkxer("Csycon_rook", infot, nout, lerr, ok);
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
        //        Csytrf_rk
        //
        srnamt = "Csytrf_rk";
        infot = 1;
        Csytrf_rk("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("Csytrf_rk", infot, nout, lerr, ok);
        infot = 2;
        Csytrf_rk("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("Csytrf_rk", infot, nout, lerr, ok);
        infot = 4;
        Csytrf_rk("U", 2, a, 1, e, ip, w, 4, info);
        chkxer("Csytrf_rk", infot, nout, lerr, ok);
        infot = 8;
        Csytrf_rk("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("Csytrf_rk", infot, nout, lerr, ok);
        infot = 8;
        Csytrf_rk("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("Csytrf_rk", infot, nout, lerr, ok);
        //
        //        Csytf2_rk
        //
        srnamt = "Csytf2_rk";
        infot = 1;
        Csytf2_rk("/", 0, a, 1, e, ip, info);
        chkxer("Csytf2_rk", infot, nout, lerr, ok);
        infot = 2;
        Csytf2_rk("U", -1, a, 1, e, ip, info);
        chkxer("Csytf2_rk", infot, nout, lerr, ok);
        infot = 4;
        Csytf2_rk("U", 2, a, 1, e, ip, info);
        chkxer("Csytf2_rk", infot, nout, lerr, ok);
        //
        //        Csytri_3
        //
        srnamt = "Csytri_3";
        infot = 1;
        Csytri_3("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("Csytri_3", infot, nout, lerr, ok);
        infot = 2;
        Csytri_3("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("Csytri_3", infot, nout, lerr, ok);
        infot = 4;
        Csytri_3("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("Csytri_3", infot, nout, lerr, ok);
        infot = 8;
        Csytri_3("U", 0, a, 1, e, ip, w, 0, info);
        chkxer("Csytri_3", infot, nout, lerr, ok);
        infot = 8;
        Csytri_3("U", 0, a, 1, e, ip, w, -2, info);
        chkxer("Csytri_3", infot, nout, lerr, ok);
        //
        //        Csytri_3x
        //
        srnamt = "Csytri_3x";
        infot = 1;
        Csytri_3x("/", 0, a, 1, e, ip, w, 1, info);
        chkxer("Csytri_3x", infot, nout, lerr, ok);
        infot = 2;
        Csytri_3x("U", -1, a, 1, e, ip, w, 1, info);
        chkxer("Csytri_3x", infot, nout, lerr, ok);
        infot = 4;
        Csytri_3x("U", 2, a, 1, e, ip, w, 1, info);
        chkxer("Csytri_3x", infot, nout, lerr, ok);
        //
        //        Csytrs_3
        //
        srnamt = "Csytrs_3";
        infot = 1;
        Csytrs_3("/", 0, 0, a, 1, e, ip, b, 1, info);
        chkxer("Csytrs_3", infot, nout, lerr, ok);
        infot = 2;
        Csytrs_3("U", -1, 0, a, 1, e, ip, b, 1, info);
        chkxer("Csytrs_3", infot, nout, lerr, ok);
        infot = 3;
        Csytrs_3("U", 0, -1, a, 1, e, ip, b, 1, info);
        chkxer("Csytrs_3", infot, nout, lerr, ok);
        infot = 5;
        Csytrs_3("U", 2, 1, a, 1, e, ip, b, 2, info);
        chkxer("Csytrs_3", infot, nout, lerr, ok);
        infot = 9;
        Csytrs_3("U", 2, 1, a, 2, e, ip, b, 1, info);
        chkxer("Csytrs_3", infot, nout, lerr, ok);
        //
        //        Csycon_3
        //
        srnamt = "Csycon_3";
        infot = 1;
        Csycon_3("/", 0, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("Csycon_3", infot, nout, lerr, ok);
        infot = 2;
        Csycon_3("U", -1, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("Csycon_3", infot, nout, lerr, ok);
        infot = 4;
        Csycon_3("U", 2, a, 1, e, ip, anrm, rcond, w, info);
        chkxer("Csycon_3", infot, nout, lerr, ok);
        infot = 7;
        Csycon_3("U", 1, a, 1, e, ip, -1.0, rcond, w, info);
        chkxer("Csycon_3", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SP")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite packed matrix with patrial
        //        (Bunch-Kaufman) pivoting.
        //
        //        Csptrf
        //
        srnamt = "Csptrf";
        infot = 1;
        Csptrf("/", 0, a, ip, info);
        chkxer("Csptrf", infot, nout, lerr, ok);
        infot = 2;
        Csptrf("U", -1, a, ip, info);
        chkxer("Csptrf", infot, nout, lerr, ok);
        //
        //        Csptri
        //
        srnamt = "Csptri";
        infot = 1;
        Csptri("/", 0, a, ip, w, info);
        chkxer("Csptri", infot, nout, lerr, ok);
        infot = 2;
        Csptri("U", -1, a, ip, w, info);
        chkxer("Csptri", infot, nout, lerr, ok);
        //
        //        Csptrs
        //
        srnamt = "Csptrs";
        infot = 1;
        Csptrs("/", 0, 0, a, ip, b, 1, info);
        chkxer("Csptrs", infot, nout, lerr, ok);
        infot = 2;
        Csptrs("U", -1, 0, a, ip, b, 1, info);
        chkxer("Csptrs", infot, nout, lerr, ok);
        infot = 3;
        Csptrs("U", 0, -1, a, ip, b, 1, info);
        chkxer("Csptrs", infot, nout, lerr, ok);
        infot = 7;
        Csptrs("U", 2, 1, a, ip, b, 1, info);
        chkxer("Csptrs", infot, nout, lerr, ok);
        //
        //        Csprfs
        //
        srnamt = "Csprfs";
        infot = 1;
        Csprfs("/", 0, 0, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Csprfs", infot, nout, lerr, ok);
        infot = 2;
        Csprfs("U", -1, 0, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Csprfs", infot, nout, lerr, ok);
        infot = 3;
        Csprfs("U", 0, -1, a, af, ip, b, 1, x, 1, r1, r2, w, r, info);
        chkxer("Csprfs", infot, nout, lerr, ok);
        infot = 8;
        Csprfs("U", 2, 1, a, af, ip, b, 1, x, 2, r1, r2, w, r, info);
        chkxer("Csprfs", infot, nout, lerr, ok);
        infot = 10;
        Csprfs("U", 2, 1, a, af, ip, b, 2, x, 1, r1, r2, w, r, info);
        chkxer("Csprfs", infot, nout, lerr, ok);
        //
        //        Cspcon
        //
        srnamt = "Cspcon";
        infot = 1;
        Cspcon("/", 0, a, ip, anrm, rcond, w, info);
        chkxer("Cspcon", infot, nout, lerr, ok);
        infot = 2;
        Cspcon("U", -1, a, ip, anrm, rcond, w, info);
        chkxer("Cspcon", infot, nout, lerr, ok);
        infot = 5;
        Cspcon("U", 1, a, ip, -anrm, rcond, w, info);
        chkxer("Cspcon", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SA")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with Aasen's algorithm.
        //
        //        Csytrf_aa
        //
        srnamt = "Csytrf_aa";
        infot = 1;
        Csytrf_aa("/", 0, a, 1, ip, w, 1, info);
        chkxer("Csytrf_aa", infot, nout, lerr, ok);
        infot = 2;
        Csytrf_aa("U", -1, a, 1, ip, w, 1, info);
        chkxer("Csytrf_aa", infot, nout, lerr, ok);
        infot = 4;
        Csytrf_aa("U", 2, a, 1, ip, w, 4, info);
        chkxer("Csytrf_aa", infot, nout, lerr, ok);
        infot = 7;
        Csytrf_aa("U", 0, a, 1, ip, w, 0, info);
        chkxer("Csytrf_aa", infot, nout, lerr, ok);
        infot = 7;
        Csytrf_aa("U", 0, a, 1, ip, w, -2, info);
        chkxer("Csytrf_aa", infot, nout, lerr, ok);
        //
        //        Csytrs_aa
        //
        srnamt = "Csytrs_aa";
        infot = 1;
        Csytrs_aa("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Csytrs_aa", infot, nout, lerr, ok);
        infot = 2;
        Csytrs_aa("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Csytrs_aa", infot, nout, lerr, ok);
        infot = 3;
        Csytrs_aa("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("Csytrs_aa", infot, nout, lerr, ok);
        infot = 5;
        Csytrs_aa("U", 2, 1, a, 1, ip, b, 2, w, 1, info);
        chkxer("Csytrs_aa", infot, nout, lerr, ok);
        infot = 8;
        Csytrs_aa("U", 2, 1, a, 2, ip, b, 1, w, 1, info);
        chkxer("Csytrs_aa", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "S2")) {
        //
        //        Test error exits of the routines that use factorization
        //        of a symmetric indefinite matrix with Aasen's algorithm.
        //
        //        Csytrf_aa_2stage
        //
        srnamt = "Csytrf_aa_2stage";
        infot = 1;
        Csytrf_aa_2stage("/", 0, a, 1, a, 1, ip, ip, w, 1, info);
        chkxer("Csytrf_aa_2stage", infot, nout, lerr, ok);
        infot = 2;
        Csytrf_aa_2stage("U", -1, a, 1, a, 1, ip, ip, w, 1, info);
        chkxer("Csytrf_aa_2stage", infot, nout, lerr, ok);
        infot = 4;
        Csytrf_aa_2stage("U", 2, a, 1, a, 2, ip, ip, w, 1, info);
        chkxer("Csytrf_aa_2stage", infot, nout, lerr, ok);
        infot = 6;
        Csytrf_aa_2stage("U", 2, a, 2, a, 1, ip, ip, w, 1, info);
        chkxer("Csytrf_aa_2stage", infot, nout, lerr, ok);
        infot = 10;
        Csytrf_aa_2stage("U", 2, a, 2, a, 8, ip, ip, w, 0, info);
        chkxer("Csytrf_aa_2stage", infot, nout, lerr, ok);
        //
        //        CHETRS_AA_2STAGE
        //
        srnamt = "Csytrs_aa_2stage";
        infot = 1;
        Csytrs_aa_2stage("/", 0, 0, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("Csytrs_aa_2stage", infot, nout, lerr, ok);
        infot = 2;
        Csytrs_aa_2stage("U", -1, 0, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("Csytrs_aa_2stage", infot, nout, lerr, ok);
        infot = 3;
        Csytrs_aa_2stage("U", 0, -1, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("Csytrs_aa_2stage", infot, nout, lerr, ok);
        infot = 5;
        Csytrs_aa_2stage("U", 2, 1, a, 1, a, 1, ip, ip, b, 1, info);
        chkxer("Csytrs_aa_2stage", infot, nout, lerr, ok);
        infot = 7;
        Csytrs_aa_2stage("U", 2, 1, a, 2, a, 1, ip, ip, b, 1, info);
        chkxer("Csytrs_aa_2stage", infot, nout, lerr, ok);
        infot = 11;
        Csytrs_aa_2stage("U", 2, 1, a, 2, a, 8, ip, ip, b, 1, info);
        chkxer("Csytrs_aa_STAGE", infot, nout, lerr, ok);
        //
    }
    //
    //     Print a summary line.
    //
    Alaesm(path, ok, nout);
    //
    //     End of Cerrsy
    //
}
