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
#include <mplapack_eig.h>

#include <mplapack_debug.h>

void Clarhs(const char *path, const char *xtype, const char *uplo, const char *trans, INTEGER const m, INTEGER const n, INTEGER const kl, INTEGER const ku, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const ldx, COMPLEX *b, INTEGER const ldb, INTEGER *iseed, INTEGER &info) {
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters.
    //
    info = 0;
    char c1;
    c1 = path[0];
    char c2[2];
    c2[0] = path[1];
    c2[1] = path[2];
    bool tran = Mlsame(trans, "T") || Mlsame(trans, "C");
    bool notran = !tran;
    bool gen = Mlsame(&path[1], "G");
    bool qrs = Mlsame(&path[1], "Q") || Mlsame(&path[2], "Q");
    bool sym = Mlsame(&path[1], "P") || Mlsame(&path[1], "S") || Mlsame(&path[1], "H");
    bool tri = Mlsame(&path[1], "T");
    bool band = Mlsame(&path[2], "B");
    if (!Mlsame(&c1, "Zomplex precision")) {
        info = -1;
    } else if (!(Mlsame(xtype, "N") || Mlsame(xtype, "C"))) {
        info = -2;
    } else if ((sym || tri) && !(Mlsame(uplo, "U") || Mlsame(uplo, "L"))) {
        info = -3;
    } else if ((gen || qrs) && !(tran || Mlsame(trans, "N"))) {
        info = -4;
    } else if (m < 0) {
        info = -5;
    } else if (n < 0) {
        info = -6;
    } else if (band && kl < 0) {
        info = -7;
    } else if (band && ku < 0) {
        info = -8;
    } else if (nrhs < 0) {
        info = -9;
    } else if ((!band && lda < max((INTEGER)1, m)) || (band && (sym || tri) && lda < kl + 1) || (band && gen && lda < kl + ku + 1)) {
        info = -11;
    } else if ((notran && ldx < max((INTEGER)1, n)) || (tran && ldx < max((INTEGER)1, m))) {
        info = -13;
    } else if ((notran && ldb < max((INTEGER)1, m)) || (tran && ldb < max((INTEGER)1, n))) {
        info = -15;
    }
    if (info != 0) {
        Mxerbla("Clarhs", -info);
        return;
    }
    //
    //     Initialize X to NRHS random vectors unless XTYPE = 'C'.
    //
    INTEGER nx = 0;
    INTEGER mb = 0;
    if (tran) {
        nx = m;
        mb = n;
    } else {
        nx = n;
        mb = m;
    }
    INTEGER j = 0;
    if (!Mlsame(xtype, "C")) {
        for (j = 1; j <= nrhs; j = j + 1) {
            Clarnv(2, iseed, n, &x[(j - 1) * ldx]);
        }
    }
    //
    //     Multiply X by op( A ) using an appropriate
    //     matrix multiply routine.
    //
    const COMPLEX one = COMPLEX(1.0, 0.0);
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    char diag;
    if (Mlsamen(2, c2, "GE") || Mlsamen(2, c2, "QR") || Mlsamen(2, c2, "LQ") || Mlsamen(2, c2, "QL") || Mlsamen(2, c2, "RQ")) {
        //
        //        General matrix
        //
        Cgemm(trans, "N", mb, nrhs, nx, one, a, lda, x, ldx, zero, b, ldb);
        //
    } else if (Mlsamen(2, c2, "PO") || Mlsamen(2, c2, "HE")) {
        //
        //        Hermitian matrix, 2-D storage
        //
        Chemm("Left", uplo, n, nrhs, one, a, lda, x, ldx, zero, b, ldb);
        //
    } else if (Mlsamen(2, c2, "SY")) {
        //
        //        Symmetric matrix, 2-D storage
        //
        Csymm("Left", uplo, n, nrhs, one, a, lda, x, ldx, zero, b, ldb);
        //
    } else if (Mlsamen(2, c2, "GB")) {
        //
        //        General matrix, band storage
        //
        for (j = 1; j <= nrhs; j = j + 1) {
            Cgbmv(trans, m, n, kl, ku, one, a, lda, &x[(j - 1) * ldx], 1, zero, &b[(j - 1) * ldb], 1);
        }
        //
    } else if (Mlsamen(2, c2, "PB") || Mlsamen(2, c2, "HB")) {
        //
        //        Hermitian matrix, band storage
        //
        for (j = 1; j <= nrhs; j = j + 1) {
            Chbmv(uplo, n, kl, one, a, lda, &x[(j - 1) * ldx], 1, zero, &b[(j - 1) * ldb], 1);
        }
        //
    } else if (Mlsamen(2, c2, "SB")) {
        //
        //        Symmetric matrix, band storage
        //
        for (j = 1; j <= nrhs; j = j + 1) {
            Csbmv(uplo, n, kl, one, a, lda, &x[(j - 1) * ldx], 1, zero, &b[(j - 1) * ldb], 1);
        }
        //
    } else if (Mlsamen(2, c2, "PP") || Mlsamen(2, c2, "HP")) {
        //
        //        Hermitian matrix, packed storage
        //
        for (j = 1; j <= nrhs; j = j + 1) {
            Chpmv(uplo, n, one, a, &x[(j - 1) * ldx], 1, zero, &b[(j - 1) * ldb], 1);
        }
        //
    } else if (Mlsamen(2, c2, "SP")) {
        //
        //        Symmetric matrix, packed storage
        //
        for (j = 1; j <= nrhs; j = j + 1) {
            Cspmv(uplo, n, one, a, &x[(j - 1) * ldx], 1, zero, &b[(j - 1) * ldb], 1);
        }
        //
    } else if (Mlsamen(2, c2, "TR")) {
        //
        //        Triangular matrix.  Note that for triangular matrices,
        //           KU = 1 => non-unit triangular
        //           KU = 2 => unit triangular
        //
        Clacpy("Full", n, nrhs, x, ldx, b, ldb);
        if (ku == 2) {
            diag = 'U';
        } else {
            diag = 'N';
        }
        Ctrmm("Left", uplo, trans, &diag, n, nrhs, one, a, lda, b, ldb);
        //
    } else if (Mlsamen(2, c2, "TP")) {
        //
        //        Triangular matrix, packed storage
        //
        Clacpy("Full", n, nrhs, x, ldx, b, ldb);
        if (ku == 2) {
            diag = 'U';
        } else {
            diag = 'N';
        }
        for (j = 1; j <= nrhs; j = j + 1) {
            Ctpmv(uplo, trans, &diag, n, a, &b[(j - 1) * ldb], 1);
        }
        //
    } else if (Mlsamen(2, c2, "TB")) {
        //
        //        Triangular matrix, banded storage
        //
        Clacpy("Full", n, nrhs, x, ldx, b, ldb);
        if (ku == 2) {
            diag = 'U';
        } else {
            diag = 'N';
        }
        for (j = 1; j <= nrhs; j = j + 1) {
            Ctbmv(uplo, trans, &diag, n, kl, a, lda, &b[(j - 1) * ldb], 1);
        }
        //
    } else {
        //
        //        If none of the above, set INFO = -1 and return
        //
        info = -1;
        Mxerbla("Clarhs", -info);
    }
    //
    //     End of Clarhs
    //
}
