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

void Rlatb4(const char *path, INTEGER const imat, INTEGER const m, INTEGER const n, char *type, INTEGER &kl, INTEGER &ku, REAL &anorm, INTEGER &mode, REAL &cndnum, char *dist) {
    FEM_CMN_SVE(Rlatb4);
    // SAVE
    REAL &badc1 = sve.badc1;
    REAL &badc2 = sve.badc2;
    REAL &eps = sve.eps;
    bool &first = sve.first;
    REAL &large = sve.large;
    REAL &small = sve.small;
    //
    if (is_called_first_time) {
        first = true;
    }
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
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Save statement ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Set some constants for use in the subroutine.
    //
    const REAL tenth = 0.1e+0;
    const REAL one = 1.0;
    const REAL shrink = 0.25e0;
    if (first) {
        first = false;
        eps = Rlamch("Precision");
        badc2 = tenth / eps;
        badc1 = sqrt(badc2);
        small = Rlamch("Safe minimum");
        large = one / small;
        //
        //        If it looks like we're on a Cray, take the square root of
        //        SMALL and LARGE to avoid overflow and underflow problems.
        //
        Rlabad(small, large);
        small = shrink * (small / eps);
        large = one / small;
    }
    //
    char c2[2] = path[(2 - 1) + (3 - 1) * ldpath];
    //
    //     Set some parameters we don't plan to change.
    //
    dist = "S";
    mode = 3;
    //
    const REAL two = 2.0e+0;
    INTEGER mat = 0;
    if (Mlsamen(2, c2, "QR") || Mlsamen(2, c2, "LQ") || Mlsamen(2, c2, "QL") || Mlsamen(2, c2, "RQ")) {
        //
        //        xQR, xLQ, xQL, xRQ:  Set parameters to generate a general
        //                             M x N matrix.
        //
        //        Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        //        Set the lower and upper bandwidths.
        //
        if (imat == 1) {
            kl = 0;
            ku = 0;
        } else if (imat == 2) {
            kl = 0;
            ku = max(n - 1, 0);
        } else if (imat == 3) {
            kl = max(m - 1, 0);
            ku = 0;
        } else {
            kl = max(m - 1, 0);
            ku = max(n - 1, 0);
        }
        //
        //        Set the condition number and norm.
        //
        if (imat == 5) {
            cndnum = badc1;
        } else if (imat == 6) {
            cndnum = badc2;
        } else {
            cndnum = two;
        }
        //
        if (imat == 7) {
            anorm = small;
        } else if (imat == 8) {
            anorm = large;
        } else {
            anorm = one;
        }
        //
    } else if (Mlsamen(2, c2, "GE")) {
        //
        //        xGE:  Set parameters to generate a general M x N matrix.
        //
        //        Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        //        Set the lower and upper bandwidths.
        //
        if (imat == 1) {
            kl = 0;
            ku = 0;
        } else if (imat == 2) {
            kl = 0;
            ku = max(n - 1, 0);
        } else if (imat == 3) {
            kl = max(m - 1, 0);
            ku = 0;
        } else {
            kl = max(m - 1, 0);
            ku = max(n - 1, 0);
        }
        //
        //        Set the condition number and norm.
        //
        if (imat == 8) {
            cndnum = badc1;
        } else if (imat == 9) {
            cndnum = badc2;
        } else {
            cndnum = two;
        }
        //
        if (imat == 10) {
            anorm = small;
        } else if (imat == 11) {
            anorm = large;
        } else {
            anorm = one;
        }
        //
    } else if (Mlsamen(2, c2, "GB")) {
        //
        //        xGB:  Set parameters to generate a general banded matrix.
        //
        //        Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        //        Set the condition number and norm.
        //
        if (imat == 5) {
            cndnum = badc1;
        } else if (imat == 6) {
            cndnum = tenth * badc2;
        } else {
            cndnum = two;
        }
        //
        if (imat == 7) {
            anorm = small;
        } else if (imat == 8) {
            anorm = large;
        } else {
            anorm = one;
        }
        //
    } else if (Mlsamen(2, c2, "GT")) {
        //
        //        xGT:  Set parameters to generate a general tridiagonal matrix.
        //
        //        Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        //        Set the lower and upper bandwidths.
        //
        if (imat == 1) {
            kl = 0;
        } else {
            kl = 1;
        }
        ku = kl;
        //
        //        Set the condition number and norm.
        //
        if (imat == 3) {
            cndnum = badc1;
        } else if (imat == 4) {
            cndnum = badc2;
        } else {
            cndnum = two;
        }
        //
        if (imat == 5 || imat == 11) {
            anorm = small;
        } else if (imat == 6 || imat == 12) {
            anorm = large;
        } else {
            anorm = one;
        }
        //
    } else if (Mlsamen(2, c2, "PO") || Mlsamen(2, c2, "PP")) {
        //
        //        xPO, xPP: Set parameters to generate a
        //        symmetric positive definite matrix.
        //
        //        Set TYPE, the type of matrix to be generated.
        //
        type = c2[(1 - 1)];
        //
        //        Set the lower and upper bandwidths.
        //
        if (imat == 1) {
            kl = 0;
        } else {
            kl = max(n - 1, 0);
        }
        ku = kl;
        //
        //        Set the condition number and norm.
        //
        if (imat == 6) {
            cndnum = badc1;
        } else if (imat == 7) {
            cndnum = badc2;
        } else {
            cndnum = two;
        }
        //
        if (imat == 8) {
            anorm = small;
        } else if (imat == 9) {
            anorm = large;
        } else {
            anorm = one;
        }
        //
    } else if (Mlsamen(2, c2, "SY") || Mlsamen(2, c2, "SP")) {
        //
        //        xSY, xSP: Set parameters to generate a
        //        symmetric matrix.
        //
        //        Set TYPE, the type of matrix to be generated.
        //
        type = c2[(1 - 1)];
        //
        //        Set the lower and upper bandwidths.
        //
        if (imat == 1) {
            kl = 0;
        } else {
            kl = max(n - 1, 0);
        }
        ku = kl;
        //
        //        Set the condition number and norm.
        //
        if (imat == 7) {
            cndnum = badc1;
        } else if (imat == 8) {
            cndnum = badc2;
        } else {
            cndnum = two;
        }
        //
        if (imat == 9) {
            anorm = small;
        } else if (imat == 10) {
            anorm = large;
        } else {
            anorm = one;
        }
        //
    } else if (Mlsamen(2, c2, "PB")) {
        //
        //        xPB:  Set parameters to generate a symmetric band matrix.
        //
        //        Set TYPE, the type of matrix to be generated.
        //
        type = "P";
        //
        //        Set the norm and condition number.
        //
        if (imat == 5) {
            cndnum = badc1;
        } else if (imat == 6) {
            cndnum = badc2;
        } else {
            cndnum = two;
        }
        //
        if (imat == 7) {
            anorm = small;
        } else if (imat == 8) {
            anorm = large;
        } else {
            anorm = one;
        }
        //
    } else if (Mlsamen(2, c2, "PT")) {
        //
        //        xPT:  Set parameters to generate a symmetric positive definite
        //        tridiagonal matrix.
        //
        type = "P";
        if (imat == 1) {
            kl = 0;
        } else {
            kl = 1;
        }
        ku = kl;
        //
        //        Set the condition number and norm.
        //
        if (imat == 3) {
            cndnum = badc1;
        } else if (imat == 4) {
            cndnum = badc2;
        } else {
            cndnum = two;
        }
        //
        if (imat == 5 || imat == 11) {
            anorm = small;
        } else if (imat == 6 || imat == 12) {
            anorm = large;
        } else {
            anorm = one;
        }
        //
    } else if (Mlsamen(2, c2, "TR") || Mlsamen(2, c2, "TP")) {
        //
        //        xTR, xTP:  Set parameters to generate a triangular matrix
        //
        //        Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        //        Set the lower and upper bandwidths.
        //
        mat = abs(imat);
        if (mat == 1 || mat == 7) {
            kl = 0;
            ku = 0;
        } else if (imat < 0) {
            kl = max(n - 1, 0);
            ku = 0;
        } else {
            kl = 0;
            ku = max(n - 1, 0);
        }
        //
        //        Set the condition number and norm.
        //
        if (mat == 3 || mat == 9) {
            cndnum = badc1;
        } else if (mat == 4) {
            cndnum = badc2;
        } else if (mat == 10) {
            cndnum = badc2;
        } else {
            cndnum = two;
        }
        //
        if (mat == 5) {
            anorm = small;
        } else if (mat == 6) {
            anorm = large;
        } else {
            anorm = one;
        }
        //
    } else if (Mlsamen(2, c2, "TB")) {
        //
        //        xTB:  Set parameters to generate a triangular band matrix.
        //
        //        Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        //        Set the norm and condition number.
        //
        if (imat == 2 || imat == 8) {
            cndnum = badc1;
        } else if (imat == 3 || imat == 9) {
            cndnum = badc2;
        } else {
            cndnum = two;
        }
        //
        if (imat == 4) {
            anorm = small;
        } else if (imat == 5) {
            anorm = large;
        } else {
            anorm = one;
        }
    }
    if (n <= 1) {
        cndnum = one;
    }
    //
    //     End of Rlatb4
    //
}
