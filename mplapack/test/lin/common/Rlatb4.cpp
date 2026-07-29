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

// Derived from LAPACK routine DLATB4.
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

void Rlatb4(fem::str_cref path, INTEGER const imat, INTEGER const m, INTEGER const n, fem::str_ref type, INTEGER &kl, INTEGER &ku, REAL &anorm, INTEGER &mode, REAL &cndnum, fem::str_ref dist) {
    //
    // Set some constants for use in the subroutine.
    //
    const REAL tenth = 0.1;
    const REAL one = 1.0;
    const REAL shrink = 0.25;
    REAL eps = Rlamch("Precision");
    REAL badc2 = tenth / eps;
#if defined MPLAPACK_BUILD_WITH_DD || defined MPLAPACK_BUILD_WITH_BINARY128 || defined MPLAPACK_BUILD_WITH_MPFR
    const REAL badc2_cap = 1.0e24;
    badc2 = min(badc2, badc2_cap);
#elif defined MPLAPACK_BUILD_WITH_QD || defined MPLAPACK_BUILD_WITH_GMP
    const REAL badc2_cap = 1.0e30;
    badc2 = min(badc2, badc2_cap);
#endif
    REAL badc1 = sqrt(badc2);
    REAL small = Rlamch("Safe minimum");
    REAL large = one / small;
    //
#if defined MPLAPACK_BUILD_WITH_GMP
    Rlabad(small, large);
#endif
    small = shrink * (small / eps);
    large = one / small;
    //
    fem::str<2> c2 = path(2, 3);
    //
    // Set some parameters we don't plan to change.
    //
    dist = "S";
    mode = 3;
    //
    const REAL two = 2.0;
    INTEGER mat = 0;
    if (Mlsamen(2, c2.elems, "QR") || Mlsamen(2, c2.elems, "LQ") || Mlsamen(2, c2.elems, "QL") || Mlsamen(2, c2.elems, "RQ")) {
        //
        // xQR, xLQ, xQL, xRQ:  Set parameters to generate a general
        // M x N matrix.
        //
        // Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        // Set the lower and upper bandwidths.
        //
        if (imat == 1) {
            kl = 0;
            ku = 0;
        } else if (imat == 2) {
            kl = 0;
            ku = max(n - 1, (INTEGER)0);
        } else if (imat == 3) {
            kl = max(m - 1, (INTEGER)0);
            ku = 0;
        } else {
            kl = max(m - 1, (INTEGER)0);
            ku = max(n - 1, (INTEGER)0);
        }
        //
        // Set the condition number and norm.
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
    } else if (Mlsamen(2, c2.elems, "QK")) {
        //
        // xQK: truncated QR with pivoting.
        // Set parameters to generate a general
        // M x N matrix.
        //
        // Set TYPE, the type of matrix to be generated.  'N' is nonsymmetric.
        //
        type = "N";
        //
        // Set DIST, the type of distribution for the random
        // number generator. 'S' is
        //
        dist = "S";
        //
        // Set the lower and upper bandwidths.
        //
        if (imat == 2) {
            //
            // 2. Random, Diagonal, CNDNUM = 2
            //
            kl = 0;
            ku = 0;
            cndnum = two;
            anorm = one;
            mode = 3;
        } else if (imat == 3) {
            //
            // 3. Random, Upper triangular,  CNDNUM = 2
            //
            kl = 0;
            ku = max(n - 1, (INTEGER)0);
            cndnum = two;
            anorm = one;
            mode = 3;
        } else if (imat == 4) {
            //
            // 4. Random, Lower triangular,  CNDNUM = 2
            //
            kl = max(m - 1, (INTEGER)0);
            ku = 0;
            cndnum = two;
            anorm = one;
            mode = 3;
        } else {
            //
            // 5.-19. Rectangular matrix
            //
            kl = max(m - 1, (INTEGER)0);
            ku = max(n - 1, (INTEGER)0);
            //
            if (imat >= 5 && imat <= 14) {
                //
                // 5.-14. Random, CNDNUM = 2.
                //
                cndnum = two;
                anorm = one;
                mode = 3;
                //
            } else if (imat == 15) {
                //
                // 15. Random, CNDNUM = sqrt(0.1/EPS)
                //
                cndnum = badc1;
                anorm = one;
                mode = 3;
                //
            } else if (imat == 16) {
                //
                // 16. Random, CNDNUM = 0.1/EPS
                //
                cndnum = badc2;
                anorm = one;
                mode = 3;
                //
            } else if (imat == 17) {
                //
                // 17. Random, CNDNUM = 0.1/EPS,
                // one small singular value S(N)=1/CNDNUM
                //
                cndnum = badc2;
                anorm = one;
                mode = 2;
                //
            } else if (imat == 18) {
                //
                // 18. Random, scaled near underflow
                //
                cndnum = two;
                anorm = small;
                mode = 3;
                //
            } else if (imat == 19) {
                //
                // 19. Random, scaled near overflow
                //
                cndnum = two;
                anorm = large;
                mode = 3;
                //
            }
            //
        }
        //
    } else if (Mlsamen(2, c2.elems, "GE")) {
        //
        // xGE:  Set parameters to generate a general M x N matrix.
        //
        // Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        // Set the lower and upper bandwidths.
        //
        if (imat == 1) {
            kl = 0;
            ku = 0;
        } else if (imat == 2) {
            kl = 0;
            ku = max(n - 1, (INTEGER)0);
        } else if (imat == 3) {
            kl = max(m - 1, (INTEGER)0);
            ku = 0;
        } else {
            kl = max(m - 1, (INTEGER)0);
            ku = max(n - 1, (INTEGER)0);
        }
        //
        // Set the condition number and norm.
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
    } else if (Mlsamen(2, c2.elems, "GB")) {
        //
        // xGB:  Set parameters to generate a general banded matrix.
        //
        // Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        // Set the condition number and norm.
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
    } else if (Mlsamen(2, c2.elems, "GT")) {
        //
        // xGT:  Set parameters to generate a general tridiagonal matrix.
        //
        // Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        // Set the lower and upper bandwidths.
        //
        if (imat == 1) {
            kl = 0;
        } else {
            kl = 1;
        }
        ku = kl;
        //
        // Set the condition number and norm.
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
    } else if (Mlsamen(2, c2.elems, "PO") || Mlsamen(2, c2.elems, "PP")) {
        //
        // xPO, xPP: Set parameters to generate a
        // symmetric positive definite matrix.
        //
        // Set TYPE, the type of matrix to be generated.
        //
        type = c2(1, 1);
        //
        // Set the lower and upper bandwidths.
        //
        if (imat == 1) {
            kl = 0;
        } else {
            kl = max(n - 1, (INTEGER)0);
        }
        ku = kl;
        //
        // Set the condition number and norm.
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
    } else if (Mlsamen(2, c2.elems, "SY") || Mlsamen(2, c2.elems, "SP")) {
        //
        // xSY, xSP: Set parameters to generate a
        // symmetric matrix.
        //
        // Set TYPE, the type of matrix to be generated.
        //
        type = c2(1, 1);
        //
        // Set the lower and upper bandwidths.
        //
        if (imat == 1) {
            kl = 0;
        } else {
            kl = max(n - 1, (INTEGER)0);
        }
        ku = kl;
        //
        // Set the condition number and norm.
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
    } else if (Mlsamen(2, c2.elems, "PB")) {
        //
        // xPB:  Set parameters to generate a symmetric band matrix.
        //
        // Set TYPE, the type of matrix to be generated.
        //
        type = "P";
        //
        // Set the norm and condition number.
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
    } else if (Mlsamen(2, c2.elems, "PT")) {
        //
        // xPT:  Set parameters to generate a symmetric positive definite
        // tridiagonal matrix.
        //
        type = "P";
        if (imat == 1) {
            kl = 0;
        } else {
            kl = 1;
        }
        ku = kl;
        //
        // Set the condition number and norm.
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
    } else if (Mlsamen(2, c2.elems, "TR") || Mlsamen(2, c2.elems, "TP")) {
        //
        // xTR, xTP:  Set parameters to generate a triangular matrix
        //
        // Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        // Set the lower and upper bandwidths.
        //
        mat = abs(imat);
        if (mat == 1 || mat == 7) {
            kl = 0;
            ku = 0;
        } else if (imat < 0) {
            kl = max(n - 1, (INTEGER)0);
            ku = 0;
        } else {
            kl = 0;
            ku = max(n - 1, (INTEGER)0);
        }
        //
        // Set the condition number and norm.
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
    } else if (Mlsamen(2, c2.elems, "TB")) {
        //
        // xTB:  Set parameters to generate a triangular band matrix.
        //
        // Set TYPE, the type of matrix to be generated.
        //
        type = "N";
        //
        // Set the norm and condition number.
        //
        mat = abs(imat);
        if (mat == 2 || mat == 8) {
            cndnum = badc1;
        } else if (mat == 3 || mat == 9) {
            cndnum = badc2;
        } else {
            cndnum = two;
        }
        //
        if (mat == 4) {
            anorm = small;
        } else if (mat == 5) {
            anorm = large;
        } else {
            anorm = one;
        }
    }
    if (n <= 1) {
        cndnum = one;
    }
    //
    // End of Rlatb4
    //
}
