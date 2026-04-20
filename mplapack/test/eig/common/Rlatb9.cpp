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

// Derived from LAPACK routine DLATB9.
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
#include <mplapack_eig.h>

void Rlatb9(fem::str_cref path, INTEGER const imat, INTEGER const m, INTEGER const p, INTEGER const n, fem::str_ref type, INTEGER &kla, INTEGER &kua, INTEGER &klb, INTEGER &kub, REAL &anorm, REAL &bnorm, INTEGER &modea, INTEGER &modeb, REAL &cndnma, REAL &cndnmb, fem::str_ref dista, fem::str_ref distb) {
    //
    // Set some constants for use in the subroutine.
    //
    const REAL tenth = 0.1;
    const REAL one = 1.0;
    const REAL shrink = 0.25;
    REAL eps = Rlamch("Precision");
    REAL badc2 = tenth / eps;
    REAL badc1 = sqrt(badc2);
    REAL small = Rlamch("Safe minimum");
    REAL large = one / small;
    small = shrink * (small / eps);
    large = one / small;

    //
    // Set some parameters we don't plan to change.
    //
    type = "N";
    dista = "S";
    distb = "S";
    modea = 3;
    modeb = 4;
    //
    // Set the lower and upper bandwidths.
    //
    if (Mlsamen(3, path.elems(), "GRQ") || Mlsamen(3, path.elems(), "LSE") || Mlsamen(3, path.elems(), "GSV")) {
        //
        // A: M by N, B: P by N
        //
        if (imat == 1) {
            //
            // A: diagonal, B: upper triangular
            //
            kla = 0;
            kua = 0;
            klb = 0;
            kub = max(n - 1, (INTEGER)0);
            //
        } else if (imat == 2) {
            //
            // A: upper triangular, B: upper triangular
            //
            kla = 0;
            kua = max(n - 1, (INTEGER)0);
            klb = 0;
            kub = max(n - 1, (INTEGER)0);
            //
        } else if (imat == 3) {
            //
            // A: lower triangular, B: upper triangular
            //
            kla = max(m - 1, (INTEGER)0);
            kua = 0;
            klb = 0;
            kub = max(n - 1, (INTEGER)0);
            //
        } else {
            //
            // A: general dense, B: general dense
            //
            kla = max(m - 1, (INTEGER)0);
            kua = max(n - 1, (INTEGER)0);
            klb = max(p - 1, (INTEGER)0);
            kub = max(n - 1, (INTEGER)0);
            //
        }
        //
    } else if (Mlsamen(3, path.elems(), "GQR") || Mlsamen(3, path.elems(), "GLM")) {
        //
        // A: N by M, B: N by P
        //
        if (imat == 1) {
            //
            // A: diagonal, B: lower triangular
            //
            kla = 0;
            kua = 0;
            klb = max(n - 1, (INTEGER)0);
            kub = 0;
        } else if (imat == 2) {
            //
            // A: lower triangular, B: diagonal
            //
            kla = max(n - 1, (INTEGER)0);
            kua = 0;
            klb = 0;
            kub = 0;
            //
        } else if (imat == 3) {
            //
            // A: lower triangular, B: upper triangular
            //
            kla = max(n - 1, (INTEGER)0);
            kua = 0;
            klb = 0;
            kub = max(p - 1, (INTEGER)0);
            //
        } else {
            //
            // A: general dense, B: general dense
            //
            kla = max(n - 1, (INTEGER)0);
            kua = max(m - 1, (INTEGER)0);
            klb = max(n - 1, (INTEGER)0);
            kub = max(p - 1, (INTEGER)0);
        }
        //
    }
    //
    // Set the condition number and norm.
    //
    const REAL ten = 10.0;
    cndnma = ten * ten;
    cndnmb = ten;
    if (Mlsamen(3, path.elems(), "GQR") || Mlsamen(3, path.elems(), "GRQ") || Mlsamen(3, path.elems(), "GSV")) {
        if (imat == 5) {
            cndnma = badc1;
            cndnmb = badc1;
        } else if (imat == 6) {
            cndnma = badc2;
            cndnmb = badc2;
        } else if (imat == 7) {
            cndnma = badc1;
            cndnmb = badc2;
        } else if (imat == 8) {
            cndnma = badc2;
            cndnmb = badc1;
        }
    }
    //
    anorm = ten;
    bnorm = ten * ten * ten;
    if (Mlsamen(3, path.elems(), "GQR") || Mlsamen(3, path.elems(), "GRQ")) {
        if (imat == 7) {
            anorm = small;
            bnorm = large;
        } else if (imat == 8) {
            anorm = large;
            bnorm = small;
        }
    }
    //
    if (n <= 1) {
        cndnma = one;
        cndnmb = one;
    }
    //
    // End of Rlatb9
    //
}
