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

void Rlatb9(const char *path, INTEGER const imat, INTEGER const m, INTEGER const p, INTEGER const n, char *type, INTEGER &kla, INTEGER &kua, INTEGER &klb, INTEGER &kub, REAL &anorm, REAL &bnorm, INTEGER &modea, INTEGER &modeb, REAL &cndnma, REAL &cndnmb, char *dista, char *distb) {

    REAL badc1;
    REAL badc2;
    REAL eps;
    bool first;
    REAL large;
    REAL small;
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
    //
    //     Set some parameters we don't plan to change.
    //
    *type = 'N';
    *dista = 'S';
    *distb = 'S';
    modea = 3;
    modeb = 4;
    //
    //     Set the lower and upper bandwidths.
    //
    if (Mlsamen(3, path, "GRQ") || Mlsamen(3, path, "LSE") || Mlsamen(3, path, "GSV")) {
        //
        //        A: M by N, B: P by N
        //
        if (imat == 1) {
            //
            //           A: diagonal, B: upper triangular
            //
            kla = 0;
            kua = 0;
            klb = 0;
            kub = max(n - 1, 0);
            //
        } else if (imat == 2) {
            //
            //           A: upper triangular, B: upper triangular
            //
            kla = 0;
            kua = max(n - 1, 0);
            klb = 0;
            kub = max(n - 1, 0);
            //
        } else if (imat == 3) {
            //
            //           A: lower triangular, B: upper triangular
            //
            kla = max(m - 1, 0);
            kua = 0;
            klb = 0;
            kub = max(n - 1, 0);
            //
        } else {
            //
            //           A: general dense, B: general dense
            //
            kla = max(m - 1, 0);
            kua = max(n - 1, 0);
            klb = max(p - 1, 0);
            kub = max(n - 1, 0);
            //
        }
        //
    } else if (Mlsamen(3, path, "GQR") || Mlsamen(3, path, "GLM")) {
        //
        //        A: N by M, B: N by P
        //
        if (imat == 1) {
            //
            //           A: diagonal, B: lower triangular
            //
            kla = 0;
            kua = 0;
            klb = max(n - 1, 0);
            kub = 0;
        } else if (imat == 2) {
            //
            //           A: lower triangular, B: diagonal
            //
            kla = max(n - 1, 0);
            kua = 0;
            klb = 0;
            kub = 0;
            //
        } else if (imat == 3) {
            //
            //           A: lower triangular, B: upper triangular
            //
            kla = max(n - 1, 0);
            kua = 0;
            klb = 0;
            kub = max(p - 1, 0);
            //
        } else {
            //
            //           A: general dense, B: general dense
            //
            kla = max(n - 1, 0);
            kua = max(m - 1, 0);
            klb = max(n - 1, 0);
            kub = max(p - 1, 0);
        }
        //
    }
    //
    //     Set the condition number and norm.
    //
    const REAL ten = 1.0e+1;
    cndnma = ten * ten;
    cndnmb = ten;
    if (Mlsamen(3, path, "GQR") || Mlsamen(3, path, "GRQ") || Mlsamen(3, path, "GSV")) {
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
    if (Mlsamen(3, path, "GQR") || Mlsamen(3, path, "GRQ")) {
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
    //     End of Rlatb9
    //
}
