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

void Clatsy(const char *uplo, INTEGER const n, COMPLEX *x, INTEGER const ldx, INTEGER *iseed) {
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Initialize constants
    //

    REAL seventeen = 8.0;
    REAL eight = 8.0;
    REAL two = 2.0;
    REAL alpha = (1.0 + sqrt(seventeen)) / eight;
    REAL beta = alpha - 1.0 / 1000.0;
    REAL alpha3 = alpha * alpha * alpha;
    //
    //     UPLO = 'U':  Upper triangular storage
    //
    INTEGER j = 0;
    INTEGER i = 0;
    INTEGER n5 = 0;
    COMPLEX a = 0.0;
    COMPLEX b = 0.0;
    const COMPLEX eye = COMPLEX(0.0, 1.0);
    COMPLEX c = 0.0;
    COMPLEX r = 0.0;
    if (Mlsame(uplo, "U")) {
        //
        //        Fill the upper triangle of the matrix with zeros.
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= j; i = i + 1) {
                x[(i - 1) + (j - 1) * ldx] = 0.0;
            }
        }
        n5 = n / 5;
        n5 = n - 5 * n5 + 1;
        //
        for (i = n; i >= n5; i = i - 5) {
            a = alpha3 * Clarnd(5, iseed);
            b = Clarnd(5, iseed) / alpha;
            c = a - two * b * eye;
            r = c / beta;
            x[(i - 1) + (i - 1) * ldx] = a;
            x[((i - 2) - 1) + (i - 1) * ldx] = b;
            x[((i - 2) - 1) + ((i - 1) - 1) * ldx] = r;
            x[((i - 2) - 1) + ((i - 2) - 1) * ldx] = c;
            x[((i - 1) - 1) + ((i - 1) - 1) * ldx] = Clarnd(2, iseed);
            x[((i - 3) - 1) + ((i - 3) - 1) * ldx] = Clarnd(2, iseed);
            x[((i - 4) - 1) + ((i - 4) - 1) * ldx] = Clarnd(2, iseed);
            if (abs(x[((i - 3) - 1) + ((i - 3) - 1) * ldx]) > abs(x[((i - 4) - 1) + ((i - 4) - 1) * ldx])) {
                x[((i - 4) - 1) + ((i - 3) - 1) * ldx] = two * x[((i - 3) - 1) + ((i - 3) - 1) * ldx];
            } else {
                x[((i - 4) - 1) + ((i - 3) - 1) * ldx] = two * x[((i - 4) - 1) + ((i - 4) - 1) * ldx];
            }
        }
        //
        //        Clean-up for N not a multiple of 5.
        //
        i = n5 - 1;
        if (i > 2) {
            a = alpha3 * Clarnd(5, iseed);
            b = Clarnd(5, iseed) / alpha;
            c = a - two * b * eye;
            r = c / beta;
            x[(i - 1) + (i - 1) * ldx] = a;
            x[((i - 2) - 1) + (i - 1) * ldx] = b;
            x[((i - 2) - 1) + ((i - 1) - 1) * ldx] = r;
            x[((i - 2) - 1) + ((i - 2) - 1) * ldx] = c;
            x[((i - 1) - 1) + ((i - 1) - 1) * ldx] = Clarnd(2, iseed);
            i = i - 3;
        }
        if (i > 1) {
            x[(i - 1) + (i - 1) * ldx] = Clarnd(2, iseed);
            x[((i - 1) - 1) + ((i - 1) - 1) * ldx] = Clarnd(2, iseed);
            if (abs(x[(i - 1) + (i - 1) * ldx]) > abs(x[((i - 1) - 1) + ((i - 1) - 1) * ldx])) {
                x[((i - 1) - 1) + (i - 1) * ldx] = two * x[(i - 1) + (i - 1) * ldx];
            } else {
                x[((i - 1) - 1) + (i - 1) * ldx] = two * x[((i - 1) - 1) + ((i - 1) - 1) * ldx];
            }
            i = i - 2;
        } else if (i == 1) {
            x[(i - 1) + (i - 1) * ldx] = Clarnd(2, iseed);
            i = i - 1;
        }
        //
        //     UPLO = 'L':  Lower triangular storage
        //
    } else {
        //
        //        Fill the lower triangle of the matrix with zeros.
        //
        for (j = 1; j <= n; j = j + 1) {
            for (i = j; i <= n; i = i + 1) {
                x[(i - 1) + (j - 1) * ldx] = 0.0;
            }
        }
        n5 = n / 5;
        n5 = n5 * 5;
        //
        for (i = 1; i <= n5; i = i + 5) {
            a = alpha3 * Clarnd(5, iseed);
            b = Clarnd(5, iseed) / alpha;
            c = a - two * b * eye;
            r = c / beta;
            x[(i - 1) + (i - 1) * ldx] = a;
            x[((i + 2) - 1) + (i - 1) * ldx] = b;
            x[((i + 2) - 1) + ((i + 1) - 1) * ldx] = r;
            x[((i + 2) - 1) + ((i + 2) - 1) * ldx] = c;
            x[((i + 1) - 1) + ((i + 1) - 1) * ldx] = Clarnd(2, iseed);
            x[((i + 3) - 1) + ((i + 3) - 1) * ldx] = Clarnd(2, iseed);
            x[((i + 4) - 1) + ((i + 4) - 1) * ldx] = Clarnd(2, iseed);
            if (abs(x[((i + 3) - 1) + ((i + 3) - 1) * ldx]) > abs(x[((i + 4) - 1) + ((i + 4) - 1) * ldx])) {
                x[((i + 4) - 1) + ((i + 3) - 1) * ldx] = two * x[((i + 3) - 1) + ((i + 3) - 1) * ldx];
            } else {
                x[((i + 4) - 1) + ((i + 3) - 1) * ldx] = two * x[((i + 4) - 1) + ((i + 4) - 1) * ldx];
            }
        }
        //
        //        Clean-up for N not a multiple of 5.
        //
        i = n5 + 1;
        if (i < n - 1) {
            a = alpha3 * Clarnd(5, iseed);
            b = Clarnd(5, iseed) / alpha;
            c = a - two * b * eye;
            r = c / beta;
            x[(i - 1) + (i - 1) * ldx] = a;
            x[((i + 2) - 1) + (i - 1) * ldx] = b;
            x[((i + 2) - 1) + ((i + 1) - 1) * ldx] = r;
            x[((i + 2) - 1) + ((i + 2) - 1) * ldx] = c;
            x[((i + 1) - 1) + ((i + 1) - 1) * ldx] = Clarnd(2, iseed);
            i += 3;
        }
        if (i < n) {
            x[(i - 1) + (i - 1) * ldx] = Clarnd(2, iseed);
            x[((i + 1) - 1) + ((i + 1) - 1) * ldx] = Clarnd(2, iseed);
            if (abs(x[(i - 1) + (i - 1) * ldx]) > abs(x[((i + 1) - 1) + ((i + 1) - 1) * ldx])) {
                x[((i + 1) - 1) + (i - 1) * ldx] = two * x[(i - 1) + (i - 1) * ldx];
            } else {
                x[((i + 1) - 1) + (i - 1) * ldx] = two * x[((i + 1) - 1) + ((i + 1) - 1) * ldx];
            }
            i += 2;
        } else if (i == n) {
            x[(i - 1) + (i - 1) * ldx] = Clarnd(2, iseed);
            i++;
        }
    }
    //
    //     End of Clatsy
    //
}
