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

void Clatsp(const char *uplo, INTEGER const n, COMPLEX *x, INTEGER *iseed) {
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
    REAL alpha = (1.0 + sqrt(17.e0)) / 8.e0;
    REAL beta = alpha - 1.0 / 1000.0;
    REAL alpha3 = alpha * alpha * alpha;
    //
    //     Fill the matrix with zeros.
    //
    INTEGER j = 0;
    for (j = 1; j <= n * (n + 1) / 2; j = j + 1) {
        x[j - 1] = 0.0;
    }
    //
    //     UPLO = 'U':  Upper triangular storage
    //
    INTEGER n5 = 0;
    INTEGER jj = 0;
    COMPLEX a = 0.0;
    COMPLEX b = 0.0;
    const COMPLEX eye = COMPLEX(0.0, 1.0);
    COMPLEX c = 0.0;
    COMPLEX r = 0.0;
    if (uplo == "U") {
        n5 = n / 5;
        n5 = n - 5 * n5 + 1;
        //
        jj = n * (n + 1) / 2;
        for (j = n; j >= n5; j = j - 5) {
            a = alpha3 * zlarnd(5, iseed);
            b = zlarnd(5, iseed) / alpha;
            c = a - 2.0 * b * eye;
            r = c / beta;
            x[jj - 1] = a;
            x[(jj - 2) - 1] = b;
            jj = jj - j;
            x[jj - 1] = zlarnd(2, iseed);
            x[(jj - 1) - 1] = r;
            jj = jj - (j - 1);
            x[jj - 1] = c;
            jj = jj - (j - 2);
            x[jj - 1] = zlarnd(2, iseed);
            jj = jj - (j - 3);
            x[jj - 1] = zlarnd(2, iseed);
            if (abs(x[(jj + (j - 3)) - 1]) > abs(x[jj - 1])) {
                x[(jj + (j - 4)) - 1] = 2.0 * x[(jj + (j - 3)) - 1];
            } else {
                x[(jj + (j - 4)) - 1] = 2.0 * x[jj - 1];
            }
            jj = jj - (j - 4);
        }
        //
        //        Clean-up for N not a multiple of 5.
        //
        j = n5 - 1;
        if (j > 2) {
            a = alpha3 * zlarnd(5, iseed);
            b = zlarnd(5, iseed) / alpha;
            c = a - 2.0 * b * eye;
            r = c / beta;
            x[jj - 1] = a;
            x[(jj - 2) - 1] = b;
            jj = jj - j;
            x[jj - 1] = zlarnd(2, iseed);
            x[(jj - 1) - 1] = r;
            jj = jj - (j - 1);
            x[jj - 1] = c;
            jj = jj - (j - 2);
            j = j - 3;
        }
        if (j > 1) {
            x[jj - 1] = zlarnd(2, iseed);
            x[(jj - j) - 1] = zlarnd(2, iseed);
            if (abs(x[jj - 1]) > abs(x[(jj - j) - 1])) {
                x[(jj - 1) - 1] = 2.0 * x[jj - 1];
            } else {
                x[(jj - 1) - 1] = 2.0 * x[(jj - j) - 1];
            }
            jj = jj - j - (j - 1);
            j = j - 2;
        } else if (j == 1) {
            x[jj - 1] = zlarnd(2, iseed);
            j = j - 1;
        }
        //
        //     UPLO = 'L':  Lower triangular storage
        //
    } else {
        n5 = n / 5;
        n5 = n5 * 5;
        //
        jj = 1;
        for (j = 1; j <= n5; j = j + 5) {
            a = alpha3 * zlarnd(5, iseed);
            b = zlarnd(5, iseed) / alpha;
            c = a - 2.0 * b * eye;
            r = c / beta;
            x[jj - 1] = a;
            x[(jj + 2) - 1] = b;
            jj += (n - j + 1);
            x[jj - 1] = zlarnd(2, iseed);
            x[(jj + 1) - 1] = r;
            jj += (n - j);
            x[jj - 1] = c;
            jj += (n - j - 1);
            x[jj - 1] = zlarnd(2, iseed);
            jj += (n - j - 2);
            x[jj - 1] = zlarnd(2, iseed);
            if (abs(x[(jj - (n - j - 2)) - 1]) > abs(x[jj - 1])) {
                x[(jj - (n - j - 2) + 1) - 1] = 2.0 * x[(jj - (n - j - 2)) - 1];
            } else {
                x[(jj - (n - j - 2) + 1) - 1] = 2.0 * x[jj - 1];
            }
            jj += (n - j - 3);
        }
        //
        //        Clean-up for N not a multiple of 5.
        //
        j = n5 + 1;
        if (j < n - 1) {
            a = alpha3 * zlarnd(5, iseed);
            b = zlarnd(5, iseed) / alpha;
            c = a - 2.0 * b * eye;
            r = c / beta;
            x[jj - 1] = a;
            x[(jj + 2) - 1] = b;
            jj += (n - j + 1);
            x[jj - 1] = zlarnd(2, iseed);
            x[(jj + 1) - 1] = r;
            jj += (n - j);
            x[jj - 1] = c;
            jj += (n - j - 1);
            j += 3;
        }
        if (j < n) {
            x[jj - 1] = zlarnd(2, iseed);
            x[(jj + (n - j + 1)) - 1] = zlarnd(2, iseed);
            if (abs(x[jj - 1]) > abs(x[(jj + (n - j + 1)) - 1])) {
                x[(jj + 1) - 1] = 2.0 * x[jj - 1];
            } else {
                x[(jj + 1) - 1] = 2.0 * x[(jj + (n - j + 1)) - 1];
            }
            jj += (n - j + 1) + (n - j);
            j += 2;
        } else if (j == n) {
            x[jj - 1] = zlarnd(2, iseed);
            jj += (n - j + 1);
            j++;
        }
    }
    //
    //     End of Clatsp
    //
}
