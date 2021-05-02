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

void Rlahilb(INTEGER const n, INTEGER const nrhs, REAL *a, INTEGER const lda, REAL *x, INTEGER const ldx, REAL *b, INTEGER const ldb, REAL *work, INTEGER &info) {
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //     .. Local Scalars ..
    //     ..
    //     .. Parameters ..
    //                  exact.
    //                  a small componentwise relative error.
    //
    //     ..
    //     .. External Functions
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments
    //
    info = 0;
    const INTEGER nmax_approx = 11;
    if (n < 0 || n > nmax_approx) {
        info = -1;
    } else if (nrhs < 0) {
        info = -2;
    } else if (lda < n) {
        info = -4;
    } else if (ldx < n) {
        info = -6;
    } else if (ldb < n) {
        info = -8;
    }
    if (info < 0) {
        Mxerbla("Rlahilb", -info);
        return;
    }
    const INTEGER nmax_exact = 6;
    if (n > nmax_exact) {
        info = 1;
    }
    //
    //     Compute M = the LCM of the integers [1, 2*N-1].  The largest
    //     reasonable N is small enough that integers suffice (up to N = 11).
    INTEGER m = 1;
    INTEGER i = 0;
    INTEGER tm = 0;
    INTEGER ti = 0;
    INTEGER r = 0;
    for (i = 2; i <= (2 * n - 1); i = i + 1) {
        tm = m;
        ti = i;
        r = mod(tm, ti);
        while (r != 0) {
            tm = ti;
            ti = r;
            r = mod(tm, ti);
        }
        m = (m / ti) * i;
    }
    //
    //     Generate the scaled Hilbert matrix in A
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        for (i = 1; i <= n; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = m.real() / (i + j - 1);
        }
    }
    //
    //     Generate matrix B as simply the first NRHS columns of M * the
    //     identity.
    Rlaset("Full", n, nrhs, 0.0, m.real(), b, ldb);
    //
    //     Generate the true solutions in X.  Because B = the first NRHS
    //     columns of M*I, the true solutions are just the first NRHS columns
    //     of the inverse Hilbert matrix.
    work[1 - 1] = n;
    for (j = 2; j <= n; j = j + 1) {
        work[j - 1] = (((work[(j - 1) - 1] / (j - 1)) * (j - 1 - n)) / (j - 1)) * (n + j - 1);
    }
    //
    for (j = 1; j <= nrhs; j = j + 1) {
        for (i = 1; i <= n; i = i + 1) {
            x[(i - 1) + (j - 1) * ldx] = (work[i - 1] * work[j - 1]) / (i + j - 1);
        }
    }
    //
}
