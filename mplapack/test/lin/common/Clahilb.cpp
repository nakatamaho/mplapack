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

void Clahilb(INTEGER const n, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *x, INTEGER const ldx, COMPLEX *b, INTEGER const ldb, REAL *work, INTEGER &info, const char *path) {
    FEM_CMN_SVE(Clahilb);
    a([lda * n]);
    x([ldx * nrhs]);
    b([ldb * nrhs]);
    work([n]);
    // SAVE
    COMPLEX *d1(sve.d1, [8]);
    COMPLEX *d2(sve.d2, [8]);
    COMPLEX *invd1(sve.invd1, [8]);
    COMPLEX *invd2(sve.invd2, [8]);
    //
    if (is_called_first_time) {
        data((values, cmplx(-1, 0), cmplx(0, 1), cmplx(-1, -1), cmplx(0, -1), cmplx(1, 0), cmplx(-1, 1), cmplx(1, 1), cmplx(1, -1))), d1;
        data((values, cmplx(-1, 0), cmplx(0, -1), cmplx(-1, 1), cmplx(0, 1), cmplx(1, 0), cmplx(-1, -1), cmplx(1, -1), cmplx(1, 1))), d2;
        data((values, cmplx(-1, 0), cmplx(0, -1), cmplx(-.5f, .5f), cmplx(0, 1), cmplx(1, 0), cmplx(-.5f, -.5f), cmplx(.5f, -.5f), cmplx(.5f, .5f))), invd1;
        data((values, cmplx(-1, 0), cmplx(0, 1), cmplx(-.5f, -.5f), cmplx(0, -1), cmplx(1, 0), cmplx(-.5f, .5f), cmplx(.5f, .5f), cmplx(.5f, -.5f))), invd2;
    }
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
    //     ??? complex uses how many bits ???
    //
    //     d's are generated from random permutation of those eight elements.
    //
    //     ..
    //     .. External Functions
    //     ..
    //     .. Executable Statements ..
    char c2[2] = path[(2 - 1) + (3 - 1) * ldpath];
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
        Mxerbla("Clahilb", -info);
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
    //     If we are testing SY routines,
    //        take D1_i = D2_i, else, D1_i = D2_i*
    INTEGER j = 0;
    const INTEGER size_d = 8;
    if (Mlsamen(2, c2, "SY")) {
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = d1[(mod(j - 1) + ((size_d) + 1) - 1) * ldd1] * (castREAL(m) / (i + j - 1)) * d1[(mod(i - 1) + ((size_d) + 1) - 1) * ldd1];
            }
        }
    } else {
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                a[(i - 1) + (j - 1) * lda] = d1[(mod(j - 1) + ((size_d) + 1) - 1) * ldd1] * (castREAL(m) / (i + j - 1)) * d2[(mod(i - 1) + ((size_d) + 1) - 1) * ldd2];
            }
        }
    }
    //
    //     Generate matrix B as simply the first NRHS columns of M * the
    //     identity.
    COMPLEX tmp = castREAL(m);
    Claset("Full", n, nrhs, (0.0, 0.0), tmp, b, ldb);
    //
    //     Generate the true solutions in X.  Because B = the first NRHS
    //     columns of M*I, the true solutions are just the first NRHS columns
    //     of the inverse Hilbert matrix.
    work[1 - 1] = n;
    for (j = 2; j <= n; j = j + 1) {
        work[j - 1] = (((work[(j - 1) - 1] / (j - 1)) * (j - 1 - n)) / (j - 1)) * (n + j - 1);
    }
    //
    //     If we are testing SY routines,
    //           take D1_i = D2_i, else, D1_i = D2_i*
    if (Mlsamen(2, c2, "SY")) {
        for (j = 1; j <= nrhs; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                x[(i - 1) + (j - 1) * ldx] = invd1[(mod(j - 1) + ((size_d) + 1) - 1) * ldinvd1] * ((work[i - 1] * work[j - 1]) / (i + j - 1)) * invd1[(mod(i - 1) + ((size_d) + 1) - 1) * ldinvd1];
            }
        }
    } else {
        for (j = 1; j <= nrhs; j = j + 1) {
            for (i = 1; i <= n; i = i + 1) {
                x[(i - 1) + (j - 1) * ldx] = invd2[(mod(j - 1) + ((size_d) + 1) - 1) * ldinvd2] * ((work[i - 1] * work[j - 1]) / (i + j - 1)) * invd1[(mod(i - 1) + ((size_d) + 1) - 1) * ldinvd1];
            }
        }
    }
}
