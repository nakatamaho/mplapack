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

void Rget52(bool const left, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *e, INTEGER const lde, REAL *alphar, REAL *alphai, REAL *beta, REAL *work, REAL *result) {
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
    const REAL zero = 0.0;
    result[1 - 1] = zero;
    result[2 - 1] = zero;
    if (n <= 0) {
        return;
    }
    //
    REAL safmin = Rlamch("Safe minimum");
    const REAL one = 1.0;
    REAL safmax = one / safmin;
    REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
    //
    char trans;
    char normab;
    if (left) {
        trans = 'T';
        normab = 'I';
    } else {
        trans = 'N';
        normab = 'O';
    }
    //
    //     Norm of A, B, and E:
    //
    REAL anorm = max({Rlange(&normab, n, n, a, lda, work), safmin});
    REAL bnorm = max({Rlange(&normab, n, n, b, ldb, work), safmin});
    REAL enorm = max({Rlange("O", n, n, e, lde, work), ulp});
    REAL alfmax = safmax / max(one, bnorm);
    REAL betmax = safmax / max(one, anorm);
    //
    //     Compute error matrix.
    //     Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B| |b(i) A| )
    //
    bool ilcplx = false;
    INTEGER jvec = 0;
    REAL salfr = 0.0;
    REAL salfi = 0.0;
    REAL sbeta = 0.0;
    REAL abmax = 0.0;
    REAL scale = 0.0;
    REAL acoef = 0.0;
    REAL bcoefr = 0.0;
    const REAL ten = 10.0;
    REAL bcoefi = 0.0;
    for (jvec = 1; jvec <= n; jvec = jvec + 1) {
        if (ilcplx) {
            //
            //           2nd Eigenvalue/-vector of pair -- do nothing
            //
            ilcplx = false;
        } else {
            salfr = alphar[jvec - 1];
            salfi = alphai[jvec - 1];
            sbeta = beta[jvec - 1];
            if (salfi == zero) {
                //
                //              Real eigenvalue and -vector
                //
                abmax = max(abs(salfr), abs(sbeta));
                if (abs(salfr) > alfmax || abs(sbeta) > betmax || abmax < one) {
                    scale = one / max(abmax, safmin);
                    salfr = scale * salfr;
                    sbeta = scale * sbeta;
                }
                scale = one / max({abs(salfr) * bnorm, abs(sbeta) * anorm, safmin});
                acoef = scale * sbeta;
                bcoefr = scale * salfr;
                Rgemv(&trans, n, n, acoef, a, lda, &e[(jvec - 1) * lde], 1, zero, &work[(n * (jvec - 1) + 1) - 1], 1);
                Rgemv(&trans, n, n, -bcoefr, b, lda, &e[(jvec - 1) * lde], 1, one, &work[(n * (jvec - 1) + 1) - 1], 1);
            } else {
                //
                //              Complex conjugate pair
                //
                ilcplx = true;
                if (jvec == n) {
                    result[1 - 1] = ten / ulp;
                    return;
                }
                abmax = max(abs(salfr) + abs(salfi), abs(sbeta));
                if (abs(salfr) + abs(salfi) > alfmax || abs(sbeta) > betmax || abmax < one) {
                    scale = one / max(abmax, safmin);
                    salfr = scale * salfr;
                    salfi = scale * salfi;
                    sbeta = scale * sbeta;
                }
                scale = one / max({(abs(salfr) + abs(salfi)) * bnorm, abs(sbeta) * anorm, safmin});
                acoef = scale * sbeta;
                bcoefr = scale * salfr;
                bcoefi = scale * salfi;
                if (left) {
                    bcoefi = -bcoefi;
                }
                //
                Rgemv(&trans, n, n, acoef, a, lda, &e[(jvec - 1) * lde], 1, zero, &work[(n * (jvec - 1) + 1) - 1], 1);
                Rgemv(&trans, n, n, -bcoefr, b, lda, &e[(jvec - 1) * lde], 1, one, &work[(n * (jvec - 1) + 1) - 1], 1);
                Rgemv(&trans, n, n, bcoefi, b, lda, &e[((jvec + 1) - 1) * lde], 1, one, &work[(n * (jvec - 1) + 1) - 1], 1);
                //
                Rgemv(&trans, n, n, acoef, a, lda, &e[((jvec + 1) - 1) * lde], 1, zero, &work[(n * jvec + 1) - 1], 1);
                Rgemv(&trans, n, n, -bcoefi, b, lda, &e[(jvec - 1) * lde], 1, one, &work[(n * jvec + 1) - 1], 1);
                Rgemv(&trans, n, n, -bcoefr, b, lda, &e[((jvec + 1) - 1) * lde], 1, one, &work[(n * jvec + 1) - 1], 1);
            }
        }
    }
    //
    REAL errnrm = Rlange("One", n, n, work, n, &work[(pow2(n) + 1) - 1]) / enorm;
    //
    //     Compute RESULT(1)
    //
    result[1 - 1] = errnrm / ulp;
    //
    //     Normalization of E:
    //
    REAL enrmer = zero;
    ilcplx = false;
    REAL temp1 = 0.0;
    INTEGER j = 0;
    for (jvec = 1; jvec <= n; jvec = jvec + 1) {
        if (ilcplx) {
            ilcplx = false;
        } else {
            temp1 = zero;
            if (alphai[jvec - 1] == zero) {
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max(temp1, abs(e[(j - 1) + (jvec - 1) * lde]));
                }
                enrmer = max(enrmer, temp1 - one);
            } else {
                ilcplx = true;
                for (j = 1; j <= n; j = j + 1) {
                    temp1 = max(temp1, abs(e[(j - 1) + (jvec - 1) * lde]) + abs(e[(j - 1) + ((jvec + 1) - 1) * lde]));
                }
                enrmer = max(enrmer, temp1 - one);
            }
        }
    }
    //
    //     Compute RESULT(2) : the normalization error in E.
    //
    result[2 - 1] = enrmer / (castREAL(n) * ulp);
    //
    //     End of Rget52
    //
}
