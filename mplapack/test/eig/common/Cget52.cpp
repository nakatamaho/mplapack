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

inline REAL abs1(COMPLEX x) { return abs(x.real()) + abs(x.imag()); }

void Cget52(bool const left, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *e, INTEGER const lde, COMPLEX *alpha, COMPLEX *beta, COMPLEX *work, REAL *rwork, REAL *result) {
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
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    COMPLEX x = 0.0;
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
        trans = 'C';
        normab = 'I';
    } else {
        trans = 'N';
        normab = 'O';
    }
    //
    //     Norm of A, B, and E:
    //
    REAL anorm = max({Clange(&normab, n, n, a, lda, rwork), safmin});
    REAL bnorm = max({Clange(&normab, n, n, b, ldb, rwork), safmin});
    REAL enorm = max({Clange("O", n, n, e, lde, rwork), ulp});
    REAL alfmax = safmax / max(one, bnorm);
    REAL betmax = safmax / max(one, anorm);
    //
    //     Compute error matrix.
    //     Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B| |b(i) A| )
    //
    INTEGER jvec = 0;
    COMPLEX alphai = 0.0;
    COMPLEX betai = 0.0;
    REAL abmax = 0.0;
    REAL scale = 0.0;
    COMPLEX acoeff = 0.0;
    COMPLEX bcoeff = 0.0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    for (jvec = 1; jvec <= n; jvec = jvec + 1) {
        alphai = alpha[jvec - 1];
        betai = beta[jvec - 1];
        abmax = max(abs1(alphai), abs1(betai));
        if (abs1(alphai) > alfmax || abs1(betai) > betmax || abmax < one) {
            scale = one / max(abmax, safmin);
            alphai = scale * alphai;
            betai = scale * betai;
        }
        scale = one / max({abs1(alphai) * bnorm, abs1(betai) * anorm, safmin});
        acoeff = scale * betai;
        bcoeff = scale * alphai;
        if (left) {
            acoeff = conj(acoeff);
            bcoeff = conj(bcoeff);
        }
        Cgemv(&trans, n, n, acoeff, a, lda, &e[(jvec - 1) * lde], 1, czero, &work[(n * (jvec - 1) + 1) - 1], 1);
        Cgemv(&trans, n, n, -bcoeff, b, lda, &e[(jvec - 1) * lde], 1, cone, &work[(n * (jvec - 1) + 1) - 1], 1);
    }
    //
    REAL errnrm = Clange("One", n, n, work, n, rwork) / enorm;
    //
    //     Compute RESULT(1)
    //
    result[1 - 1] = errnrm / ulp;
    //
    //     Normalization of E:
    //
    REAL enrmer = zero;
    REAL temp1 = 0.0;
    INTEGER j = 0;
    for (jvec = 1; jvec <= n; jvec = jvec + 1) {
        temp1 = zero;
        for (j = 1; j <= n; j = j + 1) {
            temp1 = max(temp1, abs1(e[(j - 1) + (jvec - 1) * lde]));
        }
        enrmer = max(enrmer, temp1 - one);
    }
    //
    //     Compute RESULT(2) : the normalization error in E.
    //
    result[2 - 1] = enrmer / (castREAL(n) * ulp);
    //
    //     End of Cget52
    //
}
