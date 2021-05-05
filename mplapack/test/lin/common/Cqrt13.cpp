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

void Cqrt13(INTEGER const scale, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL &norma, INTEGER *iseed) {
    a([lda * star]);
    iseed([4]);
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
    //     .. Local Arrays ..
    //     ..
    //     .. Executable Statements ..
    //
    if (m <= 0 || n <= 0) {
        return;
    }
    //
    //     benign matrix
    //
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &a[(j - 1) * lda]);
        if (j <= m) {
            a[(j - 1) + (j - 1) * lda] += COMPLEX(sign(RCasum(m, &a[(j - 1) * lda], 1), &a[(j - 1) + (j - 1) * lda].real()));
        }
    }
    //
    //     scaled versions
    //
    REAL dummy[1];
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    INTEGER info = 0;
    if (scale != 1) {
        norma = Clange("Max", m, n, a, lda, dummy);
        smlnum = Rlamch("Safe minimum");
        bignum = one / smlnum;
        Rlabad(smlnum, bignum);
        smlnum = smlnum / Rlamch("Epsilon");
        bignum = one / smlnum;
        //
        if (scale == 2) {
            //
            //           matrix scaled up
            //
            Clascl("General", 0, 0, norma, bignum, m, n, a, lda, info);
        } else if (scale == 3) {
            //
            //           matrix scaled down
            //
            Clascl("General", 0, 0, norma, smlnum, m, n, a, lda, info);
        }
    }
    //
    norma = Clange("One-norm", m, n, a, lda, dummy);
    //
    //     End of Cqrt13
    //
}
