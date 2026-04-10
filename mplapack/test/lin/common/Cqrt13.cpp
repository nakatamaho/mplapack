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

// Derived from LAPACK routine ZQRT13.
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

void Cqrt13(INTEGER const scale, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, REAL &norma, INTEGER (&iseed)[4]) {
    //
    if (m <= 0 || n <= 0) {
        return;
    }
    //
    // benign matrix
    //
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &a[(j - 1) * lda]);
        if (j <= m) {
            a[(j - 1) + (j - 1) * lda] += COMPLEX(sign(RCasum(m, &a[(j - 1) * lda], 1), a[(j - 1) + (j - 1) * lda].real()));
        }
    }
    //
    // scaled versions
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
        smlnum = smlnum / Rlamch("Epsilon");
        bignum = one / smlnum;
#if defined ___MPLAPACK_BUILD_WITH_BINARY80___ || defined ___MPLAPACK_BUILD_WITH_BINARY128___
        // Historical DLABAD-style range reduction for very wide exponent ranges.
        //
        // Current LAPACK DLABAD is a no-op, but its former logic reduced SMALL and
        // LARGE by square roots when LOG10(LARGE) > 2000. That condition is false
        // for ordinary IEEE binary32/binary64 and true for wide-range formats such
        // as binary80/binary128. Apply the same policy here so that standard
        // backends keep LAPACK's original test scaling while wide-range backends
        // avoid excessively extreme near-underflow/near-overflow matrices.
        if (log10(bignum) > 2000.0) {
            smlnum = sqrt(smlnum);
            bignum = sqrt(bignum);
        }
#endif
        //
        if (scale == 2) {
            //
            // matrix scaled up
            //
            Clascl("General", 0, 0, norma, bignum, m, n, a, lda, info);
        } else if (scale == 3) {
            //
            // matrix scaled down
            //
            Clascl("General", 0, 0, norma, smlnum, m, n, a, lda, info);
        }
    }
    //
    norma = Clange("One-norm", m, n, a, lda, dummy);
    //
    // End of Cqrt13
    //
}
