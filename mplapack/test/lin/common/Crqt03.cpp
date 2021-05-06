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

void Crqt03(INTEGER const m, INTEGER const n, INTEGER const k, COMPLEX *af, COMPLEX *c, COMPLEX *cc, COMPLEX *q, INTEGER const lda, COMPLEX *tau, COMPLEX *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    //     .. Local Arrays ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER ldaf = lda;
    INTEGER ldc = lda;
    INTEGER ldcc = lda;
    INTEGER ldq = lda;
    INTEGER iseed[] = {1988, 1989, 1990, 1991};
    REAL eps = Rlamch("Epsilon");
    INTEGER minmn = min(m, n);
    //
    //     Quick return if possible
    //
    const REAL zero = 0.0;
    if (minmn == 0) {
        result[1 - 1] = zero;
        result[2 - 1] = zero;
        result[3 - 1] = zero;
        result[4 - 1] = zero;
        return;
    }
    //
    //     Copy the last k rows of the factorization to the array Q
    //
    const COMPLEX rogue = COMPLEX(-1.0e+10, -1.0e+10);
    Claset("Full", n, n, rogue, rogue, q, lda);
    if (k > 0 && n > k) {
        Clacpy("Full", k, n - k, &af[((m - k + 1) - 1)], lda, &q[((n - k + 1) - 1)], lda);
    }
    if (k > 1) {
        Clacpy("Lower", k - 1, k - 1, &af[((m - k + 2) - 1) + ((n - k + 1) - 1) * ldaf], lda, &q[((n - k + 2) - 1) + ((n - k + 1) - 1) * ldq], lda);
    }
    //
    //     Generate the n-by-n matrix Q
    //
    INTEGER info = 0;
    Cungrq(n, n, k, q, lda, &tau[(minmn - k + 1) - 1], work, lwork, info);
    //
    INTEGER iside = 0;
    char side;
    INTEGER mc = 0;
    INTEGER nc = 0;
    INTEGER j = 0;
    REAL cnorm = 0.0;
    const REAL one = 1.0;
    INTEGER itrans = 0;
    char trans;
    REAL resid = 0.0;
    for (iside = 1; iside <= 2; iside = iside + 1) {
        if (iside == 1) {
            side = 'L';
            mc = n;
            nc = m;
        } else {
            side = 'R';
            mc = m;
            nc = n;
        }
        //
        //        Generate MC by NC matrix C
        //
        for (j = 1; j <= nc; j = j + 1) {
            Clarnv(2, iseed, mc, &c[(j - 1) * ldc]);
        }
        cnorm = Clange("1", mc, nc, c, lda, rwork);
        if (cnorm == zero) {
            cnorm = one;
        }
        //
        for (itrans = 1; itrans <= 2; itrans = itrans + 1) {
            if (itrans == 1) {
                trans = 'N';
            } else {
                trans = 'C';
            }
            //
            //           Copy C
            //
            Clacpy("Full", mc, nc, c, lda, cc, lda);
            //
            //           Apply Q or Q' to C
            //
            if (k > 0) {
                Cunmrq(&side, &trans, mc, nc, k, &af[((m - k + 1) - 1)], lda, &tau[(minmn - k + 1) - 1], cc, lda, work, lwork, info);
            }
            //
            //           Form explicit product and subtract
            //
            if (Mlsame(&side, "L")) {
                Cgemm(&trans, "No transpose", mc, nc, mc, COMPLEX(-one), q, lda, c, lda, COMPLEX(one), cc, lda);
            } else {
                Cgemm("No transpose", &trans, mc, nc, nc, COMPLEX(-one), c, lda, q, lda, COMPLEX(one), cc, lda);
            }
            //
            //           Compute error in the difference
            //
            resid = Clange("1", mc, nc, cc, lda, rwork);
            result[((iside - 1) * 2 + itrans) - 1] = resid / castREAL((max((INTEGER)1, n)) * cnorm * eps);
            //
        }
    }
    //
    //     End of Crqt03
    //
}
