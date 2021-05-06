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

void Rlqt03(INTEGER const m, INTEGER const n, INTEGER const k, REAL *af, REAL *c, REAL *cc, REAL *q, INTEGER const lda, REAL *tau, REAL *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    //
    //     Copy the first k rows of the factorization to the array Q
    //
    const REAL rogue = -1.0e+10;
    Rlaset("Full", n, n, rogue, rogue, q, lda);
    Rlacpy("Upper", k, n - 1, &af[(2 - 1) * ldaf], lda, &q[(2 - 1) * ldq], lda);
    //
    //     Generate the n-by-n matrix Q
    //
    INTEGER info = 0;
    Rorglq(n, n, k, q, lda, tau, work, lwork, info);
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
            Rlarnv(2, iseed, mc, &c[(j - 1) * ldc]);
        }
        cnorm = Rlange("1", mc, nc, c, lda, rwork);
        if (cnorm == 0.0) {
            cnorm = one;
        }
        //
        for (itrans = 1; itrans <= 2; itrans = itrans + 1) {
            if (itrans == 1) {
                trans = 'N';
            } else {
                trans = 'T';
            }
            //
            //           Copy C
            //
            Rlacpy("Full", mc, nc, c, lda, cc, lda);
            //
            //           Apply Q or Q' to C
            //
            Rormlq(&side, &trans, mc, nc, k, af, lda, tau, cc, lda, work, lwork, info);
            //
            //           Form explicit product and subtract
            //
            if (Mlsame(&side, "L")) {
                Rgemm(&trans, "No transpose", mc, nc, mc, -one, q, lda, c, lda, one, cc, lda);
            } else {
                Rgemm("No transpose", &trans, mc, nc, nc, -one, c, lda, q, lda, one, cc, lda);
            }
            //
            //           Compute error in the difference
            //
            resid = Rlange("1", mc, nc, cc, lda, rwork);
            result[((iside - 1) * 2 + itrans) - 1] = resid / castREAL((max((INTEGER)1, n)) * cnorm * eps);
            //
        }
    }
    //
    //     End of Rlqt03
    //
}
