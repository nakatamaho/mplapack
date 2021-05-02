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

void Cqrt15(INTEGER const scale, INTEGER const rksel, INTEGER const m, INTEGER const n, INTEGER const nrhs, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL *s, INTEGER &rank, REAL &norma, REAL &normb, INTEGER *iseed, COMPLEX *work, INTEGER const lwork) {
    INTEGER mn = 0;
    REAL smlnum = 0.0;
    const REAL one = 1.0;
    REAL bignum = 0.0;
    REAL eps = 0.0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    REAL temp = 0.0;
    const REAL svmin = 0.1e+0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const REAL two = 2.0e+0;
    INTEGER info = 0;
    arr_1d<1, REAL> dummy;
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    mn = min(m, n);
    if (lwork < max({m + mn, mn * nrhs, 2 * n + m})) {
        Mxerbla("Cqrt15", 16);
        return;
    }
    //
    smlnum = Rlamch("Safe minimum");
    bignum = one / smlnum;
    Rlabad(smlnum, bignum);
    eps = Rlamch("Epsilon");
    smlnum = (smlnum / eps) / eps;
    bignum = one / smlnum;
    //
    //     Determine rank and (unscaled) singular values
    //
    if (rksel == 1) {
        rank = mn;
    } else if (rksel == 2) {
        rank = (3 * mn) / 4;
        for (j = rank + 1; j <= mn; j = j + 1) {
            s[j - 1] = zero;
        }
    } else {
        Mxerbla("Cqrt15", 2);
    }
    //
    if (rank > 0) {
        //
        //        Nontrivial case
        //
        s[1 - 1] = one;
        for (j = 2; j <= rank; j = j + 1) {
        statement_20:
            temp = dlarnd(1, iseed);
            if (temp > svmin) {
                s[j - 1] = abs(temp);
            } else {
                goto statement_20;
            }
        }
        Rlaord("Decreasing", rank, s, 1);
        //
        //        Generate 'rank' columns of a random orthogonal matrix in A
        //
        Clarnv(2, iseed, m, work);
        CRscal(m, one / RCnrm2(m, work, 1), work, 1);
        Claset("Full", m, rank, czero, cone, a, lda);
        Clarf("Left", m, rank, work, 1, COMPLEX(two), a, lda, &work[(m + 1) - 1]);
        //
        //        workspace used: m+mn
        //
        //        Generate consistent rhs in the range space of A
        //
        Clarnv(2, iseed, rank * nrhs, work);
        Cgemm("No transpose", "No transpose", m, nrhs, rank, cone, a, lda, work, rank, czero, b, ldb);
        //
        //        work space used: <= mn *nrhs
        //
        //        generate (unscaled) matrix A
        //
        for (j = 1; j <= rank; j = j + 1) {
            CRscal(m, s[j - 1], &a[(j - 1) * lda], 1);
        }
        if (rank < n) {
            Claset("Full", m, n - rank, czero, czero, &a[((rank + 1) - 1) * lda], lda);
        }
        zlaror("Right", "No initialization", m, n, a, lda, iseed, work, info);
        //
    } else {
        //
        //        work space used 2*n+m
        //
        //        Generate null matrix and rhs
        //
        for (j = 1; j <= mn; j = j + 1) {
            s[j - 1] = zero;
        }
        Claset("Full", m, n, czero, czero, a, lda);
        Claset("Full", m, nrhs, czero, czero, b, ldb);
        //
    }
    //
    //     Scale the matrix
    //
    if (scale != 1) {
        norma = Clange("Max", m, n, a, lda, dummy);
        if (norma != zero) {
            if (scale == 2) {
                //
                //              matrix scaled up
                //
                Clascl("General", 0, 0, norma, bignum, m, n, a, lda, info);
                Rlascl("General", 0, 0, norma, bignum, mn, 1, s, mn, info);
                Clascl("General", 0, 0, norma, bignum, m, nrhs, b, ldb, info);
            } else if (scale == 3) {
                //
                //              matrix scaled down
                //
                Clascl("General", 0, 0, norma, smlnum, m, n, a, lda, info);
                Rlascl("General", 0, 0, norma, smlnum, mn, 1, s, mn, info);
                Clascl("General", 0, 0, norma, smlnum, m, nrhs, b, ldb, info);
            } else {
                Mxerbla("Cqrt15", 1);
                return;
            }
        }
    }
    //
    norma = Rasum(mn, s, 1);
    normb = Clange("One-norm", m, nrhs, b, ldb, dummy);
    //
    //     End of Cqrt15
    //
}
