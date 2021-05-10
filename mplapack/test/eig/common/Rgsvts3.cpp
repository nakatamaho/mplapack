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

void Rgsvts3(INTEGER const m, INTEGER const p, INTEGER const n, REAL *a, REAL *af, INTEGER const lda, REAL *b, REAL *bf, INTEGER const ldb, REAL *u, INTEGER const ldu, REAL *v, INTEGER const ldv, REAL *q, INTEGER const ldq, REAL *alpha, REAL *beta, REAL *r, INTEGER const ldr, INTEGER *iwork, REAL *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    INTEGER ldaf = lda;
    INTEGER ldbf = ldb;
    REAL ulp = Rlamch("Precision");
    const REAL one = 1.0;
    REAL ulpinv = one / ulp;
    REAL unfl = Rlamch("Safe minimum");
    //
    //     Copy the matrix A to the array AF.
    //
    Rlacpy("Full", m, n, a, lda, af, lda);
    Rlacpy("Full", p, n, b, ldb, bf, ldb);
    //
    REAL anorm = max({Rlange("1", m, n, a, lda, rwork), unfl});
    REAL bnorm = max({Rlange("1", p, n, b, ldb, rwork), unfl});
    //
    //     Factorize the matrices A and B in the arrays AF and BF.
    //
    INTEGER k = 0;
    INTEGER l = 0;
    INTEGER info = 0;
    Rggsvd3("U", "V", "Q", m, n, p, k, l, af, lda, bf, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info);
    //
    //     Copy R
    //
    INTEGER i = 0;
    INTEGER j = 0;
    for (i = 1; i <= min(k + l, m); i = i + 1) {
        for (j = i; j <= k + l; j = j + 1) {
            r[(i - 1) + (j - 1) * ldr] = af[(i - 1) + ((n - k - l + j) - 1) * ldaf];
        }
    }
    //
    if (m - k - l < 0) {
        for (i = m + 1; i <= k + l; i = i + 1) {
            for (j = i; j <= k + l; j = j + 1) {
                r[(i - 1) + (j - 1) * ldr] = bf[((i - k) - 1) + ((n - k - l + j) - 1) * ldbf];
            }
        }
    }
    //
    //     Compute A:= U'*A*Q - D1*R
    //
    const REAL zero = 0.0;
    Rgemm("No transpose", "No transpose", m, n, n, one, a, lda, q, ldq, zero, work, lda);
    //
    Rgemm("Transpose", "No transpose", m, n, m, one, u, ldu, work, lda, zero, a, lda);
    //
    for (i = 1; i <= k; i = i + 1) {
        for (j = i; j <= k + l; j = j + 1) {
            a[(i - 1) + ((n - k - l + j) - 1) * lda] = a[(i - 1) + ((n - k - l + j) - 1) * lda] - r[(i - 1) + (j - 1) * ldr];
        }
    }
    //
    for (i = k + 1; i <= min(k + l, m); i = i + 1) {
        for (j = i; j <= k + l; j = j + 1) {
            a[(i - 1) + ((n - k - l + j) - 1) * lda] = a[(i - 1) + ((n - k - l + j) - 1) * lda] - alpha[i - 1] * r[(i - 1) + (j - 1) * ldr];
        }
    }
    //
    //     Compute norm( U'*A*Q - D1*R ) / ( MAX(1,M,N)*norm(A)*ULP ) .
    //
    REAL resid = Rlange("1", m, n, a, lda, rwork);
    //
    if (anorm > zero) {
        result[1 - 1] = ((resid / castREAL(max({(INTEGER)1, m, n}))) / anorm) / ulp;
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute B := V'*B*Q - D2*R
    //
    Rgemm("No transpose", "No transpose", p, n, n, one, b, ldb, q, ldq, zero, work, ldb);
    //
    Rgemm("Transpose", "No transpose", p, n, p, one, v, ldv, work, ldb, zero, b, ldb);
    //
    for (i = 1; i <= l; i = i + 1) {
        for (j = i; j <= l; j = j + 1) {
            b[(i - 1) + ((n - l + j) - 1) * ldb] = b[(i - 1) + ((n - l + j) - 1) * ldb] - beta[(k + i) - 1] * r[((k + i) - 1) + ((k + j) - 1) * ldr];
        }
    }
    //
    //     Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .
    //
    resid = Rlange("1", p, n, b, ldb, rwork);
    if (bnorm > zero) {
        result[2 - 1] = ((resid / castREAL(max({(INTEGER)1, p, n}))) / bnorm) / ulp;
    } else {
        result[2 - 1] = zero;
    }
    //
    //     Compute I - U'*U
    //
    Rlaset("Full", m, m, zero, one, work, ldq);
    Rsyrk("Upper", "Transpose", m, m, -one, u, ldu, one, work, ldu);
    //
    //     Compute norm( I - U'*U ) / ( M * ULP ) .
    //
    resid = Rlansy("1", "Upper", m, work, ldu, rwork);
    result[3 - 1] = (resid / castREAL(max((INTEGER)1, m))) / ulp;
    //
    //     Compute I - V'*V
    //
    Rlaset("Full", p, p, zero, one, work, ldv);
    Rsyrk("Upper", "Transpose", p, p, -one, v, ldv, one, work, ldv);
    //
    //     Compute norm( I - V'*V ) / ( P * ULP ) .
    //
    resid = Rlansy("1", "Upper", p, work, ldv, rwork);
    result[4 - 1] = (resid / castREAL(max((INTEGER)1, p))) / ulp;
    //
    //     Compute I - Q'*Q
    //
    Rlaset("Full", n, n, zero, one, work, ldq);
    Rsyrk("Upper", "Transpose", n, n, -one, q, ldq, one, work, ldq);
    //
    //     Compute norm( I - Q'*Q ) / ( N * ULP ) .
    //
    resid = Rlansy("1", "Upper", n, work, ldq, rwork);
    result[5 - 1] = (resid / castREAL(max((INTEGER)1, n))) / ulp;
    //
    //     Check sorting
    //
    Rcopy(n, alpha, 1, work, 1);
    REAL temp = 0.0;
    for (i = k + 1; i <= min(k + l, m); i = i + 1) {
        j = iwork[i - 1];
        if (i != j) {
            temp = work[i - 1];
            work[i - 1] = work[j - 1];
            work[j - 1] = temp;
        }
    }
    //
    result[6 - 1] = zero;
    for (i = k + 1; i <= min(k + l, m) - 1; i = i + 1) {
        if (work[i - 1] < work[(i + 1) - 1]) {
            result[6 - 1] = ulpinv;
        }
    }
    //
    //     End of Rgsvts3
    //
}
