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

void Cgsvts3(INTEGER const m, INTEGER const p, INTEGER const n, COMPLEX *a, COMPLEX *af, INTEGER const lda, COMPLEX *b, COMPLEX *bf, INTEGER const ldb, COMPLEX *u, INTEGER const ldu, COMPLEX *v, INTEGER const ldv, COMPLEX *q, INTEGER const ldq, REAL *alpha, REAL *beta, COMPLEX *r, INTEGER const ldr, INTEGER *iwork, COMPLEX *work, INTEGER const lwork, REAL *rwork, REAL *result) {
    a([lda * star]);
    af([lda * star]);
    b([ldb * star]);
    bf([ldb * star]);
    u([ldu * star]);
    v([ldv * star]);
    q([ldq * star]);
    r([ldr * star]);
    work([lwork]);
    result([6]);
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
    REAL ulp = Rlamch("Precision");
    const REAL one = 1.0;
    REAL ulpinv = one / ulp;
    REAL unfl = Rlamch("Safe minimum");
    //
    //     Copy the matrix A to the array AF.
    //
    Clacpy("Full", m, n, a, lda, af, lda);
    Clacpy("Full", p, n, b, ldb, bf, ldb);
    //
    REAL anorm = max({Clange("1", m, n, a, lda, rwork), unfl});
    REAL bnorm = max({Clange("1", p, n, b, ldb, rwork), unfl});
    //
    //     Factorize the matrices A and B in the arrays AF and BF.
    //
    INTEGER k = 0;
    INTEGER l = 0;
    INTEGER info = 0;
    Cggsvd3("U", "V", "Q", m, n, p, k, l, af, lda, bf, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, rwork, iwork, info);
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
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    Cgemm("No transpose", "No transpose", m, n, n, cone, a, lda, q, ldq, czero, work, lda);
    //
    Cgemm("Conjugate transpose", "No transpose", m, n, m, cone, u, ldu, work, lda, czero, a, lda);
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
    REAL resid = Clange("1", m, n, a, lda, rwork);
    const REAL zero = 0.0;
    if (anorm > zero) {
        result[1 - 1] = ((resid / (max({(INTEGER)1, m, n})).real()) / anorm) / ulp;
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute B := V'*B*Q - D2*R
    //
    Cgemm("No transpose", "No transpose", p, n, n, cone, b, ldb, q, ldq, czero, work, ldb);
    //
    Cgemm("Conjugate transpose", "No transpose", p, n, p, cone, v, ldv, work, ldb, czero, b, ldb);
    //
    for (i = 1; i <= l; i = i + 1) {
        for (j = i; j <= l; j = j + 1) {
            b[(i - 1) + ((n - l + j) - 1) * ldb] = b[(i - 1) + ((n - l + j) - 1) * ldb] - beta[(k + i) - 1] * r[((k + i) - 1) + ((k + j) - 1) * ldr];
        }
    }
    //
    //     Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .
    //
    resid = Clange("1", p, n, b, ldb, rwork);
    if (bnorm > zero) {
        result[2 - 1] = ((resid / (max({(INTEGER)1, p, n})).real()) / bnorm) / ulp;
    } else {
        result[2 - 1] = zero;
    }
    //
    //     Compute I - U'*U
    //
    Claset("Full", m, m, czero, cone, work, ldq);
    Cherk("Upper", "Conjugate transpose", m, m, -one, u, ldu, one, work, ldu);
    //
    //     Compute norm( I - U'*U ) / ( M * ULP ) .
    //
    resid = Clanhe("1", "Upper", m, work, ldu, rwork);
    result[3 - 1] = (resid / (max((INTEGER)1, m)).real()) / ulp;
    //
    //     Compute I - V'*V
    //
    Claset("Full", p, p, czero, cone, work, ldv);
    Cherk("Upper", "Conjugate transpose", p, p, -one, v, ldv, one, work, ldv);
    //
    //     Compute norm( I - V'*V ) / ( P * ULP ) .
    //
    resid = Clanhe("1", "Upper", p, work, ldv, rwork);
    result[4 - 1] = (resid / (max((INTEGER)1, p)).real()) / ulp;
    //
    //     Compute I - Q'*Q
    //
    Claset("Full", n, n, czero, cone, work, ldq);
    Cherk("Upper", "Conjugate transpose", n, n, -one, q, ldq, one, work, ldq);
    //
    //     Compute norm( I - Q'*Q ) / ( N * ULP ) .
    //
    resid = Clanhe("1", "Upper", n, work, ldq, rwork);
    result[5 - 1] = (resid / (max((INTEGER)1, n)).real()) / ulp;
    //
    //     Check sorting
    //
    Rcopy(n, alpha, 1, rwork, 1);
    REAL temp = 0.0;
    for (i = k + 1; i <= min(k + l, m); i = i + 1) {
        j = iwork[i - 1];
        if (i != j) {
            temp = rwork[i - 1];
            rwork[i - 1] = rwork[j - 1];
            rwork[j - 1] = temp;
        }
    }
    //
    result[6 - 1] = zero;
    for (i = k + 1; i <= min(k + l, m) - 1; i = i + 1) {
        if (rwork[i - 1] < rwork[(i + 1) - 1]) {
            result[6 - 1] = ulpinv;
        }
    }
    //
    //     End of Cgsvts3
    //
}
