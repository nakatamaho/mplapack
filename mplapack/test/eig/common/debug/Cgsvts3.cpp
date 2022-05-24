/*
 * Copyright (c) 2021-2022
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
#include <lapacke.h>

void Cgsvts3(INTEGER const m, INTEGER const p, INTEGER const n, COMPLEX *a, COMPLEX *af, INTEGER const lda, COMPLEX *b, COMPLEX *bf, INTEGER const ldb, COMPLEX *u, INTEGER const ldu, COMPLEX *v, INTEGER const ldv, COMPLEX *q, INTEGER const ldq, REAL *alpha, REAL *beta, COMPLEX *r, INTEGER const ldr, INTEGER *iwork, COMPLEX *work, INTEGER const lwork, REAL *rwork, REAL *result) {
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
    printf("a="); printmat(m, n, af, lda); printf("\n");
    printf("b="); printmat(p, n, bf, ldb); printf("\n");
    Cggsvd3("U", "V", "Q", m, n, p, k, l, af, lda, bf, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, rwork, iwork, info);
    printf("u="); printmat(m, m, u, ldu); printf("\n");
    printf("v="); printmat(p, p, v, ldv); printf("\n");
    printf("q="); printmat(n, n, q, ldq); printf("\n");
    {
        __complex__ double *a_d = new __complex__ double[max(m * n, (INTEGER)1)];
        __complex__ double *b_d = new __complex__ double[max(p * n, (INTEGER)1)];
        __complex__ double *u_d = new __complex__ double[max(m * m, (INTEGER)1)];
        __complex__ double *v_d = new __complex__ double[max(p * p, (INTEGER)1)];
        __complex__ double *q_d = new __complex__ double[max(n * n, (INTEGER)1)];
        int *iwork_d = new int[max(n, (INTEGER)1)];
        int lda_d = m;
        int ldb_d = p;
        int ldu_d = m;
        int ldv_d = p;
        int ldq_d = n;
        int k_d, l_d;
        double *alpha_d = new double[max(n, (INTEGER)1)];
        double *beta_d = new double[max(n, (INTEGER)1)];
        double dtmp_r, dtmp_i;
        for (int pp = 0; pp < m; pp++) {
            for (int qq = 0; qq < n; qq++) {
                dtmp_r = cast2double(a[pp + qq * lda].real());
                dtmp_i = cast2double(a[pp + qq * lda].imag());
                __real__ a_d[pp + qq * lda_d] = dtmp_r;
                __imag__ a_d[pp + qq * lda_d] = dtmp_i;
            }
        }
        for (int pp = 0; pp < p; pp++) {
            for (int qq = 0; qq < n; qq++) {
                dtmp_r = cast2double(b[pp + qq * ldb].real());
                dtmp_i = cast2double(b[pp + qq * ldb].imag());
                __real__ b_d[pp + qq * ldb_d] = dtmp_r;
                __imag__ b_d[pp + qq * ldb_d] = dtmp_i;
	    }
        }
	printf("\n");
        printf("a_d="); printmat(m, n, a_d, lda_d); printf("\n");
        printf("b_d="); printmat(p, n, b_d, ldb_d); printf("\n");
        LAPACKE_zggsvd3(LAPACK_COL_MAJOR, 'U', 'V', 'Q', (int)m, (int)n, (int)p, &k_d, &l_d, a_d, lda_d, b_d, ldb_d, alpha_d, beta_d, u_d, ldu_d, v_d, ldv_d, q_d, ldq_d, iwork_d);
        printf("u_d="); printmat(m, m, u_d, ldu_d); printf("\n");
        printf("v_d="); printmat(p, p, v_d, ldv_d); printf("\n");
        printf("q_d="); printmat(n, n, q_d, ldq_d); printf("\n");
        delete[] beta_d;
        delete[] alpha_d;
        delete[] iwork_d;
        delete[] q_d;
        delete[] u_d;
        delete[] v_d;
        delete[] b_d;
        delete[] a_d;
    }

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
        result[1 - 1] = ((resid / castREAL(max({(INTEGER)1, m, n}))) / anorm) / ulp;
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
        result[2 - 1] = ((resid / castREAL(max({(INTEGER)1, p, n}))) / bnorm) / ulp;
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
    result[3 - 1] = (resid / castREAL(max((INTEGER)1, m))) / ulp;
    //
    //     Compute I - V'*V
    //
    Claset("Full", p, p, czero, cone, work, ldv);
    Cherk("Upper", "Conjugate transpose", p, p, -one, v, ldv, one, work, ldv);
    //
    //     Compute norm( I - V'*V ) / ( P * ULP ) .
    //
    resid = Clanhe("1", "Upper", p, work, ldv, rwork);
    result[4 - 1] = (resid / castREAL(max((INTEGER)1, p))) / ulp;
    //
    //     Compute I - Q'*Q
    //
    Claset("Full", n, n, czero, cone, work, ldq);
    Cherk("Upper", "Conjugate transpose", n, n, -one, q, ldq, one, work, ldq);
    //
    //     Compute norm( I - Q'*Q ) / ( N * ULP ) .
    //
    resid = Clanhe("1", "Upper", n, work, ldq, rwork);
    result[5 - 1] = (resid / castREAL(max((INTEGER)1, n))) / ulp;
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
