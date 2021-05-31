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

void Cunhr_col01(INTEGER const m, INTEGER const n, INTEGER const mb1, INTEGER const nb1, INTEGER const nb2, REAL *result) {
    //
    //     TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
    //
    INTEGER iseed[4] = {1988, 1989, 1990, 1991};
    bool testzeros = false;
    //
    REAL eps = Rlamch("Epsilon");
    INTEGER k = min(m, n);
    INTEGER l = max({m, n, (INTEGER)1});
    //
    //     Dynamically allocate local arrays
    //
    //     Put random numbers into A and copy to AF
    //
    INTEGER j = 0;
    COMPLEX *a = new COMPLEX[m * n];
    INTEGER lda = m;
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &a[(j - 1) * lda]);
    }
    if (testzeros) {
        if (m >= 4) {
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(2, iseed, m / 2, &a[((m / 4) - 1) + (j - 1) * lda]);
            }
        }
    }
    COMPLEX *af = new COMPLEX[m * n];
    INTEGER ldaf = m;
    Clacpy("Full", m, n, a, m, af, m);
    //
    //     Number of row blocks in Clatsqr
    //
    INTEGER nrb = max((INTEGER)1, ceil(castREAL(m - n) / castREAL(mb1 - n)));
    //
    //     Begin determine LWORK for the array WORK and allocate memory.
    //
    //     Clatsqr requires NB1 to be bounded by N.
    //
    INTEGER nb1_ub = min(nb1, n);
    //
    //     Cgemqrt requires NB2 to be bounded by N.
    //
    INTEGER nb2_ub = min(nb2, n);
    //
    COMPLEX *t1 = new COMPLEX[nb1 * n * nrb];
    COMPLEX workquery[1];
    INTEGER info = 0;
    Clatsqr(m, n, mb1, nb1_ub, af, m, t1, nb1, workquery, -1, info);
    INTEGER lwork = castINTEGER(workquery[1 - 1].real());
    Cungtsqr(m, n, mb1, nb1, af, m, t1, nb1, workquery, -1, info);
    //
    lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
    //
    //     In Cgemqrt, WORK is N*NB2_UB if SIDE = 'L',
    //                or  M*NB2_UB if SIDE = 'R'.
    //
    lwork = max({lwork, nb2_ub * n, nb2_ub * m});
    //
    //
    //     End allocate memory for WORK.
    //
    //     Begin Householder reconstruction routines
    //
    //     Factor the matrix A in the array AF.
    //
    COMPLEX *work = new COMPLEX[lwork];
    Clatsqr(m, n, mb1, nb1_ub, af, m, t1, nb1, work, lwork, info);
    //
    //     Copy the factor R into the array R.
    //
    COMPLEX *r = new COMPLEX[m * l];
    INTEGER ldr = m;
    Clacpy("U", n, n, af, m, r, m);
    //
    //     Reconstruct the orthogonal matrix Q.
    //
    Cungtsqr(m, n, mb1, nb1, af, m, t1, nb1, work, lwork, info);
    //
    //     Perform the Householder reconstruction, the result is stored
    //     the arrays AF and T2.
    //
    COMPLEX *t2 = new COMPLEX[nb2 * n];
    COMPLEX *diag = new COMPLEX[n];
    Cunhr_col(m, n, nb2, af, m, t2, nb2, diag, info);
    //
    //     Compute the factor R_hr corresponding to the Householder
    //     reconstructed Q_hr and place it in the upper triangle of AF to
    //     match the Q storage format in Cgeqrt. R_hr = R_tsqr * S,
    //     this means changing the sign of I-th row of the matrix R_tsqr
    //     according to sign of of I-th diagonal element DIAG(I) of the
    //     matrix S.
    //
    Clacpy("U", n, n, r, m, af, m);
    //
    INTEGER i = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    for (i = 1; i <= n; i = i + 1) {
        if (diag[i - 1] == -cone) {
            Cscal(n + 1 - i, -cone, &af[(i - 1) + (i - 1) * ldaf], m);
        }
    }
    //
    //     End Householder reconstruction routines.
    //
    //     Generate the m-by-m matrix Q
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    COMPLEX *q = new COMPLEX[l * l];
    INTEGER ldq = l;
    Claset("Full", m, m, czero, cone, q, m);
    //
    Cgemqrt("L", "N", m, m, k, nb2_ub, af, m, t2, nb2, q, m, work, info);
    //
    //     Copy R
    //
    Claset("Full", m, n, czero, czero, r, m);
    //
    Clacpy("Upper", m, n, af, m, r, m);
    //
    //     TEST 1
    //     Compute |R - (Q**H)*A| / ( eps * m * |A| ) and store in RESULT(1)
    //
    Cgemm("C", "N", m, n, m, -cone, q, m, a, m, cone, r, m);
    //
    REAL *rwork = new REAL[l];
    REAL anorm = Clange("1", m, n, a, m, rwork);
    REAL resid = Clange("1", m, n, r, m, rwork);
    const REAL zero = 0.0;
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * max((INTEGER)1, m) * anorm);
    } else {
        result[1 - 1] = zero;
    }
    //
    //     TEST 2
    //     Compute |I - (Q**H)*Q| / ( eps * m ) and store in RESULT(2)
    //
    Claset("Full", m, m, czero, cone, r, m);
    Cherk("U", "C", m, m, -cone.real(), q, m, cone.real(), r, m);
    resid = Clansy("1", "Upper", m, r, m, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, m));
    //
    //     Generate random m-by-n matrix C
    //
    COMPLEX *c = new COMPLEX[m * n];
    INTEGER ldc = m;
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &c[(j - 1) * ldc]);
    }
    REAL cnorm = Clange("1", m, n, c, m, rwork);
    COMPLEX *cf = new COMPLEX[m * n];
    INTEGER ldcf = m;
    Clacpy("Full", m, n, c, m, cf, m);
    //
    //     Apply Q to C as Q*C = CF
    //
    Cgemqrt("L", "N", m, n, k, nb2_ub, af, m, t2, nb2, cf, m, work, info);
    //
    //     TEST 3
    //     Compute |CF - Q*C| / ( eps *  m * |C| )
    //
    Cgemm("N", "N", m, n, m, -cone, q, m, c, m, cone, cf, m);
    resid = Clange("1", m, n, cf, m, rwork);
    if (cnorm > zero) {
        result[3 - 1] = resid / (eps * max((INTEGER)1, m) * cnorm);
    } else {
        result[3 - 1] = zero;
    }
    //
    //     Copy C into CF again
    //
    Clacpy("Full", m, n, c, m, cf, m);
    //
    //     Apply Q to C as (Q**H)*C = CF
    //
    Cgemqrt("L", "C", m, n, k, nb2_ub, af, m, t2, nb2, cf, m, work, info);
    //
    //     TEST 4
    //     Compute |CF - (Q**H)*C| / ( eps * m * |C|)
    //
    Cgemm("C", "N", m, n, m, -cone, q, m, c, m, cone, cf, m);
    resid = Clange("1", m, n, cf, m, rwork);
    if (cnorm > zero) {
        result[4 - 1] = resid / (eps * max((INTEGER)1, m) * cnorm);
    } else {
        result[4 - 1] = zero;
    }
    //
    //     Generate random n-by-m matrix D and a copy DF
    //
    COMPLEX *d = new COMPLEX[n * m];
    INTEGER ldd = n;
    for (j = 1; j <= m; j = j + 1) {
        Clarnv(2, iseed, n, &d[(j - 1) * ldd]);
    }
    REAL dnorm = Clange("1", n, m, d, n, rwork);
    COMPLEX df[n * m];
    Clacpy("Full", n, m, d, n, df, n);
    //
    //     Apply Q to D as D*Q = DF
    //
    Cgemqrt("R", "N", n, m, k, nb2_ub, af, m, t2, nb2, df, n, work, info);
    //
    //     TEST 5
    //     Compute |DF - D*Q| / ( eps * m * |D| )
    //
    Cgemm("N", "N", n, m, m, -cone, d, n, q, m, cone, df, n);
    resid = Clange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    //     Copy D into DF again
    //
    Clacpy("Full", n, m, d, n, df, n);
    //
    //     Apply Q to D as D*QT = DF
    //
    Cgemqrt("R", "C", n, m, k, nb2_ub, af, m, t2, nb2, df, n, work, info);
    //
    //     TEST 6
    //     Compute |DF - D*(Q**H)| / ( eps * m * |D| )
    //
    Cgemm("N", "C", n, m, m, -cone, d, n, q, m, cone, df, n);
    resid = Clange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    //     Deallocate all arrays
    delete[] a;
    delete[] af;
    delete[] t1;
    delete[] work;
    delete[] rwork;
    delete[] r;
    delete[] t2;
    delete[] diag;
    delete[] q;
    delete[] c;
    delete[] cf;
    delete[] d;
    //
}
