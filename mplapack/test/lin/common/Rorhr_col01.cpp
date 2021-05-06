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

void Rorhr_col01(INTEGER const m, INTEGER const n, INTEGER const mb1, INTEGER const nb1, INTEGER const nb2, REAL *result) {
    FEM_CMN_SVE(Rorhr_col01);
    result([6]);
    // COMMON srmnamc
    //
    // SAVE
    INTEGER *iseed(sve.iseed, [4]);
    //
    if (is_called_first_time) {
        static const INTEGER values[] = {1988, 1989, 1990, 1991};
        data_of_type<int>(FEM_VALUES_AND_SIZE), iseed;
    }
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     .. Return values ..
    //
    //  =====================================================================
    //
    //     ..
    //     .. Local allocatable arrays
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
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Data statements ..
    //
    //     TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
    //
    bool testzeros = false;
    //
    REAL eps = Rlamch("Epsilon");
    INTEGER k = min(m, n);
    INTEGER l = max({m, n, 1});
    //
    //     Dynamically allocate local arrays
    //
    //     Put random numbers into A and copy to AF
    //
    INTEGER j = 0;
    REAL a[m * n];
    for (j = 1; j <= n; j = j + 1) {
        Rlarnv(2, iseed, m, &a[(j - 1) * lda]);
    }
    if (testzeros) {
        if (m >= 4) {
            for (j = 1; j <= n; j = j + 1) {
                Rlarnv(2, iseed, m / 2, &a[((m / 4) - 1) + (j - 1) * lda]);
            }
        }
    }
    REAL af[m * n];
    Rlacpy("Full", m, n, a, m, af, m);
    //
    //     Number of row blocks in Rlatsqr
    //
    INTEGER nrb = max((INTEGER)1, ceiling((m - n).real() / (mb1 - n).real()));
    //
    //     Begin determine LWORK for the array WORK and allocate memory.
    //
    //     Rlatsqr requires NB1 to be bounded by N.
    //
    INTEGER nb1_ub = min(nb1, n);
    //
    //     Rgemqrt requires NB2 to be bounded by N.
    //
    INTEGER nb2_ub = min(nb2, n);
    //
    REAL t1[nb1 * n * nrb];
    REAL workquery[1];
    INTEGER info = 0;
    Rlatsqr(m, n, mb1, nb1_ub, af, m, t1, nb1, workquery, -1, info);
    INTEGER lwork = int(workquery[1 - 1]);
    Rorgtsqr(m, n, mb1, nb1, af, m, t1, nb1, workquery, -1, info);
    //
    lwork = max(lwork, int(workquery[1 - 1]));
    //
    //     In Rgemqrt, WORK is N*NB2_UB if SIDE = 'L',
    //                or  M*NB2_UB if SIDE = 'R'.
    //
    lwork = max({lwork, nb2_ub * n, nb2_ub * m});
    //
    //     End allocate memory for WORK.
    //
    //     Begin Householder reconstruction routines
    //
    //     Factor the matrix A in the array AF.
    //
    REAL work[lwork];
    Rlatsqr(m, n, mb1, nb1_ub, af, m, t1, nb1, work, lwork, info);
    //
    //     Copy the factor R into the array R.
    //
    REAL r[m * l];
    Rlacpy("U", n, n, af, m, r, m);
    //
    //     Reconstruct the orthogonal matrix Q.
    //
    Rorgtsqr(m, n, mb1, nb1, af, m, t1, nb1, work, lwork, info);
    //
    //     Perform the Householder reconstruction, the result is stored
    //     the arrays AF and T2.
    //
    REAL t2[nb2 * n];
    REAL diag[n];
    Rorhr_col(m, n, nb2, af, m, t2, nb2, diag, info);
    //
    //     Compute the factor R_hr corresponding to the Householder
    //     reconstructed Q_hr and place it in the upper triangle of AF to
    //     match the Q storage format in Rgeqrt. R_hr = R_tsqr * S,
    //     this means changing the sign of I-th row of the matrix R_tsqr
    //     according to sign of of I-th diagonal element DIAG(I) of the
    //     matrix S.
    //
    Rlacpy("U", n, n, r, m, af, m);
    //
    INTEGER i = 0;
    const REAL one = 1.0;
    for (i = 1; i <= n; i = i + 1) {
        if (diag[i - 1] == -one) {
            Rscal(n + 1 - i, -one, af[(i - 1) + (i - 1) * ldaf], m);
        }
    }
    //
    //     End Householder reconstruction routines.
    //
    //     Generate the m-by-m matrix Q
    //
    const REAL zero = 0.0;
    REAL q[l * l];
    Rlaset("Full", m, m, zero, one, q, m);
    //
    Rgemqrt("L", "N", m, m, k, nb2_ub, af, m, t2, nb2, q, m, work, info);
    //
    //     Copy R
    //
    Rlaset("Full", m, n, zero, zero, r, m);
    //
    Rlacpy("Upper", m, n, af, m, r, m);
    //
    //     TEST 1
    //     Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)
    //
    Rgemm("T", "N", m, n, m, -one, q, m, a, m, one, r, m);
    //
    REAL rwork[l];
    REAL anorm = Rlange("1", m, n, a, m, rwork);
    REAL resid = Rlange("1", m, n, r, m, rwork);
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * max((INTEGER)1, m) * anorm);
    } else {
        result[1 - 1] = zero;
    }
    //
    //     TEST 2
    //     Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)
    //
    Rlaset("Full", m, m, zero, one, r, m);
    Rsyrk("U", "T", m, m, -one, q, m, one, r, m);
    resid = Rlansy("1", "Upper", m, r, m, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, m));
    //
    //     Generate random m-by-n matrix C
    //
    REAL c[m * n];
    for (j = 1; j <= n; j = j + 1) {
        Rlarnv(2, iseed, m, &c[(j - 1) * ldc]);
    }
    REAL cnorm = Rlange("1", m, n, c, m, rwork);
    REAL cf[m * n];
    Rlacpy("Full", m, n, c, m, cf, m);
    //
    //     Apply Q to C as Q*C = CF
    //
    Rgemqrt("L", "N", m, n, k, nb2_ub, af, m, t2, nb2, cf, m, work, info);
    //
    //     TEST 3
    //     Compute |CF - Q*C| / ( eps *  m * |C| )
    //
    Rgemm("N", "N", m, n, m, -one, q, m, c, m, one, cf, m);
    resid = Rlange("1", m, n, cf, m, rwork);
    if (cnorm > zero) {
        result[3 - 1] = resid / (eps * max((INTEGER)1, m) * cnorm);
    } else {
        result[3 - 1] = zero;
    }
    //
    //     Copy C into CF again
    //
    Rlacpy("Full", m, n, c, m, cf, m);
    //
    //     Apply Q to C as (Q**T)*C = CF
    //
    Rgemqrt("L", "T", m, n, k, nb2_ub, af, m, t2, nb2, cf, m, work, info);
    //
    //     TEST 4
    //     Compute |CF - (Q**T)*C| / ( eps * m * |C|)
    //
    Rgemm("T", "N", m, n, m, -one, q, m, c, m, one, cf, m);
    resid = Rlange("1", m, n, cf, m, rwork);
    if (cnorm > zero) {
        result[4 - 1] = resid / (eps * max((INTEGER)1, m) * cnorm);
    } else {
        result[4 - 1] = zero;
    }
    //
    //     Generate random n-by-m matrix D and a copy DF
    //
    REAL d[n * m];
    for (j = 1; j <= m; j = j + 1) {
        Rlarnv(2, iseed, n, &d[(j - 1) * ldd]);
    }
    REAL dnorm = Rlange("1", n, m, d, n, rwork);
    REAL df[n * m];
    Rlacpy("Full", n, m, d, n, df, n);
    //
    //     Apply Q to D as D*Q = DF
    //
    Rgemqrt("R", "N", n, m, k, nb2_ub, af, m, t2, nb2, df, n, work, info);
    //
    //     TEST 5
    //     Compute |DF - D*Q| / ( eps * m * |D| )
    //
    Rgemm("N", "N", n, m, m, -one, d, n, q, m, one, df, n);
    resid = Rlange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    //     Copy D into DF again
    //
    Rlacpy("Full", n, m, d, n, df, n);
    //
    //     Apply Q to D as D*QT = DF
    //
    Rgemqrt("R", "T", n, m, k, nb2_ub, af, m, t2, nb2, df, n, work, info);
    //
    //     TEST 6
    //     Compute |DF - D*(Q**T)| / ( eps * m * |D| )
    //
    Rgemm("N", "T", n, m, m, -one, d, n, q, m, one, df, n);
    resid = Rlange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    //     Deallocate all arrays
    //
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t1,t2,diag,c,d,cf,d"
                        "f)");
    //
    //     End of Rorhr_col01
    //
}
