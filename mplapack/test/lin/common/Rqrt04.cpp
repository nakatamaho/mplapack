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

void Rqrt04(INTEGER const m, INTEGER const n, INTEGER const nb, REAL *result) {
    FEM_CMN_SVE(Rqrt04);
    result([6]);
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Data statements ..
    //
    REAL eps = Rlamch("Epsilon");
    INTEGER k = min(m, n);
    INTEGER l = max(m, n);
    INTEGER lwork = max((INTEGER)2, l) * max((INTEGER)2, l) * nb;
    //
    //     Dynamically allocate local arrays
    //
    //     Put random numbers into A and copy to AF
    //
    INTEGER ldt = nb;
    INTEGER j = 0;
    REAL a[m * n];
    for (j = 1; j <= n; j = j + 1) {
        Rlarnv(2, iseed, m, &a[(j - 1) * lda]);
    }
    REAL af[m * n];
    Rlacpy("Full", m, n, a, m, af, m);
    //
    //     Factor the matrix A in the array AF.
    //
    REAL t[nb * n];
    REAL work[lwork];
    INTEGER info = 0;
    Rgeqrt(m, n, nb, af, m, t, ldt, work, info);
    //
    //     Generate the m-by-m matrix Q
    //
    const REAL zero = 0.0f;
    const REAL one = 1.0f;
    REAL q[m * m];
    Rlaset("Full", m, m, zero, one, q, m);
    Rgemqrt("R", "N", m, m, k, nb, af, m, t, ldt, q, m, work, info);
    //
    //     Copy R
    //
    REAL r[m * l];
    Rlaset("Full", m, n, zero, zero, r, m);
    Rlacpy("Upper", m, n, af, m, r, m);
    //
    //     Compute |R - Q'*A| / |A| and store in RESULT(1)
    //
    Rgemm("T", "N", m, n, m, -one, q, m, a, m, one, r, m);
    REAL rwork[l];
    REAL anorm = Rlange("1", m, n, a, m, rwork);
    REAL resid = Rlange("1", m, n, r, m, rwork);
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * max((INTEGER)1, m) * anorm);
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute |I - Q'*Q| and store in RESULT(2)
    //
    Rlaset("Full", m, m, zero, one, r, m);
    Rsyrk("U", "C", m, m, -one, q, m, one, r, m);
    resid = Rlansy("1", "Upper", m, r, m, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, m));
    //
    //     Generate random m-by-n matrix C and a copy CF
    //
    REAL c[m * n];
    for (j = 1; j <= n; j = j + 1) {
        Rlarnv(2, iseed, m, &c[(j - 1) * ldc]);
    }
    REAL cnorm = Rlange("1", m, n, c, m, rwork);
    REAL cf[m * n];
    Rlacpy("Full", m, n, c, m, cf, m);
    //
    //     Apply Q to C as Q*C
    //
    Rgemqrt("L", "N", m, n, k, nb, af, m, t, nb, cf, m, work, info);
    //
    //     Compute |Q*C - Q*C| / |C|
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
    //     Apply Q to C as QT*C
    //
    Rgemqrt("L", "T", m, n, k, nb, af, m, t, nb, cf, m, work, info);
    //
    //     Compute |QT*C - QT*C| / |C|
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
    //     Apply Q to D as D*Q
    //
    Rgemqrt("R", "N", n, m, k, nb, af, m, t, nb, df, n, work, info);
    //
    //     Compute |D*Q - D*Q| / |D|
    //
    Rgemm("N", "N", n, m, m, -one, d, n, q, m, one, df, n);
    resid = Rlange("1", n, m, df, n, rwork);
    if (cnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    //     Copy D into DF again
    //
    Rlacpy("Full", n, m, d, n, df, n);
    //
    //     Apply Q to D as D*QT
    //
    Rgemqrt("R", "T", n, m, k, nb, af, m, t, nb, df, n, work, info);
    //
    //     Compute |D*QT - D*QT| / |D|
    //
    Rgemm("N", "T", n, m, m, -one, d, n, q, m, one, df, n);
    resid = Rlange("1", n, m, df, n, rwork);
    if (cnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    //     Deallocate all arrays
    //
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t,c,d,cf,df)");
    //
}
