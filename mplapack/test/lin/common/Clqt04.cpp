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

void Clqt04(INTEGER const m, INTEGER const n, INTEGER const nb, REAL *result) {
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
    INTEGER ll = max(m, n);
    INTEGER lwork = max((INTEGER)2, ll) * max((INTEGER)2, ll) * nb;
    //
    //     Dynamically allocate local arrays
    //
    //     Put random numbers into A and copy to AF
    //
    INTEGER ldt = nb;
    INTEGER j = 0;
    COMPLEX *a = new COMPLEX[m * n];
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &a[(j - 1) * lda]);
    }
    COMPLEX *af = new COMPLEX[m * n];
    Clacpy("Full", m, n, a, m, af, m);
    //
    //     Factor the matrix A in the array AF.
    //
    COMPLEX *t = new COMPLEX[nb * n];
    COMPLEX *work = new COMPLEX[lwork];
    INTEGER info = 0;
    Cgelqt(m, n, nb, af, m, t, ldt, work, info);
    //
    //     Generate the n-by-n matrix Q
    //
    const COMPLEX czero = COMPLEX(0.0f, 0.0f);
    const COMPLEX one = COMPLEX(1.0f, 0.0f);
    COMPLEX *q = new COMPLEX[n * n];
    Claset("Full", n, n, czero, one, q, n);
    Cgemlqt("R", "N", n, n, k, nb, af, m, t, ldt, q, n, work, info);
    //
    //     Copy L
    //
    COMPLEX *l = new COMPLEX[ll * n];
    Claset("Full", ll, n, czero, czero, l, ll);
    Clacpy("Lower", m, n, af, m, l, ll);
    //
    //     Compute |L - A*Q'| / |A| and store in RESULT(1)
    //
    Cgemm("N", "C", m, n, n, -one, a, m, q, n, one, l, ll);
    COMPLEX *rwork = new COMPLEX[ll];
    REAL anorm = Clange("1", m, n, a, m, rwork);
    REAL resid = Clange("1", m, n, l, ll, rwork);
    const REAL zero = 0.0f;
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * max((INTEGER)1, m) * anorm);
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute |I - Q'*Q| and store in RESULT(2)
    //
    Claset("Full", n, n, czero, one, l, ll);
    Cherk("U", "C", n, n, dreal[-one - 1], q, n, dreal[one - 1], l, ll);
    resid = Clansy("1", "Upper", n, l, ll, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, n));
    //
    //     Generate random m-by-n matrix C and a copy CF
    //
    COMPLEX *d = new COMPLEX[n * m];
    for (j = 1; j <= m; j = j + 1) {
        Clarnv(2, iseed, n, &d[(j - 1) * ldd]);
    }
    REAL dnorm = Clange("1", n, m, d, n, rwork);
    COMPLEX *df = new COMPLEX[n * m];
    Clacpy("Full", n, m, d, n, df, n);
    //
    //     Apply Q to C as Q*C
    //
    Cgemlqt("L", "N", n, m, k, nb, af, m, t, nb, df, n, work, info);
    //
    //     Compute |Q*D - Q*D| / |D|
    //
    Cgemm("N", "N", n, m, n, -one, q, n, d, n, one, df, n);
    resid = Clange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[3 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[3 - 1] = zero;
    }
    //
    //     Copy D into DF again
    //
    Clacpy("Full", n, m, d, n, df, n);
    //
    //     Apply Q to D as QT*D
    //
    Cgemlqt("L", "C", n, m, k, nb, af, m, t, nb, df, n, work, info);
    //
    //     Compute |QT*D - QT*D| / |D|
    //
    Cgemm("C", "N", n, m, n, -one, q, n, d, n, one, df, n);
    resid = Clange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[4 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[4 - 1] = zero;
    }
    //
    //     Generate random n-by-m matrix D and a copy DF
    //
    COMPLEX c[m * n];
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &c[(j - 1) * ldc]);
    }
    REAL cnorm = Clange("1", m, n, c, m, rwork);
    COMPLEX cf[m * n];
    Clacpy("Full", m, n, c, m, cf, m);
    //
    //     Apply Q to C as C*Q
    //
    Cgemlqt("R", "N", m, n, k, nb, af, m, t, nb, cf, m, work, info);
    //
    //     Compute |C*Q - C*Q| / |C|
    //
    Cgemm("N", "N", m, n, n, -one, c, m, q, n, one, cf, m);
    resid = Clange("1", n, m, df, n, rwork);
    if (cnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    //     Copy C into CF again
    //
    Clacpy("Full", m, n, c, m, cf, m);
    //
    //     Apply Q to D as D*QT
    //
    Cgemlqt("R", "C", m, n, k, nb, af, m, t, nb, cf, m, work, info);
    //
    //     Compute |C*QT - C*QT| / |C|
    //
    Cgemm("N", "C", m, n, n, -one, c, m, q, n, one, cf, m);
    resid = Clange("1", m, n, cf, m, rwork);
    if (cnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    //     Deallocate all arrays
    //
    //
    delete[] a;
    delete[] af;
    delete[] t;
    delete[] work;
    delete[] q;
    delete[] l;
    delete[] rwork;
    delete[] d;
    delete[] df;
}
