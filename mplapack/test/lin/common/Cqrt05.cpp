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

void Cqrt05(INTEGER const m, INTEGER const n, INTEGER const l, INTEGER const nb, REAL *result) {
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
    //     .. Data statements ..
    //
    INTEGER iseed[] = {1988, 1989, 1990, 1991};
    REAL eps = Rlamch("Epsilon");
    INTEGER k = n;
    INTEGER m2 = m + n;
    INTEGER np1 = 0;
    if (m > 0) {
        np1 = n + 1;
    } else {
        np1 = 1;
    }
    INTEGER lwork = m2 * m2 * nb;
    //
    //     Dynamically allocate all arrays
    //
    //     Put random stuff into A
    //
    const COMPLEX czero = COMPLEX(0.0f, 0.0f);
    COMPLEX *a = new COMPLEX[m2 * n];
    INTEGER lda = m2;
    Claset("Full", m2, n, czero, czero, a, m2);
    COMPLEX *t = new COMPLEX[nb * n];
    INTEGER ldt = nb;
    Claset("Full", nb, n, czero, czero, t, nb);
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, j, &a[(j - 1) * lda]);
    }
    if (m > 0) {
        for (j = 1; j <= n; j = j + 1) {
            Clarnv(2, iseed, m - l, &a[(min(n + m, n + 1) - 1) + (j - 1) * lda]);
        }
    }
    if (l > 0) {
        for (j = 1; j <= n; j = j + 1) {
            Clarnv(2, iseed, min(j, l), &a[(min(n + m, n + m - l + 1) - 1) + (j - 1) * lda]);
        }
    }
    //
    //     Copy the matrix A to the array AF.
    //
    COMPLEX *af = new COMPLEX[m2 * n];
    INTEGER ldaf = m2;
    Clacpy("Full", m2, n, a, m2, af, m2);
    //
    //     Factor the matrix A in the array AF.
    //
    COMPLEX *work = new COMPLEX[lwork];
    INTEGER info = 0;
    Ctpqrt(m, n, l, nb, af, m2, &af[(np1 - 1)], m2, t, ldt, work, info);
    //
    //     Generate the (M+N)-by-(M+N) matrix Q by applying H to I
    //
    const COMPLEX one = COMPLEX(1.0f, 0.0f);
    COMPLEX *q = new COMPLEX[m2 * m2];
    INTEGER ldq = m2;
    Claset("Full", m2, m2, czero, one, q, m2);
    Cgemqrt("R", "N", m2, m2, k, nb, af, m2, t, ldt, q, m2, work, info);
    //
    //     Copy R
    //
    COMPLEX *r = new COMPLEX[m2 * m2];
    INTEGER ldr = m2;
    Claset("Full", m2, n, czero, czero, r, m2);
    Clacpy("Upper", m2, n, af, m2, r, m2);
    //
    //     Compute |R - Q'*A| / |A| and store in RESULT(1)
    //
    Cgemm("C", "N", m2, n, m2, -one, q, m2, a, m2, one, r, m2);
    REAL *rwork = new REAL[m2];
    REAL anorm = Clange("1", m2, n, a, m2, rwork);
    REAL resid = Clange("1", m2, n, r, m2, rwork);
    const REAL zero = 0.0f;
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * anorm * max((INTEGER)1, m2));
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute |I - Q'*Q| and store in RESULT(2)
    //
    Claset("Full", m2, m2, czero, one, r, m2);
    Cherk("U", "C", m2, m2, -one.real(), q, m2, one.real(), r, m2);
    resid = Clansy("1", "Upper", m2, r, m2, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, m2));
    //
    //     Generate random m-by-n matrix C and a copy CF
    //
    COMPLEX *c = new COMPLEX[m2 * n];
    INTEGER ldc = m2;
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m2, &c[(j - 1) * ldc]);
    }
    REAL cnorm = Clange("1", m2, n, c, m2, rwork);
    COMPLEX *cf = new COMPLEX[m2 * n];
    INTEGER ldcf = m2;
    Clacpy("Full", m2, n, c, m2, cf, m2);
    //
    //     Apply Q to C as Q*C
    //
    Ctpmqrt("L", "N", m, n, k, l, nb, &af[(np1 - 1)], m2, t, ldt, cf, m2, &cf[(np1 - 1)], m2, work, info);
    //
    //     Compute |Q*C - Q*C| / |C|
    //
    Cgemm("N", "N", m2, n, m2, -one, q, m2, c, m2, one, cf, m2);
    resid = Clange("1", m2, n, cf, m2, rwork);
    if (cnorm > zero) {
        result[3 - 1] = resid / (eps * max((INTEGER)1, m2) * cnorm);
    } else {
        result[3 - 1] = zero;
    }
    //
    //     Copy C into CF again
    //
    Clacpy("Full", m2, n, c, m2, cf, m2);
    //
    //     Apply Q to C as QT*C
    //
    Ctpmqrt("L", "C", m, n, k, l, nb, &af[(np1 - 1)], m2, t, ldt, cf, m2, &cf[(np1 - 1)], m2, work, info);
    //
    //     Compute |QT*C - QT*C| / |C|
    //
    Cgemm("C", "N", m2, n, m2, -one, q, m2, c, m2, one, cf, m2);
    resid = Clange("1", m2, n, cf, m2, rwork);
    if (cnorm > zero) {
        result[4 - 1] = resid / (eps * max((INTEGER)1, m2) * cnorm);
    } else {
        result[4 - 1] = zero;
    }
    //
    //     Generate random n-by-m matrix D and a copy DF
    //
    COMPLEX *d = new COMPLEX[n * m2];
    INTEGER ldd = n;

    for (j = 1; j <= m2; j = j + 1) {
        Clarnv(2, iseed, n, &d[(j - 1) * ldd]);
    }
    REAL dnorm = Clange("1", n, m2, d, n, rwork);
    COMPLEX *df = new COMPLEX[n * m2];
    INTEGER lddf = n;
    Clacpy("Full", n, m2, d, n, df, n);
    //
    //     Apply Q to D as D*Q
    //
    Ctpmqrt("R", "N", n, m, n, l, nb, &af[(np1 - 1)], m2, t, ldt, df, n, &df[(np1 - 1) * lddf], n, work, info);
    //
    //     Compute |D*Q - D*Q| / |D|
    //
    Cgemm("N", "N", n, m2, m2, -one, d, n, q, m2, one, df, n);
    resid = Clange("1", n, m2, df, n, rwork);
    if (cnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, m2) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    //     Copy D into DF again
    //
    Clacpy("Full", n, m2, d, n, df, n);
    //
    //     Apply Q to D as D*QT
    //
    Ctpmqrt("R", "C", n, m, n, l, nb, &af[(np1 - 1)], m2, t, ldt, df, n, &df[(np1 - 1) * lddf], n, work, info);
    //
    //     Compute |D*QT - D*QT| / |D|
    //
    Cgemm("N", "C", n, m2, m2, -one, d, n, q, m2, one, df, n);
    resid = Clange("1", n, m2, df, n, rwork);
    if (cnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, m2) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    //     Deallocate all arrays
    //
    delete[] a;
    delete[] af;
    delete[] q;
    delete[] r;
    delete[] rwork;
    delete[] work;
    delete[] t;
    delete[] c;
    delete[] d;
    delete[] cf;
    delete[] df;
}
