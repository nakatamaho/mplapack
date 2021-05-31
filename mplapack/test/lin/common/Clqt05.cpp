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

void Clqt05(INTEGER const m, INTEGER const n, INTEGER const l, INTEGER const nb, REAL *result) {
    //
    INTEGER iseed[] = {1988, 1989, 1990, 1991};
    REAL eps = Rlamch("Epsilon");
    INTEGER k = m;
    INTEGER n2 = m + n;
    INTEGER np1 = 0;
    if (n > 0) {
        np1 = m + 1;
    } else {
        np1 = 1;
    }
    INTEGER lwork = n2 * n2 * nb;
    //
    //     Dynamically allocate all arrays
    //
    //     Put random stuff into A
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    COMPLEX *a = new COMPLEX[m * n2];
    INTEGER lda = m;
    Claset("Full", m, n2, czero, czero, a, m);
    COMPLEX *t = new COMPLEX[nb * m];
    INTEGER ldt = nb;
    Claset("Full", nb, m, czero, czero, t, nb);
    INTEGER j = 0;
    for (j = 1; j <= m; j = j + 1) {
        Clarnv(2, iseed, m - j + 1, &a[(j - 1) + (j - 1) * lda]);
    }
    if (n > 0) {
        for (j = 1; j <= n - l; j = j + 1) {
            Clarnv(2, iseed, m, &a[(min(n + m, m + 1) + j - 1 - 1) * lda]);
        }
    }
    if (l > 0) {
        for (j = 1; j <= l; j = j + 1) {
            Clarnv(2, iseed, m - j + 1, &a[(j - 1) + (min(n + m, n + m - l + 1) + j - 1 - 1) * lda]);
        }
    }
    //
    //     Copy the matrix A to the array AF.
    //
    COMPLEX *af = new COMPLEX[m * n2];
    INTEGER ldaf = m;
    Clacpy("Full", m, n2, a, m, af, m);
    //
    //     Factor the matrix A in the array AF.
    //
    COMPLEX *work = new COMPLEX[lwork];
    INTEGER info = 0;
    Ctplqt(m, n, l, nb, af, m, &af[(np1 - 1) * ldaf], m, t, ldt, work, info);
    //
    //     Generate the (M+N)-by-(M+N) matrix Q by applying H to I
    //
    const COMPLEX one = COMPLEX(1.0, 0.0);
    COMPLEX *q = new COMPLEX[n2 * n2];
    INTEGER ldq = n2;
    Claset("Full", n2, n2, czero, one, q, n2);
    Cgemlqt("L", "N", n2, n2, k, nb, af, m, t, ldt, q, n2, work, info);
    //
    //     Copy L
    //
    COMPLEX *r = new COMPLEX[n2 * n2];
    INTEGER ldr = n2;
    Claset("Full", n2, n2, czero, czero, r, n2);
    Clacpy("Lower", m, n2, af, m, r, n2);
    //
    //     Compute |L - A*Q*C| / |A| and store in RESULT(1)
    //
    Cgemm("N", "C", m, n2, n2, -one, a, m, q, n2, one, r, n2);
    REAL *rwork = new REAL[n2];
    REAL anorm = Clange("1", m, n2, a, m, rwork);
    REAL resid = Clange("1", m, n2, r, n2, rwork);
    const REAL zero = 0.0f;
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * anorm * max((INTEGER)1, n2));
    } else {
        result[1 - 1] = zero;
    }
    //
    //     Compute |I - Q*Q'| and store in RESULT(2)
    //
    Claset("Full", n2, n2, czero, one, r, n2);
    Cherk("U", "N", n2, n2, (-one).real(), q, n2, one.real(), r, n2);
    resid = Clansy("1", "Upper", n2, r, n2, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, n2));
    //
    //     Generate random m-by-n matrix C and a copy CF
    //
    COMPLEX *c = new COMPLEX[n2 * m];
    INTEGER ldc = n2;
    Claset("Full", n2, m, czero, one, c, n2);
    for (j = 1; j <= m; j = j + 1) {
        Clarnv(2, iseed, n2, &c[(j - 1) * ldc]);
    }
    REAL cnorm = Clange("1", n2, m, c, n2, rwork);
    COMPLEX *cf = new COMPLEX[n2 * m];
    INTEGER ldcf = n2;
    Clacpy("Full", n2, m, c, n2, cf, n2);
    //
    //     Apply Q to C as Q*C
    //
    Ctpmlqt("L", "N", n, m, k, l, nb, &af[(np1 - 1) * ldaf], m, t, ldt, cf, n2, &cf[(np1 - 1)], n2, work, info);
    //
    //     Compute |Q*C - Q*C| / |C|
    //
    Cgemm("N", "N", n2, m, n2, -one, q, n2, c, n2, one, cf, n2);
    resid = Clange("1", n2, m, cf, n2, rwork);
    if (cnorm > zero) {
        result[3 - 1] = resid / (eps * max((INTEGER)1, n2) * cnorm);
    } else {
        result[3 - 1] = zero;
    }
    //
    //     Copy C into CF again
    //
    Clacpy("Full", n2, m, c, n2, cf, n2);
    //
    //     Apply Q to C as QT*C
    //
    Ctpmlqt("L", "C", n, m, k, l, nb, &af[(np1 - 1) * ldaf], m, t, ldt, cf, n2, &cf[(np1 - 1)], n2, work, info);
    //
    //     Compute |QT*C - QT*C| / |C|
    //
    Cgemm("C", "N", n2, m, n2, -one, q, n2, c, n2, one, cf, n2);
    resid = Clange("1", n2, m, cf, n2, rwork);
    //
    if (cnorm > zero) {
        result[4 - 1] = resid / (eps * max((INTEGER)1, n2) * cnorm);
    } else {
        result[4 - 1] = zero;
    }
    //
    //     Generate random m-by-n matrix D and a copy DF
    //
    COMPLEX *d = new COMPLEX[m * n2];
    INTEGER ldd = m;
    for (j = 1; j <= n2; j = j + 1) {
        Clarnv(2, iseed, m, &d[(j - 1) * ldd]);
    }
    REAL dnorm = Clange("1", m, n2, d, m, rwork);
    COMPLEX *df = new COMPLEX[m * n2];
    INTEGER lddf = m;
    Clacpy("Full", m, n2, d, m, df, m);
    //
    //     Apply Q to D as D*Q
    //
    Ctpmlqt("R", "N", m, n, k, l, nb, &af[(np1 - 1) * ldaf], m, t, ldt, df, m, &df[(np1 - 1) * lddf], m, work, info);
    //
    //     Compute |D*Q - D*Q| / |D|
    //
    Cgemm("N", "N", m, n2, n2, -one, d, m, q, n2, one, df, m);
    resid = Clange("1", m, n2, df, m, rwork);
    if (cnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, n2) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    //     Copy D into DF again
    //
    Clacpy("Full", m, n2, d, m, df, m);
    //
    //     Apply Q to D as D*QT
    //
    Ctpmlqt("R", "C", m, n, k, l, nb, &af[(np1 - 1) * ldaf], m, t, ldt, df, m, &df[(np1 - 1) * lddf], m, work, info);
    //
    //     Compute |D*QT - D*QT| / |D|
    //
    Cgemm("N", "C", m, n2, n2, -one, d, m, q, n2, one, df, m);
    resid = Clange("1", m, n2, df, m, rwork);
    if (cnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, n2) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    //     Deallocate all arrays
    //
    delete[] t;
    delete[] af;
    delete[] work;
    delete[] q;
    delete[] r;
    delete[] c;
    delete[] cf;
    delete[] d;
    delete[] df;
}
