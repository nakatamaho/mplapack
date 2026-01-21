/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ZLQT05.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>
#include <memory>

void Clqt05(INTEGER const m, INTEGER const n, INTEGER const l, INTEGER const nb, REAL *result) {
    static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
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
    INTEGER ldt = nb;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    std::unique_ptr<COMPLEX[]> __a_storage(new COMPLEX[m * n2]);
    COMPLEX *a = __a_storage.get();
    Claset("Full", m, n2, czero, czero, a, m);
    std::unique_ptr<COMPLEX[]> __t_storage(new COMPLEX[nb * m]);
    COMPLEX *t = __t_storage.get();
    Claset("Full", nb, m, czero, czero, t, nb);
    INTEGER j = 0;
    for (j = 1; j <= m; j = j + 1) {
        Clarnv(2, iseed, m - j + 1, &a[(j - 1) + (j - 1) * m]);
    }
    if (n > 0) {
        for (j = 1; j <= n - l; j = j + 1) {
            Clarnv(2, iseed, m, &a[((min(n + m, m + 1) + j - 1) - 1) * m]);
        }
    }
    if (l > 0) {
        for (j = 1; j <= l; j = j + 1) {
            Clarnv(2, iseed, m - j + 1, &a[(j - 1) + ((min(n + m, n + m - l + 1) + j - 1) - 1) * m]);
        }
    }
    std::unique_ptr<COMPLEX[]> __af_storage(new COMPLEX[m * n2]);
    COMPLEX *af = __af_storage.get();
    Clacpy("Full", m, n2, a, m, af, m);
    std::unique_ptr<COMPLEX[]> __work_storage(new COMPLEX[lwork]);
    COMPLEX *work = __work_storage.get();
    INTEGER info = 0;
    Ctplqt(m, n, l, nb, af, m, &af[(np1 - 1) * m], m, t, ldt, work, info);
    const COMPLEX one = COMPLEX(1.0, 0.0);
    std::unique_ptr<COMPLEX[]> __q_storage(new COMPLEX[n2 * n2]);
    COMPLEX *q = __q_storage.get();
    Claset("Full", n2, n2, czero, one, q, n2);
    Cgemlqt("L", "N", n2, n2, k, nb, af, m, t, ldt, q, n2, work, info);
    std::unique_ptr<COMPLEX[]> __r_storage(new COMPLEX[n2 * n2]);
    COMPLEX *r = __r_storage.get();
    Claset("Full", n2, n2, czero, czero, r, n2);
    Clacpy("Lower", m, n2, af, m, r, n2);
    // Compute |L - A*Q*C| / |A| and store in RESULT(1)
    //
    Cgemm("N", "C", m, n2, n2, -one, a, m, q, n2, one, r, n2);
    std::unique_ptr<REAL[]> __rwork_storage(new REAL[n2]);
    REAL *rwork = __rwork_storage.get();
    REAL anorm = Clange("1", m, n2, a, m, rwork);
    REAL resid = Clange("1", m, n2, r, n2, rwork);
    const REAL zero = 0.0;
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * anorm * max((INTEGER)1, n2));
    } else {
        result[1 - 1] = zero;
    }
    //
    // Compute |I - Q*Q'| and store in RESULT(2)
    //
    Claset("Full", n2, n2, czero, one, r, n2);
    Cherk("U", "N", n2, n2, (-one).real(), q, n2, one.real(), r, n2);
    resid = Clansy("1", "Upper", n2, r, n2, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, n2));
    //
    // Generate random m-by-n matrix C and a copy CF
    //
    std::unique_ptr<COMPLEX[]> __c_storage(new COMPLEX[n2 * m]);
    COMPLEX *c = __c_storage.get();
    Claset("Full", n2, m, czero, one, c, n2);
    for (j = 1; j <= m; j = j + 1) {
        Clarnv(2, iseed, n2, &c[(j - 1) * n2]);
    }
    REAL cnorm = Clange("1", n2, m, c, n2, rwork);
    std::unique_ptr<COMPLEX[]> __cf_storage(new COMPLEX[n2 * m]);
    COMPLEX *cf = __cf_storage.get();
    Clacpy("Full", n2, m, c, n2, cf, n2);
    //
    // Apply Q to C as Q*C
    //
    Ctpmlqt("L", "N", n, m, k, l, nb, &af[(np1 - 1) * m], m, t, ldt, cf, n2, &cf[(np1 - 1)], n2, work, info);
    //
    // Compute |Q*C - Q*C| / |C|
    //
    Cgemm("N", "N", n2, m, n2, -one, q, n2, c, n2, one, cf, n2);
    resid = Clange("1", n2, m, cf, n2, rwork);
    if (cnorm > zero) {
        result[3 - 1] = resid / (eps * max((INTEGER)1, n2) * cnorm);
    } else {
        result[3 - 1] = zero;
    }
    //
    // Copy C into CF again
    //
    Clacpy("Full", n2, m, c, n2, cf, n2);
    //
    // Apply Q to C as QT*C
    //
    Ctpmlqt("L", "C", n, m, k, l, nb, &af[(np1 - 1) * m], m, t, ldt, cf, n2, &cf[(np1 - 1)], n2, work, info);
    //
    // Compute |QT*C - QT*C| / |C|
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
    // Generate random m-by-n matrix D and a copy DF
    //
    std::unique_ptr<COMPLEX[]> __d_storage(new COMPLEX[m * n2]);
    COMPLEX *d = __d_storage.get();
    for (j = 1; j <= n2; j = j + 1) {
        Clarnv(2, iseed, m, &d[(j - 1) * m]);
    }
    REAL dnorm = Clange("1", m, n2, d, m, rwork);
    std::unique_ptr<COMPLEX[]> __df_storage(new COMPLEX[m * n2]);
    COMPLEX *df = __df_storage.get();
    Clacpy("Full", m, n2, d, m, df, m);
    //
    // Apply Q to D as D*Q
    //
    Ctpmlqt("R", "N", m, n, k, l, nb, &af[(np1 - 1) * m], m, t, ldt, df, m, &df[(np1 - 1) * m], m, work, info);
    //
    // Compute |D*Q - D*Q| / |D|
    //
    Cgemm("N", "N", m, n2, n2, -one, d, m, q, n2, one, df, m);
    resid = Clange("1", m, n2, df, m, rwork);
    if (cnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, n2) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    // Copy D into DF again
    //
    Clacpy("Full", m, n2, d, m, df, m);
    //
    // Apply Q to D as D*QT
    //
    Ctpmlqt("R", "C", m, n, k, l, nb, &af[(np1 - 1) * m], m, t, ldt, df, m, &df[(np1 - 1) * m], m, work, info);
    //
    // Compute |D*QT - D*QT| / |D|
    //
    Cgemm("N", "C", m, n2, n2, -one, d, m, q, n2, one, df, m);
    resid = Clange("1", m, n2, df, m, rwork);
    if (cnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, n2) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    // Deallocate all arrays
    //
}
