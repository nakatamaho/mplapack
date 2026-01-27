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

// Derived from LAPACK routine DLQT05.
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

void Rlqt05(INTEGER const m, INTEGER const n, INTEGER const l, INTEGER const nb, REAL *result) {
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
    std::unique_ptr<REAL[]> a_storage;
    REAL *a = nullptr;
    a_storage = std::make_unique<REAL[]>(max((INTEGER)1, m * n2));
    a = a_storage.get();
    std::unique_ptr<REAL[]> af_storage;
    REAL *af = nullptr;
    af_storage = std::make_unique<REAL[]>(max((INTEGER)1, m * n2));
    af = af_storage.get();
    std::unique_ptr<REAL[]> q_storage;
    REAL *q = nullptr;
    q_storage = std::make_unique<REAL[]>(max((INTEGER)1, n2 * n2));
    q = q_storage.get();
    std::unique_ptr<REAL[]> r_storage;
    REAL *r = nullptr;
    r_storage = std::make_unique<REAL[]>(max((INTEGER)1, n2 * n2));
    r = r_storage.get();
    std::unique_ptr<REAL[]> rwork_storage;
    REAL *rwork = nullptr;
    rwork_storage = std::make_unique<REAL[]>(max((INTEGER)1, n2));
    rwork = rwork_storage.get();
    std::unique_ptr<REAL[]> work_storage;
    REAL *work = nullptr;
    work_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
    work = work_storage.get();
    std::unique_ptr<REAL[]> t_storage;
    REAL *t = nullptr;
    t_storage = std::make_unique<REAL[]>(max((INTEGER)1, nb * m));
    t = t_storage.get();
    std::unique_ptr<REAL[]> c_storage;
    REAL *c = nullptr;
    c_storage = std::make_unique<REAL[]>(max((INTEGER)1, n2 * m));
    c = c_storage.get();
    std::unique_ptr<REAL[]> cf_storage;
    REAL *cf = nullptr;
    cf_storage = std::make_unique<REAL[]>(max((INTEGER)1, n2 * m));
    cf = cf_storage.get();
    std::unique_ptr<REAL[]> d_storage;
    REAL *d = nullptr;
    d_storage = std::make_unique<REAL[]>(max((INTEGER)1, m * n2));
    d = d_storage.get();
    std::unique_ptr<REAL[]> df_storage;
    REAL *df = nullptr;
    df_storage = std::make_unique<REAL[]>(max((INTEGER)1, m * n2));
    df = df_storage.get();
    INTEGER ldt = nb;
    const REAL zero = 0.0;
    Rlaset("Full", m, n2, zero, zero, a, m);
    Rlaset("Full", nb, m, zero, zero, t, nb);
    INTEGER j = 0;
    for (j = 1; j <= m; j = j + 1) {
        Rlarnv(2, iseed, m - j + 1, &a[(j - 1) + (j - 1) * m]);
    }
    if (n > 0) {
        for (j = 1; j <= n - l; j = j + 1) {
            Rlarnv(2, iseed, m, &a[((min(n + m, m + 1) + j - 1) - 1) * m]);
        }
    }
    if (l > 0) {
        for (j = 1; j <= l; j = j + 1) {
            Rlarnv(2, iseed, m - j + 1, &a[(j - 1) + ((min(n + m, n + m - l + 1) + j - 1) - 1) * m]);
        }
    }
    Rlacpy("Full", m, n2, a, m, af, m);
    INTEGER info = 0;
    Rtplqt(m, n, l, nb, af, m, &af[(np1 - 1) * m], m, t, ldt, work, info);
    const REAL one = 1.0;
    Rlaset("Full", n2, n2, zero, one, q, n2);
    Rgemlqt("L", "N", n2, n2, k, nb, af, m, t, ldt, q, n2, work, info);
    Rlaset("Full", n2, n2, zero, zero, r, n2);
    Rlacpy("Lower", m, n2, af, m, r, n2);
    // Compute |L - A*Q*T| / |A| and store in RESULT(1)
    //
    Rgemm("N", "T", m, n2, n2, -one, a, m, q, n2, one, r, n2);
    REAL anorm = Rlange("1", m, n2, a, m, rwork);
    REAL resid = Rlange("1", m, n2, r, n2, rwork);
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * anorm * max((INTEGER)1, n2));
    } else {
        result[1 - 1] = zero;
    }
    //
    // Compute |I - Q*Q'| and store in RESULT(2)
    //
    Rlaset("Full", n2, n2, zero, one, r, n2);
    Rsyrk("U", "N", n2, n2, -one, q, n2, one, r, n2);
    resid = Rlansy("1", "Upper", n2, r, n2, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, n2));
    //
    // Generate random m-by-n matrix C and a copy CF
    //
    Rlaset("Full", n2, m, zero, one, c, n2);
    for (j = 1; j <= m; j = j + 1) {
        Rlarnv(2, iseed, n2, &c[(j - 1) * n2]);
    }
    REAL cnorm = Rlange("1", n2, m, c, n2, rwork);
    Rlacpy("Full", n2, m, c, n2, cf, n2);
    //
    // Apply Q to C as Q*C
    //
    Rtpmlqt("L", "N", n, m, k, l, nb, &af[(np1 - 1) * m], m, t, ldt, cf, n2, &cf[(np1 - 1)], n2, work, info);
    //
    // Compute |Q*C - Q*C| / |C|
    //
    Rgemm("N", "N", n2, m, n2, -one, q, n2, c, n2, one, cf, n2);
    resid = Rlange("1", n2, m, cf, n2, rwork);
    if (cnorm > zero) {
        result[3 - 1] = resid / (eps * max((INTEGER)1, n2) * cnorm);
    } else {
        result[3 - 1] = zero;
    }
    //
    // Copy C into CF again
    //
    Rlacpy("Full", n2, m, c, n2, cf, n2);
    //
    // Apply Q to C as QT*C
    //
    Rtpmlqt("L", "T", n, m, k, l, nb, &af[(np1 - 1) * m], m, t, ldt, cf, n2, &cf[(np1 - 1)], n2, work, info);
    //
    // Compute |QT*C - QT*C| / |C|
    //
    Rgemm("T", "N", n2, m, n2, -one, q, n2, c, n2, one, cf, n2);
    resid = Rlange("1", n2, m, cf, n2, rwork);
    //
    if (cnorm > zero) {
        result[4 - 1] = resid / (eps * max((INTEGER)1, n2) * cnorm);
    } else {
        result[4 - 1] = zero;
    }
    //
    // Generate random m-by-n matrix D and a copy DF
    //
    for (j = 1; j <= n2; j = j + 1) {
        Rlarnv(2, iseed, m, &d[(j - 1) * m]);
    }
    REAL dnorm = Rlange("1", m, n2, d, m, rwork);
    Rlacpy("Full", m, n2, d, m, df, m);
    //
    // Apply Q to D as D*Q
    //
    Rtpmlqt("R", "N", m, n, k, l, nb, &af[(np1 - 1) * m], m, t, ldt, df, m, &df[(np1 - 1) * m], m, work, info);
    //
    // Compute |D*Q - D*Q| / |D|
    //
    Rgemm("N", "N", m, n2, n2, -one, d, m, q, n2, one, df, m);
    resid = Rlange("1", m, n2, df, m, rwork);
    if (cnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, n2) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    // Copy D into DF again
    //
    Rlacpy("Full", m, n2, d, m, df, m);
    //
    // Apply Q to D as D*QT
    //
    Rtpmlqt("R", "T", m, n, k, l, nb, &af[(np1 - 1) * m], m, t, ldt, df, m, &df[(np1 - 1) * m], m, work, info);
    //
    // Compute |D*QT - D*QT| / |D|
    //
    Rgemm("N", "T", m, n2, n2, -one, d, m, q, n2, one, df, m);
    resid = Rlange("1", m, n2, df, m, rwork);
    if (cnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, n2) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    // Deallocate all arrays
    //
}
