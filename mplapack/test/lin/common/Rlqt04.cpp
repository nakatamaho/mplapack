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

// Derived from LAPACK routine DLQT04.
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

void Rlqt04(INTEGER const m, INTEGER const n, INTEGER const nb, REAL *result) {
    static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
    REAL eps = Rlamch("Epsilon");
    INTEGER k = min(m, n);
    INTEGER ll = max(m, n);
    INTEGER lwork = max((INTEGER)2, ll) * max((INTEGER)2, ll) * nb;
    std::unique_ptr<REAL[]> a_storage;
    REAL *a = nullptr;
    a_storage = std::make_unique<REAL[]>(max((INTEGER)1, m * n));
    a = a_storage.get();
    std::unique_ptr<REAL[]> af_storage;
    REAL *af = nullptr;
    af_storage = std::make_unique<REAL[]>(max((INTEGER)1, m * n));
    af = af_storage.get();
    std::unique_ptr<REAL[]> q_storage;
    REAL *q = nullptr;
    q_storage = std::make_unique<REAL[]>(max((INTEGER)1, n * n));
    q = q_storage.get();
    std::unique_ptr<REAL[]> l_storage;
    REAL *l = nullptr;
    l_storage = std::make_unique<REAL[]>(max((INTEGER)1, ll * n));
    l = l_storage.get();
    std::unique_ptr<REAL[]> rwork_storage;
    REAL *rwork = nullptr;
    rwork_storage = std::make_unique<REAL[]>(max((INTEGER)1, ll));
    rwork = rwork_storage.get();
    std::unique_ptr<REAL[]> work_storage;
    REAL *work = nullptr;
    work_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
    work = work_storage.get();
    std::unique_ptr<REAL[]> t_storage;
    REAL *t = nullptr;
    t_storage = std::make_unique<REAL[]>(max((INTEGER)1, nb * n));
    t = t_storage.get();
    std::unique_ptr<REAL[]> c_storage;
    REAL *c = nullptr;
    c_storage = std::make_unique<REAL[]>(max((INTEGER)1, m * n));
    c = c_storage.get();
    std::unique_ptr<REAL[]> cf_storage;
    REAL *cf = nullptr;
    cf_storage = std::make_unique<REAL[]>(max((INTEGER)1, m * n));
    cf = cf_storage.get();
    std::unique_ptr<REAL[]> d_storage;
    REAL *d = nullptr;
    d_storage = std::make_unique<REAL[]>(max((INTEGER)1, n * m));
    d = d_storage.get();
    std::unique_ptr<REAL[]> df_storage;
    REAL *df = nullptr;
    df_storage = std::make_unique<REAL[]>(max((INTEGER)1, n * m));
    df = df_storage.get();
    INTEGER ldt = nb;
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        Rlarnv(2, iseed, m, &a[(j - 1) * m]);
    }
    Rlacpy("Full", m, n, a, m, af, m);
    INTEGER info = 0;
    Rgelqt(m, n, nb, af, m, t, ldt, work, info);
    const REAL zero = 0.0;
    const REAL one = 1.0;
    Rlaset("Full", n, n, zero, one, q, n);
    Rgemlqt("R", "N", n, n, k, nb, af, m, t, ldt, q, n, work, info);
    Rlaset("Full", m, n, zero, zero, l, ll);
    Rlacpy("Lower", m, n, af, m, l, ll);
    // Compute |L - A*Q'| / |A| and store in RESULT(1)
    //
    Rgemm("N", "T", m, n, n, -one, a, m, q, n, one, l, ll);
    REAL anorm = Rlange("1", m, n, a, m, rwork);
    REAL resid = Rlange("1", m, n, l, ll, rwork);
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * max((INTEGER)1, m) * anorm);
    } else {
        result[1 - 1] = zero;
    }
    //
    // Compute |I - Q'*Q| and store in RESULT(2)
    //
    Rlaset("Full", n, n, zero, one, l, ll);
    Rsyrk("U", "C", n, n, -one, q, n, one, l, ll);
    resid = Rlansy("1", "Upper", n, l, ll, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, n));
    //
    // Generate random m-by-n matrix C and a copy CF
    //
    for (j = 1; j <= m; j = j + 1) {
        Rlarnv(2, iseed, n, &d[(j - 1) * n]);
    }
    REAL dnorm = Rlange("1", n, m, d, n, rwork);
    Rlacpy("Full", n, m, d, n, df, n);
    //
    // Apply Q to C as Q*C
    //
    Rgemlqt("L", "N", n, m, k, nb, af, m, t, nb, df, n, work, info);
    //
    // Compute |Q*D - Q*D| / |D|
    //
    Rgemm("N", "N", n, m, n, -one, q, n, d, n, one, df, n);
    resid = Rlange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[3 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[3 - 1] = zero;
    }
    //
    // Copy D into DF again
    //
    Rlacpy("Full", n, m, d, n, df, n);
    //
    // Apply Q to D as QT*D
    //
    Rgemlqt("L", "T", n, m, k, nb, af, m, t, nb, df, n, work, info);
    //
    // Compute |QT*D - QT*D| / |D|
    //
    Rgemm("T", "N", n, m, n, -one, q, n, d, n, one, df, n);
    resid = Rlange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[4 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[4 - 1] = zero;
    }
    //
    // Generate random n-by-m matrix D and a copy DF
    //
    for (j = 1; j <= n; j = j + 1) {
        Rlarnv(2, iseed, m, &c[(j - 1) * m]);
    }
    REAL cnorm = Rlange("1", m, n, c, m, rwork);
    Rlacpy("Full", m, n, c, m, cf, m);
    //
    // Apply Q to C as C*Q
    //
    Rgemlqt("R", "N", m, n, k, nb, af, m, t, nb, cf, m, work, info);
    //
    // Compute |C*Q - C*Q| / |C|
    //
    Rgemm("N", "N", m, n, n, -one, c, m, q, n, one, cf, m);
    resid = Rlange("1", n, m, df, n, rwork);
    if (cnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    // Copy C into CF again
    //
    Rlacpy("Full", m, n, c, m, cf, m);
    //
    // Apply Q to D as D*QT
    //
    Rgemlqt("R", "T", m, n, k, nb, af, m, t, nb, cf, m, work, info);
    //
    // Compute |C*QT - C*QT| / |C|
    //
    Rgemm("N", "T", m, n, n, -one, c, m, q, n, one, cf, m);
    resid = Rlange("1", m, n, cf, m, rwork);
    if (cnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    // Deallocate all arrays
    //
}
