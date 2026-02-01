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

// Derived from LAPACK routine ZQRT04.
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

void Cqrt04(INTEGER const m, INTEGER const n, INTEGER const nb, REAL *result) {
    static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
    REAL eps = Rlamch("Epsilon");
    INTEGER k = min(m, n);
    INTEGER l = max(m, n);
    INTEGER lwork = max((INTEGER)2, l) * max((INTEGER)2, l) * nb;
    std::unique_ptr<COMPLEX[]> a_storage;
    COMPLEX *a = nullptr;
    a_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m * n));
    a = a_storage.get();
    std::unique_ptr<COMPLEX[]> af_storage;
    COMPLEX *af = nullptr;
    af_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m * n));
    af = af_storage.get();
    std::unique_ptr<COMPLEX[]> q_storage;
    COMPLEX *q = nullptr;
    q_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m * m));
    q = q_storage.get();
    std::unique_ptr<COMPLEX[]> r_storage;
    COMPLEX *r = nullptr;
    r_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m * l));
    r = r_storage.get();
    std::unique_ptr<REAL[]> rwork_storage;
    REAL *rwork = nullptr;
    rwork_storage = std::make_unique<REAL[]>(max((INTEGER)1, l));
    rwork = rwork_storage.get();
    std::unique_ptr<COMPLEX[]> work_storage;
    COMPLEX *work = nullptr;
    work_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lwork));
    work = work_storage.get();
    std::unique_ptr<COMPLEX[]> t_storage;
    COMPLEX *t = nullptr;
    t_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, nb * n));
    t = t_storage.get();
    std::unique_ptr<COMPLEX[]> c_storage;
    COMPLEX *c = nullptr;
    c_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m * n));
    c = c_storage.get();
    std::unique_ptr<COMPLEX[]> cf_storage;
    COMPLEX *cf = nullptr;
    cf_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m * n));
    cf = cf_storage.get();
    std::unique_ptr<COMPLEX[]> d_storage;
    COMPLEX *d = nullptr;
    d_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, n * m));
    d = d_storage.get();
    std::unique_ptr<COMPLEX[]> df_storage;
    COMPLEX *df = nullptr;
    df_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, n * m));
    df = df_storage.get();
    INTEGER ldt = nb;
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &a[(j - 1) * m]);
    }
    Clacpy("Full", m, n, a, m, af, m);
    INTEGER info = 0;
    Cgeqrt(m, n, nb, af, m, t, ldt, work, info);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX one = COMPLEX(1.0, 0.0);
    Claset("Full", m, m, czero, one, q, m);
    Cgemqrt("R", "N", m, m, k, nb, af, m, t, ldt, q, m, work, info);
    Claset("Full", m, n, czero, czero, r, m);
    Clacpy("Upper", m, n, af, m, r, m);
    // Compute |R - Q'*A| / |A| and store in RESULT(1)
    //
    Cgemm("C", "N", m, n, m, -one, q, m, a, m, one, r, m);
    REAL anorm = Clange("1", m, n, a, m, rwork);
    REAL resid = Clange("1", m, n, r, m, rwork);
    const REAL zero = 0.0;
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * max((INTEGER)1, m) * anorm);
    } else {
        result[1 - 1] = zero;
    }
    //
    // Compute |I - Q'*Q| and store in RESULT(2)
    //
    Claset("Full", m, m, czero, one, r, m);
    Cherk("U", "C", m, m, (-one).real(), q, m, one.real(), r, m);
    resid = Clansy("1", "Upper", m, r, m, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, m));
    //
    // Generate random m-by-n matrix C and a copy CF
    //
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &c[(j - 1) * m]);
    }
    REAL cnorm = Clange("1", m, n, c, m, rwork);
    Clacpy("Full", m, n, c, m, cf, m);
    //
    // Apply Q to C as Q*C
    //
    Cgemqrt("L", "N", m, n, k, nb, af, m, t, nb, cf, m, work, info);
    //
    // Compute |Q*C - Q*C| / |C|
    //
    Cgemm("N", "N", m, n, m, -one, q, m, c, m, one, cf, m);
    resid = Clange("1", m, n, cf, m, rwork);
    if (cnorm > zero) {
        result[3 - 1] = resid / (eps * max((INTEGER)1, m) * cnorm);
    } else {
        result[3 - 1] = zero;
    }
    //
    // Copy C into CF again
    //
    Clacpy("Full", m, n, c, m, cf, m);
    //
    // Apply Q to C as QT*C
    //
    Cgemqrt("L", "C", m, n, k, nb, af, m, t, nb, cf, m, work, info);
    //
    // Compute |QT*C - QT*C| / |C|
    //
    Cgemm("C", "N", m, n, m, -one, q, m, c, m, one, cf, m);
    resid = Clange("1", m, n, cf, m, rwork);
    if (cnorm > zero) {
        result[4 - 1] = resid / (eps * max((INTEGER)1, m) * cnorm);
    } else {
        result[4 - 1] = zero;
    }
    //
    // Generate random n-by-m matrix D and a copy DF
    //
    for (j = 1; j <= m; j = j + 1) {
        Clarnv(2, iseed, n, &d[(j - 1) * n]);
    }
    REAL dnorm = Clange("1", n, m, d, n, rwork);
    Clacpy("Full", n, m, d, n, df, n);
    //
    // Apply Q to D as D*Q
    //
    Cgemqrt("R", "N", n, m, k, nb, af, m, t, nb, df, n, work, info);
    //
    // Compute |D*Q - D*Q| / |D|
    //
    Cgemm("N", "N", n, m, m, -one, d, n, q, m, one, df, n);
    resid = Clange("1", n, m, df, n, rwork);
    if (cnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    // Copy D into DF again
    //
    Clacpy("Full", n, m, d, n, df, n);
    //
    // Apply Q to D as D*QT
    //
    Cgemqrt("R", "C", n, m, k, nb, af, m, t, nb, df, n, work, info);
    //
    // Compute |D*QT - D*QT| / |D|
    //
    Cgemm("N", "C", n, m, m, -one, d, n, q, m, one, df, n);
    resid = Clange("1", n, m, df, n, rwork);
    if (cnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    // Deallocate all arrays
    //
}
