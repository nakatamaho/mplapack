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

// Derived from LAPACK routine ZUNHR_COL02.
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

void Cunhr_col02(INTEGER const m, INTEGER const n, INTEGER const mb1, INTEGER const nb1, INTEGER const nb2, REAL *result) {
    static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
    // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
    //
    bool testzeros = false;
    //
    REAL eps = Rlamch("Epsilon");
    INTEGER k = min(m, n);
    INTEGER l = max(m, n, (INTEGER)1);
    //
    // Dynamically allocate local arrays
    //
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
    q_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, l * l));
    q = q_storage.get();
    std::unique_ptr<COMPLEX[]> r_storage;
    COMPLEX *r = nullptr;
    r_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m * l));
    r = r_storage.get();
    std::unique_ptr<REAL[]> rwork_storage;
    REAL *rwork = nullptr;
    rwork_storage = std::make_unique<REAL[]>(max((INTEGER)1, l));
    rwork = rwork_storage.get();
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
    //
    // Put random numbers into A and copy to AF
    //
    INTEGER j = 0;
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &a[(j - 1) * m]);
    }
    if (testzeros) {
        if (m >= 4) {
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(2, iseed, m / 2, &a[((m / 4) - 1) + (j - 1) * m]);
            }
        }
    }
    Clacpy("Full", m, n, a, m, af, m);
    //
    // Number of row blocks in Clatsqr
    //
    INTEGER nrb = max((INTEGER)1, iceil(castREAL(m - n) / castREAL(mb1 - n)));
    //
    std::unique_ptr<COMPLEX[]> t1_storage;
    COMPLEX *t1 = nullptr;
    t1_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, nb1 * (n * nrb)));
    t1 = t1_storage.get();
    std::unique_ptr<COMPLEX[]> t2_storage;
    COMPLEX *t2 = nullptr;
    t2_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, nb2 * n));
    t2 = t2_storage.get();
    std::unique_ptr<COMPLEX[]> diag_storage;
    COMPLEX *diag = nullptr;
    diag_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, n));
    diag = diag_storage.get();
    //
    // Begin determine LWORK for the array WORK and allocate memory.
    //
    // Cgemqrt requires NB2 to be bounded by N.
    //
    INTEGER nb2_ub = min(nb2, n);
    //
    COMPLEX workquery[1];
    INTEGER info = 0;
    Cgetsqrhrt(m, n, mb1, nb1, nb2, af, m, t2, nb2, workquery, -1, info);
    //
    INTEGER lwork = castINTEGER(workquery[1 - 1].real());
    //
    // In Cgemqrt, WORK is N*NB2_UB if SIDE = 'L',
    // or  M*NB2_UB if SIDE = 'R'.
    //
    lwork = max(lwork, nb2_ub * n, nb2_ub * m);
    //
    std::unique_ptr<COMPLEX[]> work_storage;
    COMPLEX *work = nullptr;
    work_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lwork));
    work = work_storage.get();
    //
    // End allocate memory for WORK.
    //
    // Begin Householder reconstruction routines
    //
    // Factor the matrix A in the array AF.
    //
    srnamt = "Cgetsqrhrt";
    Cgetsqrhrt(m, n, mb1, nb1, nb2, af, m, t2, nb2, work, lwork, info);
    //
    // End Householder reconstruction routines.
    //
    // Generate the m-by-m matrix Q
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    Claset("Full", m, m, czero, cone, q, m);
    //
    srnamt = "Cgemqrt";
    Cgemqrt("L", "N", m, m, k, nb2_ub, af, m, t2, nb2, q, m, work, info);
    //
    // Copy R
    //
    Claset("Full", m, n, czero, czero, r, m);
    //
    Clacpy("Upper", m, n, af, m, r, m);
    //
    // TEST 1
    // Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)
    //
    Cgemm("C", "N", m, n, m, -cone, q, m, a, m, cone, r, m);
    //
    REAL anorm = Clange("1", m, n, a, m, rwork);
    REAL resid = Clange("1", m, n, r, m, rwork);
    const REAL zero = 0.0;
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * max((INTEGER)1, m) * anorm);
    } else {
        result[1 - 1] = zero;
    }
    //
    // TEST 2
    // Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)
    //
    Claset("Full", m, m, czero, cone, r, m);
    Cherk("U", "C", m, m, -(cone).real(), q, m, cone.real(), r, m);
    resid = Clansy("1", "Upper", m, r, m, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, m));
    //
    // Generate random m-by-n matrix C
    //
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &c[(j - 1) * m]);
    }
    REAL cnorm = Clange("1", m, n, c, m, rwork);
    Clacpy("Full", m, n, c, m, cf, m);
    //
    // Apply Q to C as Q*C = CF
    //
    srnamt = "Cgemqrt";
    Cgemqrt("L", "N", m, n, k, nb2_ub, af, m, t2, nb2, cf, m, work, info);
    //
    // TEST 3
    // Compute |CF - Q*C| / ( eps *  m * |C| )
    //
    Cgemm("N", "N", m, n, m, -cone, q, m, c, m, cone, cf, m);
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
    // Apply Q to C as (Q**T)*C = CF
    //
    srnamt = "Cgemqrt";
    Cgemqrt("L", "C", m, n, k, nb2_ub, af, m, t2, nb2, cf, m, work, info);
    //
    // TEST 4
    // Compute |CF - (Q**T)*C| / ( eps * m * |C|)
    //
    Cgemm("C", "N", m, n, m, -cone, q, m, c, m, cone, cf, m);
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
    // Apply Q to D as D*Q = DF
    //
    srnamt = "Cgemqrt";
    Cgemqrt("R", "N", n, m, k, nb2_ub, af, m, t2, nb2, df, n, work, info);
    //
    // TEST 5
    // Compute |DF - D*Q| / ( eps * m * |D| )
    //
    Cgemm("N", "N", n, m, m, -cone, d, n, q, m, cone, df, n);
    resid = Clange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[5 - 1] = zero;
    }
    //
    // Copy D into DF again
    //
    Clacpy("Full", n, m, d, n, df, n);
    //
    // Apply Q to D as D*QT = DF
    //
    srnamt = "Cgemqrt";
    Cgemqrt("R", "C", n, m, k, nb2_ub, af, m, t2, nb2, df, n, work, info);
    //
    // TEST 6
    // Compute |DF - D*(Q**T)| / ( eps * m * |D| )
    //
    Cgemm("N", "C", n, m, m, -cone, d, n, q, m, cone, df, n);
    resid = Clange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    // Deallocate all arrays
    //
    // End of Cunhr_col02
    //
}
