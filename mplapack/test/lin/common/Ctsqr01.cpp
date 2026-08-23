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

// Derived from LAPACK routine ZTSQR01.
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

void Ctsqr01(fem::str_cref tssw, INTEGER const m, INTEGER const n, INTEGER const mb, INTEGER const nb, REAL *result) {
    static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
    // TEST TALL SKINNY OR SHORT WIDE
    //
    bool ts = Mlsame(tssw.elems(), "TS");
    //
    // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
    //
    bool testzeros = false;
    //
    REAL eps = Rlamch("Epsilon");
    INTEGER k = min(m, n);
    INTEGER l = max(m, n, (INTEGER)1);
    INTEGER mnb = max(mb, nb);
    INTEGER lwork = max((INTEGER)3, l) * mnb;
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
    std::unique_ptr<COMPLEX[]> lq_storage;
    COMPLEX *lq = nullptr;
    lq_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, l * n));
    lq = lq_storage.get();
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
    COMPLEX tquery[5];
    COMPLEX workquery[1];
    INTEGER info = 0;
    INTEGER tsize = 0;
    std::unique_ptr<COMPLEX[]> t_storage;
    COMPLEX *t = nullptr;
    std::unique_ptr<COMPLEX[]> work_storage;
    COMPLEX *work = nullptr;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX one = COMPLEX(1.0, 0.0);
    REAL anorm = 0.0;
    REAL resid = 0.0;
    const REAL zero = 0.0;
    REAL cnorm = 0.0;
    REAL dnorm = 0.0;
    if (ts) {
        //
        // Factor the matrix A in the array AF.
        //
        Cgeqr(m, n, af, m, tquery, -1, workquery, -1, info);
        tsize = castINTEGER(tquery[1 - 1].real());
        lwork = castINTEGER(workquery[1 - 1].real());
        Cgemqr("L", "N", m, m, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
        Cgemqr("L", "N", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
        Cgemqr("L", "C", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
        Cgemqr("R", "N", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
        Cgemqr("R", "C", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
        t_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, tsize));
        t = t_storage.get();
        work_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lwork));
        work = work_storage.get();
        srnamt = "Cgeqr";
        Cgeqr(m, n, af, m, t, tsize, work, lwork, info);
        //
        // Generate the m-by-m matrix Q
        //
        Claset("Full", m, m, czero, one, q, m);
        srnamt = "Cgemqr";
        Cgemqr("L", "N", m, m, k, af, m, t, tsize, q, m, work, lwork, info);
        //
        // Copy R
        //
        Claset("Full", m, n, czero, czero, r, m);
        Clacpy("Upper", m, n, af, m, r, m);
        //
        // Compute |R - Q'*A| / |A| and store in RESULT(1)
        //
        Cgemm("C", "N", m, n, m, -one, q, m, a, m, one, r, m);
        anorm = Clange("1", m, n, a, m, rwork);
        resid = Clange("1", m, n, r, m, rwork);
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
        cnorm = Clange("1", m, n, c, m, rwork);
        Clacpy("Full", m, n, c, m, cf, m);
        //
        // Apply Q to C as Q*C
        //
        srnamt = "Cgemqr";
        Cgemqr("L", "N", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
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
        srnamt = "Cgemqr";
        Cgemqr("L", "C", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
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
        dnorm = Clange("1", n, m, d, n, rwork);
        Clacpy("Full", n, m, d, n, df, n);
        //
        // Apply Q to D as D*Q
        //
        srnamt = "Cgemqr";
        Cgemqr("R", "N", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        // Compute |D*Q - D*Q| / |D|
        //
        Cgemm("N", "N", n, m, m, -one, d, n, q, m, one, df, n);
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
        // Apply Q to D as D*QT
        //
        Cgemqr("R", "C", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
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
        // Short and wide
        //
    } else {
        Cgelq(m, n, af, m, tquery, -1, workquery, -1, info);
        tsize = castINTEGER(tquery[1 - 1].real());
        lwork = castINTEGER(workquery[1 - 1].real());
        Cgemlq("R", "N", n, n, k, af, m, tquery, tsize, q, n, workquery, -1, info);
        lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
        Cgemlq("L", "N", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
        Cgemlq("L", "C", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
        Cgemlq("R", "N", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
        Cgemlq("R", "C", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, castINTEGER(workquery[1 - 1].real()));
        t_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, tsize));
        t = t_storage.get();
        work_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lwork));
        work = work_storage.get();
        srnamt = "Cgelq";
        Cgelq(m, n, af, m, t, tsize, work, lwork, info);
        //
        // Generate the n-by-n matrix Q
        //
        Claset("Full", n, n, czero, one, q, n);
        srnamt = "Cgemlq";
        Cgemlq("R", "N", n, n, k, af, m, t, tsize, q, n, work, lwork, info);
        //
        // Copy R
        //
        Claset("Full", m, n, czero, czero, lq, l);
        Clacpy("Lower", m, n, af, m, lq, l);
        //
        // Compute |L - A*Q'| / |A| and store in RESULT(1)
        //
        Cgemm("N", "C", m, n, n, -one, a, m, q, n, one, lq, l);
        anorm = Clange("1", m, n, a, m, rwork);
        resid = Clange("1", m, n, lq, l, rwork);
        if (anorm > zero) {
            result[1 - 1] = resid / (eps * max((INTEGER)1, n) * anorm);
        } else {
            result[1 - 1] = zero;
        }
        //
        // Compute |I - Q'*Q| and store in RESULT(2)
        //
        Claset("Full", n, n, czero, one, lq, l);
        Cherk("U", "C", n, n, (-one).real(), q, n, one.real(), lq, l);
        resid = Clansy("1", "Upper", n, lq, l, rwork);
        result[2 - 1] = resid / (eps * max((INTEGER)1, n));
        //
        // Generate random m-by-n matrix C and a copy CF
        //
        for (j = 1; j <= m; j = j + 1) {
            Clarnv(2, iseed, n, &d[(j - 1) * n]);
        }
        dnorm = Clange("1", n, m, d, n, rwork);
        Clacpy("Full", n, m, d, n, df, n);
        //
        // Apply Q to C as Q*C
        //
        Cgemlq("L", "N", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        // Compute |Q*D - Q*D| / |D|
        //
        Cgemm("N", "N", n, m, n, -one, q, n, d, n, one, df, n);
        resid = Clange("1", n, m, df, n, rwork);
        if (dnorm > zero) {
            result[3 - 1] = resid / (eps * max((INTEGER)1, n) * dnorm);
        } else {
            result[3 - 1] = zero;
        }
        //
        // Copy D into DF again
        //
        Clacpy("Full", n, m, d, n, df, n);
        //
        // Apply Q to D as QT*D
        //
        Cgemlq("L", "C", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        // Compute |QT*D - QT*D| / |D|
        //
        Cgemm("C", "N", n, m, n, -one, q, n, d, n, one, df, n);
        resid = Clange("1", n, m, df, n, rwork);
        if (dnorm > zero) {
            result[4 - 1] = resid / (eps * max((INTEGER)1, n) * dnorm);
        } else {
            result[4 - 1] = zero;
        }
        //
        // Generate random n-by-m matrix D and a copy DF
        //
        for (j = 1; j <= n; j = j + 1) {
            Clarnv(2, iseed, m, &c[(j - 1) * m]);
        }
        cnorm = Clange("1", m, n, c, m, rwork);
        Clacpy("Full", m, n, c, m, cf, m);
        //
        // Apply Q to C as C*Q
        //
        Cgemlq("R", "N", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
        //
        // Compute |C*Q - C*Q| / |C|
        //
        Cgemm("N", "N", m, n, n, -one, c, m, q, n, one, cf, m);
        resid = Clange("1", m, n, cf, m, rwork);
        if (cnorm > zero) {
            result[5 - 1] = resid / (eps * max((INTEGER)1, n) * cnorm);
        } else {
            result[5 - 1] = zero;
        }
        //
        // Copy C into CF again
        //
        Clacpy("Full", m, n, c, m, cf, m);
        //
        // Apply Q to D as D*QT
        //
        Cgemlq("R", "C", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
        //
        // Compute |C*QT - C*QT| / |C|
        //
        Cgemm("N", "C", m, n, n, -one, c, m, q, n, one, cf, m);
        resid = Clange("1", m, n, cf, m, rwork);
        if (cnorm > zero) {
            result[6 - 1] = resid / (eps * max((INTEGER)1, n) * cnorm);
        } else {
            result[6 - 1] = zero;
        }
        //
    }
    //
    // Deallocate all arrays
    //
}
