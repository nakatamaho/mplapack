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

// Derived from LAPACK routine ZLQT04.
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

void Clqt04(INTEGER const m, INTEGER const n, INTEGER const nb, REAL *result) {
    common cmn;
    static INTEGER iseed[4] = {1988, 1989, 1990, 1991};
    REAL eps = Rlamch("Epsilon");
    INTEGER k = min(m, n);
    INTEGER ll = max(m, n);
    INTEGER lwork = max((INTEGER)2, ll) * max((INTEGER)2, ll) * nb;
    INTEGER ldt = nb;
    INTEGER j = 0;
    std::unique_ptr<COMPLEX[]> __a_storage(new COMPLEX[m * n]);
    COMPLEX *a = __a_storage.get();
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &a[(j - 1) * m]);
    }
    std::unique_ptr<COMPLEX[]> __af_storage(new COMPLEX[m * n]);
    COMPLEX *af = __af_storage.get();
    Clacpy("Full", m, n, a, m, af, m);
    std::unique_ptr<COMPLEX[]> __t_storage(new COMPLEX[nb * n]);
    COMPLEX *t = __t_storage.get();
    std::unique_ptr<COMPLEX[]> __work_storage(new COMPLEX[lwork]);
    COMPLEX *work = __work_storage.get();
    INTEGER info = 0;
    Cgelqt(m, n, nb, af, m, t, ldt, work, info);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const COMPLEX one = COMPLEX(1.0, 0.0);
    std::unique_ptr<COMPLEX[]> __q_storage(new COMPLEX[n * n]);
    COMPLEX *q = __q_storage.get();
    Claset("Full", n, n, czero, one, q, n);
    Cgemlqt("R", "N", n, n, k, nb, af, m, t, ldt, q, n, work, info);
    std::unique_ptr<COMPLEX[]> __l_storage(new COMPLEX[ll * n]);
    COMPLEX *l = __l_storage.get();
    Claset("Full", ll, n, czero, czero, l, ll);
    Clacpy("Lower", m, n, af, m, l, ll);
    // Compute |L - A*Q'| / |A| and store in RESULT(1)
    //
    Cgemm("N", "C", m, n, n, -one, a, m, q, n, one, l, ll);
    std::unique_ptr<COMPLEX[]> __rwork_storage(new COMPLEX[ll]);
    COMPLEX *rwork = __rwork_storage.get();
    REAL anorm = Clange("1", m, n, a, m, rwork);
    REAL resid = Clange("1", m, n, l, ll, rwork);
    const REAL zero = 0.0;
    if (anorm > zero) {
        result[1 - 1] = resid / (eps * max((INTEGER)1, m) * anorm);
    } else {
        result[1 - 1] = zero;
    }
    //
    // Compute |I - Q'*Q| and store in RESULT(2)
    //
    Claset("Full", n, n, czero, one, l, ll);
    Cherk("U", "C", n, n, dreal(-one), q, n, dreal(one), l, ll);
    resid = Clansy("1", "Upper", n, l, ll, rwork);
    result[2 - 1] = resid / (eps * max((INTEGER)1, n));
    //
    // Generate random m-by-n matrix C and a copy CF
    //
    std::unique_ptr<COMPLEX[]> __d_storage(new COMPLEX[n * m]);
    COMPLEX *d = __d_storage.get();
    for (j = 1; j <= m; j = j + 1) {
        Clarnv(2, iseed, n, &d[(j - 1) * n]);
    }
    REAL dnorm = Clange("1", n, m, d, n, rwork);
    std::unique_ptr<COMPLEX[]> __df_storage(new COMPLEX[n * m]);
    COMPLEX *df = __df_storage.get();
    Clacpy("Full", n, m, d, n, df, n);
    //
    // Apply Q to C as Q*C
    //
    Cgemlqt("L", "N", n, m, k, nb, af, m, t, nb, df, n, work, info);
    //
    // Compute |Q*D - Q*D| / |D|
    //
    Cgemm("N", "N", n, m, n, -one, q, n, d, n, one, df, n);
    resid = Clange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[3 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
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
    Cgemlqt("L", "C", n, m, k, nb, af, m, t, nb, df, n, work, info);
    //
    // Compute |QT*D - QT*D| / |D|
    //
    Cgemm("C", "N", n, m, n, -one, q, n, d, n, one, df, n);
    resid = Clange("1", n, m, df, n, rwork);
    if (dnorm > zero) {
        result[4 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[4 - 1] = zero;
    }
    //
    // Generate random n-by-m matrix D and a copy DF
    //
    std::unique_ptr<COMPLEX[]> __c_storage(new COMPLEX[m * n]);
    COMPLEX *c = __c_storage.get();
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &c[(j - 1) * m]);
    }
    REAL cnorm = Clange("1", m, n, c, m, rwork);
    std::unique_ptr<COMPLEX[]> __cf_storage(new COMPLEX[m * n]);
    COMPLEX *cf = __cf_storage.get();
    Clacpy("Full", m, n, c, m, cf, m);
    //
    // Apply Q to C as C*Q
    //
    Cgemlqt("R", "N", m, n, k, nb, af, m, t, nb, cf, m, work, info);
    //
    // Compute |C*Q - C*Q| / |C|
    //
    Cgemm("N", "N", m, n, n, -one, c, m, q, n, one, cf, m);
    resid = Clange("1", n, m, df, n, rwork);
    if (cnorm > zero) {
        result[5 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
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
    Cgemlqt("R", "C", m, n, k, nb, af, m, t, nb, cf, m, work, info);
    //
    // Compute |C*QT - C*QT| / |C|
    //
    Cgemm("N", "C", m, n, n, -one, c, m, q, n, one, cf, m);
    resid = Clange("1", m, n, cf, m, rwork);
    if (cnorm > zero) {
        result[6 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
    } else {
        result[6 - 1] = zero;
    }
    //
    // Deallocate all arrays
    //
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,l,rwork,work,t,c,d,cf,df)");
    //
}
