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

void Rtsqr01(const char *tssw, INTEGER const m, INTEGER const n, INTEGER const mb, INTEGER const nb, REAL *result) {
    FEM_CMN_SVE(Rtsqr01);
    result([6]);
    // COMMON srnamc
    //
    // SAVE
    INTEGER *iseed(sve.iseed, [4]);
    //
    if (is_called_first_time) {
        static const INTEGER values[] = {1988, 1989, 1990, 1991};
        data_of_type<int>(FEM_VALUES_AND_SIZE), iseed;
    }
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
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Data statements ..
    //
    //     TEST TALL SKINNY OR SHORT WIDE
    //
    bool ts = Mlsame(tssw, "TS");
    //
    //     TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
    //
    bool testzeros = false;
    //
    REAL eps = Rlamch("Epsilon");
    INTEGER k = min(m, n);
    INTEGER l = max({m, n, 1});
    INTEGER mnb = max(mb, nb);
    INTEGER lwork = max(3, l) * mnb;
    //
    //     Dynamically allocate local arrays
    //
    //     Put random numbers into A and copy to AF
    //
    INTEGER j = 0;
    REAL a[m * n];
    for (j = 1; j <= n; j = j + 1) {
        Rlarnv(2, iseed, m, &a[(j - 1) * lda]);
    }
    if (testzeros) {
        if (m >= 4) {
            for (j = 1; j <= n; j = j + 1) {
                Rlarnv(2, iseed, m / 2, &a[((m / 4) - 1) + (j - 1) * lda]);
            }
        }
    }
    REAL af[m * n];
    Rlacpy("Full", m, n, a, m, af, m);
    //
    REAL tquery[5];
    REAL workquery[1];
    INTEGER info = 0;
    INTEGER tsize = 0;
    REAL cf[m * n];
    REAL df[n * m];
    REAL t[tsize];
    REAL work[lwork];
    const REAL zero = 0.0f;
    const REAL one = 1.0f;
    REAL q[l * l];
    REAL r[m * l];
    REAL rwork[l];
    REAL anorm = 0.0;
    REAL resid = 0.0;
    REAL c[m * n];
    REAL cnorm = 0.0;
    REAL d[n * m];
    REAL dnorm = 0.0;
    REAL lq[l * n];
    if (ts) {
        //
        //     Factor the matrix A in the array AF.
        //
        Rgeqr(m, n, af, m, tquery, -1, workquery, -1, info);
        tsize = int(tquery[1 - 1]);
        lwork = int(workquery[1 - 1]);
        Rgemqr("L", "N", m, m, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Rgemqr("L", "N", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Rgemqr("L", "T", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Rgemqr("R", "N", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Rgemqr("R", "T", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Rgeqr(m, n, af, m, t, tsize, work, lwork, info);
        //
        //     Generate the m-by-m matrix Q
        //
        Rlaset("Full", m, m, zero, one, q, m);
        Rgemqr("L", "N", m, m, k, af, m, t, tsize, q, m, work, lwork, info);
        //
        //     Copy R
        //
        Rlaset("Full", m, n, zero, zero, r, m);
        Rlacpy("Upper", m, n, af, m, r, m);
        //
        //     Compute |R - Q'*A| / |A| and store in RESULT(1)
        //
        Rgemm("T", "N", m, n, m, -one, q, m, a, m, one, r, m);
        anorm = Rlange("1", m, n, a, m, rwork);
        resid = Rlange("1", m, n, r, m, rwork);
        if (anorm > zero) {
            result[1 - 1] = resid / (eps * max((INTEGER)1, m) * anorm);
        } else {
            result[1 - 1] = zero;
        }
        //
        //     Compute |I - Q'*Q| and store in RESULT(2)
        //
        Rlaset("Full", m, m, zero, one, r, m);
        Rsyrk("U", "C", m, m, -one, q, m, one, r, m);
        resid = Rlansy("1", "Upper", m, r, m, rwork);
        result[2 - 1] = resid / (eps * max((INTEGER)1, m));
        //
        //     Generate random m-by-n matrix C and a copy CF
        //
        for (j = 1; j <= n; j = j + 1) {
            Rlarnv(2, iseed, m, &c[(j - 1) * ldc]);
        }
        cnorm = Rlange("1", m, n, c, m, rwork);
        Rlacpy("Full", m, n, c, m, cf, m);
        //
        //     Apply Q to C as Q*C
        //
        Rgemqr("L", "N", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
        //
        //     Compute |Q*C - Q*C| / |C|
        //
        Rgemm("N", "N", m, n, m, -one, q, m, c, m, one, cf, m);
        resid = Rlange("1", m, n, cf, m, rwork);
        if (cnorm > zero) {
            result[3 - 1] = resid / (eps * max((INTEGER)1, m) * cnorm);
        } else {
            result[3 - 1] = zero;
        }
        //
        //     Copy C into CF again
        //
        Rlacpy("Full", m, n, c, m, cf, m);
        //
        //     Apply Q to C as QT*C
        //
        Rgemqr("L", "T", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
        //
        //     Compute |QT*C - QT*C| / |C|
        //
        Rgemm("T", "N", m, n, m, -one, q, m, c, m, one, cf, m);
        resid = Rlange("1", m, n, cf, m, rwork);
        if (cnorm > zero) {
            result[4 - 1] = resid / (eps * max((INTEGER)1, m) * cnorm);
        } else {
            result[4 - 1] = zero;
        }
        //
        //     Generate random n-by-m matrix D and a copy DF
        //
        for (j = 1; j <= m; j = j + 1) {
            Rlarnv(2, iseed, n, &d[(j - 1) * ldd]);
        }
        dnorm = Rlange("1", n, m, d, n, rwork);
        Rlacpy("Full", n, m, d, n, df, n);
        //
        //     Apply Q to D as D*Q
        //
        Rgemqr("R", "N", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        //     Compute |D*Q - D*Q| / |D|
        //
        Rgemm("N", "N", n, m, m, -one, d, n, q, m, one, df, n);
        resid = Rlange("1", n, m, df, n, rwork);
        if (dnorm > zero) {
            result[5 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
        } else {
            result[5 - 1] = zero;
        }
        //
        //     Copy D into DF again
        //
        Rlacpy("Full", n, m, d, n, df, n);
        //
        //     Apply Q to D as D*QT
        //
        Rgemqr("R", "T", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        //     Compute |D*QT - D*QT| / |D|
        //
        Rgemm("N", "T", n, m, m, -one, d, n, q, m, one, df, n);
        resid = Rlange("1", n, m, df, n, rwork);
        if (cnorm > zero) {
            result[6 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
        } else {
            result[6 - 1] = zero;
        }
        //
        //     Short and wide
        //
    } else {
        Rgelq(m, n, af, m, tquery, -1, workquery, -1, info);
        tsize = int(tquery[1 - 1]);
        lwork = int(workquery[1 - 1]);
        Rgemlq("R", "N", n, n, k, af, m, tquery, tsize, q, n, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Rgemlq("L", "N", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Rgemlq("L", "T", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Rgemlq("R", "N", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Rgemlq("R", "T", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        FEM_THROW_UNHANDLED("executable allocate: allocate(t(tsize))");
        FEM_THROW_UNHANDLED("executable allocate: allocate(work(lwork))");
        Rgelq(m, n, af, m, t, tsize, work, lwork, info);
        //
        //     Generate the n-by-n matrix Q
        //
        Rlaset("Full", n, n, zero, one, q, n);
        Rgemlq("R", "N", n, n, k, af, m, t, tsize, q, n, work, lwork, info);
        //
        //     Copy R
        //
        Rlaset("Full", m, n, zero, zero, lq, l);
        Rlacpy("Lower", m, n, af, m, lq, l);
        //
        //     Compute |L - A*Q'| / |A| and store in RESULT(1)
        //
        Rgemm("N", "T", m, n, n, -one, a, m, q, n, one, lq, l);
        anorm = Rlange("1", m, n, a, m, rwork);
        resid = Rlange("1", m, n, lq, l, rwork);
        if (anorm > zero) {
            result[1 - 1] = resid / (eps * max((INTEGER)1, n) * anorm);
        } else {
            result[1 - 1] = zero;
        }
        //
        //     Compute |I - Q'*Q| and store in RESULT(2)
        //
        Rlaset("Full", n, n, zero, one, lq, l);
        Rsyrk("U", "C", n, n, -one, q, n, one, lq, l);
        resid = Rlansy("1", "Upper", n, lq, l, rwork);
        result[2 - 1] = resid / (eps * max((INTEGER)1, n));
        //
        //     Generate random m-by-n matrix C and a copy CF
        //
        for (j = 1; j <= m; j = j + 1) {
            Rlarnv(2, iseed, n, &d[(j - 1) * ldd]);
        }
        dnorm = Rlange("1", n, m, d, n, rwork);
        Rlacpy("Full", n, m, d, n, df, n);
        //
        //     Apply Q to C as Q*C
        //
        Rgemlq("L", "N", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        //     Compute |Q*D - Q*D| / |D|
        //
        Rgemm("N", "N", n, m, n, -one, q, n, d, n, one, df, n);
        resid = Rlange("1", n, m, df, n, rwork);
        if (dnorm > zero) {
            result[3 - 1] = resid / (eps * max((INTEGER)1, n) * dnorm);
        } else {
            result[3 - 1] = zero;
        }
        //
        //     Copy D into DF again
        //
        Rlacpy("Full", n, m, d, n, df, n);
        //
        //     Apply Q to D as QT*D
        //
        Rgemlq("L", "T", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        //     Compute |QT*D - QT*D| / |D|
        //
        Rgemm("T", "N", n, m, n, -one, q, n, d, n, one, df, n);
        resid = Rlange("1", n, m, df, n, rwork);
        if (dnorm > zero) {
            result[4 - 1] = resid / (eps * max((INTEGER)1, n) * dnorm);
        } else {
            result[4 - 1] = zero;
        }
        //
        //     Generate random n-by-m matrix D and a copy DF
        //
        for (j = 1; j <= n; j = j + 1) {
            Rlarnv(2, iseed, m, &c[(j - 1) * ldc]);
        }
        cnorm = Rlange("1", m, n, c, m, rwork);
        Rlacpy("Full", m, n, c, m, cf, m);
        //
        //     Apply Q to C as C*Q
        //
        Rgemlq("R", "N", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
        //
        //     Compute |C*Q - C*Q| / |C|
        //
        Rgemm("N", "N", m, n, n, -one, c, m, q, n, one, cf, m);
        resid = Rlange("1", n, m, df, n, rwork);
        if (cnorm > zero) {
            result[5 - 1] = resid / (eps * max((INTEGER)1, n) * cnorm);
        } else {
            result[5 - 1] = zero;
        }
        //
        //     Copy C into CF again
        //
        Rlacpy("Full", m, n, c, m, cf, m);
        //
        //     Apply Q to D as D*QT
        //
        Rgemlq("R", "T", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
        //
        //     Compute |C*QT - C*QT| / |C|
        //
        Rgemm("N", "T", m, n, n, -one, c, m, q, n, one, cf, m);
        resid = Rlange("1", m, n, cf, m, rwork);
        if (cnorm > zero) {
            result[6 - 1] = resid / (eps * max((INTEGER)1, n) * cnorm);
        } else {
            result[6 - 1] = zero;
        }
        //
    }
    //
    //     Deallocate all arrays
    //
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t,c,d,cf,df)");
    //
}
