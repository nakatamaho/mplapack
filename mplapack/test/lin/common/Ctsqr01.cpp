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

void Ctsqr01(const char *tssw, INTEGER &m, INTEGER &n, INTEGER const mb, INTEGER const nb, REAL *result) {
    FEM_CMN_SVE(Ctsqr01);
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
    m = 10;
    n = 10;
    INTEGER l = 10;
    INTEGER lwork = 10;
    //
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
    l = max({m, n, 1});
    INTEGER mnb = max(mb, nb);
    lwork = max(3, l) * mnb;
    //
    //     Dynamically allocate local arrays
    //
    //     Put random numbers into A and copy to AF
    //
    INTEGER j = 0;
    COMPLEX a[m * n];
    for (j = 1; j <= n; j = j + 1) {
        Clarnv(2, iseed, m, &a[(j - 1) * lda]);
    }
    if (testzeros) {
        if (m >= 4) {
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(2, iseed, m / 2, &a[((m / 4) - 1) + (j - 1) * lda]);
            }
        }
    }
    COMPLEX af[m * n];
    Clacpy("Full", m, n, a, m, af, m);
    //
    COMPLEX tquery[5];
    COMPLEX workquery[1];
    INTEGER info = 0;
    INTEGER tsize = 0;
    COMPLEX cf[m * n];
    COMPLEX df[n * m];
    COMPLEX t[m * n];
    COMPLEX work[lwork];
    const COMPLEX czero = COMPLEX(0.0f, 0.0f);
    const COMPLEX one = COMPLEX(1.0f, 0.0f);
    COMPLEX q[l * l];
    COMPLEX r[m * l];
    COMPLEX rwork[l];
    REAL anorm = 0.0;
    REAL resid = 0.0;
    const REAL zero = 0.0f;
    COMPLEX c[m * n];
    REAL cnorm = 0.0;
    COMPLEX d[n * m];
    REAL dnorm = 0.0;
    COMPLEX lq[l * n];
    if (ts) {
        //
        //     Factor the matrix A in the array AF.
        //
        Cgeqr(m, n, af, m, tquery, -1, workquery, -1, info);
        tsize = int(tquery[1 - 1]);
        lwork = int(workquery[1 - 1]);
        Cgemqr("L", "N", m, m, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Cgemqr("L", "N", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Cgemqr("L", "C", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Cgemqr("R", "N", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Cgemqr("R", "C", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Cgeqr(m, n, af, m, t, tsize, work, lwork, info);
        //
        //     Generate the m-by-m matrix Q
        //
        Claset("Full", m, m, czero, one, q, m);
        Cgemqr("L", "N", m, m, k, af, m, t, tsize, q, m, work, lwork, info);
        //
        //     Copy R
        //
        Claset("Full", m, n, czero, czero, r, m);
        Clacpy("Upper", m, n, af, m, r, m);
        //
        //     Compute |R - Q'*A| / |A| and store in RESULT(1)
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
        //     Compute |I - Q'*Q| and store in RESULT(2)
        //
        Claset("Full", m, m, czero, one, r, m);
        Cherk("U", "C", m, m, dreal[-one - 1], q, m, dreal[one - 1], r, m);
        resid = Clansy("1", "Upper", m, r, m, rwork);
        result[2 - 1] = resid / (eps * max((INTEGER)1, m));
        //
        //     Generate random m-by-n matrix C and a copy CF
        //
        for (j = 1; j <= n; j = j + 1) {
            Clarnv(2, iseed, m, &c[(j - 1) * ldc]);
        }
        cnorm = Clange("1", m, n, c, m, rwork);
        Clacpy("Full", m, n, c, m, cf, m);
        //
        //     Apply Q to C as Q*C
        //
        Cgemqr("L", "N", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
        //
        //     Compute |Q*C - Q*C| / |C|
        //
        Cgemm("N", "N", m, n, m, -one, q, m, c, m, one, cf, m);
        resid = Clange("1", m, n, cf, m, rwork);
        if (cnorm > zero) {
            result[3 - 1] = resid / (eps * max((INTEGER)1, m) * cnorm);
        } else {
            result[3 - 1] = zero;
        }
        //
        //     Copy C into CF again
        //
        Clacpy("Full", m, n, c, m, cf, m);
        //
        //     Apply Q to C as QT*C
        //
        Cgemqr("L", "C", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
        //
        //     Compute |QT*C - QT*C| / |C|
        //
        Cgemm("C", "N", m, n, m, -one, q, m, c, m, one, cf, m);
        resid = Clange("1", m, n, cf, m, rwork);
        if (cnorm > zero) {
            result[4 - 1] = resid / (eps * max((INTEGER)1, m) * cnorm);
        } else {
            result[4 - 1] = zero;
        }
        //
        //     Generate random n-by-m matrix D and a copy DF
        //
        for (j = 1; j <= m; j = j + 1) {
            Clarnv(2, iseed, n, &d[(j - 1) * ldd]);
        }
        dnorm = Clange("1", n, m, d, n, rwork);
        Clacpy("Full", n, m, d, n, df, n);
        //
        //     Apply Q to D as D*Q
        //
        Cgemqr("R", "N", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        //     Compute |D*Q - D*Q| / |D|
        //
        Cgemm("N", "N", n, m, m, -one, d, n, q, m, one, df, n);
        resid = Clange("1", n, m, df, n, rwork);
        if (dnorm > zero) {
            result[5 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
        } else {
            result[5 - 1] = zero;
        }
        //
        //     Copy D into DF again
        //
        Clacpy("Full", n, m, d, n, df, n);
        //
        //     Apply Q to D as D*QT
        //
        Cgemqr("R", "C", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        //     Compute |D*QT - D*QT| / |D|
        //
        Cgemm("N", "C", n, m, m, -one, d, n, q, m, one, df, n);
        resid = Clange("1", n, m, df, n, rwork);
        if (cnorm > zero) {
            result[6 - 1] = resid / (eps * max((INTEGER)1, m) * dnorm);
        } else {
            result[6 - 1] = zero;
        }
        //
        //     Short and wide
        //
    } else {
        Cgelq(m, n, af, m, tquery, -1, workquery, -1, info);
        tsize = int(tquery[1 - 1]);
        lwork = int(workquery[1 - 1]);
        Cgemlq("R", "N", n, n, k, af, m, tquery, tsize, q, n, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Cgemlq("L", "N", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Cgemlq("L", "C", n, m, k, af, m, tquery, tsize, df, n, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Cgemlq("R", "N", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Cgemlq("R", "C", m, n, k, af, m, tquery, tsize, cf, m, workquery, -1, info);
        lwork = max(lwork, int(workquery[1 - 1]));
        Cgelq(m, n, af, m, t, tsize, work, lwork, info);
        //
        //     Generate the n-by-n matrix Q
        //
        Claset("Full", n, n, czero, one, q, n);
        Cgemlq("R", "N", n, n, k, af, m, t, tsize, q, n, work, lwork, info);
        //
        //     Copy R
        //
        Claset("Full", m, n, czero, czero, lq, l);
        Clacpy("Lower", m, n, af, m, lq, l);
        //
        //     Compute |L - A*Q'| / |A| and store in RESULT(1)
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
        //     Compute |I - Q'*Q| and store in RESULT(2)
        //
        Claset("Full", n, n, czero, one, lq, l);
        Cherk("U", "C", n, n, dreal[-one - 1], q, n, dreal[one - 1], lq, l);
        resid = Clansy("1", "Upper", n, lq, l, rwork);
        result[2 - 1] = resid / (eps * max((INTEGER)1, n));
        //
        //     Generate random m-by-n matrix C and a copy CF
        //
        for (j = 1; j <= m; j = j + 1) {
            Clarnv(2, iseed, n, &d[(j - 1) * ldd]);
        }
        dnorm = Clange("1", n, m, d, n, rwork);
        Clacpy("Full", n, m, d, n, df, n);
        //
        //     Apply Q to C as Q*C
        //
        Cgemlq("L", "N", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        //     Compute |Q*D - Q*D| / |D|
        //
        Cgemm("N", "N", n, m, n, -one, q, n, d, n, one, df, n);
        resid = Clange("1", n, m, df, n, rwork);
        if (dnorm > zero) {
            result[3 - 1] = resid / (eps * max((INTEGER)1, n) * dnorm);
        } else {
            result[3 - 1] = zero;
        }
        //
        //     Copy D into DF again
        //
        Clacpy("Full", n, m, d, n, df, n);
        //
        //     Apply Q to D as QT*D
        //
        Cgemlq("L", "C", n, m, k, af, m, t, tsize, df, n, work, lwork, info);
        //
        //     Compute |QT*D - QT*D| / |D|
        //
        Cgemm("C", "N", n, m, n, -one, q, n, d, n, one, df, n);
        resid = Clange("1", n, m, df, n, rwork);
        if (dnorm > zero) {
            result[4 - 1] = resid / (eps * max((INTEGER)1, n) * dnorm);
        } else {
            result[4 - 1] = zero;
        }
        //
        //     Generate random n-by-m matrix D and a copy DF
        //
        for (j = 1; j <= n; j = j + 1) {
            Clarnv(2, iseed, m, &c[(j - 1) * ldc]);
        }
        cnorm = Clange("1", m, n, c, m, rwork);
        Clacpy("Full", m, n, c, m, cf, m);
        //
        //     Apply Q to C as C*Q
        //
        Cgemlq("R", "N", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
        //
        //     Compute |C*Q - C*Q| / |C|
        //
        Cgemm("N", "N", m, n, n, -one, c, m, q, n, one, cf, m);
        resid = Clange("1", n, m, df, n, rwork);
        if (cnorm > zero) {
            result[5 - 1] = resid / (eps * max((INTEGER)1, n) * cnorm);
        } else {
            result[5 - 1] = zero;
        }
        //
        //     Copy C into CF again
        //
        Clacpy("Full", m, n, c, m, cf, m);
        //
        //     Apply Q to D as D*QT
        //
        Cgemlq("R", "C", m, n, k, af, m, t, tsize, cf, m, work, lwork, info);
        //
        //     Compute |C*QT - C*QT| / |C|
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
    //     Deallocate all arrays
    //
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(a,af,q,r,rwork,work,t,c,d,cf,df)");
    //
}
