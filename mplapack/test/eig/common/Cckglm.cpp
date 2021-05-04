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
#include <mplapack_eig.h>

#include <mplapack_debug.h>

void Cckglm(INTEGER const nn, INTEGER *nval, INTEGER *mval, INTEGER *pval, INTEGER const nmats, INTEGER *iseed, REAL const thresh, INTEGER const nmax, COMPLEX *a, COMPLEX *af, COMPLEX *b, COMPLEX *bf, COMPLEX *x, COMPLEX *work, REAL *rwork, INTEGER const nin, INTEGER const nout, INTEGER &info) {
    common_write write(cmn);
    char path[3];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    bool firstt = false;
    const INTEGER ntypes = 8;
    arr_1d<ntypes, bool> dotype;
    INTEGER lda = 0;
    INTEGER ldb = 0;
    INTEGER lwork = 0;
    INTEGER ik = 0;
    INTEGER m = 0;
    INTEGER p = 0;
    INTEGER n = 0;
    INTEGER imat = 0;
    char type[1];
    INTEGER kla = 0;
    INTEGER kua = 0;
    INTEGER klb = 0;
    INTEGER kub = 0;
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    INTEGER modea = 0;
    INTEGER modeb = 0;
    REAL cndnma = 0.0;
    REAL cndnmb = 0.0;
    char dista[1];
    char distb[1];
    INTEGER iinfo = 0;
    INTEGER i = 0;
    REAL resid = 0.0;
    static const char *format_9999 = "(' ZLATMS in Cckglm INFO = ',i5)";
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Initialize constants.
    //
    path[(3 - 1) * ldpath] = "GLM";
    info = 0;
    nrun = 0;
    nfail = 0;
    firstt = true;
    Alareq(path, nmats, dotype, ntypes, nin, nout);
    lda = nmax;
    ldb = nmax;
    lwork = nmax * nmax;
    //
    //     Check for valid input values.
    //
    for (ik = 1; ik <= nn; ik = ik + 1) {
        m = mval[ik - 1];
        p = pval[ik - 1];
        n = nval[ik - 1];
        if (m > n || n > m + p) {
            if (firstt) {
                write(nout, star);
                firstt = false;
            }
            write(nout, "(' *** Invalid input  for GLM:  M = ',i6,', P = ',i6,', N = ',i6,';',"
                        "/,'     must satisfy M <= N <= M+P  ',"
                        "'(this set of values will be skipped)')"),
                m, p, n;
        }
    }
    firstt = true;
    //
    //     Do for each value of M in MVAL.
    //
    for (ik = 1; ik <= nn; ik = ik + 1) {
        m = mval[ik - 1];
        p = pval[ik - 1];
        n = nval[ik - 1];
        if (m > n || n > m + p) {
            goto statement_40;
        }
        //
        for (imat = 1; imat <= ntypes; imat = imat + 1) {
            //
            //           Do the tests only if DOTYPE( IMAT ) is true.
            //
            if (!dotype[imat - 1]) {
                goto statement_30;
            }
            //
            //           Set up parameters with Rlatb9 and generate test
            //           matrices A and B with ZLATMS.
            //
            Rlatb9(path, imat, m, p, n, type, kla, kua, klb, kub, anorm, bnorm, modea, modeb, cndnma, cndnmb, dista, distb);
            //
            zlatms(n, m, dista, iseed, type, rwork, modea, cndnma, anorm, kla, kua, "No packing", a, lda, work, iinfo);
            if (iinfo != 0) {
                write(nout, format_9999), iinfo;
                info = abs(iinfo);
                goto statement_30;
            }
            //
            zlatms(n, p, distb, iseed, type, rwork, modeb, cndnmb, bnorm, klb, kub, "No packing", b, ldb, work, iinfo);
            if (iinfo != 0) {
                write(nout, format_9999), iinfo;
                info = abs(iinfo);
                goto statement_30;
            }
            //
            //           Generate random left hand side vector of GLM
            //
            for (i = 1; i <= n; i = i + 1) {
                x[i - 1] = zlarnd(2, iseed);
            }
            //
            Cglmts(n, m, p, a, af, lda, b, bf, ldb, x, &x[(nmax + 1) - 1], &x[(2 * nmax + 1) - 1], &x[(3 * nmax + 1) - 1], work, lwork, rwork, resid);
            //
            //           Print information about the tests that did not
            //           pass the threshold.
            //
            if (resid >= thresh) {
                if (nfail == 0 && firstt) {
                    firstt = false;
                    Alahdg(nout, path);
                }
                write(nout, "(' N=',i4,' M=',i4,', P=',i4,', type ',i2,', test ',i2,', ratio=',"
                            "g13.6)"),
                    n, m, p, imat, 1, resid;
                nfail++;
            }
            nrun++;
        //
        statement_30:;
        }
    statement_40:;
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, 0);
    //
    //     End of Cckglm
    //
}
