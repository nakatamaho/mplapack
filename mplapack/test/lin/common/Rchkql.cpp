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

#include <mplapack_debug.h>

void Rchkql(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nn, INTEGER *nval, INTEGER const nnb, INTEGER *nbval, INTEGER *nxval, INTEGER const nrhs, REAL const thresh, bool const tsterr, INTEGER const nmax, REAL *a, REAL *af, REAL *aq, REAL *al, REAL *ac, REAL *b, REAL *x, REAL *xact, REAL *tau, REAL *work, REAL *rwork, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    //
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    char path[3];
    char buf[1024];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    INTEGER lda = 0;
    INTEGER lwork = 0;
    INTEGER im = 0;
    INTEGER m = 0;
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER minmn = 0;
    INTEGER imat = 0;
    const INTEGER ntypes = 8;
    char type;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist;
    INTEGER info = 0;
    INTEGER kval[4];
    INTEGER nk = 0;
    INTEGER ik = 0;
    INTEGER k = 0;
    INTEGER inb = 0;
    INTEGER nb = 0;
    INTEGER nx = 0;
    const INTEGER ntests = 7;
    const REAL zero = 0.0;
    REAL result[ntests];
    INTEGER nt = 0;
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Initialize constants and the random number seed.
    //
    path[0] = 'R';
    path[1] = 'Q';
    path[2] = 'L';
    nrun = 0;
    nfail = 0;
    nerrs = 0;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    //
    //     Test the error exits
    //
    if (tsterr) {
        Rerrql(path, nout);
    }
    infot = 0;
    xlaenv(2, 2);
    //
    lda = nmax;
    lwork = nmax * max(nmax, nrhs);
    //
    //     Do for each value of M in MVAL.
    //
    for (im = 1; im <= nm; im = im + 1) {
        m = mval[im - 1];
        //
        //        Do for each value of N in NVAL.
        //
        for (in = 1; in <= nn; in = in + 1) {
            n = nval[in - 1];
            minmn = min(m, n);
            for (imat = 1; imat <= ntypes; imat = imat + 1) {
                //
                //              Do the tests only if DOTYPE( IMAT ) is true.
                //
                if (!dotype[imat - 1]) {
                    goto statement_50;
                }
                //
                //              Set up parameters with Rlatb4 and generate a test matrix
                //              with Rlatms.
                //
                Rlatb4(path, imat, m, n, &type, kl, ku, anorm, mode, cndnum, &dist);
                //
                strncpy(srnamt, "Rlatms", srnamt_len);
                Rlatms(m, n, &dist, iseed, &type, rwork, mode, cndnum, anorm, kl, ku, "No packing", a, lda, work, info);
                //
                //              Check error code from Rlatms.
                //
                if (info != 0) {
                    Alaerh(path, "Rlatms", info, 0, " ", m, n, -1, -1, -1, imat, nfail, nerrs, nout);
                    goto statement_50;
                }
                //
                //              Set some values for K: the first value must be MINMN,
                //              corresponding to the call of Rqlt01; other values are
                //              used in the calls of Rqlt02, and must not exceed MINMN.
                //
                kval[1 - 1] = minmn;
                kval[2 - 1] = 0;
                kval[3 - 1] = 1;
                kval[4 - 1] = minmn / 2;
                if (minmn == 0) {
                    nk = 1;
                } else if (minmn == 1) {
                    nk = 2;
                } else if (minmn <= 3) {
                    nk = 3;
                } else {
                    nk = 4;
                }
                //
                //              Do for each value of K in KVAL
                //
                for (ik = 1; ik <= nk; ik = ik + 1) {
                    k = kval[ik - 1];
                    //
                    //                 Do for each pair of values (NB,NX) in NBVAL and NXVAL.
                    //
                    for (inb = 1; inb <= nnb; inb = inb + 1) {
                        nb = nbval[inb - 1];
                        xlaenv(1, nb);
                        nx = nxval[inb - 1];
                        xlaenv(3, nx);
                        for (i = 1; i <= ntests; i = i + 1) {
                            result[i - 1] = zero;
                        }
                        nt = 2;
                        if (ik == 1) {
                            //
                            //                       Test Rgeqlf
                            //
                            Rqlt01(m, n, a, af, aq, al, lda, tau, work, lwork, rwork, &result[1 - 1]);
                        } else if (m >= n) {
                            //
                            //                       Test Rorgql, using factorization
                            //                       returned by Rqlt01
                            //
                            Rqlt02(m, n, k, a, af, aq, al, lda, tau, work, lwork, rwork, &result[1 - 1]);
                        }
                        if (m >= k) {
                            //
                            //                       Test Rormql, using factorization returned
                            //                       by Rqlt01
                            //
                            Rqlt03(m, n, k, af, ac, al, aq, lda, tau, work, lwork, rwork, &result[3 - 1]);
                            nt += 4;
                            //
                            //                       If M>=N and K=N, call Rgeqls to solve a system
                            //                       with NRHS right hand sides and compute the
                            //                       residual.
                            //
                            if (k == n && inb == 1) {
                                //
                                //                          Generate a solution and set the right
                                //                          hand side.
                                //
                                strncpy(srnamt, "Rlarhs", srnamt_len); 
                                Rlarhs(path, "New", "Full", "No transpose", m, n, 0, 0, nrhs, a, lda, xact, lda, b, lda, iseed, info);
                                //
                                Rlacpy("Full", m, nrhs, b, lda, x, lda);
                                strncpy(srnamt, "Rgeqls", srnamt_len);  
                                Rgeqls(m, n, nrhs, af, lda, tau, x, lda, work, lwork, info);
                                //
                                //                          Check error code from Rgeqls.
                                //
                                if (info != 0) {
                                    Alaerh(path, "Rgeqls", info, 0, " ", m, n, nrhs, -1, nb, imat, nfail, nerrs, nout);
                                }
                                //
                                Rget02("No transpose", m, n, nrhs, a, lda, &x[(m - n + 1) - 1], lda, b, lda, rwork, result[7 - 1]);
                                nt++;
                            }
                        }
                        //
                        //                    Print information about the tests that did not
                        //                    pass the threshold.
                        //
                        for (i = 1; i <= nt; i = i + 1) {
                            if (result[i - 1] >= thresh) {
                                if (nfail == 0 && nerrs == 0) {
                                    Alahd(nout, path);
                                }
                                sprintnum_short(buf, result[i - 1]);
                                write(nout, "(' M=',i5,', N=',i5,', K=',i5,', NB=',i4,', NX=',i5,"
                                            "', type ',i2,', test(',i2,')=',a)"),
                                    m, n, k, nb, nx, imat, i, buf;
                                nfail++;
                            }
                        }
                        nrun += nt;
                    }
                }
            statement_50:;
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Rchkql
    //
}
