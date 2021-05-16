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

void Rchkpt(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, REAL const thresh, bool const tsterr, REAL *a, REAL *d, REAL *e, REAL *b, REAL *x, REAL *xact, REAL *work, REAL *rwork, INTEGER const nout) {
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
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    const INTEGER ntypes = 12;
    INTEGER nimat = 0;
    INTEGER imat = 0;
    char type;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cond = 0.0;
    char dist;
    bool zerot = false;
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER ia = 0;
    INTEGER ix = 0;
    REAL dmax = 0.0;
    REAL z[3];
    const REAL zero = 0.0;
    REAL rcondc = 0.0;
    const INTEGER ntests = 7;
    REAL result[ntests];
    REAL ainvnm = 0.0;
    INTEGER j = 0;
    const REAL one = 1.0;
    INTEGER irhs = 0;
    INTEGER nrhs = 0;
    INTEGER k = 0;
    REAL rcond = 0.0;
    static const char *format_9999 = "(' N =',i5,', type ',i2,', test ',i2,', ratio = ',a)";
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
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    path[0] = 'R';
    path[1] = 'P';
    path[2] = 'T';
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
        Rerrgt(path, nout);
    }
    infot = 0;
    //
    for (in = 1; in <= nn; in = in + 1) {
        //
        //        Do for each value of N in NVAL.
        //
        n = nval[in - 1];
        lda = max((INTEGER)1, n);
        nimat = ntypes;
        if (n <= 0) {
            nimat = 1;
        }
        //
        for (imat = 1; imat <= nimat; imat = imat + 1) {
            //
            //           Do the tests only if DOTYPE( IMAT ) is true.
            //
            if (n > 0 && !dotype[imat - 1]) {
                goto statement_100;
            }
            //
            //           Set up parameters with Rlatb4.
            //
            Rlatb4(path, imat, n, n, &type, kl, ku, anorm, mode, cond, &dist);
            //
            zerot = imat >= 8 && imat <= 10;
            if (imat <= 6) {
                //
                //              Type 1-6:  generate a symmetric tridiagonal matrix of
                //              known condition number in lower triangular band storage.
                //
                strncpy(srnamt, "Rlatms", srnamt_len);
                Rlatms(n, n, &dist, iseed, &type, rwork, mode, cond, anorm, kl, ku, "B", a, 2, work, info);
                //
                //              Check the error code from Rlatms.
                //
                if (info != 0) {
                    Alaerh(path, "Rlatms", info, 0, " ", n, n, kl, ku, -1, imat, nfail, nerrs, nout);
                    goto statement_100;
                }
                izero = 0;
                //
                //              Copy the matrix to D and E.
                //
                ia = 1;
                for (i = 1; i <= n - 1; i = i + 1) {
                    d[i - 1] = a[ia - 1];
                    e[i - 1] = a[(ia + 1) - 1];
                    ia += 2;
                }
                if (n > 0) {
                    d[n - 1] = a[ia - 1];
                }
            } else {
                //
                //              Type 7-12:  generate a diagonally dominant matrix with
                //              unknown condition number in the vectors D and E.
                //
                if (!zerot || !dotype[7 - 1]) {
                    //
                    //                 Let D and E have values from [-1,1].
                    //
                    Rlarnv(2, iseed, n, d);
                    Rlarnv(2, iseed, n - 1, e);
                    //
                    //                 Make the tridiagonal matrix diagonally dominant.
                    //
                    if (n == 1) {
                        d[1 - 1] = abs(d[1 - 1]);
                    } else {
                        d[1 - 1] = abs(d[1 - 1]) + abs(e[1 - 1]);
                        d[n - 1] = abs(d[n - 1]) + abs(e[(n - 1) - 1]);
                        for (i = 2; i <= n - 1; i = i + 1) {
                            d[i - 1] = abs(d[i - 1]) + abs(e[i - 1]) + abs(e[(i - 1) - 1]);
                        }
                    }
                    //
                    //                 Scale D and E so the maximum element is ANORM.
                    //
                    ix = iRamax(n, d, 1);
                    dmax = d[ix - 1];
                    Rscal(n, anorm / dmax, d, 1);
                    Rscal(n - 1, anorm / dmax, e, 1);
                    //
                } else if (izero > 0) {
                    //
                    //                 Reuse the last matrix by copying back the zeroed out
                    //                 elements.
                    //
                    if (izero == 1) {
                        d[1 - 1] = z[2 - 1];
                        if (n > 1) {
                            e[1 - 1] = z[3 - 1];
                        }
                    } else if (izero == n) {
                        e[(n - 1) - 1] = z[1 - 1];
                        d[n - 1] = z[2 - 1];
                    } else {
                        e[(izero - 1) - 1] = z[1 - 1];
                        d[izero - 1] = z[2 - 1];
                        e[izero - 1] = z[3 - 1];
                    }
                }
                //
                //              For types 8-10, set one row and column of the matrix to
                //              zero.
                //
                izero = 0;
                if (imat == 8) {
                    izero = 1;
                    z[2 - 1] = d[1 - 1];
                    d[1 - 1] = zero;
                    if (n > 1) {
                        z[3 - 1] = e[1 - 1];
                        e[1 - 1] = zero;
                    }
                } else if (imat == 9) {
                    izero = n;
                    if (n > 1) {
                        z[1 - 1] = e[(n - 1) - 1];
                        e[(n - 1) - 1] = zero;
                    }
                    z[2 - 1] = d[n - 1];
                    d[n - 1] = zero;
                } else if (imat == 10) {
                    izero = (n + 1) / 2;
                    if (izero > 1) {
                        z[1 - 1] = e[(izero - 1) - 1];
                        e[(izero - 1) - 1] = zero;
                        z[3 - 1] = e[izero - 1];
                        e[izero - 1] = zero;
                    }
                    z[2 - 1] = d[izero - 1];
                    d[izero - 1] = zero;
                }
            }
            //
            Rcopy(n, d, 1, &d[(n + 1) - 1], 1);
            if (n > 1) {
                Rcopy(n - 1, e, 1, &e[(n + 1) - 1], 1);
            }
            //
            //+    TEST 1
            //           Factor A as L*D*L' and compute the ratio
            //              norm(L*D*L' - A) / (n * norm(A) * EPS )
            //
            Rpttrf(n, &d[(n + 1) - 1], &e[(n + 1) - 1], info);
            //
            //           Check error code from Rpttrf.
            //
            if (info != izero) {
                Alaerh(path, "Rpttrf", info, izero, " ", n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                goto statement_100;
            }
            //
            if (info > 0) {
                rcondc = zero;
                goto statement_90;
            }
            //
            Rptt01(n, d, e, &d[(n + 1) - 1], &e[(n + 1) - 1], work, result[1 - 1]);
            //
            //           Print the test ratio if greater than or equal to THRESH.
            //
            if (result[1 - 1] >= thresh) {
                if (nfail == 0 && nerrs == 0) {
                    Alahd(nout, path);
                }
                sprintnum_short(buf, result[1 - 1]);
                write(nout, format_9999), n, imat, 1, buf;
                nfail++;
            }
            nrun++;
            //
            //           Compute RCONDC = 1 / (norm(A) * norm(inv(A))
            //
            //           Compute norm(A).
            //
            anorm = Rlanst("1", n, d, e);
            //
            //           Use Rpttrs to solve for one column at a time of inv(A),
            //           computing the maximum column sum as we go.
            //
            ainvnm = zero;
            for (i = 1; i <= n; i = i + 1) {
                for (j = 1; j <= n; j = j + 1) {
                    x[j - 1] = zero;
                }
                x[i - 1] = one;
                Rpttrs(n, 1, &d[(n + 1) - 1], &e[(n + 1) - 1], x, lda, info);
                ainvnm = max({ainvnm, Rasum(n, x, 1)});
            }
            rcondc = one / max(one, anorm * ainvnm);
            //
            for (irhs = 1; irhs <= nns; irhs = irhs + 1) {
                nrhs = nsval[irhs - 1];
                //
                //           Generate NRHS random solution vectors.
                //
                ix = 1;
                for (j = 1; j <= nrhs; j = j + 1) {
                    Rlarnv(2, iseed, n, &xact[ix - 1]);
                    ix += lda;
                }
                //
                //           Set the right hand side.
                //
                Rlaptm(n, nrhs, one, d, e, xact, lda, zero, b, lda);
                //
                //+    TEST 2
                //           Solve A*x = b and compute the residual.
                //
                Rlacpy("Full", n, nrhs, b, lda, x, lda);
                Rpttrs(n, nrhs, &d[(n + 1) - 1], &e[(n + 1) - 1], x, lda, info);
                //
                //           Check error code from Rpttrs.
                //
                if (info != 0) {
                    Alaerh(path, "Rpttrs", info, 0, " ", n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                }
                //
                Rlacpy("Full", n, nrhs, b, lda, work, lda);
                Rptt02(n, nrhs, d, e, x, lda, work, lda, result[2 - 1]);
                //
                //+    TEST 3
                //           Check solution from generated exact solution.
                //
                Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                //
                //+    TESTS 4, 5, and 6
                //           Use iterative refinement to improve the solution.
                //
                strncpy(srnamt, "Rptrfs", srnamt_len);
                Rptrfs(n, nrhs, d, e, &d[(n + 1) - 1], &e[(n + 1) - 1], b, lda, x, lda, rwork, &rwork[(nrhs + 1) - 1], work, info);
                //
                //           Check error code from Rptrfs.
                //
                if (info != 0) {
                    Alaerh(path, "Rptrfs", info, 0, " ", n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                }
                //
                Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[4 - 1]);
                Rptt05(n, nrhs, d, e, b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], &result[5 - 1]);
                //
                //           Print information about the tests that did not pass the
                //           threshold.
                //
                for (k = 2; k <= 6; k = k + 1) {
                    if (result[k - 1] >= thresh) {
                        if (nfail == 0 && nerrs == 0) {
                            Alahd(nout, path);
                        }
                        sprintnum_short(buf, result[k - 1]);
                        write(nout, "(' N =',i5,', NRHS=',i3,', type ',i2,', test(',i2,') = ',"
                                    "a)"),
                            n, nrhs, imat, k, buf;
                        nfail++;
                    }
                }
                nrun += 5;
            }
        //
        //+    TEST 7
        //           Estimate the reciprocal of the condition number of the
        //           matrix.
        //
        statement_90:
            strncpy(srnamt, "Rptcon", srnamt_len);
            Rptcon(n, &d[(n + 1) - 1], &e[(n + 1) - 1], anorm, rcond, rwork, info);
            //
            //           Check error code from Rptcon.
            //
            if (info != 0) {
                Alaerh(path, "Rptcon", info, 0, " ", n, n, -1, -1, -1, imat, nfail, nerrs, nout);
            }
            //
            result[7 - 1] = Rget06(rcond, rcondc);
            //
            //           Print the test ratio if greater than or equal to THRESH.
            //
            if (result[7 - 1] >= thresh) {
                if (nfail == 0 && nerrs == 0) {
                    Alahd(nout, path);
                }
                sprintnum_short(buf, result[7 - 1]);
                write(nout, format_9999), n, imat, 7, buf;
                nfail++;
            }
            nrun++;
        statement_100:;
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Rchkpt
    //
}
