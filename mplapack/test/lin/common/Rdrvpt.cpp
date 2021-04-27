/*
 * Copyright (c) 2008-2021
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

void Rdrvpt(common &cmn, bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nrhs, REAL const thresh, bool const tsterr, REAL *a, REAL *d, REAL *e, REAL *b, REAL *x, REAL *xact, REAL *work, REAL *rwork, INTEGER const nout) {
    FEM_CMN_SVE(Rdrvpt);
    common_write write(cmn);
    str<32> &srnamt = cmn.srnamt;
    //
    if (is_called_first_time) {
        static const INTEGER values[] = {0, 0, 0, 1};
        data_of_type<int>(FEM_VALUES_AND_SIZE), iseedy;
    }
    str<3> path = char0;
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    arr_1d<4, int> iseed(fill0);
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    const INTEGER ntypes = 12;
    INTEGER nimat = 0;
    INTEGER imat = 0;
    char type = char0;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cond = 0.0;
    char dist = char0;
    bool zerot = false;
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER ia = 0;
    INTEGER ix = 0;
    REAL dmax = 0.0;
    arr_1d<3, REAL> z(fill0);
    const REAL zero = 0.0;
    INTEGER j = 0;
    const REAL one = 1.0;
    INTEGER ifact = 0;
    char fact = char0;
    REAL rcondc = 0.0;
    REAL ainvnm = 0.0;
    INTEGER nt = 0;
    const INTEGER ntests = 6;
    arr_1d<ntests, REAL> result(fill0);
    INTEGER k = 0;
    REAL rcond = 0.0;
    INTEGER k1 = 0;
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
    path[(1 - 1)] = "Double precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "PT";
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
        Rerrvx(path, nout);
    }
    cmn.infot = 0;
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
                goto statement_110;
            }
            //
            //           Set up parameters with Rlatb4.
            //
            Rlatb4(path, imat, n, n, type, kl, ku, anorm, mode, cond, dist);
            //
            zerot = imat >= 8 && imat <= 10;
            if (imat <= 6) {
                //
                //              Type 1-6:  generate a symmetric tridiagonal matrix of
                //              known condition number in lower triangular band storage.
                //
                srnamt = "DLATMS";
                dlatms(n, n, dist, iseed, type, rwork, mode, cond, anorm, kl, ku, "B", a, 2, work, info);
                //
                //              Check the error code from DLATMS.
                //
                if (info != 0) {
                    Alaerh(path, "DLATMS", info, 0, " ", n, n, kl, ku, -1, imat, nfail, nerrs, nout);
                    goto statement_110;
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
                    if (n > 1) {
                        Rscal(n - 1, anorm / dmax, e, 1);
                    }
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
                        z[3 - 1] = e[izero - 1];
                        e[(izero - 1) - 1] = zero;
                        e[izero - 1] = zero;
                    }
                    z[2 - 1] = d[izero - 1];
                    d[izero - 1] = zero;
                }
            }
            //
            //           Generate NRHS random solution vectors.
            //
            ix = 1;
            for (j = 1; j <= nrhs; j = j + 1) {
                Rlarnv(2, iseed, n, xact[ix - 1]);
                ix += lda;
            }
            //
            //           Set the right hand side.
            //
            Rlaptm(n, nrhs, one, d, e, xact, lda, zero, b, lda);
            //
            for (ifact = 1; ifact <= 2; ifact = ifact + 1) {
                if (ifact == 1) {
                    fact = "F";
                } else {
                    fact = "N";
                }
                //
                //              Compute the condition number for comparison with
                //              the value returned by Rptsvx.
                //
                if (zerot) {
                    if (ifact == 1) {
                        goto statement_100;
                    }
                    rcondc = zero;
                    //
                } else if (ifact == 1) {
                    //
                    //                 Compute the 1-norm of A.
                    //
                    anorm = Rlanst("1", n, d, e);
                    //
                    Rcopy(n, d, 1, &d[(n + 1) - 1], 1);
                    if (n > 1) {
                        Rcopy(n - 1, e, 1, &e[(n + 1) - 1], 1);
                    }
                    //
                    //                 Factor the matrix A.
                    //
                    Rpttrf(n, &d[(n + 1) - 1], &e[(n + 1) - 1], info);
                    //
                    //                 Use Rpttrs to solve for one column at a time of
                    //                 inv(A), computing the maximum column sum as we go.
                    //
                    ainvnm = zero;
                    for (i = 1; i <= n; i = i + 1) {
                        for (j = 1; j <= n; j = j + 1) {
                            x[j - 1] = zero;
                        }
                        x[i - 1] = one;
                        Rpttrs(n, 1, &d[(n + 1) - 1], &e[(n + 1) - 1], x, lda, info);
                        ainvnm = max(ainvnm, Rasum(n, x, 1));
                    }
                    //
                    //                 Compute the 1-norm condition number of A.
                    //
                    if (anorm <= zero || ainvnm <= zero) {
                        rcondc = one;
                    } else {
                        rcondc = (one / anorm) / ainvnm;
                    }
                }
                //
                if (ifact == 2) {
                    //
                    //                 --- Test Rptsv --
                    //
                    Rcopy(n, d, 1, &d[(n + 1) - 1], 1);
                    if (n > 1) {
                        Rcopy(n - 1, e, 1, &e[(n + 1) - 1], 1);
                    }
                    Rlacpy("Full", n, nrhs, b, lda, x, lda);
                    //
                    //                 Factor A as L*D*L' and solve the system A*X = B.
                    //
                    srnamt = "Rptsv ";
                    Rptsv(n, nrhs, &d[(n + 1) - 1], &e[(n + 1) - 1], x, lda, info);
                    //
                    //                 Check error code from Rptsv .
                    //
                    if (info != izero) {
                        Alaerh(path, "Rptsv ", info, izero, " ", n, n, 1, 1, nrhs, imat, nfail, nerrs, nout);
                    }
                    nt = 0;
                    if (izero == 0) {
                        //
                        //                    Check the factorization by computing the ratio
                        //                       norm(L*D*L' - A) / (n * norm(A) * EPS )
                        //
                        Rptt01(n, d, e, &d[(n + 1) - 1], &e[(n + 1) - 1], work, result[1 - 1]);
                        //
                        //                    Compute the residual in the solution.
                        //
                        Rlacpy("Full", n, nrhs, b, lda, work, lda);
                        Rptt02(n, nrhs, d, e, x, lda, work, lda, result[2 - 1]);
                        //
                        //                    Check solution from generated exact solution.
                        //
                        Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                        nt = 3;
                    }
                    //
                    //                 Print information about the tests that did not pass
                    //                 the threshold.
                    //
                    for (k = 1; k <= nt; k = k + 1) {
                        if (result[k - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Aladhd(nout, path);
                            }
                            write(nout, "(1x,a,', N =',i5,', type ',i2,', test ',i2,', ratio = ',"
                                        "g12.5)"),
                                "Rptsv ", n, imat, k, result(k);
                            nfail++;
                        }
                    }
                    nrun += nt;
                }
                //
                //              --- Test Rptsvx ---
                //
                if (ifact > 1) {
                    //
                    //                 Initialize D( N+1:2*N ) and E( N+1:2*N ) to zero.
                    //
                    for (i = 1; i <= n - 1; i = i + 1) {
                        d[(n + i) - 1] = zero;
                        e[(n + i) - 1] = zero;
                    }
                    if (n > 0) {
                        d[(n + n) - 1] = zero;
                    }
                }
                //
                Rlaset("Full", n, nrhs, zero, zero, x, lda);
                //
                //              Solve the system and compute the condition number and
                //              error bounds using Rptsvx.
                //
                srnamt = "Rptsvx";
                Rptsvx(fact, n, nrhs, d, e, &d[(n + 1) - 1], &e[(n + 1) - 1], b, lda, x, lda, rcond, rwork, &rwork[(nrhs + 1) - 1], work, info);
                //
                //              Check the error code from Rptsvx.
                //
                if (info != izero) {
                    Alaerh(path, "Rptsvx", info, izero, fact, n, n, 1, 1, nrhs, imat, nfail, nerrs, nout);
                }
                if (izero == 0) {
                    if (ifact == 2) {
                        //
                        //                    Check the factorization by computing the ratio
                        //                       norm(L*D*L' - A) / (n * norm(A) * EPS )
                        //
                        k1 = 1;
                        Rptt01(n, d, e, &d[(n + 1) - 1], &e[(n + 1) - 1], work, result[1 - 1]);
                    } else {
                        k1 = 2;
                    }
                    //
                    //                 Compute the residual in the solution.
                    //
                    Rlacpy("Full", n, nrhs, b, lda, work, lda);
                    Rptt02(n, nrhs, d, e, x, lda, work, lda, result[2 - 1]);
                    //
                    //                 Check solution from generated exact solution.
                    //
                    Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                    //
                    //                 Check error bounds from iterative refinement.
                    //
                    Rptt05(n, nrhs, d, e, b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], result[4 - 1]);
                } else {
                    k1 = 6;
                }
                //
                //              Check the reciprocal of the condition number.
                //
                result[6 - 1] = Rget06[(rcond - 1) + (rcondc - 1) * ldRget06];
                //
                //              Print information about the tests that did not pass
                //              the threshold.
                //
                for (k = k1; k <= 6; k = k + 1) {
                    if (result[k - 1] >= thresh) {
                        if (nfail == 0 && nerrs == 0) {
                            Aladhd(nout, path);
                        }
                        write(nout, "(1x,a,', FACT=''',a1,''', N =',i5,', type ',i2,', test ',i2,"
                                    "', ratio = ',g12.5)"),
                            "Rptsvx", fact, n, imat, k, result(k);
                        nfail++;
                    }
                }
                nrun += 7 - k1;
            statement_100:;
            }
        statement_110:;
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasvm(path, nout, nfail, nrun, nerrs);
    //
    //     End of Rdrvpt
    //
}
