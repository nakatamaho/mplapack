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

void Cchkhp(common &cmn, bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, REAL const thresh, bool const tsterr, INTEGER const /* nmax */, COMPLEX *a, COMPLEX *afac, COMPLEX *ainv, COMPLEX *b, COMPLEX *x, COMPLEX *xact, COMPLEX *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
    FEM_CMN_SVE(Cchkhp);
    common_write write(cmn);
    str<32> &srnamt = cmn.srnamt;
    //
    if (is_called_first_time) {
        {
            static const INTEGER values[] = {1988, 1989, 1990, 1991};
            data_of_type<int>(FEM_VALUES_AND_SIZE), iseedy;
        }
        {
            static const char *values[] = {"U", "L"};
            data_of_type_str(FEM_VALUES_AND_SIZE), uplos;
        }
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
    char xtype = char0;
    const INTEGER ntypes = 10;
    INTEGER nimat = 0;
    INTEGER izero = 0;
    INTEGER imat = 0;
    bool zerot = false;
    INTEGER iuplo = 0;
    char uplo = char0;
    char packit = char0;
    char type = char0;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist = char0;
    INTEGER info = 0;
    INTEGER ioff = 0;
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER i2 = 0;
    INTEGER i1 = 0;
    INTEGER npp = 0;
    INTEGER k = 0;
    bool trfcon = false;
    const INTEGER ntests = 8;
    arr_1d<ntests, REAL> result(fill0);
    INTEGER nt = 0;
    REAL rcondc = 0.0;
    INTEGER irhs = 0;
    INTEGER nrhs = 0;
    REAL rcond = 0.0;
    static const char *format_9999 = "(' UPLO = ''',a1,''', N =',i5,', type ',i2,', test ',i2,', ratio =',"
                                     "g12.5)";
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
    //     Initialize constants and the random number seed.
    //
    path[(1 - 1)] = "Zomplex precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "HP";
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
        Cerrsy(path, nout);
    }
    cmn.infot = 0;
    //
    //     Do for each value of N in NVAL
    //
    for (in = 1; in <= nn; in = in + 1) {
        n = nval[in - 1];
        lda = max(n, 1);
        xtype = "N";
        nimat = ntypes;
        if (n <= 0) {
            nimat = 1;
        }
        //
        izero = 0;
        for (imat = 1; imat <= nimat; imat = imat + 1) {
            //
            //           Do the tests only if DOTYPE( IMAT ) is true.
            //
            if (!dotype[imat - 1]) {
                goto statement_160;
            }
            //
            //           Skip types 3, 4, 5, or 6 if the matrix size is too small.
            //
            zerot = imat >= 3 && imat <= 6;
            if (zerot && n < imat - 2) {
                goto statement_160;
            }
            //
            //           Do first for UPLO = 'U', then for UPLO = 'L'
            //
            for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                uplo = uplos[iuplo - 1];
                if (Mlsame(uplo, "U")) {
                    packit = "C";
                } else {
                    packit = "R";
                }
                //
                //              Set up parameters with Clatb4 and generate a test matrix
                //              with ZLATMS.
                //
                Clatb4(path, imat, n, n, type, kl, ku, anorm, mode, cndnum, dist);
                //
                srnamt = "ZLATMS";
                zlatms(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, kl, ku, packit, a, lda, work, info);
                //
                //              Check error code from ZLATMS.
                //
                if (info != 0) {
                    Alaerh(path, "ZLATMS", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                    goto statement_150;
                }
                //
                //              For types 3-6, zero one or more rows and columns of
                //              the matrix to test that INFO is returned correctly.
                //
                if (zerot) {
                    if (imat == 3) {
                        izero = 1;
                    } else if (imat == 4) {
                        izero = n;
                    } else {
                        izero = n / 2 + 1;
                    }
                    //
                    if (imat < 6) {
                        //
                        //                    Set row and column IZERO to zero.
                        //
                        if (iuplo == 1) {
                            ioff = (izero - 1) * izero / 2;
                            for (i = 1; i <= izero - 1; i = i + 1) {
                                a[(ioff + i) - 1] = zero;
                            }
                            ioff += izero;
                            for (i = izero; i <= n; i = i + 1) {
                                a[ioff - 1] = zero;
                                ioff += i;
                            }
                        } else {
                            ioff = izero;
                            for (i = 1; i <= izero - 1; i = i + 1) {
                                a[ioff - 1] = zero;
                                ioff += n - i;
                            }
                            ioff = ioff - izero;
                            for (i = izero; i <= n; i = i + 1) {
                                a[(ioff + i) - 1] = zero;
                            }
                        }
                    } else {
                        ioff = 0;
                        if (iuplo == 1) {
                            //
                            //                       Set the first IZERO rows and columns to zero.
                            //
                            for (j = 1; j <= n; j = j + 1) {
                                i2 = min(j, izero);
                                for (i = 1; i <= i2; i = i + 1) {
                                    a[(ioff + i) - 1] = zero;
                                }
                                ioff += j;
                            }
                        } else {
                            //
                            //                       Set the last IZERO rows and columns to zero.
                            //
                            for (j = 1; j <= n; j = j + 1) {
                                i1 = max(j, izero);
                                for (i = i1; i <= n; i = i + 1) {
                                    a[(ioff + i) - 1] = zero;
                                }
                                ioff += n - j;
                            }
                        }
                    }
                } else {
                    izero = 0;
                }
                //
                //              Set the imaginary part of the diagonals.
                //
                if (iuplo == 1) {
                    Claipd(n, a, 2, 1);
                } else {
                    Claipd(n, a, n, -1);
                }
                //
                //              Compute the L*D*L' or U*D*U' factorization of the matrix.
                //
                npp = n * (n + 1) / 2;
                Ccopy(npp, a, 1, afac, 1);
                srnamt = "Chptrf";
                Chptrf(uplo, n, afac, iwork, info);
                //
                //              Adjust the expected value of INFO to account for
                //              pivoting.
                //
                k = izero;
                if (k > 0) {
                statement_100:
                    if (iwork[k - 1] < 0) {
                        if (iwork[k - 1] != -k) {
                            k = -iwork[k - 1];
                            goto statement_100;
                        }
                    } else if (iwork[k - 1] != k) {
                        k = iwork[k - 1];
                        goto statement_100;
                    }
                }
                //
                //              Check error code from Chptrf.
                //
                if (info != k) {
                    Alaerh(path, "Chptrf", info, k, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                }
                if (info != 0) {
                    trfcon = true;
                } else {
                    trfcon = false;
                }
                //
                //+    TEST 1
                //              Reconstruct matrix from factors and compute residual.
                //
                Chpt01(uplo, n, a, afac, iwork, ainv, lda, rwork, result[1 - 1]);
                nt = 1;
                //
                //+    TEST 2
                //              Form the inverse and compute the residual.
                //
                if (!trfcon) {
                    Ccopy(npp, afac, 1, ainv, 1);
                    srnamt = "Chptri";
                    Chptri(uplo, n, ainv, iwork, work, info);
                    //
                    //              Check error code from Chptri.
                    //
                    if (info != 0) {
                        Alaerh(path, "Chptri", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                    }
                    //
                    Cppt03(uplo, n, a, ainv, work, lda, rwork, rcondc, result[2 - 1]);
                    nt = 2;
                }
                //
                //              Print information about the tests that did not pass
                //              the threshold.
                //
                for (k = 1; k <= nt; k = k + 1) {
                    if (result[k - 1] >= thresh) {
                        if (nfail == 0 && nerrs == 0) {
                            Alahd(nout, path);
                        }
                        write(nout, format_9999), uplo, n, imat, k, result(k);
                        nfail++;
                    }
                }
                nrun += nt;
                //
                //              Do only the condition estimate if INFO is not 0.
                //
                if (trfcon) {
                    rcondc = zero;
                    goto statement_140;
                }
                //
                for (irhs = 1; irhs <= nns; irhs = irhs + 1) {
                    nrhs = nsval[irhs - 1];
                    //
                    //+    TEST 3
                    //              Solve and compute residual for  A * X = B.
                    //
                    srnamt = "Clarhs";
                    Clarhs(path, xtype, uplo, " ", n, n, kl, ku, nrhs, a, lda, xact, lda, b, lda, iseed, info);
                    xtype = "C";
                    Clacpy("Full", n, nrhs, b, lda, x, lda);
                    //
                    srnamt = "Chptrs";
                    Chptrs(uplo, n, nrhs, afac, iwork, x, lda, info);
                    //
                    //              Check error code from Chptrs.
                    //
                    if (info != 0) {
                        Alaerh(path, "Chptrs", info, 0, uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                    }
                    //
                    Clacpy("Full", n, nrhs, b, lda, work, lda);
                    Cppt02(uplo, n, nrhs, a, x, lda, work, lda, rwork, result[3 - 1]);
                    //
                    //+    TEST 4
                    //              Check solution from generated exact solution.
                    //
                    Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[4 - 1]);
                    //
                    //+    TESTS 5, 6, and 7
                    //              Use iterative refinement to improve the solution.
                    //
                    srnamt = "ChprFS";
                    Chprfs(uplo, n, nrhs, a, afac, iwork, b, lda, x, lda, rwork, &rwork[(nrhs + 1) - 1], work, &rwork[(2 * nrhs + 1) - 1], info);
                    //
                    //              Check error code from ChprFS.
                    //
                    if (info != 0) {
                        Alaerh(path, "ChprFS", info, 0, uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                    }
                    //
                    Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[5 - 1]);
                    Cppt05(uplo, n, nrhs, a, b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], result[6 - 1]);
                    //
                    //                 Print information about the tests that did not pass
                    //                 the threshold.
                    //
                    for (k = 3; k <= 7; k = k + 1) {
                        if (result[k - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Alahd(nout, path);
                            }
                            write(nout, "(' UPLO = ''',a1,''', N =',i5,', NRHS=',i3,', type ',i2,"
                                        "', test(',i2,') =',g12.5)"),
                                uplo, n, nrhs, imat, k, result(k);
                            nfail++;
                        }
                    }
                    nrun += 5;
                }
            //
            //+    TEST 8
            //              Get an estimate of RCOND = 1/CNDNUM.
            //
            statement_140:
                anorm = Clanhp("1", uplo, n, a, rwork);
                srnamt = "Chpcon";
                Chpcon(uplo, n, afac, iwork, anorm, rcond, work, info);
                //
                //              Check error code from Chpcon.
                //
                if (info != 0) {
                    Alaerh(path, "Chpcon", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                }
                //
                result[8 - 1] = Rget06[(rcond - 1) + (rcondc - 1) * ldRget06];
                //
                //              Print the test ratio if it is .GE. THRESH.
                //
                if (result[8 - 1] >= thresh) {
                    if (nfail == 0 && nerrs == 0) {
                        Alahd(nout, path);
                    }
                    write(nout, format_9999), uplo, n, imat, 8, result(8);
                    nfail++;
                }
                nrun++;
            statement_150:;
            }
        statement_160:;
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Cchkhp
    //
}
