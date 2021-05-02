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

void Cchkpp(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, REAL const thresh, bool const tsterr, INTEGER const  /* nmax */, COMPLEX *a, COMPLEX *afac, COMPLEX *ainv, COMPLEX *b, COMPLEX *x, COMPLEX *xact, COMPLEX *work, REAL *rwork, INTEGER const nout) {
    FEM_CMN_SVE(Cchkpp);
    common_write write(cmn);
    char &srnamt = cmn.srnamt;
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
        {
            static const char *values[] = {"C", "R"};
            data_of_type_str(FEM_VALUES_AND_SIZE), packs;
        }
    }
    char path[3];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    arr_1d<4, int> iseed;
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    char xtype[1];
    const INTEGER ntypes = 9;
    INTEGER nimat = 0;
    INTEGER imat = 0;
    bool zerot = false;
    INTEGER iuplo = 0;
    char uplo[1];
    char packit[1];
    char type[1];
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist[1];
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER ioff = 0;
    const REAL zero = 0.0;
    INTEGER npp = 0;
    const INTEGER ntests = 8;
    arr_1d<ntests, REAL> result;
    REAL rcondc = 0.0;
    INTEGER k = 0;
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
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Initialize constants and the random number seed.
    //
    path[(1 - 1)] = "Zomplex precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "PP";
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
        Cerrpo(path, nout);
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
        for (imat = 1; imat <= nimat; imat = imat + 1) {
            //
            //           Do the tests only if DOTYPE( IMAT ) is true.
            //
            if (!dotype[imat - 1]) {
                goto statement_100;
            }
            //
            //           Skip types 3, 4, or 5 if the matrix size is too small.
            //
            zerot = imat >= 3 && imat <= 5;
            if (zerot && n < imat - 2) {
                goto statement_100;
            }
            //
            //           Do first for UPLO = 'U', then for UPLO = 'L'
            //
            for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                uplo = uplos[iuplo - 1];
                packit = packs[iuplo - 1];
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
                    goto statement_90;
                }
                //
                //              For types 3-5, zero one row and column of the matrix to
                //              test that INFO is returned correctly.
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
                    //                 Set row and column IZERO of A to 0.
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
                //              Compute the L*L' or U'*U factorization of the matrix.
                //
                npp = n * (n + 1) / 2;
                Ccopy(npp, a, 1, afac, 1);
                srnamt = "Cpptrf";
                Cpptrf(uplo, n, afac, info);
                //
                //              Check error code from Cpptrf.
                //
                if (info != izero) {
                    Alaerh(path, "Cpptrf", info, izero, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                    goto statement_90;
                }
                //
                //              Skip the tests if INFO is not 0.
                //
                if (info != 0) {
                    goto statement_90;
                }
                //
                //+    TEST 1
                //              Reconstruct matrix from factors and compute residual.
                //
                Ccopy(npp, afac, 1, ainv, 1);
                Cppt01(uplo, n, a, ainv, rwork, result[1 - 1]);
                //
                //+    TEST 2
                //              Form the inverse and compute the residual.
                //
                Ccopy(npp, afac, 1, ainv, 1);
                srnamt = "Cpptri";
                Cpptri(uplo, n, ainv, info);
                //
                //              Check error code from Cpptri.
                //
                if (info != 0) {
                    Alaerh(path, "Cpptri", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                }
                //
                Cppt03(uplo, n, a, ainv, work, lda, rwork, rcondc, result[2 - 1]);
                //
                //              Print information about the tests that did not pass
                //              the threshold.
                //
                for (k = 1; k <= 2; k = k + 1) {
                    if (result[k - 1] >= thresh) {
                        if (nfail == 0 && nerrs == 0) {
                            Alahd(nout, path);
                        }
                        write(nout, format_9999), uplo, n, imat, k, result(k);
                        nfail++;
                    }
                }
                nrun += 2;
                //
                for (irhs = 1; irhs <= nns; irhs = irhs + 1) {
                    nrhs = nsval[irhs - 1];
                    //
                    //+    TEST 3
                    //              Solve and compute residual for  A * X = B.
                    //
                    srnamt = "Clarhs";
                    Clarhs(path, xtype, uplo, " ", n, n, kl, ku, nrhs, a, lda, xact, lda, b, lda, iseed, info);
                    Clacpy("Full", n, nrhs, b, lda, x, lda);
                    //
                    srnamt = "Cpptrs";
                    Cpptrs(uplo, n, nrhs, afac, x, lda, info);
                    //
                    //              Check error code from Cpptrs.
                    //
                    if (info != 0) {
                        Alaerh(path, "Cpptrs", info, 0, uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
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
                    srnamt = "Cpprfs";
                    Cpprfs(uplo, n, nrhs, a, afac, b, lda, x, lda, rwork, &rwork[(nrhs + 1) - 1], work, &rwork[(2 * nrhs + 1) - 1], info);
                    //
                    //              Check error code from Cpprfs.
                    //
                    if (info != 0) {
                        Alaerh(path, "Cpprfs", info, 0, uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
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
                anorm = Clanhp("1", uplo, n, a, rwork);
                srnamt = "Cppcon";
                Cppcon(uplo, n, afac, anorm, rcond, work, rwork, info);
                //
                //              Check error code from Cppcon.
                //
                if (info != 0) {
                    Alaerh(path, "Cppcon", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                }
                //
                result[8 - 1] = Rget06[(rcond - 1) + (rcondc - 1) * ldRget06];
                //
                //              Print the test ratio if greater than or equal to THRESH.
                //
                if (result[8 - 1] >= thresh) {
                    if (nfail == 0 && nerrs == 0) {
                        Alahd(nout, path);
                    }
                    write(nout, format_9999), uplo, n, imat, 8, result(8);
                    nfail++;
                }
                nrun++;
            //
            statement_90:;
            }
        statement_100:;
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Cchkpp
    //
}
