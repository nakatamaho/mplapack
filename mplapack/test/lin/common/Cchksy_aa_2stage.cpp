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

void Cchksy_aa_2stage(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nnb, INTEGER *nbval, INTEGER const nns, INTEGER *nsval, REAL const thresh, bool const tsterr, INTEGER const nmax, COMPLEX *a, COMPLEX *afac, COMPLEX *ainv, COMPLEX *b, COMPLEX *x, COMPLEX *xact, COMPLEX *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    char uplos[] = {'U', 'L'};
    char path[3];
    char matpath[3];
    char buf[1024];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    char xtype;
    const INTEGER ntypes = 10;
    INTEGER nimat = 0;
    INTEGER izero = 0;
    INTEGER imat = 0;
    bool zerot = false;
    INTEGER iuplo = 0;
    char uplo;
    char type;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist;
    INTEGER info = 0;
    INTEGER ioff = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER j = 0;
    INTEGER i2 = 0;
    INTEGER i1 = 0;
    INTEGER inb = 0;
    INTEGER nb = 0;
    INTEGER lwork = 0;
    INTEGER k = 0;
    INTEGER nt = 0;
    const INTEGER ntests = 9;
    REAL result[ntests];
    INTEGER irhs = 0;
    INTEGER nrhs = 0;
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
    //     Test path
    //
    path[0] = 'C';
    path[1] = 'S';
    path[2] = '2';
    //     Path to generate matrices
    //
    matpath[0] = 'C';
    matpath[1] = 'S';
    matpath[2] = 'Y';
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
    //
    //     Set the minimum block size for which the block routine should
    //     be used, which will be later returned by iMlaenv
    xlaenv(2, 2);
    //
    //
    //     Do for each value of N in NVAL
    //
    for (in = 1; in <= nn; in = in + 1) {
        n = nval[in - 1];
        if (n > nmax) {
            nfail++;
            write(nout, "(' Invalid input value: ',a4,'=',i6,'; must be <=',i6)"), "M ", n, nmax;
            goto statement_180;
        }
        lda = max(n, 1);
        xtype = 'N';
        nimat = ntypes;
        if (n <= 0) {
            nimat = 1;
        }
        //
        izero = 0;
        //
        //        Do for each value of matrix type IMAT
        //
        for (imat = 1; imat <= nimat; imat = imat + 1) {
            //
            //           Do the tests only if DOTYPE( IMAT ) is true.
            //
            if (!dotype[imat - 1]) {
                goto statement_170;
            }
            //
            //           Skip types 3, 4, 5, or 6 if the matrix size is too small.
            //
            zerot = imat >= 3 && imat <= 6;
            if (zerot && n < imat - 2) {
                goto statement_170;
            }
            //
            //           Do first for UPLO = 'U', then for UPLO = 'L'
            //
            for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                uplo = uplos[iuplo - 1];
                //
                //              Begin generate the test matrix A.
                //
                //              Set up parameters with Clatb4 for the matrix generator
                //              based on the type of matrix to be generated.
                //
                Clatb4(matpath, imat, n, n, &type, kl, ku, anorm, mode, cndnum, &dist);
                //
                //              Generate a matrix with Clatms.
                //
                Clatms(n, n, &dist, iseed, &type, rwork, mode, cndnum, anorm, kl, ku, &uplo, a, lda, work, info);
                //
                //              Check error code from Clatms and handle error.
                //
                if (info != 0) {
                    Alaerh(path, "Clatms", info, 0, &uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                    //
                    //                    Skip all tests for this generated matrix
                    //
                    goto statement_160;
                }
                //
                //              For matrix types 3-6, zero one or more rows and
                //              columns of the matrix to test that INFO is returned
                //              correctly.
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
                            ioff = (izero - 1) * lda;
                            for (i = 1; i <= izero - 1; i = i + 1) {
                                a[(ioff + i) - 1] = czero;
                            }
                            ioff += izero;
                            for (i = izero; i <= n; i = i + 1) {
                                a[ioff - 1] = czero;
                                ioff += lda;
                            }
                        } else {
                            ioff = izero;
                            for (i = 1; i <= izero - 1; i = i + 1) {
                                a[ioff - 1] = czero;
                                ioff += lda;
                            }
                            ioff = ioff - izero;
                            for (i = izero; i <= n; i = i + 1) {
                                a[(ioff + i) - 1] = czero;
                            }
                        }
                    } else {
                        if (iuplo == 1) {
                            //
                            //                       Set the first IZERO rows and columns to zero.
                            //
                            ioff = 0;
                            for (j = 1; j <= n; j = j + 1) {
                                i2 = min(j, izero);
                                for (i = 1; i <= i2; i = i + 1) {
                                    a[(ioff + i) - 1] = czero;
                                }
                                ioff += lda;
                            }
                            izero = 1;
                        } else {
                            //
                            //                       Set the last IZERO rows and columns to zero.
                            //
                            ioff = 0;
                            for (j = 1; j <= n; j = j + 1) {
                                i1 = max(j, izero);
                                for (i = i1; i <= n; i = i + 1) {
                                    a[(ioff + i) - 1] = czero;
                                }
                                ioff += lda;
                            }
                        }
                    }
                } else {
                    izero = 0;
                }
                //
                //              End generate the test matrix A.
                //
                //              Do for each value of NB in NBVAL
                //
                for (inb = 1; inb <= nnb; inb = inb + 1) {
                    //
                    //                 Set the optimal blocksize, which will be later
                    //                 returned by iMlaenv.
                    //
                    nb = nbval[inb - 1];
                    xlaenv(1, nb);
                    //
                    //                 Copy the test matrix A into matrix AFAC which
                    //                 will be factorized in place. This is needed to
                    //                 preserve the test matrix A for subsequent tests.
                    //
                    Clacpy(&uplo, n, n, a, lda, afac, lda);
                    //
                    //                 Compute the L*D*L**T or U*D*U**T factorization of the
                    //                 matrix. IWORK stores details of the interchanges and
                    //                 the block structure of D. AINV is a work array for
                    //                 block factorization, LWORK is the length of AINV.
                    //
                    lwork = min(n * nb, 3 * nmax * nmax);
                    Csytrf_aa_2stage(&uplo, n, afac, lda, ainv, (3 * nb + 1) * n, iwork, &iwork[(1 + n) - 1], work, lwork, info);
                    //
                    //                 Adjust the expected value of INFO to account for
                    //                 pivoting.
                    //
                    if (izero > 0) {
                        j = 1;
                        k = izero;
                    statement_100:
                        if (j == k) {
                            k = iwork[j - 1];
                        } else if (iwork[j - 1] == k) {
                            k = j;
                        }
                        if (j < k) {
                            j++;
                            goto statement_100;
                        }
                    } else {
                        k = 0;
                    }
                    //
                    //                 Check error code from Csytrf and handle error.
                    //
                    if (info != k) {
                        Alaerh(path, "Csytrf_aa_2stage", info, k, &uplo, n, n, -1, -1, nb, imat, nfail, nerrs, nout);
                    }
                    //
                    //+    TEST 1
                    //                 Reconstruct matrix from factors and compute residual.
                    //
                    //                  CALL Csyt01_aa( UPLO, N, A, LDA, AFAC, LDA, IWORK,
                    //     $                            AINV, LDA, RWORK, RESULT( 1 ) )
                    //                  NT = 1
                    nt = 0;
                    //
                    //                 Print information about the tests that did not pass
                    //                 the threshold.
                    //
                    for (k = 1; k <= nt; k = k + 1) {
                        if (result[k - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Alahd(nout, path);
                            }
                            sprintnum_short(buf, result[k - 1]);
                            write(nout, "(' UPLO = ''',a1,''', N =',i5,', NB =',i4,', type ',i2,"
                                        "', test ',i2,', ratio =',a)"),
                                uplo, n, nb, imat, k, buf;
                            nfail++;
                        }
                    }
                    nrun += nt;
                    //
                    //                 Skip solver test if INFO is not 0.
                    //
                    if (info != 0) {
                        goto statement_140;
                    }
                    //
                    //                 Do for each value of NRHS in NSVAL.
                    //
                    for (irhs = 1; irhs <= nns; irhs = irhs + 1) {
                        nrhs = nsval[irhs - 1];
                        //
                        //+    TEST 2 (Using TRS)
                        //                 Solve and compute residual for  A * X = B.
                        //
                        //                    Choose a set of NRHS random solution vectors
                        //                    stored in XACT and set up the right hand side B
                        //
                        Clarhs(matpath, &xtype, &uplo, " ", n, n, kl, ku, nrhs, a, lda, xact, lda, b, lda, iseed, info);
                        Clacpy("Full", n, nrhs, b, lda, x, lda);
                        //
                        lwork = max((INTEGER)1, 3 * n - 2);
                        Csytrs_aa_2stage(&uplo, n, nrhs, afac, lda, ainv, (3 * nb + 1) * n, iwork, &iwork[(1 + n) - 1], x, lda, info);
                        //
                        //                    Check error code from Csytrs and handle error.
                        //
                        if (info != 0) {
                            if (izero == 0) {
                                Alaerh(path, "Csytrs_aa_2stage", info, 0, &uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                            }
                        } else {
                            Clacpy("Full", n, nrhs, b, lda, work, lda);
                            //
                            //                       Compute the residual for the solution
                            //
                            Csyt02(&uplo, n, nrhs, a, lda, x, lda, work, lda, rwork, result[2 - 1]);
                            //
                            //                       Print information about the tests that did not pass
                            //                       the threshold.
                            //
                            for (k = 2; k <= 2; k = k + 1) {
                                if (result[k - 1] >= thresh) {
                                    if (nfail == 0 && nerrs == 0) {
                                        Alahd(nout, path);
                                    }
                                    sprintnum_short(buf, result[k - 1]);
                                    write(nout, "(' UPLO = ''',a1,''', N =',i5,', NRHS=',i3,', type ',i2,"
                                                "', test(',i2,') =',a)"),
                                        uplo, n, nrhs, imat, k, buf;
                                    nfail++;
                                }
                            }
                        }
                        nrun++;
                        //
                        //                 End do for each value of NRHS in NSVAL.
                        //
                    }
                statement_140:;
                }
            statement_160:;
            }
        statement_170:;
        }
    statement_180:;
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Cchksy_aa_2stage
    //
}
