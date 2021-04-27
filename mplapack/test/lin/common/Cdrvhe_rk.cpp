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

void Cdrvhe_rk(common &cmn, bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nrhs, REAL const thresh, bool const tsterr, INTEGER const nmax, COMPLEX *a, COMPLEX *afac, COMPLEX *e, COMPLEX *ainv, COMPLEX *b, COMPLEX *x, COMPLEX *xact, COMPLEX *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
    FEM_CMN_SVE(Cdrvhe_rk);
    common_write write(cmn);
    str<32> &srnamt = cmn.srnamt;
    //
    const INTEGER nfact = 2;
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
            static const char *values[] = {"F", "N"};
            data_of_type_str(FEM_VALUES_AND_SIZE), facts;
        }
    }
    str<3> path = char0;
    str<3> matpath = char0;
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    arr_1d<4, int> iseed(fill0);
    INTEGER lwork = 0;
    INTEGER nb = 0;
    INTEGER nbmin = 0;
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    char xtype = char0;
    const INTEGER ntypes = 10;
    INTEGER nimat = 0;
    INTEGER imat = 0;
    bool zerot = false;
    INTEGER iuplo = 0;
    char uplo = char0;
    char type = char0;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist = char0;
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER ioff = 0;
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER i2 = 0;
    INTEGER i1 = 0;
    INTEGER ifact = 0;
    char fact = char0;
    REAL rcondc = 0.0;
    REAL ainvnm = 0.0;
    const REAL one = 1.0;
    INTEGER k = 0;
    const INTEGER ntests = 3;
    arr_1d<ntests, REAL> result(fill0);
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
    //
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
    //     Test path
    //
    path[(1 - 1)] = "Zomplex precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "HK";
    //
    //     Path to generate matrices
    //
    matpath[(1 - 1)] = "Zomplex precision";
    matpath[(2 - 1) + (3 - 1) * ldmatpath] = "HE";
    //
    nrun = 0;
    nfail = 0;
    nerrs = 0;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    lwork = max((INTEGER)2 * nmax, nmax * nrhs);
    //
    //     Test the error exits
    //
    if (tsterr) {
        Cerrvx(path, nout);
    }
    cmn.infot = 0;
    //
    //     Set the block size and minimum block size for which the block
    //     routine should be used, which will be later returned by ILAENV.
    //
    nb = 1;
    nbmin = 2;
    xlaenv(1, nb);
    xlaenv(2, nbmin);
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
                //                 Begin generate the test matrix A.
                //
                //                 Set up parameters with Clatb4 for the matrix generator
                //                 based on the type of matrix to be generated.
                //
                Clatb4(matpath, imat, n, n, type, kl, ku, anorm, mode, cndnum, dist);
                //
                //                 Generate a matrix with ZLATMS.
                //
                srnamt = "ZLATMS";
                zlatms(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, kl, ku, uplo, a, lda, work, info);
                //
                //                 Check error code from ZLATMS and handle error.
                //
                if (info != 0) {
                    Alaerh(path, "ZLATMS", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                    goto statement_160;
                }
                //
                //                 For types 3-6, zero one or more rows and columns of
                //                 the matrix to test that INFO is returned correctly.
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
                        //                       Set row and column IZERO to zero.
                        //
                        if (iuplo == 1) {
                            ioff = (izero - 1) * lda;
                            for (i = 1; i <= izero - 1; i = i + 1) {
                                a[(ioff + i) - 1] = zero;
                            }
                            ioff += izero;
                            for (i = izero; i <= n; i = i + 1) {
                                a[ioff - 1] = zero;
                                ioff += lda;
                            }
                        } else {
                            ioff = izero;
                            for (i = 1; i <= izero - 1; i = i + 1) {
                                a[ioff - 1] = zero;
                                ioff += lda;
                            }
                            ioff = ioff - izero;
                            for (i = izero; i <= n; i = i + 1) {
                                a[(ioff + i) - 1] = zero;
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
                                    a[(ioff + i) - 1] = zero;
                                }
                                ioff += lda;
                            }
                        } else {
                            //
                            //                       Set the first IZERO rows and columns to zero.
                            //
                            ioff = 0;
                            for (j = 1; j <= n; j = j + 1) {
                                i1 = max(j, izero);
                                for (i = i1; i <= n; i = i + 1) {
                                    a[(ioff + i) - 1] = zero;
                                }
                                ioff += lda;
                            }
                        }
                    }
                } else {
                    izero = 0;
                }
                //
                //                 End generate the test matrix A.
                //
                for (ifact = 1; ifact <= nfact; ifact = ifact + 1) {
                    //
                    //                 Do first for FACT = 'F', then for other values.
                    //
                    fact = facts[ifact - 1];
                    //
                    //                 Compute the condition number
                    //
                    if (zerot) {
                        if (ifact == 1) {
                            goto statement_150;
                        }
                        rcondc = zero;
                        //
                    } else if (ifact == 1) {
                        //
                        //                    Compute the 1-norm of A.
                        //
                        anorm = Clanhe("1", uplo, n, a, lda, rwork);
                        //
                        //                    Factor the matrix A.
                        //
                        Clacpy(uplo, n, n, a, lda, afac, lda);
                        Chetrf_rk(uplo, n, afac, lda, e, iwork, work, lwork, info);
                        //
                        //                    Compute inv(A) and take its norm.
                        //
                        Clacpy(uplo, n, n, afac, lda, ainv, lda);
                        lwork = (n + nb + 1) * (nb + 3);
                        //
                        //                    We need to compute the inverse to compute
                        //                    RCONDC that is used later in TEST3.
                        //
                        Chetri_3(uplo, n, ainv, lda, e, iwork, work, lwork, info);
                        ainvnm = Clanhe("1", uplo, n, ainv, lda, rwork);
                        //
                        //                    Compute the 1-norm condition number of A.
                        //
                        if (anorm <= zero || ainvnm <= zero) {
                            rcondc = one;
                        } else {
                            rcondc = (one / anorm) / ainvnm;
                        }
                    }
                    //
                    //                 Form an exact solution and set the right hand side.
                    //
                    srnamt = "Clarhs";
                    Clarhs(matpath, xtype, uplo, " ", n, n, kl, ku, nrhs, a, lda, xact, lda, b, lda, iseed, info);
                    xtype = "C";
                    //
                    //                 --- Test Chesv_rk  ---
                    //
                    if (ifact == 2) {
                        Clacpy(uplo, n, n, a, lda, afac, lda);
                        Clacpy("Full", n, nrhs, b, lda, x, lda);
                        //
                        //                    Factor the matrix and solve the system using
                        //                    Chesv_rk.
                        //
                        srnamt = "Chesv_rk";
                        Chesv_rk(uplo, n, nrhs, afac, lda, e, iwork, x, lda, work, lwork, info);
                        //
                        //                    Adjust the expected value of INFO to account for
                        //                    pivoting.
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
                        //                    Check error code from Chesv_rk and handle error.
                        //
                        if (info != k) {
                            Alaerh(path, "Chesv_rk", info, k, uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                            goto statement_120;
                        } else if (info != 0) {
                            goto statement_120;
                        }
                        //
                        //+    TEST 1      Reconstruct matrix from factors and compute
                        //                 residual.
                        //
                        Chet01_3(uplo, n, a, lda, afac, lda, e, iwork, ainv, lda, rwork, result[1 - 1]);
                        //
                        //+    TEST 2      Compute residual of the computed solution.
                        //
                        Clacpy("Full", n, nrhs, b, lda, work, lda);
                        Cpot02(uplo, n, nrhs, a, lda, x, lda, work, lda, rwork, result[2 - 1]);
                        //
                        //+    TEST 3
                        //                 Check solution from generated exact solution.
                        //
                        Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                        nt = 3;
                        //
                        //                    Print information about the tests that did not pass
                        //                    the threshold.
                        //
                        for (k = 1; k <= nt; k = k + 1) {
                            if (result[k - 1] >= thresh) {
                                if (nfail == 0 && nerrs == 0) {
                                    Aladhd(nout, path);
                                }
                                write(nout, "(1x,a,', UPLO=''',a1,''', N =',i5,', type ',i2,', test ',"
                                            "i2,', ratio =',g12.5)"),
                                    "Chesv_rk", uplo, n, imat, k, result(k);
                                nfail++;
                            }
                        }
                        nrun += nt;
                    statement_120:;
                    }
                //
                statement_150:;
                }
            //
            statement_160:;
            }
        statement_170:;
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasvm(path, nout, nfail, nrun, nerrs);
    //
    //     End of Cdrvhe_rk
    //
}
