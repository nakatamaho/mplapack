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

void Rdrvsy(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nrhs, REAL const thresh, bool const tsterr, INTEGER const nmax, REAL *a, REAL *afac, REAL *ainv, REAL *b, REAL *x, REAL *xact, REAL *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
    common_write write(cmn);
    //
    const INTEGER nfact = 2;
    str_arr_ref<1> facts(sve.facts, [nfact]);
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    str_arr_ref<1> uplos(sve.uplos, [2]);
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
    char path[3];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    INTEGER lwork = 0;
    INTEGER nb = 0;
    INTEGER nbmin = 0;
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    char xtype;
    const INTEGER ntypes = 10;
    INTEGER nimat = 0;
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
    INTEGER izero = 0;
    INTEGER ioff = 0;
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER i2 = 0;
    INTEGER i1 = 0;
    INTEGER ifact = 0;
    char fact;
    REAL rcondc = 0.0;
    REAL ainvnm = 0.0;
    const REAL one = 1.0;
    INTEGER k = 0;
    const INTEGER ntests = 6;
    REAL result[ntests];
    INTEGER nt = 0;
    REAL rcond = 0.0;
    INTEGER k1 = 0;
    INTEGER n_err_bnds = 0;
    char equed;
    REAL rpvgrw_svxx = 0.0;
    REAL berr[nrhs];
    REAL errbnds_n[nrhs * 3];
    REAL errbnds_c[nrhs * 3];
    static const char *format_9998 = "(1x,a,', FACT=''',a1,''', UPLO=''',a1,''', N =',i5,', type ',i2,"
                                     "', test ',i2,', ratio =',g12.5)";
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
    path[(1 - 1)] = "Double precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "SY";
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
        Rerrvx(path, nout);
    }
    cmn.infot = 0;
    //
    //     Set the block size and minimum block size for testing.
    //
    nb = 1;
    nbmin = 2;
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
                //              Set up parameters with Rlatb4 and generate a test matrix
                //              with Rlatms.
                //
                Rlatb4(path, imat, n, n, type, kl, ku, anorm, mode, cndnum, dist);
                //
                Rlatms(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, kl, ku, uplo, a, lda, work, info);
                //
                //              Check error code from Rlatms.
                //
                if (info != 0) {
                    Alaerh(path, "Rlatms", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                    goto statement_160;
                }
                //
                //              For types 3-6, zero one or more rows and columns of the
                //              matrix to test that INFO is returned correctly.
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
                                ioff += lda;
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
                                ioff += lda;
                            }
                        }
                    }
                } else {
                    izero = 0;
                }
                //
                for (ifact = 1; ifact <= nfact; ifact = ifact + 1) {
                    //
                    //                 Do first for FACT = 'F', then for other values.
                    //
                    fact = facts[ifact - 1];
                    //
                    //                 Compute the condition number for comparison with
                    //                 the value returned by Rsysvx.
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
                        anorm = Rlansy("1", uplo, n, a, lda, rwork);
                        //
                        //                    Factor the matrix A.
                        //
                        Rlacpy(uplo, n, n, a, lda, afac, lda);
                        Rsytrf(uplo, n, afac, lda, iwork, work, lwork, info);
                        //
                        //                    Compute inv(A) and take its norm.
                        //
                        Rlacpy(uplo, n, n, afac, lda, ainv, lda);
                        lwork = (n + nb + 1) * (nb + 3);
                        Rsytri2(uplo, n, ainv, lda, iwork, work, lwork, info);
                        ainvnm = Rlansy("1", uplo, n, ainv, lda, rwork);
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
                    Rlarhs(path, xtype, uplo, " ", n, n, kl, ku, nrhs, a, lda, xact, lda, b, lda, iseed, info);
                    xtype = "C";
                    //
                    //                 --- Test Rsysv  ---
                    //
                    if (ifact == 2) {
                        Rlacpy(uplo, n, n, a, lda, afac, lda);
                        Rlacpy("Full", n, nrhs, b, lda, x, lda);
                        //
                        //                    Factor the matrix and solve the system using Rsysv.
                        //
                        Rsysv(uplo, n, nrhs, afac, lda, iwork, x, lda, work, lwork, info);
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
                        //                    Check error code from Rsysv .
                        //
                        if (info != k) {
                            Alaerh(path, "Rsysv ", info, k, uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                            goto statement_120;
                        } else if (info != 0) {
                            goto statement_120;
                        }
                        //
                        //                    Reconstruct matrix from factors and compute
                        //                    residual.
                        //
                        Rsyt01(uplo, n, a, lda, afac, lda, iwork, ainv, lda, rwork, result[1 - 1]);
                        //
                        //                    Compute residual of the computed solution.
                        //
                        Rlacpy("Full", n, nrhs, b, lda, work, lda);
                        Rpot02(uplo, n, nrhs, a, lda, x, lda, work, lda, rwork, result[2 - 1]);
                        //
                        //                    Check solution from generated exact solution.
                        //
                        Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
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
                                    "Rsysv ", uplo, n, imat, k, result(k);
                                nfail++;
                            }
                        }
                        nrun += nt;
                    statement_120:;
                    }
                    //
                    //                 --- Test Rsysvx ---
                    //
                    if (ifact == 2) {
                        Rlaset(uplo, n, n, zero, zero, afac, lda);
                    }
                    Rlaset("Full", n, nrhs, zero, zero, x, lda);
                    //
                    //                 Solve the system and compute the condition number and
                    //                 error bounds using Rsysvx.
                    //
                    Rsysvx(fact, uplo, n, nrhs, a, lda, afac, lda, iwork, b, lda, x, lda, rcond, rwork, &rwork[(nrhs + 1) - 1], work, lwork, &iwork[(n + 1) - 1], info);
                    //
                    //                 Adjust the expected value of INFO to account for
                    //                 pivoting.
                    //
                    k = izero;
                    if (k > 0) {
                    statement_130:
                        if (iwork[k - 1] < 0) {
                            if (iwork[k - 1] != -k) {
                                k = -iwork[k - 1];
                                goto statement_130;
                            }
                        } else if (iwork[k - 1] != k) {
                            k = iwork[k - 1];
                            goto statement_130;
                        }
                    }
                    //
                    //                 Check the error code from Rsysvx.
                    //
                    if (info != k) {
                        Alaerh(path, "Rsysvx", info, k, fact + uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                        goto statement_150;
                    }
                    //
                    if (info == 0) {
                        if (ifact >= 2) {
                            //
                            //                       Reconstruct matrix from factors and compute
                            //                       residual.
                            //
                            Rsyt01(uplo, n, a, lda, afac, lda, iwork, ainv, lda, &rwork[(2 * nrhs + 1) - 1], result[1 - 1]);
                            k1 = 1;
                        } else {
                            k1 = 2;
                        }
                        //
                        //                    Compute residual of the computed solution.
                        //
                        Rlacpy("Full", n, nrhs, b, lda, work, lda);
                        Rpot02(uplo, n, nrhs, a, lda, x, lda, work, lda, &rwork[(2 * nrhs + 1) - 1], result[2 - 1]);
                        //
                        //                    Check solution from generated exact solution.
                        //
                        Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                        //
                        //                    Check the error bounds from iterative refinement.
                        //
                        Rpot05(uplo, n, nrhs, a, lda, b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], result[4 - 1]);
                    } else {
                        k1 = 6;
                    }
                    //
                    //                 Compare RCOND from Rsysvx with the computed value
                    //                 in RCONDC.
                    //
                    result[6 - 1] = Rget06[(rcond - 1) + (rcondc - 1) * ldRget06];
                    //
                    //                 Print information about the tests that did not pass
                    //                 the threshold.
                    //
                    for (k = k1; k <= 6; k = k + 1) {
                        if (result[k - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Aladhd(nout, path);
                            }
                            write(nout, format_9998), "Rsysvx", fact, uplo, n, imat, k, result(k);
                            nfail++;
                        }
                    }
                    nrun += 7 - k1;
                    //
                    //                 --- Test Rsysvxx ---
                    //
                    //                 Restore the matrices A and B.
                    //
                    if (ifact == 2) {
                        Rlaset(uplo, n, n, zero, zero, afac, lda);
                    }
                    Rlaset("Full", n, nrhs, zero, zero, x, lda);
                    //
                    //                 Solve the system and compute the condition number
                    //                 and error bounds using Rsysvxx.
                    //
                    n_err_bnds = 3;
                    equed = "N";
                    Rsysvxx(fact, uplo, n, nrhs, a, lda, afac, lda, iwork, equed, &work[(n + 1) - 1], b, lda, x, lda, rcond, rpvgrw_svxx, berr, n_err_bnds, errbnds_n, errbnds_c, 0, zero, work, &iwork[(n + 1) - 1], info);
                    //
                    //                 Adjust the expected value of INFO to account for
                    //                 pivoting.
                    //
                    k = izero;
                    if (k > 0) {
                    statement_135:
                        if (iwork[k - 1] < 0) {
                            if (iwork[k - 1] != -k) {
                                k = -iwork[k - 1];
                                goto statement_135;
                            }
                        } else if (iwork[k - 1] != k) {
                            k = iwork[k - 1];
                            goto statement_135;
                        }
                    }
                    //
                    //                 Check the error code from Rsysvxx.
                    //
                    if (info != k && info <= n) {
                        Alaerh(path, "Rsysvxx", info, k, fact + uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                        goto statement_150;
                    }
                    //
                    if (info == 0) {
                        if (ifact >= 2) {
                            //
                            //                 Reconstruct matrix from factors and compute
                            //                 residual.
                            //
                            Rsyt01(uplo, n, a, lda, afac, lda, iwork, ainv, lda, &rwork[(2 * nrhs + 1) - 1], result[1 - 1]);
                            k1 = 1;
                        } else {
                            k1 = 2;
                        }
                        //
                        //                 Compute residual of the computed solution.
                        //
                        Rlacpy("Full", n, nrhs, b, lda, work, lda);
                        Rpot02(uplo, n, nrhs, a, lda, x, lda, work, lda, &rwork[(2 * nrhs + 1) - 1], result[2 - 1]);
                        //
                        //                 Check solution from generated exact solution.
                        //
                        Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                        //
                        //                 Check the error bounds from iterative refinement.
                        //
                        Rpot05(uplo, n, nrhs, a, lda, b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], result[4 - 1]);
                    } else {
                        k1 = 6;
                    }
                    //
                    //                 Compare RCOND from Rsysvxx with the computed value
                    //                 in RCONDC.
                    //
                    result[6 - 1] = Rget06[(rcond - 1) + (rcondc - 1) * ldRget06];
                    //
                    //                 Print information about the tests that did not pass
                    //                 the threshold.
                    //
                    for (k = k1; k <= 6; k = k + 1) {
                        if (result[k - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Aladhd(nout, path);
                            }
                            write(nout, format_9998), "Rsysvxx", fact, uplo, n, imat, k, result(k);
                            nfail++;
                        }
                    }
                    nrun += 7 - k1;
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
    //     Test Error Bounds from Rsysvxx
    //
    Rebchvxx(thresh, path);
    //
    //     End of Rdrvsy
    //
}
