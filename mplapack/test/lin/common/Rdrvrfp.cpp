/*
 * Copyright (c) 2021-2022
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

void Rdrvrfp(INTEGER const nout, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, INTEGER const nnt, INTEGER *ntval, REAL const thresh, REAL *a, REAL *asav, REAL *afac, REAL *ainv, REAL *b, REAL *bsav, REAL *xact, REAL *x, REAL *arf, REAL *arfinv, REAL *d_work_Rlatms, REAL *d_work_Rpot01, REAL *d_temp_Rpot02, REAL *d_temp_Rpot03, REAL *d_work_Rlansy, REAL *d_work_Rpot02, REAL *d_work_Rpot03) {
    common cmn;
    common_write write(cmn);
    //
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    char uplos[] = {'U', 'L'};
    char forms[] = {'N', 'T'};
    char buf[1024];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    INTEGER iin = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    INTEGER ldb = 0;
    INTEGER iis = 0;
    INTEGER nrhs = 0;
    INTEGER iit = 0;
    INTEGER imat = 0;
    INTEGER iuplo = 0;
    char uplo;
    INTEGER iform = 0;
    char cform;
    char ctype;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist;
    INTEGER info = 0;
    bool zerot = false;
    INTEGER izero = 0;
    INTEGER ioff = 0;
    const REAL zero = 0.0;
    REAL rcondc = 0.0;
    REAL ainvnm = 0.0;
    const REAL one = 1.0;
    const INTEGER ntests = 4;
    REAL result[ntests];
    INTEGER nt = 0;
    INTEGER k = 0;
    //
    //     Initialize constants and the random number seed.
    //
    nrun = 0;
    nfail = 0;
    nerrs = 0;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    //
    for (iin = 1; iin <= nn; iin = iin + 1) {
        //
        n = nval[iin - 1];
        lda = max(n, (INTEGER)1);
        ldb = max(n, (INTEGER)1);
        //
        for (iis = 1; iis <= nns; iis = iis + 1) {
            //
            nrhs = nsval[iis - 1];
            //
            for (iit = 1; iit <= nnt; iit = iit + 1) {
                //
                imat = ntval[iit - 1];
                //
                //              If N.EQ.0, only consider the first type
                //
                if (n == 0 && iit >= 1) {
                    goto statement_120;
                }
                //
                //              Skip types 3, 4, or 5 if the matrix size is too small.
                //
                if (imat == 4 && n <= 1) {
                    goto statement_120;
                }
                if (imat == 5 && n <= 2) {
                    goto statement_120;
                }
                //
                //              Do first for UPLO = 'U', then for UPLO = 'L'
                //
                for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                    uplo = uplos[iuplo - 1];
                    //
                    //                 Do first for CFORM = 'N', then for CFORM = 'C'
                    //
                    for (iform = 1; iform <= 2; iform = iform + 1) {
                        cform = forms[iform - 1];
                        //
                        //                    Set up parameters with Rlatb4 and generate a test
                        //                    matrix with Rlatms.
                        //
                        Rlatb4("RPO", imat, n, n, &ctype, kl, ku, anorm, mode, cndnum, &dist);
                        //
                        Rlatms(n, n, &dist, iseed, &ctype, d_work_Rlatms, mode, cndnum, anorm, kl, ku, &uplo, a, lda, d_work_Rlatms, info);
                        //
                        //                    Check error code from Rlatms.
                        //
                        if (info != 0) {
                            Alaerh("RPF", "Rlatms", info, 0, &uplo, n, n, -1, -1, -1, iit, nfail, nerrs, nout);
                            goto statement_100;
                        }
                        //
                        //                    For types 3-5, zero one row and column of the matrix to
                        //                    test that INFO is returned correctly.
                        //
                        zerot = imat >= 3 && imat <= 5;
                        if (zerot) {
                            if (iit == 3) {
                                izero = 1;
                            } else if (iit == 4) {
                                izero = n;
                            } else {
                                izero = n / 2 + 1;
                            }
                            ioff = (izero - 1) * lda;
                            //
                            //                       Set row and column IZERO of A to 0.
                            //
                            if (iuplo == 1) {
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
                            izero = 0;
                        }
                        //
                        //                    Save a copy of the matrix A in ASAV.
                        //
                        Rlacpy(&uplo, n, n, a, lda, asav, lda);
                        //
                        //                    Compute the condition number of A (RCONDC).
                        //
                        if (zerot) {
                            rcondc = zero;
                        } else {
                            //
                            //                       Compute the 1-norm of A.
                            //
                            anorm = Rlansy("1", &uplo, n, a, lda, d_work_Rlansy);
                            //
                            //                       Factor the matrix A.
                            //
                            Rpotrf(&uplo, n, a, lda, info);
                            //
                            //                       Form the inverse of A.
                            //
                            Rpotri(&uplo, n, a, lda, info);
                            //
                            if (n != 0) {
                                //
                                //                          Compute the 1-norm condition number of A.
                                //
                                ainvnm = Rlansy("1", &uplo, n, a, lda, d_work_Rlansy);
                                rcondc = (one / anorm) / ainvnm;
                                //
                                //                          Restore the matrix A.
                                //
                                Rlacpy(&uplo, n, n, asav, lda, a, lda);
                            }
                            //
                        }
                        //
                        //                    Form an exact solution and set the right hand side.
                        //
                        Rlarhs("RPO", "N", &uplo, " ", n, n, kl, ku, nrhs, a, lda, xact, lda, b, lda, iseed, info);
                        Rlacpy("Full", n, nrhs, b, lda, bsav, lda);
                        //
                        //                    Compute the L*L' or U'*U factorization of the
                        //                    matrix and solve the system.
                        //
                        Rlacpy(&uplo, n, n, a, lda, afac, lda);
                        Rlacpy("Full", n, nrhs, b, ldb, x, ldb);
                        //
                        Rtrttf(&cform, &uplo, n, afac, lda, arf, info);
                        Rpftrf(&cform, &uplo, n, arf, info);
                        //
                        //                    Check error code from Rpftrf.
                        //
                        if (info != izero) {
                            //
                            //                       LANGOU: there is a small hick here: IZERO should
                            //                       always be INFO however if INFO is ZERO, Alaerh does not
                            //                       complain.
                            //
                            Alaerh("RPF", "Rpfsv", info, izero, &uplo, n, n, -1, -1, nrhs, iit, nfail, nerrs, nout);
                            goto statement_100;
                        }
                        //
                        //                    Skip the tests if INFO is not 0.
                        //
                        if (info != 0) {
                            goto statement_100;
                        }
                        //
                        Rpftrs(&cform, &uplo, n, nrhs, arf, x, ldb, info);
                        //
                        Rtfttr(&cform, &uplo, n, arf, afac, lda, info);
                        //
                        //                    Reconstruct matrix from factors and compute
                        //                    residual.
                        //
                        Rlacpy(&uplo, n, n, afac, lda, asav, lda);
                        Rpot01(&uplo, n, a, lda, afac, lda, d_work_Rpot01, result[1 - 1]);
                        Rlacpy(&uplo, n, n, asav, lda, afac, lda);
                        //
                        //                    Form the inverse and compute the residual.
                        //
                        if (mod(n, 2) == 0) {
                            Rlacpy("A", n + 1, n / 2, arf, n + 1, arfinv, n + 1);
                        } else {
                            Rlacpy("A", n, (n + 1) / 2, arf, n, arfinv, n);
                        }
                        //
                        Rpftri(&cform, &uplo, n, arfinv, info);
                        //
                        Rtfttr(&cform, &uplo, n, arfinv, ainv, lda, info);
                        //
                        //                    Check error code from Rpftri.
                        //
                        if (info != 0) {
                            Alaerh("RPO", "Rpftri", info, 0, &uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                        }
                        //
                        Rpot03(&uplo, n, a, lda, ainv, lda, d_temp_Rpot03, lda, d_work_Rpot03, rcondc, result[2 - 1]);
                        //
                        //                    Compute residual of the computed solution.
                        //
                        Rlacpy("Full", n, nrhs, b, lda, d_temp_Rpot02, lda);
                        Rpot02(&uplo, n, nrhs, a, lda, x, lda, d_temp_Rpot02, lda, d_work_Rpot02, result[3 - 1]);
                        //
                        //                    Check solution from generated exact solution.
                        //
                        Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[4 - 1]);
                        nt = 4;
                        //
                        //                    Print information about the tests that did not
                        //                    pass the threshold.
                        //
                        for (k = 1; k <= nt; k = k + 1) {
                            if (result[k - 1] >= thresh) {
                                if (nfail == 0 && nerrs == 0) {
                                    Aladhd(nout, "DPF");
                                }
                                sprintnum_short(buf, result[k - 1]);
                                write(nout, "(1x,a6,', UPLO=''',a1,''', N =',i5,', type ',i1,', test(',"
                                            "i1,')=',a)"),
                                    "Rpfsv", &uplo, n, iit, k, buf;
                                nfail++;
                            }
                        }
                        nrun += nt;
                    statement_100:;
                    }
                }
            statement_120:;
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasvm("RPF", nout, nfail, nrun, nerrs);
    //
    //     End of Rdrvrfp
    //
}
