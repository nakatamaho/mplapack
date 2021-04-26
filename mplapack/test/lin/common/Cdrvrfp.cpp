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

void Cdrvrfp(common &cmn, INTEGER const nout, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, INTEGER const nnt, INTEGER *ntval, REAL const thresh, COMPLEX *a, COMPLEX *asav, COMPLEX *afac, COMPLEX *ainv, COMPLEX *b, COMPLEX *bsav, COMPLEX *xact, COMPLEX *x, COMPLEX *arf, COMPLEX *arfinv, COMPLEX *z_work_zlatms, COMPLEX *z_work_Cpot02, COMPLEX *z_work_Cpot03, REAL *d_work_zlatms, REAL *d_work_zlanhe, REAL *d_work_Cpot01, REAL *d_work_Cpot02, REAL *d_work_Cpot03) {
    FEM_CMN_SVE(Cdrvrfp);
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
        {
            static const char *values[] = {"N", "C"};
            data_of_type_str(FEM_VALUES_AND_SIZE), forms;
        }
    }
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    arr_1d<4, int> iseed(fill0);
    INTEGER iin = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    INTEGER ldb = 0;
    INTEGER iis = 0;
    INTEGER nrhs = 0;
    INTEGER iit = 0;
    INTEGER imat = 0;
    INTEGER iuplo = 0;
    char uplo = char0;
    INTEGER iform = 0;
    char cform = char0;
    char ctype = char0;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist = char0;
    INTEGER info = 0;
    bool zerot = false;
    INTEGER izero = 0;
    INTEGER ioff = 0;
    const REAL zero = 0.0;
    REAL rcondc = 0.0;
    REAL ainvnm = 0.0;
    const REAL one = 1.0;
    const INTEGER ntests = 4;
    arr_1d<ntests, REAL> result(fill0);
    INTEGER nt = 0;
    INTEGER k = 0;
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
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
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
        lda = max(n, 1);
        ldb = max(n, 1);
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
                        //                    Set up parameters with Clatb4 and generate a test
                        //                    matrix with ZLATMS.
                        //
                        Clatb4("ZPO", imat, n, n, ctype, kl, ku, anorm, mode, cndnum, dist);
                        //
                        srnamt = "ZLATMS";
                        zlatms(n, n, dist, iseed, ctype, d_work_zlatms, mode, cndnum, anorm, kl, ku, uplo, a, lda, z_work_zlatms, info);
                        //
                        //                    Check error code from ZLATMS.
                        //
                        if (info != 0) {
                            Alaerh("ZPF", "ZLATMS", info, 0, uplo, n, n, -1, -1, -1, iit, nfail, nerrs, nout);
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
                        //                    Set the imaginary part of the diagonals.
                        //
                        Claipd(n, a, lda + 1, 0);
                        //
                        //                    Save a copy of the matrix A in ASAV.
                        //
                        zlacpy(uplo, n, n, a, lda, asav, lda);
                        //
                        //                    Compute the condition number of A (RCONDC).
                        //
                        if (zerot) {
                            rcondc = zero;
                        } else {
                            //
                            //                       Compute the 1-norm of A.
                            //
                            anorm = zlanhe("1", uplo, n, a, lda, d_work_zlanhe);
                            //
                            //                       Factor the matrix A.
                            //
                            zpotrf(uplo, n, a, lda, info);
                            //
                            //                       Form the inverse of A.
                            //
                            zpotri(uplo, n, a, lda, info);
                            //
                            if (n != 0) {
                                //
                                //                          Compute the 1-norm condition number of A.
                                //
                                ainvnm = zlanhe("1", uplo, n, a, lda, d_work_zlanhe);
                                rcondc = (one / anorm) / ainvnm;
                                //
                                //                          Restore the matrix A.
                                //
                                zlacpy(uplo, n, n, asav, lda, a, lda);
                            }
                            //
                        }
                        //
                        //                    Form an exact solution and set the right hand side.
                        //
                        srnamt = "Clarhs";
                        Clarhs("ZPO", "N", uplo, " ", n, n, kl, ku, nrhs, a, lda, xact, lda, b, lda, iseed, info);
                        zlacpy("Full", n, nrhs, b, lda, bsav, lda);
                        //
                        //                    Compute the L*L' or U'*U factorization of the
                        //                    matrix and solve the system.
                        //
                        zlacpy(uplo, n, n, a, lda, afac, lda);
                        zlacpy("Full", n, nrhs, b, ldb, x, ldb);
                        //
                        srnamt = "ZTRTTF";
                        ztrttf(cform, uplo, n, afac, lda, arf, info);
                        srnamt = "ZPFTRF";
                        zpftrf(cform, uplo, n, arf, info);
                        //
                        //                    Check error code from ZPFTRF.
                        //
                        if (info != izero) {
                            //
                            //                       LANGOU: there is a small hick here: IZERO should
                            //                       always be INFO however if INFO is ZERO, Alaerh does not
                            //                       complain.
                            //
                            Alaerh("ZPF", "ZPFSV ", info, izero, uplo, n, n, -1, -1, nrhs, iit, nfail, nerrs, nout);
                            goto statement_100;
                        }
                        //
                        //                     Skip the tests if INFO is not 0.
                        //
                        if (info != 0) {
                            goto statement_100;
                        }
                        //
                        srnamt = "ZPFTRS";
                        zpftrs(cform, uplo, n, nrhs, arf, x, ldb, info);
                        //
                        srnamt = "ZTFTTR";
                        ztfttr(cform, uplo, n, arf, afac, lda, info);
                        //
                        //                    Reconstruct matrix from factors and compute
                        //                    residual.
                        //
                        zlacpy(uplo, n, n, afac, lda, asav, lda);
                        Cpot01(uplo, n, a, lda, afac, lda, d_work_Cpot01, result[1 - 1]);
                        zlacpy(uplo, n, n, asav, lda, afac, lda);
                        //
                        //                    Form the inverse and compute the residual.
                        //
                        if (mod(n, 2) == 0) {
                            zlacpy("A", n + 1, n / 2, arf, n + 1, arfinv, n + 1);
                        } else {
                            zlacpy("A", n, (n + 1) / 2, arf, n, arfinv, n);
                        }
                        //
                        srnamt = "ZPFTRI";
                        zpftri(cform, uplo, n, arfinv, info);
                        //
                        srnamt = "ZTFTTR";
                        ztfttr(cform, uplo, n, arfinv, ainv, lda, info);
                        //
                        //                    Check error code from ZPFTRI.
                        //
                        if (info != 0) {
                            Alaerh("ZPO", "ZPFTRI", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                        }
                        //
                        Cpot03(uplo, n, a, lda, ainv, lda, z_work_Cpot03, lda, d_work_Cpot03, rcondc, result[2 - 1]);
                        //
                        //                    Compute residual of the computed solution.
                        //
                        zlacpy("Full", n, nrhs, b, lda, z_work_Cpot02, lda);
                        Cpot02(uplo, n, nrhs, a, lda, x, lda, z_work_Cpot02, lda, d_work_Cpot02, result[3 - 1]);
                        //
                        //                    Check solution from generated exact solution.
                        //
                        Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[4 - 1]);
                        nt = 4;
                        //
                        //                    PrINTEGER information about the tests that did not
                        //                    pass the threshold.
                        //
                        for (k = 1; k <= nt; k = k + 1) {
                            if (result[k - 1] >= thresh) {
                                if (nfail == 0 && nerrs == 0) {
                                    Aladhd(nout, "ZPF");
                                }
                                write(nout, "(1x,a6,', UPLO=''',a1,''', N =',i5,', type ',i1,', test(',"
                                            "i1,')=',g12.5)"),
                                    "ZPFSV ", uplo, n, iit, k, result(k);
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
    //     PrINTEGER a summary of the results.
    //
    Alasvm("ZPF", nout, nfail, nrun, nerrs);
    //
    //     End of Cdrvrfp
    //
}
