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

void Cchktb(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, REAL const thresh, bool const tsterr, INTEGER const  /* nmax */, COMPLEX *ab, COMPLEX *ainv, COMPLEX *b, COMPLEX *x, COMPLEX *xact, COMPLEX *work, REAL *rwork, INTEGER const nout) {
    FEM_CMN_SVE(Cchktb);
    common_write write(cmn);
    char &srnamt = cmn.srnamt;
    //
    const INTEGER ntran = 3;
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
            static const char *values[] = {"N", "T", "C"};
            data_of_type_str(FEM_VALUES_AND_SIZE), transs;
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
    const INTEGER ntype1 = 9;
    INTEGER nimat = 0;
    const INTEGER ntypes = 17;
    INTEGER nimat2 = 0;
    INTEGER nk = 0;
    INTEGER ik = 0;
    INTEGER kd = 0;
    INTEGER ldab = 0;
    INTEGER imat = 0;
    INTEGER iuplo = 0;
    char uplo[1];
    char diag[1];
    INTEGER info = 0;
    INTEGER idiag = 0;
    const REAL zero = 0.0;
    const REAL one = 1.0;
    INTEGER j = 0;
    REAL anorm = 0.0;
    REAL ainvnm = 0.0;
    REAL rcondo = 0.0;
    REAL rcondi = 0.0;
    INTEGER irhs = 0;
    INTEGER nrhs = 0;
    INTEGER itran = 0;
    char trans[1];
    char norm[1];
    REAL rcondc = 0.0;
    const INTEGER ntests = 8;
    arr_1d<ntests, REAL> result;
    INTEGER k = 0;
    REAL rcond = 0.0;
    REAL scale = 0.0;
    static const char *format_9997 = "(1x,a,'( ''',a1,''', ''',a1,''', ''',a1,''', ''',a1,''',',i5,',',i5,"
                                     "', ...  ),  type ',i2,', test(',i1,')=',g12.5)";
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
    path[(2 - 1) + (3 - 1) * ldpath] = "TB";
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
        Cerrtr(path, nout);
    }
    cmn.infot = 0;
    //
    for (in = 1; in <= nn; in = in + 1) {
        //
        //        Do for each value of N in NVAL
        //
        n = nval[in - 1];
        lda = max((INTEGER)1, n);
        xtype = "N";
        nimat = ntype1;
        nimat2 = ntypes;
        if (n <= 0) {
            nimat = 1;
            nimat2 = ntype1 + 1;
        }
        //
        nk = min(n + 1, 4);
        for (ik = 1; ik <= nk; ik = ik + 1) {
            //
            //           Do for KD = 0, N, (3N-1)/4, and (N+1)/4. This order makes
            //           it easier to skip redundant values for small values of N.
            //
            if (ik == 1) {
                kd = 0;
            } else if (ik == 2) {
                kd = max(n, 0);
            } else if (ik == 3) {
                kd = (3 * n - 1) / 4;
            } else if (ik == 4) {
                kd = (n + 1) / 4;
            }
            ldab = kd + 1;
            //
            for (imat = 1; imat <= nimat; imat = imat + 1) {
                //
                //              Do the tests only if DOTYPE( IMAT ) is true.
                //
                if (!dotype[imat - 1]) {
                    goto statement_90;
                }
                //
                for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                    //
                    //                 Do first for UPLO = 'U', then for UPLO = 'L'
                    //
                    uplo = uplos[iuplo - 1];
                    //
                    //                 Call Clattb to generate a triangular test matrix.
                    //
                    srnamt = "Clattb";
                    Clattb(imat, uplo, "No transpose", diag, iseed, n, kd, ab, ldab, x, work, rwork, info);
                    //
                    //                 Set IDIAG = 1 for non-unit matrices, 2 for unit.
                    //
                    if (Mlsame(diag, "N")) {
                        idiag = 1;
                    } else {
                        idiag = 2;
                    }
                    //
                    //                 Form the inverse of A so we can get a good estimate
                    //                 of RCONDC = 1/(norm(A) * norm(inv(A))).
                    //
                    Claset("Full", n, n, COMPLEX(zero), COMPLEX(one), ainv, lda);
                    if (Mlsame(uplo, "U")) {
                        for (j = 1; j <= n; j = j + 1) {
                            Ctbsv(uplo, "No transpose", diag, j, kd, ab, ldab, ainv[((j - 1) * lda + 1) - 1], 1);
                        }
                    } else {
                        for (j = 1; j <= n; j = j + 1) {
                            Ctbsv(uplo, "No transpose", diag, n - j + 1, kd, &ab[((j - 1) * ldab + 1) - 1], ldab, ainv[((j - 1) * lda + j) - 1], 1);
                        }
                    }
                    //
                    //                 Compute the 1-norm condition number of A.
                    //
                    anorm = Clantb("1", uplo, diag, n, kd, ab, ldab, rwork);
                    ainvnm = Clantr("1", uplo, diag, n, n, ainv, lda, rwork);
                    if (anorm <= zero || ainvnm <= zero) {
                        rcondo = one;
                    } else {
                        rcondo = (one / anorm) / ainvnm;
                    }
                    //
                    //                 Compute the infinity-norm condition number of A.
                    //
                    anorm = Clantb("I", uplo, diag, n, kd, ab, ldab, rwork);
                    ainvnm = Clantr("I", uplo, diag, n, n, ainv, lda, rwork);
                    if (anorm <= zero || ainvnm <= zero) {
                        rcondi = one;
                    } else {
                        rcondi = (one / anorm) / ainvnm;
                    }
                    //
                    for (irhs = 1; irhs <= nns; irhs = irhs + 1) {
                        nrhs = nsval[irhs - 1];
                        xtype = "N";
                        //
                        for (itran = 1; itran <= ntran; itran = itran + 1) {
                            //
                            //                    Do for op(A) = A, A**T, or A**H.
                            //
                            trans = transs[itran - 1];
                            if (itran == 1) {
                                norm = "O";
                                rcondc = rcondo;
                            } else {
                                norm = 'I';
                                rcondc = rcondi;
                            }
                            //
                            //+    TEST 1
                            //                    Solve and compute residual for op(A)*x = b.
                            //
                            srnamt = "Clarhs";
                            Clarhs(path, xtype, uplo, trans, n, n, kd, idiag, nrhs, ab, ldab, xact, lda, b, lda, iseed, info);
                            xtype = "C";
                            Clacpy("Full", n, nrhs, b, lda, x, lda);
                            //
                            srnamt = "Ctbtrs";
                            Ctbtrs(uplo, trans, diag, n, kd, nrhs, ab, ldab, x, lda, info);
                            //
                            //                    Check error code from Ctbtrs.
                            //
                            if (info != 0) {
                                Alaerh(path, "Ctbtrs", info, 0, uplo + trans + diag, n, n, kd, kd, nrhs, imat, nfail, nerrs, nout);
                            }
                            //
                            Ctbt02(uplo, trans, diag, n, kd, nrhs, ab, ldab, x, lda, b, lda, work, rwork, result[1 - 1]);
                            //
                            //+    TEST 2
                            //                    Check solution from generated exact solution.
                            //
                            Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[2 - 1]);
                            //
                            //+    TESTS 3, 4, and 5
                            //                    Use iterative refinement to improve the solution
                            //                    and compute error bounds.
                            //
                            srnamt = "Ctbrfs";
                            Ctbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, lda, x, lda, rwork, &rwork[(nrhs + 1) - 1], work, &rwork[(2 * nrhs + 1) - 1], info);
                            //
                            //                    Check error code from Ctbrfs.
                            //
                            if (info != 0) {
                                Alaerh(path, "Ctbrfs", info, 0, uplo + trans + diag, n, n, kd, kd, nrhs, imat, nfail, nerrs, nout);
                            }
                            //
                            Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                            Ctbt05(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], result[4 - 1]);
                            //
                            //                       Print information about the tests that did not
                            //                       pass the threshold.
                            //
                            for (k = 1; k <= 5; k = k + 1) {
                                if (result[k - 1] >= thresh) {
                                    if (nfail == 0 && nerrs == 0) {
                                        Alahd(nout, path);
                                    }
                                    write(nout, "(' UPLO=''',a1,''', TRANS=''',a1,''',      DIAG=''',a1,"
                                                "''', N=',i5,', KD=',i5,', NRHS=',i5,', type ',i2,"
                                                "', test(',i2,')=',g12.5)"),
                                        uplo, trans, diag, n, kd, nrhs, imat, k, result(k);
                                    nfail++;
                                }
                            }
                            nrun += 5;
                        }
                    }
                    //
                    //+    TEST 6
                    //                    Get an estimate of RCOND = 1/CNDNUM.
                    //
                    for (itran = 1; itran <= 2; itran = itran + 1) {
                        if (itran == 1) {
                            norm = "O";
                            rcondc = rcondo;
                        } else {
                            norm = 'I';
                            rcondc = rcondi;
                        }
                        srnamt = "Ctbcon";
                        Ctbcon(norm, uplo, diag, n, kd, ab, ldab, rcond, work, rwork, info);
                        //
                        //                    Check error code from Ctbcon.
                        //
                        if (info != 0) {
                            Alaerh(path, "Ctbcon", info, 0, norm + uplo + diag, n, n, kd, kd, -1, imat, nfail, nerrs, nout);
                        }
                        //
                        Ctbt06(rcond, rcondc, uplo, diag, n, kd, ab, ldab, rwork, result[6 - 1]);
                        //
                        //                    Print the test ratio if it is .GE. THRESH.
                        //
                        if (result[6 - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Alahd(nout, path);
                            }
                            write(nout, "(1x,a,'( ''',a1,''', ''',a1,''', ''',a1,''',',i5,',',i5,"
                                        "',  ... ), type ',i2,', test(',i2,')=',g12.5)"),
                                "Ctbcon", norm, uplo, diag, n, kd, imat, 6, result(6);
                            nfail++;
                        }
                        nrun++;
                    }
                }
            statement_90:;
            }
            //
            //           Use pathological test matrices to test Clatbs.
            //
            for (imat = ntype1 + 1; imat <= nimat2; imat = imat + 1) {
                //
                //              Do the tests only if DOTYPE( IMAT ) is true.
                //
                if (!dotype[imat - 1]) {
                    goto statement_120;
                }
                //
                for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                    //
                    //                 Do first for UPLO = 'U', then for UPLO = 'L'
                    //
                    uplo = uplos[iuplo - 1];
                    for (itran = 1; itran <= ntran; itran = itran + 1) {
                        //
                        //                    Do for op(A) = A, A**T, and A**H.
                        //
                        trans = transs[itran - 1];
                        //
                        //                    Call Clattb to generate a triangular test matrix.
                        //
                        srnamt = "Clattb";
                        Clattb(imat, uplo, trans, diag, iseed, n, kd, ab, ldab, x, work, rwork, info);
                        //
                        //+    TEST 7
                        //                    Solve the system op(A)*x = b
                        //
                        srnamt = "Clatbs";
                        Ccopy(n, x, 1, b, 1);
                        Clatbs(uplo, trans, diag, "N", n, kd, ab, ldab, b, scale, rwork, info);
                        //
                        //                    Check error code from Clatbs.
                        //
                        if (info != 0) {
                            Alaerh(path, "Clatbs", info, 0, uplo + trans + diag + const char *("N"), n, n, kd, kd, -1, imat, nfail, nerrs, nout);
                        }
                        //
                        Ctbt03(uplo, trans, diag, n, kd, 1, ab, ldab, scale, rwork, one, b, lda, x, lda, work, result[7 - 1]);
                        //
                        //+    TEST 8
                        //                    Solve op(A)*x = b again with NORMIN = 'Y'.
                        //
                        Ccopy(n, x, 1, b, 1);
                        Clatbs(uplo, trans, diag, "Y", n, kd, ab, ldab, b, scale, rwork, info);
                        //
                        //                    Check error code from Clatbs.
                        //
                        if (info != 0) {
                            Alaerh(path, "Clatbs", info, 0, uplo + trans + diag + const char *("Y"), n, n, kd, kd, -1, imat, nfail, nerrs, nout);
                        }
                        //
                        Ctbt03(uplo, trans, diag, n, kd, 1, ab, ldab, scale, rwork, one, b, lda, x, lda, work, result[8 - 1]);
                        //
                        //                    Print information about the tests that did not pass
                        //                    the threshold.
                        //
                        if (result[7 - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Alahd(nout, path);
                            }
                            write(nout, format_9997), "Clatbs", uplo, trans, diag, "N", n, kd, imat, 7, result(7);
                            nfail++;
                        }
                        if (result[8 - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Alahd(nout, path);
                            }
                            write(nout, format_9997), "Clatbs", uplo, trans, diag, "Y", n, kd, imat, 8, result(8);
                            nfail++;
                        }
                        nrun += 2;
                    }
                }
            statement_120:;
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Cchktb
    //
}
