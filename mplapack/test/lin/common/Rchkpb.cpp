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

void Rchkpb(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nnb, INTEGER *nbval, INTEGER const nns, INTEGER *nsval, REAL const thresh, bool const tsterr, INTEGER const  /* nmax */, REAL *a, REAL *afac, REAL *ainv, REAL *b, REAL *x, REAL *xact, REAL *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
    FEM_CMN_SVE(Rchkpb);
    common_write write(cmn);
    //
    INTEGER *iseedy(sve.iseedy, [4]);
    if (is_called_first_time) {
        static const INTEGER values[] = {1988, 1989, 1990, 1991};
        data_of_type<int>(FEM_VALUES_AND_SIZE), iseedy;
    }
    char path[3];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    const INTEGER nbw = 4;
    INTEGER kdval[nbw];
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    char xtype;
    INTEGER nkd = 0;
    const INTEGER ntypes = 8;
    INTEGER nimat = 0;
    INTEGER ikd = 0;
    INTEGER kd = 0;
    INTEGER ldab = 0;
    INTEGER iuplo = 0;
    INTEGER koff = 0;
    char uplo;
    char packit;
    INTEGER imat = 0;
    bool zerot = false;
    char type;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist;
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER iw = 0;
    INTEGER ioff = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    const REAL zero = 0.0;
    INTEGER inb = 0;
    INTEGER nb = 0;
    const INTEGER ntests = 7;
    REAL result[ntests];
    const REAL one = 1.0;
    REAL ainvnm = 0.0;
    REAL rcondc = 0.0;
    INTEGER irhs = 0;
    INTEGER nrhs = 0;
    INTEGER k = 0;
    REAL rcond = 0.0;
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
    path[(1 - 1)] = "Double precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "PB";
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
        Rerrpo(path, nout);
    }
    cmn.infot = 0;
    xlaenv(2, 2);
    kdval[1 - 1] = 0;
    //
    //     Do for each value of N in NVAL
    //
    for (in = 1; in <= nn; in = in + 1) {
        n = nval[in - 1];
        lda = max(n, 1);
        xtype = "N";
        //
        //        Set limits on the number of loop iterations.
        //
        nkd = max({(INTEGER)1, min(n, 4)});
        nimat = ntypes;
        if (n == 0) {
            nimat = 1;
        }
        //
        kdval[2 - 1] = n + (n + 1) / 4;
        kdval[3 - 1] = (3 * n - 1) / 4;
        kdval[4 - 1] = (n + 1) / 4;
        //
        for (ikd = 1; ikd <= nkd; ikd = ikd + 1) {
            //
            //           Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
            //           makes it easier to skip redundant values for small values
            //           of N.
            //
            kd = kdval[ikd - 1];
            ldab = kd + 1;
            //
            //           Do first for UPLO = 'U', then for UPLO = 'L'
            //
            for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                koff = 1;
                if (iuplo == 1) {
                    uplo = "U";
                    koff = max((INTEGER)1, kd + 2 - n);
                    packit = "Q";
                } else {
                    uplo = "L";
                    packit = "B";
                }
                //
                for (imat = 1; imat <= nimat; imat = imat + 1) {
                    //
                    //                 Do the tests only if DOTYPE( IMAT ) is true.
                    //
                    if (!dotype[imat - 1]) {
                        goto statement_60;
                    }
                    //
                    //                 Skip types 2, 3, or 4 if the matrix size is too small.
                    //
                    zerot = imat >= 2 && imat <= 4;
                    if (zerot && n < imat - 1) {
                        goto statement_60;
                    }
                    //
                    if (!zerot || !dotype[1 - 1]) {
                        //
                        //                    Set up parameters with Rlatb4 and generate a test
                        //                    matrix with Rlatms.
                        //
                        Rlatb4(path, imat, n, n, type, kl, ku, anorm, mode, cndnum, dist);
                        //
                        Rlatms(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, kd, kd, packit, &a[koff - 1], ldab, work, info);
                        //
                        //                    Check error code from Rlatms.
                        //
                        if (info != 0) {
                            Alaerh(path, "Rlatms", info, 0, uplo, n, n, kd, kd, -1, imat, nfail, nerrs, nout);
                            goto statement_60;
                        }
                    } else if (izero > 0) {
                        //
                        //                    Use the same matrix for types 3 and 4 as for type
                        //                    2 by copying back the zeroed out column,
                        //
                        iw = 2 * lda + 1;
                        if (iuplo == 1) {
                            ioff = (izero - 1) * ldab + kd + 1;
                            Rcopy(izero - i1, &work[iw - 1], 1, &a[(ioff - izero + i1) - 1], 1);
                            iw += izero - i1;
                            Rcopy(i2 - izero + 1, &work[iw - 1], 1, &a[ioff - 1], max(ldab - 1, 1));
                        } else {
                            ioff = (i1 - 1) * ldab + 1;
                            Rcopy(izero - i1, &work[iw - 1], 1, &a[(ioff + izero - i1) - 1], max(ldab - 1, 1));
                            ioff = (izero - 1) * ldab + 1;
                            iw += izero - i1;
                            Rcopy(i2 - izero + 1, &work[iw - 1], 1, &a[ioff - 1], 1);
                        }
                    }
                    //
                    //                 For types 2-4, zero one row and column of the matrix
                    //                 to test that INFO is returned correctly.
                    //
                    izero = 0;
                    if (zerot) {
                        if (imat == 2) {
                            izero = 1;
                        } else if (imat == 3) {
                            izero = n;
                        } else {
                            izero = n / 2 + 1;
                        }
                        //
                        //                    Save the zeroed out row and column in WORK(*,3)
                        //
                        iw = 2 * lda;
                        for (i = 1; i <= min(2 * kd + 1, n); i = i + 1) {
                            work[(iw + i) - 1] = zero;
                        }
                        iw++;
                        i1 = max(izero - kd, 1);
                        i2 = min(izero + kd, n);
                        //
                        if (iuplo == 1) {
                            ioff = (izero - 1) * ldab + kd + 1;
                            Rswap(izero - i1, &a[(ioff - izero + i1) - 1], 1, &work[iw - 1], 1);
                            iw += izero - i1;
                            Rswap(i2 - izero + 1, &a[ioff - 1], max(ldab - 1, 1), &work[iw - 1], 1);
                        } else {
                            ioff = (i1 - 1) * ldab + 1;
                            Rswap(izero - i1, &a[(ioff + izero - i1) - 1], max(ldab - 1, 1), &work[iw - 1], 1);
                            ioff = (izero - 1) * ldab + 1;
                            iw += izero - i1;
                            Rswap(i2 - izero + 1, &a[ioff - 1], 1, &work[iw - 1], 1);
                        }
                    }
                    //
                    //                 Do for each value of NB in NBVAL
                    //
                    for (inb = 1; inb <= nnb; inb = inb + 1) {
                        nb = nbval[inb - 1];
                        xlaenv(1, nb);
                        //
                        //                    Compute the L*L' or U'*U factorization of the band
                        //                    matrix.
                        //
                        Rlacpy("Full", kd + 1, n, a, ldab, afac, ldab);
                        Rpbtrf(uplo, n, kd, afac, ldab, info);
                        //
                        //                    Check error code from Rpbtrf.
                        //
                        if (info != izero) {
                            Alaerh(path, "Rpbtrf", info, izero, uplo, n, n, kd, kd, nb, imat, nfail, nerrs, nout);
                            goto statement_50;
                        }
                        //
                        //                    Skip the tests if INFO is not 0.
                        //
                        if (info != 0) {
                            goto statement_50;
                        }
                        //
                        //+    TEST 1
                        //                    Reconstruct matrix from factors and compute
                        //                    residual.
                        //
                        Rlacpy("Full", kd + 1, n, afac, ldab, ainv, ldab);
                        Rpbt01(uplo, n, kd, a, ldab, ainv, ldab, rwork, result[1 - 1]);
                        //
                        //                    Print the test ratio if it is .GE. THRESH.
                        //
                        if (result[1 - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Alahd(nout, path);
                            }
                            write(nout, "(' UPLO=''',a1,''', N=',i5,', KD=',i5,', NB=',i4,', type ',"
                                        "i2,', test ',i2,', ratio= ',g12.5)"),
                                uplo, n, kd, nb, imat, 1, result(1);
                            nfail++;
                        }
                        nrun++;
                        //
                        //                    Only do other tests if this is the first blocksize.
                        //
                        if (inb > 1) {
                            goto statement_50;
                        }
                        //
                        //                    Form the inverse of A so we can get a good estimate
                        //                    of RCONDC = 1/(norm(A) * norm(inv(A))).
                        //
                        Rlaset("Full", n, n, zero, one, ainv, lda);
                        Rpbtrs(uplo, n, kd, n, afac, ldab, ainv, lda, info);
                        //
                        //                    Compute RCONDC = 1/(norm(A) * norm(inv(A))).
                        //
                        anorm = Rlansb("1", uplo, n, kd, a, ldab, rwork);
                        ainvnm = Rlange("1", n, n, ainv, lda, rwork);
                        if (anorm <= zero || ainvnm <= zero) {
                            rcondc = one;
                        } else {
                            rcondc = (one / anorm) / ainvnm;
                        }
                        //
                        for (irhs = 1; irhs <= nns; irhs = irhs + 1) {
                            nrhs = nsval[irhs - 1];
                            //
                            //+    TEST 2
                            //                    Solve and compute residual for A * X = B.
                            //
                            Rlarhs(path, xtype, uplo, " ", n, n, kd, kd, nrhs, a, ldab, xact, lda, b, lda, iseed, info);
                            Rlacpy("Full", n, nrhs, b, lda, x, lda);
                            //
                            Rpbtrs(uplo, n, kd, nrhs, afac, ldab, x, lda, info);
                            //
                            //                    Check error code from Rpbtrs.
                            //
                            if (info != 0) {
                                Alaerh(path, "Rpbtrs", info, 0, uplo, n, n, kd, kd, nrhs, imat, nfail, nerrs, nout);
                            }
                            //
                            Rlacpy("Full", n, nrhs, b, lda, work, lda);
                            Rpbt02(uplo, n, kd, nrhs, a, ldab, x, lda, work, lda, rwork, result[2 - 1]);
                            //
                            //+    TEST 3
                            //                    Check solution from generated exact solution.
                            //
                            Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                            //
                            //+    TESTS 4, 5, and 6
                            //                    Use iterative refinement to improve the solution.
                            //
                            Rpbrfs(uplo, n, kd, nrhs, a, ldab, afac, ldab, b, lda, x, lda, rwork, &rwork[(nrhs + 1) - 1], work, iwork, info);
                            //
                            //                    Check error code from Rpbrfs.
                            //
                            if (info != 0) {
                                Alaerh(path, "Rpbrfs", info, 0, uplo, n, n, kd, kd, nrhs, imat, nfail, nerrs, nout);
                            }
                            //
                            Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[4 - 1]);
                            Rpbt05(uplo, n, kd, nrhs, a, ldab, b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], result[5 - 1]);
                            //
                            //                       Print information about the tests that did not
                            //                       pass the threshold.
                            //
                            for (k = 2; k <= 6; k = k + 1) {
                                if (result[k - 1] >= thresh) {
                                    if (nfail == 0 && nerrs == 0) {
                                        Alahd(nout, path);
                                    }
                                    write(nout, "(' UPLO=''',a1,''', N=',i5,', KD=',i5,', NRHS=',i3,"
                                                "', type ',i2,', test(',i2,') = ',g12.5)"),
                                        uplo, n, kd, nrhs, imat, k, result(k);
                                    nfail++;
                                }
                            }
                            nrun += 5;
                        }
                        //
                        //+    TEST 7
                        //                    Get an estimate of RCOND = 1/CNDNUM.
                        //
                        Rpbcon(uplo, n, kd, afac, ldab, anorm, rcond, work, iwork, info);
                        //
                        //                    Check error code from Rpbcon.
                        //
                        if (info != 0) {
                            Alaerh(path, "Rpbcon", info, 0, uplo, n, n, kd, kd, -1, imat, nfail, nerrs, nout);
                        }
                        //
                        result[7 - 1] = Rget06[(rcond - 1) + (rcondc - 1) * ldRget06];
                        //
                        //                    Print the test ratio if it is .GE. THRESH.
                        //
                        if (result[7 - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Alahd(nout, path);
                            }
                            write(nout, "(' UPLO=''',a1,''', N=',i5,', KD=',i5,',',10x,' type ',i2,"
                                        "', test(',i2,') = ',g12.5)"),
                                uplo, n, kd, imat, 7, result(7);
                            nfail++;
                        }
                        nrun++;
                    statement_50:;
                    }
                statement_60:;
                }
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Rchkpb
    //
}
