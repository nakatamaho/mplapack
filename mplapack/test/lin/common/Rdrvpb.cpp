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

void Rdrvpb(common &cmn, bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nrhs, REAL const thresh, bool const tsterr, INTEGER const /* nmax */, REAL *a, REAL *afac, REAL *asav, REAL *b, REAL *bsav, REAL *x, REAL *xact, REAL *s, REAL *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
    FEM_CMN_SVE(Rdrvpb);
    common_write write(cmn);
    str<32> &srnamt = cmn.srnamt;
    //
    if (is_called_first_time) {
        {
            static const INTEGER values[] = {1988, 1989, 1990, 1991};
            data_of_type<int>(FEM_VALUES_AND_SIZE), iseedy;
        }
        {
            static const char *values[] = {"F", "N", "E"};
            data_of_type_str(FEM_VALUES_AND_SIZE), facts;
        }
        {
            static const char *values[] = {"N", "Y"};
            data_of_type_str(FEM_VALUES_AND_SIZE), equeds;
        }
    }
    str<3> path = char0;
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    arr_1d<4, int> iseed(fill0);
    const INTEGER nbw = 4;
    arr_1d<nbw, int> kdval(fill0);
    INTEGER nb = 0;
    INTEGER nbmin = 0;
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    char xtype = char0;
    INTEGER nkd = 0;
    const INTEGER ntypes = 8;
    INTEGER nimat = 0;
    INTEGER ikd = 0;
    INTEGER kd = 0;
    INTEGER ldab = 0;
    INTEGER iuplo = 0;
    INTEGER koff = 0;
    char uplo = char0;
    char packit = char0;
    INTEGER imat = 0;
    bool zerot = false;
    char type = char0;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist = char0;
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER iw = 0;
    INTEGER ioff = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    const REAL zero = 0.0;
    INTEGER iequed = 0;
    char equed = char0;
    INTEGER nfact = 0;
    INTEGER ifact = 0;
    char fact = char0;
    bool prefac = false;
    bool nofact = false;
    bool equil = false;
    REAL rcondc = 0.0;
    REAL scond = 0.0;
    REAL amax = 0.0;
    REAL roldc = 0.0;
    const REAL one = 1.0;
    REAL ainvnm = 0.0;
    const INTEGER ntests = 6;
    arr_1d<ntests, REAL> result(fill0);
    INTEGER nt = 0;
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
        Rerrvx(path, nout);
    }
    cmn.infot = 0;
    kdval[1 - 1] = 0;
    //
    //     Set the block size and minimum block size for testing.
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
        //
        //        Set limits on the number of loop iterations.
        //
        nkd = max((INTEGER)1, min(n, 4));
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
                    packit = "Q";
                    koff = max((INTEGER)1, kd + 2 - n);
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
                        goto statement_80;
                    }
                    //
                    //                 Skip types 2, 3, or 4 if the matrix size is too small.
                    //
                    zerot = imat >= 2 && imat <= 4;
                    if (zerot && n < imat - 1) {
                        goto statement_80;
                    }
                    //
                    if (!zerot || !dotype[1 - 1]) {
                        //
                        //                    Set up parameters with Rlatb4 and generate a test
                        //                    matrix with DLATMS.
                        //
                        Rlatb4(path, imat, n, n, type, kl, ku, anorm, mode, cndnum, dist);
                        //
                        srnamt = "DLATMS";
                        dlatms(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, kd, kd, packit, &a[koff - 1], ldab, work, info);
                        //
                        //                    Check error code from DLATMS.
                        //
                        if (info != 0) {
                            Alaerh(path, "DLATMS", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                            goto statement_80;
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
                    //                 Save a copy of the matrix A in ASAV.
                    //
                    Rlacpy("Full", kd + 1, n, a, ldab, asav, ldab);
                    //
                    for (iequed = 1; iequed <= 2; iequed = iequed + 1) {
                        equed = equeds[iequed - 1];
                        if (iequed == 1) {
                            nfact = 3;
                        } else {
                            nfact = 1;
                        }
                        //
                        for (ifact = 1; ifact <= nfact; ifact = ifact + 1) {
                            fact = facts[ifact - 1];
                            prefac = Mlsame(fact, "F");
                            nofact = Mlsame(fact, "N");
                            equil = Mlsame(fact, "E");
                            //
                            if (zerot) {
                                if (prefac) {
                                    goto statement_60;
                                }
                                rcondc = zero;
                                //
                            } else if (!Mlsame(fact, "N")) {
                                //
                                //                          Compute the condition number for comparison
                                //                          with the value returned by Rpbsvx (FACT =
                                //                          'N' reuses the condition number from the
                                //                          previous iteration with FACT = 'F').
                                //
                                Rlacpy("Full", kd + 1, n, asav, ldab, afac, ldab);
                                if (equil || iequed > 1) {
                                    //
                                    //                             Compute row and column scale factors to
                                    //                             equilibrate the matrix A.
                                    //
                                    Rpbequ(uplo, n, kd, afac, ldab, s, scond, amax, info);
                                    if (info == 0 && n > 0) {
                                        if (iequed > 1) {
                                            scond = zero;
                                        }
                                        //
                                        //                                Equilibrate the matrix.
                                        //
                                        Rlaqsb(uplo, n, kd, afac, ldab, s, scond, amax, equed);
                                    }
                                }
                                //
                                //                          Save the condition number of the
                                //                          non-equilibrated system for use in Rget04.
                                //
                                if (equil) {
                                    roldc = rcondc;
                                }
                                //
                                //                          Compute the 1-norm of A.
                                //
                                anorm = Rlansb("1", uplo, n, kd, afac, ldab, rwork);
                                //
                                //                          Factor the matrix A.
                                //
                                Rpbtrf(uplo, n, kd, afac, ldab, info);
                                //
                                //                          Form the inverse of A.
                                //
                                Rlaset("Full", n, n, zero, one, a, lda);
                                srnamt = "Rpbtrs";
                                Rpbtrs(uplo, n, kd, n, afac, ldab, a, lda, info);
                                //
                                //                          Compute the 1-norm condition number of A.
                                //
                                ainvnm = Rlange("1", n, n, a, lda, rwork);
                                if (anorm <= zero || ainvnm <= zero) {
                                    rcondc = one;
                                } else {
                                    rcondc = (one / anorm) / ainvnm;
                                }
                            }
                            //
                            //                       Restore the matrix A.
                            //
                            Rlacpy("Full", kd + 1, n, asav, ldab, a, ldab);
                            //
                            //                       Form an exact solution and set the right hand
                            //                       side.
                            //
                            srnamt = "Rlarhs";
                            Rlarhs(path, xtype, uplo, " ", n, n, kd, kd, nrhs, a, ldab, xact, lda, b, lda, iseed, info);
                            xtype = "C";
                            Rlacpy("Full", n, nrhs, b, lda, bsav, lda);
                            //
                            if (nofact) {
                                //
                                //                          --- Test Rpbsv  ---
                                //
                                //                          Compute the L*L' or U'*U factorization of the
                                //                          matrix and solve the system.
                                //
                                Rlacpy("Full", kd + 1, n, a, ldab, afac, ldab);
                                Rlacpy("Full", n, nrhs, b, lda, x, lda);
                                //
                                srnamt = "Rpbsv ";
                                Rpbsv(uplo, n, kd, nrhs, afac, ldab, x, lda, info);
                                //
                                //                          Check error code from Rpbsv .
                                //
                                if (info != izero) {
                                    Alaerh(path, "Rpbsv ", info, izero, uplo, n, n, kd, kd, nrhs, imat, nfail, nerrs, nout);
                                    goto statement_40;
                                } else if (info != 0) {
                                    goto statement_40;
                                }
                                //
                                //                          Reconstruct matrix from factors and compute
                                //                          residual.
                                //
                                Rpbt01(uplo, n, kd, a, ldab, afac, ldab, rwork, result[1 - 1]);
                                //
                                //                          Compute residual of the computed solution.
                                //
                                Rlacpy("Full", n, nrhs, b, lda, work, lda);
                                Rpbt02(uplo, n, kd, nrhs, a, ldab, x, lda, work, lda, rwork, result[2 - 1]);
                                //
                                //                          Check solution from generated exact solution.
                                //
                                Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                                nt = 3;
                                //
                                //                          Print information about the tests that did
                                //                          not pass the threshold.
                                //
                                for (k = 1; k <= nt; k = k + 1) {
                                    if (result[k - 1] >= thresh) {
                                        if (nfail == 0 && nerrs == 0) {
                                            Aladhd(nout, path);
                                        }
                                        write(nout, "(1x,a,', UPLO=''',a1,''', N =',i5,', KD =',i5,"
                                                    "', type ',i1,', test(',i1,')=',g12.5)"),
                                            "Rpbsv ", uplo, n, kd, imat, k, result(k);
                                        nfail++;
                                    }
                                }
                                nrun += nt;
                            statement_40:;
                            }
                            //
                            //                       --- Test Rpbsvx ---
                            //
                            if (!prefac) {
                                Rlaset("Full", kd + 1, n, zero, zero, afac, ldab);
                            }
                            Rlaset("Full", n, nrhs, zero, zero, x, lda);
                            if (iequed > 1 && n > 0) {
                                //
                                //                          Equilibrate the matrix if FACT='F' and
                                //                          EQUED='Y'
                                //
                                Rlaqsb(uplo, n, kd, a, ldab, s, scond, amax, equed);
                            }
                            //
                            //                       Solve the system and compute the condition
                            //                       number and error bounds using Rpbsvx.
                            //
                            srnamt = "Rpbsvx";
                            Rpbsvx(fact, uplo, n, kd, nrhs, a, ldab, afac, ldab, equed, s, b, lda, x, lda, rcond, rwork, &rwork[(nrhs + 1) - 1], work, iwork, info);
                            //
                            //                       Check the error code from Rpbsvx.
                            //
                            if (info != izero) {
                                Alaerh(path, "Rpbsvx", info, izero, fact + uplo, n, n, kd, kd, nrhs, imat, nfail, nerrs, nout);
                                goto statement_60;
                            }
                            //
                            if (info == 0) {
                                if (!prefac) {
                                    //
                                    //                             Reconstruct matrix from factors and
                                    //                             compute residual.
                                    //
                                    Rpbt01(uplo, n, kd, a, ldab, afac, ldab, &rwork[(2 * nrhs + 1) - 1], result[1 - 1]);
                                    k1 = 1;
                                } else {
                                    k1 = 2;
                                }
                                //
                                //                          Compute residual of the computed solution.
                                //
                                Rlacpy("Full", n, nrhs, bsav, lda, work, lda);
                                Rpbt02(uplo, n, kd, nrhs, asav, ldab, x, lda, work, lda, &rwork[(2 * nrhs + 1) - 1], result[2 - 1]);
                                //
                                //                          Check solution from generated exact solution.
                                //
                                if (nofact || (prefac && Mlsame(equed, "N"))) {
                                    Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                                } else {
                                    Rget04(n, nrhs, x, lda, xact, lda, roldc, result[3 - 1]);
                                }
                                //
                                //                          Check the error bounds from iterative
                                //                          refinement.
                                //
                                Rpbt05(uplo, n, kd, nrhs, asav, ldab, b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], result[4 - 1]);
                            } else {
                                k1 = 6;
                            }
                            //
                            //                       Compare RCOND from Rpbsvx with the computed
                            //                       value in RCONDC.
                            //
                            result[6 - 1] = Rget06[(rcond - 1) + (rcondc - 1) * ldRget06];
                            //
                            //                       Print information about the tests that did not
                            //                       pass the threshold.
                            //
                            for (k = k1; k <= 6; k = k + 1) {
                                if (result[k - 1] >= thresh) {
                                    if (nfail == 0 && nerrs == 0) {
                                        Aladhd(nout, path);
                                    }
                                    if (prefac) {
                                        write(nout, "(1x,a,'( ''',a1,''', ''',a1,''', ',i5,', ',i5,"
                                                    "', ... ), EQUED=''',a1,''', type ',i1,', test(',i1,"
                                                    "')=',g12.5)"),
                                            "Rpbsvx", fact, uplo, n, kd, equed, imat, k, result(k);
                                    } else {
                                        write(nout, "(1x,a,'( ''',a1,''', ''',a1,''', ',i5,', ',i5,"
                                                    "', ... ), type ',i1,', test(',i1,')=',g12.5)"),
                                            "Rpbsvx", fact, uplo, n, kd, imat, k, result(k);
                                    }
                                    nfail++;
                                }
                            }
                            nrun += 7 - k1;
                        statement_60:;
                        }
                    }
                statement_80:;
                }
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasvm(path, nout, nfail, nrun, nerrs);
    //
    //     End of Rdrvpb
    //
}
