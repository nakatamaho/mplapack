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

void Cdrvpp(common &cmn, bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nrhs, REAL const thresh, bool const tsterr, INTEGER const /* nmax */, COMPLEX *a, COMPLEX *afac, COMPLEX *asav, COMPLEX *b, COMPLEX *bsav, COMPLEX *x, COMPLEX *xact, REAL *s, COMPLEX *work, REAL *rwork, INTEGER const nout) {
    FEM_CMN_SVE(Cdrvpp);
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
            static const char *values[] = {"F", "N", "E"};
            data_of_type_str(FEM_VALUES_AND_SIZE), facts;
        }
        {
            static const char *values[] = {"C", "R"};
            data_of_type_str(FEM_VALUES_AND_SIZE), packs;
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
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    INTEGER npp = 0;
    char xtype = char0;
    const INTEGER ntypes = 9;
    INTEGER nimat = 0;
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
    const REAL one = 1.0;
    REAL rcondc = 0.0;
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER ioff = 0;
    const REAL zero = 0.0;
    INTEGER iequed = 0;
    char equed = char0;
    INTEGER nfact = 0;
    INTEGER ifact = 0;
    char fact = char0;
    bool prefac = false;
    bool nofact = false;
    bool equil = false;
    REAL scond = 0.0;
    REAL amax = 0.0;
    REAL roldc = 0.0;
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
        Cerrvx(path, nout);
    }
    cmn.infot = 0;
    //
    //     Do for each value of N in NVAL
    //
    for (in = 1; in <= nn; in = in + 1) {
        n = nval[in - 1];
        lda = max(n, 1);
        npp = n * (n + 1) / 2;
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
                goto statement_130;
            }
            //
            //           Skip types 3, 4, or 5 if the matrix size is too small.
            //
            zerot = imat >= 3 && imat <= 5;
            if (zerot && n < imat - 2) {
                goto statement_130;
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
                rcondc = one / cndnum;
                //
                srnamt = "ZLATMS";
                zlatms(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, kl, ku, packit, a, lda, work, info);
                //
                //              Check error code from ZLATMS.
                //
                if (info != 0) {
                    Alaerh(path, "ZLATMS", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                    goto statement_120;
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
                //              Save a copy of the matrix A in ASAV.
                //
                Ccopy(npp, a, 1, asav, 1);
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
                                goto statement_100;
                            }
                            rcondc = zero;
                            //
                        } else if (!Mlsame(fact, "N")) {
                            //
                            //                       Compute the condition number for comparison with
                            //                       the value returned by Cppsvx (FACT = 'N' reuses
                            //                       the condition number from the previous iteration
                            //                          with FACT = 'F').
                            //
                            Ccopy(npp, asav, 1, afac, 1);
                            if (equil || iequed > 1) {
                                //
                                //                          Compute row and column scale factors to
                                //                          equilibrate the matrix A.
                                //
                                Cppequ(uplo, n, afac, s, scond, amax, info);
                                if (info == 0 && n > 0) {
                                    if (iequed > 1) {
                                        scond = zero;
                                    }
                                    //
                                    //                             Equilibrate the matrix.
                                    //
                                    Claqhp(uplo, n, afac, s, scond, amax, equed);
                                }
                            }
                            //
                            //                       Save the condition number of the
                            //                       non-equilibrated system for use in Cget04.
                            //
                            if (equil) {
                                roldc = rcondc;
                            }
                            //
                            //                       Compute the 1-norm of A.
                            //
                            anorm = Clanhp("1", uplo, n, afac, rwork);
                            //
                            //                       Factor the matrix A.
                            //
                            Cpptrf(uplo, n, afac, info);
                            //
                            //                       Form the inverse of A.
                            //
                            Ccopy(npp, afac, 1, a, 1);
                            Cpptri(uplo, n, a, info);
                            //
                            //                       Compute the 1-norm condition number of A.
                            //
                            ainvnm = Clanhp("1", uplo, n, a, rwork);
                            if (anorm <= zero || ainvnm <= zero) {
                                rcondc = one;
                            } else {
                                rcondc = (one / anorm) / ainvnm;
                            }
                        }
                        //
                        //                    Restore the matrix A.
                        //
                        Ccopy(npp, asav, 1, a, 1);
                        //
                        //                    Form an exact solution and set the right hand side.
                        //
                        srnamt = "Clarhs";
                        Clarhs(path, xtype, uplo, " ", n, n, kl, ku, nrhs, a, lda, xact, lda, b, lda, iseed, info);
                        xtype = "C";
                        Clacpy("Full", n, nrhs, b, lda, bsav, lda);
                        //
                        if (nofact) {
                            //
                            //                       --- Test Cppsv  ---
                            //
                            //                       Compute the L*L' or U'*U factorization of the
                            //                       matrix and solve the system.
                            //
                            Ccopy(npp, a, 1, afac, 1);
                            Clacpy("Full", n, nrhs, b, lda, x, lda);
                            //
                            srnamt = "Cppsv ";
                            Cppsv(uplo, n, nrhs, afac, x, lda, info);
                            //
                            //                       Check error code from Cppsv .
                            //
                            if (info != izero) {
                                Alaerh(path, "Cppsv ", info, izero, uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                                goto statement_70;
                            } else if (info != 0) {
                                goto statement_70;
                            }
                            //
                            //                       Reconstruct matrix from factors and compute
                            //                       residual.
                            //
                            Cppt01(uplo, n, a, afac, rwork, result[1 - 1]);
                            //
                            //                       Compute residual of the computed solution.
                            //
                            Clacpy("Full", n, nrhs, b, lda, work, lda);
                            Cppt02(uplo, n, nrhs, a, x, lda, work, lda, rwork, result[2 - 1]);
                            //
                            //                       Check solution from generated exact solution.
                            //
                            Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                            nt = 3;
                            //
                            //                       Print information about the tests that did not
                            //                       pass the threshold.
                            //
                            for (k = 1; k <= nt; k = k + 1) {
                                if (result[k - 1] >= thresh) {
                                    if (nfail == 0 && nerrs == 0) {
                                        Aladhd(nout, path);
                                    }
                                    write(nout, "(1x,a,', UPLO=''',a1,''', N =',i5,', type ',i1,', test(',"
                                                "i1,')=',g12.5)"),
                                        "Cppsv ", uplo, n, imat, k, result(k);
                                    nfail++;
                                }
                            }
                            nrun += nt;
                        statement_70:;
                        }
                        //
                        //                    --- Test Cppsvx ---
                        //
                        if (!prefac && npp > 0) {
                            Claset("Full", npp, 1, COMPLEX(zero), COMPLEX(zero), afac, npp);
                        }
                        Claset("Full", n, nrhs, COMPLEX(zero), COMPLEX(zero), x, lda);
                        if (iequed > 1 && n > 0) {
                            //
                            //                       Equilibrate the matrix if FACT='F' and
                            //                       EQUED='Y'.
                            //
                            Claqhp(uplo, n, a, s, scond, amax, equed);
                        }
                        //
                        //                    Solve the system and compute the condition number
                        //                    and error bounds using Cppsvx.
                        //
                        srnamt = "Cppsvx";
                        Cppsvx(fact, uplo, n, nrhs, a, afac, equed, s, b, lda, x, lda, rcond, rwork, &rwork[(nrhs + 1) - 1], work, &rwork[(2 * nrhs + 1) - 1], info);
                        //
                        //                    Check the error code from Cppsvx.
                        //
                        if (info != izero) {
                            Alaerh(path, "Cppsvx", info, izero, fact + uplo, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                            goto statement_90;
                        }
                        //
                        if (info == 0) {
                            if (!prefac) {
                                //
                                //                          Reconstruct matrix from factors and compute
                                //                          residual.
                                //
                                Cppt01(uplo, n, a, afac, &rwork[(2 * nrhs + 1) - 1], result[1 - 1]);
                                k1 = 1;
                            } else {
                                k1 = 2;
                            }
                            //
                            //                       Compute residual of the computed solution.
                            //
                            Clacpy("Full", n, nrhs, bsav, lda, work, lda);
                            Cppt02(uplo, n, nrhs, asav, x, lda, work, lda, &rwork[(2 * nrhs + 1) - 1], result[2 - 1]);
                            //
                            //                       Check solution from generated exact solution.
                            //
                            if (nofact || (prefac && Mlsame(equed, "N"))) {
                                Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                            } else {
                                Cget04(n, nrhs, x, lda, xact, lda, roldc, result[3 - 1]);
                            }
                            //
                            //                       Check the error bounds from iterative
                            //                       refinement.
                            //
                            Cppt05(uplo, n, nrhs, asav, b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], result[4 - 1]);
                        } else {
                            k1 = 6;
                        }
                        //
                        //                    Compare RCOND from Cppsvx with the computed value
                        //                    in RCONDC.
                        //
                        result[6 - 1] = Rget06[(rcond - 1) + (rcondc - 1) * ldRget06];
                        //
                        //                    Print information about the tests that did not pass
                        //                    the threshold.
                        //
                        for (k = k1; k <= 6; k = k + 1) {
                            if (result[k - 1] >= thresh) {
                                if (nfail == 0 && nerrs == 0) {
                                    Aladhd(nout, path);
                                }
                                if (prefac) {
                                    write(nout, "(1x,a,', FACT=''',a1,''', UPLO=''',a1,''', N=',i5,"
                                                "', EQUED=''',a1,''', type ',i1,', test(',i1,')=',g12.5)"),
                                        "Cppsvx", fact, uplo, n, equed, imat, k, result(k);
                                } else {
                                    write(nout, "(1x,a,', FACT=''',a1,''', UPLO=''',a1,''', N=',i5,"
                                                "', type ',i1,', test(',i1,')=',g12.5)"),
                                        "Cppsvx", fact, uplo, n, imat, k, result(k);
                                }
                                nfail++;
                            }
                        }
                        nrun += 7 - k1;
                    statement_90:
                    statement_100:;
                    }
                }
            statement_120:;
            }
        statement_130:;
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasvm(path, nout, nfail, nrun, nerrs);
    //
    //     End of Cdrvpp
    //
}
