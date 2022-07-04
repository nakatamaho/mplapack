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

void Cdrvpb(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nrhs, REAL const thresh, bool const tsterr, INTEGER const /* nmax */, COMPLEX *a, COMPLEX *afac, COMPLEX *asav, COMPLEX *b, COMPLEX *bsav, COMPLEX *x, COMPLEX *xact, REAL *s, COMPLEX *work, REAL *rwork, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    //
    char facts[] = {'F', 'N', 'E'};
    char equeds[] = {'N', 'Y'};
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    char path[4] = {};
    char buf[1024];
    char fact_uplo[3];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    const INTEGER nbw = 4;
    INTEGER kdval[nbw];
    INTEGER nb = 0;
    INTEGER nbmin = 0;
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    char xtype[1];
    INTEGER nkd = 0;
    const INTEGER ntypes = 8;
    INTEGER nimat = 0;
    INTEGER ikd = 0;
    INTEGER kd = 0;
    INTEGER ldab = 0;
    INTEGER iuplo = 0;
    INTEGER koff = 0;
    char uplo[1];
    char packit[1];
    INTEGER imat = 0;
    bool zerot = false;
    char type[1];
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist[1];
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER iw = 0;
    INTEGER ioff = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    const REAL zero = 0.0;
    INTEGER iequed = 0;
    char equed[1];
    INTEGER nfact = 0;
    INTEGER ifact = 0;
    char fact[1];
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
    REAL result[ntests];
    INTEGER nt = 0;
    INTEGER k = 0;
    REAL rcond = 0.0;
    INTEGER k1 = 0;
    //
    //     Initialize constants and the random number seed.
    //
    path[0] = 'C';
    path[1] = 'P';
    path[2] = 'B';
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
    infot = 0;
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
        lda = max(n, (INTEGER)1);
        xtype[0] = 'N';
        //
        //        Set limits on the number of loop iterations.
        //
        nkd = max((INTEGER)1, min(n, (INTEGER)4));
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
                    uplo[0] = 'U';
                    packit[0] = 'Q';
                    koff = max((INTEGER)1, kd + 2 - n);
                } else {
                    uplo[0] = 'L';
                    packit[0] = 'B';
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
                        //                    Set up parameters with Clatb4 and generate a test
                        //                    matrix with Clatms.
                        //
                        Clatb4(path, imat, n, n, type, kl, ku, anorm, mode, cndnum, dist);
                        //
                        Clatms(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, kd, kd, packit, &a[koff - 1], ldab, work, info);
                        //
                        //                    Check error code from Clatms.
                        //
                        if (info != 0) {
                            Alaerh(path, "Clatms", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
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
                            Ccopy(izero - i1, &work[iw - 1], 1, &a[(ioff - izero + i1) - 1], 1);
                            iw += izero - i1;
                            Ccopy(i2 - izero + 1, &work[iw - 1], 1, &a[ioff - 1], max(ldab - 1, (INTEGER)1));
                        } else {
                            ioff = (i1 - 1) * ldab + 1;
                            Ccopy(izero - i1, &work[iw - 1], 1, &a[(ioff + izero - i1) - 1], max(ldab - 1, (INTEGER)1));
                            ioff = (izero - 1) * ldab + 1;
                            iw += izero - i1;
                            Ccopy(i2 - izero + 1, &work[iw - 1], 1, &a[ioff - 1], 1);
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
                        i1 = max(izero - kd, (INTEGER)1);
                        i2 = min(izero + kd, n);
                        //
                        if (iuplo == 1) {
                            ioff = (izero - 1) * ldab + kd + 1;
                            Cswap(izero - i1, &a[(ioff - izero + i1) - 1], 1, &work[iw - 1], 1);
                            iw += izero - i1;
                            Cswap(i2 - izero + 1, &a[ioff - 1], max(ldab - 1, (INTEGER)1), &work[iw - 1], 1);
                        } else {
                            ioff = (i1 - 1) * ldab + 1;
                            Cswap(izero - i1, &a[(ioff + izero - i1) - 1], max(ldab - 1, (INTEGER)1), &work[iw - 1], 1);
                            ioff = (izero - 1) * ldab + 1;
                            iw += izero - i1;
                            Cswap(i2 - izero + 1, &a[ioff - 1], 1, &work[iw - 1], 1);
                        }
                    }
                    //
                    //                 Set the imaginary part of the diagonals.
                    //
                    if (iuplo == 1) {
                        Claipd(n, &a[(kd + 1) - 1], ldab, 0);
                    } else {
                        Claipd(n, &a[1 - 1], ldab, 0);
                    }
                    //
                    //                 Save a copy of the matrix A in ASAV.
                    //
                    Clacpy("Full", kd + 1, n, a, ldab, asav, ldab);
                    //
                    for (iequed = 1; iequed <= 2; iequed = iequed + 1) {
                        equed[0] = equeds[iequed - 1];
                        if (iequed == 1) {
                            nfact = 3;
                        } else {
                            nfact = 1;
                        }
                        //
                        for (ifact = 1; ifact <= nfact; ifact = ifact + 1) {
                            fact[0] = facts[ifact - 1];
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
                                //                          with the value returned by Cpbsvx (FACT =
                                //                          'N' reuses the condition number from the
                                //                          previous iteration with FACT = 'F').
                                //
                                Clacpy("Full", kd + 1, n, asav, ldab, afac, ldab);
                                if (equil || iequed > 1) {
                                    //
                                    //                             Compute row and column scale factors to
                                    //                             equilibrate the matrix A.
                                    //
                                    Cpbequ(uplo, n, kd, afac, ldab, s, scond, amax, info);
                                    if (info == 0 && n > 0) {
                                        if (iequed > 1) {
                                            scond = zero;
                                        }
                                        //
                                        //                                Equilibrate the matrix.
                                        //
                                        Claqhb(uplo, n, kd, afac, ldab, s, scond, amax, equed);
                                    }
                                }
                                //
                                //                          Save the condition number of the
                                //                          non-equilibrated system for use in Cget04.
                                //
                                if (equil) {
                                    roldc = rcondc;
                                }
                                //
                                //                          Compute the 1-norm of A.
                                //
                                anorm = Clanhb("1", uplo, n, kd, afac, ldab, rwork);
                                //
                                //                          Factor the matrix A.
                                //
                                Cpbtrf(uplo, n, kd, afac, ldab, info);
                                //
                                //                          Form the inverse of A.
                                //
                                Claset("Full", n, n, COMPLEX(zero), COMPLEX(one), a, lda);
                                Cpbtrs(uplo, n, kd, n, afac, ldab, a, lda, info);
                                //
                                //                          Compute the 1-norm condition number of A.
                                //
                                ainvnm = Clange("1", n, n, a, lda, rwork);
                                if (anorm <= zero || ainvnm <= zero) {
                                    rcondc = one;
                                } else {
                                    rcondc = (one / anorm) / ainvnm;
                                }
                            }
                            //
                            //                       Restore the matrix A.
                            //
                            Clacpy("Full", kd + 1, n, asav, ldab, a, ldab);
                            //
                            //                       Form an exact solution and set the right hand
                            //                       side.
                            //
                            Clarhs(path, xtype, uplo, " ", n, n, kd, kd, nrhs, a, ldab, xact, lda, b, lda, iseed, info);
                            xtype[0] = 'C';
                            Clacpy("Full", n, nrhs, b, lda, bsav, lda);
                            //
                            if (nofact) {
                                //
                                //                          --- Test Cpbsv  ---
                                //
                                //                          Compute the L*L' or U'*U factorization of the
                                //                          matrix and solve the system.
                                //
                                Clacpy("Full", kd + 1, n, a, ldab, afac, ldab);
                                Clacpy("Full", n, nrhs, b, lda, x, lda);
                                //
                                Cpbsv(uplo, n, kd, nrhs, afac, ldab, x, lda, info);
                                //
                                //                          Check error code from Cpbsv .
                                //
                                if (info != izero) {
                                    Alaerh(path, "Cpbsv ", info, izero, uplo, n, n, kd, kd, nrhs, imat, nfail, nerrs, nout);
                                    goto statement_40;
                                } else if (info != 0) {
                                    goto statement_40;
                                }
                                //
                                //                          Reconstruct matrix from factors and compute
                                //                          residual.
                                //
                                Cpbt01(uplo, n, kd, a, ldab, afac, ldab, rwork, result[1 - 1]);
                                //
                                //                          Compute residual of the computed solution.
                                //
                                Clacpy("Full", n, nrhs, b, lda, work, lda);
                                Cpbt02(uplo, n, kd, nrhs, a, ldab, x, lda, work, lda, rwork, result[2 - 1]);
                                //
                                //                          Check solution from generated exact solution.
                                //
                                Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
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
                                        sprintnum_short(buf, result[k - 1]);
                                        write(nout, "(1x,a,', UPLO=''',a1,''', N =',i5,', KD =',i5,"
                                                    "', type ',i1,', test(',i1,')=',a)"),
                                            "Cpbsv ", uplo, n, kd, imat, k, buf;
                                        nfail++;
                                    }
                                }
                                nrun += nt;
                            statement_40:;
                            }
                            //
                            //                       --- Test Cpbsvx ---
                            //
                            if (!prefac) {
                                Claset("Full", kd + 1, n, COMPLEX(zero), COMPLEX(zero), afac, ldab);
                            }
                            Claset("Full", n, nrhs, COMPLEX(zero), COMPLEX(zero), x, lda);
                            if (iequed > 1 && n > 0) {
                                //
                                //                          Equilibrate the matrix if FACT='F' and
                                //                          EQUED='Y'
                                //
                                Claqhb(uplo, n, kd, a, ldab, s, scond, amax, equed);
                            }
                            //
                            //                       Solve the system and compute the condition
                            //                       number and error bounds using Cpbsvx.
                            //
                            Cpbsvx(fact, uplo, n, kd, nrhs, a, ldab, afac, ldab, equed, s, b, lda, x, lda, rcond, rwork, &rwork[(nrhs + 1) - 1], work, &rwork[(2 * nrhs + 1) - 1], info);
                            //
                            //                       Check the error code from Cpbsvx.
                            //
                            if (info != izero) {
                                fact_uplo[0] = fact[0];
                                fact_uplo[1] = uplo[0];
                                fact_uplo[2] = '\0';
                                Alaerh(path, "Cpbsvx", info, izero, fact_uplo, n, n, kd, kd, nrhs, imat, nfail, nerrs, nout);
                                goto statement_60;
                            }
                            //
                            if (info == 0) {
                                if (!prefac) {
                                    //
                                    //                             Reconstruct matrix from factors and

                                    //
                                    Cpbt01(uplo, n, kd, a, ldab, afac, ldab, &rwork[(2 * nrhs + 1) - 1], result[1 - 1]);
                                    k1 = 1;
                                } else {
                                    k1 = 2;
                                }
                                //
                                //                          Compute residual of the computed solution.
                                //
                                Clacpy("Full", n, nrhs, bsav, lda, work, lda);
                                Cpbt02(uplo, n, kd, nrhs, asav, ldab, x, lda, work, lda, &rwork[(2 * nrhs + 1) - 1], result[2 - 1]);
                                //
                                //                          Check solution from generated exact solution.
                                //
                                if (nofact || (prefac && Mlsame(equed, "N"))) {
                                    Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                                } else {
                                    Cget04(n, nrhs, x, lda, xact, lda, roldc, result[3 - 1]);
                                }
                                //
                                //                          Check the error bounds from iterative
                                //                          refinement.
                                //
                                Cpbt05(uplo, n, kd, nrhs, asav, ldab, b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], &result[4 - 1]);
                            } else {
                                k1 = 6;
                            }
                            //
                            //                       Compare RCOND from Cpbsvx with the computed
                            //                       value in RCONDC.
                            //
                            result[6 - 1] = Rget06(rcond, rcondc);
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
                                        sprintnum_short(buf, result[k - 1]);
                                        write(nout, "(1x,a,'( ''',a1,''', ''',a1,''', ',i5,', ',i5,"
                                                    "', ... ), EQUED=''',a1,''', type ',i1,', test(',i1,"
                                                    "')=',a)"),
                                            "Cpbsvx", fact, uplo, n, kd, equed, imat, k, buf;
                                    } else {
                                        sprintnum_short(buf, result[k - 1]);
                                        write(nout, "(1x,a,'( ''',a1,''', ''',a1,''', ',i5,', ',i5,"
                                                    "', ... ), type ',i1,', test(',i1,')=',a)"),
                                            "Cpbsvx", fact, uplo, n, kd, imat, k, buf;
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
    //     End of Cdrvpb
    //
}
