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

void Cdrvgb(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nrhs, REAL const thresh, bool const tsterr, COMPLEX *a, INTEGER const la, COMPLEX *afb, INTEGER const lafb, COMPLEX *asav, COMPLEX *b, COMPLEX *bsav, COMPLEX *x, COMPLEX *xact, REAL *s, COMPLEX *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
    common_write write(cmn);
    //
    str_arr_ref<1> equeds(sve.equeds, [4]);
    str_arr_ref<1> facts(sve.facts, [3]);
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    const INTEGER ntran = 3;
    str_arr_ref<1> transs(sve.transs, [ntran]);
    if (is_called_first_time) {
        {
            static const INTEGER values[] = {1988, 1989, 1990, 1991};
            data_of_type<int>(FEM_VALUES_AND_SIZE), iseedy;
        }
        {
            static const char *values[] = {"N", "T", "C"};
            data_of_type_str(FEM_VALUES_AND_SIZE), transs;
        }
        {
            static const char *values[] = {"F", "N", "E"};
            data_of_type_str(FEM_VALUES_AND_SIZE), facts;
        }
        {
            static const char *values[] = {"N", "R", "C", "B"};
            data_of_type_str(FEM_VALUES_AND_SIZE), equeds;
        }
    }
    char path[3];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    INTEGER nb = 0;
    INTEGER nbmin = 0;
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER ldb = 0;
    char xtype;
    INTEGER nkl = 0;
    INTEGER nku = 0;
    const INTEGER ntypes = 8;
    INTEGER nimat = 0;
    INTEGER ikl = 0;
    INTEGER kl = 0;
    INTEGER iku = 0;
    INTEGER ku = 0;
    INTEGER lda = 0;
    INTEGER ldafb = 0;
    INTEGER imat = 0;
    bool zerot = false;
    char type;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist;
    const REAL one = 1.0;
    REAL rcondc = 0.0;
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER ioff = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER iequed = 0;
    char equed;
    INTEGER nfact = 0;
    INTEGER ifact = 0;
    char fact;
    bool prefac = false;
    bool nofact = false;
    bool equil = false;
    REAL rcondo = 0.0;
    REAL rcondi = 0.0;
    REAL rowcnd = 0.0;
    REAL colcnd = 0.0;
    REAL amax = 0.0;
    REAL roldo = 0.0;
    REAL roldi = 0.0;
    REAL anormo = 0.0;
    REAL anormi = 0.0;
    REAL ainvnm = 0.0;
    INTEGER itran = 0;
    char trans;
    const INTEGER ntests = 7;
    REAL result[ntests];
    INTEGER nt = 0;
    INTEGER k = 0;
    REAL rcond = 0.0;
    REAL anrmpv = 0.0;
    REAL rdum[1];
    REAL rpvgrw = 0.0;
    INTEGER k1 = 0;
    bool trfcon = false;
    REAL roldc = 0.0;
    INTEGER n_err_bnds = 0;
    REAL rpvgrw_svxx = 0.0;
    REAL berr[nrhs];
    REAL errbnds_n[nrhs * 3];
    REAL errbnds_c[nrhs * 3];
    static const char *format_9995 = "(1x,a,'( ''',a1,''',''',a1,''',',i5,',',i5,',',i5,',...), EQUED=''',a1,"
                                     "''', type ',i1,', test(',i1,')=',g12.5)";
    static const char *format_9996 = "(1x,a,'( ''',a1,''',''',a1,''',',i5,',',i5,',',i5,',...), type ',i1,"
                                     "', test(',i1,')=',g12.5)";
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
    path[(1 - 1)] = "Zomplex precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "GB";
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
        ldb = max(n, 1);
        xtype = "N";
        //
        //        Set limits on the number of loop iterations.
        //
        nkl = max({(INTEGER)1, min(n, 4)});
        if (n == 0) {
            nkl = 1;
        }
        nku = nkl;
        nimat = ntypes;
        if (n <= 0) {
            nimat = 1;
        }
        //
        for (ikl = 1; ikl <= nkl; ikl = ikl + 1) {
            //
            //           Do for KL = 0, N-1, (3N-1)/4, and (N+1)/4. This order makes
            //           it easier to skip redundant values for small values of N.
            //
            if (ikl == 1) {
                kl = 0;
            } else if (ikl == 2) {
                kl = max(n - 1, 0);
            } else if (ikl == 3) {
                kl = (3 * n - 1) / 4;
            } else if (ikl == 4) {
                kl = (n + 1) / 4;
            }
            for (iku = 1; iku <= nku; iku = iku + 1) {
                //
                //              Do for KU = 0, N-1, (3N-1)/4, and (N+1)/4. This order
                //              makes it easier to skip redundant values for small
                //              values of N.
                //
                if (iku == 1) {
                    ku = 0;
                } else if (iku == 2) {
                    ku = max(n - 1, 0);
                } else if (iku == 3) {
                    ku = (3 * n - 1) / 4;
                } else if (iku == 4) {
                    ku = (n + 1) / 4;
                }
                //
                //              Check that A and AFB are big enough to generate this
                //              matrix.
                //
                lda = kl + ku + 1;
                ldafb = 2 * kl + ku + 1;
                if (lda * n > la || ldafb * n > lafb) {
                    if (nfail == 0 && nerrs == 0) {
                        Aladhd(nout, path);
                    }
                    if (lda * n > la) {
                        write(nout, "(' *** In Cdrvgb, LA=',i5,' is too small for N=',i5,', KU=',i5,"
                                    "', KL=',i5,/,' ==> Increase LA to at least ',i5)"),
                            la, n, kl, ku, n *(kl + ku + 1);
                        nerrs++;
                    }
                    if (ldafb * n > lafb) {
                        write(nout, "(' *** In Cdrvgb, LAFB=',i5,' is too small for N=',i5,', KU=',"
                                    "i5,', KL=',i5,/,' ==> Increase LAFB to at least ',i5)"),
                            lafb, n, kl, ku, n *(2 * kl + ku + 1);
                        nerrs++;
                    }
                    goto statement_130;
                }
                //
                for (imat = 1; imat <= nimat; imat = imat + 1) {
                    //
                    //                 Do the tests only if DOTYPE( IMAT ) is true.
                    //
                    if (!dotype[imat - 1]) {
                        goto statement_120;
                    }
                    //
                    //                 Skip types 2, 3, or 4 if the matrix is too small.
                    //
                    zerot = imat >= 2 && imat <= 4;
                    if (zerot && n < imat - 1) {
                        goto statement_120;
                    }
                    //
                    //                 Set up parameters with Clatb4 and generate a
                    //                 test matrix with Clatms.
                    //
                    Clatb4(path, imat, n, n, type, kl, ku, anorm, mode, cndnum, dist);
                    rcondc = one / cndnum;
                    //
                    Clatms(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, kl, ku, "Z", a, lda, work, info);
                    //
                    //                 Check the error code from Clatms.
                    //
                    if (info != 0) {
                        Alaerh(path, "Clatms", info, 0, " ", n, n, kl, ku, -1, imat, nfail, nerrs, nout);
                        goto statement_120;
                    }
                    //
                    //                 For types 2, 3, and 4, zero one or more columns of
                    //                 the matrix to test that INFO is returned correctly.
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
                        ioff = (izero - 1) * lda;
                        if (imat < 4) {
                            i1 = max((INTEGER)1, ku + 2 - izero);
                            i2 = min(kl + ku + 1, ku + 1 + (n - izero));
                            for (i = i1; i <= i2; i = i + 1) {
                                a[(ioff + i) - 1] = zero;
                            }
                        } else {
                            for (j = izero; j <= n; j = j + 1) {
                                for (i = max((INTEGER)1, ku + 2 - j); i <= min(kl + ku + 1, ku + 1 + (n - j)); i = i + 1) {
                                    a[(ioff + i) - 1] = zero;
                                }
                                ioff += lda;
                            }
                        }
                    }
                    //
                    //                 Save a copy of the matrix A in ASAV.
                    //
                    Clacpy("Full", kl + ku + 1, n, a, lda, asav, lda);
                    //
                    for (iequed = 1; iequed <= 4; iequed = iequed + 1) {
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
                                rcondo = zero;
                                rcondi = zero;
                                //
                            } else if (!nofact) {
                                //
                                //                          Compute the condition number for comparison
                                //                          with the value returned by Rgesvx (FACT =
                                //                          'N' reuses the condition number from the
                                //                          previous iteration with FACT = 'F').
                                //
                                Clacpy("Full", kl + ku + 1, n, asav, lda, afb[(kl + 1) - 1], ldafb);
                                if (equil || iequed > 1) {
                                    //
                                    //                             Compute row and column scale factors to
                                    //                             equilibrate the matrix A.
                                    //
                                    Cgbequ(n, n, kl, ku, afb[(kl + 1) - 1], ldafb, s, s[(n + 1) - 1], rowcnd, colcnd, amax, info);
                                    if (info == 0 && n > 0) {
                                        if (Mlsame(equed, "R")) {
                                            rowcnd = zero;
                                            colcnd = one;
                                        } else if (Mlsame(equed, "C")) {
                                            rowcnd = one;
                                            colcnd = zero;
                                        } else if (Mlsame(equed, "B")) {
                                            rowcnd = zero;
                                            colcnd = zero;
                                        }
                                        //
                                        //                                Equilibrate the matrix.
                                        //
                                        Claqgb(n, n, kl, ku, afb[(kl + 1) - 1], ldafb, s, s[(n + 1) - 1], rowcnd, colcnd, amax, equed);
                                    }
                                }
                                //
                                //                          Save the condition number of the
                                //                          non-equilibrated system for use in Cget04.
                                //
                                if (equil) {
                                    roldo = rcondo;
                                    roldi = rcondi;
                                }
                                //
                                //                          Compute the 1-norm and infinity-norm of A.
                                //
                                anormo = Clangb("1", n, kl, ku, afb[(kl + 1) - 1], ldafb, rwork);
                                anormi = Clangb("I", n, kl, ku, afb[(kl + 1) - 1], ldafb, rwork);
                                //
                                //                          Factor the matrix A.
                                //
                                Cgbtrf(n, n, kl, ku, afb, ldafb, iwork, info);
                                //
                                //                          Form the inverse of A.
                                //
                                Claset("Full", n, n, COMPLEX(zero), COMPLEX(one), work, ldb);
                                Cgbtrs("No transpose", n, kl, ku, n, afb, ldafb, iwork, work, ldb, info);
                                //
                                //                          Compute the 1-norm condition number of A.
                                //
                                ainvnm = Clange("1", n, n, work, ldb, rwork);
                                if (anormo <= zero || ainvnm <= zero) {
                                    rcondo = one;
                                } else {
                                    rcondo = (one / anormo) / ainvnm;
                                }
                                //
                                //                          Compute the infinity-norm condition number
                                //                          of A.
                                //
                                ainvnm = Clange("I", n, n, work, ldb, rwork);
                                if (anormi <= zero || ainvnm <= zero) {
                                    rcondi = one;
                                } else {
                                    rcondi = (one / anormi) / ainvnm;
                                }
                            }
                            //
                            for (itran = 1; itran <= ntran; itran = itran + 1) {
                                //
                                //                          Do for each value of TRANS.
                                //
                                trans = transs[itran - 1];
                                if (itran == 1) {
                                    rcondc = rcondo;
                                } else {
                                    rcondc = rcondi;
                                }
                                //
                                //                          Restore the matrix A.
                                //
                                Clacpy("Full", kl + ku + 1, n, asav, lda, a, lda);
                                //
                                //                          Form an exact solution and set the right hand
                                //                          side.
                                //
                                Clarhs(path, xtype, "Full", trans, n, n, kl, ku, nrhs, a, lda, xact, ldb, b, ldb, iseed, info);
                                xtype = "C";
                                Clacpy("Full", n, nrhs, b, ldb, bsav, ldb);
                                //
                                if (nofact && itran == 1) {
                                    //
                                    //                             --- Test Cgbsv  ---
                                    //
                                    //                             Compute the LU factorization of the matrix
                                    //                             and solve the system.
                                    //
                                    Clacpy("Full", kl + ku + 1, n, a, lda, afb[(kl + 1) - 1], ldafb);
                                    Clacpy("Full", n, nrhs, b, ldb, x, ldb);
                                    //
                                    Cgbsv(n, kl, ku, nrhs, afb, ldafb, iwork, x, ldb, info);
                                    //
                                    //                             Check error code from Cgbsv .
                                    //
                                    if (info != izero) {
                                        Alaerh(path, "Cgbsv ", info, izero, " ", n, n, kl, ku, nrhs, imat, nfail, nerrs, nout);
                                    }
                                    //
                                    //                             Reconstruct matrix from factors and
                                    //                             compute residual.
                                    //
                                    Cgbt01(n, n, kl, ku, a, lda, afb, ldafb, iwork, work, result[1 - 1]);
                                    nt = 1;
                                    if (izero == 0) {
                                        //
                                        //                                Compute residual of the computed
                                        //                                solution.
                                        //
                                        Clacpy("Full", n, nrhs, b, ldb, work, ldb);
                                        Cgbt02("No transpose", n, n, kl, ku, nrhs, a, lda, x, ldb, work, ldb, result[2 - 1]);
                                        //
                                        //                                Check solution from generated exact
                                        //                                solution.
                                        //
                                        Cget04(n, nrhs, x, ldb, xact, ldb, rcondc, result[3 - 1]);
                                        nt = 3;
                                    }
                                    //
                                    //                             Print information about the tests that did
                                    //                             not pass the threshold.
                                    //
                                    for (k = 1; k <= nt; k = k + 1) {
                                        if (result[k - 1] >= thresh) {
                                            if (nfail == 0 && nerrs == 0) {
                                                Aladhd(nout, path);
                                            }
                                            write(nout, "(1x,a,', N=',i5,', KL=',i5,', KU=',i5,', type ',i1,"
                                                        "', test(',i1,')=',g12.5)"),
                                                "Cgbsv ", n, kl, ku, imat, k, result(k);
                                            nfail++;
                                        }
                                    }
                                    nrun += nt;
                                }
                                //
                                //                          --- Test Cgbsvx ---
                                //
                                if (!prefac) {
                                    Claset("Full", 2 * kl + ku + 1, n, COMPLEX(zero), COMPLEX(zero), afb, ldafb);
                                }
                                Claset("Full", n, nrhs, COMPLEX(zero), COMPLEX(zero), x, ldb);
                                if (iequed > 1 && n > 0) {
                                    //
                                    //                             Equilibrate the matrix if FACT = 'F' and
                                    //                             EQUED = 'R', 'C', or 'B'.
                                    //
                                    Claqgb(n, n, kl, ku, a, lda, s, s[(n + 1) - 1], rowcnd, colcnd, amax, equed);
                                }
                                //
                                //                          Solve the system and compute the condition
                                //                          number and error bounds using Cgbsvx.
                                //
                                Cgbsvx(fact, trans, n, kl, ku, nrhs, a, lda, afb, ldafb, iwork, equed, s, s[(ldb + 1) - 1], b, ldb, x, ldb, rcond, rwork, &rwork[(nrhs + 1) - 1], work, &rwork[(2 * nrhs + 1) - 1], info);
                                //
                                //                          Check the error code from Cgbsvx.
                                //
                                if (info != izero) {
                                    Alaerh(path, "Cgbsvx", info, izero, fact + trans, n, n, kl, ku, nrhs, imat, nfail, nerrs, nout);
                                }
                                //
                                //                          Compare RWORK(2*NRHS+1) from Cgbsvx with the
                                //                          computed reciprocal pivot growth RPVGRW
                                //
                                if (info != 0) {
                                    anrmpv = zero;
                                    for (j = 1; j <= info; j = j + 1) {
                                        for (i = max(ku + 2 - j, 1); i <= min(n + ku + 1 - j, kl + ku + 1); i = i + 1) {
                                            anrmpv = max(anrmpv, abs(a[(i + (j - 1) * lda) - 1]));
                                        }
                                    }
                                    rpvgrw = Clantb("M", "U", "N", info, min(info - 1, kl + ku), afb[(max((kl + ku + 2 - info)) - 1) * ldafb], ldafb, rdum);
                                    if (rpvgrw == zero) {
                                        rpvgrw = one;
                                    } else {
                                        rpvgrw = anrmpv / rpvgrw;
                                    }
                                } else {
                                    rpvgrw = Clantb("M", "U", "N", n, kl + ku, afb, ldafb, rdum);
                                    if (rpvgrw == zero) {
                                        rpvgrw = one;
                                    } else {
                                        rpvgrw = Clangb("M", n, kl, ku, a, lda, rdum) / rpvgrw;
                                    }
                                }
                                result[7 - 1] = abs(rpvgrw - rwork[(2 * nrhs + 1) - 1]) / max(rwork[(2 * nrhs + 1) - 1], rpvgrw) / Rlamch("E");
                                //
                                if (!prefac) {
                                    //
                                    //                             Reconstruct matrix from factors and
                                    //                             compute residual.
                                    //
                                    Cgbt01(n, n, kl, ku, a, lda, afb, ldafb, iwork, work, result[1 - 1]);
                                    k1 = 1;
                                } else {
                                    k1 = 2;
                                }
                                //
                                if (info == 0) {
                                    trfcon = false;
                                    //
                                    //                             Compute residual of the computed solution.
                                    //
                                    Clacpy("Full", n, nrhs, bsav, ldb, work, ldb);
                                    Cgbt02(trans, n, n, kl, ku, nrhs, asav, lda, x, ldb, work, ldb, result[2 - 1]);
                                    //
                                    //                             Check solution from generated exact
                                    //                             solution.
                                    //
                                    if (nofact || (prefac && Mlsame(equed, "N"))) {
                                        Cget04(n, nrhs, x, ldb, xact, ldb, rcondc, result[3 - 1]);
                                    } else {
                                        if (itran == 1) {
                                            roldc = roldo;
                                        } else {
                                            roldc = roldi;
                                        }
                                        Cget04(n, nrhs, x, ldb, xact, ldb, roldc, result[3 - 1]);
                                    }
                                    //
                                    //                             Check the error bounds from iterative
                                    //                             refinement.
                                    //
                                    Cgbt05(trans, n, kl, ku, nrhs, asav, lda, bsav, ldb, x, ldb, xact, ldb, rwork, &rwork[(nrhs + 1) - 1], result[4 - 1]);
                                } else {
                                    trfcon = true;
                                }
                                //
                                //                          Compare RCOND from Cgbsvx with the computed
                                //                          value in RCONDC.
                                //
                                result[6 - 1] = Rget06[(rcond - 1) + (rcondc - 1) * ldRget06];
                                //
                                //                          Print information about the tests that did
                                //                          not pass the threshold.
                                //
                                if (!trfcon) {
                                    for (k = k1; k <= ntests; k = k + 1) {
                                        if (result[k - 1] >= thresh) {
                                            if (nfail == 0 && nerrs == 0) {
                                                Aladhd(nout, path);
                                            }
                                            if (prefac) {
                                                write(nout, format_9995), "Cgbsvx", fact, trans, n, kl, ku, equed, imat, k, result(k);
                                            } else {
                                                write(nout, format_9996), "Cgbsvx", fact, trans, n, kl, ku, imat, k, result(k);
                                            }
                                            nfail++;
                                        }
                                    }
                                    nrun += 7 - k1;
                                } else {
                                    if (result[1 - 1] >= thresh && !prefac) {
                                        if (nfail == 0 && nerrs == 0) {
                                            Aladhd(nout, path);
                                        }
                                        if (prefac) {
                                            write(nout, format_9995), "Cgbsvx", fact, trans, n, kl, ku, equed, imat, 1, result(1);
                                        } else {
                                            write(nout, format_9996), "Cgbsvx", fact, trans, n, kl, ku, imat, 1, result(1);
                                        }
                                        nfail++;
                                        nrun++;
                                    }
                                    if (result[6 - 1] >= thresh) {
                                        if (nfail == 0 && nerrs == 0) {
                                            Aladhd(nout, path);
                                        }
                                        if (prefac) {
                                            write(nout, format_9995), "Cgbsvx", fact, trans, n, kl, ku, equed, imat, 6, result(6);
                                        } else {
                                            write(nout, format_9996), "Cgbsvx", fact, trans, n, kl, ku, imat, 6, result(6);
                                        }
                                        nfail++;
                                        nrun++;
                                    }
                                    if (result[7 - 1] >= thresh) {
                                        if (nfail == 0 && nerrs == 0) {
                                            Aladhd(nout, path);
                                        }
                                        if (prefac) {
                                            write(nout, format_9995), "Cgbsvx", fact, trans, n, kl, ku, equed, imat, 7, result(7);
                                        } else {
                                            write(nout, format_9996), "Cgbsvx", fact, trans, n, kl, ku, imat, 7, result(7);
                                        }
                                        nfail++;
                                        nrun++;
                                    }
                                }
                                //
                                //                    --- Test Cgbsvxx ---
                                //
                                //                    Restore the matrices A and B.
                                //
                                //                     write(*,*) 'begin Cgbsvxx testing'
                                //
                                Clacpy("Full", kl + ku + 1, n, asav, lda, a, lda);
                                Clacpy("Full", n, nrhs, bsav, ldb, b, ldb);
                                //
                                if (!prefac) {
                                    Claset("Full", 2 * kl + ku + 1, n, COMPLEX(zero), COMPLEX(zero), afb, ldafb);
                                }
                                Claset("Full", n, nrhs, COMPLEX(zero), COMPLEX(zero), x, ldb);
                                if (iequed > 1 && n > 0) {
                                    //
                                    //                       Equilibrate the matrix if FACT = 'F' and
                                    //                       EQUED = 'R', 'C', or 'B'.
                                    //
                                    Claqgb(n, n, kl, ku, a, lda, s, s[(n + 1) - 1], rowcnd, colcnd, amax, equed);
                                }
                                //
                                //                    Solve the system and compute the condition number
                                //                    and error bounds using Cgbsvxx.
                                //
                                n_err_bnds = 3;
                                Cgbsvxx(fact, trans, n, kl, ku, nrhs, a, lda, afb, ldafb, iwork, equed, s, s[(n + 1) - 1], b, ldb, x, ldb, rcond, rpvgrw_svxx, berr, n_err_bnds, errbnds_n, errbnds_c, 0, zero, work, rwork, info);
                                //
                                //                    Check the error code from Cgbsvxx.
                                //
                                if (info == n + 1) {
                                    goto statement_90;
                                }
                                if (info != izero) {
                                    Alaerh(path, "Cgbsvxx", info, izero, fact + trans, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                                    goto statement_90;
                                }
                                //
                                //                    Compare rpvgrw_svxx from Cgesvxx with the computed
                                //                    reciprocal pivot growth factor RPVGRW
                                //
                                if (info > 0 && info < n + 1) {
                                    rpvgrw = Cla_gbrpvgrw(n, kl, ku, info, a, lda, afb, ldafb);
                                } else {
                                    rpvgrw = Cla_gbrpvgrw(n, kl, ku, n, a, lda, afb, ldafb);
                                }
                                //
                                result[7 - 1] = abs(rpvgrw - rpvgrw_svxx) / max(rpvgrw_svxx, rpvgrw) / Rlamch("E");
                                //
                                if (!prefac) {
                                    //
                                    //                       Reconstruct matrix from factors and compute
                                    //                       residual.
                                    //
                                    Cgbt01(n, n, kl, ku, a, lda, afb, ldafb, iwork, &work[(2 * nrhs + 1) - 1], result[1 - 1]);
                                    k1 = 1;
                                } else {
                                    k1 = 2;
                                }
                                //
                                if (info == 0) {
                                    trfcon = false;
                                    //
                                    //                       Compute residual of the computed solution.
                                    //
                                    Clacpy("Full", n, nrhs, bsav, ldb, work, ldb);
                                    Cgbt02(trans, n, n, kl, ku, nrhs, asav, lda, x, ldb, work, ldb, result[2 - 1]);
                                    //
                                    //                       Check solution from generated exact solution.
                                    //
                                    if (nofact || (prefac && Mlsame(equed, "N"))) {
                                        Cget04(n, nrhs, x, ldb, xact, ldb, rcondc, result[3 - 1]);
                                    } else {
                                        if (itran == 1) {
                                            roldc = roldo;
                                        } else {
                                            roldc = roldi;
                                        }
                                        Cget04(n, nrhs, x, ldb, xact, ldb, roldc, result[3 - 1]);
                                    }
                                } else {
                                    trfcon = true;
                                }
                                //
                                //                    Compare RCOND from Cgbsvxx with the computed value
                                //                    in RCONDC.
                                //
                                result[6 - 1] = Rget06[(rcond - 1) + (rcondc - 1) * ldRget06];
                                //
                                //                    Print information about the tests that did not pass
                                //                    the threshold.
                                //
                                if (!trfcon) {
                                    for (k = k1; k <= ntests; k = k + 1) {
                                        if (result[k - 1] >= thresh) {
                                            if (nfail == 0 && nerrs == 0) {
                                                Aladhd(nout, path);
                                            }
                                            if (prefac) {
                                                write(nout, format_9995), "Cgbsvxx", fact, trans, n, kl, ku, equed, imat, k, result(k);
                                            } else {
                                                write(nout, format_9996), "Cgbsvxx", fact, trans, n, kl, ku, imat, k, result(k);
                                            }
                                            nfail++;
                                        }
                                    }
                                    nrun += 7 - k1;
                                } else {
                                    if (result[1 - 1] >= thresh && !prefac) {
                                        if (nfail == 0 && nerrs == 0) {
                                            Aladhd(nout, path);
                                        }
                                        if (prefac) {
                                            write(nout, format_9995), "Cgbsvxx", fact, trans, n, kl, ku, equed, imat, 1, result(1);
                                        } else {
                                            write(nout, format_9996), "Cgbsvxx", fact, trans, n, kl, ku, imat, 1, result(1);
                                        }
                                        nfail++;
                                        nrun++;
                                    }
                                    if (result[6 - 1] >= thresh) {
                                        if (nfail == 0 && nerrs == 0) {
                                            Aladhd(nout, path);
                                        }
                                        if (prefac) {
                                            write(nout, format_9995), "Cgbsvxx", fact, trans, n, kl, ku, equed, imat, 6, result(6);
                                        } else {
                                            write(nout, format_9996), "Cgbsvxx", fact, trans, n, kl, ku, imat, 6, result(6);
                                        }
                                        nfail++;
                                        nrun++;
                                    }
                                    if (result[7 - 1] >= thresh) {
                                        if (nfail == 0 && nerrs == 0) {
                                            Aladhd(nout, path);
                                        }
                                        if (prefac) {
                                            write(nout, format_9995), "Cgbsvxx", fact, trans, n, kl, ku, equed, imat, 7, result(7);
                                        } else {
                                            write(nout, format_9996), "Cgbsvxx", fact, trans, n, kl, ku, imat, 7, result(7);
                                        }
                                        nfail++;
                                        nrun++;
                                    }
                                    //
                                }
                            //
                            statement_90:;
                            }
                        statement_100:;
                        }
                    }
                statement_120:;
                }
            statement_130:;
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasvm(path, nout, nfail, nrun, nerrs);
    //
    //     Test Error Bounds from Cgbsvxx
    //
    Cebchvxx(thresh, path);
    //
    //     End of Cdrvgb
    //
}
