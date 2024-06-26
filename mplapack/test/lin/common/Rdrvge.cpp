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

void Rdrvge(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nrhs, REAL const thresh, bool const tsterr, INTEGER const nmax, REAL *a, REAL *afac, REAL *asav, REAL *b, REAL *bsav, REAL *x, REAL *xact, REAL *s, REAL *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
    //
    common cmn;
    common_write write(cmn);
    const INTEGER ntran = 3;
    char transs[ntran] = {'N', 'T', 'C'};
    char facts[] = {'F', 'N', 'E'};
    char equeds[] = {'N', 'R', 'C', 'B'};
    char path[4] = {};
    char fact_trans[3];
    char buf[1024];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    INTEGER nb = 0;
    INTEGER nbmin = 0;
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    char xtype;
    const INTEGER ntypes = 11;
    INTEGER nimat = 0;
    INTEGER imat = 0;
    bool zerot = false;
    char type;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist;
    const REAL one = 1.0;
    REAL rcondc = 0.0;
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER ioff = 0;
    const REAL zero = 0.0;
    INTEGER iequed = 0;
    char equed[1];
    INTEGER nfact = 0;
    INTEGER ifact = 0;
    char fact[1];
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
    INTEGER lwork = 0;
    REAL ainvnm = 0.0;
    INTEGER itran = 0;
    char trans[1];
    const INTEGER ntests = 7;
    REAL result[ntests];
    INTEGER nt = 0;
    INTEGER k = 0;
    REAL rcond = 0.0;
    REAL rpvgrw = 0.0;
    INTEGER k1 = 0;
    bool trfcon = false;
    REAL roldc = 0.0;
    static const char *format_9997 = "(1x,a,', FACT=''',a1,''', TRANS=''',a1,''', N=',i5,', EQUED=''',a1,"
                                     "''', type ',i2,', test(',i1,')=',a)";
    static const char *format_9998 = "(1x,a,', FACT=''',a1,''', TRANS=''',a1,''', N=',i5,', type ',i2,"
                                     "', test(',i1,')=',a)";
    //
    path[0] = 'R';
    path[1] = 'G';
    path[2] = 'E';
    path[3] = '\0';
    nrun = 0;
    nfail = 0;
    nerrs = 0;
    //
    //     Test the error exits
    //
    if (tsterr) {
        Rerrvx(path, nout);
    }
    infot = 0;
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
        xtype = 'N';
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
                goto statement_80;
            }
            //
            //           Skip types 5, 6, or 7 if the matrix size is too small.
            //
            zerot = imat >= 5 && imat <= 7;
            if (zerot && n < imat - 4) {
                goto statement_80;
            }
            //
            //           Set up parameters with Rlatb4 and generate a test matrix
            //           with Rlatms.
            //
            Rlatb4(path, imat, n, n, &type, kl, ku, anorm, mode, cndnum, &dist);
            rcondc = one / cndnum;
            //
            strncpy(srnamt, "Rlatms", srnamt_len);
            Rlatms(n, n, &dist, iseed, &type, rwork, mode, cndnum, anorm, kl, ku, "No packing", a, lda, work, info);
            //
            //           Check error code from Rlatms.
            //
            if (info != 0) {
                Alaerh(path, "Rlatms", info, 0, " ", n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                goto statement_80;
            }
            //
            //           For types 5-7, zero one or more columns of the matrix to
            //           test that INFO is returned correctly.
            //
            if (zerot) {
                if (imat == 5) {
                    izero = 1;
                } else if (imat == 6) {
                    izero = n;
                } else {
                    izero = n / 2 + 1;
                }
                ioff = (izero - 1) * lda;
                if (imat < 7) {
                    for (i = 1; i <= n; i = i + 1) {
                        a[(ioff + i) - 1] = zero;
                    }
                } else {
                    Rlaset("Full", n, n - izero + 1, zero, zero, &a[(ioff + 1) - 1], lda);
                }
            } else {
                izero = 0;
            }
            //
            //           Save a copy of the matrix A in ASAV.
            //
            Rlacpy("Full", n, n, a, lda, asav, lda);
            //
            for (iequed = 1; iequed <= 4; iequed = iequed + 1) {
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
                        rcondo = zero;
                        rcondi = zero;
                        //
                    } else if (!nofact) {
                        //
                        //                    Compute the condition number for comparison with
                        //                    the value returned by Rgesvx (FACT = 'N' reuses
                        //                    the condition number from the previous iteration
                        //                    with FACT = 'F').
                        //
                        Rlacpy("Full", n, n, asav, lda, afac, lda);
                        if (equil || iequed > 1) {
                            //
                            //                       Compute row and column scale factors to
                            //                       equilibrate the matrix A.
                            //
                            Rgeequ(n, n, afac, lda, s, &s[(n + 1) - 1], rowcnd, colcnd, amax, info);
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
                                //                          Equilibrate the matrix.
                                //
                                Rlaqge(n, n, afac, lda, s, &s[(n + 1) - 1], rowcnd, colcnd, amax, equed);
                            }
                        }
                        //
                        //                    Save the condition number of the non-equilibrated
                        //                    system for use in Rget04.
                        //
                        if (equil) {
                            roldo = rcondo;
                            roldi = rcondi;
                        }
                        //
                        //                    Compute the 1-norm and infinity-norm of A.
                        //
                        anormo = Rlange("1", n, n, afac, lda, rwork);
                        anormi = Rlange("I", n, n, afac, lda, rwork);
                        //
                        //                    Factor the matrix A.
                        //
                        strncpy(srnamt, "Rgetrf", srnamt_len);
                        Rgetrf(n, n, afac, lda, iwork, info);
                        //
                        //                    Form the inverse of A.
                        //
                        Rlacpy("Full", n, n, afac, lda, a, lda);
                        lwork = nmax * max((INTEGER)3, nrhs);
                        strncpy(srnamt, "Rgetri", srnamt_len);
                        Rgetri(n, a, lda, iwork, work, lwork, info);
                        //
                        //                    Compute the 1-norm condition number of A.
                        //
                        ainvnm = Rlange("1", n, n, a, lda, rwork);
                        if (anormo <= zero || ainvnm <= zero) {
                            rcondo = one;
                        } else {
                            rcondo = (one / anormo) / ainvnm;
                        }
                        //
                        //                    Compute the infinity-norm condition number of A.
                        //
                        ainvnm = Rlange("I", n, n, a, lda, rwork);
                        if (anormi <= zero || ainvnm <= zero) {
                            rcondi = one;
                        } else {
                            rcondi = (one / anormi) / ainvnm;
                        }
                    }
                    //
                    for (itran = 1; itran <= ntran; itran = itran + 1) {
                        //
                        //                    Do for each value of TRANS.
                        //
                        trans[0] = transs[itran - 1];
                        if (itran == 1) {
                            rcondc = rcondo;
                        } else {
                            rcondc = rcondi;
                        }
                        //
                        //                    Restore the matrix A.
                        //
                        Rlacpy("Full", n, n, asav, lda, a, lda);
                        //
                        //                    Form an exact solution and set the right hand side.
                        //
                        strncpy(srnamt, "Rlarhs", srnamt_len);
                        Rlarhs(path, &xtype, "Full", trans, n, n, kl, ku, nrhs, a, lda, xact, lda, b, lda, iseed, info);
                        xtype = 'C';
                        Rlacpy("Full", n, nrhs, b, lda, bsav, lda);
                        //
                        if (nofact && itran == 1) {
                            //
                            //                       --- Test Rgesv  ---
                            //
                            //                       Compute the LU factorization of the matrix and
                            //                       solve the system.
                            //
                            Rlacpy("Full", n, n, a, lda, afac, lda);
                            Rlacpy("Full", n, nrhs, b, lda, x, lda);
                            //
                            strncpy(srnamt, "Rgesv", srnamt_len);
                            Rgesv(n, nrhs, afac, lda, iwork, x, lda, info);
                            //
                            //                       Check error code from Rgesv .
                            //
                            if (info != izero) {
                                Alaerh(path, "Rgesv", info, izero, " ", n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                            }
                            //
                            //                       Reconstruct matrix from factors and compute
                            //                       residual.
                            //
                            Rget01(n, n, a, lda, afac, lda, iwork, rwork, result[1 - 1]);
                            nt = 1;
                            if (izero == 0) {
                                //
                                //                          Compute residual of the computed solution.
                                //
                                Rlacpy("Full", n, nrhs, b, lda, work, lda);
                                Rget02("No transpose", n, n, nrhs, a, lda, x, lda, work, lda, rwork, result[2 - 1]);
                                //
                                //                          Check solution from generated exact solution.
                                //
                                Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                                nt = 3;
                            }
                            //
                            //                       Print information about the tests that did not
                            //                       pass the threshold.
                            //
                            for (k = 1; k <= nt; k = k + 1) {
                                if (result[k - 1] >= thresh) {
                                    if (nfail == 0 && nerrs == 0) {
                                        Aladhd(nout, path);
                                    }
                                    sprintnum_short(buf, result[k - 1]);
                                    write(nout, "(1x,a,', N =',i5,', type ',i2,', test(',i2,') =',a)"), "Rgesv", n, imat, k, buf;
                                    nfail++;
                                }
                            }
                            nrun += nt;
                        }
                        //
                        //                    --- Test Rgesvx ---
                        //
                        if (!prefac) {
                            Rlaset("Full", n, n, zero, zero, afac, lda);
                        }
                        Rlaset("Full", n, nrhs, zero, zero, x, lda);
                        if (iequed > 1 && n > 0) {
                            //
                            //                       Equilibrate the matrix if FACT = 'F' and
                            //                       EQUED = 'R', 'C', or 'B'.
                            //
                            Rlaqge(n, n, a, lda, s, &s[(n + 1) - 1], rowcnd, colcnd, amax, equed);
                        }
                        //
                        //                    Solve the system and compute the condition number
                        //                    and error bounds using Rgesvx.
                        //
                        strncpy(srnamt, "Rgesvx", srnamt_len);
                        Rgesvx(fact, trans, n, nrhs, a, lda, afac, lda, iwork, equed, s, &s[(n + 1) - 1], b, lda, x, lda, rcond, rwork, &rwork[(nrhs + 1) - 1], work, &iwork[(n + 1) - 1], info);
                        //
                        //                    Check the error code from Rgesvx.
                        //
                        if (info != izero) {
                            fact_trans[0] = fact[0];
                            fact_trans[1] = trans[0];
                            fact_trans[2] = '\0';
                            Alaerh(path, "Rgesvx", info, izero, fact_trans, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                        }
                        //
                        //                    Compare WORK(1) from Rgesvx with the computed
                        //                    reciprocal pivot growth factor RPVGRW
                        //
                        if (info != 0 && info <= n) {
                            rpvgrw = Rlantr("M", "U", "N", info, info, afac, lda, work);
                            if (rpvgrw == zero) {
                                rpvgrw = one;
                            } else {
                                rpvgrw = Rlange("M", n, info, a, lda, work) / rpvgrw;
                            }
                        } else {
                            rpvgrw = Rlantr("M", "U", "N", n, n, afac, lda, work);
                            if (rpvgrw == zero) {
                                rpvgrw = one;
                            } else {
                                rpvgrw = Rlange("M", n, n, a, lda, work) / rpvgrw;
                            }
                        }
                        result[7 - 1] = abs(rpvgrw - work[1 - 1]) / max(work[1 - 1], rpvgrw) / Rlamch("E");
                        //
                        if (!prefac) {
                            //
                            //                       Reconstruct matrix from factors and compute
                            //                       residual.
                            //
                            Rget01(n, n, a, lda, afac, lda, iwork, &rwork[(2 * nrhs + 1) - 1], result[1 - 1]);
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
                            Rlacpy("Full", n, nrhs, bsav, lda, work, lda);
                            Rget02(trans, n, n, nrhs, asav, lda, x, lda, work, lda, &rwork[(2 * nrhs + 1) - 1], result[2 - 1]);
                            //
                            //                       Check solution from generated exact solution.
                            //
                            if (nofact || (prefac && Mlsame(equed, "N"))) {
                                Rget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                            } else {
                                if (itran == 1) {
                                    roldc = roldo;
                                } else {
                                    roldc = roldi;
                                }
                                Rget04(n, nrhs, x, lda, xact, lda, roldc, result[3 - 1]);
                            }
                            //
                            //                       Check the error bounds from iterative
                            //                       refinement.
                            //
                            Rget07(trans, n, nrhs, asav, lda, b, lda, x, lda, xact, lda, rwork, true, &rwork[(nrhs + 1) - 1], &result[4 - 1]);
                        } else {
                            trfcon = true;
                        }
                        //
                        //                    Compare RCOND from Rgesvx with the computed value
                        //                    in RCONDC.
                        //
                        result[6 - 1] = Rget06(rcond, rcondc);
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
                                        sprintnum_short(buf, result[k - 1]);
                                        write(nout, format_9997), "Rgesvx", fact, trans, n, equed, imat, k, buf;
                                    } else {
                                        sprintnum_short(buf, result[k - 1]);
                                        write(nout, format_9998), "Rgesvx", fact, trans, n, imat, k, buf;
                                    }
                                    nfail++;
                                }
                            }
                            nrun += ntests - k1 + 1;
                        } else {
                            if (result[1 - 1] >= thresh && !prefac) {
                                if (nfail == 0 && nerrs == 0) {
                                    Aladhd(nout, path);
                                }
                                if (prefac) {
                                    sprintnum_short(buf, result[1 - 1]);
                                    write(nout, format_9997), "Rgesvx", fact, trans, n, equed, imat, 1, buf;
                                } else {
                                    sprintnum_short(buf, result[1 - 1]);
                                    write(nout, format_9998), "Rgesvx", fact, trans, n, imat, 1, buf;
                                }
                                nfail++;
                                nrun++;
                            }
                            if (result[6 - 1] >= thresh) {
                                if (nfail == 0 && nerrs == 0) {
                                    Aladhd(nout, path);
                                }
                                if (prefac) {
                                    sprintnum_short(buf, result[6 - 1]);
                                    write(nout, format_9997), "Rgesvx", fact, trans, n, equed, imat, 6, buf;
                                } else {
                                    sprintnum_short(buf, result[6 - 1]);
                                    write(nout, format_9998), "Rgesvx", fact, trans, n, imat, 6, buf;
                                }
                                nfail++;
                                nrun++;
                            }
                            if (result[7 - 1] >= thresh) {
                                if (nfail == 0 && nerrs == 0) {
                                    Aladhd(nout, path);
                                }
                                if (prefac) {
                                    sprintnum_short(buf, result[7 - 1]);
                                    write(nout, format_9997), "Rgesvx", fact, trans, n, equed, imat, 7, buf;
                                } else {
                                    sprintnum_short(buf, result[7 - 1]);
                                    write(nout, format_9998), "Rgesvx", fact, trans, n, imat, 7, buf;
                                }
                                nfail++;
                                nrun++;
                            }
                            //
                        }
                        //
                    }
                statement_60:;
                }
            }
        statement_80:;
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasvm(path, nout, nfail, nrun, nerrs);
    //
    //     End of Rdrvge
    //
}
