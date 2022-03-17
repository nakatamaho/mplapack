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

void Rchkgb(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nn, INTEGER *nval, INTEGER const nnb, INTEGER *nbval, INTEGER const nns, INTEGER *nsval, REAL const thresh, bool const tsterr, REAL *a, INTEGER const la, REAL *afac, INTEGER const lafac, REAL *b, REAL *x, REAL *xact, REAL *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    //
    const INTEGER ntran = 3;
    char transs[ntran] = {'N', 'T', 'C'};
    char path[4];
    char buf[1024];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    const INTEGER nbw = 4;
    INTEGER klval[nbw];
    INTEGER kuval[nbw];
    INTEGER im = 0;
    INTEGER m = 0;
    INTEGER in = 0;
    INTEGER n = 0;
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
    INTEGER ldafac = 0;
    INTEGER imat = 0;
    bool zerot = false;
    char type;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist;
    INTEGER koff = 0;
    const REAL zero = 0.0;
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER i2 = 0;
    INTEGER i1 = 0;
    INTEGER ioff = 0;
    INTEGER j = 0;
    INTEGER inb = 0;
    INTEGER nb = 0;
    bool trfcon = false;
    const INTEGER ntests = 7;
    REAL result[ntests];
    REAL anormo = 0.0;
    REAL anormi = 0.0;
    INTEGER ldb = 0;
    const REAL one = 1.0;
    REAL ainvnm = 0.0;
    REAL rcondo = 0.0;
    REAL rcondi = 0.0;
    INTEGER irhs = 0;
    INTEGER nrhs = 0;
    INTEGER itran = 0;
    char trans;
    REAL rcondc = 0.0;
    char norm;
    INTEGER k = 0;
    REAL rcond = 0.0;
    //
    //     Initialize constants and the random number seed.
    //
    path[0] = 'R';
    path[1] = 'G';
    path[2] = 'B';
    path[3] = '\0';
    nrun = 0;
    nfail = 0;
    nerrs = 0;
    //
    //     Test the error exits
    //
    if (tsterr) {
        Rerrge(path, nout);
    }
    infot = 0;
    xlaenv(2, 2);
    //
    //     Initialize the first value for the lower and upper bandwidths.
    //
    klval[1 - 1] = 0;
    kuval[1 - 1] = 0;
    //
    //     Do for each value of M in MVAL
    //
    for (im = 1; im <= nm; im = im + 1) {
        m = mval[im - 1];
        //
        //        Set values to use for the lower bandwidth.
        //
        klval[2 - 1] = m + (m + 1) / 4;
        //
        //        KLVAL( 2 ) = MAX( M-1, 0 )
        //
        klval[3 - 1] = (3 * m - 1) / 4;
        klval[4 - 1] = (m + 1) / 4;
        //
        //        Do for each value of N in NVAL
        //
        for (in = 1; in <= nn; in = in + 1) {
            n = nval[in - 1];
            xtype = 'N';
            //
            //           Set values to use for the upper bandwidth.
            //
            kuval[2 - 1] = n + (n + 1) / 4;
            //
            //           KUVAL( 2 ) = MAX( N-1, 0 )
            //
            kuval[3 - 1] = (3 * n - 1) / 4;
            kuval[4 - 1] = (n + 1) / 4;
            //
            //           Set limits on the number of loop iterations.
            //
            nkl = min(m + 1, (INTEGER)4);
            if (n == 0) {
                nkl = 2;
            }
            nku = min(n + 1, (INTEGER)4);
            if (m == 0) {
                nku = 2;
            }
            nimat = ntypes;
            if (m <= 0 || n <= 0) {
                nimat = 1;
            }
            //
            for (ikl = 1; ikl <= nkl; ikl = ikl + 1) {
                //
                //              Do for KL = 0, (5*M+1)/4, (3M-1)/4, and (M+1)/4. This
                //              order makes it easier to skip redundant values for small
                //              values of M.
                //
                kl = klval[ikl - 1];
                for (iku = 1; iku <= nku; iku = iku + 1) {
                    //
                    //                 Do for KU = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This
                    //                 order makes it easier to skip redundant values for
                    //                 small values of N.
                    //
                    ku = kuval[iku - 1];
                    //
                    //                 Check that A and AFAC are big enough to generate this
                    //                 matrix.
                    //
                    lda = kl + ku + 1;
                    ldafac = 2 * kl + ku + 1;
                    if ((lda * n) > la || (ldafac * n) > lafac) {
                        if (nfail == 0 && nerrs == 0) {
                            Alahd(nout, path);
                        }
                        if (n * (kl + ku + 1) > la) {
                            write(nout, "(' *** In Rchkgb, LA=',i5,' is too small for M=',i5,', N=',"
                                        "i5,', KL=',i4,', KU=',i4,/,' ==> Increase LA to at least ',"
                                        "i5)"),
                                la, m, n, kl, ku, n *(kl + ku + 1);
                            nerrs++;
                        }
                        if (n * (2 * kl + ku + 1) > lafac) {
                            write(nout, "(' *** In Rchkgb, LAFAC=',i5,' is too small for M=',i5,"
                                        "', N=',i5,', KL=',i4,', KU=',i4,/,"
                                        "' ==> Increase LAFAC to at least ',i5)"),
                                lafac, m, n, kl, ku, n *(2 * kl + ku + 1);
                            nerrs++;
                        }
                        goto statement_130;
                    }
                    //
                    for (imat = 1; imat <= nimat; imat = imat + 1) {
                        //
                        //                    Do the tests only if DOTYPE( IMAT ) is true.
                        //
                        if (!dotype[imat - 1]) {
                            goto statement_120;
                        }
                        //
                        //                    Skip types 2, 3, or 4 if the matrix size is too
                        //                    small.
                        //
                        zerot = imat >= 2 && imat <= 4;
                        if (zerot && n < imat - 1) {
                            goto statement_120;
                        }
                        //
                        if (!zerot || !dotype[1 - 1]) {
                            //
                            //                       Set up parameters with Rlatb4 and generate a
                            //                       test matrix with Rlatms.
                            //
                            Rlatb4(path, imat, m, n, &type, kl, ku, anorm, mode, cndnum, &dist);
                            //
                            koff = max((INTEGER)1, ku + 2 - n);
                            for (i = 1; i <= koff - 1; i = i + 1) {
                                a[i - 1] = zero;
                            }
                            strncpy(srnamt, "Rlatms", srnamt_len);
                            Rlatms(m, n, &dist, iseed, &type, rwork, mode, cndnum, anorm, kl, ku, "Z", &a[koff - 1], lda, work, info);
                            //
                            //                       Check the error code from Rlatms.
                            //
                            if (info != 0) {
                                Alaerh(path, "Rlatms", info, 0, " ", m, n, kl, ku, -1, imat, nfail, nerrs, nout);
                                goto statement_120;
                            }
                        } else if (izero > 0) {
                            //
                            //                       Use the same matrix for types 3 and 4 as for
                            //                       type 2 by copying back the zeroed out column.
                            //
                            Rcopy(i2 - i1 + 1, b, 1, &a[(ioff + i1) - 1], 1);
                        }
                        //
                        //                    For types 2, 3, and 4, zero one or more columns of
                        //                    the matrix to test that INFO is returned correctly.
                        //
                        izero = 0;
                        if (zerot) {
                            if (imat == 2) {
                                izero = 1;
                            } else if (imat == 3) {
                                izero = min(m, n);
                            } else {
                                izero = min(m, n) / 2 + 1;
                            }
                            ioff = (izero - 1) * lda;
                            if (imat < 4) {
                                //
                                //                          Store the column to be zeroed out in B.
                                //
                                i1 = max((INTEGER)1, ku + 2 - izero);
                                i2 = min(kl + ku + 1, ku + 1 + (m - izero));
                                Rcopy(i2 - i1 + 1, &a[(ioff + i1) - 1], 1, b, 1);
                                //
                                for (i = i1; i <= i2; i = i + 1) {
                                    a[(ioff + i) - 1] = zero;
                                }
                            } else {
                                for (j = izero; j <= n; j = j + 1) {
                                    for (i = max((INTEGER)1, ku + 2 - j); i <= min(kl + ku + 1, ku + 1 + (m - j)); i = i + 1) {
                                        a[(ioff + i) - 1] = zero;
                                    }
                                    ioff += lda;
                                }
                            }
                        }
                        //
                        //                    These lines, if used in place of the calls in the
                        //                    loop over INB, cause the code to bomb on a Sun
                        //                    SPARCstation.
                        //
                        //                     ANORMO = Rlangb( 'O', N, KL, KU, A, LDA, RWORK )
                        //                     ANORMI = Rlangb( 'I', N, KL, KU, A, LDA, RWORK )
                        //
                        //                    Do for each blocksize in NBVAL
                        //
                        for (inb = 1; inb <= nnb; inb = inb + 1) {
                            nb = nbval[inb - 1];
                            xlaenv(1, nb);
                            //
                            //                       Compute the LU factorization of the band matrix.
                            //
                            if (m > 0 && n > 0) {
                                Rlacpy("Full", kl + ku + 1, n, a, lda, &afac[(kl + 1) - 1], ldafac);
                            }
                            strncpy(srnamt, "Rgbtrf", srnamt_len);
                            Rgbtrf(m, n, kl, ku, afac, ldafac, iwork, info);
                            //
                            //                       Check error code from Rgbtrf.
                            //
                            if (info != izero) {
                                Alaerh(path, "Rgbtrf", info, izero, " ", m, n, kl, ku, nb, imat, nfail, nerrs, nout);
                            }
                            trfcon = false;
                            //
                            //+    TEST 1
                            //                       Reconstruct matrix from factors and compute
                            //                       residual.
                            //
                            Rgbt01(m, n, kl, ku, a, lda, afac, ldafac, iwork, work, result[1 - 1]);
                            //
                            //                       Print information about the tests so far that
                            //                       did not pass the threshold.
                            //
                            if (result[1 - 1] >= thresh) {
                                if (nfail == 0 && nerrs == 0) {
                                    Alahd(nout, path);
                                }
                                sprintnum_short(buf, result[1 - 1]);
                                write(nout, "(' M =',i5,', N =',i5,', KL=',i5,', KU=',i5,', NB =',i4,"
                                            "', type ',i1,', test(',i1,')=',a)"),
                                    m, n, kl, ku, nb, imat, 1, buf;
                                nfail++;
                            }
                            nrun++;
                            //
                            //                       Skip the remaining tests if this is not the
                            //                       first block size or if M .ne. N.
                            //
                            if (inb > 1 || m != n) {
                                goto statement_110;
                            }
                            //
                            anormo = Rlangb("O", n, kl, ku, a, lda, rwork);
                            anormi = Rlangb("I", n, kl, ku, a, lda, rwork);
                            //
                            if (info == 0) {
                                //
                                //                          Form the inverse of A so we can get a good
                                //                          estimate of CNDNUM = norm(A) * norm(inv(A)).
                                //
                                ldb = max((INTEGER)1, n);
                                Rlaset("Full", n, n, zero, one, work, ldb);
                                strncpy(srnamt, "Rgbtrs", srnamt_len);
                                Rgbtrs("No transpose", n, kl, ku, n, afac, ldafac, iwork, work, ldb, info);
                                //
                                //                          Compute the 1-norm condition number of A.
                                //
                                ainvnm = Rlange("O", n, n, work, ldb, rwork);
                                if (anormo <= zero || ainvnm <= zero) {
                                    rcondo = one;
                                } else {
                                    rcondo = (one / anormo) / ainvnm;
                                }
                                //
                                //                          Compute the infinity-norm condition number of
                                //                          A.
                                //
                                ainvnm = Rlange("I", n, n, work, ldb, rwork);
                                if (anormi <= zero || ainvnm <= zero) {
                                    rcondi = one;
                                } else {
                                    rcondi = (one / anormi) / ainvnm;
                                }
                            } else {
                                //
                                //                          Do only the condition estimate if INFO.NE.0.
                                //
                                trfcon = true;
                                rcondo = zero;
                                rcondi = zero;
                            }
                            //
                            //                       Skip the solve tests if the matrix is singular.
                            //
                            if (trfcon) {
                                goto statement_90;
                            }
                            //
                            for (irhs = 1; irhs <= nns; irhs = irhs + 1) {
                                nrhs = nsval[irhs - 1];
                                xtype = 'N';
                                //
                                for (itran = 1; itran <= ntran; itran = itran + 1) {
                                    trans = transs[itran - 1];
                                    if (itran == 1) {
                                        rcondc = rcondo;
                                        norm = 'O';
                                    } else {
                                        rcondc = rcondi;
                                        norm = 'I';
                                    }
                                    //
                                    //+    TEST 2:
                                    //                             Solve and compute residual for A * X = B.
                                    //
                                    strncpy(srnamt, "Rlarhs", srnamt_len);
                                    Rlarhs(path, &xtype, " ", &trans, n, n, kl, ku, nrhs, a, lda, xact, ldb, b, ldb, iseed, info);
                                    xtype = 'C';
                                    Rlacpy("Full", n, nrhs, b, ldb, x, ldb);
                                    //
                                    strncpy(srnamt, "Rgbtrs", srnamt_len);
                                    Rgbtrs(&trans, n, kl, ku, nrhs, afac, ldafac, iwork, x, ldb, info);
                                    //
                                    //                             Check error code from Rgbtrs.
                                    //
                                    if (info != 0) {
                                        Alaerh(path, "Rgbtrs", info, 0, &trans, n, n, kl, ku, -1, imat, nfail, nerrs, nout);
                                    }
                                    //
                                    Rlacpy("Full", n, nrhs, b, ldb, work, ldb);
                                    Rgbt02(&trans, m, n, kl, ku, nrhs, a, lda, x, ldb, work, ldb, result[2 - 1]);
                                    //
                                    //+    TEST 3:
                                    //                             Check solution from generated exact
                                    //                             solution.
                                    //
                                    Rget04(n, nrhs, x, ldb, xact, ldb, rcondc, result[3 - 1]);
                                    //
                                    //+    TESTS 4, 5, 6:
                                    //                             Use iterative refinement to improve the
                                    //                             solution.
                                    //
                                    strncpy(srnamt, "Rgbrfs", srnamt_len);
                                    Rgbrfs(&trans, n, kl, ku, nrhs, a, lda, afac, ldafac, iwork, b, ldb, x, ldb, rwork, &rwork[(nrhs + 1) - 1], work, &iwork[(n + 1) - 1], info);
                                    //
                                    //                             Check error code from Rgbrfs.
                                    //
                                    if (info != 0) {
                                        Alaerh(path, "Rgbrfs", info, 0, &trans, n, n, kl, ku, nrhs, imat, nfail, nerrs, nout);
                                    }
                                    //
                                    Rget04(n, nrhs, x, ldb, xact, ldb, rcondc, result[4 - 1]);
                                    Rgbt05(&trans, n, kl, ku, nrhs, a, lda, b, ldb, x, ldb, xact, ldb, rwork, &rwork[(nrhs + 1) - 1], &result[5 - 1]);
                                    for (k = 2; k <= 6; k = k + 1) {
                                        if (result[k - 1] >= thresh) {
                                            if (nfail == 0 && nerrs == 0) {
                                                Alahd(nout, path);
                                            }
                                            sprintnum_short(buf, result[k - 1]);
                                            write(nout, "(' TRANS=''',a1,''', N=',i5,', KL=',i5,', KU=',i5,"
                                                        "', NRHS=',i3,', type ',i1,', test(',i1,')=',a)"),
                                                trans, n, kl, ku, nrhs, imat, k, buf;
                                            nfail++;
                                        }
                                    }
                                    nrun += 5;
                                }
                            }
                        //
                        //+    TEST 7:
                        //                          Get an estimate of RCOND = 1/CNDNUM.
                        //
                        statement_90:
                            for (itran = 1; itran <= 2; itran = itran + 1) {
                                if (itran == 1) {
                                    anorm = anormo;
                                    rcondc = rcondo;
                                    norm = 'O';
                                } else {
                                    anorm = anormi;
                                    rcondc = rcondi;
                                    norm = 'I';
                                }
                                strncpy(srnamt, "Rgbcon", srnamt_len);
                                Rgbcon(&norm, n, kl, ku, afac, ldafac, iwork, anorm, rcond, work, &iwork[(n + 1) - 1], info);
                                //
                                //                             Check error code from Rgbcon.
                                //
                                if (info != 0) {
                                    Alaerh(path, "Rgbcon", info, 0, &norm, n, n, kl, ku, -1, imat, nfail, nerrs, nout);
                                }
                                //
                                result[7 - 1] = Rget06(rcond, rcondc);
                                //
                                //                          Print information about the tests that did
                                //                          not pass the threshold.
                                //
                                if (result[7 - 1] >= thresh) {
                                    if (nfail == 0 && nerrs == 0) {
                                        Alahd(nout, path);
                                    }
                                    sprintnum_short(buf, result[7 - 1]);
                                    write(nout, "(' NORM =''',a1,''', N=',i5,', KL=',i5,', KU=',i5,',',"
                                                "10x,' type ',i1,', test(',i1,')=',a)"),
                                        norm, n, kl, ku, imat, 7, buf;
                                    nfail++;
                                }
                                nrun++;
                            }
                        //
                        statement_110:;
                        }
                    statement_120:;
                    }
                statement_130:;
                }
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Rchkgb
    //
}
