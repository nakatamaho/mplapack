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

void Cchkgt(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, REAL const thresh, bool const tsterr, COMPLEX *a, COMPLEX *af, COMPLEX *b, COMPLEX *x, COMPLEX *xact, COMPLEX *work, REAL *rwork, INTEGER *iwork, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    char transs[] = {'N', 'T', 'C'};
    char path[3];
    char buf[1024];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER m = 0;
    INTEGER lda = 0;
    const INTEGER ntypes = 12;
    INTEGER nimat = 0;
    INTEGER imat = 0;
    char type[1];
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cond = 0.0;
    char dist[1];
    bool zerot = false;
    INTEGER koff = 0;
    INTEGER info = 0;
    INTEGER izero = 0;
    const REAL one = 1.0;
    COMPLEX z[3];
    const REAL zero = 0.0;
    bool trfcon = false;
    const INTEGER ntests = 7;
    REAL result[ntests];
    INTEGER itran = 0;
    char trans[1];
    char norm[1];
    REAL ainvnm = 0.0;
    INTEGER j = 0;
    REAL rcondc = 0.0;
    REAL rcondo = 0.0;
    REAL rcondi = 0.0;
    REAL rcond = 0.0;
    INTEGER irhs = 0;
    INTEGER nrhs = 0;
    INTEGER ix = 0;
    INTEGER k = 0;
    //
    path[0] = 'C';
    path[1] = 'G';
    path[2] = 'T';
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
        Cerrge(path, nout);
    }
    infot = 0;
    //
    for (in = 1; in <= nn; in = in + 1) {
        //
        //        Do for each value of N in NVAL.
        //
        n = nval[in - 1];
        m = max(n - 1, 0);
        lda = max((INTEGER)1, n);
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
                goto statement_100;
            }
            //
            //           Set up parameters with Clatb4.
            //
            Clatb4(path, imat, n, n, type, kl, ku, anorm, mode, cond, dist);
            //
            zerot = imat >= 8 && imat <= 10;
            if (imat <= 6) {
                //
                //              Types 1-6:  generate matrices of known condition number.
                //
                koff = max({(INTEGER)2 - ku, 3 - max((INTEGER)1, n)});
                strncpy(srnamt, "Clatms", srnamt_len);
                Clatms(n, n, dist, iseed, type, rwork, mode, cond, anorm, kl, ku, "Z", &af[koff - 1], 3, work, info);
                //
                //              Check the error code from Clatms.
                //
                if (info != 0) {
                    Alaerh(path, "Clatms", info, 0, " ", n, n, kl, ku, -1, imat, nfail, nerrs, nout);
                    goto statement_100;
                }
                izero = 0;
                //
                if (n > 1) {
                    Ccopy(n - 1, &af[4 - 1], 3, a, 1);
                    Ccopy(n - 1, &af[3 - 1], 3, &a[(n + m + 1) - 1], 1);
                }
                Ccopy(n, &af[2 - 1], 3, &a[(m + 1) - 1], 1);
            } else {
                //
                //              Types 7-12:  generate tridiagonal matrices with
                //              unknown condition numbers.
                //
                if (!zerot || !dotype[7 - 1]) {
                    //
                    //                 Generate a matrix with elements whose real and
                    //                 imaginary parts are from [-1,1].
                    //
                    Clarnv(2, iseed, n + 2 * m, a);
                    if (anorm != one) {
                        CRscal(n + 2 * m, anorm, a, 1);
                    }
                } else if (izero > 0) {
                    //
                    //                 Reuse the last matrix by copying back the zeroed out
                    //                 elements.
                    //
                    if (izero == 1) {
                        a[n - 1] = z[2 - 1];
                        if (n > 1) {
                            a[1 - 1] = z[3 - 1];
                        }
                    } else if (izero == n) {
                        a[(3 * n - 2) - 1] = z[1 - 1];
                        a[(2 * n - 1) - 1] = z[2 - 1];
                    } else {
                        a[(2 * n - 2 + izero) - 1] = z[1 - 1];
                        a[(n - 1 + izero) - 1] = z[2 - 1];
                        a[izero - 1] = z[3 - 1];
                    }
                }
                //
                //              If IMAT > 7, set one column of the matrix to 0.
                //
                if (!zerot) {
                    izero = 0;
                } else if (imat == 8) {
                    izero = 1;
                    z[2 - 1] = a[n - 1];
                    a[n - 1] = zero;
                    if (n > 1) {
                        z[3 - 1] = a[1 - 1];
                        a[1 - 1] = zero;
                    }
                } else if (imat == 9) {
                    izero = n;
                    z[1 - 1] = a[(3 * n - 2) - 1];
                    z[2 - 1] = a[(2 * n - 1) - 1];
                    a[(3 * n - 2) - 1] = zero;
                    a[(2 * n - 1) - 1] = zero;
                } else {
                    izero = (n + 1) / 2;
                    for (i = izero; i <= n - 1; i = i + 1) {
                        a[(2 * n - 2 + i) - 1] = zero;
                        a[(n - 1 + i) - 1] = zero;
                        a[i - 1] = zero;
                    }
                    a[(3 * n - 2) - 1] = zero;
                    a[(2 * n - 1) - 1] = zero;
                }
            }
            //
            //+    TEST 1
            //           Factor A as L*U and compute the ratio
            //              norm(L*U - A) / (n * norm(A) * EPS )
            //
            Ccopy(n + 2 * m, a, 1, af, 1);
            strncpy(srnamt, "Cgttrf", srnamt_len);
            Cgttrf(n, af, &af[(m + 1) - 1], &af[(n + m + 1) - 1], &af[(n + 2 * m + 1) - 1], iwork, info);
            //
            //           Check error code from Cgttrf.
            //
            if (info != izero) {
                Alaerh(path, "Cgttrf", info, izero, " ", n, n, 1, 1, -1, imat, nfail, nerrs, nout);
            }
            trfcon = info != 0;
            //
            Cgtt01(n, a, &a[(m + 1) - 1], &a[(n + m + 1) - 1], af, &af[(m + 1) - 1], &af[(n + m + 1) - 1], &af[(n + 2 * m + 1) - 1], iwork, work, lda, rwork, result[1 - 1]);
            //
            //           Print the test ratio if it is .GE. THRESH.
            //
            if (result[1 - 1] >= thresh) {
                if (nfail == 0 && nerrs == 0) {
                    Alahd(nout, path);
                }
                sprintnum_short(buf, result[1 - 1]);
                write(nout, "(12x,'N =',i5,',',10x,' type ',i2,', test(',i2,') = ',a)"), n, imat, 1, buf;
                nfail++;
            }
            nrun++;
            //
            for (itran = 1; itran <= 2; itran = itran + 1) {
                trans[0] = transs[itran - 1];
                if (itran == 1) {
                    norm[0] = 'O';
                } else {
                    norm[0] = 'I';
                }
                anorm = Clangt(norm, n, a, &a[(m + 1) - 1], &a[(n + m + 1) - 1]);
                //
                if (!trfcon) {
                    //
                    //                 Use Cgttrs to solve for one column at a time of
                    //                 inv(A), computing the maximum column sum as we go.
                    //
                    ainvnm = zero;
                    for (i = 1; i <= n; i = i + 1) {
                        for (j = 1; j <= n; j = j + 1) {
                            x[j - 1] = zero;
                        }
                        x[i - 1] = one;
                        Cgttrs(trans, n, 1, af, &af[(m + 1) - 1], &af[(n + m + 1) - 1], &af[(n + 2 * m + 1) - 1], iwork, x, lda, info);
                        ainvnm = max({ainvnm, RCasum(n, x, 1)});
                    }
                    //
                    //                 Compute RCONDC = 1 / (norm(A) * norm(inv(A))
                    //
                    if (anorm <= zero || ainvnm <= zero) {
                        rcondc = one;
                    } else {
                        rcondc = (one / anorm) / ainvnm;
                    }
                    if (itran == 1) {
                        rcondo = rcondc;
                    } else {
                        rcondi = rcondc;
                    }
                } else {
                    rcondc = zero;
                }
                //
                //+    TEST 7
                //              Estimate the reciprocal of the condition number of the
                //              matrix.
                //
                strncpy(srnamt, "Cgtcon", srnamt_len);
                Cgtcon(norm, n, af, &af[(m + 1) - 1], &af[(n + m + 1) - 1], &af[(n + 2 * m + 1) - 1], iwork, anorm, rcond, work, info);
                //
                //              Check error code from Cgtcon.
                //
                if (info != 0) {
                    Alaerh(path, "Cgtcon", info, 0, norm, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                }
                //
                result[7 - 1] = Rget06(rcond, rcondc);
                //
                //              Print the test ratio if it is .GE. THRESH.
                //
                if (result[7 - 1] >= thresh) {
                    if (nfail == 0 && nerrs == 0) {
                        Alahd(nout, path);
                    }
                    sprintnum_short(buf, result[7 - 1]);
                    write(nout, "(' NORM =''',a1,''', N =',i5,',',10x,' type ',i2,', test(',i2,"
                                "') = ',a)"),
                        norm, n, imat, 7, buf;
                    nfail++;
                }
                nrun++;
            }
            //
            //           Skip the remaining tests if the matrix is singular.
            //
            if (trfcon) {
                goto statement_100;
            }
            //
            for (irhs = 1; irhs <= nns; irhs = irhs + 1) {
                nrhs = nsval[irhs - 1];
                //
                //              Generate NRHS random solution vectors.
                //
                ix = 1;
                for (j = 1; j <= nrhs; j = j + 1) {
                    Clarnv(2, iseed, n, &xact[ix - 1]);
                    ix += lda;
                }
                //
                for (itran = 1; itran <= 3; itran = itran + 1) {
                    trans[0] = transs[itran - 1];
                    if (itran == 1) {
                        rcondc = rcondo;
                    } else {
                        rcondc = rcondi;
                    }
                    //
                    //                 Set the right hand side.
                    //
                    Clagtm(trans, n, nrhs, one, a, &a[(m + 1) - 1], &a[(n + m + 1) - 1], xact, lda, zero, b, lda);
                    //
                    //+    TEST 2
                    //              Solve op(A) * X = B and compute the residual.
                    //
                    Clacpy("Full", n, nrhs, b, lda, x, lda);
                    strncpy(srnamt, "Cgttrs", srnamt_len);
                    Cgttrs(trans, n, nrhs, af, &af[(m + 1) - 1], &af[(n + m + 1) - 1], &af[(n + 2 * m + 1) - 1], iwork, x, lda, info);
                    //
                    //              Check error code from Cgttrs.
                    //
                    if (info != 0) {
                        Alaerh(path, "Cgttrs", info, 0, trans, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                    }
                    //
                    Clacpy("Full", n, nrhs, b, lda, work, lda);
                    Cgtt02(trans, n, nrhs, a, &a[(m + 1) - 1], &a[(n + m + 1) - 1], x, lda, work, lda, result[2 - 1]);
                    //
                    //+    TEST 3
                    //              Check solution from generated exact solution.
                    //
                    Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[3 - 1]);
                    //
                    //+    TESTS 4, 5, and 6
                    //              Use iterative refinement to improve the solution.
                    //
                    strncpy(srnamt, "Cgtrfs", srnamt_len);
                    Cgtrfs(trans, n, nrhs, a, &a[(m + 1) - 1], &a[(n + m + 1) - 1], af, &af[(m + 1) - 1], &af[(n + m + 1) - 1], &af[(n + 2 * m + 1) - 1], iwork, b, lda, x, lda, rwork, &rwork[(nrhs + 1) - 1], work, &rwork[(2 * nrhs + 1) - 1], info);
                    //
                    //              Check error code from Cgtrfs.
                    //
                    if (info != 0) {
                        Alaerh(path, "Cgtrfs", info, 0, trans, n, n, -1, -1, nrhs, imat, nfail, nerrs, nout);
                    }
                    //
                    Cget04(n, nrhs, x, lda, xact, lda, rcondc, result[4 - 1]);
                    Cgtt05(trans, n, nrhs, a, &a[(m + 1) - 1], &a[(n + m + 1) - 1], b, lda, x, lda, xact, lda, rwork, &rwork[(nrhs + 1) - 1], &result[5 - 1]);
                    //
                    //              Print information about the tests that did not pass the
                    //              threshold.
                    //
                    for (k = 2; k <= 6; k = k + 1) {
                        if (result[k - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Alahd(nout, path);
                            }
                            sprintnum_short(buf, result[k - 1]);
                            write(nout, "(' TRANS=''',a1,''', N =',i5,', NRHS=',i3,', type ',i2,"
                                        "', test(',i2,') = ',a)"),
                                trans, n, nrhs, imat, k, buf;
                            nfail++;
                        }
                    }
                    nrun += 5;
                }
            }
        statement_100:;
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Cchkgt
    //
}
