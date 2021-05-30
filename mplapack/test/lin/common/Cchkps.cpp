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

void Cchkps(bool *dotype, INTEGER const nn, INTEGER *nval, INTEGER const nnb, INTEGER *nbval, INTEGER const nrank, INTEGER *rankval, REAL const thresh, bool const tsterr, INTEGER const /* nmax */, COMPLEX *a, COMPLEX *afac, COMPLEX *perm, INTEGER *piv, COMPLEX *work, REAL *rwork, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    //
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    const char uplos[] = {'U', 'L'};
    char path[3];
    char buf[1024];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    INTEGER in = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    const INTEGER ntypes = 9;
    INTEGER nimat = 0;
    INTEGER izero = 0;
    INTEGER imat = 0;
    INTEGER irank = 0;
    INTEGER rank = 0;
    INTEGER iuplo = 0;
    char uplo[1];
    char type[1];
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist[1];
    INTEGER info = 0;
    INTEGER inb = 0;
    INTEGER nb = 0;
    const REAL one = 1.0f;
    REAL tol = 0.0;
    INTEGER comprank = 0;
    REAL result = 0.0;
    INTEGER rankdiff = 0;
    //
    path[0] = 'Z';
    path[1] = 'P';
    path[2] = 'S';
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
        Cerrps(path, nout);
    }
    infot = 0;
    //
    //     Do for each value of N in NVAL
    //
    for (in = 1; in <= nn; in = in + 1) {
        n = nval[in - 1];
        lda = max(n, 1);
        nimat = ntypes;
        if (n <= 0) {
            nimat = 1;
        }
        //
        izero = 0;
        for (imat = 1; imat <= nimat; imat = imat + 1) {
            //
            //           Do the tests only if DOTYPE( IMAT ) is true.
            //
            if (!dotype[imat - 1]) {
                goto statement_140;
            }
            //
            //              Do for each value of RANK in RANKVAL
            //
            for (irank = 1; irank <= nrank; irank = irank + 1) {
                //
                //              Only repeat test 3 to 5 for different ranks
                //              Other tests use full rank
                //
                if ((imat < 3 || imat > 5) && irank > 1) {
                    goto statement_130;
                }
                //
                rank = ceil((n * castREAL(rankval[irank - 1])) / 100.0);
                //
                //           Do first for UPLO = 'U', then for UPLO = 'L'
                //
                for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                    uplo[0] = uplos[iuplo - 1];
                    //
                    //              Set up parameters with Clatb5 and generate a test matrix
                    //              with Clatmt.
                    //
                    Clatb5(path, imat, n, type, kl, ku, anorm, mode, cndnum, dist);
                    //
                    strncpy(srnamt, "Clatmt", srnamt_len);
                    Clatmt(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, rank, kl, ku, uplo, a, lda, work, info);
                    //
                    //              Check error code from Clatmt.
                    //
                    if (info != 0) {
                        Alaerh(path, "Clatmt", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                        goto statement_120;
                    }
                    //
                    //              Do for each value of NB in NBVAL
                    //
                    for (inb = 1; inb <= nnb; inb = inb + 1) {
                        nb = nbval[inb - 1];
                        xlaenv(1, nb);
                        //
                        //                 Compute the pivoted L*L' or U'*U factorization
                        //                 of the matrix.
                        //
                        Clacpy(uplo, n, n, a, lda, afac, lda);
                        strncpy(srnamt, "Cpstrf", srnamt_len);
                        //
                        //                 Use default tolerance
                        //
                        tol = -one;
                        Cpstrf(uplo, n, afac, lda, piv, comprank, tol, rwork, info);
                        //
                        //                 Check error code from Cpstrf.
                        //
                        if ((info < izero) || (info != izero && rank == n) || (info <= izero && rank < n)) {
                            Alaerh(path, "Cpstrf", info, izero, uplo, n, n, -1, -1, nb, imat, nfail, nerrs, nout);
                            goto statement_110;
                        }
                        //
                        //                 Skip the test if INFO is not 0.
                        //
                        if (info != 0) {
                            goto statement_110;
                        }
                        //
                        //                 Reconstruct matrix from factors and compute residual.
                        //
                        //                 PERM holds permuted L*L^T or U^T*U
                        //
                        Cpst01(uplo, n, a, lda, afac, lda, perm, lda, piv, rwork, result, comprank);
                        //
                        //                 Print information about the tests that did not pass
                        //                 the threshold or where computed rank was not RANK.
                        //
                        if (n == 0) {
                            comprank = 0;
                        }
                        rankdiff = rank - comprank;
                        if (result >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                Alahd(nout, path);
                            }
                            sprintnum_short(buf, result);
                            write(nout, "(' UPLO = ''',a1,''', N =',i5,', RANK =',i3,', Diff =',i5,"
                                        "', NB =',i4,', type ',i2,', Ratio =',a)"),
                                uplo, n, rank, rankdiff, nb, imat, buf;
                            nfail++;
                        }
                        nrun++;
                    statement_110:;
                    }
                //
                statement_120:;
                }
            statement_130:;
            }
        statement_140:;
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Cchkps
    //
}
