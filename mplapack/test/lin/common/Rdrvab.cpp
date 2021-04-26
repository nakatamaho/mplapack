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

void Rdrvab(common &cmn, bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nns, INTEGER *nsval, REAL const thresh, INTEGER const  /* nmax */, REAL *a, REAL *afac, REAL *b, REAL *x, REAL *work, REAL *rwork, arr_cref<float> swork, INTEGER *iwork, INTEGER const nout) {
    FEM_CMN_SVE(Rdrvab);
    common_write write(cmn);
    str<32> &srnamt = cmn.srnamt;
    //
    if (is_called_first_time) {
        static const INTEGER values[] = {2006, 2007, 2008, 2009};
        data_of_type<int>(FEM_VALUES_AND_SIZE), iseedy;
    }
    INTEGER kase = 0;
    str<3> path = char0;
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    arr_1d<4, int> iseed(fill0);
    INTEGER im = 0;
    INTEGER m = 0;
    INTEGER lda = 0;
    INTEGER n = 0;
    const INTEGER ntypes = 11;
    INTEGER nimat = 0;
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
    INTEGER ioff = 0;
    const REAL zero = 0.0;
    INTEGER irhs = 0;
    INTEGER nrhs = 0;
    char xtype = char0;
    char trans = char0;
    INTEGER iter = 0;
    const INTEGER ntests = 1;
    arr_1d<ntests, REAL> result(fill0);
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
    //     .. Local Variables ..
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
    kase = 0;
    path[(1 - 1)] = "Double precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "GE";
    nrun = 0;
    nfail = 0;
    nerrs = 0;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    //
    cmn.infot = 0;
    //
    //     Do for each value of M in MVAL
    //
    for (im = 1; im <= nm; im = im + 1) {
        m = mval[im - 1];
        lda = max((INTEGER)1, m);
        //
        n = m;
        nimat = ntypes;
        if (m <= 0 || n <= 0) {
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
            //           Skip types 5, 6, or 7 if the matrix size is too small.
            //
            zerot = imat >= 5 && imat <= 7;
            if (zerot && n < imat - 4) {
                goto statement_100;
            }
            //
            //           Set up parameters with Rlatb4 and generate a test matrix
            //           with DLATMS.
            //
            Rlatb4(path, imat, m, n, type, kl, ku, anorm, mode, cndnum, dist);
            //
            srnamt = "DLATMS";
            dlatms(m, n, dist, iseed, type, rwork, mode, cndnum, anorm, kl, ku, "No packing", a, lda, work, info);
            //
            //           Check error code from DLATMS.
            //
            if (info != 0) {
                Alaerh(path, "DLATMS", info, 0, " ", m, n, -1, -1, -1, imat, nfail, nerrs, nout);
                goto statement_100;
            }
            //
            //           For types 5-7, zero one or more columns of the matrix to
            //           test that INFO is returned correctly.
            //
            if (zerot) {
                if (imat == 5) {
                    izero = 1;
                } else if (imat == 6) {
                    izero = min(m, n);
                } else {
                    izero = min(m, n) / 2 + 1;
                }
                ioff = (izero - 1) * lda;
                if (imat < 7) {
                    for (i = 1; i <= m; i = i + 1) {
                        a[(ioff + i) - 1] = zero;
                    }
                } else {
                    dlaset("Full", m, n - izero + 1, zero, zero, &a[(ioff + 1) - 1], lda);
                }
            } else {
                izero = 0;
            }
            //
            for (irhs = 1; irhs <= nns; irhs = irhs + 1) {
                nrhs = nsval[irhs - 1];
                xtype = "N";
                trans = "N";
                //
                srnamt = "Rlarhs";
                Rlarhs(path, xtype, " ", trans, n, n, kl, ku, nrhs, a, lda, x, lda, b, lda, iseed, info);
                //
                srnamt = "DSGESV";
                //
                kase++;
                //
                dlacpy("Full", m, n, a, lda, afac, lda);
                //
                dsgesv(n, nrhs, a, lda, iwork, b, lda, x, lda, work, swork, iter, info);
                //
                if (iter < 0) {
                    dlacpy("Full", m, n, afac, lda, a, lda);
                }
                //
                //              Check error code from DSGESV. This should be the same as
                //              the one of DGETRF.
                //
                if (info != izero) {
                    //
                    if (nfail == 0 && nerrs == 0) {
                        Alahd(nout, path);
                    }
                    nerrs++;
                    //
                    if (info != izero && izero != 0) {
                        write(nout, "(' *** ',a6,' returned with INFO =',i5,' instead of ',i5,/,"
                                    "' ==> M =',i5,', type ',i2)"),
                            "DSGESV", info, izero, m, imat;
                    } else {
                        write(nout, "(' *** Error code from ',a6,'=',i5,' for M=',i5,', type ',i2)"), "DSGESV", info, m, imat;
                    }
                }
                //
                //              Skip the remaining test if the matrix is singular.
                //
                if (info != 0) {
                    goto statement_100;
                }
                //
                //              Check the quality of the solution
                //
                dlacpy("Full", n, nrhs, b, lda, work, lda);
                //
                Rget08(trans, n, n, nrhs, a, lda, x, lda, work, lda, rwork, result[1 - 1]);
                //
                //              Check if the test passes the tesing.
                //              PrINTEGER information about the tests that did not
                //              pass the testing.
                //
                //              If iterative refinement has been used and claimed to
                //              be successful (ITER>0), we want
                //                NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS*SRQT(N)) < 1
                //
                //              If REAL precision has been used (ITER<0), we want
                //                NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS) < THRES
                //              (Cf. the linear solver testing routines)
                //
                if ((thresh <= 0.00f) || ((iter >= 0) && (n > 0) && (result[1 - 1] >= sqrt(n.real()))) || ((iter < 0) && (result[1 - 1] >= thresh))) {
                    //
                    if (nfail == 0 && nerrs == 0) {
                        write(nout, "(/,1x,a3,':  General dense matrices')"), "DGE";
                        write(nout, "(' Matrix types:')");
                        write(nout, "(4x,'1. Diagonal',24x,'7. Last n/2 columns zero',/,4x,"
                                    "'2. Upper triangular',16x,'8. Random, CNDNUM = sqrt(0.1/EPS)',"
                                    "/,4x,'3. Lower triangular',16x,'9. Random, CNDNUM = 0.1/EPS',/,"
                                    "4x,'4. Random, CNDNUM = 2',13x,'10. Scaled near underflow',/,"
                                    "4x,'5. First column zero',14x,'11. Scaled near overflow',/,4x,"
                                    "'6. Last column zero')");
                        write(nout, "(' Test ratios:')");
                        write(nout, "(3x,i2,': norm_1( B - A * X )  / ',"
                                    "'( norm_1(A) * norm_1(X) * EPS * SQRT(N) ) > 1 if ITERREF',/,"
                                    "4x,'or norm_1( B - A * X )  / ',"
                                    "'( norm_1(A) * norm_1(X) * EPS ) > THRES if DGETRF')"),
                            1;
                        write(nout, "(' Messages:')");
                    }
                    //
                    write(nout, "(' TRANS=''',a1,''', N =',i5,', NRHS=',i3,', type ',i2,', test(',"
                                "i2,') =',g12.5)"),
                        trans, n, nrhs, imat, 1, result(1);
                    nfail++;
                }
                nrun++;
            }
        statement_100:;
        }
    }
    //
    //     PrINTEGER a summary of the results.
    //
    if (nfail > 0) {
        write(nout, "(1x,a6,': ',i6,' out of ',i6,' tests failed to pass the threshold')"), "DSGESV", nfail, nrun;
    } else {
        write(nout, "(/,1x,'All tests for ',a6,' routines passed the threshold ( ',i6,"
                    "' tests run)')"),
            "DSGESV", nrun;
    }
    if (nerrs > 0) {
        write(nout, "(6x,i6,' error messages recorded')"), nerrs;
    }
    //
    //     SUBNAM, INFO, INFOE, M, IMAT
    //
    //     SUBNAM, INFO, M, IMAT
    //
    //     End of Rdrvab
    //
}
