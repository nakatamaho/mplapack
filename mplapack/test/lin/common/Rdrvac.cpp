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

void Rdrvac(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nns, INTEGER *nsval, REAL const thresh, INTEGER const /* nmax */, REAL *a, REAL *afac, REAL *b, REAL *x, REAL *work, REAL *rwork, arr_cref<float> swork, INTEGER const nout) {
    FEM_CMN_SVE(Rdrvac);
    common_write write(cmn);
    //
    INTEGER *iseedy(sve.iseedy, [4]);
    str_arr_ref<1> uplos(sve.uplos, [2]);
    if (is_called_first_time) {
        {
            static const INTEGER values[] = {1988, 1989, 1990, 1991};
            data_of_type<int>(FEM_VALUES_AND_SIZE), iseedy;
        }
        {
            static const char *values[] = {"U", "L"};
            data_of_type_str(FEM_VALUES_AND_SIZE), uplos;
        }
    }
    INTEGER kase = 0;
    char path[3];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    INTEGER im = 0;
    INTEGER n = 0;
    INTEGER lda = 0;
    const INTEGER ntypes = 9;
    INTEGER nimat = 0;
    INTEGER imat = 0;
    bool zerot = false;
    INTEGER iuplo = 0;
    char uplo;
    char type;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist;
    INTEGER info = 0;
    INTEGER izero = 0;
    INTEGER ioff = 0;
    const REAL zero = 0.0;
    INTEGER irhs = 0;
    INTEGER nrhs = 0;
    char xtype;
    INTEGER iter = 0;
    const INTEGER ntests = 1;
    REAL result[ntests];
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
    kase = 0;
    path[(1 - 1)] = "Double precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "PO";
    nrun = 0;
    nfail = 0;
    nerrs = 0;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    //
    cmn.infot = 0;
    //
    //     Do for each value of N in MVAL
    //
    for (im = 1; im <= nm; im = im + 1) {
        n = mval[im - 1];
        lda = max(n, 1);
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
                goto statement_110;
            }
            //
            //           Skip types 3, 4, or 5 if the matrix size is too small.
            //
            zerot = imat >= 3 && imat <= 5;
            if (zerot && n < imat - 2) {
                goto statement_110;
            }
            //
            //           Do first for UPLO = 'U', then for UPLO = 'L'
            //
            for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                uplo = uplos[iuplo - 1];
                //
                //              Set up parameters with Rlatb4 and generate a test matrix
                //              with Rlatms.
                //
                Rlatb4(path, imat, n, n, type, kl, ku, anorm, mode, cndnum, dist);
                //
                Rlatms(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, kl, ku, uplo, a, lda, work, info);
                //
                //              Check error code from Rlatms.
                //
                if (info != 0) {
                    Alaerh(path, "Rlatms", info, 0, uplo, n, n, -1, -1, -1, imat, nfail, nerrs, nout);
                    goto statement_100;
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
                    ioff = (izero - 1) * lda;
                    //
                    //                 Set row and column IZERO of A to 0.
                    //
                    if (iuplo == 1) {
                        for (i = 1; i <= izero - 1; i = i + 1) {
                            a[(ioff + i) - 1] = zero;
                        }
                        ioff += izero;
                        for (i = izero; i <= n; i = i + 1) {
                            a[ioff - 1] = zero;
                            ioff += lda;
                        }
                    } else {
                        ioff = izero;
                        for (i = 1; i <= izero - 1; i = i + 1) {
                            a[ioff - 1] = zero;
                            ioff += lda;
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
                for (irhs = 1; irhs <= nns; irhs = irhs + 1) {
                    nrhs = nsval[irhs - 1];
                    xtype = "N";
                    //
                    //                 Form an exact solution and set the right hand side.
                    //
                    Rlarhs(path, xtype, uplo, " ", n, n, kl, ku, nrhs, a, lda, x, lda, b, lda, iseed, info);
                    //
                    //                 Compute the L*L' or U'*U factorization of the
                    //                 matrix and solve the system.
                    //
                    kase++;
                    //
                    Rlacpy("All", n, n, a, lda, afac, lda);
                    //
                    Rsposv(uplo, n, nrhs, afac, lda, b, lda, x, lda, work, swork, iter, info);
                    //
                    if (iter < 0) {
                        Rlacpy("All", n, n, a, lda, afac, lda);
                    }
                    //
                    //                 Check error code from Rsposv .
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
                                        "' ==> N =',i5,', type ',i2)"),
                                "Rsposv", info, izero, n, imat;
                        } else {
                            write(nout, "(' *** Error code from ',a6,'=',i5,' for M=',i5,', type ',"
                                        "i2)"),
                                "Rsposv", info, n, imat;
                        }
                    }
                    //
                    //                 Skip the remaining test if the matrix is singular.
                    //
                    if (info != 0) {
                        goto statement_110;
                    }
                    //
                    //                 Check the quality of the solution
                    //
                    Rlacpy("All", n, nrhs, b, lda, work, lda);
                    //
                    Rpot06(uplo, n, nrhs, a, lda, x, lda, work, lda, rwork, result[1 - 1]);
                    //
                    //                 Check if the test passes the tesing.
                    //                 Print information about the tests that did not
                    //                 pass the testing.
                    //
                    //                 If iterative refinement has been used and claimed to
                    //                 be successful (ITER>0), we want
                    //                 NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS*SRQT(N)) < 1
                    //
                    //                 If REAL precision has been used (ITER<0), we want
                    //                 NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS) < THRES
                    //                 (Cf. the linear solver testing routines)
                    //
                    if ((thresh <= 0.00f) || ((iter >= 0) && (n > 0) && (result[1 - 1] >= sqrt(n.real()))) || ((iter < 0) && (result[1 - 1] >= thresh))) {
                        //
                        if (nfail == 0 && nerrs == 0) {
                            write(nout, "(/,1x,a3,':  positive definite dense matrices')"), "DPO";
                            write(nout, "(' Matrix types:')");
                            write(nout, "(4x,'1. Diagonal',24x,'7. Last n/2 columns zero',/,4x,"
                                        "'2. Upper triangular',16x,"
                                        "'8. Random, CNDNUM = sqrt(0.1/EPS)',/,4x,"
                                        "'3. Lower triangular',16x,'9. Random, CNDNUM = 0.1/EPS',/,4x,"
                                        "'4. Random, CNDNUM = 2',13x,'10. Scaled near underflow',/,4x,"
                                        "'5. First column zero',14x,'11. Scaled near overflow',/,4x,"
                                        "'6. Last column zero')");
                            write(nout, "(' Test ratios:')");
                            write(nout, "(3x,i2,': norm_1( B - A * X )  / ',"
                                        "'( norm_1(A) * norm_1(X) * EPS * SQRT(N) ) > 1 if ITERREF',/,"
                                        "4x,'or norm_1( B - A * X )  / ',"
                                        "'( norm_1(A) * norm_1(X) * EPS ) > THRES if Rpotrf')"),
                                1;
                            write(nout, "(' Messages:')");
                        }
                        //
                        write(nout, "(' UPLO=''',a1,''', N =',i5,', NRHS=',i3,', type ',i2,"
                                    "', test(',i2,') =',g12.5)"),
                            uplo, n, nrhs, imat, 1, result(1);
                        //
                        nfail++;
                        //
                    }
                    //
                    nrun++;
                    //
                }
            statement_100:;
            }
        statement_110:;
        }
    }
    //
    //     Print a summary of the results.
    //
    if (nfail > 0) {
        write(nout, "(1x,a6,': ',i6,' out of ',i6,' tests failed to pass the threshold')"), "Rsposv", nfail, nrun;
    } else {
        write(nout, "(/,1x,'All tests for ',a6,' routines passed the threshold ( ',i6,"
                    "' tests run)')"),
            "Rsposv", nrun;
    }
    if (nerrs > 0) {
        write(nout, "(6x,i6,' error messages recorded')"), nerrs;
    }
    //
    //     SUBNAM, INFO, INFOE, N, IMAT
    //
    //     SUBNAM, INFO, N, IMAT
    //
    //     End of Rdrvac
    //
}
