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

void Rdrvrf1(INTEGER const nout, INTEGER const nn, INTEGER *nval, REAL const thresh, REAL *a, INTEGER const lda, REAL *arf, REAL *work) {
    FEM_CMN_SVE(Rdrvrf1);
    common_write write(cmn);
    char &srnamt = cmn.srnamt;
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
            static const char *values[] = {"N", "T"};
            data_of_type_str(FEM_VALUES_AND_SIZE), forms;
        }
        {
            static const char *values[] = {"M", "1", "I", "F"};
            data_of_type_str(FEM_VALUES_AND_SIZE), norms;
        }
    }
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER info = 0;
    INTEGER i = 0;
    arr_1d<4, int> iseed;
    REAL eps = 0.0;
    REAL small = 0.0;
    const REAL one = 1.0;
    REAL large = 0.0;
    INTEGER iin = 0;
    INTEGER n = 0;
    INTEGER iit = 0;
    INTEGER j = 0;
    INTEGER iuplo = 0;
    char uplo[1];
    INTEGER iform = 0;
    char cform[1];
    INTEGER inorm = 0;
    char norm[1];
    REAL normarf = 0.0;
    REAL norma = 0.0;
    const INTEGER ntests = 1;
    arr_1d<ntests, REAL> result;
    static const char *format_9999 = "(1x,' *** Error(s) or Failure(s) while testing Rlansf         ***')";
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
    //     ..
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
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Initialize constants and the random number seed.
    //
    nrun = 0;
    nfail = 0;
    nerrs = 0;
    info = 0;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    //
    eps = Rlamch("Precision");
    small = Rlamch("Safe minimum");
    large = one / small;
    small = small * lda * lda;
    large = large / lda / lda;
    //
    for (iin = 1; iin <= nn; iin = iin + 1) {
        //
        n = nval[iin - 1];
        //
        for (iit = 1; iit <= 3; iit = iit + 1) {
            //           Nothing to do for N=0
            if (n == 0) {
                break;
            }
            //
            //           IIT = 1 : random matrix
            //           IIT = 2 : random matrix scaled near underflow
            //           IIT = 3 : random matrix scaled near overflow
            //
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= n; i = i + 1) {
                    a[(i - 1) + (j - 1) * lda] = dlarnd(2, iseed);
                }
            }
            //
            if (iit == 2) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= n; i = i + 1) {
                        a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * large;
                    }
                }
            }
            //
            if (iit == 3) {
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= n; i = i + 1) {
                        a[(i - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda] * small;
                    }
                }
            }
            //
            //           Do first for UPLO = 'U', then for UPLO = 'L'
            //
            for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                //
                uplo = uplos[iuplo - 1];
                //
                //              Do first for CFORM = 'N', then for CFORM = 'C'
                //
                for (iform = 1; iform <= 2; iform = iform + 1) {
                    //
                    cform = forms[iform - 1];
                    //
                    srnamt = "Rtrttf";
                    Rtrttf(cform, uplo, n, a, lda, arf, info);
                    //
                    //                 Check error code from Rtrttf
                    //
                    if (info != 0) {
                        if (nfail == 0 && nerrs == 0) {
                            write(nout, star);
                            write(nout, format_9999);
                        }
                        write(nout, "(1x,'     Error in ',a6,' with UPLO=''',a1,''', FORM=''',a1,"
                                    "''', N=',i5)"),
                            srnamt, uplo, cform, n;
                        nerrs++;
                        goto statement_100;
                    }
                    //
                    for (inorm = 1; inorm <= 4; inorm = inorm + 1) {
                        //
                        //                    Check all four norms: 'M', '1', 'I', 'F'
                        //
                        norm = norms[inorm - 1];
                        normarf = Rlansf(norm, cform, uplo, n, arf, work);
                        norma = Rlansy(norm, uplo, n, a, lda, work);
                        //
                        result[1 - 1] = (norma - normarf) / norma / eps;
                        nrun++;
                        //
                        if (result[1 - 1] >= thresh) {
                            if (nfail == 0 && nerrs == 0) {
                                write(nout, star);
                                write(nout, format_9999);
                            }
                            write(nout, "(1x,'     Failure in ',a6,' N=',i5,' TYPE=',i5,' UPLO=''',a1,"
                                        "''', FORM =''',a1,''', NORM=''',a1,''', test=',g12.5)"),
                                "Rlansf", n, iit, uplo, cform, norm, result(1);
                            nfail++;
                        }
                    }
                statement_100:;
                }
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    if (nfail == 0) {
        write(nout, "(1x,'All tests for ',a6,' auxiliary routine passed the ',"
                    "'threshold ( ',i5,' tests run)')"),
            "Rlansf", nrun;
    } else {
        write(nout, "(1x,a6,' auxiliary routine: ',i5,' out of ',i5,"
                    "' tests failed to pass the threshold')"),
            "Rlansf", nfail, nrun;
    }
    if (nerrs != 0) {
        write(nout, "(26x,i5,' error message recorded (',a6,')')"), nerrs, "Rlansf";
    }
    //
    //     End of Rdrvrf1
    //
}
