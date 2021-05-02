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

void Rdrvrf2(INTEGER const nout, INTEGER const nn, INTEGER *nval, REAL *a, INTEGER const lda, REAL *arf, REAL *ap, REAL *asav) {
    FEM_CMN_SVE(Rdrvrf2);
    common_write write(cmn);
    // COMMON srnamc
    char &srnamt = cmn.srnamt;
    //
    // SAVE
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
    }
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
    INTEGER nrun = 0;
    INTEGER nerrs = 0;
    INTEGER info = 0;
    INTEGER i = 0;
    arr_1d<4, int> iseed;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    //
    INTEGER iin = 0;
    INTEGER n = 0;
    INTEGER iuplo = 0;
    char uplo[1];
    bool lower = false;
    INTEGER iform = 0;
    char cform[1];
    INTEGER j = 0;
    bool ok1 = false;
    bool ok2 = false;
    for (iin = 1; iin <= nn; iin = iin + 1) {
        //
        n = nval[iin - 1];
        //
        //        Do first for UPLO = 'U', then for UPLO = 'L'
        //
        for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
            //
            uplo = uplos[iuplo - 1];
            lower = true;
            if (iuplo == 1) {
                lower = false;
            }
            //
            //           Do first for CFORM = 'N', then for CFORM = 'T'
            //
            for (iform = 1; iform <= 2; iform = iform + 1) {
                //
                cform = forms[iform - 1];
                //
                nrun++;
                //
                for (j = 1; j <= n; j = j + 1) {
                    for (i = 1; i <= n; i = i + 1) {
                        a[(i - 1) + (j - 1) * lda] = dlarnd(2, iseed);
                    }
                }
                //
                srnamt = "Rtrttf";
                Rtrttf(cform, uplo, n, a, lda, arf, info);
                //
                srnamt = "Rtfttp";
                Rtfttp(cform, uplo, n, arf, ap, info);
                //
                srnamt = "Rtpttr";
                Rtpttr(uplo, n, ap, asav, lda, info);
                //
                ok1 = true;
                if (lower) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            if (a[(i - 1) + (j - 1) * lda] != asav[(i - 1) + (j - 1) * ldasav]) {
                                ok1 = false;
                            }
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            if (a[(i - 1) + (j - 1) * lda] != asav[(i - 1) + (j - 1) * ldasav]) {
                                ok1 = false;
                            }
                        }
                    }
                }
                //
                nrun++;
                //
                srnamt = "Rtrttp";
                Rtrttp(uplo, n, a, lda, ap, info);
                //
                srnamt = "Rtpttf";
                Rtpttf(cform, uplo, n, ap, arf, info);
                //
                srnamt = "Rtfttr";
                Rtfttr(cform, uplo, n, arf, asav, lda, info);
                //
                ok2 = true;
                if (lower) {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = j; i <= n; i = i + 1) {
                            if (a[(i - 1) + (j - 1) * lda] != asav[(i - 1) + (j - 1) * ldasav]) {
                                ok2 = false;
                            }
                        }
                    }
                } else {
                    for (j = 1; j <= n; j = j + 1) {
                        for (i = 1; i <= j; i = i + 1) {
                            if (a[(i - 1) + (j - 1) * lda] != asav[(i - 1) + (j - 1) * ldasav]) {
                                ok2 = false;
                            }
                        }
                    }
                }
                //
                if ((!ok1) || (!ok2)) {
                    if (nerrs == 0) {
                        write(nout, star);
                        write(nout, "(1x,' *** Error(s) while testing the RFP conversion',"
                                    "' routines ***')");
                    }
                    write(nout, "(1x,'     Error in RFP,conversion routines N=',i5,' UPLO=''',a1,"
                                "''', FORM =''',a1,'''')"),
                        n, uplo, cform;
                    nerrs++;
                }
                //
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    if (nerrs == 0) {
        write(nout, "(1x,'All tests for the RFP conversion routines passed ( ',i5,"
                    "' tests run)')"),
            nrun;
    } else {
        write(nout, "(1x,'RFP conversion routines: ',i5,' out of ',i5,"
                    "' error message recorded')"),
            nerrs, nrun;
    }
    //
    //     End of Rdrvrf2
    //
}
