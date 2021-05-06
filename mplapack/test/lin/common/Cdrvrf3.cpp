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

void Cdrvrf3(INTEGER const nout, INTEGER const nn, INTEGER *nval, REAL const thresh, COMPLEX *a, INTEGER const lda, COMPLEX *arf, COMPLEX *b1, COMPLEX *b2, REAL *d_work_Clange, COMPLEX *z_work_Cgeqrf, COMPLEX *tau) {
    nval([nn]);
    a([lda * star]);
    b1([lda * star]);
    b2([lda * star]);
    common_write write(cmn);
    // COMMON srnamc
    //
    // SAVE
    str_arr_ref<1> diags(sve.diags, [2]);
    str_arr_ref<1> forms(sve.forms, [2]);
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    str_arr_ref<1> sides(sve.sides, [2]);
    str_arr_ref<1> transs(sve.transs, [2]);
    str_arr_ref<1> uplos(sve.uplos, [2]);
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
            static const char *values[] = {"N", "C"};
            data_of_type_str(FEM_VALUES_AND_SIZE), forms;
        }
        {
            static const char *values[] = {"L", "R"};
            data_of_type_str(FEM_VALUES_AND_SIZE), sides;
        }
        {
            static const char *values[] = {"N", "C"};
            data_of_type_str(FEM_VALUES_AND_SIZE), transs;
        }
        {
            static const char *values[] = {"N", "U"};
            data_of_type_str(FEM_VALUES_AND_SIZE), diags;
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
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER info = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    REAL eps = Rlamch("Precision");
    //
    INTEGER iim = 0;
    INTEGER m = 0;
    INTEGER iin = 0;
    INTEGER n = 0;
    INTEGER iform = 0;
    char cform;
    INTEGER iuplo = 0;
    char uplo;
    INTEGER iside = 0;
    char side;
    INTEGER itrans = 0;
    char trans;
    INTEGER idiag = 0;
    char diag;
    INTEGER ialpha = 0;
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    COMPLEX alpha = 0.0;
    const COMPLEX one = COMPLEX(1.0, 0.0);
    INTEGER na = 0;
    INTEGER j = 0;
    const INTEGER ntests = 1;
    REAL result[ntests];
    for (iim = 1; iim <= nn; iim = iim + 1) {
        //
        m = nval[iim - 1];
        //
        for (iin = 1; iin <= nn; iin = iin + 1) {
            //
            n = nval[iin - 1];
            //
            for (iform = 1; iform <= 2; iform = iform + 1) {
                //
                cform = forms[iform - 1];
                //
                for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                    //
                    uplo = uplos[iuplo - 1];
                    //
                    for (iside = 1; iside <= 2; iside = iside + 1) {
                        //
                        side = sides[iside - 1];
                        //
                        for (itrans = 1; itrans <= 2; itrans = itrans + 1) {
                            //
                            trans = transs[itrans - 1];
                            //
                            for (idiag = 1; idiag <= 2; idiag = idiag + 1) {
                                //
                                diag = diags[idiag - 1];
                                //
                                for (ialpha = 1; ialpha <= 3; ialpha = ialpha + 1) {
                                    //
                                    if (ialpha == 1) {
                                        alpha = zero;
                                    } else if (ialpha == 2) {
                                        alpha = one;
                                    } else {
                                        alpha = Clarnd(4, iseed);
                                    }
                                    //
                                    //                             All the parameters are set:
                                    //                                CFORM, SIDE, UPLO, TRANS, DIAG, M, N,
                                    //                                and ALPHA
                                    //                             READY TO TEST!
                                    //
                                    nrun++;
                                    //
                                    if (iside == 1) {
                                        //
                                        //                                The case ISIDE.EQ.1 is when SIDE.EQ.'L'
                                        //                                -> A is M-by-M ( B is M-by-N )
                                        //
                                        na = m;
                                        //
                                    } else {
                                        //
                                        //                                The case ISIDE.EQ.2 is when SIDE.EQ.'R'
                                        //                                -> A is N-by-N ( B is M-by-N )
                                        //
                                        na = n;
                                        //
                                    }
                                    //
                                    //                             Generate A our NA--by--NA triangular
                                    //                             matrix.
                                    //                             Our test is based on forward error so we
                                    //                             do want A to be well conditioned! To get
                                    //                             a well-conditioned triangular matrix, we
                                    //                             take the R factor of the QR/LQ factorization
                                    //                             of a random matrix.
                                    //
                                    for (j = 1; j <= na; j = j + 1) {
                                        for (i = 1; i <= na; i = i + 1) {
                                            a[(i - 1) + (j - 1) * lda] = Clarnd(4, iseed);
                                        }
                                    }
                                    //
                                    if (iuplo == 1) {
                                        //
                                        //                                The case IUPLO.EQ.1 is when SIDE.EQ.'U'
                                        //                                -> QR factorization.
                                        //
                                        Cgeqrf(na, na, a, lda, tau, z_work_Cgeqrf, lda, info);
                                    } else {
                                        //
                                        //                                The case IUPLO.EQ.2 is when SIDE.EQ.'L'
                                        //                                -> QL factorization.
                                        //
                                        Cgelqf(na, na, a, lda, tau, z_work_Cgeqrf, lda, info);
                                    }
                                    //
                                    //                             After the QR factorization, the diagonal
                                    //                             of A is made of real numbers, we multiply
                                    //                             by a random complex number of absolute
                                    //                             value 1.0E+00.
                                    //
                                    for (j = 1; j <= na; j = j + 1) {
                                        a[(j - 1) + (j - 1) * lda] = a[(j - 1) + (j - 1) * lda] * Clarnd(5, iseed);
                                    }
                                    //
                                    //                             Store a copy of A in RFP format (in ARF).
                                    //
                                    Ctrttf(cform, uplo, na, a, lda, arf, info);
                                    //
                                    //                             Generate B1 our M--by--N right-hand side
                                    //                             and store a copy in B2.
                                    //
                                    for (j = 1; j <= n; j = j + 1) {
                                        for (i = 1; i <= m; i = i + 1) {
                                            b1[(i - 1) + (j - 1) * ldb1] = Clarnd(4, iseed);
                                            b2[(i - 1) + (j - 1) * ldb2] = b1[(i - 1) + (j - 1) * ldb1];
                                        }
                                    }
                                    //
                                    //                             Solve op( A ) X = B or X op( A ) = B
                                    //                             with Ctrsm
                                    //
                                    Ctrsm(side, uplo, trans, diag, m, n, alpha, a, lda, b1, lda);
                                    //
                                    //                             Solve op( A ) X = B or X op( A ) = B
                                    //                             with Ctfsm
                                    //
                                    Ctfsm(cform, side, uplo, trans, diag, m, n, alpha, arf, b2, lda);
                                    //
                                    //                             Check that the result agrees.
                                    //
                                    for (j = 1; j <= n; j = j + 1) {
                                        for (i = 1; i <= m; i = i + 1) {
                                            b1[(i - 1) + (j - 1) * ldb1] = b2[(i - 1) + (j - 1) * ldb2] - b1[(i - 1) + (j - 1) * ldb1];
                                        }
                                    }
                                    //
                                    result[1 - 1] = Clange("I", m, n, b1, lda, d_work_Clange);
                                    //
                                    result[1 - 1] = result[1 - 1] / sqrt(eps) / max({max(m, n), 1});
                                    //
                                    if (result[1 - 1] >= thresh) {
                                        if (nfail == 0) {
                                            write(nout, star);
                                            write(nout, "(1x,' *** Error(s) or Failure(s) while testing Ctfsm "
                                                        "        ***')");
                                        }
                                        write(nout, "(1x,'     Failure in ',a5,', CFORM=''',a1,''',',"
                                                    "' SIDE=''',a1,''',',' UPLO=''',a1,''',',' TRANS=''',a1,"
                                                    "''',',' DIAG=''',a1,''',',' M=',i3,', N =',i3,"
                                                    "', test=',g12.5)"),
                                            "Ctfsm", cform, side, uplo, trans, diag, m, n, result(1);
                                        nfail++;
                                    }
                                    //
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    if (nfail == 0) {
        write(nout, "(1x,'All tests for ',a5,' auxiliary routine passed the ',"
                    "'threshold ( ',i5,' tests run)')"),
            "Ctfsm", nrun;
    } else {
        write(nout, "(1x,a6,' auxiliary routine:',i5,' out of ',i5,"
                    "' tests failed to pass the threshold')"),
            "Ctfsm", nfail, nrun;
    }
    //
    //     End of Cdrvrf3
    //
}
