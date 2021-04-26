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

void Rdrvrf4(common &cmn, INTEGER const nout, INTEGER const nn, INTEGER *nval, REAL const thresh, REAL *c1, REAL *c2, INTEGER const ldc, REAL *crf, REAL *a, INTEGER const lda, REAL *d_work_dlange) {
    FEM_CMN_SVE(Rdrvrf4);
    common_write write(cmn);
    // COMMON srnamc
    str<32> &srnamt = cmn.srnamt;
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
        {
            static const char *values[] = {"N", "T"};
            data_of_type_str(FEM_VALUES_AND_SIZE), transs;
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
    arr_1d<4, int> iseed(fill0);
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    REAL eps = Rlamch("Precision");
    //
    INTEGER iin = 0;
    INTEGER n = 0;
    INTEGER iik = 0;
    INTEGER k = 0;
    INTEGER iform = 0;
    char cform = char0;
    INTEGER iuplo = 0;
    char uplo = char0;
    INTEGER itrans = 0;
    char trans = char0;
    INTEGER ialpha = 0;
    const REAL zero = 0.0;
    REAL alpha = 0.0;
    REAL beta = 0.0;
    const REAL one = 1.0;
    INTEGER j = 0;
    REAL norma = 0.0;
    REAL normc = 0.0;
    const INTEGER ntests = 1;
    arr_1d<ntests, REAL> result(fill0);
    for (iin = 1; iin <= nn; iin = iin + 1) {
        //
        n = nval[iin - 1];
        //
        for (iik = 1; iik <= nn; iik = iik + 1) {
            //
            k = nval[iin - 1];
            //
            for (iform = 1; iform <= 2; iform = iform + 1) {
                //
                cform = forms[iform - 1];
                //
                for (iuplo = 1; iuplo <= 2; iuplo = iuplo + 1) {
                    //
                    uplo = uplos[iuplo - 1];
                    //
                    for (itrans = 1; itrans <= 2; itrans = itrans + 1) {
                        //
                        trans = transs[itrans - 1];
                        //
                        for (ialpha = 1; ialpha <= 4; ialpha = ialpha + 1) {
                            //
                            if (ialpha == 1) {
                                alpha = zero;
                                beta = zero;
                            } else if (ialpha == 2) {
                                alpha = one;
                                beta = zero;
                            } else if (ialpha == 3) {
                                alpha = zero;
                                beta = one;
                            } else {
                                alpha = dlarnd(2, iseed);
                                beta = dlarnd(2, iseed);
                            }
                            //
                            //                       All the parameters are set:
                            //                          CFORM, UPLO, TRANS, M, N,
                            //                          ALPHA, and BETA
                            //                       READY TO TEST!
                            //
                            nrun++;
                            //
                            if (itrans == 1) {
                                //
                                //                          In this case we are NOTRANS, so A is N-by-K
                                //
                                for (j = 1; j <= k; j = j + 1) {
                                    for (i = 1; i <= n; i = i + 1) {
                                        a[(i - 1) + (j - 1) * lda] = dlarnd(2, iseed);
                                    }
                                }
                                //
                                norma = dlange("I", n, k, a, lda, d_work_dlange);
                                //
                            } else {
                                //
                                //                          In this case we are TRANS, so A is K-by-N
                                //
                                for (j = 1; j <= n; j = j + 1) {
                                    for (i = 1; i <= k; i = i + 1) {
                                        a[(i - 1) + (j - 1) * lda] = dlarnd(2, iseed);
                                    }
                                }
                                //
                                norma = dlange("I", k, n, a, lda, d_work_dlange);
                                //
                            }
                            //
                            //                       Generate C1 our N--by--N symmetric matrix.
                            //                       Make sure C2 has the same upper/lower part,
                            //                       (the one that we do not touch), so
                            //                       copy the initial C1 in C2 in it.
                            //
                            for (j = 1; j <= n; j = j + 1) {
                                for (i = 1; i <= n; i = i + 1) {
                                    c1[(i - 1) + (j - 1) * ldc1] = dlarnd(2, iseed);
                                    c2[(i - 1) + (j - 1) * ldc2] = c1[(i - 1) + (j - 1) * ldc1];
                                }
                            }
                            //
                            //                       (See comment later on for why we use DLANGE and
                            //                       not DLANSY for C1.)
                            //
                            normc = dlange("I", n, n, c1, ldc, d_work_dlange);
                            //
                            srnamt = "DTRTTF";
                            dtrttf(cform, uplo, n, c1, ldc, crf, info);
                            //
                            //                       call Rsyrk the BLAS routine -> gives C1
                            //
                            srnamt = "Rsyrk ";
                            Rsyrk(uplo, trans, n, k, alpha, a, lda, beta, c1, ldc);
                            //
                            //                       call dsfrk the RFP routine -> gives CRF
                            //
                            srnamt = "DSFRK ";
                            dsfrk(cform, uplo, trans, n, k, alpha, a, lda, beta, crf);
                            //
                            //                       convert CRF in full format -> gives C2
                            //
                            srnamt = "DTFTTR";
                            dtfttr(cform, uplo, n, crf, c2, ldc, info);
                            //
                            //                       compare C1 and C2
                            //
                            for (j = 1; j <= n; j = j + 1) {
                                for (i = 1; i <= n; i = i + 1) {
                                    c1[(i - 1) + (j - 1) * ldc1] = c1[(i - 1) + (j - 1) * ldc1] - c2[(i - 1) + (j - 1) * ldc2];
                                }
                            }
                            //
                            //                       Yes, C1 is symmetric so we could call DLANSY,
                            //                       but we want to check the upper part that is
                            //                       supposed to be unchanged and the diagonal that
                            //                       is supposed to be real -> DLANGE
                            //
                            result[1 - 1] = dlange("I", n, n, c1, ldc, d_work_dlange);
                            result[1 - 1] = result[1 - 1] / max(abs(alpha) * norma + abs(beta), one) / max(n, 1) / eps;
                            //
                            if (result[1 - 1] >= thresh) {
                                if (nfail == 0) {
                                    write(nout, star);
                                    write(nout, "(1x,' *** Error(s) or Failure(s) while testing DSFRK     "
                                                "    ***')");
                                }
                                write(nout, "(1x,'     Failure in ',a5,', CFORM=''',a1,''',',' UPLO=''',"
                                            "a1,''',',' TRANS=''',a1,''',',' N=',i3,', K =',i3,"
                                            "', test=',g12.5)"),
                                    "DSFRK", cform, uplo, trans, n, k, result(1);
                                nfail++;
                            }
                            //
                        }
                    }
                }
            }
        }
    }
    //
    //     PrINTEGER a summary of the results.
    //
    if (nfail == 0) {
        write(nout, "(1x,'All tests for ',a5,' auxiliary routine passed the ',"
                    "'threshold ( ',i5,' tests run)')"),
            "DSFRK", nrun;
    } else {
        write(nout, "(1x,a6,' auxiliary routine: ',i5,' out of ',i5,"
                    "' tests failed to pass the threshold')"),
            "DSFRK", nfail, nrun;
    }
    //
    //     End of Rdrvrf4
    //
}
