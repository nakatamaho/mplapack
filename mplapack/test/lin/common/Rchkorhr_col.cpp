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

void Rchkorhr_col(REAL const thresh, bool const tsterr, INTEGER const nm, INTEGER *mval, INTEGER const nn, INTEGER *nval, INTEGER const nnb, INTEGER *nbval, INTEGER const nout) {
    common_write write(cmn);
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
    //
    //     .. Local Arrays ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Initialize constants
    //
    char path[1] = "D";
    path[(2 - 1) + (3 - 1) * ldpath] = "HH";
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    //
    //     Test the error exits
    //
    if (tsterr) {
        Rerrorhr_col(path, nout);
    }
    cmn.infot = 0;
    //
    //     Do for each value of M in MVAL.
    //
    INTEGER i = 0;
    INTEGER m = 0;
    INTEGER j = 0;
    INTEGER n = 0;
    INTEGER imb1 = 0;
    INTEGER mb1 = 0;
    INTEGER inb1 = 0;
    INTEGER nb1 = 0;
    INTEGER inb2 = 0;
    INTEGER nb2 = 0;
    const INTEGER ntests = 6;
    arr_1d<ntests, REAL> result;
    INTEGER t = 0;
    for (i = 1; i <= nm; i = i + 1) {
        m = mval[i - 1];
        //
        //        Do for each value of N in NVAL.
        //
        for (j = 1; j <= nn; j = j + 1) {
            n = nval[j - 1];
            //
            //           Only for M >= N
            //
            if (min(m, n) > 0 && m >= n) {
                //
                //              Do for each possible value of MB1
                //
                for (imb1 = 1; imb1 <= nnb; imb1 = imb1 + 1) {
                    mb1 = nbval[imb1 - 1];
                    //
                    //                 Only for MB1 > N
                    //
                    if (mb1 > n) {
                        //
                        //                    Do for each possible value of NB1
                        //
                        for (inb1 = 1; inb1 <= nnb; inb1 = inb1 + 1) {
                            nb1 = nbval[inb1 - 1];
                            //
                            //                       Do for each possible value of NB2
                            //
                            for (inb2 = 1; inb2 <= nnb; inb2 = inb2 + 1) {
                                nb2 = nbval[inb2 - 1];
                                //
                                if (nb1 > 0 && nb2 > 0) {
                                    //
                                    //                             Test Rorhr_col
                                    //
                                    Rorhr_col01(m, n, mb1, nb1, nb2, result);
                                    //
                                    //                             Print information about the tests that did
                                    //                             not pass the threshold.
                                    //
                                    for (t = 1; t <= ntests; t = t + 1) {
                                        if (result[t - 1] >= thresh) {
                                            if (nfail == 0 && nerrs == 0) {
                                                Alahd(nout, path);
                                            }
                                            write(nout, "('Rorgtsqr and Rorhr_col: M=',i5,', N=',i5,', MB1=',"
                                                        "i5,', NB1=',i5,', NB2=',i5,' test(',i2,')=',g12.5)"),
                                                m, n, mb1, nb1, nb2, t, result(t);
                                            nfail++;
                                        }
                                    }
                                    nrun += ntests;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //
    //     Do for each value of M in MVAL.
    //
    for (i = 1; i <= nm; i = i + 1) {
        m = mval[i - 1];
        //
        //        Do for each value of N in NVAL.
        //
        for (j = 1; j <= nn; j = j + 1) {
            n = nval[j - 1];
            //
            //           Only for M >= N
            //
            if (min(m, n) > 0 && m >= n) {
                //
                //              Do for each possible value of MB1
                //
                for (imb1 = 1; imb1 <= nnb; imb1 = imb1 + 1) {
                    mb1 = nbval[imb1 - 1];
                    //
                    //                 Only for MB1 > N
                    //
                    if (mb1 > n) {
                        //
                        //                    Do for each possible value of NB1
                        //
                        for (inb1 = 1; inb1 <= nnb; inb1 = inb1 + 1) {
                            nb1 = nbval[inb1 - 1];
                            //
                            //                       Do for each possible value of NB2
                            //
                            for (inb2 = 1; inb2 <= nnb; inb2 = inb2 + 1) {
                                nb2 = nbval[inb2 - 1];
                                //
                                if (nb1 > 0 && nb2 > 0) {
                                    //
                                    //                             Test Rorhr_col
                                    //
                                    Rorhr_col02(m, n, mb1, nb1, nb2, result);
                                    //
                                    //                             Print information about the tests that did
                                    //                             not pass the threshold.
                                    //
                                    for (t = 1; t <= ntests; t = t + 1) {
                                        if (result[t - 1] >= thresh) {
                                            if (nfail == 0 && nerrs == 0) {
                                                Alahd(nout, path);
                                            }
                                            write(nout, "('Rorgtsqr_row and Rorhr_col: M=',i5,', N=',i5,"
                                                        "', MB1=',i5,', NB1=',i5,', NB2=',i5,' test(',i2,')=',"
                                                        "g12.5)"),
                                                m, n, mb1, nb1, nb2, t, result(t);
                                            nfail++;
                                        }
                                    }
                                    nrun += ntests;
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
    Alasum(path, nout, nfail, nrun, nerrs);
    //
    //     End of Rchkorhr_col
    //
}
