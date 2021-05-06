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

void Rchktsqr(REAL const thresh, bool const tsterr, INTEGER const nm, INTEGER *mval, INTEGER const nn, INTEGER *nval, INTEGER const nnb, INTEGER *nbval, INTEGER const nout) {
    common cmn;
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
    char path[3];
    path[0] = 'D';
    path[1] = 'T';
    path[2] = 'S';
    char buf[1024];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    //
    //     Test the error exits
    //
    if (tsterr) {
        Rerrtsqr(path, nout);
    }
    //
    //     Do for each value of M in MVAL.
    //
    INTEGER i = 0;
    INTEGER m = 0;
    INTEGER j = 0;
    INTEGER n = 0;
    INTEGER inb = 0;
    INTEGER mb = 0;
    INTEGER imb = 0;
    INTEGER nb = 0;
    const INTEGER ntests = 6;
    REAL result[ntests];
    INTEGER t = 0;
    for (i = 1; i <= nm; i = i + 1) {
        m = mval[i - 1];
        //
        //        Do for each value of N in NVAL.
        //
        for (j = 1; j <= nn; j = j + 1) {
            n = nval[j - 1];
            if (min(m, n) != 0) {
                for (inb = 1; inb <= nnb; inb = inb + 1) {
                    mb = nbval[inb - 1];
                    for (imb = 1; imb <= nnb; imb = imb + 1) {
                        nb = nbval[imb - 1];
                        //
                        //                 Test Rgeqr and Rgemqr
                        //
                        Rtsqr01("TS", m, n, mb, nb, result);
                        //
                        //                 Print information about the tests that did not
                        //                 pass the threshold.
                        //
                        for (t = 1; t <= ntests; t = t + 1) {
                            if (result[t - 1] >= thresh) {
                                if (nfail == 0 && nerrs == 0) {
                                    Alahd(nout, path);
                                }
                                sprintnum_short(buf, result[t - 1]);
                                write(nout, "('TS: M=',i5,', N=',i5,', MB=',i5,', NB=',i5,' test(',i2,"
                                            "')=',a)"),
                                    m, n, mb, nb, t, buf;
                                nfail++;
                            }
                        }
                        nrun += ntests;
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
            if (min(m, n) != 0) {
                for (inb = 1; inb <= nnb; inb = inb + 1) {
                    mb = nbval[inb - 1];
                    for (imb = 1; imb <= nnb; imb = imb + 1) {
                        nb = nbval[imb - 1];
                        //
                        //                 Test Rgeqr and Rgemqr
                        //
                        Rtsqr01("SW", m, n, mb, nb, result);
                        //
                        //                 Print information about the tests that did not
                        //                 pass the threshold.
                        //
                        for (t = 1; t <= ntests; t = t + 1) {
                            if (result[t - 1] >= thresh) {
                                if (nfail == 0 && nerrs == 0) {
                                    Alahd(nout, path);
                                }
                                sprintnum_short(buf, result[t - 1]);
                                write(nout, "('SW: M=',i5,', N=',i5,', MB=',i5,', NB=',i5,' test(',i2,"
                                            "')=',a)"),
                                    m, n, mb, nb, t, buf;
                                nfail++;
                            }
                        }
                        nrun += ntests;
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
    //     End of Rchkqrt
    //
}
