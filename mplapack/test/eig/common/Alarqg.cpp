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
#include <mplapack_eig.h>

#include <mplapack_debug.h>

void Alarqg(const char *path, INTEGER const nmats, bool *dotype, INTEGER const ntypes, INTEGER const nin, INTEGER const nout) {
    FEM_CMN_SVE(Alarqg);
    common_read read(cmn);
    common_write write(cmn);
    char &intstr = sve.intstr;
    if (is_called_first_time) {
        intstr = "0123456789";
    }
    INTEGER i = 0;
    bool firstt = false;
    char line[80];
    INTEGER lenp = 0;
    INTEGER j = 0;
    INTEGER nreq[100];
    INTEGER i1 = 0;
    char c1;
    INTEGER k = 0;
    INTEGER ic = 0;
    INTEGER nt = 0;
    static const char *format_9994 = "(' ==> Specify ',i4,' matrix types on this line or ',"
                                     "'adjust NTYPES on previous line')";
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
    // ======================================================================
    //
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    if (nmats >= ntypes) {
        //
        //        Test everything if NMATS >= NTYPES.
        //
        for (i = 1; i <= ntypes; i = i + 1) {
            dotype[i - 1] = true;
        }
    } else {
        for (i = 1; i <= ntypes; i = i + 1) {
            dotype[i - 1] = false;
        }
        firstt = true;
        //
        //        Read a line of matrix types if 0 < NMATS < NTYPES.
        //
        if (nmats > 0) {
            try {
                read(nin, "(a80)"), line;
            } catch (read_end const ) {
                goto statement_90;
            }
            lenp = len[line - 1];
            i = 0;
            for (j = 1; j <= nmats; j = j + 1) {
                nreq[j - 1] = 0;
                i1 = 0;
            statement_30:
                i++;
                if (i > lenp) {
                    if (j == nmats && i1 > 0) {
                        goto statement_60;
                    } else {
                        write(nout, "(/,/,' *** Not enough matrix types on input line',/,a79)"), line;
                        write(nout, format_9994), nmats;
                        goto statement_80;
                    }
                }
                if (line[(i - 1) + (i - 1) * ldline] != " " && line[(i - 1) + (i - 1) * ldline] != ",") {
                    i1 = i;
                    c1 = line[(i1 - 1) + (i1 - 1) * ldline];
                    //
                    //              Check that a valid integer was read
                    //
                    for (k = 1; k <= 10; k = k + 1) {
                        if (c1 == intstr[(k - 1) + (k - 1) * ldintstr]) {
                            ic = k - 1;
                            goto statement_50;
                        }
                    }
                    write(nout, "(/,/,' *** Invalid integer value in column ',i2,' of input',"
                                "' line:',/,a79)"),
                        i, line;
                    write(nout, format_9994), nmats;
                    goto statement_80;
                statement_50:
                    nreq[j - 1] = 10 * nreq[j - 1] + ic;
                    goto statement_30;
                } else if (i1 > 0) {
                    goto statement_60;
                } else {
                    goto statement_30;
                }
            statement_60:;
            }
        }
        for (i = 1; i <= nmats; i = i + 1) {
            nt = nreq[i - 1];
            if (nt > 0 && nt <= ntypes) {
                if (dotype[nt - 1]) {
                    if (firstt) {
                        write(nout, star);
                    }
                    firstt = false;
                    write(nout, "(' *** Warning:  duplicate request of matrix type ',i2,' for ',"
                                "a3)"),
                        nt, path;
                }
                dotype[nt - 1] = true;
            } else {
                write(nout, "(' *** Invalid type request for ',a3,', type  ',i4,"
                            "': must satisfy  1 <= type <= ',i2)"),
                    path, nt, ntypes;
            }
        }
    statement_80:;
    }
    return;
//
statement_90:
    write(nout, "(/,' *** End of file reached when trying to read matrix ','types for ',"
                "a3,/,' *** Check that you are requesting the',"
                "' right number of types for each path',/)"),
        path;
    write(nout, star);
    FEM_STOP(0);
    //
    //     End of Alarqg
    //
}
