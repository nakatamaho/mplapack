/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ALARQG.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_eig.h>

void Alarqg(fem::str_cref path, INTEGER const nmats, bool *dotype, INTEGER const ntypes, INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    static fem::str<10> intstr = "0123456789";
    INTEGER i = 0;
    bool firstt = false;
    fem::str<80> line;
    INTEGER lenp = 0;
    INTEGER j = 0;
    INTEGER nreq[100];
    INTEGER i1 = 0;
    fem::str<1> c1;
    INTEGER k = 0;
    INTEGER ic = 0;
    INTEGER nt = 0;
    static const char *format_9999 = "(' *** Invalid type request for ',a3,', type  ',i4,"
                                     "': must satisfy  1 <= type <= ',i2)";
    static const char *format_9998 = "(/,' *** End of file reached when trying to read matrix ','types for ',"
                                     "a3,/,' *** Check that you are requesting the',"
                                     "' right number of types for each path',/)";
    static const char *format_9997 = "(' *** Warning:  duplicate request of matrix type ',i2,' for ',a3)";
    static const char *format_9996 = "(/,/,' *** Invalid integer value in column ',i2,' of input',' line:',/,"
                                     "a79)";
    static const char *format_9995 = "(/,/,' *** Not enough matrix types on input line',/,a79)";
    static const char *format_9994 = "(' ==> Specify ',i4,' matrix types on this line or ',"
                                     "'adjust NTYPES on previous line')";
    //
    if (nmats >= ntypes) {
        //
        // Test everything if NMATS >= NTYPES.
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
        // Read a line of matrix types if 0 < NMATS < NTYPES.
        //
        if (nmats > 0) {
            try {
                read(nin, "(a80)"), line;
            } catch (fem::read_end const &) {
                goto statement_90;
            }
            lenp = fem::len(line);
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
                        write(nout, format_9995), line;
                        write(nout, format_9994), nmats;
                        goto statement_80;
                    }
                }
                if (line(i, i) != " " && line(i, i) != ",") {
                    i1 = i;
                    c1 = line(i1, i1);
                    //
                    // Check that a valid integer was read
                    //
                    for (k = 1; k <= 10; k = k + 1) {
                        if (c1 == intstr(k, k)) {
                            ic = k - 1;
                            goto statement_50;
                        }
                    }
                    write(nout, format_9996), i, line;
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
                    write(nout, format_9997), nt, path;
                }
                dotype[nt - 1] = true;
            } else {
                write(nout, format_9999), path, nt, ntypes;
            }
        }
    statement_80:;
    }
    return;
//
statement_90:
    write(nout, format_9998), path;
    write(nout, star);
    FEM_STOP(0);
    //
    // End of Alarqg
    //
}
