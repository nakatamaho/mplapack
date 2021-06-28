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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

void Alareq(const char *path, INTEGER const nmats, bool *dotype, INTEGER const ntypes, INTEGER const nin, INTEGER const nout) {
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    INTEGER i = 0;
    bool firstt = false;
    INTEGER lenp = 0;
    INTEGER j = 0;
    INTEGER nreq[100];
    INTEGER i1 = 0;
    char c1;
    INTEGER k = 0;
    INTEGER ic = 0;
    INTEGER nt = 0;
    string str;
    istringstream iss;
    char line[1024];
    double dtmp;
    int itmp;

    static const char *format_9994 = "(' ==> Specify ',i4,' matrix types on this line or ',"
                                     "'adjust NTYPES on previous line')";
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
            getline(cin, str);
            for (j = 1; j <= nmats; j = j + 1) {
                nreq[j - 1] = 0;
            }
            iss.clear();
            iss.str(str);
            for (j = 1; j <= nmats; j = j + 1) {
                iss >> itmp;
                nreq[j - 1] = itmp;
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
        //
    }
    return;
statement_90:
    write(nout, "(/,' *** End of file reached when trying to read matrix ','types for ',"
                "a3,/,' *** Check that you are requesting the',"
                "' right number of types for each path',/)"),
        path;
    write(nout, star);
    //
    //     End of Alareq
    //
}
