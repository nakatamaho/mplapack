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

void Rlafts(const char *type, INTEGER const m, INTEGER const n, INTEGER const imat, INTEGER const ntests, REAL *result, INTEGER *iseed, REAL const thresh, INTEGER const iounit, INTEGER &ie) {
    common cmn;
    common_write write(cmn);
    char buf[1024];
    //
    INTEGER k = 0;
    if (m == n) {
        //
        //     Output for square matrices:
        //
        for (k = 1; k <= ntests; k = k + 1) {
            if (result[k - 1] >= thresh) {
                //
                //           If this is the first test to fail, call Rlahd2
                //           to print a header to the data file.
                //
                if (ie == 0) {
                    Rlahd2(iounit, type);
                }
                ie++;
                if (result[k - 1] < 10000.0) {
                    sprintnum_short(buf, result[k - 1]);
                    write(iounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),' result ',"
                                  "i3,' is ',0p,a)"),
                        n, imat, iseed[0], iseed[1], iseed[2], iseed[3], k, buf;
                } else {
                    sprintnum_short(buf, result[k - 1]);
                    write(iounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),' result ',"
                                  "i3,' is ',1p,a)"),
                        n, imat, iseed[0], iseed[1], iseed[2], iseed[3], k, buf;
                }
            }
        }
    } else {
        //
        //     Output for rectangular matrices
        //
        for (k = 1; k <= ntests; k = k + 1) {
            if (result[k - 1] >= thresh) {
                //
                //              If this is the first test to fail, call Rlahd2
                //              to print a header to the data file.
                //
                if (ie == 0) {
                    Rlahd2(iounit, type);
                }
                ie++;
                if (result[k - 1] < 10000.0) {
                    sprintnum_short(buf, result[k - 1]);
                    write(iounit, "(1x,i5,' x',i5,' matrix, type=',i2,', seed=',3(i4,','),i4,"
                                  "': result ',i3,' is ',0p,a)"),
                        m, n, imat, iseed[0], iseed[1], iseed[2], iseed[3], k, buf;
                } else {
                    sprintnum_short(buf, result[k - 1]);
                    write(iounit, "(1x,i5,' x',i5,' matrix, type=',i2,', seed=',3(i4,','),i4,"
                                  "': result ',i3,' is ',1p,a)"),
                        m, n, imat, iseed[0], iseed[1], iseed[2], iseed[3], k, buf;
                }
            }
        }
        //
    }
    //
    //     End of Rlafts
    //
}
