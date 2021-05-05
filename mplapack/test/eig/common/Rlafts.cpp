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
    iseed([4]);
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
    //     .. Local Scalars ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
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
                //           to prINTEGER a header to the data file.
                //
                if (ie == 0) {
                    Rlahd2(iounit, type);
                }
                ie++;
                if (result[k - 1] < 10000.0) {
                    write(iounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),' result ',"
                                  "i3,' is',0p,f8.2)"),
                        n, imat, iseed, k, result(k);
                } else {
                    write(iounit, "(' Matrix order=',i5,', type=',i2,', seed=',4(i4,','),' result ',"
                                  "i3,' is',1p,d10.3)"),
                        n, imat, iseed, k, result(k);
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
                //              to prINTEGER a header to the data file.
                //
                if (ie == 0) {
                    Rlahd2(iounit, type);
                }
                ie++;
                if (result[k - 1] < 10000.0) {
                    write(iounit, "(1x,i5,' x',i5,' matrix, type=',i2,', s','eed=',3(i4,','),i4,"
                                  "': result ',i3,' is',0p,f8.2)"),
                        m, n, imat, iseed, k, result(k);
                } else {
                    write(iounit, "(1x,i5,' x',i5,' matrix, type=',i2,', s','eed=',3(i4,','),i4,"
                                  "': result ',i3,' is',1p,d10.3)"),
                        m, n, imat, iseed, k, result(k);
                }
            }
        }
        //
    }
    //
    //     End of Rlafts
    //
}
