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
#include <mplapack_common.h>

#include <mplapack_matgen.h>
#include <mplapack_eig.h>

#include <mplapack_debug.h>

bool Rlctsx(common &cmn, REAL const /* ar */, REAL const /* ai */, REAL const /* beta */) {
    bool return_value = false;
    INTEGER mplusn = cmn.mplusn;
    INTEGER m = cmn.m;
    INTEGER i = cmn.i;
    INTEGER n = cmn.n;
    bool fs = cmn.fs;
    //
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Save statement ..
    //     ..
    //     .. Executable Statements ..
    //
    if (fs) {
        i++;
        if (i <= m) {
            return_value = false;
        } else {
            return_value = true;
        }
        if (i == mplusn) {
            fs = false;
            i = 0;
        }
    } else {
        i++;
        if (i <= cmn.n) {
            return_value = true;
        } else {
            return_value = false;
        }
        if (i == mplusn) {
            fs = true;
            i = 0;
        }
    }
    //
    //       IF( AR/BETA.GT.0.0 )THEN
    //          Rlctsx = .TRUE.
    //       ELSE
    //          Rlctsx = .FALSE.
    //       END IF
    //
    return return_value;
    //
    //     End of Rlctsx
    //
}
