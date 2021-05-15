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

INTEGER iMlaenv(INTEGER const ispec, const char *name, const char * /* opts */, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const /* n4 */) {
    INTEGER return_value = 0;
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
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Arrays in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Save statement ..
    //     ..
    //     .. Executable Statements ..
    //
    std::string str = name;
    char *subname = new char[str.size() + 1];
    std::strcpy(subname, str.c_str());
    if (ispec >= 1 && ispec <= 5) {
        //
        //        Return a value from the common block.
        //
        if (strncmp(subname, "GEQR ", 5)) {
            if (n3 == 2) {
                return_value = iparms[2 - 1];
            } else {
                return_value = iparms[1 - 1];
            }
        } else if (strncmp(subname, "GELQ ", 5)) {
            if (n3 == 2) {
                return_value = iparms[2 - 1];
            } else {
                return_value = iparms[1 - 1];
            }
        } else {
            return_value = iparms[ispec - 1];
        }
        //
    } else if (ispec == 6) {
        //
        //        Compute SVD crossover point.
        //
        return_value = castINTEGER(castREAL(min(n1, n2 * 1.6e0f)));
        //
    } else if (ispec >= 7 && ispec <= 9) {
        //
        //        Return a value from the common block.
        //
        return_value = iparms[ispec - 1];
        //
    } else if (ispec == 10) {
        //
        //        IEEE NaN arithmetic can be trusted not to trap
        //
        //        iMlaenv = 0
        return_value = 1;
        if (return_value == 1) {
            return_value = iMieeeck(1, 0.0, 1.0);
        }
        //
    } else if (ispec == 11) {
        //
        //        Infinity arithmetic can be trusted not to trap
        //
        //        iMlaenv = 0
        return_value = 1;
        if (return_value == 1) {
            return_value = iMieeeck(1, 0.0, 1.0);
        }
        //
    } else {
        //
        //        Invalid value for ISPEC
        //
        return_value = -1;
    }
    //
    delete[] subname;
    return return_value;
    //
    //     End of iMlaenv
    //
}

INTEGER iMlaenv2stage(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4) {
    INTEGER return_value = 0;
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Local variables ..
    //     .. External Functions ..
    //     ..
    //     .. Arrays in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Save statement ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER iispec = 0;
    if ((ispec >= 1) && (ispec <= 5)) {
        //
        //     1 <= ISPEC <= 5: 2stage eigenvalues SVD routines.
        //
        if (ispec == 1) {
            return_value = iparms[1 - 1];
        } else {
            iispec = 16 + ispec;
            return_value = iMparam2stage(iispec, name, opts, n1, n2, n3, n4);
        }
        //
    } else {
        //
        //        Invalid value for ISPEC
        //
        return_value = -1;
    }
    //
    return return_value;
}
