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

// Derived from LAPACK routine XERBLA.
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

#define __MPLAPACK_MXERBLA__

#include <string>
#include <mplapack_matgen.h>
#include <mplapack_lin.h>
#include <mplapack_debug.h>

void Mxerbla(const char *srname, int info) {
    common cmn;
    common_write write(cmn);
    //
    // Trim trailing null bytes and spaces from srnamt.
    // fem::str<N>::operator std::string() copies all N bytes via std::string(elems, StrLen).
    // fem::str<N>::str() zero-fills elems with memset, so unset srnamt is all '\0'.
    // We must strip both '\0' and ' ' to get a clean C++ string.
    std::string srnamt_trimmed = static_cast<std::string>(srnamt);
    {
        size_t end = srnamt_trimmed.find_last_not_of(std::string("\0 ", 2));
        srnamt_trimmed = (end == std::string::npos) ? "" : srnamt_trimmed.substr(0, end + 1);
    }
    //
    lerr = true;
    if (info != infot) {
        if (infot != 0) {
            // srnamt was set: report which routine was expected and what info was received
            write(nout, "(' *** XERBLA was called from ',a,' with INFO = ',i6,' instead of ',"
                        "i2,' ***')"),
                srnamt_trimmed, info, infot;
        } else {
            // srnamt was not set: this is an unexpected XERBLA call from srname itself
            write(nout, "(' *** On entry to ',a,' parameter number ',i6,"
                        "' had an illegal value ***')"),
                srname, info;
        }
        ok = false;
    }
    // Check that XERBLA was called by the expected routine.
    // Guard with infot != 0: when infot == 0, srnamt was not set by the test
    // harness (all '\0' after memset), so comparing against srname is meaningless.
    if (infot == 0) {
        // srnamt was not set by the test harness: unexpected XERBLA call
        write(nout, "(' *** XERBLA was called with srname= ',a,' but no srname was expected (srnamt not set) ***')"), srname;
        ok = false;
    } else if (!(srnamt == srname)) {
        write(nout, "(' *** XERBLA was called with srname= ',a,' instead of ',a,' ***')"), srname, srnamt_trimmed;
        ok = false;
    }
    //
    // End of Mxerbla
    //
}
