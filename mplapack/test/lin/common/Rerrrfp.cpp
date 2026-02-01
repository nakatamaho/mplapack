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

// Derived from LAPACK routine DERRRFP.
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
#include <mplapack_lin.h>

void Rerrrfp(INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    INTEGER lda = 1;
    INTEGER ldb = 1;
    //
    static const char *format_9999 = "(1x,'DOUBLE PRECISION RFP routines passed the tests of ',"
                                     "'the error exits')";
    static const char *format_9998 = "(' *** RFP routines failed the tests of the error ','exits ***')";
    //
    nout = nunit;
    ok = true;
    REAL a[1 * 1];
    a[0] = 1.0;
    REAL b[1 * 1];
    b[0] = 1.0;
    REAL alpha = 1.0;
    REAL beta = 1.0;
    //
    srnamt = "Rpftrf";
    infot = 1;
    INTEGER info = 0;
    Rpftrf("/", "U", 0, a, info);
    Chkxer("Rpftrf", infot, nout, lerr, ok);
    infot = 2;
    Rpftrf("N", "/", 0, a, info);
    Chkxer("Rpftrf", infot, nout, lerr, ok);
    infot = 3;
    Rpftrf("N", "U", -1, a, info);
    Chkxer("Rpftrf", infot, nout, lerr, ok);
    //
    srnamt = "Rpftrs";
    infot = 1;
    Rpftrs("/", "U", 0, 0, a, b, 1, info);
    Chkxer("Rpftrs", infot, nout, lerr, ok);
    infot = 2;
    Rpftrs("N", "/", 0, 0, a, b, 1, info);
    Chkxer("Rpftrs", infot, nout, lerr, ok);
    infot = 3;
    Rpftrs("N", "U", -1, 0, a, b, 1, info);
    Chkxer("Rpftrs", infot, nout, lerr, ok);
    infot = 4;
    Rpftrs("N", "U", 0, -1, a, b, 1, info);
    Chkxer("Rpftrs", infot, nout, lerr, ok);
    infot = 7;
    Rpftrs("N", "U", 0, 0, a, b, 0, info);
    Chkxer("Rpftrs", infot, nout, lerr, ok);
    //
    srnamt = "Rpftri";
    infot = 1;
    Rpftri("/", "U", 0, a, info);
    Chkxer("Rpftri", infot, nout, lerr, ok);
    infot = 2;
    Rpftri("N", "/", 0, a, info);
    Chkxer("Rpftri", infot, nout, lerr, ok);
    infot = 3;
    Rpftri("N", "U", -1, a, info);
    Chkxer("Rpftri", infot, nout, lerr, ok);
    //
    srnamt = "Rtfsm";
    infot = 1;
    Rtfsm("/", "L", "U", "T", "U", 0, 0, alpha, a, b, 1);
    Chkxer("Rtfsm", infot, nout, lerr, ok);
    infot = 2;
    Rtfsm("N", "/", "U", "T", "U", 0, 0, alpha, a, b, 1);
    Chkxer("Rtfsm", infot, nout, lerr, ok);
    infot = 3;
    Rtfsm("N", "L", "/", "T", "U", 0, 0, alpha, a, b, 1);
    Chkxer("Rtfsm", infot, nout, lerr, ok);
    infot = 4;
    Rtfsm("N", "L", "U", "/", "U", 0, 0, alpha, a, b, 1);
    Chkxer("Rtfsm", infot, nout, lerr, ok);
    infot = 5;
    Rtfsm("N", "L", "U", "T", "/", 0, 0, alpha, a, b, 1);
    Chkxer("Rtfsm", infot, nout, lerr, ok);
    infot = 6;
    Rtfsm("N", "L", "U", "T", "U", -1, 0, alpha, a, b, 1);
    Chkxer("Rtfsm", infot, nout, lerr, ok);
    infot = 7;
    Rtfsm("N", "L", "U", "T", "U", 0, -1, alpha, a, b, 1);
    Chkxer("Rtfsm", infot, nout, lerr, ok);
    infot = 11;
    Rtfsm("N", "L", "U", "T", "U", 0, 0, alpha, a, b, 0);
    Chkxer("Rtfsm", infot, nout, lerr, ok);
    //
    srnamt = "Rtftri";
    infot = 1;
    Rtftri("/", "L", "N", 0, a, info);
    Chkxer("Rtftri", infot, nout, lerr, ok);
    infot = 2;
    Rtftri("N", "/", "N", 0, a, info);
    Chkxer("Rtftri", infot, nout, lerr, ok);
    infot = 3;
    Rtftri("N", "L", "/", 0, a, info);
    Chkxer("Rtftri", infot, nout, lerr, ok);
    infot = 4;
    Rtftri("N", "L", "N", -1, a, info);
    Chkxer("Rtftri", infot, nout, lerr, ok);
    //
    srnamt = "Rtfttr";
    infot = 1;
    Rtfttr("/", "U", 0, a, b, 1, info);
    Chkxer("Rtfttr", infot, nout, lerr, ok);
    infot = 2;
    Rtfttr("N", "/", 0, a, b, 1, info);
    Chkxer("Rtfttr", infot, nout, lerr, ok);
    infot = 3;
    Rtfttr("N", "U", -1, a, b, 1, info);
    Chkxer("Rtfttr", infot, nout, lerr, ok);
    infot = 6;
    Rtfttr("N", "U", 0, a, b, 0, info);
    Chkxer("Rtfttr", infot, nout, lerr, ok);
    //
    srnamt = "Rtrttf";
    infot = 1;
    Rtrttf("/", "U", 0, a, 1, b, info);
    Chkxer("Rtrttf", infot, nout, lerr, ok);
    infot = 2;
    Rtrttf("N", "/", 0, a, 1, b, info);
    Chkxer("Rtrttf", infot, nout, lerr, ok);
    infot = 3;
    Rtrttf("N", "U", -1, a, 1, b, info);
    Chkxer("Rtrttf", infot, nout, lerr, ok);
    infot = 5;
    Rtrttf("N", "U", 0, a, 0, b, info);
    Chkxer("Rtrttf", infot, nout, lerr, ok);
    //
    srnamt = "Rtfttp";
    infot = 1;
    Rtfttp("/", "U", 0, a, b, info);
    Chkxer("Rtfttp", infot, nout, lerr, ok);
    infot = 2;
    Rtfttp("N", "/", 0, a, b, info);
    Chkxer("Rtfttp", infot, nout, lerr, ok);
    infot = 3;
    Rtfttp("N", "U", -1, a, b, info);
    Chkxer("Rtfttp", infot, nout, lerr, ok);
    //
    srnamt = "Rtpttf";
    infot = 1;
    Rtpttf("/", "U", 0, a, b, info);
    Chkxer("Rtpttf", infot, nout, lerr, ok);
    infot = 2;
    Rtpttf("N", "/", 0, a, b, info);
    Chkxer("Rtpttf", infot, nout, lerr, ok);
    infot = 3;
    Rtpttf("N", "U", -1, a, b, info);
    Chkxer("Rtpttf", infot, nout, lerr, ok);
    //
    srnamt = "Rtrttp";
    infot = 1;
    Rtrttp("/", 0, a, 1, b, info);
    Chkxer("Rtrttp", infot, nout, lerr, ok);
    infot = 2;
    Rtrttp("U", -1, a, 1, b, info);
    Chkxer("Rtrttp", infot, nout, lerr, ok);
    infot = 4;
    Rtrttp("U", 0, a, 0, b, info);
    Chkxer("Rtrttp", infot, nout, lerr, ok);
    //
    srnamt = "Rtpttr";
    infot = 1;
    Rtpttr("/", 0, a, b, 1, info);
    Chkxer("Rtpttr", infot, nout, lerr, ok);
    infot = 2;
    Rtpttr("U", -1, a, b, 1, info);
    Chkxer("Rtpttr", infot, nout, lerr, ok);
    infot = 5;
    Rtpttr("U", 0, a, b, 0, info);
    Chkxer("Rtpttr", infot, nout, lerr, ok);
    //
    srnamt = "Rsfrk";
    infot = 1;
    Rsfrk("/", "U", "N", 0, 0, alpha, a, 1, beta, b);
    Chkxer("Rsfrk", infot, nout, lerr, ok);
    infot = 2;
    Rsfrk("N", "/", "N", 0, 0, alpha, a, 1, beta, b);
    Chkxer("Rsfrk", infot, nout, lerr, ok);
    infot = 3;
    Rsfrk("N", "U", "/", 0, 0, alpha, a, 1, beta, b);
    Chkxer("Rsfrk", infot, nout, lerr, ok);
    infot = 4;
    Rsfrk("N", "U", "N", -1, 0, alpha, a, 1, beta, b);
    Chkxer("Rsfrk", infot, nout, lerr, ok);
    infot = 5;
    Rsfrk("N", "U", "N", 0, -1, alpha, a, 1, beta, b);
    Chkxer("Rsfrk", infot, nout, lerr, ok);
    infot = 8;
    Rsfrk("N", "U", "N", 0, 0, alpha, a, 0, beta, b);
    Chkxer("Rsfrk", infot, nout, lerr, ok);
    //
    // Print a summary line.
    //
    if (ok) {
        write(nout, format_9999);
    } else {
        write(nout, format_9998);
    }
    //
    // End of Rerrrfp
    //
}
