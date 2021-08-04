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

inline REAL abs1(COMPLEX zdum) { return abs(zdum.real()) + abs(zdum.imag()); }

void Cget04(INTEGER const n, INTEGER const nrhs, COMPLEX *x, INTEGER const ldx, COMPLEX *xact, INTEGER const ldxact, REAL const rcond, REAL &resid) {
    COMPLEX zdum = 0.0;
    //
    //     Quick exit if N = 0 or NRHS = 0.
    //
    const REAL zero = 0.0;
    if (n <= 0 || nrhs <= 0) {
        resid = zero;
        return;
    }
    //
    //     Exit with RESID = 1/EPS if RCOND is invalid.
    //
    REAL eps = Rlamch("Epsilon");
    if (rcond < zero) {
        resid = 1.0 / eps;
        return;
    }
    //
    //     Compute the maximum of
    //        norm(X - XACT) / ( norm(XACT) * EPS )
    //     over all the vectors X and XACT .
    //
    resid = zero;
    INTEGER j = 0;
    INTEGER ix = 0;
    REAL xnorm = 0.0;
    REAL diffnm = 0.0;
    INTEGER i = 0;
    for (j = 1; j <= nrhs; j = j + 1) {
        ix = iCamax(n, &xact[(j - 1) * ldxact], 1);
        xnorm = abs1(xact[(ix - 1) + (j - 1) * ldxact]);
        diffnm = zero;
        for (i = 1; i <= n; i = i + 1) {
            diffnm = max(diffnm, abs1(x[(i - 1) + (j - 1) * ldx] - xact[(i - 1) + (j - 1) * ldxact]));
        }
        if (xnorm <= zero) {
            if (diffnm > zero) {
                resid = 1.0 / eps;
            }
        } else {
            resid = max(resid, REAL((diffnm / xnorm) * rcond));
        }
    }
    if (resid * eps < 1.0) {
        resid = resid / eps;
    }
    //
    //     End of Cget04
    //
}
