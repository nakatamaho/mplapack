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

// Derived from LAPACK routine IEEECK.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

INTEGER iMieeeck(INTEGER const ispec, REAL const zero, REAL const one) {
    INTEGER return_value = 0;
#if defined ___MPLAPACK_BUILD_WITH_GMP___
    // GMP uses arbitrary-precision integers/rationals internally and does not
    // implement IEEE 754 semantics: division by zero does not yield +Inf/-Inf,
    // and NaN/signed-zero behavior is undefined.  The runtime checks below
    // would invoke undefined behavior on GMP arithmetic, so return 0 here.
    return 0;
#endif
#if defined ___MPLAPACK_BUILD_WITH_DD___
    // DD (double-double) arithmetic does not comply with IEEE 754: it lacks
    // proper handling of infinities, NaN, and signed zero.  The runtime checks
    // below would invoke undefined behavior on DD arithmetic, so return 0 here.
    return 0;
#endif
#if defined ___MPLAPACK_BUILD_WITH_QD___
    // QD (quad-double) arithmetic does not comply with IEEE 754: it lacks
    // proper handling of infinities, NaN, and signed zero.  The runtime checks
    // below would invoke undefined behavior on QD arithmetic, so return 0 here.
    return 0;
#endif
    return_value = 1;
    //
    REAL posinf = one / zero;
    if (!Risinf(posinf) || posinf <= zero) {
        return_value = 0;
        return return_value;
    }
    //
    REAL neginf = -one / zero;
    if (!Risinf(neginf) || neginf >= zero) {
        return_value = 0;
        return return_value;
    }
    //
    REAL negzro = one / (neginf + one);
    if (negzro != zero) {
        return_value = 0;
        return return_value;
    }
    //
    neginf = one / negzro;
    if (!Risinf(neginf) || neginf >= zero) {
        return_value = 0;
        return return_value;
    }
    //
    REAL newzro = negzro + zero;
    if (newzro != zero) {
        return_value = 0;
        return return_value;
    }
    //
    posinf = one / newzro;
    if (!Risinf(posinf) || posinf <= zero) {
        return_value = 0;
        return return_value;
    }
    //
    neginf = neginf * posinf;
    if (!Risinf(neginf) || neginf >= zero) {
        return_value = 0;
        return return_value;
    }
    //
    posinf = posinf * posinf;
    if (!Risinf(posinf) || posinf <= zero) {
        return_value = 0;
        return return_value;
    }
    //
    // Return if we were only asked to check infinity arithmetic
    //
    if (ispec == 0) {
        return return_value;
    }
    //
    REAL nan1 = posinf + neginf;
    //
    REAL nan2 = posinf / neginf;
    //
    REAL nan3 = posinf / posinf;
    //
    REAL nan4 = posinf * zero;
    //
    REAL nan5 = neginf * negzro;
    //
    REAL nan6 = nan5 * zero;
    //
    if (!Risnan(nan1)) {
        return_value = 0;
        return return_value;
    }
    //
    if (!Risnan(nan2)) {
        return_value = 0;
        return return_value;
    }
    //
    if (!Risnan(nan3)) {
        return_value = 0;
        return return_value;
    }
    //
    if (!Risnan(nan4)) {
        return_value = 0;
        return return_value;
    }
    //
    if (!Risnan(nan5)) {
        return_value = 0;
        return return_value;
    }
    //
    if (!Risnan(nan6)) {
        return_value = 0;
        return return_value;
    }
    //
    return return_value;
}
