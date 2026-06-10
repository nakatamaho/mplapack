/*
 * Copyright (c) 2008-2026
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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include <mpblas.h>
#include <mplapack.h>

static void fail(const char *what) {
    printf("*** Testing Rnaninf failed: %s ***\n", what);
    exit(1);
}

static void check(bool cond, const char *what) {
    if (!cond)
        fail(what);
}

int main() {
    printf("*** Testing Rnaninf start ***\n");

#if defined ___MPLAPACK_BUILD_WITH_GMP___
    REAL finite_value = 1.0;

    check(!Risnan(finite_value), "GMP finite value recognized as NaN");
    check(!Risinf(finite_value), "GMP finite value recognized as Inf");
#else
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
    REAL nan_value(NAN, 71, MPFR_RNDN);
    REAL pos_inf_value(INFINITY, 71, MPFR_RNDN);
    REAL neg_inf_value(-INFINITY, 71, MPFR_RNDN);
    REAL finite_value(1.0, 71, MPFR_RNDN);

    check(mpfr_nan_p(mpfr_ptr(nan_value)) != 0, "mpfr_nan_p did not recognize NaN test value");
#else
    REAL nan_value = std::numeric_limits<double>::quiet_NaN();
    REAL pos_inf_value = std::numeric_limits<double>::infinity();
    REAL neg_inf_value = -std::numeric_limits<double>::infinity();
    REAL finite_value = 1.0;
#endif

    check(Risnan(nan_value), "Risnan did not recognize NaN");
    check(!Risnan(finite_value), "Risnan recognized finite value as NaN");
    check(!Risnan(pos_inf_value), "Risnan recognized +Inf as NaN");
    check(!Risnan(neg_inf_value), "Risnan recognized -Inf as NaN");

    check(Risinf(pos_inf_value), "Risinf did not recognize +Inf");
    check(Risinf(neg_inf_value), "Risinf did not recognize -Inf");
    check(!Risinf(finite_value), "Risinf recognized finite value as Inf");
    check(!Risinf(nan_value), "Risinf recognized NaN as Inf");
#endif

    printf("*** Testing Rnaninf successful ***\n");
    return 0;
}
