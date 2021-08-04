/*
 * Copyright (c) 2012-2021
 *	Nakata, Maho
 * 	All rights reserved.
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
#include <stdio.h>

#if defined ___MPLAPACK_BUILD_WITH_GMP___
void __attribute__((constructor)) mplapack_initialize_gmp(void);
void __attribute__((destructor)) mplapack_finalize_gmp(void);
void mplapack_initialize_gmp(void) {
    char *p = getenv("MPLAPACK_GMP_PRECISION");
    if (p) {
        mpf_set_default_prec(atoi(p));
    }
}
void mplapack_finalize_gmp(void) {
    // no finalization needed
}
#endif


#if defined ___MPLAPACK_BUILD_WITH_MPFR___
void __attribute__((constructor)) mplapack_initialize_mpfr(void);
void mplapack_initialize_mpfr(void) {
    char *p = getenv("MPLAPACK_MPFR_PRECISION");
    if (p) {
        mpreal::set_default_prec(atoi(p));
        mpcomplex::set_default_prec(atoi(p));
    }
}
void __attribute__((destructor)) mplapack_finalize_mpfr(void);
void mplapack_finalize_mpfr(void) {}
#endif


#if defined ___MPLAPACK_BUILD_WITH_QD___ || defined ___MPLAPACK_BUILD_WITH_DD___
void __attribute__((constructor)) mplapack_initialize_qd(void);
void __attribute__((destructor)) mplapack_finalize_qd(void);
static unsigned int oldcw;
void mplapack_initialize_qd(void) { fpu_fix_start(&oldcw); }
void mplapack_finalize_qd(void) { fpu_fix_end(&oldcw); }
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
void __attribute__((constructor)) mplapack_initialize_double(void);
void __attribute__((destructor)) mplapack_finalize_double(void);
void mplapack_initialize_double(void) {
    // no initializization needed
}

void mplapack_finalize_double(void) {
    // no finalization needed
}
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
void __attribute__((constructor)) mplapack_initialize__Float64x(void);
void __attribute__((destructor)) mplapack_finalize__Float64x(void);
void mplapack_initialize__Float64x(void) {
    // no initializization needed
}

void mplapack_finalize__Float64x(void) {
    // no finalization needed
}
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___
void __attribute__((constructor)) mplapack_initialize_binary128(void);
void __attribute__((destructor)) mplapack_finalize_binary128(void);
void mplapack_initialize__Float128(void) {
    // no initializization needed
}

void mplapack_finalize__Float128(void) {
    // no finalization needed
}
#endif
