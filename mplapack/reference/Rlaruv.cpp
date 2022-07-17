/*
 * Copyright (c) 2008-2021
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
#include <random>
#include <time.h>

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
gmp_randstate_t ___random_mplapack_mpfr_state;
void __attribute__((constructor)) ___mplapack_Rlaruv_mpfr_initialize(void) {
    gmp_randinit_default(___random_mplapack_mpfr_state);
    gmp_randseed_ui(___random_mplapack_mpfr_state, (unsigned long int)time(NULL));
}
void __attribute__((destructor)) ___mplapack_Rlaruv_mpfr_finalize(void) {
    gmp_randclear(___random_mplapack_mpfr_state);
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
gmp_randstate_t ___random_mplapack_gmp_state;
gmp_randclass ___random_mplapack_gmp(gmp_randinit_default);
void __attribute__((constructor)) ___mplapack_Rlaruv_gmp_initialize(void) {
    gmp_randinit_default(___random_mplapack_gmp_state);
    gmp_randseed_ui(___random_mplapack_gmp_state, (unsigned long int)time(NULL));
}
void __attribute__((destructor)) ___mplapack_Rlaruv_gmp_finalize(void) {
    gmp_randclear(___random_mplapack_gmp_state);
}
#endif

void Rlaruv(INTEGER *iseed, INTEGER const n, REAL *x) {

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
    for (int i = 0; i < n; i++)
        x[i] = urandom(___random_mplapack_mpfr_state);
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
    for (int i = 0; i < n; i++)
        x[i] = ___random_mplapack_gmp.get_f();
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < n; i++) {
        x[i].x[0] = dist(mt);
        x[i].x[1] = dist(mt) * 0x1p-53;
    }
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < n; i++) {
        x[i].x[0] = dist(mt);
        x[i].x[1] = dist(mt) * 0x1p-53;
        x[i].x[2] = dist(mt) * 0x1p-106;
        x[i].x[3] = dist(mt) * 0x1p-159;
    }
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < n; i++)
        x[i] = dist(mt) + dist(mt) * 0x1p-53 + dist(mt) * 0x1p-106;
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> dist(0, 1.0);
    for (int i = 0; i < n; i++)
        x[i] = dist(mt) + dist(mt) * 0x1p-64;
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < n; i++)
        x[i] = dist(mt);
#endif
}
