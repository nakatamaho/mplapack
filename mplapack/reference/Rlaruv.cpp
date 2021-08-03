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

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
gmp_randstate_t random_mplapack_mpfr_state;
void __attribute__((constructor)) mplapack_Rlaruv_mpfr_initialize(void);
void mplapack_Rlaruv_mpfr_initialize(void) { gmp_randinit_default(random_mplapack_mpfr_state); } // this is gmp_randinit_mt
void __attribute__((destructor)) mplapack_Rlaruv_mpfr_finalize(void);
void mplapack_Rlaruv_mpfr_finalize(void) { gmp_randclear(random_mplapack_mpfr_state); mpfr_free_cache(); } // this is gmp_randinit_mt
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
gmp_randstate_t random_mplapack_gmp_state;
gmp_randclass random_mplapack_gmp(gmp_randinit_default);
void __attribute__((constructor)) mplapack_Rlaruv_gmp_initialize(void);
void mplapack_Rlaruv_gmp_initialize(void) { random_mplapack_gmp.seed((unsigned long int)time(NULL)); } // XXX better initializaition req'ed
#endif

void Rlaruv(INTEGER *iseed, INTEGER const n, REAL *x) {

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
    for (int i = 0; i < n; i++)
        x[i] = urandom(random_mplapack_mpfr_state);
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
    for (int i = 0; i < n; i++)
        x[i] = random_mplapack_gmp.get_f();
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
