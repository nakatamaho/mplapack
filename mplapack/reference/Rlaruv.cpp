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
#include <cstdint>
#include <random>

#if defined ___MPLAPACK_BUILD_WITH_MPFR___ || defined ___MPLAPACK_BUILD_WITH_GMP___
#include <ctime>
#endif

// MPFR backend: global non-deterministic state
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
gmp_randstate_t ___random_mplapack_mpfr_state;
void __attribute__((constructor)) ___mplapack_Rlaruv_mpfr_initialize(void) {
    gmp_randinit_default(___random_mplapack_mpfr_state);
    gmp_randseed_ui(___random_mplapack_mpfr_state, static_cast<unsigned long int>(std::time(nullptr)));
}
void __attribute__((destructor)) ___mplapack_Rlaruv_mpfr_finalize(void) { gmp_randclear(___random_mplapack_mpfr_state); }
#endif

// GMP backend: global non-deterministic state
// Note: C++ guarantees globals in the same TU are initialized in definition order.
// gmp_randclass is constructed first, then the seeder initializes it.
#if defined ___MPLAPACK_BUILD_WITH_GMP___
gmp_randclass ___random_mplapack_gmp(gmp_randinit_default);
namespace {
struct ___mplapack_Rlaruv_gmp_seeder {
    ___mplapack_Rlaruv_gmp_seeder() { ___random_mplapack_gmp.seed(static_cast<unsigned long int>(std::time(nullptr))); }
} ___mplapack_Rlaruv_gmp_seeder_instance;
} // namespace
#endif

// For DD, QD, binary128, binary80, and double backends:
// Non-deterministic engine, seeded once via std::random_device.
#if defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH_QD___ || defined ___MPLAPACK_BUILD_WITH_BINARY128___ || defined ___MPLAPACK_BUILD_WITH_BINARY80___ || defined ___MPLAPACK_BUILD_WITH_DOUBLE___
namespace {
inline double nondeterministic_rand() {
    static std::mt19937_64 mt(std::random_device{}());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(mt);
}
} // namespace
#endif

// Utility: convert LAPACK-style iseed[4] to a 64-bit seed.
// LAPACK convention: iseed[0..2] in [0, 4095], iseed[3] odd and in [0, 4095].
// iseed == {0,0,0,0} selects non-deterministic mode.
namespace {
inline bool iseed_is_nondeterministic(INTEGER *iseed) { return (iseed[0] == 0 && iseed[1] == 0 && iseed[2] == 0 && iseed[3] == 0); }

inline uint64_t iseed_to_seed64(INTEGER *iseed) {
    uint64_t s = static_cast<uint64_t>(iseed[0] & 0xFFF);
    s |= static_cast<uint64_t>(iseed[1] & 0xFFF) << 12;
    s |= static_cast<uint64_t>(iseed[2] & 0xFFF) << 24;
    s |= static_cast<uint64_t>(iseed[3] & 0xFFF) << 36;
    // splitmix64 finalizer to improve avalanche for small iseed values
    s ^= s >> 17;
    s *= 0xbf58476d1ce4e5b9ULL;
    s ^= s >> 31;
    s *= 0x94d049bb133111ebULL;
    s ^= s >> 32;
    return s;
}

// Update iseed so that consecutive calls with the same initial iseed
// produce different sequences. Values stay in LAPACK convention.
inline void advance_iseed(INTEGER *iseed, INTEGER n) {
    std::mt19937_64 tmp(iseed_to_seed64(iseed));
    tmp.discard(static_cast<unsigned long long>(n));
    auto next = tmp();
    iseed[0] = static_cast<INTEGER>((next >> 0) & 0xFFF);
    iseed[1] = static_cast<INTEGER>((next >> 12) & 0xFFF);
    iseed[2] = static_cast<INTEGER>((next >> 24) & 0xFFF);
    iseed[3] = static_cast<INTEGER>(((next >> 36) & 0xFFF) | 1); // ensure odd
}
} // namespace

void Rlaruv(INTEGER *iseed, INTEGER const n, REAL *x) {
    bool nondet = iseed_is_nondeterministic(iseed);

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
    if (nondet) {
        for (int i = 0; i < n; i++)
            x[i] = urandom(___random_mplapack_mpfr_state);
    } else {
        gmp_randstate_t det_state;
        gmp_randinit_default(det_state);
        uint64_t s64 = iseed_to_seed64(iseed);
        mpz_t seed_z;
        mpz_init(seed_z);
        mpz_import(seed_z, 1, 1, sizeof(uint64_t), 0, 0, &s64);
        gmp_randseed(det_state, seed_z);
        mpz_clear(seed_z);
        for (int i = 0; i < n; i++)
            x[i] = urandom(det_state);
        gmp_randclear(det_state);
        advance_iseed(iseed, n);
    }
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
    if (nondet) {
        for (int i = 0; i < n; i++)
            x[i] = ___random_mplapack_gmp.get_f();
    } else {
        gmp_randclass det_rng(gmp_randinit_default);
        uint64_t s64 = iseed_to_seed64(iseed);
        mpz_class seed_z;
        mpz_import(seed_z.get_mpz_t(), 1, 1, sizeof(uint64_t), 0, 0, &s64);
        det_rng.seed(seed_z);
        for (int i = 0; i < n; i++)
            x[i] = det_rng.get_f();
        advance_iseed(iseed, n);
    }
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
    if (nondet) {
        for (int i = 0; i < n; i++) {
            x[i].x[0] = nondeterministic_rand();
            x[i].x[1] = nondeterministic_rand() * 0x1p-53;
        }
    } else {
        std::mt19937_64 mt(iseed_to_seed64(iseed));
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < n; i++) {
            x[i].x[0] = dist(mt);
            x[i].x[1] = dist(mt) * 0x1p-53;
        }
        advance_iseed(iseed, n);
    }
#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
    if (nondet) {
        for (int i = 0; i < n; i++) {
            x[i].x[0] = nondeterministic_rand();
            x[i].x[1] = nondeterministic_rand() * 0x1p-53;
            x[i].x[2] = nondeterministic_rand() * 0x1p-106;
            x[i].x[3] = nondeterministic_rand() * 0x1p-159;
        }
    } else {
        std::mt19937_64 mt(iseed_to_seed64(iseed));
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < n; i++) {
            x[i].x[0] = dist(mt);
            x[i].x[1] = dist(mt) * 0x1p-53;
            x[i].x[2] = dist(mt) * 0x1p-106;
            x[i].x[3] = dist(mt) * 0x1p-159;
        }
        advance_iseed(iseed, n);
    }
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
    if (nondet) {
        for (int i = 0; i < n; i++)
            x[i] = nondeterministic_rand() + nondeterministic_rand() * 0x1p-53 + nondeterministic_rand() * 0x1p-106;
    } else {
        std::mt19937_64 mt(iseed_to_seed64(iseed));
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < n; i++)
            x[i] = dist(mt) + dist(mt) * 0x1p-53 + dist(mt) * 0x1p-106;
        advance_iseed(iseed, n);
    }
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
    if (nondet) {
        for (int i = 0; i < n; i++)
            x[i] = nondeterministic_rand() + nondeterministic_rand() * 0x1p-64;
    } else {
        std::mt19937_64 mt(iseed_to_seed64(iseed));
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < n; i++)
            x[i] = dist(mt) + dist(mt) * 0x1p-64;
        advance_iseed(iseed, n);
    }
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
    if (nondet) {
        for (int i = 0; i < n; i++)
            x[i] = nondeterministic_rand();
    } else {
        std::mt19937_64 mt(iseed_to_seed64(iseed));
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        for (int i = 0; i < n; i++)
            x[i] = dist(mt);
        advance_iseed(iseed, n);
    }
#endif
}
