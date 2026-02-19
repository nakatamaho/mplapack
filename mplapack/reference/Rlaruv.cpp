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
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <atomic>
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

// Utility: convert iseed[4] to a 64-bit seed.
// Each element contributes 16 bits (packed at bit offsets 0/16/32/48).
// iseed == {-1,-1,-1,-1} is reserved and rejected in deterministic mode.
namespace {
inline bool rlaruv_nondeterministic_enabled() {
    const char *v = std::getenv("MPLAPACK_RLARUV_NONDET");
    return (v != nullptr && v[0] == '1' && v[1] == '\0');
}

inline void rlaruv_print_nondet_banner_once() {
    static std::atomic<bool> printed{false};
    if (!printed.exchange(true)) {
        std::fprintf(stderr, "Rlaruv: non-deterministic mode enabled by MPLAPACK_RLARUV_NONDET=1\n");
        std::fflush(stderr);
    }
}

inline bool iseed_is_all_minus_one(const INTEGER *iseed) {
    return (iseed[0] == -1 && iseed[1] == -1 && iseed[2] == -1 && iseed[3] == -1);
}

inline uint64_t iseed_to_seed64(const INTEGER *iseed) {
    uint64_t s = static_cast<uint64_t>(static_cast<uint16_t>(iseed[0]));
    s |= static_cast<uint64_t>(static_cast<uint16_t>(iseed[1])) << 16;
    s |= static_cast<uint64_t>(static_cast<uint16_t>(iseed[2])) << 32;
    s |= static_cast<uint64_t>(static_cast<uint16_t>(iseed[3])) << 48;
    // splitmix64 finalizer to improve avalanche for small iseed values
    s ^= s >> 17;
    s *= 0xbf58476d1ce4e5b9ULL;
    s ^= s >> 31;
    s *= 0x94d049bb133111ebULL;
    s ^= s >> 32;
    return s;
}

// Update iseed so that consecutive calls with the same initial iseed
// produce different sequences. Each element is kept within 16 bits.
inline void advance_iseed(INTEGER *iseed, INTEGER n) {
    std::mt19937_64 tmp(iseed_to_seed64(iseed));
    tmp.discard(static_cast<unsigned long long>(n));
    uint64_t next = tmp();
    iseed[0] = static_cast<INTEGER>((next >> 0) & 0xFFFFu);
    iseed[1] = static_cast<INTEGER>((next >> 16) & 0xFFFFu);
    iseed[2] = static_cast<INTEGER>((next >> 32) & 0xFFFFu);
    iseed[3] = static_cast<INTEGER>((next >> 48) & 0xFFFFu);
}


// Deterministic real generation from mt19937_64 without std::uniform_real_distribution.
// Build x in [0,1) as a fixed-point fraction constructed from raw RNG bits.
// This avoids implementation-defined rounding differences across libstdc++/libc++ and libm.
template <typename T> inline T fixed_point_from_u32x4(std::mt19937_64 &mt) {
    // 128-bit fraction: (u0 + u1*2^-32 + u2*2^-64 + u3*2^-96) * 2^-32
    // Equivalent: sum_{i=1..4} u_i * 2^(-32*i)
    uint64_t a = mt();
    uint64_t b = mt();
    uint32_t u0 = static_cast<uint32_t>(a & 0xFFFFFFFFu);
    uint32_t u1 = static_cast<uint32_t>((a >> 32) & 0xFFFFFFFFu);
    uint32_t u2 = static_cast<uint32_t>(b & 0xFFFFFFFFu);
    uint32_t u3 = static_cast<uint32_t>((b >> 32) & 0xFFFFFFFFu);

    const T s32 = T(0x1p-32); // exact power of two
    T x = static_cast<T>(u0) * s32;
    x += static_cast<T>(u1) * (s32 * s32);
    x += static_cast<T>(u2) * (s32 * s32 * s32);
    x += static_cast<T>(u3) * (s32 * s32 * s32 * s32);
    return x;
}

template <typename T> inline T fixed_point_63bits(std::mt19937_64 &mt) {
    // 63-bit fraction: hi32 * 2^-32 + lo31 * 2^-63
    uint64_t u = mt();
    u >>= 1; // keep 63 bits
    uint32_t hi32 = static_cast<uint32_t>((u >> 31) & 0xFFFFFFFFu);
    uint32_t lo31 = static_cast<uint32_t>(u & 0x7FFFFFFFu);
    T x = static_cast<T>(hi32) * T(0x1p-32);
    x += static_cast<T>(lo31) * T(0x1p-63);
    return x;
}
} // namespace

void Rlaruv(INTEGER *iseed, INTEGER const n, REAL *x) {
    bool nondet = rlaruv_nondeterministic_enabled();
    if (nondet) {
        rlaruv_print_nondet_banner_once();
    } else {
        if (iseed_is_all_minus_one(iseed)) {
            Mxerbla("Rlaruv", 1);
            return;
        }
    }


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
        for (int i = 0; i < n; i++)
            x[i] = fixed_point_from_u32x4<REAL>(mt);
        advance_iseed(iseed, n);
    }
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
    if (nondet) {
        for (int i = 0; i < n; i++)
            x[i] = nondeterministic_rand() + nondeterministic_rand() * 0x1p-64;
    } else {
        std::mt19937_64 mt(iseed_to_seed64(iseed));
        for (int i = 0; i < n; i++)
            x[i] = fixed_point_63bits<REAL>(mt);
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
