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
#include <atomic>
#include <random>

// For DD, QD, binary128, binary80, and double backends:
// Non-deterministic engine, seeded once via std::random_device.
#if defined MPLAPACK_BUILD_WITH_DD || defined MPLAPACK_BUILD_WITH_QD || defined MPLAPACK_BUILD_WITH_BINARY128 || defined MPLAPACK_BUILD_WITH_BINARY80 || defined MPLAPACK_BUILD_WITH_DOUBLE
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

inline uint64_t rlaruv_nondeterministic_seed64() {
    std::random_device source;
    const uint64_t high = static_cast<uint64_t>(source());
    const uint64_t low = static_cast<uint64_t>(source());
    return (high << 32) ^ low;
}

#if defined MPLAPACK_BUILD_WITH_MPFR
inline mpfrxx::mpfr_randclass &rlaruv_mpfr_nondeterministic_state() {
    thread_local mpfrxx::mpfr_randclass state(gmp_randinit_default);
    thread_local const bool seeded = (state.seed_u64(rlaruv_nondeterministic_seed64()), true);
    (void)seeded;
    return state;
}
#endif

#if defined MPLAPACK_BUILD_WITH_GMP
inline gmpxx::gmp_randclass &rlaruv_gmp_nondeterministic_state() {
    thread_local gmpxx::gmp_randclass state(gmp_randinit_default);
    thread_local const bool seeded = (state.seed_u64(rlaruv_nondeterministic_seed64()), true);
    (void)seeded;
    return state;
}
#endif

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

#if defined MPLAPACK_BUILD_WITH_MPFR
    if (nondet) {
        mpfrxx::mpfr_randclass &rng = rlaruv_mpfr_nondeterministic_state();
        for (int i = 0; i < n; i++)
            x[i] = rng.get_fr();
    } else {
        mpfrxx::mpfr_randclass rng(gmp_randinit_default);
        rng.seed_u64(iseed_to_seed64(iseed));
        for (int i = 0; i < n; i++)
            x[i] = rng.get_fr();
        advance_iseed(iseed, n);
    }
#endif

#if defined MPLAPACK_BUILD_WITH_GMP
    if (nondet) {
        gmpxx::gmp_randclass &rng = rlaruv_gmp_nondeterministic_state();
        for (int i = 0; i < n; i++)
            x[i] = rng.get_f();
    } else {
        gmpxx::gmp_randclass rng(gmp_randinit_default);
        rng.seed_u64(iseed_to_seed64(iseed));
        for (int i = 0; i < n; i++)
            x[i] = rng.get_f();
        advance_iseed(iseed, n);
    }
#endif

#if defined MPLAPACK_BUILD_WITH_DD
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

#if defined MPLAPACK_BUILD_WITH_QD
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

#if defined MPLAPACK_BUILD_WITH_BINARY128
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

#if defined MPLAPACK_BUILD_WITH_BINARY80
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

#if defined MPLAPACK_BUILD_WITH_DOUBLE
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
