/*
 * Copyright (c) 2008-2026
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

#include <iostream>
#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_compare_debug.h>
#include <blas.h>
#include <lapack.h>

#define VERBOSE_TEST

/* dlamch result on FreeBSD 6/i386 + Lapack 3.1.1
  Epsilon                      =   1.110223024625157E-016
  Safe minimum                 =   2.225073858507201E-308
  Base                         =    2.00000000000000
  Precision                    =   2.220446049250313E-016
  Number of digits in mantissa =    53.0000000000000
  Rounding mode                =    1.00000000000000
  Minimum exponent             =   -1021.00000000000
  Underflow threshold          =   2.225073858507201E-308
  Largest exponent             =    1024.00000000000
  Overflow threshold           =   1.797693134862316E+308
*/

// -----------------------------------------------------------------------------
// Rlamch test: Comprehensive test for machine constants
//
// Test strategy:
//   E (eps):   Unit roundoff. Verify 1 + eps/2 == 1 and 1 + eps > 1.
//   S (sfmin): Safe minimum. Verify 1/sfmin is finite and sfmin >= max(rmin, 1/rmax).
//   B (base):  Must be 2.
//   P (prec):  Must equal eps * base.
//   N (digits):Must match mpf_get_prec().
//   R (round): Must be 1 (rounding occurs).
//   M (emin):  Verify consistency with rmin = base^(emin-1).
//   U (rmin):  Underflow threshold. Verify rmin > 0 and basic scaling.
//   L (emax):  Verify consistency with rmax = (1-eps) * base^emax.
//   O (rmax):  Overflow threshold. Verify rmax > 0 and 1/rmax > 0.
// -----------------------------------------------------------------------------

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
#include <climits>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <type_traits>
#include <utility>
#include <mpfr.h>

namespace {

// Fail fast with context (used by MPFR-only tests).
[[noreturn]] static void mpfr_test_fail(const char *what) {
    printf("*** Testing Mutils (MPFR) failed: %s ***\n", what);
    exit(1);
}

static void mpfr_assert(bool cond, const char *what) {
    if (!cond)
        mpfr_test_fail(what);
}

static void mpfr_assert_case(bool cond, const char *tag, const char *what) {
    if (cond)
        return;
    char buf[256];
    snprintf(buf, sizeof(buf), "%s: %s", tag, what);
    mpfr_test_fail(buf);
}

static bool rnd_rounds_up_for_positive(mpfr_rnd_t rnd) {
    if (rnd == MPFR_RNDU)
        return true;
#ifdef MPFR_RNDA
    if (rnd == MPFR_RNDA)
        return true;
#endif
    return false;
}

static long checked_exp_to_long(mpfr_exp_t e, const char *what) {
    // mpfr_mul_2si() / mpfr_set_si() accept long exponents.
    // If mpfr_exp_t is wider than long, guard against truncation.
    if (sizeof(mpfr_exp_t) > sizeof(long)) {
        if (e > static_cast<mpfr_exp_t>(LONG_MAX) || e < static_cast<mpfr_exp_t>(LONG_MIN)) {
            char buf[256];
            snprintf(buf, sizeof(buf), "%s (mpfr_exp_t out of range for long)", what);
            mpfr_test_fail(buf);
        }
    }
    return static_cast<long>(e);
}

static unsigned long checked_prec_to_ulong(mpfr_prec_t p, const char *what) {
    // mpfr_div_2ui() uses an unsigned long shift count.
    if (sizeof(mpfr_prec_t) > sizeof(unsigned long)) {
        if (p > static_cast<mpfr_prec_t>(std::numeric_limits<unsigned long>::max())) {
            char buf[256];
            snprintf(buf, sizeof(buf), "%s (mpfr_prec_t out of range for unsigned long)", what);
            mpfr_test_fail(buf);
        }
    }
    return static_cast<unsigned long>(p);
}

static mpfr_prec_t current_real_default_prec() {
    // Capture the precision of newly created REAL values. In some builds, MPLAPACK's
    // mpfr::mpreal default precision may differ from mpfr_get_default_prec().
    REAL probe = 1.0;
    return mpfr_get_prec(mpfr_ptr(probe));
}

template <typename T, typename = void> struct has_get_default_rnd : std::false_type {};

template <typename T> struct has_get_default_rnd<T, std::void_t<decltype(T::get_default_rnd())>> : std::true_type {};

template <typename T, typename = void> struct has_set_default_rnd : std::false_type {};

template <typename T> struct has_set_default_rnd<T, std::void_t<decltype(T::set_default_rnd(std::declval<mpfr_rnd_t>()))>> : std::true_type {};

static mpfr_rnd_t current_real_default_rnd() {
    // Prefer mpfr::mpreal's own default rounding mode if available.
    if constexpr (has_get_default_rnd<mpfr::mpreal>::value) {
        return mpfr::mpreal::default_rnd;
    } else {
        // Fallback: assume mpfr::mpreal follows MPFR's global default rounding mode.
        return mpfr_get_default_rounding_mode();
    }
}

static void set_real_default_rnd(mpfr_rnd_t rnd) {
    mpfr::mpreal::default_rnd = rnd;
    mpfr_set_default_rounding_mode(rnd);
}

struct MpfrEnvSnapshot {
    mpfr_prec_t mpfr_prec;
    mpfr_prec_t real_prec;
    mpfr_exp_t emin;
    mpfr_exp_t emax;
    mpfr_rnd_t mpfr_rnd;
    mpfr_rnd_t real_rnd;
};

static MpfrEnvSnapshot mpfr_env_capture() {
    MpfrEnvSnapshot s;
    s.mpfr_prec = mpreal::default_prec; // XXX mpreal.h mpreal::default_prec is not in sync with mpfr_get_default_prec();
                                        // s.mpfr_prec = mpfr_get_default_prec();
    s.real_prec = current_real_default_prec();
    s.emin = mpfr_get_emin();
    s.emax = mpfr_get_emax();
    s.mpfr_rnd = mpfr_get_default_rounding_mode();
    s.real_rnd = current_real_default_rnd();
    return s;
}

static void mpfr_env_apply(const MpfrEnvSnapshot &s) {
    // Apply rounding modes first (they affect subsequent computations).
    mpfr_set_default_rounding_mode(s.mpfr_rnd);
    set_real_default_rnd(s.real_rnd);

    // Apply exponent range.
    if (mpfr_set_emin(s.emin) != 0) {
        mpfr_test_fail("mpfr_set_emin() rejected the requested value");
    }
    if (mpfr_set_emax(s.emax) != 0) {
        mpfr_test_fail("mpfr_set_emax() rejected the requested value");
    }

    // Apply default precision for both MPFR (C API) and mpfr::mpreal (C++ wrapper).
    mpfr_set_default_prec(s.mpfr_prec);
    mpreal::default_prec = s.mpfr_prec;
}

struct MpfrEnvGuard {
    MpfrEnvSnapshot saved;
    MpfrEnvGuard() : saved(mpfr_env_capture()) {}
    ~MpfrEnvGuard() { mpfr_env_apply(saved); }
};

static REAL unit_roundoff_from_prec(mpfr_prec_t prec, mpfr_rnd_t rnd) {
    // Unit roundoff for base=2: 2^(-prec).
    // Use an exact shift to avoid any platform-dependent exp2() path.
    REAL x = 1.0;
    mpfr_div_2ui(mpfr_ptr(x), mpfr_ptr(x), checked_prec_to_ulong(prec, "prec out of range for mpfr_div_2ui"), rnd);
    return x;
}

static REAL pow2_from_exp(mpfr_exp_t exp2, mpfr_rnd_t rnd) {
    // Return 2^exp2 exactly, if representable under the current exponent range.
    REAL x = 0.0;
    mpfr_clear_flags();
    mpfr_set_ui_2exp(mpfr_ptr(x), 1ul, exp2, rnd);
    if (!mpfr_number_p(mpfr_ptr(x)) || mpfr_overflow_p() || mpfr_underflow_p()) {
        mpfr_test_fail("pow2_from_exp under/overflowed or produced non-number");
    }
    return x;
}

struct LamchExpected {
    REAL E; // eps
    REAL S; // sfmin
    REAL B; // base
    REAL P; // precision
    REAL N; // digits
    REAL R; // rounding
    REAL M; // min exponent
    REAL U; // rmin
    REAL L; // max exponent
    REAL O; // rmax
};

static LamchExpected compute_expected(mpfr_prec_t real_prec, mpfr_exp_t emin, mpfr_exp_t emax, mpfr_rnd_t rnd) {
    const REAL one = 1.0;
    const REAL two = 2.0;

    LamchExpected ex{};
    ex.B = two;

    ex.E = unit_roundoff_from_prec(real_prec, rnd);
    ex.P = ex.E * two;

    // N, M, L are exact integers in these test cases.
    mpfr_set_ui(mpfr_ptr(ex.N), checked_prec_to_ulong(real_prec, "prec out of range for N"), MPFR_RNDN); // Exact integer; rounding irrelevant.
    mpfr_set_si(mpfr_ptr(ex.M), checked_exp_to_long(emin, "emin out of range for M"), MPFR_RNDN);        // Exact integer; rounding irrelevant.
    mpfr_set_si(mpfr_ptr(ex.L), checked_exp_to_long(emax, "emax out of range for L"), MPFR_RNDN);        // Exact integer; rounding irrelevant.

    // DLAMCH('R'): 1.0 if rounding occurs in addition, 0.0 if chopping.
    // Treat MPFR_RNDZ (round toward zero) as chopping; other modes as rounding.
    ex.R = (rnd == MPFR_RNDZ) ? REAL(0.0) : REAL(1.0);

    // U = 2^(emin-1)
    if (emin <= mpfr_get_emin_min()) {
        mpfr_test_fail("emin too small to form 2^(emin-1)");
    }
    ex.U = pow2_from_exp(emin - 1, rnd);

    // O = (1 - E) * 2^emax computed without forming 2^emax as a separate value.
    ex.O = one - ex.E;
    mpfr_clear_flags();
    mpfr_mul_2si(mpfr_ptr(ex.O), mpfr_ptr(ex.O), checked_exp_to_long(emax, "emax out of range for mpfr_mul_2si"), rnd);
    if (!mpfr_number_p(mpfr_ptr(ex.O)) || mpfr_overflow_p() || mpfr_underflow_p()) {
        mpfr_test_fail("rmax computation under/overflowed or produced non-number");
    }

    // S: DLAMCH safe minimum rule.
    // sfmin = max( U, (1/O) * (1 + E) )
    const REAL small = one / ex.O;
    const REAL candidate = small * (one + ex.E);
    ex.S = (candidate >= ex.U) ? candidate : ex.U;

    return ex;
}

static void check_operational_eps_prec(const char *tag, mpfr_rnd_t rnd, const REAL &E, const REAL &P) {
    const REAL one = 1.0;
    const REAL half = 0.5;

    // P should be the spacing (ULP) at 1: nextabove(1) - 1.
    REAL next = one;
    mpfr_nextabove(mpfr_ptr(next));
    const REAL ulp = next - one;
    mpfr_assert_case(ulp == P, tag, "P is not equal to ulp at 1 (nextabove(1)-1)");

    // Operational behavior around 1 depends on the rounding direction.
    // For round-to-nearest, downward, and toward-zero modes: half-ulp does not change 1.
    // For upward (and away-from-zero for positive values): any positive increment rounds to the next value.
    const REAL half_ulp = P * half;
    const REAL a = one + half_ulp;
    const REAL b = one + P;

    if (rnd_rounds_up_for_positive(rnd)) {
        mpfr_assert_case(a == b, tag, "round-up: expected fl(1 + P/2) == 1 + P");
    } else {
        mpfr_assert_case(a == one, tag, "expected fl(1 + P/2) == 1");
    }
    mpfr_assert_case(b > one, tag, "expected fl(1 + P) > 1");

    // Local consistency: P must equal 2*E.
    mpfr_assert_case(E * REAL(2.0) == P, tag, "expected P == 2*E");
}

static void check_range_checks(const char *tag, const REAL &U, const REAL &O, const REAL &S) {
    const REAL zero = 0.0;
    const REAL one = 1.0;
    const REAL two = 2.0;

    mpfr_assert_case(U > zero && U < one, tag, "U (rmin) is not in (0,1)");
    mpfr_assert_case(O > one, tag, "O (rmax) is not > 1");
    mpfr_assert_case(S > zero && S < one, tag, "S (sfmin) is not in (0,1)");

    // Scaling check near zero (avoid being stuck at 0).
    mpfr_assert_case((U * two) > U, tag, "U * 2 <= U (stuck at zero?)");

    // NOTE: In an MPFR model with IEEE-like exponent limits, MPFR has no subnormals.
    // Therefore 1/O may underflow to +0 (especially in RNDN/RNDZ/RNDD) even when O is valid.
    mpfr_assert_case((one / O) >= zero, tag, "1/O is negative or NaN");
    mpfr_assert_case((one / S) > zero, tag, "1/S is not positive");
}

static void check_sfmin_inequalities(const char *tag, const REAL &S, const REAL &U, const REAL &O) {
    const REAL one = 1.0;

    const REAL invS = one / S;
    const REAL invO = one / O;

    // Netlib-style sfmin inequalities (operational checks under the active rounding mode).
    mpfr_assert_case(invS <= O, tag, "1/S > O (violates safe minimum contract)");
    mpfr_assert_case(S >= U, tag, "S < U (violates sfmin >= rmin)");
    mpfr_assert_case(S >= invO, tag, "S < 1/O (violates sfmin >= 1/rmax)");
}

static void check_cross_consistency_rmin_rmax(const char *tag, mpfr_exp_t emin, mpfr_exp_t emax, mpfr_rnd_t rnd, const REAL &E, const REAL &U, const REAL &O) {
    const REAL one = 1.0;

    // In this MPFR model:
    //   U = 2^(emin-1)
    //   O = (1-E)*2^emax
    // Therefore:
    //   U*O = (1-E)*2^(emin+emax-1)
    const REAL got_prod = U * O;

    REAL expected = one - E;
    mpfr_clear_flags();
    mpfr_mul_2si(mpfr_ptr(expected), mpfr_ptr(expected), checked_exp_to_long(emin + emax - 1, "emin+emax out of range for rmin*rmax cross-check"), rnd);
    if (!mpfr_number_p(mpfr_ptr(expected)) || mpfr_overflow_p() || mpfr_underflow_p()) {
        mpfr_test_fail("rmin*rmax cross-check under/overflowed or produced non-number");
    }

    mpfr_assert_case(got_prod == expected, tag, "U*O cross-check failed (inconsistent rmin/rmax model)");
}

static void assert_equal_real(const char *tag, const char *name, const REAL &got, const REAL &expected) {
    // Option 1: Use operator== if available
    if (got == expected)
        return;

    printf("*** Testing Mutils (MPFR) failed: %s mismatch in %s ***\n", tag, name);
    printf("    got      = ");
    printnum(got);
    printf("\n");
    printf("    expected = ");
    printnum(expected);
    printf("\n");
    exit(1);
}

static void check_lamch_mpfr_values(const char *tag, const MpfrEnvSnapshot &cfg, bool print_values) {
    // Sanity: verify we are running under the intended MPFR environment.
    mpfr_assert(mpfr_get_default_prec() == cfg.mpfr_prec, "mpfr default precision did not take effect");
    mpfr_assert(current_real_default_prec() == cfg.real_prec, "REAL default precision did not take effect");
    mpfr_assert(mpfr_get_emin() == cfg.emin, "emin did not take effect");
    mpfr_assert(mpfr_get_emax() == cfg.emax, "emax did not take effect");
    mpfr_assert(mpfr_get_default_rounding_mode() == cfg.mpfr_rnd, "mpfr rounding mode did not take effect");
    mpfr_assert(current_real_default_rnd() == cfg.real_rnd, "REAL rounding mode did not take effect");

    const LamchExpected ex = compute_expected(cfg.real_prec, cfg.emin, cfg.emax, cfg.real_rnd);

    const REAL gotE = Rlamch_mpfr("E");
    const REAL gotS = Rlamch_mpfr("S");
    const REAL gotB = Rlamch_mpfr("B");
    const REAL gotP = Rlamch_mpfr("P");
    const REAL gotN = Rlamch_mpfr("N");
    const REAL gotR = Rlamch_mpfr("R");
    const REAL gotM = Rlamch_mpfr("M");
    const REAL gotU = Rlamch_mpfr("U");
    const REAL gotL = Rlamch_mpfr("L");
    const REAL gotO = Rlamch_mpfr("O");
    const REAL gotZ = Rlamch_mpfr("Z");

    // Re-validate: Rlamch must not mutate the global MPFR environment.
    mpfr_assert(mpfr_get_default_prec() == cfg.mpfr_prec, "Rlamch mutated mpfr default precision");
    mpfr_assert(current_real_default_prec() == cfg.real_prec, "Rlamch mutated REAL default precision");
    mpfr_assert(mpfr_get_emin() == cfg.emin, "Rlamch mutated emin");
    mpfr_assert(mpfr_get_emax() == cfg.emax, "Rlamch mutated emax");
    mpfr_assert(mpfr_get_default_rounding_mode() == cfg.mpfr_rnd, "Rlamch mutated mpfr rounding mode");
    mpfr_assert(current_real_default_rnd() == cfg.real_rnd, "Rlamch mutated REAL rounding mode");

    assert_equal_real(tag, "E", gotE, ex.E);
    assert_equal_real(tag, "B", gotB, ex.B);
    assert_equal_real(tag, "P", gotP, ex.P);
    assert_equal_real(tag, "N", gotN, ex.N);
    assert_equal_real(tag, "R", gotR, ex.R);
    assert_equal_real(tag, "M", gotM, ex.M);
    assert_equal_real(tag, "U", gotU, ex.U);
    assert_equal_real(tag, "L", gotL, ex.L);
    assert_equal_real(tag, "O", gotO, ex.O);
    assert_equal_real(tag, "S", gotS, ex.S);
    mpfr_assert_case(gotZ == REAL(0.0), tag, "Z (dummy) is not 0");

    // Operational property checks
    check_operational_eps_prec(tag, cfg.real_rnd, gotE, gotP);
    check_range_checks(tag, gotU, gotO, gotS);
    check_sfmin_inequalities(tag, gotS, gotU, gotO);
    check_cross_consistency_rmin_rmax(tag, cfg.emin, cfg.emax, cfg.real_rnd, gotE, gotU, gotO);

    if (print_values) {
        printf("Rlamch E: Epsilon                      ");
        printnum(gotE);
        printf("\n");
        printf("Rlamch S: Safe minimum                 ");
        printnum(gotS);
        printf("\n");
        printf("Rlamch B: Base                         ");
        printnum(gotB);
        printf("\n");
        printf("Rlamch P: Precision                    ");
        printnum(gotP);
        printf("\n");
        printf("Rlamch N: Number of digits in mantissa ");
        printnum(gotN);
        printf("\n");
        printf("Rlamch R: Rounding mode                ");
        printnum(gotR);
        printf("\n");
        printf("Rlamch M: Minimum exponent             ");
        printnum(gotM);
        printf("\n");
        printf("Rlamch U: Underflow threshold          ");
        printnum(gotU);
        printf("\n");
        printf("Rlamch L: Largest exponent             ");
        printnum(gotL);
        printf("\n");
        printf("Rlamch O: Overflow threshold           ");
        printnum(gotO);
        printf("\n");
    }
}

static void run_mpfr_env_test(const char *tag, const MpfrEnvSnapshot &cfg, bool print_values) {
    MpfrEnvGuard guard; // Save current environment; restore on scope exit.
    mpfr_env_apply(cfg);
    check_lamch_mpfr_values(tag, cfg, print_values);
}

} // namespace

void Rlamch_mpfr_test() {
#if defined VERBOSE_TEST
    const bool print_values = true;
#else
    const bool print_values = false;
#endif

    struct MpfrEnvCase {
        const char *tag;
        bool use_current;
        mpfr_prec_t prec;
        mpfr_exp_t emin;
        mpfr_exp_t emax;
    };

    struct MpfrRndCase {
        const char *name;
        mpfr_rnd_t rnd;
    };

    // Capture the true current environment once (including REAL defaults).
    const MpfrEnvSnapshot cur = mpfr_env_capture();

    // Test matrix for MPFR environment emulation.
    // "current" starts from the captured process-wide settings.
    static const MpfrEnvCase env_cases[] = {
        {"current", true, 0, 0, 0},

        // Emulations: IEEE-like formats.
        {"binary64", false, (mpfr_prec_t)53, (mpfr_exp_t)-1021, (mpfr_exp_t)1024},
        {"binary80x", false, (mpfr_prec_t)64, (mpfr_exp_t)-16381, (mpfr_exp_t)16384},
        {"binary128", false, (mpfr_prec_t)113, (mpfr_exp_t)-16381, (mpfr_exp_t)16384},

        // Vary precision only (keep binary64 exponent range).
        {"binary64_prec60", false, (mpfr_prec_t)60, (mpfr_exp_t)-1021, (mpfr_exp_t)1024},

        // Vary emax only (force the sfmin branch that depends on (1/O)*(1+E)).
        {"binary64_emax200", false, (mpfr_prec_t)53, (mpfr_exp_t)-1021, (mpfr_exp_t)200},

        // Vary emin only (raise the minimum exponent).
        {"binary64_emin-100", false, (mpfr_prec_t)53, (mpfr_exp_t)-100, (mpfr_exp_t)1024},
    };

    static const MpfrRndCase rnd_cases[] = {
        {"RNDN", MPFR_RNDN}, {"RNDZ", MPFR_RNDZ}, {"RNDU", MPFR_RNDU}, {"RNDD", MPFR_RNDD},
#ifdef MPFR_RNDA
        {"RNDA", MPFR_RNDA},
#endif
    };

    for (const auto &ec: env_cases) {
        MpfrEnvSnapshot base;
        if (ec.use_current) {
            base = cur;
        } else {
            base.mpfr_prec = ec.prec;
            base.real_prec = ec.prec;
            base.emin = ec.emin;
            base.emax = ec.emax;
            // Rounding modes set per rnd_cases below.
            base.mpfr_rnd = cur.mpfr_rnd;
            base.real_rnd = cur.real_rnd;
        }

        for (const auto &rc: rnd_cases) {
            MpfrEnvSnapshot cfg = base;
            cfg.mpfr_rnd = rc.rnd;
            cfg.real_rnd = rc.rnd;

            char tagbuf[128];
            snprintf(tagbuf, sizeof(tagbuf), "%s/%s", ec.tag, rc.name);

            // Print the target environment settings for each case.
            printf("[MPFR] case=%s mpfr_prec=%llu real_prec=%llu emin=%lld emax=%lld rnd=%s\n", tagbuf, (unsigned long long)cfg.mpfr_prec, (unsigned long long)cfg.real_prec, (long long)cfg.emin, (long long)cfg.emax, rc.name);

            run_mpfr_env_test(tagbuf, cfg, print_values);
        }
    }
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
#include <limits>
#include <cstdio>
#include <cstdlib>
// -----------------------------------------------------------------------------
// Helper: choose a conservative exponent limit for GMP mpf.
// This matches the production implementation policy.
// -----------------------------------------------------------------------------
static inline mp_bitcnt_t get_max_safe_exponent() {
    const mp_exp_t max_e = std::numeric_limits<mp_exp_t>::max();
    return static_cast<mp_bitcnt_t>(max_e / 2);
}
// -----------------------------------------------------------------------------
// Rlamch_gmp_test: GMP mpf precision semantics
//
// GMP mpf arithmetic truncates results to the current precision.
// In GMP, the effective precision is rounded up to a limb boundary, and
// mpf_get_prec() reports that effective precision in bits.
//
// With that definition of "prec", the spacing near 1.0 is 2^(-prec).
// Therefore we expect:
//   - Rlamch("R") == 0.0 (no round-to-nearest).
//   - Rlamch("E") == 2^(-prec), so that fl(1 + E/2) == 1 and fl(1 + E) > 1.
//   - Rlamch("P") == base * E (for base=2, P == 2*E).
//
// The exponent-related constants (M, L, U, O, S) are design-defined via
// get_max_safe_exponent() because GMP mpf overflow behavior is not IEEE-like.
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Rlamch_gmp_test: single-case validation (no rounding-mode sweep, no precision sweep)
//
// This follows the same overall structure as the MPFR test:
//   1) Capture the current environment.
//   2) Compute model-expected Lamch values.
//   3) Compare all Lamch outputs against the model.
//   4) Run a small set of operational sanity checks.
// -----------------------------------------------------------------------------

struct GmpEnvSnapshot {
    mp_bitcnt_t default_prec;
};

static inline GmpEnvSnapshot gmp_env_capture() {
    GmpEnvSnapshot s{};
    s.default_prec = mpf_get_default_prec();
    return s;
}
static void gmp_test_fail(const char *tag, const char *what) {
    // Provide enough context to diagnose precision-dependent failures.
    const REAL probe = 1.0;
    const mp_bitcnt_t real_prec = mpf_get_prec(probe.get_mpf_t());
    const mp_bitcnt_t def_prec = mpf_get_default_prec();

    printf("*** Testing Rlamch (GMP) failed: %s (%s) ***\n", tag, what);
    printf("    mpf_default_prec = %lu bits\n", static_cast<unsigned long>(def_prec));
    printf("    REAL prec        = %lu bits\n", static_cast<unsigned long>(real_prec));
    exit(1);
}

static void gmp_assert_case(bool cond, const char *tag, const char *what) {
    if (!cond)
        gmp_test_fail(tag, what);
}

static void gmp_assert_equal_real(const char *tag, const char *name, const REAL &got, const REAL &expected) {
    if (got == expected)
        return;

    printf("*** Testing Rlamch (GMP) failed: %s mismatch in %s ***\n", tag, name);
    printf("    got      = ");
    printnum(got);
    printf("\n");
    printf("    expected = ");
    printnum(expected);
    printf("\n");
    exit(1);
}

struct LamchExpectedGmp {
    REAL E; // "E"
    REAL S; // "S"
    REAL B; // "B"
    REAL P; // "P"
    REAL N; // "N"
    REAL R; // "R"
    REAL M; // "M"
    REAL U; // "U"
    REAL L; // "L"
    REAL O; // "O"
};

static LamchExpectedGmp compute_expected_gmp(mp_bitcnt_t prec_bits, mp_bitcnt_t emax_design) {
    const REAL zero = 0.0;
    const REAL one = 1.0;
    const REAL two = 2.0;

    (void)zero; // silence unused warnings in some build modes

    LamchExpectedGmp ex{};
    ex.B = two;

    // The GMP backend uses truncation (chopping). This test uses the same model
    // as the production implementation:
    //   N = mpf_get_prec(REAL(1))
    //   E = 2^(-N)
    //   P = 2*E
    ex.N = REAL(prec_bits);

    ex.E = one;
    mpf_div_2exp(ex.E.get_mpf_t(), one.get_mpf_t(), prec_bits);
    ex.P = ex.E * two;

    ex.R = REAL(0.0);

    // Exponent range is design-defined for GMP (not IEEE hardware constants).
    ex.L = REAL(emax_design);
    ex.M = one - ex.L;

    // U = 2^(-L)
    ex.U = one;
    mpf_div_2exp(ex.U.get_mpf_t(), one.get_mpf_t(), emax_design);

    // O = (1 - E) * 2^(L)
    ex.O = one - ex.E;
    mpf_mul_2exp(ex.O.get_mpf_t(), ex.O.get_mpf_t(), emax_design);

    // S = max(U, (1/O) * (1 + E))  (DLAMCH-style safe minimum rule)
    const REAL small = one / ex.O;
    const REAL candidate = small * (one + ex.E);
    ex.S = (candidate >= ex.U) ? candidate : ex.U;

    return ex;
}

static void check_lamch_gmp_values(const char *tag, mp_bitcnt_t prec_bits, bool print_values) {
    const REAL zero = 0.0;
    const REAL one = 1.0;
    const REAL two = 2.0;
    const REAL half = 0.5;

    const GmpEnvSnapshot env = gmp_env_capture();

    const mp_bitcnt_t Eexp = get_max_safe_exponent();
    const LamchExpectedGmp ex = compute_expected_gmp(prec_bits, Eexp);

    const REAL gotE = Rlamch_gmp("E");
    const REAL gotS = Rlamch_gmp("S");
    const REAL gotB = Rlamch_gmp("B");
    const REAL gotP = Rlamch_gmp("P");
    const REAL gotN = Rlamch_gmp("N");
    const REAL gotR = Rlamch_gmp("R");
    const REAL gotM = Rlamch_gmp("M");
    const REAL gotU = Rlamch_gmp("U");
    const REAL gotL = Rlamch_gmp("L");
    const REAL gotO = Rlamch_gmp("O");
    const REAL gotZ = Rlamch_gmp("Z");

    // Re-validate: Rlamch must not mutate the GMP global environment.
    gmp_assert_case(mpf_get_default_prec() == env.default_prec, tag, "Rlamch mutated mpf default precision");
    // Exact-value checks (model-defined).
    gmp_assert_equal_real(tag, "B", gotB, ex.B);
    gmp_assert_equal_real(tag, "N", gotN, ex.N);
    gmp_assert_equal_real(tag, "R", gotR, ex.R);

    gmp_assert_equal_real(tag, "E", gotE, ex.E);
    gmp_assert_equal_real(tag, "P", gotP, ex.P);

    gmp_assert_equal_real(tag, "L", gotL, ex.L);
    gmp_assert_equal_real(tag, "M", gotM, ex.M);
    gmp_assert_equal_real(tag, "U", gotU, ex.U);
    gmp_assert_equal_real(tag, "O", gotO, ex.O);
    gmp_assert_equal_real(tag, "S", gotS, ex.S);

    gmp_assert_case(gotZ == zero, tag, "Z (dummy) is not 0");

    // Operational sanity checks.
    gmp_assert_case(gotE > zero, tag, "E is not positive");
    gmp_assert_case(gotU > zero && gotU < one, tag, "U (rmin) is not in (0,1)");
    gmp_assert_case(gotO > one, tag, "O (rmax) is not > 1");
    gmp_assert_case(gotS > zero && gotS < one, tag, "S (sfmin) is not in (0,1)");
    gmp_assert_case((gotU * two) > gotU, tag, "U * 2 <= U (stuck at zero?)");

    // Chopping behavior around 1 (definition of E for this backend).
    // fl(1 + E/2) == 1, and fl(1 + E) > 1.
    const REAL a = one + gotE * half;
    const REAL b = one + gotE;
    gmp_assert_case(a == one, tag, "expected fl(1 + E/2) == 1 (chopping)");
    gmp_assert_case(b > one, tag, "expected fl(1 + E) > 1");

    // Safe minimum contract.
    const REAL invO = one / gotO;
    const REAL invS = one / gotS;
    gmp_assert_case(invS <= gotO, tag, "1/S > O (violates safe minimum contract)");
    gmp_assert_case(gotS >= gotU, tag, "S < U (violates sfmin >= rmin)");
    gmp_assert_case(gotS >= invO, tag, "S < 1/O (violates sfmin >= 1/rmax)");
    gmp_assert_case(invS > zero, tag, "1/S is not positive");
    gmp_assert_case(invO >= zero, tag, "1/O is negative or NaN");

    // Cross-consistency in this GMP model:
    //   O = (1 - E) * 2^(L)  with L == Eexp (design constant).
    // Therefore O * 2^(-L) == 1 - E (use an exact exponent shift).
    REAL shifted = gotO;
    mpf_div_2exp(shifted.get_mpf_t(), shifted.get_mpf_t(), Eexp);
    gmp_assert_case(shifted == (one - gotE), tag, "O * 2^(-L) cross-check failed (expected 1 - E)");

    if (print_values) {
        printf("Rlamch E: Epsilon                      ");
        printnum(gotE);
        printf("\n");
        printf("Rlamch S: Safe minimum                 ");
        printnum(gotS);
        printf("\n");
        printf("Rlamch B: Base                         ");
        printnum(gotB);
        printf("\n");
        printf("Rlamch P: Precision                    ");
        printnum(gotP);
        printf("\n");
        printf("Rlamch N: Number of digits in mantissa ");
        printnum(gotN);
        printf("\n");
        printf("Rlamch R: Rounding mode                ");
        printnum(gotR);
        printf("\n");
        printf("Rlamch M: Minimum exponent             ");
        printnum(gotM);
        printf("\n");
        printf("Rlamch U: Underflow threshold          ");
        printnum(gotU);
        printf("\n");
        printf("Rlamch L: Largest exponent             ");
        printnum(gotL);
        printf("\n");
        printf("Rlamch O: Overflow threshold           ");
        printnum(gotO);
        printf("\n");
    }
}

void Rlamch_gmp_test() {
#if defined VERBOSE_TEST
    const bool print_values = true;
#else
    const bool print_values = false;
#endif

    // Current precision used by REAL (in bits).
    const REAL prec_probe = 1.0;
    const mp_bitcnt_t prec_bits = mpf_get_prec(prec_probe.get_mpf_t());

    check_lamch_gmp_values("current", prec_bits, print_values);
    // Additional precision cases requested: 4096, 128, 64 bits.
    // Note: mpf_set_default_prec() sets the default for newly created mpf values.
    const mp_bitcnt_t requested_precisions[] = {4096, 128, 64};

    for (size_t i = 0; i < (sizeof(requested_precisions) / sizeof(requested_precisions[0])); ++i) {
        printf("[GMP] precision=%ld\n", requested_precisions[i]);
        const mp_bitcnt_t req = requested_precisions[i];
        const mp_bitcnt_t saved = mpf_get_default_prec();
        mpf_set_default_prec(req);

        // Use the effective precision actually used by REAL/mpf objects.
        const REAL probe = 1.0;
        const mp_bitcnt_t eff = mpf_get_prec(probe.get_mpf_t());

        char tagbuf[64];
        snprintf(tagbuf, sizeof(tagbuf), "prec=%llu (eff=%llu)", (unsigned long long)req, (unsigned long long)eff);

        check_lamch_gmp_values(tagbuf, eff, print_values);

        mpf_set_default_prec(saved);
    }
}
#endif // ___MPLAPACK_BUILD_WITH_GMP___

#if defined ___MPLAPACK_BUILD_WITH_QD___
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace {

// -----------------------------------------------------------------------------
// qd_real machine constants
// -----------------------------------------------------------------------------
// This test matches the implementation style used by Rlamch_qd:
//   - E is taken from qd_real::_eps
//   - U and O are taken from qd_real::_min_normalized / qd_real::_max
//   - N is taken from std::numeric_limits<qd_real>::digits
//   - S follows the Netlib-style "safe minimum" logic using U, O, and E
//
// IMPORTANT:
// qd_real is a quad-double type (sum of 4 doubles). Arithmetic does NOT round
// to a fixed p-bit floating-point format, so classic half-ULP tie tests like
// fl(1 + P/2) == 1 do not apply here.
// -----------------------------------------------------------------------------

constexpr int QD_MANTISSA_BITS = std::numeric_limits<qd_real>::digits;
constexpr int QD_EMAX = std::numeric_limits<double>::max_exponent; // expected 1024

[[noreturn]] static void qd_test_fail(const char *what) {
    printf("*** Testing Mutils (qd) failed: %s ***\n", what);
    exit(1);
}

static void qd_assert_case(bool cond, const char *tag, const char *what) {
    if (cond)
        return;
    char buf[256];
    snprintf(buf, sizeof(buf), "%s: %s", tag, what);
    qd_test_fail(buf);
}

static void assert_equal_qd(const char *tag, const char *name, const qd_real &got, const qd_real &expected) {
    if (got == expected)
        return;

    printf("*** Testing Mutils (qd) failed: %s mismatch in %s ***\n", tag, name);
    printf("    got      = ");
    printnum(got);
    printf("\n");
    printf("    expected = ");
    printnum(expected);
    printf("\n");
    exit(1);
}

struct LamchExpectedQD {
    qd_real E;
    qd_real S;
    qd_real B;
    qd_real P;
    qd_real N;
    qd_real R;
    qd_real M;
    qd_real U;
    qd_real L;
    qd_real O;
    qd_real Z;
};

static LamchExpectedQD compute_expected_qd() {
    const qd_real zero(0.0);
    const qd_real one(1.0);
    const qd_real two(2.0);

    LamchExpectedQD ex{};
    ex.B = two;
    ex.E = qd_real::_eps;
    ex.P = ex.E * ex.B;
    ex.N = qd_real(static_cast<double>(QD_MANTISSA_BITS));
    ex.R = one;

    int emin = 0;
    (void)std::frexp(qd_real::_min_normalized, &emin);
    ex.M = qd_real(static_cast<double>(emin));

    ex.U = qd_real::_min_normalized;
    ex.L = qd_real(static_cast<double>(QD_EMAX));
    ex.O = qd_real::_max;

    qd_real sfmin = ex.U;
    const qd_real small = one / ex.O;
    if (small >= sfmin) {
        sfmin = small * (one + ex.E);
    }
    ex.S = sfmin;
    ex.Z = zero;
    return ex;
}

static void check_operational_eps_prec_qd(const char *tag, const qd_real &E, const qd_real &P) {
    const qd_real one(1.0);
    const qd_real two(2.0);
    qd_assert_case(E * two == P, tag, "expected P == 2*E (base 2)");
    qd_assert_case(one + E > one, tag, "expected 1 + E > 1");
    qd_assert_case(one + P > one, tag, "expected 1 + P > 1");
}

static void check_range_checks_qd(const char *tag, const qd_real &U, const qd_real &O, const qd_real &S) {
    const qd_real zero(0.0);
    const qd_real one(1.0);
    const qd_real two(2.0);

    qd_assert_case(U > zero && U < one, tag, "U (rmin) is not in (0,1)");
    qd_assert_case(O > one, tag, "O (rmax) is not > 1");
    qd_assert_case(S > zero && S < one, tag, "S (sfmin) is not in (0,1)");
    qd_assert_case((U * two) > U, tag, "U * 2 <= U (stuck at zero?)");

    const double invS = to_double(one / S);
    qd_assert_case(invS > 0.0, tag, "1/S is not positive");
    qd_assert_case(std::isfinite(invS), tag, "1/S is not finite");
}

static void check_sfmin_inequalities_qd(const char *tag, const qd_real &S, const qd_real &U, const qd_real &O) {
    const qd_real one(1.0);
    const qd_real invS = one / S;
    const qd_real invO = one / O;
    qd_assert_case(invS <= O, tag, "1/S > O (violates safe minimum contract)");
    qd_assert_case(S >= U, tag, "S < U (violates sfmin >= rmin)");
    qd_assert_case(S >= invO, tag, "S < 1/O (violates sfmin >= 1/rmax)");
}

static long qd_to_long_checked(const qd_real &x, const char *tag, const char *msg) {
    const double d = to_double(x);
    const long v = static_cast<long>(d);
    qd_assert_case(qd_real(static_cast<double>(v)) == x, tag, msg);
    return v;
}

static void check_cross_consistency_rmin_rmax_qd(const char *tag, const qd_real &M, const qd_real &U, const qd_real &L, const qd_real &O) {
    const qd_real two(2.0);
    qd_assert_case(U > qd_real(0.0), tag, "U is not positive");
    qd_assert_case(U / two < U, tag, "U/2 >= U (not smaller)");
    qd_assert_case(std::isfinite(to_double(O)), tag, "O converts to non-finite double");

    const long m = qd_to_long_checked(M, tag, "M is not an exact integer");
    const long l = qd_to_long_checked(L, tag, "L is not an exact integer");

    const double dU = to_double(U);
    qd_assert_case(dU > 0.0, tag, "U converts to non-positive double");
    int expU = 0;
    (void)std::frexp(dU, &expU);
    qd_assert_case(expU == m, tag, "frexp exponent of U inconsistent with M");

    const double dO = to_double(O);
    qd_assert_case(std::isfinite(dO) && dO > 0.0, tag, "O converts to non-finite/non-positive double");
    int expO = 0;
    (void)std::frexp(dO, &expO);
    qd_assert_case(expO == l, tag, "frexp exponent of O inconsistent with L");
}

static void check_lamch_qd_values(const char *tag, bool print_values) {
    const LamchExpectedQD ex = compute_expected_qd();

    const qd_real gotE = Rlamch_qd("E");
    const qd_real gotS = Rlamch_qd("S");
    const qd_real gotB = Rlamch_qd("B");
    const qd_real gotP = Rlamch_qd("P");
    const qd_real gotN = Rlamch_qd("N");
    const qd_real gotR = Rlamch_qd("R");
    const qd_real gotM = Rlamch_qd("M");
    const qd_real gotU = Rlamch_qd("U");
    const qd_real gotL = Rlamch_qd("L");
    const qd_real gotO = Rlamch_qd("O");
    const qd_real gotZ = Rlamch_qd("Z");

    assert_equal_qd(tag, "E", gotE, ex.E);
    assert_equal_qd(tag, "B", gotB, ex.B);
    assert_equal_qd(tag, "P", gotP, ex.P);
    assert_equal_qd(tag, "N", gotN, ex.N);
    assert_equal_qd(tag, "R", gotR, ex.R);
    assert_equal_qd(tag, "M", gotM, ex.M);
    assert_equal_qd(tag, "U", gotU, ex.U);
    assert_equal_qd(tag, "L", gotL, ex.L);
    assert_equal_qd(tag, "O", gotO, ex.O);
    assert_equal_qd(tag, "S", gotS, ex.S);
    qd_assert_case(gotZ == ex.Z, tag, "Z (dummy) is not 0");

    check_operational_eps_prec_qd(tag, gotE, gotP);
    check_range_checks_qd(tag, gotU, gotO, gotS);
    check_sfmin_inequalities_qd(tag, gotS, gotU, gotO);
    check_cross_consistency_rmin_rmax_qd(tag, gotM, gotU, gotL, gotO);

    if (print_values) {
        printf("Rlamch E: Epsilon                      ");
        printnum(gotE);
        printf("\n");
        printf("Rlamch S: Safe minimum                 ");
        printnum(gotS);
        printf("\n");
        printf("Rlamch B: Base                         ");
        printnum(gotB);
        printf("\n");
        printf("Rlamch P: Precision                    ");
        printnum(gotP);
        printf("\n");
        printf("Rlamch N: Number of digits in mantissa ");
        printnum(gotN);
        printf("\n");
        printf("Rlamch R: Rounding mode                ");
        printnum(gotR);
        printf("\n");
        printf("Rlamch M: Minimum exponent             ");
        printnum(gotM);
        printf("\n");
        printf("Rlamch U: Underflow threshold          ");
        printnum(gotU);
        printf("\n");
        printf("Rlamch L: Largest exponent             ");
        printnum(gotL);
        printf("\n");
        printf("Rlamch O: Overflow threshold           ");
        printnum(gotO);
        printf("\n");
    }
}

} // namespace

void Rlamch_qd_test() {
#if defined VERBOSE_TEST
    const bool print_values = true;
#else
    const bool print_values = false;
#endif

    const char *tag = "qd_real";
    check_lamch_qd_values(tag, print_values);
}

#endif // ___MPLAPACK_BUILD_WITH_QD___

#if defined ___MPLAPACK_BUILD_WITH_DD___

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace {

// -----------------------------------------------------------------------------
// dd_real machine constants
// -----------------------------------------------------------------------------
// This test matches the implementation style shown in RlamchE_dd/RlamchS_dd:
//   - E is taken from dd_real::_eps
//   - U and O are taken from dd_real::_min_normalized / dd_real::_max
//   - N is taken from std::numeric_limits<dd_real>::digits
// - The exponent range matches IEEE-754 binary64 (double)
//
// IMPORTANT:
// Do NOT derive O via std::ldexp(1.0, 1024) (it overflows in double).
// Use dd_real::_max (as the implementation does).
// -----------------------------------------------------------------------------

constexpr int DD_MANTISSA_BITS = std::numeric_limits<dd_real>::digits;
constexpr int DD_EMIN = std::numeric_limits<double>::min_exponent; // -1021
constexpr int DD_EMAX = std::numeric_limits<double>::max_exponent; // 1024

// Fail/assert utilities
[[noreturn]] static void dd_test_fail(const char *what) {
    printf("*** Testing Mutils (dd) failed: %s ***\n", what);
    exit(1);
}

static void dd_assert(bool cond, const char *what) {
    if (!cond)
        dd_test_fail(what);
}

static void dd_assert_case(bool cond, const char *tag, const char *what) {
    if (cond)
        return;
    char buf[256];
    snprintf(buf, sizeof(buf), "%s: %s", tag, what);
    dd_test_fail(buf);
}

static void assert_equal_dd(const char *tag, const char *name, const dd_real &got, const dd_real &expected) {
    if (got == expected)
        return;

    printf("*** Testing Mutils (dd) failed: %s mismatch in %s ***\n", tag, name);
    printf("    got      = ");
    printnum(got);
    printf("\n");
    printf("    expected = ");
    printnum(expected);
    printf("\n");
    exit(1);
}

// Expected values structure
struct LamchExpectedDD {
    dd_real E; // eps (unit roundoff)
    dd_real S; // sfmin (safe minimum)
    dd_real B; // base
    dd_real P; // precision (ulp at 1)
    dd_real N; // digits (mantissa bits)
    dd_real R; // rounding mode indicator
    dd_real M; // min exponent
    dd_real U; // rmin (underflow threshold)
    dd_real L; // max exponent
    dd_real O; // rmax (overflow threshold)
    dd_real Z; // dummy (always 0)
};

// Compute expected values for dd_real
static LamchExpectedDD compute_expected_dd() {
    const dd_real zero(0.0);
    const dd_real one(1.0);
    const dd_real two(2.0);

    LamchExpectedDD ex{};

    // Match Rlamch_dd implementation.
    ex.B = two;
    ex.E = dd_real::_eps;
    ex.P = ex.E * ex.B;
    ex.N = dd_real(static_cast<double>(DD_MANTISSA_BITS));
    ex.R = one;

    // M: replicate RlamchM_dd (frexp exponent of the minimum normalized value).
    int emin = 0;
    (void)std::frexp(dd_real::_min_normalized, &emin);
    ex.M = dd_real(static_cast<double>(emin));

    // U/O/L: direct constants used by RlamchU_dd/RlamchO_dd/RlamchL_dd.
    ex.U = dd_real::_min_normalized;
    ex.L = dd_real(static_cast<double>(DD_EMAX));
    ex.O = dd_real::_max;

    // S: replicate RlamchS_dd (Netlib-style safe minimum logic).
    dd_real sfmin = ex.U;
    const dd_real small = one / ex.O;
    if (small >= sfmin) {
        sfmin = small * (one + ex.E);
    }
    ex.S = sfmin;

    // Z: Dummy (always 0)
    ex.Z = zero;

    return ex;
}

// Operational checks
static void check_operational_eps_prec_dd(const char *tag, const dd_real &E, const dd_real &P) {
    const dd_real one(1.0);
    const dd_real two(2.0);

    // dd_real is a double-double type. Addition does NOT round to a fixed p-bit
    // floating-point format, so the classic half-ULP tie test (fl(1 + P/2) == 1)
    // does not apply.
    //
    // Here we only test invariants consistent with the implementation:
    //   P == 2*E (base 2), and E matches the nominal precision implied by digits.
    dd_assert_case(E * two == P, tag, "expected P == 2*E (base 2)");

    // Sanity: these should strictly increase (dd_real can represent tiny increments near 1).
    dd_assert_case(one + E > one, tag, "expected 1 + E > 1 (dd_real has no half-ULP tie behavior)");
    dd_assert_case(one + P > one, tag, "expected 1 + P > 1");

    // Do not assume E == 2^(-digits) for double-double types.
    // Many dd_real implementations define digits and epsilon-like constants independently.
    // Exact E is already validated elsewhere against the implementation's constant.
}

static void check_range_checks_dd(const char *tag, const dd_real &U, const dd_real &O, const dd_real &S) {
    const dd_real zero(0.0);
    const dd_real one(1.0);
    const dd_real two(2.0);

    dd_assert_case(U > zero && U < one, tag, "U (rmin) is not in (0,1)");
    dd_assert_case(O > one, tag, "O (rmax) is not > 1");
    dd_assert_case(S > zero && S < one, tag, "S (sfmin) is not in (0,1)");

    // Scaling check near zero (avoid being stuck at 0)
    dd_assert_case((U * two) > U, tag, "U * 2 <= U (stuck at zero?)");

    // Reciprocal checks
    dd_assert_case((one / O) >= zero, tag, "1/O is negative or NaN");
    dd_assert_case((one / S) > zero, tag, "1/S is not positive");
    dd_assert_case(isfinite(one / S), tag, "1/S is not finite");
}

static void check_sfmin_inequalities_dd(const char *tag, const dd_real &S, const dd_real &U, const dd_real &O) {
    const dd_real one(1.0);

    const dd_real invS = one / S;
    const dd_real invO = one / O;

    // Netlib-style sfmin inequalities
    dd_assert_case(invS <= O, tag, "1/S > O (violates safe minimum contract)");
    dd_assert_case(S >= U, tag, "S < U (violates sfmin >= rmin)");
    dd_assert_case(S >= invO, tag, "S < 1/O (violates sfmin >= 1/rmax)");
}

static long dd_to_long_checked(const dd_real &x, const char *what) {
    const double d = to_double(x);
    const long v = static_cast<long>(d);
    // Ensure x is exactly an integer representable as double (expected for M/L in this test)
    dd_assert_case(dd_real(static_cast<double>(v)) == x, what, "expected an exact integer value");
    return v;
}

static void check_cross_consistency_rmin_rmax_dd(const char *tag, const dd_real &M, const dd_real &U, const dd_real &L, const dd_real &O) {
    const dd_real two(2.0);

    dd_assert_case(U > dd_real(0.0), tag, "U is not positive");
    dd_assert_case(U / two < U, tag, "U/2 >= U (not smaller)");
    dd_assert_case(isfinite(O), tag, "O is not finite");

    // Match RlamchM_dd / RlamchL_dd semantics:
    // They are derived from frexp() / known exponent constants in the implementation.
    const long m = dd_to_long_checked(M, "M");
    const long l = dd_to_long_checked(L, "L");

    const double dU = to_double(U);
    dd_assert_case(dU > 0.0, tag, "U converts to non-positive double");
    int expU = 0;
    (void)std::frexp(dU, &expU);
    dd_assert_case(expU == m, tag, "frexp exponent of U inconsistent with M");

    const double dO = to_double(O);
    dd_assert_case(std::isfinite(dO) && dO > 0.0, tag, "O converts to non-finite/non-positive double");
    int expO = 0;
    (void)std::frexp(dO, &expO);
    dd_assert_case(expO == l, tag, "frexp exponent of O inconsistent with L");
}

static void check_lamch_dd_values(const char *tag, bool print_values) {
    const LamchExpectedDD ex = compute_expected_dd();

    // Fetch actual values from Rlamch_dd
    const dd_real gotE = Rlamch_dd("E");
    const dd_real gotS = Rlamch_dd("S");
    const dd_real gotB = Rlamch_dd("B");
    const dd_real gotP = Rlamch_dd("P");
    const dd_real gotN = Rlamch_dd("N");
    const dd_real gotR = Rlamch_dd("R");
    const dd_real gotM = Rlamch_dd("M");
    const dd_real gotU = Rlamch_dd("U");
    const dd_real gotL = Rlamch_dd("L");
    const dd_real gotO = Rlamch_dd("O");
    const dd_real gotZ = Rlamch_dd("Z");

    // Exact-value checks (match implementation-defined constants)
    assert_equal_dd(tag, "E", gotE, ex.E);
    assert_equal_dd(tag, "B", gotB, ex.B);
    assert_equal_dd(tag, "P", gotP, ex.P);
    assert_equal_dd(tag, "N", gotN, ex.N);
    assert_equal_dd(tag, "R", gotR, ex.R);
    assert_equal_dd(tag, "M", gotM, ex.M);
    assert_equal_dd(tag, "U", gotU, ex.U);
    assert_equal_dd(tag, "L", gotL, ex.L);
    assert_equal_dd(tag, "O", gotO, ex.O);
    assert_equal_dd(tag, "S", gotS, ex.S);
    dd_assert_case(gotZ == ex.Z, tag, "Z (dummy) is not 0");

    // Operational property checks
    check_operational_eps_prec_dd(tag, gotE, gotP);
    check_range_checks_dd(tag, gotU, gotO, gotS);
    check_sfmin_inequalities_dd(tag, gotS, gotU, gotO);
    check_cross_consistency_rmin_rmax_dd(tag, gotM, gotU, gotL, gotO);

    if (print_values) {
        printf("Rlamch E: Epsilon                      ");
        printnum(gotE);
        printf("\n");
        printf("Rlamch S: Safe minimum                 ");
        printnum(gotS);
        printf("\n");
        printf("Rlamch B: Base                         ");
        printnum(gotB);
        printf("\n");
        printf("Rlamch P: Precision                    ");
        printnum(gotP);
        printf("\n");
        printf("Rlamch N: Number of digits in mantissa ");
        printnum(gotN);
        printf("\n");
        printf("Rlamch R: Rounding mode                ");
        printnum(gotR);
        printf("\n");
        printf("Rlamch M: Minimum exponent             ");
        printnum(gotM);
        printf("\n");
        printf("Rlamch U: Underflow threshold          ");
        printnum(gotU);
        printf("\n");
        printf("Rlamch L: Largest exponent             ");
        printnum(gotL);
        printf("\n");
        printf("Rlamch O: Overflow threshold           ");
        printnum(gotO);
        printf("\n");
    }
}

} // namespace

void Rlamch_dd_test() {
#if defined VERBOSE_TEST
    const bool print_values = true;
#else
    const bool print_values = false;
#endif

    const char *tag = "dd_real";
    check_lamch_dd_values(tag, print_values);
}

#endif // ___MPLAPACK_BUILD_WITH_DD___

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace {
[[noreturn]] static void double_test_fail(const char *what) {
    printf("*** Testing Mutils (double) failed: %s ***\n", what);
    exit(1);
}

static void double_assert(bool cond, const char *what) {
    if (!cond)
        double_test_fail(what);
}

static void double_assert_case(bool cond, const char *tag, const char *what) {
    if (cond)
        return;
    char buf[256];
    snprintf(buf, sizeof(buf), "%s: %s", tag, what);
    double_test_fail(buf);
}

static void assert_equal_double(const char *tag, const char *name, double got, double expected) {
    if (got == expected)
        return;

    printf("*** Testing Mutils (double) failed: %s mismatch in %s ***\n", tag, name);
    printf("    got      = ");
    printnum(got);
    printf("\n");
    printf("    expected = ");
    printnum(expected);
    printf("\n");
    exit(1);
}

// Expected values structure
struct LamchExpectedDouble {
    double E; // eps (unit roundoff)
    double S; // sfmin (safe minimum)
    double B; // base
    double P; // precision (ulp at 1)
    double N; // digits (mantissa bits)
    double R; // rounding mode indicator
    double M; // min exponent
    double U; // rmin (underflow threshold)
    double L; // max exponent
    double O; // rmax (overflow threshold)
    double Z; // dummy (always 0)
};

// Compute expected values from std::numeric_limits<double>
static LamchExpectedDouble compute_expected_double() {
    const double zero = 0.0;
    const double one = 1.0;
    const double two = 2.0;

    LamchExpectedDouble ex{};

    // B: Base (radix)
    ex.B = static_cast<double>(std::numeric_limits<double>::radix);

    // N: Number of digits in mantissa
    ex.N = static_cast<double>(std::numeric_limits<double>::digits);

    // E: Unit roundoff = epsilon / 2
    // std::numeric_limits<double>::epsilon() returns ulp(1), i.e., nextafter(1,2) - 1
    // LAPACK convention: E = epsilon / 2 (half an ulp)
    const double ulp1 = std::numeric_limits<double>::epsilon();
    ex.E = ulp1 * 0.5;

    // P: Precision = epsilon (full ulp at 1)
    ex.P = ulp1;

    // M: Minimum exponent (DLAMCH style: such that base^(M-1) = rmin)
    ex.M = static_cast<double>(std::numeric_limits<double>::min_exponent);

    // L: Maximum exponent (DLAMCH style: such that (1-E)*base^L = rmax)
    ex.L = static_cast<double>(std::numeric_limits<double>::max_exponent);

    // U: Underflow threshold = smallest normalized positive number
    ex.U = std::numeric_limits<double>::min();

    // O: Overflow threshold = largest finite number
    ex.O = std::numeric_limits<double>::max();

    // S: Safe minimum as in Netlib DLAMCH
    //    sfmin = max(rmin, (1/rmax)*(1+eps))
    const double small = one / ex.O;
    const double candidate = small * (one + ex.E);
    ex.S = (candidate >= ex.U) ? candidate : ex.U;

    // R: Rounding indicator
    //    1.0 if rounding occurs in addition (IEEE round-to-nearest)
    //    0.0 if truncation/chopping
    ex.R = one;

    // Z: Dummy (always 0)
    ex.Z = zero;

    return ex;
}

// Operational checks
static void check_operational_eps_prec_double(const char *tag, double E, double P) {
    const double one = 1.0;
    const double half = 0.5;

    // P should be ulp(1) = nextafter(1, 2) - 1
    const double next = std::nextafter(one, 2.0);
    const double ulp = next - one;
    double_assert_case(ulp == P, tag, "P is not equal to ulp at 1 (nextafter(1,2)-1)");

    // Operational behavior: IEEE round-to-nearest
    // fl(1 + P/2) == 1 (half-ulp rounds down)
    // fl(1 + P) > 1
    const double half_ulp = P * half;
    const double a = one + half_ulp;
    const double b = one + P;

    double_assert_case(a == one, tag, "expected fl(1 + P/2) == 1 (round-to-nearest)");
    double_assert_case(b > one, tag, "expected fl(1 + P) > 1");

    // Local consistency: P must equal 2*E
    double_assert_case(E * 2.0 == P, tag, "expected P == 2*E");
}

static void check_range_checks_double(const char *tag, double U, double O, double S) {
    const double zero = 0.0;
    const double one = 1.0;
    const double two = 2.0;

    double_assert_case(U > zero && U < one, tag, "U (rmin) is not in (0,1)");
    double_assert_case(O > one, tag, "O (rmax) is not > 1");
    double_assert_case(S > zero && S < one, tag, "S (sfmin) is not in (0,1)");

    // Scaling check near zero (avoid being stuck at 0)
    double_assert_case((U * two) > U, tag, "U * 2 <= U (stuck at zero?)");

    // Reciprocal checks
    double_assert_case((one / O) >= zero, tag, "1/O is negative or NaN");
    double_assert_case((one / S) > zero, tag, "1/S is not positive");
    double_assert_case(std::isfinite(one / S), tag, "1/S is not finite");
}

static void check_sfmin_inequalities_double(const char *tag, double S, double U, double O) {
    const double one = 1.0;

    const double invS = one / S;
    const double invO = one / O;

    // Netlib-style sfmin inequalities
    double_assert_case(invS <= O, tag, "1/S > O (violates safe minimum contract)");
    double_assert_case(S >= U, tag, "S < U (violates sfmin >= rmin)");
    double_assert_case(S >= invO, tag, "S < 1/O (violates sfmin >= 1/rmax)");
}

static void check_cross_consistency_rmin_rmax_double(const char *tag, double E, double U, double O) {
    // For IEEE-754 binary64:
    //   U = 2^(emin-1) = 2^(-1022)
    //   O = (1 - 2^(-53)) * 2^1024
    //
    // Cross-check: verify U and O are consistent with their definitions
    // U * 2 should not overflow, and O / 2 should not underflow

    const double two = 2.0;

    // Verify U is the smallest normalized positive
    double_assert_case(U > 0.0, tag, "U is not positive");
    double_assert_case(U / two < U, tag, "U/2 >= U (not smallest normalized)");

    // Verify O is the largest finite
    double_assert_case(std::isfinite(O), tag, "O is not finite");
    double_assert_case(!std::isfinite(O * two), tag, "O*2 is finite (not largest)");

    // Verify exponent relationship via log2
    // For binary64: U = 2^(-1022), O  2^1024
    // log2(U) should be close to M-1 = -1022
    // log2(O) should be close to L = 1024
    const double log2U = std::log2(U);
    const double log2O = std::log2(O);

    double_assert_case(std::abs(log2U - (-1022.0)) < 1.0, tag, "log2(U) inconsistent with expected min exponent");
    double_assert_case(std::abs(log2O - 1024.0) < 1.0, tag, "log2(O) inconsistent with expected max exponent");
}

static void check_lamch_double_values(const char *tag, bool print_values) {
    const LamchExpectedDouble ex = compute_expected_double();

    // Fetch actual values from Rlamch_double
    const double gotE = Rlamch_double("E");
    const double gotS = Rlamch_double("S");
    const double gotB = Rlamch_double("B");
    const double gotP = Rlamch_double("P");
    const double gotN = Rlamch_double("N");
    const double gotR = Rlamch_double("R");
    const double gotM = Rlamch_double("M");
    const double gotU = Rlamch_double("U");
    const double gotL = Rlamch_double("L");
    const double gotO = Rlamch_double("O");
    const double gotZ = Rlamch_double("Z");

    // Exact-value checks
    assert_equal_double(tag, "E", gotE, ex.E);
    assert_equal_double(tag, "B", gotB, ex.B);
    assert_equal_double(tag, "P", gotP, ex.P);
    assert_equal_double(tag, "N", gotN, ex.N);
    assert_equal_double(tag, "R", gotR, ex.R);
    assert_equal_double(tag, "M", gotM, ex.M);
    assert_equal_double(tag, "U", gotU, ex.U);
    assert_equal_double(tag, "L", gotL, ex.L);
    assert_equal_double(tag, "O", gotO, ex.O);
    assert_equal_double(tag, "S", gotS, ex.S);
    double_assert_case(gotZ == ex.Z, tag, "Z (dummy) is not 0");

    // Operational property checks
    check_operational_eps_prec_double(tag, gotE, gotP);
    check_range_checks_double(tag, gotU, gotO, gotS);
    check_sfmin_inequalities_double(tag, gotS, gotU, gotO);
    check_cross_consistency_rmin_rmax_double(tag, gotE, gotU, gotO);

    if (print_values) {
        printf("Rlamch E: Epsilon                      ");
        printnum(gotE);
        printf("\n");
        printf("Rlamch S: Safe minimum                 ");
        printnum(gotS);
        printf("\n");
        printf("Rlamch B: Base                         ");
        printnum(gotB);
        printf("\n");
        printf("Rlamch P: Precision                    ");
        printnum(gotP);
        printf("\n");
        printf("Rlamch N: Number of digits in mantissa ");
        printnum(gotN);
        printf("\n");
        printf("Rlamch R: Rounding mode                ");
        printnum(gotR);
        printf("\n");
        printf("Rlamch M: Minimum exponent             ");
        printnum(gotM);
        printf("\n");
        printf("Rlamch U: Underflow threshold          ");
        printnum(gotU);
        printf("\n");
        printf("Rlamch L: Largest exponent             ");
        printnum(gotL);
        printf("\n");
        printf("Rlamch O: Overflow threshold           ");
        printnum(gotO);
        printf("\n");
    }
}

} // namespace

void Rlamch_double_test() {
#if defined VERBOSE_TEST
    const bool print_values = true;
#else
    const bool print_values = false;
#endif

    const char *tag = "binary64";
    check_lamch_double_values(tag, print_values);
}

#endif // ___MPLAPACK_BUILD_WITH_DOUBLE___

#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___
#include <cfenv> // std::fegetround, FE_TONEAREST

void Rlamch__Float128_test() {
#if defined VERBOSE_TEST
    const bool print_values = true;
#else
    const bool print_values = false;
#endif

    const char *tag = "binary128";

    auto fail = [&](const char *what) {
        printf("*** Testing Mutils (_Float128) failed: %s ***\n", what);
        exit(1);
    };

    auto assert_case = [&](bool cond, const char *what) {
        if (cond)
            return;
        char buf[256];
        snprintf(buf, sizeof(buf), "%s: %s", tag, what);
        fail(buf);
    };

    auto assert_equal = [&](const char *name, _Float128 got, _Float128 expected) {
        if (got == expected)
            return;

        printf("*** Testing Mutils (_Float128) failed: %s mismatch in %s ***\n", tag, name);
        printf("    got      = ");
        printnum(got);
        printf("\n");
        printf("    expected = ");
        printnum(expected);
        printf("\n");
        exit(1);
    };

    // libm entry points for binary128 helpers.
    _Float128 ldexpf128(_Float128, int);
    _Float128 nextafterf128(_Float128, _Float128);

    const _Float128 zero = (_Float128)0.0;
    const _Float128 one = (_Float128)1.0;
    const _Float128 two = (_Float128)2.0;

    // Expected values for IEEE-754 binary128 based on compiler-provided parameters.
    const int p = (int)__FLT128_MANT_DIG__;
    const int emin = (int)__FLT128_MIN_EXP__;
    const int emax = (int)__FLT128_MAX_EXP__;

    const _Float128 exB = two;
    const _Float128 exN = (_Float128)p;

    // ulp(1) = 2^(1-p), unit roundoff E = ulp(1)/2.
    const _Float128 exP = ldexpf128(one, 1 - p);
    const _Float128 exE = exP / two;

    // DLAMCH-style exponents: rmin = 2^(emin-1), rmax = (1-E)*2^emax.
    const _Float128 exM = (_Float128)emin;
    const _Float128 exL = (_Float128)emax;
    const _Float128 exU = ldexpf128(one, emin - 1);
    const _Float128 exO = ldexpf128(one - ldexpf128(one, -p), emax);

    // Safe minimum: max(rmin, (1/rmax)*(1+E)).
    const _Float128 small = one / exO;
    const _Float128 candidate = small * (one + exE);
    const _Float128 exS = (candidate >= exU) ? candidate : exU;

    const _Float128 exR = one;
    const _Float128 exZ = zero;

    // Fetch actual values from Rlamch__Float128.
    const _Float128 gotE = Rlamch__Float128("E");
    const _Float128 gotS = Rlamch__Float128("S");
    const _Float128 gotB = Rlamch__Float128("B");
    const _Float128 gotP = Rlamch__Float128("P");
    const _Float128 gotN = Rlamch__Float128("N");
    const _Float128 gotR = Rlamch__Float128("R");
    const _Float128 gotM = Rlamch__Float128("M");
    const _Float128 gotU = Rlamch__Float128("U");
    const _Float128 gotL = Rlamch__Float128("L");
    const _Float128 gotO = Rlamch__Float128("O");
    const _Float128 gotZ = Rlamch__Float128("Z");

    // Exact-value checks.
    assert_equal("E", gotE, exE);
    assert_equal("B", gotB, exB);
    assert_equal("P", gotP, exP);
    assert_equal("N", gotN, exN);
    assert_equal("R", gotR, exR);
    assert_equal("M", gotM, exM);
    assert_equal("U", gotU, exU);
    assert_equal("L", gotL, exL);
    assert_equal("O", gotO, exO);
    assert_equal("S", gotS, exS);
    assert_case(gotZ == exZ, "Z (dummy) is not 0");

    // Operational property checks.
    const _Float128 next = nextafterf128(one, two);
    assert_case((next - one) == gotP, "P is not equal to ulp(1) = nextafter(1,2)-1");
    assert_case((gotE * two) == gotP, "expected P == 2*E");

    // Strong, rounding-mode-independent check for 1 + P
    volatile _Float128 b = one + gotP;
    assert_case(b == next, "expected fl(1 + P) == nextafter(1,2)");
    assert_case(b > one, "expected fl(1 + 2E) > 1"); // this is effectivly covered by b==next.

    // Tie check: only assert when FE_TONEAREST
    if (std::fegetround() == FE_TONEAREST) {
        volatile _Float128 a = one + gotE;
        assert_case(a == one, "expected fl(1 + E) == 1 under ties-to-even");
    }

    // Range and reciprocal sanity checks.
    assert_case(gotU > zero && gotU < one, "U (rmin) is not in (0,1)");
    assert_case(gotO > one, "O (rmax) is not > 1");
    assert_case(gotS > zero && gotS < one, "S (sfmin) is not in (0,1)");
    assert_case((gotU * two) > gotU, "U*2 <= U (stuck at zero?)");
    assert_case((one / gotS) > zero, "1/S is not positive");
    assert_case(__builtin_isfinite(one / gotS), "1/S is not finite");

    // Netlib-style sfmin inequalities.
    assert_case((one / gotS) <= gotO, "1/S > O (violates safe minimum contract)");
    assert_case(gotS >= gotU, "S < U (violates sfmin >= rmin)");
    assert_case(gotS >= (one / gotO), "S < 1/O (violates sfmin >= 1/rmax)");

    // Cross-check: U*O == (1-E)*2^(emin+emax-1).
    const _Float128 got_prod = gotU * gotO;
    const _Float128 expected_prod = ldexpf128(one - gotE, emin + emax - 1);
    assert_case(got_prod == expected_prod, "U*O cross-check failed (inconsistent rmin/rmax model)");

    if (print_values) {
        printf("Rlamch E: Epsilon                      ");
        printnum(gotE);
        printf("\n");
        printf("Rlamch S: Safe minimum                 ");
        printnum(gotS);
        printf("\n");
        printf("Rlamch B: Base                         ");
        printnum(gotB);
        printf("\n");
        printf("Rlamch P: Precision                    ");
        printnum(gotP);
        printf("\n");
        printf("Rlamch N: Number of digits in mantissa ");
        printnum(gotN);
        printf("\n");
        printf("Rlamch R: Rounding mode                ");
        printnum(gotR);
        printf("\n");
        printf("Rlamch M: Minimum exponent             ");
        printnum(gotM);
        printf("\n");
        printf("Rlamch U: Underflow threshold          ");
        printnum(gotU);
        printf("\n");
        printf("Rlamch L: Largest exponent             ");
        printnum(gotL);
        printf("\n");
        printf("Rlamch O: Overflow threshold           ");
        printnum(gotO);
        printf("\n");
    }
}
#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
#include <cfenv> // std::fegetround, FE_TONEAREST

void Rlamch__Float64x_test() {
#if defined VERBOSE_TEST
    const bool print_values = true;
#else
    const bool print_values = false;
#endif

    const char *tag = "Float64x";

    auto fail = [&](const char *what) {
        printf("*** Testing Mutils (_Float64x) failed: %s ***\n", what);
        exit(1);
    };

    auto assert_case = [&](bool cond, const char *what) {
        if (cond)
            return;
        char buf[256];
        snprintf(buf, sizeof(buf), "%s: %s", tag, what);
        fail(buf);
    };

    auto assert_equal = [&](const char *name, _Float64x got, _Float64x expected) {
        if (got == expected)
            return;

        printf("*** Testing Mutils (_Float64x) failed: %s mismatch in %s ***\n", tag, name);
        printf("    got      = ");
        printnum(got);
        printf("\n");
        printf("    expected = ");
        printnum(expected);
        printf("\n");
        exit(1);
    };

    const _Float64x zero = (_Float64x)0.0;
    const _Float64x one = (_Float64x)1.0;
    const _Float64x two = (_Float64x)2.0;

    // Expected values for _Float64x based on compiler-provided parameters.
    const int p = (int)__FLT64X_MANT_DIG__;
    const int emin = (int)__FLT64X_MIN_EXP__;
    const int emax = (int)__FLT64X_MAX_EXP__;

    const _Float64x exB = two;
    const _Float64x exN = (_Float64x)p;

    // ulp(1) = 2^(1-p), unit roundoff E = ulp(1)/2.
    _Float64x exP = one;
    for (int i = 0; i < p - 1; ++i) {
        exP /= two;
    }
    const _Float64x exE = exP / two;

    // rmin = 2^(emin-1)
    _Float64x exU = one;
    // emin is negative for IEEE binary formats.
    for (int i = 0; i < (-emin + 1); ++i) {
        exU /= two;
    }

    // rmax = (1 - 2^(-p)) * 2^emax  (computed without overflowing intermediates)
    _Float64x two_to_minus_p = one;
    for (int i = 0; i < p; ++i) {
        two_to_minus_p /= two;
    }
    _Float64x exO = one - two_to_minus_p;
    for (int i = 0; i < emax; ++i) {
        exO *= two;
    }

    const _Float64x exM = (_Float64x)emin;
    const _Float64x exL = (_Float64x)emax;

    // Safe minimum: max(rmin, (1/rmax)*(1+E)).
    const _Float64x small = one / exO;
    const _Float64x candidate = small * (one + exE);
    const _Float64x exS = (candidate >= exU) ? candidate : exU;

    const _Float64x exR = one;
    const _Float64x exZ = zero;

    // Fetch actual values from Rlamch__Float64x.
    const _Float64x gotE = Rlamch__Float64x("E");
    const _Float64x gotS = Rlamch__Float64x("S");
    const _Float64x gotB = Rlamch__Float64x("B");
    const _Float64x gotP = Rlamch__Float64x("P");
    const _Float64x gotN = Rlamch__Float64x("N");
    const _Float64x gotR = Rlamch__Float64x("R");
    const _Float64x gotM = Rlamch__Float64x("M");
    const _Float64x gotU = Rlamch__Float64x("U");
    const _Float64x gotL = Rlamch__Float64x("L");
    const _Float64x gotO = Rlamch__Float64x("O");
    const _Float64x gotZ = Rlamch__Float64x("Z");

    // Exact-value checks.
    assert_equal("E", gotE, exE);
    assert_equal("B", gotB, exB);
    assert_equal("P", gotP, exP);
    assert_equal("N", gotN, exN);
    assert_equal("R", gotR, exR);
    assert_equal("M", gotM, exM);
    assert_equal("U", gotU, exU);
    assert_equal("L", gotL, exL);
    assert_equal("O", gotO, exO);
    assert_equal("S", gotS, exS);
    assert_case(gotZ == exZ, "Z (dummy) is not 0");

    // Operational property checks (_Float64x)
    const _Float64x next = nextafterf64x(one, two);
    // const _Float64x next = (_Float64x)nextafterl((long double)one, (long double)two);

    // P == ulp(1) == nextafter(1,2) - 1
    assert_case((next - one) == gotP, "P is not equal to ulp(1) = nextafterf64x(1,2)-1");

    // Consistency: P == 2*E
    assert_case((gotE * two) == gotP, "expected P == 2*E");

    // Strong, rounding-mode-independent check: fl(1 + P) == nextafter(1,2)
    volatile _Float64x b = one + gotP;
    assert_case(b == next, "expected fl(1 + P) == nextafterf64x(1,2)");
    assert_case(b > one, "expected fl(1 + 2E) > 1");

    // Tie check: only assert under round-to-nearest (ties-to-even)
    if (std::fegetround() == FE_TONEAREST) {
        volatile _Float64x a = one + gotE; // 1 + ulp/2 (midpoint)
        assert_case(a == one, "expected fl(1 + E) == 1 under FE_TONEAREST (ties-to-even)");
    }
    // Range and reciprocal sanity checks.
    assert_case(gotU > zero && gotU < one, "U (rmin) is not in (0,1)");
    assert_case(gotO > one, "O (rmax) is not > 1");
    assert_case(gotS > zero && gotS < one, "S (sfmin) is not in (0,1)");
    assert_case((gotU * two) > gotU, "U*2 <= U (stuck at zero?)");
    assert_case((one / gotS) > zero, "1/S is not positive");
    assert_case(__builtin_isfinite(one / gotS), "1/S is not finite");

    // Netlib-style sfmin inequalities.
    assert_case((one / gotS) <= gotO, "1/S > O (violates safe minimum contract)");
    assert_case(gotS >= gotU, "S < U (violates sfmin >= rmin)");
    assert_case(gotS >= (one / gotO), "S < 1/O (violates sfmin >= 1/rmax)");

    // Cross-check: U*O == (1-E)*2^(emin+emax-1).
    const int exp_sum = emin + emax - 1;
    if (exp_sum >= -128 && exp_sum <= 128) {
        _Float64x pow2 = one;
        if (exp_sum > 0) {
            for (int i = 0; i < exp_sum; ++i)
                pow2 *= two;
        } else if (exp_sum < 0) {
            for (int i = 0; i < -exp_sum; ++i)
                pow2 /= two;
        }
        const _Float64x got_prod = gotU * gotO;
        const _Float64x expected_prod = (one - exE) * pow2;
        assert_case(got_prod == expected_prod, "U*O cross-check failed (inconsistent rmin/rmax model)");
    }

    if (print_values) {
        printf("Rlamch E: Epsilon                      ");
        printnum(gotE);
        printf("\n");
        printf("Rlamch S: Safe minimum                 ");
        printnum(gotS);
        printf("\n");
        printf("Rlamch B: Base                         ");
        printnum(gotB);
        printf("\n");
        printf("Rlamch P: Precision                    ");
        printnum(gotP);
        printf("\n");
        printf("Rlamch N: Number of digits in mantissa ");
        printnum(gotN);
        printf("\n");
        printf("Rlamch R: Rounding mode                ");
        printnum(gotR);
        printf("\n");
        printf("Rlamch M: Minimum exponent             ");
        printnum(gotM);
        printf("\n");
        printf("Rlamch U: Underflow threshold          ");
        printnum(gotU);
        printf("\n");
        printf("Rlamch L: Largest exponent             ");
        printnum(gotL);
        printf("\n");
        printf("Rlamch O: Overflow threshold           ");
        printnum(gotO);
        printf("\n");
    }
}
#endif

int main(int argc, char *argv[]) {
    printf("*** Testing Rlamch start ***\n");
#if defined ___MPLAPACK_BUILD_WITH_MPFR___
    Rlamch_mpfr_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH_GMP___
    Rlamch_gmp_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH_QD___
    Rlamch_qd_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH_DD___
    Rlamch_dd_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
    Rlamch_double_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
    Rlamch__Float64x_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___
    Rlamch__Float128_test();
#endif
    printf("*** Testing Rlamch successful ***\n");
    return (0);
}
