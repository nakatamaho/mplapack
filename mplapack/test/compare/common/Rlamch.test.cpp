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
#include <cstdarg>
#include <cstdio>
#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_compare_debug.h>
#include <mplapack_arithmetic_params.h>
#include <blas.h>
#include <lapack.h>

#define VERBOSE_TEST

// ---------------------------------------------------------------------------
// Dual output: write to stdout and to Rlamch.txt simultaneously.
// ---------------------------------------------------------------------------
static FILE *g_lamch_file = nullptr;

static void dual_printf(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    if (g_lamch_file) {
        va_start(ap, fmt);
        vfprintf(g_lamch_file, fmt, ap);
        va_end(ap);
    }
}

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
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        sprintnum(_spbuf, gotE);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), gotE);
        dual_printf("Rlamch E: Epsilon                      %80s    %80s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotS);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), gotS);
        dual_printf("Rlamch S: Safe minimum                 %80s    %80s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotB);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), gotB);
        dual_printf("Rlamch B: Base                         %80s    %80s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotP);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), gotP);
        dual_printf("Rlamch P: Precision                    %80s    %80s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotN);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), gotN);
        dual_printf("Rlamch N: Number of digits in mantissa %80s    %80s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotR);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), gotR);
        dual_printf("Rlamch R: Rounding mode                %80s    %80s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotM);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), gotM);
        dual_printf("Rlamch M: Minimum exponent             %80s    %80s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotU);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), gotU);
        dual_printf("Rlamch U: Underflow threshold          %80s    %80s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotL);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), gotL);
        dual_printf("Rlamch L: Largest exponent             %80s    %80s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotO);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), gotO);
        dual_printf("Rlamch O: Overflow threshold           %80s    %80s\n", _spbuf, _sphexbuf);
    }
}

template <typename BlueQ> static int classify_blue_mpfr_value(const BlueQ &q, const REAL &x) {
    REAL ax = x;
    if (ax < REAL(0.0))
        ax = -ax;
    if (ax > q.tbig)
        return +1;
    if (ax < q.tsml)
        return -1;
    return 0;
}

static void check_arithmetic_params_mpfr(const char *tag, bool print_values) {
    using mplapack::arithmetic_int;
    const auto p = mplapack::get_arithmetic_params<REAL>();
    const auto q = mplapack::get_blue_scaling_params<REAL>();
    const auto q2 = mplapack::make_blue_scaling_params(p);

    assert_equal_real(tag, "params.E", p.eps, Rlamch_mpfr("E"));
    assert_equal_real(tag, "params.S", p.sfmin, Rlamch_mpfr("S"));
    assert_equal_real(tag, "params.B", p.base, Rlamch_mpfr("B"));
    assert_equal_real(tag, "params.P", p.prec, Rlamch_mpfr("P"));
    assert_equal_real(tag, "params.R", p.rnd, Rlamch_mpfr("R"));
    assert_equal_real(tag, "params.U", p.rmin, Rlamch_mpfr("U"));
    assert_equal_real(tag, "params.O", p.rmax, Rlamch_mpfr("O"));

    assert_equal_real(tag, "params.N", mplapack::detail::to_rlamch_real<REAL>(p.t), Rlamch_mpfr("N"));
    assert_equal_real(tag, "params.M", mplapack::detail::to_rlamch_real<REAL>(p.emin), Rlamch_mpfr("M"));
    assert_equal_real(tag, "params.L", mplapack::detail::to_rlamch_real<REAL>(p.emax), Rlamch_mpfr("L"));

    assert_equal_real(tag, "params.prec_consistency", p.prec, p.eps * p.base);
    assert_equal_real(tag, "params.safmin", p.safmin, mplapack::detail::compute_safmin<REAL>(p.emin, p.emax));
    assert_equal_real(tag, "params.safmax", p.safmax, mplapack::detail::compute_safmax<REAL>(p.emin, p.emax));

    mpfr_assert_case(q.exp_tsml == q2.exp_tsml, tag, "ArithmeticParams->Blue builder mismatch: exp_tsml");
    mpfr_assert_case(q.exp_tbig == q2.exp_tbig, tag, "ArithmeticParams->Blue builder mismatch: exp_tbig");
    mpfr_assert_case(q.exp_ssml == q2.exp_ssml, tag, "ArithmeticParams->Blue builder mismatch: exp_ssml");
    mpfr_assert_case(q.exp_sbig == q2.exp_sbig, tag, "ArithmeticParams->Blue builder mismatch: exp_sbig");

    assert_equal_real(tag, "ArithmeticParams->Blue tsml", q.tsml, q2.tsml);
    assert_equal_real(tag, "ArithmeticParams->Blue tbig", q.tbig, q2.tbig);

    const arithmetic_int ex_exp_tsml = q2.exp_tsml;
    const arithmetic_int ex_exp_tbig = q2.exp_tbig;
    const arithmetic_int ex_exp_ssml = q2.exp_ssml;
    const arithmetic_int ex_exp_sbig = q2.exp_sbig;

    const arithmetic_int emin = p.emin;
    const arithmetic_int emax = p.emax;

    const bool ordering_valid = (ex_exp_tsml < 0) && (ex_exp_tbig > 0) && (ex_exp_ssml > ex_exp_tbig) && (ex_exp_sbig < ex_exp_tsml);

    const bool finite_valid = (q2.exp_ssml >= emin) && (q2.exp_ssml <= emax) && (q2.exp_sbig >= emin) && (q2.exp_sbig <= emax);

    // Stress / pathological MPFR environments may violate the usual Blue
    // ordering chain or finite-representability assumptions.
    const bool strict_blue_valid = ordering_valid && finite_valid;

    if (strict_blue_valid) {
        assert_equal_real(tag, "ArithmeticParams->Blue ssml", q.ssml, q2.ssml);
        assert_equal_real(tag, "ArithmeticParams->Blue sbig", q.sbig, q2.sbig);
    }

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("[params/%s] t=%lld emin=%lld emax=%lld\n", tag, (long long)p.t, (long long)p.emin, (long long)p.emax);

        sprintnum(_spbuf, p.safmin);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), p.safmin);
        dual_printf("[params/%s] safmin %80s    %80s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, p.safmax);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), p.safmax);
        dual_printf("[params/%s] safmax %80s    %80s\n", tag, _spbuf, _sphexbuf);

        dual_printf("[params/%s] builder exp_tsml=%lld exp_tbig=%lld exp_ssml=%lld exp_sbig=%lld\n", tag, (long long)q2.exp_tsml, (long long)q2.exp_tbig, (long long)q2.exp_ssml, (long long)q2.exp_sbig);

        sprintnum(_spbuf, q2.tsml);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), q2.tsml);
        dual_printf("[params/%s] builder tsml   %80s    %80s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.tbig);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), q2.tbig);
        dual_printf("[params/%s] builder tbig   %80s    %80s\n", tag, _spbuf, _sphexbuf);

        if (strict_blue_valid) {
            sprintnum(_spbuf, q2.ssml);
            sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), q2.ssml);
            dual_printf("[params/%s] builder ssml   %80s    %80s\n", tag, _spbuf, _sphexbuf);

            sprintnum(_spbuf, q2.sbig);
            sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), q2.sbig);
            dual_printf("[params/%s] builder sbig   %80s    %80s\n", tag, _spbuf, _sphexbuf);
        } else {
            dual_printf("[params/%s] builder ssml   [informational-only case omitted]\n", tag);
            dual_printf("[params/%s] builder sbig   [informational-only case omitted]\n", tag);
        }
    }
}

template <typename BlueQ> static void check_blue_threshold_boundaries_mpfr(const char *tag, const BlueQ &q, bool ordering_valid) {
    if (!ordering_valid)
        return;

    const REAL zero(0.0), one(1.0);

    REAL below_tsml = q.tsml;
    REAL above_tsml = q.tsml;
    REAL below_tbig = q.tbig;
    REAL above_tbig = q.tbig;
    mpfr_nextbelow(mpfr_ptr(below_tsml));
    mpfr_nextabove(mpfr_ptr(above_tsml));
    mpfr_nextbelow(mpfr_ptr(below_tbig));
    mpfr_nextabove(mpfr_ptr(above_tbig));

    mpfr_assert_case(below_tsml < q.tsml, tag, "BlueScale boundary: nextbelow(tsml) is not < tsml");
    mpfr_assert_case(above_tsml > q.tsml, tag, "BlueScale boundary: nextabove(tsml) is not > tsml");
    mpfr_assert_case(below_tbig < q.tbig, tag, "BlueScale boundary: nextbelow(tbig) is not < tbig");
    mpfr_assert_case(above_tbig > q.tbig, tag, "BlueScale boundary: nextabove(tbig) is not > tbig");

    mpfr_assert_case(classify_blue_mpfr_value(q, below_tsml) == -1, tag, "BlueScale boundary: nextbelow(tsml) must classify as small");
    mpfr_assert_case(classify_blue_mpfr_value(q, q.tsml) == 0, tag, "BlueScale boundary: tsml must classify as medium");
    mpfr_assert_case(classify_blue_mpfr_value(q, above_tsml) == 0, tag, "BlueScale boundary: nextabove(tsml) must classify as medium");
    mpfr_assert_case(classify_blue_mpfr_value(q, below_tbig) == 0, tag, "BlueScale boundary: nextbelow(tbig) must classify as medium");
    mpfr_assert_case(classify_blue_mpfr_value(q, q.tbig) == 0, tag, "BlueScale boundary: tbig must classify as medium");
    mpfr_assert_case(classify_blue_mpfr_value(q, above_tbig) == +1, tag, "BlueScale boundary: nextabove(tbig) must classify as big");

    mpfr_assert_case(classify_blue_mpfr_value(q, REAL(-1.0) * below_tsml) == -1, tag, "BlueScale boundary: -nextbelow(tsml) must classify as small");
    mpfr_assert_case(classify_blue_mpfr_value(q, REAL(-1.0) * above_tbig) == +1, tag, "BlueScale boundary: -nextabove(tbig) must classify as big");
    mpfr_assert_case(classify_blue_mpfr_value(q, zero) == -1, tag, "BlueScale boundary: 0 must classify as small");
    mpfr_assert_case(classify_blue_mpfr_value(q, one) == 0, tag, "BlueScale boundary: 1 must classify as medium");

    REAL scaled_small = below_tsml * q.ssml;
    REAL scaled_big = above_tbig * q.sbig;
    mpfr_assert_case(mpfr_number_p(mpfr_ptr(scaled_small)) != 0, tag, "BlueScale boundary: scaled nextbelow(tsml) is not finite");
    mpfr_assert_case(mpfr_number_p(mpfr_ptr(scaled_big)) != 0, tag, "BlueScale boundary: scaled nextabove(tbig) is not finite");
    mpfr_assert_case(classify_blue_mpfr_value(q, scaled_small) == 0, tag, "BlueScale boundary: nextbelow(tsml) * ssml must classify as medium");
    mpfr_assert_case(classify_blue_mpfr_value(q, scaled_big) == 0, tag, "BlueScale boundary: nextabove(tbig) * sbig must classify as medium");
}

// ---------------------------------------------------------------------------
// Blue scaling parameter check for MPFR.
// Must be called with the target MPFR environment already active.
// ---------------------------------------------------------------------------
static void check_blue_scaling_mpfr(const char *tag, bool print_values) {
    using mplapack::arithmetic_int;

    const auto q = mplapack::get_blue_scaling_params<REAL>();

    const mpfr_prec_t real_prec = current_real_default_prec();
    const mpfr_exp_t active_emin = mpfr_get_emin();
    const mpfr_exp_t active_emax = mpfr_get_emax();

    const arithmetic_int digits = static_cast<arithmetic_int>(real_prec);
    const arithmetic_int emin = static_cast<arithmetic_int>(active_emin);
    const arithmetic_int emax = static_cast<arithmetic_int>(active_emax);

    const arithmetic_int ex_exp_tsml = mplapack::detail::ceildiv2(emin - 1);
    const arithmetic_int ex_exp_tbig = mplapack::detail::floordiv2(emax - digits + 1);
    const arithmetic_int ex_exp_ssml = -mplapack::detail::floordiv2(emin - digits);
    const arithmetic_int ex_exp_sbig = -mplapack::detail::ceildiv2(emax + digits - 1);

    mpfr_assert_case(q.exp_tsml == ex_exp_tsml, tag, "BlueScale: exp_tsml mismatch");
    mpfr_assert_case(q.exp_tbig == ex_exp_tbig, tag, "BlueScale: exp_tbig mismatch");
    mpfr_assert_case(q.exp_ssml == ex_exp_ssml, tag, "BlueScale: exp_ssml mismatch");
    mpfr_assert_case(q.exp_sbig == ex_exp_sbig, tag, "BlueScale: exp_sbig mismatch");

    const mpfr_rnd_t rnd = current_real_default_rnd();

    // Always-valid side: tsml / tbig can be checked directly.
    const REAL ex_tsml = pow2_from_exp(static_cast<mpfr_exp_t>(q.exp_tsml), rnd);
    const REAL ex_tbig = pow2_from_exp(static_cast<mpfr_exp_t>(q.exp_tbig), rnd);

    assert_equal_real(tag, "BlueScale tsml", q.tsml, ex_tsml);
    assert_equal_real(tag, "BlueScale tbig", q.tbig, ex_tbig);

    const REAL zero(0.0), one(1.0);
    mpfr_assert_case(q.tsml > zero, tag, "BlueScale: tsml must be positive");
    mpfr_assert_case(q.tsml < one, tag, "BlueScale: tsml must be < 1");
    mpfr_assert_case(q.tbig > one, tag, "BlueScale: tbig must be > 1");

    // Strict Blue checks only make sense when the usual ordering chain can hold:
    //   0 < sbig < tsml < 1 < tbig < ssml
    const bool ordering_valid = (ex_exp_tsml < 0) && (ex_exp_tbig > 0) && (ex_exp_ssml > ex_exp_tbig) && (ex_exp_sbig < ex_exp_tsml);

    // In addition, ssml/sbig must be representable as finite MPFR numbers
    // under the active exponent range.
    const bool finite_valid = (q.exp_ssml >= static_cast<arithmetic_int>(active_emin)) && (q.exp_ssml <= static_cast<arithmetic_int>(active_emax)) && (q.exp_sbig >= static_cast<arithmetic_int>(active_emin)) && (q.exp_sbig <= static_cast<arithmetic_int>(active_emax));

    const bool strict_blue_valid = ordering_valid && finite_valid;

    if (!strict_blue_valid) {
        dual_printf("BlueScale note [%s]: informational-only MPFR Blue case "
                    "(exp_tsml=%lld exp_tbig=%lld exp_ssml=%lld exp_sbig=%lld "
                    "emin=%lld emax=%lld ordering_valid=%d finite_valid=%d) "
                    "― strict ssml/sbig checks skipped\n",
                    tag, (long long)ex_exp_tsml, (long long)ex_exp_tbig, (long long)ex_exp_ssml, (long long)ex_exp_sbig, (long long)emin, (long long)emax, ordering_valid ? 1 : 0, finite_valid ? 1 : 0);
    }

    if (strict_blue_valid) {
        REAL ex_ssml(1.0), ex_sbig(1.0);
        mpfr_mul_2si(mpfr_ptr(ex_ssml), mpfr_ptr(ex_ssml), checked_exp_to_long(static_cast<mpfr_exp_t>(q.exp_ssml), "BlueScale: exp_ssml overflows long"), rnd);
        mpfr_mul_2si(mpfr_ptr(ex_sbig), mpfr_ptr(ex_sbig), checked_exp_to_long(static_cast<mpfr_exp_t>(q.exp_sbig), "BlueScale: exp_sbig overflows long"), rnd);

        assert_equal_real(tag, "BlueScale ssml", q.ssml, ex_ssml);
        assert_equal_real(tag, "BlueScale sbig", q.sbig, ex_sbig);

        mpfr_assert_case(q.ssml > q.tbig, tag, "BlueScale: ssml must be > tbig");
        mpfr_assert_case(q.sbig > zero, tag, "BlueScale: sbig must be positive");
        mpfr_assert_case(q.sbig < q.tsml, tag, "BlueScale: sbig must be < tsml");

        REAL prod_ts = q.tsml * q.ssml;
        mpfr_assert_case(mpfr_number_p(mpfr_ptr(prod_ts)) != 0, tag, "BlueScale: tsml*ssml is not a number");
        mpfr_assert_case(prod_ts > zero, tag, "BlueScale: tsml*ssml must be positive");

        REAL prod_bs = q.tbig * q.sbig;
        mpfr_assert_case(mpfr_number_p(mpfr_ptr(prod_bs)) != 0, tag, "BlueScale: tbig*sbig is not a number");
        mpfr_assert_case(prod_bs > zero, tag, "BlueScale: tbig*sbig must be positive");

        REAL tsml_sq = q.tsml * q.tsml;
        mpfr_assert_case(mpfr_number_p(mpfr_ptr(tsml_sq)) != 0, tag, "BlueScale: tsml^2 is not a number");
        mpfr_assert_case(tsml_sq > zero, tag, "BlueScale: tsml^2 underflowed to zero");
    }

    // tbig^2 is still a useful check even in informational mode.
    REAL tbig_sq = q.tbig * q.tbig;
    mpfr_assert_case(mpfr_number_p(mpfr_ptr(tbig_sq)) != 0, tag, "BlueScale: tbig^2 overflowed (not a finite number)");

    check_blue_threshold_boundaries_mpfr(tag, q, strict_blue_valid);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("BlueScale exp_tsml: %lld\n", (long long)q.exp_tsml);
        dual_printf("BlueScale exp_tbig: %lld\n", (long long)q.exp_tbig);
        dual_printf("BlueScale exp_ssml: %lld\n", (long long)q.exp_ssml);
        dual_printf("BlueScale exp_sbig: %lld\n", (long long)q.exp_sbig);

        sprintnum(_spbuf, q.tsml);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), q.tsml);
        dual_printf("BlueScale tsml:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.tbig);
        sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), q.tbig);
        dual_printf("BlueScale tbig:     %40s    %40s\n", _spbuf, _sphexbuf);

        if (strict_blue_valid) {
            sprintnum(_spbuf, q.ssml);
            sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), q.ssml);
            dual_printf("BlueScale ssml:     %40s    %40s\n", _spbuf, _sphexbuf);

            sprintnum(_spbuf, q.sbig);
            sprinthex_mpfr_fixed(_sphexbuf, sizeof(_sphexbuf), q.sbig);
            dual_printf("BlueScale sbig:     %40s    %40s\n", _spbuf, _sphexbuf);
        } else {
            dual_printf("BlueScale ssml:     [informational-only case omitted]\n");
            dual_printf("BlueScale sbig:     [informational-only case omitted]\n");
        }
    }
}

static void run_mpfr_env_test(const char *tag, const MpfrEnvSnapshot &cfg, bool print_values) {
    MpfrEnvGuard guard; // Save current environment; restore on scope exit.
    mpfr_env_apply(cfg);
    check_arithmetic_params_mpfr(tag, print_values);
    check_lamch_mpfr_values(tag, cfg, print_values);
    check_blue_scaling_mpfr(tag, print_values);
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
        {"binary512", false, (mpfr_prec_t)513, (mpfr_exp_t)-65533, (mpfr_exp_t)65536},

        // Emulations: IEEE-like formats.
        {"binary64", false, (mpfr_prec_t)53, (mpfr_exp_t)-1021, (mpfr_exp_t)1024},
        {"binary80", false, (mpfr_prec_t)64, (mpfr_exp_t)-16381, (mpfr_exp_t)16384},
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
            dual_printf("[MPFR] case=%s mpfr_prec=%llu real_prec=%llu emin=%lld emax=%lld rnd=%s\n", tagbuf, (unsigned long long)cfg.mpfr_prec, (unsigned long long)cfg.real_prec, (long long)cfg.emin, (long long)cfg.emax, rc.name);

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

static void check_arithmetic_params_gmp(const char *tag, mp_bitcnt_t prec_bits, bool print_values) {
    // p.sfmin corresponds to Rlamch("S") = safe minimum.
    // p.rmin corresponds to Rlamch("U") = underflow threshold.
    // p.safmin/p.safmax are internal finite scaling-side constants used by ArithmeticParams.

    const auto p = mplapack::get_arithmetic_params<REAL>();
    const auto q = mplapack::get_blue_scaling_params<REAL>();
    const auto q2 = mplapack::make_blue_scaling_params(p);

    gmp_assert_equal_real(tag, "params.E", p.eps, Rlamch_gmp("E"));
    gmp_assert_equal_real(tag, "params.S", p.sfmin, Rlamch_gmp("S"));
    gmp_assert_equal_real(tag, "params.B", p.base, Rlamch_gmp("B"));
    gmp_assert_equal_real(tag, "params.P", p.prec, Rlamch_gmp("P"));
    gmp_assert_equal_real(tag, "params.R", p.rnd, Rlamch_gmp("R"));
    gmp_assert_equal_real(tag, "params.U", p.rmin, Rlamch_gmp("U"));
    gmp_assert_equal_real(tag, "params.O", p.rmax, Rlamch_gmp("O"));

    gmp_assert_equal_real(tag, "params.N", mplapack::detail::to_rlamch_real<REAL>(p.t), Rlamch_gmp("N"));
    gmp_assert_equal_real(tag, "params.M", mplapack::detail::to_rlamch_real<REAL>(p.emin), Rlamch_gmp("M"));
    gmp_assert_equal_real(tag, "params.L", mplapack::detail::to_rlamch_real<REAL>(p.emax), Rlamch_gmp("L"));

    gmp_assert_case(p.t == static_cast<mplapack::arithmetic_int>(prec_bits), tag, "ArithmeticParams.t mismatch");
    gmp_assert_equal_real(tag, "params.prec_consistency", p.prec, p.eps * p.base);
    gmp_assert_equal_real(tag, "params.safmin", p.safmin, mplapack::detail::compute_safmin<REAL>(p.emin, p.emax));
    gmp_assert_equal_real(tag, "params.safmax", p.safmax, mplapack::detail::compute_safmax<REAL>(p.emin, p.emax));

    gmp_assert_case(p.sfmin >= p.rmin, tag, "safe minimum must be >= underflow threshold");
    gmp_assert_case(p.base > REAL(1.0), tag, "base must be > 1");
    gmp_assert_case(p.eps > REAL(0.0), tag, "epsilon must be positive");
    gmp_assert_case(p.prec > REAL(0.0), tag, "precision must be positive");
    gmp_assert_case(p.rmin > REAL(0.0), tag, "underflow threshold must be positive");
    gmp_assert_case(p.rmax > REAL(1.0), tag, "overflow threshold must be > 1");
    gmp_assert_case(p.emin < p.emax, tag, "minimum exponent must be less than maximum exponent");

    gmp_assert_case(q.exp_tsml == q2.exp_tsml, tag, "ArithmeticParams->Blue builder mismatch: exp_tsml");
    gmp_assert_case(q.exp_tbig == q2.exp_tbig, tag, "ArithmeticParams->Blue builder mismatch: exp_tbig");
    gmp_assert_case(q.exp_ssml == q2.exp_ssml, tag, "ArithmeticParams->Blue builder mismatch: exp_ssml");
    gmp_assert_case(q.exp_sbig == q2.exp_sbig, tag, "ArithmeticParams->Blue builder mismatch: exp_sbig");

    gmp_assert_equal_real(tag, "ArithmeticParams->Blue tsml", q.tsml, q2.tsml);
    gmp_assert_equal_real(tag, "ArithmeticParams->Blue tbig", q.tbig, q2.tbig);
    gmp_assert_equal_real(tag, "ArithmeticParams->Blue ssml", q.ssml, q2.ssml);
    gmp_assert_equal_real(tag, "ArithmeticParams->Blue sbig", q.sbig, q2.sbig);

    gmp_assert_case(q.tsml > REAL(0.0), tag, "Blue tsml must be positive");
    gmp_assert_case(q.tbig > REAL(1.0), tag, "Blue tbig must be > 1");
    gmp_assert_case(q.ssml > REAL(1.0), tag, "Blue ssml must be > 1");
    gmp_assert_case(q.sbig > REAL(0.0), tag, "Blue sbig must be positive");
    gmp_assert_case(q.sbig < REAL(1.0), tag, "Blue sbig must be < 1");
    gmp_assert_case(q.tsml < q.tbig, tag, "Blue thresholds must satisfy tsml < tbig");
    gmp_assert_case(q.exp_tsml < q.exp_tbig, tag, "Blue exponents must satisfy exp_tsml < exp_tbig");
    gmp_assert_case(q.exp_sbig < q.exp_ssml, tag, "Blue exponents must satisfy exp_sbig < exp_ssml");

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("[params/%s] t=%lld emin=%lld emax=%lld\n",
                    tag, (long long)p.t, (long long)p.emin, (long long)p.emax);

        sprintnum(_spbuf, p.safmin);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), p.safmin);
        dual_printf("[params/%s] safmin %130s    %130s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, p.safmax);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), p.safmax);
        dual_printf("[params/%s] safmax %130s    %130s\n", tag, _spbuf, _sphexbuf);

        dual_printf("[params/%s] builder exp_tsml=%lld exp_tbig=%lld exp_ssml=%lld exp_sbig=%lld\n",
                    tag,
                    (long long)q2.exp_tsml, (long long)q2.exp_tbig,
                    (long long)q2.exp_ssml, (long long)q2.exp_sbig);

        sprintnum(_spbuf, q2.tsml);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), q2.tsml);
        dual_printf("[params/%s] builder tsml   %130s    %130s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.tbig);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), q2.tbig);
        dual_printf("[params/%s] builder tbig   %130s    %130s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.ssml);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), q2.ssml);
        dual_printf("[params/%s] builder ssml   %130s    %130s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.sbig);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), q2.sbig);
        dual_printf("[params/%s] builder sbig   %130s    %130s\n", tag, _spbuf, _sphexbuf);
    }
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

template <typename BlueQ> static int classify_blue_gmp_value(const BlueQ &q, const REAL &x) {
    REAL ax = x;
    if (ax < REAL(0.0))
        ax = -ax;
    if (ax > q.tbig)
        return +1;
    if (ax < q.tsml)
        return -1;
    return 0;
}

template <typename BlueQ> static void check_blue_threshold_boundaries_gmp(const char *tag, const BlueQ &q, const REAL &delta, bool ordering_valid) {
    if (!ordering_valid)
        return;

    const REAL zero(0.0), one(1.0), minus_one(-1.0);

    const REAL below_tsml = q.tsml * (one - delta);
    const REAL above_tsml = q.tsml * (one + delta);
    const REAL below_tbig = q.tbig * (one - delta);
    const REAL above_tbig = q.tbig * (one + delta);

    gmp_assert_case(below_tsml < q.tsml, tag, "BlueScale boundary: below-tsml probe is not < tsml");
    gmp_assert_case(above_tsml > q.tsml, tag, "BlueScale boundary: above-tsml probe is not > tsml");
    gmp_assert_case(below_tbig < q.tbig, tag, "BlueScale boundary: below-tbig probe is not < tbig");
    gmp_assert_case(above_tbig > q.tbig, tag, "BlueScale boundary: above-tbig probe is not > tbig");

    gmp_assert_case(classify_blue_gmp_value(q, below_tsml) == -1, tag, "BlueScale boundary: below-tsml probe must classify as small");
    gmp_assert_case(classify_blue_gmp_value(q, q.tsml) == 0, tag, "BlueScale boundary: tsml must classify as medium");
    gmp_assert_case(classify_blue_gmp_value(q, above_tsml) == 0, tag, "BlueScale boundary: above-tsml probe must classify as medium");
    gmp_assert_case(classify_blue_gmp_value(q, below_tbig) == 0, tag, "BlueScale boundary: below-tbig probe must classify as medium");
    gmp_assert_case(classify_blue_gmp_value(q, q.tbig) == 0, tag, "BlueScale boundary: tbig must classify as medium");
    gmp_assert_case(classify_blue_gmp_value(q, above_tbig) == +1, tag, "BlueScale boundary: above-tbig probe must classify as big");

    gmp_assert_case(classify_blue_gmp_value(q, minus_one * below_tsml) == -1, tag, "BlueScale boundary: negative below-tsml probe must classify as small");
    gmp_assert_case(classify_blue_gmp_value(q, minus_one * above_tbig) == +1, tag, "BlueScale boundary: negative above-tbig probe must classify as big");
    gmp_assert_case(classify_blue_gmp_value(q, zero) == -1, tag, "BlueScale boundary: 0 must classify as small");
    gmp_assert_case(classify_blue_gmp_value(q, one) == 0, tag, "BlueScale boundary: 1 must classify as medium");

    const REAL scaled_small = below_tsml * q.ssml;
    const REAL scaled_big = above_tbig * q.sbig;
    gmp_assert_case(classify_blue_gmp_value(q, scaled_small) == 0, tag, "BlueScale boundary: below-tsml probe * ssml must classify as medium");
    gmp_assert_case(classify_blue_gmp_value(q, scaled_big) == 0, tag, "BlueScale boundary: above-tbig probe * sbig must classify as medium");
}

static void check_blue_scaling_gmp(const char *tag, mp_bitcnt_t prec_bits, bool print_values) {
    using mplapack::arithmetic_int;

    const auto q = mplapack::get_blue_scaling_params<REAL>();

    const mp_bitcnt_t Eexp = get_max_safe_exponent();
    const arithmetic_int emin = static_cast<arithmetic_int>(1) - static_cast<arithmetic_int>(Eexp);
    const arithmetic_int emax = static_cast<arithmetic_int>(Eexp);
    const arithmetic_int digits = static_cast<arithmetic_int>(prec_bits);

    const arithmetic_int ex_exp_tsml = mplapack::detail::ceildiv2(emin - 1);
    const arithmetic_int ex_exp_tbig = mplapack::detail::floordiv2(emax - digits + 1);
    const arithmetic_int ex_exp_ssml = -mplapack::detail::floordiv2(emin - digits);
    const arithmetic_int ex_exp_sbig = -mplapack::detail::ceildiv2(emax + digits - 1);

    gmp_assert_case(q.exp_tsml == ex_exp_tsml, tag, "BlueScale: exp_tsml mismatch");
    gmp_assert_case(q.exp_tbig == ex_exp_tbig, tag, "BlueScale: exp_tbig mismatch");
    gmp_assert_case(q.exp_ssml == ex_exp_ssml, tag, "BlueScale: exp_ssml mismatch");
    gmp_assert_case(q.exp_sbig == ex_exp_sbig, tag, "BlueScale: exp_sbig mismatch");

    const REAL one(1.0), zero(0.0);
    REAL ex_tsml(1.0), ex_tbig(1.0), ex_ssml(1.0), ex_sbig(1.0);
    if (q.exp_tsml >= 0)
        mpf_mul_2exp(ex_tsml.get_mpf_t(), ex_tsml.get_mpf_t(), static_cast<mp_bitcnt_t>(q.exp_tsml));
    else
        mpf_div_2exp(ex_tsml.get_mpf_t(), ex_tsml.get_mpf_t(), static_cast<mp_bitcnt_t>(-q.exp_tsml));
    if (q.exp_tbig >= 0)
        mpf_mul_2exp(ex_tbig.get_mpf_t(), ex_tbig.get_mpf_t(), static_cast<mp_bitcnt_t>(q.exp_tbig));
    else
        mpf_div_2exp(ex_tbig.get_mpf_t(), ex_tbig.get_mpf_t(), static_cast<mp_bitcnt_t>(-q.exp_tbig));
    if (q.exp_ssml >= 0)
        mpf_mul_2exp(ex_ssml.get_mpf_t(), ex_ssml.get_mpf_t(), static_cast<mp_bitcnt_t>(q.exp_ssml));
    else
        mpf_div_2exp(ex_ssml.get_mpf_t(), ex_ssml.get_mpf_t(), static_cast<mp_bitcnt_t>(-q.exp_ssml));
    if (q.exp_sbig >= 0)
        mpf_mul_2exp(ex_sbig.get_mpf_t(), ex_sbig.get_mpf_t(), static_cast<mp_bitcnt_t>(q.exp_sbig));
    else
        mpf_div_2exp(ex_sbig.get_mpf_t(), ex_sbig.get_mpf_t(), static_cast<mp_bitcnt_t>(-q.exp_sbig));

    gmp_assert_equal_real(tag, "BlueScale tsml", q.tsml, ex_tsml);
    gmp_assert_equal_real(tag, "BlueScale tbig", q.tbig, ex_tbig);
    gmp_assert_equal_real(tag, "BlueScale ssml", q.ssml, ex_ssml);
    gmp_assert_equal_real(tag, "BlueScale sbig", q.sbig, ex_sbig);

    const bool ordering_valid = (ex_exp_tsml < 0) && (ex_exp_tbig > 0) && (ex_exp_ssml > ex_exp_tbig) && (ex_exp_sbig < ex_exp_tsml);
    if (!ordering_valid) {
        dual_printf("BlueScale note [%s]: ordering does not hold "
                    "(exp_tsml=%lld exp_tbig=%lld exp_ssml=%lld exp_sbig=%lld "
                    "emin=%lld emax=%lld) ― non-standard GMP env, "
                    "ssml/sbig checks skipped\n",
                    tag, (long long)ex_exp_tsml, (long long)ex_exp_tbig, (long long)ex_exp_ssml, (long long)ex_exp_sbig, (long long)emin, (long long)emax);
    }

    gmp_assert_equal_real(tag, "BlueScale tsml", q.tsml, ex_tsml);
    gmp_assert_equal_real(tag, "BlueScale tbig", q.tbig, ex_tbig);
    if (ordering_valid) {
        gmp_assert_equal_real(tag, "BlueScale ssml", q.ssml, ex_ssml);
        gmp_assert_equal_real(tag, "BlueScale sbig", q.sbig, ex_sbig);
    }

    gmp_assert_case(q.tsml > zero, tag, "BlueScale: tsml must be positive");
    gmp_assert_case(q.tsml < one, tag, "BlueScale: tsml must be < 1");
    gmp_assert_case(q.tbig > one, tag, "BlueScale: tbig must be > 1");
    if (ordering_valid) {
        gmp_assert_case(q.ssml > q.tbig, tag, "BlueScale: ssml must be > tbig");
        gmp_assert_case(q.sbig > zero, tag, "BlueScale: sbig must be positive");
        gmp_assert_case(q.sbig < q.tsml, tag, "BlueScale: sbig must be < tsml");
    }

    if (ordering_valid) {
        const REAL prod_ts = q.tsml * q.ssml;
        gmp_assert_case(prod_ts > zero, tag, "BlueScale: tsml*ssml must be positive");
        const REAL prod_bs = q.tbig * q.sbig;
        gmp_assert_case(prod_bs > zero, tag, "BlueScale: tbig*sbig must be positive");
        const REAL tsml_sq = q.tsml * q.tsml;
        gmp_assert_case(tsml_sq > zero, tag, "BlueScale: tsml^2 underflowed to zero");
    }

    const REAL delta = Rlamch_gmp("P");
    check_blue_threshold_boundaries_gmp(tag, q, delta, ordering_valid);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("BlueScale exp_tsml: %lld\n", (long long)q.exp_tsml);
        dual_printf("BlueScale exp_tbig: %lld\n", (long long)q.exp_tbig);
        dual_printf("BlueScale exp_ssml: %lld\n", (long long)q.exp_ssml);
        dual_printf("BlueScale exp_sbig: %lld\n", (long long)q.exp_sbig);

        sprintnum(_spbuf, q.tsml);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), q.tsml);
        dual_printf("BlueScale tsml:     %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.tbig);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), q.tbig);
        dual_printf("BlueScale tbig:     %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.ssml);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), q.ssml);
        dual_printf("BlueScale ssml:     %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.sbig);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), q.sbig);
        dual_printf("BlueScale sbig:     %130s    %130s\n", _spbuf, _sphexbuf);
    }
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
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        sprintnum(_spbuf, gotE);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), gotE);
        dual_printf("Rlamch E: Epsilon                      %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotS);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), gotS);
        dual_printf("Rlamch S: Safe minimum                 %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotB);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), gotB);
        dual_printf("Rlamch B: Base                         %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotP);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), gotP);
        dual_printf("Rlamch P: Precision                    %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotN);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), gotN);
        dual_printf("Rlamch N: Number of digits in mantissa %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotR);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), gotR);
        dual_printf("Rlamch R: Rounding mode                %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotM);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), gotM);
        dual_printf("Rlamch M: Minimum exponent             %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotU);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), gotU);
        dual_printf("Rlamch U: Underflow threshold          %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotL);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), gotL);
        dual_printf("Rlamch L: Largest exponent             %130s    %130s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotO);
        sprinthex_gmp_fixed(_sphexbuf, sizeof(_sphexbuf), gotO);
        dual_printf("Rlamch O: Overflow threshold           %130s    %130s\n", _spbuf, _sphexbuf);
    }

    check_blue_scaling_gmp(tag, prec_bits, print_values);
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

    check_arithmetic_params_gmp("current", prec_bits, print_values);
    check_lamch_gmp_values("current", prec_bits, print_values);
    // Additional precision cases requested: 4096, 128, 64 bits.
    // Note: mpf_set_default_prec() sets the default for newly created mpf values.
    const mp_bitcnt_t requested_precisions[] = {4096, 128, 64};

    for (size_t i = 0; i < (sizeof(requested_precisions) / sizeof(requested_precisions[0])); ++i) {
        dual_printf("[GMP] precision=%ld\n", requested_precisions[i]);
        const mp_bitcnt_t req = requested_precisions[i];
        const mp_bitcnt_t saved = mpf_get_default_prec();
        mpf_set_default_prec(req);

        // Use the effective precision actually used by REAL/mpf objects.
        const REAL probe = 1.0;
        const mp_bitcnt_t eff = mpf_get_prec(probe.get_mpf_t());

        char tagbuf[64];
        snprintf(tagbuf, sizeof(tagbuf), "prec=%llu (eff=%llu)", (unsigned long long)req, (unsigned long long)eff);

        check_arithmetic_params_gmp(tagbuf, eff, print_values);
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
// This test matches the current arithmetic_params-based implementation:
// - E and P are fixed to MPLAPACK's canonical QD literals
// - U and O are taken from qd_real::_min_normalized / qd_real::_max
// - N is taken from std::numeric_limits<qd_real>::digits
// - S follows the Netlib-style "safe minimum" logic using U, O, and E
//
// IMPORTANT:
// qd_real is a quad-double type (sum of 4 doubles). Arithmetic does not round
// to a fixed p-bit floating-point format, so classic half-ULP tie tests such as
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

static void check_arithmetic_params_qd(const char *tag, bool print_values) {
    const auto p = mplapack::get_arithmetic_params<qd_real>();
    const auto q = mplapack::get_blue_scaling_params<qd_real>();
    const auto q2 = mplapack::make_blue_scaling_params(p);

    assert_equal_qd(tag, "params.E", p.eps, Rlamch_qd("E"));
    assert_equal_qd(tag, "params.S", p.sfmin, Rlamch_qd("S"));
    assert_equal_qd(tag, "params.B", p.base, Rlamch_qd("B"));
    assert_equal_qd(tag, "params.P", p.prec, Rlamch_qd("P"));
    assert_equal_qd(tag, "params.R", p.rnd, Rlamch_qd("R"));
    assert_equal_qd(tag, "params.U", p.rmin, Rlamch_qd("U"));
    assert_equal_qd(tag, "params.O", p.rmax, Rlamch_qd("O"));

    assert_equal_qd(tag, "params.N", mplapack::detail::to_rlamch_real<qd_real>(p.t), Rlamch_qd("N"));
    assert_equal_qd(tag, "params.M", mplapack::detail::to_rlamch_real<qd_real>(p.emin), Rlamch_qd("M"));
    assert_equal_qd(tag, "params.L", mplapack::detail::to_rlamch_real<qd_real>(p.emax), Rlamch_qd("L"));

    assert_equal_qd(tag, "params.prec_consistency", p.prec, p.eps * p.base);
    assert_equal_qd(tag, "params.safmin", p.safmin, mplapack::detail::compute_safmin<qd_real>(p.emin, p.emax));
    assert_equal_qd(tag, "params.safmax", p.safmax, mplapack::detail::compute_safmax<qd_real>(p.emin, p.emax));

    qd_assert_case(q.exp_tsml == q2.exp_tsml, tag, "ArithmeticParams->Blue builder mismatch: exp_tsml");
    qd_assert_case(q.exp_tbig == q2.exp_tbig, tag, "ArithmeticParams->Blue builder mismatch: exp_tbig");
    qd_assert_case(q.exp_ssml == q2.exp_ssml, tag, "ArithmeticParams->Blue builder mismatch: exp_ssml");
    qd_assert_case(q.exp_sbig == q2.exp_sbig, tag, "ArithmeticParams->Blue builder mismatch: exp_sbig");

    assert_equal_qd(tag, "ArithmeticParams->Blue tsml", q.tsml, q2.tsml);
    assert_equal_qd(tag, "ArithmeticParams->Blue tbig", q.tbig, q2.tbig);
    assert_equal_qd(tag, "ArithmeticParams->Blue ssml", q.ssml, q2.ssml);
    assert_equal_qd(tag, "ArithmeticParams->Blue sbig", q.sbig, q2.sbig);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("[params/%s] t=%lld emin=%lld emax=%lld\n",
                    tag, (long long)p.t, (long long)p.emin, (long long)p.emax);

        sprintnum(_spbuf, p.safmin);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), p.safmin);
        dual_printf("[params/%s] safmin %70s    %70s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, p.safmax);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), p.safmax);
        dual_printf("[params/%s] safmax %70s    %70s\n", tag, _spbuf, _sphexbuf);

        dual_printf("[params/%s] builder exp_tsml=%lld exp_tbig=%lld exp_ssml=%lld exp_sbig=%lld\n",
                    tag,
                    (long long)q2.exp_tsml, (long long)q2.exp_tbig,
                    (long long)q2.exp_ssml, (long long)q2.exp_sbig);

        sprintnum(_spbuf, q2.tsml);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), q2.tsml);
        dual_printf("[params/%s] builder tsml   %70s    %70s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.tbig);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), q2.tbig);
        dual_printf("[params/%s] builder tbig   %70s    %70s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.ssml);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), q2.ssml);
        dual_printf("[params/%s] builder ssml   %70s    %70s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.sbig);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), q2.sbig);
        dual_printf("[params/%s] builder sbig   %70s    %70s\n", tag, _spbuf, _sphexbuf);
    }
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
    ex.E = qd_real(+0x1.fffffffffffffp-209, +0x0.0000000000000p+0000, +0x0.0000000000000p+0000, +0x0.0000000000000p+0000);
    ex.P = qd_real(+0x1.fffffffffffffp-208, +0x0.0000000000000p+0000, +0x0.0000000000000p+0000, +0x0.0000000000000p+0000);
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
    const qd_real invS = one / S; // S = min_normalized, 1/S is in normal range: safe

    qd_assert_case(invS <= O, tag, "1/S > O (violates safe minimum contract)");
    qd_assert_case(S >= U, tag, "S < U (violates sfmin >= rmin)");

    // 1/O is in the subnormal double range; qd_real division is unreliable there.
    // Use double arithmetic, which handles subnormals correctly.
    const double invO_d = 1.0 / to_double(O);
    qd_assert_case(to_double(S) >= invO_d, tag, "S < 1/O (violates sfmin >= 1/rmax)");
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

template <typename BlueQ> static int classify_blue_qd_value(const BlueQ &q, const qd_real &x) {
    qd_real ax = x;
    if (ax < qd_real(0.0))
        ax = -ax;
    if (ax > q.tbig)
        return +1;
    if (ax < q.tsml)
        return -1;
    return 0;
}

template <typename BlueQ> static void check_blue_threshold_boundaries_qd(const char *tag, const BlueQ &q, const qd_real &delta) {
    const qd_real zero(0.0), one(1.0), minus_one(-1.0);

    const qd_real below_tsml = q.tsml * (one - delta);
    const qd_real above_tsml = q.tsml * (one + delta);
    const qd_real below_tbig = q.tbig * (one - delta);
    const qd_real above_tbig = q.tbig * (one + delta);

    qd_assert_case(below_tsml < q.tsml, tag, "BlueScale boundary: below-tsml probe is not < tsml");
    qd_assert_case(above_tsml > q.tsml, tag, "BlueScale boundary: above-tsml probe is not > tsml");
    qd_assert_case(below_tbig < q.tbig, tag, "BlueScale boundary: below-tbig probe is not < tbig");
    qd_assert_case(above_tbig > q.tbig, tag, "BlueScale boundary: above-tbig probe is not > tbig");

    qd_assert_case(classify_blue_qd_value(q, below_tsml) == -1, tag, "BlueScale boundary: below-tsml probe must classify as small");
    qd_assert_case(classify_blue_qd_value(q, q.tsml) == 0, tag, "BlueScale boundary: tsml must classify as medium");
    qd_assert_case(classify_blue_qd_value(q, above_tsml) == 0, tag, "BlueScale boundary: above-tsml probe must classify as medium");
    qd_assert_case(classify_blue_qd_value(q, below_tbig) == 0, tag, "BlueScale boundary: below-tbig probe must classify as medium");
    qd_assert_case(classify_blue_qd_value(q, q.tbig) == 0, tag, "BlueScale boundary: tbig must classify as medium");
    qd_assert_case(classify_blue_qd_value(q, above_tbig) == +1, tag, "BlueScale boundary: above-tbig probe must classify as big");

    qd_assert_case(classify_blue_qd_value(q, minus_one * below_tsml) == -1, tag, "BlueScale boundary: negative below-tsml probe must classify as small");
    qd_assert_case(classify_blue_qd_value(q, minus_one * above_tbig) == +1, tag, "BlueScale boundary: negative above-tbig probe must classify as big");
    qd_assert_case(classify_blue_qd_value(q, zero) == -1, tag, "BlueScale boundary: 0 must classify as small");
    qd_assert_case(classify_blue_qd_value(q, one) == 0, tag, "BlueScale boundary: 1 must classify as medium");

    const qd_real scaled_small = below_tsml * q.ssml;
    const qd_real scaled_big = above_tbig * q.sbig;
    qd_assert_case(classify_blue_qd_value(q, scaled_small) == 0, tag, "BlueScale boundary: below-tsml probe * ssml must classify as medium");
    qd_assert_case(classify_blue_qd_value(q, scaled_big) == 0, tag, "BlueScale boundary: above-tbig probe * sbig must classify as medium");
}

static void check_blue_scaling_qd(const char *tag, bool print_values) {
    using mplapack::arithmetic_int;

    const auto q = mplapack::get_blue_scaling_params<qd_real>();

    int emin_int = 0;
    (void)std::frexp(qd_real::_min_normalized, &emin_int);
    const arithmetic_int emin = static_cast<arithmetic_int>(emin_int);
    const arithmetic_int emax = static_cast<arithmetic_int>(std::numeric_limits<double>::max_exponent);
    const arithmetic_int digits = static_cast<arithmetic_int>(std::numeric_limits<qd_real>::digits);

    qd_assert_case(q.exp_tsml == mplapack::detail::ceildiv2(emin - 1), tag, "BlueScale: exp_tsml mismatch");
    qd_assert_case(q.exp_tbig == mplapack::detail::floordiv2(emax - digits + 1), tag, "BlueScale: exp_tbig mismatch");
    qd_assert_case(q.exp_ssml == -mplapack::detail::floordiv2(emin - digits), tag, "BlueScale: exp_ssml mismatch");
    qd_assert_case(q.exp_sbig == -mplapack::detail::ceildiv2(emax + digits - 1), tag, "BlueScale: exp_sbig mismatch");

    assert_equal_qd(tag, "BlueScale tsml", q.tsml, ::ldexp(qd_real(1.0), static_cast<int>(q.exp_tsml)));
    assert_equal_qd(tag, "BlueScale tbig", q.tbig, ::ldexp(qd_real(1.0), static_cast<int>(q.exp_tbig)));
    assert_equal_qd(tag, "BlueScale ssml", q.ssml, ::ldexp(qd_real(1.0), static_cast<int>(q.exp_ssml)));
    assert_equal_qd(tag, "BlueScale sbig", q.sbig, ::ldexp(qd_real(1.0), static_cast<int>(q.exp_sbig)));

    const qd_real zero(0.0), one(1.0);
    qd_assert_case(q.tsml > zero, tag, "BlueScale: tsml must be positive");
    qd_assert_case(q.tsml < one, tag, "BlueScale: tsml must be < 1");
    qd_assert_case(q.tbig > one, tag, "BlueScale: tbig must be > 1");
    qd_assert_case(q.ssml > q.tbig, tag, "BlueScale: ssml must be > tbig");
    qd_assert_case(q.sbig > zero, tag, "BlueScale: sbig must be positive");
    qd_assert_case(q.sbig < q.tsml, tag, "BlueScale: sbig must be < tsml");

    const qd_real prod_ts = q.tsml * q.ssml;
    qd_assert_case(std::isfinite(to_double(prod_ts)) && to_double(prod_ts) > 0.0, tag, "BlueScale: tsml*ssml is not a positive finite value");
    const qd_real prod_bs = q.tbig * q.sbig;
    qd_assert_case(std::isfinite(to_double(prod_bs)) && to_double(prod_bs) > 0.0, tag, "BlueScale: tbig*sbig is not a positive finite value");
    const qd_real tsml_sq = q.tsml * q.tsml;
    qd_assert_case(to_double(tsml_sq) > 0.0, tag, "BlueScale: tsml^2 underflowed to zero");
    const qd_real tbig_sq = q.tbig * q.tbig;
    qd_assert_case(std::isfinite(to_double(tbig_sq)), tag, "BlueScale: tbig^2 overflowed");
    // Boundary classification probes specialized for exact Blue constants.
    // q.tsml and q.tbig are exact pure powers of two here, so their lower limbs are zero.
    // Use a 1-double-step probe in the lowest qd slot instead of relative (1 +/- delta).
    const double qd_low_up = std::nextafter(0.0, +std::numeric_limits<double>::infinity());
    const double qd_low_dn = std::nextafter(0.0, -std::numeric_limits<double>::infinity());

    const qd_real qd_probe_up(0.0, 0.0, 0.0, qd_low_up);
    const qd_real qd_probe_dn(0.0, 0.0, 0.0, qd_low_dn);

    const qd_real below_tsml = q.tsml + qd_probe_dn;
    const qd_real at_tsml = q.tsml;
    const qd_real above_tsml = q.tsml + qd_probe_up;

    const qd_real below_tbig = q.tbig + qd_probe_dn;
    const qd_real at_tbig = q.tbig;
    const qd_real above_tbig = q.tbig + qd_probe_up;

    const auto classify_blue = [&](const qd_real &x) -> int {
        const qd_real ax = (x < zero) ? -x : x;
        if (ax > q.tbig)
            return +1; // big
        if (ax < q.tsml)
            return -1; // small
        return 0;      // medium
    };

    qd_assert_case(classify_blue(below_tsml) == -1, tag, "BlueScale boundary: below_tsml must classify as small");
    qd_assert_case(classify_blue(at_tsml) == 0, tag, "BlueScale boundary: tsml must classify as medium");
    qd_assert_case(classify_blue(above_tsml) == 0, tag, "BlueScale boundary: above_tsml must classify as medium");

    qd_assert_case(classify_blue(below_tbig) == 0, tag, "BlueScale boundary: below_tbig must classify as medium");
    qd_assert_case(classify_blue(at_tbig) == 0, tag, "BlueScale boundary: tbig must classify as medium");
    qd_assert_case(classify_blue(above_tbig) == +1, tag, "BlueScale boundary: above_tbig must classify as big");

    qd_assert_case(classify_blue(-below_tsml) == -1, tag, "BlueScale boundary: negative below_tsml must classify as small");
    qd_assert_case(classify_blue(-above_tbig) == +1, tag, "BlueScale boundary: negative above_tbig must classify as big");

    qd_assert_case(classify_blue(qd_real(0.0)) == -1, tag, "BlueScale boundary: zero must classify as small");
    qd_assert_case(classify_blue(qd_real(1.0)) == 0, tag, "BlueScale boundary: one must classify as medium");

    const qd_real rescaled_small = below_tsml * q.ssml;
    const qd_real rescaled_big = above_tbig * q.sbig;
    qd_assert_case(classify_blue(rescaled_small) == 0, tag, "BlueScale boundary: below_tsml * ssml must classify as medium");
    qd_assert_case(classify_blue(rescaled_big) == 0, tag, "BlueScale boundary: above_tbig * sbig must classify as medium");

    const qd_real delta = Rlamch_qd("P");
    check_blue_threshold_boundaries_qd(tag, q, delta);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("BlueScale exp_tsml: %lld\n", (long long)q.exp_tsml);
        dual_printf("BlueScale exp_tbig: %lld\n", (long long)q.exp_tbig);
        dual_printf("BlueScale exp_ssml: %lld\n", (long long)q.exp_ssml);
        dual_printf("BlueScale exp_sbig: %lld\n", (long long)q.exp_sbig);

        sprintnum(_spbuf, q.tsml);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), q.tsml);
        dual_printf("BlueScale tsml:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.tbig);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), q.tbig);
        dual_printf("BlueScale tbig:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.ssml);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), q.ssml);
        dual_printf("BlueScale ssml:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.sbig);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), q.sbig);
        dual_printf("BlueScale sbig:     %40s    %40s\n", _spbuf, _sphexbuf);
    }
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
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        sprintnum(_spbuf, gotE);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), gotE);
        dual_printf("Rlamch E: Epsilon                      %70s    %70s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotS);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), gotS);
        dual_printf("Rlamch S: Safe minimum                 %70s    %70s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotB);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), gotB);
        dual_printf("Rlamch B: Base                         %70s    %70s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotP);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), gotP);
        dual_printf("Rlamch P: Precision                    %70s    %70s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotN);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), gotN);
        dual_printf("Rlamch N: Number of digits in mantissa %70s    %70s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotR);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), gotR);
        dual_printf("Rlamch R: Rounding mode                %70s    %70s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotM);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), gotM);
        dual_printf("Rlamch M: Minimum exponent             %70s    %70s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotU);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), gotU);
        dual_printf("Rlamch U: Underflow threshold          %70s    %70s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotL);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), gotL);
        dual_printf("Rlamch L: Largest exponent             %70s    %70s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotO);
        sprinthex_qd(_sphexbuf, sizeof(_sphexbuf), gotO);
        dual_printf("Rlamch O: Overflow threshold           %70s    %70s\n", _spbuf, _sphexbuf);
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
    check_arithmetic_params_qd(tag, print_values);
    check_lamch_qd_values(tag, print_values);
    check_blue_scaling_qd(tag, print_values);
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
// This test matches the current arithmetic_params-based implementation:
// - E and P are fixed to MPLAPACK's canonical DD literals
// - U and O are taken from dd_real::_min_normalized / dd_real::_max
// - N is taken from std::numeric_limits<dd_real>::digits
// - The exponent range matches IEEE-754 binary64 (double)
//
// IMPORTANT:
// Do NOT derive O via std::ldexp(1.0, 1024) (it overflows in double).
// Use dd_real::_max (as the implementation does).
// -----------------------------------------------------------------------------

constexpr int DD_MANTISSA_BITS = std::numeric_limits<dd_real>::digits;
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

static void check_arithmetic_params_dd(const char *tag, bool print_values) {
    const auto p = mplapack::get_arithmetic_params<dd_real>();
    const auto q = mplapack::get_blue_scaling_params<dd_real>();
    const auto q2 = mplapack::make_blue_scaling_params(p);

    assert_equal_dd(tag, "params.E", p.eps, Rlamch_dd("E"));
    assert_equal_dd(tag, "params.S", p.sfmin, Rlamch_dd("S"));
    assert_equal_dd(tag, "params.B", p.base, Rlamch_dd("B"));
    assert_equal_dd(tag, "params.P", p.prec, Rlamch_dd("P"));
    assert_equal_dd(tag, "params.R", p.rnd, Rlamch_dd("R"));
    assert_equal_dd(tag, "params.U", p.rmin, Rlamch_dd("U"));
    assert_equal_dd(tag, "params.O", p.rmax, Rlamch_dd("O"));

    assert_equal_dd(tag, "params.N", mplapack::detail::to_rlamch_real<dd_real>(p.t), Rlamch_dd("N"));
    assert_equal_dd(tag, "params.M", mplapack::detail::to_rlamch_real<dd_real>(p.emin), Rlamch_dd("M"));
    assert_equal_dd(tag, "params.L", mplapack::detail::to_rlamch_real<dd_real>(p.emax), Rlamch_dd("L"));

    assert_equal_dd(tag, "params.prec_consistency", p.prec, p.eps * p.base);
    assert_equal_dd(tag, "params.safmin", p.safmin, mplapack::detail::compute_safmin<dd_real>(p.emin, p.emax));
    assert_equal_dd(tag, "params.safmax", p.safmax, mplapack::detail::compute_safmax<dd_real>(p.emin, p.emax));

    dd_assert_case(q.exp_tsml == q2.exp_tsml, tag, "ArithmeticParams->Blue builder mismatch: exp_tsml");
    dd_assert_case(q.exp_tbig == q2.exp_tbig, tag, "ArithmeticParams->Blue builder mismatch: exp_tbig");
    dd_assert_case(q.exp_ssml == q2.exp_ssml, tag, "ArithmeticParams->Blue builder mismatch: exp_ssml");
    dd_assert_case(q.exp_sbig == q2.exp_sbig, tag, "ArithmeticParams->Blue builder mismatch: exp_sbig");

    assert_equal_dd(tag, "ArithmeticParams->Blue tsml", q.tsml, q2.tsml);
    assert_equal_dd(tag, "ArithmeticParams->Blue tbig", q.tbig, q2.tbig);
    assert_equal_dd(tag, "ArithmeticParams->Blue ssml", q.ssml, q2.ssml);
    assert_equal_dd(tag, "ArithmeticParams->Blue sbig", q.sbig, q2.sbig);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("[params/%s] t=%lld emin=%lld emax=%lld\n",
                    tag, (long long)p.t, (long long)p.emin, (long long)p.emax);

        sprintnum(_spbuf, p.safmin);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), p.safmin);
        dual_printf("[params/%s] safmin %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, p.safmax);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), p.safmax);
        dual_printf("[params/%s] safmax %40s    %40s\n", tag, _spbuf, _sphexbuf);

        dual_printf("[params/%s] builder exp_tsml=%lld exp_tbig=%lld exp_ssml=%lld exp_sbig=%lld\n",
                    tag,
                    (long long)q2.exp_tsml, (long long)q2.exp_tbig,
                    (long long)q2.exp_ssml, (long long)q2.exp_sbig);

        sprintnum(_spbuf, q2.tsml);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), q2.tsml);
        dual_printf("[params/%s] builder tsml   %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.tbig);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), q2.tbig);
        dual_printf("[params/%s] builder tbig   %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.ssml);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), q2.ssml);
        dual_printf("[params/%s] builder ssml   %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.sbig);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), q2.sbig);
        dual_printf("[params/%s] builder sbig   %40s    %40s\n", tag, _spbuf, _sphexbuf);
    }
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
    ex.E = dd_real(+0x1.fffffffffffffp-105, +0x0.0000000000000p+0000);
    ex.P = dd_real(+0x1.fffffffffffffp-104, +0x0.0000000000000p+0000);
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
    dd_assert_case((one / S) > zero, tag, "1/S is not positive");
    dd_assert_case(isfinite(one / S), tag, "1/S is not finite");
}

static void check_sfmin_inequalities_dd(const char *tag, const dd_real &S, const dd_real &U, const dd_real &O) {
    const dd_real one(1.0);

    const dd_real invS = one / S; // S = min_normalized, 1/S is in normal range: safe

    // Netlib-style sfmin inequalities
    dd_assert_case(invS <= O, tag, "1/S > O (violates safe minimum contract)");
    dd_assert_case(S >= U, tag, "S < U (violates sfmin >= rmin)");

    // 1/O is in the subnormal double range (~5.56e-309); dd_real division is
    // unreliable there. Use double arithmetic, which handles subnormals correctly.
    const double invO_d = 1.0 / to_double(O);
    dd_assert_case(to_double(S) >= invO_d, tag, "S < 1/O (violates sfmin >= 1/rmax)");
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

template <typename BlueQ> static int classify_blue_dd_value(const BlueQ &q, const dd_real &x) {
    dd_real ax = x;
    if (ax < dd_real(0.0))
        ax = -ax;
    if (ax > q.tbig)
        return +1;
    if (ax < q.tsml)
        return -1;
    return 0;
}

template <typename BlueQ> static void check_blue_threshold_boundaries_dd(const char *tag, const BlueQ &q, const dd_real &delta) {
    const dd_real zero(0.0), one(1.0), minus_one(-1.0);

    const dd_real below_tsml = q.tsml * (one - delta);
    const dd_real above_tsml = q.tsml * (one + delta);
    const dd_real below_tbig = q.tbig * (one - delta);
    const dd_real above_tbig = q.tbig * (one + delta);

    dd_assert_case(below_tsml < q.tsml, tag, "BlueScale boundary: below-tsml probe is not < tsml");
    dd_assert_case(above_tsml > q.tsml, tag, "BlueScale boundary: above-tsml probe is not > tsml");
    dd_assert_case(below_tbig < q.tbig, tag, "BlueScale boundary: below-tbig probe is not < tbig");
    dd_assert_case(above_tbig > q.tbig, tag, "BlueScale boundary: above-tbig probe is not > tbig");

    dd_assert_case(classify_blue_dd_value(q, below_tsml) == -1, tag, "BlueScale boundary: below-tsml probe must classify as small");
    dd_assert_case(classify_blue_dd_value(q, q.tsml) == 0, tag, "BlueScale boundary: tsml must classify as medium");
    dd_assert_case(classify_blue_dd_value(q, above_tsml) == 0, tag, "BlueScale boundary: above-tsml probe must classify as medium");
    dd_assert_case(classify_blue_dd_value(q, below_tbig) == 0, tag, "BlueScale boundary: below-tbig probe must classify as medium");
    dd_assert_case(classify_blue_dd_value(q, q.tbig) == 0, tag, "BlueScale boundary: tbig must classify as medium");
    dd_assert_case(classify_blue_dd_value(q, above_tbig) == +1, tag, "BlueScale boundary: above-tbig probe must classify as big");

    dd_assert_case(classify_blue_dd_value(q, minus_one * below_tsml) == -1, tag, "BlueScale boundary: negative below-tsml probe must classify as small");
    dd_assert_case(classify_blue_dd_value(q, minus_one * above_tbig) == +1, tag, "BlueScale boundary: negative above-tbig probe must classify as big");
    dd_assert_case(classify_blue_dd_value(q, zero) == -1, tag, "BlueScale boundary: 0 must classify as small");
    dd_assert_case(classify_blue_dd_value(q, one) == 0, tag, "BlueScale boundary: 1 must classify as medium");

    const dd_real scaled_small = below_tsml * q.ssml;
    const dd_real scaled_big = above_tbig * q.sbig;
    dd_assert_case(classify_blue_dd_value(q, scaled_small) == 0, tag, "BlueScale boundary: below-tsml probe * ssml must classify as medium");
    dd_assert_case(classify_blue_dd_value(q, scaled_big) == 0, tag, "BlueScale boundary: above-tbig probe * sbig must classify as medium");
}

static void check_blue_scaling_dd(const char *tag, bool print_values) {
    using mplapack::arithmetic_int;

    const auto q = mplapack::get_blue_scaling_params<dd_real>();

    int emin_int = 0;
    (void)std::frexp(dd_real::_min_normalized, &emin_int);
    const arithmetic_int emin = static_cast<arithmetic_int>(emin_int);
    const arithmetic_int emax = static_cast<arithmetic_int>(std::numeric_limits<double>::max_exponent);
    const arithmetic_int digits = static_cast<arithmetic_int>(std::numeric_limits<dd_real>::digits);

    dd_assert_case(q.exp_tsml == mplapack::detail::ceildiv2(emin - 1), tag, "BlueScale: exp_tsml mismatch");
    dd_assert_case(q.exp_tbig == mplapack::detail::floordiv2(emax - digits + 1), tag, "BlueScale: exp_tbig mismatch");
    dd_assert_case(q.exp_ssml == -mplapack::detail::floordiv2(emin - digits), tag, "BlueScale: exp_ssml mismatch");
    dd_assert_case(q.exp_sbig == -mplapack::detail::ceildiv2(emax + digits - 1), tag, "BlueScale: exp_sbig mismatch");

    assert_equal_dd(tag, "BlueScale tsml", q.tsml, ::ldexp(dd_real(1.0), static_cast<int>(q.exp_tsml)));
    assert_equal_dd(tag, "BlueScale tbig", q.tbig, ::ldexp(dd_real(1.0), static_cast<int>(q.exp_tbig)));
    assert_equal_dd(tag, "BlueScale ssml", q.ssml, ::ldexp(dd_real(1.0), static_cast<int>(q.exp_ssml)));
    assert_equal_dd(tag, "BlueScale sbig", q.sbig, ::ldexp(dd_real(1.0), static_cast<int>(q.exp_sbig)));

    const dd_real zero(0.0), one(1.0);
    dd_assert_case(q.tsml > zero, tag, "BlueScale: tsml must be positive");
    dd_assert_case(q.tsml < one, tag, "BlueScale: tsml must be < 1");
    dd_assert_case(q.tbig > one, tag, "BlueScale: tbig must be > 1");
    dd_assert_case(q.ssml > q.tbig, tag, "BlueScale: ssml must be > tbig");
    dd_assert_case(q.sbig > zero, tag, "BlueScale: sbig must be positive");
    dd_assert_case(q.sbig < q.tsml, tag, "BlueScale: sbig must be < tsml");

    const dd_real prod_ts = q.tsml * q.ssml;
    dd_assert_case(isfinite(prod_ts) && prod_ts > zero, tag, "BlueScale: tsml*ssml is not a positive finite value");
    const dd_real prod_bs = q.tbig * q.sbig;
    dd_assert_case(isfinite(prod_bs) && prod_bs > zero, tag, "BlueScale: tbig*sbig is not a positive finite value");
    const dd_real tsml_sq = q.tsml * q.tsml;
    dd_assert_case(tsml_sq > zero, tag, "BlueScale: tsml^2 underflowed to zero");
    const dd_real tbig_sq = q.tbig * q.tbig;
    dd_assert_case(isfinite(tbig_sq), tag, "BlueScale: tbig^2 overflowed");

    // Boundary classification probes specialized for exact Blue constants.
    // q.tsml and q.tbig are exact pure powers of two here, so their lower limbs are zero.
    // Use a 1-double-step probe in the lowest dd slot instead of relative (1 +/- delta).
    const double dd_low_up = std::nextafter(0.0, +std::numeric_limits<double>::infinity());
    const double dd_low_dn = std::nextafter(0.0, -std::numeric_limits<double>::infinity());

    const dd_real dd_probe_up(0.0, dd_low_up);
    const dd_real dd_probe_dn(0.0, dd_low_dn);

    const dd_real below_tsml = q.tsml + dd_probe_dn;
    const dd_real at_tsml = q.tsml;
    const dd_real above_tsml = q.tsml + dd_probe_up;

    const dd_real below_tbig = q.tbig + dd_probe_dn;
    const dd_real at_tbig = q.tbig;
    const dd_real above_tbig = q.tbig + dd_probe_up;

    const auto classify_blue = [&](const dd_real &x) -> int {
        const dd_real ax = (x < zero) ? -x : x;
        if (ax > q.tbig)
            return +1; // big
        if (ax < q.tsml)
            return -1; // small
        return 0;      // medium
    };

    dd_assert_case(classify_blue(below_tsml) == -1, tag, "BlueScale boundary: below_tsml must classify as small");
    dd_assert_case(classify_blue(at_tsml) == 0, tag, "BlueScale boundary: tsml must classify as medium");
    dd_assert_case(classify_blue(above_tsml) == 0, tag, "BlueScale boundary: above_tsml must classify as medium");

    dd_assert_case(classify_blue(below_tbig) == 0, tag, "BlueScale boundary: below_tbig must classify as medium");
    dd_assert_case(classify_blue(at_tbig) == 0, tag, "BlueScale boundary: tbig must classify as medium");
    dd_assert_case(classify_blue(above_tbig) == +1, tag, "BlueScale boundary: above_tbig must classify as big");

    dd_assert_case(classify_blue(-below_tsml) == -1, tag, "BlueScale boundary: negative below_tsml must classify as small");
    dd_assert_case(classify_blue(-above_tbig) == +1, tag, "BlueScale boundary: negative above_tbig must classify as big");

    dd_assert_case(classify_blue(dd_real(0.0)) == -1, tag, "BlueScale boundary: zero must classify as small");
    dd_assert_case(classify_blue(dd_real(1.0)) == 0, tag, "BlueScale boundary: one must classify as medium");

    const dd_real rescaled_small = below_tsml * q.ssml;
    const dd_real rescaled_big = above_tbig * q.sbig;
    dd_assert_case(classify_blue(rescaled_small) == 0, tag, "BlueScale boundary: below_tsml * ssml must classify as medium");
    dd_assert_case(classify_blue(rescaled_big) == 0, tag, "BlueScale boundary: above_tbig * sbig must classify as medium");

    const dd_real delta = Rlamch_dd("P");
    check_blue_threshold_boundaries_dd(tag, q, delta);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("BlueScale exp_tsml: %lld\n", (long long)q.exp_tsml);
        dual_printf("BlueScale exp_tbig: %lld\n", (long long)q.exp_tbig);
        dual_printf("BlueScale exp_ssml: %lld\n", (long long)q.exp_ssml);
        dual_printf("BlueScale exp_sbig: %lld\n", (long long)q.exp_sbig);

        sprintnum(_spbuf, q.tsml);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), q.tsml);
        dual_printf("BlueScale tsml:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.tbig);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), q.tbig);
        dual_printf("BlueScale tbig:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.ssml);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), q.ssml);
        dual_printf("BlueScale ssml:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.sbig);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), q.sbig);
        dual_printf("BlueScale sbig:     %40s    %40s\n", _spbuf, _sphexbuf);
    }
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
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        sprintnum(_spbuf, gotE);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), gotE);
        dual_printf("Rlamch E: Epsilon                      %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotS);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), gotS);
        dual_printf("Rlamch S: Safe minimum                 %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotB);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), gotB);
        dual_printf("Rlamch B: Base                         %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotP);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), gotP);
        dual_printf("Rlamch P: Precision                    %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotN);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), gotN);
        dual_printf("Rlamch N: Number of digits in mantissa %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotR);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), gotR);
        dual_printf("Rlamch R: Rounding mode                %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotM);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), gotM);
        dual_printf("Rlamch M: Minimum exponent             %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotU);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), gotU);
        dual_printf("Rlamch U: Underflow threshold          %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotL);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), gotL);
        dual_printf("Rlamch L: Largest exponent             %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotO);
        sprinthex_dd(_sphexbuf, sizeof(_sphexbuf), gotO);
        dual_printf("Rlamch O: Overflow threshold           %40s    %40s\n", _spbuf, _sphexbuf);
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
    check_arithmetic_params_dd(tag, print_values);
    check_lamch_dd_values(tag, print_values);
    check_blue_scaling_dd(tag, print_values);
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

static void check_arithmetic_params_double(const char *tag, bool print_values) {
    const auto p = mplapack::get_arithmetic_params<double>();
    const auto q = mplapack::get_blue_scaling_params<double>();
    const auto q2 = mplapack::make_blue_scaling_params(p);

    assert_equal_double(tag, "params.E", p.eps, Rlamch_double("E"));
    assert_equal_double(tag, "params.S", p.sfmin, Rlamch_double("S"));
    assert_equal_double(tag, "params.B", p.base, Rlamch_double("B"));
    assert_equal_double(tag, "params.P", p.prec, Rlamch_double("P"));
    assert_equal_double(tag, "params.R", p.rnd, Rlamch_double("R"));
    assert_equal_double(tag, "params.U", p.rmin, Rlamch_double("U"));
    assert_equal_double(tag, "params.O", p.rmax, Rlamch_double("O"));

    assert_equal_double(tag, "params.N", mplapack::detail::to_rlamch_real<double>(p.t), Rlamch_double("N"));
    assert_equal_double(tag, "params.M", mplapack::detail::to_rlamch_real<double>(p.emin), Rlamch_double("M"));
    assert_equal_double(tag, "params.L", mplapack::detail::to_rlamch_real<double>(p.emax), Rlamch_double("L"));

    assert_equal_double(tag, "params.prec_consistency", p.prec, p.eps * p.base);
    assert_equal_double(tag, "params.safmin", p.safmin, mplapack::detail::compute_safmin<double>(p.emin, p.emax));
    assert_equal_double(tag, "params.safmax", p.safmax, mplapack::detail::compute_safmax<double>(p.emin, p.emax));

    double_assert_case(q.exp_tsml == q2.exp_tsml, tag, "ArithmeticParams->Blue builder mismatch: exp_tsml");
    double_assert_case(q.exp_tbig == q2.exp_tbig, tag, "ArithmeticParams->Blue builder mismatch: exp_tbig");
    double_assert_case(q.exp_ssml == q2.exp_ssml, tag, "ArithmeticParams->Blue builder mismatch: exp_ssml");
    double_assert_case(q.exp_sbig == q2.exp_sbig, tag, "ArithmeticParams->Blue builder mismatch: exp_sbig");

    assert_equal_double(tag, "ArithmeticParams->Blue tsml", q.tsml, q2.tsml);
    assert_equal_double(tag, "ArithmeticParams->Blue tbig", q.tbig, q2.tbig);
    assert_equal_double(tag, "ArithmeticParams->Blue ssml", q.ssml, q2.ssml);
    assert_equal_double(tag, "ArithmeticParams->Blue sbig", q.sbig, q2.sbig);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("[params/%s] t=%lld emin=%lld emax=%lld\n",
                    tag, (long long)p.t, (long long)p.emin, (long long)p.emax);

        sprintnum(_spbuf, p.safmin);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), p.safmin);
        dual_printf("[params/%s] safmin %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, p.safmax);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), p.safmax);
        dual_printf("[params/%s] safmax %40s    %40s\n", tag, _spbuf, _sphexbuf);

        dual_printf("[params/%s] builder exp_tsml=%lld exp_tbig=%lld exp_ssml=%lld exp_sbig=%lld\n",
                    tag,
                    (long long)q2.exp_tsml, (long long)q2.exp_tbig,
                    (long long)q2.exp_ssml, (long long)q2.exp_sbig);

        sprintnum(_spbuf, q2.tsml);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), q2.tsml);
        dual_printf("[params/%s] builder tsml   %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.tbig);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), q2.tbig);
        dual_printf("[params/%s] builder tbig   %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.ssml);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), q2.ssml);
        dual_printf("[params/%s] builder ssml   %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.sbig);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), q2.sbig);
        dual_printf("[params/%s] builder sbig   %40s    %40s\n", tag, _spbuf, _sphexbuf);
    }
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

template <typename BlueQ> static int classify_blue_double_value(const BlueQ &q, double x) {
    const double ax = std::abs(x);
    if (ax > q.tbig)
        return +1;
    if (ax < q.tsml)
        return -1;
    return 0;
}

template <typename BlueQ> static void check_blue_threshold_boundaries_double(const char *tag, const BlueQ &q) {
    const double zero = 0.0;
    const double one = 1.0;
    const double below_tsml = std::nextafter(q.tsml, zero);
    const double above_tsml = std::nextafter(q.tsml, one);
    const double below_tbig = std::nextafter(q.tbig, zero);
    const double above_tbig = std::nextafter(q.tbig, std::numeric_limits<double>::infinity());

    double_assert_case(below_tsml < q.tsml, tag, "BlueScale boundary: nextbelow(tsml) is not < tsml");
    double_assert_case(above_tsml > q.tsml, tag, "BlueScale boundary: nextabove(tsml) is not > tsml");
    double_assert_case(below_tbig < q.tbig, tag, "BlueScale boundary: nextbelow(tbig) is not < tbig");
    double_assert_case(above_tbig > q.tbig, tag, "BlueScale boundary: nextabove(tbig) is not > tbig");

    double_assert_case(classify_blue_double_value(q, below_tsml) == -1, tag, "BlueScale boundary: nextbelow(tsml) must classify as small");
    double_assert_case(classify_blue_double_value(q, q.tsml) == 0, tag, "BlueScale boundary: tsml must classify as medium");
    double_assert_case(classify_blue_double_value(q, above_tsml) == 0, tag, "BlueScale boundary: nextabove(tsml) must classify as medium");
    double_assert_case(classify_blue_double_value(q, below_tbig) == 0, tag, "BlueScale boundary: nextbelow(tbig) must classify as medium");
    double_assert_case(classify_blue_double_value(q, q.tbig) == 0, tag, "BlueScale boundary: tbig must classify as medium");
    double_assert_case(classify_blue_double_value(q, above_tbig) == +1, tag, "BlueScale boundary: nextabove(tbig) must classify as big");

    double_assert_case(classify_blue_double_value(q, -below_tsml) == -1, tag, "BlueScale boundary: -nextbelow(tsml) must classify as small");
    double_assert_case(classify_blue_double_value(q, -above_tbig) == +1, tag, "BlueScale boundary: -nextabove(tbig) must classify as big");
    double_assert_case(classify_blue_double_value(q, zero) == -1, tag, "BlueScale boundary: 0 must classify as small");
    double_assert_case(classify_blue_double_value(q, one) == 0, tag, "BlueScale boundary: 1 must classify as medium");

    const double scaled_small = below_tsml * q.ssml;
    const double scaled_big = above_tbig * q.sbig;
    double_assert_case(classify_blue_double_value(q, scaled_small) == 0, tag, "BlueScale boundary: nextbelow(tsml) * ssml must classify as medium");
    double_assert_case(classify_blue_double_value(q, scaled_big) == 0, tag, "BlueScale boundary: nextabove(tbig) * sbig must classify as medium");
}

static void check_blue_scaling_double(const char *tag, bool print_values) {
    using mplapack::arithmetic_int;

    const auto q = mplapack::get_blue_scaling_params<double>();

    constexpr arithmetic_int emin = static_cast<arithmetic_int>(DBL_MIN_EXP);
    constexpr arithmetic_int emax = static_cast<arithmetic_int>(DBL_MAX_EXP);
    constexpr arithmetic_int digits = static_cast<arithmetic_int>(DBL_MANT_DIG);

    double_assert_case(q.exp_tsml == mplapack::detail::ceildiv2(emin - 1), tag, "BlueScale: exp_tsml mismatch");
    double_assert_case(q.exp_tbig == mplapack::detail::floordiv2(emax - digits + 1), tag, "BlueScale: exp_tbig mismatch");
    double_assert_case(q.exp_ssml == -mplapack::detail::floordiv2(emin - digits), tag, "BlueScale: exp_ssml mismatch");
    double_assert_case(q.exp_sbig == -mplapack::detail::ceildiv2(emax + digits - 1), tag, "BlueScale: exp_sbig mismatch");

    double_assert_case(q.exp_tsml == -511, tag, "BlueScale: exp_tsml must be -511 for binary64");
    double_assert_case(q.exp_tbig == 486, tag, "BlueScale: exp_tbig must be  486 for binary64");
    double_assert_case(q.exp_ssml == 537, tag, "BlueScale: exp_ssml must be  537 for binary64");
    double_assert_case(q.exp_sbig == -538, tag, "BlueScale: exp_sbig must be -538 for binary64");

    double_assert_case(q.tsml == std::ldexp(1.0, static_cast<int>(q.exp_tsml)), tag, "BlueScale: tsml value mismatch");
    double_assert_case(q.tbig == std::ldexp(1.0, static_cast<int>(q.exp_tbig)), tag, "BlueScale: tbig value mismatch");
    double_assert_case(q.ssml == std::ldexp(1.0, static_cast<int>(q.exp_ssml)), tag, "BlueScale: ssml value mismatch");
    double_assert_case(q.sbig == std::ldexp(1.0, static_cast<int>(q.exp_sbig)), tag, "BlueScale: sbig value mismatch");

    double_assert_case(q.tsml > 0.0, tag, "BlueScale: tsml must be positive");
    double_assert_case(q.tsml < 1.0, tag, "BlueScale: tsml must be < 1");
    double_assert_case(q.tbig > 1.0, tag, "BlueScale: tbig must be > 1");
    double_assert_case(q.ssml > q.tbig, tag, "BlueScale: ssml must be > tbig");
    double_assert_case(q.sbig > 0.0, tag, "BlueScale: sbig must be positive");
    double_assert_case(q.sbig < q.tsml, tag, "BlueScale: sbig must be < tsml");

    const double prod_ts = q.tsml * q.ssml;
    double_assert_case(std::isfinite(prod_ts) && prod_ts > 0.0, tag, "BlueScale: tsml*ssml is not a positive finite value");
    const double prod_bs = q.tbig * q.sbig;
    double_assert_case(std::isfinite(prod_bs) && prod_bs > 0.0, tag, "BlueScale: tbig*sbig is not a positive finite value");
    const double tsml_sq = q.tsml * q.tsml;
    double_assert_case(tsml_sq >= DBL_MIN, tag, "BlueScale: tsml^2 must be >= DBL_MIN");
    const double tbig_sq = q.tbig * q.tbig;
    double_assert_case(std::isfinite(tbig_sq) && tbig_sq <= DBL_MAX, tag, "BlueScale: tbig^2 overflowed");

    check_blue_threshold_boundaries_double(tag, q);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("BlueScale exp_tsml: %lld\n", (long long)q.exp_tsml);
        dual_printf("BlueScale exp_tbig: %lld\n", (long long)q.exp_tbig);
        dual_printf("BlueScale exp_ssml: %lld\n", (long long)q.exp_ssml);
        dual_printf("BlueScale exp_sbig: %lld\n", (long long)q.exp_sbig);

        sprintnum(_spbuf, q.tsml);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), q.tsml);
        dual_printf("BlueScale tsml:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.tbig);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), q.tbig);
        dual_printf("BlueScale tbig:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.ssml);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), q.ssml);
        dual_printf("BlueScale ssml:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.sbig);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), q.sbig);
        dual_printf("BlueScale sbig:     %40s    %40s\n", _spbuf, _sphexbuf);
    }
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
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        sprintnum(_spbuf, gotE);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), gotE);
        dual_printf("Rlamch E: Epsilon                      %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotS);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), gotS);
        dual_printf("Rlamch S: Safe minimum                 %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotB);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), gotB);
        dual_printf("Rlamch B: Base                         %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotP);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), gotP);
        dual_printf("Rlamch P: Precision                    %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotN);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), gotN);
        dual_printf("Rlamch N: Number of digits in mantissa %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotR);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), gotR);
        dual_printf("Rlamch R: Rounding mode                %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotM);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), gotM);
        dual_printf("Rlamch M: Minimum exponent             %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotU);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), gotU);
        dual_printf("Rlamch U: Underflow threshold          %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotL);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), gotL);
        dual_printf("Rlamch L: Largest exponent             %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotO);
        sprinthex_double(_sphexbuf, sizeof(_sphexbuf), gotO);
        dual_printf("Rlamch O: Overflow threshold           %40s    %40s\n", _spbuf, _sphexbuf);
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
    check_arithmetic_params_double(tag, print_values);
    check_lamch_double_values(tag, print_values);
    check_blue_scaling_double(tag, print_values);
}

#endif // ___MPLAPACK_BUILD_WITH_DOUBLE___

#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
#include <cfenv> // std::fegetround, FE_TONEAREST
#include <limits>

template <typename BlueQ> static int classify_blue_binary128_value(const BlueQ &q, mplapack_binary128_t x) {
    const mplapack_binary128_t zero = (mplapack_binary128_t)0.0;
    if (x < zero)
        x = -x;
    if (x > q.tbig)
        return +1;
    if (x < q.tsml)
        return -1;
    return 0;
}

template <typename BlueQ> static void check_blue_threshold_boundaries_binary128(const char *tag, const BlueQ &q) {
    const mplapack_binary128_t zero = (mplapack_binary128_t)0.0;
    const mplapack_binary128_t one = (mplapack_binary128_t)1.0;
    const mplapack_binary128_t minus_one = (mplapack_binary128_t)-1.0;
    const mplapack_binary128_t delta = Rlamch_binary128("P");

    const mplapack_binary128_t below_tsml = q.tsml * (one - delta);
    const mplapack_binary128_t above_tsml = q.tsml * (one + delta);
    const mplapack_binary128_t below_tbig = q.tbig * (one - delta);
    const mplapack_binary128_t above_tbig = q.tbig * (one + delta);

    auto assert_case = [&](bool cond, const char *what) {
        if (cond)
            return;
        char buf[256];
        snprintf(buf, sizeof(buf), "%s: %s", tag, what);
        printf("*** Testing Mutils (binary128) BlueScale failed: %s ***\n", buf);
        exit(1);
    };

    assert_case(below_tsml < q.tsml, "BlueScale boundary: below-tsml probe is not < tsml");
    assert_case(above_tsml > q.tsml, "BlueScale boundary: above-tsml probe is not > tsml");
    assert_case(below_tbig < q.tbig, "BlueScale boundary: below-tbig probe is not < tbig");
    assert_case(above_tbig > q.tbig, "BlueScale boundary: above-tbig probe is not > tbig");

    assert_case(classify_blue_binary128_value(q, below_tsml) == -1, "BlueScale boundary: below-tsml probe must classify as small");
    assert_case(classify_blue_binary128_value(q, q.tsml) == 0, "BlueScale boundary: tsml must classify as medium");
    assert_case(classify_blue_binary128_value(q, above_tsml) == 0, "BlueScale boundary: above-tsml probe must classify as medium");
    assert_case(classify_blue_binary128_value(q, below_tbig) == 0, "BlueScale boundary: below-tbig probe must classify as medium");
    assert_case(classify_blue_binary128_value(q, q.tbig) == 0, "BlueScale boundary: tbig must classify as medium");
    assert_case(classify_blue_binary128_value(q, above_tbig) == +1, "BlueScale boundary: above-tbig probe must classify as big");

    assert_case(classify_blue_binary128_value(q, minus_one * below_tsml) == -1, "BlueScale boundary: negative below-tsml probe must classify as small");
    assert_case(classify_blue_binary128_value(q, minus_one * above_tbig) == +1, "BlueScale boundary: negative above-tbig probe must classify as big");
    assert_case(classify_blue_binary128_value(q, zero) == -1, "BlueScale boundary: 0 must classify as small");
    assert_case(classify_blue_binary128_value(q, one) == 0, "BlueScale boundary: 1 must classify as medium");

    const mplapack_binary128_t scaled_small = below_tsml * q.ssml;
    const mplapack_binary128_t scaled_big = above_tbig * q.sbig;
    assert_case(classify_blue_binary128_value(q, scaled_small) == 0, "BlueScale boundary: below-tsml probe * ssml must classify as medium");
    assert_case(classify_blue_binary128_value(q, scaled_big) == 0, "BlueScale boundary: above-tbig probe * sbig must classify as medium");
}

static void check_blue_scaling_binary128(const char *tag, int emin, int emax, int p, bool print_values) {
    using mplapack::arithmetic_int;
    const mplapack_binary128_t zero = (mplapack_binary128_t)0.0;
    const mplapack_binary128_t one = (mplapack_binary128_t)1.0;

    auto fail = [&](const char *what) {
        printf("*** Testing Mutils (binary128) BlueScale failed: %s: %s ***\n", tag, what);
        exit(1);
    };
    auto assert_case = [&](bool cond, const char *what) {
        if (!cond) {
            char buf[256];
            snprintf(buf, sizeof(buf), "%s: %s", tag, what);
            fail(buf);
        }
    };

    const auto q = mplapack::get_blue_scaling_params<mplapack_binary128_t>();

    const arithmetic_int ai_emin = static_cast<arithmetic_int>(emin);
    const arithmetic_int ai_emax = static_cast<arithmetic_int>(emax);
    const arithmetic_int ai_digits = static_cast<arithmetic_int>(p);

    assert_case(q.exp_tsml == mplapack::detail::ceildiv2(ai_emin - 1), "BlueScale: exp_tsml mismatch");
    assert_case(q.exp_tbig == mplapack::detail::floordiv2(ai_emax - ai_digits + 1), "BlueScale: exp_tbig mismatch");
    assert_case(q.exp_ssml == -mplapack::detail::floordiv2(ai_emin - ai_digits), "BlueScale: exp_ssml mismatch");
    assert_case(q.exp_sbig == -mplapack::detail::ceildiv2(ai_emax + ai_digits - 1), "BlueScale: exp_sbig mismatch");

    assert_case(q.exp_tsml == -8191LL, "BlueScale: exp_tsml must be -8191 for binary128");
    assert_case(q.exp_tbig == 8136LL, "BlueScale: exp_tbig must be  8136 for binary128");
    assert_case(q.exp_ssml == 8247LL, "BlueScale: exp_ssml must be  8247 for binary128");
    assert_case(q.exp_sbig == -8248LL, "BlueScale: exp_sbig must be -8248 for binary128");

    assert_case(q.tsml == ldexp(one, static_cast<int>(q.exp_tsml)), "BlueScale: tsml value mismatch");
    assert_case(q.tbig == ldexp(one, static_cast<int>(q.exp_tbig)), "BlueScale: tbig value mismatch");
    assert_case(q.ssml == ldexp(one, static_cast<int>(q.exp_ssml)), "BlueScale: ssml value mismatch");
    assert_case(q.sbig == ldexp(one, static_cast<int>(q.exp_sbig)), "BlueScale: sbig value mismatch");

    assert_case(q.tsml > zero, "BlueScale: tsml must be positive");
    assert_case(q.tsml < one, "BlueScale: tsml must be < 1");
    assert_case(q.tbig > one, "BlueScale: tbig must be > 1");
    assert_case(q.ssml > q.tbig, "BlueScale: ssml must be > tbig");
    assert_case(q.sbig > zero, "BlueScale: sbig must be positive");
    assert_case(q.tsml < q.tbig, "BlueScale: tsml must be < tbig");
    assert_case(q.sbig < q.tsml, "BlueScale: sbig must be < tsml");

    const mplapack_binary128_t prod_ts = q.tsml * q.ssml;
    assert_case(__builtin_isfinite(prod_ts) && prod_ts > zero, "BlueScale: tsml*ssml must be a positive finite value");
    const mplapack_binary128_t prod_bs = q.tbig * q.sbig;
    assert_case(__builtin_isfinite(prod_bs) && prod_bs > zero, "BlueScale: tbig*sbig must be a positive finite value");
    const mplapack_binary128_t tsml_sq = q.tsml * q.tsml;
    assert_case(tsml_sq > zero, "BlueScale: tsml^2 underflowed to zero");
    const mplapack_binary128_t tbig_sq = q.tbig * q.tbig;
    assert_case(__builtin_isfinite(tbig_sq), "BlueScale: tbig^2 overflowed");

    check_blue_threshold_boundaries_binary128(tag, q);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("BlueScale exp_tsml: %lld\n", (long long)q.exp_tsml);
        dual_printf("BlueScale exp_tbig: %lld\n", (long long)q.exp_tbig);
        dual_printf("BlueScale exp_ssml: %lld\n", (long long)q.exp_ssml);
        dual_printf("BlueScale exp_sbig: %lld\n", (long long)q.exp_sbig);

        sprintnum(_spbuf, q.tsml);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), q.tsml);
        dual_printf("BlueScale tsml:     %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.tbig);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), q.tbig);
        dual_printf("BlueScale tbig:     %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.ssml);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), q.ssml);
        dual_printf("BlueScale ssml:     %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.sbig);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), q.sbig);
        dual_printf("BlueScale sbig:     %50s    %50s\n", _spbuf, _sphexbuf);
    }
}

static void check_arithmetic_params_binary128(const char *tag, bool print_values) {
    auto fail = [&](const char *what) {
        printf("*** Testing Mutils (binary128) failed: %s: %s ***\n", tag, what);
        exit(1);
    };

    auto assert_case = [&](bool cond, const char *what) {
        if (!cond) {
            fail(what);
        }
    };

    auto assert_equal = [&](const char *name, mplapack_binary128_t got, mplapack_binary128_t expected) {
        if (got == expected)
            return;

        printf("*** Testing Mutils (binary128) failed: %s mismatch in %s ***\n", tag, name);
        printf("    got      = ");
        printnum(got);
        printf("\n");
        printf("    expected = ");
        printnum(expected);
        printf("\n");
        exit(1);
    };

    const auto p = mplapack::get_arithmetic_params<mplapack_binary128_t>();
    const auto q = mplapack::get_blue_scaling_params<mplapack_binary128_t>();
    const auto q2 = mplapack::make_blue_scaling_params(p);

    assert_equal("params.E", p.eps, Rlamch_binary128("E"));
    assert_equal("params.S", p.sfmin, Rlamch_binary128("S"));
    assert_equal("params.B", p.base, Rlamch_binary128("B"));
    assert_equal("params.P", p.prec, Rlamch_binary128("P"));
    assert_equal("params.R", p.rnd, Rlamch_binary128("R"));
    assert_equal("params.U", p.rmin, Rlamch_binary128("U"));
    assert_equal("params.O", p.rmax, Rlamch_binary128("O"));

    assert_equal("params.N", mplapack::detail::to_rlamch_real<mplapack_binary128_t>(p.t), Rlamch_binary128("N"));
    assert_equal("params.M", mplapack::detail::to_rlamch_real<mplapack_binary128_t>(p.emin), Rlamch_binary128("M"));
    assert_equal("params.L", mplapack::detail::to_rlamch_real<mplapack_binary128_t>(p.emax), Rlamch_binary128("L"));

    assert_equal("params.prec_consistency", p.prec, p.eps * p.base);
    assert_equal("params.safmin", p.safmin, mplapack::detail::compute_safmin<mplapack_binary128_t>(p.emin, p.emax));
    assert_equal("params.safmax", p.safmax, mplapack::detail::compute_safmax<mplapack_binary128_t>(p.emin, p.emax));

    assert_case(q.exp_tsml == q2.exp_tsml, "ArithmeticParams->Blue builder mismatch: exp_tsml");
    assert_case(q.exp_tbig == q2.exp_tbig, "ArithmeticParams->Blue builder mismatch: exp_tbig");
    assert_case(q.exp_ssml == q2.exp_ssml, "ArithmeticParams->Blue builder mismatch: exp_ssml");
    assert_case(q.exp_sbig == q2.exp_sbig, "ArithmeticParams->Blue builder mismatch: exp_sbig");

    assert_equal("ArithmeticParams->Blue tsml", q.tsml, q2.tsml);
    assert_equal("ArithmeticParams->Blue tbig", q.tbig, q2.tbig);
    assert_equal("ArithmeticParams->Blue ssml", q.ssml, q2.ssml);
    assert_equal("ArithmeticParams->Blue sbig", q.sbig, q2.sbig);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("[params/%s] t=%lld emin=%lld emax=%lld\n",
                    tag, (long long)p.t, (long long)p.emin, (long long)p.emax);

        sprintnum(_spbuf, p.safmin);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), p.safmin);
        dual_printf("[params/%s] safmin %50s    %50s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, p.safmax);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), p.safmax);
        dual_printf("[params/%s] safmax %50s    %50s\n", tag, _spbuf, _sphexbuf);

        dual_printf("[params/%s] builder exp_tsml=%lld exp_tbig=%lld exp_ssml=%lld exp_sbig=%lld\n",
                    tag,
                    (long long)q2.exp_tsml, (long long)q2.exp_tbig,
                    (long long)q2.exp_ssml, (long long)q2.exp_sbig);

        sprintnum(_spbuf, q2.tsml);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), q2.tsml);
        dual_printf("[params/%s] builder tsml   %50s    %50s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.tbig);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), q2.tbig);
        dual_printf("[params/%s] builder tbig   %50s    %50s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.ssml);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), q2.ssml);
        dual_printf("[params/%s] builder ssml   %50s    %50s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.sbig);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), q2.sbig);
        dual_printf("[params/%s] builder sbig   %50s    %50s\n", tag, _spbuf, _sphexbuf);
    }
}

void Rlamch_binary128_test() {
#if defined VERBOSE_TEST
    const bool print_values = true;
#else
    const bool print_values = false;
#endif

    const char *tag = "binary128";
    check_arithmetic_params_binary128(tag, print_values);

    auto fail = [&](const char *what) {
        printf("*** Testing Mutils (binary128) failed: %s ***\n", what);
        exit(1);
    };

    auto assert_case = [&](bool cond, const char *what) {
        if (cond)
            return;
        char buf[256];
        snprintf(buf, sizeof(buf), "%s: %s", tag, what);
        fail(buf);
    };

    auto assert_equal = [&](const char *name, mplapack_binary128_t got, mplapack_binary128_t expected) {
        if (got == expected)
            return;

        printf("*** Testing Mutils (binary128) failed: %s mismatch in %s ***\n", tag, name);
        printf("    got      = ");
        printnum(got);
        printf("\n");
        printf("    expected = ");
        printnum(expected);
        printf("\n");
        exit(1);
    };

    const mplapack_binary128_t zero = (mplapack_binary128_t)0.0;
    const mplapack_binary128_t one = (mplapack_binary128_t)1.0;
    const mplapack_binary128_t two = (mplapack_binary128_t)2.0;

// Expected values for IEEE-754 binary128 based on compiler-provided parameters.
#if MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_FLOAT128
// Using _Float128
#if defined(__FLT128_MANT_DIG__)
    const int p = (int)__FLT128_MANT_DIG__;
    const int emin = (int)__FLT128_MIN_EXP__;
    const int emax = (int)__FLT128_MAX_EXP__;
#else
#error "_Float128 mode selected but __FLT128_MANT_DIG__ not available"
#endif
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_QUADMATH
// Using __float128 (GCC quadmath extension)
#if defined(__FLT128_MANT_DIG__)
    const int p = (int)__FLT128_MANT_DIG__;
    const int emin = (int)__FLT128_MIN_EXP__;
    const int emax = (int)__FLT128_MAX_EXP__;
#else
    const int p = (int)FLT128_MANT_DIG; // defined in quadmath.h
    const int emin = (int)FLT128_MIN_EXP;
    const int emax = (int)FLT128_MAX_EXP;
#endif
#elif MPLAPACK_BINARY128_MODE == MPLAPACK_BINARY128_MODE_LDBL
// Using long double (must be binary128)
#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 113
    const int p = (int)__LDBL_MANT_DIG__;
    const int emin = (int)__LDBL_MIN_EXP__;
    const int emax = (int)__LDBL_MAX_EXP__;
#else
#error "long double mode selected but long double is not binary128 (113-bit mantissa required)"
#endif
#else
#error "Invalid or disabled MPLAPACK_BINARY128_MODE"
#endif

    const mplapack_binary128_t exB = two;
    const mplapack_binary128_t exN = (mplapack_binary128_t)p;

    // ulp(1) = 2^(1-p), unit roundoff E = ulp(1)/2.
    const mplapack_binary128_t exP = ldexp(one, 1 - p);
    const mplapack_binary128_t exE = exP / two;

    // DLAMCH-style exponents: rmin = 2^(emin-1), rmax = (1-E)*2^emax.
    const mplapack_binary128_t exM = (mplapack_binary128_t)emin;
    const mplapack_binary128_t exL = (mplapack_binary128_t)emax;
    const mplapack_binary128_t exU = ldexp(one, emin - 1);
    const mplapack_binary128_t exO = ldexp(one - ldexp(one, -p), emax);

    // Safe minimum: max(rmin, (1/rmax)*(1+E)).
    const mplapack_binary128_t small = one / exO;
    const mplapack_binary128_t candidate = small * (one + exE);
    const mplapack_binary128_t exS = (candidate >= exU) ? candidate : exU;

    const mplapack_binary128_t exR = one;
    const mplapack_binary128_t exZ = zero;

    // Fetch actual values from Rlamch_mplapack_binary128_t.
    const mplapack_binary128_t gotE = Rlamch_binary128("E");
    const mplapack_binary128_t gotS = Rlamch_binary128("S");
    const mplapack_binary128_t gotB = Rlamch_binary128("B");
    const mplapack_binary128_t gotP = Rlamch_binary128("P");
    const mplapack_binary128_t gotN = Rlamch_binary128("N");
    const mplapack_binary128_t gotR = Rlamch_binary128("R");
    const mplapack_binary128_t gotM = Rlamch_binary128("M");
    const mplapack_binary128_t gotU = Rlamch_binary128("U");
    const mplapack_binary128_t gotL = Rlamch_binary128("L");
    const mplapack_binary128_t gotO = Rlamch_binary128("O");
    const mplapack_binary128_t gotZ = Rlamch_binary128("Z");

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
    const mplapack_binary128_t next = nextafter(one, two);
    assert_case((next - one) == gotP, "P is not equal to ulp(1) = nextafter(1,2)-1");
    assert_case((gotE * two) == gotP, "expected P == 2*E");

    // Strong, rounding-mode-independent check for 1 + P
    volatile mplapack_binary128_t b = one + gotP;
    assert_case(b == next, "expected fl(1 + P) == nextafter(1,2)");
    assert_case(b > one, "expected fl(1 + 2E) > 1"); // this is effectivly covered by b==next.

    // Tie check: only assert when FE_TONEAREST
    if (std::fegetround() == FE_TONEAREST) {
        volatile mplapack_binary128_t a = one + gotE;
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
    const mplapack_binary128_t got_prod = gotU * gotO;
    const mplapack_binary128_t expected_prod = ldexp(one - gotE, emin + emax - 1);
    assert_case(got_prod == expected_prod, "U*O cross-check failed (inconsistent rmin/rmax model)");

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        sprintnum(_spbuf, gotE);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), gotE);
        dual_printf("Rlamch E: Epsilon                      %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotS);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), gotS);
        dual_printf("Rlamch S: Safe minimum                 %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotB);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), gotB);
        dual_printf("Rlamch B: Base                         %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotP);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), gotP);
        dual_printf("Rlamch P: Precision                    %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotN);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), gotN);
        dual_printf("Rlamch N: Number of digits in mantissa %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotR);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), gotR);
        dual_printf("Rlamch R: Rounding mode                %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotM);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), gotM);
        dual_printf("Rlamch M: Minimum exponent             %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotU);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), gotU);
        dual_printf("Rlamch U: Underflow threshold          %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotL);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), gotL);
        dual_printf("Rlamch L: Largest exponent             %50s    %50s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotO);
        sprinthex_binary128_fixed(_sphexbuf, sizeof(_sphexbuf), gotO);
        dual_printf("Rlamch O: Overflow threshold           %50s    %50s\n", _spbuf, _sphexbuf);
    }

    check_blue_scaling_binary128(tag, emin, emax, p, print_values);
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
#include <cfenv> // std::fegetround, FE_TONEAREST
#include <cmath> // scalbnl

template <typename BlueQ> static int classify_blue_binary80_value(const BlueQ &q, mplapack_binary80_t x) {
    const mplapack_binary80_t zero = (mplapack_binary80_t)0.0;
    if (x < zero)
        x = -x;
    if (x > q.tbig)
        return +1;
    if (x < q.tsml)
        return -1;
    return 0;
}

template <typename BlueQ> static void check_blue_threshold_boundaries_binary80(const char *tag, const BlueQ &q) {
    const mplapack_binary80_t zero = (mplapack_binary80_t)0.0;
    const mplapack_binary80_t one = (mplapack_binary80_t)1.0;
    const mplapack_binary80_t minus_one = (mplapack_binary80_t)-1.0;
    const mplapack_binary80_t delta = Rlamch_binary80("P");

    const mplapack_binary80_t below_tsml = q.tsml * (one - delta);
    const mplapack_binary80_t above_tsml = q.tsml * (one + delta);
    const mplapack_binary80_t below_tbig = q.tbig * (one - delta);
    const mplapack_binary80_t above_tbig = q.tbig * (one + delta);

    auto assert_case = [&](bool cond, const char *what) {
        if (cond)
            return;
        char buf[256];
        snprintf(buf, sizeof(buf), "%s: %s", tag, what);
        printf("*** Testing Mutils (mplapack_binary80_t) BlueScale failed: %s ***\n", buf);
        exit(1);
    };

    assert_case(below_tsml < q.tsml, "BlueScale boundary: below-tsml probe is not < tsml");
    assert_case(above_tsml > q.tsml, "BlueScale boundary: above-tsml probe is not > tsml");
    assert_case(below_tbig < q.tbig, "BlueScale boundary: below-tbig probe is not < tbig");
    assert_case(above_tbig > q.tbig, "BlueScale boundary: above-tbig probe is not > tbig");

    assert_case(classify_blue_binary80_value(q, below_tsml) == -1, "BlueScale boundary: below-tsml probe must classify as small");
    assert_case(classify_blue_binary80_value(q, q.tsml) == 0, "BlueScale boundary: tsml must classify as medium");
    assert_case(classify_blue_binary80_value(q, above_tsml) == 0, "BlueScale boundary: above-tsml probe must classify as medium");
    assert_case(classify_blue_binary80_value(q, below_tbig) == 0, "BlueScale boundary: below-tbig probe must classify as medium");
    assert_case(classify_blue_binary80_value(q, q.tbig) == 0, "BlueScale boundary: tbig must classify as medium");
    assert_case(classify_blue_binary80_value(q, above_tbig) == +1, "BlueScale boundary: above-tbig probe must classify as big");

    assert_case(classify_blue_binary80_value(q, minus_one * below_tsml) == -1, "BlueScale boundary: negative below-tsml probe must classify as small");
    assert_case(classify_blue_binary80_value(q, minus_one * above_tbig) == +1, "BlueScale boundary: negative above-tbig probe must classify as big");
    assert_case(classify_blue_binary80_value(q, zero) == -1, "BlueScale boundary: 0 must classify as small");
    assert_case(classify_blue_binary80_value(q, one) == 0, "BlueScale boundary: 1 must classify as medium");

    const mplapack_binary80_t scaled_small = below_tsml * q.ssml;
    const mplapack_binary80_t scaled_big = above_tbig * q.sbig;
    assert_case(classify_blue_binary80_value(q, scaled_small) == 0, "BlueScale boundary: below-tsml probe * ssml must classify as medium");
    assert_case(classify_blue_binary80_value(q, scaled_big) == 0, "BlueScale boundary: above-tbig probe * sbig must classify as medium");
}

static void check_blue_scaling_binary80(const char *tag, int emin, int emax, int p, bool print_values) {
    using mplapack::arithmetic_int;
    const mplapack_binary80_t zero = (mplapack_binary80_t)0.0;
    const mplapack_binary80_t one = (mplapack_binary80_t)1.0;

    auto fail = [&](const char *what) {
        printf("*** Testing Mutils (mplapack_binary80_t) BlueScale failed: %s: %s ***\n", tag, what);
        exit(1);
    };
    auto assert_case = [&](bool cond, const char *what) {
        if (!cond) {
            char buf[256];
            snprintf(buf, sizeof(buf), "%s: %s", tag, what);
            fail(buf);
        }
    };

    const auto q = mplapack::get_blue_scaling_params<mplapack_binary80_t>();

    const arithmetic_int ai_emin = static_cast<arithmetic_int>(emin);
    const arithmetic_int ai_emax = static_cast<arithmetic_int>(emax);
    const arithmetic_int ai_digits = static_cast<arithmetic_int>(p);

    assert_case(q.exp_tsml == mplapack::detail::ceildiv2(ai_emin - 1), "BlueScale: exp_tsml mismatch");
    assert_case(q.exp_tbig == mplapack::detail::floordiv2(ai_emax - ai_digits + 1), "BlueScale: exp_tbig mismatch");
    assert_case(q.exp_ssml == -mplapack::detail::floordiv2(ai_emin - ai_digits), "BlueScale: exp_ssml mismatch");
    assert_case(q.exp_sbig == -mplapack::detail::ceildiv2(ai_emax + ai_digits - 1), "BlueScale: exp_sbig mismatch");

    assert_case(q.exp_tsml == -8191LL, "BlueScale: exp_tsml must be -8191 for binary80");
    assert_case(q.exp_tbig == 8160LL, "BlueScale: exp_tbig must be  8160 for binary80");
    assert_case(q.exp_ssml == 8223LL, "BlueScale: exp_ssml must be  8223 for binary80");
    assert_case(q.exp_sbig == -8224LL, "BlueScale: exp_sbig must be -8224 for binary80");

    assert_case(q.tsml == scalbnl(one, static_cast<int>(q.exp_tsml)), "BlueScale: tsml value mismatch");
    assert_case(q.tbig == scalbnl(one, static_cast<int>(q.exp_tbig)), "BlueScale: tbig value mismatch");
    assert_case(q.ssml == scalbnl(one, static_cast<int>(q.exp_ssml)), "BlueScale: ssml value mismatch");
    assert_case(q.sbig == scalbnl(one, static_cast<int>(q.exp_sbig)), "BlueScale: sbig value mismatch");

    assert_case(q.tsml > zero, "BlueScale: tsml must be positive");
    assert_case(q.tsml < one, "BlueScale: tsml must be < 1");
    assert_case(q.tbig > one, "BlueScale: tbig must be > 1");
    assert_case(q.ssml > q.tbig, "BlueScale: ssml must be > tbig");
    assert_case(q.tsml < q.tbig, "BlueScale: tsml must be < tbig");
    assert_case(q.exp_sbig < q.exp_tsml, "BlueScale: exp_sbig must be < exp_tsml");
    assert_case(q.exp_tsml < q.exp_tbig, "BlueScale: exp_tsml must be < exp_tbig");
    assert_case(q.exp_tbig < q.exp_ssml, "BlueScale: exp_tbig must be < exp_ssml");
    assert_case(q.sbig > zero, "BlueScale: sbig must be positive");
    assert_case(q.sbig < q.tsml, "BlueScale: sbig must be < tsml");

    const mplapack_binary80_t prod_ts = q.tsml * q.ssml;
    assert_case(__builtin_isfinite(prod_ts) && prod_ts > zero, "BlueScale: tsml*ssml must be a positive finite value");
    const mplapack_binary80_t prod_bs = q.tbig * q.sbig;
    assert_case(__builtin_isfinite(prod_bs) && prod_bs > zero, "BlueScale: tbig*sbig must be a positive finite value");
    const mplapack_binary80_t tsml_sq = q.tsml * q.tsml;
    assert_case(tsml_sq > zero, "BlueScale: tsml^2 underflowed to zero");
    const mplapack_binary80_t tbig_sq = q.tbig * q.tbig;
    assert_case(__builtin_isfinite(tbig_sq), "BlueScale: tbig^2 overflowed");

    check_blue_threshold_boundaries_binary80(tag, q);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("BlueScale exp_tsml: %lld\n", (long long)q.exp_tsml);
        dual_printf("BlueScale exp_tbig: %lld\n", (long long)q.exp_tbig);
        dual_printf("BlueScale exp_ssml: %lld\n", (long long)q.exp_ssml);
        dual_printf("BlueScale exp_sbig: %lld\n", (long long)q.exp_sbig);

        sprintnum(_spbuf, q.tsml);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), q.tsml);
        dual_printf("BlueScale tsml:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.tbig);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), q.tbig);
        dual_printf("BlueScale tbig:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.ssml);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), q.ssml);
        dual_printf("BlueScale ssml:     %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, q.sbig);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), q.sbig);
        dual_printf("BlueScale sbig:     %40s    %40s\n", _spbuf, _sphexbuf);
    }
}

static void check_arithmetic_params_binary80(const char *tag, bool print_values) {
    auto fail = [&](const char *what) {
        printf("*** Testing Mutils (mplapack_binary80_t) failed: %s: %s ***\n", tag, what);
        exit(1);
    };

    auto assert_case = [&](bool cond, const char *what) {
        if (!cond) {
            fail(what);
        }
    };

    auto assert_equal = [&](const char *name, mplapack_binary80_t got, mplapack_binary80_t expected) {
        if (got == expected)
            return;

        printf("*** Testing Mutils (mplapack_binary80_t) failed: %s mismatch in %s ***\n", tag, name);
        printf("    got      = ");
        printnum(got);
        printf("\n");
        printf("    expected = ");
        printnum(expected);
        printf("\n");
        exit(1);
    };

    const auto p = mplapack::get_arithmetic_params<mplapack_binary80_t>();
    const auto q = mplapack::get_blue_scaling_params<mplapack_binary80_t>();
    const auto q2 = mplapack::make_blue_scaling_params(p);

    assert_equal("params.E", p.eps, Rlamch_binary80("E"));
    assert_equal("params.S", p.sfmin, Rlamch_binary80("S"));
    assert_equal("params.B", p.base, Rlamch_binary80("B"));
    assert_equal("params.P", p.prec, Rlamch_binary80("P"));
    assert_equal("params.R", p.rnd, Rlamch_binary80("R"));
    assert_equal("params.U", p.rmin, Rlamch_binary80("U"));
    assert_equal("params.O", p.rmax, Rlamch_binary80("O"));

    assert_equal("params.N", mplapack::detail::to_rlamch_real<mplapack_binary80_t>(p.t), Rlamch_binary80("N"));
    assert_equal("params.M", mplapack::detail::to_rlamch_real<mplapack_binary80_t>(p.emin), Rlamch_binary80("M"));
    assert_equal("params.L", mplapack::detail::to_rlamch_real<mplapack_binary80_t>(p.emax), Rlamch_binary80("L"));

    assert_equal("params.prec_consistency", p.prec, p.eps * p.base);
    assert_equal("params.safmin", p.safmin, mplapack::detail::compute_safmin<mplapack_binary80_t>(p.emin, p.emax));
    assert_equal("params.safmax", p.safmax, mplapack::detail::compute_safmax<mplapack_binary80_t>(p.emin, p.emax));

    assert_case(q.exp_tsml == q2.exp_tsml, "ArithmeticParams->Blue builder mismatch: exp_tsml");
    assert_case(q.exp_tbig == q2.exp_tbig, "ArithmeticParams->Blue builder mismatch: exp_tbig");
    assert_case(q.exp_ssml == q2.exp_ssml, "ArithmeticParams->Blue builder mismatch: exp_ssml");
    assert_case(q.exp_sbig == q2.exp_sbig, "ArithmeticParams->Blue builder mismatch: exp_sbig");

    assert_equal("ArithmeticParams->Blue tsml", q.tsml, q2.tsml);
    assert_equal("ArithmeticParams->Blue tbig", q.tbig, q2.tbig);
    assert_equal("ArithmeticParams->Blue ssml", q.ssml, q2.ssml);
    assert_equal("ArithmeticParams->Blue sbig", q.sbig, q2.sbig);

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        dual_printf("[params/%s] t=%lld emin=%lld emax=%lld\n",
                    tag, (long long)p.t, (long long)p.emin, (long long)p.emax);

        sprintnum(_spbuf, p.safmin);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), p.safmin);
        dual_printf("[params/%s] safmin %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, p.safmax);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), p.safmax);
        dual_printf("[params/%s] safmax %40s    %40s\n", tag, _spbuf, _sphexbuf);

        dual_printf("[params/%s] builder exp_tsml=%lld exp_tbig=%lld exp_ssml=%lld exp_sbig=%lld\n",
                    tag,
                    (long long)q2.exp_tsml, (long long)q2.exp_tbig,
                    (long long)q2.exp_ssml, (long long)q2.exp_sbig);

        sprintnum(_spbuf, q2.tsml);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), q2.tsml);
        dual_printf("[params/%s] builder tsml   %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.tbig);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), q2.tbig);
        dual_printf("[params/%s] builder tbig   %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.ssml);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), q2.ssml);
        dual_printf("[params/%s] builder ssml   %40s    %40s\n", tag, _spbuf, _sphexbuf);

        sprintnum(_spbuf, q2.sbig);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), q2.sbig);
        dual_printf("[params/%s] builder sbig   %40s    %40s\n", tag, _spbuf, _sphexbuf);
    }
}

void Rlamch_binary80_test() {
#if defined VERBOSE_TEST
    const bool print_values = true;
#else
    const bool print_values = false;
#endif
    const char *tag = "binary80";
    check_arithmetic_params_binary80(tag, print_values);

    auto fail = [&](const char *what) {
        printf("*** Testing Mutils (mplapack_binary80_t) failed: %s ***\n", what);
        exit(1);
    };

    auto assert_case = [&](bool cond, const char *what) {
        if (cond)
            return;
        char buf[256];
        snprintf(buf, sizeof(buf), "%s: %s", tag, what);
        fail(buf);
    };

    auto assert_equal = [&](const char *name, mplapack_binary80_t got, mplapack_binary80_t expected) {
        if (got == expected)
            return;

        printf("*** Testing Mutils (mplapack_binary80_t) failed: %s mismatch in %s ***\n", tag, name);
        printf("    got      = ");
        printnum(got);
        printf("\n");
        printf("    expected = ");
        printnum(expected);
        printf("\n");
        exit(1);
    };

    const mplapack_binary80_t zero = (mplapack_binary80_t)0.0;
    const mplapack_binary80_t one = (mplapack_binary80_t)1.0;
    const mplapack_binary80_t two = (mplapack_binary80_t)2.0;

    // Expected values for mplapack_binary80_t based on compiler-provided parameters.
#if MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_FLOAT64X
    // Using _Float64x
#if defined(__FLT64X_MANT_DIG__)
    const int p = (int)__FLT64X_MANT_DIG__;
    const int emin = (int)__FLT64X_MIN_EXP__;
    const int emax = (int)__FLT64X_MAX_EXP__;
#else
#error "_Float64x mode selected but __FLT64X_MANT_DIG__ not available"
#endif
#elif MPLAPACK_BINARY80_MODE == MPLAPACK_BINARY80_MODE_LDBL80
    // Using long double (must be binary80)
#if defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 64
    const int p = (int)__LDBL_MANT_DIG__;
    const int emin = (int)__LDBL_MIN_EXP__;
    const int emax = (int)__LDBL_MAX_EXP__;
#else
#error "long double mode selected but long double is not binary80 (64-bit mantissa required)"
#endif
#else
#error "Invalid MPLAPACK_BINARY80_MODE specified"
#endif
    const mplapack_binary80_t exB = two;
    const mplapack_binary80_t exN = (mplapack_binary80_t)p;

    // ulp(1) = 2^(1-p), unit roundoff E = ulp(1)/2.
    mplapack_binary80_t exP = one;
    for (int i = 0; i < p - 1; ++i) {
        exP /= two;
    }
    const mplapack_binary80_t exE = exP / two;

    // rmin = 2^(emin-1)
    mplapack_binary80_t exU = one;
    // emin is negative for IEEE binary formats.
    for (int i = 0; i < (-emin + 1); ++i) {
        exU /= two;
    }

    // rmax = (1 - 2^(-p)) * 2^emax  (computed without overflowing intermediates)
    mplapack_binary80_t two_to_minus_p = one;
    for (int i = 0; i < p; ++i) {
        two_to_minus_p /= two;
    }
    mplapack_binary80_t exO = one - two_to_minus_p;
    for (int i = 0; i < emax; ++i) {
        exO *= two;
    }

    const mplapack_binary80_t exM = (mplapack_binary80_t)emin;
    const mplapack_binary80_t exL = (mplapack_binary80_t)emax;

    // Safe minimum: max(rmin, (1/rmax)*(1+E)).
    const mplapack_binary80_t small = one / exO;
    const mplapack_binary80_t candidate = small * (one + exE);
    const mplapack_binary80_t exS = (candidate >= exU) ? candidate : exU;

    const mplapack_binary80_t exR = one;
    const mplapack_binary80_t exZ = zero;

    // Fetch actual values from Rlamch_mplapack_binary80_t.
    const mplapack_binary80_t gotE = Rlamch_binary80("E");
    const mplapack_binary80_t gotS = Rlamch_binary80("S");
    const mplapack_binary80_t gotB = Rlamch_binary80("B");
    const mplapack_binary80_t gotP = Rlamch_binary80("P");
    const mplapack_binary80_t gotN = Rlamch_binary80("N");
    const mplapack_binary80_t gotR = Rlamch_binary80("R");
    const mplapack_binary80_t gotM = Rlamch_binary80("M");
    const mplapack_binary80_t gotU = Rlamch_binary80("U");
    const mplapack_binary80_t gotL = Rlamch_binary80("L");
    const mplapack_binary80_t gotO = Rlamch_binary80("O");
    const mplapack_binary80_t gotZ = Rlamch_binary80("Z");

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

    // Operational property checks (mplapack_binary80_t)
    const mplapack_binary80_t next = nextafter(one, two);

    // P == ulp(1) == nextafter(1,2) - 1
    assert_case((next - one) == gotP, "P is not equal to ulp(1) = nextafter(1,2)-1");

    // Consistency: P == 2*E
    assert_case((gotE * two) == gotP, "expected P == 2*E");

    // Strong, rounding-mode-independent check: fl(1 + P) == nextafter(1,2)
    volatile mplapack_binary80_t b = one + gotP;
    assert_case(b == next, "expected fl(1 + P) == nextafter(1,2)");
    assert_case(b > one, "expected fl(1 + 2E) > 1");

    // Tie check: only assert under round-to-nearest (ties-to-even)
    if (std::fegetround() == FE_TONEAREST) {
        volatile mplapack_binary80_t a = one + gotE; // 1 + ulp/2 (midpoint)
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
        mplapack_binary80_t pow2 = one;
        if (exp_sum > 0) {
            for (int i = 0; i < exp_sum; ++i)
                pow2 *= two;
        } else if (exp_sum < 0) {
            for (int i = 0; i < -exp_sum; ++i)
                pow2 /= two;
        }
        const mplapack_binary80_t got_prod = gotU * gotO;
        const mplapack_binary80_t expected_prod = (one - exE) * pow2;
        assert_case(got_prod == expected_prod, "U*O cross-check failed (inconsistent rmin/rmax model)");
    }

    if (print_values) {
        char _spbuf[__MPLAPACK_BUFLEN__];
        char _sphexbuf[__MPLAPACK_BUFLEN__];

        sprintnum(_spbuf, gotE);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), gotE);
        dual_printf("Rlamch E: Epsilon                      %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotS);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), gotS);
        dual_printf("Rlamch S: Safe minimum                 %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotB);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), gotB);
        dual_printf("Rlamch B: Base                         %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotP);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), gotP);
        dual_printf("Rlamch P: Precision                    %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotN);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), gotN);
        dual_printf("Rlamch N: Number of digits in mantissa %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotR);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), gotR);
        dual_printf("Rlamch R: Rounding mode                %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotM);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), gotM);
        dual_printf("Rlamch M: Minimum exponent             %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotU);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), gotU);
        dual_printf("Rlamch U: Underflow threshold          %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotL);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), gotL);
        dual_printf("Rlamch L: Largest exponent             %40s    %40s\n", _spbuf, _sphexbuf);

        sprintnum(_spbuf, gotO);
        sprinthex_binary80_fixed(_sphexbuf, sizeof(_sphexbuf), gotO);
        dual_printf("Rlamch O: Overflow threshold           %40s    %40s\n", _spbuf, _sphexbuf);
    }
    check_blue_scaling_binary80(tag, emin, emax, p, print_values);
}
#endif

int main(int argc, char *argv[]) {
    g_lamch_file = fopen("Rlamch.txt", "wb"); // binary mode: suppress \r\n conversion on MinGW/Windows
    if (!g_lamch_file) {
        fprintf(stderr, "Warning: could not open Rlamch.txt for writing\n");
    }

    dual_printf("*** Testing Rlamch start ***\n");
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
#if defined ___MPLAPACK_BUILD_WITH_BINARY80___
    Rlamch_binary80_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH_BINARY128___
    Rlamch_binary128_test();
#endif
    dual_printf("*** Testing Rlamch successful ***\n");

    if (g_lamch_file) {
        fclose(g_lamch_file);
        g_lamch_file = nullptr;
    }
    return (0);
}
