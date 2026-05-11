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

#if defined(___MPLAPACK_BUILD_WITH_MPFR___)
#error "primary GMP test must not be compiled with MPFR macro set"
#endif

#if !defined(___MPLAPACK_BUILD_WITH_GMP___)
#error "this test requires the GMP backend"
#endif

#include <cstdio>
#include <cstdlib>

#include <mpblas_gmp.h>
#include <mplapack_gmp_transcendents.h>

using mplapack_gmp_transcendents::pi;
using mplapack_gmp_transcendents::log_two;
using mplapack_gmp_transcendents::compute_log1p;
using mplapack_gmp_transcendents::compute_log;
using mplapack_gmp_transcendents::compute_exp;
using mplapack_gmp_transcendents::compute_expm1;
using mplapack_gmp_transcendents::compute_sin;
using mplapack_gmp_transcendents::compute_cos;
using mplapack_gmp_transcendents::compute_sincos;
using mplapack_gmp_transcendents::compute_atan;
using mplapack_gmp_transcendents::compute_atan2;
using mplapack_gmp_transcendents::compute_pow;
using mplapack_gmp_transcendents::div;
using mplapack_gmp_transcendents::div_ui;
using mplapack_gmp_transcendents::make_ui;
using mplapack_gmp_transcendents::quo_rem;
using mplapack_gmp_transcendents::set_prec_copy;
using mplapack_gmp_transcendents::sqr;

static mpf_class parse_decimal_literal(const char *text, mp_bitcnt_t precision) {
    mpf_class value(0, precision);
    if (mpf_set_str(value.get_mpf_t(), text, 10) != 0) {
        std::fprintf(stderr, "failed to parse decimal literal: %s\n", text);
        std::exit(1);
    }
    return value;
}

static mpf_class pi_ulp(mp_bitcnt_t precision) {
    mpf_class ulp(1, precision + 8);
    if (precision > 2) {
        mpf_div_2exp(ulp.get_mpf_t(), ulp.get_mpf_t(), precision - 2);
    } else if (precision < 2) {
        mpf_mul_2exp(ulp.get_mpf_t(), ulp.get_mpf_t(), 2 - precision);
    }
    return ulp;
}

static mpf_class ulp_for_value(const mpf_class &value, mp_bitcnt_t precision) {
    mp_exp_t exponent = 0;
    mpf_get_d_2exp(&exponent, value.get_mpf_t());

    mpf_class ulp(1, precision + 8);
    if (exponent >= static_cast<mp_exp_t>(precision)) {
        mpf_mul_2exp(ulp.get_mpf_t(), ulp.get_mpf_t(), exponent - static_cast<mp_exp_t>(precision));
    } else if (precision < 2) {
        mpf_mul_2exp(ulp.get_mpf_t(), ulp.get_mpf_t(), 2 - precision);
    } else {
        mpf_div_2exp(ulp.get_mpf_t(), ulp.get_mpf_t(), static_cast<mp_bitcnt_t>(static_cast<mp_exp_t>(precision) - exponent));
    }
    return ulp;
}

static void test_pi_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5BPi%2C+200%5D
    // Query: N[Pi, 200]
    // Retrieval date: 2026-04-22
    const char *pi_literal =
        "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899"
        "8628034825342117067982148086513282306647093844609550582231725359408128481117450284";

    const mpf_class reference_hi = parse_decimal_literal(pi_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class value = pi(target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = pi_ulp(target_precision);

    if (diff > ulp) {
        std::printf("pi literal test failed\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("pi 256-bit literal check passed\n");
}

static void test_pi_precision_doubling_pair(mp_bitcnt_t low_precision, mp_bitcnt_t high_precision) {
    const mpf_class low_value = pi(low_precision);
    const mpf_class high_value = pi(high_precision);
    mpf_class rounded_high(0, low_precision);
    mpf_set(rounded_high.get_mpf_t(), high_value.get_mpf_t());

    const mpf_class diff = abs(low_value - rounded_high);
    const mpf_class ulp = pi_ulp(low_precision);

    if (diff > ulp) {
        std::printf("pi precision-doubling test failed for %lu -> %lu bits\n",
            static_cast<unsigned long>(low_precision),
            static_cast<unsigned long>(high_precision));
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("pi %lu-bit self-check passed (vs rounded %lu-bit value)\n",
        static_cast<unsigned long>(low_precision),
        static_cast<unsigned long>(high_precision));
}

static void test_log_two_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5BLog%5B2%5D%2C+200%5D
    // Query: N[Log[2], 200]
    // Retrieval date: 2026-04-22
    const char *log_two_literal =
        "0.69314718055994530941723212145817656807550013436025525412068000949339362196969471"
        "56058633269964186875420014810205706857336855202357581305570326707516350759619307";

    const mpf_class reference_hi = parse_decimal_literal(log_two_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class value = log_two(target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, target_precision);

    if (diff > ulp) {
        std::printf("log(2) literal test failed\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("log(2) 256-bit literal check passed\n");
}

static void test_log_two_precision_doubling_pair(mp_bitcnt_t low_precision, mp_bitcnt_t high_precision) {
    const mpf_class low_value = log_two(low_precision);
    const mpf_class high_value = log_two(high_precision);
    mpf_class rounded_high(0, low_precision);
    mpf_set(rounded_high.get_mpf_t(), high_value.get_mpf_t());

    const mpf_class diff = abs(low_value - rounded_high);
    const mpf_class ulp = ulp_for_value(rounded_high, low_precision);

    if (diff > ulp) {
        std::printf("log(2) precision-doubling test failed for %lu -> %lu bits\n",
            static_cast<unsigned long>(low_precision),
            static_cast<unsigned long>(high_precision));
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("log(2) %lu-bit self-check passed (vs rounded %lu-bit value)\n",
        static_cast<unsigned long>(low_precision),
        static_cast<unsigned long>(high_precision));
}

static void test_log1p_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5BLog%5B1%2B1%2F10%5D%2C+200%5D
    // Query: N[Log[1 + 1/10], 200]
    // Retrieval date: 2026-04-22
    // Cross-checked locally with Python Decimal.ln() at 260 digits.
    const char *log1p_literal =
        "0.09531017980432486004395212328076509222060536530864419918523980816300101423588423"
        "283905750291303649307274794184585174988884604369351298063868901502170232637556873";

    const mpf_class reference_hi = parse_decimal_literal(log1p_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class x = div(make_ui(1, target_precision), make_ui(10, target_precision), target_precision);
    const mpf_class value = compute_log1p(x, target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, target_precision);

    if (diff > ulp) {
        std::printf("log1p literal test failed\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("log1p 256-bit literal check passed\n");
}

static void test_log1p_precision_doubling_pair(mp_bitcnt_t low_precision, mp_bitcnt_t high_precision) {
    const mpf_class x = div(make_ui(1, high_precision), make_ui(10, high_precision), high_precision);
    const mpf_class low_x = set_prec_copy(x, low_precision);
    const mpf_class low_value = compute_log1p(low_x, low_precision);
    const mpf_class high_value = compute_log1p(x, high_precision);
    mpf_class rounded_high(0, low_precision);
    mpf_set(rounded_high.get_mpf_t(), high_value.get_mpf_t());

    const mpf_class diff = abs(low_value - rounded_high);
    const mpf_class ulp = ulp_for_value(rounded_high, low_precision);

    if (diff > ulp) {
        std::printf("log1p precision-doubling test failed for %lu -> %lu bits\n",
            static_cast<unsigned long>(low_precision),
            static_cast<unsigned long>(high_precision));
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("log1p %lu-bit self-check passed (vs rounded %lu-bit value)\n",
        static_cast<unsigned long>(low_precision),
        static_cast<unsigned long>(high_precision));
}

static void test_log1p_near_zero() {
    const mp_bitcnt_t precision = 256;
    mpf_class x = make_ui(1, precision);
    mpf_div_2exp(x.get_mpf_t(), x.get_mpf_t(), 120);
    const mpf_class value = compute_log1p(x, precision);
    const mpf_class diff = abs(value - x);

    mpf_class bound = sqr(x, precision);
    if (diff > bound) {
        std::printf("log1p near-zero seam test failed\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nbound = ");
        printnum(bound);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("log1p near-zero seam check passed\n");
}

static void test_log_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5BLog%5B25%5D%2C+200%5D
    // Query: N[Log[25], 200]
    // Retrieval date: 2026-04-22
    // Cross-checked locally with Python Decimal.ln() at 260 digits.
    const char *log_literal =
        "3.21887582486820074920151866645237527905120270853703544382529578294835797541531552"
        "926026775618635922159993260604343112579944801045864935239926723323492741145510435";

    const mpf_class reference_hi = parse_decimal_literal(log_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class x = make_ui(25, target_precision);
    const mpf_class value = compute_log(x, target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, target_precision);

    if (diff > ulp) {
        std::printf("log literal test failed\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("log 256-bit literal check passed\n");
}

static void test_log_precision_doubling_pair(mp_bitcnt_t low_precision, mp_bitcnt_t high_precision) {
    const mpf_class x = make_ui(25, high_precision);
    const mpf_class low_x = set_prec_copy(x, low_precision);
    const mpf_class low_value = compute_log(low_x, low_precision);
    const mpf_class high_value = compute_log(x, high_precision);
    mpf_class rounded_high(0, low_precision);
    mpf_set(rounded_high.get_mpf_t(), high_value.get_mpf_t());

    const mpf_class diff = abs(low_value - rounded_high);
    const mpf_class ulp = ulp_for_value(rounded_high, low_precision);

    if (diff > ulp) {
        std::printf("log precision-doubling test failed for %lu -> %lu bits\n",
            static_cast<unsigned long>(low_precision),
            static_cast<unsigned long>(high_precision));
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("log %lu-bit self-check passed (vs rounded %lu-bit value)\n",
        static_cast<unsigned long>(low_precision),
        static_cast<unsigned long>(high_precision));
}

static void test_log_near_one_branch() {
    const mp_bitcnt_t precision = 256;
    mpf_class x = make_ui(1, precision);
    mpf_class delta = make_ui(1, precision);
    mpf_div_2exp(delta.get_mpf_t(), delta.get_mpf_t(), 80);
    x = mplapack_gmp_transcendents::add(x, delta, precision);

    const mpf_class value = compute_log(x, precision);
    const mpf_class reference = compute_log1p(delta, precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, precision);

    if (diff > ulp) {
        std::printf("log near-one branch test failed\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("log near-one branch check passed\n");
}

static void test_log_extreme_exponents() {
    const mp_bitcnt_t precision = 256;
    const mp_bitcnt_t work = 512;
    const mp_bitcnt_t shift = 600;

    mpf_class x = make_ui(1, work);
    mpf_mul_2exp(x.get_mpf_t(), x.get_mpf_t(), shift);

    mpf_class expected = log_two(work);
    mpf_mul_ui(expected.get_mpf_t(), expected.get_mpf_t(), static_cast<unsigned long>(shift));
    mpf_class rounded_expected(0, precision);
    mpf_set(rounded_expected.get_mpf_t(), expected.get_mpf_t());

    mpf_class value = compute_log(x, precision);
    mpf_class bound = ulp_for_value(rounded_expected, precision);
    mpf_mul_ui(bound.get_mpf_t(), bound.get_mpf_t(), 4ul);
    if (abs(value - rounded_expected) > bound) {
        std::printf("log extreme positive exponent test failed\n");
        std::exit(1);
    }

    x = make_ui(1, work);
    mpf_div_2exp(x.get_mpf_t(), x.get_mpf_t(), shift);
    expected = log_two(work);
    mpf_mul_ui(expected.get_mpf_t(), expected.get_mpf_t(), static_cast<unsigned long>(shift));
    expected = mplapack_gmp_transcendents::sub(make_ui(0, work), expected, work);
    mpf_set(rounded_expected.get_mpf_t(), expected.get_mpf_t());

    value = compute_log(x, precision);
    bound = ulp_for_value(rounded_expected, precision);
    mpf_mul_ui(bound.get_mpf_t(), bound.get_mpf_t(), 4ul);
    if (abs(value - rounded_expected) > bound) {
        std::printf("log extreme negative exponent test failed\n");
        std::exit(1);
    }

    std::printf("log extreme exponent checks passed\n");
}

static void test_exp_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5BExp%5B1%5D%2C+200%5D
    // Query: N[Exp[1], 200]
    // Retrieval date: 2026-04-22
    const char *exp_literal =
        "2.71828182845904523536028747135266249775724709369995957496696762772407663035354759"
        "457138217852516642742746639193200305992181741359662904357290033429526059563073814";

    const mpf_class reference_hi = parse_decimal_literal(exp_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class x = make_ui(1, target_precision);
    const mpf_class value = compute_exp(x, target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, target_precision);

    if (diff > ulp) {
        std::printf("exp literal test failed\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("exp 256-bit literal check passed\n");
}

static void test_exp_precision_doubling_pair(mp_bitcnt_t low_precision, mp_bitcnt_t high_precision) {
    const mpf_class x = make_ui(1, high_precision);
    const mpf_class low_x = set_prec_copy(x, low_precision);
    const mpf_class low_value = compute_exp(low_x, low_precision);
    const mpf_class high_value = compute_exp(x, high_precision);
    mpf_class rounded_high(0, low_precision);
    mpf_set(rounded_high.get_mpf_t(), high_value.get_mpf_t());

    const mpf_class diff = abs(low_value - rounded_high);
    const mpf_class ulp = ulp_for_value(rounded_high, low_precision);

    if (diff > ulp) {
        std::printf("exp precision-doubling test failed for %lu -> %lu bits\n",
            static_cast<unsigned long>(low_precision),
            static_cast<unsigned long>(high_precision));
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("exp %lu-bit self-check passed (vs rounded %lu-bit value)\n",
        static_cast<unsigned long>(low_precision),
        static_cast<unsigned long>(high_precision));
}

static void test_exp_near_zero() {
    const mp_bitcnt_t precision = 256;
    mpf_class x = make_ui(1, precision);
    mpf_div_2exp(x.get_mpf_t(), x.get_mpf_t(), 120);
    const mpf_class value = compute_exp(x, precision);
    const mpf_class one = make_ui(1, precision);
    const mpf_class linear = mplapack_gmp_transcendents::add(one, x, precision);
    const mpf_class diff = abs(value - linear);

    mpf_class bound = sqr(x, precision);
    if (diff > bound) {
        std::printf("exp near-zero seam test failed\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nbound = ");
        printnum(bound);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("exp near-zero seam check passed\n");
}

static void test_expm1_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5BExp%5B1%2F10%5D-1%2C+200%5D
    // Query: N[Exp[1/10] - 1, 200]
    // Retrieval date: 2026-04-22
    const char *expm1_literal =
        "0.10517091807564762481170782649024666822454719473751871879286328944096796674765462"
        "180826680334383576123364162622389881639224377083652885920639130690370248999245655";

    const mpf_class reference_hi = parse_decimal_literal(expm1_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class x = div(make_ui(1, target_precision), make_ui(10, target_precision), target_precision);
    const mpf_class value = compute_expm1(x, target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, target_precision);

    if (diff > ulp) {
        std::printf("expm1 literal test failed\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("expm1 256-bit literal check passed\n");
}

static void test_expm1_precision_doubling_pair(mp_bitcnt_t low_precision, mp_bitcnt_t high_precision) {
    const mpf_class x = div(make_ui(1, high_precision), make_ui(10, high_precision), high_precision);
    const mpf_class low_x = set_prec_copy(x, low_precision);
    const mpf_class low_value = compute_expm1(low_x, low_precision);
    const mpf_class high_value = compute_expm1(x, high_precision);
    mpf_class rounded_high(0, low_precision);
    mpf_set(rounded_high.get_mpf_t(), high_value.get_mpf_t());

    const mpf_class diff = abs(low_value - rounded_high);
    const mpf_class ulp = ulp_for_value(rounded_high, low_precision);

    if (diff > ulp) {
        std::printf("expm1 precision-doubling test failed for %lu -> %lu bits\n",
            static_cast<unsigned long>(low_precision),
            static_cast<unsigned long>(high_precision));
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("expm1 %lu-bit self-check passed (vs rounded %lu-bit value)\n",
        static_cast<unsigned long>(low_precision),
        static_cast<unsigned long>(high_precision));
}

static void test_expm1_near_zero() {
    const mp_bitcnt_t precision = 256;
    mpf_class x = make_ui(1, precision);
    mpf_div_2exp(x.get_mpf_t(), x.get_mpf_t(), 120);
    const mpf_class value = compute_expm1(x, precision);
    const mpf_class diff = abs(value - x);

    mpf_class bound = sqr(x, precision);
    if (diff > bound) {
        std::printf("expm1 near-zero seam test failed\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nbound = ");
        printnum(bound);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("expm1 near-zero seam check passed\n");
}

static void test_quo_rem_rounding() {
    const mp_bitcnt_t precision = 256;
    const mpf_class one = make_ui(1, precision);
    const mpf_class ten = make_ui(10, precision);

    {
        const mpf_class x = div(make_ui(17, precision), ten, precision);
        const auto qr = quo_rem(x, one, precision);
        if (qr.quotient != 2 || abs(qr.remainder - parse_decimal_literal("-0.3", precision)) > ulp_for_value(one, precision)) {
            std::printf("quo_rem positive nearest test failed\n");
            std::exit(1);
        }
    }
    {
        const mpf_class x = parse_decimal_literal("-1.7", precision);
        const auto qr = quo_rem(x, one, precision);
        if (qr.quotient != -2 || abs(qr.remainder - parse_decimal_literal("0.3", precision)) > ulp_for_value(one, precision)) {
            std::printf("quo_rem negative nearest test failed\n");
            std::exit(1);
        }
    }
    {
        const auto qr = quo_rem(div(make_ui(5, precision), make_ui(2, precision), precision), one, precision);
        if (qr.quotient != 2 || abs(qr.remainder - parse_decimal_literal("0.5", precision)) > ulp_for_value(one, precision)) {
            std::printf("quo_rem positive tie-even test failed\n");
            std::exit(1);
        }
    }
    {
        const auto qr = quo_rem(div(make_ui(3, precision), make_ui(2, precision), precision), one, precision);
        if (qr.quotient != 2 || abs(qr.remainder - parse_decimal_literal("-0.5", precision)) > ulp_for_value(one, precision)) {
            std::printf("quo_rem positive tie-up-to-even test failed\n");
            std::exit(1);
        }
    }
    {
        const auto qr = quo_rem(parse_decimal_literal("-2.5", precision), one, precision);
        if (qr.quotient != -2 || abs(qr.remainder - parse_decimal_literal("-0.5", precision)) > ulp_for_value(one, precision)) {
            std::printf("quo_rem negative tie-even test failed\n");
            std::exit(1);
        }
    }
    {
        const auto qr = quo_rem(parse_decimal_literal("-1.5", precision), one, precision);
        if (qr.quotient != -2 || abs(qr.remainder - parse_decimal_literal("0.5", precision)) > ulp_for_value(one, precision)) {
            std::printf("quo_rem negative tie-down-to-even test failed\n");
            std::exit(1);
        }
    }
    {
        const mpf_class x = parse_decimal_literal("7.25", precision);
        const mpf_class y = parse_decimal_literal("0.7", precision);
        const auto qr = quo_rem(x, y, precision);
        const mpf_class half_y = div(abs(y), make_ui(2, precision), precision);
        if (abs(qr.remainder) > half_y) {
            std::printf("quo_rem remainder bound test failed\n");
            std::exit(1);
        }
    }
    {
        const mp_bitcnt_t target_precision = 160;
        const mp_bitcnt_t source_precision = 1024;
        const mp_bitcnt_t shift = 600;

        mpf_class x = make_ui(1, source_precision);
        mpf_mul_2exp(x.get_mpf_t(), x.get_mpf_t(), shift);
        const mpf_class expected_remainder = div(make_ui(1, source_precision), make_ui(3, source_precision), source_precision);
        x = mplapack_gmp_transcendents::add(x, expected_remainder, source_precision);

        const auto qr = quo_rem(x, make_ui(1, source_precision), target_precision);
        mpz_class expected_quotient = 1;
        mpz_mul_2exp(expected_quotient.get_mpz_t(), expected_quotient.get_mpz_t(), shift);

        const mpf_class rounded_expected = set_prec_copy(expected_remainder, target_precision);
        const mpf_class diff = abs(qr.remainder - rounded_expected);
        mpf_class bound = ulp_for_value(rounded_expected, target_precision);
        mpf_mul_2exp(bound.get_mpf_t(), bound.get_mpf_t(), 4);
        if (qr.quotient != expected_quotient || diff > bound) {
            std::printf("quo_rem large cancellation test failed\n");
            std::printf("diff = ");
            printnum(diff);
            std::printf("\nbound = ");
            printnum(bound);
            std::printf("\n");
            std::exit(1);
        }
    }
    std::printf("quo_rem round-to-nearest-even checks passed\n");
}

static void test_sin_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5BSin%5B1%5D%2C+200%5D
    // Query: N[Sin[1], 200]
    // Retrieval date: 2026-04-22
    // Cross-checked locally with Python Decimal Taylor summation at 260 digits.
    const char *sin_literal =
        "0.84147098480789650665250232163029899962256306079837106567275170999191040439123966"
        "894863974354305269585434903790792067429325911892099189888119341032772921240948079";

    const mpf_class reference_hi = parse_decimal_literal(sin_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class x = make_ui(1, target_precision);
    const mpf_class value = compute_sin(x, target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, target_precision);

    if (diff > ulp) {
        std::printf("sin literal test failed\n");
        std::printf("value = ");
        printnum(value);
        std::printf("\nreference = ");
        printnum(reference);
        std::printf("\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("sin 256-bit literal check passed\n");
}

static void test_cos_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5BCos%5B1%5D%2C+200%5D
    // Query: N[Cos[1], 200]
    // Retrieval date: 2026-04-22
    // Cross-checked locally with Python Decimal Taylor summation at 260 digits.
    const char *cos_literal =
        "0.54030230586813971740093660744297660373231042061792222767009725538110039477447176"
        "451795185608718308934357173116003008909786063376002166345640651226541731858471797";

    const mpf_class reference_hi = parse_decimal_literal(cos_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class x = make_ui(1, target_precision);
    const mpf_class value = compute_cos(x, target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, target_precision);

    if (diff > ulp) {
        std::printf("cos literal test failed\n");
        std::printf("value = ");
        printnum(value);
        std::printf("\nreference = ");
        printnum(reference);
        std::printf("\n");
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("cos 256-bit literal check passed\n");
}

static void test_sin_precision_doubling_pair(mp_bitcnt_t low_precision, mp_bitcnt_t high_precision) {
    const mpf_class x = make_ui(1, high_precision);
    const mpf_class low_x = set_prec_copy(x, low_precision);
    const mpf_class low_value = compute_sin(low_x, low_precision);
    const mpf_class high_value = compute_sin(x, high_precision);
    mpf_class rounded_high(0, low_precision);
    mpf_set(rounded_high.get_mpf_t(), high_value.get_mpf_t());

    const mpf_class diff = abs(low_value - rounded_high);
    const mpf_class ulp = ulp_for_value(rounded_high, low_precision);
    if (diff > ulp) {
        std::printf("sin precision-doubling test failed for %lu -> %lu bits\n",
            static_cast<unsigned long>(low_precision),
            static_cast<unsigned long>(high_precision));
        std::exit(1);
    }
    std::printf("sin %lu-bit self-check passed (vs rounded %lu-bit value)\n",
        static_cast<unsigned long>(low_precision),
        static_cast<unsigned long>(high_precision));
}

static void test_cos_precision_doubling_pair(mp_bitcnt_t low_precision, mp_bitcnt_t high_precision) {
    const mpf_class x = make_ui(1, high_precision);
    const mpf_class low_x = set_prec_copy(x, low_precision);
    const mpf_class low_value = compute_cos(low_x, low_precision);
    const mpf_class high_value = compute_cos(x, high_precision);
    mpf_class rounded_high(0, low_precision);
    mpf_set(rounded_high.get_mpf_t(), high_value.get_mpf_t());

    const mpf_class diff = abs(low_value - rounded_high);
    const mpf_class ulp = ulp_for_value(rounded_high, low_precision);
    if (diff > ulp) {
        std::printf("cos precision-doubling test failed for %lu -> %lu bits\n",
            static_cast<unsigned long>(low_precision),
            static_cast<unsigned long>(high_precision));
        std::exit(1);
    }
    std::printf("cos %lu-bit self-check passed (vs rounded %lu-bit value)\n",
        static_cast<unsigned long>(low_precision),
        static_cast<unsigned long>(high_precision));
}

static void test_sincos_identity() {
    const mp_bitcnt_t precision = 256;
    const mpf_class x = div(make_ui(7, precision), make_ui(10, precision), precision);
    const auto values = compute_sincos(x, precision);
    const mpf_class sum = mplapack_gmp_transcendents::add(sqr(values.sin_value, precision), sqr(values.cos_value, precision), precision);
    const mpf_class diff = abs(sum - make_ui(1, precision));
    const mpf_class ulp = ulp_for_value(make_ui(1, precision), precision);
    if (diff > ulp) {
        std::printf("sincos identity test failed\n");
        std::exit(1);
    }
    std::printf("sincos identity check passed\n");
}

static void test_sincos_reduction_seam() {
    const mp_bitcnt_t precision = 256;
    mpf_class x = pi(precision);
    mpf_div_2exp(x.get_mpf_t(), x.get_mpf_t(), 1);
    mpf_class delta = make_ui(1, precision);
    mpf_div_2exp(delta.get_mpf_t(), delta.get_mpf_t(), 80);
    x = mplapack_gmp_transcendents::sub(x, delta, precision);

    const auto values = compute_sincos(x, precision);
    const mpf_class sin_diff = abs(values.sin_value - make_ui(1, precision));
    if (sin_diff > delta) {
        std::printf("sincos reduction seam test failed\n");
        std::exit(1);
    }
    std::printf("sincos reduction seam check passed\n");
}

static void test_sincos_large_argument_reduction() {
    const mp_bitcnt_t precision = 256;
    const mp_bitcnt_t source_precision = 1536;
    const mp_bitcnt_t shift = 600;

    mpf_class x = pi(source_precision);
    mpf_mul_2exp(x.get_mpf_t(), x.get_mpf_t(), shift);

    const auto values = compute_sincos(x, precision);
    mpf_class bound = make_ui(1, precision);
    mpf_div_2exp(bound.get_mpf_t(), bound.get_mpf_t(), 180);

    if (abs(values.sin_value) > bound || abs(values.cos_value - make_ui(1, precision)) > bound) {
        std::printf("sincos large-argument reduction test failed\n");
        std::printf("sin = ");
        printnum(values.sin_value);
        std::printf("\ncos = ");
        printnum(values.cos_value);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("sincos large-argument reduction check passed\n");
}

static void check_sincos_taylor_denominator_step(const char *name, unsigned long den1, unsigned long den2) {
    const mp_bitcnt_t precision = 512;
    const mp_bitcnt_t reference_precision = 1024;

    const mpf_class numerator = div(make_ui(7, precision), make_ui(13, precision), precision);
    mpf_class split = div_ui(numerator, den1, precision);
    split = div_ui(split, den2, precision);

    const mpz_class exact_denominator = mpz_class(den1) * mpz_class(den2);
    mpf_class denominator(0, reference_precision);
    mpf_set_z(denominator.get_mpf_t(), exact_denominator.get_mpz_t());
    const mpf_class numerator_hi = set_prec_copy(numerator, reference_precision);
    const mpf_class reference_hi = div(numerator_hi, denominator, reference_precision);
    const mpf_class reference = set_prec_copy(reference_hi, precision);

    const mpf_class diff = abs(split - reference);
    mpf_class bound = ulp_for_value(reference, precision);
    mpf_mul_2exp(bound.get_mpf_t(), bound.get_mpf_t(), 4);
    if (diff > bound) {
        std::printf("sincos Taylor denominator split test failed for %s\n", name);
        std::printf("diff = ");
        printnum(diff);
        std::printf("\nbound = ");
        printnum(bound);
        std::printf("\n");
        std::exit(1);
    }

    if (sizeof(unsigned long) == 4) {
        const unsigned long wrapped_denominator = den1 * den2;
        const mpf_class wrapped = div(numerator, make_ui(wrapped_denominator, precision), precision);
        mpf_class wrapped_gap_bound = bound;
        mpf_mul_2exp(wrapped_gap_bound.get_mpf_t(), wrapped_gap_bound.get_mpf_t(), 12);
        if (abs(wrapped - reference) <= wrapped_gap_bound) {
            std::printf("sincos Taylor denominator overflow guard did not exercise %s\n", name);
            std::exit(1);
        }
    }
}

static void test_sincos_taylor_no_denominator_overflow() {
    const unsigned long sin_k = 32768ul;
    check_sincos_taylor_denominator_step("sin", 2ul * sin_k, 2ul * sin_k + 1ul);

    const unsigned long cos_k = 32769ul;
    check_sincos_taylor_denominator_step("cos", 2ul * cos_k - 1ul, 2ul * cos_k);

    std::printf("sincos Taylor denominator overflow check passed\n");
}

static void test_atan_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5BArcTan%5B1%5D%2C+200%5D
    // Query: N[ArcTan[1], 200]
    // Retrieval date: 2026-04-22
    const char *atan_literal =
        "0.78539816339744830961566084581987572104929234984377645524373614807695410157155224"
        "965700870633552926699553702162832057666177346115238764555793133985203212027936257";

    const mpf_class reference_hi = parse_decimal_literal(atan_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class x = make_ui(1, target_precision);
    const mpf_class value = compute_atan(x, target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, target_precision);

    if (diff > ulp) {
        std::printf("atan literal test failed\n");
        std::printf("value = ");
        printnum(value);
        std::printf("\nreference = ");
        printnum(reference);
        std::printf("\ndiff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("atan 256-bit literal check passed\n");
}

static void test_atan_precision_doubling_pair(mp_bitcnt_t low_precision, mp_bitcnt_t high_precision) {
    const mpf_class x = make_ui(1, high_precision);
    const mpf_class low_x = set_prec_copy(x, low_precision);
    const mpf_class low_value = compute_atan(low_x, low_precision);
    const mpf_class high_value = compute_atan(x, high_precision);
    mpf_class rounded_high(0, low_precision);
    mpf_set(rounded_high.get_mpf_t(), high_value.get_mpf_t());

    const mpf_class diff = abs(low_value - rounded_high);
    const mpf_class ulp = ulp_for_value(rounded_high, low_precision);
    if (diff > ulp) {
        std::printf("atan precision-doubling test failed for %lu -> %lu bits\n",
            static_cast<unsigned long>(low_precision),
            static_cast<unsigned long>(high_precision));
        std::exit(1);
    }
    std::printf("atan %lu-bit self-check passed (vs rounded %lu-bit value)\n",
        static_cast<unsigned long>(low_precision),
        static_cast<unsigned long>(high_precision));
}

static void test_atan_odd_symmetry() {
    const mp_bitcnt_t precision = 256;
    const mpf_class x = parse_decimal_literal("0.7", precision);
    const mpf_class value = compute_atan(x, precision);
    const mpf_class neg_value = compute_atan(parse_decimal_literal("-0.7", precision), precision);
    const mpf_class diff = abs(mplapack_gmp_transcendents::add(value, neg_value, precision));
    const mpf_class ulp = ulp_for_value(value, precision);
    if (diff > ulp) {
        std::printf("atan odd symmetry test failed\n");
        std::exit(1);
    }
    std::printf("atan odd symmetry check passed\n");
}

static void test_atan2_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5BArcTan%5B1%2C+1%5D%2C+200%5D
    // Query: N[ArcTan[1, 1], 200]
    // Retrieval date: 2026-04-22
    const char *atan2_literal =
        "0.78539816339744830961566084581987572104929234984377645524373614807695410157155224"
        "965700870633552926699553702162832057666177346115238764555793133985203212027936257";

    const mpf_class reference_hi = parse_decimal_literal(atan2_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class one = make_ui(1, target_precision);
    const mpf_class value = compute_atan2(one, one, target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, target_precision);
    if (diff > ulp) {
        std::printf("atan2 literal test failed\n");
        std::printf("value = ");
        printnum(value);
        std::printf("\nreference = ");
        printnum(reference);
        std::printf("\ndiff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("atan2 256-bit literal check passed\n");
}

static void test_atan2_axis_convention() {
    const mp_bitcnt_t precision = 256;
    const mpf_class zero = make_ui(0, precision);
    const mpf_class one = make_ui(1, precision);
    const mpf_class minus_one = parse_decimal_literal("-1", precision);
    mpf_class pio2 = pi(precision);
    mpf_div_2exp(pio2.get_mpf_t(), pio2.get_mpf_t(), 1);
    const mpf_class pi_value = pi(precision);
    const mpf_class ulp = ulp_for_value(pi_value, precision);

    if (abs(compute_atan2(zero, one, precision) - zero) > ulp) {
        std::printf("atan2 axis convention test failed: atan2(0,+1)\n");
        std::exit(1);
    }
    if (abs(compute_atan2(zero, minus_one, precision) - pi_value) > ulp) {
        std::printf("atan2 axis convention test failed: atan2(0,-1)\n");
        std::exit(1);
    }
    if (abs(compute_atan2(one, zero, precision) - pio2) > ulp) {
        std::printf("atan2 axis convention test failed: atan2(+1,0)\n");
        std::exit(1);
    }
    if (abs(compute_atan2(minus_one, zero, precision) + pio2) > ulp) {
        std::printf("atan2 axis convention test failed: atan2(-1,0)\n");
        std::exit(1);
    }
    if (abs(compute_atan2(zero, zero, precision) - zero) > ulp) {
        std::printf("atan2 axis convention test failed: atan2(0,0)\n");
        std::exit(1);
    }
    std::printf("atan2 axis convention check passed\n");
}

static void test_pow_reference_literal() {
    const mp_bitcnt_t target_precision = 256;
    const mp_bitcnt_t literal_precision = 768;

    // Wolfram Alpha:
    // https://www.wolframalpha.com/input?i=N%5B2%5E%281%2F2%29%2C+200%5D
    // Query: N[2^(1/2), 200]
    // Retrieval date: 2026-04-22
    const char *pow_literal =
        "1.41421356237309504880168872420969807856967187537694807317667973799073247846210703"
        "885038753432764157273501384623091229702492483605585073721264412149709993583141322";

    const mpf_class reference_hi = parse_decimal_literal(pow_literal, literal_precision);
    mpf_class reference(0, target_precision);
    mpf_set(reference.get_mpf_t(), reference_hi.get_mpf_t());

    const mpf_class x = make_ui(2, target_precision);
    const mpf_class y = div(make_ui(1, target_precision), make_ui(2, target_precision), target_precision);
    const mpf_class value = compute_pow(x, y, target_precision);
    const mpf_class diff = abs(value - reference);
    const mpf_class ulp = ulp_for_value(reference, target_precision);
    if (diff > ulp) {
        std::printf("pow literal test failed\n");
        std::printf("value = ");
        printnum(value);
        std::printf("\nreference = ");
        printnum(reference);
        std::printf("\ndiff = ");
        printnum(diff);
        std::printf("\nulp  = ");
        printnum(ulp);
        std::printf("\n");
        std::exit(1);
    }
    std::printf("pow 256-bit literal check passed\n");
}

static void test_pow_integer_exact() {
    const mp_bitcnt_t precision = 256;
    const mpf_class x = make_ui(2, precision);
    const mpf_class y = make_ui(10, precision);
    const mpf_class value = compute_pow(x, y, precision);
    const mpf_class expected = make_ui(1024, precision);
    if (value != expected) {
        std::printf("pow integer exact test failed\n");
        std::exit(1);
    }
    std::printf("pow integer exact check passed\n");
}

static void test_pow_negative_integer_exponent() {
    const mp_bitcnt_t precision = 256;
    const mpf_class x = parse_decimal_literal("-2", precision);
    const mpf_class y = make_ui(3, precision);
    const mpf_class value = compute_pow(x, y, precision);
    const mpf_class expected = parse_decimal_literal("-8", precision);
    if (value != expected) {
        std::printf("pow negative-base integer exponent test failed\n");
        std::exit(1);
    }
    std::printf("pow negative-base integer exponent check passed\n");
}

static void test_pow_precision_doubling_pair(mp_bitcnt_t low_precision, mp_bitcnt_t high_precision) {
    const mpf_class x = make_ui(2, high_precision);
    const mpf_class y = div(make_ui(1, high_precision), make_ui(2, high_precision), high_precision);
    const mpf_class low_x = set_prec_copy(x, low_precision);
    const mpf_class low_y = set_prec_copy(y, low_precision);
    const mpf_class low_value = compute_pow(low_x, low_y, low_precision);
    const mpf_class high_value = compute_pow(x, y, high_precision);
    mpf_class rounded_high(0, low_precision);
    mpf_set(rounded_high.get_mpf_t(), high_value.get_mpf_t());
    const mpf_class diff = abs(low_value - rounded_high);
    const mpf_class ulp = ulp_for_value(rounded_high, low_precision);
    if (diff > ulp) {
        std::printf("pow precision-doubling test failed for %lu -> %lu bits\n",
            static_cast<unsigned long>(low_precision),
            static_cast<unsigned long>(high_precision));
        std::exit(1);
    }
    std::printf("pow %lu-bit self-check passed (vs rounded %lu-bit value)\n",
        static_cast<unsigned long>(low_precision),
        static_cast<unsigned long>(high_precision));
}

int main() {
    std::printf("*** Testing GMP transcendents step0/10 start ***\n");
    test_pi_reference_literal();
    test_pi_precision_doubling_pair(128, 256);
    test_pi_precision_doubling_pair(512, 1024);
    test_pi_precision_doubling_pair(1024, 4096);
    test_pi_precision_doubling_pair(2048, 4096);
    test_log_two_reference_literal();
    test_log_two_precision_doubling_pair(128, 256);
    test_log_two_precision_doubling_pair(512, 1024);
    test_log_two_precision_doubling_pair(1024, 4096);
    test_log_two_precision_doubling_pair(2048, 4096);
    test_log1p_reference_literal();
    test_log1p_precision_doubling_pair(128, 256);
    test_log1p_precision_doubling_pair(512, 1024);
    test_log1p_precision_doubling_pair(1024, 4096);
    test_log1p_precision_doubling_pair(2048, 4096);
    test_log1p_near_zero();
    test_log_reference_literal();
    test_log_precision_doubling_pair(128, 256);
    test_log_precision_doubling_pair(512, 1024);
    test_log_precision_doubling_pair(1024, 4096);
    test_log_precision_doubling_pair(2048, 4096);
    test_log_near_one_branch();
    test_log_extreme_exponents();
    test_exp_reference_literal();
    test_exp_precision_doubling_pair(128, 256);
    test_exp_precision_doubling_pair(512, 1024);
    test_exp_precision_doubling_pair(1024, 4096);
    test_exp_precision_doubling_pair(2048, 4096);
    test_exp_near_zero();
    test_expm1_reference_literal();
    test_expm1_precision_doubling_pair(128, 256);
    test_expm1_precision_doubling_pair(512, 1024);
    test_expm1_precision_doubling_pair(1024, 4096);
    test_expm1_precision_doubling_pair(2048, 4096);
    test_expm1_near_zero();
    test_quo_rem_rounding();
    test_sin_reference_literal();
    test_cos_reference_literal();
    test_sin_precision_doubling_pair(128, 256);
    test_sin_precision_doubling_pair(512, 1024);
    test_sin_precision_doubling_pair(1024, 4096);
    test_cos_precision_doubling_pair(128, 256);
    test_cos_precision_doubling_pair(512, 1024);
    test_cos_precision_doubling_pair(1024, 4096);
    test_sincos_identity();
    test_sincos_reduction_seam();
    test_sincos_large_argument_reduction();
    test_sincos_taylor_no_denominator_overflow();
    test_atan_reference_literal();
    test_atan_precision_doubling_pair(128, 256);
    test_atan_precision_doubling_pair(512, 1024);
    test_atan_precision_doubling_pair(1024, 4096);
    test_atan_odd_symmetry();
    test_atan2_reference_literal();
    test_atan2_axis_convention();
    test_pow_reference_literal();
    test_pow_integer_exact();
    test_pow_negative_integer_exponent();
    test_pow_precision_doubling_pair(128, 256);
    test_pow_precision_doubling_pair(512, 1024);
    test_pow_precision_doubling_pair(1024, 4096);
    std::printf("*** Testing GMP transcendents step0/10 successful ***\n");
    return 0;
}
