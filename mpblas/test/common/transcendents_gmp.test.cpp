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
using mplapack_gmp_transcendents::div;
using mplapack_gmp_transcendents::make_ui;
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

int main() {
    std::printf("*** Testing GMP transcendents step0/4 start ***\n");
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
    std::printf("*** Testing GMP transcendents step0/4 successful ***\n");
    return 0;
}
