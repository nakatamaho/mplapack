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

int main() {
    std::printf("*** Testing GMP transcendents step0/1 start ***\n");
    test_pi_reference_literal();
    test_pi_precision_doubling_pair(128, 256);
    test_pi_precision_doubling_pair(512, 1024);
    test_pi_precision_doubling_pair(1024, 4096);
    test_pi_precision_doubling_pair(2048, 4096);
    std::printf("*** Testing GMP transcendents step0/1 successful ***\n");
    return 0;
}
