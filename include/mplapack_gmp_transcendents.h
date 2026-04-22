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

#ifndef _MPLAPACK_GMP_TRANSCENDENTS_H_
#define _MPLAPACK_GMP_TRANSCENDENTS_H_

#include <algorithm>
#include <mutex>
#include <stdexcept>

#include "gmpxx.h"

namespace mplapack_gmp_transcendents {

using precision_type = mp_bitcnt_t;

inline precision_type normalize_target_precision(precision_type target_precision) {
    return std::max<precision_type>(target_precision, 32);
}

inline unsigned long ceil_log2_precision(precision_type value) {
    if (value <= 1) {
        return 0;
    }
    unsigned long result = 0;
    precision_type power = 1;
    while (power < value) {
        power <<= 1;
        result++;
    }
    return result;
}

inline precision_type guard_bits_for_pi(precision_type) {
    //
    // Faithful rounding needs roughly 1 ulp at the target precision.
    // The AGM iteration is self-correcting, but the final division and
    // repeated square roots still need a non-trivial guard margin.
    // Use a fixed +96-bit budget for step 1; benchmark/tune later.
    //
    return 96;
}

inline precision_type working_precision_for_pi(precision_type target_precision) {
    return normalize_target_precision(target_precision) + guard_bits_for_pi(target_precision);
}

inline mpf_class make_ui(unsigned long value, precision_type precision) {
    mpf_class result(0, precision);
    mpf_set_ui(result.get_mpf_t(), value);
    return result;
}

inline mpf_class set_prec_copy(const mpf_class &value, precision_type precision) {
    mpf_class result(0, precision);
    mpf_set(result.get_mpf_t(), value.get_mpf_t());
    return result;
}

inline mpf_class add(const mpf_class &a, const mpf_class &b, precision_type precision) {
    mpf_class result(0, precision);
    mpf_add(result.get_mpf_t(), a.get_mpf_t(), b.get_mpf_t());
    return result;
}

inline mpf_class sub(const mpf_class &a, const mpf_class &b, precision_type precision) {
    mpf_class result(0, precision);
    mpf_sub(result.get_mpf_t(), a.get_mpf_t(), b.get_mpf_t());
    return result;
}

inline mpf_class mul(const mpf_class &a, const mpf_class &b, precision_type precision) {
    mpf_class result(0, precision);
    mpf_mul(result.get_mpf_t(), a.get_mpf_t(), b.get_mpf_t());
    return result;
}

inline mpf_class sqr(const mpf_class &a, precision_type precision) {
    return mul(a, a, precision);
}

inline mpf_class div(const mpf_class &a, const mpf_class &b, precision_type precision) {
    mpf_class result(0, precision);
    mpf_div(result.get_mpf_t(), a.get_mpf_t(), b.get_mpf_t());
    return result;
}

inline mpf_class mul_ui(const mpf_class &a, unsigned long value, precision_type precision) {
    mpf_class result(0, precision);
    mpf_mul_ui(result.get_mpf_t(), a.get_mpf_t(), value);
    return result;
}

struct quo_rem_result {
    mpz_class quotient;
    mpf_class remainder = make_ui(0, 32);
};

struct sincos_result {
    mpf_class sin_value = make_ui(0, 32);
    mpf_class cos_value = make_ui(1, 32);
};


struct trig_constant_cache_state {
    std::mutex mutex;
    precision_type cached_precision = 0;
    mpf_class pi_value = make_ui(0, 32);
    mpf_class pi_over_two_value = make_ui(0, 32);
    mpf_class two_over_pi_value = make_ui(0, 32);
    bool initialized = false;
};

inline mpf_class sqrt_prec(const mpf_class &a, precision_type precision) {
    mpf_class result(0, precision);
    mpf_sqrt(result.get_mpf_t(), a.get_mpf_t());
    return result;
}

inline mpf_class average(const mpf_class &a, const mpf_class &b, precision_type precision) {
    mpf_class result = add(a, b, precision);
    mpf_div_2exp(result.get_mpf_t(), result.get_mpf_t(), 1);
    return result;
}

inline mpf_class inv_sqrt_ui(unsigned long value, precision_type precision) {
    mpf_class denominator(0, precision);
    mpf_set_ui(denominator.get_mpf_t(), value);
    mpf_sqrt(denominator.get_mpf_t(), denominator.get_mpf_t());
    mpf_class result(0, precision);
    mpf_ui_div(result.get_mpf_t(), 1ul, denominator.get_mpf_t());
    return result;
}

inline trig_constant_cache_state &trig_constant_cache() {
    static trig_constant_cache_state cache;
    return cache;
}

inline quo_rem_result quo_rem(const mpf_class &x_input, const mpf_class &y_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = target + 64;
    const mpf_class x = set_prec_copy(x_input, work);
    const mpf_class y = set_prec_copy(y_input, work);
    if (y == make_ui(0, work)) {
        throw std::domain_error("quo_rem division by zero");
    }

    const mpf_class q = div(x, y, work);
    mpz_class k = 0;
    mpz_set_f(k.get_mpz_t(), q.get_mpf_t());

    const mpf_class integer_part(k, work);
    const mpf_class frac = sub(q, integer_part, work);
    const mpf_class half = div(make_ui(1, work), make_ui(2, work), work);
    const mpf_class neg_half = sub(make_ui(0, work), half, work);

    const int frac_vs_half = mpf_cmp(frac.get_mpf_t(), half.get_mpf_t());
    const int frac_vs_neg_half = mpf_cmp(frac.get_mpf_t(), neg_half.get_mpf_t());
    if (frac_vs_half > 0) {
        k += 1;
    } else if (frac_vs_half == 0) {
        if (mpz_odd_p(k.get_mpz_t())) {
            k += 1;
        }
    } else if (frac_vs_neg_half < 0) {
        k -= 1;
    } else if (frac_vs_neg_half == 0) {
        if (mpz_odd_p(k.get_mpz_t())) {
            k -= 1;
        }
    }

    quo_rem_result result;
    result.quotient = k;
    result.remainder = set_prec_copy(sub(x, mul(y, mpf_class(k, work), work), work), target);
    return result;
}

inline mpf_class compute_pi_gauss_legendre(precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_pi(target);

    mpf_class zero = make_ui(0, work);
    mpf_class one = make_ui(1, work);
    mpf_class two = make_ui(2, work);
    mpf_class four = make_ui(4, work);
    mpf_class quarter = set_prec_copy(one, work);
    mpf_div_2exp(quarter.get_mpf_t(), quarter.get_mpf_t(), 2);

    mpf_class a = one;
    mpf_class b = inv_sqrt_ui(2, work);
    mpf_class t = quarter;
    mpf_class p = one;
    mpf_class epsilon = make_ui(1, work);
    mpf_div_2exp(epsilon.get_mpf_t(), epsilon.get_mpf_t(), work);
    mpf_class a_next = zero;
    mpf_class b_next = zero;
    mpf_class t_next = zero;
    mpf_class pi_previous = zero;
    mpf_class pi_current = zero;
    mpf_class tmp = zero;

    while (true) {
        a_next = average(a, b, work);
        b_next = sqrt_prec(mul(a, b, work), work);
        const mpf_class diff = sub(a, a_next, work);
        t_next = sub(t, mul(p, sqr(diff, work), work), work);

        a = a_next;
        b = b_next;
        t = t_next;
        p = mul(p, two, work);

        const mpf_class sum = add(a, b, work);
        const mpf_class numerator = sqr(sum, work);
        const mpf_class denominator = mul(t, four, work);

        pi_previous = pi_current;
        pi_current = div(numerator, denominator, work);
        tmp = abs(sub(pi_current, pi_previous, work));
        if (tmp < epsilon) {
            break;
        }
    }

    return set_prec_copy(pi_current, target);
}

struct pi_cache_state {
    std::mutex mutex;
    precision_type cached_precision = 0;
    mpf_class cached_value = make_ui(0, 32);
    bool initialized = false;
};

inline pi_cache_state &pi_cache() {
    static pi_cache_state cache;
    return cache;
}

inline mpf_class pi(precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    pi_cache_state &cache = pi_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    if (!cache.initialized || cache.cached_precision < target) {
        cache.cached_value = compute_pi_gauss_legendre(target);
        cache.cached_precision = target;
        cache.initialized = true;
    }
    return set_prec_copy(cache.cached_value, target);
}

inline precision_type guard_bits_for_log_two(precision_type) {
    return 128;
}

inline precision_type working_precision_for_log_two(precision_type target_precision) {
    return normalize_target_precision(target_precision) + guard_bits_for_log_two(target_precision);
}

inline mpf_class theta_series_threshold(precision_type precision) {
    mpf_class threshold = make_ui(1, precision);
    mpf_div_2exp(threshold.get_mpf_t(), threshold.get_mpf_t(), precision);
    return threshold;
}

inline mpf_class agm_converged(const mpf_class &a_initial, const mpf_class &b_initial, precision_type precision) {
    mpf_class a = set_prec_copy(a_initial, precision);
    mpf_class b = set_prec_copy(b_initial, precision);
    mpf_class epsilon = make_ui(1, precision);
    mpf_div_2exp(epsilon.get_mpf_t(), epsilon.get_mpf_t(), precision);

    while (true) {
        const mpf_class a_next = average(a, b, precision);
        const mpf_class b_next = sqrt_prec(mul(a, b, precision), precision);
        if (abs(sub(a_next, b_next, precision)) < epsilon) {
            return average(a_next, b_next, precision);
        }
        a = a_next;
        b = b_next;
    }
}

inline mpf_class theta3_from_power_of_two_q(mp_bitcnt_t q_exponent, precision_type precision) {
    const mpf_class threshold = theta_series_threshold(precision);
    mpf_class sum = make_ui(1, precision);
    mpf_class term = make_ui(1, precision);
    mpf_div_2exp(term.get_mpf_t(), term.get_mpf_t(), q_exponent);

    for (unsigned long k = 1;; ++k) {
        const mpf_class contribution = mul_ui(term, 2ul, precision);
        if (contribution < threshold) {
            break;
        }
        sum = add(sum, contribution, precision);
        mpf_div_2exp(term.get_mpf_t(), term.get_mpf_t(), q_exponent * (2ul * k + 1ul));
    }
    return sum;
}

inline mpf_class theta2_from_power_of_two_q(mp_bitcnt_t q_exponent, precision_type precision) {
    const mpf_class threshold = theta_series_threshold(precision);
    mpf_class q = make_ui(1, precision);
    mpf_div_2exp(q.get_mpf_t(), q.get_mpf_t(), q_exponent);
    mpf_class term = sqrt_prec(sqrt_prec(q, precision), precision);
    mpf_class sum = make_ui(0, precision);

    for (unsigned long k = 0;; ++k) {
        const mpf_class contribution = mul_ui(term, 2ul, precision);
        if (contribution < threshold) {
            break;
        }
        sum = add(sum, contribution, precision);
        mpf_div_2exp(term.get_mpf_t(), term.get_mpf_t(), q_exponent * (2ul * k + 2ul));
    }
    return sum;
}

inline mpf_class compute_log_two_theta_agm(precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_log_two(target);
    const mp_bitcnt_t q_exponent = (work / 2) - 2;

    const mpf_class theta2 = theta2_from_power_of_two_q(q_exponent, work);
    const mpf_class theta3 = theta3_from_power_of_two_q(q_exponent, work);
    const mpf_class agm_value = agm_converged(sqr(theta2, work), sqr(theta3, work), work);
    const mpf_class q_scale = make_ui(static_cast<unsigned long>(q_exponent), work);
    const mpf_class denominator = mul(q_scale, agm_value, work);
    return set_prec_copy(div(pi(work), denominator, work), target);
}

struct log_two_cache_state {
    std::mutex mutex;
    precision_type cached_precision = 0;
    mpf_class cached_value = make_ui(0, 32);
    bool initialized = false;
};

inline log_two_cache_state &log_two_cache() {
    static log_two_cache_state cache;
    return cache;
}

inline mpf_class log_two(precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    log_two_cache_state &cache = log_two_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    if (!cache.initialized || cache.cached_precision < target) {
        cache.cached_value = compute_log_two_theta_agm(target);
        cache.cached_precision = target;
        cache.initialized = true;
    }
    return set_prec_copy(cache.cached_value, target);
}

inline precision_type guard_bits_for_log1p(precision_type) {
    return 96;
}

inline precision_type working_precision_for_log1p(precision_type target_precision) {
    return normalize_target_precision(target_precision) + guard_bits_for_log1p(target_precision);
}

inline mpf_class log1p_taylor_small(const mpf_class &x, precision_type precision) {
    mpf_class epsilon = make_ui(1, precision);
    mpf_div_2exp(epsilon.get_mpf_t(), epsilon.get_mpf_t(), precision);

    mpf_class sum = set_prec_copy(x, precision);
    mpf_class power = set_prec_copy(x, precision);

    for (unsigned long k = 2;; ++k) {
        power = mul(power, x, precision);
        mpf_class term = div(power, make_ui(k, precision), precision);
        if ((k & 1ul) == 0ul) {
            term = sub(make_ui(0, precision), term, precision);
        }
        sum = add(sum, term, precision);
        if (abs(term) < epsilon) {
            break;
        }
    }
    return sum;
}

inline mpf_class log1p_atanh_series(const mpf_class &x, precision_type precision) {
    //
    // log(1 + x) = 2 * atanh(x / (2 + x))
    //            = 2 * sum_{k>=0} y^(2k+1)/(2k+1),  y = x / (2 + x)
    //
    // This avoids the cancellation that breaks a direct log(x) path near 1.
    //
    mpf_class epsilon = make_ui(1, precision);
    mpf_div_2exp(epsilon.get_mpf_t(), epsilon.get_mpf_t(), precision);

    const mpf_class two = make_ui(2, precision);
    const mpf_class y = div(x, add(two, x, precision), precision);
    const mpf_class y2 = sqr(y, precision);

    mpf_class sum = set_prec_copy(y, precision);
    mpf_class term = set_prec_copy(y, precision);

    for (unsigned long k = 1;; ++k) {
        term = mul(term, y2, precision);
        const unsigned long denominator_ui = 2ul * k + 1ul;
        const mpf_class contribution = div(term, make_ui(denominator_ui, precision), precision);
        sum = add(sum, contribution, precision);
        if (abs(contribution) < epsilon) {
            break;
        }
    }
    return mul_ui(sum, 2ul, precision);
}

inline mpf_class compute_log1p(const mpf_class &x_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_log1p(target);
    const mpf_class x = set_prec_copy(x_input, work);
    const mpf_class one = make_ui(1, work);
    const mpf_class one_plus_x = add(one, x, work);

    if (one_plus_x == make_ui(0, work)) {
        throw std::domain_error("log1p(x) pole at x = -1");
    }
    if (one_plus_x < make_ui(0, work)) {
        throw std::domain_error("log1p(x) is undefined for x < -1");
    }
    if (x == make_ui(0, work)) {
        return make_ui(0, target);
    }

    mpf_class small_threshold = make_ui(1, work);
    mpf_div_2exp(small_threshold.get_mpf_t(), small_threshold.get_mpf_t(), work / 2);

    mpf_class result = make_ui(0, work);
    if (abs(x) <= small_threshold) {
        result = log1p_taylor_small(x, work);
    } else {
        result = log1p_atanh_series(x, work);
    }
    return set_prec_copy(result, target);
}

inline precision_type guard_bits_for_log(precision_type) {
    return 96;
}

inline precision_type log_cancellation_bits(const mpf_class &x_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type probe_precision = target + 32;
    const mpf_class one = make_ui(1, probe_precision);
    const mpf_class x = set_prec_copy(x_input, probe_precision);
    const mpf_class delta = abs(sub(x, one, probe_precision));

    if (delta == make_ui(0, probe_precision)) {
        return target;
    }

    mp_exp_t delta_exponent = 0;
    mpf_get_d_2exp(&delta_exponent, delta.get_mpf_t());
    const long numerator_bits = static_cast<long>(ceil_log2_precision(target)) + 1;
    const long estimate = numerator_bits - delta_exponent + 2;
    return estimate > 0 ? static_cast<precision_type>(estimate) : 0;
}

inline precision_type working_precision_for_log(const mpf_class &x_input, precision_type target_precision) {
    return normalize_target_precision(target_precision) + guard_bits_for_log(target_precision) + log_cancellation_bits(x_input, target_precision);
}

inline mpf_class mul_signed_exp(const mpf_class &value, mp_exp_t multiplier, precision_type precision) {
    if (multiplier == 0) {
        return make_ui(0, precision);
    }
    const unsigned long magnitude_ui = static_cast<unsigned long>(multiplier > 0 ? multiplier : -multiplier);
    const mpf_class magnitude = mul_ui(value, magnitude_ui, precision);
    if (multiplier < 0) {
        return sub(make_ui(0, precision), magnitude, precision);
    }
    return magnitude;
}

inline mpf_class compute_log(const mpf_class &x_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_log(x_input, target);
    const mpf_class x = set_prec_copy(x_input, work);
    const mpf_class zero = make_ui(0, work);
    const mpf_class one = make_ui(1, work);
    const mpf_class two = make_ui(2, work);

    if (x == zero) {
        throw std::domain_error("log(x) pole at x = 0");
    }
    if (x < zero) {
        throw std::domain_error("log(x) is undefined for x < 0");
    }
    if (x == one) {
        return make_ui(0, target);
    }

    const mpf_class delta = sub(x, one, work);
    mpf_class near_one_threshold = make_ui(1, work);
    mpf_div_2exp(near_one_threshold.get_mpf_t(), near_one_threshold.get_mpf_t(), 1);
    if (abs(delta) < near_one_threshold) {
        return set_prec_copy(compute_log1p(delta, work), target);
    }

    //
    // Variable-scale AGM reduction:
    //   choose m so that s = x * 2^m > 2^(p/2 + safety)
    //   then evaluate on b = 4 / s
    //   evaluate
    //     log(x) ~= pi / (2 * AGM(1, 4 / (x * 2^m))) - m * log(2)
    //
    // The near-1 branch stays on log1p() to avoid cancellation.
    //
    mp_exp_t x_exponent = 0;
    mpf_get_d_2exp(&x_exponent, x.get_mpf_t());
    const mp_exp_t desired_exponent = static_cast<mp_exp_t>(work / 2 + 16);
    const mp_exp_t exponent = desired_exponent - x_exponent + 1;

    mpf_class s = set_prec_copy(x, work);
    if (exponent >= 0) {
        mpf_mul_2exp(s.get_mpf_t(), s.get_mpf_t(), static_cast<mp_bitcnt_t>(exponent));
    } else {
        mpf_div_2exp(s.get_mpf_t(), s.get_mpf_t(), static_cast<mp_bitcnt_t>(-exponent));
    }

    const mpf_class b = div(make_ui(4, work), s, work);
    const mpf_class agm_value = agm_converged(one, b, work);
    const mpf_class leading = div(pi(work), mul(two, agm_value, work), work);
    const mpf_class correction = mul_signed_exp(log_two(work), exponent, work);
    return set_prec_copy(sub(leading, correction, work), target);
}

inline precision_type guard_bits_for_exp(precision_type) {
    return 96;
}

inline precision_type working_precision_for_exp(precision_type target_precision) {
    return normalize_target_precision(target_precision) + guard_bits_for_exp(target_precision);
}

inline mp_exp_t round_to_nearest_mp_exp(const mpf_class &value, precision_type precision) {
    const mpf_class zero = make_ui(0, precision);
    mpf_class half = make_ui(1, precision);
    mpf_div_2exp(half.get_mpf_t(), half.get_mpf_t(), 1);
    mpf_class adjusted = set_prec_copy(value, precision);
    if (adjusted >= zero) {
        adjusted = add(adjusted, half, precision);
    } else {
        adjusted = sub(adjusted, half, precision);
    }

    mpz_class rounded_integer;
    mpz_set_f(rounded_integer.get_mpz_t(), adjusted.get_mpf_t());
    if (!rounded_integer.fits_slong_p()) {
        throw std::overflow_error("exp(x) scaling exponent is too large");
    }
    return static_cast<mp_exp_t>(rounded_integer.get_si());
}

inline mpf_class exp_taylor_reduced(const mpf_class &x, precision_type precision) {
    mpf_class epsilon = make_ui(1, precision);
    mpf_div_2exp(epsilon.get_mpf_t(), epsilon.get_mpf_t(), precision);

    mpf_class sum = make_ui(1, precision);
    mpf_class term = make_ui(1, precision);
    for (unsigned long n = 1;; ++n) {
        term = div(mul(term, x, precision), make_ui(n, precision), precision);
        sum = add(sum, term, precision);
        if (abs(term) < epsilon) {
            break;
        }
    }
    return sum;
}

inline mpf_class compute_exp(const mpf_class &x_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_exp(target);
    const mpf_class x = set_prec_copy(x_input, work);
    const mpf_class zero = make_ui(0, work);

    if (x == zero) {
        return make_ui(1, target);
    }

    const mpf_class log2_value = log_two(work);
    const mpf_class quotient = div(x, log2_value, work);
    const mp_exp_t k = round_to_nearest_mp_exp(quotient, work);
    const mpf_class reduced = sub(x, mul_signed_exp(log2_value, k, work), work);

    mpf_class result = exp_taylor_reduced(reduced, work);
    if (k >= 0) {
        mpf_mul_2exp(result.get_mpf_t(), result.get_mpf_t(), static_cast<mp_bitcnt_t>(k));
    } else {
        mpf_div_2exp(result.get_mpf_t(), result.get_mpf_t(), static_cast<mp_bitcnt_t>(-k));
    }
    return set_prec_copy(result, target);
}

inline precision_type guard_bits_for_expm1(precision_type) {
    return 96;
}

inline precision_type working_precision_for_expm1(precision_type target_precision) {
    return normalize_target_precision(target_precision) + guard_bits_for_expm1(target_precision);
}

inline mpf_class expm1_taylor_small(const mpf_class &x, precision_type precision) {
    mpf_class epsilon = make_ui(1, precision);
    mpf_div_2exp(epsilon.get_mpf_t(), epsilon.get_mpf_t(), precision);

    mpf_class sum = set_prec_copy(x, precision);
    mpf_class term = set_prec_copy(x, precision);
    for (unsigned long n = 2;; ++n) {
        term = div(mul(term, x, precision), make_ui(n, precision), precision);
        sum = add(sum, term, precision);
        if (abs(term) < epsilon) {
            break;
        }
    }
    return sum;
}

inline mpf_class compute_expm1(const mpf_class &x_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_expm1(target);
    const mpf_class x = set_prec_copy(x_input, work);

    if (x == make_ui(0, work)) {
        return make_ui(0, target);
    }

    mpf_class small_threshold = make_ui(1, work);
    mpf_div_2exp(small_threshold.get_mpf_t(), small_threshold.get_mpf_t(), work / 2);

    mpf_class result = make_ui(0, work);
    if (abs(x) <= small_threshold) {
        result = expm1_taylor_small(x, work);
    } else {
        result = sub(compute_exp(x, work), make_ui(1, work), work);
    }
    return set_prec_copy(result, target);
}

inline precision_type guard_bits_for_trig(precision_type) {
    return 128;
}

inline precision_type working_precision_for_trig(precision_type target_precision) {
    return normalize_target_precision(target_precision) + guard_bits_for_trig(target_precision);
}

inline precision_type trig_constant_precision(precision_type target_precision) {
    return (2 * normalize_target_precision(target_precision)) + 64;
}

inline sincos_result sincos_taylor_small(const mpf_class &x, precision_type precision) {
    mpf_class epsilon = make_ui(1, precision);
    mpf_div_2exp(epsilon.get_mpf_t(), epsilon.get_mpf_t(), precision);

    sincos_result result;
    result.sin_value = set_prec_copy(x, precision);
    result.cos_value = make_ui(1, precision);

    mpf_class x2 = sqr(x, precision);
    mpf_class sin_term = set_prec_copy(x, precision);
    mpf_class cos_term = make_ui(1, precision);

    for (unsigned long k = 1;; ++k) {
        const unsigned long sin_den1 = 2ul * k;
        const unsigned long sin_den2 = 2ul * k + 1ul;
        sin_term = div(mul(sin_term, x2, precision), make_ui(sin_den1 * sin_den2, precision), precision);
        sin_term = sub(make_ui(0, precision), sin_term, precision);
        result.sin_value = add(result.sin_value, sin_term, precision);

        const unsigned long cos_den1 = 2ul * k - 1ul;
        const unsigned long cos_den2 = 2ul * k;
        cos_term = div(mul(cos_term, x2, precision), make_ui(cos_den1 * cos_den2, precision), precision);
        cos_term = sub(make_ui(0, precision), cos_term, precision);
        result.cos_value = add(result.cos_value, cos_term, precision);

        if (abs(sin_term) < epsilon && abs(cos_term) < epsilon) {
            break;
        }
    }
    return result;
}

inline void ensure_trig_constants(precision_type target_precision) {
    const precision_type cache_precision = trig_constant_precision(target_precision);
    trig_constant_cache_state &cache = trig_constant_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    if (!cache.initialized || cache.cached_precision < cache_precision) {
        cache.pi_value = compute_pi_gauss_legendre(cache_precision);
        cache.pi_over_two_value = set_prec_copy(cache.pi_value, cache_precision);
        mpf_div_2exp(cache.pi_over_two_value.get_mpf_t(), cache.pi_over_two_value.get_mpf_t(), 1);
        cache.two_over_pi_value = div(make_ui(2, cache_precision), cache.pi_value, cache_precision);
        cache.cached_precision = cache_precision;
        cache.initialized = true;
    }
}

inline mpf_class trig_pi_over_two(precision_type target_precision) {
    ensure_trig_constants(target_precision);
    trig_constant_cache_state &cache = trig_constant_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    return set_prec_copy(cache.pi_over_two_value, trig_constant_precision(target_precision));
}

inline mpf_class trig_two_over_pi(precision_type target_precision) {
    ensure_trig_constants(target_precision);
    trig_constant_cache_state &cache = trig_constant_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    return set_prec_copy(cache.two_over_pi_value, trig_constant_precision(target_precision));
}

inline sincos_result compute_sincos(const mpf_class &x_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_trig(target);
    const mpf_class x = set_prec_copy(x_input, work);
    const mpf_class zero = make_ui(0, work);
    const precision_type const_precision = trig_constant_precision(target);
    const mpf_class pio2 = set_prec_copy(trig_pi_over_two(target), const_precision);
    const mpf_class two_over_pi = set_prec_copy(trig_two_over_pi(target), const_precision);

    const mpf_class scaled_x = set_prec_copy(x_input, const_precision);
    const mpf_class q = mul(scaled_x, two_over_pi, const_precision);
    mpz_class k = 0;
    mpz_set_f(k.get_mpz_t(), q.get_mpf_t());
    const mpf_class integer_part(k, const_precision);
    const mpf_class frac = sub(q, integer_part, const_precision);
    const mpf_class half = div(make_ui(1, const_precision), make_ui(2, const_precision), const_precision);
    const mpf_class neg_half = sub(make_ui(0, const_precision), half, const_precision);
    const int frac_vs_half = mpf_cmp(frac.get_mpf_t(), half.get_mpf_t());
    const int frac_vs_neg_half = mpf_cmp(frac.get_mpf_t(), neg_half.get_mpf_t());
    if (frac_vs_half > 0) {
        k += 1;
    } else if (frac_vs_half == 0) {
        if (mpz_odd_p(k.get_mpz_t())) {
            k += 1;
        }
    } else if (frac_vs_neg_half < 0) {
        k -= 1;
    } else if (frac_vs_neg_half == 0) {
        if (mpz_odd_p(k.get_mpz_t())) {
            k -= 1;
        }
    }

    const mpf_class remainder_hi = sub(scaled_x, mul(pio2, mpf_class(k, const_precision), const_precision), const_precision);
    const mpf_class reduced_argument = set_prec_copy(remainder_hi, work);
    const sincos_result base = sincos_taylor_small(reduced_argument, work);

    const unsigned long quadrant = mpz_fdiv_ui(k.get_mpz_t(), 4ul);
    sincos_result result;
    switch (quadrant) {
    case 0:
        result.sin_value = set_prec_copy(base.sin_value, work);
        result.cos_value = set_prec_copy(base.cos_value, work);
        break;
    case 1:
        result.sin_value = set_prec_copy(base.cos_value, work);
        result.cos_value = set_prec_copy(sub(zero, base.sin_value, work), work);
        break;
    case 2:
        result.sin_value = set_prec_copy(sub(zero, base.sin_value, work), work);
        result.cos_value = set_prec_copy(sub(zero, base.cos_value, work), work);
        break;
    default:
        result.sin_value = set_prec_copy(sub(zero, base.cos_value, work), work);
        result.cos_value = set_prec_copy(base.sin_value, work);
        break;
    }

    result.sin_value = set_prec_copy(result.sin_value, target);
    result.cos_value = set_prec_copy(result.cos_value, target);
    return result;
}

inline mpf_class compute_sin(const mpf_class &x_input, precision_type target_precision) {
    return compute_sincos(x_input, target_precision).sin_value;
}

inline mpf_class compute_cos(const mpf_class &x_input, precision_type target_precision) {
    return compute_sincos(x_input, target_precision).cos_value;
}

inline precision_type guard_bits_for_atan(precision_type) {
    return 96;
}

inline precision_type working_precision_for_atan(precision_type target_precision) {
    return normalize_target_precision(target_precision) + guard_bits_for_atan(target_precision);
}

inline mpf_class atan_taylor_small(const mpf_class &x, precision_type precision) {
    mpf_class epsilon = make_ui(1, precision);
    mpf_div_2exp(epsilon.get_mpf_t(), epsilon.get_mpf_t(), precision);

    const mpf_class x2 = sqr(x, precision);
    mpf_class sum = set_prec_copy(x, precision);
    mpf_class power = set_prec_copy(x, precision);

    for (unsigned long k = 1;; ++k) {
        power = mul(power, x2, precision);
        mpf_class contribution = div(power, make_ui(2ul * k + 1ul, precision), precision);
        if ((k & 1ul) == 1ul) {
            contribution = sub(make_ui(0, precision), contribution, precision);
        }
        sum = add(sum, contribution, precision);
        if (abs(contribution) < epsilon) {
            break;
        }
    }
    return sum;
}

inline mpf_class compute_atan(const mpf_class &x_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_atan(target);
    const mpf_class x = set_prec_copy(x_input, work);
    const mpf_class zero = make_ui(0, work);
    const mpf_class one = make_ui(1, work);

    if (x == zero) {
        return make_ui(0, target);
    }

    const bool negate = x < zero;
    const mpf_class ax = negate ? sub(zero, x, work) : x;

    mpf_class y = set_prec_copy(ax, work);
    unsigned long reductions = 0;
    const mpf_class threshold = div(one, make_ui(2, work), work);
    while (y > threshold) {
        const mpf_class sqrt_term = sqrt_prec(add(one, sqr(y, work), work), work);
        y = div(y, add(one, sqrt_term, work), work);
        reductions++;
    }

    mpf_class result = atan_taylor_small(y, work);
    for (unsigned long i = 0; i < reductions; ++i) {
        result = mul_ui(result, 2ul, work);
    }

    if (negate) {
        result = sub(zero, result, work);
    }
    return set_prec_copy(result, target);
}

inline mpf_class compute_atan2(const mpf_class &y_input, const mpf_class &x_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_atan(target);
    const mpf_class y = set_prec_copy(y_input, work);
    const mpf_class x = set_prec_copy(x_input, work);
    const mpf_class zero = make_ui(0, work);

    // GMP mpf_class does not expose signed zero. Adopt a simple fallback
    // convention:
    //   atan2(0, x>0) = 0
    //   atan2(0, x<0) = pi
    //   atan2(y>0, 0) = pi/2
    //   atan2(y<0, 0) = -pi/2
    //   atan2(0, 0)   = 0
    if (y == zero) {
        if (x > zero) {
            return make_ui(0, target);
        }
        if (x < zero) {
            return set_prec_copy(pi(target + 2), target);
        }
        return make_ui(0, target);
    }

    if (x == zero) {
        mpf_class pio2 = pi(target + 2);
        mpf_div_2exp(pio2.get_mpf_t(), pio2.get_mpf_t(), 1);
        if (y > zero) {
            return set_prec_copy(pio2, target);
        }
        return set_prec_copy(sub(make_ui(0, target + 2), pio2, target + 2), target);
    }

    const mpf_class ratio = div(y, x, work);
    mpf_class result = compute_atan(ratio, work);
    const mpf_class pi_value = pi(work);
    if (x < zero) {
        if (y > zero) {
            result = add(result, pi_value, work);
        } else {
            result = sub(result, pi_value, work);
        }
    }
    return set_prec_copy(result, target);
}

} // namespace mplapack_gmp_transcendents

#endif
