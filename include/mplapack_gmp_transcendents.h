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

#include "gmpxx.h"

namespace mplapack_gmp_transcendents {

using precision_type = mp_bitcnt_t;

inline precision_type normalize_target_precision(precision_type target_precision) {
    return std::max<precision_type>(target_precision, 32);
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

} // namespace mplapack_gmp_transcendents

#endif
