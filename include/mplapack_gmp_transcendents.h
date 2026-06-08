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
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "gmpxx.h"

namespace mplapack_gmp_transcendents {

using precision_type = mp_bitcnt_t;

inline constexpr precision_type minimum_target_precision = 32;
inline constexpr precision_type log_cancellation_probe_guard_bits = 32;

inline precision_type normalize_target_precision(precision_type target_precision) {
    return std::max<precision_type>(target_precision, minimum_target_precision);
}

inline unsigned long ceil_log2_precision(precision_type value) {
    if (value <= 1) {
        return 0;
    }
    --value;
    unsigned long bits = 0;
    while (value != 0) {
        value >>= 1;
        ++bits;
    }
    return bits;
}

inline std::uintmax_t mpf_mp_exp_negative_magnitude(mp_exp_t value) noexcept {
    return static_cast<std::uintmax_t>(-(value + 1)) + std::uintmax_t{1};
}

inline void ensure_mpf_shift_result_exponent_fits(mpf_srcptr value, mp_bitcnt_t bits, bool left_shift, const char *message) {
    if (mpf_sgn(value) == 0 || bits == 0) {
        return;
    }

    mp_exp_t exponent = 0;
    mpf_get_d_2exp(&exponent, value);
    std::uintmax_t allowed = 0;
    if (left_shift) {
        if (exponent >= 0) {
            allowed = static_cast<std::uintmax_t>(std::numeric_limits<mp_exp_t>::max() - exponent);
        } else {
            allowed = static_cast<std::uintmax_t>(std::numeric_limits<mp_exp_t>::max()) + mpf_mp_exp_negative_magnitude(exponent);
        }
    } else {
        if (exponent >= 0) {
            allowed = static_cast<std::uintmax_t>(exponent) + mpf_mp_exp_negative_magnitude(std::numeric_limits<mp_exp_t>::min());
        } else {
            const std::uintmax_t exponent_magnitude = mpf_mp_exp_negative_magnitude(exponent);
            const std::uintmax_t min_magnitude = mpf_mp_exp_negative_magnitude(std::numeric_limits<mp_exp_t>::min());
            allowed = min_magnitude - exponent_magnitude;
        }
    }

    if (static_cast<std::uintmax_t>(bits) > allowed) {
        throw std::overflow_error(message);
    }
}

inline std::uint64_t iteration_limit_for_precision(precision_type precision) {
    constexpr std::uint64_t minimum_iterations = 64;
    constexpr std::uint64_t precision_multiplier = 16;
    const std::uintmax_t precision_value = static_cast<std::uintmax_t>(precision);
    const std::uintmax_t max_value = static_cast<std::uintmax_t>(std::numeric_limits<std::uint64_t>::max());
    if (precision_value > max_value / precision_multiplier) {
        throw std::overflow_error("MPF iteration limit exceeds uint64_t");
    }
    return std::max(minimum_iterations, static_cast<std::uint64_t>(precision_value * precision_multiplier));
}

[[noreturn]] inline void throw_iteration_limit_exceeded(const char *function_name) {
    throw std::runtime_error(std::string(function_name) + " failed to converge within iteration limit");
}

inline std::uint64_t checked_taylor_counter_product(std::uint64_t lhs, std::uint64_t rhs) {
    if (lhs != 0 && rhs > std::numeric_limits<std::uint64_t>::max() / lhs) {
        throw std::overflow_error("MPF Taylor denominator exceeds uint64_t");
    }
    return lhs * rhs;
}

inline std::uint64_t checked_taylor_counter_add(std::uint64_t lhs, std::uint64_t rhs) {
    if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
        throw std::overflow_error("MPF Taylor denominator exceeds uint64_t");
    }
    return lhs + rhs;
}

inline std::uint64_t checked_taylor_odd_denominator(std::uint64_t k) {
    return checked_taylor_counter_add(checked_taylor_counter_product(std::uint64_t{2}, k), std::uint64_t{1});
}

inline mp_bitcnt_t checked_taylor_shift_count(mp_bitcnt_t scale, std::uint64_t factor) {
    const std::uint64_t shift = checked_taylor_counter_product(static_cast<std::uint64_t>(scale), factor);
    if (shift > static_cast<std::uint64_t>(std::numeric_limits<mp_bitcnt_t>::max())) {
        throw std::overflow_error("MPF Taylor shift count exceeds mp_bitcnt_t");
    }
    return static_cast<mp_bitcnt_t>(shift);
}

inline void increment_taylor_counter(std::uint64_t &counter) {
    if (counter == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("MPF Taylor counter exceeds uint64_t");
    }
    ++counter;
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

inline mpf_class make_u64(std::uint64_t value, precision_type precision) {
    mpf_class result(0, precision);
    mpz_t integer;
    mpz_init(integer);
    mpz_import(integer, 1, -1, sizeof(value), 0, 0, &value);
    mpf_set_z(result.get_mpf_t(), integer);
    mpz_clear(integer);
    return result;
}

inline mpf_class make_decimal(const char *value, precision_type precision) {
    mpf_class result(0, precision);
    if (mpf_set_str(result.get_mpf_t(), value, 10) != 0) {
        throw std::invalid_argument("invalid MPF decimal literal");
    }
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

inline mpf_class negate(const mpf_class &a, precision_type precision) {
    mpf_class result(0, precision);
    mpf_neg(result.get_mpf_t(), a.get_mpf_t());
    return result;
}

inline mpf_class neg_prec(const mpf_class &a, precision_type precision) {
    return negate(a, precision);
}

inline mpf_class abs_prec(const mpf_class &a, precision_type precision) {
    mpf_class result(0, precision);
    mpf_abs(result.get_mpf_t(), a.get_mpf_t());
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
    if (mpf_sgn(b.get_mpf_t()) == 0) {
        throw std::domain_error("mpf division by zero");
    }
    mpf_class result(0, precision);
    mpf_div(result.get_mpf_t(), a.get_mpf_t(), b.get_mpf_t());
    return result;
}

inline mpf_class div_ui(const mpf_class &a, unsigned long value, precision_type precision) {
    if (value == 0ul) {
        throw std::domain_error("mpf division by zero");
    }
    mpf_class result(0, precision);
    mpf_div_ui(result.get_mpf_t(), a.get_mpf_t(), value);
    return result;
}

inline mpf_class mul_ui(const mpf_class &a, unsigned long value, precision_type precision) {
    mpf_class result(0, precision);
    mpf_mul_ui(result.get_mpf_t(), a.get_mpf_t(), value);
    return result;
}

inline mp_bitcnt_t checked_mp_exp_magnitude(mp_exp_t value) {
    const std::uintmax_t bit_limit = static_cast<std::uintmax_t>(std::numeric_limits<mp_bitcnt_t>::max());
    std::uintmax_t magnitude = 0;
    if (value >= 0) {
        magnitude = static_cast<std::uintmax_t>(value);
    } else {
        magnitude = static_cast<std::uintmax_t>(-(value + 1));
        ++magnitude;
    }
    if (magnitude > bit_limit) {
        throw std::overflow_error("MPF exponent shift exceeds mp_bitcnt_t");
    }
    return static_cast<mp_bitcnt_t>(magnitude);
}

inline mp_bitcnt_t mp_exp_magnitude(mp_exp_t value) {
    return checked_mp_exp_magnitude(value);
}

inline mp_exp_t checked_mp_bitcnt_to_mp_exp(mp_bitcnt_t value, const char *message) {
    if (static_cast<std::uintmax_t>(value) > static_cast<std::uintmax_t>(std::numeric_limits<mp_exp_t>::max())) {
        throw std::overflow_error(message);
    }
    return static_cast<mp_exp_t>(value);
}

inline void mpz_set_mp_exp(mpz_t target, mp_exp_t value) {
    if constexpr (std::numeric_limits<mp_exp_t>::min() >= std::numeric_limits<long>::min() &&
                  std::numeric_limits<mp_exp_t>::max() <= std::numeric_limits<long>::max()) {
        mpz_set_si(target, static_cast<long>(value));
    } else {
        using unsigned_exp_t = std::make_unsigned_t<mp_exp_t>;
        unsigned_exp_t magnitude = 0;
        if (value >= 0) {
            magnitude = static_cast<unsigned_exp_t>(value);
        } else {
            magnitude = static_cast<unsigned_exp_t>(-(value + 1));
            ++magnitude;
        }
        mpz_import(target, 1, -1, sizeof(magnitude), 0, 0, &magnitude);
        if (value < 0) {
            mpz_neg(target, target);
        }
    }
}

inline mp_exp_t checked_log_scaling_exponent(mp_exp_t desired_exponent, mp_exp_t x_exponent) {
    mpz_t value;
    mpz_t x_value;
    mpz_init(value);
    mpz_init(x_value);
    mpz_set_mp_exp(value, desired_exponent);
    mpz_set_mp_exp(x_value, x_exponent);
    mpz_sub(value, value, x_value);
    mpz_add_ui(value, value, 1);
    if (!mpz_fits_slong_p(value)) {
        mpz_clear(x_value);
        mpz_clear(value);
        throw std::overflow_error("log(x) scaling exponent is too large");
    }
    const long result = mpz_get_si(value);
    mpz_clear(x_value);
    mpz_clear(value);
    if constexpr (sizeof(mp_exp_t) < sizeof(long)) {
        if (result < static_cast<long>(std::numeric_limits<mp_exp_t>::min()) ||
            result > static_cast<long>(std::numeric_limits<mp_exp_t>::max())) {
            throw std::overflow_error("log(x) scaling exponent is too large");
        }
    }
    return static_cast<mp_exp_t>(result);
}

inline void ensure_exp_scaling_exponent_fits(const mpf_class &value, mp_exp_t shift) {
    mp_exp_t exponent = 0;
    mpf_get_d_2exp(&exponent, value.get_mpf_t());
    if (shift > 0) {
        if (exponent > std::numeric_limits<mp_exp_t>::max() - shift) {
            throw std::overflow_error("exp(x) result exponent is too large");
        }
    } else if (shift < 0) {
        if (exponent < std::numeric_limits<mp_exp_t>::min() - shift) {
            throw std::overflow_error("exp(x) result exponent is too small");
        }
    }
}

inline void ensure_mpf_mul_exponent_fits(const mpf_class &lhs, const mpf_class &rhs) {
    if (mpf_sgn(lhs.get_mpf_t()) == 0 || mpf_sgn(rhs.get_mpf_t()) == 0) {
        return;
    }

    mp_exp_t lhs_exponent = 0;
    mp_exp_t rhs_exponent = 0;
    mpf_get_d_2exp(&lhs_exponent, lhs.get_mpf_t());
    mpf_get_d_2exp(&rhs_exponent, rhs.get_mpf_t());

    mpz_t sum;
    mpz_t rhs_value;
    mpz_init(sum);
    mpz_init(rhs_value);
    mpz_set_mp_exp(sum, lhs_exponent);
    mpz_set_mp_exp(rhs_value, rhs_exponent);
    mpz_add(sum, sum, rhs_value);
    mpz_clear(rhs_value);

    mpz_t min_value;
    mpz_t max_value;
    mpz_init(min_value);
    mpz_init(max_value);
    mpz_set_mp_exp(min_value, std::numeric_limits<mp_exp_t>::min());
    mpz_set_mp_exp(max_value, std::numeric_limits<mp_exp_t>::max());
    const bool too_large = mpz_cmp(sum, max_value) > 0;
    mpz_sub_ui(sum, sum, 1);
    const bool too_small = mpz_cmp(sum, min_value) < 0;
    mpz_clear(max_value);
    mpz_clear(min_value);
    mpz_clear(sum);

    if (too_large) {
        throw std::overflow_error("mpf multiplication result exponent is too large");
    }
    if (too_small) {
        throw std::overflow_error("mpf multiplication result exponent is too small");
    }
}

inline precision_type positive_exponent_bits(const mpf_class &value) {
    if (mpf_sgn(value.get_mpf_t()) == 0) {
        return 0;
    }
    mp_exp_t exponent = 0;
    mpf_get_d_2exp(&exponent, value.get_mpf_t());
    return exponent > 0 ? static_cast<precision_type>(exponent) : 0;
}

struct quo_rem_result {
    mpz_class quotient;
    mpf_class remainder = make_ui(0, 32);
};

struct sincos_result {
    mpf_class sin_value = make_ui(0, 32);
    mpf_class cos_value = make_ui(1, 32);
};

struct trig_constants_result {
    mpf_class pi_over_two_value = make_ui(0, 32);
    mpf_class two_over_pi_value = make_ui(0, 32);
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
    if (mpf_sgn(a.get_mpf_t()) < 0) {
        throw std::domain_error("mpf square root of negative value");
    }
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
    if (value == 0ul) {
        throw std::domain_error("mpf inverse square root of zero");
    }
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

// Computes the nearest-integer quotient k and remainder r such that
// x_input = k * y_input + r, with ties rounded to even k.
//
// Work precision policy: target + exponent_bits(x_input) + 64. This protects
// the final subtraction x - k*y from cancellation when |x_input| is large.
// The input value itself must already carry enough precision for its magnitude;
// quo_rem cannot recover bits that were not present in x_input or y_input.
inline quo_rem_result quo_rem(const mpf_class &x_input, const mpf_class &y_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = target + positive_exponent_bits(x_input) + 64;
    const mpf_class x = set_prec_copy(x_input, work);
    const mpf_class y = set_prec_copy(y_input, work);
    if (mpf_sgn(y.get_mpf_t()) == 0) {
        throw std::domain_error("quo_rem division by zero");
    }

    const mpf_class q = div(x, y, work);
    mpz_class k = 0;
    mpz_set_f(k.get_mpz_t(), q.get_mpf_t());

    const mpf_class integer_part(k, work);
    const mpf_class frac = sub(q, integer_part, work);
    const mpf_class half = div_ui(make_ui(1, work), 2ul, work);
    const mpf_class neg_half = negate(half, work);

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
    mpf_class pi_previous = zero;
    mpf_class pi_current = zero;

    const std::uint64_t max_iterations = iteration_limit_for_precision(work);
    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        mpf_class a_next = average(a, b, work);
        mpf_class b_next = sqrt_prec(mul(a, b, work), work);
        const mpf_class diff = sub(a, a_next, work);
        mpf_class t_next = sub(t, mul(p, sqr(diff, work), work), work);

        a = a_next;
        b = b_next;
        t = t_next;
        p = mul(p, two, work);

        const mpf_class sum = add(a, b, work);
        const mpf_class numerator = sqr(sum, work);
        const mpf_class denominator = mul(t, four, work);

        pi_previous = pi_current;
        pi_current = div(numerator, denominator, work);
        if (abs_prec(sub(pi_current, pi_previous, work), work) < epsilon) {
            return set_prec_copy(pi_current, target);
        }
    }

    throw_iteration_limit_exceeded("compute_pi_gauss_legendre");
}

inline const char *pi_decimal_literal() {
    return "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"
           "82148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964"
           "42881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127372"
           "45870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330"
           "57270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949129833"
           "67336244065664308602139494639522473719070217986094370277053921717629317675238467481846766940513200056"
           "81271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199"
           "56112129021960864034418159813629774771309960518707211349999998372978049951059731732816096318595024459"
           "45534690830264252230825334468503526193118817101000313783875288658753320838142061717766914730359825349"
           "0428755468731159562863882353787593751957781857780532171226806613001927876611195909216420199";
}

inline bool has_hardcoded_pi(precision_type target_precision) {
    constexpr precision_type decimal_digits = 1000;
    constexpr precision_type fractional_digits = decimal_digits - 1;
    constexpr precision_type conservative_bits = (fractional_digits * 332) / 100;
    return target_precision <= conservative_bits;
}

inline mpf_class hardcoded_pi(precision_type target_precision) {
    return make_decimal(pi_decimal_literal(), target_precision);
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
    if (has_hardcoded_pi(target)) {
        return hardcoded_pi(target);
    }

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

inline const char *log_two_decimal_literal() {
    return "0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875"
           "42001481020570685733685520235758130557032670751635075961930727570828371435190307038623891673471123350"
           "11536449795523912047517268157493206515552473413952588295045300709532636664265410423915781495204374043"
           "03855008019441706416715186447128399681717845469570262716310645461502572074024816377733896385506952606"
           "68341137273873722928956493547025762652098859693201965058554764703306793654432547632744951250406069438"
           "14710468994650622016772042452452961268794654619316517468139267250410380254625965686914419287160829380"
           "31727143677826548775664850856740776484514644399404614226031930967354025744460703080960850474866385231"
           "38181676751438667476647890881437141985494231519973548803751658612753529166100071053558249879414729509"
           "29311389715599820565439287170007218085761025236889213244971389320378439353088774825970171559107088236"
           "83627589842589185353024363421436706118923678919237231467232172053401649256872747782344535348";
}

inline bool has_hardcoded_log_two(precision_type target_precision) {
    constexpr precision_type decimal_digits = 1000;
    constexpr precision_type fractional_digits = decimal_digits - 1;
    constexpr precision_type conservative_bits = (fractional_digits * 332) / 100;
    return target_precision <= conservative_bits;
}

inline mpf_class hardcoded_log_two(precision_type target_precision) {
    return make_decimal(log_two_decimal_literal(), target_precision);
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

    const std::uint64_t max_iterations = iteration_limit_for_precision(precision);
    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        const mpf_class a_next = average(a, b, precision);
        const mpf_class b_next = sqrt_prec(mul(a, b, precision), precision);
        if (abs_prec(sub(a_next, b_next, precision), precision) < epsilon) {
            return average(a_next, b_next, precision);
        }
        a = a_next;
        b = b_next;
    }
    throw_iteration_limit_exceeded("agm_converged");
}

inline mpf_class theta3_from_power_of_two_q(mp_bitcnt_t q_exponent, precision_type precision) {
    const mpf_class threshold = theta_series_threshold(precision);
    mpf_class sum = make_ui(1, precision);
    mpf_class term = make_ui(1, precision);
    mpf_div_2exp(term.get_mpf_t(), term.get_mpf_t(), q_exponent);

    std::uint64_t k = 1;
    const std::uint64_t max_iterations = iteration_limit_for_precision(precision);
    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        const mpf_class contribution = mul_ui(term, 2ul, precision);
        if (contribution < threshold) {
            return sum;
        }
        sum = add(sum, contribution, precision);
        mpf_div_2exp(term.get_mpf_t(), term.get_mpf_t(), checked_taylor_shift_count(q_exponent, checked_taylor_odd_denominator(k)));
        increment_taylor_counter(k);
    }
    throw_iteration_limit_exceeded("theta3_from_power_of_two_q");
}

inline mpf_class theta2_from_power_of_two_q(mp_bitcnt_t q_exponent, precision_type precision) {
    const mpf_class threshold = theta_series_threshold(precision);
    mpf_class q = make_ui(1, precision);
    mpf_div_2exp(q.get_mpf_t(), q.get_mpf_t(), q_exponent);
    mpf_class term = sqrt_prec(sqrt_prec(q, precision), precision);
    mpf_class sum = make_ui(0, precision);

    std::uint64_t k = 0;
    const std::uint64_t max_iterations = iteration_limit_for_precision(precision);
    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        const mpf_class contribution = mul_ui(term, 2ul, precision);
        if (contribution < threshold) {
            return sum;
        }
        sum = add(sum, contribution, precision);
        const std::uint64_t factor = checked_taylor_counter_add(checked_taylor_counter_product(std::uint64_t{2}, k), std::uint64_t{2});
        mpf_div_2exp(term.get_mpf_t(), term.get_mpf_t(), checked_taylor_shift_count(q_exponent, factor));
        increment_taylor_counter(k);
    }
    throw_iteration_limit_exceeded("theta2_from_power_of_two_q");
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

// Cache dependency note:
// log_two() computes a cache miss while holding log_two_cache.mutex, and the
// miss path calls pi(work). This is safe only while pi_cache remains a leaf:
// compute_pi_gauss_legendre() must not call log_two() or trig_constants().
// trig_constants() also avoids pi_cache and computes its pi value directly, so
// there is no current cache-lock cycle. If these dependencies change, replace
// the miss path with an in-progress/condition-variable scheme rather than
// adding a reverse cache call under an existing lock.
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
    if (has_hardcoded_log_two(target)) {
        return hardcoded_log_two(target);
    }

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
    return 160;
}

inline precision_type working_precision_for_log1p(precision_type target_precision) {
    return normalize_target_precision(target_precision) + guard_bits_for_log1p(target_precision);
}

inline mpf_class log1p_taylor_small(const mpf_class &x, precision_type precision) {
    mpf_class epsilon = make_ui(1, precision);
    mpf_div_2exp(epsilon.get_mpf_t(), epsilon.get_mpf_t(), precision);

    mpf_class sum = set_prec_copy(x, precision);
    mpf_class power = set_prec_copy(x, precision);

    std::uint64_t k = 2;
    const std::uint64_t max_iterations = iteration_limit_for_precision(precision);
    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        power = mul(power, x, precision);
        mpf_class term = div(power, make_u64(k, precision), precision);
        if ((k & std::uint64_t{1}) == std::uint64_t{0}) {
            term = negate(term, precision);
        }
        sum = add(sum, term, precision);
        if (abs_prec(term, precision) < epsilon) {
            return sum;
        }
        increment_taylor_counter(k);
    }
    throw_iteration_limit_exceeded("log1p_taylor_small");
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

    std::uint64_t k = 1;
    const std::uint64_t max_iterations = iteration_limit_for_precision(precision);
    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        term = mul(term, y2, precision);
        const std::uint64_t denominator = checked_taylor_odd_denominator(k);
        const mpf_class contribution = div(term, make_u64(denominator, precision), precision);
        sum = add(sum, contribution, precision);
        if (abs_prec(contribution, precision) < epsilon) {
            return mul_ui(sum, 2ul, precision);
        }
        increment_taylor_counter(k);
    }
    throw_iteration_limit_exceeded("log1p_atanh_series");
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
    const precision_type probe_precision = target + log_cancellation_probe_guard_bits;
    const mpf_class one = make_ui(1, probe_precision);
    const mpf_class x = set_prec_copy(x_input, probe_precision);
    const mpf_class delta = abs_prec(sub(x, one, probe_precision), probe_precision);

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
    const mpf_class factor = make_u64(static_cast<std::uint64_t>(checked_mp_exp_magnitude(multiplier)), precision);
    const mpf_class magnitude = mul(value, factor, precision);
    if (multiplier < 0) {
        return negate(magnitude, precision);
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
        return compute_log1p(delta, target);
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
    const precision_type desired_exponent_bits = work / 2;
    if (desired_exponent_bits > std::numeric_limits<precision_type>::max() - 16) {
        throw std::overflow_error("log(x) scaling exponent is too large");
    }
    const mp_exp_t desired_exponent = checked_mp_bitcnt_to_mp_exp(desired_exponent_bits + 16, "log(x) scaling exponent is too large");
    const mp_exp_t exponent = checked_log_scaling_exponent(desired_exponent, x_exponent);

    mpf_class s = set_prec_copy(x, work);
    if (exponent >= 0) {
        const mp_bitcnt_t shift = checked_mp_exp_magnitude(exponent);
        ensure_mpf_shift_result_exponent_fits(s.get_mpf_t(), shift, true, "log(x) scaling exponent is too large");
        mpf_mul_2exp(s.get_mpf_t(), s.get_mpf_t(), shift);
    } else {
        const mp_bitcnt_t shift = checked_mp_exp_magnitude(exponent);
        ensure_mpf_shift_result_exponent_fits(s.get_mpf_t(), shift, false, "log(x) scaling exponent is too small");
        mpf_div_2exp(s.get_mpf_t(), s.get_mpf_t(), shift);
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
    const mpf_class adjusted = set_prec_copy(value, precision);

    mpz_class rounded_integer;
    mpz_set_f(rounded_integer.get_mpz_t(), adjusted.get_mpf_t());

    mpf_class integer_part(0, precision);
    mpf_set_z(integer_part.get_mpf_t(), rounded_integer.get_mpz_t());
    const mpf_class frac = sub(adjusted, integer_part, precision);

    mpf_class half = make_ui(1, precision);
    mpf_div_2exp(half.get_mpf_t(), half.get_mpf_t(), 1);
    const mpf_class neg_half = negate(half, precision);
    const int frac_vs_half = mpf_cmp(frac.get_mpf_t(), half.get_mpf_t());
    const int frac_vs_neg_half = mpf_cmp(frac.get_mpf_t(), neg_half.get_mpf_t());
    if (frac_vs_half > 0 || (frac_vs_half == 0 && mpz_odd_p(rounded_integer.get_mpz_t()))) {
        rounded_integer += 1;
    } else if (frac_vs_neg_half < 0 || (frac_vs_neg_half == 0 && mpz_odd_p(rounded_integer.get_mpz_t()))) {
        rounded_integer -= 1;
    }

    if (!rounded_integer.fits_slong_p()) {
        throw std::overflow_error("exp(x) scaling exponent is too large");
    }
    const long rounded = rounded_integer.get_si();
    if constexpr (sizeof(mp_exp_t) < sizeof(long)) {
        if (rounded < static_cast<long>(std::numeric_limits<mp_exp_t>::min()) ||
            rounded > static_cast<long>(std::numeric_limits<mp_exp_t>::max())) {
            throw std::overflow_error("exp(x) scaling exponent is too large");
        }
    }
    return static_cast<mp_exp_t>(rounded);
}

inline mpf_class exp_taylor_reduced(const mpf_class &x, precision_type precision) {
    mpf_class epsilon = make_ui(1, precision);
    mpf_div_2exp(epsilon.get_mpf_t(), epsilon.get_mpf_t(), precision);

    mpf_class sum = make_ui(1, precision);
    mpf_class term = make_ui(1, precision);
    std::uint64_t n = 1;
    const std::uint64_t max_iterations = iteration_limit_for_precision(precision);
    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        term = div(mul(term, x, precision), make_u64(n, precision), precision);
        sum = add(sum, term, precision);
        if (abs_prec(term, precision) < epsilon) {
            return sum;
        }
        increment_taylor_counter(n);
    }
    throw_iteration_limit_exceeded("exp_taylor_reduced");
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
    ensure_exp_scaling_exponent_fits(result, k);
    if (k >= 0) {
        mpf_mul_2exp(result.get_mpf_t(), result.get_mpf_t(), checked_mp_exp_magnitude(k));
    } else {
        mpf_div_2exp(result.get_mpf_t(), result.get_mpf_t(), checked_mp_exp_magnitude(k));
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
    std::uint64_t n = 2;
    const std::uint64_t max_iterations = iteration_limit_for_precision(precision);
    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        term = div(mul(term, x, precision), make_u64(n, precision), precision);
        sum = add(sum, term, precision);
        if (abs_prec(term, precision) < epsilon) {
            return sum;
        }
        increment_taylor_counter(n);
    }
    throw_iteration_limit_exceeded("expm1_taylor_small");
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
    const precision_type target = normalize_target_precision(target_precision);
    if (target > (std::numeric_limits<precision_type>::max() - 64) / 2) {
        throw std::overflow_error("MPF trigonometric constant precision exceeds mp_bitcnt_t");
    }
    return (2 * target) + 64;
}

inline precision_type checked_trig_precision_add(precision_type lhs, precision_type rhs) {
    if (rhs > std::numeric_limits<precision_type>::max() - lhs) {
        throw std::overflow_error("MPF trigonometric argument reduction precision exceeds mp_bitcnt_t");
    }
    return lhs + rhs;
}

inline precision_type trig_argument_exponent_bits(const mpf_class &x_input) {
    return positive_exponent_bits(x_input);
}

inline precision_type trig_constant_precision_for_argument(const mpf_class &x_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    precision_type precision = trig_constant_precision(target);
    if (mpf_sgn(x_input.get_mpf_t()) == 0) {
        return precision;
    }

    mp_exp_t x_exponent = 0;
    mpf_get_d_2exp(&x_exponent, x_input.get_mpf_t());
    if (x_exponent <= 0) {
        return precision;
    }

    precision_type required = checked_trig_precision_add(target, guard_bits_for_trig(target));
    required = checked_trig_precision_add(required, checked_mp_exp_magnitude(x_exponent));
    return std::max(precision, required);
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

    std::uint64_t k = 1;
    const std::uint64_t max_iterations = iteration_limit_for_precision(precision);
    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        const std::uint64_t sin_den1 = checked_taylor_counter_product(std::uint64_t{2}, k);
        const std::uint64_t sin_den2 = checked_taylor_counter_add(sin_den1, std::uint64_t{1});
        const std::uint64_t sin_den = checked_taylor_counter_product(sin_den1, sin_den2);
        sin_term = div(mul(sin_term, x2, precision), make_u64(sin_den, precision), precision);
        sin_term = negate(sin_term, precision);
        result.sin_value = add(result.sin_value, sin_term, precision);

        const std::uint64_t cos_den1 = sin_den1 - 1u;
        const std::uint64_t cos_den2 = sin_den1;
        const std::uint64_t cos_den = checked_taylor_counter_product(cos_den1, cos_den2);
        cos_term = div(mul(cos_term, x2, precision), make_u64(cos_den, precision), precision);
        cos_term = negate(cos_term, precision);
        result.cos_value = add(result.cos_value, cos_term, precision);

        if (abs_prec(sin_term, precision) < epsilon && abs_prec(cos_term, precision) < epsilon) {
            return result;
        }
        increment_taylor_counter(k);
    }
    throw_iteration_limit_exceeded("sincos_taylor_small");
}

inline void refresh_trig_constants_locked(trig_constant_cache_state &cache, precision_type cache_precision) {
    cache_precision = normalize_target_precision(cache_precision);
    if (!cache.initialized || cache.cached_precision < cache_precision) {
        cache.pi_value = compute_pi_gauss_legendre(cache_precision);
        cache.pi_over_two_value = set_prec_copy(cache.pi_value, cache_precision);
        mpf_div_2exp(cache.pi_over_two_value.get_mpf_t(), cache.pi_over_two_value.get_mpf_t(), 1);
        cache.two_over_pi_value = div(make_ui(2, cache_precision), cache.pi_value, cache_precision);
        cache.cached_precision = cache_precision;
        cache.initialized = true;
    }
}

inline trig_constants_result trig_constants(precision_type cache_precision) {
    trig_constant_cache_state &cache = trig_constant_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    refresh_trig_constants_locked(cache, cache_precision);
    trig_constants_result result;
    result.pi_over_two_value = set_prec_copy(cache.pi_over_two_value, cache_precision);
    result.two_over_pi_value = set_prec_copy(cache.two_over_pi_value, cache_precision);
    return result;
}

inline sincos_result compute_sincos(const mpf_class &x_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_trig(target);
    const precision_type const_precision = trig_constant_precision_for_argument(x_input, target);
    const trig_constants_result constants = trig_constants(const_precision);
    const mpf_class pio2 = constants.pi_over_two_value;
    const mpf_class two_over_pi = constants.two_over_pi_value;

    const mpf_class scaled_x = set_prec_copy(x_input, const_precision);
    const mpf_class q = mul(scaled_x, two_over_pi, const_precision);
    mpz_class k = 0;
    mpz_set_f(k.get_mpz_t(), q.get_mpf_t());

    mpf_class integer_part(0, const_precision);
    mpf_set_z(integer_part.get_mpf_t(), k.get_mpz_t());
    const mpf_class frac = sub(q, integer_part, const_precision);
    const mpf_class half = div_ui(make_ui(1, const_precision), 2ul, const_precision);
    const mpf_class neg_half = negate(half, const_precision);
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

    mpf_class k_mpf(0, const_precision);
    mpf_set_z(k_mpf.get_mpf_t(), k.get_mpz_t());
    const mpf_class remainder_hi = sub(scaled_x, mul(pio2, k_mpf, const_precision), const_precision);
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
        result.cos_value = set_prec_copy(negate(base.sin_value, work), work);
        break;
    case 2:
        result.sin_value = set_prec_copy(negate(base.sin_value, work), work);
        result.cos_value = set_prec_copy(negate(base.cos_value, work), work);
        break;
    default:
        result.sin_value = set_prec_copy(negate(base.cos_value, work), work);
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

    std::uint64_t k = 1;
    const std::uint64_t max_iterations = iteration_limit_for_precision(precision);
    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        power = mul(power, x2, precision);
        mpf_class contribution = div(power, make_u64(checked_taylor_odd_denominator(k), precision), precision);
        if ((k & std::uint64_t{1}) == std::uint64_t{1}) {
            contribution = negate(contribution, precision);
        }
        sum = add(sum, contribution, precision);
        if (abs_prec(contribution, precision) < epsilon) {
            return sum;
        }
        increment_taylor_counter(k);
    }
    throw_iteration_limit_exceeded("atan_taylor_small");
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

    const bool negative = x < zero;
    const mpf_class ax = negative ? negate(x, work) : x;

    mpf_class y = set_prec_copy(ax, work);
    std::uint64_t reductions = 0;
    const std::uint64_t max_reductions = iteration_limit_for_precision(work);
    const mpf_class threshold = div_ui(one, 2ul, work);
    while (y > threshold) {
        if (reductions >= max_reductions) {
            throw_iteration_limit_exceeded("compute_atan argument reduction");
        }
        const mpf_class sqrt_term = sqrt_prec(add(one, sqr(y, work), work), work);
        y = div(y, add(one, sqrt_term, work), work);
        ++reductions;
    }

    mpf_class result = atan_taylor_small(y, work);
    for (std::uint64_t i = 0; i < reductions; ++i) {
        result = mul_ui(result, 2ul, work);
    }

    if (negative) {
        result = negate(result, work);
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
        return set_prec_copy(negate(pio2, target + 2), target);
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

inline mpf_class e(precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_exp(target) + 8;
    return set_prec_copy(compute_exp(make_ui(1, work), work), target);
}

struct log_ten_cache_state {
    std::mutex mutex;
    precision_type cached_precision = 0;
    mpf_class cached_value = make_ui(0, 32);
    bool initialized = false;
};

inline log_ten_cache_state &log_ten_cache() {
    static log_ten_cache_state cache;
    return cache;
}

inline mpf_class log_ten(precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = target + guard_bits_for_log(target) + 8;

    log_ten_cache_state &cache = log_ten_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    if (!cache.initialized || cache.cached_precision < work) {
        cache.cached_value = compute_log(make_ui(10, work), work);
        cache.cached_precision = work;
        cache.initialized = true;
    }
    return set_prec_copy(cache.cached_value, target);
}

inline mpf_class inv_log_two(precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_log_two(target) + 8;
    return set_prec_copy(div(make_ui(1, work), log_two(work), work), target);
}

inline mpf_class pi_over_two(precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    mpf_class result = pi(target + 8);
    mpf_div_2exp(result.get_mpf_t(), result.get_mpf_t(), 1);
    return set_prec_copy(result, target);
}

inline mpf_class pi_over_four(precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    mpf_class result = pi(target + 8);
    mpf_div_2exp(result.get_mpf_t(), result.get_mpf_t(), 2);
    return set_prec_copy(result, target);
}

inline mpf_class two_pi(precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    mpf_class result = pi(target + 8);
    mpf_mul_2exp(result.get_mpf_t(), result.get_mpf_t(), 1);
    return set_prec_copy(result, target);
}

inline precision_type guard_bits_for_pow(precision_type) {
    return 96;
}

inline precision_type working_precision_for_pow(precision_type target_precision) {
    return normalize_target_precision(target_precision) + guard_bits_for_pow(target_precision);
}

inline bool mpf_is_exact_integer(const mpf_class &x, mpz_class &integer_value) {
    mpz_set_f(integer_value.get_mpz_t(), x.get_mpf_t());
    mpf_class rounded(0, x.get_prec());
    mpf_set_z(rounded.get_mpf_t(), integer_value.get_mpz_t());
    return rounded == x;
}

inline mpf_class pow_integer_unsigned(const mpf_class &base_input, const mpz_class &exponent, precision_type precision) {
    mpf_class result = make_ui(1, precision);
    mpf_class base = set_prec_copy(base_input, precision);
    mpz_class e = exponent;

    while (e > 0) {
        if (mpz_odd_p(e.get_mpz_t())) {
            ensure_mpf_mul_exponent_fits(result, base);
            result = mul(result, base, precision);
        }
        e >>= 1;
        if (e > 0) {
            ensure_mpf_mul_exponent_fits(base, base);
            base = sqr(base, precision);
        }
    }
    return result;
}

inline mpf_class compute_pow(const mpf_class &x_input, const mpf_class &y_input, precision_type target_precision) {
    const precision_type target = normalize_target_precision(target_precision);
    const precision_type work = working_precision_for_pow(target);
    const mpf_class x = set_prec_copy(x_input, work);
    const mpf_class y = set_prec_copy(y_input, work);
    const mpf_class zero = make_ui(0, work);
    const mpf_class one = make_ui(1, work);

    if (y == zero) {
        if (x == zero) {
            throw std::domain_error("pow(0, 0) is undefined");
        }
        return make_ui(1, target);
    }

    if (x == zero) {
        if (y > zero) {
            return make_ui(0, target);
        }
        throw std::domain_error("pow(0, y) is undefined for y <= 0");
    }

    mpz_class integer_exponent;
    if (mpf_is_exact_integer(y, integer_exponent)) {
        const bool negative_exponent = integer_exponent < 0;
        if (negative_exponent) {
            integer_exponent = -integer_exponent;
        }
        mpf_class magnitude = pow_integer_unsigned(x, integer_exponent, work);
        if (negative_exponent) {
            magnitude = div(one, magnitude, work);
        }
        return set_prec_copy(magnitude, target);
    }

    if (x < zero) {
        throw std::domain_error("pow(x, y) is undefined for x < 0 and non-integer y");
    }

    const mpf_class exponent_product = mul(y, compute_log(x, work), work);
    return set_prec_copy(compute_exp(exponent_product, work), target);
}

} // namespace mplapack_gmp_transcendents

#endif
