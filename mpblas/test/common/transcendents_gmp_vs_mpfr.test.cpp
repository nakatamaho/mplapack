#ifdef ___MPLAPACK_BUILD_WITH_MPFR___

#include <cstdio>
#include <cstdlib>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <mplapack_gmp_transcendents.h>

using mplapack_gmp_transcendents::compute_atan;
using mplapack_gmp_transcendents::compute_atan2;
using mplapack_gmp_transcendents::compute_cos;
using mplapack_gmp_transcendents::compute_exp;
using mplapack_gmp_transcendents::compute_expm1;
using mplapack_gmp_transcendents::compute_log;
using mplapack_gmp_transcendents::compute_log1p;
using mplapack_gmp_transcendents::compute_pow;
using mplapack_gmp_transcendents::compute_sin;
using mplapack_gmp_transcendents::div;
using mplapack_gmp_transcendents::guard_bits_for_expm1;
using mplapack_gmp_transcendents::log_two;
using mplapack_gmp_transcendents::make_ui;
using mplapack_gmp_transcendents::pi;

namespace {

struct mpfr_holder {
    mpfr_t value;
    explicit mpfr_holder(mpfr_prec_t precision) {
        mpfr_init2(value, precision);
    }
    ~mpfr_holder() {
        mpfr_clear(value);
    }
    mpfr_holder(const mpfr_holder &) = delete;
    mpfr_holder &operator=(const mpfr_holder &) = delete;
};

mpf_class parse_decimal_literal(const char *text, mp_bitcnt_t precision) {
    mpf_class value(0, precision);
    if (mpf_set_str(value.get_mpf_t(), text, 10) != 0) {
        std::fprintf(stderr, "failed to parse decimal literal: %s\n", text);
        std::exit(1);
    }
    return value;
}

void set_from_mpf(mpfr_t dst, const mpf_class &src) {
    mpfr_set_f(dst, src.get_mpf_t(), MPFR_RNDN);
}

void set_from_ui(mpfr_t dst, unsigned long value) {
    mpfr_set_ui(dst, value, MPFR_RNDN);
}

double ulp_distance(const mpf_class &gmp_value, mpfr_t mpfr_value, mp_bitcnt_t precision) {
    mpfr_holder gmp_as_mpfr(precision + 64);
    mpfr_holder diff(precision + 64);
    mpfr_holder ulp(precision + 64);

    set_from_mpf(gmp_as_mpfr.value, gmp_value);
    mpfr_sub(diff.value, gmp_as_mpfr.value, mpfr_value, MPFR_RNDN);
    mpfr_abs(diff.value, diff.value, MPFR_RNDN);

    if (mpfr_zero_p(mpfr_value)) {
        mpfr_set_ui(ulp.value, 1, MPFR_RNDN);
        mpfr_div_2si(ulp.value, ulp.value, precision, MPFR_RNDN);
    } else {
        const mpfr_exp_t exponent = mpfr_get_exp(mpfr_value);
        mpfr_set_ui(ulp.value, 1, MPFR_RNDN);
        mpfr_mul_2si(ulp.value, ulp.value, exponent - static_cast<mpfr_exp_t>(precision), MPFR_RNDN);
    }

    mpfr_div(diff.value, diff.value, ulp.value, MPFR_RNDN);
    return mpfr_get_d(diff.value, MPFR_RNDN);
}

void require_ulp(const char *label, const mpf_class &gmp_value, mpfr_t mpfr_value, mp_bitcnt_t precision, double tolerance = 1.0) {
    const double distance = ulp_distance(gmp_value, mpfr_value, precision);
    if (distance > tolerance) {
        std::fprintf(stderr, "%s failed at %lu bits: ulp=%0.17g\n", label, static_cast<unsigned long>(precision), distance);
        std::exit(1);
    }
}

mpf_class make_fraction(unsigned long numer, unsigned long denom, mp_bitcnt_t precision) {
    return div(make_ui(numer, precision), make_ui(denom, precision), precision);
}

mpf_class make_power_of_two_negative(mp_bitcnt_t exponent, mp_bitcnt_t precision) {
    mpf_class value = make_ui(1, precision);
    mpf_div_2exp(value.get_mpf_t(), value.get_mpf_t(), exponent);
    return value;
}

void test_pi_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        mpfr_holder ref(p);
        mpfr_const_pi(ref.value, MPFR_RNDN);
        require_ulp("pi_vs_mpfr", pi(p), ref.value, p);
    }
    std::printf("pi vs mpfr checks passed\n");
}

void test_log2_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        mpfr_holder ref(p);
        mpfr_const_log2(ref.value, MPFR_RNDN);
        require_ulp("log2_vs_mpfr", log_two(p), ref.value, p);
    }
    std::printf("log2 vs mpfr checks passed\n");
}

void test_log1p_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 1, MPFR_RNDN);
            mpfr_div_ui(x.value, x.value, 10, MPFR_RNDN);
            mpfr_log1p(ref.value, x.value, MPFR_RNDN);
            require_ulp("log1p_vs_mpfr_pos_tenth", compute_log1p(make_fraction(1, 10, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_si(x.value, -1, MPFR_RNDN);
            mpfr_div_ui(x.value, x.value, 10, MPFR_RNDN);
            mpfr_log1p(ref.value, x.value, MPFR_RNDN);
            require_ulp("log1p_vs_mpfr_neg_tenth", compute_log1p(-make_fraction(1, 10, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui_2exp(x.value, 1, -static_cast<mpfr_exp_t>(p / 2), MPFR_RNDN);
            mpfr_log1p(ref.value, x.value, MPFR_RNDN);
            require_ulp("log1p_vs_mpfr_small", compute_log1p(make_power_of_two_negative(p / 2, p), p), ref.value, p);
        }
    }
    std::printf("log1p vs mpfr checks passed\n");
}

void test_log_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 25, MPFR_RNDN);
            mpfr_log(ref.value, x.value, MPFR_RNDN);
            require_ulp("log_vs_mpfr_25", compute_log(make_ui(25, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 11, MPFR_RNDN);
            mpfr_div_ui(x.value, x.value, 10, MPFR_RNDN);
            mpfr_log(ref.value, x.value, MPFR_RNDN);
            require_ulp("log_vs_mpfr_near_one", compute_log(make_fraction(11, 10, p), p), ref.value, p, 4.0);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 3, MPFR_RNDN);
            mpfr_div_ui(x.value, x.value, 4, MPFR_RNDN);
            mpfr_log(ref.value, x.value, MPFR_RNDN);
            require_ulp("log_vs_mpfr_three_quarters", compute_log(make_fraction(3, 4, p), p), ref.value, p);
        }
    }
    std::printf("log vs mpfr checks passed\n");
}

void test_exp_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 1, MPFR_RNDN);
            mpfr_exp(ref.value, x.value, MPFR_RNDN);
            require_ulp("exp_vs_mpfr_one", compute_exp(make_ui(1, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_si(x.value, -1, MPFR_RNDN);
            mpfr_exp(ref.value, x.value, MPFR_RNDN);
            require_ulp("exp_vs_mpfr_minus_one", compute_exp(-make_ui(1, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 10, MPFR_RNDN);
            mpfr_exp(ref.value, x.value, MPFR_RNDN);
            require_ulp("exp_vs_mpfr_ten", compute_exp(make_ui(10, p), p), ref.value, p);
        }
    }
    std::printf("exp vs mpfr checks passed\n");
}

void test_expm1_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 1, MPFR_RNDN);
            mpfr_div_ui(x.value, x.value, 10, MPFR_RNDN);
            mpfr_expm1(ref.value, x.value, MPFR_RNDN);
            require_ulp("expm1_vs_mpfr_pos_tenth", compute_expm1(make_fraction(1, 10, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_si(x.value, -1, MPFR_RNDN);
            mpfr_div_ui(x.value, x.value, 10, MPFR_RNDN);
            mpfr_expm1(ref.value, x.value, MPFR_RNDN);
            require_ulp("expm1_vs_mpfr_neg_tenth", compute_expm1(-make_fraction(1, 10, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui_2exp(x.value, 1, -static_cast<mpfr_exp_t>(p / 2), MPFR_RNDN);
            mpfr_expm1(ref.value, x.value, MPFR_RNDN);
            require_ulp("expm1_vs_mpfr_small", compute_expm1(make_power_of_two_negative(p / 2, p), p), ref.value, p);
        }
        {
            const mp_bitcnt_t old_work = p + guard_bits_for_expm1(p);
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui_2exp(x.value, 3, -static_cast<mpfr_exp_t>(old_work / 2), MPFR_RNDN);
            mpfr_expm1(ref.value, x.value, MPFR_RNDN);

            mpf_class gx = make_ui(3, p);
            mpf_div_2exp(gx.get_mpf_t(), gx.get_mpf_t(), old_work / 2);
            require_ulp("expm1_vs_mpfr_cancellation_boundary", compute_expm1(gx, p), ref.value, p);
        }
    }
    std::printf("expm1 vs mpfr checks passed\n");
}

void test_sin_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 1, MPFR_RNDN);
            mpfr_sin(ref.value, x.value, MPFR_RNDN);
            require_ulp("sin_vs_mpfr_one", compute_sin(make_ui(1, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            const mpf_class gx = div(pi(p), make_ui(6, p), p);
            set_from_mpf(x.value, gx);
            mpfr_sin(ref.value, x.value, MPFR_RNDN);
            require_ulp("sin_vs_mpfr_pi_over_six", compute_sin(gx, p), ref.value, p, 2.0);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 100, MPFR_RNDN);
            mpfr_sin(ref.value, x.value, MPFR_RNDN);
            require_ulp("sin_vs_mpfr_hundred", compute_sin(make_ui(100, p), p), ref.value, p);
        }
    }
    std::printf("sin vs mpfr checks passed\n");
}

void test_cos_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 1, MPFR_RNDN);
            mpfr_cos(ref.value, x.value, MPFR_RNDN);
            require_ulp("cos_vs_mpfr_one", compute_cos(make_ui(1, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            const mpf_class gx = div(pi(p), make_ui(3, p), p);
            set_from_mpf(x.value, gx);
            mpfr_cos(ref.value, x.value, MPFR_RNDN);
            require_ulp("cos_vs_mpfr_pi_over_three", compute_cos(gx, p), ref.value, p, 2.0);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 100, MPFR_RNDN);
            mpfr_cos(ref.value, x.value, MPFR_RNDN);
            require_ulp("cos_vs_mpfr_hundred", compute_cos(make_ui(100, p), p), ref.value, p);
        }
    }
    std::printf("cos vs mpfr checks passed\n");
}

void test_atan_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 1, MPFR_RNDN);
            mpfr_atan(ref.value, x.value, MPFR_RNDN);
            require_ulp("atan_vs_mpfr_one", compute_atan(make_ui(1, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 10, MPFR_RNDN);
            mpfr_atan(ref.value, x.value, MPFR_RNDN);
            require_ulp("atan_vs_mpfr_ten", compute_atan(make_ui(10, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 1, MPFR_RNDN);
            mpfr_div_ui(x.value, x.value, 10, MPFR_RNDN);
            mpfr_atan(ref.value, x.value, MPFR_RNDN);
            require_ulp("atan_vs_mpfr_tenth", compute_atan(make_fraction(1, 10, p), p), ref.value, p);
        }
    }
    std::printf("atan vs mpfr checks passed\n");
}

void test_atan2_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        {
            mpfr_holder y(p);
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(y.value, 1, MPFR_RNDN);
            mpfr_set_ui(x.value, 1, MPFR_RNDN);
            mpfr_atan2(ref.value, y.value, x.value, MPFR_RNDN);
            require_ulp("atan2_vs_mpfr_q1", compute_atan2(make_ui(1, p), make_ui(1, p), p), ref.value, p);
        }
        {
            mpfr_holder y(p);
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_ui(y.value, 1, MPFR_RNDN);
            mpfr_set_si(x.value, -1, MPFR_RNDN);
            mpfr_atan2(ref.value, y.value, x.value, MPFR_RNDN);
            require_ulp("atan2_vs_mpfr_q2", compute_atan2(make_ui(1, p), -make_ui(1, p), p), ref.value, p);
        }
        {
            mpfr_holder y(p);
            mpfr_holder x(p);
            mpfr_holder ref(p);
            mpfr_set_si(y.value, -1, MPFR_RNDN);
            mpfr_set_ui(x.value, 1, MPFR_RNDN);
            mpfr_atan2(ref.value, y.value, x.value, MPFR_RNDN);
            require_ulp("atan2_vs_mpfr_q4", compute_atan2(-make_ui(1, p), make_ui(1, p), p), ref.value, p);
        }
    }
    std::printf("atan2 vs mpfr checks passed\n");
}

void test_pow_vs_mpfr(const std::vector<mp_bitcnt_t> &precisions) {
    for (auto p : precisions) {
        {
            mpfr_holder x(p);
            mpfr_holder y(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 2, MPFR_RNDN);
            mpfr_set_ui(y.value, 1, MPFR_RNDN);
            mpfr_div_ui(y.value, y.value, 2, MPFR_RNDN);
            mpfr_pow(ref.value, x.value, y.value, MPFR_RNDN);
            require_ulp("pow_vs_mpfr_sqrt2", compute_pow(make_ui(2, p), make_fraction(1, 2, p), p), ref.value, p);
        }
        {
            mpfr_holder x(p);
            mpfr_holder y(p);
            mpfr_holder ref(p);
            mpfr_set_ui(x.value, 3, MPFR_RNDN);
            mpfr_set_ui(y.value, 5, MPFR_RNDN);
            mpfr_div_ui(y.value, y.value, 4, MPFR_RNDN);
            mpfr_pow(ref.value, x.value, y.value, MPFR_RNDN);
            require_ulp("pow_vs_mpfr_three_five_fourths", compute_pow(make_ui(3, p), make_fraction(5, 4, p), p), ref.value, p);
        }
    }
    std::printf("pow vs mpfr checks passed\n");
}

} // namespace

int main() {
    const std::vector<mp_bitcnt_t> precisions = {64, 128, 256, 512, 1024, 2048, 4096};
    std::printf("*** Testing GMP transcendents vs MPFR start ***\n");
    test_pi_vs_mpfr(precisions);
    test_log2_vs_mpfr(precisions);
    test_log1p_vs_mpfr(precisions);
    test_log_vs_mpfr(precisions);
    test_exp_vs_mpfr(precisions);
    test_expm1_vs_mpfr(precisions);
    test_sin_vs_mpfr(precisions);
    test_cos_vs_mpfr(precisions);
    test_atan_vs_mpfr(precisions);
    test_atan2_vs_mpfr(precisions);
    test_pow_vs_mpfr(precisions);
    std::printf("*** Testing GMP transcendents vs MPFR successful ***\n");
    return 0;
}

#else
#error "transcendents_gmp_vs_mpfr.test.cpp requires ___MPLAPACK_BUILD_WITH_MPFR___"
#endif
