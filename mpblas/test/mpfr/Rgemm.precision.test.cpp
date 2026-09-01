/*
 * Copyright (c) 2026
 *      Nakata, Maho
 *      All rights reserved.
 *
 * SPDX-License-Identifier: BSD-2-Clause
 */

#include <algorithm>
#include <array>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <mpblas_mpfr.h>

namespace {

class DefaultPrecisionGuard {
  public:
    explicit DefaultPrecisionGuard(mpfr_prec_t precision) : saved_(mpfrxx::default_precision_bits()) {
        mpfrxx::set_default_precision_bits(precision);
        require(mpfrxx::default_precision_bits() == precision, "failed to set the test default precision");
    }

    DefaultPrecisionGuard(const DefaultPrecisionGuard &) = delete;
    DefaultPrecisionGuard &operator=(const DefaultPrecisionGuard &) = delete;

    ~DefaultPrecisionGuard() { mpfrxx::set_default_precision_bits(saved_); }

    static void require(bool condition, const char *message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

  private:
    mpfr_prec_t saved_;
};

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

struct Matrix {
    mplapackint rows;
    mplapackint columns;
    mplapackint leading_dimension;
    mpfr_prec_t precision;
    std::vector<mpfr_class> values;

    Matrix(mplapackint row_count, mplapackint column_count, mpfr_prec_t precision_bits) : rows(row_count), columns(column_count), leading_dimension(std::max<mplapackint>(1, row_count)), precision(precision_bits) {
        values.reserve(static_cast<std::size_t>(leading_dimension * columns));
        for (mplapackint column = 0; column < columns; ++column) {
            for (mplapackint row = 0; row < leading_dimension; ++row) {
                values.emplace_back(mpfr_class::with_precision(precision));
                mpfr_set_zero(values.back().mpfr_data(), 1);
            }
        }
    }

    mpfr_class &at(mplapackint row, mplapackint column) { return values[static_cast<std::size_t>(row + column * leading_dimension)]; }

    const mpfr_class &at(mplapackint row, mplapackint column) const { return values[static_cast<std::size_t>(row + column * leading_dimension)]; }
};

bool no_transpose(const char *trans) { return trans[0] == 'N' || trans[0] == 'n'; }

Matrix make_operand(const char *trans, mplapackint op_rows, mplapackint op_columns, mpfr_prec_t precision) { return no_transpose(trans) ? Matrix(op_rows, op_columns, precision) : Matrix(op_columns, op_rows, precision); }

mpfr_class &op_at(Matrix &matrix, const char *trans, mplapackint row, mplapackint column) { return no_transpose(trans) ? matrix.at(row, column) : matrix.at(column, row); }

const mpfr_class &op_at(const Matrix &matrix, const char *trans, mplapackint row, mplapackint column) { return no_transpose(trans) ? matrix.at(row, column) : matrix.at(column, row); }

void set_unsigned(mpfr_class &value, unsigned long integer) { mpfr_set_ui(value.mpfr_data(), integer, MPFR_RNDN); }

void set_power_of_two(mpfr_class &value, long exponent) { mpfr_set_ui_2exp(value.mpfr_data(), 1, exponent, MPFR_RNDN); }

void add_power_of_two(mpfr_class &value, long exponent) {
    mpfr_class increment = mpfr_class::with_precision(value.precision());
    set_power_of_two(increment, exponent);
    mpfr_add(value.mpfr_data(), value.mpfr_data(), increment.mpfr_data(), MPFR_RNDN);
}

mpfr_class explicit_value(mpfr_prec_t precision, unsigned long integer) {
    mpfr_class value = mpfr_class::with_precision(precision);
    set_unsigned(value, integer);
    return value;
}

void require_equal(const mpfr_class &actual, const mpfr_class &expected, const std::string &message) { require(mpfr_equal_p(actual.mpfr_data(), expected.mpfr_data()) != 0, message); }

void require_matrix_precision(const Matrix &matrix, mpfr_prec_t expected, const std::string &message) {
    for (const auto &value: matrix.values) {
        require(value.precision() == expected, message);
    }
}

void call_rgemm(const char *transa, const char *transb, mplapackint m, mplapackint n, mplapackint k, const mpfr_class &alpha, Matrix &a, Matrix &b, const mpfr_class &beta, Matrix &c) { Rgemm(transa, transb, m, n, k, alpha, a.values.data(), a.leading_dimension, b.values.data(), b.leading_dimension, beta, c.values.data(), c.leading_dimension); }

void test_identity_case(mpfr_prec_t operation_precision, long low_bit, mpfr_prec_t external_default, const char *transa, const char *transb) {
    DefaultPrecisionGuard default_guard(external_default);
    Matrix identity_a = make_operand(transa, 2, 2, operation_precision);
    Matrix identity_b = make_operand(transb, 2, 2, operation_precision);
    Matrix significant_a = make_operand(transa, 2, 2, operation_precision);
    Matrix significant_b = make_operand(transb, 2, 2, operation_precision);
    Matrix left_result(2, 2, operation_precision);
    Matrix right_result(2, 2, operation_precision);
    mpfr_class alpha = explicit_value(operation_precision, 1);
    mpfr_class beta = explicit_value(operation_precision, 0);

    set_unsigned(op_at(identity_a, transa, 0, 0), 1);
    set_unsigned(op_at(identity_a, transa, 1, 1), 1);
    set_unsigned(op_at(identity_b, transb, 0, 0), 1);
    set_unsigned(op_at(identity_b, transb, 1, 1), 1);
    set_unsigned(op_at(significant_a, transa, 0, 0), 1);
    set_unsigned(op_at(significant_a, transa, 1, 1), 1);
    set_unsigned(op_at(significant_b, transb, 0, 0), 1);
    set_unsigned(op_at(significant_b, transb, 1, 1), 1);
    add_power_of_two(op_at(significant_a, transa, 0, 0), -low_bit);
    add_power_of_two(op_at(significant_b, transb, 0, 0), -low_bit);

    call_rgemm(transa, transb, 2, 2, 2, alpha, identity_a, significant_b, beta, left_result);
    call_rgemm(transa, transb, 2, 2, 2, alpha, significant_a, identity_b, beta, right_result);

    for (mplapackint column = 0; column < 2; ++column) {
        for (mplapackint row = 0; row < 2; ++row) {
            require_equal(left_result.at(row, column), op_at(significant_b, transb, row, column), "left identity lost MPFR significance");
            require_equal(right_result.at(row, column), op_at(significant_a, transa, row, column), "right identity lost MPFR significance");
        }
    }
    require_matrix_precision(left_result, operation_precision, "left identity changed result precision");
    require_matrix_precision(right_result, operation_precision, "right identity changed result precision");
    require(mpfrxx::default_precision_bits() == external_default, "Rgemm changed the caller's default precision");
}

void test_precision_invariance_matrix() {
    const std::array<std::pair<mpfr_prec_t, long>, 5> operations{{{32, 20}, {128, 90}, {512, 400}, {1024, 700}, {2048, 1500}}};
    const std::array<std::vector<mpfr_prec_t>, 5> defaults{{{128, 2048, 4096}, {64, 512}, {128, 2048}, {64, 512, 2048}, {128, 4096}}};
    for (std::size_t index = 0; index < operations.size(); ++index) {
        for (mpfr_prec_t external_default: defaults[index]) {
            test_identity_case(operations[index].first, operations[index].second, external_default, "N", "N");
        }
    }
}

void test_transpose_branches() {
    for (const char *transa: {"N", "T"}) {
        for (const char *transb: {"N", "T"}) {
            test_identity_case(1024, 700, 64, transa, transb);
        }
    }
}

void test_general_alpha_beta() {
    constexpr mpfr_prec_t precision = 1024;
    DefaultPrecisionGuard default_guard(64);
    Matrix a(1, 1, precision);
    Matrix b(1, 1, precision);
    Matrix c(1, 1, precision);
    set_unsigned(a.at(0, 0), 1);
    set_unsigned(b.at(0, 0), 1);
    add_power_of_two(b.at(0, 0), -700);
    set_power_of_two(c.at(0, 0), -900);

    mpfr_class alpha = explicit_value(precision, 1);
    mpfr_class beta = explicit_value(precision, 1);
    add_power_of_two(alpha, -800);
    add_power_of_two(beta, -950);

    mpfr_class expected = mpfr_class::with_precision(precision);
    mpfr_class product = mpfr_class::with_precision(precision);
    mpfr_mul(expected.mpfr_data(), beta.mpfr_data(), c.at(0, 0).mpfr_data(), MPFR_RNDN);
    mpfr_mul(product.mpfr_data(), alpha.mpfr_data(), b.at(0, 0).mpfr_data(), MPFR_RNDN);
    mpfr_fma(expected.mpfr_data(), product.mpfr_data(), a.at(0, 0).mpfr_data(), expected.mpfr_data(), MPFR_RNDN);

    call_rgemm("N", "N", 1, 1, 1, alpha, a, b, beta, c);
    require_equal(c.at(0, 0), expected, "general alpha/beta path disagrees with the explicit MPFR reference");
    require(mpfrxx::default_precision_bits() == 64, "general alpha/beta path changed the default precision");
}

void test_alpha_zero_beta_paths() {
    constexpr mpfr_prec_t precision = 1024;
    DefaultPrecisionGuard default_guard(64);
    Matrix a(1, 1, precision);
    Matrix b(1, 1, precision);
    Matrix c(1, 1, precision);
    set_unsigned(a.at(0, 0), 3);
    set_unsigned(b.at(0, 0), 5);
    set_unsigned(c.at(0, 0), 1);
    add_power_of_two(c.at(0, 0), -700);

    mpfr_class alpha = explicit_value(precision, 0);
    mpfr_class beta = explicit_value(precision, 1);
    const mpfr_class original = c.at(0, 0);
    call_rgemm("N", "N", 1, 1, 1, alpha, a, b, beta, c);
    require_equal(c.at(0, 0), original, "alpha=0, beta=1 path lost MPFR significance");

    beta = explicit_value(precision, 1);
    add_power_of_two(beta, -800);
    mpfr_class expected = mpfr_class::with_precision(precision);
    mpfr_mul(expected.mpfr_data(), beta.mpfr_data(), original.mpfr_data(), MPFR_RNDN);
    call_rgemm("N", "N", 1, 1, 1, alpha, a, b, beta, c);
    require_equal(c.at(0, 0), expected, "alpha=0, general beta path disagrees with the explicit MPFR reference");
    require(mpfrxx::default_precision_bits() == 64, "alpha=0 beta paths changed the default precision");
}

void test_multiterm_accumulation() {
    constexpr mpfr_prec_t precision = 1024;
    DefaultPrecisionGuard default_guard(64);
    Matrix a = make_operand("T", 1, 3, precision);
    Matrix b = make_operand("N", 3, 1, precision);
    Matrix c(1, 1, precision);
    mpfr_class alpha = explicit_value(precision, 1);
    mpfr_class beta = explicit_value(precision, 0);

    for (mplapackint column = 0; column < 3; ++column) {
        set_unsigned(op_at(a, "T", 0, column), 1);
    }
    set_unsigned(op_at(b, "N", 0, 0), 1);
    set_power_of_two(op_at(b, "N", 1, 0), -700);
    set_power_of_two(op_at(b, "N", 2, 0), -800);

    mpfr_class expected = explicit_value(precision, 1);
    add_power_of_two(expected, -700);
    add_power_of_two(expected, -800);
    call_rgemm("T", "N", 1, 1, 3, alpha, a, b, beta, c);
    require_equal(c.at(0, 0), expected, "multi-term Rgemm accumulation lost MPFR significance");
}

void test_operation_precision_scan() {
    Matrix a(2, 2, 128);
    Matrix b(2, 2, 256);
    Matrix c(2, 2, 333);
    a.at(1, 1).set_prec(1024);
    mpfr_class alpha = explicit_value(64, 1);
    mpfr_class beta = explicit_value(512, 0);
    const mpfr_prec_t precision = mplapack_mpfr_rgemm_operation_precision(true, true, 2, 2, 2, alpha, a.values.data(), a.leading_dimension, b.values.data(), b.leading_dimension, beta, c.values.data(), c.leading_dimension);
    require(precision == 1024, "operation precision did not include a "
                               "non-leading high-precision element");
}

void test_mixed_precision_direct_call() {
    DefaultPrecisionGuard default_guard(64);
    Matrix identity(2, 2, 128);
    Matrix operand(2, 2, 1024);
    Matrix result(2, 2, 1024);
    set_unsigned(identity.at(0, 0), 1);
    set_unsigned(identity.at(1, 1), 1);
    set_unsigned(operand.at(0, 0), 1);
    add_power_of_two(operand.at(0, 0), -700);
    set_unsigned(operand.at(1, 1), 1);
    mpfr_class alpha = explicit_value(256, 1);
    mpfr_class beta = explicit_value(512, 0);

    call_rgemm("N", "N", 2, 2, 2, alpha, identity, operand, beta, result);
    for (mplapackint column = 0; column < 2; ++column) {
        for (mplapackint row = 0; row < 2; ++row) {
            require_equal(result.at(row, column), operand.at(row, column), "mixed-precision direct Rgemm call lost MPFR significance");
        }
    }
    require_matrix_precision(identity, 128, "mixed-precision call changed A precision");
    require_matrix_precision(operand, 1024, "mixed-precision call changed B precision");
    require_matrix_precision(result, 1024, "mixed-precision call changed C precision");
}

} // namespace

int main() {
    try {
        const mpfr_prec_t original_default = mpfrxx::default_precision_bits();
        test_precision_invariance_matrix();
        test_transpose_branches();
        test_general_alpha_beta();
        test_alpha_zero_beta_paths();
        test_multiterm_accumulation();
        test_operation_precision_scan();
        test_mixed_precision_direct_call();
        require(mpfrxx::default_precision_bits() == original_default, "test did not restore the original default precision");
        std::cout << "Rgemm MPFR operation-precision regression: PASS\n";
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "Rgemm MPFR operation-precision regression: FAIL: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
