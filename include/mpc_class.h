/*
 * Copyright (c) 2008-2022
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mpc_class.h,v 1.17 2010/08/07 03:15:46 nakatamaho Exp $
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
/****************************************************************
Copyright 1990 - 1997 by AT&T, Lucent Technologies and Bellcore.

Permission to use, copy, modify, and distribute this software
and its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the names of AT&T, Bell Laboratories,
Lucent or Bellcore or any of their entities not be used in
advertising or publicity pertaining to distribution of the
software without specific, written prior permission.

AT&T, Lucent and Bellcore disclaim all warranties with regard to
this software, including all implied warranties of
merchantability and fitness.  In no event shall AT&T, Lucent or
Bellcore be liable for any special, indirect or consequential
damages or any damages whatsoever resulting from loss of use,
data or profits, whether in an action of contract, negligence or
other tortious action, arising out of or in connection with the
use or performance of this software.
****************************************************************/
/*
Complex class declare for the GMP.
*/

#ifndef _MPC_CLASS_H_
#define _MPC_CLASS_H_

#include "gmpxx.h"
#include <algorithm>
#include <complex>
#include <utility>

class mpc_class {
  private:
    mpf_class re, im;

  public:
    // constructor and deconstructor
    mpc_class();
    mpc_class(const mpc_class &a);
    mpc_class(const mpf_class &a, const mpf_class &b);
    mpc_class(const mpf_class &a);
    mpc_class(const mpf_t a, const mpf_t b);
    mpc_class(const std::complex<double> a);
    mpc_class(const double &a, const double &b);
    mpc_class(const double &a);
    ~mpc_class();

    // extraction of real and imaginary parts
    inline mpf_class &real() { return re; }
    inline const mpf_class &real() const { return re; }
    inline mpf_class &imag() { return im; }
    inline const mpf_class &imag() const { return im; }
    // assignment of real part
    inline void real(const mpf_class r) { this->re = r; }
    inline void imag(const mpf_class r) { this->im = r; }

    // comparison
    friend bool operator==(const mpc_class &a, const mpc_class &b);
    friend bool operator==(const mpc_class &a, const mpf_class &b);
    friend bool operator==(const mpf_class &a, const mpc_class &b);
    friend bool operator==(const mpc_class &a, const double &b);
    friend bool operator==(const double &a, const mpc_class &b);

    friend bool operator!=(const mpc_class &a, const mpc_class &b);
    friend bool operator!=(const mpc_class &a, const mpf_class &b);
    friend bool operator!=(const mpf_class &a, const mpc_class &b);
    friend bool operator!=(const mpc_class &a, const double &b);
    friend bool operator!=(const double &a, const mpc_class &b);

    // subsututtion
    // difficult to implement; mpc_class& operator=(const mpc_class& b);
    mpc_class &operator=(const std::complex<double> b);
    mpc_class &operator=(const mpf_class b);
    mpc_class &operator=(const double b);

    // addition
    mpc_class &operator+=(const mpc_class &b);
    mpc_class &operator+=(const mpf_class &b);
    const mpc_class operator+() const;

    // subtraction
    mpc_class &operator-=(const mpc_class &b);
    mpc_class &operator-=(const mpf_class &b);
    const mpc_class operator-() const;

    // multiplication
    mpc_class &operator*=(const mpc_class &b);
    mpc_class &operator*=(const mpf_class &b);

    // division
    mpc_class &operator/=(const mpc_class &b);
    mpc_class &operator/=(const mpf_class &b);
};

const mpc_class operator+(const mpc_class &a, const mpc_class &b);
const mpc_class operator+(const mpc_class &a, const mpf_class &b);
const mpc_class operator+(const mpf_class &a, const mpc_class &b);
const mpc_class operator+(const mpc_class &a, std::complex<double> b);
const mpc_class operator+(std::complex<double> a, const mpc_class &b);

const mpc_class operator-(const mpc_class &a, const mpc_class &b);
const mpc_class operator-(const mpc_class &a, const mpf_class &b);
const mpc_class operator-(const mpf_class &a, const mpc_class &b);
const mpc_class operator-(const mpc_class &a, std::complex<double> b);
const mpc_class operator-(std::complex<double> a, const mpc_class &b);

const mpc_class operator*(const mpc_class &a, const mpc_class &b);
const mpc_class operator*(const mpc_class &a, const mpf_class &b);
const mpc_class operator*(const mpf_class &a, const mpc_class &b);
const mpc_class operator*(const mpc_class &a, std::complex<double> b);
const mpc_class operator*(std::complex<double> a, const mpc_class &b);

const mpc_class operator/(const mpc_class &a, const mpc_class &b);
const mpc_class operator/(const mpc_class &a, const mpf_class &b);
const mpc_class operator/(const mpf_class &a, const mpc_class &b);
const mpc_class operator/(const mpc_class &a, std::complex<double> b);
const mpc_class operator/(std::complex<double> a, const mpc_class &b);
mpf_class real(const mpc_class &a);
mpf_class imag(const mpc_class &a);
mpf_class norm(const mpc_class &a);
mpf_class abs(const mpc_class &a);
mpc_class sqrt(mpc_class z);
mpc_class conj(const mpc_class &a);
void swap(mpc_class &a, mpc_class &b);

// constructor
inline mpc_class::mpc_class() {
    re = 0.0;
    im = 0.0;
}

inline mpc_class::mpc_class(const mpc_class &a) {
    re = a.re;
    im = a.im;
}

inline mpc_class::mpc_class(const mpf_class &a, const mpf_class &b) {
    re = a;
    im = b;
}

inline mpc_class::mpc_class(const mpf_class &a) {
    re = a;
    im = 0.0;
}

inline mpc_class::mpc_class(const mpf_t a, const mpf_t b) {
    mpf_class tmpr(a), tmpi(b);
    re = tmpr;
    im = tmpi;
}

inline mpc_class::mpc_class(const std::complex<double> a) {
    mpf_class tmpr(a.real()), tmpi(a.imag());
    re = tmpr;
    im = tmpi;
}

inline mpc_class::mpc_class(const double &a, const double &b) {
    mpf_class tmpr(a), tmpi(b);
    re = tmpr;
    im = tmpi;
}

inline mpc_class::mpc_class(const double &a) {
    mpf_class tmp(a);
    re = tmp;
    im = 0.0;
}

inline mpc_class::~mpc_class() {}

// comparison
inline bool operator==(const mpc_class &a, const mpc_class &b) {
    if ((a.re == b.re) && (a.im == b.im)) {
        return true;
    } else
        return false;
}

inline bool operator==(const mpf_class &a, const mpc_class &b) {
    if ((a == b.re) && (b.im == 0.0)) {
        return true;
    } else
        return false;
}

inline bool operator==(const mpc_class &a, const mpf_class &b) {
    if ((a.re == b) && (a.im == 0.0)) {
        return true;
    } else
        return false;
}

inline bool operator==(const mpc_class &a, const double &b) {
    if ((a.re == b) && (a.im == 0.0)) {
        return true;
    } else
        return false;
}

inline bool operator==(const double &a, const mpc_class &b) { return b == a; }

inline bool operator!=(const mpc_class &a, const mpc_class &b) {
    if ((a.re != b.re) || (a.im != b.im)) {
        return true;
    } else
        return false;
}

inline bool operator!=(const mpf_class &a, const mpc_class &b) {
    if ((a != b.re) || (b.im != 0.0)) {
        return true;
    } else
        return false;
}

inline bool operator!=(const mpc_class &a, const mpf_class &b) {
    if ((a.re != b) || (a.im != 0.0)) {
        return true;
    } else
        return false;
}

inline bool operator!=(const mpc_class &a, const double &b) {
    if ((a.re != b) || (a.im != 0.0)) {
        return true;
    } else
        return false;
}

inline bool operator!=(const double &a, const mpc_class &b) { return b != a; }

inline mpc_class &mpc_class::operator=(const std::complex<double> b) {
    re = b.real();
    im = b.imag();
    return *this;
}

inline mpc_class &mpc_class::operator=(const mpf_class b) {
    re = b;
    im = 0.0;
    return *this;
}

inline mpc_class &mpc_class::operator=(const double b) {
    re = b;
    im = 0.0;
    return *this;
}

inline mpc_class &mpc_class::operator+=(const mpc_class &b) {
    re += b.re;
    im += b.im;
    return *this;
}

inline mpc_class &mpc_class::operator+=(const mpf_class &b) {
    re += b;
    return *this;
}

inline const mpc_class mpc_class::operator+() const { return mpc_class(*this); }

inline const mpc_class operator+(const mpc_class &a, const mpc_class &b) {
    mpc_class tmp(a);
    tmp += b;
    return tmp;
}

inline const mpc_class operator+(const mpc_class &a, const mpf_class &b) {
    mpc_class tmp(b);
    tmp += a;
    return tmp;
}

inline const mpc_class operator+(const mpf_class &a, const mpc_class &b) {
    mpc_class tmp(a);
    tmp += b;
    return tmp;
}

inline const mpc_class operator+(const mpc_class &a, std::complex<double> b) {
    return a + mpc_class(b);
}

inline const mpc_class operator+(std::complex<double> a, const mpc_class &b) {
    return mpc_class(a) + b;
}

inline mpc_class &mpc_class::operator-=(const mpc_class &b) {
    re -= b.re;
    im -= b.im;
    return *this;
}

inline mpc_class &mpc_class::operator-=(const mpf_class &b) {
    re -= b;
    return *this;
}

inline const mpc_class mpc_class::operator-() const {
    mpc_class tmp;
    tmp.re = -re;
    tmp.im = -im;
    return tmp;
}

inline const mpc_class operator-(const mpc_class &a, const mpc_class &b) {
    mpc_class tmp(a);
    return tmp -= b;
}

inline const mpc_class operator-(const mpc_class &a, const mpf_class &b) {
    mpc_class tmp(b);
    return -(tmp -= a);
}

inline const mpc_class operator-(const mpf_class &a, const mpc_class &b) {
    mpc_class tmp(a);
    return tmp -= b;
}

inline const mpc_class operator-(const mpc_class &a, std::complex<double> b) {
    return a - mpc_class(b);
}

inline const mpc_class operator-(std::complex<double> a, const mpc_class &b) {
    return mpc_class(a) - b;
}

// not so bad as overflow might not occur with GMP; exponet range is extraordinarly large.
inline mpc_class &mpc_class::operator*=(const mpc_class &b) {
    mpc_class tmp(*this);
    re = tmp.re * b.re - tmp.im * b.im;
    im = tmp.re * b.im + tmp.im * b.re;
    return (*this);
}

inline mpc_class &mpc_class::operator*=(const mpf_class &b) {
    mpc_class tmp(*this);
    re = tmp.re * b;
    im = tmp.im * b;
    return (*this);
}

inline const mpc_class operator*(const mpc_class &a, const mpc_class &b) {
    mpc_class tmp(a);
    tmp *= b;
    return tmp;
}

inline const mpc_class operator*(const mpf_class &a, const mpc_class &b) {
    mpc_class tmp;
    tmp.real(a * b.real());
    tmp.imag(a * b.imag());
    return tmp;
}

inline const mpc_class operator*(const mpc_class &a, const mpf_class &b) {
    mpc_class tmp;

    tmp.real(a.real() * b);
    tmp.imag(a.imag() * b);
    return tmp;
}

inline const mpc_class operator*(const mpc_class &a, std::complex<double> b) {
    return a * mpc_class(b);
}

inline const mpc_class operator*(std::complex<double> a, const mpc_class &b) {
    return mpc_class(a) * b;
}

inline mpc_class &mpc_class::operator/=(const mpc_class &b) {
    mpc_class tmp(*this);
    mpf_class abr, abi, ratio, den;

    if ((abr = b.re) < 0.)
        abr = -abr;
    if ((abi = b.im) < 0.)
        abi = -abi;
    if (abr <= abi) {
        if (abi == 0) {
            if (tmp.im != 0 || tmp.re != 0)
                abi = 1.;
            tmp.im = tmp.re = abi / abr;
            return (*this);
        }
        ratio = b.re / b.im;
        den = b.im * (1.0 + ratio * ratio);
        re = (tmp.re * ratio + tmp.im) / den;
        im = (tmp.im * ratio - tmp.re) / den;
    } else {
        ratio = b.im / b.re;
        den = b.re * (1.0 + ratio * ratio);
        re = (tmp.re + tmp.im * ratio) / den;
        im = (tmp.im - tmp.re * ratio) / den;
    }
    return (*this);
}

inline mpc_class &mpc_class::operator/=(const mpf_class &b) {
    mpc_class tmp(*this);
    re = (tmp.re / b);
    im = (tmp.im / b);
    return (*this);
}

inline const mpc_class operator/(const mpc_class &a, const mpc_class &b) {
    mpc_class tmp(a);
    tmp /= b;
    return tmp;
}

inline const mpc_class operator/(const mpc_class &a, const mpf_class &b) {
    mpc_class tmp;
    tmp.real(a.real() / b);
    tmp.imag(a.imag() / b);
    return tmp;
}

inline const mpc_class operator/(const mpf_class &a, const mpc_class &b) {
    mpc_class tmp;

    tmp.real((a * b.real()) / (b.real() * b.real() + b.imag() * b.imag()));
    tmp.imag(-(a * b.imag()) / (b.real() * b.real() + b.imag() * b.imag()));
    return tmp;
}

inline const mpc_class operator/(const mpc_class &a, std::complex<double> b) {
    return a / mpc_class(b);
}

inline const mpc_class operator/(std::complex<double> a, const mpc_class &b) {
    return mpc_class(a) / b;
}

namespace mplapack_gmp_mpc_detail {

inline mp_bitcnt_t precision_for(const mpc_class &value) {
    return std::max(value.real().get_prec(), value.imag().get_prec());
}

inline mpf_class zero(mp_bitcnt_t precision) {
    return mpf_class(0, precision);
}

inline mpf_class one(mp_bitcnt_t precision) {
    mpf_class result(0, precision);
    mpf_set_ui(result.get_mpf_t(), 1ul);
    return result;
}

inline mpf_class half(mp_bitcnt_t precision) {
    mpf_class result = one(precision);
    mpf_div_2exp(result.get_mpf_t(), result.get_mpf_t(), 1);
    return result;
}

inline mpf_class set_prec_copy(const mpf_class &value, mp_bitcnt_t precision) {
    mpf_class result(0, precision);
    mpf_set(result.get_mpf_t(), value.get_mpf_t());
    return result;
}

inline mpf_class abs_prec(const mpf_class &value, mp_bitcnt_t precision) {
    mpf_class result(0, precision);
    mpf_abs(result.get_mpf_t(), value.get_mpf_t());
    return result;
}

inline mpf_class neg_prec(const mpf_class &value, mp_bitcnt_t precision) {
    mpf_class result(0, precision);
    mpf_neg(result.get_mpf_t(), value.get_mpf_t());
    return result;
}

inline mpf_class add(const mpf_class &lhs, const mpf_class &rhs, mp_bitcnt_t precision) {
    mpf_class result(0, precision);
    mpf_add(result.get_mpf_t(), lhs.get_mpf_t(), rhs.get_mpf_t());
    return result;
}

inline mpf_class sub(const mpf_class &lhs, const mpf_class &rhs, mp_bitcnt_t precision) {
    mpf_class result(0, precision);
    mpf_sub(result.get_mpf_t(), lhs.get_mpf_t(), rhs.get_mpf_t());
    return result;
}

inline mpf_class mul(const mpf_class &lhs, const mpf_class &rhs, mp_bitcnt_t precision) {
    mpf_class result(0, precision);
    mpf_mul(result.get_mpf_t(), lhs.get_mpf_t(), rhs.get_mpf_t());
    return result;
}

inline mpf_class div(const mpf_class &lhs, const mpf_class &rhs, mp_bitcnt_t precision) {
    mpf_class result(0, precision);
    mpf_div(result.get_mpf_t(), lhs.get_mpf_t(), rhs.get_mpf_t());
    return result;
}

inline mpf_class sqr(const mpf_class &value, mp_bitcnt_t precision) {
    return mul(value, value, precision);
}

inline mpf_class sqrt_prec(const mpf_class &value, mp_bitcnt_t precision) {
    mpf_class result(0, precision);
    mpf_sqrt(result.get_mpf_t(), value.get_mpf_t());
    return result;
}

inline mpf_class scaled_hypot_abs(const mpf_class &lhs, const mpf_class &rhs, mp_bitcnt_t precision) {
    mpf_class lhs_abs = abs_prec(lhs, precision);
    mpf_class rhs_abs = abs_prec(rhs, precision);

    const mpf_class *high = &lhs_abs;
    const mpf_class *low = &rhs_abs;
    if (mpf_cmp(lhs_abs.get_mpf_t(), rhs_abs.get_mpf_t()) < 0) {
        high = &rhs_abs;
        low = &lhs_abs;
    }

    if (mpf_sgn(high->get_mpf_t()) == 0) {
        return zero(precision);
    }

    mpf_class ratio = div(*low, *high, precision);
    mpf_mul(ratio.get_mpf_t(), ratio.get_mpf_t(), ratio.get_mpf_t());
    mpf_add_ui(ratio.get_mpf_t(), ratio.get_mpf_t(), 1ul);
    mpf_sqrt(ratio.get_mpf_t(), ratio.get_mpf_t());

    mpf_class result(0, precision);
    mpf_mul(result.get_mpf_t(), high->get_mpf_t(), ratio.get_mpf_t());
    return result;
}

} // namespace mplapack_gmp_mpc_detail

inline mpf_class real(const mpc_class &a) { return a.real(); }

inline mpf_class imag(const mpc_class &a) { return a.imag(); }

inline mpf_class norm(const mpc_class &a) {
    const mp_bitcnt_t precision = mplapack_gmp_mpc_detail::precision_for(a);
    mpf_class real_abs = mplapack_gmp_mpc_detail::abs_prec(a.real(), precision);
    mpf_class imag_abs = mplapack_gmp_mpc_detail::abs_prec(a.imag(), precision);

    const mpf_class *high = &real_abs;
    const mpf_class *low = &imag_abs;
    if (mpf_cmp(real_abs.get_mpf_t(), imag_abs.get_mpf_t()) < 0) {
        high = &imag_abs;
        low = &real_abs;
    }

    if (mpf_sgn(high->get_mpf_t()) == 0) {
        return mplapack_gmp_mpc_detail::zero(precision);
    }

    mpf_class ratio = mplapack_gmp_mpc_detail::div(*low, *high, precision);
    mpf_mul(ratio.get_mpf_t(), ratio.get_mpf_t(), ratio.get_mpf_t());
    mpf_add_ui(ratio.get_mpf_t(), ratio.get_mpf_t(), 1ul);

    mpf_class result = mplapack_gmp_mpc_detail::sqr(*high, precision);
    mpf_mul(result.get_mpf_t(), result.get_mpf_t(), ratio.get_mpf_t());
    return result;
}

inline mpf_class abs(const mpc_class &a) {
    const mp_bitcnt_t precision = mplapack_gmp_mpc_detail::precision_for(a);
    return mplapack_gmp_mpc_detail::scaled_hypot_abs(a.real(), a.imag(), precision);
}

inline mpc_class sqrt(mpc_class z) {
    const mp_bitcnt_t precision = mplapack_gmp_mpc_detail::precision_for(z);
    const mpf_class zero = mplapack_gmp_mpc_detail::zero(precision);
    const mpf_class real_part_input = mplapack_gmp_mpc_detail::set_prec_copy(z.real(), precision);
    const mpf_class imag_part_input = mplapack_gmp_mpc_detail::set_prec_copy(z.imag(), precision);

    if (mpf_sgn(imag_part_input.get_mpf_t()) == 0) {
        if (mpf_sgn(real_part_input.get_mpf_t()) < 0) {
            return mpc_class(zero, mplapack_gmp_mpc_detail::sqrt_prec(mplapack_gmp_mpc_detail::neg_prec(real_part_input, precision), precision));
        }
        return mpc_class(mplapack_gmp_mpc_detail::sqrt_prec(real_part_input, precision), zero);
    }

    const mpf_class half = mplapack_gmp_mpc_detail::half(precision);
    const mpf_class magnitude = mplapack_gmp_mpc_detail::scaled_hypot_abs(real_part_input, imag_part_input, precision);

    if (mpf_sgn(real_part_input.get_mpf_t()) >= 0) {
        const mpf_class real_argument = mplapack_gmp_mpc_detail::mul(mplapack_gmp_mpc_detail::add(magnitude, real_part_input, precision), half, precision);
        const mpf_class real_part = mplapack_gmp_mpc_detail::sqrt_prec(real_argument, precision);
        const mpf_class imag_part = mplapack_gmp_mpc_detail::div(imag_part_input, mplapack_gmp_mpc_detail::add(real_part, real_part, precision), precision);
        return mpc_class(real_part, imag_part);
    }

    const mpf_class imag_argument = mplapack_gmp_mpc_detail::mul(mplapack_gmp_mpc_detail::sub(magnitude, real_part_input, precision), half, precision);
    mpf_class imag_magnitude = mplapack_gmp_mpc_detail::sqrt_prec(imag_argument, precision);
    const mpf_class real_part = mplapack_gmp_mpc_detail::div(mplapack_gmp_mpc_detail::abs_prec(imag_part_input, precision), mplapack_gmp_mpc_detail::add(imag_magnitude, imag_magnitude, precision), precision);
    if (mpf_sgn(imag_part_input.get_mpf_t()) < 0) {
        imag_magnitude = mplapack_gmp_mpc_detail::neg_prec(imag_magnitude, precision);
    }
    return mpc_class(real_part, imag_magnitude);
}

inline mpc_class conj(const mpc_class &ctmp) {
    mpc_class ctmp2;

    ctmp2.real(ctmp.real());
    ctmp2.imag(-ctmp.imag());
    return ctmp2;
}

inline void swap(mpc_class &a, mpc_class &b) {
    using std::swap;
    swap(a.real(), b.real());
    swap(a.imag(), b.imag());
}

#endif
