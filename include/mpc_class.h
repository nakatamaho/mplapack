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
#include <complex>

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

    friend bool operator!=(const mpc_class &a, const mpc_class &b);
    friend bool operator!=(const mpc_class &a, const mpf_class &b);
    friend bool operator!=(const mpf_class &a, const mpc_class &b);
    friend bool operator!=(const mpc_class &a, const double &b);

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
    if ((a.re != b) || (a.im == 0.0)) {
        return true;
    } else
        return false;
}

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
    re = (tmp.re * b);
    im = (tmp.im * b);
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

inline mpf_class abs(mpc_class ctemp) {
    mpf_class temp;
    if (ctemp.real() < 0)
        ctemp.real() = -ctemp.real();
    if (ctemp.imag() < 0)
        ctemp.imag() = -ctemp.imag();
    if (ctemp.imag() > ctemp.real()) {
        temp = ctemp.real();
        ctemp.real() = ctemp.imag();
        ctemp.imag() = temp;
    }
    if ((ctemp.real() + ctemp.imag()) == ctemp.real())
        return (ctemp.real());

    temp = ctemp.imag() / ctemp.real();
    temp = ctemp.real() * sqrt(1.0 + temp * temp); /*overflow!!*/

    return temp;
}

inline mpc_class sqrt(mpc_class z) {
    mpf_class mag;
    mpc_class r;

    mag = abs(z);
    if (abs(mag) == 0.0) {
        r.real(0.0), r.imag(0.0);
    } else if (z.real() > 0.0) {
        r.real(sqrt(0.5 * (mag + z.real())));
        r.imag(z.imag() / (2.0 * r.real()));
    } else {
        r.imag(sqrt(0.5 * (mag - z.real())));
        if (z.imag() < 0.0)
            r.imag(-r.imag());
        r.real(z.imag() / (2.0 * r.imag()));
    }
    return r;
}

inline mpc_class conj(mpc_class ctmp) {
    mpc_class ctmp2;

    ctmp2.real(ctmp.real());
    ctmp2.imag(-ctmp.imag());
    return ctmp2;
}

#endif
