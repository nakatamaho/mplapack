/*
 * Copyright (c) 2008-2010
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
/*
Complex class declare for the GMP.
*/

#ifndef _MPC_CLASS_H_
#define _MPC_CLASS_H_

#include <complex>
#include "gmpxx.h"

class mpc_class {
  private:
    mpf_class re, im;

  public:
//constructor
    mpc_class();
    mpc_class(const mpc_class & a);
    mpc_class(const mpf_class & a, const mpf_class & b);
    mpc_class(const mpf_class & a);
    mpc_class(const mpf_t a, const mpf_t b);

//extraction of real and imaginary parts
    const mpf_class & real() const
    {
	return re;
    }
    const mpf_class & imag() const
    {
	return im;
    }
    void real(const mpf_class r)
    {
	this->re = r;
    }
    void imag(const mpf_class r)
    {
	this->im = r;
    }

//comparison
    friend bool operator==(const mpc_class & a, const mpc_class & b);
    friend bool operator==(const mpc_class & a, const mpf_class & b);
    friend bool operator==(const mpf_class & a, const mpc_class & b);

    friend bool operator!=(const mpc_class & a, const mpc_class & b);
    friend bool operator!=(const mpc_class & a, const mpf_class & b);
    friend bool operator!=(const mpf_class & a, const mpc_class & b);

//subsututtion
//difficult to implement; mpc_class& operator=(const mpc_class& b);
    mpc_class & operator=(const std::complex < double >b);
    mpc_class & operator=(const mpf_class b);
    mpc_class & operator=(const double b);

//addition
    mpc_class & operator+=(const mpc_class & b);
    mpc_class & operator+=(const mpf_class & b);
    const mpc_class operator+() const;

//subtraction
    mpc_class & operator-=(const mpc_class & b);
    mpc_class & operator-=(const mpf_class & b);
    const mpc_class operator-() const;

//multiplication
    mpc_class & operator*=(const mpc_class & b);
    mpc_class & operator*=(const mpf_class & b);

//division
    mpc_class & operator/=(const mpc_class & b);
    mpc_class & operator/=(const mpf_class & b);
};

const mpc_class operator+(const mpc_class & a, const mpc_class & b);
const mpc_class operator+(const mpc_class & a, const mpf_class & b);
const mpc_class operator+(const mpf_class & a, const mpc_class & b);
const mpc_class operator+(const mpc_class & a, std::complex < double >&b);
const mpc_class operator+(std::complex < double >&a, const mpc_class & b);

const mpc_class operator-(const mpc_class & a, const mpc_class & b);
const mpc_class operator-(const mpc_class & a, const mpf_class & b);
const mpc_class operator-(const mpf_class & a, const mpc_class & b);
const mpc_class operator-(const mpc_class & a, std::complex < double >&b);
const mpc_class operator-(std::complex < double >&a, const mpc_class & b);

const mpc_class operator*(const mpc_class & a, const mpc_class & b);
const mpc_class operator*(const mpc_class & a, const mpf_class & b);
const mpc_class operator*(const mpf_class & a, const mpc_class & b);
const mpc_class operator*(const mpc_class & a, std::complex < double >&b);
const mpc_class operator*(std::complex < double >&a, const mpc_class & b);

const mpc_class operator/(const mpc_class & a, const mpc_class & b);
const mpc_class operator/(const mpc_class & a, const mpf_class & b);
const mpc_class operator/(const mpf_class & a, const mpc_class & b);
const mpc_class operator/(const mpc_class & a, std::complex < double >&b);
const mpc_class operator/(std::complex < double >&a, const mpc_class & b);

//constructor
inline mpc_class::mpc_class()
{
    re = 0.0;
    im = 0.0;
}

inline mpc_class::mpc_class(const mpc_class & a)
{
    re = a.re;
    im = a.im;
}

inline mpc_class::mpc_class(const mpf_class & a, const mpf_class & b)
{
    re = a;
    im = b;
}

inline mpc_class::mpc_class(const mpf_class & a)
{
    re = a;
    im = 0.0;
}

inline mpc_class::mpc_class(const mpf_t a, const mpf_t b)
{
   mpf_class tmpr(a), tmpi(b);
   re = tmpr; im = tmpi;
}

//comparison
inline bool operator==(const mpc_class & a, const mpc_class & b)
{
    if ((a.re == b.re) && (a.im == b.im)) {
	return true;
    } else
	return false;
}

inline bool operator==(const mpf_class & a, const mpc_class & b)
{
    if ((a == b.re) && (b.im == 0.0)) {
	return true;
    } else
	return false;
}

inline bool operator==(const mpc_class & a, const mpf_class & b)
{
    if ((a.re == b) && (a.im == 0.0)) {
	return true;
    } else
	return false;
}

inline bool operator!=(const mpc_class & a, const mpc_class & b)
{
    if ((a.re != b.re) || (a.im != b.im)) {
	return true;
    } else
	return false;
}

inline bool operator!=(const mpf_class & a, const mpc_class & b)
{
    if ((a != b.re) || (b.im != 0.0)) {
	return true;
    } else
	return false;
}

inline bool operator!=(const mpc_class & a, const mpf_class & b)
{
    if ((a.re != b) || (a.im != 0.0)) {
	return true;
    } else
	return false;
}

inline mpc_class & mpc_class::operator=(const std::complex < double >b)
{
    re = b.real();
    im = b.imag();
    return *this;
}

inline mpc_class & mpc_class::operator=(const mpf_class b)
{
    re = b;
    im = 0.0;
    return *this;
}

inline mpc_class & mpc_class::operator=(const double b)
{
    re = b;
    im = 0.0;
    return *this;
}


inline mpc_class & mpc_class::operator+=(const mpc_class & b)
{
    re += b.re;
    im += b.im;
    return *this;
}

inline mpc_class & mpc_class::operator+=(const mpf_class & b)
{
    re += b;
    return *this;
}

inline const mpc_class mpc_class::operator+() const
{
    return mpc_class(*this);
}

inline const mpc_class operator+(const mpc_class & a, const mpc_class & b)
{
    mpc_class tmp(a);
    tmp += b;
    return tmp;
}

inline const mpc_class operator+(const mpc_class & a, const mpf_class & b)
{
    mpc_class tmp(b);
    tmp += a;
    return tmp;
}

inline const mpc_class operator+(const mpf_class & a, const mpc_class & b)
{
    mpc_class tmp(a);
    tmp += b;
    return tmp;
}

inline mpc_class & mpc_class::operator-=(const mpc_class & b)
{
    re -= b.re;
    im -= b.im;
    return *this;
}

inline mpc_class & mpc_class::operator-=(const mpf_class & b)
{
    re -= b;
    return *this;
}

inline const mpc_class mpc_class::operator-() const
{
    mpc_class tmp;
    tmp.re = -re;
    tmp.im = -im;
    return tmp;
}

inline const mpc_class operator-(const mpc_class & a, const mpc_class & b)
{
    mpc_class tmp(a);
    return tmp -= b;
}

inline const mpc_class operator-(const mpc_class & a, const mpf_class & b)
{
    mpc_class tmp(b);
    return -(tmp -= a);
}

inline const mpc_class operator-(const mpf_class & a, const mpc_class & b)
{
    mpc_class tmp(a);
    return tmp -= b;
}

//not so bad as overflow might not occur with GMP; exponet range is extraordinarly large.
inline mpc_class & mpc_class::operator*=(const mpc_class & b)
{
    mpc_class tmp(*this);
    re = tmp.re * b.re - tmp.im * b.im;
    im = tmp.re * b.im + tmp.im * b.re;
    return (*this);
}

inline mpc_class & mpc_class::operator*=(const mpf_class & b)
{
    mpc_class tmp(*this);
    re = tmp.re * b;
    im = tmp.im * b;
    return (*this);
}

inline const mpc_class operator*(const mpc_class & a, const mpc_class & b)
{
    mpc_class tmp(a);
    tmp *= b;
    return tmp;
}

inline const mpc_class operator*(const mpf_class & a, const mpc_class & b)
{
    mpc_class tmp;
    tmp.real(a * b.real());
    tmp.imag(a * b.imag());
    return tmp;
}

inline const mpc_class operator*(const mpc_class & a, const mpf_class & b)
{
    mpc_class tmp;

    tmp.real(a.real() * b);
    tmp.imag(a.imag() * b);
    return tmp;
}

//not so bad as overflow might not occur with GMP; exponet range is extraordinarly large.
inline mpc_class & mpc_class::operator/=(const mpc_class & b)
{
    mpc_class tmp(*this);
    re = (tmp.re * b.re + tmp.im * b.im) / (b.re * b.re + b.im * b.im);
    im = (tmp.im * b.re - tmp.re * b.im) / (b.re * b.re + b.im * b.im);
    return (*this);
}

inline mpc_class & mpc_class::operator/=(const mpf_class & b)
{
    mpc_class tmp(*this);
    re = (tmp.re * b);
    im = (tmp.im * b);
    return (*this);
}

inline const mpc_class operator/(const mpc_class & a, const mpc_class & b)
{
    mpc_class tmp(a);
    tmp /= b;
    return tmp;
}

inline const mpc_class operator/(const mpc_class & a, const mpf_class & b)
{
    mpc_class tmp;
    tmp.real(a.real() / b);
    tmp.imag(a.imag() / b);
    return tmp;

}

inline const mpc_class operator/(const mpf_class & a, const mpc_class & b)
{
    mpc_class tmp;

    tmp.real((a * b.real()) / (b.real() * b.real() + b.imag() * b.imag()));
    tmp.imag(-(a * b.imag()) / (b.real() * b.real() + b.imag() * b.imag()));
    return tmp;
}


//not so bad as overflow might not occur with GMP; exponet range is extraordinarly large.
inline mpf_class abs(mpc_class ctmp)
{
    return sqrt(ctmp.real() * ctmp.real() + ctmp.imag() * ctmp.imag());
}

inline mpc_class sqrt(mpc_class z)
{
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

inline mpc_class conj(mpc_class ctmp)
{
    mpc_class ctmp2;

    ctmp2.real(ctmp.real());
    ctmp2.imag(-ctmp.imag());
    return ctmp2;
}

#endif
