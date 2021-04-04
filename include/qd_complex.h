/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: qd_complex.h,v 1.11 2010/08/07 03:15:46 nakatamaho Exp $
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
Complex class declare for the QD (double-double-double-double)
Above copyright (AT&T) only holds for complex division
*/

#ifndef _QD_COMPLEX_H_
#define _QD_COMPLEX_H_

#include <complex>

class qd_complex {
  private:
    qd_real re, im;

  public:
//constructor and deconstructor
    qd_complex();
    qd_complex(const qd_complex & a);
    qd_complex(const qd_real & a, const qd_real & b);
    qd_complex(const qd_real & a);
    qd_complex(const std::complex<double> & a);
    qd_complex(const double & a, const double & b);
    qd_complex(const double & a);

    ~qd_complex();
//extraction of real and imaginary parts
    inline qd_real & real()
    {
	return re;
    }
    inline const qd_real & real() const
    {
	return re;
    }
    inline qd_real & imag()
    {
	return im;
    }
    inline const qd_real & imag() const
    {
	return im;
    }
    inline void real(const qd_real r)
    {
	this->re = r;
    }
    inline void imag(const qd_real r)
    {
	this->im = r;
    }

//comparison
    friend bool operator==(const qd_complex & a, const qd_complex & b);
    friend bool operator==(const qd_complex & a, const qd_real & b);
    friend bool operator==(const qd_real & a, const qd_complex & b);

    friend bool operator!=(const qd_complex & a, const qd_complex & b);
    friend bool operator!=(const qd_complex & a, const qd_real & b);
    friend bool operator!=(const qd_real & a, const qd_complex & b);

//subsututtion
//difficult to implement; qd_complex& operator=(const qd_complex& b);
    qd_complex & operator=(const std::complex < double >b);
    qd_complex & operator=(const qd_real b);
    qd_complex & operator=(const double b);

//addition
    qd_complex & operator+=(const qd_complex & b);
    qd_complex & operator+=(const qd_real & b);
    const qd_complex operator+() const;

//subtraction
    qd_complex & operator-=(const qd_complex & b);
    qd_complex & operator-=(const qd_real & b);
    const qd_complex operator-() const;

//multiplication
    qd_complex & operator*=(const qd_complex & b);
    qd_complex & operator*=(const qd_real & b);

//division
    qd_complex & operator/=(const qd_complex & b);
    qd_complex & operator/=(const qd_real & b);

//operators
    operator std::complex<double> () const;
};

const qd_complex operator+(const qd_complex & a, const qd_complex & b);
const qd_complex operator+(const qd_complex & a, const qd_real & b);
const qd_complex operator+(const qd_real & a, const qd_complex & b);
const qd_complex operator+(const qd_complex & a, const std::complex < double >&b);
const qd_complex operator+(const std::complex < double >&a, const qd_complex & b);

const qd_complex operator-(const qd_complex & a, const qd_complex & b);
const qd_complex operator-(const qd_complex & a, const qd_real & b);
const qd_complex operator-(const qd_real & a, const qd_complex & b);
const qd_complex operator-(const qd_complex & a, std::complex < double >&b);
const qd_complex operator-(std::complex < double >&a, const qd_complex & b);

const qd_complex operator*(const qd_complex & a, const qd_complex & b);
const qd_complex operator*(const qd_complex & a, const qd_real & b);
const qd_complex operator*(const qd_real & a, const qd_complex & b);
const qd_complex operator*(const qd_complex & a, std::complex < double >&b);
const qd_complex operator*(std::complex < double >&a, const qd_complex & b);

const qd_complex operator/(const qd_complex & a, const qd_complex & b);
const qd_complex operator/(const qd_complex & a, const qd_real & b);
const qd_complex operator/(const qd_real & a, const qd_complex & b);
const qd_complex operator/(const qd_complex & a, std::complex < double >&b);
const qd_complex operator/(std::complex < double >&a, const qd_complex & b);

//constructor
inline qd_complex::qd_complex()
{
    re = 0.0;
    im = 0.0;
}

inline qd_complex::qd_complex(const qd_complex & a)
{
    re = a.re;
    im = a.im;
}

inline qd_complex::qd_complex(const qd_real & a, const qd_real & b)
{
    re = a;
    im = b;
}

inline qd_complex::qd_complex(const qd_real & a)
{
    re = a;
    im = 0.0;
}

inline qd_complex::qd_complex(const double& a, const double &b)
{
    re = a;
    im = b;
}

inline qd_complex::qd_complex(const double& a)
{
    re = a;
    im = 0.0;
}

inline qd_complex::qd_complex(const std::complex<double>& a)
{
    re = a.real();
    im = a.imag();
}

inline qd_complex::~qd_complex() { }

//comparison
inline bool operator==(const qd_complex & a, const qd_complex & b)
{
    if ((a.re == b.re) && (a.im == b.im)) {
	return true;
    } else
	return false;
}

inline bool operator==(const qd_real & a, const qd_complex & b)
{
    if ((a == b.re) && (b.im == 0.0)) {
	return true;
    } else
	return false;
}

inline bool operator==(const qd_complex & a, const qd_real & b)
{
    if ((a.re == b) && (a.im == 0.0)) {
	return true;
    } else
	return false;
}

inline bool operator!=(const qd_complex & a, const qd_complex & b)
{
    if ((a.re != b.re) || (a.im != b.im)) {
	return true;
    } else
	return false;
}

inline bool operator!=(const qd_real & a, const qd_complex & b)
{
    if ((a != b.re) || (b.im != 0.0)) {
	return true;
    } else
	return false;
}

inline bool operator!=(const qd_complex & a, const qd_real & b)
{
    if ((a.re != b) || (a.im != 0.0)) {
	return true;
    } else
	return false;
}

inline qd_complex & qd_complex::operator=(const std::complex < double >b)
{
    re = b.real();
    im = b.imag();
    return *this;
}

inline qd_complex & qd_complex::operator=(const qd_real b)
{
    re = b;
    im = 0.0;
    return *this;
}

inline qd_complex & qd_complex::operator=(const double b)
{
    re = b;
    im = 0.0;
    return *this;
}

inline qd_complex & qd_complex::operator+=(const qd_complex & b)
{
    re += b.re;
    im += b.im;
    return *this;
}

inline qd_complex & qd_complex::operator+=(const qd_real & b)
{
    re += b;
    return *this;
}

inline const qd_complex qd_complex::operator+() const
{
    return qd_complex(*this);
}

inline const qd_complex operator+(const qd_complex & a, const qd_complex & b)
{
    qd_complex tmp(a);
    tmp += b;
    return tmp;
}

inline const qd_complex operator+(const qd_complex & a, const std::complex<double> & b)
{
    qd_complex tmp(b);
    tmp += a;
    return tmp;
}

inline const qd_complex operator+(const std::complex<double> & a, const qd_complex &b)
{
    qd_complex tmp(a);
    tmp += b;
    return tmp;
}

inline const qd_complex operator+(const qd_complex & a, const qd_real & b)
{
    qd_complex tmp(b);
    tmp += a;
    return tmp;
}

inline const qd_complex operator+(const qd_real & a, const qd_complex & b)
{
    qd_complex tmp(a);
    tmp += b;
    return tmp;
}

inline qd_complex & qd_complex::operator-=(const qd_complex & b)
{
    re -= b.re;
    im -= b.im;
    return *this;
}

inline qd_complex & qd_complex::operator-=(const qd_real & b)
{
    re -= b;
    return *this;
}

inline const qd_complex qd_complex::operator-() const
{
    qd_complex tmp;
    tmp.re = -re;
    tmp.im = -im;
    return tmp;
}

inline const qd_complex operator-(const qd_complex & a, const qd_complex & b)
{
    qd_complex tmp(a);
    return tmp -= b;
}

//XXX can overflow
inline qd_complex & qd_complex::operator*=(const qd_complex & b)
{
    qd_complex tmp(*this);
    re = tmp.re * b.re - tmp.im * b.im;
    im = tmp.re * b.im + tmp.im * b.re;
    return (*this);
}

inline qd_complex & qd_complex::operator*=(const qd_real & b)
{
    qd_complex tmp(*this);
    re = tmp.re * b;
    im = tmp.im * b;
    return (*this);
}

inline const qd_complex operator*(const qd_complex & a, const qd_complex & b)
{
    qd_complex tmp(a);
    tmp *= b;
    return tmp;
}

inline qd_complex & qd_complex::operator/=(const qd_complex & b)
{
    qd_complex tmp(*this);
    qd_real abr, abi, ratio, den;

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

inline qd_complex & qd_complex::operator/=(const qd_real & b)
{
    qd_complex tmp(*this);
    re = (tmp.re * b);
    im = (tmp.im * b);
    return (*this);
}

inline qd_complex::operator std::complex<double> () const
{
  std::complex<double> p;
  p.real(re.x[0]);
  p.imag(im.x[0]);
  return p;
}

inline const qd_complex operator/(const qd_complex & a, const qd_complex & b)
{
    qd_complex tmp(a);
    tmp /= b;
    return tmp;
}
inline const qd_complex operator/(const qd_complex & a, const qd_real & b)
{
   qd_complex tmp;
   tmp.real(a.real() / b);
   tmp.imag(a.imag() / b);
   return tmp;
}
inline const qd_complex operator/(const qd_real & a, const qd_complex & b)
{
    qd_complex tmp(a);
    tmp /= b;
    return tmp;
}

inline const qd_complex operator*(const qd_complex & a, const qd_real & b)
{
  qd_complex p;
  p.real(a.real() * b);
  p.imag(a.imag() * b);
  return p;
}
inline const qd_complex operator*(const qd_real & a, const qd_complex & b)
{
  qd_complex p;
  p.real(a * b.real());
  p.imag(a * b.imag());
  return p;
}

inline const qd_complex operator-(const qd_complex & a, const qd_real & b)
{
  qd_complex p;
  p.real(a.real() - b);
  p.imag(a.imag());
  return p;
}
inline const qd_complex operator-(const qd_real & a, const qd_complex & b)
{
  qd_complex p;
  p.real(a - b.real());
  p.imag(-b.imag());
  return p;
}

//XXX can overflow
inline qd_real abs(qd_complex ctmp)
{
    return sqrt(ctmp.real() * ctmp.real() + ctmp.imag() * ctmp.imag());
}

inline qd_complex sqrt(qd_complex z)
{
    qd_real mag;
    qd_complex r;

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

inline qd_complex conj(qd_complex ctmp)
{
    qd_complex ctmp2;
    ctmp2.real(ctmp.real());
    ctmp2.imag(-ctmp.imag());
    return ctmp2;
}

#endif
