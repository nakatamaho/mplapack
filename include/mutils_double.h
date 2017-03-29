/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mutils_double.h,v 1.11 2010/08/07 03:15:46 nakatamaho Exp $
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

#ifndef _MUTILS_DOUBLE_H_
#define _MUTILS_DOUBLE_H_

#if !(_XOPEN_SOURCE >= 600) && !(_ISOC99_SOURCE)
inline double log2(double x)
{
   return log(x) / log(2.0);
}
#endif

#include <complex>
double pi(double dummy);

double sign(double a, double b);
double cast2double(double a);
long nint(double a);
std::complex<double> Real2Complex(double a, double b);
double Cabs1(std::complex<double> zdum);

//implementation of sign transfer function.
inline double sign(double a, double b)
{
  double mtmp;
  mtmp=std::abs(a);
  if (b<0.0) {
    mtmp=-mtmp;
  }
  return mtmp;
}

inline double
cast2double(double a)
{
    return a;
}

inline long nint(double a)
{
    long i;
    double tmp;
    a = a + 0.5;
    tmp = floor(a);
    i = (long)tmp;
    return i;
}

inline std::complex<double> Real2Complex(double a, double b)
{
    std::complex<double> ret;
    ret.real() = a;
    ret.imag() = b;
    return ret;
}

inline double Cabs1(std::complex<double> zdum)
{
    double ret;
    ret = std::abs(zdum.real()) + std::abs(zdum.imag());
    return ret;
}

#endif
