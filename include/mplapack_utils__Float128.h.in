/*
 * Copyright (c) 2012
 *	Nakata, Maho
 * 	All rights reserved.
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

#ifndef _MUTILS__FLOAT128_H_
#define _MUTILS__FLOAT128_H_

#ifdef __cplusplus
extern "C" {
#endif
#include <quadmath.h>
#ifdef __cplusplus
}
#endif

namespace std {
#ifndef _GLIBCXX_BITS_STD_ABS_H
inline _Float128 abs(const _Float128 &a) {
   return fabsq(a);
}
#endif

inline _Float128 sqrt(const _Float128 &a) {
   return sqrtq(a);
}
inline _Float128 log(const _Float128 &a) {
   return logq(a);
}
inline _Float128 log2(const _Float128 &a) {
   return log2q(a);
}
inline _Float128 pow(const _Float128 &a, const _Float128 &b) {
   return powq(a, b);
}

inline _Float128 exp(const _Float128 &a) {
   return expq(a);
}

inline _Float128 sin(const _Float128 &a) {
   return sinq(a);
}

inline _Float128 cos(const _Float128 &a) {
   return cosq(a);
}

inline _Float128 log10(const _Float128 &a) {
   return log10q(a);
}

}

#include <complex>

namespace std {

inline _Float128 abs(const std::complex<_Float128> &a) {
  return sqrtq(a.real() * a.real() + a.imag() * a.imag());
}


inline std::complex<_Float128> sqrt(const std::complex<_Float128> &a) {
   __complex128 b,tmp;
   std::complex<_Float128> c;
   __real__(b) = (a.real()); __imag__(b) = (a.imag());
   tmp = csqrtq (b);
   c.real( __real__(tmp)); c.imag(__imag__(tmp));
   return c;
}

inline std::complex<_Float128> exp(const std::complex<_Float128> &a) {
   __complex128 b,tmp;
   std::complex<_Float128> c;
   __real__(b) = (a.real()); __imag__(b) = (a.imag());
   tmp = cexpq (b);
   c.real(__real__(tmp)); c.imag(__imag__(tmp));
   return c;
}

inline std::complex<_Float128> log(const std::complex<_Float128> &a) {
   __complex128 b,tmp;
   std::complex<_Float128> c;
   __real__(b) = (a.real()); __imag__(b) = (a.imag());
   tmp = clogq (b);
   c.real(__real__(tmp)); c.imag(__imag__(tmp));
   return c;
}

}

#if !(_XOPEN_SOURCE >= 600) && !(_ISOC99_SOURCE)
inline _Float128 log2(_Float128 x)
{
   return log(x) / log(2.0);
}
#endif

_Float128 pi(_Float128 dummy);
_Float128 sign(_Float128 a, _Float128 b);
double cast2double(_Float128 a);
long nint(_Float128 a);
std::complex<_Float128> Real2Complex(_Float128 a, _Float128 b);
_Float128 Cabs1(std::complex<_Float128> zdum);

//implementation of sign transfer function.
inline _Float128 sign(_Float128 a, _Float128 b)
{
  _Float128 mtmp;
  mtmp=std::abs(a);
  if (b<0.0) {
    mtmp=-mtmp;
  }
  return mtmp;
}

inline double
cast2double(_Float128 a)
{
    return (double)a;
}

inline long nint(_Float128 a)
{
    long i;
    _Float128 tmp;
    a = a + 0.5;
    tmp = floorq(a);
    i = (long)tmp;
    return i;
}

inline std::complex<_Float128> Real2Complex(_Float128 a, _Float128 b)
{
    std::complex<_Float128> ret;
    ret.real(a);
    ret.imag(b);
    return ret;
}

inline _Float128 Cabs1(std::complex<_Float128> zdum)
{
    _Float128 ret;
    ret = std::abs(zdum.real()) + std::abs(zdum.imag());
    return ret;
}

#endif
