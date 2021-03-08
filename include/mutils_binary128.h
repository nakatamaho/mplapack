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

#ifndef _MUTILS_BINARY128_H_
#define _MUTILS_BINARY128_H_

#ifdef __cplusplus
extern "C" {
#endif
#include <quadmath.h>
#ifdef __cplusplus
}
#endif

namespace std {
#ifndef _GLIBCXX_BITS_STD_ABS_H
inline binary128 abs(const binary128 &a) {
   return fabsq(a);
}
#endif

inline binary128 sqrt(const binary128 &a) {
   return sqrtq(a);
}
inline binary128 log(const binary128 &a) {
   return logq(a);
}
inline binary128 log2(const binary128 &a) {
   return log2q(a);
}
inline binary128 pow(const binary128 &a, const binary128 &b) {
   return powq(a, b);
}

inline binary128 exp(const binary128 &a) {
   return expq(a);
}

inline binary128 sin(const binary128 &a) {
   return sinq(a);
}

inline binary128 cos(const binary128 &a) {
   return cosq(a);
}

inline binary128 log10(const binary128 &a) {
   return log10q(a);
}

}

#include <complex>

namespace std {

inline binary128 abs(const std::complex<binary128> &a) {
  return sqrtq(a.real() * a.real() + a.imag() * a.imag());
}


inline std::complex<binary128> sqrt(const std::complex<binary128> &a) {
   __complex128 b,tmp;
   std::complex<binary128> c;
   __real__(b) = (a.real()); __imag__(b) = (a.imag());
   tmp = csqrtq (b);
   c.real( __real__(tmp)); c.imag(__imag__(tmp));
   return c;
}

inline std::complex<binary128> exp(const std::complex<binary128> &a) {
   __complex128 b,tmp;
   std::complex<binary128> c;
   __real__(b) = (a.real()); __imag__(b) = (a.imag());
   tmp = cexpq (b);
   c.real(__real__(tmp)); c.imag(__imag__(tmp));
   return c;
}

inline std::complex<binary128> log(const std::complex<binary128> &a) {
   __complex128 b,tmp;
   std::complex<binary128> c;
   __real__(b) = (a.real()); __imag__(b) = (a.imag());
   tmp = clogq (b);
   c.real(__real__(tmp)); c.imag(__imag__(tmp));
   return c;
}

}

#if !(_XOPEN_SOURCE >= 600) && !(_ISOC99_SOURCE)
inline binary128 log2(binary128 x)
{
   return log(x) / log(2.0);
}
#endif

binary128 pi(binary128 dummy);
binary128 sign(binary128 a, binary128 b);
double cast2double(binary128 a);
long nint(binary128 a);
std::complex<binary128> Real2Complex(binary128 a, binary128 b);
binary128 Cabs1(std::complex<binary128> zdum);

//implementation of sign transfer function.
inline binary128 sign(binary128 a, binary128 b)
{
  binary128 mtmp;
  mtmp=std::abs(a);
  if (b<0.0) {
    mtmp=-mtmp;
  }
  return mtmp;
}

inline double
cast2double(binary128 a)
{
    return (double)a;
}

inline long nint(binary128 a)
{
    long i;
    binary128 tmp;
    a = a + 0.5;
    tmp = floorq(a);
    i = (long)tmp;
    return i;
}

inline std::complex<binary128> Real2Complex(binary128 a, binary128 b)
{
    std::complex<binary128> ret;
    ret.real(a);
    ret.imag(b);
    return ret;
}

inline binary128 Cabs1(std::complex<binary128> zdum)
{
    binary128 ret;
    ret = std::abs(zdum.real()) + std::abs(zdum.imag());
    return ret;
}

#endif
