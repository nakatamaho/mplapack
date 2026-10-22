/*
 * Copyright (c) 2026
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

/*
 * Device versions of dd_real/td_real/qd_real and of libQD3's qd3_complex,
 * for running the text of the reference BLAS routines on the GPU
 * (Rblas_dev_reference.h): every operator performs the operations of the
 * libQD3 operator it replaces (with QD_BF), so the device code computes the
 * reference results bit for bit (NaN signs aside).  Only the operators the
 * reference BLAS use are provided; the operations of a real with a double
 * are missing on purpose (they are separate libQD3 functions).
 */

#ifndef _MPLAPACK_RBLAS_DEV_TYPES_H_
#define _MPLAPACK_RBLAS_DEV_TYPES_H_

#include "Rgemm_cuda_ops_common.h"

namespace mplapack_blas_dev {

#define MPLAPACK_DEV_HD MPLAPACK_GEMM_HD

template <class Ops> struct dreal {
    typedef typename Ops::F F;
    F v;
    MPLAPACK_DEV_HD dreal() : v(Ops::from_d(0.0)) {}
    MPLAPACK_DEV_HD dreal(double d) : v(Ops::from_d(d)) {}
    MPLAPACK_DEV_HD explicit dreal(const F &f) : v(f) {}
    MPLAPACK_DEV_HD dreal operator-() const { return dreal(Ops::neg(v)); }
    MPLAPACK_DEV_HD dreal operator+() const { return *this; }
    MPLAPACK_DEV_HD dreal &operator+=(const dreal &a) { v = Ops::add(v, a.v); return *this; }
    MPLAPACK_DEV_HD dreal &operator-=(const dreal &a) { v = Ops::sub(v, a.v); return *this; }
    MPLAPACK_DEV_HD dreal &operator*=(const dreal &a) { v = Ops::mul(v, a.v); return *this; }
    MPLAPACK_DEV_HD dreal &operator/=(const dreal &a) { v = Ops::div(v, a.v); return *this; }
};

template <class Ops> MPLAPACK_DEV_HD dreal<Ops> operator+(const dreal<Ops> &a, const dreal<Ops> &b) { return dreal<Ops>(Ops::add(a.v, b.v)); }
template <class Ops> MPLAPACK_DEV_HD dreal<Ops> operator-(const dreal<Ops> &a, const dreal<Ops> &b) { return dreal<Ops>(Ops::sub(a.v, b.v)); }
template <class Ops> MPLAPACK_DEV_HD dreal<Ops> operator*(const dreal<Ops> &a, const dreal<Ops> &b) { return dreal<Ops>(Ops::mul(a.v, b.v)); }
template <class Ops> MPLAPACK_DEV_HD dreal<Ops> operator/(const dreal<Ops> &a, const dreal<Ops> &b) { return dreal<Ops>(Ops::div(a.v, b.v)); }

// comparisons: qd::*comparison_* (NaN: false; Inf: leading words; else the sign of a - b)
template <class Ops> MPLAPACK_DEV_HD int dreal_cmp(const dreal<Ops> &a, const dreal<Ops> &b, bool &unordered)
{
    unordered = Ops::is_nan(a.v) || Ops::is_nan(b.v);
    if (unordered)
        return 0;
    if (MPLAPACK_GEMM_ISINF(a.v.x[0]) || MPLAPACK_GEMM_ISINF(b.v.x[0]))
        return a.v.x[0] < b.v.x[0] ? -1 : (a.v.x[0] > b.v.x[0] ? 1 : 0);
    const typename Ops::F d = Ops::sub(a.v, b.v);
    return Ops::is_zero(d) ? 0 : (d.x[0] < 0.0 ? -1 : 1);
}
template <class Ops> MPLAPACK_DEV_HD bool operator==(const dreal<Ops> &a, const dreal<Ops> &b) { return Ops::eq(a.v, b.v); }
template <class Ops> MPLAPACK_DEV_HD bool operator!=(const dreal<Ops> &a, const dreal<Ops> &b) { return !Ops::eq(a.v, b.v); }
template <class Ops> MPLAPACK_DEV_HD bool operator<(const dreal<Ops> &a, const dreal<Ops> &b) { bool u; const int c = dreal_cmp(a, b, u); return !u && c < 0; }
template <class Ops> MPLAPACK_DEV_HD bool operator>(const dreal<Ops> &a, const dreal<Ops> &b) { bool u; const int c = dreal_cmp(a, b, u); return !u && c > 0; }
template <class Ops> MPLAPACK_DEV_HD bool operator<=(const dreal<Ops> &a, const dreal<Ops> &b) { bool u; const int c = dreal_cmp(a, b, u); return !u && c <= 0; }
template <class Ops> MPLAPACK_DEV_HD bool operator>=(const dreal<Ops> &a, const dreal<Ops> &b) { bool u; const int c = dreal_cmp(a, b, u); return !u && c >= 0; }

// a == double: libQD3 compares with the double converted to the type
template <class Ops> MPLAPACK_DEV_HD bool operator==(const dreal<Ops> &a, double b) { return a == dreal<Ops>(b); }
template <class Ops> MPLAPACK_DEV_HD bool operator!=(const dreal<Ops> &a, double b) { return a != dreal<Ops>(b); }

// abs: (a[0] < 0.0) ? -a : a
template <class Ops> MPLAPACK_DEV_HD dreal<Ops> abs(const dreal<Ops> &a) { return a.v.x[0] < 0.0 ? -a : a; }

// qd3_complex<Real>
template <class R> struct dcomplex {
    R re_, im_;
    MPLAPACK_DEV_HD dcomplex() : re_(0.0), im_(0.0) {}
    MPLAPACK_DEV_HD dcomplex(const R &re) : re_(re), im_(0.0) {}
    MPLAPACK_DEV_HD dcomplex(const R &re, const R &im) : re_(re), im_(im) {}
    MPLAPACK_DEV_HD dcomplex(double re) : re_(re), im_(0.0) {}
    MPLAPACK_DEV_HD dcomplex(double re, double im) : re_(re), im_(im) {}
    MPLAPACK_DEV_HD R &real() { return re_; }
    MPLAPACK_DEV_HD const R &real() const { return re_; }
    MPLAPACK_DEV_HD R &imag() { return im_; }
    MPLAPACK_DEV_HD const R &imag() const { return im_; }
    MPLAPACK_DEV_HD dcomplex &operator+=(const dcomplex &r) { re_ += r.re_; im_ += r.im_; return *this; }
    MPLAPACK_DEV_HD dcomplex &operator+=(const R &r) { re_ += r; return *this; }
    MPLAPACK_DEV_HD dcomplex &operator-=(const dcomplex &r) { re_ -= r.re_; im_ -= r.im_; return *this; }
    MPLAPACK_DEV_HD dcomplex &operator-=(const R &r) { re_ -= r; return *this; }
    MPLAPACK_DEV_HD dcomplex &operator*=(const dcomplex &r);
    MPLAPACK_DEV_HD dcomplex &operator*=(const R &r) { re_ *= r; im_ *= r; return *this; }
    MPLAPACK_DEV_HD dcomplex &operator/=(const dcomplex &r);
    MPLAPACK_DEV_HD dcomplex &operator/=(const R &r) { re_ /= r; im_ /= r; return *this; }
    MPLAPACK_DEV_HD dcomplex operator+() const { return *this; }
    MPLAPACK_DEV_HD dcomplex operator-() const { return dcomplex(-re_, -im_); }
};

template <class R> MPLAPACK_DEV_HD bool operator==(const dcomplex<R> &a, const dcomplex<R> &b) { return a.real() == b.real() && a.imag() == b.imag(); }
template <class R> MPLAPACK_DEV_HD bool operator==(const dcomplex<R> &a, const R &b) { return a.real() == b && a.imag() == R(0.0); }
template <class R> MPLAPACK_DEV_HD bool operator==(const R &a, const dcomplex<R> &b) { return b == a; }
template <class R> MPLAPACK_DEV_HD bool operator!=(const dcomplex<R> &a, const dcomplex<R> &b) { return !(a == b); }
template <class R> MPLAPACK_DEV_HD bool operator!=(const dcomplex<R> &a, const R &b) { return !(a == b); }
template <class R> MPLAPACK_DEV_HD bool operator!=(const R &a, const dcomplex<R> &b) { return !(a == b); }

template <class R> MPLAPACK_DEV_HD dcomplex<R> operator+(dcomplex<R> a, const dcomplex<R> &b) { a += b; return a; }
template <class R> MPLAPACK_DEV_HD dcomplex<R> operator+(dcomplex<R> a, const R &b) { a += b; return a; }
template <class R> MPLAPACK_DEV_HD dcomplex<R> operator+(const R &a, dcomplex<R> b) { b += a; return b; }
template <class R> MPLAPACK_DEV_HD dcomplex<R> operator-(dcomplex<R> a, const dcomplex<R> &b) { a -= b; return a; }
template <class R> MPLAPACK_DEV_HD dcomplex<R> operator-(dcomplex<R> a, const R &b) { a -= b; return a; }
template <class R> MPLAPACK_DEV_HD dcomplex<R> operator-(const R &a, const dcomplex<R> &b) { return dcomplex<R>(a - b.real(), -b.imag()); }
template <class R> MPLAPACK_DEV_HD dcomplex<R> operator*(const dcomplex<R> &a, const dcomplex<R> &b)
{
    return dcomplex<R>(a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real());
}
template <class R> MPLAPACK_DEV_HD dcomplex<R> operator*(dcomplex<R> a, const R &b) { a *= b; return a; }
template <class R> MPLAPACK_DEV_HD dcomplex<R> operator*(const R &a, dcomplex<R> b) { b *= a; return b; }
template <class R> MPLAPACK_DEV_HD dcomplex<R> operator/(const dcomplex<R> &a, const dcomplex<R> &b)
{
    const R abr = abs(b.real());
    const R abi = abs(b.imag());
    if (abr >= abi) {
        const R r = b.imag() / b.real();
        const R den = b.real() + b.imag() * r;
        return dcomplex<R>((a.real() + a.imag() * r) / den, (a.imag() - a.real() * r) / den);
    }
    const R r = b.real() / b.imag();
    const R den = b.real() * r + b.imag();
    return dcomplex<R>((a.real() * r + a.imag()) / den, (a.imag() * r - a.real()) / den);
}
template <class R> MPLAPACK_DEV_HD dcomplex<R> operator/(dcomplex<R> a, const R &b) { a /= b; return a; }
template <class R> MPLAPACK_DEV_HD dcomplex<R> &dcomplex<R>::operator*=(const dcomplex<R> &r) { *this = *this * r; return *this; }
template <class R> MPLAPACK_DEV_HD dcomplex<R> &dcomplex<R>::operator/=(const dcomplex<R> &r) { *this = *this / r; return *this; }
template <class R> MPLAPACK_DEV_HD dcomplex<R> conj(const dcomplex<R> &z) { return dcomplex<R>(z.real(), -z.imag()); }

// mpblas/reference/RCabs1.cpp
template <class R> MPLAPACK_DEV_HD R RCabs1(const dcomplex<R> &z)
{
    R return_value = 0.0;
    return_value = abs(z.real()) + abs(z.imag());
    return return_value;
}

// The functions the reference text calls (Mlsame, Mxerbla, max, min).
MPLAPACK_DEV_HD bool lsame(const char *a, const char *b)
{
    char x = a[0], y = b[0];
    if (x >= 'a' && x <= 'z')
        x = (char)(x - 'a' + 'A');
    if (y >= 'a' && y <= 'z')
        y = (char)(y - 'a' + 'A');
    return x == y;
}
MPLAPACK_DEV_HD void xerbla(const char *, int) {} // the host checks the arguments
template <class T> MPLAPACK_DEV_HD T max(T a, T b) { return a > b ? a : b; }
template <class T> MPLAPACK_DEV_HD T min(T a, T b) { return a < b ? a : b; }

} // namespace mplapack_blas_dev

#endif
