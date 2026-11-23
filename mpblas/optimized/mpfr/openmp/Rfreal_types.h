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

// fr<PB> and fc<PB>: stand-ins for mpfr_class and mpc_class of precision PB
// whose operators give the values of the MPFR/MPC operations of the
// libraries' expressions bit for bit (round to nearest, MPC_RNDNN),
// signed zeros included:
//   fr  +, -, *, / (mpfr_add, mpfr_sub, mpfr_mul, mpfr_div), unary -, abs,
//       sqrt, comparisons;
//   fc  +, - (per part), * (mpc_mul: correctly rounded parts; the zero rules
//       of its mul_real / mul_imag paths), conj, real(), imag(); a real
//       operand is (x, +0) (mpc_set_fr), as in the mpc_class expressions.
// Arithmetic on cu_freal<PB> (../lgpl/R*_host.h, LGPL-3.0-or-later) with
// unbounded exponents: a nonzero result outside [emin, emax] of
// mplapack_freal_tls (Rfreal_state.h) -- where MPFR overflows or underflows
// -- and a complex division (not correctly rounded here) are reported, and
// the call that made them is run by MPFR instead.
//
// The includer defines MPFREAL_NS (its own namespace) and has included
// ../lgpl/Rfreal_host.h, Rfround_host.h and Rfcomplex_host.h.

#include <cmath>
#include "Rfreal_state.h"

#if defined(__GNUC__) || defined(__clang__)
#define MPLAPACK_FREAL_NOINLINE __attribute__((noinline))
#else
#define MPLAPACK_FREAL_NOINLINE
#endif

namespace MPFREAL_NS {

template <int PB> struct fr {
    typedef cu_freal<PB> F;
    F v; // a zero is v.is_zero() with its sign in v.sign
    fr() { v.set_zero(); }
    fr(double d)
    {
        if (d == 0.0) {
            v.set_zero();
            v.sign = std::signbit(d) ? -1 : 1;
        } else
            v = F::from_double(d);
    }
    fr &operator=(double d) { return *this = fr(d); }
    bool is_zero() const { return v.is_zero(); }
    fr &operator+=(const fr &b);
    fr &operator-=(const fr &b);
    fr &operator*=(const fr &b);
    fr &operator/=(const fr &b);
};

template <int PB> inline void fr_check(const fr<PB> &r)
{
    if (!r.v.is_zero() && (r.v.exp > mplapack_freal_tls.emax || r.v.exp < mplapack_freal_tls.emin))
        mplapack_freal_fail();
}
template <int PB> inline fr<PB> fr_zero(int sign)
{
    fr<PB> r;
    r.v.sign = sign;
    return r;
}

template <int PB> MPLAPACK_FREAL_NOINLINE fr<PB> fr_mul(const fr<PB> &a, const fr<PB> &b)
{
    if (a.is_zero() || b.is_zero())
        return fr_zero<PB>(a.v.sign * b.v.sign);
    fr<PB> r;
    r.v = cu_fmul<PB>(a.v, b.v);
    fr_check(r);
    return r;
}
// a + s b, s = +1 or -1 (a sum of two zeros is -0 only if both are -0)
template <int PB> MPLAPACK_FREAL_NOINLINE fr<PB> fr_addsub(const fr<PB> &a, const fr<PB> &b, int s)
{
    if (b.is_zero()) {
        if (a.is_zero())
            return fr_zero<PB>((a.v.sign < 0 && s * b.v.sign < 0) ? -1 : 1);
        return a;
    }
    fr<PB> r;
    if (a.is_zero()) {
        r = b;
        r.v.sign *= s;
        return r;
    }
    if (s > 0)
        r.v = cu_fadd<PB>(a.v, b.v);
    else
        r.v = cu_fsub<PB>(a.v, b.v);
    fr_check(r); // an exact zero is +0, as MPFR's
    return r;
}
template <int PB> MPLAPACK_FREAL_NOINLINE fr<PB> fr_div(const fr<PB> &a, const fr<PB> &b)
{
    if (b.is_zero()) { // MPFR: inf or nan
        mplapack_freal_fail();
        return fr<PB>();
    }
    if (a.is_zero())
        return fr_zero<PB>(a.v.sign * b.v.sign);
    fr<PB> r;
    if (cu_fdiv_rn<PB>(r.v, a.v, b.v) != 0) {
        mplapack_freal_fail();
        return fr<PB>();
    }
    fr_check(r);
    return r;
}
template <int PB> MPLAPACK_FREAL_NOINLINE fr<PB> fr_sqrt(const fr<PB> &a)
{
    if (a.is_zero()) // mpfr_sqrt(-0) = -0
        return a;
    fr<PB> r;
    if (a.v.sign < 0 || cu_sqrt_rn<PB>(r.v, a.v) != 0) {
        mplapack_freal_fail();
        return fr<PB>();
    }
    fr_check(r);
    return r;
}

template <int PB> inline fr<PB> operator*(const fr<PB> &a, const fr<PB> &b) { return fr_mul(a, b); }
template <int PB> inline fr<PB> operator+(const fr<PB> &a, const fr<PB> &b) { return fr_addsub(a, b, 1); }
template <int PB> inline fr<PB> operator-(const fr<PB> &a, const fr<PB> &b) { return fr_addsub(a, b, -1); }
template <int PB> inline fr<PB> operator/(const fr<PB> &a, const fr<PB> &b) { return fr_div(a, b); }
template <int PB> inline fr<PB> operator-(const fr<PB> &a)
{
    fr<PB> r = a;
    r.v.sign = -r.v.sign;
    return r;
}
template <int PB> inline fr<PB> operator+(const fr<PB> &a) { return a; }
template <int PB> inline fr<PB> abs(const fr<PB> &a)
{
    fr<PB> r = a;
    r.v.sign = 1;
    return r;
}
template <int PB> inline fr<PB> sqrt(const fr<PB> &a) { return fr_sqrt(a); }
template <int PB> inline fr<PB> &fr<PB>::operator+=(const fr<PB> &b) { return *this = *this + b; }
template <int PB> inline fr<PB> &fr<PB>::operator-=(const fr<PB> &b) { return *this = *this - b; }
template <int PB> inline fr<PB> &fr<PB>::operator*=(const fr<PB> &b) { return *this = *this * b; }
template <int PB> inline fr<PB> &fr<PB>::operator/=(const fr<PB> &b) { return *this = *this / b; }

// comparisons (no NaN here; -0 == +0)
template <int PB> inline int fr_cmp(const fr<PB> &a, const fr<PB> &b)
{
    if (a.is_zero() && b.is_zero())
        return 0;
    return cu_cmp<PB>(a.v, b.v);
}
template <int PB> inline bool operator==(const fr<PB> &a, const fr<PB> &b)
{
    if (a.is_zero() || b.is_zero())
        return a.is_zero() && b.is_zero();
    if (a.v.sign != b.v.sign || a.v.exp != b.v.exp)
        return false;
    for (int i = 0; i < cu_freal<PB>::N; i++)
        if (a.v.m[i] != b.v.m[i])
            return false;
    return true;
}
template <int PB> inline bool operator!=(const fr<PB> &a, const fr<PB> &b) { return !(a == b); }
template <int PB> inline bool operator<(const fr<PB> &a, const fr<PB> &b) { return fr_cmp(a, b) < 0; }
template <int PB> inline bool operator>(const fr<PB> &a, const fr<PB> &b) { return fr_cmp(a, b) > 0; }
template <int PB> inline bool operator<=(const fr<PB> &a, const fr<PB> &b) { return fr_cmp(a, b) <= 0; }
template <int PB> inline bool operator>=(const fr<PB> &a, const fr<PB> &b) { return fr_cmp(a, b) >= 0; }
template <int PB> inline bool operator==(const fr<PB> &a, double b) { return b == 0.0 ? a.is_zero() : a == fr<PB>(b); }
template <int PB> inline bool operator!=(const fr<PB> &a, double b) { return !(a == b); }

template <int PB> struct fc {
    fr<PB> re, im;
    fc() {}
    fc(double r) : re(r), im() {}
    fc(double r, double i) : re(r), im(i) {}
    fc(const fr<PB> &r) : re(r), im() {}
    fc(const fr<PB> &r, const fr<PB> &i) : re(r), im(i) {}
    fc &operator=(const fr<PB> &r)
    {
        re = r;
        im = fr<PB>();
        return *this;
    }
    fc &operator=(double r) { return *this = fc(r); }
    const fr<PB> &real() const { return re; }
    const fr<PB> &imag() const { return im; }
    fc &operator+=(const fc &b);
    fc &operator-=(const fc &b);
    fc &operator*=(const fc &b);
    fc &operator/=(const fc &b);
};

// p q - r s (sign = -1) or p q + r s (sign = +1), rounded once, with the
// zero rules of the exact formula: the products with a zero factor are
// signed zeros, a nonzero product is rounded alone.
template <int PB> inline fr<PB> fc_part(const fr<PB> &p, const fr<PB> &q, const fr<PB> &r, const fr<PB> &s, int sign)
{
    const bool z1 = p.is_zero() || q.is_zero(), z2 = r.is_zero() || s.is_zero();
    if (z1 && z2) {
        const int s1 = p.v.sign * q.v.sign, s2 = sign * r.v.sign * s.v.sign;
        return fr_zero<PB>((s1 < 0 && s2 < 0) ? -1 : 1);
    }
    if (z2)
        return fr_mul(p, q);
    if (z1) {
        fr<PB> t = fr_mul(r, s);
        t.v.sign *= sign;
        return t;
    }
    fr<PB> t;
    t.v = sign < 0 ? cu_fmms<PB>(p.v, q.v, r.v, s.v) : cu_fmma<PB>(p.v, q.v, r.v, s.v);
    fr_check(t);
    return t;
}
template <int PB> MPLAPACK_FREAL_NOINLINE fc<PB> fc_mul(const fc<PB> &a, const fc<PB> &b)
{
    fc<PB> z;
    if (!a.re.is_zero() && !a.im.is_zero() && !b.re.is_zero() && !b.im.is_zero()) {
        cu_fcomplex<PB> x(a.re.v, a.im.v), y(b.re.v, b.im.v), w = cu_cmul<PB>(x, y);
        z.re.v = w.re, z.im.v = w.im;
        fr_check(z.re), fr_check(z.im);
        return z;
    }
    z.re = fc_part(a.re, b.re, a.im, b.im, -1);
    z.im = fc_part(a.re, b.im, a.im, b.re, 1);
    return z;
}
template <int PB> inline fc<PB> operator+(const fc<PB> &a, const fc<PB> &b) { return fc<PB>(a.re + b.re, a.im + b.im); }
template <int PB> inline fc<PB> operator-(const fc<PB> &a, const fc<PB> &b) { return fc<PB>(a.re - b.re, a.im - b.im); }
template <int PB> inline fc<PB> operator*(const fc<PB> &a, const fc<PB> &b) { return fc_mul(a, b); }
template <int PB> inline fc<PB> operator/(const fc<PB> &, const fc<PB> &)
{
    mplapack_freal_fail(); // mpc_div is not here
    return fc<PB>();
}
template <int PB> inline fc<PB> operator-(const fc<PB> &a) { return fc<PB>(-a.re, -a.im); }
template <int PB> inline fc<PB> operator+(const fc<PB> &a) { return a; }
template <int PB> inline fc<PB> conj(const fc<PB> &a) { return fc<PB>(a.re, -a.im); }
template <int PB> inline fr<PB> real(const fc<PB> &a) { return a.re; }
template <int PB> inline fr<PB> imag(const fc<PB> &a) { return a.im; }
// a real operand in a complex expression is (x, +0)
template <int PB> inline fc<PB> operator+(const fr<PB> &a, const fc<PB> &b) { return fc<PB>(a) + b; }
template <int PB> inline fc<PB> operator+(const fc<PB> &a, const fr<PB> &b) { return a + fc<PB>(b); }
template <int PB> inline fc<PB> operator-(const fr<PB> &a, const fc<PB> &b) { return fc<PB>(a) - b; }
template <int PB> inline fc<PB> operator-(const fc<PB> &a, const fr<PB> &b) { return a - fc<PB>(b); }
template <int PB> inline fc<PB> operator*(const fr<PB> &a, const fc<PB> &b) { return fc<PB>(a) * b; }
template <int PB> inline fc<PB> operator*(const fc<PB> &a, const fr<PB> &b) { return a * fc<PB>(b); }
template <int PB> inline fc<PB> operator/(const fc<PB> &a, const fr<PB> &b) { return a / fc<PB>(b); }
template <int PB> inline fc<PB> &fc<PB>::operator+=(const fc<PB> &b) { return *this = *this + b; }
template <int PB> inline fc<PB> &fc<PB>::operator-=(const fc<PB> &b) { return *this = *this - b; }
template <int PB> inline fc<PB> &fc<PB>::operator*=(const fc<PB> &b) { return *this = *this * b; }
template <int PB> inline fc<PB> &fc<PB>::operator/=(const fc<PB> &b) { return *this = *this / b; }
template <int PB> inline bool operator==(const fc<PB> &a, const fc<PB> &b) { return a.re == b.re && a.im == b.im; }
template <int PB> inline bool operator!=(const fc<PB> &a, const fc<PB> &b) { return !(a == b); }

} // namespace MPFREAL_NS
