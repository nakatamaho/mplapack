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
 * Device types for running the text of the reference complex BLAS
 * (../../common/blas_dev/) with mpfr_class/mpc_class of a fixed precision
 * PB: mreal<PB> on cu_fp::cu_freal<PB>, mcomplex<PB> on its components.
 * gmpfrxx_mkII evaluates every operator of the reference text as one MPFR or
 * MPC call rounded to nearest at the default precision, never fused:
 *
 *   mpfr_class + - *       mpfr_add/sub/mul       cu_freal + - *
 *   mpc_class + -          mpc_add/sub            componentwise
 *   mpc_class * mpc_class  mpc_mul                cu_cmul (correctly rounded
 *                                                 components, as MPC)
 *   mpc_class op mpfr_class the real operand becomes (r, +0), then mpc_*:
 *                          for *, componentwise products
 *   (expr).real()          the complex value, then its real part
 *   negation, conj         exact
 *   ==                     mpc_cmp == 0 (both components equal)
 *
 * so the results equal the reference ones bit for bit, except that
 * cu_freal zeros are +0.  The exponent range of cu_freal is unbounded: the
 * host keeps the operands inside a bound where MPFR never overflows or
 * underflows (Cblas_mpfr_bridge_cuda.cpp).  There is no correctly rounded
 * complex division: the division operators below must not run, and the
 * host never sends Ctrsm with diag = "N" (the only routine that divides).
 * This header is compiled by nvcc and the host compiler; it must not
 * include gmpfrxx_mkII or MPFR.
 */

#ifndef _MPLAPACK_CBLAS_MPFR_DEV_TYPES_H_
#define _MPLAPACK_CBLAS_MPFR_DEV_TYPES_H_

#include <mpc_cuda/cu_freal.cuh>
#include <mpc_cuda/cu_fcomplex.cuh>
#include "../../common/Rblas_dev_base.h"

namespace mplapack_blas_dev {

template <int PB> struct mpfr_fixed {}; // the tag of cgpu_types

// The arithmetic as calls: inlined into the reference text it makes the
// device code too large to compile in reasonable time.
#if defined(__CUDACC__)
#define MPLAPACK_MPFR_DEV_OP __host__ __device__ __noinline__
#else
#define MPLAPACK_MPFR_DEV_OP inline
#endif
template <int PB> MPLAPACK_MPFR_DEV_OP cu_fp::cu_freal<PB> mp_add(const cu_fp::cu_freal<PB> &a, const cu_fp::cu_freal<PB> &b) { return cu_fp::cu_fadd<PB>(a, b); }
template <int PB> MPLAPACK_MPFR_DEV_OP cu_fp::cu_freal<PB> mp_sub(const cu_fp::cu_freal<PB> &a, const cu_fp::cu_freal<PB> &b) { return cu_fp::cu_fsub<PB>(a, b); }
template <int PB> MPLAPACK_MPFR_DEV_OP cu_fp::cu_freal<PB> mp_mul(const cu_fp::cu_freal<PB> &a, const cu_fp::cu_freal<PB> &b) { return cu_fp::cu_fmul<PB>(a, b); }
template <int PB> MPLAPACK_MPFR_DEV_OP int mp_cmp(const cu_fp::cu_freal<PB> &a, const cu_fp::cu_freal<PB> &b) { return cu_fp::cu_cmp<PB>(a, b); }
template <int PB> MPLAPACK_MPFR_DEV_OP cu_fp::cu_freal<PB> mp_from_d(double d) { return cu_fp::cu_freal<PB>::from_double(d); }
// ar*br - ai*bi and ar*bi + ai*br, each rounded once (cu_cmul)
template <int PB> MPLAPACK_MPFR_DEV_OP cu_fp::cu_freal<PB> mp_cmul_re(const cu_fp::cu_freal<PB> &ar, const cu_fp::cu_freal<PB> &ai, const cu_fp::cu_freal<PB> &br, const cu_fp::cu_freal<PB> &bi) { return cu_fp::cu_fmms<PB>(ar, br, ai, bi); }
template <int PB> MPLAPACK_MPFR_DEV_OP cu_fp::cu_freal<PB> mp_cmul_im(const cu_fp::cu_freal<PB> &ar, const cu_fp::cu_freal<PB> &ai, const cu_fp::cu_freal<PB> &br, const cu_fp::cu_freal<PB> &bi) { return cu_fp::cu_fmma<PB>(ar, bi, ai, br); }

template <int PB> struct mreal {
    typedef cu_fp::cu_freal<PB> F;
    F v;
    MPLAPACK_DEV_HD mreal() { v.set_zero(); }
    MPLAPACK_DEV_HD mreal(double d) : v(mp_from_d<PB>(d)) {}
    MPLAPACK_DEV_HD explicit mreal(const F &f) : v(f) {}
    MPLAPACK_DEV_HD mreal operator-() const { return mreal(cu_fp::cu_neg<PB>(v)); }
    MPLAPACK_DEV_HD mreal operator+() const { return *this; }
    MPLAPACK_DEV_HD mreal &operator+=(const mreal &a) { v = mp_add<PB>(v, a.v); return *this; }
    MPLAPACK_DEV_HD mreal &operator-=(const mreal &a) { v = mp_sub<PB>(v, a.v); return *this; }
    MPLAPACK_DEV_HD mreal &operator*=(const mreal &a) { v = mp_mul<PB>(v, a.v); return *this; }
};

template <int PB> MPLAPACK_DEV_HD mreal<PB> operator+(const mreal<PB> &a, const mreal<PB> &b) { return mreal<PB>(mp_add<PB>(a.v, b.v)); }
template <int PB> MPLAPACK_DEV_HD mreal<PB> operator-(const mreal<PB> &a, const mreal<PB> &b) { return mreal<PB>(mp_sub<PB>(a.v, b.v)); }
template <int PB> MPLAPACK_DEV_HD mreal<PB> operator*(const mreal<PB> &a, const mreal<PB> &b) { return mreal<PB>(mp_mul<PB>(a.v, b.v)); }

template <int PB> MPLAPACK_DEV_HD int mreal_cmp(const mreal<PB> &a, const mreal<PB> &b) { return mp_cmp<PB>(a.v, b.v); }
template <int PB> MPLAPACK_DEV_HD bool operator==(const mreal<PB> &a, const mreal<PB> &b) { return mreal_cmp(a, b) == 0; }
template <int PB> MPLAPACK_DEV_HD bool operator!=(const mreal<PB> &a, const mreal<PB> &b) { return mreal_cmp(a, b) != 0; }
template <int PB> MPLAPACK_DEV_HD bool operator<(const mreal<PB> &a, const mreal<PB> &b) { return mreal_cmp(a, b) < 0; }
template <int PB> MPLAPACK_DEV_HD bool operator>(const mreal<PB> &a, const mreal<PB> &b) { return mreal_cmp(a, b) > 0; }
template <int PB> MPLAPACK_DEV_HD bool operator<=(const mreal<PB> &a, const mreal<PB> &b) { return mreal_cmp(a, b) <= 0; }
template <int PB> MPLAPACK_DEV_HD bool operator>=(const mreal<PB> &a, const mreal<PB> &b) { return mreal_cmp(a, b) >= 0; }
template <int PB> MPLAPACK_DEV_HD bool operator==(const mreal<PB> &a, double b) { return a == mreal<PB>(b); }
template <int PB> MPLAPACK_DEV_HD bool operator!=(const mreal<PB> &a, double b) { return a != mreal<PB>(b); }
template <int PB> MPLAPACK_DEV_HD mreal<PB> abs(const mreal<PB> &a) { return mreal<PB>(cu_fp::cu_abs<PB>(a.v)); }

template <int PB> struct mcomplex {
    typedef mreal<PB> R;
    R re_, im_;
    MPLAPACK_DEV_HD mcomplex() {}
    MPLAPACK_DEV_HD mcomplex(const R &re) : re_(re) {}
    MPLAPACK_DEV_HD mcomplex(const R &re, const R &im) : re_(re), im_(im) {}
    MPLAPACK_DEV_HD mcomplex(double re) : re_(re) {}
    MPLAPACK_DEV_HD mcomplex(double re, double im) : re_(re), im_(im) {}
    MPLAPACK_DEV_HD R &real() { return re_; }
    MPLAPACK_DEV_HD const R &real() const { return re_; }
    MPLAPACK_DEV_HD R &imag() { return im_; }
    MPLAPACK_DEV_HD const R &imag() const { return im_; }
    MPLAPACK_DEV_HD mcomplex &operator+=(const mcomplex &r) { re_ += r.re_; im_ += r.im_; return *this; }
    MPLAPACK_DEV_HD mcomplex &operator+=(const R &r) { re_ += r; return *this; }
    MPLAPACK_DEV_HD mcomplex &operator-=(const mcomplex &r) { re_ -= r.re_; im_ -= r.im_; return *this; }
    MPLAPACK_DEV_HD mcomplex &operator-=(const R &r) { re_ -= r; return *this; }
    MPLAPACK_DEV_HD mcomplex &operator*=(const mcomplex &r);
    MPLAPACK_DEV_HD mcomplex &operator*=(const R &r) { re_ *= r; im_ *= r; return *this; }
    MPLAPACK_DEV_HD mcomplex &operator/=(const mcomplex &r);
    MPLAPACK_DEV_HD mcomplex operator+() const { return *this; }
    MPLAPACK_DEV_HD mcomplex operator-() const { return mcomplex(-re_, -im_); }
};

template <int PB> MPLAPACK_DEV_HD bool operator==(const mcomplex<PB> &a, const mcomplex<PB> &b) { return a.real() == b.real() && a.imag() == b.imag(); }
template <int PB> MPLAPACK_DEV_HD bool operator==(const mcomplex<PB> &a, const mreal<PB> &b) { return a.real() == b && a.imag() == mreal<PB>(); }
template <int PB> MPLAPACK_DEV_HD bool operator==(const mreal<PB> &a, const mcomplex<PB> &b) { return b == a; }
template <int PB> MPLAPACK_DEV_HD bool operator!=(const mcomplex<PB> &a, const mcomplex<PB> &b) { return !(a == b); }
template <int PB> MPLAPACK_DEV_HD bool operator!=(const mcomplex<PB> &a, const mreal<PB> &b) { return !(a == b); }
template <int PB> MPLAPACK_DEV_HD bool operator!=(const mreal<PB> &a, const mcomplex<PB> &b) { return !(a == b); }

template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator+(mcomplex<PB> a, const mcomplex<PB> &b) { a += b; return a; }
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator+(mcomplex<PB> a, const mreal<PB> &b) { a += b; return a; }
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator+(const mreal<PB> &a, mcomplex<PB> b) { return mcomplex<PB>(a + b.real(), b.imag()); }
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator-(mcomplex<PB> a, const mcomplex<PB> &b) { a -= b; return a; }
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator-(mcomplex<PB> a, const mreal<PB> &b) { a -= b; return a; }
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator-(const mreal<PB> &a, const mcomplex<PB> &b) { return mcomplex<PB>(a - b.real(), -b.imag()); }
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator*(const mcomplex<PB> &a, const mcomplex<PB> &b)
{
    return mcomplex<PB>(mreal<PB>(mp_cmul_re<PB>(a.re_.v, a.im_.v, b.re_.v, b.im_.v)), mreal<PB>(mp_cmul_im<PB>(a.re_.v, a.im_.v, b.re_.v, b.im_.v)));
}
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator*(mcomplex<PB> a, const mreal<PB> &b) { a *= b; return a; }
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator*(const mreal<PB> &a, const mcomplex<PB> &b) { return mcomplex<PB>(a * b.real(), a * b.imag()); }
// Not computed (see above); the result is 0.
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator/(const mcomplex<PB> &, const mcomplex<PB> &) { return mcomplex<PB>(); }
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> &mcomplex<PB>::operator*=(const mcomplex<PB> &r) { *this = *this * r; return *this; }
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> &mcomplex<PB>::operator/=(const mcomplex<PB> &r) { *this = *this / r; return *this; }
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> conj(const mcomplex<PB> &z) { return mcomplex<PB>(z.real(), -z.imag()); }

// mpblas/reference/RCabs1.cpp
template <int PB> MPLAPACK_DEV_HD mreal<PB> RCabs1(const mcomplex<PB> &z)
{
    mreal<PB> return_value = 0.0;
    return_value = abs(z.real()) + abs(z.imag());
    return return_value;
}

template <int PB> struct cgpu_types<mpfr_fixed<PB> > {
    typedef mreal<PB> DR;
    typedef mcomplex<PB> DC;
};

} // namespace mplapack_blas_dev

#endif
