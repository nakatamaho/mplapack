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
 *   mpc_class / mpc_class  mpc_div                mpc_div (libmpc_cuda's cu_mpc_div)
 *   negation, conj         exact
 *   ==                     mpc_cmp == 0 (both components equal)
 *
 * so the results equal the reference ones bit for bit, except that
 * cu_freal zeros are +0.  The exponent range of cu_freal is unbounded: the
 * host keeps the operands inside a bound where MPFR never overflows or
 * underflows (Cblas_mpfr_bridge_cuda.cpp), and sends no zero divisor.
 * This header is compiled by nvcc and the host compiler; it does not
 * include gmpfrxx_mkII (the host emulation includes the system MPC).
 */

#ifndef _MPLAPACK_CBLAS_MPFR_DEV_TYPES_H_
#define _MPLAPACK_CBLAS_MPFR_DEV_TYPES_H_

#include <mpc_cuda/cu_freal.cuh>
#include <mpc_cuda/cu_fcomplex.cuh>
#include "../../common/Rblas_dev_base.h"
#include "../../common/Rblas_gpu_operand_common.h"

// Complex division: mpc_div (MPC_RNDNN), as gmpfrxx_mkII performs it, on the
// device with libmpc_cuda's cu_mpc_div and, in the host emulation, with the
// system MPC.  The operands are views of the cu_freal limbs.
#if defined(MPLAPACK_MPFR_CUDA_HOST_EMULATION)
#include <mpc.h>
#define MPLAPACK_MPFR_CDIV_HOST 1
#elif defined(__CUDACC__)
#include <mpc_cuda.cuh>
#define MPLAPACK_MPFR_CDIV_DEVICE 1
#endif

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

// r = a / b (complex) rounded to nearest per component; an infinite or NaN
// quotient (the host never sends a zero divisor) gives 0.
#define MPLAPACK_MPFR_CDIV_BODY(P, U, LIMB)                                                         \
    LIMB l[6][cu_fp::cu_freal<PB>::N];                                                              \
    P##mpc_t za, zb, zr;                                                                            \
    const cu_fp::cu_freal<PB> *in[4] = {&ar, &ai, &br, &bi};                                        \
    P##mpfr_ptr cs[4] = {P##mpc_realref(za), P##mpc_imagref(za), P##mpc_realref(zb), P##mpc_imagref(zb)}; \
    for (int q = 0; q < 4; q++) {                                                                   \
        if (in[q]->is_zero()) {                                                                     \
            P##mpfr_custom_init_set(cs[q], U##MPFR_ZERO_KIND, 0, PB, l[q]);                         \
        } else {                                                                                    \
            for (int t = 0; t < cu_fp::cu_freal<PB>::N; t++)                                        \
                l[q][t] = in[q]->m[t];                                                              \
            P##mpfr_custom_init_set(cs[q], in[q]->sign > 0 ? U##MPFR_REGULAR_KIND : -U##MPFR_REGULAR_KIND, in[q]->exp, PB, l[q]); \
        }                                                                                           \
    }                                                                                               \
    P##mpfr_custom_init_set(P##mpc_realref(zr), U##MPFR_ZERO_KIND, 0, PB, l[4]);                   \
    P##mpfr_custom_init_set(P##mpc_imagref(zr), U##MPFR_ZERO_KIND, 0, PB, l[5]);                   \
    P##mpc_div(zr, za, zb, U##MPC_RNDNN);                                                           \
    cu_fp::cu_freal<PB> *out[2] = {&rr, &ri};                                                       \
    P##mpfr_ptr os[2] = {P##mpc_realref(zr), P##mpc_imagref(zr)};                                   \
    for (int q = 0; q < 2; q++) {                                                                   \
        const int kind = P##mpfr_custom_get_kind(os[q]);                                            \
        out[q]->set_zero();                                                                         \
        if (kind == U##MPFR_REGULAR_KIND || kind == -U##MPFR_REGULAR_KIND) {                        \
            const LIMB *s = (const LIMB *)P##mpfr_custom_get_significand(os[q]);                    \
            out[q]->sign = kind > 0 ? 1 : -1;                                                       \
            out[q]->exp = (long)P##mpfr_custom_get_exp(os[q]);                                      \
            for (int t = 0; t < cu_fp::cu_freal<PB>::N; t++)                                        \
                out[q]->m[t] = s[t];                                                                \
        }                                                                                           \
    }
template <int PB> MPLAPACK_MPFR_DEV_OP void mp_cdiv(const cu_fp::cu_freal<PB> &ar, const cu_fp::cu_freal<PB> &ai, const cu_fp::cu_freal<PB> &br, const cu_fp::cu_freal<PB> &bi, cu_fp::cu_freal<PB> &rr, cu_fp::cu_freal<PB> &ri)
{
#if defined(MPLAPACK_MPFR_CDIV_HOST) && !defined(__CUDA_ARCH__)
    MPLAPACK_MPFR_CDIV_BODY(, , mp_limb_t)
#elif defined(MPLAPACK_MPFR_CDIV_DEVICE)
    MPLAPACK_MPFR_CDIV_BODY(cu_, CU_, cu_mp_limb_t)
#else
    (void)ar, (void)ai, (void)br, (void)bi;
    rr.set_zero();
    ri.set_zero();
#endif
}

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
template <int PB> MPLAPACK_DEV_HD mcomplex<PB> operator/(const mcomplex<PB> &a, const mcomplex<PB> &b)
{
    mcomplex<PB> r;
    mp_cdiv<PB>(a.re_.v, a.im_.v, b.re_.v, b.im_.v, r.re_.v, r.im_.v);
    return r;
}
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

// ---- runtime precision: the MPFR/MPC functions themselves ----
//
// mpfr_rt<NL>: mpfr_class/mpc_class values of the precision p of the call
// (p <= 64 NL bits), stored as kind, exponent and NL limbs; every operator
// calls the MPFR or MPC function gmpfrxx_mkII calls (libmpc_cuda's cu_mpfr_*
// and cu_mpc_* on the device, the system MPFR and MPC in the host
// emulation), rounded to nearest at p, so the results are the reference
// ones bit for bit, signed zeros included.  The host keeps the operands
// finite and inside the exponent bound (no overflow, underflow or division
// by zero).  p is mplapack_mpfr_rt_prec (Cblas_mpfr_device_cuda.cu).

template <int NL> struct mpfr_rt {}; // the tag of cgpu_types

// the limb counts of the runtime-precision types (p <= 64 NL bits): up to
// 4096 bits, and with MPLAPACK_MPFR_CUDA_RT_LARGE (their device code takes
// long to compile) 8192 and 16384 bits
#define MPLAPACK_MPFR_CUDA_RT_TIERS(X) X(2) X(4) X(8) X(16) X(32) X(64)
#define MPLAPACK_MPFR_CUDA_RT_LARGE_TIERS(X) X(128) X(256)

// Defined in Cblas_mpfr_device_cuda.cu: sets the precision, then cgpu_run.
template <class Job> int cgpu_run_rt(Job job, const mplapack_blas::gpu_operand<typename Job::E> *ops, int nops, long prec);

#if defined(MPLAPACK_MPFR_CDIV_HOST)
#define MPLAPACK_RT(x) x
#define MPLAPACK_RTU(x) x
typedef mp_limb_t mplapack_rt_limb;
#elif defined(MPLAPACK_MPFR_CDIV_DEVICE)
#define MPLAPACK_RT(x) cu_##x
#define MPLAPACK_RTU(x) CU_##x
typedef cu_mp_limb_t mplapack_rt_limb;
#else
typedef unsigned long long mplapack_rt_limb;
#endif
// the operators run where the MPFR functions are: on the device, or on the
// host in the emulation
#if (defined(MPLAPACK_MPFR_CDIV_HOST) && !defined(__CUDA_ARCH__)) || (defined(MPLAPACK_MPFR_CDIV_DEVICE) && defined(__CUDA_ARCH__))
#define MPLAPACK_RT_ACTIVE 1
#endif

extern long mplapack_mpfr_rt_prec_h;
#if defined(MPLAPACK_MPFR_CDIV_DEVICE)
extern __device__ long mplapack_mpfr_rt_prec_d;
#endif
MPLAPACK_DEV_HD long rt_prec()
{
#if defined(MPLAPACK_MPFR_CDIV_DEVICE) && defined(__CUDA_ARCH__)
    return mplapack_mpfr_rt_prec_d;
#else
    return mplapack_mpfr_rt_prec_h;
#endif
}

template <int NL> struct rreal {
    int kind; // MPFR custom kind (zero or regular, with the sign)
    long exp;
    mplapack_rt_limb m[NL];
    MPLAPACK_DEV_HD rreal();
    MPLAPACK_DEV_HD rreal(double d);
    MPLAPACK_DEV_HD rreal operator-() const;
    MPLAPACK_DEV_HD rreal operator+() const { return *this; }
    MPLAPACK_DEV_HD rreal &operator+=(const rreal &a);
    MPLAPACK_DEV_HD rreal &operator-=(const rreal &a);
    MPLAPACK_DEV_HD rreal &operator*=(const rreal &a);
};

#ifdef MPLAPACK_RT_ACTIVE
#define MPLAPACK_RT_VIEW(v, x) MPLAPACK_RT(mpfr_custom_init_set)(v, (x).kind, (x).exp, rt_prec(), const_cast<mplapack_rt_limb *>((x).m))
#define MPLAPACK_RT_DECL(v) MPLAPACK_RT(mpfr_t) v
// r := the value of v (whose significand is r.m)
#define MPLAPACK_RT_TAKE(r, v) ((r).kind = MPLAPACK_RT(mpfr_custom_get_kind)(v), (r).exp = (long)MPLAPACK_RT(mpfr_custom_get_exp)(v))
#endif

template <int NL> MPLAPACK_MPFR_DEV_OP void rt_set_d(rreal<NL> &r, double d)
{
#ifdef MPLAPACK_RT_ACTIVE
    MPLAPACK_RT_DECL(v);
    MPLAPACK_RT(mpfr_custom_init_set)(v, MPLAPACK_RTU(MPFR_ZERO_KIND), 0, rt_prec(), r.m);
    MPLAPACK_RT(mpfr_set_d)(v, d, MPLAPACK_RTU(MPFR_RNDN));
    MPLAPACK_RT_TAKE(r, v);
#else
    (void)d;
    r.kind = 2, r.exp = 0;
#endif
}
// op: 0 add, 1 sub, 2 mul
template <int NL> MPLAPACK_MPFR_DEV_OP void rt_op(rreal<NL> &r, const rreal<NL> &a, const rreal<NL> &b, int op)
{
#ifdef MPLAPACK_RT_ACTIVE
    MPLAPACK_RT_DECL(va);
    MPLAPACK_RT_DECL(vb);
    MPLAPACK_RT_DECL(vr);
    rreal<NL> t;
    MPLAPACK_RT_VIEW(va, a);
    MPLAPACK_RT_VIEW(vb, b);
    MPLAPACK_RT(mpfr_custom_init_set)(vr, MPLAPACK_RTU(MPFR_ZERO_KIND), 0, rt_prec(), t.m);
    if (op == 0)
        MPLAPACK_RT(mpfr_add)(vr, va, vb, MPLAPACK_RTU(MPFR_RNDN));
    else if (op == 1)
        MPLAPACK_RT(mpfr_sub)(vr, va, vb, MPLAPACK_RTU(MPFR_RNDN));
    else
        MPLAPACK_RT(mpfr_mul)(vr, va, vb, MPLAPACK_RTU(MPFR_RNDN));
    MPLAPACK_RT_TAKE(t, vr);
    r = t;
#else
    (void)a, (void)b, (void)op;
    r.kind = 2, r.exp = 0;
#endif
}
template <int NL> MPLAPACK_MPFR_DEV_OP int rt_cmp(const rreal<NL> &a, const rreal<NL> &b)
{
#ifdef MPLAPACK_RT_ACTIVE
    MPLAPACK_RT_DECL(va);
    MPLAPACK_RT_DECL(vb);
    MPLAPACK_RT_VIEW(va, a);
    MPLAPACK_RT_VIEW(vb, b);
    return MPLAPACK_RT(mpfr_cmp)(va, vb);
#else
    (void)a, (void)b;
    return 0;
#endif
}

template <int NL> MPLAPACK_DEV_HD rreal<NL>::rreal() : kind(2), exp(0) {} // +0 (MPFR_ZERO_KIND)
template <int NL> MPLAPACK_DEV_HD rreal<NL>::rreal(double d) { rt_set_d(*this, d); }
template <int NL> MPLAPACK_DEV_HD rreal<NL> rreal<NL>::operator-() const { rreal r = *this; r.kind = -r.kind; return r; } // mpfr_neg: exact
template <int NL> MPLAPACK_DEV_HD rreal<NL> &rreal<NL>::operator+=(const rreal &a) { rt_op(*this, *this, a, 0); return *this; }
template <int NL> MPLAPACK_DEV_HD rreal<NL> &rreal<NL>::operator-=(const rreal &a) { rt_op(*this, *this, a, 1); return *this; }
template <int NL> MPLAPACK_DEV_HD rreal<NL> &rreal<NL>::operator*=(const rreal &a) { rt_op(*this, *this, a, 2); return *this; }
template <int NL> MPLAPACK_DEV_HD rreal<NL> operator+(const rreal<NL> &a, const rreal<NL> &b) { rreal<NL> r; rt_op(r, a, b, 0); return r; }
template <int NL> MPLAPACK_DEV_HD rreal<NL> operator-(const rreal<NL> &a, const rreal<NL> &b) { rreal<NL> r; rt_op(r, a, b, 1); return r; }
template <int NL> MPLAPACK_DEV_HD rreal<NL> operator*(const rreal<NL> &a, const rreal<NL> &b) { rreal<NL> r; rt_op(r, a, b, 2); return r; }
template <int NL> MPLAPACK_DEV_HD bool operator==(const rreal<NL> &a, const rreal<NL> &b) { return rt_cmp(a, b) == 0; }
template <int NL> MPLAPACK_DEV_HD bool operator!=(const rreal<NL> &a, const rreal<NL> &b) { return rt_cmp(a, b) != 0; }
template <int NL> MPLAPACK_DEV_HD bool operator<(const rreal<NL> &a, const rreal<NL> &b) { return rt_cmp(a, b) < 0; }
template <int NL> MPLAPACK_DEV_HD bool operator>(const rreal<NL> &a, const rreal<NL> &b) { return rt_cmp(a, b) > 0; }
template <int NL> MPLAPACK_DEV_HD bool operator<=(const rreal<NL> &a, const rreal<NL> &b) { return rt_cmp(a, b) <= 0; }
template <int NL> MPLAPACK_DEV_HD bool operator>=(const rreal<NL> &a, const rreal<NL> &b) { return rt_cmp(a, b) >= 0; }
template <int NL> MPLAPACK_DEV_HD bool operator==(const rreal<NL> &a, double b) { return a == rreal<NL>(b); }
template <int NL> MPLAPACK_DEV_HD bool operator!=(const rreal<NL> &a, double b) { return a != rreal<NL>(b); }
template <int NL> MPLAPACK_DEV_HD rreal<NL> abs(const rreal<NL> &a) { rreal<NL> r = a; r.kind = r.kind < 0 ? -r.kind : r.kind; return r; } // mpfr_abs

template <int NL> struct rcomplex {
    typedef rreal<NL> R;
    R re_, im_;
    MPLAPACK_DEV_HD rcomplex() {}
    MPLAPACK_DEV_HD rcomplex(const R &re) : re_(re) {}
    MPLAPACK_DEV_HD rcomplex(const R &re, const R &im) : re_(re), im_(im) {}
    MPLAPACK_DEV_HD rcomplex(double re) : re_(re) {}
    MPLAPACK_DEV_HD rcomplex(double re, double im) : re_(re), im_(im) {}
    MPLAPACK_DEV_HD R &real() { return re_; }
    MPLAPACK_DEV_HD const R &real() const { return re_; }
    MPLAPACK_DEV_HD R &imag() { return im_; }
    MPLAPACK_DEV_HD const R &imag() const { return im_; }
    MPLAPACK_DEV_HD rcomplex &operator+=(const rcomplex &r);
    MPLAPACK_DEV_HD rcomplex &operator+=(const R &r);
    MPLAPACK_DEV_HD rcomplex &operator-=(const rcomplex &r);
    MPLAPACK_DEV_HD rcomplex &operator-=(const R &r);
    MPLAPACK_DEV_HD rcomplex &operator*=(const rcomplex &r);
    MPLAPACK_DEV_HD rcomplex &operator*=(const R &r);
    MPLAPACK_DEV_HD rcomplex &operator/=(const rcomplex &r);
    MPLAPACK_DEV_HD rcomplex operator+() const { return *this; }
    MPLAPACK_DEV_HD rcomplex operator-() const { return rcomplex(-re_, -im_); } // mpc_neg
};

// op: 0 mpc_add, 1 mpc_sub, 2 mpc_mul, 3 mpc_div; a real operand of
// gmpfrxx_mkII becomes (r, +0), as here
template <int NL> MPLAPACK_MPFR_DEV_OP void rt_cop(rcomplex<NL> &r, const rcomplex<NL> &a, const rcomplex<NL> &b, int op)
{
#ifdef MPLAPACK_RT_ACTIVE
    MPLAPACK_RT(mpc_t) za, zb, zr;
    rcomplex<NL> t;
    MPLAPACK_RT_VIEW(MPLAPACK_RT(mpc_realref)(za), a.re_);
    MPLAPACK_RT_VIEW(MPLAPACK_RT(mpc_imagref)(za), a.im_);
    MPLAPACK_RT_VIEW(MPLAPACK_RT(mpc_realref)(zb), b.re_);
    MPLAPACK_RT_VIEW(MPLAPACK_RT(mpc_imagref)(zb), b.im_);
    MPLAPACK_RT(mpfr_custom_init_set)(MPLAPACK_RT(mpc_realref)(zr), MPLAPACK_RTU(MPFR_ZERO_KIND), 0, rt_prec(), t.re_.m);
    MPLAPACK_RT(mpfr_custom_init_set)(MPLAPACK_RT(mpc_imagref)(zr), MPLAPACK_RTU(MPFR_ZERO_KIND), 0, rt_prec(), t.im_.m);
    if (op == 0)
        MPLAPACK_RT(mpc_add)(zr, za, zb, MPLAPACK_RTU(MPC_RNDNN));
    else if (op == 1)
        MPLAPACK_RT(mpc_sub)(zr, za, zb, MPLAPACK_RTU(MPC_RNDNN));
    else if (op == 2)
        MPLAPACK_RT(mpc_mul)(zr, za, zb, MPLAPACK_RTU(MPC_RNDNN));
    else
        MPLAPACK_RT(mpc_div)(zr, za, zb, MPLAPACK_RTU(MPC_RNDNN));
    MPLAPACK_RT_TAKE(t.re_, MPLAPACK_RT(mpc_realref)(zr));
    MPLAPACK_RT_TAKE(t.im_, MPLAPACK_RT(mpc_imagref)(zr));
    r = t;
#else
    (void)a, (void)b, (void)op;
    r = rcomplex<NL>();
#endif
}
template <int NL> MPLAPACK_MPFR_DEV_OP int rt_ccmp(const rcomplex<NL> &a, const rcomplex<NL> &b)
{
#ifdef MPLAPACK_RT_ACTIVE
    MPLAPACK_RT(mpc_t) za, zb;
    MPLAPACK_RT_VIEW(MPLAPACK_RT(mpc_realref)(za), a.re_);
    MPLAPACK_RT_VIEW(MPLAPACK_RT(mpc_imagref)(za), a.im_);
    MPLAPACK_RT_VIEW(MPLAPACK_RT(mpc_realref)(zb), b.re_);
    MPLAPACK_RT_VIEW(MPLAPACK_RT(mpc_imagref)(zb), b.im_);
    return MPLAPACK_RT(mpc_cmp)(za, zb);
#else
    (void)a, (void)b;
    return 0;
#endif
}
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> operator+(const rcomplex<NL> &a, const rcomplex<NL> &b) { rcomplex<NL> r; rt_cop(r, a, b, 0); return r; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> operator-(const rcomplex<NL> &a, const rcomplex<NL> &b) { rcomplex<NL> r; rt_cop(r, a, b, 1); return r; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> operator*(const rcomplex<NL> &a, const rcomplex<NL> &b) { rcomplex<NL> r; rt_cop(r, a, b, 2); return r; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> operator/(const rcomplex<NL> &a, const rcomplex<NL> &b) { rcomplex<NL> r; rt_cop(r, a, b, 3); return r; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> operator+(const rcomplex<NL> &a, const rreal<NL> &b) { return a + rcomplex<NL>(b); }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> operator+(const rreal<NL> &a, const rcomplex<NL> &b) { return rcomplex<NL>(a) + b; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> operator-(const rcomplex<NL> &a, const rreal<NL> &b) { return a - rcomplex<NL>(b); }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> operator-(const rreal<NL> &a, const rcomplex<NL> &b) { return rcomplex<NL>(a) - b; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> operator*(const rcomplex<NL> &a, const rreal<NL> &b) { return a * rcomplex<NL>(b); }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> operator*(const rreal<NL> &a, const rcomplex<NL> &b) { return rcomplex<NL>(a) * b; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> &rcomplex<NL>::operator+=(const rcomplex &r) { *this = *this + r; return *this; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> &rcomplex<NL>::operator+=(const R &r) { *this = *this + r; return *this; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> &rcomplex<NL>::operator-=(const rcomplex &r) { *this = *this - r; return *this; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> &rcomplex<NL>::operator-=(const R &r) { *this = *this - r; return *this; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> &rcomplex<NL>::operator*=(const rcomplex &r) { *this = *this * r; return *this; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> &rcomplex<NL>::operator*=(const R &r) { *this = *this * r; return *this; }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> &rcomplex<NL>::operator/=(const rcomplex &r) { *this = *this / r; return *this; }
template <int NL> MPLAPACK_DEV_HD bool operator==(const rcomplex<NL> &a, const rcomplex<NL> &b) { return rt_ccmp(a, b) == 0; }
template <int NL> MPLAPACK_DEV_HD bool operator==(const rcomplex<NL> &a, const rreal<NL> &b) { return a.real() == b && a.imag() == rreal<NL>(); }
template <int NL> MPLAPACK_DEV_HD bool operator==(const rreal<NL> &a, const rcomplex<NL> &b) { return b == a; }
template <int NL> MPLAPACK_DEV_HD bool operator!=(const rcomplex<NL> &a, const rcomplex<NL> &b) { return !(a == b); }
template <int NL> MPLAPACK_DEV_HD bool operator!=(const rcomplex<NL> &a, const rreal<NL> &b) { return !(a == b); }
template <int NL> MPLAPACK_DEV_HD bool operator!=(const rreal<NL> &a, const rcomplex<NL> &b) { return !(a == b); }
template <int NL> MPLAPACK_DEV_HD rcomplex<NL> conj(const rcomplex<NL> &z) { return rcomplex<NL>(z.real(), -z.imag()); } // mpc_conj
template <int NL> MPLAPACK_DEV_HD rreal<NL> RCabs1(const rcomplex<NL> &z)
{
    rreal<NL> return_value = 0.0;
    return_value = abs(z.real()) + abs(z.imag());
    return return_value;
}

template <int NL> struct cgpu_types<mpfr_rt<NL> > {
    typedef rreal<NL> DR;
    typedef rcomplex<NL> DC;
};

template <int PB> struct cgpu_types<mpfr_fixed<PB> > {
    typedef mreal<PB> DR;
    typedef mcomplex<PB> DC;
};

} // namespace mplapack_blas_dev

#endif
