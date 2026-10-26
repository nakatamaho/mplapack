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
 * Double-, triple- and quad-double arithmetic on plain doubles, usable on
 * the host and in CUDA device code.  Every function performs the operations
 * of the corresponding libQD3 inline function in the same order:
 *
 *   dd_bf_add  dd_real::bf_add (Kouya Algorithm 6, DWBFAdd)
 *   dd_bf_mul  dd_real::bf_mul (Algorithm 8, DWBFMul)
 *   td_bf_add  td_real::bf_add (TWBFAdd)    td_bf_mul  td_real::bf_mul (TWBFMul)
 *   qd_bf_add  qd_real::bf_add (QWBFAdd)    qd_bf_mul  qd_real::bf_mul (QWBFMul)
 *   dd_add_d / td_add_d / qd_add_d          operator+(x_real, double)
 *   dd_mul_d / td_mul_d / qd_mul_d          operator*(x_real, double)
 *   dd_div / td_div / qd_div                operator/(x_real, x_real) (accurate_div)
 *
 * and, for the opt-in sloppy GPU Rgemm only (<P>_CUDA_SLOPPY_ADD/_MUL):
 *
 *   dd_sloppy_add / td_sloppy_add / qd_sloppy_add   x_real::sloppy_add
 *   td_sloppy_mul  operator*(td_real, td_real) without QD_BF_MUL
 *   qd_sloppy_mul  qd_real::sloppy_mul
 *   dd_sloppy_mul  SLOPPY_MUL dd_mul of the former MPLAPACK CUDA kernels
 *
 * MPLAPACK computes dd, td and qd with the branch-free (BF) arithmetic of
 * libQD3: libQD3 must be built with QD_BF (QD_BF_ADD and QD_BF_MUL), as
 * external/qd does, so that operator+, operator- (a + (-b)) and operator*
 * of dd_real, td_real and qd_real are these functions and the results are
 * bit-identical to the C++ types.  The error-free product uses fma when
 * libQD3 does (QD_FMS defined in qd/qd_config.h), Dekker's splitting
 * otherwise.  Code using these functions must be compiled without
 * floating-point contraction (-ffp-contract=off, nvcc -fmad=false) and
 * without value-unsafe optimizations.
 */

#ifndef _MPLAPACK_RGEMM_QD_OPS_COMMON_H_
#define _MPLAPACK_RGEMM_QD_OPS_COMMON_H_

#include <cmath>
#include <limits>
#include <qd/qd_config.h>

#if defined(__CUDACC__)
#define MPLAPACK_GEMM_HD __host__ __device__ __forceinline__
#elif defined(__GNUC__)
#define MPLAPACK_GEMM_HD inline __attribute__((always_inline))
#else
#define MPLAPACK_GEMM_HD inline
#endif

#if defined(__CUDA_ARCH__)
#define MPLAPACK_GEMM_FMA(a, b, c) fma(a, b, c)
#define MPLAPACK_GEMM_FABS(a) fabs(a)
#define MPLAPACK_GEMM_ISINF(a) isinf(a)
#else
#define MPLAPACK_GEMM_FMA(a, b, c) std::fma(a, b, c)
#define MPLAPACK_GEMM_FABS(a) std::fabs(a)
#define MPLAPACK_GEMM_ISINF(a) std::isinf(a)
#endif

// The configuration of libQD3 these functions reproduce: BF addition and
// multiplication (QD_BF implies both).
#if defined(QD_BF) || (defined(QD_BF_ADD) && defined(QD_BF_MUL))
#define MPLAPACK_GEMM_QD_OPS_MATCH_LIBQD 1
#else
#error "MPLAPACK needs libQD3 built with branch-free arithmetic (QD_BF; CMake -DQD_ENABLE_BF=ON)"
#endif

namespace mplapack_gemm {

struct dd2 {
    double x[2];
};
struct qd4 {
    double x[4];
};
struct td3 {
    double x[3];
};

MPLAPACK_GEMM_HD double quick_two_sum(double a, double b, double &err) {
    double s = a + b;
    err = b - (s - a);
    return s;
}

MPLAPACK_GEMM_HD double two_sum(double a, double b, double &err) {
    double s = a + b;
    double bb = s - a;
    err = (a - (s - bb)) + (b - bb);
    return s;
}

MPLAPACK_GEMM_HD double two_diff(double a, double b, double &err) {
    double s = a - b;
    double bb = s - a;
    err = (a - (s - bb)) - (b + bb);
    return s;
}

#ifndef QD_FMS
MPLAPACK_GEMM_HD void split(double a, double &hi, double &lo) {
    const double thresh = 6.69692879491417e+299; // 2^996
    const double splitter = 134217729.0;          // 2^27 + 1
    double temp;
    if (a > thresh || a < -thresh) {
        a *= 3.7252902984619140625e-09; // 2^-28
        temp = splitter * a;
        hi = temp - (temp - a);
        lo = a - hi;
        hi *= 268435456.0; // 2^28
        lo *= 268435456.0;
    } else {
        temp = splitter * a;
        hi = temp - (temp - a);
        lo = a - hi;
    }
}
#endif

MPLAPACK_GEMM_HD double two_prod(double a, double b, double &err) {
    double p = a * b;
#ifdef QD_FMS
    err = MPLAPACK_GEMM_FMA(a, b, -p);
#else
    double a_hi, a_lo, b_hi, b_lo;
    split(a, a_hi, a_lo);
    split(b, b_hi, b_lo);
    err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
#endif
    return p;
}

// ---- double-double ----

// Kouya Algorithm 6, DWBFAdd (dd_real::bf_add)
MPLAPACK_GEMM_HD void dd_bf_add_w(double a0, double a1, double b0, double b1, double &r0, double &r1) {
    double g1, g1e, g2, g2e, g3, g3e, g4, g5;
    g1 = two_sum(a0, b0, g1e);
    g2 = two_sum(a1, b1, g2e);
    g3 = quick_two_sum(g1, g2, g3e);
    g4 = g1e + g2e;
    g5 = g4 + g3e;
    r0 = quick_two_sum(g3, g5, r1);
}

// Kouya Algorithm 8, DWBFMul (dd_real::bf_mul)
MPLAPACK_GEMM_HD void dd_bf_mul_w(double a0, double a1, double b0, double b1, double &r0, double &r1) {
    double p00, pe00, p01, p10, g1, g2;
    p00 = two_prod(a0, b0, pe00);
    p01 = a0 * b1;
    p10 = a1 * b0;
    g1 = p01 + p10;
    g2 = pe00 + g1;
    r0 = quick_two_sum(p00, g2, r1);
}

MPLAPACK_GEMM_HD dd2 dd_bf_add(const dd2 &a, const dd2 &b) {
    dd2 r;
    dd_bf_add_w(a.x[0], a.x[1], b.x[0], b.x[1], r.x[0], r.x[1]);
    return r;
}

MPLAPACK_GEMM_HD dd2 dd_bf_mul(const dd2 &a, const dd2 &b) {
    dd2 r;
    dd_bf_mul_w(a.x[0], a.x[1], b.x[0], b.x[1], r.x[0], r.x[1]);
    return r;
}

// dd + double (operator+(const dd_real &, double))
MPLAPACK_GEMM_HD dd2 dd_add_d(const dd2 &a, double b) {
    double s1, s2;
    s1 = two_sum(a.x[0], b, s2);
    s2 += a.x[1];
    dd2 r;
    r.x[0] = quick_two_sum(s1, s2, r.x[1]);
    return r;
}

// ---- quad-double ----

MPLAPACK_GEMM_HD void renorm(double &c0, double &c1, double &c2, double &c3, double &c4) {
    double s0, s1, s2 = 0.0, s3 = 0.0;
    if (MPLAPACK_GEMM_ISINF(c0))
        return;
    s0 = quick_two_sum(c3, c4, c4);
    s0 = quick_two_sum(c2, s0, c3);
    s0 = quick_two_sum(c1, s0, c2);
    c0 = quick_two_sum(c0, s0, c1);
    s0 = c0;
    s1 = c1;
    if (s1 != 0.0) {
        s1 = quick_two_sum(s1, c2, s2);
        if (s2 != 0.0) {
            s2 = quick_two_sum(s2, c3, s3);
            if (s3 != 0.0)
                s3 += c4;
            else
                s2 = quick_two_sum(s2, c4, s3);
        } else {
            s1 = quick_two_sum(s1, c3, s2);
            if (s2 != 0.0)
                s2 = quick_two_sum(s2, c4, s3);
            else
                s1 = quick_two_sum(s1, c4, s2);
        }
    } else {
        s0 = quick_two_sum(s0, c2, s1);
        if (s1 != 0.0) {
            s1 = quick_two_sum(s1, c3, s2);
            if (s2 != 0.0)
                s2 = quick_two_sum(s2, c4, s3);
            else
                s1 = quick_two_sum(s1, c4, s2);
        } else {
            s0 = quick_two_sum(s0, c3, s1);
            if (s1 != 0.0)
                s1 = quick_two_sum(s1, c4, s2);
            else
                s0 = quick_two_sum(s0, c4, s1);
        }
    }
    c0 = s0;
    c1 = s1;
    c2 = s2;
    c3 = s3;
}

// qd + double (operator+(const qd_real &, double))
MPLAPACK_GEMM_HD qd4 qd_add_d(const qd4 &a, double b) {
    double c0, c1, c2, c3;
    double e;
    c0 = two_sum(a.x[0], b, e);
    c1 = two_sum(a.x[1], e, e);
    c2 = two_sum(a.x[2], e, e);
    c3 = two_sum(a.x[3], e, e);
    renorm(c0, c1, c2, c3, e);
    qd4 r;
    r.x[0] = c0;
    r.x[1] = c1;
    r.x[2] = c2;
    r.x[3] = c3;
    return r;
}

// Kouya, branch-free quad-word addition (QWBFAdd)
MPLAPACK_GEMM_HD void qd_bf_add_w(double a_0, double a_1, double a_2, double a_3, double b_0, double b_1, double b_2, double b_3, double &r0, double &r1, double &r2, double &r3) {
    double a1, b1, c1, d1, e1, f1, g1, h1;
    double a2, b2, c2, d2, e2, f2, g2, b3, g3;
    double c3, d3, e3, f3, a4, c4, d4, e4, b5, d5, e5;
    double b6, c6, d6, e6, a7, b7, c7, d7, e8, b8, c8;
    double d9, b10, c10, d10, c11;

    a1 = two_sum(a_0, b_0, b1);
    c1 = two_sum(a_1, b_1, d1);
    e1 = two_sum(a_2, b_2, f1);
    g1 = two_sum(a_3, b_3, h1);
    a2 = quick_two_sum(a1, c1, c2);
    b2 = b1 + h1;
    d2 = two_sum(d1, e1, e2);
    f2 = two_sum(f1, g1, g2);
    b3 = two_sum(b2, g2, g3);
    c3 = quick_two_sum(c2, d2, d3);
    e3 = two_sum(e2, f2, f3);
    a4 = quick_two_sum(a2, c3, c4);
    d4 = quick_two_sum(d3, e3, e4);
    b5 = two_sum(b3, d4, d5);
    e5 = e4 + f3;
    b6 = two_sum(b5, c4, c6);
    d6 = two_sum(d5, e5, e6);
    a7 = quick_two_sum(a4, b6, b7);
    c7 = quick_two_sum(c6, d6, d7);
    e8 = e6 + g3;
    b8 = quick_two_sum(b7, c7, c8);
    d9 = d7 + e8;
    r0 = quick_two_sum(a7, b8, b10);
    c10 = quick_two_sum(c8, d9, d10);
    r1 = quick_two_sum(b10, c10, c11);
    r2 = quick_two_sum(c11, d10, r3);
}

MPLAPACK_GEMM_HD qd4 qd_bf_add(const qd4 &a, const qd4 &b) {
    qd4 r;
    qd_bf_add_w(a.x[0], a.x[1], a.x[2], a.x[3], b.x[0], b.x[1], b.x[2], b.x[3], r.x[0], r.x[1], r.x[2], r.x[3]);
    return r;
}

// Kouya, branch-free quad-word multiplication (QDBFMul)
MPLAPACK_GEMM_HD void qd_bf_mul_w(double a_0, double a_1, double a_2, double a_3, double b_0, double b_1, double b_2, double b_3, double &r0, double &r1, double &r2, double &r3) {
    double a0, b0, c0, e0, d0, f0, g0, j0, h0, k0, i0, l0;
    double m0, n0, o0, p0, c1, d1, e1, f1, g1, i1;
    double j1, m1, n1, b2, c2, e2, h2, f2, i2, m2;
    double a3, b3, c3, d3, e3, g3, f3, h3, c4, e4, d4, f4;
    double d5, c6, d6, b7, c7, d7, b8, c8, d8, c9;

    a0 = two_prod(a_0, b_0, b0);
    c0 = two_prod(a_0, b_1, e0);
    d0 = two_prod(a_1, b_0, f0);
    g0 = two_prod(a_0, b_2, j0);
    h0 = two_prod(a_1, b_1, k0);
    i0 = two_prod(a_2, b_0, l0);
    m0 = a_0 * b_3;
    n0 = a_1 * b_2;
    o0 = a_2 * b_1;
    p0 = a_3 * b_0;
    c1 = two_sum(c0, d0, d1);
    e1 = two_sum(e0, f0, f1);
    g1 = two_sum(g0, i0, i1);
    j1 = j0 + l0;
    m1 = m0 + p0;
    n1 = n0 + o0;
    b2 = two_sum(b0, c1, c2);
    e2 = two_sum(e1, h0, h2);
    f2 = f1 + j1;
    i2 = i1 + k0;
    m2 = m1 + n1;
    a3 = quick_two_sum(a0, b2, b3);
    c3 = quick_two_sum(c2, d1, d3);
    e3 = two_sum(e2, g1, g3);
    f3 = f2 + m2;
    h3 = h2 + i2;
    c4 = two_sum(c3, e3, e4);
    d4 = d3 + h3;
    f4 = f3 + g3;
    d5 = d4 + e4;
    c6 = two_sum(c4, d5, d6);
    b7 = two_sum(b3, c6, c7);
    d7 = d6 + f4;
    r0 = quick_two_sum(a3, b7, b8);
    c8 = two_sum(c7, d7, d8);
    r1 = two_sum(b8, c8, c9);
    r2 = quick_two_sum(c9, d8, r3);
    }

MPLAPACK_GEMM_HD qd4 qd_bf_mul(const qd4 &a, const qd4 &b) {
    qd4 r;
    qd_bf_mul_w(a.x[0], a.x[1], a.x[2], a.x[3], b.x[0], b.x[1], b.x[2], b.x[3], r.x[0], r.x[1], r.x[2], r.x[3]);
    return r;
}

// ---- triple-double (td:: namespace of libqd's td_inline.h) ----

MPLAPACK_GEMM_HD void td_renorm(double &c0, double &c1, double &c2, double &c3) {
    double s0, s1, s2 = 0.0;
    if (MPLAPACK_GEMM_ISINF(c0))
        return;
    s0 = quick_two_sum(c2, c3, c3);
    s0 = quick_two_sum(c1, s0, c2);
    c0 = quick_two_sum(c0, s0, c1);
    s0 = c0;
    s1 = c1;
    if (s1 != 0.0) {
        s1 = quick_two_sum(s1, c2, s2);
        if (s2 != 0.0) {
            s2 += c3;
        } else {
            s1 = quick_two_sum(s1, c3, s2);
        }
    } else {
        s0 = quick_two_sum(s0, c2, s1);
        if (s1 != 0.0) {
            s1 = quick_two_sum(s1, c3, s2);
        } else {
            s0 = quick_two_sum(s0, c3, s1);
        }
    }
    c0 = s0;
    c1 = s1;
    c2 = s2;
}

MPLAPACK_GEMM_HD td3 td_make(double a, double b, double c) {
    td3 r;
    r.x[0] = a;
    r.x[1] = b;
    r.x[2] = c;
    return r;
}

// operator+(const td_real &, double)
MPLAPACK_GEMM_HD td3 td_add_d(const td3 &a, double b) {
    double c0, c1, c2, e;
    c0 = two_sum(a.x[0], b, e);
    c1 = two_sum(a.x[1], e, e);
    c2 = two_sum(a.x[2], e, e);
    td_renorm(c0, c1, c2, e);
    return td_make(c0, c1, c2);
}

// Kouya, branch-free triple-word addition (TWBFAdd)
MPLAPACK_GEMM_HD void td_bf_add_w(double a_0, double a_1, double a_2, double b_0, double b_1, double b_2, double &r0, double &r1, double &r2) {
    double a1, b1, c1, d1, e1, f1;
    double a2, b2, c2, d2, e2, a3, d3;
    double b3, c3, c4, c5, d5, b6, c6, b7, c7;

    a1 = two_sum(a_0, b_0, b1);
    c1 = two_sum(a_1, b_1, d1);
    e1 = two_sum(a_2, b_2, f1);
    a2 = quick_two_sum(a1, c1, c2);
    b2 = b1 + f1;
    d2 = two_sum(d1, e1, e2);
    a3 = quick_two_sum(a2, d2, d3);
    b3 = two_sum(b2, c2, c3);
    c4 = c3 + e2;
    c5 = two_sum(c4, d3, d5);
    b6 = two_sum(b3, c5, c6);
    r0 = quick_two_sum(a3, b6, b7);
    c7 = c6 + d5;
    r1 = quick_two_sum(b7, c7, r2);
}

// Kouya, branch-free triple-word multiplication (TWBFMul)
MPLAPACK_GEMM_HD void td_bf_mul_w(double a_0, double a_1, double a_2, double b_0, double b_1, double b_2, double &r0, double &r1, double &r2) {
    double a0, b0, c0, e0, d0, f0, g0, h0, i0;
    double c1, d1, e1, g1, b2, c2, g2, a3, b3;
    double c3, e3, c4, b5, c5, b6;

    a0 = two_prod(a_0, b_0, b0);
    c0 = two_prod(a_0, b_1, e0);
    d0 = two_prod(a_1, b_0, f0);
    g0 = a_0 * b_2;
    h0 = a_1 * b_1;
    i0 = a_2 * b_0;
    c1 = two_sum(c0, d0, d1);
    e1 = e0 + f0;
    g1 = g0 + i0;
    b2 = two_sum(b0, c1, c2);
    g2 = g1 + h0;
    a3 = quick_two_sum(a0, b2, b3);
    c3 = c2 + d1;
    e3 = e1 + g2;
    c4 = c3 + e3;
    b5 = quick_two_sum(b3, c4, c5);
    r0 = quick_two_sum(a3, b5, b6);
    r1 = quick_two_sum(b6, c5, r2);
}

MPLAPACK_GEMM_HD td3 td_bf_add(const td3 &a, const td3 &b) {
    td3 r;
    td_bf_add_w(a.x[0], a.x[1], a.x[2], b.x[0], b.x[1], b.x[2], r.x[0], r.x[1], r.x[2]);
    return r;
}

MPLAPACK_GEMM_HD td3 td_bf_mul(const td3 &a, const td3 &b) {
    td3 r;
    td_bf_mul_w(a.x[0], a.x[1], a.x[2], b.x[0], b.x[1], b.x[2], r.x[0], r.x[1], r.x[2]);
    return r;
}


// ---- division (operator/ of libQD3 without QD_SLOPPY_DIV: accurate_div,
// with its special values and the rescaling of large numerators; the
// additions and subtractions inside are the BF ones) ----

MPLAPACK_GEMM_HD double qd_ops_nan() {
#if defined(__CUDA_ARCH__)
    return __longlong_as_double(0x7ff8000000000000LL);
#else
    return std::numeric_limits<double>::quiet_NaN();
#endif
}
MPLAPACK_GEMM_HD double qd_ops_inf() {
#if defined(__CUDA_ARCH__)
    return __longlong_as_double(0x7ff0000000000000LL);
#else
    return std::numeric_limits<double>::infinity();
#endif
}
MPLAPACK_GEMM_HD bool qd_ops_isnan(double a) {
#if defined(__CUDA_ARCH__)
    return isnan(a);
#else
    return std::isnan(a);
#endif
}
MPLAPACK_GEMM_HD bool qd_ops_signbit(double a) {
#if defined(__CUDA_ARCH__)
    return signbit(a);
#else
    return std::signbit(a);
#endif
}
MPLAPACK_GEMM_HD double qd_ops_copysign(double a, double b) {
#if defined(__CUDA_ARCH__)
    return copysign(a, b);
#else
    return std::copysign(a, b);
#endif
}

// Special cases of x_div_special: 0 = none, 1 = NaN, 2 = +-Inf, 3 = +-0;
// the sign is negative if the signs of the leading words differ.
MPLAPACK_GEMM_HD int qd_ops_div_special(bool a_nan, bool b_nan, bool a_inf, bool b_inf, bool a_zero, bool b_zero) {
    if (a_nan || b_nan || (a_inf && b_inf) || (a_zero && b_zero))
        return 1;
    if (b_zero || a_inf)
        return 2;
    if (a_zero || b_inf)
        return 3;
    return 0;
}

// operator*(const dd_real &, double)
MPLAPACK_GEMM_HD dd2 dd_mul_d(const dd2 &a, double b) {
    double p1, p2;
    p1 = two_prod(a.x[0], b, p2);
    p2 += (a.x[1] * b);
    dd2 r;
    r.x[0] = quick_two_sum(p1, p2, r.x[1]);
    return r;
}

MPLAPACK_GEMM_HD dd2 dd_fill(double v) {
    dd2 r;
    r.x[0] = r.x[1] = v;
    return r;
}

MPLAPACK_GEMM_HD dd2 dd_neg(const dd2 &a) {
    dd2 r;
    r.x[0] = -a.x[0];
    r.x[1] = -a.x[1];
    return r;
}

// dd_real::accurate_div
MPLAPACK_GEMM_HD dd2 dd_div(const dd2 &a, const dd2 &b) {
    const bool negative = qd_ops_signbit(a.x[0]) != qd_ops_signbit(b.x[0]);
    const int sp = qd_ops_div_special(qd_ops_isnan(a.x[0]) || qd_ops_isnan(a.x[1]), qd_ops_isnan(b.x[0]) || qd_ops_isnan(b.x[1]),
                                      MPLAPACK_GEMM_ISINF(a.x[0]), MPLAPACK_GEMM_ISINF(b.x[0]), a.x[0] == 0.0, b.x[0] == 0.0);
    if (sp == 1)
        return dd_fill(qd_ops_nan());
    if (sp == 2)
        return negative ? dd_fill(-qd_ops_inf()) : dd_fill(qd_ops_inf());
    if (sp == 3) {
        dd2 z;
        z.x[0] = qd_ops_copysign(0.0, negative ? -1.0 : 1.0);
        z.x[1] = 0.0;
        return z;
    }
    const bool rescale = MPLAPACK_GEMM_FABS(a.x[0]) > 0x1p+969;
    dd2 aa = a;
    if (rescale)
        aa.x[0] = a.x[0] * 0x1p-53, aa.x[1] = a.x[1] * 0x1p-53;
    double q1 = aa.x[0] / b.x[0];
    if (MPLAPACK_GEMM_ISINF(q1))
        return qd_ops_signbit(q1) ? dd_fill(-qd_ops_inf()) : dd_fill(qd_ops_inf());
    dd2 r = dd_bf_add(aa, dd_neg(dd_mul_d(b, q1)));
    double q2 = r.x[0] / b.x[0];
    r = dd_bf_add(r, dd_neg(dd_mul_d(b, q2)));
    double q3 = r.x[0] / b.x[0];
    q1 = quick_two_sum(q1, q2, q2);
    dd2 t;
    t.x[0] = q1;
    t.x[1] = q2;
    r = dd_add_d(t, q3);
    if (rescale)
        r.x[0] = r.x[0] * 0x1p+53, r.x[1] = r.x[1] * 0x1p+53;
    if (MPLAPACK_GEMM_ISINF(r.x[0]))
        return qd_ops_signbit(r.x[0]) ? dd_fill(-qd_ops_inf()) : dd_fill(qd_ops_inf());
    return r;
}

// operator*(const td_real &, double)
MPLAPACK_GEMM_HD td3 td_mul_d(const td3 &a, double b) {
    double p0, p1, p2, q0, q1, q2;
    p0 = two_prod(a.x[0], b, q0);
    p1 = two_prod(a.x[1], b, q1);
    p2 = two_prod(a.x[2], b, q2);
    p1 = two_sum(p1, q0, q0);
    p2 = two_sum(p2, q1, q1);
    q0 += q1 + q2;
    td_renorm(p0, p1, p2, q0);
    return td_make(p0, p1, p2);
}

MPLAPACK_GEMM_HD td3 td_neg(const td3 &a) { return td_make(-a.x[0], -a.x[1], -a.x[2]); }

// td_real::accurate_div
MPLAPACK_GEMM_HD td3 td_div(const td3 &a, const td3 &b) {
    const bool negative = qd_ops_signbit(a.x[0]) != qd_ops_signbit(b.x[0]);
    const int sp = qd_ops_div_special(qd_ops_isnan(a.x[0]) || qd_ops_isnan(a.x[1]) || qd_ops_isnan(a.x[2]),
                                      qd_ops_isnan(b.x[0]) || qd_ops_isnan(b.x[1]) || qd_ops_isnan(b.x[2]),
                                      MPLAPACK_GEMM_ISINF(a.x[0]), MPLAPACK_GEMM_ISINF(b.x[0]),
                                      a.x[0] == 0.0 && a.x[1] == 0.0 && a.x[2] == 0.0, b.x[0] == 0.0 && b.x[1] == 0.0 && b.x[2] == 0.0);
    const double inf = qd_ops_inf();
    if (sp == 1)
        return td_make(qd_ops_nan(), qd_ops_nan(), qd_ops_nan());
    if (sp == 2)
        return negative ? td_make(-inf, -inf, -inf) : td_make(inf, inf, inf);
    if (sp == 3)
        return td_make(qd_ops_copysign(0.0, negative ? -1.0 : 1.0), 0.0, 0.0);
    const bool rescale = MPLAPACK_GEMM_FABS(a.x[0]) > 0x1p+969;
    const td3 aa = rescale ? td_make(a.x[0] * 0x1p-53, a.x[1] * 0x1p-53, a.x[2] * 0x1p-53) : a;
    double q0 = aa.x[0] / b.x[0];
    if (MPLAPACK_GEMM_ISINF(q0))
        return qd_ops_signbit(q0) ? td_make(-inf, -inf, -inf) : td_make(inf, inf, inf);
    td3 r = td_bf_add(aa, td_neg(td_mul_d(b, q0)));
    double q1 = r.x[0] / b.x[0];
    r = td_bf_add(r, td_neg(td_mul_d(b, q1)));
    double q2 = r.x[0] / b.x[0];
    r = td_bf_add(r, td_neg(td_mul_d(b, q2)));
    double q3 = r.x[0] / b.x[0];
    td_renorm(q0, q1, q2, q3);
    td3 res = td_make(q0, q1, q2);
    if (rescale)
        res = td_make(res.x[0] * 0x1p+53, res.x[1] * 0x1p+53, res.x[2] * 0x1p+53);
    if (MPLAPACK_GEMM_ISINF(res.x[0]))
        return qd_ops_signbit(res.x[0]) ? td_make(-inf, -inf, -inf) : td_make(inf, inf, inf);
    return res;
}

MPLAPACK_GEMM_HD void qd_three_sum(double &a, double &b, double &c) {
    double t1, t2, t3;
    t1 = two_sum(a, b, t2);
    a = two_sum(c, t1, t3);
    b = two_sum(t2, t3, c);
}

MPLAPACK_GEMM_HD void qd_three_sum2(double &a, double &b, double &c) {
    double t1, t2, t3;
    t1 = two_sum(a, b, t2);
    a = two_sum(c, t1, t3);
    b = t2 + t3;
}

MPLAPACK_GEMM_HD qd4 qd_make(double a, double b, double c, double d) {
    qd4 r;
    r.x[0] = a;
    r.x[1] = b;
    r.x[2] = c;
    r.x[3] = d;
    return r;
}

// operator*(const qd_real &, double)
MPLAPACK_GEMM_HD qd4 qd_mul_d(const qd4 &a, double b) {
    double p0, p1, p2, p3, q0, q1, q2, s0, s1, s2, s3, s4;
    p0 = two_prod(a.x[0], b, q0);
    p1 = two_prod(a.x[1], b, q1);
    p2 = two_prod(a.x[2], b, q2);
    p3 = a.x[3] * b;
    s0 = p0;
    s1 = two_sum(q0, p1, s2);
    qd_three_sum(s2, q1, p2);
    qd_three_sum2(q1, q2, p3);
    s3 = q1;
    s4 = q2 + p2;
    renorm(s0, s1, s2, s3, s4);
    return qd_make(s0, s1, s2, s3);
}

MPLAPACK_GEMM_HD qd4 qd_neg(const qd4 &a) { return qd_make(-a.x[0], -a.x[1], -a.x[2], -a.x[3]); }

// qd_real::accurate_div
MPLAPACK_GEMM_HD qd4 qd_div(const qd4 &a, const qd4 &b) {
    const bool negative = qd_ops_signbit(a.x[0]) != qd_ops_signbit(b.x[0]);
    const int sp = qd_ops_div_special(qd_ops_isnan(a.x[0]) || qd_ops_isnan(a.x[1]) || qd_ops_isnan(a.x[2]) || qd_ops_isnan(a.x[3]),
                                      qd_ops_isnan(b.x[0]) || qd_ops_isnan(b.x[1]) || qd_ops_isnan(b.x[2]) || qd_ops_isnan(b.x[3]),
                                      MPLAPACK_GEMM_ISINF(a.x[0]), MPLAPACK_GEMM_ISINF(b.x[0]), a.x[0] == 0.0, b.x[0] == 0.0);
    const double inf = qd_ops_inf(), nan = qd_ops_nan();
    if (sp == 1)
        return qd_make(nan, nan, nan, nan);
    if (sp == 2)
        return negative ? qd_make(-inf, -inf, -inf, -inf) : qd_make(inf, inf, inf, inf);
    if (sp == 3)
        return qd_make(qd_ops_copysign(0.0, negative ? -1.0 : 1.0), 0.0, 0.0, 0.0);
    const bool rescale = MPLAPACK_GEMM_FABS(a.x[0]) > 0x1p+969;
    const qd4 aa = rescale ? qd_make(a.x[0] * 0x1p-53, a.x[1] * 0x1p-53, a.x[2] * 0x1p-53, a.x[3] * 0x1p-53) : a;
    double q0 = aa.x[0] / b.x[0];
    if (MPLAPACK_GEMM_ISINF(q0))
        return qd_ops_signbit(q0) ? qd_make(-inf, -inf, -inf, -inf) : qd_make(inf, inf, inf, inf);
    qd4 r = qd_bf_add(aa, qd_neg(qd_mul_d(b, q0)));
    double q1 = r.x[0] / b.x[0];
    r = qd_bf_add(r, qd_neg(qd_mul_d(b, q1)));
    double q2 = r.x[0] / b.x[0];
    r = qd_bf_add(r, qd_neg(qd_mul_d(b, q2)));
    double q3 = r.x[0] / b.x[0];
    r = qd_bf_add(r, qd_neg(qd_mul_d(b, q3)));
    double q4 = r.x[0] / b.x[0];
    renorm(q0, q1, q2, q3, q4);
    qd4 res = qd_make(q0, q1, q2, q3);
    if (rescale)
        res = qd_make(res.x[0] * 0x1p+53, res.x[1] * 0x1p+53, res.x[2] * 0x1p+53, res.x[3] * 0x1p+53);
    if (MPLAPACK_GEMM_ISINF(res.x[0]))
        return qd_ops_signbit(res.x[0]) ? qd_make(-inf, -inf, -inf, -inf) : qd_make(inf, inf, inf, inf);
    return res;
}


// ---- sloppy arithmetic (GPU Rgemm with <P>_CUDA_SLOPPY_ADD / _MUL) ----
// Not used by default: the results differ from the C++ types.

// dd_real::sloppy_add (Cray-style error bound)
MPLAPACK_GEMM_HD dd2 dd_sloppy_add(const dd2 &a, const dd2 &b) {
    double s, e;
    s = two_sum(a.x[0], b.x[0], e);
    e += (a.x[1] + b.x[1]);
    dd2 r;
    r.x[0] = quick_two_sum(s, e, r.x[1]);
    return r;
}

// The SLOPPY_MUL dd_mul of the former MPLAPACK CUDA kernels (libQD3 has no
// sloppy dd multiplication): the low-order products are added to the fma
// of the leading ones, without renormalization.
MPLAPACK_GEMM_HD dd2 dd_sloppy_mul(const dd2 &a, const dd2 &b) {
    double m1, m2, t, ph, e;
    m1 = a.x[0] * b.x[1];
    m2 = a.x[1] * b.x[0];
    t = m1 + m2;
    ph = MPLAPACK_GEMM_FMA(a.x[0], b.x[0], t);
    e = MPLAPACK_GEMM_FMA(a.x[0], b.x[0], -ph);
    dd2 r;
    r.x[0] = ph;
    r.x[1] = e + t;
    return r;
}

// td::three_sum
MPLAPACK_GEMM_HD void td_three_sum(double &a, double &b, double &c) {
    double t1, t2, t3;
    t1 = two_sum(a, b, t2);
    a = two_sum(c, t1, t3);
    b = two_sum(t2, t3, c);
}

// td::renorm(c0, c1, c2, c3, c4)
MPLAPACK_GEMM_HD void td_renorm5(double &c0, double &c1, double &c2, double &c3, double &c4) {
    double s0, s1, s2 = 0.0, s3 = 0.0;
    if (MPLAPACK_GEMM_ISINF(c0))
        return;
    s0 = quick_two_sum(c3, c4, c4);
    s0 = quick_two_sum(c2, s0, c3);
    s0 = quick_two_sum(c1, s0, c2);
    c0 = quick_two_sum(c0, s0, c1);
    s0 = c0;
    s1 = c1;
    if (s1 != 0.0) {
        s1 = quick_two_sum(s1, c2, s2);
        if (s2 != 0.0) {
            s2 = quick_two_sum(s2, c3, s3);
            if (s3 != 0.0)
                s3 += c4;
            else
                s2 = quick_two_sum(s2, c4, s3);
        } else {
            s1 = quick_two_sum(s1, c3, s2);
            if (s2 != 0.0)
                s2 = quick_two_sum(s2, c4, s3);
            else
                s1 = quick_two_sum(s1, c4, s2);
        }
    } else {
        s0 = quick_two_sum(s0, c2, s1);
        if (s1 != 0.0) {
            s1 = quick_two_sum(s1, c3, s2);
            if (s2 != 0.0)
                s2 = quick_two_sum(s2, c4, s3);
            else
                s1 = quick_two_sum(s1, c4, s2);
        } else {
            s0 = quick_two_sum(s0, c3, s1);
            if (s1 != 0.0)
                s1 = quick_two_sum(s1, c4, s2);
            else
                s0 = quick_two_sum(s0, c4, s1);
        }
    }
    c0 = s0;
    c1 = s1;
    c2 = s2;
}

// td_real::sloppy_add
MPLAPACK_GEMM_HD td3 td_sloppy_add(const td3 &a, const td3 &b) {
    double s0, s1, s2, t0, t1, t2;
    s0 = two_sum(a.x[0], b.x[0], t0);
    s1 = two_sum(a.x[1], b.x[1], t1);
    s2 = two_sum(a.x[2], b.x[2], t2);
    s1 = two_sum(s1, t0, t0);
    s2 += (t0 + t1);
    t0 = t2;
    td_renorm(s0, s1, s2, t0);
    return td_make(s0, s1, s2);
}

// operator*(td_real, td_real) of libQD3 built without QD_BF_MUL
MPLAPACK_GEMM_HD td3 td_sloppy_mul(const td3 &a, const td3 &b) {
    double p0, p1, p2, p3, p4, p5;
    double q0, q1, q2, q3, q4, q5;
    double t0, t1, s0, s1, s2;
    p0 = two_prod(a.x[0], b.x[0], q0);
    p1 = two_prod(a.x[0], b.x[1], q1);
    p2 = two_prod(a.x[1], b.x[0], q2);
    p3 = two_prod(a.x[0], b.x[2], q3);
    p4 = two_prod(a.x[1], b.x[1], q4);
    p5 = two_prod(a.x[2], b.x[0], q5);
    td_three_sum(p1, p2, q0);
    td_three_sum(p2, q1, q2);
    td_three_sum(p3, p4, p5);
    s0 = two_sum(p2, p3, t0);
    s1 = two_sum(q1, p4, t1);
    s2 = q2 + p5;
    s1 = two_sum(s1, t0, t0);
    s2 += t0 + t1;
    s1 += a.x[1] * b.x[2] + a.x[2] * b.x[1] + q0 + q3 + q4 + q5;
    td_renorm5(p0, p1, s0, s1, s2);
    return td_make(p0, p1, s0);
}

// qd_real::sloppy_add (the version reorganized to minimize data dependency)
MPLAPACK_GEMM_HD qd4 qd_sloppy_add(const qd4 &a, const qd4 &b) {
    double s0, s1, s2, s3, t0, t1, t2, t3;
    double v0, v1, v2, v3, u0, u1, u2, u3, w0, w1, w2, w3;
    s0 = a.x[0] + b.x[0];
    s1 = a.x[1] + b.x[1];
    s2 = a.x[2] + b.x[2];
    s3 = a.x[3] + b.x[3];
    v0 = s0 - a.x[0];
    v1 = s1 - a.x[1];
    v2 = s2 - a.x[2];
    v3 = s3 - a.x[3];
    u0 = s0 - v0;
    u1 = s1 - v1;
    u2 = s2 - v2;
    u3 = s3 - v3;
    w0 = a.x[0] - u0;
    w1 = a.x[1] - u1;
    w2 = a.x[2] - u2;
    w3 = a.x[3] - u3;
    u0 = b.x[0] - v0;
    u1 = b.x[1] - v1;
    u2 = b.x[2] - v2;
    u3 = b.x[3] - v3;
    t0 = w0 + u0;
    t1 = w1 + u1;
    t2 = w2 + u2;
    t3 = w3 + u3;
    s1 = two_sum(s1, t0, t0);
    qd_three_sum(s2, t0, t1);
    qd_three_sum2(s3, t0, t2);
    t0 = t0 + t1 + t3;
    renorm(s0, s1, s2, s3, t0);
    return qd_make(s0, s1, s2, s3);
}

// qd_real::sloppy_mul
MPLAPACK_GEMM_HD qd4 qd_sloppy_mul(const qd4 &a, const qd4 &b) {
    double p0, p1, p2, p3, p4, p5;
    double q0, q1, q2, q3, q4, q5;
    double t0, t1, s0, s1, s2;
    p0 = two_prod(a.x[0], b.x[0], q0);
    p1 = two_prod(a.x[0], b.x[1], q1);
    p2 = two_prod(a.x[1], b.x[0], q2);
    p3 = two_prod(a.x[0], b.x[2], q3);
    p4 = two_prod(a.x[1], b.x[1], q4);
    p5 = two_prod(a.x[2], b.x[0], q5);
    qd_three_sum(p1, p2, q0);
    qd_three_sum(p2, q1, q2);
    qd_three_sum(p3, p4, p5);
    s0 = two_sum(p2, p3, t0);
    s1 = two_sum(q1, p4, t1);
    s2 = q2 + p5;
    s1 = two_sum(s1, t0, t0);
    s2 += (t0 + t1);
    s1 += a.x[0] * b.x[3] + a.x[1] * b.x[2] + a.x[2] * b.x[1] + a.x[3] * b.x[0] + q0 + q3 + q4 + q5;
    renorm(p0, p1, s0, s1, s2);
    return qd_make(p0, p1, s0, s1);
}

} // namespace mplapack_gemm

#endif
