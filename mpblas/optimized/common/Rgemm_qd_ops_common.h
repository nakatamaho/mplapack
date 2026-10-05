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
 * Double-double and quad-double arithmetic on plain doubles, usable on the
 * host and in CUDA device code.  Every function performs the operations of
 * the corresponding libqd (QD3) inline function in the same order:
 *
 *   dd_add  dd_real::ieee_add        dd_mul  operator*(dd_real, dd_real)
 *   qd_add  qd_real::ieee_add        qd_mul  qd_real::accurate_mul
 *   qd_bf_add / qd_bf_mul            qd_real::bf_add / qd_real::bf_mul
 *                                    (Kouya, branch-free)
 *
 * so that the results are bit-identical to libqd configured with
 * QD_IEEE_ADD and without QD_SLOPPY_MUL / QD_BF_*, which is how MPLAPACK
 * builds it (external/qd).  The error-free product uses fma when libqd does
 * (QD_FMS defined in qd/qd_config.h), Dekker's splitting otherwise.  Code
 * using these functions must be compiled without floating-point contraction
 * (-ffp-contract=off, nvcc -fmad=false) and without value-unsafe
 * optimizations.
 */

#ifndef _MPLAPACK_RGEMM_QD_OPS_COMMON_H_
#define _MPLAPACK_RGEMM_QD_OPS_COMMON_H_

#include <cmath>
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

// The configurations of libqd these functions reproduce.
#if defined(QD_IEEE_ADD) && !defined(QD_SLOPPY_MUL) && !defined(QD_BF_ADD) && !defined(QD_BF_MUL) && !defined(QD_BF)
#define MPLAPACK_GEMM_QD_OPS_MATCH_LIBQD 1
#else
#define MPLAPACK_GEMM_QD_OPS_MATCH_LIBQD 0
#endif

namespace mplapack_gemm {

struct dd2 {
    double x[2];
};
struct qd4 {
    double x[4];
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

MPLAPACK_GEMM_HD dd2 dd_add(const dd2 &a, const dd2 &b) {
    double s1, s2, t1, t2;
    s1 = two_sum(a.x[0], b.x[0], s2);
    t1 = two_sum(a.x[1], b.x[1], t2);
    s2 += t1;
    s1 = quick_two_sum(s1, s2, s2);
    s2 += t2;
    dd2 r;
    r.x[0] = quick_two_sum(s1, s2, r.x[1]);
    return r;
}

MPLAPACK_GEMM_HD dd2 dd_mul(const dd2 &a, const dd2 &b) {
    double p1, p2;
    p1 = two_prod(a.x[0], b.x[0], p2);
    p2 += (a.x[0] * b.x[1] + a.x[1] * b.x[0]);
    dd2 r;
    r.x[0] = quick_two_sum(p1, p2, r.x[1]);
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

// The same on separate words (vectorizes better than the struct forms).
MPLAPACK_GEMM_HD void dd_mul_w(double a0, double a1, double b0, double b1, double &r0, double &r1) {
    double p1, p2;
    p1 = two_prod(a0, b0, p2);
    p2 += (a0 * b1 + a1 * b0);
    r0 = quick_two_sum(p1, p2, r1);
}

MPLAPACK_GEMM_HD void dd_add_w(double a0, double a1, double b0, double b1, double &r0, double &r1) {
    double s1, s2, t1, t2;
    s1 = two_sum(a0, b0, s2);
    t1 = two_sum(a1, b1, t2);
    s2 += t1;
    s1 = quick_two_sum(s1, s2, s2);
    s2 += t2;
    r0 = quick_two_sum(s1, s2, r1);
}

// ---- quad-double ----

MPLAPACK_GEMM_HD void renorm(double &c0, double &c1, double &c2, double &c3) {
    double s0, s1, s2 = 0.0, s3 = 0.0;
    if (MPLAPACK_GEMM_ISINF(c0))
        return;
    s0 = quick_two_sum(c2, c3, c3);
    s0 = quick_two_sum(c1, s0, c2);
    c0 = quick_two_sum(c0, s0, c1);
    s0 = c0;
    s1 = c1;
    if (s1 != 0.0) {
        s1 = quick_two_sum(s1, c2, s2);
        if (s2 != 0.0)
            s2 = quick_two_sum(s2, c3, s3);
        else
            s1 = quick_two_sum(s1, c3, s2);
    } else {
        s0 = quick_two_sum(s0, c2, s1);
        if (s1 != 0.0)
            s1 = quick_two_sum(s1, c3, s2);
        else
            s0 = quick_two_sum(s0, c3, s1);
    }
    c0 = s0;
    c1 = s1;
    c2 = s2;
    c3 = s3;
}

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

MPLAPACK_GEMM_HD void three_sum(double &a, double &b, double &c) {
    double t1, t2, t3;
    t1 = two_sum(a, b, t2);
    a = two_sum(c, t1, t3);
    b = two_sum(t2, t3, c);
}

MPLAPACK_GEMM_HD double quick_three_accum(double &a, double &b, double c) {
    double s;
    bool za, zb;
    s = two_sum(b, c, b);
    s = two_sum(a, s, a);
    za = (a != 0.0);
    zb = (b != 0.0);
    if (za && zb)
        return s;
    if (!zb) {
        b = a;
        a = s;
    } else {
        a = s;
    }
    return 0.0;
}

MPLAPACK_GEMM_HD qd4 qd_add(const qd4 &a, const qd4 &b) {
    int i, j, k;
    double s, t;
    double u, v; // double-length accumulator
    double x[4] = {0.0, 0.0, 0.0, 0.0};

    i = j = k = 0;
    if (MPLAPACK_GEMM_FABS(a.x[i]) > MPLAPACK_GEMM_FABS(b.x[j]))
        u = a.x[i++];
    else
        u = b.x[j++];
    if (MPLAPACK_GEMM_FABS(a.x[i]) > MPLAPACK_GEMM_FABS(b.x[j]))
        v = a.x[i++];
    else
        v = b.x[j++];

    u = quick_two_sum(u, v, v);

    while (k < 4) {
        if (i >= 4 && j >= 4) {
            x[k] = u;
            if (k < 3)
                x[++k] = v;
            break;
        }
        if (i >= 4)
            t = b.x[j++];
        else if (j >= 4)
            t = a.x[i++];
        else if (MPLAPACK_GEMM_FABS(a.x[i]) > MPLAPACK_GEMM_FABS(b.x[j])) {
            t = a.x[i++];
        } else
            t = b.x[j++];

        s = quick_three_accum(u, v, t);
        if (s != 0.0) {
            x[k++] = s;
        }
    }

    // add the rest
    for (k = i; k < 4; k++)
        x[3] += a.x[k];
    for (k = j; k < 4; k++)
        x[3] += b.x[k];

    renorm(x[0], x[1], x[2], x[3]);
    qd4 r;
    r.x[0] = x[0];
    r.x[1] = x[1];
    r.x[2] = x[2];
    r.x[3] = x[3];
    return r;
}

MPLAPACK_GEMM_HD qd4 qd_mul(const qd4 &qa, const qd4 &qb) {
    const double *a = qa.x, *b = qb.x;
    double p0, p1, p2, p3, p4, p5;
    double q0, q1, q2, q3, q4, q5;
    double p6, p7, p8, p9;
    double q6, q7, q8, q9;
    double r0, r1;
    double t0, t1;
    double s0, s1, s2;

    p0 = two_prod(a[0], b[0], q0);

    p1 = two_prod(a[0], b[1], q1);
    p2 = two_prod(a[1], b[0], q2);

    p3 = two_prod(a[0], b[2], q3);
    p4 = two_prod(a[1], b[1], q4);
    p5 = two_prod(a[2], b[0], q5);

    // Start Accumulation
    three_sum(p1, p2, q0);

    // Six-Three Sum of p2, q1, q2, p3, p4, p5.
    three_sum(p2, q1, q2);
    three_sum(p3, p4, p5);
    // compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5).
    s0 = two_sum(p2, p3, t0);
    s1 = two_sum(q1, p4, t1);
    s2 = q2 + p5;
    s1 = two_sum(s1, t0, t0);
    s2 += (t0 + t1);

    // O(eps^3) order terms
    p6 = two_prod(a[0], b[3], q6);
    p7 = two_prod(a[1], b[2], q7);
    p8 = two_prod(a[2], b[1], q8);
    p9 = two_prod(a[3], b[0], q9);

    // Nine-Two-Sum of q0, s1, q3, q4, q5, p6, p7, p8, p9.
    q0 = two_sum(q0, q3, q3);
    q4 = two_sum(q4, q5, q5);
    p6 = two_sum(p6, p7, p7);
    p8 = two_sum(p8, p9, p9);
    // Compute (t0, t1) = (q0, q3) + (q4, q5).
    t0 = two_sum(q0, q4, t1);
    t1 += (q3 + q5);
    // Compute (r0, r1) = (p6, p7) + (p8, p9).
    r0 = two_sum(p6, p8, r1);
    r1 += (p7 + p9);
    // Compute (q3, q4) = (t0, t1) + (r0, r1).
    q3 = two_sum(t0, r0, q4);
    q4 += (t1 + r1);
    // Compute (t0, t1) = (q3, q4) + s1.
    t0 = two_sum(q3, s1, t1);
    t1 += q4;

    // O(eps^4) terms -- Nine-One-Sum
    t1 += a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + q6 + q7 + q8 + q9 + s2;

    renorm(p0, p1, s0, t0, t1);
    qd4 r;
    r.x[0] = p0;
    r.x[1] = p1;
    r.x[2] = s0;
    r.x[3] = t0;
    return r;
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

} // namespace mplapack_gemm

#endif
